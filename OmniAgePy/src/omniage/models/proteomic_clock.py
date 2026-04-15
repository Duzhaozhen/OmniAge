import os
import json
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.interpolate import interp1d

class WyssCorayOrganAge:
    """
    Implementation of the Wyss-Coray Organ-Specific Aging Clocks (Nature, 2023).
    
    This suite uses plasma proteomics (SomaScan) to estimate biological age 
    across 11 major organs and systems.
    """

    METADATA = {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "omic_type": "Proteomics (SomaScanV4/V4.1/V5)",
        "prediction": "Biological Age (Years)",
        "source": "https://doi.org/10.1038/s41586-023-06802-1"
    }
    
    def __init__(self, model_data_path=None):
        """
        Initializes the ProteomicOrganAge suite.
        
        :param model_data_path: Path to the JSON model file. If None, it automatically 
                                resolves to the default packaged JSON.
        """
        if model_data_path is None:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            parent_dir = os.path.dirname(current_dir)
            model_data_path = os.path.join(parent_dir, "data", "proteomic_clock", "proteomic_wysscoray_organage_models.json")
            
        self.models = self._load_models(model_data_path)

        self.name = "WyssCorayOrganAge"
        self.metadata = self.METADATA
        self.models = self._load_models(model_data_path)
    def info(self):
        """
        Print formatted model metadata to the console.
        Consistent with the omniage model API.
        """
        print(f"[{self.name}] Model Information:")
        print("-" * 40)
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"{label:12}: {value}")
        print("-" * 40)
    

    def _load_models(self, path):
        with open(path, 'r') as f:
            return json.load(f)

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the LASSO weights for all 25 organ models and consolidates them into a single DataFrame.
        
        :return: A DataFrame formatted for downstream frameworks: ['Organ', 'probe', 'coef']
        """
        all_coefs = []
        for organ_name, cfg in self.models["organs"].items():
            df_weights = pd.DataFrame(cfg["lasso_weights"])
            df_weights = df_weights.rename(columns={"Feature": "probe", "Weight": "coef"})
            
            intercept_row = pd.DataFrame([{
                "probe": "Intercept", 
                "coef": cfg["lasso_intercept"]
            }])
            
            df_organ = pd.concat([intercept_row, df_weights], ignore_index=True)
            df_organ["Organ"] = organ_name
            all_coefs.append(df_organ)
            
        df_final = pd.concat(all_coefs, ignore_index=True)
        return df_final[["Organ", "probe", "coef"]]

    def predict(self, prot_df, metadata, assay_version="v4.1", reference="cohort", verbose=True):
        """
        Core prediction engine: Executes data ingestion, validation, and multi-organ age estimation in a single step.
        
        :param prot_df: Protein expression matrix. Supports both standard ML format (Samples x Proteins) 
                        and Bioconductor format (Proteins x Samples).
        :param metadata: Sample metadata. Must include a 'Sex' or 'Sex_F' column. 'Age' is optional but required for AgeGap calculations.
        :param assay_version: SomaScan assay version ("v4.1", "v4", or "v5"). Default is "v4.1".
        :param reference: Reference baseline for Age Gap ("original" or "cohort"). Default is "cohort".
        :param verbose: Boolean indicating whether to print progress logs.
        """
        p_exp = prot_df.copy()
        m_dat = metadata.copy()

        # ==========================================
        # 1. Intelligent Dimension Detection & Auto-Transpose
        # ==========================================
        if p_exp.index.equals(m_dat.index):
            # Mode A: ML Standard (Rows = Samples, Columns = Features)
            pass 
        elif p_exp.columns.equals(m_dat.index):
            # Mode B: Bioconductor Standard (Rows = Features, Columns = Samples)
            if verbose: 
                print("💡 Auto-Detection: Detected Bioconductor-style input (Rows=Proteins, Columns=Samples). Matrix has been automatically transposed.")
            p_exp = p_exp.T
        else:
            raise ValueError(
                "Sample alignment error: Neither the rows nor the columns of 'prot_df' "
                "perfectly match the indices of 'metadata'. Please ensure exact matching of Sample IDs."
            )

        # ==========================================
        # 2. Metadata Validation & Pre-processing
        # ==========================================
        if "Sex_F" not in m_dat.columns:
            if "Sex" in m_dat.columns:
                m_dat['Sex_F'] = np.where(m_dat['Sex'] == 'Female', 1, 0)
            else:
                raise ValueError("Metadata must contain either a 'Sex' (Female/Male) or 'Sex_F' (1/0) column.")

        has_age = "Age" in m_dat.columns
        if not has_age and verbose:
            print("Note: 'Age' column not found in metadata. Only PredictedAge will be calculated. AgeGap metrics will be set to NA.")

        if assay_version not in ["v4.1", "v4", "v5"]:
            raise ValueError("Invalid assay_version: must be 'v4.1', 'v4', or 'v5'.")
        if reference not in ["original", "cohort"]:
            raise ValueError("Invalid reference baseline: must be 'original' or 'cohort'.")

        # ==========================================
        # 3. Protein Scale Adjustment & Log10 Transformation
        # ==========================================
        if assay_version != "v4.1":
            scales = pd.DataFrame(self.models["scale_factors"][assay_version])
            common_prots = list(set(p_exp.columns) & set(scales["SeqId"]))
            if common_prots:
                scales_sub = scales.set_index("SeqId").loc[common_prots, "ScaleFactor"]
                p_exp[common_prots] = p_exp[common_prots] * scales_sub

        if p_exp.values.mean() < 500:
            warnings.warn("Data scale warning: The global mean of the input protein matrix is low. Ensure data is provided as raw RFU.", UserWarning)
        
        p_exp_log = np.log10(p_exp)

        # ==========================================
        # 4. Multi-organ Prediction Loop
        # ==========================================
        res_dict = {}
        for organ, cfg in self.models["organs"].items():
            if verbose: print(f"Processing organ model: {organ}...")
            
            weights_df = pd.DataFrame(cfg["lasso_weights"])
            target_prots = weights_df[weights_df["Feature"] != "Sex_F"]["Feature"].tolist()
            
            if not all(p in p_exp_log.columns for p in target_prots):
                if verbose: print(f"  -> Skipped '{organ}': Insufficient protein features found in input data.")
                continue
                
            scaler_df = pd.DataFrame(cfg["prot_scaler"]).set_index("SeqId").loc[target_prots]
            z_prot = (p_exp_log[target_prots] - scaler_df["mean"].values) / scaler_df["sd"].values
            
            input_m = pd.concat([m_dat[["Sex_F"]], z_prot], axis=1)
            beta = weights_df.set_index("Feature").loc[input_m.columns, "Weight"].values
            pred_age = cfg["lasso_intercept"] + np.dot(input_m.values, beta)
            
            # ==========================================
            # 5. Expected Age and Age Gap Calculation
            # ==========================================
            if has_age:
                if reference == "original":
                    lowess_df = pd.DataFrame(cfg["lowess_table"]).dropna()
                    interp_func = interp1d(lowess_df["Chronological_Age"], 
                                           lowess_df["Expected_Predicted_Age"], 
                                           kind='linear', bounds_error=False, fill_value=np.nan)
                    exp_age = interp_func(m_dat["Age"])
                    
                    raw_gap = pred_age - exp_age
                    z_gap = raw_gap / cfg["agegap_scaler"]["sd"][0]
                    
                elif reference == "cohort":
                    if len(m_dat) < 30:
                        raise ValueError(f"Cohort normalization mode for '{organ}' requires a minimum of 30 samples with valid Age data.")
                    
                    lowess_fit = sm.nonparametric.lowess(pred_age, m_dat["Age"], frac=2/3, it=5)
                    lowess_func = interp1d(lowess_fit[:, 0], lowess_fit[:, 1], 
                                           kind='linear', bounds_error=False, fill_value='extrapolate')
                    exp_age = lowess_func(m_dat["Age"])
                    raw_gap = pred_age - exp_age
                    
                    if np.isnan(raw_gap).all():
                        z_gap = np.nan
                    else:
                        mu = np.nanmean(raw_gap)
                        sigma = np.nanstd(raw_gap, ddof=0)
                        z_gap = (raw_gap - mu) / sigma
            else:
                exp_age = np.nan
                raw_gap = np.nan
                z_gap = np.nan

            res_dict[organ] = pd.DataFrame({
                'SampleID': m_dat.index,
                'Organ': organ,
                'ChronologicalAge': m_dat['Age'] if has_age else np.nan,
                'PredictedAge': pred_age,
                'ExpectedAge': exp_age,
                'AgeGap': raw_gap,
                'AgeGapZ': z_gap
            })

        if verbose: print("✅ Age estimation completed successfully.")
        return res_dict



class GladyshevOrganAge:
    """
    Implementation of the Gladyshev Organ-Specific Biological Age and Mortality Clocks.
    
    This suite provides models for predicting biological age (Gen1) and mortality risk 
    (Gen2) across multiple organs using proteomic data, specifically optimized 
    for Olink platforms.
    """

    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "Proteomics(Olink 3072/1536)",
        "prediction": "Age(Years) & Mortality Risk",
        "source": "https://doi.org/10.1016/j.cmet.2024.10.005"
    }
   

    def __init__(self, model_data_path=None):
        """
        Initialize the GladyshevOrganAge suite.

        Args:
            model_data_path (str, optional): Path to the JSON model file. 
                If None, defaults to the internal package data directory.
        """
        if model_data_path is None:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            parent_dir = os.path.dirname(current_dir)
            # Default path for the converted JSON model weights
            model_data_path = os.path.join(
                parent_dir, "data", "proteomic_clock", "proteomic_gladyshev_organage_models.json"
            )
            
        self.models = self._load_models(model_data_path)

        self.name = "GladyshevOrganAge"
        self.metadata = self.METADATA
        self.models = self._load_models(model_data_path)

    def info(self):
        """
        Print formatted model metadata to the console.
        Consistent with other clocks in the omniage package.
        """
        print(f"[{self.name}] Model Metadata:")
        print("-" * 40) 
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').title() 
            print(f"{label:12}: {value}")
        print("-" * 40)

    def get_coefs(self) -> pd.DataFrame:
        """
        Extract all model weights for the Gladyshev Organ clocks.

        Flattens the nested model dictionary into a comprehensive long-format DataFrame.

        Returns:
            pd.DataFrame: A DataFrame with columns 
                ['Platform', 'ModelType', 'Organ', 'probe', 'fold_1', ..., 'fold_5'].
        """
        all_records = []
        
        # Iterate through Platforms (full / reduced)
        for platform, type_dict in self.models["coef"].items():
            # Iterate through Model Types (gen1: biological age / gen2: mortality)
            for model_type, organ_dict in type_dict.items():
                # Iterate through Organs (e.g., Adipose, Heart, Lung)
                for organ, weights_list in organ_dict.items():
                    # weights_list is a list of dicts: [{'feature_name': '...', 'fold_1': ...}]
                    for row in weights_list:
                        # Use 'probe' to maintain compatibility with other package clocks
                        feature = row.get("feature_name", "Intercept")
                        
                        record = {
                            "Platform": platform,
                            "ModelType": model_type,
                            "Organ": organ,
                            "probe": feature
                        }
                        
                        # Extract weights for all available folds
                        for k, v in row.items():
                            if k != "feature_name":
                                record[k] = v
                                
                        all_records.append(record)
                        
        df_final = pd.DataFrame(all_records)
        return df_final

    def _load_models(self, path):
        """Internal: Load model parameters from a JSON file."""
        with open(path, 'r') as f:
            return json.load(f)

    def predict(self, prot_df, model_type="gen2", platform="reduced", 
                fold=1, to_years=False, standardize=True, verbose=True):
        """
        Predict organ-specific biological age or mortality risk.

        Args:
            prot_df (pd.DataFrame): Protein expression matrix. Supports ML format 
                (Samples x Proteins) or Bioconductor format (Proteins x Samples).
            model_type (str): "gen1" (biological age) or "gen2" (mortality risk). 
                Defaults to "gen2".
            platform (str): "full" (Olink 3072) or "reduced" (Olink 1536). 
                Defaults to "reduced".
            fold (Union[int, str]): Integer (1-5) or "ensemble" to average results. 
                Defaults to 1.
            to_years (bool): If True, converts Gen2 relative hazard to years.
            standardize (bool): Align variance to the UK Biobank (UKB) cohort. 
                Highly recommended for Olink data.
            verbose (bool): If True, prints progress and auto-detection logs.

        Returns:
            pd.DataFrame: Predicted scores/ages for each organ per sample.
        """
        if model_type not in ["gen1", "gen2"]:
            raise ValueError("Invalid model_type: must be 'gen1' or 'gen2'.")
        if platform not in ["full", "reduced"]:
            raise ValueError("Invalid platform: must be 'full' or 'reduced'.")

        p_exp = prot_df.copy()
        
        # Extract model-specific parameters
        ukb_sds = self.models["sds"][platform]
        coef_dict = self.models["coef"][platform][model_type]
        
        # 1. Intelligent Dimension Detection & Auto-Transpose
        # Detect orientation based on feature name overlap with model features
        col_overlap = len(set(p_exp.columns) & set(ukb_sds.keys()))
        idx_overlap = len(set(p_exp.index) & set(ukb_sds.keys()))
        
        if idx_overlap > col_overlap and idx_overlap > 50:
            if verbose: 
                print("[Info] Detected Bioconductor-style input (Proteins x Samples). "
                      "Automatically transposing matrix.")
            p_exp = p_exp.T

        # 2. Feature Intersection
        common_prots = list(set(p_exp.columns) & set(ukb_sds.keys()))
        if not common_prots:
            raise ValueError("Feature mismatch: No proteins from input data found in model features.")
            
        p_exp = p_exp[common_prots]

        # 3. Variance Alignment (UK Biobank Rescaling)
        if standardize:
            if verbose: print("[Action] Aligning variance with UK Biobank reference cohort...")
            # Use ddof=1 to align with R's default sd() behavior
            user_sds = p_exp.std(axis=0, ddof=1)
            user_sds[user_sds == 0] = 1.0  # Avoid division by zero
            
            # Rescale user values: (Value / User_SD) * UKB_SD
            ukb_sds_series = pd.Series(ukb_sds)[common_prots]
            p_exp = (p_exp / user_sds) * ukb_sds_series

        # 4. Multi-organ Prediction Loop
        res_df = pd.DataFrame(index=p_exp.index)
        
        for organ, df_weights_list in coef_dict.items():
            weights_df = pd.DataFrame(df_weights_list)
            
            # Select weights based on fold or ensemble average
            if str(fold).lower() == "ensemble":
                fold_cols = [c for c in weights_df.columns if c.startswith("fold_")]
                final_coefs = weights_df.set_index("feature_name")[fold_cols].mean(axis=1)
            else:
                col_name = f"fold_{fold}"
                if col_name not in weights_df.columns:
                    raise ValueError(f"Fold '{col_name}' not found for organ '{organ}'.")
                final_coefs = weights_df.set_index("feature_name")[col_name]
                
            # Final protein alignment for the specific organ model
            target_prots = list(set(final_coefs.index) & set(common_prots))
            matched_weights = final_coefs[target_prots]
            
            # Dot product for score calculation
            score = p_exp[target_prots].dot(matched_weights)
            
            # Handle Intercept for Gen1 (Chronological Age) models
            if model_type == "gen1" and "Intercept" in final_coefs.index:
                score += final_coefs["Intercept"]
                
            # Apply Gompertz transform for Gen2 (Mortality) models if requested
            if model_type == "gen2" and to_years:
                gompertz = self.models.get("mortality_transform")
                if not gompertz or "slope" not in gompertz:
                    raise ValueError("Gompertz parameters missing in model metadata.")
                
                # Formula: (RelativeHazard - H0) / Slope - Intercept
                h0 = gompertz["avg_rel_log_mort_hazard"]
                score = (-h0 + score) / gompertz["slope"] - gompertz["intercept"]
                
            res_df[organ] = score
            
        if verbose: print("✅ Gladyshev age/risk estimation completed successfully.")
        return res_df


