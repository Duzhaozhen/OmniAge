import os
import json
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.interpolate import interp1d as _interp1d

class ProteomicOrganAge:
    """
    Implements the Proteomic Organ Age models to estimate biological age 
    and age gaps across 25 different organs using proteomic data.
    """

    def __init__(self, model_data_path=None):
        """
        Initialize the ProteomicOrganAge object.

        Args:
            model_data_path (str, optional): Path to the model configuration JSON. 
                If None, defaults to the internal package data directory.
        """
        if model_data_path is None:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            parent_dir = os.path.dirname(current_dir)
            model_data_path = os.path.join(parent_dir, "data", "proteomic_clock", "proteomic_organage_models.json")

        self.models = self._load_models(model_data_path)
        self.df_prot = None
        self.metadata = None
        self.results = None

    def _load_models(self, path):
        """Load model weights and scaling factors from a JSON file."""
        with open(path, 'r') as f:
            return json.load(f)

    def add_data(self, df_prot: pd.DataFrame, metadata: pd.DataFrame):
        """
        Inject and validate proteomic and clinical metadata.

        Args:
            df_prot: DataFrame where rows are samples and columns are protein probes (RFU).
            metadata: DataFrame containing sample info. Must include 'Age' and 'Sex'.

        Raises:
            ValueError: If required columns are missing or indices are misaligned.
        """
        # Validate required columns
        required_cols = ["Age", "Sex"]
        if not all(col in metadata.columns for col in required_cols):
            raise ValueError(f"Metadata is missing required columns: {required_cols}")

        # Strict sample alignment check (index must match in content and order)
        if not df_prot.index.equals(metadata.index):
            raise ValueError(
                "Sample alignment error: The index of 'df_prot' must exactly match "
                "the index of 'metadata' in both content and order."
            )

        self.df_prot = df_prot.copy()
        self.metadata = metadata.copy()
        
        # Binarize sex (Female = 1, Male = 0)
        self.metadata['Sex_F'] = np.where(self.metadata['Sex'] == 'Female', 1, 0)
        
        return self

    def get_coefs(self) -> pd.DataFrame:
        """
        Extract coefficients for all 25 organ models.

        Returns:
            pd.DataFrame: A long-format DataFrame with columns ['Organ', 'probe', 'coef'].
        """
        all_coefs = []
        
        for organ_name, cfg in self.models["organs"].items():
            # Extract protein weights
            df_weights = pd.DataFrame(cfg["lasso_weights"])
            df_weights = df_weights.rename(columns={"Feature": "probe", "Weight": "coef"})
            
            # Prepend intercept
            intercept_row = pd.DataFrame([{
                "probe": "Intercept", 
                "coef": cfg["lasso_intercept"]
            }])
            
            df_organ = pd.concat([intercept_row, df_weights], ignore_index=True)
            df_organ["Organ"] = organ_name
            all_coefs.append(df_organ)
            
        df_final = pd.concat(all_coefs, ignore_index=True)
        return df_final[["Organ", "probe", "coef"]]

    def predict(self, *args, **kwargs):
        """
        Dummy predict method for compatibility with model registration hooks.
        Actual age estimation is performed via .estimate_ages().
        """
        pass

    def estimate_ages(self, assay_version="v4.1", reference="cohort", verbose=True):
        """
        Predict organ-specific ages and calculate age gaps.

        Args:
            assay_version (str): SomaScan version ('v4', 'v4.1', or 'v5').
            reference (str): 'original' (uses fixed LOWESS from training) or 
                             'cohort' (re-calculates LOWESS on the current data).
            verbose (bool): If True, prints processing progress.

        Returns:
            dict: A dictionary of DataFrames, keyed by organ name.
        """
        if self.df_prot is None or self.metadata is None:
            raise ValueError("Data not found. Call .add_data() before estimation.")
            
        if assay_version not in ["v4.1", "v4", "v5"]:
            raise ValueError("Invalid assay_version. Choose from 'v4.1', 'v4', or 'v5'.")
        if reference not in ["original", "cohort"]:
            raise ValueError("Invalid reference. Choose 'original' or 'cohort'.")

        p_exp = self.df_prot.copy()
        m_dat = self.metadata

        # Step 1: Version Scaling
        if assay_version != "v4.1":
            scales = pd.DataFrame(self.models["scale_factors"][assay_version])
            common_prots = list(set(p_exp.columns) & set(scales["SeqId"]))
            
            if common_prots:
                scales_sub = scales.set_index("SeqId").loc[common_prots, "ScaleFactor"]
                p_exp[common_prots] = p_exp[common_prots] * scales_sub

        # Step 2: Log Transformation
        if p_exp.values.mean() < 500:
            warnings.warn("Data scale warning: Low global mean detected. Ensure input is raw RFU.", UserWarning)
        
        p_exp_log = np.log10(p_exp)

        # Step 3: Multi-organ Prediction Loop
        res_dict = {}
        for organ, cfg in self.models["organs"].items():
            if verbose: print(f"Processing organ: {organ}...")
            
            weights_df = pd.DataFrame(cfg["lasso_weights"])
            target_prots = weights_df[weights_df["Feature"] != "Sex_F"]["Feature"].tolist()
            
            # Check for missing protein features
            if not all(p in p_exp_log.columns for p in target_prots):
                if verbose: print(f"Skipped '{organ}': Features missing in input.")
                continue
                
            # Z-score normalization using JSON-defined parameters
            scaler_df = pd.DataFrame(cfg["prot_scaler"]).set_index("SeqId").loc[target_prots]
            z_prot = (p_exp_log[target_prots] - scaler_df["mean"].values) / scaler_df["sd"].values
            
            # Construct feature matrix [Sex_F, Proteins...]
            input_m = pd.concat([m_dat[["Sex_F"]], z_prot], axis=1)
            beta = weights_df.set_index("Feature").loc[input_m.columns, "Weight"].values
            
            # Linear combination: Intercept + (X * Beta)
            pred_age = cfg["lasso_intercept"] + np.dot(input_m.values, beta)
            
            # Step 4: Calculate Expected Age and Age Gaps
            if reference == "original":
                # Use pre-calculated LOWESS table from the training cohort
                lowess_df = pd.DataFrame(cfg["lowess_table"]).dropna()
                interp_func = _interp1d(lowess_df["Chronological_Age"], 
                                       lowess_df["Expected_Predicted_Age"], 
                                       kind='linear', bounds_error=False, fill_value=np.nan)
                exp_age = interp_func(m_dat["Age"])
                
                raw_gap = pred_age - exp_age
                z_gap = raw_gap / cfg["agegap_scaler"]["sd"]
                
            elif reference == "cohort":
                if len(m_dat) < 30:
                    raise ValueError(f"Cohort mode for '{organ}' requires N >= 30 for stable LOWESS fit.")
                
                # Fit LOWESS to the current cohort (equivalent to R's lowess)
                lowess_fit = sm.nonparametric.lowess(pred_age, m_dat["Age"], frac=2/3, it=5)
                lowess_func = _interp1d(lowess_fit[:, 0], lowess_fit[:, 1], 
                                       kind='linear', bounds_error=False, fill_value='extrapolate')
                exp_age = lowess_func(m_dat["Age"])
                
                raw_gap = pred_age - exp_age
                
                # Standardize gaps based on current cohort distribution
                if np.isnan(raw_gap).all():
                    z_gap = np.full(len(raw_gap), np.nan)
                else:
                    mu = np.nanmean(raw_gap)
                    sigma = np.nanstd(raw_gap, ddof=0) # Population standard deviation
                    z_gap = (raw_gap - mu) / sigma

            # Assemble results
            res_dict[organ] = pd.DataFrame({
                'SampleID': m_dat.index,
                'Organ': organ,
                'ChronologicalAge': m_dat['Age'],
                'PredictedAge': pred_age,
                'ExpectedAge': exp_age,
                'AgeGap': raw_gap,
                'AgeGapZ': z_gap
            })

        self.results = res_dict
        if verbose: print("Age estimation completed successfully.")
        return res_dict