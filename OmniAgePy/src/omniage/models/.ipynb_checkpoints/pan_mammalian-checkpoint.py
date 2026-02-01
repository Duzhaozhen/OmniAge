import pandas as pd
import numpy as np
import os
from typing import Optional, Literal, List, Dict
from .base import BaseLinearClock

class PanMammalianClock:
    """
    Universal Pan-Mammalian Epigenetic Clocks (Lu et al., 2023).
    
    This class implements the multi-species epigenetic clocks applicable across 
    different mammalian tissues. It supports three specific variations:
    1. Universal: 3 clocks applicable to all tissues.
    2. Blood: 2 clocks optimized for blood tissue.
    3. Skin: 2 clocks optimized for skin tissue.
    
    These clocks require species-specific life history traits (e.g., gestation time, 
    max lifespan) which are loaded from the AnAge database included in the package.
    
    Reference:
    Lu, A.T. et al. Universal DNA methylation age across mammalian tissues. 
    Nature Aging (2023).
    """
    
    def __init__(self, tissue_type: Literal["Universal", "Blood", "Skin"] = "Universal", metadata: Dict = None,data_dir: str = None):
        """
        Initialize the Pan-Mammalian Clock.
        
        Args:
            tissue_type: The specific clock variant ('Universal', 'Blood', or 'Skin').
            data_dir: Path to the data directory. If None, defaults to package internal path.
        """
        self.tissue_type = tissue_type
        self.metadata = metadata or {}
        self.models = [] # Stores coefficient DataFrames for each sub-clock
        
        # 1. Automatically locate resource path
        if data_dir is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data")
        
        self.data_dir = data_dir
        
        # 2. Load AnAge Database (Life History Traits)
        anage_path = os.path.join(data_dir, "pan_mammalian", "anage_data.csv")
        if os.path.exists(anage_path):
            self.anage_df = pd.read_csv(anage_path)
        else:
            print(f"[Warning] AnAge database not found at {anage_path}. Predictions will fail.")
            self.anage_df = None

        # 3. Load Model Coefficients
        # Universal has 3 clocks; Blood/Skin have 2 clocks (technically Clock 2 & 3 variations)
        if tissue_type == "Universal":
            clock_indices = [1, 2, 3]
            prefix = "PanMammalian_Universal"
        elif tissue_type == "Blood":
            clock_indices = [1, 2] # Represents specific Blood clocks (analogous to Univ Clock 2 & 3)
            prefix = "PanMammalian_Blood"
        elif tissue_type == "Skin":
            clock_indices = [1, 2] # Represents specific Skin clocks
            prefix = "PanMammalian_Skin"
        else:
            raise ValueError("tissue_type must be 'Universal', 'Blood', or 'Skin'")

        for k in clock_indices:
            path = os.path.join(data_dir, "pan_mammalian", f"{prefix}_Clock{k}_coefs.csv")
            if os.path.exists(path):
                df = pd.read_csv(path)
                self.models.append(df)
            else:
                print(f"[Error] Model file not found: {path}")

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves coefficients for all sub-clocks (Clock 1, 2, 3) in this model.
        """
        if not self.models:
            return pd.DataFrame(columns=['model_name', 'probe', 'coef'])

        dfs = []
        
        # Re-derive indices logic matching __init__
        if self.tissue_type == "Universal":
            clock_indices = [1, 2, 3]
        else:
            clock_indices = [1, 2]

        prefix = f"PanMammalian_{self.tissue_type}"

        for i, df in enumerate(self.models):
            res_df = df.copy()
            
            # Construct name: e.g., PanMammalian_Universal_Clock1
            clock_num = clock_indices[i]
            model_name = f"{prefix}_Clock{clock_num}"
            
            res_df['model_name'] = model_name
            
            # Ensure column order
            if 'probe' in res_df.columns and 'coef' in res_df.columns:
                res_df = res_df[['model_name', 'probe', 'coef']]
            
            dfs.append(res_df)

        return pd.concat(dfs, ignore_index=True)
                
    def info(self):
        """Prints metadata information about the clock."""
        print(f"[Pan-Mammalian Clock: {self.tissue_type}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    # --- Mathematical Helper Functions (Ported from R) ---
    
    @staticmethod
    def _F2_antitrans_clock2(y, y_maxAge, y_gestation, const=1):
        """
        Inverse transformation for Clock 2.
        Converts relative age metric back to chronological years.
        """
        x0 = const * np.exp(-np.exp(-1 * y))
        x1 = x0 * (y_maxAge + y_gestation)
        x = x1 - y_gestation
        return x

    @staticmethod
    def _F1_logli(age1, m1, m2=None, c1=1):
        """Log-linear transformation function."""
        if m2 is None: m2 = m1
        
        condition = age1 >= m1
        res_true = (age1 - m1) / m2
        res_false = c1 * np.log((age1 - m1) / m2 / c1 + 1)
        
        return np.where(condition, res_true, res_false)

    @staticmethod
    def _F2_revtrsf_clock3(y_pred, m1, m2=None, c1=1):
        """
        Reverse transformation for Clock 3 (Relative Adult Age).
        """
        if m2 is None: m2 = m1
        
        condition = y_pred < 0
        res_true = (np.exp(y_pred / c1) - 1) * m2 * c1 + m1
        res_false = y_pred * m2 + m1
        
        return np.where(condition, res_true, res_false)

    @staticmethod
    def _F3_loglifn(dat1, b1=1, max_tage=4, c1=5, c2=0.38, c0=0):
        """
        Calculates intermediate life-history variables (a_Logli) required for transformation.
        Requires columns: 'GestationTimeInYears', 'averagedMaturity.yrs'.
        """
        df = dat1.copy()
        a2 = (df['GestationTimeInYears'] + c0) / df['averagedMaturity.yrs']
        a_Logli = c1 * (a2 ** c2)
        df['a_Logli'] = a_Logli
        return df

    def predict(self, beta_df: pd.DataFrame, sample_info: pd.DataFrame = None, verbose: bool = True) -> pd.DataFrame:
        """
        Execute prediction for the pan-mammalian clocks.
        
        Args:
            beta_df: Methylation matrix (Rows=CpGs, Columns=Samples).
            sample_info: DataFrame containing metadata. REQUIRED.
            verbose: If True, print warnings about missing species.
            
        Returns:
            pd.DataFrame: Predicted ages.
        """
        # --- 0. Enhanced Input Validation ---
        # Define a detailed error message template
        error_msg_template = """
[Input Error] sample_info is missing, invalid, or incomplete.

The Pan-Mammalian Clock requires species-specific life history traits (max lifespan, gestation, etc.) 
to calculate age. You must provide a 'sample_info' DataFrame to map samples to species.

REQUIRED FORMAT:
---------------------------------------------------------
Type: pandas.DataFrame
Columns:
  1. 'Sample'           : Must match column names in your beta_df.
  2. 'SpeciesLatinName' : Scientific name (e.g., 'Homo sapiens', 'Mus musculus').

EXAMPLE:
   Sample      SpeciesLatinName
0  GSM123456   Homo sapiens
1  GSM123457   Mus musculus
---------------------------------------------------------
Current Status: {status}
"""

        # Check 1: Is None
        if sample_info is None:
            raise ValueError(error_msg_template.format(status="sample_info is None."))

        # Check 2: Is DataFrame
        if not isinstance(sample_info, pd.DataFrame):
            raise TypeError(error_msg_template.format(status=f"sample_info is type '{type(sample_info).__name__}', expected 'pandas.DataFrame'."))

        # Check 3: Contains required columns
        req_cols = {'Sample', 'SpeciesLatinName'}
        missing_cols = req_cols - set(sample_info.columns)
        if missing_cols:
            raise ValueError(error_msg_template.format(status=f"Missing columns: {missing_cols}"))

        # --- Original Logic Follows ---

        if self.anage_df is None or not self.models:
            print("[Error] Missing dependencies (AnAge DB or Model Coefs). Returning empty.")
            return pd.DataFrame()

        # Ensure sample_info matches beta_df columns
        valid_samples = [s for s in sample_info['Sample'] if s in beta_df.columns]
        
        # Check 4: Sample matching check (Alert user if no matches found)
        if not valid_samples:
            raise ValueError(f"""
[Data Error] No matching samples found.
- Samples in sample_info : {len(sample_info)} (First 3: {sample_info['Sample'].head(3).tolist()})
- Samples in beta_df     : {len(beta_df.columns)} (First 3: {beta_df.columns[:3].tolist()})

Please check if 'Sample' column strictly matches the column names in beta_df.
""")
        
        # Filter metadata for valid samples
        info = sample_info[sample_info['Sample'].isin(valid_samples)].copy()
        
        # 2. Merge with AnAge Database
        info = pd.merge(info, self.anage_df, on='SpeciesLatinName', how='left')
        
        # Check for missing species data
        if info['maxAge'].isnull().any():
            missing_species = info[info['maxAge'].isnull()]['SpeciesLatinName'].unique()
            if verbose:
                print(f"[Warning] The following species were not found in AnAge DB and will be dropped: {missing_species}")
                print("Tip: Check spelling. Examples: 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus'.")
            info = info.dropna(subset=['maxAge'])
            
        if info.empty:
            print("[Warning] No valid samples remaining after species filtering.")
            return pd.DataFrame()

        # 3. Calculate HighmaxAge (Upper bound for age normalization)
        MYMAX = 1.3
        info['HighmaxAge'] = MYMAX * info['maxAge']
        
        # Specific overrides for Human and Mouse (as defined in original paper)
        mask_human = info['SpeciesLatinName'] == 'Homo sapiens'
        mask_mouse = info['SpeciesLatinName'] == 'Mus musculus'
        info.loc[mask_human, 'HighmaxAge'] = info.loc[mask_human, 'maxAge']
        info.loc[mask_mouse, 'HighmaxAge'] = info.loc[mask_mouse, 'maxAge']

        # 4. Prepare Result Container
        results = info[['Sample', 'SpeciesLatinName']].copy()

        # 5. Raw Prediction Loop
        raw_preds = [] 
        
        for k, model_df in enumerate(self.models):
            model_probes = model_df['probe'].values
            model_coefs = model_df['coef'].values
            
            # Handle Intercept
            intercept_mask = model_df['probe'] == 'Intercept'
            intercept_val = 0.0
            if intercept_mask.any():
                intercept_val = model_df.loc[intercept_mask, 'coef'].values[0]
                valid_mask = ~intercept_mask
                model_probes = model_probes[valid_mask]
                model_coefs = model_coefs[valid_mask]

            # Subset beta_df to relevant samples and transpose to (Samples x CpGs)
            X_subset = beta_df.loc[:, info['Sample']].T 
            
            # Reindex to match model probes, filling missing CpGs with 0
            X_aligned = X_subset.reindex(columns=model_probes, fill_value=0.0).fillna(0.0)
            
            # Calculate Linear Predictor (X * beta + intercept)
            y_pred = X_aligned.dot(model_coefs) + intercept_val
            raw_preds.append(y_pred.values)

        # 6. Inverse Transformation (Convert raw scores to Age)
        if self.tissue_type == "Universal":
            # --- Clock 1: Log-linear Age ---
            y1 = raw_preds[0]
            results['DNAmAgePanMammalianClock1'] = np.exp(y1) - 2
            
            # --- Clock 2: Relative Age ---
            y2 = raw_preds[1]
            results['DNAmRelativeAge'] = np.exp(-np.exp(-1 * y2))
            
            # Transform relative age to chronological years
            results['DNAmAgePanMammalianClock2'] = self._F2_antitrans_clock2(
                y2, info['HighmaxAge'], info['GestationTimeInYears']
            )
            
            # --- Clock 3: Relative Adult Age ---
            y3 = raw_preds[2]
            info_aug = self._F3_loglifn(info) 
            m1 = info_aug['a_Logli']
            
            results['DNAmRelativeAdultAge'] = self._F2_revtrsf_clock3(y3, m1)
            
            # Final Age 3 Calculation
            results['DNAmAgePanMammalianClock3'] = (
                results['DNAmRelativeAdultAge'] * (info['averagedMaturity.yrs'] + info['GestationTimeInYears']) - 
                info['GestationTimeInYears']
            )

        elif self.tissue_type in ["Blood", "Skin"]:
            suffix = f"_{self.tissue_type}"
            
            # --- Clock 2 Variation ---
            y2 = raw_preds[0]
            results[f'DNAmRelativeAge{suffix}'] = np.exp(-np.exp(-1 * y2))
            results[f'DNAmAgePanMammalian{self.tissue_type}2'] = self._F2_antitrans_clock2(
                y2, info['HighmaxAge'], info['GestationTimeInYears']
            )
            
            # --- Clock 3 Variation ---
            y3 = raw_preds[1]
            info_aug = self._F3_loglifn(info)
            m1 = info_aug['a_Logli']
            
            results[f'DNAmRelativeAdultAge{suffix}'] = self._F2_revtrsf_clock3(y3, m1)
            results[f'DNAmAgePanMammalian{self.tissue_type}3'] = (
                results[f'DNAmRelativeAdultAge{suffix}'] * (info['averagedMaturity.yrs'] + info['GestationTimeInYears']) - 
                info['GestationTimeInYears']
            )

        return results

# --- Convenience Wrappers for Specific Tissue Types ---

class PanMammalianUniversal(PanMammalianClock):
    """Universal Pan-Mammalian Clocks (3 variants)."""
    METADATA = {
        "year": 2023,
        "species": "Pan-Mammalian",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(Mammal40k)",
        "prediction": "Chronological Age & Relative Age",
        "source": "https://doi.org/10.1038/s43587-023-00462-6"
    }
    def __init__(self):
        super().__init__("Universal", metadata=self.METADATA)

class PanMammalianBlood(PanMammalianClock):
    """Blood-specific Pan-Mammalian Clocks (2 variants)."""
    METADATA = {
        "year": 2023,
        "species": "Pan-Mammalian",
        "tissue": "Blood",
        "omic type": "DNAm(Mammal40k)",
        "prediction": "Chronological Age & Relative Age",
        "source": "https://doi.org/10.1038/s43587-023-00462-6"
    }
    def __init__(self):
        super().__init__("Blood", metadata=self.METADATA)

class PanMammalianSkin(PanMammalianClock):
    """Skin-specific Pan-Mammalian Clocks (2 variants)."""
    METADATA = {
        "year": 2023,
        "species": "Pan-Mammalian",
        "tissue": "Skin",
        "omic type": "DNAm(Mammal40k)",
        "prediction": "Chronological Age & Relative Age",
        "source": "https://doi.org/10.1038/s43587-023-00462-6"
    }
    def __init__(self):
        super().__init__("Skin", metadata=self.METADATA)