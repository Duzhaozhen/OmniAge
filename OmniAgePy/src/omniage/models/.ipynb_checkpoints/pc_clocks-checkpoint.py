import pandas as pd
import numpy as np
import os
from pathlib import Path
from typing import Optional, Dict, Union
from ..utils import check_and_download_zenodo

PC_CLOCKS_URL = "https://zenodo.org/records/18287372/files/PCClocks.tar.gz?download=1"

class BasePCClock:
    """
    Base class for Principal Component (PC) based epigenetic clocks.
    
    This class handles data ingestion, resource loading, high-performance preprocessing,
    and PCA projection.
    
    [Performance Note]: 
    This implementation utilizes an "Intersection-Concat" strategy and enforces 
    `float32` precision to maximize throughput on large methylation datasets.

    Attributes:
        name (str): Name of the clock.
        metadata (dict): Metadata info.
        center (pd.Series): PCA centering vector.
        rotation (pd.DataFrame): PCA rotation matrix (CpGs x PCs).
        impute_ref (pd.Series): Reference means for missing value imputation.
    """
    
    def __init__(self, clock_name: str, metadata: Dict = None):
        """
        Initialize the PC Clock and load resources.

        Args:
            clock_name (str): Directory name of the clock in the data folder.
            metadata (Dict, optional): Metadata dictionary.
        """
        self.name = clock_name
        self.metadata = metadata or {} 
        current_file_dir = Path(__file__).parent.resolve()
        self.base_path = current_file_dir.parent / "data" / "PCClocks"
        self.clock_dir = self.base_path / clock_name

        check_and_download_zenodo(
            target_dir=self.base_path, 
            download_url=PC_CLOCKS_URL
        )
        
        if not self.clock_dir.exists():
            raise FileNotFoundError(f"Data for {clock_name} not found at {self.clock_dir}")

        self._load_assets()
    
    def info(self):
        """Prints summary information about the clock model."""
        print(f"[{self.name}] Model information:")
        if not self.metadata:
            print("No metadata available.")
            return
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")
            
    def _load_assets(self):
        """
        Loads necessary Parquet assets (centers, rotations, imputation means).
        """
        # Load Global Imputation Means
        impute_path = self.base_path / "PC_Impute_Means.parquet"
        self.impute_ref = pd.read_parquet(impute_path).set_index("probe")["mean"].astype(np.float32)
        
        # Load PCA Center and Rotation Matrix
        self.center = pd.read_parquet(self.clock_dir / "center.parquet").set_index("probe")["mean"].astype(np.float32)
        self.rotation = pd.read_parquet(self.clock_dir / "rotation.parquet").set_index("probe").astype(np.float32)

    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        [High-Performance Preprocessing] 
        Performs CpG alignment, transposition, and missing value imputation.
        
        Args:
            beta_df (pd.DataFrame): Raw Methylation Data.
                - **Rows**: CpGs (Probes)
                - **Columns**: Samples
                - **Values**: Beta values (0-1)

        Returns:
            pd.DataFrame: A **Sample-by-CpG** matrix aligned with the clock's rotation matrix.
            (Note: The dimensions are transposed compared to input).
        """
        # 1. Identify Required vs. Missing CpGs
        required_cpgs = self.rotation.index
        
        # Fast Intersection: Filter indices before data extraction to reduce memory overhead
        common_cpgs = beta_df.index.intersection(required_cpgs)
        missing_cpgs = required_cpgs.difference(common_cpgs)
        
        # 2. Extract Existing Data & Transpose
        # Input (CpG x Sample) -> Output (Sample x CpG)
        X_existing = beta_df.loc[common_cpgs].T.astype(np.float32)
        
        # ============================================================
        # 3. Layer 1 Imputation: Handle Sporadic NAs in Existing Data
        if X_existing.isna().values.any():
            nan_cols = X_existing.columns[X_existing.isna().any()]
            col_means = X_existing[nan_cols].mean(axis=0, skipna=True)
            col_means = col_means.fillna(self.impute_ref)
            X_existing[nan_cols] = X_existing[nan_cols].fillna(col_means)

        # 4. Layer 2 Imputation: Handle Completely Missing Columns
        if len(missing_cpgs) > 0:
            fill_vals = self.impute_ref.loc[missing_cpgs].values.astype(np.float32)
            n_samples = X_existing.shape[0]
            missing_data = np.tile(fill_vals, (n_samples, 1))
            
            X_missing = pd.DataFrame(
                missing_data, 
                index=X_existing.index, 
                columns=missing_cpgs
            )
            X = pd.concat([X_existing, X_missing], axis=1)
        else:
            X = X_existing

        # 5. Final Reorder & Safety Fallback
        X = X[required_cpgs]
        X = X.fillna(0.0)
            
        return X

    def get_pcs(self, X: pd.DataFrame) -> pd.DataFrame:
        """Calculates Principal Components (PCs) using vectorized linear algebra."""
        center_vec = self.center.loc[X.columns].values
        X_values = X.values 
        X_centered = X_values - center_vec
        rotation_values = self.rotation.values
        pcs = X_centered @ rotation_values
        return pd.DataFrame(pcs, index=X.index, columns=self.rotation.columns)


# --- 3. Standard PC Clock Implementation ---
class StandardPCClock(BasePCClock):
    """
    Implementation for standard linear PC clocks (e.g., Horvath, Hannum, PhenoAge).
    Calculates age based on a linear combination of PCs + Intercept.
    """
    def __init__(self, clock_name: str, do_anti_trafo: bool, metadata: Dict = None):
        super().__init__(clock_name, metadata=metadata)
        self.do_anti_trafo = do_anti_trafo
        
        model_df = pd.read_parquet(self.clock_dir / "model.parquet")
        self.coefs = pd.Series(model_df['coef'].values.astype(np.float32), index=model_df['pc'])
        
        with open(self.clock_dir / "intercept.txt", "r") as f:
            self.intercept = float(f.read().strip())

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the list of CpGs required by this PC clock.
        
        Returns:
            pd.DataFrame: Columns ['model_name', 'probe', 'coef'].
            Note: 'coef' is NaN because weights are applied to PCs, not CpGs directly.
        """
        required_cpgs = self.rotation.index
        df = pd.DataFrame({
            'probe': required_cpgs,
            'coef': np.nan 
        })
        df['model_name'] = self.name
        return df[['model_name', 'probe', 'coef']]    

    def anti_trafo(self, x: np.ndarray, adult_age: float = 20) -> np.ndarray:
        """Applies inverse log-linear transformation (Horvath-style)."""
        return np.where(x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age)

    def predict(self, beta_df: pd.DataFrame, verbose: bool = False, **kwargs) -> pd.Series:
        """
        Predict Age using the PC Clock.

        Args:
            beta_df (pd.DataFrame): Input methylation data (Rows: CpGs, Columns: Samples).
            verbose (bool): Unused, kept for compatibility.
            **kwargs: Unused arguments (e.g. 'ages', 'sex') ignored by standard clocks.

        Returns:
            pd.Series: Predicted ages, indexed by sample ID.
        """
        X = self.preprocess(beta_df)
        pcs_df = self.get_pcs(X)
        common_pcs = self.coefs.index.intersection(pcs_df.columns)
        raw_score = pcs_df[common_pcs] @ self.coefs[common_pcs] + self.intercept
        
        if self.do_anti_trafo:
            return pd.Series(self.anti_trafo(raw_score.values), index=X.index)
        return raw_score


# --- 4. GrimAge Implementation (Optimized) ---
class PCGrimAge1Impl(BasePCClock):
    """
    Implementation of the PCGrimAge clock logic.
    Predicts composite age/mortality risk based on PC-estimated surrogates.
    """
    def __init__(self, metadata: Dict = None):
        super().__init__("PCGrimAge1", metadata=metadata)
        
        surr_df = pd.read_parquet(self.clock_dir / "surrogate_weights.parquet")
        self.surr_intercepts = surr_df[surr_df['term'] == 'Intercept'].set_index('target')['coef'].astype(np.float32)
        self.surr_weights = surr_df[surr_df['term'] != 'Intercept'].copy()
        self.surr_weights['coef'] = self.surr_weights['coef'].astype(np.float32)
        
        final_df = pd.read_parquet(self.clock_dir / "final_model.parquet")
        self.final_intercept = final_df[final_df['term'] == 'Intercept']['coef'].iloc[0]
        self.final_weights = final_df[final_df['term'] != 'Intercept'].set_index('term')['coef'].astype(np.float32)
        
    def get_coefs(self) -> pd.DataFrame:
        """Retrieves required features (CpGs + Metadata) for PCGrimAge."""
        required_cpgs = self.rotation.index
        df_cpgs = pd.DataFrame({'probe': required_cpgs, 'coef': np.nan})
        df_meta = pd.DataFrame([{'probe': 'Age', 'coef': np.nan}, {'probe': 'Female', 'coef': np.nan}])
        final_df = pd.concat([df_meta, df_cpgs], ignore_index=True)
        final_df['model_name'] = self.name
        return final_df[['model_name', 'probe', 'coef']]
        
    def predict(self, beta_df: pd.DataFrame, ages: pd.Series, sex: pd.Series, verbose: bool = False, **kwargs) -> pd.DataFrame:
        """
        Predict PCGrimAge and its surrogate biomarkers.

        Args:
            beta_df (pd.DataFrame): Methylation data (Rows: CpGs, Columns: Samples).
            ages (pd.Series): Chronological age for each sample.
            sex (pd.Series): Sex for each sample. 
                Supported values (case-insensitive):
                - Female: 'F', 'Female', 'Woman', 'W'
                - Male: 'M', 'Male', 'Man'
            verbose (bool): Unused.

        Returns:
            pd.DataFrame: A DataFrame containing:
                - Surrogate Scores (e.g., DNAmPACKYRS, DNAmADM...)
                - **PCGrimAge1**: The final predicted biological age.
        """
        if ages is None or sex is None:
            raise ValueError("[PCGrimAge1] Missing required 'ages' or 'sex' arguments.")

        # --- Step 1: Preprocessing & PCA ---
        X = self.preprocess(beta_df)
        pcs_df = self.get_pcs(X)
        samples = X.index
        
        # --- Step 2: Metadata Alignment ---
        ages = ages.reindex(samples).astype(float)
        sex = sex.reindex(samples)

        # --- Step 3: Sex Normalization ---
        sex_cleaned = sex.astype(str).str.strip().str.lower()
        sex_map_to_int = {
            'f': 1, 'female': 1, 'woman': 1, 'w': 1, 
            'm': 0, 'male': 0, 'man': 0,
        }
        is_female = sex_cleaned.map(sex_map_to_int)
        
        if is_female.isna().any():
            invalid_entries = sex[is_female.isna()].unique()
            raise ValueError(f"[PCGrimAge1] Invalid sex values: {list(invalid_entries)}. Expected: {list(sex_map_to_int.keys())}")

        # --- Step 4: Surrogate Calculation ---
        pcs_df['Age'] = ages
        pcs_df['Female'] = is_female
        
        targets = self.surr_weights['target'].unique()
        surrogate_scores = pd.DataFrame(index=samples, columns=targets)
        
        for target in targets:
            w_df = self.surr_weights[self.surr_weights['target'] == target]
            w_series = pd.Series(w_df['coef'].values, index=w_df['term'])
            intercept = self.surr_intercepts.get(target, 0.0)
            
            common_features = w_series.index.intersection(pcs_df.columns)
            val = pcs_df[common_features] @ w_series[common_features] + intercept
            surrogate_scores[target] = val
            
        # --- Step 5: Final Prediction ---
        rename_map = {
            'PCPACKYRS': 'DNAmPACKYRS', 'PCADM': 'DNAmADM', 'PCB2M': 'DNAmB2M',
            'PCCystatinC': 'DNAmCystatinC', 'PCGDF15': 'DNAmGDF15', 'PCLeptin': 'DNAmLeptin',
            'PCPAI1': 'DNAmPAI1', 'PCTIMP1': 'DNAmTIMP1'
        }
        final_features = surrogate_scores.rename(columns=rename_map)
        final_features['Age'] = ages
        final_features['Female'] = is_female
        
        aligned_features = final_features.reindex(columns=self.final_weights.index, fill_value=0)
        pred_age = aligned_features @ self.final_weights + self.final_intercept
        
        surrogate_scores[self.name] = pred_age
        
        return surrogate_scores


# --- 5. Export Wrapper Classes ---

class PCHorvath2013(StandardPCClock):
    """
    PC-Horvath (2013) Clock.
    
    A multi-tissue predictor of chronological age, adapted to use Principal Components
    to improve reliability across different datasets.
    """
    METADATA = {
        "year": 2022, "species": "Human", "tissue": "Multi-tissue",
        "omic type": "DNAm(450k)", "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1038/s43587-022-00248-2", "original_clock": "Horvath2013"
    }
    def __init__(self):
        """Initialize PC-Horvath2013 (No arguments required)."""
        super().__init__("PCHorvath2013", do_anti_trafo=True, metadata=self.METADATA)

class PCHorvath2018(StandardPCClock):
    """
    PC-Horvath (2018) Skin & Blood Clock.
    
    Optimized for Skin and Blood tissues, utilizing PC projections for robustness.
    """
    METADATA = {
        "year": 2022, "species": "Human", "tissue": "Skin/Blood",
        "omic type": "DNAm(450k)", "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1038/s43587-022-00248-2", "original_clock": "Horvath2018"
    }
    def __init__(self):
        """Initialize PC-Horvath2018 (No arguments required)."""
        super().__init__("PCHorvath2018", do_anti_trafo=True, metadata=self.METADATA)

class PCHannum(StandardPCClock):
    """
    PC-Hannum Clock.
    
    A blood-based chronological age predictor.
    """
    METADATA = {
        "year": 2022, "species": "Human", "tissue": "Blood",
        "omic type": "DNAm(450k)", "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1038/s43587-022-00248-2", "original_clock": "Hannum2013"
    }
    def __init__(self):
        """Initialize PC-Hannum (No arguments required)."""
        super().__init__("PCHannum", do_anti_trafo=False, metadata=self.METADATA)

class PCPhenoAge(StandardPCClock):
    """
    PC-PhenoAge.
    
    Predicts 'Phenotypic Age', a marker of mortality risk and morbidity.
    """
    METADATA = {
        "year": 2022, "species": "Human", "tissue": "Blood",
        "omic type": "DNAm(450k)", "prediction": "Mortality",
        "source": "https://doi.org/10.1038/s43587-022-00248-2", "original_clock": "PhenoAge"
    }
    def __init__(self):
        """Initialize PC-PhenoAge (No arguments required)."""
        super().__init__("PCPhenoAge", do_anti_trafo=False, metadata=self.METADATA)

class PCDNAmTL(StandardPCClock):
    """
    PC-DNAmTL (Telomere Length).
    
    Predicts telomere length based on DNA methylation patterns.
    """
    METADATA = {
        "year": 2022, "species": "Human", "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)", "prediction": "Telomere Length",
        "source": "https://doi.org/10.1038/s43587-022-00248-2", "original_clock": "DNAmTL"
    }
    def __init__(self):
        """Initialize PC-DNAmTL (No arguments required)."""
        super().__init__("PCDNAmTL", do_anti_trafo=False, metadata=self.METADATA)

class PCGrimAge1(PCGrimAge1Impl):
    """
    PC-GrimAge (Version 1).
    
    A robust predictor of lifespan and healthspan. 
    Requires Chronological Age and Sex as additional inputs.
    """
    METADATA = {
        "year": 2022, "species": "Human", "tissue": "Blood",
        "omic type": "DNAm(450k)", "prediction": "Mortality",
        "source": "https://doi.org/10.1038/s43587-022-00248-2", "original_clock": "GrimAge1"
    }
    def __init__(self):
        """Initialize PC-GrimAge1 (No arguments required)."""
        super().__init__(metadata=self.METADATA)