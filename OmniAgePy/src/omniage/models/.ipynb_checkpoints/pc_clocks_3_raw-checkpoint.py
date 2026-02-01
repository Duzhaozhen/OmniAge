import pandas as pd
import numpy as np
import os
from pathlib import Path
from typing import Optional, Dict

class BasePCClock:
    """
    Base class for Principal Component (PC) based epigenetic clocks.
    
    This class handles data ingestion, resource loading, high-performance preprocessing,
    and PCA projection.
    
    [Performance Note]: 
    This implementation utilizes an "Intersection-Concat" strategy (mirroring efficient R logic)
    and enforces `float32` precision to maximize throughput and minimize memory usage 
    on large methylation datasets (e.g., EPIC array).
    """
    def __init__(self, clock_name: str, metadata: Dict = None):
        self.name = clock_name
        self.metadata = metadata or {} 
        current_file_dir = Path(__file__).parent.resolve()
        self.base_path = current_file_dir.parent / "data" / "PCClocks"
        self.clock_dir = self.base_path / clock_name

        if not self.clock_dir.exists():
            raise FileNotFoundError(f"Data for {clock_name} not found at {self.clock_dir}")

        self._load_assets()
    
    def info(self):
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
        
        [Optimization]: 
        Forces all numerical data to `float32`. This reduces memory footprint by 50% 
        and significantly accelerates CPU-bound matrix multiplication operations compared 
        to standard `float64`.
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
        Performs CpG alignment, transposition, and a two-stage missing value imputation.
        
        Strategy:
        1. Intersection: Only process CpGs present in the input data (avoids unnecessary allocation).
        2. Layer 1 Imputation: Fill sporadic NAs in existing columns using Cohort Means.
        3. Layer 2 Imputation: Construct missing columns using Reference Means via efficient concatenation.
        
        Returns:
            pd.DataFrame: A sample-by-CpG matrix aligned with the clock's rotation matrix.
        """
        # 1. Identify Required vs. Missing CpGs
        required_cpgs = self.rotation.index
        
        # Fast Intersection: Filter indices before data extraction to reduce memory overhead
        common_cpgs = beta_df.index.intersection(required_cpgs)
        missing_cpgs = required_cpgs.difference(common_cpgs)
        
        # 2. Extract Existing Data & Transpose
        # Enforce float32 immediately to optimize downstream operations
        X_existing = beta_df.loc[common_cpgs].T.astype(np.float32)
        
        # ============================================================
        # 3. Layer 1 Imputation: Handle Sporadic NAs in Existing Data
        #    (Method: Cohort Mean Imputation with Reference Fallback)
        # ============================================================
        if X_existing.isna().values.any():
            # Identify columns containing NaNs (Boolean indexing is highly efficient)
            nan_cols = X_existing.columns[X_existing.isna().any()]
            
            # Calculate Cohort Mean for these specific columns
            col_means = X_existing[nan_cols].mean(axis=0, skipna=True)
            
            # Fallback: If a column is entirely NaN, use the Reference Mean
            col_means = col_means.fillna(self.impute_ref)
            
            # Apply imputation only to affected columns
            X_existing[nan_cols] = X_existing[nan_cols].fillna(col_means)

        # ============================================================
        # 4. Layer 2 Imputation: Handle Completely Missing Columns
        #    (Method: Vectorized construction & Concatenation)
        # ============================================================
        if len(missing_cpgs) > 0:
            # Retrieve reference means for the missing CpGs
            fill_vals = self.impute_ref.loc[missing_cpgs].values.astype(np.float32)
            
            # Use np.tile for O(1) matrix generation (instant broadcasting)
            # This is significantly faster than iteratively filling a DataFrame
            n_samples = X_existing.shape[0]
            missing_data = np.tile(fill_vals, (n_samples, 1))
            
            # Construct DataFrame for missing parts
            X_missing = pd.DataFrame(
                missing_data, 
                index=X_existing.index, 
                columns=missing_cpgs
            )
            
            # Concatenate existing data with imputed missing data
            X = pd.concat([X_existing, X_missing], axis=1)
        else:
            X = X_existing

        # 5. Final Reorder
        # Ensure columns align strictly with the rotation matrix index
        X = X[required_cpgs]
        
        # 6. Safety Fallback
        # Fill any remaining NaNs with 0.0 (though logically unlikely)
        X = X.fillna(0.0)
            
        return X

    def get_pcs(self, X: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates Principal Components (PCs) using vectorized linear algebra.
        
        Args:
            X (pd.DataFrame): Preprocessed Beta matrix (Samples x CpGs).
            
        Returns:
            pd.DataFrame: Projected PC scores (Samples x PCs).
        """
        # Align center vector to input columns
        center_vec = self.center.loc[X.columns].values
        
        # Use Numpy for broadcast subtraction (Centering)
        X_values = X.values 
        X_centered = X_values - center_vec
        
        # Matrix Multiplication (Projection)
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
        
        # Load Model Coefficients (float32)
        model_df = pd.read_parquet(self.clock_dir / "model.parquet")
        self.coefs = pd.Series(model_df['coef'].values.astype(np.float32), index=model_df['pc'])
        
        # Load Intercept
        with open(self.clock_dir / "intercept.txt", "r") as f:
            self.intercept = float(f.read().strip())

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the list of CpGs required by this PC clock.
        
        Since this is a PC-based clock, exact linear coefficients for CpGs are 
        dynamic (based on projection). This method returns the list of input CpGs 
        defined in the rotation matrix.
        
        Returns:
            pd.DataFrame: columns ['model_name', 'probe', 'coef']
            Note: 'coef' is set to NaN because we are skipping back-projection.
        """
        # 1. Get required CpGs from Rotation Matrix Index
        # This is O(1) - instant access
        required_cpgs = self.rotation.index
        
        # 2. Construct DataFrame
        df = pd.DataFrame({
            'probe': required_cpgs,
            'coef': np.nan  # Placeholder to indicate "Used" but not "Weighted"
        })
        
        
        df['model_name'] = self.name
        
        return df[['model_name', 'probe', 'coef']]   

    
    def anti_trafo(self, x: np.ndarray, adult_age: float = 20) -> np.ndarray:
        """
        Applies the inverse transformation to convert linear scores back to biological age units.
        Used primarily by Horvath-style clocks.
        """
        return np.where(x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age)

    def predict(self, beta_df: pd.DataFrame, verbose: bool = False, **kwargs) -> pd.Series:
        # Preprocess and project data
        X = self.preprocess(beta_df)
        pcs_df = self.get_pcs(X)
        
        # Align coefficients
        common_pcs = self.coefs.index.intersection(pcs_df.columns)
        
        # Calculate raw score (Linear Combination)
        raw_score = pcs_df[common_pcs] @ self.coefs[common_pcs] + self.intercept
        
        # Apply anti-transformation if required
        if self.do_anti_trafo:
            return pd.Series(self.anti_trafo(raw_score.values), index=X.index)
        return raw_score

# --- 4. GrimAge Implementation (Optimized) ---
class PCGrimAge1Impl(BasePCClock):
    """
    Implementation of the PCGrimAge clock.
    
    This complex clock involves a multi-stage process:
    1. Estimating surrogate biomarkers (Surrogates) from DNAm PCs, Age, and Sex.
    2. Predicting composite age from these surrogates.
    """
    def __init__(self, metadata: Dict = None):
        super().__init__("PCGrimAge1", metadata=metadata)
        
        # Load Surrogate Weights (PC -> Surrogate)
        surr_df = pd.read_parquet(self.clock_dir / "surrogate_weights.parquet")
        self.surr_intercepts = surr_df[surr_df['term'] == 'Intercept'].set_index('target')['coef'].astype(np.float32)
        self.surr_weights = surr_df[surr_df['term'] != 'Intercept'].copy()
        self.surr_weights['coef'] = self.surr_weights['coef'].astype(np.float32)
        
        # Load Final Model Weights (Surrogates -> GrimAge)
        final_df = pd.read_parquet(self.clock_dir / "final_model.parquet")
        self.final_intercept = final_df[final_df['term'] == 'Intercept']['coef'].iloc[0]
        self.final_weights = final_df[final_df['term'] != 'Intercept'].set_index('term')['coef'].astype(np.float32)
        
    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the list of features (CpGs + Metadata) required by PCGrimAge.
        """
        # 1. Get required CpGs from Rotation Matrix
        required_cpgs = self.rotation.index
        
        df_cpgs = pd.DataFrame({
            'probe': required_cpgs,
            'coef': np.nan
        })
        
        # 2. Add Metadata Features (Age and Sex are required inputs)
        # This helps users know that these variables are part of the model
        df_meta = pd.DataFrame([
            {'probe': 'Age', 'coef': np.nan},
            {'probe': 'Female', 'coef': np.nan} # Representing Sex
        ])
        
        # 3. Combine
        final_df = pd.concat([df_meta, df_cpgs], ignore_index=True)
        final_df['model_name'] = self.name
        
        return final_df[['model_name', 'probe', 'coef']]
        
    def predict(self, beta_df: pd.DataFrame, ages: pd.Series, sex: pd.Series, verbose: bool = False, **kwargs) -> pd.DataFrame:
        
        # Validate required metadata inputs
        if ages is None or sex is None:
            raise ValueError("[PCGrimAge1] Missing required 'ages' or 'sex' arguments.")

        # --- Step 1: High-Throughput Preprocessing ---
        X = self.preprocess(beta_df)
        pcs_df = self.get_pcs(X)
        samples = X.index
        
        # --- Step 2: Metadata Alignment & Type Casting ---
        # Ensure Age is float for mathematical operations
        ages = ages.reindex(samples).astype(float)
        sex = sex.reindex(samples)

        # --- Step 3: Vectorized Sex Normalization ---
        # Standardize sex input to numeric binary format (Female=1, Male=0)
        sex_cleaned = sex.astype(str).str.strip().str.lower()
        sex_map_to_int = {
            'f': 1, 'female': 1, 'woman': 1, 'w': 1, 
            'm': 0, 'male': 0, 'man': 0,
        }
        is_female = sex_cleaned.map(sex_map_to_int)
        
        # Check for unrecognized sex labels
        if is_female.isna().any():
            invalid_entries = sex[is_female.isna()].unique()
            raise ValueError(
                f"[PCGrimAge1] Input 'sex' contains unrecognized values: {list(invalid_entries)}. "
                f"Supported variants (case-insensitive): {list(sex_map_to_int.keys())}"
            )

        # --- Step 4: Optimized Surrogate Calculation ---
        # Augment the PC matrix with Age and Sex (In-place modification to avoid copy)
        pcs_df['Age'] = ages
        pcs_df['Female'] = is_female
        
        # Pre-allocate results container
        targets = self.surr_weights['target'].unique()
        surrogate_scores = pd.DataFrame(index=samples, columns=targets)
        
        for target in targets:
            # Extract weights for the specific surrogate
            w_df = self.surr_weights[self.surr_weights['target'] == target]
            w_series = pd.Series(w_df['coef'].values, index=w_df['term'])
            intercept = self.surr_intercepts.get(target, 0.0)
            
            # [Optimization]: Use Intersection Indexing
            # Only compute dot product on features present in both the weights and the PC matrix.
            # Implicitly assumes missing features contribute 0.0.
            common_features = w_series.index.intersection(pcs_df.columns)
            
            # Vectorized Dot Product
            val = pcs_df[common_features] @ w_series[common_features] + intercept
            surrogate_scores[target] = val
            
        # --- Step 5: Feature Aggregation & Nomenclature Mapping ---
        # Prepare feature matrix for the final model
        final_features = surrogate_scores 
        final_features['Age'] = ages
        final_features['Female'] = is_female
        
        # Map calculated PC-based names to the DNAm-based names expected by the final model
        rename_map = {
            'PCPACKYRS': 'DNAmPACKYRS',
            'PCADM': 'DNAmADM',
            'PCB2M': 'DNAmB2M',
            'PCCystatinC': 'DNAmCystatinC',
            'PCGDF15': 'DNAmGDF15',
            'PCLeptin': 'DNAmLeptin',
            'PCPAI1': 'DNAmPAI1',
            'PCTIMP1': 'DNAmTIMP1'
        }
        # Note: This operation creates a copy due to renaming
        final_features_renamed = final_features.rename(columns=rename_map)

        # --- Step 6: Final Composite Age Prediction ---
        # Align features with the final linear model
        aligned_features = final_features_renamed.reindex(columns=self.final_weights.index, fill_value=0)
        
        pred_age = aligned_features @ self.final_weights + self.final_intercept
        
        # Append final GrimAge score to the results
        surrogate_scores[self.name] = pred_age
        
        return surrogate_scores

# --- 5. Export Wrapper Classes ---
# These classes hardcode parameters for specific clocks to allow parameter-less instantiation.

class PCHorvath2013(StandardPCClock):
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1038/s43587-022-00248-2",
        "original_clock": "Horvath2013"
    }
    def __init__(self):
        super().__init__("PCHorvath2013", do_anti_trafo=True, metadata=self.METADATA)

class PCHorvath2018(StandardPCClock):
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Skin/Blood",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1038/s43587-022-00248-2",
        "original_clock": "Horvath2018"
    }
    def __init__(self):
        super().__init__("PCHorvath2018", do_anti_trafo=True, metadata=self.METADATA)

class PCHannum(StandardPCClock):
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1038/s43587-022-00248-2",
        "original_clock": "Hannum2013"
    }
    def __init__(self):
        super().__init__("PCHannum", do_anti_trafo=False, metadata=self.METADATA)

class PCPhenoAge(StandardPCClock):
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Mortality",
        "source": "https://doi.org/10.1038/s43587-022-00248-2",
        "original_clock": "PhenoAge"
    }
    def __init__(self):
        super().__init__("PCPhenoAge", do_anti_trafo=False, metadata=self.METADATA)

class PCDNAmTL(StandardPCClock):
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "Telomere Length",
        "source": "https://doi.org/10.1038/s43587-022-00248-2",
        "original_clock": "DNAmTL"
    }
    def __init__(self):
        super().__init__("PCDNAmTL", do_anti_trafo=False, metadata=self.METADATA)

class PCGrimAge1(PCGrimAge1Impl):
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Mortality",
        "source": "https://doi.org/10.1038/s43587-022-00248-2",
        "original_clock": "GrimAge1"
    }

    def __init__(self):
        super().__init__(metadata=self.METADATA)