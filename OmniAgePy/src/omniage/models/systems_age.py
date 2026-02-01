import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict
from ..utils import check_and_download_zenodo
SYSTEMS_AGE_URL = "https://zenodo.org/records/18287372/files/SystemsAge.tar.gz?download=1"

class SystemsAge:
    """
    Implements the Systems Age epigenetic clock (Sehgal et al. 2025).
    
    This is a multi-stage, PC-based clock that calculates biological age for **11 distinct 
    physiological systems** and aggregates them into a single composite 'Systems Age'.
    
    **Algorithm Pipeline:**
    1.  **Level 1**: Project CpG methylation into DNAm Principal Components (PCs).
    2.  **Level 2**: Project DNAm PCs into System PCs using system-specific vectors.
    3.  **Level 3**: Calculate raw scores for 11 systems (Blood, Brain, Heart, etc.).
    4.  **Level 4**: Calculate Composite 'SystemsAge' via PCA of system scores.
    5.  **Final**: Normalize all scores to biological age (years) using Z-score transformation.

    Attributes:
        name (str): Model name.
        metadata (dict): Model metadata.
        base_dir (Path): Path to model assets.
    
    References:
        Sehgal, R. et al. Systems Age: a single blood methylation test to quantify aging heterogeneity across 11 physiological systems.
        Nature Aging (2025). https://doi.org/10.1038/s43587-025-00958-3
    """
    
    
    METADATA = {
        "year": 2025,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm (450k)",
        "prediction": "Mortality",
        "source": "https://doi.org/10.1038/s43587-025-00958-3"
    }

    def __init__(self):
        """Initialize the SystemsAge model and load all projection matrices."""
        self.name = "SystemsAge"
        self.metadata = self.METADATA
        
        # 1. Locate Data Directory
        # Assumes data is in ../data/SystemsAge relative to this file
        current_file_dir = Path(__file__).parent.resolve()
        self.base_dir = current_file_dir.parent / "data" / "SystemsAge"

        check_and_download_zenodo(
            target_dir=self.base_dir, 
            download_url=SYSTEMS_AGE_URL
        )
        
        if not self.base_dir.exists():
            raise FileNotFoundError(f"SystemsAge data not found at {self.base_dir}")
            
        self._load_assets()

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the list of CpGs required by the SystemsAge model.
        
        Note: SystemsAge is a multi-stage PC-based model (CpG -> DNAmPC -> SystemPC -> Score).
        Therefore, it implies dynamic weighting rather than static linear coefficients.
        This method returns the input feature space (CpGs used in the primary PCA).
        
        Returns:
            pd.DataFrame: columns ['model_name', 'probe', 'coef']
            ('coef' is set to NaN as per PC clock convention)
        """
        # 1. Identify Input Features
        # The first layer of the model is the DNAm PCA Projection.
        # Any CpG in this rotation matrix is a required input.
        required_cpgs = self.dnam_pca_rot.index
        
        # 2. Construct DataFrame
        df = pd.DataFrame({
            'probe': required_cpgs,
            'coef': np.nan  # Placeholder
        })
        
        # 3. Add Model Name
        df['model_name'] = self.name
        
        return df[['model_name', 'probe', 'coef']]
    

    def info(self):
        """Prints metadata information about the clock."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def _load_assets(self):
        """
        Loads all required Parquet assets (weights, rotations, scalers).
        
        [Optimization]: 
        Forces `float32` precision for all numerical matrices. This reduces memory 
        bandwidth requirements and accelerates vector operations.
        """
        
        # A. Imputation (Reference Means)
        self.impute_ref = pd.read_parquet(self.base_dir / "impute_means.parquet").set_index("probe")["mean"].astype(np.float32)
        
        # B. DNAm PCA (Level 1 Projection: CpG -> DNAmPC)
        self.dnam_pca_center = pd.read_parquet(self.base_dir / "dnam_pca_center.parquet").set_index("probe")["mean"].astype(np.float32)
        self.dnam_pca_rot = pd.read_parquet(self.base_dir / "dnam_pca_rotation.parquet").set_index("probe").astype(np.float32)
        
        # C. System Vector Coefficients (Level 2 Projection: DNAmPC -> SystemPC)
        df_sys_vec = pd.read_parquet(self.base_dir / "system_vector_coefs.parquet")
        if 'pc' in df_sys_vec.columns:
            self.sys_vec_coefs = df_sys_vec.set_index("pc")
        else:
            self.sys_vec_coefs = df_sys_vec
        self.sys_vec_coefs = self.sys_vec_coefs.astype(np.float32)
            
        # D. System Score Coefs (Level 3: SystemPC -> Raw Score)
        self.sys_score_coefs = pd.read_parquet(self.base_dir / "system_score_coefs.parquet").set_index("term")["coef"].astype(np.float32)
        
        # E. Chronological Age Prediction Model
        age_coefs_df = pd.read_parquet(self.base_dir / "age_pred_coefs.parquet")
        self.age_intercept = age_coefs_df.loc[age_coefs_df['term'] == 'Intercept', 'coef'].iloc[0]
        self.age_linear_weights = age_coefs_df[age_coefs_df['term'] != 'Intercept'].set_index("term")["coef"].astype(np.float32)
        
        # Quadratic Model Params
        age_params_df = pd.read_parquet(self.base_dir / "age_model_params.parquet").set_index("param")
        self.age_model_const = age_params_df.iloc[0, 0]
        self.age_model_lin = age_params_df.iloc[1, 0]
        self.age_model_quad = age_params_df.iloc[2, 0]

        # F. Systems PCA (Level 4: System Scores -> Composite SystemsAge)
        self.sys_pca_center = pd.read_parquet(self.base_dir / "systems_pca_center.parquet").set_index("term")["mean"].astype(np.float32)
        self.sys_pca_scale = pd.read_parquet(self.base_dir / "systems_pca_scale.parquet").set_index("term")["scale"].astype(np.float32)
        self.sys_pca_rot = pd.read_parquet(self.base_dir / "systems_pca_rotation.parquet").set_index("term").astype(np.float32)
        
        final_coefs_df = pd.read_parquet(self.base_dir / "final_coefs.parquet")
        self.final_coefs = pd.Series(final_coefs_df['coef'].values, index=final_coefs_df['pc']).astype(np.float32)

        # G. Final Transformation Coefficients
        self.trans_coefs = pd.read_parquet(self.base_dir / "transformation_coefs.parquet").set_index("system").astype(np.float32)

    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        [High-Performance Preprocessing] 
        Replicates the optimized R logic: Intersection -> Split -> Impute -> Concat.
        
        Strategy:
        1. Intersection: Only process CpGs present in the input to save memory.
        2. Existing Data: Impute sporadic NaNs using the cohort mean.
        3. Missing Columns: Construct using reference means via O(1) Numpy tiling.
        4. Concatenation: Merge and reorder.
        """
        # 1. Identify Required CpGs
        required_cpgs = self.dnam_pca_rot.index 
        
        # 2. Set Operations: Identify Present vs. Absent Features
        common_cpgs = beta_df.index.intersection(required_cpgs)
        missing_cpgs = required_cpgs.difference(common_cpgs)
        
        # 3. Extract Existing Data (Transposing & Casting to float32)
        X_existing = beta_df.loc[common_cpgs].T.astype(np.float32)
        
        # ============================================================
        # 4. Layer 1 Imputation: Existing Data (Sporadic Missingness)
        # ============================================================
        if X_existing.isna().values.any():
            nan_cols = X_existing.columns[X_existing.isna().any()]
            col_means = X_existing[nan_cols].mean(axis=0, skipna=True)
            col_means = col_means.fillna(self.impute_ref)
            X_existing[nan_cols] = X_existing[nan_cols].fillna(col_means)
            
        # ============================================================
        # 5. Layer 2 Imputation: Missing Columns (Structural Missingness)
        # ============================================================
        if len(missing_cpgs) > 0:
            fill_vals = self.impute_ref.loc[missing_cpgs].values
            n_samples = X_existing.shape[0]
            
            # Broadcasting via np.tile (Instant O(1) operation)
            missing_data = np.tile(fill_vals, (n_samples, 1))
            
            X_missing = pd.DataFrame(
                missing_data, 
                index=X_existing.index, 
                columns=missing_cpgs,
                dtype=np.float32
            )
            
            # ============================================================
            # 6. Final Assembly
            # ============================================================
            X = pd.concat([X_existing, X_missing], axis=1)
        else:
            X = X_existing

        # 7. Final Reorder
        X = X[required_cpgs]
        
        # 8. Fallback
        X = X.fillna(0.0)
        
        return X

    def predict(self, beta_df: pd.DataFrame, verbose: bool = False, **kwargs) -> pd.DataFrame:
        """
        Calculates Systems Age for all physiological systems.

        Args:
            beta_df (pd.DataFrame): Methylation data (Rows: CpGs, Columns: Samples).
            verbose (bool): Whether to print progress messages.

        Returns:
            pd.DataFrame: A DataFrame (Samples x Systems) containing biological ages (in Years).
            
            **Output Columns:**
            - ``SystemsAge``: The overall composite biological age.
            - ``Blood``, ``Brain``, ``Heart``, ``Lung``, ``Kidney``, ``Liver``: Organ-specific ages.
            - ``Inflammation``, ``Immune``, ``Metabolic``, ``Hormone``, ``MusculoSkeletal``: System-specific ages.
            - ``Age_prediction``: Predicted chronological age (intermediate step).
        """
        
        # --- Step 1: Preprocessing ---
        X = self.preprocess(beta_df)
        samples = X.index
        
        # --- Step 2: DNAm PCA (Level 1 Projection) ---
        X_centered = X - self.dnam_pca_center
        dnam_pcs = X_centered @ self.dnam_pca_rot
        
        # --- Step 3: DNAm System PCs (Level 2 Projection) ---
        # Strictly slice first 4017 PCs as per R implementation
        dnam_pcs_subset = dnam_pcs.iloc[:, :4017]
        sys_vec_subset = self.sys_vec_coefs.iloc[:4017]
        
        dnam_sys_pcs = dnam_pcs_subset @ sys_vec_subset
        
        # --- Step 4: Calculate System Scores ---
        system_groups = ["Blood", "Brain", "Cytokine", "Heart", "Hormone", "Immune", 
                         "Kidney", "Liver", "Metab", "Lung", "MusculoSkeletal"]
        
        system_scores = pd.DataFrame(index=samples, columns=system_groups)
        
        for group in system_groups:
            relevant_cols = [c for c in dnam_sys_pcs.columns if group in c]
            
            sub_matrix = dnam_sys_pcs[relevant_cols]
            sub_coefs = self.sys_score_coefs.loc[relevant_cols]
            
            if len(relevant_cols) == 1:
                score = sub_matrix.iloc[:, 0] * -1.0
            else:
                score = sub_matrix @ sub_coefs
            
            system_scores[group] = score

        # --- Step 5: Predicted Chronological Age ---
        common_age_pcs = dnam_pcs.columns.intersection(self.age_linear_weights.index)
        
        age_pred_raw = dnam_pcs[common_age_pcs] @ self.age_linear_weights.loc[common_age_pcs] + self.age_intercept
        
        # Quadratic Correction
        age_final = (age_pred_raw * self.age_model_lin) + \
                    (age_pred_raw**2 * self.age_model_quad) + \
                    self.age_model_const
        
        age_final_years = age_final / 12.0
        system_scores["Age_prediction"] = age_final_years
        
        rename_map = {"Cytokine": "Inflammation", "Metab": "Metabolic"}
        system_scores = system_scores.rename(columns=rename_map)
        
        # --- Step 6: Composite SystemsAge (Level 3 Projection) ---
        ordered_cols = ["Blood", "Brain", "Inflammation", "Heart", "Hormone", "Immune", 
                        "Kidney", "Liver", "Metabolic", "Lung", "MusculoSkeletal", "Age_prediction"]
        
        scores_for_pca = system_scores[ordered_cols]
        scores_centered = scores_for_pca - self.sys_pca_center
        
        if hasattr(self, 'sys_pca_scale') and not self.sys_pca_scale.empty:
             scores_scaled = scores_centered / self.sys_pca_scale
        else:
             scores_scaled = scores_centered

        sys_pca_res = scores_scaled @ self.sys_pca_rot
        
        final_pcs = sys_pca_res.columns.intersection(self.final_coefs.index)
        composite_age = sys_pca_res[final_pcs] @ self.final_coefs.loc[final_pcs]
        
        system_scores["SystemsAge"] = composite_age

        # --- Step 7: Final Scaling (Normalization: Z-score -> Years) ---
        final_columns = ordered_cols + ["SystemsAge"]
        system_scores = system_scores[final_columns]

        for i, col_name in enumerate(final_columns):
            y = system_scores[col_name]
            
            # Use iloc to access coefficients by position
            coef_row = self.trans_coefs.iloc[i]
            v1, v2, v3, v4 = coef_row.iloc[0], coef_row.iloc[1], coef_row.iloc[2], coef_row.iloc[3]
            
            # Transformation
            val = (((y - v1) / v2) * v4) + v3
            val = val / 12.0
            
            system_scores[col_name] = val
            
        return system_scores