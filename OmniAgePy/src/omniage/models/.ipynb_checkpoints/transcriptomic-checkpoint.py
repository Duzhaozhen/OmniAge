import pandas as pd
import numpy as np
import os
from typing import List, Optional, Dict, Tuple, Union
from ..utils import get_data_path
import scipy.sparse
import scanpy as sc
import warnings
import anndata

def scImmuAging_generate_pseudocells(adata: anndata.AnnData, size: int = 15, n_repeats: int = 100) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Implements the Pseudocell Bootstrapping logic used in scImmuAging.
    
    This function generates pseudocells by randomly sampling single cells from each donor.
    It serves as a preprocessing step for transcriptomic clocks requiring robust, 
    aggregated expression profiles.
    
    Logic:
    For each unique donor (identified by 'donor_id' and 'age'):
      1. Repeat the process `n_repeats` times.
      2. Randomly sample `size` cells (with replacement if the cell count < `size`).
      3. Compute the mean expression profile of these sampled cells.
      
    Args:
        adata (AnnData): Input single-cell object. Must contain 'donor_id' and 'age' in .obs.
        size (int): Number of cells to aggregate per pseudocell.
        n_repeats (int): Number of pseudocells to generate per donor.
        
    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: 
            - X_pseudo: Pseudocell expression matrix (Rows=Pseudocells, Cols=Genes).
            - meta_pseudo: Corresponding metadata (donor_id, age).
    """
    pseudocells = []
    metadata = []
    
    # Metadata Validation
    if "donor_id" not in adata.obs.columns:
        raise ValueError("AnnData.obs must contain 'donor_id' column.")
    if "age" not in adata.obs.columns:
        raise ValueError("AnnData.obs must contain 'age' column.")

    # Group by Donor for efficient processing
    grouped = adata.obs.groupby(['donor_id', 'age'])
    
    for (donor, age), group_indices in grouped.indices.items():
        n_cells = len(group_indices)
        replace = n_cells < size
        
        # Extract expression matrix for the current donor
        # Optimize handling for sparse vs dense matrices
        if hasattr(adata.X, "toarray"):
            donor_X = adata.X[group_indices].toarray()
        else:
            donor_X = adata.X[group_indices]
            
        donor_df = pd.DataFrame(donor_X, columns=adata.var_names)
        
        # Bootstrap Loop
        donor_pseudos = []
        for _ in range(n_repeats):
            # Sample indices
            sample_idx = np.random.choice(n_cells, size=size, replace=replace)
            # Calculate mean expression (Pseudocell profile)
            pseudo_expr = donor_df.iloc[sample_idx].mean(axis=0)
            donor_pseudos.append(pseudo_expr)
        
        # Collect results
        pseudocells.extend(donor_pseudos)
        metadata.extend([{'donor_id': donor, 'age': age}] * n_repeats)
        
    # Return empty if no pseudocells generated
    if not pseudocells:
        return pd.DataFrame(), pd.DataFrame()

    X_pseudo = pd.DataFrame(pseudocells)
    meta_pseudo = pd.DataFrame(metadata)
    
    return X_pseudo, meta_pseudo


class scImmuAging:
    """
    Implements scImmuAging (Li et al. 2025).
    
    A single-cell transcriptomic clock that predicts biological age for specific 
    immune cell types.
    
    Mechanism:
    1. Preprocessing: Generates 'Pseudocells' via bootstrapping.
    2. Prediction: Applies cell-type-specific ElasticNet/Lasso linear models.
    3. Aggregation: Averages pseudocell predictions to produce a donor-level age.
    
    Supported Cell Types: CD4T, CD8T, MONO, NK, B.
    
    References:
        Li W, et al. Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and vaccination. Nat Aging (2025).
        https://doi.org/10.1038/s43587-025-00819-z
    """
    METADATA = {
        "year": 2025,
        "species": "Human",
        "tissue": "PBMC",
        "omic type": "Transcriptomics(scRNA-seq)",
        "prediction": "Chronological Age(Years)",
        "source": "https://doi.org/10.1038/s43587-025-00819-z"
    }
    def __init__(self):
        self.name = "scImmuAging"
        self.metadata = self.METADATA
        self.data_dir = get_data_path("scImmuAging")
        
        if not os.path.exists(self.data_dir):
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")
            
        # Preload models for all available cell types
        self.models = {}
        
        if os.path.exists(self.data_dir):
            for filename in os.listdir(self.data_dir):
                if filename.endswith(".csv"):
                    ct_name = filename.replace(".csv", "")
                    try:
                        df = pd.read_csv(os.path.join(self.data_dir, filename))
                        
                        # Parse coefficients and intercept
                        intercept_mask = df['gene'] == 'Intercept'
                        if intercept_mask.any():
                            intercept = float(df[intercept_mask]['coef'].iloc[0])
                            weights = df[~intercept_mask].set_index('gene')['coef']
                        else:
                            intercept = 0.0
                            weights = df.set_index('gene')['coef']
                        
                        self.models[ct_name] = {
                            "intercept": intercept,
                            "weights": weights
                        }
                    except Exception as e:
                        print(f"[scImmuAging] Warning: Failed to load model {filename}: {e}")
    
    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves coefficients for all cell-type-specific models.
        
        Returns:
            pd.DataFrame: columns ['model_name', 'gene', 'coef']
            The 'model_name' column will be formatted as 'scImmuAging_{CellType}' (e.g., scImmuAging_CD4T).
        """
        if not self.models:
            return pd.DataFrame(columns=['model_name', 'gene', 'coef'])

        dfs = []
        
        for ct_name, params in self.models.items():
            # 1. Extract weights (Genes)
            weights_series = params['weights']
            df_genes = pd.DataFrame({
                'gene': weights_series.index,
                'coef': weights_series.values
            })
            
            # 2. Extract Intercept
            # We add it as a row with gene='(Intercept)'
            df_intercept = pd.DataFrame([{
                'gene': '(Intercept)',
                'coef': params['intercept']
            }])
            
            # 3. Combine
            df_ct = pd.concat([df_intercept, df_genes], ignore_index=True)
            
            # 4. Add Model Name (e.g., scImmuAging_CD4T)
            df_ct['model_name'] = f"{self.name}_{ct_name}"
            
            dfs.append(df_ct)
            
        # Combine all cell types
        final_df = pd.concat(dfs, ignore_index=True)
        
        return final_df[['model_name', 'gene', 'coef']]
        
    def info(self):
        """Prints metadata information about the scImmuAging clock."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def predict(
        self, 
        adata: anndata.AnnData, 
        cell_types: List[str] = ["CD4T", "CD8T", "MONO", "NK", "B"], 
        verbose: bool = True
    ) -> Dict[str, Dict[str, pd.DataFrame]]:
        """
        Executes the prediction pipeline: Subset -> Pseudocell -> Predict -> Aggregate.
        
        Args:
            adata (AnnData): **Log-Normalized** expression object.
                - **.var_names**: Must use Gene Symbols (e.g., 'RPS4Y1', 'CD74') to match the model.
                - **.obs Requirements**:
                    - ``'donor_id'``: For donor aggregation.
                    - ``'age'``: For reference/plotting.
                    - ``'celltype'``: For subsetting. Values must match ``cell_types`` list.
            cell_types (List[str], optional): Cell populations to predict. 
                Defaults to ["CD4T", "CD8T", "MONO", "NK", "B"].
            verbose (bool, optional): If True, prints progress.
            
        Returns:
            Dict: A nested dictionary containing results per cell type:
                {
                    "CD4T": {
                        "BootstrapCell": pd.DataFrame (Pseudocell-level predictions),
                        "Donor": pd.DataFrame (Aggregated donor-level predictions)
                    },
                    "MONO": { ... }
                    ...
                }
        """
        results = {}
        
        if "celltype" not in adata.obs.columns:
            raise ValueError("AnnData.obs must contain 'celltype' column.")

        for ct in cell_types:
            if ct not in self.models:
                if verbose: print(f"[scImmuAging] Warning: No model found for cell type '{ct}'. Skipping.")
                continue
                
            if verbose: print(f"\n--- Processing cell type: {ct} ---")
            
            # 1. Subset Data
            subset_idx = adata.obs['celltype'] == ct
            if subset_idx.sum() == 0:
                if verbose: print(f"No cells found for {ct}.")
                continue
                
            adata_sub = adata[subset_idx]
            
            # 2. Pre-process (Generate Pseudocells)
            if verbose: print("Generating pseudocells (Bootstrap)...")
            X_pseudo, meta_pseudo = scImmuAging_generate_pseudocells(adata_sub)
            
            if X_pseudo.empty:
                continue

            # 3. Predict (Linear Model)
            if verbose: print("Predicting age...")
            model_params = self.models[ct]
            weights = model_params['weights']
            intercept = model_params['intercept']
            
            # Feature alignment (fill missing genes with 0)
            X_aligned = X_pseudo.reindex(columns=weights.index, fill_value=0)
            
            # Calculate raw predictions
            raw_preds = X_aligned.dot(weights) + intercept
            
            # Compile BootstrapCell results
            bootstrap_res = meta_pseudo.copy()
            bootstrap_res['Prediction'] = raw_preds.values
            
            # 4. Aggregate per Donor
            if verbose: print("Aggregating per donor...")
            donor_res = bootstrap_res.groupby('donor_id').agg({
                'age': 'first',          # Age is constant per donor
                'Prediction': 'mean'     # Mean of pseudocell predictions
            }).rename(columns={'Prediction': 'predicted'}).reset_index()
            
            # Rounding (matching original R implementation)
            donor_res['predicted'] = donor_res['predicted'].round()
            
            # Store results
            results[ct] = {
                "BootstrapCell": bootstrap_res,
                "Donor": donor_res
            }
            
            if verbose: print(f"Prediction complete for {ct}!")
            
        return results


class BrainCTClock:
    """
    Human Brain Cell-Type-Specific Aging Clocks (Muralidharan et al. 2025).
    
    This class predicts biological age from single-cell transcriptomic data using 
    cell-type-specific elastic net models. It supports robust prediction via 
    bootstrapping and imputation of missing genes.

    **Operational Modes:**
    1.  **SC (Single-Cell)**: 
        Predicts age for every single cell individually. High resolution but noisy.
    2.  **Pseudobulk**: 
        Aggregates all cells of a specific type per donor into one profile before prediction. 
        Most robust for donor-level aging assessment.
    3.  **Bootstrap**: 
        Resamples cells (with replacement) to create multiple "pseudocells" per donor. 
        Provides a distribution of predicted ages and confidence intervals.

    Attributes:
        name (str): Model name.
        metadata (dict): Model metadata.
        models (dict): Loaded coefficients for each cell type (5 folds per type).
        imputation_data (dict): Reference values for imputing missing genes.

    References:
        Muralidharan, et al. Human Brain Cell-Type-Specific Aging Clocks Based on Single-Nuclei Transcriptomics. 
        Advanced Science (2025). https://doi.org/10.1002/advs.202506109
    """
    METADATA = {
        "year": 2025,
        "species": "Human",
        "tissue": "Brain",
        "omic type": "Transcriptomics (scRNA-seq)",
        "prediction": "Chronological Age(Years)",
        "source": "https://doi.org/10.1002/advs.202506109"
    }

    # Cell count thresholds for bootstrapping
    BOOTSTRAP_THRESHOLDS = {
        "Oligodendrocytes": 200,
        "Astrocytes": 50,
        "Microglia": 50,
        "OPCs": 50,
        "Excitatory Neurons": 100,
        "Inhibitory Neurons": 100
    }

    def __init__(self, data_dir: str = None):
        """
        Initialize the BrainCTClock by loading model coefficients and imputation data.
        """
        self.name = "BrainCTClock"
        self.metadata = self.METADATA
        self.models = {} 
        self.imputation_data = {}
        
        # --- 1. Automatic Path Resolution ---
        if data_dir is None:
            # Locate package data directory
            current_file_path = os.path.abspath(__file__)
            models_dir = os.path.dirname(current_file_path)
            package_root = os.path.dirname(models_dir)
            data_dir = os.path.join(package_root, "data", "Brain_CT_Clock")
        
        if os.path.exists(data_dir):
            self._load_models(data_dir)
        else:
            print(f"[Error] BrainCT data directory not found: {data_dir}")

    
    
    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the coefficients for every fold of each model.
        
        Returns:
            pd.DataFrame: columns ['model_name', 'gene', 'coef']
            (Note: Changed from 'probe' to 'gene' to reflect transcriptomic nature)
        """
        if not self.models:
            return pd.DataFrame(columns=['model_name', 'gene', 'coef'])

        results = []

    
        for model_key, folds_list in self.models.items():
            if not folds_list:
                continue

            for i, fold_df in enumerate(folds_list):
                df = fold_df.copy()


                df = df.rename(columns={
                    "feature_name": "gene",
                    "coefficient": "coef"
                })

                fold_id = i + 1
                df["model_name"] = f"{self.name}_{model_key}_fold{fold_id}"

                if 'model_name' in df.columns and 'gene' in df.columns and 'coef' in df.columns:
                    results.append(df[['model_name', 'gene', 'coef']])

        if not results:
             return pd.DataFrame(columns=['model_name', 'gene', 'coef'])

        final_df = pd.concat(results, ignore_index=True)
        return final_df

    
    def info(self):
        """Prints metadata information about the BrainCTClock."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")
    
    def _load_models(self, data_dir):
        """Helper to load and restructure model coefficients from CSV."""
        coef_path = os.path.join(data_dir, "brain_ct_coefs.csv")
        impute_path = os.path.join(data_dir, "brain_ct_imputation.csv")

        # --- Load Imputation Data ---
        if os.path.exists(impute_path):
            df_imp_all = pd.read_csv(impute_path)
            if "model_key" in df_imp_all.columns:
                for model_key, group in df_imp_all.groupby("model_key"):
                    self.imputation_data[model_key] = group[["feature_name", "imputation_value"]].reset_index(drop=True)
            else:
                print(f"[Error] 'model_key' column missing in {impute_path}")
        else:
            print(f"[Warning] Imputation file not found: {impute_path}")

        # --- Load Coefficients ---
        if os.path.exists(coef_path):
            df_coef_all = pd.read_csv(coef_path)
            
            if "model_key" in df_coef_all.columns and "fold" in df_coef_all.columns:
                for model_key, model_group in df_coef_all.groupby("model_key"):
                    fold_list = []
                    unique_folds = model_group["fold"].unique()
                    
                    # Smart sort for folds (e.g., fold_1, fold_2, fold_10)
                    try:
                        sorted_folds = sorted(unique_folds, key=lambda x: int(str(x).split('_')[-1]) if '_' in str(x) else x)
                    except:
                        sorted_folds = sorted(unique_folds)

                    for fold in sorted_folds:
                        fold_df = model_group[model_group["fold"] == fold]
                        clean_df = fold_df[["feature_name", "coefficient"]].reset_index(drop=True)
                        fold_list.append(clean_df)
                    
                    self.models[model_key] = fold_list
            else:
                print(f"[Error] 'model_key' or 'fold' column missing in {coef_path}")
        else:
            print(f"[Error] Coefficient file not found: {coef_path}")

    def predict(
        self, 
        adata: anndata.AnnData, 
        cell_types: List[str], 
        model_name: Union[str, List[str]] = "all"
    ) -> Dict[str, pd.DataFrame]:
        """
        Executes the BrainCT prediction pipeline.
        
        This method automatically handles:
        1.  **Gene Matching**: Case-insensitive matching of Gene Symbols.
        2.  **Imputation**: Filling missing genes with training set mean values.
        3.  **Ensembling**: Averaging predictions across 5 cross-validation folds.

        Args:
            adata (anndata.AnnData): Input data object.
                - **.X**: Should be **Log-Normalized** expression data.
                - **.var_names**: Must be **Gene Symbols** (e.g., 'APOE', 'TREM2').
                - **.obs Requirements**:
                    - ``'donor_id'``: Unique donor identifier.
                    - ``'age'``: Chronological age.
                    - ``'celltype'``: Cell type labels matching ``cell_types`` list.
            cell_types (List[str]): List of cell types to predict (e.g., ['Microglia', 'Astrocytes']).
            model_name (Union[str, List[str]], optional): Mode(s) to run. 
                Options: ``['SC', 'Pseudobulk', 'Bootstrap']``. Defaults to "all".

        Returns:
            Dict[str, pd.DataFrame]: A dictionary where keys are the mode names (e.g., 'SC', 'Pseudobulk').
            Each value is a DataFrame containing:
                - ``predictions``: Predicted biological age.
                - ``ages``: True chronological age (from metadata).
                - ``donors``: Donor ID.
                - ``celltype``: The cell type associated with the prediction.
                - ``sample_type``: The mode used (SC/Pseudobulk/Bootstrap).
        
        Warning:
            If your ``adata.var_names`` are Ensembl IDs (e.g., ENSG000...), prediction 
            will fail or rely solely on the intercept. Please map to Symbols first.
        """
        # --- 1. Standard Metadata Column Names ---
        donor_col = "donor_id"
        age_col = "age"
        cell_type_col = "celltype"

        # --- 2. Validation ---
        required_cols = [donor_col, age_col, cell_type_col]
        missing = [col for col in required_cols if col not in adata.obs.columns]
        
        if missing:
            raise ValueError(
                f"Input AnnData is missing required metadata columns: {missing}.\n"
                f"Please rename your columns to match: 'donor_id', 'age', 'celltype'."
            )

        # --- 3. Parse Model Types ---
        valid_models = ["SC", "Pseudobulk", "Bootstrap"]
        if isinstance(model_name, str):
            models_to_run = valid_models if model_name == "all" else [model_name]
        else:
            models_to_run = model_name

        invalid = [m for m in models_to_run if m not in valid_models]
        if invalid:
            raise ValueError(f"Invalid model_name(s): {invalid}. Valid options: {valid_models}")

        all_results = {}
        print("Starting Brain_CT_clock prediction...")

        # --- 4. Execution Loop ---
        for current_model_type in models_to_run:
            print(f"Running mode: {current_model_type}")
            
            res = self._run_prediction_pipeline(
                sample_type=current_model_type,
                adata=adata,
                common_celltypes=cell_types,
                donor_col=donor_col,
                age_col=age_col,
                cell_type_col=cell_type_col
            )
            all_results[current_model_type] = res

        print("--- All Brain_CT_clock predictions complete! ---")
        return all_results

    def _run_prediction_pipeline(self, sample_type, adata, common_celltypes, 
                                 donor_col, age_col, cell_type_col):
        """
        Executes the prediction pipeline (5-fold averaging) for a specific model type.
        """
        final_results_list = []

        for cell_type in common_celltypes:
            clock_key = f"{sample_type}_{cell_type}"
            
            if clock_key not in self.models:
                print(f"Skipping {cell_type}: Model key {clock_key} not found.")
                continue
            
            model_folds = self.models[clock_key]
            impute_df = self.imputation_data.get(clock_key, None)

            # 2. Data Preparation
            df_base = self._prepare_data(adata, cell_type, sample_type, 
                                         donor_col, age_col, cell_type_col)

            if df_base.empty:
                continue

            # 3. 5-Fold Loop Prediction
            fold_preds = []
            for i, model_fold in enumerate(model_folds):
                pred_res = self._predict_core(
                    data=df_base.drop(columns=[donor_col, age_col, cell_type_col], errors='ignore'), 
                    impute_data=impute_df,
                    model=model_fold,
                    sample_type=sample_type
                )
                fold_preds.append(pred_res['predictions'].values)

            # 4. Average Predictions across folds
            avg_preds = np.column_stack(fold_preds).mean(axis=1)

            # 5. Construct Result DataFrame
            # Ensures correct row alignment using df_base
            res_df = pd.DataFrame({
                "predictions": avg_preds,
                "ages": df_base[age_col].values,      
                "donors": df_base[donor_col].values,  
                "sample_type": sample_type,
                "celltype": cell_type
            })
            final_results_list.append(res_df)

        if not final_results_list:
            return pd.DataFrame()

        return pd.concat(final_results_list, ignore_index=True)

    def _prepare_data(self, adata, cell_type, sample_type, 
                      donor_col, age_col, cell_type_col):
        """
        Handles data extraction logic for SC, Pseudobulk, and Bootstrap modes.
        """
        if cell_type_col not in adata.obs:
             return pd.DataFrame()
             
        subset = adata[adata.obs[cell_type_col] == cell_type].copy()
        if subset.n_obs == 0:
            return pd.DataFrame()

        if hasattr(subset.X, "toarray"):
            X = subset.X.toarray()
        else:
            X = subset.X
            
        df = pd.DataFrame(X, index=subset.obs_names, columns=subset.var_names)
        
        # Attach metadata
        meta_cols = [donor_col, age_col, cell_type_col]
        for col in meta_cols:
            df[col] = subset.obs[col].values

        gene_cols = subset.var_names.tolist()

        if sample_type == "SC":
            return df

        elif sample_type == "Pseudobulk":
            # Group by donor, age, celltype and compute mean expression
            df_ct = df.groupby([donor_col, age_col, cell_type_col], observed=True)[gene_cols].mean().reset_index()
            return df_ct

        elif sample_type == "Bootstrap":
            return self._bootstrap_sampling(df, cell_type, gene_cols, 
                                            donor_col, age_col, cell_type_col)
        
        return pd.DataFrame()

    def _bootstrap_sampling(self, df, cell_type, gene_cols, 
                            donor_col, age_col, cell_type_col):
        """
        Implements bootstrap sampling logic.
        Samples 100 replicates per donor.
        """
        threshold = self.BOOTSTRAP_THRESHOLDS.get(cell_type, 50)
        donors = df[donor_col].unique()
        boot_rows = []

        np.random.seed(42)

        for donor in donors:
            df_donor = df[df[donor_col] == donor]
            n_cells = len(df_donor)
            
            expr_matrix = df_donor[gene_cols].values

            if n_cells > threshold:
                for _ in range(100):
                    indices = np.random.choice(n_cells, size=threshold, replace=False)
                    sample_mean = expr_matrix[indices, :].mean(axis=0)
                    boot_rows.append(sample_mean)
            else:
                true_mean = expr_matrix.mean(axis=0)
                for _ in range(100):
                    boot_rows.append(true_mean)
        
        df_boot = pd.DataFrame(boot_rows, columns=gene_cols)
        
        # Reconstruct metadata for bootstrap samples
        meta_rows = []
        for donor in donors:
            df_donor = df[df[donor_col] == donor]
            d_age = df_donor[age_col].iloc[0]
            d_type = df_donor[cell_type_col].iloc[0]
            for _ in range(100):
                meta_rows.append({
                    donor_col: donor,
                    age_col: d_age,
                    cell_type_col: d_type
                })
        
        meta_df = pd.DataFrame(meta_rows)
        df_final = pd.concat([meta_df, df_boot], axis=1)
        return df_final
    
    def _predict_core(self, data, impute_data, model, sample_type):
        """
        Core prediction logic using matrix multiplication.
        Robust to case sensitivity and intercept-only models.
        """
        # Force uppercase gene names for consistency
        data.columns = data.columns.astype(str).str.upper().str.strip()
        
        model = model.copy()
        model['feature_name'] = model['feature_name'].astype(str).str.upper().str.strip()
        
        # 1. Extract and separate Intercept
        if "INTERCEPT" in model['feature_name'].values:
            intercept_val = model.loc[model['feature_name'] == 'INTERCEPT', 'coefficient'].values[0]
        else:
            intercept_val = 0.0
            
        coef_df = model[model['feature_name'] != 'INTERCEPT']
        model_genes = coef_df['feature_name'].values
        
        # 2. Check Gene Matching
        if len(model_genes) > 0:
            common_genes = set(data.columns) & set(model_genes)
            if len(common_genes) == 0:
                print(f"[CRITICAL WARNING] 0 Genes matched! Prediction will reflect intercept only.")
        
        # 3. Align Data (Reindex)
        expr_data = data.reindex(columns=model_genes)
        
        # 4. Imputation
        if len(model_genes) > 0:
            missing_genes = expr_data.columns[expr_data.isna().any()].tolist()
            if missing_genes and impute_data is not None:
                impute_data = impute_data.copy()
                impute_data['feature_name'] = impute_data['feature_name'].astype(str).str.upper().str.strip()
                impute_dict = dict(zip(impute_data['feature_name'], impute_data['imputation_value']))
                
                fill_values = {g: impute_dict.get(g, 0.0) for g in missing_genes}
                expr_data = expr_data.fillna(fill_values)
        
        # Fill remaining NaNs with 0
        expr_data = expr_data.fillna(0.0)
        
        # 5. Calculation
        if len(model_genes) > 0:
            coef_series = coef_df.set_index('feature_name').loc[model_genes, 'coefficient']
            raw_preds = expr_data.dot(coef_series.values) + intercept_val
        else:
            # Fallback for Intercept-only models
            raw_preds = np.full(len(data), intercept_val)
        
        return pd.DataFrame({
            "predictions": raw_preds,
            "sample_type": sample_type
        })


class PASTA_Clock:
    """
    PASTA Transcriptomic Clock (Python Implementation).
    
    This class implements the "Propensity Adjustment by Subsampling Transcriptomic Age" (PASTA)
    clock and its variants. These clocks predict biological age based on **bulk RNA-seq** data.

    

    **Available Models:**
    1.  **PASTA**: The standard model using rank-based scoring. Robust to normalization differences.
    2.  **REG**: A regularized version (Lasso/ElasticNet) for potentially higher precision.
    3.  **CT46**: A simplified version using a core set of 46 transcripts.

    Attributes:
        model_type (str): The variant being used ("PASTA", "REG", "CT46").
        metadata (dict): Model metadata (species, tissue, source).
        coefs (pd.Series): Model coefficients (Gene weights).
        intercept (float): Model intercept.
        scaling_factor (float): Final scaling factor applied to the score.

    Reference:
        Salignon, J. et al. Pasta, a versatile transcriptomic clock, maps the chemical and genetic determinants of aging and rejuvenation.
        bioRxiv (2025). https://doi.org/10.1101/2025.06.04.657785
    """
    METADATA_VARIANTS = {
        "PASTA": {
            "year": 2025,
            "species": "Human",
            "tissue": "Multi-tissue",
            "omic type": "Transcriptomics (Bulk RNA-seq)",
            "prediction": "Age Score",
            "source": "https://doi.org/10.1101/2025.06.04.657785",
            "description": "Propensity Adjustment by Subsampling Transcriptomic Age (Standard Model)."
        },
        "REG": {
            "year": 2025,
            "species": "Human",
            "tissue": "Multi-tissue",
            "omic type": "Transcriptomics (Bulk RNA-seq)",
            "prediction": "Chronological Age (Years)",
            "source": "https://doi.org/10.1101/2025.06.04.657785",
            "description": "Regularized version of the PASTA clock."
        },
        "CT46": {
            "year": 2025,
            "species": "Human",
            "tissue": "Multi-tissue",
            "omic type": "Transcriptomics (Bulk RNA-seq)",
            "prediction": "Age Score",
            "source": "https://doi.org/10.1101/2025.06.04.657785"
        }
    }

    
    def __init__(self, model_type: str = "PASTA", data_dir: str = None):
        """
        Initialize the PASTA Clock.
        
        Args:
            model_type: "PASTA", "REG", or "CT46".
            data_dir: Path to directory containing coefficients and metadata.
        """
        self.model_type = model_type
        
        if model_type in self.METADATA_VARIANTS:
            self.metadata = self.METADATA_VARIANTS[model_type]
            self.name = model_type 
        else:
            raise ValueError(f"Invalid model_type. Choose from {list(self.METADATA_VARIANTS.keys())}")
            
        self.coefs = None
        self.intercept = 0.0
        self.scaling_factor = 1.0
        
        # 1. Path Resolution
        if data_dir is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data", "PASTA")
            
        file_map = {
            "PASTA": "PASTA_PASTA",
            "REG":   "PASTA_REG",
            "CT46":  "PASTA_CT46"
        }
        
        if model_type not in file_map:
            raise ValueError(f"Invalid model_type. Choose from {list(file_map.keys())}")
            
        file_base = file_map[model_type]
        coef_path = os.path.join(data_dir, f"{file_base}_coefs.csv")
        meta_path = os.path.join(data_dir, f"{file_base}_meta.csv")

        # 2. Load Coefficients
        if os.path.exists(coef_path):
            self._load_coefficients(coef_path)
        else:
            print(f"[Error] Coefficient file not found: {coef_path}")

        # 3. Load Metadata (Intercept & Scaling)
        if os.path.exists(meta_path):
            self._load_metadata(meta_path)
        else:
            print(f"[Warning] Metadata file not found: {meta_path}. Using default intercept=0, scale=1.")

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the model coefficients including the intercept.
        
        Returns:
            pd.DataFrame: columns ['model_name', 'gene', 'coef']
        """
        if self.coefs is None:
            return pd.DataFrame(columns=['model_name', 'gene', 'coef'])

        # 1. Extract gene weights
        df_genes = pd.DataFrame({
            'gene': self.coefs.index,
            'coef': self.coefs.values
        })

        # 2. Extract the intercept 
        df_intercept = pd.DataFrame([{
            'gene': '(Intercept)',
            'coef': self.intercept
        }])

        # 3. merge
        final_df = pd.concat([df_intercept, df_genes], ignore_index=True)

        final_df['model_name'] = self.name

        return final_df[['model_name', 'gene', 'coef']]
    
    def info(self):
        """Prints metadata information about the PASTA clock."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def _load_coefficients(self, path):
        """Helper to load and standardize coefficients."""
        try:
            df = pd.read_csv(path)
            df.columns = df.columns.str.lower().str.strip()
            
            # Extract gene and coefficient
            if "gene" in df.columns and "coefficient" in df.columns:
                self.coefs = df.set_index("gene")["coefficient"]
            else:
                # Fallback: Guess columns by position
                print("[Warning] 'gene'/'coefficient' columns not found. Guessing by position.")
                obj_cols = df.select_dtypes(include=['object']).columns
                num_cols = df.select_dtypes(include=['number']).columns
                if len(obj_cols) > 0 and len(num_cols) > 0:
                    self.coefs = df.set_index(obj_cols[0])[num_cols[-1]]
                else:
                    raise ValueError("Cannot parse coefficient file structure.")
            
            self.coefs = self.coefs.astype(float)
            
        except Exception as e:
            print(f"[Error] Failed to load coefficients: {e}")

    def _load_metadata(self, path):
        """Helper to load intercept and scaling factor."""
        try:
            df = pd.read_csv(path)
            if "intercept" in df.columns:
                self.intercept = float(df["intercept"].values[0])
            
            if "scaling_factor" in df.columns:
                self.scaling_factor = float(df["scaling_factor"].values[0])
                
        except Exception as e:
            print(f"[Error] Failed to load metadata: {e}")

    def predict(self, adata: anndata.AnnData, rank_norm: bool = True) -> pd.DataFrame:
        """
        Calculates the PASTA Age Score.

        The prediction pipeline includes:
        1.  **Feature Alignment**: Reindexes input data to match model genes.
        2.  **Imputation**: Fills missing genes with the global sample median.
        3.  **Rank Normalization** (Crucial): Converts expression to ranks to handle batch effects.
        4.  **Scoring**: Computes weighted sum and applies scaling.

        Args:
            adata (anndata.AnnData): Input expression object.
                - **.var_names**: Must be **Gene Symbols** (e.g., 'GAPDH', 'ACTB'). 
                  *Warning*: If Ensembl IDs are provided, feature alignment will fail.
                - **.X**: Can be Raw Counts, TPM, or FPKM. Rank normalization handles the scale.
            rank_norm (bool, optional): Whether to apply rank normalization. 
                Defaults to True (Strongly Recommended for PASTA).

        Returns:
            pd.DataFrame: A DataFrame containing the predicted age scores, indexed by sample ID.

        """
        if self.coefs is None:
            return pd.DataFrame()

        # --- Step 1: Prepare Expression Matrix ---
        if scipy.sparse.issparse(adata.X):
            X = adata.X.toarray()
        else:
            X = adata.X
        
        # Samples as rows, Genes as columns
        exp_df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)

        # --- Step 2: Gene Filtering & Imputation ---
        model_genes = self.coefs.index
        
        # Reindex to match model genes (introduces NaNs for missing genes)
        exp_df = exp_df.reindex(columns=model_genes)
        
        # Impute missing values with global median
        global_median = exp_df.stack().median()
        exp_df = exp_df.fillna(global_median)

        # --- Step 3: Rank Normalization ---
        # Normalize expression ranks within each sample (axis=1)
        if rank_norm:
            exp_df = exp_df.rank(axis=1, method='average')

        # --- Step 4: Calculate Score ---
        exp_df = exp_df.astype(float)
        
        # Linear dot product
        raw_score = exp_df.dot(self.coefs)
        
        # Add intercept
        score = raw_score + self.intercept
        
        # --- Step 5: Scaling ---
        score = score * self.scaling_factor
        
        return score.to_frame(name=f"Predicted_{self.model_type}_Age")