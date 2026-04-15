import pandas as pd
import numpy as np
import os
from typing import List, Optional, Dict, Tuple, Union
from ..utils import get_data_path
import scipy.sparse
import scanpy as sc
import warnings
import anndata



def generate_pseudocells(
    adata: anndata.AnnData, 
    size: int = 15, 
    n_repeats: int = 100, 
    aggr_method: str = "mean"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generates pseudocells by randomly sampling and aggregating single cells from each donor.
    
    Args:
        adata (AnnData): Input single-cell object. Must contain 'donor_id' and 'age' in .obs.
        size (int): Number of cells to aggregate per pseudocell.
        n_repeats (int): Number of pseudocells to generate per donor.
        aggr_method (str): "mean" (for scImmuAging) or "sum" (for Buckley SVZ clocks).
        
    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: X_pseudo (expression), meta_pseudo (metadata).
    """
    if aggr_method not in ["mean", "sum"]:
        raise ValueError("aggr_method must be 'mean' or 'sum'")

    pseudocells = []
    metadata = []
    
    if "donor_id" not in adata.obs.columns or "age" not in adata.obs.columns:
        raise ValueError("AnnData.obs must contain 'donor_id' and 'age' columns.")

    grouped = adata.obs.groupby(['donor_id', 'age'])
    
    for (donor, age), group_indices in grouped.indices.items():
        n_cells = len(group_indices)
        replace = n_cells < size
        
        # Optimize handling for sparse vs dense matrices
        if hasattr(adata.X, "toarray"):
            donor_X = adata.X[group_indices].toarray()
        else:
            donor_X = adata.X[group_indices]
            
        donor_df = pd.DataFrame(donor_X, columns=adata.var_names)
        
        donor_pseudos = []
        for _ in range(n_repeats):
            sample_idx = np.random.choice(n_cells, size=size, replace=replace)
            
            # Use requested aggregation method
            if aggr_method == "sum":
                pseudo_expr = donor_df.iloc[sample_idx].sum(axis=0)
            else:
                pseudo_expr = donor_df.iloc[sample_idx].mean(axis=0)
                
            donor_pseudos.append(pseudo_expr)
        
        pseudocells.extend(donor_pseudos)
        metadata.extend([{'donor_id': donor, 'age': age}] * n_repeats)
        
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
            X_pseudo, meta_pseudo = generate_pseudocells(
                adata_sub,
                size=15,
                n_repeats=100,
                aggr_method="mean"
                )

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




class BuckleyMouseSVZ:
    """
    Implements the SVZ Bootstrap Clocks (Buckley et al. 2023).
    
    Estimates chronological or biological age from scRNA-seq data derived 
    from the mouse subventricular zone (SVZ).
    
    Mechanism:
    1. Preprocessing: Generates Pseudocells via bootstrapping (aggrMethod = "sum").
    2. Normalization: CP10k + log1p on pseudocell counts.
    3. Prediction: Matrix multiplication with specific Elastic Net weights.
    4. Biological Post-processing: Transforms proliferation fractions to age.
    
    References:
        Buckley MT, et al. Cell-type-specific aging clocks to quantify aging and 
        rejuvenation in neurogenic regions of the brain. Nat Aging (2023).
    """
    METADATA = {
        "year": 2023,
        "species": "Mouse",
        "tissue": "Subventricular Zone (SVZ)",
        "omic type": "Transcriptomics(scRNA-seq)",
        "prediction": "Chronological Age(Months) or Biological Age",
        "source": "https://doi.org/10.1038/s43587-022-00335-4"
    }

    def __init__(self):
        self.name = "BuckleyMouseSVZ"
        self.metadata = self.METADATA
        self.data_dir = get_data_path("BuckleyMouseSVZ")
        
        # Structure: self.models["chronological"]["Microglia"] = {'intercept': 0.0, 'weights': pd.Series}
        self.models = {
            "chronological": {},
            "biological": {}
        }
        
        if not os.path.exists(self.data_dir):
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")

        # Example expected file naming: "aNSC_NPC_Biological_Bootstrap.csv" or "Neuroblast_Biological_Bootstrap.csv"
        for filename in os.listdir(self.data_dir):
            if filename.endswith(".csv"):
                filename_clean = filename.replace(".csv", "")
                
                if "Chronological_Bootstrap" in filename_clean:
                    clock_type = "chronological"
                    ct_name = filename_clean.replace("_Chronological_Bootstrap", "")
                elif "Biological_Bootstrap" in filename_clean:
                    clock_type = "biological"
                    ct_name = filename_clean.replace("_Biological_Bootstrap", "")
                else:
                    continue 
                    
                if clock_type in self.models:
                    try:
                        df = pd.read_csv(os.path.join(self.data_dir, filename))
                        
                        intercept_mask = df['gene'].str.lower().isin(['intercept', '(intercept)'])
                        if intercept_mask.any():
                            intercept = float(df[intercept_mask]['coef'].iloc[0])
                            weights = df[~intercept_mask].set_index('gene')['coef']
                        else:
                            intercept = 0.0
                            weights = df.set_index('gene')['coef']
                            
                        self.models[clock_type][ct_name] = {
                            "intercept": intercept,
                            "weights": weights
                        }
                    except Exception as e:
                        print(f"[{self.name}] Warning: Failed to load {filename}: {e}")
        

    def info(self):
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            print(f"  - {key.replace('_', ' ').title()}: {value}")

    def predict(
        self, 
        adata: anndata.AnnData, 
        cell_types: List[str] = ["aNSC_NPC", "Astrocyte_qNSC", "Endothelial", "Microglia", "Neuroblast", "Oligodendro"], 
        clock_type: str = "chronological",
        verbose: bool = True
    ) -> Dict[str, Dict[str, pd.DataFrame]]:
        """
        Execute age prediction strictly following the Buckley et al. (2023) methodology.

        This method performs internal data orchestration, including bootstrap-based 
        pseudocell generation and specific CP10k + log1p normalization.

        Note:
            - Input data (adata.X) MUST contain RAW COUNTS. 
            - Biological age is derived from proliferation fractions and transformed 
              via: Age = 35 - (100 * Predicted_Fraction).

        Required metadata (adata.obs):
            - 'celltype': Must match the names in the 'cell_types' list.
            - 'donor_id': Used for aggregating bootstrap results per donor.
            - 'age': Reference chronological age.

        Args:
            adata: AnnData object with raw expression counts and required obs columns.
            cell_types: List of cell types to evaluate. Defaults to the six major 
                types identified in the SVZ neurogenic niche.
            clock_type: Type of prediction to run:
                - 'chronological': Predicts age in months.
                - 'biological': Predicts a biological score (proliferation-based).
            verbose: If True, prints detailed progress logs for each processing step.

        Returns:
            Dict: A nested dictionary where keys are cell types. Each value contains:
                - 'BootstrapCell': Raw predictions for each generated pseudocell.
                - 'Donor': Aggregated mean predictions per donor/sample.
        """
        results = {}
        clock_type = clock_type.lower()
        
        if clock_type not in ["chronological", "biological"]:
            raise ValueError("clock_type must be 'chronological' or 'biological'")
            
        if "celltype" not in adata.obs.columns:
            raise ValueError("AnnData.obs must contain 'celltype' column.")

        for ct in cell_types:
            if ct not in self.models.get(clock_type, {}):
                if verbose: print(f"[{self.name}] Warning: '{ct}' model not found for {clock_type} clock. Skipping.")
                continue
                
            if verbose: print(f"\n--- Processing {clock_type} clock for: {ct} ---")
            
            # 1. Subset Data
            subset_idx = adata.obs['celltype'] == ct
            if subset_idx.sum() == 0:
                if verbose: print(f"No cells found for {ct}.")
                continue
            adata_sub = adata[subset_idx]
            
            # 2. Bootstrap Pseudocells (AggrMethod = "sum")
            if verbose: print("Generating pseudocells (Summation)...")
            X_pseudo, meta_pseudo = generate_pseudocells(adata_sub, size=15, n_repeats=100, aggr_method="sum")
            
            if X_pseudo.empty:
                continue

            # 3. Normalization: CP10k + log1p
            if verbose: print("Applying CP10k + log1p normalization...")
            row_sums = X_pseudo.sum(axis=1)
            row_sums[row_sums == 0] = 1 # Prevent division by zero
            
            # X / row_sums * 10000 -> log1p
            X_normed = X_pseudo.div(row_sums, axis=0) * 10000
            X_logged = np.log1p(X_normed)

            # 4. Predict
            if verbose: print("Predicting age...")
            model_params = self.models[clock_type][ct]
            weights = model_params['weights']
            intercept = model_params['intercept']
            
            # Feature alignment
            X_aligned = X_logged.reindex(columns=weights.index, fill_value=0)
            raw_preds = X_aligned.dot(weights) + intercept
            
            # 5. Biological Age Post-treatment
            if clock_type == "biological":
                # 35 - (100 * fraction)
                raw_preds = 35 - (100 * raw_preds)
                
            bootstrap_res = meta_pseudo.copy()
            bootstrap_res['Prediction'] = raw_preds.values
            bootstrap_res['CellType'] = ct
            
            # 6. Aggregate per Donor
            donor_res = bootstrap_res.groupby('donor_id').agg({
                'age': 'first',
                'Prediction': 'mean',
                'CellType': 'first'
            }).rename(columns={'Prediction': 'predicted'}).reset_index()
            
            # Round chronological age to align with R implementation output style
            if clock_type == "chronological":
                donor_res['predicted'] = donor_res['predicted'].round(2)
                
            results[ct] = {
                "BootstrapCell": bootstrap_res,
                "Donor": donor_res
            }
            
            if verbose: print(f"Prediction complete for {ct}!")
            
        return results



class BaseSingleCellClock:
    """
    Base class for single-cell transcriptomic aging clocks.
    
    Provides a standardized framework for data loading, cross-validation handling, 
    ensemble prediction, and complex sampling methods (Single-cell, Pseudobulk, 
    and Bootstrap).
    """
    
    def __init__(self, name: str, metadata: dict, data_dir: str, 
                 coef_filename: str, impute_filename: str,
                 bootstrap_thresholds: Union[int, dict] = 50,
                 bootstrap_replace: bool = True):
        """
        Initialize the base aging clock model.
        
        Args:
            name: Model identifier (e.g., 'BrainCTClock').
            metadata: Dictionary containing model provenance and specifications.
            data_dir: Path to the directory containing weights and imputation data.
            coef_filename: Filename for the model coefficients (CSV).
            impute_filename: Filename for the feature imputation values (CSV).
            bootstrap_thresholds: Cell count threshold for bootstrap sampling. 
                Can be an integer or a dictionary mapping cell types to thresholds.
            bootstrap_replace: Whether to sample with replacement during bootstrapping.
        """
        self.name = name
        self.metadata = metadata
        self.models = {} 
        self.imputation_data = {}
        
        self._bootstrap_thresholds = bootstrap_thresholds
        self._bootstrap_replace = bootstrap_replace
        
        if data_dir and os.path.exists(data_dir):
            self._load_models(data_dir, coef_filename, impute_filename)
        else:
            print(f"[Error] {self.name} data directory not found: {data_dir}")

    def info(self):
        """Print formatted model metadata to the console."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def get_coefs(self) -> pd.DataFrame:
        """
        Consolidate all model coefficients into a single DataFrame.

        Returns:
            pd.DataFrame: A long-format DataFrame with columns 
                ['model_name', 'gene', 'coef'].
        """
        if not self.models:
            return pd.DataFrame(columns=['model_name', 'gene', 'coef'])

        results = []
        for model_key, folds_list in self.models.items():
            if not folds_list: continue
            for i, fold_df in enumerate(folds_list):
                df = fold_df.copy()
                df = df.rename(columns={"feature_name": "gene", "coefficient": "coef"})
                df["model_name"] = f"{self.name}_{model_key}_fold{i+1}"
                if all(col in df.columns for col in ['model_name', 'gene', 'coef']):
                    results.append(df[['model_name', 'gene', 'coef']])

        return pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=['model_name', 'gene', 'coef'])

    def _load_models(self, data_dir: str, coef_filename: str, impute_filename: str):
        """Internal: Load coefficients and imputation parameters from CSV files."""
        coef_path = os.path.join(data_dir, coef_filename)
        impute_path = os.path.join(data_dir, impute_filename)

        # Load imputation values
        if os.path.exists(impute_path):
            df_imp_all = pd.read_csv(impute_path)
            if "model_key" in df_imp_all.columns:
                for model_key, group in df_imp_all.groupby("model_key"):
                    self.imputation_data[model_key] = group[["feature_name", "imputation_value"]].reset_index(drop=True)
        else:
            print(f"[Warning] Imputation file not found for {self.name}: {impute_path}")

        # Load cross-validation model coefficients
        if os.path.exists(coef_path):
            df_coef_all = pd.read_csv(coef_path)
            if "model_key" in df_coef_all.columns and "fold" in df_coef_all.columns:
                for model_key, model_group in df_coef_all.groupby("model_key"):
                    fold_list = []
                    unique_folds = model_group["fold"].unique()
                    # Sort folds numerically (handles 'fold_1' vs 'fold_10' correctly)
                    try:
                        sorted_folds = sorted(unique_folds, key=lambda x: int(str(x).split('_')[-1]) if '_' in str(x) else x)
                    except (ValueError, IndexError):
                        sorted_folds = sorted(unique_folds)

                    for fold in sorted_folds:
                        fold_df = model_group[model_group["fold"] == fold]
                        clean_df = fold_df[["feature_name", "coefficient"]].reset_index(drop=True)
                        fold_list.append(clean_df)
                    
                    self.models[model_key] = fold_list
        else:
            print(f"[Error] Coefficient file not found for {self.name}: {coef_path}")

    def predict(self, adata: anndata.AnnData, cell_types: List[str], model_name: Union[str, List[str]] = "all") -> Dict[str, pd.DataFrame]:
        """
        Execute the prediction pipeline for designated model modes.

        Args:
            adata: AnnData object containing expression data and metadata.
            cell_types: List of cell types to process.
            model_name: Prediction mode(s). Options: 'SC', 'Pseudobulk', 'Bootstrap', or 'all'.

        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping mode names to result DataFrames.
        """
        donor_col, age_col, cell_type_col = "donor_id", "age", "celltype"

        # Validate metadata columns
        missing = [col for col in [donor_col, age_col, cell_type_col] if col not in adata.obs.columns]
        if missing:
            raise ValueError(f"Input AnnData is missing required metadata columns: {missing}.")

        valid_models = ["SC", "Pseudobulk", "Bootstrap"]
        models_to_run = valid_models if model_name == "all" else ([model_name] if isinstance(model_name, str) else model_name)

        # Validate requested modes
        invalid = [m for m in models_to_run if m not in valid_models]
        if invalid:
            raise ValueError(f"Invalid model_name(s): {invalid}. Valid options: {valid_models}")

        all_results = {}
        print(f"Starting {self.name} prediction pipeline...")

        for current_model_type in models_to_run:
            print(f"Running mode: {current_model_type}")
            res = self._run_prediction_pipeline(
                sample_type=current_model_type, adata=adata, common_celltypes=cell_types,
                donor_col=donor_col, age_col=age_col, cell_type_col=cell_type_col
            )
            all_results[current_model_type] = res

        print(f"--- All {self.name} predictions complete! ---")
        return all_results

    def _run_prediction_pipeline(self, sample_type, adata, common_celltypes, donor_col, age_col, cell_type_col):
        """Internal: Handles CV ensemble averaging and per-celltype processing."""
        final_results_list = []

        for cell_type in common_celltypes:
            # Format key to match stored CSV keys (replacing spaces with underscores)
            clock_key = f"{sample_type}_{cell_type.replace(' ', '_')}"
            
            if clock_key not in self.models:
                print(f"Skipping {cell_type}: Model key {clock_key} not found.")
                continue
            
            model_folds = self.models[clock_key]
            impute_df = self.imputation_data.get(clock_key, None)

            # Data extraction and aggregation
            df_base = self._prepare_data(adata, cell_type, sample_type, donor_col, age_col, cell_type_col)
            if df_base.empty:
                continue

            # Accumulate predictions across all CV folds
            fold_preds = []
            for model_fold in model_folds:
                pred_res = self._predict_core(
                    data=df_base.drop(columns=[donor_col, age_col, cell_type_col], errors='ignore'), 
                    impute_data=impute_df, model=model_fold, sample_type=sample_type
                )
                fold_preds.append(pred_res['predictions'].values)

            # Average fold predictions (Ensemble)
            avg_preds = np.column_stack(fold_preds).mean(axis=1)

            res_df = pd.DataFrame({
                "predictions": avg_preds,
                "ages": df_base[age_col].values,      
                "donors": df_base[donor_col].values,  
                "sample_type": sample_type,
                "celltype": cell_type
            })
            final_results_list.append(res_df)

        return pd.concat(final_results_list, ignore_index=True) if final_results_list else pd.DataFrame()

    def _prepare_data(self, adata, cell_type, sample_type, donor_col, age_col, cell_type_col):
        """Internal: Extracts expression and applies sampling/aggregation logic."""
        if cell_type_col not in adata.obs: return pd.DataFrame()
             
        subset = adata[adata.obs[cell_type_col] == cell_type].copy()
        if subset.n_obs == 0: return pd.DataFrame()

        # Handle sparse vs dense matrices
        X = subset.X.toarray() if hasattr(subset.X, "toarray") else subset.X
        df = pd.DataFrame(X, index=subset.obs_names, columns=subset.var_names)
        
        # Attach metadata for grouping
        for col in [donor_col, age_col, cell_type_col]:
            df[col] = subset.obs[col].values

        gene_cols = subset.var_names.tolist()

        if sample_type == "SC": 
            return df
        elif sample_type == "Pseudobulk":
            # Per-donor mean aggregation
            return df.groupby([donor_col, age_col, cell_type_col], observed=True)[gene_cols].mean().reset_index()
        elif sample_type == "Bootstrap":
            return self._bootstrap_sampling(df, cell_type, gene_cols, donor_col, age_col, cell_type_col)
        
        return pd.DataFrame()

    def _bootstrap_sampling(self, df, cell_type, gene_cols, donor_col, age_col, cell_type_col):
        """Internal: Performs donor-level bootstrap sampling to simulate ensembled cell states."""
        # Resolve threshold for specific cell type
        if isinstance(self._bootstrap_thresholds, dict):
            threshold = self._bootstrap_thresholds.get(cell_type, 50)
        else:
            threshold = self._bootstrap_thresholds

        donors = df[donor_col].unique()
        boot_rows, meta_rows = [], []
        np.random.seed(42)

        for donor in donors:
            df_donor = df[df[donor_col] == donor]
            n_cells = len(df_donor)
            expr_matrix = df_donor[gene_cols].values
            d_age = df_donor[age_col].iloc[0]
            d_type = df_donor[cell_type_col].iloc[0]

            # Generate 100 bootstrap samples per donor
            if n_cells > threshold:
                for _ in range(100):
                    indices = np.random.choice(n_cells, size=threshold, replace=self._bootstrap_replace)
                    boot_rows.append(expr_matrix[indices, :].mean(axis=0))
            else:
                # If cells < threshold, fall back to true mean
                true_mean = expr_matrix.mean(axis=0)
                boot_rows.extend([true_mean] * 100)
            
            meta_rows.extend([{donor_col: donor, age_col: d_age, cell_type_col: d_type}] * 100)
        
        return pd.concat([pd.DataFrame(meta_rows), pd.DataFrame(boot_rows, columns=gene_cols)], axis=1)

    def _predict_core(self, data, impute_data, model, sample_type):
        """Internal: Core matrix multiplication logic with feature alignment and imputation."""
        # Standardize feature names to uppercase/stripped strings
        data.columns = data.columns.astype(str).str.upper().str.strip()
        model = model.copy()
        model['feature_name'] = model['feature_name'].astype(str).str.upper().str.strip()
        
        # Extract Intercept
        if "INTERCEPT" in model['feature_name'].values or "(INTERCEPT)" in model['feature_name'].values:
            intercept_idx = model['feature_name'].isin(['INTERCEPT', '(INTERCEPT)'])
            intercept_val = model.loc[intercept_idx, 'coefficient'].values[0]
            coef_df = model[~intercept_idx]
        else:
            intercept_val = 0.0
            coef_df = model
            
        model_genes = coef_df['feature_name'].values
        
        # Alignment check
        if len(model_genes) > 0:
            common_genes = set(data.columns) & set(model_genes)
            if len(common_genes) == 0:
                print(f"[Warning] Zero genes matched for model. Prediction relies solely on intercept.")
                
        # Reindex data to match model features
        expr_data = data.reindex(columns=model_genes)
        
        # Apply imputation for missing genes
        if len(model_genes) > 0:
            missing_genes = expr_data.columns[expr_data.isna().any()].tolist()
            if missing_genes and impute_data is not None:
                impute_data = impute_data.copy()
                impute_data['feature_name'] = impute_data['feature_name'].astype(str).str.upper().str.strip()
                impute_dict = dict(zip(impute_data['feature_name'], impute_data['imputation_value']))
                fill_values = {g: impute_dict.get(g, 0.0) for g in missing_genes}
                expr_data = expr_data.fillna(fill_values)
        
        # Final fallback for any remaining NaNs
        expr_data = expr_data.fillna(0.0)
        
        # Matrix multiplication: y = Xb + intercept
        if len(model_genes) > 0:
            coef_series = coef_df.set_index('feature_name').loc[model_genes, 'coefficient']
            raw_preds = expr_data.dot(coef_series.values) + intercept_val
        else:
            raw_preds = np.full(len(data), intercept_val)
        
        return pd.DataFrame({"predictions": raw_preds, "sample_type": sample_type})


class BrainCTClock(BaseSingleCellClock):
    """
    Human Brain Cell-Type-Specific Aging Clocks (Muralidharan et al. 2025).
    """
    METADATA = {
        "year": 2025,
        "species": "Human",
        "tissue": "Brain",
        "omic type": "Transcriptomics (scRNA-seq)",
        "prediction": "Chronological Age(Years)",
        "source": "https://doi.org/10.1002/advs.202506109"
    }

    BOOTSTRAP_THRESHOLDS = {
        "Oligodendrocytes": 200, "Astrocytes": 50, "Microglia": 50,
        "OPCs": 50, "Excitatory Neurons": 100, "Inhibitory Neurons": 100
    }

    def __init__(self, data_dir: str = None):
        if data_dir is None:
            current_file_path = os.path.abspath(__file__)
            models_dir = os.path.dirname(current_file_path)
            package_root = os.path.dirname(models_dir)
            data_dir = os.path.join(package_root, "data", "Brain_CT_Clock")
            
        super().__init__(
            name="BrainCTClock",
            metadata=self.METADATA,
            data_dir=data_dir,
            coef_filename="brain_ct_coefs.csv",
            impute_filename="brain_ct_imputation.csv",
            bootstrap_thresholds=self.BOOTSTRAP_THRESHOLDS,
            bootstrap_replace=False 
        )


class ScAgePolyakClock(BaseSingleCellClock):
    """
    Human PBMC Cell-Type-Specific Aging Clocks (Zakar-Polyák et al. 2024).
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood (PBMC)",
        "omic type": "Transcriptomics (scRNA-seq)",
        "prediction": "Chronological Age (Years)",
        "source": "https://doi.org/10.1038/s42003-024-07094-5"
    }

    def __init__(self, data_dir: str = None):
        if data_dir is None:
            current_file_path = os.path.abspath(__file__)
            models_dir = os.path.dirname(current_file_path)
            package_root = os.path.dirname(models_dir)
            data_dir = os.path.join(package_root, "data", "ScAge_Polyak")
            
        super().__init__(
            name="ScAgePolyakClock",
            metadata=self.METADATA,
            data_dir=data_dir,
            coef_filename="scage_polyak_coefs.csv", 
            impute_filename="scage_polyak_imputation.csv"
        )

    def predict(self, adata: anndata.AnnData, cell_types: List[str], model_name: Union[str, List[str]] = "SC") -> Dict[str, pd.DataFrame]:
        """
        Predict cell-type-specific transcriptomic age using the ScAgePolyak model.

        Note:
            This clock is optimized and strictly validated for Single-Cell (SC) mode. 
            Other modes like 'Pseudobulk' or 'Bootstrap' are disabled for this specific class.

        Required input structure (adata.obs):
            - 'donor_id': Unique identifier for each sample/donor.
            - 'age': Chronological age of the donor.
            - 'celltype': Cell type annotations matching the 'cell_types' list.
        
        Required input data (adata.X):
            - Log-normalized gene expression (typically ln(CP10k + 1) or as specified 
              in Zakar-Polyak et al. 2024).

        Args:
            adata: AnnData object containing single-cell gene expression and metadata.
            cell_types: A list of cell type names to process.
            model_name: Forced to 'SC' (Single-Cell).

        Returns:
            A dictionary where each key is a cell type and the value is a DataFrame 
            containing 'predictions', 'ages', and 'donors'.
        """
        if model_name != "SC" and model_name != ["SC"]:
            print(
                f"[Info] {self.name} is specifically designed for single-cell level "
                f"calculations. Defaulting to 'SC' mode."
            )
        
        return super().predict(adata, cell_types, model_name="SC")



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