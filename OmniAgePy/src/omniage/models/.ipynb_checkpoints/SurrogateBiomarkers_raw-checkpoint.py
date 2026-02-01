import pandas as pd
import numpy as np
import os
import glob
from ..utils import load_clock_coefs
from .base import BaseLinearClock

class CompCRP:
    """
    DNAm-based C-Reactive Protein (CRP) Score (Teschendorff et al.).
    
    Computes a DNA methylation surrogate proxy for CRP levels. 
    The calculation is based on the correlation between normalized methylation Z-scores 
    and the sign of pre-determined weights (rather than a direct linear combination).
    
    Includes two signatures:
    - CRP: Unadjusted signature (based on 1765 CpGs).
    - intCRP: Adjusted for immune cell composition (based on 62 CpGs).
    """
    
    def __init__(self):
        self.signatures = {}
        
        # Load both models
        for model_name in ["CRP", "intCRP"]:
            try:
                # Expected filename: CompCRP_CRP.csv
                df = load_clock_coefs(f"CompCRP_{model_name}")
                
                # Extract 'probe' and 'coef' columns
                # While the core calculation uses sign(coef), we store the raw 
                # coefficients for completeness and potential future use.
                self.signatures[model_name] = df.set_index("probe")["coef"]
            except Exception as e:
                print(f"[Error] Failed to load CompCRP {model_name}: {e}")

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves coefficients for both CRP and intCRP models.
        Returns a DataFrame with columns: ['model_name', 'probe', 'coef'].
        """
        dfs = []
        for model_name, series in self.signatures.items():
            df = series.reset_index()
            df.columns = ['probe', 'coef']
            
            df['model_name'] = model_name
            dfs.append(df)
            
        if not dfs:
            return pd.DataFrame(columns=['model_name', 'probe', 'coef'])
            
        return pd.concat(dfs, ignore_index=True)

    def predict(self, beta_df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        """
        Calculates CRP scores.
        
        Implementation logic follows the original R methodology:
        1. Subset to common CpGs available in the input data.
        2. Perform Row-wise Z-score normalization (standardize each CpG across all samples).
        3. Calculate the Pearson correlation between each sample's Z-scores and 
           the sign vector of the model coefficients.
        
        Args:
            beta_df: Methylation data (Rows=CpGs, Columns=Samples).
            verbose: If True, prints CpG coverage statistics.
            
        Returns:
            pd.DataFrame: A DataFrame (Rows=Samples, Columns=Scores) containing 
                          'CompCRP_CRP' and 'CompCRP_intCRP'.
        """
        results = []
        
        for name, coefs in self.signatures.items():
            # 1. Intersect CpGs
            common_cpgs = beta_df.index.intersection(coefs.index)
            
            if verbose:
                print(f"[CompCRP] {name}: Found {len(common_cpgs)} / {len(coefs)} CpGs")
                
            if len(common_cpgs) < 2:
                print(f"[Warning] Too few CpGs for {name} ({len(common_cpgs)}). Skipping.")
                continue
                
            # Subset data and coefficients
            # R equivalent: tmp.m <- nbeta.m[map2.idx,]
            X_subset = beta_df.loc[common_cpgs]
            w_subset = coefs.loc[common_cpgs]
            
            # 2. Row-wise Z-score Normalization
            # R equivalent: z.m <- (tmp.m - rowMeans(tmp.m)) / apply(tmp.m, 1, sd)
            # In Pandas, calculate stats along axis=1 (columns/samples), then broadcast along axis=0 (rows)
            row_means = X_subset.mean(axis=1)
            row_stds = X_subset.std(axis=1)
            
            # Handle division by zero (for constant rows)
            row_stds[row_stds == 0] = 1.0 
            
            # Compute Z-score matrix
            Z_matrix = X_subset.sub(row_means, axis=0).div(row_stds, axis=0)
            
            # 3. Calculate Correlation with Sign of Coefficients
            # R equivalent: score.lv[[i]] <- as.vector(cor(z.m, sign(sigCRP.v[map1.idx])))
            # We calculate the correlation of each sample (column) in Z_matrix with the target sign vector.
            
            target_vector = np.sign(w_subset)
            
            # corrwith(axis=0) computes correlation column-wise
            scores = Z_matrix.corrwith(target_vector, axis=0)
            
            scores.name = f"CompCRP_{name}"
            results.append(scores)
            
        if not results:
            return pd.DataFrame()
            
        return pd.concat(results, axis=1)



class CompCHIP:
    """
    CHIP-related Methylation Scores (Clonal Hematopoiesis of Indeterminate Potential).
    Computes scores for various CHIP signatures (e.g., DNMT3A, TET2).
    """

    # --- Added Metadata ---
    METADATA = {
        "year": 2024, 
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "CHIP score",
        "source": "https://doi.org/10.1038/s41467-025-59333-w"
    }

    def __init__(self, data_dir: str = None):
        self.name = "CompCHIP"
        self.metadata = self.METADATA # Store metadata
        self.signatures = {}
        
        # 1. Automatically locate resource path
        if data_dir is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data")
            
        # 2. Dynamically search for CompCHIP coefficient files
        pattern = os.path.join(data_dir, "CompCHIP_*.csv")
        files = glob.glob(pattern)
        
        if not files:
            print(f"[Warning] No CHIP model files found at {pattern}")
        
        for file_path in files:
            try:
                basename = os.path.basename(file_path)
                sig_name = basename[9:-10] # Extract name from "CompCHIP_{NAME}_coefs.csv"
                df = pd.read_csv(file_path)
                self.signatures[sig_name] = df.set_index("probe")["coef"]
            except Exception as e:
                print(f"[Error] Failed to load CHIP signature from {file_path}: {e}")
                
    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves coefficients for all loaded CHIP signatures.
        Returns a DataFrame with columns: ['model_name', 'probe', 'coef'].
        """
        dfs = []
        for model_name, series in self.signatures.items():
            # Convert Series back to DataFrame
            df = series.reset_index()
            df.columns = ['probe', 'coef']
            
            # Add identifier column
            df['model_name'] = model_name
            dfs.append(df)
            
        if not dfs:
            return pd.DataFrame(columns=['model_name', 'probe', 'coef'])
            
        return pd.concat(dfs, ignore_index=True)
    # --- Added Info Method ---
    def info(self):
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def predict(self, beta_df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        # ... (Rest of predict method remains unchanged) ...
        # [Content truncated for brevity, insert original predict logic here]
        results = []
        for name, coefs in self.signatures.items():
            common_cpgs = beta_df.index.intersection(coefs.index)
            if verbose: print(f"[CompCHIP] {name}: Found {len(common_cpgs)} / {len(coefs)} CpGs")
            if len(common_cpgs) < 2: continue
            
            X_subset = beta_df.loc[common_cpgs]
            w_subset = coefs.loc[common_cpgs]
            
            row_means = X_subset.mean(axis=1)
            row_stds = X_subset.std(axis=1)
            row_stds[row_stds == 0] = 1.0
            Z_matrix = X_subset.sub(row_means, axis=0).div(row_stds, axis=0)
            
            signs = np.sign(w_subset)
            pos_probes = signs[signs == 1].index
            neg_probes = signs[signs == -1].index
            
            if len(pos_probes) > 0: score_p = Z_matrix.loc[pos_probes].mean(axis=0)
            else: score_p = pd.Series(0, index=Z_matrix.columns)
                
            if len(neg_probes) > 0: score_n = Z_matrix.loc[neg_probes].mean(axis=0)
            else: score_n = pd.Series(0, index=Z_matrix.columns)
            
            final_score = score_p - score_n
            final_score.name = f"CompCHIP_{name}"
            results.append(final_score)
            
        if not results: return pd.DataFrame()
        return pd.concat(results, axis=1)

class EpiScores:
    """
    EpiScores: 109 validated epigenetic scores for the circulating proteome (Gadd et al., 2022).
    
    This class implements a robust two-stage imputation and scoring system:
    1. Sample-level NA Imputation: Missing values in the input data are filled using 
       the row-wise (CpG-wise) mean of the current dataset.
    2. Missing CpG Imputation: Probes required by the model but completely missing 
       from the input dataset are filled using reference mean beta values from the 
       original training cohort.
    """
    
    def __init__(self, data_dir: str = None):
        # 1. Locate resource path
        if data_dir is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data")
            
        file_path = os.path.join(data_dir, "EpiScores_All.csv")
        
        if not os.path.exists(file_path):
            print(f"[Error] EpiScores coefficient file not found: {file_path}")
            self.weights_matrix = None
            self.ref_means = None
            return

        # 2. Load and Restructure Data
        # Expected CSV Columns: probe, coef, mean_beta, trait
        try:
            df = pd.read_csv(file_path)
            
            # Construct Weights Matrix: Rows=CpGs, Columns=Traits
            # Pivoting allows us to calculate all 109 scores in a single matrix operation
            self.weights_matrix = df.pivot(index='probe', columns='trait', values='coef').fillna(0.0)
            
            # Construct Reference Means Map: Index=CpG, Value=Mean_Beta
            # Deduplicate probes since the same CpG can appear in multiple traits
            self.ref_means = df.drop_duplicates(subset='probe').set_index('probe')['mean_beta']
            
            print(f"[EpiScores] Loaded {self.weights_matrix.shape[1]} protein scores covering {self.weights_matrix.shape[0]} CpGs.")
            
        except Exception as e:
            print(f"[Error] Failed to process EpiScores file: {e}")
            self.weights_matrix = None

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves coefficients for all 109 EpiScores.
        Returns a DataFrame with columns: ['model_name', 'probe', 'coef', 'mean_beta'].
        """
        if self.raw_coefs is None:
            return pd.DataFrame()
            
        df = self.raw_coefs.copy()
        
        # Standardize column names to match the interface
        # trait -> model_name
        rename_map = {'trait': 'model_name'}
        
        # Ensure 'probe' and 'coef' are correct (usually they are already correct in CSV)
        if 'var' in df.columns: rename_map['var'] = 'probe'
        if 'beta' in df.columns: rename_map['beta'] = 'coef'
        
        df = df.rename(columns=rename_map)
        
        # Select relevant columns
        cols = ['model_name', 'probe', 'coef']
        # Optionally keep 'mean_beta' if present, as it's useful context for this model
        if 'mean_beta' in df.columns:
            cols.append('mean_beta')
            
        return df[cols]
   

    def predict(self, beta_df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        """
        Computes the 109 EpiScores.
        
        Args:
            beta_df: Methylation data (Rows=CpGs, Columns=Samples).
            verbose: If True, prints imputation progress.
            
        Returns:
            pd.DataFrame: A DataFrame (Rows=Samples, Columns=109 Traits) containing scores.
        """
        if self.weights_matrix is None:
            return pd.DataFrame()

        # --- Stage 1: Sample-level NA Imputation ---
        # Logic: "Any existing NA values... are imputed to the row-wise mean"
        if beta_df.isnull().any().any():
            if verbose: print("[EpiScores] Imputing internal NAs with row means...")
            # Efficient row-mean imputation: Transpose -> fillna -> Transpose
            row_means = beta_df.mean(axis=1)
            beta_df_imputed = beta_df.T.fillna(row_means).T
        else:
            beta_df_imputed = beta_df

        # --- Stage 2: Missing CpG Imputation ---
        # Identify required CpGs
        required_cpgs = self.weights_matrix.index
        
        if verbose: print("[EpiScores] Aligning CpGs and imputing missing probes from training reference...")
        
        # Align data using reindex.
        # This automatically introduces NaN rows for any CpGs missing in the input.
        beta_aligned = beta_df_imputed.reindex(required_cpgs)
        
        # Fill newly introduced NaNs (missing probes) with reference means (self.ref_means)
        if beta_aligned.isnull().any().any():
            # Align ref_means to the matrix shape and fill
            beta_aligned = beta_aligned.T.fillna(self.ref_means).T
            
            # Fallback: If NaNs persist (rare edge case where ref_mean is missing), fill with 0
            beta_aligned = beta_aligned.fillna(0.0)

        # --- Stage 3: Calculate Scores ---
        # Matrix Multiplication: (Samples x CpGs) * (CpGs x Traits) -> (Samples x Traits)
        
        # beta_aligned.T shape: (Samples, CpGs)
        # self.weights_matrix shape: (CpGs, Traits)
        scores = beta_aligned.T.dot(self.weights_matrix)
        
        # Note: No intercept is added as it is explicitly 0 in the model definition.
        
        return scores


class CompIL6(BaseLinearClock):
    """
    DNA Methylation-Based Proxy for IL-6 (Interleukin-6).
    
    IL-6 is a key cytokine and biomarker for systemic inflammation. This class 
    calculates a methylation score that serves as a proxy for serum IL-6 levels.
    """
    
    def __init__(self):
        clock_name = "CompIL6"
        
        try:
            # Expected filename: CompIL6.csv
            coef_df = load_clock_coefs(f"{clock_name}")
        except Exception as e:
            print(f"[Error] Failed to load {clock_name}: {e}")
            coef_df = pd.DataFrame(columns=['probe', 'coef'])
        
        # Initialize parent class (BaseLinearClock handles weighted sum calculation)
        # If the coefficients include an 'Intercept' row, BaseLinearClock automatically applies it.
        super().__init__(coef_df, name="IL6_Score")