import pandas as pd
import numpy as np
import os
import glob
from ..utils import load_clock_coefs
from .base import BaseLinearClock



class CompCRP:
    """
    DNAm-based C-Reactive Protein (CRP) Score.

    Computes a DNA methylation surrogate proxy for CRP levels using a correlation-based 
    scoring method.

    Attributes:
        name (str): Name of the model.
        metadata (dict): Metadata info.
        signatures (dict): Dictionary storing coefficients for 'CRP' and 'intCRP' models.

    References:
        Wielscher, M., et al. DNA methylation signature of chronic low-grade inflammation and its role in cardio-respiratory diseases. 
        Nat Commun (2022). https://doi.org/10.1038/s41467-022-29792-6
    """
    
    # --- Added Metadata ---
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "CRP score",
        "source": "https://doi.org/10.1038/s41467-022-29792-6"
    }

    def __init__(self):
        """Initialize CompCRP and load signature weights."""
        self.name = "CompCRP"
        self.metadata = self.METADATA # Store metadata
        self.signatures = {}
        
        # Load both models
        for model_name in ["CRP", "intCRP"]:
            try:
                # Expected filename: CompCRP_CRP.csv
                df = load_clock_coefs(f"CompCRP_{model_name}")
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
    
    # --- Added Info Method ---
    def info(self):
        """Prints summary information about the model."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def predict(self, beta_df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        """
        Calculate CRP Scores.

        The score is defined as the correlation between the sample's normalized 
        methylation profile and the signature weights.

        Args:
            beta_df (pd.DataFrame): Methylation data (Rows: CpGs, Columns: Samples).
            verbose (bool): Whether to print coverage stats.

        Returns:
            pd.DataFrame: Calculated scores (Columns: 'CompCRP_CRP', 'CompCRP_intCRP').
        """
        results = []
        for name, coefs in self.signatures.items():
            common_cpgs = beta_df.index.intersection(coefs.index)
            if verbose: print(f"[CompCRP] {name}: Found {len(common_cpgs)} / {len(coefs)} CpGs")
            if len(common_cpgs) < 2: continue
            
            X_subset = beta_df.loc[common_cpgs]
            w_subset = coefs.loc[common_cpgs]
            
            row_means = X_subset.mean(axis=1)
            row_stds = X_subset.std(axis=1)
            row_stds[row_stds == 0] = 1.0 
            Z_matrix = X_subset.sub(row_means, axis=0).div(row_stds, axis=0)
            
            target_vector = np.sign(w_subset)
            scores = Z_matrix.corrwith(target_vector, axis=0)
            scores.name = f"CompCRP_{name}"
            results.append(scores)
            
        if not results: return pd.DataFrame()
        return pd.concat(results, axis=1)



class CompCHIP:
    """
    CHIP-related Methylation Scores (Clonal Hematopoiesis of Indeterminate Potential).
    
    Computes scores for various CHIP signatures (e.g., DNMT3A, TET2) based on 
    differential methylation patterns associated with these mutations.

    Attributes:
        signatures (dict): Loaded coefficients for different CHIP mutations.

    References:
        Kirmani, S., et al. Epigenome-wide DNA methylation association study of CHIP provides insight into perturbed gene regulation 
        Nat Commun (2025). https://doi.org/10.1038/s41467-025-59333-w
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
        """
        Initialize CompCHIP by dynamically loading all available signature files.

        Args:
            data_dir (str, optional): Path to data directory. Defaults to package internal path.
        """
        self.name = "CompCHIP"
        self.metadata = self.METADATA # Store metadata
        self.signatures = {}
        
        # 1. Automatically locate resource path
        if data_dir is None:
            # Adjust path logic if necessary to match your folder structure
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data")
            
        # 2. Dynamically search for CompCHIP coefficient files
        pattern = os.path.join(data_dir, "CompCHIP_*.csv")
        files = glob.glob(pattern)
        
        if not files:
            # print(f"[Warning] No CHIP model files found at {pattern}")
            pass
        
        for file_path in files:
            try:
                basename = os.path.basename(file_path)
                
                # --- FIX: ROBUST PARSING LOGIC ---
                # 1. Remove extension
                filename_no_ext = os.path.splitext(basename)[0]
                
                # 2. Remove prefix "CompCHIP_"
                if filename_no_ext.startswith("CompCHIP_"):
                    sig_name = filename_no_ext[9:]
                else:
                    sig_name = filename_no_ext
                    
                # 3. Remove suffix "_coefs" if present
                if sig_name.endswith("_coefs"):
                    sig_name = sig_name[:-6]
                
                # 4. Strip extra characters
                sig_name = sig_name.strip('_')
                
                if not sig_name:
                    continue
                # ----------------------------------

                df = pd.read_csv(file_path)
                
                # Ensure correct column names
                if 'coef' not in df.columns and 'beta' in df.columns:
                    df = df.rename(columns={'beta': 'coef'})
                if 'probe' not in df.columns and 'var' in df.columns:
                    df = df.rename(columns={'var': 'probe'})

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

    def info(self):
        """Prints summary information about the model."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def predict(self, beta_df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        """
        Calculate CHIP Scores.

        Args:
            beta_df (pd.DataFrame): Methylation data.
            verbose (bool): Print warnings if CpGs are missing.

        Returns:
            pd.DataFrame: DataFrame containing scores for each CHIP signature.
        """
        results = []
        
        for name, coefs in self.signatures.items():
            common_cpgs = beta_df.index.intersection(coefs.index)
            
            if verbose:
                # print(f"[CompCHIP] {name}: Found {len(common_cpgs)} / {len(coefs)} CpGs")
                pass
                
            if len(common_cpgs) < 2:
                if verbose: print(f"[Warning] Too few CpGs for {name}. Skipping.")
                continue
                
            X_subset = beta_df.loc[common_cpgs]
            w_subset = coefs.loc[common_cpgs]
            
            row_means = X_subset.mean(axis=1)
            row_stds = X_subset.std(axis=1)
            row_stds[row_stds == 0] = 1.0
            
            Z_matrix = X_subset.sub(row_means, axis=0).div(row_stds, axis=0)
            
            signs = np.sign(w_subset)
            pos_probes = signs[signs == 1].index
            neg_probes = signs[signs == -1].index
            
            if len(pos_probes) > 0:
                score_p = Z_matrix.loc[pos_probes].mean(axis=0)
            else:
                score_p = pd.Series(0, index=Z_matrix.columns)
                
            if len(neg_probes) > 0:
                score_n = Z_matrix.loc[neg_probes].mean(axis=0)
            else:
                score_n = pd.Series(0, index=Z_matrix.columns)
            
            final_score = score_p - score_n
            final_score.name = f"CompCHIP_{name}"
            results.append(final_score)
            
        if not results:
            return pd.DataFrame()
            
        return pd.concat(results, axis=1)



class EpiScores:
    """
    EpiScores: 109 validated epigenetic scores for the circulating proteome.
    
    Predicts the levels of 109 plasma proteins based on DNA methylation data.
    Implements automatic imputation using training set reference means.

    References:
        Gadd DA, et al. Epigenetic scores for the circulating proteome as tools for disease prediction.
        eLife (2022). https://doi.org/10.7554/eLife.71802
    """
    
    # --- Added Metadata ---
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Levels of 109 Plasma Proteins",
        "source": "https://doi.org/10.7554/eLife.71802",
        "doi": "10.7554/eLife.71802"
    }

    def __init__(self, data_dir: str = None):
        """Initialize EpiScores and load the massive coefficient matrix."""
        self.name = "EpiScores"
        self.metadata = self.METADATA # Store metadata
        
        # 1. Locate resource path
        if data_dir is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data")
            
        file_path = os.path.join(data_dir, "EpiScores_All.csv")
        self.raw_coefs = None
        
        if not os.path.exists(file_path):
            print(f"[Error] EpiScores coefficient file not found: {file_path}")
            self.weights_matrix = None
            self.ref_means = None
            return

        # 2. Load and Restructure Data
        try:
            df = pd.read_csv(file_path)
            self.raw_coefs = df
            self.weights_matrix = df.pivot(index='probe', columns='trait', values='coef').fillna(0.0)
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
    # --- Added Info Method ---
    def info(self):
        """Prints summary information."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def predict(self, beta_df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        """
        Predict Protein Levels.

        Args:
            beta_df (pd.DataFrame): Methylation data (Rows: CpGs, Columns: Samples).
            verbose (bool): Whether to print status messages.

        Returns:
            pd.DataFrame: Predicted levels for 109 proteins (Rows: Samples, Columns: Proteins).
        """
        if self.weights_matrix is None: return pd.DataFrame()
        
        if beta_df.isnull().any().any():
            if verbose: print("[EpiScores] Imputing internal NAs with row means...")
            row_means = beta_df.mean(axis=1)
            beta_df_imputed = beta_df.T.fillna(row_means).T
        else:
            beta_df_imputed = beta_df

        required_cpgs = self.weights_matrix.index
        if verbose: print("[EpiScores] Aligning CpGs and imputing missing probes from training reference...")
        beta_aligned = beta_df_imputed.reindex(required_cpgs)
        
        if beta_aligned.isnull().any().any():
            beta_aligned = beta_aligned.T.fillna(self.ref_means).T
            beta_aligned = beta_aligned.fillna(0.0)

        scores = beta_aligned.T.dot(self.weights_matrix)
        return scores



class CompIL6(BaseLinearClock):
    """
    DNA Methylation-Based Proxy for IL-6 (Interleukin-6).
    
    Estimates chronic inflammation levels via an IL-6 specific methylation signature.
    """

    # --- Added Metadata ---
    METADATA = {
        "year": 2024, # Adjust based on specific paper year
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "IL-6 score",
        "source": "https://doi.org/10.1093/gerona/glab046)"
    }
    
    def __init__(self):
        """Initialize CompIL6 score."""
        clock_name = "CompIL6"
        try:
            coef_df = load_clock_coefs(f"{clock_name}")
        except Exception as e:
            print(f"[Error] Failed to load {clock_name}: {e}")
            coef_df = pd.DataFrame(columns=['probe', 'coef'])
        
        # Pass metadata to parent class
        super().__init__(coef_df, name="IL6_Score", metadata=self.METADATA)


        