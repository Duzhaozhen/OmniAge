import os
import glob
import pandas as pd
import numpy as np

class CompSmokeIndex:
    """
    DNAm-based Smoking Index (S-index).
    
    Computes a smoking index using 1,501 smoking-associated CpG sites.
    This index correlates with cancer risk across multiple tissues.

    Attributes:
        signatures (dict): Loaded coefficients for the smoking index.

    References:
        Teschendorff, Andrew E et al. Correlation of Smoking-Associated DNA 
        Methylation Changes in Buccal Cells With DNA Methylation Changes in 
        Epithelial Cancer. JAMA oncology (2015).
    """

    METADATA = {
        "year": 2015,
        "species": "Human",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(450k)",
        "prediction": "Smoking Index",
        "source": "https://doi.org/10.1001/jamaoncol.2015.1053"
    }

    def __init__(self, data_dir: str = None):
        """
        Initialize the model and load coefficients.

        Args:
            data_dir (str, optional): Path to the directory containing 'CompSmokeIndex.csv'.
        """
        self.name = "CompSmokeIndex"
        self.metadata = self.METADATA
        self.coefficients = None

        # Resolve data path
        if data_dir is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data")
        
        file_path = os.path.join(data_dir, "CompSmokeIndex.csv")

        # Load coefficients
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path)
                # Standardize column names
                df = df.rename(columns={'beta': 'coef', 'var': 'probe'})
                self.coefficients = df.set_index("probe")["coef"]
            except Exception as e:
                print(f"[Error] Failed to load coefficients: {e}")
        else:
            print(f"[Warning] Resource not found: {file_path}")

    def info(self):
        """Display model metadata."""
        print(f"[{self.name}] Configuration:")
        for key, value in self.metadata.items():
            print(f"  - {key.replace('_', ' ').capitalize()}: {value}")

    def predict(self, beta_df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        """
        Predict Smoking Index scores for given samples.

        Args:
            beta_df (pd.DataFrame): DNAm beta values (Rows: Probes, Cols: Samples).
            verbose (bool): If True, prints processing status.

        Returns:
            pd.DataFrame: S-index scores indexed by sample ID.
        """
        if self.coefficients is None:
            raise ValueError("Model coefficients not loaded.")

        # Align probes
        common_probes = beta_df.index.intersection(self.coefficients.index)
        if verbose:
            print(f"[{self.name}] Overlapping probes: {len(common_probes)}/{len(self.coefficients)}")
            
        if len(common_probes) < 2:
            return pd.DataFrame()
            
        # Subset and align
        X = beta_df.loc[common_probes]
        w = self.coefficients.loc[common_probes]
        
        # Row-wise Z-score transformation (standardize probes across samples)
        row_means = X.mean(axis=1)
        row_stds = X.std(axis=1).replace(0, 1.0)
        Z = X.sub(row_means, axis=0).div(row_stds, axis=0)
        
        # Compute weighted average (S-index)
        # S = mean(Z * coefficients)
        scores = Z.multiply(w, axis=0).mean(axis=0)
        
        return pd.DataFrame(scores, columns=[self.name])