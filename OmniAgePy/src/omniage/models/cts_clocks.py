import pandas as pd
import numpy as np
from typing import Optional, Literal
from sklearn.linear_model import LinearRegression
from .base import BaseLinearClock
from ..utils import load_clock_coefs

class CTSClock(BaseLinearClock):
    """
    Base class for Cell-Type Specific (CTS) Clocks (Tong et al. 2024).
    
    This class supports two distinct modeling strategies:
    
    1. Intrinsic Clocks (Neu-In, Glia-In, Brain): 
       - Quantify epigenetic aging independent of cell composition changes.
       - Requires 'sorted' or 'bulk' data type specification.
       - For 'bulk' data, Cell Type Fractions (CTF) are regressed out to isolate residuals.
       - Applies row-wise (per-CpG) Z-score standardization.
       
    2. Semi-intrinsic Clocks (Neu-Sin, Glia-Sin, Hep, Liver):
       - Quantify aging with partial adjustment for cell composition.
       - operate directly on raw beta values without standardization.
       
    Reference:
    Tong H, et al. Cell-type specific epigenetic clocks to quantify biological 
    age at cell-type resolution. Aging (2024).
    """
    
    def __init__(self, clock_name, mode: Literal['intrinsic', 'semi-intrinsic'], metadata: Optional[dict] = None):
        """
        Initialize the CTSClock.

        Args:
            clock_name (str): The name of the clock (e.g., 'Neu-In', 'Liver').
            mode (Literal['intrinsic', 'semi-intrinsic']): The modeling strategy.
            metadata (Optional[dict]): Additional metadata for the clock. Defaults to None.
        """
        self.mode = mode
        # Automatically load coefficients from the data/CTS directory
        coef_df = load_clock_coefs(f"CTS/{clock_name}")
        super().__init__(coef_df, name=clock_name, metadata=metadata)

    def _process_intrinsic_data(
        self, 
        beta_df: pd.DataFrame, 
        data_type: str, 
        ctf_df: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Implements the data preprocessing logic for Intrinsic clocks.
        
        Preprocessing Strategy:
        1. 'sorted': Apply row-wise Z-score normalization (across samples).
        2. 'bulk': Regress out CTF effects -> Extract Residuals -> Apply row-wise Z-score.
        """
        # --- Scenario 1: Sorted Data (Purified Cell Types) ---
        if data_type == 'sorted':
            # R equivalent: (data - rowMeans) / rowSD
            # Calculate mean and std for each CpG across all samples (axis=1)
            means = beta_df.mean(axis=1)
            stds = beta_df.std(axis=1)
            
            # Prevent division by zero for constant rows
            stds = stds.replace(0, 1)
            
            # Perform Row-wise Standardization
            # sub/div with axis=0 aligns the vector to the DataFrame index (rows)
            return beta_df.sub(means, axis=0).div(stds, axis=0)

        # --- Scenario 2: Bulk Data (Mixed Tissue) ---
        elif data_type == 'bulk':
            if ctf_df is None:
                raise ValueError(f"[{self.name}] 'ctf' (Cell Type Fractions) is required for bulk tissue predictions with Intrinsic clocks.")
            
            # Align samples (Intersection of columns and rows)
            common_samples = beta_df.columns.intersection(ctf_df.index)
            if len(common_samples) == 0:
                raise ValueError(f"[{self.name}] No common samples between methylation data columns and CTF index.")
            
            # Prepare regression matrices
            Y = beta_df[common_samples].T # Target: (Sample x CpG)
            X = ctf_df.loc[common_samples] # Predictor: (Sample x CellType)
            
            # Batch Linear Regression: Y ~ X (Regress out cell type composition)
            print(f"[{self.name}] Regressing out cell type effects for {len(common_samples)} samples...")
            
            reg = LinearRegression(fit_intercept=True, n_jobs=-1)
            reg.fit(X, Y)
            predicted = reg.predict(X)
            residuals = Y - predicted
            
            # Transpose back to (CpG x Sample)
            residuals = residuals.T
            
            # Apply Row-wise Z-score Standardization on Residuals
            means = residuals.mean(axis=1)
            stds = residuals.std(axis=1)
            stds = stds.replace(0, 1)
            
            return residuals.sub(means, axis=0).div(stds, axis=0)
            
        else:
            raise ValueError(f"Unknown data_type: {data_type}. Must be 'bulk' or 'sorted'.")

    def predict(
        self, 
        beta_df: pd.DataFrame, 
        data_type: str = 'bulk', 
        ctf: Optional[pd.DataFrame] = None, 
        verbose: bool = True
    ) -> pd.Series:
        """
        Execute age prediction using the specified CTS clock model.

        Args:
            beta_df (pd.DataFrame): Input methylation data (Rows: CpGs, Columns: Samples).
            data_type (str, optional): Type of input data ('bulk' or 'sorted'). Defaults to 'bulk'.
            ctf (Optional[pd.DataFrame], optional): Cell Type Fractions, required for 'bulk' intrinsic clocks. 
                (Rows: Samples, Columns: Cell Types). Defaults to None.
            verbose (bool, optional): Whether to print progress messages. Defaults to True.

        Returns:
            pd.Series: Predicted biological age values, indexed by sample IDs.
            
        Examples:
            # Scenario 1: Intrinsic Clock with Bulk Tissue (Requires CTF)
            >>> clock = NeuIn()
            >>> age = clock.predict(beta_bulk, data_type='bulk', ctf=neuron_proportions)

            # Scenario 2: Intrinsic Clock with Sorted/Purified Cells
            >>> clock = NeuIn()
            >>> age = clock.predict(beta_sorted, data_type='sorted')

            # Scenario 3: Semi-Intrinsic Clocks (e.g., Neu-Sin, Liver)
            >>> clock = NeuSin()
            >>> age = clock.predict(beta_values)
        """
        # --- 1. Pre-processing ---
        if self.mode == 'intrinsic':
            if verbose: print(f"[{self.name}] Performing Intrinsic preprocessing ({data_type})...")
            processed_df = self._process_intrinsic_data(beta_df, data_type, ctf)
        else:
            # Semi-intrinsic clocks use raw beta values without adjustment
            processed_df = beta_df
            
        # --- 2. Linear Prediction ---
        # Delegate to parent class for feature alignment and weighted summation
        return super().predict(processed_df, verbose=verbose)

# --- Specific Implementations ---

# 1. Intrinsic Clocks 
# --- Specific Implementations ---

# 1. Intrinsic Clocks 
class NeuIn(CTSClock):
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Brain(Neurons)",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.206184"
    }
    def __init__(self): 
        """
        Initialize the Neu-In (Neuron Intrinsic) clock.
        
        This model is pre-configured with:
        - Mode: Intrinsic (Requires cell type fractions 'ctf' if applied to bulk data)
        - Tissue: Brain (Neurons)
        """
        super().__init__("Neu-In", mode='intrinsic', metadata=self.METADATA)

class GliaIn(CTSClock):
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Brain(Glia)",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.206184"
    }
    def __init__(self): 
        """
        Initialize the Glia-In (Glia Intrinsic) clock.
        
        This model is pre-configured with:
        - Mode: Intrinsic (Requires cell type fractions 'ctf' if applied to bulk data)
        - Tissue: Brain (Glia)
        """
        super().__init__("Glia-In", mode='intrinsic', metadata=self.METADATA)

class Brain(CTSClock):
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Brain(Bulk)",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.206184"
    }
    def __init__(self): 
        """
        Initialize the Brain (Bulk Intrinsic) clock.
        
        This model is pre-configured with:
        - Mode: Intrinsic (Requires cell type fractions 'ctf' if applied to bulk data)
        - Tissue: Brain (Bulk)
        """
        super().__init__("Brain", mode='intrinsic', metadata=self.METADATA)

# 2. Semi-Intrinsic Clocks 
class NeuSin(CTSClock):
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Brain(Neurons)",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.206184"
    }
    def __init__(self): 
        """
        Initialize the Neu-Sin (Neuron Semi-intrinsic) clock.
        
        This model is pre-configured with:
        - Mode: Semi-intrinsic (Does NOT use cell type fractions)
        - Tissue: Brain (Neurons)
        """
        super().__init__("Neu-Sin", mode='semi-intrinsic', metadata=self.METADATA)

class GliaSin(CTSClock):
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Brain(Glia)",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.206184"
    }
    def __init__(self): 
        """
        Initialize the Glia-Sin (Glia Semi-intrinsic) clock.
        
        This model is pre-configured with:
        - Mode: Semi-intrinsic (Does NOT use cell type fractions)
        - Tissue: Brain (Glia)
        """
        super().__init__("Glia-Sin", mode='semi-intrinsic', metadata=self.METADATA)

class Hep(CTSClock):
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Liver(Hepatocytes)",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.206184"
    }
    def __init__(self): 
        """
        Initialize the Hep (Hepatocyte Semi-intrinsic) clock.
        
        This model is pre-configured with:
        - Mode: Semi-intrinsic (Does NOT use cell type fractions)
        - Tissue: Liver (Hepatocytes)
        """
        super().__init__("Hep", mode='semi-intrinsic', metadata=self.METADATA)

class Liver(CTSClock):
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Liver(Bulk)",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.206184"
    }
    def __init__(self): 
        """
        Initialize the Liver (Bulk Semi-intrinsic) clock.
        
        This model is pre-configured with:
        - Mode: Semi-intrinsic (Does NOT use cell type fractions)
        - Tissue: Liver (Bulk)
        """
        super().__init__("Liver", mode='semi-intrinsic', metadata=self.METADATA)