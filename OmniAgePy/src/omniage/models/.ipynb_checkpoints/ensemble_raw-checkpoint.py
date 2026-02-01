import os
import pandas as pd
import glob
from typing import Optional, Union, List
from .base import BaseLinearClock

class EnsembleAge:
    """
    EnsembleAge: Enhancing epigenetic age assessment with a multi-clock framework.
    
    This class manages the loading and prediction of multiple sub-clocks (e.g., 40 individual 
    linear models for the 'Dynamic' version) as described by Haghani et al. (2025).
    
    Supported Versions:
    - 'HumanMouse': Cross-species predictors.
    - 'Static': Mouse-specific static estimators.
    - 'Dynamic': Mouse-specific dynamic estimators (typically 40 sub-clocks).
    """
    
    def __init__(self, version: str = "HumanMouse", data_dir: str = None):
        """
        Initialize the EnsembleAge framework.
        
        Args:
            version (str): The specific ensemble version ('HumanMouse', 'Static', or 'Dynamic').
            data_dir (str, optional): Custom path to the data directory. 
                                      If None, defaults to the package's internal 'data' folder.
        """
        self.version = version
        self.sub_clocks = [] # List to store initialized BaseLinearClock instances
        
        # 1. Automatically locate resource path
        if data_dir is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data")
        
        # 2. Search for matching model files
        # Expected filename format: EnsembleAge_{version}_{sub_name}_coefs.csv
        pattern = f"EnsembleAge_{version}_*_coefs.csv"
        search_path = os.path.join(data_dir, "EnsembleAge", pattern)
        
        files = glob.glob(search_path)
        
        if not files:
            print(f"[Warning] No model files found for EnsembleAge version '{version}' at {search_path}")
            return
            
        print(f"[EnsembleAge] Found {len(files)} sub-clocks for version '{version}'. Loading...")
        
        # 3. Load each sub-clock
        for file_path in files:
            # Parse sub-clock name from filename
            # Example filename: EnsembleAge_Dynamic_Clock1_coefs.csv
            basename = os.path.basename(file_path)
            
            # Remove prefix ("EnsembleAge_{version}_") and suffix ("_coefs.csv")
            # Prefix length: len("EnsembleAge_") + len(version) + 1 (for the underscore)
            prefix_len = len("EnsembleAge_") + len(version) + 1
            suffix_len = len("_coefs.csv")
            
            sub_name = basename[prefix_len:-suffix_len]
            full_name = f"{version}_{sub_name}"
            
            try:
                # Read coefficients
                df = pd.read_csv(file_path)
                
                # Create a linear clock instance
                # BaseLinearClock automatically handles intercept extraction (row where probe='Intercept') 
                # and coefficient alignment.
                clock = BaseLinearClock(coef_df=df, name=full_name)
                self.sub_clocks.append(clock)
                
            except Exception as e:
                print(f"[Error] Failed to load sub-clock {basename}: {e}")

    def predict(self, beta_df: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
        """
        Execute prediction for all loaded sub-clocks.
        
        Args:
            beta_df: A DataFrame of methylation values (Rows=CpGs, Columns=Samples).
            verbose (bool): If True, prints progress for each sub-clock.
            
        Returns:
            pd.DataFrame: A DataFrame (Rows=Samples, Columns=Clock_Names) containing 
                          predictions for all sub-clocks.
        """
        if not self.sub_clocks:
            return pd.DataFrame()
            
        results = []
        
        # Iterate through all sub-models to perform prediction
        for clock in self.sub_clocks:
            # clock.predict returns a Series (index=Samples)
            pred = clock.predict(beta_df, verbose=verbose)
            pred.name = clock.name # Set the Series name to the clock name
            results.append(pred)
            
        if not results:
            return pd.DataFrame()
            
        # Concatenate results: Sample x Clock_Name
        final_df = pd.concat(results, axis=1)
        
        return final_df

    def calculate_dynamic_score(
        self, 
        beta_df: pd.DataFrame, 
        group_col: Optional[pd.Series] = None
    ) -> pd.DataFrame:
        """
        Helper method for 'Dynamic' version logic.
        
        In the R implementation, this function might handle complex logic like:
        1. Calculating all constituent scores.
        2. (If group_col is provided) Identifying responsive clocks.
        3. Returning aggregated metrics (e.g., median).
        
        Currently, this method delegates to `predict` and returns the raw dataframe 
        of all sub-clock scores, consistent with the base R function behavior.
        
        Args:
            beta_df: Methylation data.
            group_col: (Optional) Series indicating sample groups for responsive analysis.
            
        Returns:
            pd.DataFrame: Predictions for all sub-clocks.
        """
        return self.predict(beta_df)