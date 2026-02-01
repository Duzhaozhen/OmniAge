import pandas as pd
import numpy as np
from typing import Dict, Union

class BaseLinearClock:
    """
    Base class for linear epigenetic aging clocks.
    
    This class implements the standard weighted sum model for epigenetic clocks.
    It handles feature alignment (intersection of CpGs) and handles intercept logic.

    Attributes:
        name (str): The name of the clock.
        metadata (dict): Metadata associated with the clock (e.g., species, tissue).
        intercept (float): The model intercept (offset).
        weights (pd.Series): The model coefficients (weights), indexed by probe names.
    """
    def __init__(self, coef_df: pd.DataFrame, name: str,metadata:Dict=None):
        """
        Initialize the linear clock model.

        Args:
            coef_df (pd.DataFrame): A DataFrame containing the model coefficients.
                Must contain columns:
                - 'probe': The feature names (e.g., CpG IDs).
                - 'coef': The weight values.
                Rows with 'probe' in ['Intercept', '(Intercept)'] will be treated as the offset.
            name (str): Identifier for the clock.
            metadata (Optional[Dict], optional): Extra info (year, paper, etc.). Defaults to None.

        Raises:
            ValueError: If 'probe' or 'coef' columns are missing in coef_df.
        """
        self.name = name
        self.metadata = metadata or {}
        # --- 0. Safety Check ---
        required_columns = {'probe', 'coef'}
        if not required_columns.issubset(coef_df.columns):
            raise ValueError(f"[{name}] Input dataframe must contain columns: {required_columns}. "
                             f"Found: {list(coef_df.columns)}")

        # --- 1. Parsing Coefficients ---
        # Intercept
        intercept_mask = coef_df['probe'].isin(['Intercept', '(Intercept)','intercept'])
        intercept_row = coef_df[intercept_mask]
        
        if not intercept_row.empty:
            self.intercept = float(intercept_row.iloc[0]['coef'])
            weights_df = coef_df[~intercept_mask]
        else:
            self.intercept = 0.0
            weights_df = coef_df

        # --- 2. Store as Series ---
        self.weights = pd.Series(
            weights_df['coef'].values, 
            index=weights_df['probe'].values
        )
        
    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the model coefficients (weights + intercept) as a pandas DataFrame.
        
        Returns:
            pd.DataFrame: A DataFrame with columns ['probe', 'coef'].
                          The intercept is included as the first row with probe name '(Intercept)'.
        """
        # 1. Convert weights Series back to DataFrame
        # reset_index() turns the index (probe names) into a column
        df = self.weights.reset_index()
        df.columns = ['probe', 'coef']
        
        # 2. Create the intercept row
        # We standardize the name to '(Intercept)' regardless of input format
        intercept_df = pd.DataFrame([['(Intercept)', self.intercept]], columns=['probe', 'coef'])
        
        # 3. Combine intercept + weights
        final_df = pd.concat([intercept_df, df], ignore_index=True)
        
        return final_df
        
    def info(self):
        """Prints summary information about the clock model."""
        print(f"[{self.name}] Model information:")
        if not self.metadata:
            print("No metadata available.")
            return
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")


    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        Preprocessing hook.
        By default, this returns the data unchanged.
        """
        return beta_df

    def postprocess(self, linear_predictor):
        """Hook for post-processing (e.g., anti-log transform). Default is identity."""
        return linear_predictor

    def predict(self, beta_df: pd.DataFrame, verbose=True) -> pd.Series:
        """
        Executes the age prediction pipeline: Preprocess -> Linear Sum -> Postprocess.

        Args:
            beta_df (pd.DataFrame): Methylation data matrix.
                - Rows: Features (CpGs).
                - Columns: Samples.
            verbose (bool, optional): If True, prints feature coverage stats. Defaults to True.

        Returns:
            pd.Series: Predicted biological age values, indexed by sample IDs from beta_df.
        """
        # --- Step 0: Apply Pre-processing ---
        beta_df_processed = self.preprocess(beta_df)

        # --- Step 1: Feature Alignment ---
        common_features = self.weights.index.intersection(beta_df_processed.index)
        
        if verbose:
            n_model = len(self.weights)
            n_found = len(common_features)
            print(f"[{self.name}] Number of represented CpGs (max={n_model})={n_found}")

        if len(common_features) == 0:
            return pd.Series([self.intercept] * beta_df_processed.shape[1], index=beta_df_processed.columns)

        # --- Step 2: Linear Calculation ---
        X_subset = beta_df_processed.loc[common_features, :].T 
        w_subset = self.weights.loc[common_features]
        
        linear_predictor = X_subset.dot(w_subset) + self.intercept
        
        # --- Step 3: Post-processing ---
        final_age = self.postprocess(linear_predictor)
        
        return final_age
