import pandas as pd
from typing import Optional, Union, List, Dict
import os
from .base import BaseLinearClock
from ..utils import load_clock_coefs



class PetkovichMouse(BaseLinearClock):
    """
    Petkovich Blood Mouse DNA Methylation Age Predictor.
    Attributes:
        intercept (float): Model offset.
        weights (pd.Series): Coefficients for the 90 CpGs.

    References:
        Petkovich et al. Using DNA Methylation Profiling to Evaluate Biological Age and Longevity Interventions.
        Cell Metab. (2017). https://doi.org/10.1016/j.cmet.2017.03.016
    """
    METADATA = {
        "year": 2017,
        "species": "Mouse",
        "tissue": "Blood",
        "omic type": "RRBS",
        "prediction": "chronological age(months)",
        "source": "https://doi.org/10.1016/j.cmet.2017.03.016"
    }

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the PetkovichMouse model.
        
        Args:
            coef_df: Optional DataFrame containing model coefficients (CpG IDs and weights).
                If None, the default 'PetkovichMouse.csv' will be loaded from the package data.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            # Attempts to locate 'PetkovichMouse.csv' within the internal data directory
            coef_df = load_clock_coefs("PetkovichMouse")
            
        # 2. Initialize the base class
        super().__init__(coef_df, name="PetkovichMouse",metadata=self.METADATA)

    def postprocess(self, linear_predictor):
        """
        Applies a convertion from the output of an ElasticNet to mouse age in months.
        """
        a = 0.1666
        b = 0.4185
        c = -1.712
        age = ((linear_predictor - c) / a) ** (1 / b)
        age = age / 30.5  # days to months
        return age

    
class ThompsonMouse(BaseLinearClock):
    """
    Thompson Mouse Multi-Tissue DNA Methylation Age Predictor.
    Attributes:
        intercept (float): Model offset.
        weights (pd.Series): Coefficients for the 595 CpGs.

    References:
        PetkovichMouse et al. A multi-tissue full lifespan epigenetic clock for mice.
        Aging. (2018). https://doi.org/10.18632/aging.101590
    """
    METADATA = {
        "year": 2018,
        "species": "Mouse",
        "tissue": "Multi-tissue",
        "omic type": "RRBS",
        "prediction": "chronological age(months)",
        "source": "https://doi.org/10.18632/aging.101590"
    }

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the ThompsonMouse model.
        
        Args:
            coef_df: Optional DataFrame containing model coefficients (CpG IDs and weights).
                If None, the default 'ThompsonMouse.csv' will be loaded from the package data.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            # Attempts to locate 'ThompsonMouse.csv' within the internal data directory
            coef_df = load_clock_coefs("ThompsonMouse")
            
        # 2. Initialize the base class
        super().__init__(coef_df, name="ThompsonMouse",metadata=self.METADATA)




class MeerMouse(BaseLinearClock):
    """
    Meer Whole Lifespan Multi-Tissue Mouse DNA Methylation Age Predictor.
    Attributes:
        intercept (float): Model offset.
        weights (pd.Series): Coefficients for the 435 CpGs.

    References:
        Meer et al. A whole lifespan mouse multi-tissue DNA methylation clock.
        Elife (2018). https://doi.org/10.7554/eLife.40675
    """
    METADATA = {
        "year": 2018,
        "species": "Mouse",
        "tissue": "Multi-tissue",
        "omic type": "RRBS",
        "prediction": "chronological age(months)",
        "source": "https://doi.org/10.7554/eLife.40675"
    }

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the MeerMouse model.
        
        Args:
            coef_df: Optional DataFrame containing model coefficients (CpG IDs and weights).
                If None, the default 'MeerMouse.csv' will be loaded from the package data.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            # Attempts to locate 'MeerMouse.csv' within the internal data directory
            coef_df = load_clock_coefs("MeerMouse")
            
        # 2. Initialize the base class
        super().__init__(coef_df, name="MeerMouse",metadata=self.METADATA)

    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        Validates that input data are Beta values (0-1) and scales them to percentages (0-100).

        Args:
            beta_df (pd.DataFrame): Input methylation matrix (Features x Samples).

        Returns:
            pd.DataFrame: Scaled data (0-100) for internal calculation.

        Raises:
            ValueError: If input values are outside the expected [0, 1] range.
        """
        # 1. Calculate global range (ignoring NaNs)
        max_val = beta_df.max().max()
        min_val = beta_df.min().min()

        # 2. Strict validation: Ensure data is within [0, 1]
        # A small margin (0.01) is allowed for floating-point noise
        if max_val > 1.01 or min_val < -0.01:
            raise ValueError(
                f"[{self.name}] Input data values must be between 0 and 1 (Beta values). "
                f"Found range: [{min_val:.3f}, {max_val:.3f}]. "
                "Please ensure your input is not already in percentage (0-100) or M-values."
            )

        # 3. Scale to 0-100 as required by the Meer model coefficients
        return beta_df * 100

    def postprocess(self, linear_predictor):
        """
        Applies a convertion from the output of an ElasticNet to mouse age in months.
        """
        return linear_predictor/30.5





