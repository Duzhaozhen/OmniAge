import pandas as pd
from .base import BaseLinearClock
from ..utils import load_clock_coefs

class BohlinGA(BaseLinearClock):
    """
    Bohlin Gestational Age Clock (Cord Blood).

    This clock predicts gestational age based on CpGs described in Bohlin et al. (2016).
    
    **Unit Conversion Note:**
    The original model calculates age in **days**. This implementation automatically 
    divides the result by 7 to return **Gestational Age in Weeks**.

    Attributes:
        intercept (float): The model intercept (automatically loaded).
        weights (pd.Series): The model coefficients (automatically loaded).

    References:
        Bohlin J, et al. Prediction of gestational age based on genome-wide 
        differentially methylated regions. Genome Biol (2016).
        https://doi.org/10.1186/s13059-016-1063-4
    """
    METADATA = {
        "year": 2016,
        "species": "Human",
        "tissue": "Cord blood",
        "omic type": "DNAm(450k)",
        "prediction": "Gestational Age(weeks)",
        "source": "https://doi.org/10.1186/s13059-016-1063-4"
    }
    
    def __init__(self):
        """
        Initialize the Bohlin GA clock.
        
        This model is pre-configured and requires no arguments.
        Coefficients are loaded automatically from 'Bohlin_GA'.
        """
        clock_name = "Bohlin_GA"
        try:
            coef_df = load_clock_coefs("Bohlin_GA")
        except Exception as e:
            print(f"[Error] Failed to load {clock_name}: {e}")
            coef_df = pd.DataFrame(columns=['probe', 'coef'])
        
        super().__init__(coef_df, name=clock_name, metadata=self.METADATA)

    def postprocess(self, linear_predictor: pd.Series) -> pd.Series:
        """
        Converts the linear prediction from days to weeks.
        
        Formula: Age_weeks = Age_days / 7.0
        """
        return linear_predictor / 7.0

class EPICGA(BaseLinearClock):
    """
    EPIC Gestational Age Clock (Haftorn et al. 2021).

    Predicts gestational age using CpGs compatible with the Illumina EPIC array.
    
    **Unit Conversion Note:**
    The original model calculates age in **days**. This implementation automatically 
    divides the result by 7 to return **Gestational Age in Weeks**.

    References:
        Haftorn KL, et al. An EPIC predictor of gestational age and its application 
        to newborns conceived by assisted reproductive technologies. 
        Clin Epigenetics (2021). https://doi.org/10.1186/s13148-021-01055-z
    """
    METADATA = {
        "year": 2021,
        "species": "Human",
        "tissue": "Cord blood",
        "omic type": "DNAm(EPIC)",
        "prediction": "Gestational Age(weeks)",
        "source": "https://doi.org/10.1186/s13148-021-01055-z"
    }
    
    def __init__(self):
        """Initialize the EPIC GA clock (Pre-configured)."""
        clock_name = "EPIC_GA"
        try:
            coef_df = load_clock_coefs("EPIC_GA")
        except Exception as e:
            print(f"[Error] Failed to load {clock_name}: {e}")
            coef_df = pd.DataFrame(columns=['probe', 'coef'])
        
        super().__init__(coef_df, name=clock_name, metadata=self.METADATA)

    def postprocess(self, linear_predictor: pd.Series) -> pd.Series:
        """
        Converts the linear prediction from days to weeks.
        
        Formula: Age_weeks = Age_days / 7.0
        """
        return linear_predictor / 7.0

class KnightGA(BaseLinearClock):
    """
    Knight Gestational Age Clock (Knight et al. 2016).

    Predicts gestational age. 
    
    References:
        Knight AK, et al. An epigenetic clock for gestational age at birth 
        based on blood methylation data. Genome Biol (2016).
        https://doi.org/10.1186/s13059-016-1068-z
    """
    METADATA = {
        "year": 2016,
        "species": "Human",
        "tissue": "Cord blood or blood spots",
        "omic type": "DNAm(27K/450K)",
        "prediction": "Gestational Age(weeks)",
        "source": "https://doi.org/10.1186/s13059-016-1068-z"
    }
    
    def __init__(self):
        """Initialize the Knight GA clock (Pre-configured)."""
        clock_name = "Knight_GA"
        try:
            coef_df = load_clock_coefs("Knight_GA")
        except Exception as e:
            print(f"[Error] Failed to load {clock_name}: {e}")
            coef_df = pd.DataFrame(columns=['probe', 'coef'])
        
        super().__init__(coef_df, name=clock_name, metadata=self.METADATA)

class MayneGA(BaseLinearClock):
    """
    Mayne Placental Gestational Age Clock (Mayne et al. 2017).

    A clock specifically designed for **Placental** tissue. 
    Do not apply this to blood samples.

    References:
        Mayne BT, et al. Accelerated placental aging in early onset preeclampsia. 
        Epigenomics (2017). https://doi.org/10.2217/epi-2016-0103
    """
    METADATA = {
        "year": 2017,
        "species": "Human",
        "tissue": "Placenta",
        "omic type": "DNAm(27K/450K)",
        "prediction": "Gestational Age(weeks)",
        "source": "https://doi.org/10.2217/epi-2016-0103"
    }
    
    def __init__(self):
        """Initialize the Mayne Placental clock (Pre-configured)."""
        clock_name = "Mayne_GA"
        try:
            coef_df = load_clock_coefs(f"Mayne_GA")
        except Exception as e:
            print(f"[Error] Failed to load {clock_name}: {e}")
            coef_df = pd.DataFrame(columns=['probe', 'coef'])
        
        super().__init__(coef_df, name=clock_name, metadata=self.METADATA)


# --- Lee Clocks Hierarchy ---

class BaseLeeGA(BaseLinearClock):
    """
    Base class for the Lee et al. (2019) Placental Clocks family.
    
    Provides infrastructure to load the three variants of the Lee clock:
    Control, Robust, and Refined Robust.
    
    References:
        Lee Y, et al. Placental epigenetic clocks: estimating gestational age using placental DNA methylation levels
        Aging (Albany NY) (2019). https://doi.org/10.18632/aging.102049
    """
    METADATA = {
        "year": 2019,
        "species": "Human",
        "tissue": "Placenta",
        "omic type": "DNAm(450K/EPIC)",
        "prediction": "Gestational Age(weeks)",
        "source": "https://doi.org/10.18632/aging.102049"
    }
    
    def __init__(self, variant: str):
        """
        Initialize a specific variant of the Lee clock.

        Args:
            variant (str): Internal variant key ('LeeControl', 'LeeRobust', 'LeeRefinedRobust').
        """
        # variant: 'LeeControl', 'LeeRobust', 'LeeRefinedRobust'
        file_map = {
            'LeeControl' : 'Lee_Control', 
            'LeeRobust' : 'Lee_Robust', 
            'LeeRefinedRobust' : 'Lee_RefinedRobust'
        }
        
        try:
            coef_df = load_clock_coefs(file_map.get(variant, variant))
        except Exception as e:
            print(f"[Error] Failed to load {variant}: {e}")
            coef_df = pd.DataFrame(columns=['probe', 'coef'])
            
        super().__init__(coef_df, name=variant, metadata=self.METADATA)

class LeeControl(BaseLeeGA):
    """
    Lee Placental Clock - Control Model (546 CpGs).
    
    Trained specifically on uncomplicated/control pregnancies.
    Best for typical placental samples.
    """
    def __init__(self):
        """Initialize the Lee Control model."""
        super().__init__("LeeControl")

class LeeRobust(BaseLeeGA):
    """
    Lee Placental Clock - Robust Model (558 CpGs).
    
    Trained on a mix of pathological and control samples.
    Recommended for samples with potential pregnancy complications.
    """
    def __init__(self):
        """Initialize the Lee Robust model."""
        super().__init__("LeeRobust")

class LeeRefinedRobust(BaseLeeGA):
    """
    Lee Placental Clock - Refined Robust Model (395 CpGs).
    
    A refined version of the Robust model with fewer CpGs (RPC-like method),
    offering potentially higher stability.
    """
    def __init__(self):
        """Initialize the Lee Refined Robust model."""
        super().__init__("LeeRefinedRobust")