import os
import pandas as pd
import glob
from typing import Optional, Union, List, Dict
from .base import BaseLinearClock

# =========================================================================
# 1. Core Base Class
# =========================================================================
class EnsembleAge:
    """
    Base class for the EnsembleAge framework (Haghani et al. 2025).
    
    This framework aggregates multiple linear sub-models (sub-clocks), each trained 
    on specific perturbations or datasets. Unlike standard clocks that return a 
    single value, this class returns a matrix of predictions.

    Attributes:
        version (str): The specific version of EnsembleAge ('HumanMouse', 'Static', 'Dynamic').
        metadata (dict): Metadata associated with the version.
        sub_clocks (List[BaseLinearClock]): A list of loaded sub-model instances.
    
    References:
        Haghani, A., et al. EnsembleAge: enhancing epigenetic age assessment with a multi-clock framework.
        GeroScience (2025). https://doi.org/10.1007/s11357-025-01808-1
    """

    METADATA_VERSIONS = {
        "HumanMouse": {
            "year": 2025, "species": "Human/Mouse", "tissue": "Multi-tissue",
            "omic type": "DNAm(Mammal40k//EPIC/450k)", "prediction": "Age-to-Lifespan Ratio",
            "source": "https://doi.org/10.1007/s11357-025-01808-1"
        },
        "Static": {
            "year": 2025, "species": "Mouse", "tissue": "Multi-tissue",
            "omic type": "DNAm(Mammal40k/Mammal320k)", "prediction": "Age-to-Lifespan Ratio",
            "source": "https://doi.org/10.1007/s11357-025-01808-1"
        },
        "Dynamic": {
            "year": 2025, "species": "Mouse", "tissue": "Multi-tissue",
            "omic type": "DNAm(Mammal40k/Mammal320k)", "prediction": "Age-to-Lifespan Ratio",
            "source": "https://doi.org/10.1007/s11357-025-01808-1"
        }
    }

    def __init__(self, version: str, data_dir: str = None):
        """
        Initialize the EnsembleAge model by loading multiple sub-clocks from CSV files.

        Args:
            version (str): One of 'HumanMouse', 'Static', or 'Dynamic'.
            data_dir (str, optional): Directory containing coefficient files. 
                If None, defaults to the package's internal 'data/EnsembleAge' directory.
        
        Raises:
            ValueError: If the version is unknown.
            FileNotFoundError: If no coefficient files matching the pattern are found.
        """
        if version not in self.METADATA_VERSIONS:
            raise ValueError(f"Unknown version '{version}'. Supported: {list(self.METADATA_VERSIONS.keys())}")

        self.version = version
        self.metadata = self.METADATA_VERSIONS[version]
        self.sub_clocks = [] 
        
        if data_dir is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            data_dir = os.path.join(base_dir, "data")
        
        # Pattern: EnsembleAge_{version}_*_coefs.csv
        pattern = f"EnsembleAge_{version}_*_coefs.csv"
        search_path = os.path.join(data_dir, "EnsembleAge", pattern)
        files = glob.glob(search_path)
        
        if not files:
            print(f"[Warning] No model files found for EnsembleAge version '{version}' at {search_path}")
            return
            
        print(f"[EnsembleAge] Found {len(files)} sub-clocks for version '{version}'. Loading...")
        
        for file_path in files:
            basename = os.path.basename(file_path)
            try:
                # --- 文件名解析逻辑 ---
                filename_no_ext = os.path.splitext(basename)[0]
                
                if filename_no_ext.endswith("_coefs"):
                    clean_name = filename_no_ext[:-6] 
                else:
                    clean_name = filename_no_ext
                
                prefix = f"EnsembleAge_{version}_"
                if clean_name.startswith(prefix):
                    sub_clock_id = clean_name[len(prefix):] 
                else:
                    sub_clock_id = clean_name

                # Full name for column header
                full_name = f"{version}_{sub_clock_id}"
                
                df = pd.read_csv(file_path)
                sub_meta = self.metadata.copy()
                sub_meta["sub_clock_id"] = sub_clock_id
                
                clock = BaseLinearClock(coef_df=df, name=full_name, metadata=sub_meta)
                self.sub_clocks.append(clock)
                
            except Exception as e:
                print(f"[Error] Failed to load sub-clock {basename}: {e}")

    def get_coefs(self) -> pd.DataFrame:
        """
        Aggregates coefficients from all loaded sub-clocks into a single DataFrame.
        
        Returns:
            pd.DataFrame: A long-format DataFrame with columns ['model_name', 'probe', 'coef'].
            Useful for inspecting which CpGs are used across different sub-models.
        """
        if not self.sub_clocks:
            return pd.DataFrame(columns=['model_name', 'probe', 'coef'])

        dfs = []
        for clock in self.sub_clocks:
            dfs.append(clock.get_coefs())
        return pd.concat(dfs, ignore_index=True)

    def info(self):
        """Prints summary information about the ensemble and loaded sub-clocks."""
        print(f"[EnsembleAge: {self.version}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")
        print(f"  - Loaded sub-clocks: {len(self.sub_clocks)}")

    def get_coefs(self) -> pd.DataFrame:
        """
        Aggregates coefficients from all loaded sub-clocks.
        Ensures 'model_name' is included to distinguish between versions.
        """
        if not self.sub_clocks:
            return pd.DataFrame(columns=['model_name', 'probe', 'coef'])

        dfs = []
        for clock in self.sub_clocks:
            # 1. Obtain the coefficients of the sub-clock
            df = clock.get_coefs().copy()
            
            # 2. Force the addition of the "model_name" column
            df['model_name'] = clock.name
            
            # 3. Sort out the order
            if 'model_name' in df.columns and 'probe' in df.columns and 'coef' in df.columns:
                df = df[['model_name', 'probe', 'coef']]
            
            dfs.append(df)
            
        # 4. Merge all Dataframes
        return pd.concat(dfs, ignore_index=True)

   
    def predict(self, beta_df: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
        """
        Runs prediction on all sub-clocks.
        
        Args:
            beta_df (pd.DataFrame): Input methylation data (Rows: CpGs, Columns: Samples).
            verbose (bool): Whether to print progress for each sub-clock.

        Returns:
            pd.DataFrame: A matrix of predictions.
                - Rows: Samples (same index as beta_df columns).
                - Columns: Sub-clock names (e.g., 'Dynamic_Bmal1', 'Dynamic_GhrKO').
        
        Examples:
            >>> model = EnsembleAgeDynamic()
            >>> results = model.predict(beta_values)
        """
        if not self.sub_clocks:
            return pd.DataFrame()
            
        results = []
        # Make predictions for each sub-clock one by one
        for clock in self.sub_clocks:
            pred = clock.predict(beta_df, verbose=verbose)
            pred.name = clock.name 
            results.append(pred)
            
        if not results:
            return pd.DataFrame()
            
        final_df = pd.concat(results, axis=1)
            
        return final_df

# =========================================================================
# 2. Specialized Subclasses
# =========================================================================
class EnsembleAgeHumanMouse(EnsembleAge):
    """
    EnsembleAge (Human/Mouse) version.
    
    This ensemble is trained on datasets shared between humans and mice,
    capturing conserved aging signatures.
    """
    def __init__(self):
        super().__init__(version="HumanMouse")

class EnsembleAgeStatic(EnsembleAge):
    """
    EnsembleAge (Static) version.
    
    This ensemble consists of clocks trained on static biological differences
    (e.g., comparisons between young and old controls), reflecting baseline aging.
    """
    def __init__(self):
        super().__init__(version="Static")

class EnsembleAgeDynamic(EnsembleAge):
    """
    EnsembleAge (Dynamic) version.
    
    This ensemble consists of clocks trained on perturbation experiments 
    (e.g., caloric restriction, genetic knockouts), specifically designed 
    to be sensitive to anti-aging interventions.
    """
    def __init__(self):
        super().__init__(version="Dynamic")

    def calculate_dynamic_score(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        Wrapper for predict(), specifically for the Dynamic version.
        (Retained for backward compatibility or clarity).
        """
        return self.predict(beta_df)