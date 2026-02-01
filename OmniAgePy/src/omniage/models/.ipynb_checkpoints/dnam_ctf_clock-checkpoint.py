import os
import pandas as pd
import numpy as np
from typing import Union
import xml.etree.ElementTree as ET
try:
    from pypmml import Model
except ImportError:
    Model = None

class DNAmCTFClock:
    """
    DNA Methylation Cell-Type Fraction (CTF) Aging Clock.

    This clock predicts biological age based on **Cell Type Fractions (CTF)** derived 
    from DNA methylation data (e.g., using EpiDISH or similar deconvolution methods).
    
    It utilizes a pre-trained PMML model (Predictive Model Markup Language).

    Attributes:
        name (str): The name of the clock.
        metadata (dict): Metadata associated with the clock.
        model (pypmml.Model): The loaded PMML model instance.
        required_features (list): List of cell types required by the model.

    Note:
        This class requires **Java** to be installed and available in the system PATH, 
        as `pypmml` relies on a Java backend.
    """
    # --- Added Metadata ---
    METADATA = {
        "year": 2026, 
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "Chronological Age(Years)",
        "source": ""
    }

    def __init__(self, model_path: str = None):
        """
        Initialize the DNAm CTF Clock using a pre-trained PMML model.

        Args:
            model_path (str, optional): Path to the .pmml model file. 
                If None, automatically searches for 'dnam_ctf_clock.pmml' in the package data directory.

        Raises:
            ImportError: If the 'pypmml' package is not installed or Java is missing.
            FileNotFoundError: If the model file cannot be located.
        """
        self.name = "DNAm_CTF_Clock"
        self.metadata = self.METADATA # Store metadata
        
        # 1. Dependency Check
        if Model is None:
            raise ImportError(
                "The 'pypmml' package is required for this model but is not installed. "
                "Please run 'poetry add pypmml' and ensure Java is installed."
            )

        self.model = None
        self.required_features = None

        # 2. Automatically locate PMML model file
        if model_path is None:
            # Get the parent directory of the current script (.../src/omniage)
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            # Assumes the R-exported .pmml file is located in src/omniage/data/
            model_path = os.path.join(base_dir, "data", "dnam_ctf_clock.pmml")

        gz_path = model_path + ".gz"

        if not os.path.exists(model_path) and os.path.exists(gz_path):
            print(f"[{self.name}] Decompressing model weights...")
            try:
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(model_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            except Exception as e:
                raise RuntimeError(f"Failed to decompress {gz_path}: {e}")
        
        self.model_path = model_path
        # 3. Load the model
        if os.path.exists(model_path):
            try:
                # Load using pypmml
                self.model = Model.load(model_path)
                
                # Retrieve input field names from PMML for subsequent validation
                self.required_features = self.model.inputNames
            except Exception as e:
                print(f"[Error] Failed to load PMML model. Is Java configured correctly? Details: {e}")
        else:
            print(f"[Error] Model file not found: {model_path}")

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the list of input features (Cell Types) used by the PMML model.
        
        Note:
            Since PMML models may be non-linear (e.g., Random Forest, Neural Net), 
            this method returns `NaN` for the 'coef' column. It is primarily used 
            to inspect which cell types are required.

        Returns:
            pd.DataFrame: A DataFrame with columns ['model_name', 'probe', 'coef'].
        """
        if not os.path.exists(self.model_path):
            return pd.DataFrame()

        try:
            # Parse XML directly to find inputs
            tree = ET.parse(self.model_path)
            root = tree.getroot()
            
            features = []
            
            # Iterate elements to find MiningSchema -> MiningField
            for elem in root.iter():
                if elem.tag.endswith("MiningSchema"):
                    for field in elem:
                        if field.tag.endswith("MiningField"):
                            usage = field.get("usageType", "active")
                            if usage == "active":
                                features.append(field.get("name"))
                    if features:
                        break
            
            if not features:
                return pd.DataFrame()
            
            df = pd.DataFrame({
                'probe': features,
                'coef': np.nan 
            })
            
            df['model_name'] = self.name
            
            return df[['model_name', 'probe', 'coef']]

        except Exception as e:
            print(f"[{self.name}] Error parsing PMML schema: {e}")
            return pd.DataFrame()
    
    def info(self):
        """Prints summary information and metadata about the model."""
        print(f"[{self.name}] Model information:")
        for key, value in self.metadata.items():
            label = key.replace('_', ' ').capitalize()
            print(f"  - {label}: {value}")

    def predict(self, ctf_data: Union[pd.DataFrame, pd.Series]) -> pd.Series:
        """
        Predict age based on cell type fractions.

        Args:
            ctf_data (Union[pd.DataFrame, pd.Series]): Input data containing cell type fractions.
                - If DataFrame: Rows should be Samples, Columns should be Cell Types (e.g., 'CD8Tnv', 'Mono').
                - If Series: Should represent a single sample.

        Returns:
            pd.Series: Predicted biological age values, indexed by sample IDs.

        Raises:
            RuntimeError: If the model is not loaded (e.g., Java error or missing file).
            ValueError: If `ctf_data` contains NaNs or is missing required cell type columns.
            TypeError: If input is not a DataFrame or Series.

        Examples:
            >>> clock = DNAmCTFClock()
            
            # 1. Prepare Cell Type Fractions (e.g., estimated via EpiDISH)
            # Note: The columns MUST match the model's required features exactly.
            # Usually these fractions sum to 1.0 per sample.
            >>> ctf_df = pd.DataFrame({
            ...     'CD4Tnv':  [0.10, 0.09], # Naive CD4+ T
            ...     'CD4Tmem': [0.15, 0.14], # Memory CD4+ T
            ...     'CD8Tnv':  [0.05, 0.06], # Naive CD8+ T
            ...     'CD8Tmem': [0.08, 0.09], # Memory CD8+ T
            ...     'Bnv':     [0.03, 0.02], # Naive B cells
            ...     'Bmem':    [0.02, 0.03], # Memory B cells
            ...     'Treg':    [0.02, 0.01], # Regulatory T cells
            ...     'NK':      [0.05, 0.04], # Natural Killer cells
            ...     'Mono':    [0.10, 0.11], # Monocytes
            ...     'Neu':     [0.35, 0.36], # Neutrophils
            ...     'Eos':     [0.03, 0.04], # Eosinophils
            ...     'Baso':    [0.02, 0.01]  # Basophils
            ... }, index=['Sample1', 'Sample2'])
            
            >>> age = clock.predict(ctf_df)
        """
        if self.model is None:
            raise RuntimeError("Model is not loaded. Cannot perform prediction.")

        # --- Data Preprocessing ---
        # Convert single-sample Series to DataFrame
        if isinstance(ctf_data, pd.Series):
            ctf_data = ctf_data.to_frame().T
        
        if not isinstance(ctf_data, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Feature Validation ---
        if self.required_features:
            # Check for missing required cell type columns
            missing_cols = [c for c in self.required_features if c not in ctf_data.columns]
            if missing_cols:
                raise ValueError(
                    f"[DNAm_CTF_Clock] Input data is missing required cell types:\n"
                    f"{', '.join(missing_cols)}"
                )
            # Retain only required columns and ensure correct order
            data_for_pred = ctf_data[self.required_features]
        else:
            data_for_pred = ctf_data

        # --- NaN Check ---
        if data_for_pred.isna().any().any():
             raise ValueError("[DNAm_CTF_Clock] Input data contains NA values. Please impute missing values first.")

        # --- Execute Prediction ---
        try:
            # pypmml.predict returns a DataFrame
            raw_result = self.model.predict(data_for_pred)
            
            # --- Result Extraction Logic ---
            # PMML regression output columns are typically named 'Predicted_<TargetName>'
            # Attempt to locate the column containing 'predicted' (case-insensitive)
            pred_cols = [c for c in raw_result.columns if "predicted" in c.lower()]
            
            if pred_cols:
                # Select the first identified prediction column
                pred_values = raw_result[pred_cols[0]].values
            else:
                # Fallback: Select the first column if no naming convention matches
                pred_values = raw_result.iloc[:, 0].values

            return pd.Series(pred_values, index=data_for_pred.index, name="Predicted_CTF_Age")
            
        except Exception as e:
            raise RuntimeError(f"Prediction failed inside Java/PMML engine: {e}")