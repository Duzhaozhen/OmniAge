import inspect 
import pandas as pd
import numpy as np
from typing import Union, List, Optional, Dict
import omniage.models as models  
from tqdm import tqdm  


# --- 1. Clock Group Definitions (Configuration) ---
CLOCK_GROUPS = {
    "chronological": [
        "Horvath2013", "Hannum", "Lin", "VidalBralo", "ZhangClock",
        "Horvath2018", "Bernabeu_cAge", "CorticalClock", "PedBE",
        "CentenarianClock_40", "CentenarianClock_100", "Retro_age_V1", "Retro_age_V2",
        "PCHorvath2013","PCHorvath2018","PCHannum",
    ],  
    "biological": [
        "Zhang10", "PhenoAge", "DunedinPACE", 
        "GrimAge1", "GrimAge2","PCGrimAge1", "DNAmFitAge", "IC_Clock","SystemsAge"
    ], 
    "dnamtl": ["DNAmTL","PCDNAmTL"], 
    "cts_clock": ["NeuIn","NeuSin", "GliaIn","GliaSin","Hep"],
    "mitotic": [
        "EpiTOC1", "EpiTOC2", "EpiTOC3", "StemTOCvitro", "StemTOC", 
        "RepliTali", "HypoClock", "EpiCMIT"
    ],
    "causal": ["CausalAge", "DamAge", "AdaptAge"],
    "stochastic":["StocH", "StocZ", "StocP"],
    "trait": ["McCartneyBMI", "McCartneySmoking", "McCartneyAlcohol", "McCartneyEducation", 
                  "McCartneyTotalCholesterol", "McCartneyHDL", "McCartneyLDL", "McCartneyTotalHDLRatio", 
                  "McCartneyWHR", "McCartneyBodyFat"],
    "gestational":["BohlinGA", "EPICGA", "KnightGA","MayneGA","LeeControl","LeeRobust","LeeRefinedRobust"],
    "surrogate_biomarkers":["CompCRP","CompCHIP","EpiSocres"],
    "cross_species": ["PanMammalianUniversal", "PanMammalianBlood", "PanMammalianSkin","EnsembleAge"],
    "pc_clocks": ["PCHorvath2013", "PCHorvath2018", "PCHannum",
                  "PCPhenoAge", "PCDNAmTL","PCGrimAge1"],
}

### cellular_aging
CLOCK_GROUPS["cellular_aging"] = sorted(list(set(
    CLOCK_GROUPS["mitotic"] + CLOCK_GROUPS["dnamtl"]
)))
"""
# "all_epigenetic" 
excluded_groups = {"McCartney", "Gestational", "cellular_aging"} # 也可以把 cellular_aging 排出，因为它只是个合集
all_epi_set = set()

for group_name, clock_list in CLOCK_GROUPS.items():
    if group_name not in excluded_groups:
        all_epi_set.update(clock_list)

CLOCK_GROUPS["all_epigenetic"] = sorted(list(all_epi_set))
"""

def cal_epimarker(
    beta_df: pd.DataFrame,
    clocks: Union[str, List[str]] = "all",
    ages: Optional[Union[pd.Series, list]] = None,
    sex: Optional[Union[pd.Series, list]] = None,
    ctf: Optional[pd.DataFrame] = None,        
    data_type: Optional[str] = None,           
    verbose: bool = True,
    return_dict: bool = True 
) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Computes epigenetic age estimates or biomarker scores for a given DNA methylation matrix.

    This function serves as the primary interface for the package. It dynamically identifies
    available clock models, resolves group aliases (e.g., 'mitotic'), and aggregates predictions
    from multiple models into a single results DataFrame. It also intelligently handles 
    optional arguments (like 'ages') by inspecting model requirements.

    Parameters
    ----------
    beta_df : pd.DataFrame
        A dataframe of DNA methylation beta values.
        Structure: Rows should correspond to CpGs (probes), and columns to Samples.
    
    clocks : Union[str, List[str]], optional
        Specifies which clocks to calculate. 
        - "all": Calculates all available models in the library.
        - "group_name" (e.g., "mitotic"): Calculates all models defined in that group.
        - "ClockName" (e.g., "Horvath2013"): Calculates a specific single model.
        - List[str]: A mixed list of group names and/or specific clock names.
        Default is "all".

    ages : Union[pd.Series, list], optional
        A sequence of chronological ages corresponding to the samples in `beta_df`.
        This is required for certain mitotic clocks (e.g., EpiTOC2, EpiTOC3) to calculate 
        intrinsic stem cell division rates (irS). 
        The order must match the columns (samples) of `beta_df`.

    verbose : bool, optional
        If True, displays progress bars and logging information during execution. 
        Default is True.

    Returns
    -------
    pd.DataFrame
        A combined DataFrame containing the predicted values for all requested clocks.
        Rows correspond to samples (matching `beta_df` columns), and columns correspond 
        to the output metrics of the clocks.
    """
    # --- 1. Build Whitelist from CLOCK_GROUPS ---
    allowed_clocks = set()
    for group_list in CLOCK_GROUPS.values():
        allowed_clocks.update(group_list)
        
    # --- 2. Dynamic Model Discovery ---
    available_clocks = {}
    for name in dir(models):
        obj = getattr(models, name)
        if isinstance(obj, type) and not name.startswith("_") and hasattr(obj, "predict"):
            if name == "BaseLinearClock":
                continue
            if name in allowed_clocks:
                available_clocks[name] = obj
            else:
                pass

    # --- 3. Argument Parsing and Resolution ---
    target_clocks = []

    if clocks == "all":
        target_clocks = list(available_clocks.keys())
    elif isinstance(clocks, str):
        if clocks in CLOCK_GROUPS:
            target_clocks = CLOCK_GROUPS[clocks]
        else:
            target_clocks = [clocks]
    elif isinstance(clocks, list):
        for item in clocks:
            if item in CLOCK_GROUPS:
                target_clocks.extend(CLOCK_GROUPS[item])
            else:
                target_clocks.append(item)
    
    target_clocks = list(dict.fromkeys(target_clocks))

    # --- 4. Validation ---
    valid_clocks = []
    for clk in target_clocks:
        if clk in available_clocks:
            valid_clocks.append(clk)
        else:
            if verbose:
                print(f"Warning: Clock '{clk}' not found in library. Skipping.")

    if "DNAmFitAge" in valid_clocks:
        valid_clocks.remove("DNAmFitAge")
        valid_clocks.append("DNAmFitAge")
        
        # Automatically check whether GrimAge dependencies are included
        grim_variants = ["GrimAge1", "DNAmGrimAge", "GrimAge"]
        has_grimage = any(g in valid_clocks for g in grim_variants)
        
        if not has_grimage:
            # 尝试自动添加 GrimAge1 (如果库里有的话)
            if "GrimAge1" in available_clocks:
                if verbose: print("Note: Auto-adding 'GrimAge1' (required by DNAmFitAge).")
                valid_clocks.insert(0, "GrimAge1")
            elif "DNAmGrimAge" in available_clocks:
                if verbose: print("Note: Auto-adding 'DNAmGrimAge' (required by DNAmFitAge).")
                valid_clocks.insert(0, "DNAmGrimAge")
            else:
                if verbose: print("Warning: DNAmFitAge requested but no GrimAge model found. It may fail.")

    
    if not valid_clocks:
        raise ValueError("No valid clocks selected. Please check the input names.")

    if verbose:
        print(f"Calculating {len(valid_clocks)} clocks: {', '.join(valid_clocks)}")

    # --- Input Standardization & Validation ---

    # 1. Standardize Ages (Must be Numeric)
    if ages is not None:
        # Convert list to Series if necessary to allow pandas operations
        if isinstance(ages, list):
            ages = pd.Series(ages)
        
        # Ensure consistent index with beta_df if not already set
        if isinstance(ages, pd.Series) and not ages.index.equals(beta_df.columns):
             # If indices don't match, assume the order is correct and reset index to match samples
             ages = pd.Series(ages.values, index=beta_df.columns)

        # Validate numeric type
        if not pd.api.types.is_numeric_dtype(ages):
            try:
                # Attempt to coerce strings like "55" to 55
                ages = pd.to_numeric(ages, errors='raise')
            except ValueError:
                raise TypeError(
                    "Input 'ages' must be numeric (int/float). "
                    "Detected non-numeric values that could not be coerced."
                )

    # 2. Standardize Sex (Normalize to "Female"/"Male")
    if sex is not None:
        # Convert list to Series
        if isinstance(sex, list):
            sex = pd.Series(sex)
            
        # Ensure consistent index with beta_df
        if isinstance(sex, pd.Series) and not sex.index.equals(beta_df.columns):
             sex = pd.Series(sex.values, index=beta_df.columns)

        # Pre-check: Numeric values are dangerous (0/1 coding varies), throw error or force string
        if pd.api.types.is_numeric_dtype(sex):
             raise TypeError(
                "Input 'sex' provided as numbers. Please convert to strings ('Female'/'Male') "
                "before calling this function to avoid 0/1 coding ambiguity."
            )

        # --- Normalization Logic ---
        # 1. Convert to string, strip whitespace, convert to lower case
        sex_cleaned = sex.astype(str).str.strip().str.lower()
        
        # 2. Define mapping dictionary (Maps lower-case variations to Standard Format)
        # Note: We map to "Female"/"Male" because downstream models (e.g. GrimAge) 
        # usually check `if x == "Female"`.
        sex_map = {
            'f': 'Female', 'female': 'Female', 'woman': 'Female', 'w': 'Female',
            'm': 'Male', 'male': 'Male', 'man': 'Male',
        }
        
        # 3. Apply mapping
        sex_normalized = sex_cleaned.map(sex_map)
        
        # 4. Check for unmapped values (NaN)
        if sex_normalized.isna().any():
            invalid_values = sex[sex_normalized.isna()].unique()
            raise ValueError(
                f"Input 'sex' contains unrecognized values: {invalid_values}. "
                f"Allowed variants: F, Female, M, Male (case-insensitive)."
            )
            
        # 5. Overwrite the variable with the clean version
        sex = sex_normalized
        
    # 3. Standardize CTF (New)
    if ctf is not None:
        if not ctf.index.equals(beta_df.columns):
            common = beta_df.columns.intersection(ctf.index)
            if len(common) < len(beta_df.columns):
                if verbose: print("Warning: 'ctf' index does not fully match beta_df columns.")

    
    # --- 5. Execution Loop ---
    results_list = []
    results_cache = {}
    
    iterator = valid_clocks
    if verbose:
        try:
            from tqdm import tqdm
            iterator = tqdm(valid_clocks)
        except ImportError:
            pass

    for clock_name in iterator:
        if verbose and hasattr(iterator, "set_description"):
            iterator.set_description(f"Running {clock_name}")

        try:
            # A. Model Class Resolution
            ClockClass = available_clocks[clock_name]
            
            # B. Model Instantiation
            model = ClockClass()
            
            # C. Smart Argument Injection
            # Not all models accept 'ages'. We inspect the signature of the predict method
            # to dynamically construct the arguments dictionary.
            predict_params = {"verbose": False}
            
            # Get the signature of the model's predict method
            sig = inspect.signature(model.predict)
            
            # --- Inject Ages ---
            if "ages" in sig.parameters:
                if ages is not None:
                    predict_params["ages"] = ages
                elif sig.parameters["ages"].default == inspect.Parameter.empty:
                    if verbose: print(f"Skipping {clock_name}: Requires 'ages' input.")
                    continue 

            # --- Inject Sex ---
            if "sex" in sig.parameters:
                if sex is not None:
                    predict_params["sex"] = sex
                elif sig.parameters["sex"].default == inspect.Parameter.empty:
                    if verbose: print(f"Skipping {clock_name}: Requires 'sex' input.")
                    continue
                    
            # --- Inject CTF (Cell Type Fractions) ---
            if "ctf" in sig.parameters:
                if ctf is not None:
                    predict_params["ctf"] = ctf
                elif sig.parameters["ctf"].default == inspect.Parameter.empty:
                    if verbose: print(f"Skipping {clock_name}: Requires 'ctf' (Cell Type Fractions).")
                    continue

            # --- Inject Data Type --- 
            if "data_type" in sig.parameters:
                if data_type is not None:
                    predict_params["data_type"] = data_type
                    
            # --- Inject 'grimage' Dependency ---
            if "grimage" in sig.parameters:
                grim_val = None
                # Candidate clock name (GrimAge1 is the preferred choice)
                candidates = ["GrimAge1", "DNAmGrimAge", "GrimAge"]
                
                for cand in candidates:
                    if cand in results_cache:
                        res = results_cache[cand] # 这是一个 DataFrame
                        target_col = [c for c in res.columns if "DNAmGrimAge" in c]
                        if not target_col:
                            target_col = [c for c in res.columns if "GrimAge" in c]
                        
                        if target_col:
                            col_name = target_col[0]
                            grim_val = res[col_name] 
                            if verbose: print(f"   -> Using dependency: {cand}['{col_name}']")
                        else:
                            
                            grim_val = res.iloc[:, 0]
                        
                        break 
                
                if grim_val is not None:
                    predict_params["grimage"] = grim_val
                else:
                    if verbose: print(f"Skipping {clock_name}: Required 'grimage' input not found in previous results.")
                    continue
            
            # D. Inference Execution
            # Unpack the dynamic parameters dictionary into the function call.
            pred = model.predict(beta_df, **predict_params)
            
            # E. Output Standardization and Type Coercion
            if isinstance(pred, pd.Series):
                pred = pred.to_frame(name=clock_name)
            elif isinstance(pred, np.ndarray):
                cols = [clock_name] if pred.ndim == 1 else [f"{clock_name}_{i}" for i in range(pred.shape[1])]
                pred = pd.DataFrame(pred, index=beta_df.columns, columns=cols)
         
            results_cache[clock_name] = pred 
            results_list.append(pred)
            
        except Exception as e:
            print(f"\n[Error] {clock_name} failed execution: {str(e)}")
            pass

    # --- 6. Result Aggregation ---
    if not results_list:
        return {} if return_dict else pd.DataFrame()

    if return_dict:
        return results_cache

    final_df = pd.concat(results_list, axis=1)
    
    return final_df


