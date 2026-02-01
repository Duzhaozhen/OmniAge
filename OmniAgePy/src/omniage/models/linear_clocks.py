import numpy as np
import pandas as pd
from typing import Optional, Union, List, Dict
import os
from .base import BaseLinearClock
from ..utils import load_clock_coefs



class Horvath2013(BaseLinearClock):
    """
    Horvath 2013 Pan-Tissue Clock implementation.
    Attributes:
        intercept (float): Model offset.
        weights (pd.Series): Coefficients for the 353 CpGs.

    References:
        Horvath S. DNA methylation age of human tissues and cell types. 
        Genome Biol (2013). https://doi.org/10.1186/gb-2013-14-10-r115
    """
    METADATA = {
        "year": 2013,
        "species": "Human",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1186/gb-2013-14-10-r115"
    }

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the Horvath2013 model.
        
        Args:
            coef_df: Optional DataFrame containing model coefficients (CpG IDs and weights).
                If None, the default 'Horvath2013.csv' will be loaded from the package data.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            # Attempts to locate 'Horvath2013.csv' within the internal data directory
            coef_df = load_clock_coefs("Horvath2013")
            
        # 2. Initialize the base class
        super().__init__(coef_df, name="Horvath2013",metadata=self.METADATA)

    def postprocess(self, linear_predictor):
        """
        Implements the inverse transformation (Anti-log).
        """
        adult_age = 20
        
        # Formula for adults (linear)
        adult_transform = linear_predictor * (adult_age + 1) + adult_age
        
        # Formula for childhood (exponential)
        child_transform = np.exp(linear_predictor + np.log(adult_age + 1)) - 1
        
        # Vectorized conditional logic
        return np.where(linear_predictor > 0, adult_transform, child_transform)


class Hannum(BaseLinearClock):
    """
    Hannum 2013 Blood Clock implementation.
    
    References:
        Hannum G, et al. Genome-wide methylation profiles reveal quantitative 
        views of human aging rates. Mol Cell (2013).
        https://doi.org/10.1016/j.molcel.2012.10.016
    """
    METADATA = {
        "year": 2013,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1016/j.molcel.2012.10.016"
    }
    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the Hannum model.
        
        Args:
            coef_df: Optional DataFrame containing model coefficients.
                If None, the default 'Hannum.csv' will be loaded from the package data.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("Hannum")
            
        # 2. Invoke base class initialization
        super().__init__(coef_df, name="Hannum",metadata=self.METADATA)


class Lin(BaseLinearClock):
    """
    Lin 2016 Blood Clock
    
    References:
        Lin Q, et al. DNA methylation levels at individual age-associated CpG sites 
        can be indicative of life expectancy. Aging (2016).
        https://doi.org/10.18632/aging.100908
    """
    METADATA = {
        "year": 2016,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.100908"
    }
        
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the Lin model.
        
        Args:
            coef_df: Optional DataFrame containing model coefficients.
                If None, the default 'Lin.csv' will be loaded from the package data.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("Lin")
            
        # 2. Invoke base class initialization
        super().__init__(coef_df, name="Lin",metadata=self.METADATA)
        


class VidalBralo(BaseLinearClock):
    """
    VidalBralo clock.
    
    References:
        Vidal-Bralo L, et al. Simplified Assay for Epigenetic Age Estimation in Whole Blood of Adults.
        Front Genet (2016). https://doi.org/10.3389/fgene.2016.00126
    """
    METADATA = {
        "year": 2016,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(27k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.3389/fgene.2016.00126"
    }
    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the VidalBralo model.
        
        Args:
            coef_df: Optional DataFrame containing model coefficients.
                If None, the default 'VidalBralo.csv' will be loaded from the package data.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("VidalBralo")
            
        # 2. Invoke base class initialization
        super().__init__(coef_df, name="VidalBralo",metadata=self.METADATA)



class ZhangClock(BaseLinearClock):
    """
    Implements the Zhang (2019) Elastic Net Epigenetic Clock.

    References:
        Zhang Q, et al. Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing
        Genome Med (2019). https://doi.org/10.1186/s13073-019-0667-1
    """
    METADATA = {
        "year": 2019,
        "species": "Human",
        "tissue": "Blood and saliva",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1186/s13073-019-0667-1"
    }

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs("ZhangClock")
        super().__init__(coef_df, name="ZhangClock",metadata=self.METADATA)

    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        Override preprocess to implement Per-Sample Standardization.
        
        Logic:
        For each sample (column), calculate: (x - mean) / std
        """

        # 1. Calculate descriptive statistics for each sample
        # axis=0 computes statistics across the feature set for each column (sample)
        sample_means = beta_df.mean(axis=0)
        sample_stds = beta_df.std(axis=0)
        
        # 2. Execute Z-score transformation
        # Leverages Pandas broadcasting to perform element-wise scaling across columns
        beta_scaled = (beta_df - sample_means) / sample_stds
        
        return beta_scaled



class Horvath2018(BaseLinearClock):
    """
    Horvath 2018 "Skin & Blood" Clock
    
    References:
        Horvath S, et al. DNA methylation-based biomarkers and the epigenetic 
        clock theory of ageing. Nat Rev Genet (2018). 
        https://doi.org/10.1038/s41576-018-0004-3
    """
    METADATA = {
        "year": 2018,
        "species": "Human",
        "tissue": "Skin and blood",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.101508"
    }

    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs("Horvath2018")
        super().__init__(coef_df, name="Horvath2018",metadata=self.METADATA)

    def postprocess(self, linear_predictor):
        """
        Implements the inverse transformation (Anti-log).
        """
        adult_age = 20
        
        # Formula for adults (linear)
        adult_transform = linear_predictor * (adult_age + 1) + adult_age
        
        # Formula for childhood (exponential)
        child_transform = np.exp(linear_predictor + np.log(adult_age + 1)) - 1
        
        # Vectorized conditional logic
        return np.where(linear_predictor > 0, adult_transform, child_transform)




class Bernabeu_cAge(BaseLinearClock):
    """
    Bernabeu Chronological Age Hybrid Model (Blood).
    
    This is a sophisticated **Hybrid Linear-Quadratic** clock designed to handle the 
    non-linear nature of aging across the lifespan.

    Attributes:
        linear_model (dict): Coefficients for the main Linear-Quadratic model.
        log_model (dict): Coefficients for the backup Log-Linear model.
        ref_means (pd.Series): Reference methylation values for imputation.

    References:
        Bernabeu E, et al. Refining epigenetic prediction of chronological and biological age
        Genome Med (2023). https://doi.org/10.1186/s13073-023-01161-y
    """
    METADATA = {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(EPIC)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1186/s13073-023-01161-y"
    }
    
    def __init__(self):
        """
        Initialize the Bernabeu_cAge hybrid model.
        Loads coefficients for both linear and log sub-models and reference means.
        """
        # 1. Load all necessary data
        # Main model (linear) coefficients
        coef_linear = load_clock_coefs("Bernabeu_Linear")
        # Backup model (log) coefficients
        coef_log = load_clock_coefs("Bernabeu_Log")
        # Reference means (for imputation)
        self.ref_means = load_clock_coefs("Bernabeu_Means").set_index('probe')['ref_mean']
        
        # 2. Parse coefficients to separate linear and quadratic terms
        # These need to be stored separately to generate quadratic data during prediction
        self.linear_model = self._parse_coefficients(coef_linear)
        self.log_model = self._parse_coefficients(coef_log)
        
        # 3. Initialize parent class
        # Pass linear coefficients for compatibility, though calculation is handled manually in predict()
        super().__init__(coef_linear, name="Bernabeu_cAge",metadata=self.METADATA)

    def _parse_coefficients(self, df: pd.DataFrame) -> Dict:
        """
        Helper: Splits coefficients into intercept, linear terms, and quadratic terms.
        
        Returns:
            dict: {'intercept': float, 'linear_weights': Series, 'quadratic_weights': Series}
        """
        # Extract intercept
        intercept_mask = df['probe'].isin(['Intercept', '(Intercept)'])
        intercept = float(df[intercept_mask].iloc[0]['coef'])
        weights = df[~intercept_mask].set_index('probe')['coef']
        
        # Separate quadratic terms (names ending with _2)
        quadratic_mask = weights.index.str.endswith('_2')
        quadratic_weights = weights[quadratic_mask]
        linear_weights = weights[~quadratic_mask]
        
        # Remove '_2' suffix from quadratic term names for matching
        # e.g., 'cg0001_2' -> 'cg0001'
        quadratic_weights.index = quadratic_weights.index.str.replace('_2$', '', regex=True)
        
        return {
            'intercept': intercept,
            'linear_weights': linear_weights,
            'quadratic_weights': quadratic_weights,
            'all_cpgs': linear_weights.index.union(quadratic_weights.index)
        }

    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        Implements imputation logic:
        1. Reference imputation for completely missing probes.
        2. Mean imputation for sporadic NAs within present probes.
        """
        # Collect all required CpGs for both models
        all_required = self.linear_model['all_cpgs'].union(self.log_model['all_cpgs'])
        
        # --- 1. Impute completely missing probes ---
        missing_probes = all_required.difference(beta_df.index)
        if not missing_probes.empty:
            print(f"[{self.name}] Imputing {len(missing_probes)} missing probes.")
            # Get reference values (fill with 0.5 if missing in reference table)
            fill_vals = self.ref_means.reindex(missing_probes).fillna(0.5)
            # Create matrix and concatenate
            imputed = pd.DataFrame(
                np.tile(fill_vals.values[:, None], (1, beta_df.shape[1])),
                index=missing_probes,
                columns=beta_df.columns
            )
            beta_df = pd.concat([beta_df, imputed])
            
        # --- 2. Impute sporadic NAs (Row Mean) ---
        if beta_df.isnull().values.any():
            # Fill with row means first
            row_means = beta_df.mean(axis=1)
            beta_df = beta_df.T.fillna(row_means).T
            # If NAs remain (e.g., entire row is NA), use reference means or default to 0.5
            if beta_df.isnull().values.any():
                 beta_df = beta_df.T.fillna(self.ref_means).fillna(0.5).T
                 
        return beta_df

    def _calculate_score(self, beta_df, model_dict):
        """
        Helper: Calculates score = Intercept + (X * w_lin) + (X^2 * w_quad)
        """
        # Linear part
        common_lin = model_dict['linear_weights'].index.intersection(beta_df.index)
        X_lin = beta_df.loc[common_lin].T
        w_lin = model_dict['linear_weights'].loc[common_lin]
        score_lin = X_lin.dot(w_lin)
        
        # Quadratic part
        common_quad = model_dict['quadratic_weights'].index.intersection(beta_df.index)
        X_quad = beta_df.loc[common_quad].T ** 2  # Note the squaring here
        w_quad = model_dict['quadratic_weights'].loc[common_quad]
        score_quad = X_quad.dot(w_quad)
        
        return score_lin + score_quad + model_dict['intercept']

    def predict(self, beta_df: pd.DataFrame, verbose=True) -> pd.Series:
        """
        Executes the hybrid prediction pipeline.
        
        Logic:
        1. Predict age using the Linear-Quadratic model.
        2. Identify samples predicted as < 20 years.
        3. Re-predict those samples using the Log-Linear model.
        
        Warning:
            Input ``beta_df`` MUST be **Beta values (0-1)**. Do not use M-values, 
            as the quadratic term ($X^2$) will produce incorrect results.
        """
        # 1. Preprocessing (Imputation)
        beta_df = self.preprocess(beta_df)
        
        # 2. First Pass: Calculate age using the linear model
        age_linear = self._calculate_score(beta_df, self.linear_model)
        
        # 3. Identify samples < 20 years old
        under_20_mask = age_linear < 20
        under_20_samples = age_linear[under_20_mask].index
        
        final_predictions = age_linear.copy()
        
        if len(under_20_samples) > 0:
            if verbose:
                print(f"[{self.name}] Re-calculating for {len(under_20_samples)} samples < 20 years.")
            
            # 4. Second Pass: Calculate log model for young samples
            # Subset data
            subset_df = beta_df[under_20_samples]
            
            # Calculate log(age)
            log_age = self._calculate_score(subset_df, self.log_model)
            
            # Transform back: age = exp(log_age)
            recalculated_age = np.exp(log_age)
            
            # Update results
            final_predictions.loc[under_20_samples] = recalculated_age
            
        return final_predictions


class CorticalClock(BaseLinearClock):
    """
    Implements the DNAm Cortical Clock (Shireby et al. 2020).

    References:
        Shireby GL, et al. Recalibrating the epigenetic clock: implications for assessing biological age in the human cortex
        Brain (2020). https://doi.org/10.1093/brain/awaa334
    """
    METADATA = {
        "year": 2020,
        "species": "Human",
        "tissue": "Brain cortex",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1093/brain/awaa334"
    }
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        # 1. Automatic loading
        if coef_df is None:
            coef_df = load_clock_coefs("CorticalClock")
            
        # 2. Extract reference values (ref_mean)
        # These must be extracted and stored in self.ref_values before passing to parent
        if 'ref_mean' in coef_df.columns:
            # Extract non-intercept rows
            intercept_mask = coef_df['probe'].isin(['Intercept', '(Intercept)'])
            ref_df = coef_df[~intercept_mask]
            
            # Convert to Series: Index=CpG, Value=RefMean
            self.ref_values = pd.Series(
                ref_df['ref_mean'].values, 
                index=ref_df['probe'].values
            )
        else:
            self.ref_values = None
            print("[CorticalClock] Warning: No 'ref_mean' column found in coefficients.")

        # 3. Initialize parent class
        super().__init__(coef_df, name="CorticalClock", metadata=self.METADATA)

    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        Preprocessing pipeline:
        1. Impute completely missing probes using reference values.
        2. Impute sporadic NAs using row means.
        """
        # --- 1. Impute completely missing probes (Reference Imputation) ---
        # Identify CpGs required by the model but missing in data
        required_probes = self.weights.index
        missing_probes = required_probes.difference(beta_df.index)
        
        # If missing probes exist and reference values are available
        if not missing_probes.empty and self.ref_values is not None:
            print(f"[{self.name}] Imputing {len(missing_probes)} missing probes with reference values.")
            
            # Get reference values for missing probes
            # reindex introduces NaNs for probes not in ref_values; fillna(0) handles this edge case
            fill_values = self.ref_values.reindex(missing_probes).fillna(0)
            
            # Create imputation DataFrame (Rows=MissingCpG, Columns=Samples) using numpy broadcasting
            imputed_data = pd.DataFrame(
                np.tile(fill_values.values[:, None], (1, beta_df.shape[1])),
                index=missing_probes,
                columns=beta_df.columns
            )
            
            # Concatenate original and imputed data
            beta_df = pd.concat([beta_df, imputed_data])

        # --- 2. Impute sporadic NAs (Row Mean Imputation) ---
        # Logic: For each CpG (row), fill NAs with the mean of that row
        if beta_df.isnull().values.any():
            # Transpose, fillna with column means (original row means), then transpose back
            beta_df = beta_df.T.fillna(beta_df.mean(axis=1)).T
            
        return beta_df

    def postprocess(self, linear_predictor):
        """
        Anti-log transformation (Same as Horvath 2013).
        """
        adult_age = 20
        adult_transform = linear_predictor * (adult_age + 1) + adult_age
        child_transform = np.exp(linear_predictor + np.log(adult_age + 1)) - 1
        return np.where(linear_predictor > 0, adult_transform, child_transform)


class PedBE(BaseLinearClock):
    """
    The PedBE (Pediatric Buccal) Clock for DNAm Age in Children
    
    References:
        McEwen LM, O'Donnell KJ, McGill MG, et al.The PedBE clock accurately estimates DNA methylation age in pediatric buccal cells. 
        Proc Natl Acad Sci U S A(2020) https://doi.org/10.1073/pnas.1820843116
    """
    METADATA = {
        "year": 2019,
        "species": "Human",
        "tissue": "Buccal",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1073/pnas.1820843116"
    }
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the PedBE model.
        
        Args:
            coef_df: (Optional) A DataFrame containing coefficients.
                     If None, automatically loads 'PedBE.csv'.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("PedBE")
            
        # 2. Call parent initialization
        super().__init__(coef_df, name="PedBE",metadata=self.METADATA)

    def postprocess(self, linear_predictor):
        """
        Implements the inverse transformation (Anti-log).
        """
        adult_age = 20
        
        # Formula for adults (linear)
        adult_transform = linear_predictor * (adult_age + 1) + adult_age
        
        # Formula for childhood (exponential)
        child_transform = np.exp(linear_predictor + np.log(adult_age + 1)) - 1
        
        # Vectorized conditional logic
        return np.where(linear_predictor > 0, adult_transform, child_transform)


class CentenarianClock_40(BaseLinearClock):
    """
    CentenarianClock (trained on age 40+ cohort) implementation.

    References:
        Dec, E., et al.Centenarian clocks: epigenetic clocks for validating claims of exceptional longevity.
        GeroScience (2020) https://doi.org/10.1007/s11357-023-00731-7
    """
    METADATA = {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood and saliva",
        "omic type": "DNAm(450k and EPIC)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1007/s11357-023-00731-7"
    }
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the CentenarianClock_40 model.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("CentenarianClock_40")
            
        # 2. Call parent initialization
        super().__init__(coef_df, name="CentenarianClock_40",metadata=self.METADATA)


class CentenarianClock_100(BaseLinearClock):
    """
    CentenarianClock (trained on age 100+ cohort) implementation.

    References:
        Dec, E., et al.Centenarian clocks: epigenetic clocks for validating claims of exceptional longevity.
        GeroScience (2020) https://doi.org/10.1073/pnas.1820843116
    """
    METADATA = {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood and saliva",
        "omic type": "DNAm(450k and EPIC)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1007/s11357-023-00731-7"
    }

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the CentenarianClock_100 model.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("CentenarianClock_100")
            
        # 2. Call parent initialization
        super().__init__(coef_df, name="CentenarianClock_100",metadata=self.METADATA)


class Retro_age_V1(BaseLinearClock):
    """
    Retro_age V1 implementation.
    A retroelement-based aging clock.

    References:
         Ndhlovu LC, et al. Retro-age: A unique epigenetic biomarker of aging captured by DNA methylation states of retroelements.
         Aging Cell (2020) https://doi.org/10.1111/acel.14288
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(EPIC)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1111/acel.14288"
    }
    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the Retro_age_V1 model.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("Retro_age_V1")
            
        # 2. Call parent initialization
        super().__init__(coef_df, name="Retro_age_V1",metadata=self.METADATA)


class Retro_age_V2(BaseLinearClock):
    """
    Retro_age V2 implementation.
    A retroelement-based aging clock.

    References:
         Ndhlovu LC et al. Retro-age: A unique epigenetic biomarker of aging captured by DNA methylation states of retroelements.
         Aging Cell (2020) https://doi.org/10.1111/acel.14288
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(EPIC)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.1111/acel.14288"
    }
    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the Retro_age_V2 model.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("Retro_age_V2")
            
        # 2. Call parent initialization
        super().__init__(coef_df, name="Retro_age_V2",metadata=self.METADATA)


class Zhang10(BaseLinearClock):
    """
    Zhang 10-CpG Mortality Clock implementation.

    References:
         Zhang, et al. DNA methylation signatures in peripheral blood strongly predict all-cause mortality.
         Nat Commun (2020) https://doi.org/10.1038/ncomms14617
    """
    METADATA = {
        "year": 2017,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Mortality",
        "units": "Years",
        "source": "https://doi.org/10.1038/ncomms14617"
    }

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the Zhang10 model.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("Zhang10")
            
        # 2. Call parent initialization
        super().__init__(coef_df, name="Zhang10",metadata=self.METADATA)


class PhenoAge(BaseLinearClock):
    """
    PhenoAge (Levine et al. 2018) implementation.

    References:
         Levine et al. An epigenetic biomarker of aging for lifespan and healthspan.
         Aging (2018) https://doi.org/10.18632/aging.101414
    """
    METADATA = {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Mortality",
        "source": "https://doi.org/10.18632/aging.101414"
    }

    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize the PhenoAge model.
        """
        # 1. Automatic loading logic
        if coef_df is None:
            coef_df = load_clock_coefs("PhenoAge")
            
        # 2. Call parent initialization
        super().__init__(coef_df, name="PhenoAge",metadata=self.METADATA)


class EpiTOC1(BaseLinearClock):
    """
    Implements epiTOC1 (Yang et al. 2016).

    References:
         Yang et al. Correlation of an epigenetic mitotic clock with cancer risk.
         Genome Biol (2016) https://doi.org/10.1186/s13059-016-1064-3
    """
    METADATA = {
        "year": 2016,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Mitotic Divisions",
        "source": "https://doi.org/10.1186/s13059-016-1064-3"
    }

    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """Initialize EpiTOC1."""
        # 1. Automatic loading
        if coef_df is None:
            coef_df = load_clock_coefs("epiTOC1")
            
        # 2. Initialize parent class
        # Although unweighted, the parent class parses the probe list into self.weights.index
        super().__init__(coef_df, name="epiTOC1",metadata=self.METADATA)

    def predict(self, beta_df: pd.DataFrame, verbose=True) -> pd.Series:
        """
        Calculates the EpiTOC1 score (Average Methylation).
        
        Args:
            beta_df (pd.DataFrame): Methylation data (Rows: CpGs, Columns: Samples).
            verbose (bool): Print feature coverage stats.
            
        Returns:
            pd.Series: The EpiTOC1 score.
        """
        # 1. Preprocessing
        beta_df = self.preprocess(beta_df)
        
        # 2. Feature alignment (Retain only the required 385 CpGs)
        required_cpgs = self.weights.index
        common_features = required_cpgs.intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Number of represented CpGs (max={len(required_cpgs)})={len(common_features)}")
            
        if len(common_features) == 0:
            print(f"[{self.name}] Error: No overlapping CpGs found.")
            return pd.Series([np.nan] * beta_df.shape[1], index=beta_df.columns)

        # 3. Calculation Logic: Column Means
        # Equivalent to R: colMeans(data.m[map.idx,], na.rm=TRUE)
        
        # Subset data (Rows=CpGs, Columns=Samples)
        subset_data = beta_df.loc[common_features, :]
        
        # Calculate mean for each column (sample)
        # skipna=True corresponds to R's na.rm=TRUE
        scores = subset_data.mean(axis=0, skipna=True)
        
        return scores


class EpiTOC2(BaseLinearClock):
    """
    Implements epiTOC2 (Teschendorff et al. 2020).
    
    References:
         Teschendorff et al. Correlation of an epigenetic mitotic clock with cancer risk.
         Genome Biol (2016) https://doi.org/10.1186/s13059-016-1064-3
    """
    METADATA = {
        "year": 2020,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Mitotic Divisions",
        "source": "https://doi.org/10.1186/s13073-020-00752-3"
    }
    
    def __init__(self):
        """Initialize EpiTOC2 and load parameters (delta, beta0)."""
        # 1. Define model name
        self.name = "EpiTOC2"
        
        # 2. Automatically load coefficient data
        df = load_clock_coefs("epiTOC2")
        
        # 3. Parse epiTOC2-specific parameters (slope delta and baseline beta0)
        self.delta = df.set_index('probe')['delta']
        self.beta0 = df.set_index('probe')['beta0']
        
        # 4. Initialize parent class
        super().__init__(pd.DataFrame(columns=['probe', 'coef']), self.name, self.METADATA)
        
        # 5. Explicitly override attributes for parent-class compatibility
        self.weights = self.delta  
        self.intercept = 0.0
        
    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the specific parameters for EpiTOC2 (delta and beta0).
        Overrides BaseLinearClock.get_coefs.
        """
        # Reconstruct DataFrame
        df = pd.DataFrame({
            'probe': self.delta.index,
            'delta': self.delta.values,  
            'beta0': self.beta0.values   
        })
        return df

    def predict(self, beta_df: pd.DataFrame, ages: Optional[Union[list, pd.Series]] = None, verbose=True) -> pd.DataFrame:
        """
        Calculates epiTOC2 scores.
        
        Args:
            beta_df: CpG x Sample matrix
            ages: (Optional) Chronological age vector/series. 
                  If provided, calculates Intrinsic Rate (irS).
                  
        Returns:
            pd.DataFrame with columns: ['tnsc', 'tnsc2', 'irS', 'irS2']
        """
        # 1. Preprocessing
        beta_df = self.preprocess(beta_df)
        
        # 2. Feature alignment
        common_cpgs = self.delta.index.intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Number of represented CpGs (max={len(self.delta)})={len(common_cpgs)}")
            
        if len(common_cpgs) == 0:
            return pd.DataFrame(np.nan, index=beta_df.columns, columns=['tnsc', 'tnsc2'])

        # 3. Prepare calculation data
        # Subset and align
        beta_subset = beta_df.loc[common_cpgs, :]
        delta_subset = self.delta.loc[common_cpgs]
        beta0_subset = self.beta0.loc[common_cpgs]
        
        # 4. Core Calculation (Vectorized Implementation)
        
        # --- Calculate TNSC (Full Model) ---
        # R Formula: 2 * colMeans( diag(1/(delta*(1-beta0))) %*% (beta - beta0) )
        
        term_a = delta_subset * (1 - beta0_subset)
        
        # Pandas sub(axis=0) subtracts the vector from each column
        diff_matrix = beta_subset.sub(beta0_subset, axis=0)
        
        # div(axis=0) divides each column by the vector
        # Multiply by 2 and take column means
        tnsc = 2 * diff_matrix.div(term_a, axis=0).mean(axis=0)
        
        # --- Calculate TNSC2 (Approximation Model) ---
        # R Formula: 2 * colMeans( diag(1/delta) %*% beta )
        # Assumes beta0 approx 0
        tnsc2 = 2 * beta_subset.div(delta_subset, axis=0).mean(axis=0)
        
        # 5. Consolidate results
        results = pd.DataFrame({
            'tnsc': tnsc,
            'tnsc2': tnsc2
        })
        
        # 6. Calculate Division Rate (irS) if age is provided
        if ages is not None:
            if isinstance(ages, list):
                ages = pd.Series(ages, index=results.index)
            
            # Simple division
            results['irS'] = results['tnsc'] / ages
            results['irS2'] = results['tnsc2'] / ages
            
            # Calculate Tissue Median Rate (irT)
            irT = results['irS'].median()
            irT2 = results['irS2'].median() # Defined for completeness
            
            if verbose:
                print(f"[{self.name}] Tissue Median Rate (irT): {irT:.4f}")
                print(f"[{self.name}] Tissue Median Rate (irT2 - approximation): {irT2:.4f}")
                
        results = results.add_prefix(f"{self.name}_")
        return results


class EpiTOC3(BaseLinearClock):
    """
    Implements epiTOC3 (Teschendorff et al. 2020).
    """
    METADATA = {
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "Mitotic Divisions"
    }
    
    def __init__(self):
        """Initialize EpiTOC3."""
        # 1. Model identification
        self.name = "EpiTOC3"
        
        # 2. Automatic parameter loading
        df = load_clock_coefs("epiTOC3")
        
        # 3. Parse epiTOC3-specific parameters
        df_indexed = df.set_index('probe')
        self.delta = df_indexed['delta']
        self.beta0 = df_indexed['beta0']
        
        # 4. Parent class initialization
        super().__init__(pd.DataFrame(columns=['probe', 'coef']), self.name, self.METADATA)
        
        # 5. Explicit attribute override for framework compatibility
        self.weights = self.delta 
        self.intercept = 0.0

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the specific parameters for EpiTOC3 (delta and beta0).
        Overrides BaseLinearClock.get_coefs.
        """
        df = pd.DataFrame({
            'probe': self.delta.index,
            'delta': self.delta.values,
            'beta0': self.beta0.values
        })
        return df
        
    def predict(self, beta_df: pd.DataFrame, ages: Optional[Union[list, pd.Series]] = None, verbose=True) -> pd.DataFrame:
        """
        Calculates EpiTOC3 scores (TNSC, irS, and Average Methylation).

        Args:
            beta_df (pd.DataFrame): Methylation matrix.
            ages (list/Series, optional): Chronological age (for irS calculation).

        Returns:
            pd.DataFrame: Columns ['EpiTOC3_tnsc', 'EpiTOC3_avETOC3', 'EpiTOC3_irS', ...].
        """
        # 1. Pre-Processing 
        beta_df = self.preprocess(beta_df)
        
        # 2. Feature alignment
        common_cpgs = self.delta.index.intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Number of represented CpGs (max={len(self.delta)})={len(common_cpgs)}")
            
        if len(common_cpgs) == 0:
            cols = ['tnsc', 'tnsc2', 'avETOC3']
            return pd.DataFrame(np.nan, index=beta_df.columns, columns=cols)

        # 3. Prepare subset data
        beta_subset = beta_df.loc[common_cpgs, :]
        delta_subset = self.delta.loc[common_cpgs]
        beta0_subset = self.beta0.loc[common_cpgs]
        
        # 4. Begin calculation
        
        # --- A. avETOC3 ---
        # R: colMeans(tmp.m)
        avETOC3 = beta_subset.mean(axis=0)
        
        # --- B. TNSC (Full Model) ---
        # Formula: 2 * mean( (beta - beta0) / (delta * (1 - beta0)) )
        term_a = delta_subset * (1 - beta0_subset)
        diff_matrix = beta_subset.sub(beta0_subset, axis=0)
        tnsc = 2 * diff_matrix.div(term_a, axis=0).mean(axis=0)
        
        # --- C. TNSC2 (Approximation) ---
        # Formula: 2 * mean( beta / delta )
        tnsc2 = 2 * beta_subset.div(delta_subset, axis=0).mean(axis=0)
        
        # 5. Consolidate results
        results = pd.DataFrame({
            'tnsc': tnsc,
            'tnsc2': tnsc2,
            'avETOC3': avETOC3
        })
        
        # 6. Calculate Rates (irS, irT) - if age provided
        if ages is not None:
            if isinstance(ages, list):
                ages = pd.Series(ages, index=results.index)
            
            # Calculate sample-level rates (irS)
            results['irS'] = results['tnsc'] / ages
            results['irS2'] = results['tnsc2'] / ages
            
            # Calculate tissue-level median rate (irT)
            irT = results['irS'].median()
            irT2 = results['irS2'].median()
            
            # Store in metadata
            results.attrs['irT'] = irT
            results.attrs['irT2'] = irT2
            
            if verbose:
                print(f"[{self.name}] Median Rate (irT): {irT:.4f}")
                print(f"[{self.name}] Median Rate (irT2): {irT2:.4f}")
        
        results = results.add_prefix(f"{self.name}_")
        return results


class HypoClock(BaseLinearClock):
    """
    Implements HypoClock.

    References:
        Teschendorff AE, et al. A comparison of epigenetic mitotic-like clocks for cancer risk prediction
        Genome Med (2020). https://doi.org/10.1186/s13073-020-00752-3
    """
    
    METADATA = {
        "year": 2020,
        "species": "Human",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(450k)",
        "prediction": "Mitotic Divisions",
        "source": "https://doi.org/10.1186/s13073-020-00752-3"
    }
    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """Initialize HypoClock."""
        # 1. Automatic loading
        if coef_df is None:
            coef_df = load_clock_coefs("HypoClock")
            
        # 2. Initialize parent class
        # Parses the 678 CpGs into self.weights.index
        super().__init__(coef_df, name="HypoClock", metadata=self.METADATA)

    def predict(self, beta_df: pd.DataFrame, verbose=True) -> pd.Series:
        """
        Override predict logic: Score = 1 - Mean(Beta)
        """
        # 1. Preprocessing
        beta_df = self.preprocess(beta_df)
        
        # 2. Feature alignment
        required_cpgs = self.weights.index
        common_features = required_cpgs.intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Number of represented solo-WCGWs (max={len(required_cpgs)})={len(common_features)}")
            
        if len(common_features) == 0:
            print(f"[{self.name}] Error: No overlapping CpGs found.")
            return pd.Series([np.nan] * beta_df.shape[1], index=beta_df.columns)

        # 3. Calculation logic
        # R code: 1 - colMeans(data.m[map.idx,], na.rm=TRUE)
        
        # Subset data
        subset_data = beta_df.loc[common_features, :]
        
        # Calculate mean (column/sample-wise)
        mean_beta = subset_data.mean(axis=0, skipna=True)
        
        # Calculate final score
        hypo_score = 1 - mean_beta
        
        return hypo_score  


class RepliTali(BaseLinearClock):
    """
    Implements RepliTali (Endicott et al. 2022).

    References:
        Endicott RM, et al. Cell division drives DNA methylation loss in late-replicating domains in primary human cells
        Nat Commun (2022). https://doi.org/10.1038/s41467-022-34268-8
    """

    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(EPIC)",
        "prediction": "Mitotic Divisions",
        "source": "https://doi.org/10.1038/s41467-022-34268-8"
    }
    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        """Initialize RepliTali (Auto-loads 86 CpGs)."""
        # 1. Automatic loading
        if coef_df is None:
            coef_df = load_clock_coefs("RepliTali")
            
        # 2. Initialize parent class
        # Automatically identifies "Intercept" row and sets remaining 86 CpGs as weights
        super().__init__(coef_df, name="RepliTali",metadata=self.METADATA)
        
    # No need to override predict; parent matrix multiplication logic applies:
    # Age = Intercept + (Beta * Weights)


class EpiCMIT_Hyper(BaseLinearClock):
    """
    Implements EpiCMIT Hyper Score (Duran-Ferrer et al. 2020)
    Calculated as the mean methylation of specific hyper-methylated CpGs.

    References:
        Endicott RM, et al. Cell division drives DNA methylation loss in late-replicating domains in primary human cells
        Nat Commun (2022). https://doi.org/10.1038/s41467-022-34268-8
    """
    
    METADATA = {
        "year": 2020,
        "species": "Human",
        "tissue": "B cells",
        "omic type": "DNAm(450k)",
        "prediction": "Mitotic Divisions (Hyper)",
        "source": "https://doi.org/10.1038/s43018-020-00131-2"
    }

    def __init__(self):
        self.name = "EpiCMIT_Hyper"
        
        # 1. Load Data
        df = load_clock_coefs("EpiCMIT")
        
        # 2. Filter for Hyper probes only
        self.feature_cpgs = df[df['type'].astype(str).str.contains("hyper", case=False)]['probe'].values
        
        # 3. Parent class initialization
        super().__init__(pd.DataFrame(columns=['probe', 'coef']), self.name, self.METADATA)
        
        # 4. Compatibility attributes
        self.weights = pd.Series(1, index=self.feature_cpgs)
        self.intercept = 0.0

    def predict(self, beta_df: pd.DataFrame, verbose=True) -> pd.Series:
        """
        Calculates EpiCMIT Hyper score (Mean Methylation).
        """
        # 1. Preprocessing
        beta_df = self.preprocess(beta_df)
        
        # 2. Identify common CpGs
        common = pd.Index(self.feature_cpgs).intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Represented CpGs (max={len(self.feature_cpgs)})={len(common)}")
            
        # 3. Calculate Score: Mean of Beta Values
        if len(common) > 0:
            scores = beta_df.loc[common, :].mean(axis=0, skipna=True)
        else:
            scores = pd.Series(np.nan, index=beta_df.columns)
            
        return scores


class EpiCMIT_Hypo(BaseLinearClock):
    """
    Implements EpiCMIT Hypo Score (Duran-Ferrer et al. 2020)
    Calculated as 1 - (mean methylation of specific hypo-methylated CpGs).

    References:
        Endicott RM, et al. Cell division drives DNA methylation loss in late-replicating domains in primary human cells
        Nat Commun (2022). https://doi.org/10.1038/s41467-022-34268-8
    """
    
    METADATA = {
        "year": 2020,
        "species": "Human",
        "tissue": "B cells",
        "omic type": "DNAm(450k)",
        "prediction": "Mitotic Divisions (Hypo)",
        "source": "https://doi.org/10.1038/s43018-020-00131-2"
    }

    def __init__(self):
        self.name = "EpiCMIT_Hypo"
        
        # 1. Load Data
        df = load_clock_coefs("EpiCMIT")
        
        # 2. Filter for Hypo probes only
        self.feature_cpgs = df[df['type'].astype(str).str.contains("hypo", case=False)]['probe'].values
        
        # 3. Parent class initialization
        super().__init__(pd.DataFrame(columns=['probe', 'coef']), self.name, self.METADATA)
        
        # 4. Compatibility attributes
        self.weights = pd.Series(1, index=self.feature_cpgs)
        self.intercept = 0.0

    def predict(self, beta_df: pd.DataFrame, verbose=True) -> pd.Series:
        """
        Calculates EpiCMIT Hypo score (1 - Mean Methylation).
        """
        # 1. Preprocessing
        beta_df = self.preprocess(beta_df)
        
        # 2. Identify common CpGs
        common = pd.Index(self.feature_cpgs).intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Represented CpGs (max={len(self.feature_cpgs)})={len(common)}")
            
        # 3. Calculate Score: 1 - Mean of Beta Values
        if len(common) > 0:
            scores = 1 - beta_df.loc[common, :].mean(axis=0, skipna=True)
        else:
            scores = pd.Series(np.nan, index=beta_df.columns)
            
        return scores




class StemTOC(BaseLinearClock):
    """
    Implements StemTOC (Zhu et al. 2024).

    References:
        Zhu T, et al. An improved epigenetic counter to track mitotic age in normal and precancerous tissues
        Nat Commun (2024). https://doi.org/10.1038/s41467-024-48649-8
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "Mitotic Divisions",
        "source": "https://doi.org/10.1038/s41467-024-48649-8"
    }
    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        # 1. Automatic loading
        if coef_df is None:
            coef_df = load_clock_coefs("StemTOC")
            
        # 2. Initialize parent class
        # Parses 371 CpGs into self.weights.index
        super().__init__(coef_df, name="StemTOC", metadata=self.METADATA)

    def predict(self, beta_df: pd.DataFrame, verbose=True) -> pd.Series:
        """
        Override predict logic: Score = 95th Percentile of Beta Values
        """
        # 1. Preprocessing
        beta_df = self.preprocess(beta_df)
        
        # 2. Feature alignment
        required_cpgs = self.weights.index
        common_features = required_cpgs.intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Number of represented CpGs (max={len(required_cpgs)})={len(common_features)}")
            
        if len(common_features) == 0:
            print(f"[{self.name}] Error: No overlapping CpGs found.")
            return pd.Series([np.nan] * beta_df.shape[1], index=beta_df.columns)

        # 3. Calculation Logic
        # R code: apply(data, 2, quantile, 0.95, na.rm=T)
        # Python: quantile(q=0.95, axis=0) 
        
        # Subset data
        subset_data = beta_df.loc[common_features, :]
        
        # Calculate 0.95 quantile for each column (sample)
        # Pandas quantile implicitly ignores NaNs
        scores = subset_data.quantile(q=0.95, axis=0)
        
        return scores


class StemTOCvitro(BaseLinearClock):
    """
    Implements StemTOCvitro (Zhu et al. 2024).

    References:
        Zhu T, et al. An improved epigenetic counter to track mitotic age in normal and precancerous tissues
        Nat Commun (2024). https://doi.org/10.1038/s41467-024-48649-8
    """

    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "Mitotic Divisions",
        "source": "https://doi.org/10.1038/s41467-024-48649-8"
    }
    
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        # 1. Automatic loading
        if coef_df is None:
            coef_df = load_clock_coefs("StemTOCvitro")
            
        # 2. Initialize parent class
        # Parses 629 CpGs into self.weights.index
        super().__init__(coef_df, name="StemTOCvitro", metadata=self.METADATA)

    def predict(self, beta_df: pd.DataFrame, verbose=True) -> pd.Series:
        """
        Override predict logic: Score = 95th Percentile of Beta Values
        """
        # 1. Preprocessing
        beta_df = self.preprocess(beta_df)
        
        # 2. Feature alignment
        required_cpgs = self.weights.index
        common_features = required_cpgs.intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Number of represented CpGs (max={len(required_cpgs)})={len(common_features)}")
            
        if len(common_features) == 0:
            print(f"[{self.name}] Error: No overlapping CpGs found.")
            return pd.Series([np.nan] * beta_df.shape[1], index=beta_df.columns)

        # 3. Calculation Logic
        # R code: apply(data, 2, quantile, 0.95, na.rm=T)
        # Python: quantile(q=0.95, axis=0)
        
        # Subset data
        subset_data = beta_df.loc[common_features, :]
        
        # Calculate 0.95 quantile for each column (sample)
        scores = subset_data.quantile(q=0.95, axis=0)
        
        return scores



class DNAmTL(BaseLinearClock):
    """
    Implements DNAmTL (Lu et al. 2019).

    References:
        Lu AT, et al. DNA methylation-based estimator of telomere length.
        Aging (2019). https://doi.org/10.18632/aging.102173
    """
    METADATA = {
        "year": 2019,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "Telomere Length",
        "source": "https://doi.org/10.18632/aging.102173"
    }
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs("DNAmTL")

        super().__init__(coef_df, name="DNAmTL", metadata=self.METADATA)




class DunedinPACE(BaseLinearClock):
    """
    Implements DunedinPACE (Belsky et al. 2022).
    
    References:
        Belsky DW, et al. DunedinPACE, a DNA methylation biomarker of the pace of aging. 
        eLife (2022). https://doi.org/10.7554/elife.73420
    """
    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "prediction": "Pace of ageing",
        "source": "https://doi.org/10.7554/elife.73420"
    }
    def __init__(self):
        # 1. Load coefficients (including model_mean column)
        coef_df = load_clock_coefs("DunedinPACE")
        
        # --- Extract model means ---
        # Filter out Intercept rows to isolate probes
        cpg_df = coef_df[~coef_df['probe'].isin(['Intercept', '(Intercept)'])]
        self.model_means = cpg_df.set_index('probe')['model_mean']
        
        # 2. Load Gold Standard Means (for normalization)
        gs_df = load_clock_coefs("DunedinPACE_GoldStandard")
        self.gold_standard_means = gs_df.set_index('probe')['ref_mean']
        
        # 3. Initialize parent class
        super().__init__(coef_df, name="DunedinPACE", metadata=self.METADATA)

    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        DunedinPACE Preprocessing Pipeline:
        1. Impute missing probes using Gold Standard Means.
        2. Impute internal NAs using Model Means.
        3. Quantile Normalization against Gold Standard Target.
        """
        # --- 1. Impute completely missing probes (Gold Standard Means) ---
        # Logic: Some probes may not be in the model but are required for normalization background.
        required_gs_probes = self.gold_standard_means.index
        missing_gs_probes = required_gs_probes.difference(beta_df.index)
        
        if not missing_gs_probes.empty:
            fill_values = self.gold_standard_means.loc[missing_gs_probes]
            # Create imputation matrix using numpy broadcasting
            imputed_data = pd.DataFrame(
                np.tile(fill_values.values[:, None], (1, beta_df.shape[1])),
                index=missing_gs_probes,
                columns=beta_df.columns
            )
            beta_df = pd.concat([beta_df, imputed_data])
            
        # At this point, beta_df contains all probes required for normalization.
        # Align with Gold Standard order.
        beta_df = beta_df.loc[required_gs_probes]

        # --- 2. Impute internal NAs (Model Means) ---
        if beta_df.isnull().values.any():
            # Strategy A: First, fill with sample-wise (row) means
            row_means = beta_df.mean(axis=1)
            beta_df = beta_df.T.fillna(row_means).T
            
            # Strategy B: If NAs persist (entire row is NA), fill with Model Means
            # Note: self.model_means contains only ~173 model probes, while Gold Standard has ~20k.
            # We create a fallback series combining Model Means (priority) and Gold Standard Means.
            
            fallback_means = self.gold_standard_means.copy()
            fallback_means.update(self.model_means) # Overwrite with model means where available
            
            beta_df = beta_df.T.fillna(fallback_means).T

        # --- 3. Quantile Normalization ---
        # Sort the gold standard reference distribution
        target_sorted = np.sort(self.gold_standard_means.values)
        
        def quant_norm(col):
            # Get ranks (0 to N-1). 
            # Note: argsort().argsort() returns the rank of each element.
            rank = col.argsort().argsort()
            # Replace value with the value at the same rank in the target distribution
            return target_sorted[rank]
            
        # Apply per sample (column)
        beta_normalized = beta_df.apply(quant_norm, axis=0)
        beta_normalized.index = beta_df.index
        
        return beta_normalized


class GrimAge1(BaseLinearClock):
    """
    Implements DNAm GrimAge1 (Lu et al. 2019).

    A composite biomarker of mortality risk. GrimAge is unique because it is built 
    in two stages:
    1.  **Surrogates**: DNAm-based predictors are trained for plasma proteins (e.g., ADM, GDF15) 
        and smoking history (PACKYRS).
    2.  **Composite**: These surrogates, along with Age and Sex, are combined to predict 
        time-to-death (Mortality Risk).

    References:
        Lu AT, et al. DNA methylation GrimAge strongly predicts lifespan and healthspan. 
        Aging (Albany NY) (2019). https://doi.org/10.18632/aging.101684
    """

    METADATA = {
        "year": 2019,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Mortality",
        "source": "https://doi.org/10.18632/aging.101684"
    }
    
    def __init__(self):
        # 1. Model identification
        self.name = "GrimAge1"
        
        # 2. Load model components
        self.coefs = load_clock_coefs("GrimAge1_Coefs")
        self.gold = load_clock_coefs("GrimAge1_Gold")
        
        # 3. Identify feature requirements
        non_cpg_vars = ["Intercept", "Age", "Female"]
        all_vars = self.coefs['var'].unique()
        self.required_cpgs = [v for v in all_vars if v not in non_cpg_vars]
        
        # 4. Parent class initialization
        super().__init__(pd.DataFrame(columns=['probe', 'coef']), self.name, self.METADATA)

        # 5. Compatibility placeholders
        self.weights = pd.Series(1, index=self.required_cpgs)
        self.intercept = 0.0

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the composite coefficients for GrimAge1.
        Returns a DataFrame including the 'Y.pred' column to distinguish 
        between different surrogate models (e.g., DNAmADM, DNAmPACKYRS).
        """
        df = self.coefs.copy()
        
        # var -> probe
        # beta -> coef
        if 'var' in df.columns and 'beta' in df.columns:
            df = df.rename(columns={'var': 'probe', 'beta': 'coef'})
            
        if 'Y.pred' in df.columns:
            df = df[['Y.pred', 'probe', 'coef']]
            
        return df
    
    def predict(
        self, 
        beta_df: pd.DataFrame, 
        ages: Union[list, pd.Series], 
        sex: Union[list, pd.Series], 
        verbose: bool = True
    ) -> pd.DataFrame:
        """
        Calculates GrimAge and its components.

        Args:
            beta_df (pd.DataFrame): Methylation data (Rows: CpGs, Columns: Samples).
            ages (list/Series): **Chronological Age**. Required for the model.
            sex (list/Series): **Biological Sex**. 
                Format: 'Female'/'Male' or 'F'/'M' (Case insensitive).
            verbose (bool): Print warnings about unrecognized sex labels.

        Returns:
            pd.DataFrame: A table containing:
            - **DNAmGrimAge1**: The final biological age estimate.
            - **DNAmADM, DNAmGDF15, DNAmPACKYRS...**: Individual surrogate protein scores.
            - **raw_cox**: The raw mortality score before scaling.
        """
        # --- 1. Preprocessing and Data Preparation ---
        beta_df = self.preprocess(beta_df)
        
        # Format Age and Sex
        if isinstance(ages, list):
            ages = pd.Series(ages, index=beta_df.columns)
        if isinstance(sex, list):
            sex = pd.Series(sex, index=beta_df.columns)

        # --- Validate Input Type ---
        # Strictly reject numeric inputs to prevent ambiguity (e.g., 0=Male vs 0=Female).
        if pd.api.types.is_numeric_dtype(sex):
            raise ValueError(
                "Ambiguous sex input: Numeric values (0/1) are not allowed. "
                "Please provide clear strings, e.g., 'Female'/'F' or 'Male'/'M'."
            )
        # --- Standardize Format ---
        sex_str = sex.astype(str).str.strip().str.upper()

        # --- Binary Encoding ---
        # Matches become 1 (Female), all others default to 0 (Male).
        female_labels = ["FEMALE", "F"]
        is_female = sex_str.isin(female_labels).astype(int)

        # --- Quality Control (Verbose) ---
        if verbose:
            unknowns = sex_str[~sex_str.isin(["FEMALE", "F", "MALE", "M"])].unique()
            if len(unknowns) > 0:
                print(f"Warning: Detected unrecognized sex labels treated as Male (0): {unknowns}")
        
        # Identify common CpGs
        common_cpgs = pd.Index(self.required_cpgs).intersection(beta_df.index)
        
        if verbose:
            print(f"[{self.name}] Represented CpGs: {len(common_cpgs)}/{len(self.required_cpgs)}")
            
        # Prepare input matrix (Sample x CpG)
        # Missing CpGs are automatically filled with 0 (or preprocessed values)
        input_data = beta_df.loc[common_cpgs].T
        
        # Add metadata columns
        input_data['Intercept'] = 1
        input_data['Age'] = ages
        input_data['Female'] = is_female

        # --- 2. Phase 1: Calculate DNAm Surrogates ---
        # Filter surrogate models (excluding the final COX model)
        surrogate_coefs = self.coefs[self.coefs['Y.pred'] != 'COX']
        
        # Storage for predicted protein levels
        surrogates_df = pd.DataFrame(index=input_data.index)
        
        # Calculate each surrogate (e.g., DNAmADM, DNAmPACKYRS)
        for target_name, group_df in surrogate_coefs.groupby('Y.pred'):
            # Get required variables (CpG + Age + Sex + Intercept)
            needed_vars = group_df['var'].values
            weights = group_df.set_index('var')['beta']
            
            # Align input data (fill missing variables with 0)
            X = input_data.reindex(columns=needed_vars, fill_value=0)
            
            # Linear combination
            surrogates_df[target_name] = X.dot(weights)

        # --- 3. Phase 2: Calculate Raw COX Score ---
        # Merge calculated surrogates into input data
        combined_data = pd.concat([input_data[['Age', 'Female']], surrogates_df], axis=1)
        
        # Get COX model coefficients
        cox_coefs = self.coefs[self.coefs['Y.pred'] == 'COX']
        cox_vars = cox_coefs['var'].values
        cox_weights = cox_coefs.set_index('var')['beta']
        
        # Align features
        X_final = combined_data.reindex(columns=cox_vars, fill_value=0)
        
        # Calculate COX Score
        raw_cox = X_final.dot(cox_weights)
        
        # --- 4. Phase 3: Calibration ---
        # Formula: (COX - Mean_COX) / SD_COX * SD_Age + Mean_Age
        
        # Get calibration parameters
        param_cox = self.gold[self.gold['var'] == 'COX'].iloc[0]
        param_age = self.gold[self.gold['var'] == 'Age'].iloc[0]
        
        m_cox = param_cox['mean']
        sd_cox = param_cox['sd']
        m_age = param_age['mean']
        sd_age = param_age['sd']
        
        # Transform
        final_grimage = ((raw_cox - m_cox) / sd_cox) * sd_age + m_age
        
        # --- 5. Consolidate Output ---
        results = surrogates_df.copy()
        results['Female'] = combined_data['Female']
        results['Age'] = combined_data['Age']
        results['raw_cox'] = raw_cox
        results['DNAmGrimAge1'] = final_grimage
        
        return results


class GrimAge2(GrimAge1): # Inherits from GrimAge1
    """
    Implements DNAm GrimAge2 (Lu et al. 2022).
    Inherits calculation logic from GrimAge1 but uses updated coefficients.

    References:
        Lu AT, et al. DNA methylation GrimAge version 2. 
        Aging (Albany NY) (2022). https://doi.org/10.18632/aging.204434
    """

    METADATA = {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "Mortality",
        "source": "https://doi.org/10.18632/aging.204434"
    }
    
    def __init__(self):
        # 1. Model identification
        self.name = "GrimAge2"
        
        # 2. Load model components
        self.coefs = load_clock_coefs("GrimAge2_Coefs")
        self.gold = load_clock_coefs("GrimAge2_Gold")
        
        # 3. Identify feature requirements
        non_cpg_vars = ["Intercept", "Age", "Female"]
        all_vars = self.coefs['var'].unique()
        self.required_cpgs = [v for v in all_vars if v not in non_cpg_vars]
        
        # 4. Parent class initialization
        BaseLinearClock.__init__(self, pd.DataFrame(columns=['probe', 'coef']), self.name, self.METADATA)
        
        # 5. Compatibility placeholders
        self.weights = pd.Series(1, index=self.required_cpgs)
        self.intercept = 0.0
        
    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves the composite coefficients for GrimAge2.
        Returns a DataFrame including the 'Y.pred' column to distinguish 
        between different surrogate models (e.g., DNAmADM, DNAmPACKYRS).
        """
        df = self.coefs.copy()
        
        # var -> probe
        # beta -> coef
        if 'var' in df.columns and 'beta' in df.columns:
            df = df.rename(columns={'var': 'probe', 'beta': 'coef'})
            
        if 'Y.pred' in df.columns:
            df = df[['Y.pred', 'probe', 'coef']]
            
        return df

    def predict(
        self, 
        beta_df: pd.DataFrame, 
        ages: Union[list, pd.Series], 
        sex: Union[list, pd.Series], 
        verbose: bool = True
    ) -> pd.DataFrame:
        """
        Calculates GrimAge and its components.

        Args:
            beta_df (pd.DataFrame): Methylation data (Rows: CpGs, Columns: Samples).
            ages (list/Series): **Chronological Age**. Required for the model.
            sex (list/Series): **Biological Sex**. 
                Format: 'Female'/'Male' or 'F'/'M' (Case insensitive).
            verbose (bool): Print warnings about unrecognized sex labels.

        Returns:
            pd.DataFrame: A table containing:
            - **DNAmGrimAge2**: The final biological age estimate.
            - **DNAmADM, DNAmGDF15, DNAmPACKYRS...**: Individual surrogate protein scores.
            - **raw_cox**: The raw mortality score before scaling.
        """
        # 1. Leverage parent class logic directly
        # The parent method uses self.coefs, which is now pointing to GrimAge2 data
        results = super().predict(beta_df, ages, sex, verbose)
        
        # 2. Rename result column
        if 'DNAmGrimAge1' in results.columns:
            results = results.rename(columns={'DNAmGrimAge1': 'DNAmGrimAge2'})
            
        # 3. Apply GrimAge2-specific variable renaming (matching R implementation)
        rename_map = {
            'DNAmadm': 'DNAmADM', 
            'DNAmCystatin_C': 'DNAmCystatinC', 
            'DNAmGDF_15': 'DNAmGDF15', 
            'DNAmleptin': 'DNAmLeptin',
            'DNAmpai_1': 'DNAmPAI1', 
            'DNAmTIMP_1': 'DNAmTIMP1', 
            'DNAmlog.CRP': 'DNAmlogCRP', 
            'DNAmlog.A1C': 'DNAmlogA1C'
        }
        results = results.rename(columns=rename_map)
        
        return results


class IC_Clock(BaseLinearClock):
    """
    Implements the Intrinsic Capacity (IC) Clock (Fuentealba et al. 2025).

    References:
        Fuentealba et al. A blood-based epigenetic clock for intrinsic capacity predicts mortality and is associated with clinical, immunological and lifestyle factors.
        Nat Aging (2025). https://doi.org/10.1038/s43587-025-00883-5
    """
    METADATA = {
        "year": 2025,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(EPIC)",
        "prediction": "Intrinsic Capacity",
        "source": "https://doi.org/10.1038/s43587-025-00883-5"
    }
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs("IC_Clock")
        super().__init__(coef_df, name="IC_Clock",metadata=self.METADATA)


class DNAmFitAge(BaseLinearClock):
    """
    Implements DNAmFitAge (McGreevy et al. 2023).
    A composite biomarker combining DNA methylation, Age, Sex, and DNAmGrimAge.

    References:
        McGreevy KM, et al. DNAmFitAge: biological age indicator incorporating physical fitness. 
        Aging (Albany NY) (2023). https://doi.org/10.18632/aging.204538
    """
    METADATA = {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "prediction": "chronological age(years)",
        "source": "https://doi.org/10.18632/aging.204538"
    }
    
    def __init__(self):
        """Initialize DNAmFitAge, loading sub-models and sex-specific imputation references."""
        self.name = "DNAmFitAge"
        
        # 1. Load sub-model coefficients
        self.sub_models = load_clock_coefs("DNAmFitAge_SubModels")
        
        # 2. Load imputation reference values (sex-specific medians)
        self.fem_medians_df = load_clock_coefs("DNAmFitAge_Female_Medians")
        self.male_medians_df = load_clock_coefs("DNAmFitAge_Male_Medians")
        
        # Optimize for lookup: Convert to Series
        self.fem_medians = self.fem_medians_df.set_index('probe')['median_val']
        self.male_medians = self.male_medians_df.set_index('probe')['median_val']
        
        # 3. Define required probes
        all_probes = self.sub_models['probe'].unique()
        
        # Exclude non-CpG features (Intercept, Age) to prevent incorrect imputation logic
        self.required_cpgs = [
            p for p in all_probes 
            if p not in ["Intercept", "Age", "age"]
        ]
        # 4. Parent class initialization
        super().__init__(pd.DataFrame(columns=['probe', 'coef']), self.name, self.METADATA)
        
        # BaseLinearClock compatibility
        self.weights = pd.Series(1, index=self.required_cpgs)
        self.intercept = 0.0

    def get_coefs(self) -> pd.DataFrame:
        """
        Retrieves coefficients for all sub-models within DNAmFitAge.
        Returns a DataFrame with columns: ['model_name', 'probe', 'coef'].
        """
        # copy to avoid modifying the internal attribute
        df = self.sub_models.copy()
        
        # Ensure standard column names
        # Your load_clock_coefs likely returns what's in the CSV.
        # We want to ensure 'model_name' exists to distinguish between 
        # Gait_noAge_Females, VO2maxModel, etc.
        required_cols = {'model_name', 'probe', 'coef'}
        if not required_cols.issubset(df.columns):
            # Fallback if column names differ (adjust based on your actual CSV header)
            # For example, if your CSV has 'model' instead of 'model_name'
            rename_map = {}
            if 'model' in df.columns: rename_map['model'] = 'model_name'
            if 'var' in df.columns: rename_map['var'] = 'probe'
            if 'beta' in df.columns: rename_map['beta'] = 'coef'
            df = df.rename(columns=rename_map)
            
        return df[['model_name', 'probe', 'coef']]

    def _impute_by_sex(self, beta_df, sex_series):
        """
        Implements sex-specific median imputation for missing CpGs.
        """
        # Identify missing probes
        missing_probes = pd.Index(self.required_cpgs).difference(beta_df.index)
        
        if not missing_probes.empty:
            # Get sample indices by sex
            fem_samples = sex_series[sex_series == 1].index
            male_samples = sex_series[sex_series == 0].index
            
            # Prepare imputation matrix
            imputed_data = pd.DataFrame(index=missing_probes, columns=beta_df.columns)
            
            # Impute Females (using female medians)
            if len(fem_samples) > 0:
                fem_vals = self.fem_medians.reindex(missing_probes).fillna(0.5).values
                imputed_data[fem_samples] = np.tile(fem_vals[:, None], (1, len(fem_samples)))
                
            # Impute Males (using male medians)
            if len(male_samples) > 0:
                male_vals = self.male_medians.reindex(missing_probes).fillna(0.5).values
                imputed_data[male_samples] = np.tile(male_vals[:, None], (1, len(male_samples)))
            
            # Concatenate imputed data
            beta_df = pd.concat([beta_df, imputed_data])
            
        return beta_df.loc[self.required_cpgs]
    def predict(
        self, 
        beta_df: pd.DataFrame, 
        ages: Union[list, pd.Series], 
        sex: Union[list, pd.Series], 
        grimage: Union[list, pd.Series], 
        verbose: bool = True
    ) -> pd.DataFrame:
        """
        Calculates DNAmFitAge and its fitness components.

        Args:
            beta_df (pd.DataFrame): Methylation data (Rows: CpGs, Columns: Samples).
            ages (list/Series): Chronological age.
            sex (list/Series): Biological Sex ('Female'/'Male').
            grimage (list/Series): **Pre-calculated DNAmGrimAge**. 
                (Use ``omniage.models.GrimAge1`` to generate this).
            verbose (bool): Print warnings.

        Returns:
            pd.DataFrame: A DataFrame containing:
            - **DNAmFitAge**: The composite biological age.
            - **FitAgeAccel**: Age acceleration (Residuals vs Age).
            - **DNAmVO2max**: Predicted cardiorespiratory fitness.
            - **DNAmGrip_wAge**: Predicted grip strength (adjusted for age).
            - **DNAmGait_wAge**: Predicted gait speed (adjusted for age).
            - **DNAmFEV1_wAge**: Predicted lung function.
        """
        
        # --- 1. Data Preparation ---
        if isinstance(ages, list): ages = pd.Series(ages, index=beta_df.columns)
        if isinstance(sex, list): sex = pd.Series(sex, index=beta_df.columns)
        if isinstance(grimage, list): grimage = pd.Series(grimage, index=beta_df.columns)
        
        # --- Sex Handling & Validation ---
        if pd.api.types.is_numeric_dtype(sex):
            raise ValueError("Ambiguous sex input: Please use strings 'Female'/'F' or 'Male'/'M'.")

        sex_str = sex.astype(str).str.strip().str.upper()
        female_labels = ["FEMALE", "F"]
        is_female = sex_str.isin(female_labels).astype(int)

        if verbose:
            unknowns = sex_str[~sex_str.isin(["FEMALE", "F", "MALE", "M"])].unique()
            if len(unknowns) > 0:
                print(f"Warning: Unrecognized sex labels treated as Male (0): {unknowns}")

        # Perform sex-specific imputation
        beta_df = self._impute_by_sex(beta_df, is_female)
        
        # Transpose to (Sample x CpG)
        input_data = beta_df.T
        input_data['Intercept'] = 1
        
        # Add Age to input matrix (Crucial for *wAge models)
        # Ensure column name matches coefficients (usually "Age" or "age")
        input_data['Age'] = ages.values 
        
        # --- 2. Calculate 6 Fitness Biomarkers ---
        fitness_results = pd.DataFrame(index=input_data.index)
        
        # Map output metrics to specific sub-models
        metrics_map = {
            "DNAmGait_noAge": ["Gait_noAge_Females", "Gait_noAge_Males"],
            "DNAmGrip_noAge": ["Grip_noAge_Females", "Grip_noAge_Males"],
            "DNAmGait_wAge":  ["Gait_wAge_Females", "Gait_wAge_Males"],
            "DNAmGrip_wAge":  ["Grip_wAge_Females", "Grip_wAge_Males"],
            "DNAmFEV1_wAge":  ["FEV1_wAge_Females", "FEV1_wAge_Males"],
            "DNAmVO2max":     ["VO2maxModel", "VO2maxModel"] 
        }
        
        for metric, (fem_model, male_model) in metrics_map.items():
            fitness_results[metric] = np.nan
            
            # Predict for Females
            fem_idx = is_female[is_female == 1].index
            if len(fem_idx) > 0:
                coefs = self.sub_models[self.sub_models['model_name'] == fem_model]
                if not coefs.empty:
                    w = coefs.set_index('probe')['coef']
                    # Reindex ensures 'Age' is included if 'w' contains it
                    X = input_data.loc[fem_idx].reindex(columns=w.index, fill_value=0)
                    fitness_results.loc[fem_idx, metric] = X.dot(w)
            
            # Predict for Males
            male_idx = is_female[is_female == 0].index
            if len(male_idx) > 0:
                coefs = self.sub_models[self.sub_models['model_name'] == male_model]
                if not coefs.empty:
                    w = coefs.set_index('probe')['coef']
                    X = input_data.loc[male_idx].reindex(columns=w.index, fill_value=0)
                    fitness_results.loc[male_idx, metric] = X.dot(w)

        # --- 3. Calculate Final DNAmFitAge ---
        calc_df = fitness_results.copy()
        calc_df['DNAmGrimAge'] = grimage
        calc_df['Age'] = ages
        calc_df['Female'] = is_female
        
        fit_age = pd.Series(np.nan, index=calc_df.index)
        
        # Formula for Females (Hardcoded standardization parameters)
        fem_mask = (calc_df['Female'] == 1)
        if fem_mask.any():
            f_dat = calc_df[fem_mask]
            val = (
                0.1044232 * ((f_dat["DNAmVO2max"] - 46.825091) / (-0.13620215)) +
                0.1742083 * ((f_dat["DNAmGrip_noAge"] - 39.857718) / (-0.22074456)) +
                0.2278776 * ((f_dat["DNAmGait_noAge"] - 2.508547) / (-0.01245682)) +
                0.4934908 * ((f_dat["DNAmGrimAge"] - 7.978487) / (0.80928530))
            )
            fit_age[fem_mask] = val
            
        # Formula for Males
        male_mask = (calc_df['Female'] == 0)
        if male_mask.any():
            m_dat = calc_df[male_mask]
            val = (
                0.1390346 * ((m_dat["DNAmVO2max"] - 49.836389) / (-0.141862925)) +
                0.1787371 * ((m_dat["DNAmGrip_noAge"] - 57.514016) / (-0.253179827)) +
                0.1593873 * ((m_dat["DNAmGait_noAge"] - 2.349080) / (-0.009380061)) +
                0.5228411 * ((m_dat["DNAmGrimAge"] - 9.549733) / (0.835120557))
            )
            fit_age[male_mask] = val

        # --- 4. Calculate Age Acceleration (Residuals) ---
        valid_idx = fit_age.dropna().index
        if len(valid_idx) > 1:
            y = fit_age[valid_idx]
            x = ages[valid_idx]
            coeffs = np.polyfit(x, y, 1)
            predicted_y = np.polyval(coeffs, x)
            residuals = y - predicted_y
            
            accel = pd.Series(np.nan, index=fit_age.index)
            accel[valid_idx] = residuals
        else:
            accel = fit_age - ages 

        # --- 5. Final Output ---
        fitness_results['DNAmFitAge'] = fit_age
        fitness_results['FitAgeAccel'] = accel
        
        return fitness_results




class CausalAge(BaseLinearClock):
    """
    CausalAge: Enriched for CpGs with a causal effect on mortality.
    
    Reference:
    Ying et al. Causality-enriched epigenetic age uncouples damage and adaptation. 
    Nature Aging (2024). https://doi.org/10.1038/s43587-023-00557-0
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(EPIC)",
        "Prediction": "Chronological age(years)",
        "source": "https://doi.org/10.1038/s43587-023-00557-0" 
    }
    CLOCK_NAME = "CausalAge" 
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs(self.CLOCK_NAME)
        super().__init__(coef_df, name=self.CLOCK_NAME, metadata=self.METADATA)



class DamAge(BaseLinearClock): 
    """
    DamAge: Captures the accumulation of molecular damage (specifically associated with mortality).
    
    Reference:
    Ying et al. Causality-enriched epigenetic age uncouples damage and adaptation. 
    Nature Aging (2024). https://doi.org/10.1038/s43587-023-00557-0
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(EPIC)",
        "Prediction": "Biological Damage (Years)", 
        "source": "https://doi.org/10.1038/s43587-023-00557-0"
    }
    CLOCK_NAME = "DamAge"

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs(self.CLOCK_NAME)
        super().__init__(coef_df, name=self.CLOCK_NAME, metadata=self.METADATA)


class AdaptAge(BaseLinearClock): 
    """
    AdaptAge: Captures protective adaptations to age-related damage.

    Reference:
    Ying et al. Causality-enriched epigenetic age uncouples damage and adaptation. 
    Nature Aging (2024). https://doi.org/10.1038/s43587-023-00557-0
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(EPIC)",
        "Prediction": "Adaptive Response (years)", 
        "source": "https://doi.org/10.1038/s43587-023-00557-0"
    }
    CLOCK_NAME = "AdaptAge"

    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs(self.CLOCK_NAME)
        super().__init__(coef_df, name=self.CLOCK_NAME, metadata=self.METADATA)


class StocH(BaseLinearClock):
    """
    Stochastic analogue of the Horvath clock (Tong et al., 2024).

    Reference:
        Tong et al.  Quantifying the stochastic component of epigenetic aging.
        Nature Aging (2024). https://doi.org/10.1038/s43587-024-00600-8
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Multi-tissue",
        "omic type": "DNAm(450k/EPIC)",
        "Prediction": "Stochastic Aging Score(Horvath2013)",
        "source": "https://doi.org/10.1038/s43587-024-00600-8"
    }
    CLOCK_NAME = "StocH"
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs(self.CLOCK_NAME)
        super().__init__(coef_df, name=self.CLOCK_NAME, metadata=self.METADATA)



class StocZ(BaseLinearClock):
    """
    Stochastic analogue of the Zhang clock (Tong et al., 2024).

    Reference:
        Tong et al.  Quantifying the stochastic component of epigenetic aging.
        Nature Aging (2024). https://doi.org/10.1038/s43587-024-00600-8
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "Prediction": "Stochastic Aging Score(Zhang)",
        "source": "https://doi.org/10.1038/s43587-024-00600-8"
    }
    CLOCK_NAME = "StocZ"
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs(self.CLOCK_NAME)
        super().__init__(coef_df, name=self.CLOCK_NAME, metadata=self.METADATA)

class StocP(BaseLinearClock):
    """
    Stochastic analogue of PhenoAge (Tong et al., 2024).

    Reference:
        Tong et al.  Quantifying the stochastic component of epigenetic aging.
        Nature Aging (2024). https://doi.org/10.1038/s43587-024-00600-8
    """
    METADATA = {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k/EPIC)",
        "Prediction": "Stochastic Aging Score(PhenoAge)",
        "source": "https://doi.org/10.1038/s43587-024-00600-8"
    }
    CLOCK_NAME = "StocP"
    def __init__(self, coef_df: Optional[pd.DataFrame] = None):
        if coef_df is None:
            coef_df = load_clock_coefs(self.CLOCK_NAME)
        super().__init__(coef_df, name=self.CLOCK_NAME, metadata=self.METADATA)




class BaseMcCartneyClock(BaseLinearClock):
    """
    Base class for McCartney et al. (2018) complex trait predictors.
    
    This class manages the loading of coefficients for various phenotypic traits 
    predicted from DNA methylation data.

    References:
        McCartney DL, et al. Epigenetic prediction of complex traits and death. 
        Genome Biol (2018). https://doi.org/10.1186/s13059-018-1514-1
    """

    METADATA = {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "omic type": "DNAm(450k)",
        "Prediction": "Phenotypic Trait Score",
        "source": "https://doi.org/10.1186/s13059-018-1514-1"
    }

    def __init__(self, trait_name: str, coef_df: Optional[pd.DataFrame] = None):
        """
        Initialize a McCartney Trait Predictor.
        
        Args:
            trait_name (str): The specific trait identifier (e.g., 'BMI', 'Smoking').
            data_dir (str, optional): Custom data directory.
        """
        # Standardize clock name
        clock_name = f"McCartney_{trait_name}"
        coef_path = f"./McCartney/McCartney_{trait_name}"
        if coef_df is None:
            coef_df = load_clock_coefs(coef_path)
        # Initialize parent class
        super().__init__(coef_df, name=clock_name, metadata=self.METADATA)


# --- Specific Trait Implementations ---

class McCartneyBMI(BaseMcCartneyClock):
    """Epigenetic predictor for Body Mass Index (BMI)."""
    def __init__(self):
        super().__init__("BMI")

class McCartneySmoking(BaseMcCartneyClock):
    """Epigenetic predictor for Smoking (pack-years)."""
    def __init__(self):
        super().__init__("Smoking")

class McCartneyAlcohol(BaseMcCartneyClock):
    """Epigenetic predictor for Alcohol consumption (units/week)."""
    def __init__(self):
        super().__init__("Alcohol")

class McCartneyEducation(BaseMcCartneyClock):
    """Epigenetic predictor for Education (years)."""
    def __init__(self):
        super().__init__("Education")

class McCartneyTotalCholesterol(BaseMcCartneyClock):
    """Epigenetic predictor for Total Cholesterol."""
    def __init__(self):
        super().__init__("Total_cholesterol")

class McCartneyHDL(BaseMcCartneyClock):
    """Epigenetic predictor for HDL Cholesterol."""
    def __init__(self):
        super().__init__("HDL_cholesterol")

class McCartneyLDL(BaseMcCartneyClock):
    """Epigenetic predictor for LDL Cholesterol."""
    def __init__(self):
        super().__init__("LDL_cholesterol")

class McCartneyTotalHDLRatio(BaseMcCartneyClock):
    """Epigenetic predictor for Total/HDL Cholesterol Ratio."""
    def __init__(self):
        super().__init__("Total_HDL_ratio")

class McCartneyWHR(BaseMcCartneyClock):
    """Epigenetic predictor for Waist-to-Hip Ratio (WHR)."""
    def __init__(self):
        super().__init__("WHR")

class McCartneyBodyFat(BaseMcCartneyClock):
    """Epigenetic predictor for Body Fat Percentage."""
    def __init__(self):
        super().__init__("Body_fat_Perc")





