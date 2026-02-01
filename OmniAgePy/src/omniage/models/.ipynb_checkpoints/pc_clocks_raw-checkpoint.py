import pandas as pd
import numpy as np
import os
from pathlib import Path
from typing import Optional

class BasePCClock:
    """
    Handles data loading, Parquet reading, preprocessing, and basic PCA projection.
    [Maximum Performance Version: R-Logic Replica]
    """
    def __init__(self, clock_name: str):
        self.name = clock_name
        current_file_dir = Path(__file__).parent.resolve()
        self.base_path = current_file_dir.parent / "data" / "PCClocks"
        self.clock_dir = self.base_path / clock_name

        if not self.clock_dir.exists():
            raise FileNotFoundError(f"Data for {clock_name} not found at {self.clock_dir}")

        self._load_assets()

    def _load_assets(self):
        # Force float32 during load to save memory and speed up computation
        
        # Load Global Imputation Means
        impute_path = self.base_path / "PC_Impute_Means.parquet"
        self.impute_ref = pd.read_parquet(impute_path).set_index("probe")["mean"].astype(np.float32)
        
        # Load Center and Rotation
        self.center = pd.read_parquet(self.clock_dir / "center.parquet").set_index("probe")["mean"].astype(np.float32)
        self.rotation = pd.read_parquet(self.clock_dir / "rotation.parquet").set_index("probe").astype(np.float32)

    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        Alignment, Transpose, and Layered Imputation.
        Strategy: Intersection -> Mean Impute Existing -> Concat Missing.
        """
        # 1. Identify Required vs. Missing CpGs
        required_cpgs = self.rotation.index
        
        # Fast Intersection: Only load what we have
        common_cpgs = beta_df.index.intersection(required_cpgs)
        missing_cpgs = required_cpgs.difference(common_cpgs)
        
        # 2. Extract Existing Data & Transpose (Force float32)
        # Only manipulating actual data in memory
        X_existing = beta_df.loc[common_cpgs].T.astype(np.float32)
        
        # ============================================================
        # 3. Layer 1 Imputation: Handle Sporadic NAs in Existing Data
        #    (Corresponds to R's 'mean_impute')
        # ============================================================
        if X_existing.isna().values.any():
            # Find columns with NAs (Boolean indexing is fast)
            nan_cols = X_existing.columns[X_existing.isna().any()]
            
            # Calculate Cohort Mean for these columns
            col_means = X_existing[nan_cols].mean(axis=0, skipna=True)
            
            # Fallback: If Cohort Mean is NaN (all empty), use Reference
            col_means = col_means.fillna(self.impute_ref)
            
            # Fill (Only touches bad columns)
            X_existing[nan_cols] = X_existing[nan_cols].fillna(col_means)

        # ============================================================
        # 4. Layer 2 Imputation: Construct & Concat Missing Columns
        #    (Corresponds to R's 'needed_matrix' + 'cbind')
        # ============================================================
        if len(missing_cpgs) > 0:
            # Get reference values for missing columns
            fill_vals = self.impute_ref.loc[missing_cpgs].values.astype(np.float32)
            
            # Use np.tile to instantly generate the matrix (O(1) vs O(N))
            n_samples = X_existing.shape[0]
            missing_data = np.tile(fill_vals, (n_samples, 1))
            
            # Create DataFrame
            X_missing = pd.DataFrame(
                missing_data, 
                index=X_existing.index, 
                columns=missing_cpgs
            )
            
            # Concat existing and missing parts
            X = pd.concat([X_existing, X_missing], axis=1)
        else:
            X = X_existing

        # 5. Final Reorder (Essential for matrix multiplication alignment)
        X = X[required_cpgs]
        
        # 6. Final Fallback (Safety net)
        X = X.fillna(0.0)
            
        return X

    def get_pcs(self, X: pd.DataFrame) -> pd.DataFrame:
        """Calculate Principal Components [Optimized]"""
        # Ensure center vector is aligned
        center_vec = self.center.loc[X.columns].values
        
        # Use Numpy for broadcast subtraction and matrix multiplication (All float32)
        X_values = X.values 
        X_centered = X_values - center_vec
        
        rotation_values = self.rotation.values
        
        pcs = X_centered @ rotation_values
        
        return pd.DataFrame(pcs, index=X.index, columns=self.rotation.columns)

# --- 3. 标准 PC 时钟实现 ---
class StandardPCClock(BasePCClock):
    """
    适用于 Horvath, Hannum, PhenoAge, DNAmTL 等
    """
    def __init__(self, clock_name: str, do_anti_trafo: bool):
        super().__init__(clock_name)
        self.do_anti_trafo = do_anti_trafo
        
        # 加载模型系数 (float32)
        model_df = pd.read_parquet(self.clock_dir / "model.parquet")
        self.coefs = pd.Series(model_df['coef'].values.astype(np.float32), index=model_df['pc'])
        
        # 加载截距
        with open(self.clock_dir / "intercept.txt", "r") as f:
            self.intercept = float(f.read().strip())

    def anti_trafo(self, x: np.ndarray, adult_age: float = 20) -> np.ndarray:
        return np.where(x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age)

    def predict(self, beta_df: pd.DataFrame, verbose: bool = False, **kwargs) -> pd.Series:
        X = self.preprocess(beta_df)
        pcs_df = self.get_pcs(X)
        
        common_pcs = self.coefs.index.intersection(pcs_df.columns)
        
        # 矩阵运算
        raw_score = pcs_df[common_pcs] @ self.coefs[common_pcs] + self.intercept
        
        if self.do_anti_trafo:
            return pd.Series(self.anti_trafo(raw_score.values), index=X.index)
        return raw_score

# --- 4. GrimAge 实现 (关键适配 & 性能优化) ---
class PCGrimAge1Impl(BasePCClock):
    """
    PCGrimAge1 的具体实现逻辑
    """
    def __init__(self):
        super().__init__("PCGrimAge1")
        # 加载 Surrogate 和 Final 模型
        surr_df = pd.read_parquet(self.clock_dir / "surrogate_weights.parquet")
        self.surr_intercepts = surr_df[surr_df['term'] == 'Intercept'].set_index('target')['coef'].astype(np.float32)
        self.surr_weights = surr_df[surr_df['term'] != 'Intercept'].copy()
        # 注意：这里 coefs 转 float32
        self.surr_weights['coef'] = self.surr_weights['coef'].astype(np.float32)
        
        final_df = pd.read_parquet(self.clock_dir / "final_model.parquet")
        self.final_intercept = final_df[final_df['term'] == 'Intercept']['coef'].iloc[0]
        self.final_weights = final_df[final_df['term'] != 'Intercept'].set_index('term')['coef'].astype(np.float32)

    def predict(self, beta_df: pd.DataFrame, ages: pd.Series, sex: pd.Series, verbose: bool = False, **kwargs) -> pd.DataFrame:
        
        if ages is None or sex is None:
            raise ValueError("[PCGrimAge1] Missing required 'ages' or 'sex' arguments.")

        # --- 1. 预处理 (Fast) ---
        X = self.preprocess(beta_df)
        pcs_df = self.get_pcs(X)
        samples = X.index
        
        # --- 2. 快速对齐 & 类型转换 ---
        # 确保 Age 是 float
        ages = ages.reindex(samples).astype(float)
        sex = sex.reindex(samples)

        # --- 3. 性别标准化 (Vectorized) ---
        sex_cleaned = sex.astype(str).str.strip().str.lower()
        sex_map_to_int = {
            'f': 1, 'female': 1, 'woman': 1, 'w': 1, 
            'm': 0, 'male': 0, 'man': 0,
        }
        is_female = sex_cleaned.map(sex_map_to_int)
        
        if is_female.isna().any():
            invalid_entries = sex[is_female.isna()].unique()
            raise ValueError(
                f"[PCGrimAge1] Input 'sex' contains unrecognized values: {list(invalid_entries)}. "
                f"Supported variants (case-insensitive): {list(sex_map_to_int.keys())}"
            )

        # --- 4. 优化 Surrogate 计算 (In-place & No Reindex Loop) ---
        # 直接修改 pcs_df，避免内存拷贝 (因为 pcs_df 是这里生成的局部变量)
        pcs_df['Age'] = ages
        pcs_df['Female'] = is_female
        
        # 预分配结果 DataFrame
        targets = self.surr_weights['target'].unique()
        surrogate_scores = pd.DataFrame(index=samples, columns=targets)
        
        for target in targets:
            w_df = self.surr_weights[self.surr_weights['target'] == target]
            w_series = pd.Series(w_df['coef'].values, index=w_df['term'])
            intercept = self.surr_intercepts.get(target, 0.0)
            
            # [优化点]: 不要在循环里 reindex(fill_value=0)
            # 而是取交集。如果 pcs_df 缺某些列，在点积中自然无法贡献（等同于0）
            # 只要 input 包含了所有必要的 PCs 和 Age/Female 即可
            common_features = w_series.index.intersection(pcs_df.columns)
            
            # 向量化点积
            val = pcs_df[common_features] @ w_series[common_features] + intercept
            surrogate_scores[target] = val
            
        # --- 5. 准备最终特征 & 列名映射 ---
        # 浅拷贝即可
        final_features = surrogate_scores 
        final_features['Age'] = ages
        final_features['Female'] = is_female
        
        rename_map = {
            'PCPACKYRS': 'DNAmPACKYRS',
            'PCADM': 'DNAmADM',
            'PCB2M': 'DNAmB2M',
            'PCCystatinC': 'DNAmCystatinC',
            'PCGDF15': 'DNAmGDF15',
            'PCLeptin': 'DNAmLeptin',
            'PCPAI1': 'DNAmPAI1',
            'PCTIMP1': 'DNAmTIMP1'
        }
        # 这里会产生一次拷贝，因为改了列名
        final_features_renamed = final_features.rename(columns=rename_map)

        # --- 6. 最终预测 ---
        # 这里的 reindex 很快，因为 final_features 只有十几列
        aligned_features = final_features_renamed.reindex(columns=self.final_weights.index, fill_value=0)
        
        pred_age = aligned_features @ self.final_weights + self.final_intercept
        
        # 存入结果
        surrogate_scores[self.name] = pred_age
        
        return surrogate_scores

# --- 5. 导出类 (Wrapper Classes) ---

class PCHorvath2013(StandardPCClock):
    def __init__(self):
        super().__init__("PCHorvath2013", do_anti_trafo=True)

class PCHorvath2018(StandardPCClock):
    def __init__(self):
        super().__init__("PCHorvath2018", do_anti_trafo=True)

class PCHannum(StandardPCClock):
    def __init__(self):
        super().__init__("PCHannum", do_anti_trafo=False)

class PCPhenoAge(StandardPCClock):
    def __init__(self):
        super().__init__("PCPhenoAge", do_anti_trafo=False)

class PCDNAmTL(StandardPCClock):
    def __init__(self):
        super().__init__("PCDNAmTL", do_anti_trafo=False)

class PCGrimAge1(PCGrimAge1Impl):
    pass