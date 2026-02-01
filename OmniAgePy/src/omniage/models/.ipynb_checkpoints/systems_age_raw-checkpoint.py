import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional



class SystemsAge:
    """
    Implements the Systems Age epigenetic clock (Sehgal et al. 2025).
    Calculates 11 system-specific scores and one composite Systems Age.
    """
    def __init__(self):
        self.name = "SystemsAge"
        
        # 1. Locate Data Directory
        # Assumes data is in ../data/SystemsAge relative to this file
        current_file_dir = Path(__file__).parent.resolve()
        self.base_dir = current_file_dir.parent / "data" / "SystemsAge"
        
        if not self.base_dir.exists():
            raise FileNotFoundError(f"SystemsAge data not found at {self.base_dir}")
            
        self._load_assets()

    def _load_assets(self):
        """Load all Parquet files required for calculation."""
        
        # A. Imputation
        self.impute_ref = pd.read_parquet(self.base_dir / "impute_means.parquet").set_index("probe")["mean"]
        
        # B. DNAm PCA (CpG -> DNAmPC)
        self.dnam_pca_center = pd.read_parquet(self.base_dir / "dnam_pca_center.parquet").set_index("probe")["mean"]
        self.dnam_pca_rot = pd.read_parquet(self.base_dir / "dnam_pca_rotation.parquet").set_index("probe")
        
        # C. System Vector Coefficients (DNAmPC -> SystemPC)
        # Columns in parquet: [Blood, Brain, ..., pc]
        # We need index=pc, columns=Systems
        df_sys_vec = pd.read_parquet(self.base_dir / "system_vector_coefs.parquet")
        # Ensure 'pc' column exists and use it as index
        if 'pc' in df_sys_vec.columns:
            self.sys_vec_coefs = df_sys_vec.set_index("pc")
        else:
            # Fallback if index wasn't saved as column
            self.sys_vec_coefs = df_sys_vec
            
        # D. System Score Coefs (Scaling)
        self.sys_score_coefs = pd.read_parquet(self.base_dir / "system_score_coefs.parquet").set_index("term")["coef"]
        
        # E. Chronological Age Prediction Model
        age_coefs_df = pd.read_parquet(self.base_dir / "age_pred_coefs.parquet")
        self.age_intercept = age_coefs_df.loc[age_coefs_df['term'] == 'Intercept', 'coef'].iloc[0]
        self.age_linear_weights = age_coefs_df[age_coefs_df['term'] != 'Intercept'].set_index("term")["coef"]
        
        # Quadratic Model Params (Intercept, Linear, Quadratic)
        # Assumed structure: param=['Intercept', 'Linear', 'Quadratic'], value=[...]
        age_params_df = pd.read_parquet(self.base_dir / "age_model_params.parquet").set_index("param")
        # Map R names to standard keys if necessary, strictly following R code logic:
        # RData$Age_prediction_model[1] -> Intercept (Constant)
        # RData$Age_prediction_model[2] -> Linear Term
        # RData$Age_prediction_model[3] -> Quadratic Term
        # We assume the parquet export kept the order or names.
        # Safe access by iloc to match R indexing [1], [2], [3]
        self.age_model_const = age_params_df.iloc[0, 0]
        self.age_model_lin = age_params_df.iloc[1, 0]
        self.age_model_quad = age_params_df.iloc[2, 0]

        # F. Systems PCA (System Scores -> Composite Age)
        self.sys_pca_center = pd.read_parquet(self.base_dir / "systems_pca_center.parquet").set_index("term")["mean"]
        self.sys_pca_scale = pd.read_parquet(self.base_dir / "systems_pca_scale.parquet").set_index("term")["scale"]
        self.sys_pca_rot = pd.read_parquet(self.base_dir / "systems_pca_rotation.parquet").set_index("term")
        
        final_coefs_df = pd.read_parquet(self.base_dir / "final_coefs.parquet")
        self.final_coefs = pd.Series(final_coefs_df['coef'].values, index=final_coefs_df['pc'])

        # G. Final Transformation (Z-score to Years)
        # columns: system, V1(mean), V2(sd), V3(ref_mean), V4(ref_sd) -> based on R logic
        self.trans_coefs = pd.read_parquet(self.base_dir / "transformation_coefs.parquet").set_index("system")



    def preprocess(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """
        [R逻辑复刻版] 极致性能优化。
        核心策略：绝对不使用 reindex 进行填充。
        逻辑：
        1. 筛选出现有数据 (Common CpGs)。
        2. 仅对现有数据中有缺失值的列进行均值插补 (Mean Imputation)。
        3. 直接构造完全缺失的列矩阵 (Numpy Tile)。
        4. 最后拼接 (Concat) 并重排顺序。
        """
        # 1. 确定需要的 CpG
        required_cpgs = self.dnam_pca_rot.index # 对于 SystemsAge
        
        # 2. 集合运算：区分“有的”和“完全没的”
        # common_cpgs: 输入数据里有的
        # missing_cpgs: 输入数据里完全没有的
        common_cpgs = beta_df.index.intersection(required_cpgs)
        missing_cpgs = required_cpgs.difference(common_cpgs)
        
        # 3. 提取现有数据 (转置 + float32)
        # 这步只操作内存中实际存在的数据，不开辟多余空间
        X_existing = beta_df.loc[common_cpgs].T.astype(np.float32)
        
        # ============================================================
        # 4. 对应 R: mean_impute (处理现有数据中的 Sporadic NA)
        # ============================================================
        # 只有当确实有 NaN 时才计算，且只计算坏掉的列
        if X_existing.isna().values.any():
            # 找到有 NaN 的列名 (bool索引很快)
            nan_cols = X_existing.columns[X_existing.isna().any()]
            
            # 计算这些列的均值 (Cohort Mean)
            col_means = X_existing[nan_cols].mean(axis=0, skipna=True)
            
            # 兜底：如果 Cohort Mean 也是 NaN (说明该列全空)，用 Reference 填
            # 注意：impute_ref 已经是 float32
            col_means = col_means.fillna(self.impute_ref)
            
            # 填补 (只操作坏掉的列)
            X_existing[nan_cols] = X_existing[nan_cols].fillna(col_means)
            
        # ============================================================
        # 5. 对应 R: needed_matrix (构建完全缺失的列)
        # ============================================================
        if len(missing_cpgs) > 0:
            # 直接从 impute_ref 拿参考值
            fill_vals = self.impute_ref.loc[missing_cpgs].values
            
            # 利用 Numpy 广播直接构造矩阵，速度极快 (比 Pandas fillna 快得多)
            # R 代码对应逻辑: matrix(CpGs[needed_cpgs], nrow=..., byrow=TRUE)
            n_samples = X_existing.shape[0]
            
            # np.tile: 把一行参考值复制 n_samples 行
            missing_data = np.tile(fill_vals, (n_samples, 1))
            
            # 构建 DataFrame
            X_missing = pd.DataFrame(
                missing_data, 
                index=X_existing.index, 
                columns=missing_cpgs,
                dtype=np.float32
            )
            
            # ============================================================
            # 6. 对应 R: cbind (拼接)
            # ============================================================
            # 拼接现有部分和缺失部分
            X = pd.concat([X_existing, X_missing], axis=1)
        else:
            X = X_existing

        # 7. 最终重排 (Reorder)
        # 因为拼接后的列顺序可能是乱的 (现有在前，缺失在后)，必须按 Rotation 的顺序排好
        # 这一步虽然叫 reindex，但因为数据全齐，只是 View 操作，不涉及填充，非常快
        X = X[required_cpgs]
        
        return X

    def predict(self, beta_df: pd.DataFrame, verbose: bool = False, **kwargs) -> pd.DataFrame:
        
        # 1. Preprocess (CpGs)
        X = self.preprocess(beta_df)
        samples = X.index
        
        # 2. Calculate DNAm PCs
        # R: predict(RData$DNAmPCA, DNAm) -> (X - Center) @ Rotation
        X_centered = X - self.dnam_pca_center
        # Note: R's prcomp rotation is usually p x k. X is n x p. Result n x k.
        dnam_pcs = X_centered @ self.dnam_pca_rot
        
        # 3. Calculate DNAm System PCs
        # R Logic: DNAmPCs[, 1:4017] %*% system_vector_coefficients[1:4017, ]
        
        # 修正：不使用 intersection，而是严格切片前 4017 个
        # 假设我们导出数据时保持了顺序 (PC1, PC2...)，Parquet 通常会保留这个顺序。
        
        # 选取 DNAm PCs 的前 4017 列
        # 注意：.iloc[:, :4017] 表示取所有行，前 4017 列
        dnam_pcs_subset = dnam_pcs.iloc[:, :4017]
        
        # 选取系数矩阵的前 4017 行
        # 我们需要确保 self.sys_vec_coefs 的行顺序与 dnam_pcs 的列顺序一致
        # (通常导出时都是基于同一个 Rotation 对象的顺序)
        sys_vec_subset = self.sys_vec_coefs.iloc[:4017]
        
        # 执行矩阵乘法
        # 维度检查: (N_samples x 4017) @ (4017 x 11_Systems)
        dnam_sys_pcs = dnam_pcs_subset @ sys_vec_subset
        
        # 4. Calculate System Scores
        # Iterate through 11 systems
        system_groups = ["Blood", "Brain", "Cytokine", "Heart", "Hormone", "Immune", 
                         "Kidney", "Liver", "Metab", "Lung", "MusculoSkeletal"]
        
        system_scores = pd.DataFrame(index=samples, columns=system_groups)
        
        for group in system_groups:
            # R: grepl(group, colnames(DNAmSystemPCs))
            # Find columns in dnam_sys_pcs that contain the group name
            relevant_cols = [c for c in dnam_sys_pcs.columns if group in c]
            
            sub_matrix = dnam_sys_pcs[relevant_cols]
            
            # Get coefficients for these columns
            # R: sub_system_coefficients <- RData$system_scores_coefficients_scale[tf]
            # In our parquet, these keys are stored in 'term' column of sys_score_coefs
            sub_coefs = self.sys_score_coefs.loc[relevant_cols]
            
            # R Logic:
            # if (length == 1) score = sub * -1
            # else score = sub %*% coefs
            if len(relevant_cols) == 1:
                # Vector multiplication by -1
                score = sub_matrix.iloc[:, 0] * -1.0
            else:
                # Matrix multiplication
                score = sub_matrix @ sub_coefs
            
            system_scores[group] = score

        # 5. Generate Predicted Chronological Age
        # R: DNAmPCs %*% Predicted_age_coefficients[2:4019] + Intercept
        # We use the names from age_linear_weights (PC1...PCn) to align
        
        # Ensure we have the right PCs
        common_age_pcs = dnam_pcs.columns.intersection(self.age_linear_weights.index)
        
        # Linear prediction
        age_pred_raw = dnam_pcs[common_age_pcs] @ self.age_linear_weights.loc[common_age_pcs] + self.age_intercept
        
        # Quadratic Correction
        # R: Age_pred * Model[2] + (Age_pred^2) * Model[3] + Model[1]
        # Model[1]=Const, Model[2]=Lin, Model[3]=Quad
        age_final = (age_pred_raw * self.age_model_lin) + \
                    (age_pred_raw**2 * self.age_model_quad) + \
                    self.age_model_const
        
        # Scale to Years (R code divides by 12 at this step)
        age_final_years = age_final / 12.0
        
        # Add to scores
        system_scores["Age_prediction"] = age_final_years
        
        # Rename columns to match R final names
        # R: c("Blood", "Brain", "Inflammation", "Heart", "Hormone", "Immune", 
        #      "Kidney", "Liver", "Metabolic", "Lung", "MusculoSkeletal", "Age_prediction")
        # Note name changes: Cytokine -> Inflammation, Metab -> Metabolic
        rename_map = {
            "Cytokine": "Inflammation",
            "Metab": "Metabolic"
        }
        system_scores = system_scores.rename(columns=rename_map)
        
        # 6. Generate Overall System Index (Composite SystemsAge)
        # PCA projection of the 12 scores
        
        # R: predict(RData$systems_PCA, system_scores)
        # This implies standardizing then rotating
        # (Scores - Center) / Scale @ Rotation
        
        # Ensure column order matches PCA center/rotation
        # The R code hardcoded the column order before PCA, strictly follow it:
        ordered_cols = ["Blood", "Brain", "Inflammation", "Heart", "Hormone", "Immune", 
                        "Kidney", "Liver", "Metabolic", "Lung", "MusculoSkeletal", "Age_prediction"]
        
        scores_for_pca = system_scores[ordered_cols]
        
        # Scale
        # R's predict.prcomp does (X - center) / scale
        # Note: Check if scale is FALSE in R object. If scale vector exists, use it.
        scores_centered = scores_for_pca - self.sys_pca_center
        
        # Check if scaling is needed (if scale vector all 1s or missing, skip division)
        # Based on R `predict`, if scale is present, it divides.
        if hasattr(self, 'sys_pca_scale') and not self.sys_pca_scale.empty:
             scores_scaled = scores_centered / self.sys_pca_scale
        else:
             scores_scaled = scores_centered

        # Rotate
        sys_pca_res = scores_scaled @ self.sys_pca_rot
        
        # Linear Combination for Final Age
        # R: system_PCA %*% Systems_clock_coefficients
        # Align PCs (PC1..PC12)
        final_pcs = sys_pca_res.columns.intersection(self.final_coefs.index)
        composite_age = sys_pca_res[final_pcs] @ self.final_coefs.loc[final_pcs]
        
        system_scores["SystemsAge"] = composite_age

        
        # 7. Final Scaling (Z-score to Years)
        # R Loop over 1:13 columns
        # Formula: (((y - V1) / V2) * V4) + V3
        # Then / 12
        
        # Map R column indices to names for clarity
        # The columns are the 11 systems + Age_prediction + SystemsAge
        final_columns = ordered_cols + ["SystemsAge"]
        
        # Ensure system_scores has these columns in order
        system_scores = system_scores[final_columns]

        # 遍历这 13 列
        for i, col_name in enumerate(final_columns):
            
            # 获取该列的数据
            y = system_scores[col_name]
            
            # --- 修正点：使用 iloc 按位置获取系数 ---
            # 因为 self.trans_coefs 的 index 是 '1'...'13'，名字对不上
            # 但顺序是对应的 (第0行对应 Blood, 第1行对应 Brain...)
            # .iloc[i] 获取第 i 行的所有系数 (V1, V2, V3, V4)
            coef_row = self.trans_coefs.iloc[i]
            
            # 获取 V1, V2, V3, V4
            # 假设 Parquet 读取后的列顺序是 V1, V2, V3, V4
            # 使用 iloc[0], iloc[1]... 确保即使列名乱了也能取对值
            v1 = coef_row.iloc[0]
            v2 = coef_row.iloc[1]
            v3 = coef_row.iloc[2]
            v4 = coef_row.iloc[3]
            
            # 执行变换公式: (((y - V1) / V2) * V4) + V3
            val = (((y - v1) / v2) * v4) + v3
            
            # 转换为年 (除以 12)
            val = val / 12.0
            
            # 存回 DataFrame
            system_scores[col_name] = val
            
        return system_scores


