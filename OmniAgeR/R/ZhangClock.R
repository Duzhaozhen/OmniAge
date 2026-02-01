#' @title Predicts DNA methylation age using the Zhang clock(Elastic Net model)
#'#'
#' @description
#' This function takes a matrix of DNA methylation beta values and calculates the epigenetic age for each sample based on the elastic net model developed by Zhang et al. (2019). n.
#'
#' @param beta.m DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @details
#' This function implements the elastic net epigenetic clock from Zhang et al. (2019). A unique feature of this clock is the pre-processing step. Instead of using raw beta values, it first performs a per-sample standardization. Specifically, it calculates a Z-score for each CpG based on the mean and standard deviation of all measured CpGs within that same sample. The final age is then predicted by taking the weighted sum of these standardized values using the model's pre-defined coefficients.
#'
#' @return A numeric vector of predicted ages, with sample names preserved.
#'
#' @seealso The main function \code{\link{EpiAge}} can be used to calculate multiple clocks simultaneously.
#'
#' @references
#' Zhang Q, Vallerga CL, Walker RM, et al.
#' Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing.
#' \emph{Genome Med.} 2019
#'
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' ZhangClock.out <- ZhangClock(hannum_bmiq_m)
#'


ZhangClock <- function(beta.m) {

  # --- 1. Data Loading and Preparation ---
  data("ZhangCpG")

  # Transpose the input matrix to have samples as rows and CpGs as columns.
  # This format is required for subsequent matrix operations.
  beta.m <- t(beta.m)

  # --- 2. Check CpG Coverage ---
  # Get the list of CpGs required by the elastic net model (excluding the intercept).
  cpgs_in_model <- Zhang_coef_df $CpG[-1]
  common_cpgs <- intersect(colnames(beta.m), cpgs_in_model)
  print(paste0("[ZhangClock] Number of represented Zhang CpGs (max=",length(cpgs_in_model),")=",
             length(common_cpgs)))

  # --- 3. Per-Sample Standardization ---
  ############# for each probe, change to missing value to the mean value across all individuals #############
  # addna<-function(methy){
  #   methy[is.na(methy)]<-mean(methy,na.rm=T)
  #   return(methy)
  # }

  # print("1.2 Replacing missing values with mean value")
  # dataNona<-apply(data,2,function(x) addna(x))   ###############  replace the NA with mean value for each probe
  colCpGs <- colnames(beta.m)
  cpgs_scaled <- t(apply(beta.m, 1, scale)) # 按行(样本)进行标准化
  colnames(cpgs_scaled) <- colCpGs
  # Subset the scaled data to include only the CpGs common to both the input and the model.
  cpgs_scaled <- cpgs_scaled[, common_cpgs, drop = FALSE]

  # --- 4. Predict Age ---
  # Filter the coefficient data frame to match the available CpGs.
  coef_subset <- Zhang_coef_df [Zhang_coef_df $CpG %in% common_cpgs, ]
  coefs <- coef_subset$Coef
  names(coefs) <- coef_subset$CpG
  # Ensure the order of coefficients perfectly matches the column order of the scaled beta matrix.
  coefs <- coefs[colnames(cpgs_scaled)]
  # Extract the intercept value from the model coefficients.
  intercept_val <- Zhang_coef_df $Coef[1]

  # Calculate the predicted age using matrix multiplication and add the intercept.
  # Formula: PredictedAge = (beta_scaled * coefficients) + intercept
  predicted_age <- cpgs_scaled %*% coefs + intercept_val

  predicted_age <- as.vector(predicted_age)
  names(predicted_age) <- rownames(cpgs_scaled)
  return(predicted_age)
}
