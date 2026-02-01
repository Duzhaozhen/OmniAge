#' @title Calculate DNAm PhenoAge
#'
#' @description
#' CCalculates the DNA methylation PhenoAge, an epigenetic biomarker of aging, based on a matrix of DNA methylation beta values. The function implements the model originally developed by Levine et al. (2018), which uses a weighted linear combination of 513 specific CpG sites to predict phenotypic age.
#'
#' @details
#' This function calculates DNAm PhenoAge based on the model by Levine et al. (2018), which predicts phenotypic age by calculating a weighted sum of the beta values from 513 specific CpG sites.
#' The function automatically loads the required model coefficients, matches them with the CpGs in the input `beta.m` matrix, and computes the age.
#'
#' @param beta.m A numeric matrix of DNA methylation beta values. Rows should represent CpG sites and columns should represent individual samples.
#'
#' @return A named numeric vector containing the calculated PhenoAge for each sample.
#'
#' @references
#' Levine ME, Lu AT, Quach A, et al.
#' An epigenetic biomarker of aging for lifespan and healthspan.
#' \emph{Aging} 2018
#
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' PhenoAge.out <- PhenoAge(hannum_bmiq_m)

PhenoAge <- function(beta.m) {
  # --- 1. Data Loading and Preparation ---
  data("PhenoAgeCpG")

  # Transpose the input matrix to have samples as rows and CpGs as columns.
  # This format is required for subsequent matrix operations.
  beta.m <- t(beta.m)

  # --- 2. Check CpG Coverage ---
  # Get the list of CpGs required by the elastic net model (excluding the intercept).
  cpgs_in_model <- PhenoAge_coef_df$CpG[-1]
  common_cpgs <- intersect(colnames(beta.m), cpgs_in_model)
  print(paste0("[PhenoAge] Number of represented PhenoAge CpGs (max=",length(cpgs_in_model),")=",
             length(common_cpgs)))

  # --- 3. Subset ---
  # Subset the data to include only the CpGs common to both the input and the model.
  beta.m <- beta.m[, common_cpgs, drop = FALSE]

  # --- 4. Predict Age ---
  # Filter the coefficient data frame to match the available CpGs.
  coef_subset <- PhenoAge_coef_df[PhenoAge_coef_df $CpG %in% common_cpgs, ]
  coefs <- coef_subset$Coef
  names(coefs) <- coef_subset$CpG
  # Ensure the order of coefficients perfectly matches the column order of the scaled beta matrix.
  coefs <- coefs[colnames(beta.m)]
  # Extract the intercept value from the model coefficients.
  intercept_val <- PhenoAge_coef_df$Coef[1]

  # Calculate the predicted age using matrix multiplication and add the intercept.
  # Formula: PredictedAge = (beta_scaled * coefficients) + intercept
  predicted_age <- beta.m %*% coefs + intercept_val

  predicted_age <- as.vector(predicted_age)
  names(predicted_age) <- rownames(beta.m)
  return(predicted_age)
}
