#' @title Calculate a DNA Methylation-Based Proxy for IL-6
#'
#' @description
#' Computes a DNA methylation (DNAm) surrogate score for Interleukin-6 (IL-6) protein levels. IL-6 is a key cytokine and biomarker for systemic inflammation.
#'
#' @details
#' This function calculates the IL-6 proxy score by applying a pre-defined
#' set of coefficients to the input beta-value matrix. The coefficients
#' are loaded from the internal `IL6Coef` data object.
#'
#' @param nbeta.m normalized beta-valued DNAm data matrix with rownames the CpG identifiers. Missing values should be imputed.
#'
#' @return A numeric vector containing the calculated IL-6 proxy score for each sample. The vector is named according to the column names (sample IDs) of the input `beta.m` matrix.
#'
#'
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' CompIL6.out <- CompIL6(hannum_bmiq_m)



CompIL6 <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("IL6Coef")
  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- 0

  # Coefficients
  coefficients <- as.numeric(as.vector(IL6Coef[, 2]))
  names(coefficients) <- as.vector(IL6Coef[, 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  IL6_score <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "IL6_Score")

  return(IL6_score)
}
