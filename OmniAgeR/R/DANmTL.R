
#' @title The epigenetic age used for calculating the Leukocyte telomere length (2019).
#'
#' @description A function to calculate the the Leukocyte telomere length (2019) from a DNA methylation beta value matrix.
#'
#' @param beta.m A numeric matrix of beta values. Rows should be CpG probes and columns should be individual samples.
#'
#' @return A named vector of predicted Leukocyte telomere length.
#'
#' @export
#'
#' @references
#' Lu AT, Seeboth A, Tsai PC, et al.
#' DNA methylation-based estimator of telomere length
#' \emph{Aging} 2019
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' DNAmTL.o <- DNAmTL(hannum_bmiq_m)



DNAmTL <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("DNAmTLCoef") #

  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(DNAmTLCoef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(DNAmTLCoef[2:nrow(DNAmTLCoef), 2]))
  names(coefficients) <- as.vector(DNAmTLCoef[2:nrow(DNAmTLCoef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "DNAmTL")

  # --- Step 5: Return final age vector ---
  return(predage.v)
}
