#' @title Calculate Zhang10 DNAm Age (2017)
#'
#' @description A function to calculate the Zhang10 epigenetic clock age (2017) from a DNA methylation beta value matrix.
#'
#' @param beta.m A numeric matrix of beta values. Rows should be CpG probes and columns should be individual samples.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding to
#' the sample IDs from the input matrix's column names.
#'
#' @seealso The main function \code{\link{EpiAge}} can be used to calculate
#' multiple clocks simultaneously.
#'
#' @export
#'
#' @references
#' Zhang Y, Wilson R, Heiss J, et al.
#' DNA methylation signatures in peripheral blood strongly predict all-cause mortality.
#' \emph{Nat Commun.} 2017
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Zhang10.out <- Zhang10(hannum_bmiq_m)


Zhang10 <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("Zhang10Coef")
  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- 0

  # Coefficients
  coefficients <- as.numeric(as.vector(Zhang10Coef[1:nrow(Zhang10Coef), 2]))
  names(coefficients) <- as.vector(Zhang10Coef[1:nrow(Zhang10Coef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "Zhang10")

  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}



