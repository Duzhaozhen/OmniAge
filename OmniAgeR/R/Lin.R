#' @title Calculate Lin DNAm Age (2016)
#'
#' @description A function to calculate the Lin epigenetic clock age (2016) from a DNA methylation beta value matrix.
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
#' Lin Q, Weidner CI, Costa IG, et al.
#' DNA methylation levels at individual age-associated CpG sites can be indicative for life expectancy.
#' \emph{Aging} 2016
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Lin.out <- Lin(beta.m = hannum_bmiq_m)
#'

Lin <- function(beta.m) {
  # Step 1: Load the coefficients
  data("Lin")
  LinP.lv <- list()
  # Intercept:
  LinP.lv[[1]] <- Lin_coef$CoefficientTraining[length(Lin_coef$CoefficientTraining)]
  # Coefficients:
  coefficients <- as.numeric(Lin_coef$CoefficientTraining[-length(Lin_coef$CoefficientTraining)])
  names(coefficients) <- Lin_coef$CpGmarker[-length(Lin_coef$CoefficientTraining)]
  LinP.lv[[2]] <- coefficients
  # Step 2: predict age
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = LinP.lv,
                                        clock.name = "Lin")

  return(predage.v)
}

