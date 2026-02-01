#' @title Calculate Hannum's DNAm Age (2013)
#'
#' @description A function to calculate the Hannum epigenetic clock age (2013) from a DNA methylation beta value matrix.
#'
#' @param beta.m A numeric matrix of beta values. Rows should be CpG probes and columns should be individual samples.
#'
#' @details
#' Implements the Hannum (2013) blood-specific clock. The function calculates
#' a weighted linear predictor from 71 CpGs found in the input matrix.
#' This clock is a direct linear model without non-linear transformation.
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
#' Hannum G, Guinney J, Zhao L, et al.
#' Genome-wide methylation profiles reveal quantitative views of human aging rates.
#' \emph{Mol Cell.} 2013
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Hannum.out <- Hannum(hannum_bmiq_m)


Hannum <- function(beta.m) {
  # Step 1: Load the coefficients
  data("Hannum")
  hannumP.lv <- list()
  # Intercept:
  hannumP.lv[[1]] <- 0
  # Coefficients:
  coefficients <- as.numeric(Hannum_coef$CoefficientTraining)
  names(coefficients) <- Hannum_coef$CpGmarker
  hannumP.lv[[2]] <- coefficients
  # Step 2: Define the necessary helper function inside the main function.

  # Prediction of age function (无 'T'age, 因为它是直接预测)
  PredAge <- function(beta.m, hannumP.lv) {

    # Match required CpGs with those in the input data
    map.idx <- match(names(hannumP.lv[[2]]), rownames(beta.m))
    rep.idx <- which(!is.na(map.idx))

    print(paste0(
      "[Hannum] Number of represented Hannum CpGs (max=",
      length(hannumP.lv[[2]]),
      ")=",
      length(rep.idx)
    ))

    # Subset beta matrix to only the matched CpGs
    tmpB.m <- beta.m[map.idx[rep.idx], , drop = FALSE]

    # Calculate the weighted sum for each sample
    intercept <- hannumP.lv[[1]] # 这将是 0
    weights <- hannumP.lv[[2]][rep.idx]

    # Using matrix multiplication for efficiency
    # age.v 是 (Intercept=0) + (t(Beta) %*% Weights)
    age.v <- as.vector(intercept + t(tmpB.m) %*% weights)

    return(age.v)
  }


  # Step 3: Execute the prediction
  predMage.v <- PredAge(beta.m, hannumP.lv)

  # Step 4: Add sample names to the output vector and return it
  names(predMage.v) <- colnames(beta.m)
  return(predMage.v)
}
