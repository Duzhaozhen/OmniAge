#' @title Calculate Hannum's DNAm Age (2013)
#'
#' @description A function to calculate the Hannum epigenetic clock age (2013)
#' from a DNA methylation beta value matrix.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @details
#' Implements the Hannum (2013) blood-specific clock. The function calculates
#' a weighted linear predictor from 71 CpGs found in the input matrix.
#' This clock is a direct linear model without non-linear transformation.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding to
#' the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Hannum G, Guinney J, Zhao L, et al.
#' Genome-wide methylation profiles reveal quantitative views of human aging rates.
#' \emph{Mol Cell.} 2013
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' hannumClockOut <- hannumClock(hannumBmiqM)
#'
hannumClock <- function(betaM,
                        minCoverage = 0,
                        verbose = TRUE) {
    hannumCoef <- loadOmniAgeRdata(
        "omniager_hannum",
        verbose = verbose
    )
    predAgev <- .calLinearClock(
        betaM, hannumCoef, "hannumClock",
        minCoverage, verbose
    )
    return(predAgev)
}
