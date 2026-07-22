#' @title Calculate Hannum's DNAm Age (2013)
#'
#' @description A function to calculate the Hannum epigenetic clock age (2013)
#' from a DNA methylation beta value matrix.
#'
#' @param betaM A numeric DNA methylation beta-value matrix with CpG probe
#'   identifiers as row names and samples as columns. CpG identifiers in
#'   \code{rownames(betaM)} are required. Sample identifiers in
#'   \code{colnames(betaM)} are optional. The matched CpG values used for
#'   calculation must not contain missing values.
#'   
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @details
#' Implements the Hannum (2013) blood-specific clock. The function calculates
#' a weighted linear predictor from 71 CpGs found in the input matrix.
#' This clock is a direct linear model without non-linear transformation.
#'
#' @return A numeric vector containing one predicted age per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned.
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
                        minCoverage = 0.5,
                        verbose = TRUE) {
    
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
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
