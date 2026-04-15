#' @title Calculate Horvath's DNAm Age (2013)
#'
#' @description A function to calculate the Horvath epigenetic clock age(2013)
#' from a DNA methylation beta value matrix.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @details Implements the Horvath (2013) pan-tissue clock. The function
#' calculates a weighted linear predictor from 353 CpGs found in the input
#' matrix and then transforms this value using a non-linear function to
#' return the final DNAm age.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding to
#' the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Horvath S.
#' DNA methylation age of human tissues and cell types.
#' \emph{Genome Biol.} 2013
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' horvath2013ClockOut <- horvath2013Clock(hannumBmiqM)
horvath2013Clock <- function(betaM,
                             minCoverage = 0,
                             verbose = TRUE) {
    horvath2013Coef <- loadOmniAgeRdata(
        "omniager_horvath2013_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, horvath2013Coef, "horvath2013Clock",
        minCoverage, verbose
    )
    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}
