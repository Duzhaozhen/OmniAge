#' @title Calculate Horvath's Skin & Blood DNAm Age (2018)
#'
#' @description
#' A function to calculate the Horvath "Skin & Blood" clock age (2018)
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
#'
#' @details Implements the Horvath (2013) pan-tissue clock. The function
#' calculates a weighted linear predictor from 353 CpGs found in the input
#' matrix and then transforms this value using a non-linear function to
#' return the final DNAm age.
#'
#' @details
#' Implements the Horvath (2018) skin & blood clock. The function calculates
#' a weighted linear predictor from 391 CpGs and then transforms this value
#' using the same non-linear function as the Horvath (2013) clock to
#' return the final DNAm age.
#'
#' @return A numeric vector containing one predicted age per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned.
#'
#'
#' @export
#'
#' @references
#' Horvath S, Oshima J, Martin GM, et al.
#' Epigenetic clock for skin and blood cells applied to Hutchinson Gilford
#' Progeria Syndrome and ex vivo studies.
#' \emph{Aging} 2018
#'
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' horvath2018ClockOut <- horvath2018Clock(hannumBmiqM)
horvath2018Clock <- function(betaM,
                             minCoverage = 0.5,
                             verbose = TRUE) {
  
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    
    horvath2018Coef <- loadOmniAgeRdata(
        "omniager_horvath2018_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, horvath2018Coef, "horvath2018Clock",
        minCoverage, verbose
    )
    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}
