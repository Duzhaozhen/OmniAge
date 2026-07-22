#' @title The PedBE (Pediatric Buccal) Clock for DNAm Age in Children
#'
#' @description
#' Implements the Pediatric Buccal Epigenetic (PedBE) clock, specifically
#' developed to estimate DNA methylation age in
#' **pediatric (childhood) samples**, as described by McEwen et al. (2020).
#'
#' @details
#' This clock is specifically trained on and designed for pediatric buccal
#' epithelial (cheek swab) samples. The calculation is a two-step process:
#'
#' 1.  A linear predictor is first calculated from the beta values using the
#' 94-CpG elastic net coefficients . This value represents a *transformed* age.
#' 2.  A non-linear inverse age transformation (via the transformation function
#' from Horvath, 2013) is then applied. This converts the transformed age into
#' a final estimate of chronological age in years.
#'
#' @param betaM A numeric DNA methylation beta-value matrix with CpG probe
#'   identifiers as row names and samples as columns. CpG identifiers in
#'   \code{rownames(betaM)} are required. Sample identifiers in
#'   \code{colnames(betaM)} are optional. The matched CpG values used for
#'   calculation must not contain missing values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' 
#' 
#' @return A numeric vector containing one predicted age per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned.
#'   
#' @export
#'
#' @references
#' McEwen LM, O'Donnell KJ, McGill MG, et al.
#' The PedBE clock accurately estimates DNA methylation age in
#' pediatric buccal cells.
#' \emph{Proc Natl Acad Sci U S A.} 2020
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' pedBEClockOut <- pedBEClock(hannumBmiqM)
#'
pedBEClock <- function(betaM,
                       minCoverage = 0.5,
                       verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    pedBECoef <- loadOmniAgeRdata(
        "omniager_pedbe_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, pedBECoef, "pedBEClock",
        minCoverage, verbose
    )
    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}
