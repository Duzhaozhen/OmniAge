#' @title The epigenetic age used for calculating the Leukocyte telomere length.
#'
#' @description A function to calculate the the Leukocyte telomere length (2019)
#' from a DNA methylation beta value matrix.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
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
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' dnamTlO <- dnamTL(hannumBmiqM)
dnamTL <- function(betaM,
                   minCoverage = 0,
                   verbose = TRUE) {
    dnamTLCoef <- loadOmniAgeRdata(
        "omniager_dnamtl_coef",
        verbose = verbose
    )
    predAgev <- .calLinearClock(
        betaM, dnamTLCoef, "dnamTL",
        minCoverage, verbose
    )
    return(predAgev)
}
