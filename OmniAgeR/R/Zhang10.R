#' @title Calculate Zhang10 DNAm Age (2017)
#'
#' @description A function to calculate the Zhang10 epigenetic clock age (2017)
#' from a DNA methylation beta value matrix.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding to
#' the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Zhang Y, Wilson R, Heiss J, et al.
#' DNA methylation signatures in peripheral blood strongly predict
#' all-cause mortality.
#' \emph{Nat Commun.} 2017
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' zhang10Out <- zhang10(hannumBmiqM)
#'
zhang10 <- function(betaM,
                    minCoverage = 0,
                    verbose = TRUE) {
    zhang10Coef <- loadOmniAgeRdata(
        "omniager_zhang10_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, zhang10Coef, "zhang10",
        minCoverage, verbose
    )

    return(predAgev)
}
