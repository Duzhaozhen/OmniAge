#' @title Calculate Zhang10 DNAm Age (2017)
#'
#' @description A function to calculate the Zhang10 epigenetic clock age (2017)
#' from a DNA methylation beta value matrix.
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
#' @return A numeric vector containing one predicted DNAm age per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned.
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
                    minCoverage = 0.5,
                    verbose = TRUE) {
  
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
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
