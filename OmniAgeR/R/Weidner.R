#' @title Calculate Weidner Epigenetic Age (3-CpG Blood Clock)
#'
#' @description
#' Estimates biological age using the specific 3-CpG signature described by
#' Weidner et al. (2014). This model is a multivariate linear regression
#' originally designed for pyrosequencing data derived from blood samples,
#' but it can also be applied to microarray data.
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
#' @return A numeric vector containing one predicted chronological age per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned.
#'
#' @export
#'
#' @references
#' Weidner, C.I., Lin, Q., Koch, C.M. et al.
#' Aging of blood can be tracked by DNA methylation changes
#' at just three CpG sites.
#' \emph{Genome Biol} 2014
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' weidnerClockOut <- weidnerClock(hannumBmiqM)
weidnerClock <- function(betaM,
                         minCoverage = 0.5,
                         verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    
    weidnerCoef <- loadOmniAgeRdata(
        "omniager_weidner_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, weidnerCoef, "weidnerClock",
        minCoverage, verbose
    )
    return(predAgev)
}
