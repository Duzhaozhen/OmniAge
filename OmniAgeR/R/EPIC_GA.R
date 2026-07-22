#' @title The Gestational Age (GA) clock based on 176 Illumina EPIC CpGs
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
#' @return A numeric vector containing one predicted gestational age per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned.
#'
#'
#' @export
#'
#' @references
#' Haftorn KL, Lee Y, Denault WRP, et al.
#' An EPIC predictor of gestational age and its application to newborns
#' conceived by assisted reproductive technologies.
#' \emph{Clin Epigenetics.} 2021
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' epicGaOut <- epicGa(hannumBmiqM, minCoverage = 0)
#'
epicGa <- function(betaM,
                   minCoverage = 0.5,
                   verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    epicGACoef <- loadOmniAgeRdata(
        "omniager_epic_ga_coef",
        verbose = verbose
    )
    predAgev <- .calLinearClock(
        betaM, epicGACoef, "epicGa",
        minCoverage, verbose
    )
    ## Convert the number of days into weeks
    predAgev <- predAgev / 7
    return(predAgev)
}
