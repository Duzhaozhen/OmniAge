#' @title Wu's Epigenetic Clock for Pediatric Age Estimation
#'
#' @description
#' Estimates biological age (in years) based on DNA methylation data using the method
#' proposed by Wu et al. (2019). This clock is specifically designed for pediatric
#' cohorts and utilizes a non-linear transformation to capture rapid developmental
#' changes in early life.
#'
#' @details
#' The calculation involves two main steps:
#' \enumerate{
#'   \item Calculation of a linear predictor using weighted CpG beta values.
#'   \item Transformation of the linear predictor into biological age (years) using
#'   a specific "anti-transformation" function with a toddler age offset of 48 months.
#' }
#'
#' \strong{Data Requirements:}
#' The input \code{beta.m} must be a matrix of Beta values (0 to 1). The function
#' expects CpG probes as row names.
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
#' @return A numeric vector of predicted biological ages. The vector is
#' named using the sample IDs from the \code{rownames} of \code{betaM}.
#'
#' @export
#'
#' @references
#' Wu, Xiaohui et al.
#' DNA methylation profile is a quantitative measure of biological aging in children
#' \emph{Aging} 2019
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' wuClockOut <- wuClock(hannumBmiqM)
wuClock <- function(betaM,
                    minCoverage = 0.5,
                    verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    wuClockCoef <- loadOmniAgeRdata(
        "omniager_wu_clock_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, wuClockCoef, "wuClock",
        minCoverage, verbose
    )
    predAgev <- .antiTrafo(predAgev, 48)
    ## transform to years
    predAgev <- predAgev / 12
    return(predAgev)
}
