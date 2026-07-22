#' @title IntrinClock Age Prediction
#'
#' @description
#' Calculates the "IntrinClock" epigenetic age (intrinsic cellular age) using
#' DNA methylation data(Illumina 450K and EPIC).
#' This clock is designed to be resistant to changes in immune cell composition.
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
#' @details
#' The IntrinClock utilizes an elastic net model trained on 410 CpGs and
#' refined to 381 active predictors. It predicts age by calculating a linear
#' combination of CpG beta values and coefficients, followed by an inverse
#' Horvath transformation to convert the linear predictor to years.
#'
#' The model coefficients are stored in \code{IntrinClockCoef}, which contains:
#' \itemize{
#'   \item Intercept term
#'   \item 381 non-zero CpG coefficients (sparsity optimized)
#' }
#'
#'
#' @return A numeric vector containing one predicted age per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned
#'
#' @export
#'
#' @references
#' Tomusiak, A., Floro, A., Tiwari, R. et al.
#' Development of an epigenetic clock resistant to changes in
#' immune cell composition.
#' \emph{Commun Biol} 2024
#'
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' intrinClockO <- intrinClock(hannumBmiqM)
intrinClock <- function(betaM,
                        minCoverage = 0.5,
                        verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    
    intrinClockCoef <- loadOmniAgeRdata(
        "omniager_intrin_clock_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, intrinClockCoef, "intrinClock",
        minCoverage, verbose
    )

    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}
