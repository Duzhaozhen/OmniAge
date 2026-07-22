#' @title Calculate the Vidal-Bralo Epigenetic Clock (8-CpG Model)
#'
#' @description
#' Implements the simplified 8-CpG epigenetic clock for estimating
#' chronological age in adults, as described by Vidal-Bralo et al. (2016).
#'
#' @details
#' This function calculates the highly simplified epigenetic clock developed
#' for use in adult whole blood samples. The model is a linear predictor
#' based on 8 specific CpG sites.
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
#' Vidal-Bralo L, Lopez-Golan Y, Gonzalez A.
#' Simplified Assay for Epigenetic Age Estimation in Whole Blood of Adults.
#' \emph{Front Genet.} 2016
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' vidalBraloClockOut <- vidalBraloClock(hannumBmiqM)
vidalBraloClock <- function(betaM,
                            minCoverage = 0.5,
                            verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    vidalBraloCoef <- loadOmniAgeRdata(
        "omniager_vidalbralo_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, vidalBraloCoef, "vidalBraloClock",
        minCoverage, verbose
    )
    return(predAgev)
}
