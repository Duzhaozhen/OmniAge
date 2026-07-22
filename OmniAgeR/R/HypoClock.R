#' @title Estimate hypoClock score
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will
#' return the HypoClock score.
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
#' The hypoClock score is defined by Teschendorff (2020) and is calculated
#' as 1 minus the average beta value of 678 solo-WCGW sites.
#'
#' @return The HypoClock score of each sample.
#'
#' @references
#' Teschendorff AE.
#' A comparison of epigenetic mitotic-like clocks for cancer risk prediction.
#' \emph{Genome Med.} 2020
#'
#' @examples
#' lungInv <- loadOmniAgeRdata(
#'     "omniager_lung_inv",
#'     verbose = FALSE
#' )
#' lungInvM <- lungInv$bmiq_m
#' hypoClockOut <- hypoClock(betaM = lungInvM)
#'
#' @export
#'

hypoClock <- function(betaM, minCoverage = 0.5, verbose = TRUE) {
  
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
  
    hypoClockCpG <- loadOmniAgeRdata(
        "omniager_hypoclock_cpg",
        verbose = verbose
    )

    # Prepare reference sites
    soloCpGs <- as.character(hypoClockCpG)
    clockWeights <- setNames(rep(1, length(soloCpGs)), soloCpGs)

    # Perform coverage check
    coverageResult <- .checkCpGCoverage(
        betaM = betaM,
        allWeights = clockWeights,
        clockName = "hypoClock",
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (!coverageResult$pass) {
        scores <- rep(NA_real_, ncol(betaM))
        names(scores) <- colnames(betaM)
        return(scores)
    }

    # 5. Calculate score(1 - mean beta)
    scores <- 1 - colMeans(betaM[coverageResult$betaIdx, , drop = FALSE],
        na.rm = TRUE
    )

    return(scores)
}
