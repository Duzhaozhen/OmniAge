#' @title Estimate RepliTali score
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and
#' will return the RepliTali score.
#' @details
#' The RepliTali model is described in Endicott et al. 2022. It is based on
#' 87 population doubling associated hypomethylated CpGs.
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

#' @return The RepliTali score of each sample.
#'
#'
#' @references
#' Endicott JL, Nolte PA, Shen H, Laird PW.
#' Cell division drives DNA methylation loss in late-replicating domains
#' in primary human cells.
#' \emph{Nat Commun.} 2022
#'
#' @examples
#' lungInv <- loadOmniAgeRdata(
#'     "omniager_lung_inv",
#'     verbose = FALSE
#' )
#' lungInvM <- lungInv$bmiq_m
#' repliTaliOut <- repliTali(betaM = lungInvM, minCoverage = 0)
#' @export
#'


repliTali <- function(betaM, minCoverage = 0.5, verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    replitaliCoef <- loadOmniAgeRdata(
        "omniager_replitali_coef",
        verbose = verbose
    )
    score <- .calLinearClock(
        betaM = betaM,
        coefData = replitaliCoef,
        clockLabel = "repliTali",
        minCoverage = minCoverage,
        verbose = verbose
    )

    return(score)
}
