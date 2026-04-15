#' @title Estimate RepliTali score
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and
#' will return the RepliTali score.
#' @details
#' The RepliTali model is described in Endicott et al. 2022. It is based on
#' 87 population doubling associated hypomethylated CpGs.
#'
#' @param betaM A numeric matrix of DNAm beta values (probes as rows). Rows
#' should be Illumina 450k/EPIC CpG identifiers and columns should be samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
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


repliTali <- function(betaM, minCoverage = 0, verbose = TRUE) {
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
