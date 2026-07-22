#' @title Estimate EpiCMIT scores
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will
#' return EpiCMIT scores.
#'
#' @details
#' EpiCMIT is based on two groups of CpGs, a 184 age-associated hypermethylated
#' CpG list and a 1164 hypomethylated CpG list (Duran-Ferrer et al. 2020) .
#' The EpiCMIT function calculates the average DNAm of the two groups separately
#' and returns the EpiCMIT_hyper and EpiCMIT_hypo scores.
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
#'
#' @return A list containing the following entries.
#'
#' * hyperScore: EpiCMIT-hyper score of each sample.
#'
#' * hypoScore: EpiCMIT-hypo score of each sample.
#'
#' @references
#' Duran-Ferrer M, Clot G, Nadeu F, et al.
#' The proliferative history shapes the DNA methylome of B-cell tumors
#' and predicts clinical outcome.
#' \emph{Nat Cancer} 2020
#'
#' @examples
#' lungInv <- loadOmniAgeRdata(
#'     "omniager_lung_inv",
#'     verbose = FALSE
#' )
#' lungInvM <- lungInv$bmiq_m
#' epicmitOut <- epiCMIT(betaM = lungInvM)
#' @export
#'

epiCMIT <- function(betaM, minCoverage = 0.5, verbose = TRUE) {
    
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )  
  
    epiCMITdf <- loadOmniAgeRdata(
        "omniager_epicmit_model",
        verbose = verbose
    )

    # 1. Prepare the list of probes
    hyperProbes <- epiCMITdf[grep("hyper", epiCMITdf$epiCMIT.class), 1]
    hypoProbes <- epiCMITdf[grep("hypo", epiCMITdf$epiCMIT.class), 1]

    hyperWeights <- setNames(rep(1, length(hyperProbes)), hyperProbes)
    hypoWeights <- setNames(rep(1, length(hypoProbes)), hypoProbes)

    # 2. Call the helper function
    resHyper <- .checkCpGCoverage(
        betaM, hyperWeights, "EpiCMIT-Hyper",
        minCoverage, verbose
    )
    resHypo <- .checkCpGCoverage(
        betaM, hypoWeights, "EpiCMIT-Hypo",
        minCoverage, verbose
    )

    # 3. Predcition
    hyperScore <- if (resHyper$pass) {
        colMeans(betaM[resHyper$betaIdx, , drop = FALSE])
    } else {
        rep(NA_real_, ncol(betaM))
    }

    hypoScore <- if (resHypo$pass) {
        1 - colMeans(betaM[resHypo$betaIdx, , drop = FALSE])
    } else {
        rep(NA_real_, ncol(betaM))
    }

    return(list(
        hyperScore = hyperScore,
        hypoScore  = hypoScore
    ))
}
