#' @title  Estimate stemTOC score
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and
#' will return the stemTOC score.
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
#'
#' @details
#' The function will return the 0.95 upper quantile of the 371 stemTOC CpGs.
#' Compared to stemTOCvitro CpGs, the stemTOC CpGs are filtered for
#' significant DNA hypermethylation with chronological age
#' in large in-vivo datasets
#'
#' @return The stemTOC score of each sample.
#'
#' @references
#' Zhu, T., Tong, H., Du, Z. et al.
#' An improved epigenetic counter to track mitotic age in normal
#' and precancerous tissues.
#' \emph{Nat Commun} 2024
#'
#' @importFrom stats quantile
#' @examples
#' lungInv <- loadOmniAgeRdata(
#'     "omniager_lung_inv",
#'     verbose = FALSE
#' )
#' lungInvM <- lungInv$bmiq_m
#' stemtocOut <- stemTOC(betaM = lungInvM)
#'
#' @export
#'


stemTOC <- function(betaM, minCoverage = 0.5, verbose = TRUE) {
  
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
  
    stemTOCCpG <- loadOmniAgeRdata(
        "omniager_stemtoc_cpg",
        verbose = verbose
    )
    # Prepare the reference probe
    targetCpGs <- as.character(stemTOCCpG)
    clockWeights <- setNames(rep(1, length(targetCpGs)), targetCpGs)

    # Perform coverage check
    coverageResult <- .checkCpGCoverage(
        betaM = betaM,
        allWeights = clockWeights,
        clockName = "stemTOC",
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (!coverageResult$pass) {
        scores <- rep(NA_real_, ncol(betaM))
        names(scores) <- colnames(betaM)
        return(scores)
    }

    # 5. Calculate score
    scores <- apply(
        betaM[coverageResult$betaIdx, , drop = FALSE],
        2,
        quantile,
        probs = 0.95,
        na.rm = TRUE
    )

    return(scores)
}
