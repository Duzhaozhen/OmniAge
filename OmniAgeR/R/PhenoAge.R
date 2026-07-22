#' @title Calculate DNAm PhenoAge
#'
#' @description
#' CCalculates the DNA methylation PhenoAge, an epigenetic biomarker
#' of aging, based on a matrix of DNA methylation beta values. The
#' function implements the model originally developed by Levine et al.
#' (2018), which uses a weighted linear combination of 513 specific
#' CpG sites to predict phenotypic age.
#'
#' @details
#' This function calculates DNAm PhenoAge based on the model by
#' Levine et al. (2018), which predicts phenotypic age by
#' calculating a weighted sum of the beta values from 513 specific CpGs.
#' The function automatically loads the required model coefficients,
#' matches them with the CpGs in the input matrix, and computes the age.
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
#' @return A numeric vector containing one predicted age per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned.
#'   
#'
#' @references
#' Levine ME, Lu AT, Quach A, et al.
#' An epigenetic biomarker of aging for lifespan and healthspan.
#' \emph{Aging} 2018
#
#' @export
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' phenoAgeOut <- phenoAge(hannumBmiqM)
#'
phenoAge <- function(betaM,
                     minCoverage = 0.5,
                     verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    phenoAgeCoef <- loadOmniAgeRdata(
        "omniager_phenoage_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, phenoAgeCoef, "phenoAge",
        minCoverage, verbose
    )
    return(predAgev)
}
