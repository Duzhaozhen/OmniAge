#' @title Calculate a DNA Methylation-Based Proxy for IL-6
#'
#' @description
#' Computes a DNAm surrogate score for Interleukin-6 (IL-6)
#' protein levels.
#'
#' @details
#' This function calculates the IL-6 proxy score by applying a pre-defined
#' set of coefficients to the input beta-value matrix.
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
#' @return A numeric vector containing one prediction per
#'   sample, in the same order as the columns of \code{betaM}. If
#'   \code{betaM} has column names, these are retained as the names of the
#'   returned vector; otherwise, an unnamed numeric vector is returned.
#'
#' @references
#' Stevenson AJ et al.
#' Creating and Validating a DNA Methylation-Based Proxy for Interleukin-6.
#' \emph{J Gerontol A Biol Sci Med Sci.} 2021
#'
#' @export
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' compil6Out <- compIL6(hannumBmiqM)
compIL6 <- function(betaM,
                    minCoverage = 0.5,
                    verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    # --- Step 1: Load and parse coefficients ---
    iL6Coef <- loadOmniAgeRdata(
        "omniager_il6_coef",
        verbose = verbose
    )

    il6Score <- .calLinearClock(
        betaM, iL6Coef, "compIL6",
        minCoverage, verbose
    )

    return(il6Score)
}
