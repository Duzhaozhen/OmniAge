#' @title Calculate Lin DNAm Age (2016)
#'
#' @description A function to calculate the Lin epigenetic clock age (2016)
#' from a DNA methylation beta value matrix.
#'
#' @param betaM A matrix of beta values (CpGs in rows, samples in columns).
#' This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding to
#' the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Lin Q, Weidner CI, Costa IG, et al.
#' DNA methylation levels at individual age-associated CpG sites can be
#' indicative for life expectancy.
#' \emph{Aging} 2016
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' linClockOut <- linClock(betaM = hannumBmiqM)
linClock <- function(betaM,
                     minCoverage = 0,
                     verbose = TRUE) {
    linCoef <- loadOmniAgeRdata(
        "omniager_lin_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, linCoef, "linClock",
        minCoverage, verbose
    )
    return(predAgev)
}
