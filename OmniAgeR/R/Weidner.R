#' @title Calculate Weidner Epigenetic Age (3-CpG Blood Clock)
#'
#' @description
#' Estimates biological age using the specific 3-CpG signature described by
#' Weidner et al. (2014). This model is a multivariate linear regression
#' originally designed for pyrosequencing data derived from blood samples,
#' but it can also be applied to microarray data.
#'
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A **numeric vector** containing the predicted DNAm age
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `betaM`.
#'
#' @export
#'
#' @references
#' Weidner, C.I., Lin, Q., Koch, C.M. et al.
#' Aging of blood can be tracked by DNA methylation changes
#' at just three CpG sites.
#' \emph{Genome Biol} 2014
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' weidnerClockOut <- weidnerClock(hannumBmiqM)
weidnerClock <- function(betaM,
                         minCoverage = 0,
                         verbose = TRUE) {
    weidnerCoef <- loadOmniAgeRdata(
        "omniager_weidner_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, weidnerCoef, "weidnerClock",
        minCoverage, verbose
    )
    return(predAgev)
}
