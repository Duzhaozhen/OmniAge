#' @title Calculate HepatoXu ctDNA Methylation Scores
#' for Hepatocellular Carcinoma
#'
#' @description
#' This function implements the diagnostic prediction models for
#' Hepatocellular Carcinoma (HCC) based on circulating tumour DNA (ctDNA)
#' methylation markers as described by Xu et al. (2017)
#'
#' @details
#' The function calculates a composite score using a panel of HCC-specific
#' methylation markers identified through Random Forest and LASSO regression
#' analysis.
#'
#' For diagnosis, the model (cd-score) utilizes 10 genomic markers plus a
#' logistic regression intercept to differentiate HCC patients from healthy
#' controls or those with non-malignant liver diseases.
#'
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector containing the calculated methylation scores
#' (cd-score) for each sample.
#'
#'
#' @export
#'
#' @references
#' Xu, Rh., Wei, W., Krawczyk, M. et al.
#' Circulating tumour DNA methylation markers for diagnosis and prognosis of
#' hepatocellular carcinoma.
#' \emph{Nature Mater} 2017
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' hepatoXuRiskO <- hepatoXuRisk(hannumBmiqM)
#'
hepatoXuRisk <- function(betaM,
                         minCoverage = 0,
                         verbose = TRUE) {
    hepatoXuCoef <- loadOmniAgeRdata(
        "omniager_hepato_xu_coef",
        verbose = verbose
    )
    predAgev <- .calLinearClock(
        betaM, hepatoXuCoef, "hepatoXuRisk",
        minCoverage, verbose
    )
    return(predAgev)
}
