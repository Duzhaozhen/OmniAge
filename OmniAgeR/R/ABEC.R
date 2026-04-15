#' @title Adult Blood-based EPIC Clock (ABEC)
#'
#' @description
#' Predicts biological age using the Adult Blood-based EPIC Clock (ABEC).
#' Developed by Lee et al., this model was trained on DNA methylation (DNAm)
#' data from the Norwegian Mother, Father and Child Cohort Study (MoBa)
#' (n = 1,592, age range: 19–59 years) using the Illumina EPIC platform.
#'
#' @details
#' The function extracts the necessary CpG coefficients from the internal
#' \code{ABEC_Coef} dataset and applies them to the provided beta value matrix.
#' It uses an internal helper to handle missing probes and compute the
#' final age estimates.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#'
#' @return A numeric vector of predicted biological ages. The vector is
#' named using the sample IDs from the \code{rownames} of \code{betaM}.
#'
#' @export
#'
#' @references
#' Lee, Y., Haftorn, K.L., Denault, W.R.P. et al.
#' Blood-based epigenetic estimators of chronological age in human adults
#' using DNA methylation data from the Illumina MethylationEPIC array.
#' \emph{BMC Genomics} 2020
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' abecOut <- leeABEC(hannumBmiqM)
#'
leeABEC <- function(betaM,
                    minCoverage = 0,
                    verbose = TRUE) {
    abecCoef <- loadOmniAgeRdata("omniager_abec_coef", verbose = verbose)
    return(.calLinearClock(betaM, abecCoef, "leeABEC", minCoverage, verbose))
}


#' @title Extended Adult Blood-based EPIC Clock (eABEC)
#'
#' @description
#' Predicts biological age using the Extended Adult Blood-based EPIC Clock
#' (eABEC). This model extends the training set of ABEC by incorporating
#' public data from the Gene Expression Omnibus (GEO), resulting in a
#' broader age-span (n = 2,227, age range: 18–88 years).
#'
#' @inheritParams leeABEC
#' @inherit leeABEC return
#'
#' @details
#' Similar to \code{leeABEC}, this function utilizes the \code{eABEC_Coef}
#' dataset. It is designed for applications where a wider range of adult
#' ages is expected.
#'
#' @export
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' eabecOut <- leeExtendedABEC(hannumBmiqM)
#'
leeExtendedABEC <- function(betaM,
                            minCoverage = 0,
                            verbose = TRUE) {
    eabecCoef <- loadOmniAgeRdata("omniager_eabec_coef", verbose = verbose)
    return(.calLinearClock(betaM, eabecCoef, "leeExtendedABEC", minCoverage, verbose))
}


#' @title Common Adult Blood-based EPIC Clock (cABEC)
#'
#' @description
#' Predicts biological age using the Common Adult Blood-based EPIC Clock
#' (cABEC). This model uses the same extended training set as \code{eABEC}
#' but is restricted to CpGs common to both Illumina 450K and EPIC arrays,
#' ensuring backward compatibility and robustness across platforms.
#'
#' @inheritParams leeABEC
#' @inherit leeABEC return
#'
#' @details
#' The function uses coefficients from the \code{cABEC_Coef} dataset.
#'
#' @export
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' cabecOut <- leeCommonABEC(hannumBmiqM)
#'
leeCommonABEC <- function(betaM,
                          minCoverage = 0,
                          verbose = TRUE) {
    cabecCoef <- loadOmniAgeRdata("omniager_cabec_coef", verbose = verbose)
    return(.calLinearClock(betaM, cabecCoef, "leeCommonABEC", minCoverage, verbose))
}
