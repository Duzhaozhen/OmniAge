#' @title Calculate Stochastic Epigenetic Clocks (Tong et al., 2024)
#'
#' @description
#' Calculates DNAm-Age estimates based on one of three stochastic epigenetic
#' clocks (StocH, StocZ, StocP). These are stochastic analogues of the
#' Horvath, Zhang, and Levine/PhenoAge clocks, respectively.
#'
#' They consist of the same CpGs as the original clocks but were trained on
#' artificial cohorts generated via a stochastic process of age-related DNAm
#' change accrual. Although built using sorted monocyte samples from the
#' MESA study, their predictive ability is largely independent of immune
#' cell type or whole blood, as the DNAm data is standardized prior to
#' clock application.
#'
#' The stochastic clocks may serve as a useful complement to the original
#' clocks, providing a method to assess if specific age-acceleration patterns
#' (or decelerations) could be driven by an underlying stochastic process.
#' This may offer insights into the biological mechanisms of phenotypes
#' linked to epigenetic age acceleration.
#'
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A \code{list} containing three numeric vectors: \code{StocH},
#'   \code{StocP}, and \code{StocZ}, representing predicted DNAm ages.
#'
#' @import glmnet
#' @importFrom stats coef
#' @export
#'
#' @references
#' Tong, H., Dwaraka, V.B., Chen, Q. et al.
#' Quantifying the stochastic component of epigenetic aging.
#' \emph{Nat Aging} (2024). \doi{10.1038/s43587-024-00636-6}
#'
#' @examples
#' # Load example data
#' hannumExample <- loadOmniAgeRdata("omniager_hannum_example")
#' stochClocksOut <- stochClocks(hannumExample[[1]])
stochClocks <- function(betaM, minCoverage = 0, verbose = TRUE) {
    stocAll <- loadOmniAgeRdata("omniager_stoch_clocks")

    stocNames <- paste0("Stoc", names(stocAll))
    mageList <- list()

    for (i in seq_along(stocAll)) {
        clockLabel <- stocNames[i]
        glmObj <- stocAll[[i]]

        fullCoef <- stats::coef(glmObj)
        intercept <- fullCoef[1, ncol(fullCoef)]
        weights <- fullCoef[-1, ncol(fullCoef)]

        coefLv <- list(intercept, weights)

        mageList[[clockLabel]] <- .calculateLinearPredictor(
            betaM = betaM,
            coefLv = coefLv,
            clockName = clockLabel,
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    return(mageList)
}
