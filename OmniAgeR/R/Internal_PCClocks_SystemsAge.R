#' Preprocess and Impute DNA Methylation Data for PC Clocks
#'
#' @description
#' This internal function ensures the input methylation matrix contains all
#' required CpGs for PC-based clock projection. It performs a hierarchical
#' quality control and imputation:
#' \enumerate{
#'   \item Checks if the proportion of found CpGs meets the \code{minCoverage}
#'   threshold.
#'   \item Performs column-mean imputation for any \code{NA} values within
#'   present CpGs.
#'   \item Appends entirely missing CpGs using reference values (medians/means)
#'         provided in the \code{requiredCpGs} vector.
#' }
#'
#' @param betaM A numeric matrix of DNA methylation values.
#' @param requiredCpGs A named numeric vector of required CpG sites and their
#'   reference values.
#' @param minCoverage A numeric value (0-1) representing the minimum required
#'   proportion of CpGs.
#' @param verbose A logical flag to print progress messages.
#'
#' @return A numeric matrix with dimensions (N_samples x N_required_CpGs),
#'   where all columns match the \code{requiredCpGs} order. Returns \code{NULL}
#'   if the coverage check fails.
#'
#' @keywords internal
#' @noRd


.preprocessPcData <- function(betaM, requiredCpGs, minCoverage, verbose) {
    targetNames <- names(requiredCpGs)
    presentCpGs <- intersect(colnames(betaM), targetNames)

    # 1. Coverage Check
    coverageRatio <- length(presentCpGs) / length(targetNames)
    if (verbose) {
        message(sprintf(
            "Probe Check: Found %d / %d required CpGs (%.1f%%).",
            length(presentCpGs), length(targetNames), coverageRatio * 100
        ))
    }

    if (coverageRatio < minCoverage) {
        if (verbose) warning("Coverage below threshold. Calculation aborted.")
        return(NULL)
    }

    # 2. Extract and Mean Impute existing CpGs
    betaSub <- betaM[, presentCpGs, drop = FALSE]
    if (any(is.na(betaSub))) {
        if (verbose) message("[PCClocks] Imputing missing values within present CpGs using column means...")
        # Vectorized column mean imputation
        colMeansV <- colMeans(betaSub, na.rm = TRUE)
        naIdx <- which(is.na(betaSub), arr.ind = TRUE)
        betaSub[naIdx] <- colMeansV[naIdx[, 2]]
    }

    # 3. Append missing CpGs using Reference Medians/Means
    missingCpGs <- setdiff(targetNames, presentCpGs)
    if (length(missingCpGs) > 0) {
        if (verbose) message(sprintf("[PCClocks] Appending %d missing CpGs with reference values...", length(missingCpGs)))
        imputeMat <- matrix(rep(requiredCpGs[missingCpGs], each = nrow(betaM)),
            nrow = nrow(betaM), ncol = length(missingCpGs)
        )
        colnames(imputeMat) <- missingCpGs
        betaSub <- cbind(betaSub, imputeMat)
    }

    return(betaSub[, targetNames, drop = FALSE])
}
