#' Internal Wrapper for Linear Clock Calculation
#'
#' @description
#' A wrapper function that parses a standard coefficient table (intercept in
#' row 1, probes in rows 2+) and calls the linear predictor helper.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#' @param coefData A data frame or matrix containing the coefficients.
#'   \itemize{
#'     \item Row 1: Assumed to be the Intercept.
#'     \item Col 1: Probe IDs (CpG names).
#'     \item Col 2: Coefficient values (Weights).
#'   }
#' @param clockLabel A string, the name of the clock (for logging).
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#' @return A numeric vector of predicted ages.
#' @keywords internal
#' @noRd

.calLinearClock <- function(betaM,
                            coefData,
                            clockLabel,
                            minCoverage = 0.5,
                            verbose = TRUE) {
    # Initialize the coefficient list structure expected by .calculateLinearPredictor
    coefList <- list()

    # --- Step 1: Extract Intercept ---
    # Assumption: The intercept is located in the first row, second column
    coefList[[1]] <- as.numeric(coefData[1, 2])

    # --- Step 2: Extract Coefficients ---
    # Assumption: Rows 2 to N contain the probe weights
    idx <- 2:nrow(coefData)

    # Extract weights (Col 2) and assign probe names (Col 1)
    coefVals <- as.numeric(as.vector(coefData[idx, 2]))
    names(coefVals) <- as.vector(coefData[idx, 1])

    coefList[[2]] <- coefVals

    # --- Step 3: Compute Linear Predictor ---
    # Pass the parsed coefficients to the calculation engine
    predAge <- .calculateLinearPredictor(
        betaM,
        coefLv = coefList,
        clockName = clockLabel,
        minCoverage = minCoverage,
        verbose = verbose
    )

    return(predAge)
}


#' Check CpG coverage and data integrity for epigenetic clocks
#'
#' @description
#' An internal helper to validate the overlap between input data and
#' clock requirements.
#' It performs three primary tasks:
#' 1. Matches CpG probe IDs from the coefficients to the input matrix.
#' 2. Evaluates whether the coverage ratio meets the specified threshold.
#' 3. Enforces a strict "no-NA" policy for all required CpGs to
#' ensure calculation stability.
#'
#' @param betaM A numeric matrix of methylation beta values
#' (probes as rows, samples as columns).
#' @param allWeights A named numeric vector of CpG weights from the clock's
#' coefficient table.
#' @param clockName Character. The label of the clock for logging and
#' error reporting.
#' @param minCoverage Numeric (0-1). The minimum required proportion of CpGs
#' present in \code{betaM}.
#' @param verbose Logical. If \code{TRUE}, prints coverage
#' statistics to the console.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{pass}: Logical. \code{TRUE} if coverage meets the threshold,
#'   \code{FALSE} otherwise.
#'   \item \code{betaIdx}: Integer vector. The row indices in
#'   \code{betaM} corresponding to matched CpGs.
#'   \item \code{weightsSubset}: Numeric vector.
#'   The coefficients for the matched CpGs.
#' }
#'
#' @details
#' This function will trigger a fatal error (\code{stop}) if any \code{NA}
#' values are detected within the subset of \code{betaM} required for the clock.
#' This ensures that the linear predictor (\code{crossprod})
#' does not return \code{NA} results.
#'
#' @note
#' This is an internal function.
#'
#' @noRd
#'
.checkCpGCoverage <- function(betaM, allWeights, clockName, minCoverage, verbose) {
    requiredCpGs <- names(allWeights)
    mapIdx <- match(requiredCpGs, rownames(betaM))
    repIdx <- which(!is.na(mapIdx))

    nRequired <- length(requiredCpGs)
    nRepresented <- length(repIdx)
    currentCoverage <- if (nRequired > 0) nRepresented / nRequired else 0

    # 1. Print statistical information
    if (verbose) {
        percStr <- sprintf("%.1f%%", currentCoverage * 100)
        message(sprintf(
            "[%s] Found: %d, Required: %d (Coverage: %s)",
            clockName, nRepresented, nRequired, percStr
        ))
    }

    # 2. Coverage threshold check
    pass <- currentCoverage >= minCoverage
    if (!pass) {
        if (verbose) {
            warning(sprintf(
                "[%s] Low coverage (%.1f%%). Returning NA.",
                clockName, currentCoverage * 100
            ))
        }
        return(list(pass = FALSE))
    }

    # 3. Strict NA examination (only for the matched CpGs)
    betaIdx <- mapIdx[repIdx]

    if (anyNA(betaM[betaIdx, , drop = FALSE])) {
        stop(sprintf("[%s] NA values found in beta matrix.
                 All required CpGs must have complete data.", clockName),
            call. = FALSE
        )
    }

    return(list(
        pass = TRUE,
        betaIdx = betaIdx,
        weightsSubset = allWeights[repIdx]
    ))
}


#' @title Calculate a linear predictor for an epigenetic clock (Internal)
#' @description
#' A generic helper function to calculate the weighted linear sum
#' (predictor) for any epigenetic clock, given a beta matrix and a
#' coefficient list.
#'
#' @param betaM A numeric matrix of beta values (CpGs as rows, Samples as cols).
#' @param coefLv A list containing the clock coefficients.
#'   - `coefLv[[1]]` MUST be the numeric intercept.
#'   - `coefLv[[2]]` MUST be a named numeric vector of CpG weights.
#' @param clockName A string, the name of the clock for status messages.
#' @param minCoverage A numeric value (0-1). Threshold for CpG coverage.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages
#'
#' @return A named numeric vector of predicted ages/values for each sample.
#'
#' @keywords internal
#' @noRd
.calculateLinearPredictor <- function(betaM, coefLv, clockName = "Unknown",
                                      minCoverage = 0.5, verbose = TRUE) {
    intercept <- coefLv[[1]]
    allWeights <- coefLv[[2]]

    # Call the verification function
    coverage <- .checkCpGCoverage(betaM, allWeights, clockName, minCoverage, verbose)

    if (!coverage$pass) {
        naVec <- rep(NA_real_, ncol(betaM))
        names(naVec) <- colnames(betaM)
        return(naVec)
    }

    # Perform matrix calculation
    predAge <- as.vector(intercept +
        crossprod(
            betaM[coverage$betaIdx, , drop = FALSE],
            coverage$weightsSubset
        ))
    names(predAge) <- colnames(betaM)
    return(predAge)
}
