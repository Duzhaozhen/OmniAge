#' @title Calculate EnsembleAge Epigenetic Clocks
#'
#' @description
#' Calculates epigenetic age using the EnsembleAge multi-clock framework
#' (Haghani et al., 2025). This function computes one of the three available
#' clock versions (Static, Dynamic, or HumanMouse).
#'
#' @param betaM A numeric matrix (Rows: CpGs, Cols: Samples) or \code{SummarizedExperiment}.
#' @param clockVersion Character. One of \code{"Dynamic"}, \code{"Static"}, or \code{"HumanMouse"}.
#' @param minCoverage Numeric (0-1). Minimum required proportion of CpGs present. Default is 0.
#' @param verbose Logical. Whether to print status messages.
#'
#' @return
#' A `list` where each element is a named numeric vector of predicted values.
#' The names of the list elements correspond to the specific sub-clocks
#' calculated (e.g., `HumanMouse_HumanMouse`, `Static_Static`, etc.).
#'
#' @note
#' **Understanding the `clock_version` output:**
#' * **`"Static"` and `"HumanMouse"`:** The returned values are the final,
#'     calibrated age predictors. The `"HumanMouse"` clock returns age
#'     normalized by maximum species lifespan.
#' * **`"Dynamic"`:** This option returns a list of predictors from all
#'     40 clocks in the 'Dynamic' ensemble library. To complete the full
#'     "EnsembleAge.Dynamic" *methodology* as described in the paper,
#'     the user must perform **additional steps**:
#'     1.  Have a dataset with control and intervention groups.
#'     2.  Use the returned 40 predictors to find "responsive" clocks
#'         (e.g., those with |Z-score| > 2 between groups).
#'     3.  Calculate the final `EnsembleAge.Dynamic` value by taking the
#'         **median** of only those responsive clocks.
#'     This function provides the necessary *basis* for this calculation.
#'
#' @references
#' Haghani, A., Lu, A.T., Yan, Q. et al.
#' EnsembleAge: enhancing epigenetic age assessment with a
#' multi-clock framework.
#' \emph{GeroScience} 2025.
#'
#' @export
#' @importFrom utils data
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' # Ensure it's CpGs=rows, Samples=cols
#' # Calculate the HumanMouse clock version
#' ensembleAgeOut <- ensembleAge(hannumBmiqM, clockVersion = "HumanMouse")
#'
ensembleAge <- function(betaM, clockVersion = c("HumanMouse", "Static", "Dynamic"),
                        minCoverage = 0, verbose = TRUE) {
    clockVersion <- match.arg(clockVersion)

    # --- 1. Data Loading --
    ensembleAgeCoef <- loadOmniAgeRdata(
        "omniager_ensembleage_coef",
        verbose = verbose
    )

    if (verbose) {
        message("[EnsembleAge] Initializing EnsembleAge_", clockVersion, " calculation...")
    }

    # --- 2. Validation & SE Support ---
    if (!clockVersion %in% names(ensembleAgeCoef)) {
        stop("[EnsembleAge] Loaded coefficients do not contain version: ", clockVersion)
    }

    coefList <- ensembleAgeCoef[[clockVersion]]
    resList <- list()

    # --- 3. Iterate through sub-clocks in the Ensemble ---
    for (i in seq_along(coefList)) {
        clockSubName <- names(coefList)[i]
        coefData <- coefList[[i]]

        fullLabel <- paste0(clockVersion, "_", clockSubName)

        if (is.null(coefData) || nrow(coefData) < 2) {
            if (verbose) message(sprintf("[%s] Skipping: Invalid coefficient table.", fullLabel))
            next
        }

        # --- 4. Call Internal Wrappers ---
        resList[[fullLabel]] <- .calLinearClock(
            betaM = betaM,
            coefData = coefData,
            clockLabel = fullLabel,
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    return(resList)
}
