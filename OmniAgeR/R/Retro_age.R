#' @title Calculate the Retro-age Epigenetic Clock
#'
#' @description
#' Calculates the "Retro-age," a retroelement-based epigenetic clock for
#' chronological age, based on the models developed by Ndhlovu et al. (2024).
#' This function computes both Version 1 (V1) and Version 2 (V2) of the clock.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#'
#' @param minCoverage A numeric value between 0 and 1 (default is 0).
#' Specifies the minimum proportion of required CpGs that must be present
#' in the input matrix for the clock calculation to proceed.
#' @param verbose A logical value. If TRUE (default), the function will
#' print messages detailing the calculation steps.
#'
#' @return A list containing the predicted age for the "V1" and "V2" clocks.
#'
#' @export
#'
#' @references
#' Ndhlovu LC, Bendall ML, Dwaraka V, et al.
#' Retro-age: A unique epigenetic biomarker of aging captured by DNA methylation states of retroelements.
#' \emph{Aging Cell.} 2024
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' retroAgeRes <- retroAge(hannumBmiqM)
retroAge <- function(betaM,
                     minCoverage = 0,
                     verbose = TRUE) {
    # --- Step 1: Load Coefficients ---
    retroAgeCoef <- loadOmniAgeRdata(
        "omniager_retroage_coef",
        verbose = verbose
    )

    # Define the specific names for the these sub-clocks
    clockNames <- c("retroAgeV1", "retroAgeV2")

    # --- Step 2: Calculate Scores for Each Clock ---
    estLv <- list()

    # Loop through the list of coefficients
    for (i in seq_along(retroAgeCoef)) {
        # Call the internal helper to handle all calculation and logging
        estLv[[i]] <- .calLinearClock(
            betaM = betaM,
            coefData = retroAgeCoef[[i]],
            clockLabel = clockNames[i],
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    # Assign names to the result list
    names(estLv) <- clockNames

    return(estLv)
}
