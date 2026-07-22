#' @title Calculate the Retro-age Epigenetic Clock
#'
#' @description
#' Calculates the "Retro-age," a retroelement-based epigenetic clock for
#' chronological age, based on the models developed by Ndhlovu et al. (2024).
#' This function computes both Version 1 (V1) and Version 2 (V2) of the clock.
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
                     minCoverage = 0.5,
                     verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
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
