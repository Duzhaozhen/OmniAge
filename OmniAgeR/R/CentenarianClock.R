#' @title Calculate Centenarian Epigenetic Clocks (Eric Dec et al.)
#'
#' @description  Calculates the Centenarian epigenetic clocks
#' (ENCen40 and ENCen100) developed by Eric Dec et al. (2023). This function
#' serves as a wrapper that loads the internal clock coefficients and computes
#' the linear predictors for each clock using the helper function.
#'
#' @param betaM a matrix of methylation beta values.
#' Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A list containing the predicted scores for each Centenarian clock.
#' \itemize{
#'   \item `ENCen40`: Elastic net clock trained on individuals aged 40+.
#'   \item `ENCen100`: Elastic net clock specifically trained on centenarians (100+).
#' }
#'
#' @export
#'
#' @references
#' Dec, E., Clement, J., Cheng, K. et al.
#' Centenarian clocks: epigenetic clocks for validating claims of
#' exceptional longevity. \emph{GeroScience} 2023
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' centenarianClockOut <- centenarianClock(hannumBmiqM)
#'
centenarianClock <- function(betaM,
                             minCoverage = 0,
                             verbose = TRUE) {
    # --- Step 1: Load Coefficients ---
    CentenarianENCoef <- loadOmniAgeRdata(
        "omniager_centenarian_coef",
        verbose = verbose
    )
    # Define the specific names for the three sub-clocks
    clockNames <- c("ENCen40", "ENCen100")

    # --- Step 2: Calculate Scores for Each Clock ---
    estLv <- list()

    # Loop through the list of coefficients
    for (i in seq_along(CentenarianENCoef)) {
        # Call the internal helper to handle all calculation and logging
        estLv[[i]] <- .calLinearClock(
            betaM = betaM,
            coefData = CentenarianENCoef[[i]],
            clockLabel = clockNames[i],
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    # Assign names to the result list
    names(estLv) <- clockNames

    return(estLv)
}
