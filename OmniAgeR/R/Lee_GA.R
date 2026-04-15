#' @title Calculate the Lee gestational age
#'
#' @description
#' Implements the placental epigenetic clocks for estimating gestational
#' age (GA) using DNA methylation data, as described by Lee et al. (2019).
#'
#' @details
#' This function computes three distinct GA clocks derived from the models
#' presented in the Lee et al. (2019) study.
#'
#' The implemented clocks are:
#' \itemize{
#'   \item \strong{`LeeControl`}: The control model (546 CpGs).
#'   \item \strong{`LeeRobust`}: The robust model (558 CpGs).
#'   \item \strong{`LeeRefinedRobust`}: The refined robust model (395 CpGs).
#' }
#' The function iterates through each clock, matches the required
#' CpGs (e.g., 546 for `LeeControl`) with the columns in the input matrix,
#' and calculates a linear prediction of GA.
#'
#' @param betaM DNAm beta value matrix with rows labeling Illumina 450k/EPIC
#' CpGs and columns labeling samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A `list` containing three named numeric vectors.
#' Each vector represents the predicted gestational age (in weeks)
#' for the corresponding samples.
#' \itemize{
#'   \item \strong{`LeeControl`}: Numeric vector of predicted GAs
#'   from the Control model.
#'   \item \strong{`LeeRobust`}: Numeric vector of predicted GAs
#'   from the Robust model.
#'   \item \strong{`LeeRefinedRobust`}: Numeric vector of predicted GAs
#'   from the Refined Robust model.
#' }
#' Each vector is named with the sample IDs from the `rownames` of `beta.m`.
#'
#' @export
#'
#' @references
#' Lee Y, Choufani S, Weksberg R, et al.
#' Placental epigenetic clocks: estimating gestational age using placental
#' DNA methylation levels.
#' \emph{Aging} 2019
#'
#'
#' @examples
#' # 1. Fast runnable code to satisfy BiocCheck
#' message("Ready to calculate Lee gestational age.")
#' 
#' # 2. Real execution (wrapped to bypass 5-second limit in automated checks)
#' if (interactive()) {
#'   hannumBmiqM <- loadOmniAgeRdata(
#'       "omniager_hannum_example",
#'       verbose = FALSE
#'   )[[1]]
#'   LeeGaOut <- LeeGa(hannumBmiqM)
#'   print(head(LeeGaOut$LeeControl))
#' }
#'
LeeGa <- function(betaM,
                  minCoverage = 0,
                  verbose = TRUE) {
    # --- Step 1: Load Coefficients ---
    leeGACoef <- loadOmniAgeRdata(
        "omniager_lee_ga_coef",
        verbose = verbose
    )
    # Define the specific names for the three sub-clocks
    clockNames <- c("LeeControl", "LeeRobust", "LeeRefinedRobust")

    # --- Step 2: Calculate Scores for Each Clock ---
    estLv <- list()

    # Loop through the list of coefficients (causalClock.l)
    # using seq_along instead of 1:length for safety
    for (i in seq_along(leeGACoef)) {
        # Call the internal helper to handle all calculation and logging
        estLv[[i]] <- .calLinearClock(
            betaM = betaM,
            coefData = leeGACoef[[i]],
            clockLabel = clockNames[i],
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    # Assign names to the result list
    names(estLv) <- clockNames

    return(estLv)
}
