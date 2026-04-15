#' @title Calculate Epigenetic Surrogates for Complex Traits (McCartney et al. 2018)
#'
#' @description
#' Implements the elastic net (EN) predictors for 10 complex traits and
#' biomarkers using DNAm data, as described by McCartney et al. (2018).
#'
#' @details
#' This function computes epigenetic surrogate scores for 10 different traits
#' based on models trained on blood methylation data.
#'
#' The function iterates through and calculates scores for the following traits:
#' \itemize{
#'   \item `BMI` (Body Mass Index)
#'   \item `Smoking` (Smoking pack-years)
#'   \item `Alcohol` (Alcohol consumption units per week)
#'   \item `Education` (Years of education)
#'   \item `Total_cholesterol`
#'   \item `HDL_cholesterol`
#'   \item `LDL_cholesterol`
#'   \item `Total_HDL_ratio` (Ratio of Total to HDL cholesterol)
#'   \item `WHR` (Waist-to-Hip Ratio)
#'   \item `Body_fat_Perc` (Body Fat Percentage)
#' }
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A `list` containing 10 named elements, one for each trait. Each element is
#' a **named numeric vector** of the calculated epigenetic scores.
#' \itemize{
#'   \item \strong{`BMI`}: Numeric vector of predicted scores.
#'   \item \strong{`Smoking`}: Numeric vector of predicted scores.
#'   \item \strong{`Alcohol`}: Numeric vector of predicted scores.
#'   \item \strong{`Education`}: Numeric vector of predicted scores.
#'   \item \strong{`Total_cholesterol`}: Numeric vector of predicted scores.
#'   \item \strong{`HDL_cholesterol`}: Numeric vector of predicted scores.
#'   \item \strong{`LDL_cholesterol`}: Numeric vector of predicted scores.
#'   \item \strong{`Total_HDL_ratio`}: Numeric vector of predicted scores.
#'   \item \strong{`WHR`}: Numeric vector of predicted scores.
#'   \item \strong{`Body_fat_Perc`}: Numeric vector of predicted scores.
#' }
#' Each vector is named with the sample IDs from the `rownames` of `beta.m`.
#'
#' @export
#'
#' @references
#' McCartney DL, Hillary RF, Stevenson AJ, et al.
#' Epigenetic prediction of complex traits and death.
#' \emph{Genome Biol.} 2018
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' mcCartneyTraitOut <- mcCartneyTrait(hannumBmiqM)
#'
mcCartneyTrait <- function(betaM,
                           minCoverage = 0,
                           verbose = TRUE) {
    # --- Step 1: Load Coefficients ---
    mcCartneyTraitCoef <- loadOmniAgeRdata(
        "omniager_mccartney_trait_coef",
        verbose = verbose
    )

    # Define the specific names
    clockNames <- names(mcCartneyTraitCoef)

    # --- Step 2: Calculate Scores for Each Clock ---
    estLv <- list()

    # Loop through the list of coefficients
    # using seq_along instead of 1:length for safety
    for (i in seq_along(mcCartneyTraitCoef)) {
        tmpCoef <- mcCartneyTraitCoef[[i]]
        colnames(tmpCoef) <- c("probe", "coef")
        tmpCoef <- rbind(data.frame(probe = "(Intercept)", coef = 0), tmpCoef)

        # Call the internal helper to handle all calculation and logging
        estLv[[i]] <- .calLinearClock(
            betaM = betaM,
            coefData = tmpCoef,
            clockLabel = clockNames[i],
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    # Assign names to the result list
    names(estLv) <- clockNames

    return(estLv)
}
