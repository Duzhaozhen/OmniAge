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
#' @param beta.m A numeric matrix of DNA methylation beta values.
#'   **Rows must correspond to samples and columns to CpGs.**
#'   `rownames` (sample IDs) and `colnames` (CpG probe IDs) are required.
#'   The matrix should not contain `NA` values.
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
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' McCartney_Trait.out <- McCartney_Trait(hannum_bmiq_m)



McCartney_Trait <- function(beta.m) {

  res_list <- list()
  # --- Step 1: Load and parse coefficients ---
  data("McCartneyTraitCoef")
  for (i in seq_along(McCartneyTraitCoef)) {
    tmp_Coef  <- McCartneyTraitCoef[[i]]

    Coef_lv <- list()

    # Intercept
    Coef_lv[[1]] <- 0

    # Coefficients
    coefficients <- as.numeric(as.vector(tmp_Coef[1:nrow(tmp_Coef), 2]))
    names(coefficients) <- as.vector(tmp_Coef[1:nrow(tmp_Coef), 1])
    Coef_lv[[2]] <- coefficients

    # --- Step 2: Calculate the linear predictor ---
    # (Requires the 'calculateLinearPredictor' function)
    predage.v <- calculateLinearPredictor(beta.m,
                                          coef.lv = Coef_lv,
                                          clock.name = names(McCartneyTraitCoef)[i])

    res_list[[names(McCartneyTraitCoef)[i]]] <- predage.v
  }


  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(res_list)
}
