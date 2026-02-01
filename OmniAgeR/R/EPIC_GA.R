#' @title The Gestational Age (GA) clock based on 176 Illumina EPIC CpGs
#'
#' @param beta.m A matrix of beta values (CpGs in rows, samples in columns).  This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#'
#' @return A named vector of predicted Gestational Ages (in weeks).
#'
#' @export
#'
#' @references
#' Haftorn KL, Lee Y, Denault WRP, et al.
#' An EPIC predictor of gestational age and its application to newborns conceived by assisted reproductive technologies.
#' \emph{Clin Epigenetics.} 2021
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' EPIC_GA.out <- EPIC_GA(hannum_bmiq_m)

EPIC_GA <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("EPICGACoef")

  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(EPICGACoef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(EPICGACoef[2:nrow(EPICGACoef), 2]))
  names(coefficients) <- as.vector(EPICGACoef[2:nrow(EPICGACoef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "EPIC_GA")
  predage.v <- predage.v/7  ## Convert the number of days into weeks
  # --- Step 5: Return final age vector ---
  return(predage.v)
}
