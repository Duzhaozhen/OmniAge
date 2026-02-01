#' @title Calculate the Intrinsic Capacity (IC) Epigenetic Clock
#'
#' @description
#' Implements the DNA methylation-based predictor of Intrinsic Capacity (IC_Clock),
#' a biomarker of functional aging and mortality risk.
#'
#' @details
#' This function calculates the Intrinsic Capacity (IC) Clock, a novel
#' DNAm biomarker trained on clinical evaluations of physical and mental
#' capacities
#'
#' @param beta.m A numeric matrix of DNA methylation beta values.
#'   **Rows must correspond to samples and columns to CpGs.**
#'   `rownames` (sample IDs) and `colnames` (CpG probe IDs) are required.
#'   The matrix should not contain `NA` values.
#'
#' @return
#' A **numeric vector** containing the predicted Intrinsic Capacity (IC) score
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `beta.m`.
#'
#' @export
#'
#' @references
#' Fuentealba M, Rouch L, Guyonnet S, et al.
#' A blood-based epigenetic clock for intrinsic capacity predicts mortality and is associated with clinical, immunological and lifestyle factors.
#' \emph{Nature Aging.} 2025
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' IC_Clock.o <- IC_Clock(hannum_bmiq_m)

IC_Clock <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("IC_Clock_Coef")
  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(IC_Clock_Coef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(IC_Clock_Coef[2:nrow(IC_Clock_Coef), 2]))
  names(coefficients) <- as.vector(IC_Clock_Coef[2:nrow(IC_Clock_Coef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "IC_Clock")

  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}
