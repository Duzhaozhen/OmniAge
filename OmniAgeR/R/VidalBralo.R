#' @title Calculate the Vidal-Bralo Epigenetic Clock (8-CpG Model)
#'
#' @description
#' Implements the simplified 8-CpG epigenetic clock for estimating
#' chronological age in adults, as described by Vidal-Bralo et al. (2016).
#'
#' @details
#' This function calculates the highly simplified epigenetic clock developed
#' for use in adult whole blood samples. The model is a linear predictor
#' based on 8 specific CpG sites.
#'
#'
#' @param beta.m A numeric matrix of DNA methylation beta values.
#'   **Rows must correspond to samples and columns to CpGs.**
#'   `rownames` (sample IDs) and `colnames` (CpG probe IDs) are required.
#'   The matrix should not contain `NA` values.
#'
#' @return
#' A **numeric vector** containing the predicted chronological age (in years)
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `beta.m`.
#'
#' @export
#'
#' @references
#' Vidal-Bralo L, Lopez-Golan Y, Gonzalez A.
#' Simplified Assay for Epigenetic Age Estimation in Whole Blood of Adults.
#' \emph{Front Genet.} 2016
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' VidalBralo.out <- VidalBralo(hannum_bmiq_m)


VidalBralo <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("VidalBraloCoef")
  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(VidalBralo_Coef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(VidalBralo_Coef[2:nrow(VidalBralo_Coef), 2]))
  names(coefficients) <- as.vector(VidalBralo_Coef[2:nrow(VidalBralo_Coef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "VidalBralo")

  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}

