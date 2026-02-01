
#' @title Calculate the Retro-age Epigenetic Clock
#'
#' @description
#' Calculates the "Retro-age," a retroelement-based epigenetic clock for
#' chronological age, based on the models developed by Ndhlovu et al. (2024).
#' The function can compute either Version 1 (V1) or Version 2 (V2) of the clock.
#'
#' @param beta.m A matrix of DNA methylation beta values. **Rows must correspond
#' to samples** and **columns to CpG probes**. Rownames (sample IDs) and
#' colnames (CpG probe IDs) are required.
#' @param version A character string specifying which version of the clock to
#' use. Valid options are "V1" (EPICv1) or "V2" (EPICv1 & EPICv2 compatible,
#' recommended). Defaults to "V2".
#' @param verbose A logical value. If TRUE (default), the function will
#' print messages detailing the calculation steps.
#'
#' @return A numeric vector of predicted ages (Retro-age). The vector is
#' named with the sample IDs from the input matrix rownames.
#'
#' @export
#'
#' @references
#' Ndhlovu LC, Bendall ML, Dwaraka V, et al.
#' Retro-age: A unique epigenetic biomarker of aging captured by DNA methylation states of retroelements.
#' \emph{Aging Cell.} 2024
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Retro_age.o <- Retro_age(hannum_bmiq_m)

Retro_age <- function(beta.m,version="V2", verbose = TRUE) {

  res_list <- list()
  # --- Step 1: Load and parse coefficients ---
  data("Retro_age_Coef")
  temp_coef <- Retro_age_Coef[[version]]

  if (verbose) {
    print(paste("[Retro_age] Calculating Retro_age using coefficients for version:", version))
  }

  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(temp_coef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(temp_coef[2:nrow(temp_coef), 2]))
  names(coefficients) <- as.vector(temp_coef[2:nrow(temp_coef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "Retro_age",
                                        verbose)

  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}

