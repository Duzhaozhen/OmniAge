#' @title Calculate the Mayne Placental Gestational Age.
#'
#' @description
#' Implements the 62-CpG placental epigenetic clock for estimating gestational age (GA), as described by Mayne et al. (2017).
#'
#' @param beta.m A numeric matrix of DNA methylation beta values.
#'   **Rows must correspond to samples and columns to CpGs.**
#'   `rownames` (sample IDs) and `colnames` (CpG probe IDs) are required.
#'   The matrix should not contain `NA` values.
#'
#' @return
#' A **numeric vector** containing the predicted gestational age (in weeks) for each sample. The vector is named with the sample IDs from the `rownames` of `beta.m`.
#'
#' @export
#'
#' @references
#' Mayne BT, Leemaqz SY, Smith AK, Breen J, Roberts CT, Bianco-Miotto T.
#' Accelerated placental aging in early onset preeclampsia pregnancies identified by DNA methylation.
#' \emph{Epigenomics} 2017
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Mayne_GA.out <- Mayne_GA(hannum_bmiq_m)

Mayne_GA <- function(beta.m) {

  res_list <- list()
  # --- Step 1: Load and parse coefficients ---
  data("MayneGACoef")
  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(MayneGACoef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(MayneGACoef[2:nrow(MayneGACoef), 2]))
  names(coefficients) <- as.vector(MayneGACoef[2:nrow(MayneGACoef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "Mayne_GA")

  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}
