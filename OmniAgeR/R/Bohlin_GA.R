#' @title Calculate the Bohlin Gestational Age (Cord Blood)
#'
#' @description
#' Implements the epigenetic clock for predicting gestational age (GA) using
#' newborn cord blood, as described by Bohlin et al. (2016).
#'
#' @param beta.m a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#'
#' @return #' A **numeric vector** containing the predicted gestational age (in weeks)
#' for each sample. The vector is named with the sample IDs from the `rownames`
#' of `beta.m`.
#'
#' @export
#'
#' @references
#' Bohlin J, Håberg SE, Magnus P, et al.
#' Prediction of gestational age based on genome-wide differentially methylated regions.
#' \emph{Genome Biol.} 2016
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Bohlin_GA.out <- Bohlin_GA(hannum_bmiq_m)

Bohlin_GA <- function(beta.m) {

  res_list <- list()
  # --- Step 1: Load and parse coefficients ---
  data("BohlinGACoef")
  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(BohlinGACoef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(BohlinGACoef[2:nrow(BohlinGACoef), 2]))
  names(coefficients) <- as.vector(BohlinGACoef[2:nrow(BohlinGACoef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "Bohlin_GA")
  predage.v <- predage.v/7  ## Convert the number of days into weeks
  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}
