#' @title Calculate the Lee gestational age
#'
#' @description
#' Implements the placental epigenetic clocks for estimating gestational age (GA) using DNA methylation data, as described by Lee et al. (2019).
#'
#' @details
#' This function computes three distinct GA clocks derived from the models presented in the Lee et al. (2019) study.
#'
#' The implemented clocks are:
#' \itemize{
#'   \item \strong{`LeeControl`}: The control model (546 CpGs).
#'   \item \strong{`LeeRobust`}: The robust model (558 CpGs).
#'   \item \strong{`LeeRefinedRobust`}: The refined robust model (395 CpGs).
#' }
#' The function iterates through each clock, matches the required CpGs (e.g., 546 for `LeeControl`) with the columns in the input `beta.m` matrix, and calculates a linear prediction of GA.
#'
#' @param beta.m a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#'
#' @return
#' A `list` containing three named numeric vectors. Each vector represents the predicted gestational age (in weeks) for the corresponding samples.
#' \itemize{
#'   \item \strong{`LeeControl`}: Numeric vector of predicted GAs from the Control model.
#'   \item \strong{`LeeRobust`}: Numeric vector of predicted GAs from the Robust model.
#'   \item \strong{`LeeRefinedRobust`}: Numeric vector of predicted GAs from the Refined Robust model.
#' }
#' Each vector is named with the sample IDs from the `rownames` of `beta.m`.
#'
#' @export
#'
#' @references
#' Lee Y, Choufani S, Weksberg R, et al.
#' Placental epigenetic clocks: estimating gestational age using placental DNA methylation levels.
#' \emph{Aging} 2019
#'
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Lee_GA.out <- Lee_GA(hannum_bmiq_m)



Lee_GA <- function(beta.m) {

  res_list <- list()
  # --- Step 1: Load and parse coefficients ---
  data("LeeGACoef")
  for (i in seq_along(LeeGACoef)) {
    tmp_Coef  <- LeeGACoef[[i]]

    Coef_lv <- list()

    # Intercept
    Coef_lv[[1]] <- as.numeric(tmp_Coef[1, 2])

    # Coefficients
    coefficients <- as.numeric(as.vector(tmp_Coef[2:nrow(tmp_Coef), 2]))
    names(coefficients) <- as.vector(tmp_Coef[2:nrow(tmp_Coef), 1])
    Coef_lv[[2]] <- coefficients

    # --- Step 2: Calculate the linear predictor ---
    # (Requires the 'calculateLinearPredictor' function)
    predage.v <- calculateLinearPredictor(beta.m,
                                          coef.lv = Coef_lv,
                                          clock.name = names(LeeGACoef)[i])

    res_list[[names(LeeGACoef)[i]]] <- predage.v
  }


  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(res_list)
}
