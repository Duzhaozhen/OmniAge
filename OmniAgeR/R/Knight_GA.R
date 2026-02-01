#' @title Calculate the Knight gestational age
#'
#' @description A function to calculate the Knight gestational age
#'
#' @param beta.m a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#'
#' @return A named vector of predicted Gestational Ages (in weeks).
#'
#' @export
#'
#' @references
#' Knight AK, Craig JM, Theda C, et al.
#' An epigenetic clock for gestational age at birth based on blood methylation data.
#' \emph{Genome Biol.} 2016
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Knight_GA.out <- Knight_GA(hannum_bmiq_m)


Knight_GA <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("KnightCoef") #

  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(KnightCoef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(KnightCoef[2:nrow(KnightCoef), 2]))
  names(coefficients) <- as.vector(KnightCoef[2:nrow(KnightCoef), 1])
  Coef_lv[[2]] <- coefficients

  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "Knight_GA")

  # --- Step 5: Return final age vector ---
  return(predage.v)
}

