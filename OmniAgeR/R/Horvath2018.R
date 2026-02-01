#' @title Calculate Horvath's Skin & Blood DNAm Age (2018)
#'
#' @description
#' A function to calculate the Horvath "Skin & Blood" clock age (2018)
#' from a DNA methylation beta value matrix.
#'
#' @param beta.m A numeric matrix of beta values. Rows should be CpG probes
#' and columns should be individual samples.
#'
#' @details
#' Implements the Horvath (2018) skin & blood clock. The function calculates
#' a weighted linear predictor from 391 CpGs and then transforms this value
#' using the same non-linear function as the Horvath (2013) clock to
#' return the final DNAm age.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding
#' to the sample IDs from the input matrix's column names.
#'
#' @export
#'
#' @references
#' Horvath S, Oshima J, Martin GM, et al.
#' Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies.
#' \emph{Aging} 2018
#'
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Horvath2018.out <- Horvath2018(hannum_bmiq_m)
#'


Horvath2018 <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("Horvath_SkinBlood_CpG")
  skinP.lv <- list()

  # Intercept
  skinP.lv[[1]] <- as.numeric(Horvath2_coef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(Horvath2_coef[2:nrow(Horvath2_coef), 2]))
  names(coefficients) <- as.vector(Horvath2_coef[2:nrow(Horvath2_coef), 1])
  skinP.lv[[2]] <- coefficients


  # --- Step 2: Define the non-linear transformation function (anti.trafo) ---
  # (This is identical to the one in Horvath 2013)
  iageF <- function(ptage.v, adult.age = 20) {
    y.idx <- which(ptage.v <= 0)
    a.idx <- which(ptage.v > 0)
    mage.v <- ptage.v
    mage.v[a.idx] <- ptage.v[a.idx] * (adult.age + 1) + adult.age
    mage.v[y.idx] <- exp(ptage.v[y.idx] + log(adult.age + 1)) - 1
    return(mage.v)
  }


  # --- Step 3: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predTage.v <- calculateLinearPredictor(beta.m,
                                         coef.lv = skinP.lv,
                                         clock.name = "Horvath2018")

  # --- Step 4: Apply the non-linear transformation ---
  predMage.v <- iageF(predTage.v)

  # --- Step 5: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predMage.v)
}



