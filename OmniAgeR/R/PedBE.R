#' @title The PedBE (Pediatric Buccal) Clock for DNAm Age in Children
#'
#' @description
#' Implements the Pediatric Buccal Epigenetic (PedBE) clock, specifically developed to estimate DNA methylation age in **pediatric (childhood) samples**, as described by McEwen et al. (2020).
#'
#' @details
#' This clock is specifically trained on and designed for pediatric buccal
#' epithelial (cheek swab) samples. The calculation is a two-step process:
#'
#' 1.  A linear predictor is first calculated from the beta values using the 94-CpG elastic net coefficients . This value represents a *transformed* age.
#' 2.  A non-linear inverse age transformation (via the internal `anti.trafo` function from Horvath, 2013) is then applied. This converts the transformed age into a final estimate of chronological age in years.
#'
#' @export
#'
#' @references
#' McEwen LM, O'Donnell KJ, McGill MG, et al.
#' The PedBE clock accurately estimates DNA methylation age in pediatric buccal cells.
#' \emph{Proc Natl Acad Sci U S A.} 2020
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' PedBE.out <- PedBE(hannum_bmiq_m)


PedBE <- function(beta.m) {

  # --- Step 1: Load and parse coefficients ---
  data("PedBECoef")
  Coef_lv <- list()

  # Intercept
  Coef_lv[[1]] <- as.numeric(PedBECoef[1, 2])

  # Coefficients
  coefficients <- as.numeric(as.vector(PedBECoef[2:nrow(PedBECoef), 2]))
  names(coefficients) <- as.vector(PedBECoef[2:nrow(PedBECoef), 1])
  Coef_lv[[2]] <- coefficients


  # --- Step 2: Define the non-linear transformation function (anti.trafo) ---
  # (This is identical to the one in Horvath 2013)
  anti.trafo <- function(ptage.v, adult.age = 20) {
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
                                         coef.lv = Coef_lv,
                                         clock.name = "PedBE")

  # --- Step 4: Apply the non-linear transformation ---
  predMage.v <- anti.trafo(predTage.v)

  # --- Step 5: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predMage.v)
}





