#' @title Wu's Epigenetic Clock for Pediatric Age Estimation
#'
#' @description
#' Estimates biological age (in years) based on DNA methylation data using the method 
#' proposed by Wu et al. (2019). This clock is specifically designed for pediatric 
#' cohorts and utilizes a non-linear transformation to capture rapid developmental 
#' changes in early life.
#'
#' @details
#' The calculation involves two main steps:
#' \enumerate{
#'   \item Calculation of a linear predictor using weighted CpG beta values.
#'   \item Transformation of the linear predictor into biological age (years) using 
#'   a specific "anti-transformation" function with a toddler age offset of 48 months.
#' }
#' 
#' \strong{Data Requirements:}
#' The input \code{beta.m} must be a matrix of Beta values (0 to 1). The function 
#' expects CpG probes as row names.
#'
#' @param beta.m A numeric matrix of DNA methylation beta values.
#'   **Rows must correspond to samples and columns to CpGs.**
#'   `rownames` (sample IDs) and `colnames` (CpG probe IDs) are required.
#'   The matrix should not contain `NA` values.
#'
#' @return A numeric vector of predicted biological ages. The vector is 
#' named using the sample IDs from the \code{rownames} of \code{beta.m}.
#' 
#' @export
#'
#' @references
#' Wu, Xiaohui et al. 
#' DNA methylation profile is a quantitative measure of biological aging in children
#' \emph{Aging} 2019
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' WuClock_o <- WuClock(hannum_bmiq_m)
#' 





WuClock <- function(beta.m) {
  
  # --- Step 1: Load and parse coefficients ---
  data("WuClockCoef")
  Coef_lv <- list()
  
  # Intercept
  Coef_lv[[1]] <- as.numeric(WuClockCoef[1, 2])
  
  # Coefficients 
  coefficients <- as.numeric(as.vector(WuClockCoef[2:nrow(WuClockCoef), 2]))
  names(coefficients) <- as.vector(WuClockCoef[2:nrow(WuClockCoef), 1])
  Coef_lv[[2]] <- coefficients
  
  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predTage.v <- calculateLinearPredictor(beta.m,
                                         coef.lv = Coef_lv,
                                         clock.name = "WuClock")
  
  # --- Step 3: Define the non-linear transformation function (anti.trafo) ---
  # This is months
  iageF <- function(ptage.v, toddler.age = 48) {
    y.idx <- which(ptage.v <= 0)
    a.idx <- which(ptage.v > 0)
    mage.v <- ptage.v
    mage.v[a.idx] <- ptage.v[a.idx] * (toddler.age + 1) + toddler.age
    mage.v[y.idx] <- exp(ptage.v[y.idx] + log(toddler.age + 1)) - 1
    return(mage.v)
  }
  predage.v <- iageF(predTage.v)
  predage.v <- predage.v/12 ## transform to years
  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}