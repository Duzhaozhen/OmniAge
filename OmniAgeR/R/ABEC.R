#' @title  Adult Blood-based EPIC Clock (ABEC)
#'
#' @description
#' Predicts biological age using the Adult Blood-based EPIC Clock (ABEC). 
#' Developed by Lee et al., this model was trained on DNA methylation (DNAm) 
#' data from the Norwegian Mother, Father and Child Cohort Study (MoBa) 
#' (n = 1,592, age range: 19–59 years) using the Illumina EPIC platform.
#' 
#' @details
#' The function extracts the necessary CpG coefficients from the internal 
#' \code{ABEC_Coef} dataset and applies them to the provided beta value matrix. 
#' It relies on the \code{calculateLinearPredictor} helper function to handle 
#' missing probes and compute the final age estimates.
#'
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
#' Lee, Y., Haftorn, K.L., Denault, W.R.P. et al.
#' Blood-based epigenetic estimators of chronological age in human adults using DNA methylation data from the Illumina MethylationEPIC array. 
#' \emph{BMC Genomics} 2020
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' ABEC_o <- ABEC(hannum_bmiq_m)




ABEC <- function(beta.m) {
  
  # --- Step 1: Load and parse coefficients ---
  data("ABEC_Coef")
  Coef_lv <- list()
  
  # Intercept
  Coef_lv[[1]] <- as.numeric(ABEC_Coef[1, 2])
  
  # Coefficients
  coefficients <- as.numeric(as.vector(ABEC_Coef[2:nrow(ABEC_Coef), 2]))
  names(coefficients) <- as.vector(ABEC_Coef[2:nrow(ABEC_Coef), 1])
  Coef_lv[[2]] <- coefficients
  
  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "ABEC")
  
  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}


#' @title Extended Adult Blood-based EPIC Clock (eABEC)
#'
#' @description
#' Predicts biological age using the Extended Adult Blood-based EPIC Clock (eABEC). 
#' This model extends the training set of ABEC by incorporating public data from 
#' the Gene Expression Omnibus (GEO), resulting in a broader age-span 
#' (n = 2,227, age range: 18–88 years).
#'
#' @inheritParams ABEC
#' @inherit ABEC return
#'
#' @details
#' Similar to \code{ABEC}, this function utilizes the \code{eABEC_Coef} dataset. 
#' It is designed for applications where a wider range of adult ages is expected.
#'
#' @export
#' 
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' eABEC_o <- eABEC(hannum_bmiq_m)


eABEC <- function(beta.m) {
  
  # --- Step 1: Load and parse coefficients ---
  data("eABEC_Coef")
  Coef_lv <- list()
  
  # Intercept
  Coef_lv[[1]] <- as.numeric(eABEC_Coef[1, 2])
  
  # Coefficients
  coefficients <- as.numeric(as.vector(eABEC_Coef[2:nrow(eABEC_Coef), 2]))
  names(coefficients) <- as.vector(eABEC_Coef[2:nrow(eABEC_Coef), 1])
  Coef_lv[[2]] <- coefficients
  
  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "ABEC")
  
  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}

#' @title Common Adult Blood-based EPIC Clock (cABEC)
#'
#' @description
#' Predicts biological age using the Common Adult Blood-based EPIC Clock (cABEC). 
#' This model uses the same extended training set as \code{eABEC} but is 
#' restricted to CpGs common to both Illumina 450K and EPIC arrays, ensuring 
#' backward compatibility and robustness across platforms.
#'
#' @inheritParams ABEC
#' @inherit ABEC return
#'
#' @details
#' The function uses coefficients from the \code{cABEC_Coef} dataset.
#'
#' @export
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' cABEC_o <- cABEC(hannum_bmiq_m)

cABEC <- function(beta.m) {
  
  # --- Step 1: Load and parse coefficients ---
  data("cABEC_Coef")
  Coef_lv <- list()
  
  # Intercept
  Coef_lv[[1]] <- as.numeric(cABEC_Coef[1, 2])
  
  # Coefficients
  coefficients <- as.numeric(as.vector(cABEC_Coef[2:nrow(cABEC_Coef), 2]))
  names(coefficients) <- as.vector(cABEC_Coef[2:nrow(cABEC_Coef), 1])
  Coef_lv[[2]] <- coefficients
  
  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predage.v <- calculateLinearPredictor(beta.m,
                                        coef.lv = Coef_lv,
                                        clock.name = "cABEC_Coef")
  
  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predage.v)
}


