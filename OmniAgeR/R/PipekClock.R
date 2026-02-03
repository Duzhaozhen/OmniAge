#' @title Pipek's Multi-tissue Elastic Net Epigenetic Clock (239 CpGs)
#'
#' @description
#' Implements the "elasticNet (239)" epigenetic clock model proposed by Pipek and Csabai (2023).
#' This model was trained using elastic net regression on a large multi-tissue dataset containing
#' methylation data from Illumina 27K, 450K, and EPIC platforms. It is designed to provide
#' improved accuracy on EPIC array data compared to the original Horvath2013 clock.
#'
#' @details
#' Implements the "elasticNet (239)" model using **239 CpGs** shared across Illumina 27K, 450K, and EPIC arrays.
#' 
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed prior to input.
#'   \item **Normalization:** Recommended (e.g., BMIQ) but not strictly required.
#' }
#' 
#' Internally applies Horvath's (2013) log-linear age transformation.
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
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for methylation array data. 
#' \emph{J Math Chem} 2023
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' PipekElasticNet_o <- PipekElasticNet(hannum_bmiq_m)
#' 



PipekElasticNet <- function(beta.m) {
  
  # --- Step 1: Load and parse coefficients ---
  data("PipekElasticNet_Coef")
  Coef_lv <- list()
  
  # Intercept
  Coef_lv[[1]] <- as.numeric(PipekElasticNet_Coef[1, 2])
  
  # Coefficients 
  coefficients <- as.numeric(as.vector(PipekElasticNet_Coef[2:nrow(PipekElasticNet_Coef), 2]))
  names(coefficients) <- as.vector(PipekElasticNet_Coef[2:nrow(PipekElasticNet_Coef), 1])
  Coef_lv[[2]] <- coefficients
  
  
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
  
  transform_prediction <- function(pred, age_adult=20){
    ifelse(pred<0, (1+age_adult)*exp(pred)-1, (1+age_adult)*pred+age_adult)
  }
  # --- Step 3: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)
  predTage.v <- calculateLinearPredictor(beta.m,
                                         coef.lv = Coef_lv,
                                         clock.name = "PipekElasticNet")
  
  # --- Step 4: Apply the non-linear transformation ---
  predMage.v <- iageF(predTage.v)
  
  # --- Step 5: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predMage.v)
}




#' @title Pipek's Filtered Horvath Epigenetic Clock (272 CpGs)
#'
#' @description
#' Implements the "filtered H (272)" epigenetic clock model proposed by Pipek and Csabai (2023).
#' This model restricts feature selection to the subset of CpG sites from the original Horvath
#' pan-tissue clock that are also present on the EPIC array.
#'
#' @details
#' Implements the "filtered H (272)" model. Unlike the full ElasticNet model, this clock was trained 
#' by limiting candidate features to the intersection of original Horvath probes and the 
#' cross-platform (27K/450K/EPIC) probe set.
#' 
#' **Key Features:**
#' \itemize{
#'   \item **272 CpGs:** A subset of the original Horvath clock, re-optimized for better EPIC array compatibility.
#'   \item **Best for:** Datasets pre-filtered to Horvath clock probes but requiring updated calibration.
#' }
#' 
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed.
#'   \item **Normalization:** Recommended.
#' }
#'
#' Internally applies Horvath's (2013) log-linear age transformation.
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
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for methylation array data. 
#' \emph{J Math Chem} 2023
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' PipekFilteredh_o <- PipekFilteredh(hannum_bmiq_m)
#' 


PipekFilteredh <- function(beta.m) {
  # --- Step 1: Load and parse coefficients ---
  data("PipekFilteredh_Coef")
  Coef_lv <- list()
  
  # Intercept
  Coef_lv[[1]] <- as.numeric(PipekFilteredh_Coef[1, 2])
  
  # Coefficients 
  coefficients <- as.numeric(as.vector(PipekFilteredh_Coef[2:nrow(PipekFilteredh_Coef), 2]))
  names(coefficients) <- as.vector(PipekFilteredh_Coef[2:nrow(PipekFilteredh_Coef), 1])
  Coef_lv[[2]] <- coefficients
  
  
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
                                         coef.lv = Coef_lv,
                                         clock.name = "PipekFilteredh")
  
  # --- Step 4: Apply the non-linear transformation ---
  predMage.v <- iageF(predTage.v)
  
  # --- Step 5: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predMage.v)
}


#' @title Pipek's Retrained Horvath Epigenetic Clock (308 CpGs)
#'
#' @description
#' Implements the "retrained H (308)" epigenetic clock model proposed by Pipek and Csabai (2023).
#' This model is a direct recalibration of the original Horvath clock probes using a large 
#' multi-platform training set.
#'
#' @details
#' Implements the "retrained H (308)" model. It includes all **308 CpG sites** from the original 
#' Horvath pan-tissue clock that are present across 27K, 450K, and EPIC platforms.
#' 
#' **Key Features:**
#' \itemize{
#'   \item **No Feature Selection:** Coefficients were simply refitted to the new data without dropping probes.
#'   \item **Robust Update:** Serves as a direct update to the Horvath clock to correct for accuracy loss on EPIC arrays.
#' }
#' 
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed.
#'   \item **Normalization:** Recommended.
#' }
#'
#' Internally applies Horvath's (2013) log-linear age transformation.
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
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for methylation array data. 
#' \emph{J Math Chem} 2023
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' PipekRetrainedh_o <- PipekRetrainedh(hannum_bmiq_m)
#' 



PipekRetrainedh <- function(beta.m) {
  
  # --- Step 1: Load and parse coefficients ---
  data("PipekRetrainedh_Coef")
  Coef_lv <- list()
  
  # Intercept
  Coef_lv[[1]] <- as.numeric(PipekRetrainedh_Coef[1, 2])
  
  # Coefficients 
  coefficients <- as.numeric(as.vector(PipekRetrainedh_Coef[2:nrow(PipekRetrainedh_Coef), 2]))
  names(coefficients) <- as.vector(PipekRetrainedh_Coef[2:nrow(PipekRetrainedh_Coef), 1])
  Coef_lv[[2]] <- coefficients
  
  
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
                                         coef.lv = Coef_lv,
                                         clock.name = "PipekRetrainedh")
  
  # --- Step 4: Apply the non-linear transformation ---
  predMage.v <- iageF(predTage.v)
  
  # --- Step 5: Return final age vector ---
  # (Names are already attached by the helper function)
  return(predMage.v)
}