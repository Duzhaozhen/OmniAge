#' @title Predict Gladyshev Organ-Specific Biological Age or Mortality Risk
#'
#' @description This function computes organ-specific biological ages or mortality risk scores 
#' utilizing the proteomic aging clocks developed by Goeminne et al. (2024). It maps 
#' user-provided protein expression profiles to pre-trained Elastic Net models, offering 
#' built-in variance alignment to the UK Biobank training cohort and Gompertz-based 
#' transformation of mortality hazards into chronological age equivalents.
#'
#' @param protExp A numeric matrix or data frame of protein expression values. 
#'     Rows must represent protein identifiers (e.g., Gene Symbol/UniProt) and 
#'     columns must represent individual samples..
#' @param modelType A character string indicating the prediction objective: \code{"gen1"} 
#'     (chronological age prediction) or \code{"gen2"} (mortality risk assessment). 
#'     Default is \code{"gen2"}.
#' @param platform A character string denoting the proteomic assay platform used 
#'     for model training: \code{"full"} (Olink Explore 3072) or \code{"reduced"} 
#'     (Olink Explore 1536, suitable for longitudinal or external cohorts). Default is \code{"full"}.
#' @param fold An integer (1-5) specifying the cross-validation fold to apply, 
#'     or \code{"ensemble"} to utilize the averaged coefficients across all folds 
#'     for enhanced predictive robustness. Default is 1.
#' @param toYears A logical scalar. Applicable exclusively to \code{"gen2"} models. 
#'     If \code{TRUE}, converts the predicted relative log-mortality hazard scores 
#'     into biological age units (years) under a Gompertz mortality assumption. Default is \code{FALSE}.
#' @param standardize A logical scalar. If \code{TRUE}, performs variance alignment by 
#'     standardizing the input data and rescaling it to the standard deviations of the 
#'     UK Biobank training cohort. This step is highly recommended to mitigate batch 
#'     effects unless the data has been rigorously harmonized prior to input. Default is \code{TRUE}.
#'
#'
#' @return A data frame containing the predicted scores for each organ.
#' @export
#' 
#' @references
#' Goeminne L et al.
#' Plasma protein-based organ-specific aging and mortality models unveil 
#' diseases as accelerated aging of organismal systems 
#' \emph{Cell Metabolism} 2024
#' 
#' @export
#' @examples
#' protExample <- loadOmniAgeRdata(
#'     "omniager_proteomic_olink1536_example",
#'     verbose = FALSE
#' )
#' 
#' protExp <- protExample[[1]]
#'
#' results <- gladyshevOrganAge(protExp,modelType="gen2",platform = "reduced")
#' 
gladyshevOrganAge <- function(protExp,
                              modelType = c("gen2", "gen1"),
                              platform = c("full", "reduced"),
                              fold = 1,
                              toYears = FALSE,
                              standardize = TRUE) {
  
  # 1. Argument matching
  modelType <- match.arg(modelType)
  platform <- match.arg(platform)
  protExp <- as.matrix(t(protExp))
  # 2. Extract specific model weights and standard deviations
  modelData <- loadOmniAgeRdata("omniager_proteomic_gladyshev_organ_age_model")
  
  ukbSds <- modelData[["sds"]][[platform]]
  coefList <- modelData[["coef"]][[platform]][[modelType]]
  
  # Check Gompertz parameters if conversion is requested
  if (modelType == "gen2" && toYears) {
    gompertz <- modelData[["mortality_transform"]]
    if (is.null(gompertz) || !("slope" %in% names(gompertz))) {
      stop("Gompertz parameters missing in modelData.")
    }
  }
  
  # 3. Intersect features between user data and model
  commonProts <- intersect(colnames(protExp), names(ukbSds))
  if (length(commonProts) == 0) {
    stop("No matching proteins found between protExp and model.")
  }
  
  dataMat <- as.matrix(protExp[, commonProts, drop = FALSE])
  
  # 4. Variance alignment (rescaling to UK Biobank distribution)
  if (standardize) {
    userSds <- apply(dataMat, 2, sd, na.rm = TRUE)
    userSds[userSds == 0] <- 1  # Prevent division by zero
    
    dataMat <- scale(dataMat, center = FALSE, scale = userSds)
    dataMat <- sweep(dataMat, MARGIN = 2, STATS = ukbSds[commonProts], FUN = "*")
  }
  
  # 5. Initialize result components
  organNames <- names(coefList)
  resDf <- data.frame(row.names = rownames(protExp))
  
  # 6. Iterate and compute score per organ
  for (organ in organNames) {
    dfWeights <- coefList[[organ]]
    
    if (as.character(fold) == "ensemble") {
      # Average across all folds
      foldCols <- grep("^fold_", colnames(dfWeights), value = TRUE)
      finalCoefs <- rowMeans(dfWeights[, foldCols, drop = FALSE], na.rm = TRUE)
    } else {
      # Extract specific fold
      colName <- paste0("fold_", fold)
      if (!colName %in% colnames(dfWeights)) {
        stop(sprintf("In organ '%s', column '%s' not found.", organ, colName))
      }
      finalCoefs <- dfWeights[[colName]]
    }
    
    # Apply feature names
    names(finalCoefs) <- dfWeights[["feature_name"]]
    
    # Match target proteins
    targetProts <- intersect(names(finalCoefs), commonProts)
    matchedWeights <- finalCoefs[targetProts]
    
    # Compute base score via matrix multiplication
    score <- dataMat[, targetProts, drop = FALSE] %*% matchedWeights
    
    # Adjust intercept for chronological models (gen1)
    if (modelType == "gen1" && "Intercept" %in% names(finalCoefs)) {
      score <- score + finalCoefs["Intercept"]
    }
    
    # Convert hazard to years for mortality models (gen2)
    if (modelType == "gen2" && toYears) {
      h0 <- gompertz[["avg_rel_log_mort_hazard"]]
      score <- (-h0 + score) / gompertz[["slope"]] - gompertz[["intercept"]]
    }
    
    # Assign values with original sample names
    resDf[[organ]] <- as.vector(score)
  }
  resultList <- lapply(resDf, function(column) {
    stats::setNames(column, rownames(protExp))
  })
  
  
  return(resultList)
}