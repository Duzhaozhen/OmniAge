#' @title Predict age using SVZ Bootstrap Clocks (Supports Multiple Cell Types)
#'
#' @description
#' Estimates chronological or biological age from single-cell RNA-sequencing (scRNA-seq) data 
#' derived from the subventricular zone (SVZ). This function implements the bootstrap 
#' pseudocells methodology described by Buckley et al., integrating cell-type-specific 
#' Elastic Net models to quantify aging and rejuvenation.
#' 
#'
#' @param seuratObj A \code{Seurat} object containing the scRNA-seq data. 
#'   Must include raw counts and required metadata (\code{donorId} or \code{donor_id}, 
#'   \code{age}, and \code{celltype} column).
#' @param cellTypes Character vector. One or multiple cell types (e.g., 
#' c("aNSC_NPC", "Astrocyte_qNSC", "Endothelial", "Microglia", "Neuroblast", "Oligodendro")).
#' @param clockType Character. "chronological" or "biological".
#' @param minCoverage A numeric value (0-1). Minimum proportion of required genes.
#' @param verbose Logical flag for printing status messages.
#' 
#' @return A data frame containing donorId, age, prediction, and cellType.
#' 
#' @details
#' The prediction pipeline strictly mirrors the original methodology:
#' \enumerate{
#'   \item Extracts raw transcript counts for the specified cell types.
#'   \item Generates pseudocells via bootstrap resampling (aggregating 15 cells per pseudocell using summation).
#'   \item Applies log-normalization scaled to 10,000 transcripts per pseudocell (CP10k + log1p).
#'   \item Computes the final predicted scores via matrix multiplication against the pre-trained model weights.
#' }
#' For the \code{"biological"} clock, raw predictions (proliferation fractions) are mathematically 
#' transformed into a biological age score using the formula: \code{35 - (100 * fraction)}.
#'
#' 
#' @references
#' Buckley MT, Sun ED, George BM, et al.
#' Cell-type-specific aging clocks to quantify aging and rejuvenation in neurogenic regions of the brain.
#' \emph{Nat Aging.} 2023
#' 
#' @export
#'
#' @examples
#' # 1. Define valid cell types (Runnable code to satisfy BiocCheck)
#' valid_cell_types <- c("aNSC_NPC", "Astrocyte_qNSC", "Endothelial",
#'                        "Microglia", "Neuroblast", "Oligodendro")
#' print(valid_cell_types)
#'
#' \donttest{
#' # 2. Real pipeline execution (Wrapped in donttest because it requires
#' # downloading pre-trained models and external example datasets)
#' library(Seurat)
#' seuratObj <- loadOmniAgeRdata(
#'     "omniager_buckley_mouse_example",
#'     verbose = FALSE
#' )
#' set.seed(42)
#' predOut <- buckleyMouseSVZ(seuratObj, c("Microglia", "Oligodendro"))
#' }
#' 
#' 
buckleyMouseSVZ <- function(seuratObj, cellTypes, clockType = "chronological", 
                            minCoverage = 0, verbose = TRUE) {
  
  # 1. load model
  clockModelList <- loadOmniAgeRdata(
    "omniager_buckley_mouse_bootstrap_model",
    verbose = verbose
  )
  
  # 2. According to the user's selection, map to the correct model level
  if (tolower(clockType) == "chronological") {
    modelGroup <- "Chronological_Bootstrap"
  } else if (tolower(clockType) == "biological") {
    modelGroup <- "Biological_Bootstrap"
  } else {
    stop("clockType must be 'chronological' or 'biological'")
  }
  
  # Check whether the cell type column exists in the metadata
  if (! "celltype" %in% colnames(seuratObj@meta.data)) {
    stop(paste("Column celltype not found in Seurat metadata."))
  }
  
  # 3.  Traverse all the cell types provided by the user
  predictions_list <- lapply(cellTypes, function(cType) {
    
    if (verbose) message(sprintf("Processing cell type: %s...", cType))
    
    if (!cType %in% names(clockModelList[[modelGroup]])) {
      warning(sprintf("Model for cell type '%s' not found in %s clock. Skipping.", cType, clockType))
      return(NULL)
    }
    weightDf <- clockModelList[[modelGroup]][[cType]]
    
    cells_to_keep <- rownames(seuratObj@meta.data)[seuratObj@meta.data[["celltype"]] == cType]
    if (length(cells_to_keep) == 0) {
      warning(sprintf("No cells found for cell type '%s' in the Seurat object. Skipping.", cType))
      return(NULL)
    }
    subSeurat <- subset(seuratObj, cells = cells_to_keep)
    
    #  Extract the original counts (assayLayer = "counts")
    preprocessedData <- scImmuAgingPreprocessing(subSeurat, assayLayer = "counts") %>%
      dplyr::group_by(donorId, age) %>%
      tidyr::nest()
    
    # Bootsrap (aggrMethod = "sum")
    preprocessedData <- preprocessedData %>%
      dplyr::mutate(pseudocells = purrr::map(data, ~pseudocellScImmuAging(.x, size = 15, n = 100, aggrMethod = "sum")))
    
    preprocessedData$data <- NULL
    unNestedData <- tidyr::unnest(preprocessedData, pseudocells)
    
    # Normalization
    meta_cols <- unNestedData[, c("donorId", "age")]
    count_matrix <- as.matrix(unNestedData[, -c(1, 2)])
    
    normed <- sweep(count_matrix, MARGIN = 1, FUN = "/", STATS = rowSums(count_matrix))
    logged_matrix <- log1p(normed * 10000)
    
    final_input <- dplyr::as_tibble(cbind(meta_cols, logged_matrix))
    
    # Predcition
    clock_name_log <- paste("buckleyMouseSVZ", clockType, cType, sep = "_")
    preds <- predictSCRNAClock(
      preprocessedData = final_input,
      weightDf = weightDf,
      minCoverage = minCoverage,
      clockName = clock_name_log,
      verbose = verbose
    )
    
    # Post-treatment of biological age
    if (tolower(clockType) == "biological") {
      valid_idx <- !is.na(preds$prediction)
      preds$prediction[valid_idx] <- 35 - (100 * preds$prediction[valid_idx])
    }
    
    preds$cellType <- cType
    
    predsDonor <- ageDonor(preds)
    return(predsDonor)
  })
  
  names(predictions_list) <- cellTypes
  predictions_list <- Filter(Negate(is.null), predictions_list)
  
  if (length(predictions_list) == 0) {
    stop("No valid predictions could be generated for any of the provided cell types.")
  }
  

  return(predictions_list)
}




#' General Clock Score Calculator (Supports Dynamic Metadata & Imputation)
#'
#' @param preprocessedData A data frame/tibble output from pre-processing.
#' @param weightDf A data frame of weights with 'feature_name' and 'coefficient'.
#' @param minCoverage Minimum proportion of required genes (default 0).
#' @param imputeData Optional. A data frame with 'feature_name' and 'imputation_value'
#'   to fill in missing genes. If NULL, missing genes are filled with 0.
#' @param clockName Character for logging.
#' @param metaColNames Character vector specifying which columns are metadata. 
#'   Default is c("donorId", "age").
#' @param sampleType A character string (e.g., 'SC', 'Pseudobulk') used to
#'   tag the output data.frame.
#'   
#' @param verbose Logical flag.
#' @return Return dataframe of predicted scores safely merged with metadata.
#' @import dplyr
#' @importFrom Matrix Matrix
#' @export
predictSCRNAClock <- function(preprocessedData, weightDf, minCoverage = 0, 
                                imputeData = NULL, clockName, 
                                metaColNames = c("donorId", "age"), 
                                sampleType = "Bootstrap",
                                verbose = TRUE) {
  
  # ---------------------------------------------------------
  # 1. Parse the weight data frame
  # ---------------------------------------------------------
  interceptIdx <- tolower(weightDf$feature_name) %in% c("intercept", "(intercept)")
  
  if (any(interceptIdx)) {
    intercept <- weightDf$coefficient[interceptIdx][1]
    geneWeights <- weightDf[!interceptIdx, ]
  } else {
    intercept <- 0
    geneWeights <- weightDf
  }
  markerGenes <- geneWeights$feature_name
  
  # ---------------------------------------------------------
  # 2. Ultra-fast matrix extraction & Metadata purification
  # ---------------------------------------------------------
  availableGenes <- intersect(colnames(preprocessedData), markerGenes)
  testMat <- t(as.matrix(preprocessedData[, availableGenes, drop = FALSE]))
  
  
  missing_meta <- setdiff(metaColNames, colnames(preprocessedData))
  if (length(missing_meta) > 0) {
    if (verbose) {
      warning(sprintf("The following metadata columns were not found and will be ignored: %s", 
                      paste(missing_meta, collapse = ", ")))
    }
    validMetaCols <- intersect(metaColNames, colnames(preprocessedData))
  } else {
    validMetaCols <- metaColNames
  }
  
  if (length(validMetaCols) > 0) {
    metaData <- preprocessedData[, validMetaCols, drop = FALSE]
  } else {
    metaData <- data.frame(row.names = rownames(preprocessedData))
  }
  
  # ---------------------------------------------------------
  # 3. Coverage check
  # ---------------------------------------------------------
  fakeWeights <- rep(0, length(markerGenes))
  names(fakeWeights) <- markerGenes
  
  coverage <- .checkCpGCoverage(
    betaM = testMat,
    allWeights = fakeWeights,
    clockName = clockName,
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  # If the coverage is insufficient, simply return the metaData with "NA" included.
  if (!coverage$pass) {
    res_failed <- metaData
    res_failed$prediction <- NA_real_
    return(res_failed)
  }
  
  # ---------------------------------------------------------
  # 4. Matrix alignment and intelligent interpolation
  # ---------------------------------------------------------
  presentGenes <- rownames(testMat)[coverage$betaIdx]
  missingGenes <- setdiff(markerGenes, presentGenes)
  
  finalMat <- matrix(0, nrow = length(markerGenes), ncol = ncol(testMat))
  rownames(finalMat) <- markerGenes
  
  finalMat[presentGenes, ] <- testMat[presentGenes, ]
  
  if (length(missingGenes) > 0) {
    if (!is.null(imputeData)) {
      fillValues <- setNames(rep(0, length(missingGenes)), missingGenes)
      matchIdx <- match(missingGenes, imputeData$feature_name)
      fillValues[!is.na(matchIdx)] <- imputeData$imputation_value[matchIdx[!is.na(matchIdx)]]
    } else {
      fillValues <- rep(0, length(missingGenes))
    }
    finalMat[missingGenes, ] <- fillValues
  }
  
  # ---------------------------------------------------------
  # 5. Perform matrix multiplication
  # ---------------------------------------------------------
  weightVector <- matrix(geneWeights$coefficient, ncol = 1)
  rownames(weightVector) <- geneWeights$feature_name
  weightVector <- weightVector[rownames(finalMat), , drop = FALSE]
  
  rawPredictions <- t(finalMat) %*% weightVector
  finalPredictions <- as.numeric(rawPredictions) + intercept
  
  # ---------------------------------------------------------
  # 6. Rerurn result
  # ---------------------------------------------------------
  res_success <- metaData
  res_success$prediction <- finalPredictions
  res_success$sampleType <- sampleType
  return(res_success)
}









