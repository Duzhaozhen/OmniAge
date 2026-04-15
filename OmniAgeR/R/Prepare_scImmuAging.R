#' Pre-Processing of input dataset for scImmuAging
#'
#' @param inputObj A Seurat object. Meta data must contain "donorId" and
#' "age" columns.
#' @param cellType A string specifying the cell type.
#' @param model The pre-trained model object.
#' @param markerGenes Character vector of selected marker genes.
#' @return A nested and unnested data frame with pseudocells.
#' @import Seurat dplyr purrr tidyr
#' @export
#' @examples
#' library(Seurat)
#' library(dplyr)
#' library(tidyr)
#' library(purrr)
#'
#' # 1. Define mock marker genes
#' mock_markers <- c("GeneA", "GeneB", "GeneC")
#'
#' # 2. Create a tiny mock count matrix (3 genes, 20 cells)
#' set.seed(123)
#' mock_counts <- matrix(rpois(60, lambda = 5), nrow = 3, ncol = 20)
#' rownames(mock_counts) <- mock_markers
#' colnames(mock_counts) <- paste0("Cell_", 1:20)
#'
#' # 3. Create mock metadata
#' mock_meta <- data.frame(
#'     donorId = rep(c("Donor1", "Donor2"), each = 10),
#'     age = rep(c(30, 60), each = 10),
#'     row.names = colnames(mock_counts)
#' )
#'
#' # 4. Build Seurat object and normalize (required to create "data" layer)
#' mock_seurat <- CreateSeuratObject(
#'     counts = mock_counts,
#'     meta.data = mock_meta
#' )
#' mock_seurat <- NormalizeData(mock_seurat, verbose = FALSE)
#'
#' # 5. Run the preprocessing pipeline
#' res <- scImmuAgingPreProcess(
#'     inputObj = mock_seurat,
#'     cellType = "T_cell",
#'     model = NULL,
#'     markerGenes = mock_markers
#' )
#'
#' # 6. View the generated pseudocells
#' print(res)
scImmuAgingPreProcess <- function(inputObj, cellType, model, markerGenes) {
    newSeurat <- subset(inputObj, features = markerGenes)

    preprocessedData <- scImmuAgingPreprocessing(newSeurat) %>%
        dplyr::group_by(donorId, age) %>%
        tidyr::nest()

    preprocessedData <- preprocessedData %>%
        dplyr::mutate(pseudocells = purrr::map(data, pseudocellScImmuAging))

    preprocessedData$data <- NULL
    preprocessedData <- tidyr::unnest(preprocessedData, pseudocells)

    return(preprocessedData)
}


#' Pre-Processing of input dataset
#'
#' input is a seurat object. meta.data should have "donor_id" and "age" columns.
#' @param seuratObj seurat object with metadata "donor_id", "age"
#' @param assayLayer Character. "data" (default) for pre-normalized.
#' @return Return pre-processed input
#' @import Seurat dplyr glmnet purrr
#' @export
#' @examples
#' library(Seurat)
#'
#' # 1. Create a tiny mock count matrix (5 genes, 10 cells)
#' mock_counts <- matrix(rpois(50, lambda = 10), nrow = 5, ncol = 10)
#' rownames(mock_counts) <- paste0("Gene", 1:5)
#' colnames(mock_counts) <- paste0("Cell", 1:10)
#'
#' # 2. Create mock metadata with required columns
#' mock_meta <- data.frame(
#'     donor_id = rep(c("Donor_A", "Donor_B"), each = 5),
#'     age = rep(c(25, 65), each = 5),
#'     row.names = colnames(mock_counts)
#' )
#'
#' # 3. Create Seurat object and run basic normalization
#' # (since your function extracts the "data" layer)
#' mock_seurat <- CreateSeuratObject(counts = mock_counts, meta.data = mock_meta)
#' mock_seurat <- NormalizeData(mock_seurat, verbose = FALSE)
#'
#' # 4. Run your preprocessing function
#' processedData <- scImmuAgingPreprocessing(mock_seurat)
#'
#' # 5. View result
#' print(processedData)


scImmuAgingPreprocessing <- function(seuratObj, assayLayer = "data") {
  DefaultAssay(seuratObj) <- "RNA"
  metaData <- seuratObj@meta.data
  
  if (!"donorId" %in% colnames(metaData)) {
    if ("donor_id" %in% colnames(metaData)) {
      metaData$donorId <- metaData$donor_id
    } else {
      stop("please add 'donorId' (or 'donor_id') in your metadata!")
    }
  }
  
  if (!"age" %in% colnames(metaData)) {
    stop("please add 'age' in your metadata!")
  }
  
  metaDataSubset <- metaData[, c("donorId", "age"), drop = FALSE]

  inputMtx <- t(as.matrix(GetAssayData(seuratObj, assay = "RNA", layer = assayLayer)))
  
  combinedInput <- dplyr::as_tibble(cbind(metaDataSubset, inputMtx))
  return(combinedInput)
}




#' Pre-Processing of input dataset
#'
#' input could be interested gene list. Each row is one gene and columns are
#' "gene" and "value". "value" could be log2FC or specific value.
#' colnames of other metadata should be "subtype".
#' @param inputData the output from scImmuAgingPreProcess() function.
#' @param size how many cells would be randomly selected to generate pseudocells
#' @param n how many times to repeat random selection
#' @param replace Character or Logical. Strategy for sampling.
#'   If "dynamic" (default), it automatically sets replace to TRUE
#'   if nrow(input) <= size. Otherwise, can be set to TRUE or FALSE explicitly.
#' @param aggrMethod Character. "mean" (default) for original scImmuAging, "sum" for SVZ bootstrap clocks.  
#'   
#' @return Return ranked data.frame based on "value"
#' @import Seurat dplyr glmnet purrr
#' @export
#' @examples
#' # 1. Create a small mock dataset representing cells and genes
#' # (e.g., 20 cells as rows, 3 genes as columns)
#' set.seed(123) # Set seed for reproducibility in example
#' mock_input <- data.frame(
#'     GeneA = runif(20),
#'     GeneB = runif(20),
#'     GeneC = runif(20)
#' )
#' rownames(mock_input) <- paste0("Cell_", 1:20)
#'
#' # 2. Generate 3 pseudocells, each averaging 5 random cells
#' pseudo_res <- pseudocellScImmuAging(
#'     inputData = mock_input,
#'     size = 5,
#'     n = 3
#' )
#'
#' # 3. View the generated pseudocells
#' print(pseudo_res)


pseudocellScImmuAging <- function(inputData, size = 15, n = 100, replace = "dynamic", aggrMethod = "mean") {
  if (replace == "dynamic") {
    replace <- nrow(inputData) <= size
  }
  
  mat <- as.matrix(inputData)
  numRows <- nrow(mat)
  
  indices <- replicate(n, sample(seq_len(numRows), size = size, replace = replace))
  
  pseudoMat <- vapply(seq_len(n), function(i) {
    sub_mat <- mat[indices[, i], , drop = FALSE]
    if (aggrMethod == "sum") {
      colSums(sub_mat)
    } else {
      colMeans(sub_mat)
    }
  }, numeric(ncol(mat)))
  
  return(dplyr::as_tibble(t(pseudoMat)))
}






#' predict age for each individual
#'
#' @param preprocessedData the output from scImmuAging_PreProcess() function.
#' @param model the cell type aging clock model
#' @param markerGenes selected marker genes from the corresponding model
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required genes that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @return Return dataframe of predicted age
#' @import Seurat dplyr glmnet purrr
#' @export
#' @examples
#' library(glmnet)
#'
#' # 1. Create a tiny mock preprocessed dataset (3 donors, 4 genes)
#' mock_preprocessed <- data.frame(
#'     donorId = c("Donor1", "Donor2", "Donor3"),
#'     age = c(25, 45, 65),
#'     GeneA = c(1.2, 1.5, 2.1),
#'     GeneB = c(0.8, 0.9, 1.1),
#'     GeneC = c(3.3, 2.8, 1.5),
#'     GeneD = c(0.1, 0.5, 0.9)
#' )
#'
#' # 2. Train a valid mock cv.glmnet model
#' set.seed(123)
#' # Generate 30 rows of random data to avoid cross-validation warnings
#' X_train <- matrix(rnorm(120), nrow = 30, ncol = 4)
#' colnames(X_train) <- c("GeneA", "GeneB", "GeneC", "GeneD")
#'
#' # Create a fake linear relationship so the model learns non-zero coefficients
#' Y_train <- 40 + (5 * X_train[, "GeneA"]) - (3 * X_train[, "GeneB"]) + rnorm(30)
#' mock_model <- cv.glmnet(X_train, Y_train)
#'
#' # 3. Define the marker genes expected by the model
#' mock_markers <- c("GeneA", "GeneB", "GeneC", "GeneD")
#'
#' # 4. Run the calculator function
#' result <- scImmuAgingCalculator(
#'     preprocessedData = mock_preprocessed,
#'     model = mock_model,
#'     markerGenes = mock_markers,
#'     verbose = FALSE
#' )
#'
#' # 5. View the output
#' print(result)
scImmuAgingCalculator <- function(preprocessedData, model, markerGenes,
                                  minCoverage = 0.5, verbose = TRUE) {
    testMat <- t(as.matrix(preprocessedData[, -c(1, 2)]))


    fakeWeights <- rep(0, length(markerGenes))
    names(fakeWeights) <- markerGenes

    # 1. Check the coverage
    coverage <- .checkCpGCoverage(
        betaM = testMat,
        allWeights = fakeWeights,
        clockName = "scImmuAging",
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (!coverage$pass) {
        return(data.frame(
            donorId = preprocessedData[[1]],
            age = preprocessedData[[2]],
            prediction = NA_real_
        ))
    }

    # 2. Prepare the final matrix
    finalMat <- matrix(0, nrow = length(markerGenes), ncol = ncol(testMat))
    rownames(finalMat) <- markerGenes
    finalMat[rownames(testMat)[coverage$betaIdx], ] <- testMat[coverage$betaIdx, ]

    # 3. Prediction
    testPredictions <- predict(model, newx = t(finalMat), s = "lambda.min")

    return(data.frame(
        donorId = preprocessedData[[1]],
        age = preprocessedData[[2]],
        prediction = as.numeric(testPredictions)
    ))
}


#' predicted age for each individual
#'
#' @param predictRes the output from scImmuAgingCalculator() function
#' @return Return predicted age for each individual
#' @import Seurat dplyr glmnet purrr
#' @export
#' @examples
#' # 1. Create a mock output data.frame from scImmuAgingCalculator
#' mockPredictRes <- data.frame(
#'     donorId = c("Donor_A", "Donor_A", "Donor_A", "Donor_B", "Donor_B"),
#'     celltype = c("nCD4T", "nCD8T", "NK", "nCD4T", "nCD8T"),
#'     prediction = c(45.2, 46.8, 44.5, 61.0, 59.5)
#' )
#'
#' # 2. Run the ageDonor function to aggregate predictions by donor
#' donorAges <- ageDonor(mockPredictRes)
#'
#' # 3. View the result
#' print(donorAges)
#'
ageDonor <- function(predictRes) {
    donorDf <- predictRes %>%
        dplyr::group_by(donorId) %>%
        #dplyr::mutate(predicted = round(mean(prediction, na.rm = TRUE))) %>%
        dplyr::mutate(predicted = mean(prediction, na.rm = TRUE)) %>%
        dplyr::select(-prediction) %>%
        dplyr::distinct()

    return(as.data.frame(donorDf))
}
