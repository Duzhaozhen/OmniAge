#' Run Brain Cell Type Specific Clocks
#'
#' @description
#' A high-level wrapper function that runs the brain age prediction pipeline for
#' one or more specified sample types ('SC', 'Pseudobulk', 'Bootstrap').
#'
#' @param seuratObj The input Seurat object containing expression data and
#'   metadata (must include 'donor_id', 'age', 'celltype').
#' @param cellTypes A character vector of cell types to analyze
#'   (e.g., `c('Oligodendrocytes', 'Astrocytes')`).
#'   Available cell types are: "Oligodendrocytes", "Astrocytes", "Microglia",
#'   "OPCs", "Excitatory Neurons", and "Inhibitory Neurons".
#' @param modelNames A character string or vector specifying which models to run.
#'   \itemize{
#'     \item `"all"` (default): Runs "SC", "Pseudobulk", and "Bootstrap".
#'     \item A vector: e.g., `c("SC", "Pseudobulk")` will run only those two.
#'   }
#' @param verbose Logical. Whether to print status messages.
#' @return
#' A named list where each element corresponds to a `modelNames` that was run.
#' Each element contains a data.frame of the (5-fold averaged) predictions,
#' as returned by `runPredictionPipelineBrainCt`.
#'
#' @details
#' This function serves as the primary endpoint for the prediction pipeline.
#' It iteratively calls `runPredictionPipelineBrainCt` for each requested
#' `modelNames` (sample type) and collects the results into a single list.
#'
#' @seealso \code{\link{runPredictionPipelineBrainCt}}
#'
#' @references
#' Muralidharan C, Zakar-Polyák E, Adami A, et al.
#' Human Brain Cell-Type-Specific Aging Clocks
#' Based on Single-Nuclei Transcriptomics.
#' \emph{Adv Sci(Weinh).} 2025
#' @export
#' @examples
#' # 1. Define valid models to run (Fast, runnable example to pass BiocCheck)
#' available_models <- c("SC", "Pseudobulk", "Bootstrap")
#' print(available_models)
#' \donttest{
#' # 2. Real pipeline execution (Wrapped in donttest because it requires downloading
#' # pre-trained models and large example datasets via loadOmniAgeRdata)
#'
#' # Load the Seurat object
#' library(Seurat)
#' library(dplyr)
#' brainSeurat <- loadOmniAgeRdata(
#'     "omniager_brain_frohlich_control_example_15donors",
#'     verbose = FALSE
#' )
#'
#' # Define cell types of interest
#' cellTypes <- c("Oligodendrocytes")
#'
#' # Run the model for the specified cell types
#' clockResults <- brainCtClock(
#'     seuratObj = brainSeurat,
#'     cellTypes = cellTypes,
#'     modelNames = "SC"
#' )
#' }
brainCtClock <- function(seuratObj, cellTypes,
                         modelNames = "all", verbose = TRUE) {
    # --- 1. Define which models to run  ---
    validModels <- c("SC", "Pseudobulk", "Bootstrap")
    if (length(modelNames) == 1 && modelNames == "all") {
        modelsToRun <- validModels
    } else {
        modelsToRun <- modelNames
        if (!all(modelsToRun %in% validModels)) {
            stop(
                "Invalid modelNames.",
                " Valid: 'SC', 'Pseudobulk', 'Bootstrap', or 'all'."
            )
        }
    }

    allResults <- list()
    if (verbose) message("Starting brainCtClock for cell types...")

    # --- 2. Loop over each requested model type ---
    for (currentModelType in modelsToRun) {
        predResult <- runPredictionPipelineBrainCt(
            sampleType = currentModelType,
            seuratObj = seuratObj,
            cellTypes = cellTypes,
            verbose = verbose
        )
        allResults[[currentModelType]] <- predResult
    }

    if (verbose) message("\n--- All Brain clock predictions complete! ---")
    return(allResults)
}


#' Extract and Pre-process Data from a Seurat Object
#'
#' @description
#' Subsets a Seurat object by cell type and processes the expression data into
#' one of three formats: single-cell ('SC'), 'Pseudobulk' (averaged by donor),
#' or 'Bootstrap' (resampled within donor).
#'
#' @param seuratObj A Seurat object containing the expression data and metadata
#'   (must include 'donor_id', 'age', 'celltype').
#' @param cellType A character string specifying the cell type to subset.
#' @param sampleType A character string defining the processing method:
#'   "SC", "Pseudobulk", or "Bootstrap".
#'
#' @return
#' A data.frame (tibble) containing the processed expression data and metadata.
#' The structure depends on `sampleType`:
#' \itemize{
#'   \item{"SC": A cell-by-gene matrix with metadata.}
#'   \item{"Pseudobulk": A donor-by-gene matrix (mean expression) with metadata.}
#'   \item{"Bootstrap": A (donor * 100 replicates)-by-gene matrix with metadata.}
#' }
#' Returns an empty data.frame if no cells are found for the specified `cellType`.
#'
#' @details
#' This function defines the number of cells to sample for the 'Bootstrap' method
#' internally via the `numCellsMap` list. If a donor has fewer cells than the
#' specified number, 100 replicates of the donor's mean expression are returned.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom dplyr as_tibble group_by summarize across all_of bind_rows
#' @importFrom magrittr %>%
#' @export
#' @examples
#' library(Seurat)
#' library(dplyr)
#'
#' # 1. Create a tiny mock Seurat object (5 genes, 20 cells)
#' set.seed(123)
#' mock_counts <- matrix(rpois(100, lambda = 5), nrow = 5, ncol = 20)
#' rownames(mock_counts) <- paste0("Gene", 1:5)
#' colnames(mock_counts) <- paste0("Cell_", 1:20)
#'
#' # Metadata with 2 donors and 2 cell types
#' mock_meta <- data.frame(
#'     donor_id = rep(c("Donor1", "Donor2"), each = 10),
#'     age = rep(c(30, 60), each = 10),
#'     celltype = rep(c("Astrocytes", "Microglia"), times = 10),
#'     row.names = colnames(mock_counts)
#' )
#'
#' mock_seu <- CreateSeuratObject(counts = mock_counts, meta.data = mock_meta)
#' mock_seu <- NormalizeData(mock_seu, verbose = FALSE)
#'
#' # 2. Run getDfSeurat to generate Pseudobulk data for Astrocytes
#' res_pb <- getDfSeurat(
#'     seuratObj = mock_seu,
#'     cellType = "Astrocytes",
#'     sampleType = "Pseudobulk"
#' )
#'
#' # 3. View the generated Pseudobulk dataframe
#' print(res_pb)
#'
getDfSeurat <- function(seuratObj, cellType, sampleType) {
    numCellsMap <- list(
        "Oligodendrocytes" = 200, "Astrocytes" = 50, "Microglia" = 50,
        "OPCs" = 50, "Excitatory Neurons" = 100, "Inhibitory Neurons" = 100
    )

    # 1. Subset Seurat Object
    seuratSub <- subset(seuratObj, subset = celltype == cellType)
    if (ncol(seuratSub) == 0) {
        warning("No cells found for: ", cellType)
        return(data.frame())
    }

    # 2. Unify the column names of metadata
    metaData <- seuratSub[[]]
    metaData$cellId <- rownames(metaData)
    
    if (!"donorId" %in% colnames(metaData) && "donor_id" %in% colnames(metaData)) {
        metaData$donorId <- metaData$donor_id
    }
    
    # 3. Extract expression matrix
    exprMtx <- t(as.matrix(GetAssayData(seuratSub, assay = "RNA", layer = "data")))
    combinedData <- dplyr::as_tibble(cbind(metaData[, c("cellId","donorId", "age", "celltype")], exprMtx))
    geneCols <- colnames(exprMtx)

    # 4. Process data based on sample_type
    if (sampleType == "SC") {
        return(combinedData)
    } else if (sampleType == "Pseudobulk") {
        return(combinedData %>%
            dplyr::group_by(.data$donorId, .data$age, .data$celltype) %>%
            dplyr::summarize(dplyr::across(dplyr::all_of(geneCols), .fns = mean), .groups = "drop"))
    } else if (sampleType == "Bootstrap") {
        donors <- unique(combinedData$donorId)
        nSample <- numCellsMap[[cellType]]
        if (is.null(nSample)) {
            nSample <- 50
        }

        bootstrapList <- lapply(donors, function(d) {
            dfDonor <- combinedData[combinedData$donorId == d, ]
            numRows <- nrow(dfDonor)

            reps <- 100
            indices <- replicate(reps, sample(seq_len(numRows), size = nSample, replace = (numRows < nSample)))

            bootMat <- vapply(seq_len(reps), function(i) {
                colMeans(as.matrix(dfDonor[indices[, i], geneCols, drop = FALSE]))
            }, numeric(length(geneCols)))

            dfBoot <- dplyr::as_tibble(t(bootMat))
            dfBoot$donorId <- d
            dfBoot$age <- dfDonor$age[1]
            dfBoot$celltype <- cellType
            return(dfBoot)
        })
        return(dplyr::bind_rows(bootstrapList))
    }
}



# --- Run Prediction Flow ---

#' Apply a Single Clock Model to Expression Data (Vectorized)
#'
#' @description
#' Predicts age using a pre-trained elastic net model (coefficients and intercept)
#' on a given expression data matrix. It handles gene matching, imputation for
#' missing genes, and vectorized prediction.
#'
#' @param inputData A data.frame of expression data (samples in rows, genes in
#'   columns). Must also contain 'age' and 'donor_id' columns.
#' @param imputeData A long-format data.frame with 'feature_name' and
#'   'imputation_value' columns. Used to fill in genes present in the model
#'   but missing from `data`.
#' @param modelObj A long-format data.frame with 'feature_name' (genes + 'intercept')
#'   and 'coefficient' columns, representing one trained clock.
#' @param sampleType A character string (e.g., 'SC', 'Pseudobulk') used to
#'   tag the output data.frame.
#'
#' @return
#' A data.frame with 'predictions' (the predicted age), 'ages' (the true age),
#' 'donors', and 'sample_type'.
#'
#' @details
#' This function is the core prediction engine. It performs a matrix
#' multiplication (`expression_matrix %*% coefficients_vector + intercept`).
#' It ensures that the gene order in the expression matrix exactly matches the
#' coefficient order from the model.
#'
#' @importFrom dplyr filter left_join bind_cols
#' @importFrom tidyr pivot_wider
#' @noRd

predictBrainCtAge <- function(inputData, imputeData, modelObj, sampleType) {
  # 1. Separate coefficients and intercept
  intercept <- modelObj$coefficient[modelObj$feature_name == "intercept"]
  if (length(intercept) == 0) intercept <- 0
  
  modelGenesDf <- modelObj[modelObj$feature_name != "intercept", ]
  modelGenes <- modelGenesDf$feature_name
  
  # 2. Gene matching
  presentGenes <- intersect(modelGenes, colnames(inputData))
  missingGenes <- setdiff(modelGenes, presentGenes)
  
  exprMat <- as.matrix(inputData[, presentGenes, drop = FALSE])
  
  if (length(missingGenes) > 0) {
    fillValues <- setNames(rep(0, length(missingGenes)), missingGenes)
    matchIdx <- match(missingGenes, imputeData$feature_name)
    fillValues[!is.na(matchIdx)] <- imputeData$imputation_value[matchIdx[!is.na(matchIdx)]]
    
    fillMat <- matrix(rep(fillValues, each = nrow(exprMat)), nrow = nrow(exprMat))
    colnames(fillMat) <- missingGenes
    fullMat <- cbind(exprMat, fillMat)
  } else {
    fullMat <- exprMat
  }
  
  # 3. Strictly sort and calculate
  fullMat <- fullMat[, modelGenes, drop = FALSE]
  coefVector <- modelGenesDf$coefficient[match(modelGenes, modelGenesDf$feature_name)]
  
  predVec <- (fullMat %*% coefVector) + intercept
  
  return(data.frame(
    prediction = as.numeric(predVec),
    age = inputData$age,
    donorId = inputData$donorId,
    sampleType = sampleType
  ))
}

#' Run the Full 5-Fold Averaged Prediction Pipeline
#'
#' @description
#' Orchestrates the entire prediction process for a given sample type and set of
#' cell types. It loads the 5-fold cross-validation models, runs predictions
#' for each of the 5 folds, and returns the final averaged prediction for
#' each sample.
#'
#' @param sampleType The processing method to use: "SC", "Pseudobulk", or
#'   "Bootstrap". This is passed to `getDfSeurat`.
#' @param seuratObj The input Seurat object, passed to `getDfSeurat`.
#' @param cellTypes A character vector of cell types to process
#'   (e.g., `c('Oligodendrocytes', 'Astrocytes')`).
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @return
#' A single, combined data.frame containing the final (5-fold averaged)
#' 'predictions', 'ages', 'donors', 'sample_type', and 'celltype' for all
#' requested cell types.
#'
#' @details
#' This is the main user-facing function for this pipeline. It assumes
#' that the dataset "omniager_brain_celltype_specific_clocks_coef" (loaded
#' via \code{loadOmniAgeRdata}) is a single list containing two main components:
#' \itemize{
#'   \item \code{brain_ct_clocks_coef}: A nested list where keys
#'     (e.g., "SC_Oligodendrocytes") map to a list of 5 models (one for each fold).
#'   \item \code{Brain_CT_imputation_data_list}: A list where keys
#'     (e.g., "SC_Oligodendrocytes") map to a single imputation data.frame.
#' }
#' The function calls \code{getDfSeurat} to prepare the data, then calls
#' \code{predictBrainCtAge} five times (once for each fold). The five resulting
#' prediction vectors are then averaged (row-wise) to produce the final,
#' stable prediction.
#'
#' @importFrom dplyr bind_rows
#' @export
#' @examples
#' # 1. Define valid parameters (Runnable code to satisfy BiocCheck)
#' sample_type <- "SC"
#' target_cells <- c("Oligodendrocytes", "Astrocytes")
#' print(paste("Preparing to run pipeline for", sample_type, "samples."))
#'
#' \donttest{
#' # 2. Real pipeline execution (Wrapped in donttest because it requires
#' # downloading pre-trained models and large example datasets)
#' library(Seurat)
#'
#' # Load the example brain Seurat object
#' brainSeurat <- loadOmniAgeRdata(
#'     "omniager_brain_frohlich_control_example_15donors",
#'     verbose = FALSE
#' )
#'
#' # Run the 5-fold prediction pipeline for Oligodendrocytes
#' results <- runPredictionPipelineBrainCt(
#'     sampleType = "SC",
#'     seuratObj = brainSeurat,
#'     cellTypes = c("Oligodendrocytes"),
#'     verbose = FALSE
#' )
#'
#' print(head(results))
#' }
runPredictionPipelineBrainCt <- function(sampleType, seuratObj,
                                         cellTypes, verbose = TRUE) {
    # Load Models and Imputation Data
    brainCtClocksCoef <- loadOmniAgeRdata(
        "omniager_brain_celltype_specific_clocks_coef",
        verbose = verbose
    )
    brainCtClocksCoef <- brainCtClocksCoef[["brain_ct_clocks_coef"]]
    brainCtImputationList <- brainCtClocksCoef[["Brain_CT_imputation_data_list"]]

    finalResultsList <- list()

    for (ct in cellTypes) {
        clockKey <- paste(sampleType, ct, sep = "_")
        modelFolds <- brainCtClocksCoef[[clockKey]]
        imputeData <- brainCtImputationList[[clockKey]]

        if (is.null(modelFolds)) {
            warning("Skipping ", ct, ": Model not found for ", clockKey)
            next
        }

        # 1. Obtain preprocessed data
        dfBase <- getDfSeurat(seuratObj, ct, sampleType)
        if (nrow(dfBase) == 0) next

        # 2. 5-Fold prediction
        foldPreds <- vapply(names(modelFolds), function(f) {
            res <- predictBrainCtAge(dfBase, imputeData, modelFolds[[f]], sampleType)
            return(res$prediction)
        }, numeric(nrow(dfBase)))

        # 3. Calculate the average and construct the result
        res_df <- data.frame(
            prediction = rowMeans(foldPreds, na.rm = TRUE),
            age = dfBase$age,
            donorId = dfBase$donorId,
            sampleType = sampleType,
            celltype = ct
        )
        
        if ("cellId" %in% colnames(dfBase)) {
          res_df$cellId <- dfBase$cellId
        }
        finalResultsList[[ct]] <- res_df
        
    }

    return(dplyr::bind_rows(finalResultsList))
}
