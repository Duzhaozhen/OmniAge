#' @title Calculate PASTA-based transcriptomic age scores
#'
#' @description
#' This function computes transcriptomic age acceleration scores based on the
#' models described by Salignon et al. (2025), including the primary PASTA
#' score, a standard regression (REG) score, and the CT46 score.
#'
#' @param mat An expression matrix (numeric) with **genes as rows** and
#'   **samples as columns**.
#' @param filterGenes Logical. If \code{TRUE} (default), the matrix is
#'   subsetted to retain only the genes utilized by the pre-trained models.
#' @param rankNorm Logical. If \code{TRUE} (default), applies a rank-based
#'   inverse normal transformation (rank-normalization) to the expression data.
#' @param reg Logical. If \code{TRUE} (default), computes the REG (regression)
#'   age score.
#' @param pasta Logical. If \code{TRUE} (default), computes the PASTA age score.
#' @param ct46 Logical. If \code{TRUE} (default), computes the CT46 age score.
#' @param verbose Logical. Whether to print status messages.
#'   Default is \code{TRUE}.
#'
#' @return
#' A \code{list} where each element is a numeric vector of predicted age scores
#' for the requested model(s) (e.g., \code{res_list$PASTA}, \code{res_list$REG}).
#'
#' @export
#'
#' @references
#' Salignon, J., Tsiokou, M., Marqués, P. et al.
#' Pasta, an age-shift transcriptomic clock, maps the chemical and genetic determinants of aging and rejuvenation.
#' \emph{bioRxiv.} 2025
#'
#' @examples
#' # 1. Fast runnable example
#' print("Ready to calculate PASTA scores.")
#'
#' \donttest{
#' library(magrittr)
#' library(Seurat)
#' library(glmnet)
#'
#' seu <- loadOmniAgeRdata(
#'     "omniager_seu_gabitto_2024_filtered",
#'     verbose = FALSE
#' )
#'
#' seu$age <- seu$development_stage %>%
#'     gsub("-year.*", "", .) %>%
#'     gsub("-", " ", .) %>%
#'     gsub("80 year old and over stage", "85", .)
#'
#' seuBulk <- makePseudobulksPasta(
#'     seu,
#'     poolBy = c("cell_type", "age"),
#'     chunkSize = 512,
#'     verbose = FALSE
#' )
#'
#' lognormMatrix <- GetAssayData(seuBulk, assay = "RNA", layer = "data")
#' lognormMatrix <- as.matrix(lognormMatrix)
#'
#' seuBulkMeta <- seuBulk[[c("chunkSize", "cell_type", "age")]]
#' seuBulkMeta$age <- as.numeric(seuBulkMeta$age)
#' pastaRes <- pastaScores(lognormMatrix, filterGenes = TRUE, rankNorm = TRUE)
#' }
#'
pastaScores <- function(mat, filterGenes = TRUE, rankNorm = TRUE,
                        reg = TRUE, pasta = TRUE, ct46 = TRUE, verbose = TRUE) {
    # 1. Load model data
    pastaGenesModel <- loadOmniAgeRdata(
        "omniager_pasta_gene",
        verbose = verbose
    )

    resList <- list()

    # 2. Preprocessing
    if (filterGenes) {
        mat <- filterAgeModelGenes(mat, pastaGenesModel)
    }

    if (rankNorm) {
        mat <- applyRankNormalization(mat)
    }

    matT <- t(mat)

    # 3. Prediction
    if (reg) resList[["REG"]] <- .predictAgeScore(matT, modelType = "REG")
    if (pasta) resList[["PASTA"]] <- .predictAgeScore(matT, modelType = "PASTA")
    if (ct46) resList[["CT46"]] <- .predictAgeScore(matT, modelType = "CT46")

    return(resList)
}


#' Filter Age Model Genes from Count Matrix
#'
#' Subsets the count matrix to include only genes used in the age prediction model.
#'
#' @param mat Matrix. Count matrix.
#' @param pastaGenesModel The genes used for building the model.
#' @return Matrix. Filtered count matrix with median imputation.
#' @export
#' @examples
#' # 1. Create a mock count matrix (5 genes, 3 samples)
#' mock_mat <- matrix(1:15, nrow = 5, ncol = 3)
#' rownames(mock_mat) <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")
#' colnames(mock_mat) <- c("Sample1", "Sample2", "Sample3")
#'
#' # 2. Define the genes required by the model
#' # Note: "GeneF" is intentionally missing from the matrix to test imputation
#' model_genes <- c("GeneA", "GeneC", "GeneF")
#'
#' # 3. Run the filter function
#' filtered_mat <- filterAgeModelGenes(mat = mock_mat, pastaGenesModel = model_genes)
#'
#' # 4. View the result
#' print(filtered_mat)
filterAgeModelGenes <- function(mat, pastaGenesModel) {
    if (is.null(rownames(mat))) {
        stop("Matrix must have row names (gene symbols).")
    }

    idx <- match(pastaGenesModel, rownames(mat))
    matFiltered <- mat[idx, , drop = FALSE]

    medianVal <- stats::median(matFiltered, na.rm = TRUE)
    matFiltered[is.na(matFiltered)] <- medianVal

    rownames(matFiltered) <- pastaGenesModel
    return(matFiltered)
}

#' Apply Rank Normalization to Matrix
#'
#' Applies rank normalization across each column of the matrix.
#'
#' @param mat Matrix. Count matrix.
#' @return Matrix. Rank-normalized matrix.
#' @export
#' @examples
#' # 1. Create a mock count matrix (3 genes, 3 samples)
#' # Notice that Sample1 has a tie (10 and 10), and Sample3 has identical values
#' mock_mat <- matrix(
#'     c(
#'         10, 20, 10, # Sample 1
#'         5, 50, 15, # Sample 2
#'         8, 8, 8
#'     ), # Sample 3
#'     nrow = 3, ncol = 3
#' )
#' rownames(mock_mat) <- c("GeneA", "GeneB", "GeneC")
#' colnames(mock_mat) <- c("Sample1", "Sample2", "Sample3")
#'
#' # 2. View the original matrix
#' print(mock_mat)
#'
#' # 3. Apply rank normalization
#' norm_mat <- applyRankNormalization(mat = mock_mat)
#'
#' # 4. View the rank-normalized matrix
#' # In Sample1, the two 10s will share the average rank of 1.5
#' print(norm_mat)
applyRankNormalization <- function(mat) {
    matNormalized <- apply(mat, 2, function(x) {
        r <- rank(x, ties.method = "average")
        return(r)
    })
    return(matNormalized)
}

#' Predict Age Score from Gene Expression Matrix
#'
#' Uses pre-trained models to predict age scores based on gene expression.
#'
#' @param mat Matrix. Processed count matrix.
#' @param modelType Character. Model type ('PASTA', 'REG', or 'CT46').
#' @param verbose Logical. Whether to print status messages during data loading. Default is FALSE.
#' @return Numeric vector. Predicted age scores.
#' @noRd

.predictAgeScore <- function(mat, modelType = "PASTA", verbose = FALSE) {
    # 1. Load data according to the type
    dataName <- switch(modelType,
        "REG" = "omniager_cvfit_reg",
        "PASTA" = "omniager_cvfit_pasta",
        "CT46" = "omniager_cvfit_c46",
        stop("Invalid modelType. Choose PASTA, REG, or CT46.")
    )
    # Obtain the current model object
    curModel <- loadOmniAgeRdata(
        dataName,
        verbose = verbose
    )

    # 2. Prediction
    vAgeScores <- as.numeric(stats::predict(curModel, mat,
        s = "lambda.min",
        type = "link"
    )[, 1])

    # 3. Scaled
    if (modelType == "PASTA") {
        betaPASTA <- loadOmniAgeRdata(
            "omniager_beta_pasta",
            verbose = verbose
        )
        vAgeScores <- vAgeScores * betaPASTA
    } else if (modelType == "CT46") {
        betaC46 <- loadOmniAgeRdata(
            "omniager_beta_c46",
            verbose = verbose
        )
        vAgeScores <- vAgeScores * betaC46
    }

    return(vAgeScores)
}


#' Create Pseudobulk Samples from Seurat Object
#'
#' Aggregates single-cell expression data into pseudobulk samples based on
#' user-defined metadata variables and a specified chunk size.
#'
#' @param seu A Seurat object.
#' @param poolBy A character vector of column names from `@meta.data`. Cells are
#'   grouped by unique combinations of these variables prior to chunking.
#'   Default: `c("cell_type", "age")`.
#' @param chunkSize A numeric value. The target number of cells per pseudobulk
#'   sample. If set to 1, no aggregation is performed and the original
#'   object is returned (with metadata updated). Default: `1000`.
#' @param verbose Logical. If `TRUE`, prints a summary table detailing the
#'   number of pseudobulk samples generated per group. Default: `TRUE`.
#'
#' @return A new Seurat object where columns represent pseudobulk samples.
#'   The meta.data slot includes the original `poolBy` variables and the
#'   `chunkSize` used for aggregation.
#' @export
#' @examples
#' library(Seurat)
#'
#' # 1. Create a tiny mock count matrix (5 genes, 20 cells)
#' set.seed(123)
#' mock_counts <- matrix(rpois(100, lambda = 5), nrow = 5, ncol = 20)
#' rownames(mock_counts) <- paste0("Gene", 1:5)
#' colnames(mock_counts) <- paste0("Cell_", 1:20)
#'
#' # 2. Create mock metadata with the default 'poolBy' columns
#' mock_meta <- data.frame(
#'     cell_type = rep(c("T_cell", "B_cell"), each = 10),
#'     age = rep(c(30, 40), times = c(10, 10)),
#'     row.names = colnames(mock_counts)
#' )
#'
#' # 3. Build the Seurat object
#' seu_mock <- CreateSeuratObject(counts = mock_counts, meta.data = mock_meta)
#'
#' # 4. Run the pseudobulk function
#' # We have 10 cells per group. With chunkSize = 5, we expect exactly
#' # 2 pseudobulks for T_cell and 2 pseudobulks for B_cell.
#' seu_pb <- makePseudobulksPasta(
#'     seu = seu_mock,
#'     poolBy = c("cell_type", "age"),
#'     chunkSize = 5,
#'     verbose = FALSE
#' )
#'
#' # 5. View the metadata of the resulting pseudobulked object
#' print(seu_pb[[]])
makePseudobulksPasta <- function(seu, poolBy = c("cell_type", "age"),
                                 chunkSize = 1000, verbose = TRUE) {
    if (!all(poolBy %in% colnames(seu[[]]))) {
        stop("Required metadata columns not found in Seurat object.")
    }

    if (chunkSize == 1) {
        return(seu)
    }

    # 1. Generate grouping factors
    groupFactor <- interaction(seu[[]][, poolBy], drop = TRUE)

    # 2. Block indexing calculation
    indices <- seq_len(nrow(seu[[]]))
    chunkIds <- character(length(indices))

    splitIndices <- split(indices, groupFactor)

    for (groupName in names(splitIndices)) {
        idx <- splitIndices[[groupName]]
        n <- length(idx)
        numChunks <- ceiling(n / chunkSize)
        groupChunks <- rep(seq_len(numChunks), each = chunkSize, length.out = n)
        chunkIds[idx] <- paste(groupName, sample(groupChunks), sep = "-")
    }

    seu$tempChunkId <- make.names(chunkIds)

    # 3. Aggregation
    seuBulk <- Seurat::AggregateExpression(
        seu,
        group.by = "tempChunkId",
        return.seurat = TRUE,
        verbose = FALSE
    )

    # 4. Metadata recovery
    firstMatchIdx <- match(colnames(seuBulk), seu$tempChunkId)
    originalMeta <- seu[[]][firstMatchIdx, poolBy, drop = FALSE]

    for (col in poolBy) {
        seuBulk[[col]] <- originalMeta[[col]]
    }

    seuBulk$chunkSize <- chunkSize

    if (verbose) {
        msg <- utils::capture.output(table(seuBulk[[poolBy[1]]]))
        message(paste(msg, collapse = "\n"))
    }

    return(seuBulk)
}
