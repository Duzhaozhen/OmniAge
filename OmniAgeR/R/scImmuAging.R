#' @title
#' Run the scImmuAging prediction pipeline for multiple cell types
#'
#' @description
#' This function takes as input a log-normalized matrix of scRNA-seq data and
#' will return cell type-specific predicted age.
#'
#' @param seuratObj A Seurat object. The meta.data must contain "donor_id",
#' "age" and "celltype" columns.
#' @param cellTypes A character vector specifying the cell types for which to
#' predict age.Example: `c("CD4T", "CD8T", "MONO")`. Valid types are "CD4T",
#' "CD8T", "MONO", "NK", "B".
#' @param minCoverage Numeric (0-1). Minimum required feature coverage.
#' Default 0.
#' @param verbose Logical. Whether to print status messages.
#'
#' @details
#' This function takes a Seurat object, preprocesses the scRNA-seq data for
#' one or more specified cell types, and predicts age using the pre-trained
#' scImmuAging models. It iterates through a vector of cell types, performs
#' the full analysis for each, and returns a nested list containing all results.
#'
#' @return A list where each element is named by a cell type
#' from the input vector.
#' Each of these elements is itself a list containing two data frames:
#' \describe{
#'   \item{BootstrapCell}{A data frame with age predictions for each
#'   bootstrapped pseudocell,
#'   containing the columns `donorId`, `age`, and `Prediction`.}
#'   \item{Donor}{A data frame with the final aggregated age prediction
#'   for each donor,
#'   containing the columns `donorId`, `age`, and `predicted`.}
#' }
#'
#' @references
#' Li W, Zhang Z, Kumar S, et al.
#' Single-cell immune aging clocks reveal inter-individual heterogeneity during
#' infection and vaccination.
#' \emph{Nat Aging} 2025
#'
#' @export
#'
#' @examples
#' # 1. Define valid cell types (Runnable code to satisfy BiocCheck)
#' valid_cell_types <- c("CD4T", "CD8T", "MONO", "NK", "B")
#' print(valid_cell_types)
#'
#' \donttest{
#' # 2. Real pipeline execution (Wrapped in donttest because it requires
#' # downloading pre-trained models and external example datasets)
#' library(Seurat)
#' seuratObj <- loadOmniAgeRdata(
#'     "omniager_yazar_cd4t_cd8t_example",
#'     verbose = FALSE
#' )
#'
#' scImmuAgingOut <- scImmuAging(seuratObj, c("CD4T", "CD8T"))
#' }
scImmuAging <- function(seuratObj, cellTypes, minCoverage = 0, verbose = TRUE) {
    # 1. loading model

    scimmuagingModel <- loadOmniAgeRdata(
        "omniager_scimmuaging_model",
        verbose = verbose
    )

    allResults <- list()

    for (ct in cellTypes) {
        if (verbose) message("\n--- Processing cell type: ", ct, " ---")

        # Validate model
        if (!ct %in% names(scimmuagingModel$model_set)) {
            warning("Cell type '", ct, "' model not found. Skipping.")
            next
        }

        currentModel <- scimmuagingModel$model_set[[ct]]
        currentMarkerGenes <- scimmuagingModel$feature_set[[ct]]

        # Validate Metadata
        if (!"celltype" %in% colnames(seuratObj[[]])) {
            stop("Metadata must contain a 'celltype' column.")
        }

        # 3. Subset extraction
        seuratSub <- subset(seuratObj, subset = celltype == ct)

        # 4. Preprocessing
        preprocessedData <- scImmuAgingPreProcess(seuratSub, ct, currentModel, currentMarkerGenes)

        # 5. Predcition
        predResults <- scImmuAgingCalculator(
            preprocessedData = preprocessedData,
            model = currentModel,
            markerGenes = currentMarkerGenes,
            minCoverage = minCoverage,
            verbose = verbose
        )

        # 6. Summary
        agePerDonor <- ageDonor(predResults)

        allResults[[ct]] <- list(
            bootstrapCell = predResults,
            donor = agePerDonor
        )
    }

    return(allResults)
}
