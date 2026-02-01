#' @title
#' Run the scImmuAging prediction pipeline for multiple cell types
#'
#' @description
#' This function takes as input a log-normalized matrix of scRNA-seq data and will return cell type-specific predicted age.
#'
#' @param seurat_object A Seurat object. The meta.data must contain "donor_id", "age" and "celltype" columns.
#'
#' @param cell_types A character vector specifying the cell types for which to predict age.Example: `c("CD4T", "CD8T", "MONO")`. Valid types are "CD4T", "CD8T", "MONO", "NK", "B".
#'
#' @details
#' This function takes a Seurat object, preprocesses the scRNA-seq data for one or more specified cell types, and predicts age using the pre-trained scImmuAging models. It iterates through a vector of cell types, performs the full analysis for each, and returns a nested list containing all results.
#'
#' @return A list where each element is named by a cell type from the input vector.
#' Each of these elements is itself a list containing two data frames:
#' \describe{
#'   \item{BootstrapCell}{A data frame with age predictions for each bootstrapped pseudocell,
#'   containing the columns `donor_id`, `age`, and `Prediction`.}
#'   \item{Donor}{A data frame with the final aggregated age prediction for each donor,
#'   containing the columns `donor_id`, `age`, and `predicted`.}
#' }
#'
#' @references
#' Li W, Zhang Z, Kumar S, et al.
#' Single-cell immune aging clocks reveal inter-individual heterogeneity during infection and vaccination.
#' \emph{Nat Aging} 2025
#'
#' @export
#'
#' @examples
#' library(Seurat)
#' library(glmnet)
#' download_OmniAgeR_example("Yazar_CD4T_CD8T_example")
#' load_OmniAgeR_example("Yazar_CD4T_CD8T_example")
#' scImmuAging.out <- scImmuAging(seurat_obj,c("CD4T","CD8T"))
#'


scImmuAging <- function(seurat_object, cell_types) {
  # --- Load models and features once ---
  data("scImmuAgingFeature")
  cat("Using scImmuAging to predict age\n")

  # --- Initialize a list to store all results ---
  all_results <- list()

  # --- Loop over each requested cell type ---
  for (ct in cell_types) {
    cat(paste0("\n--- Processing cell type: ", ct, " ---\n"))

    # --- 1. Input Validation ---
    if (!ct %in% names(model_set) || !ct %in% names(feature_set)) {
      warning(paste0("Skipping: cell type '", ct, "' not found in loaded models/features."))
      next # Skip to the next cell type
    }

    # --- 2. Select Model and Genes ---
    cat(paste("Selecting model and features for:", ct, "\n"))
    model <- model_set[[ct]]
    marker_gene <- feature_set[[ct]]

    # --- 3. Pre-process the data ---
    cat("Pre-processing data and generating pseudocells...\n")
    seurat_object_sub <- subset(seurat_object,subset = celltype == ct )
    preprocessed_df <- scImmuAging_PreProcess(seurat_object_sub, ct, model, marker_gene)

    # --- 4. Predict Age ---
    cat("Calculating age predictions...\n")
    predict_res <- scImmuAgingCalculator(preprocessed_df, model, marker_gene)

    # --- 5. Aggregate results per donor ---
    cat("Aggregating predictions for each donor...\n")
    age_per_donor <- Age_Donor(predict_res)

    cat(paste("Prediction complete for", ct, "!\n"))

    # --- 6. Store results for the current cell type ---
    res_list <- list(
      BootstrapCell = as.data.frame(predict_res),
      Donor = as.data.frame(age_per_donor)
    )

    # Add the named list for the current cell type to the main results list
    all_results[[ct]] <- res_list
  }

  cat("\n--- All predictions complete! ---\n")
  return(all_results)
}
