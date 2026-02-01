#' Run Brain Cell Type Specific Clocks
#'
#' @description
#' A high-level wrapper function that runs the brain age prediction pipeline for
#' one or more specified sample types ('SC', 'Pseudobulk', 'Bootstrap').
#'
#' @param seurat_object The input Seurat object containing expression data and
#'   metadata (must include 'donor_id', 'age', 'celltype').
#' @param cell_types A character vector of cell types to analyze
#'   (e.g., `c('Oligodendrocytes', 'Astrocytes')`).
#'   Available cell types are: "Oligodendrocytes", "Astrocytes", "Microglia",
#'   "OPCs", "Excitatory Neurons", and "Inhibitory Neurons".
#' @param model_name A character string or vector specifying which models to run.
#'   \itemize{
#'     \item `"all"` (default): Runs "SC", "Pseudobulk", and "Bootstrap".
#'     \item A vector: e.g., `c("SC", "Pseudobulk")` will run only those two.
#'   }
#'
#' @return
#' A named list where each element corresponds to a `model_name` that was run.
#' Each element contains a data.frame of the (5-fold averaged) predictions,
#' as returned by `run_prediction_pipeline`.
#'
#' @details
#' This function serves as the primary endpoint for the prediction pipeline.
#' It iteratively calls `run_prediction_pipeline` for each requested `model_name`
#' (sample type) and collects the results into a single list.
#'
#' @seealso \code{\link{run_prediction_pipeline}}
#'
#' @references
#' Muralidharan C, Zakar-Polyák E, Adami A, et al.
#' Human Brain Cell-Type-Specific Aging Clocks Based on Single-Nuclei Transcriptomics.
#' \emph{Adv Sci(Weinh).} 2025
#'
#' @export
#' @examples
#' # Load the Seurat object
#' library(Seurat)
#' library(dplyr)
#' download_OmniAgeR_example("brain_frohlich_control_example_15donors")
#' load_OmniAgeR_example("brain_frohlich_control_example_15donors")
#'
#' # Define cell types of interest
#' types_to_run <- c("Oligodendrocytes")
#'
#' # Run all three models for the specified cell types
#' all_clock_results <- Brain_CT_clock(seurat_object = brain_seurat,cell_types = types_to_run, model_name = "SC")
#'

Brain_CT_clock <- function(seurat_object, cell_types, model_name = "all") {

  # --- Define which models to run ---
  if (length(model_name) == 1 && model_name == "all") {
    models_to_run <- c("SC", "Pseudobulk", "Bootstrap")
  } else {
    models_to_run <- model_name

    # Validate input
    valid_models <- c("SC", "Pseudobulk", "Bootstrap")
    if (!all(models_to_run %in% valid_models)) {
      invalid <- models_to_run[!models_to_run %in% valid_models]
      stop("Invalid model_name(s) provided: ", paste(invalid, collapse=", "),
           ". Valid names are 'SC', 'Pseudobulk', 'Bootstrap', or 'all'.")
    }
  }

  # --- Initialize a list to store all results ---
  all_results <- list()

  cat(paste("Starting Brain_CT_clock for cell types\n"))

  # --- Loop over each requested model type ---
  for (current_model_type in models_to_run) {

    # Call the main pipeline function
    pred_result <- run_prediction_pipeline(
      sample_type = current_model_type,
      seurat_obj = seurat_object,
      common_celltypes = cell_types
    )

    # Store the result in the list with the correct name
    all_results[[current_model_type]] <- pred_result
  }

  cat("\n--- All Brain_CT_clock predictions complete! ---\n")
  return(all_results)
}







#' Extract and Pre-process Data from a Seurat Object
#'
#' @description
#' Subsets a Seurat object by cell type and processes the expression data into one
#' of three formats: single-cell ('SC'), 'Pseudobulk' (averaged by donor), or
#' 'Bootstrap' (resampled within donor).
#'
#' @param seurat_obj A Seurat object containing the expression data and metadata
#'   (must include 'donor_id', 'age', 'celltype').
#' @param cell_type A character string specifying the cell type to subset.
#' @param sample_type A character string defining the processing method:
#'   "SC", "Pseudobulk", or "Bootstrap".
#'
#' @return
#' A data.frame (tibble) containing the processed expression data and metadata.
#' The structure depends on `sample_type`:
#' \itemize{
#'   \item{"SC": A cell-by-gene matrix with metadata.}
#'   \item{"Pseudobulk": A donor-by-gene matrix (mean expression) with metadata.}
#'   \item{"Bootstrap": A (donor * 100 replicates)-by-gene matrix with metadata.}
#' }
#' Returns an empty data.frame if no cells are found for the specified `cell_type`.
#'
#' @details
#' This function defines the number of cells to sample for the 'Bootstrap' method
#' internally via the `num_cells` list. If a donor has fewer cells than the
#' specified number, 100 replicates of the donor's mean expression are returned.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom dplyr as_tibble group_by summarize across all_of select everything bind_rows
#'
#' @export
get_df_seurat <- function(seurat_obj, cell_type, sample_type) {

  # Define cell counts for bootstrap sampling
  num_cells <- list(
    "Oligodendrocytes" = 200,
    "Astrocytes" = 50,
    "Microglia" = 50,
    "OPCs" = 50,
    "Excitatory Neurons" = 100,
    "Inhibitory Neurons" = 100
  )


  # 1. Subset Seurat Object
  k <- subset(seurat_obj, subset = celltype == cell_type)


  if (ncol(k) == 0) {
    warning(paste("No cells found for cell_type:", cell_type))
    # Return an empty data.frame so downstream steps can handle it safely
    return(data.frame())
  }

  # 2. Extract metadata and expression matrix
  k_meta <- k@meta.data
  k_meta <- k_meta[,c(c("donor_id", "age","celltype"))]
  input_mtx <- t(as.matrix(GetAssayData(k, assay = "RNA", layer = "data")))
  df <- dplyr::as_tibble(cbind(k_meta, input_mtx))

  # 4. Get all gene column names
  cols_to_keep <- colnames(df)[4:ncol(df)]

  # 5. Process data based on sample_type
  if (sample_type == "SC") {
    return(df)

  } else if (sample_type == "Pseudobulk") {
    df_ct <- df %>%
      group_by(donor_id, age, celltype) %>%
      summarize(
        across(all_of(cols_to_keep), .fns = mean),
        .groups = 'drop'
      )
    return(df_ct)

  } else if (sample_type == "Bootstrap") {
    donors <- unique(df$donor_id)
    n_cells_sample <- num_cells[[cell_type]]

    df_avg_list <- lapply(donors, function(donor) {
      df_aux <- df[df$donor_id == donor, ]
      set.seed(42) # Ensure reproducibility
      n_rows_aux <- nrow(df_aux)
      row_indices <- 1:n_rows_aux

      if (n_rows_aux > n_cells_sample) {
        # Sample N cells, 100 times, without replacement
        boot_list <- replicate(100, {
          sampled_indices <- sample(row_indices, size = n_cells_sample, replace = FALSE)
          df_sample_means_vec <- colMeans(df_aux[sampled_indices, cols_to_keep, drop = FALSE])
        }, simplify = FALSE)
        df_boot <- bind_rows(boot_list)
      } else {
        # If too few cells, just replicate the mean 100 times
        df_sample_means_vec <- colMeans(df_aux[, cols_to_keep, drop = FALSE])
        df_boot <- bind_rows(replicate(100, df_sample_means_vec, simplify = FALSE))
      }

      # Add metadata back
      df_boot$donor_id <- donor
      df_boot$age <- df_aux$age[1]
      df_boot$celltype <- df_aux$celltype[1]
      df_boot <- dplyr::select(df_boot, donor_id, age, celltype, dplyr::everything())
      return(df_boot)
    })

    df_avg <- bind_rows(df_avg_list)
    return(df_avg)
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
#' @param data A data.frame of expression data (samples in rows, genes in
#'   columns). Must also contain 'age' and 'donor_id' columns.
#' @param impute_data A long-format data.frame with 'feature_name' and
#'   'imputation_value' columns. Used to fill in genes present in the model
#'   but missing from `data`.
#' @param model A long-format data.frame with 'feature_name' (genes + 'intercept')
#'   and 'coefficient' columns, representing one trained clock.
#' @param sample_type A character string (e.g., 'SC', 'Pseudobulk') used to
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
#'
#' @keywords internal
predict_brain_ct_age <- function(data, impute_data, model, sample_type) {

  # 1. Prepare model (separate coefficients and intercept)
  model_intercept_df <- model[model$feature_name == "intercept", ]
  model_intercept <- if(nrow(model_intercept_df) > 0) model_intercept_df$coefficient[1] else 0

  model_genes_df <- model[model$feature_name != "intercept", ]
  model_genes <- model_genes_df$feature_name

  # If model has no genes (intercept-only), return intercept for all samples
  if (length(model_genes) == 0) {
    return(data.frame(
      predictions = rep(model_intercept, nrow(data)),
      ages = data$age,
      donors = data$donor_id,
      sample_type = sample_type
    ))
  }

  # 2. Gene matching
  genes_in_model <- intersect(model_genes, colnames(data))
  genes_to_impute <- setdiff(model_genes, genes_in_model)

  if(length(genes_to_impute) > 0) {
    # 3. Prepare imputation data
    data_for_imputation <- impute_data %>%
      dplyr::filter(feature_name %in% genes_to_impute) %>%
      tidyr::pivot_wider(
        names_from = feature_name,
        values_from = imputation_value
      )

    # Check if all required imputation genes were found in the imputation table
    missing_from_impute <- setdiff(genes_to_impute, colnames(data_for_imputation))
    if(length(missing_from_impute) > 0) {
      # If imputation table is incomplete, fill with 0
      print(paste("Warning: ", length(missing_from_impute), "genes missing from imputation file. Using 0."))
      for(g in missing_from_impute) { data_for_imputation[[g]] <- 0 }
    }

    # 4. Create the full expression matrix for prediction
    expr_data <- data[, genes_in_model, drop = FALSE]
    # Replicate imputation values for all samples
    impute_matrix <- data_for_imputation[rep(1, nrow(expr_data)), , drop = FALSE]
    full_expr_matrix <- bind_cols(expr_data, impute_matrix)

  } else {
    # No imputation needed
    expr_data <- data[, genes_in_model, drop = FALSE]
    full_expr_matrix <- expr_data
  }

  # (d) **CRITICAL**: Ensure the matrix column order exactly matches
  #     the model's gene order.
  full_expr_matrix <- full_expr_matrix[, model_genes, drop = FALSE]

  # 5. Prepare coefficient vector (also in model_genes order)
  coeff_df_ordered <- data.frame(feature_name = model_genes) %>%
    dplyr::left_join(model_genes_df, by = "feature_name")

  coeff_vector <- coeff_df_ordered$coefficient

  # 6. **Vectorized Prediction** (Matrix Multiplication)
  #    (N_samples x N_genes) %*% (N_genes x 1) -> (N_samples x 1)
  predictions_vec <- (as.matrix(full_expr_matrix) %*% coeff_vector) + model_intercept

  return(data.frame(
    predictions = predictions_vec[,1], # Convert from 1-col matrix to vector
    ages = data$age,
    donors = data$donor_id,
    sample_type = sample_type
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
#' @param sample_type The processing method to use: "SC", "Pseudobulk", or
#'   "Bootstrap". This is passed to `get_df_seurat`.
#' @param seurat_obj The input Seurat object, passed to `get_df_seurat`.
#' @param common_celltypes A character vector of cell types to process
#'   (e.g., `c('Oligodendrocytes', 'Astrocytes')`).
#'
#' @return
#' A single, combined data.frame containing the final (5-fold averaged)
#' 'predictions', 'ages', 'donors', 'sample_type', and 'celltype' for all
#' requested cell types.
#'
#' @details
#' This is the main user-facing function for this pipeline. It assumes
#' that a file named "brain_celltypeSpecific_clocks_coef_imputation_5fold.rda"
#' exists in the specified path, containing two objects:
#' \itemize{
#'   \item `brain_ct_clocks_coef`: A nested list where keys
#'     (e.g., "SC_Oligodendrocytes") map to a list of 5 models (one for each fold).
#'   \item `all_imputation_data_list`: A list where keys
#'     (e.g., "SC_Oligodendrocytes") map to a single imputation data.frame.
#' }
#' The function calls `get_df_seurat` to prepare the data, then calls
#' `predict_brain_ct_age` five times (once for each fold). The five resulting
#' prediction vectors are then averaged (row-wise) to produce the final,stable prediction.
#'
#' @importFrom dplyr bind_rows
#' @export
run_prediction_pipeline <- function(sample_type, seurat_obj,
                                    common_celltypes){

  processing_message <- paste0(sample_type," Model")

  # --- Load Models and Imputation Data ---
  data("brain_celltypeSpecific_clocks_coef")
  # This loads `brain_ct_clocks_coef` and `all_imputation_data_list`

  print(paste("--- Starting", processing_message, "Pipeline ---"))

  # Create a list to collect final (averaged) results for each cell type
  final_results_list <- list()

  for (cell_type in common_celltypes) {

    clock_key <- paste(sample_type, cell_type, sep="_")

    ## 1. Load corresponding imputation data and coefficient data
    # model_folds_list is now a LIST of 5 models
    model_folds_list <- brain_ct_clocks_coef[[ clock_key ]]
    # imputation_data is still a single data.frame
    imputation_data <- all_imputation_data_list[[ clock_key ]]

    # Check if data exists
    if (is.null(model_folds_list) || is.null(imputation_data)) {
      warning(paste("Skipping", cell_type, ": Model or imputation data not found. Key:", clock_key))
      next
    }

    print(paste("Processing", processing_message, "for:", cell_type))

    # 2. Get the base data from the Seurat object
    df_base <- get_df_seurat(seurat_obj, cell_type, sample_type)

    if (nrow(df_base) == 0) {
      warning(paste("Skipping", cell_type, ": get_df_seurat returned 0 cells/samples."))
      next
    }

    # 3. **NEW**: Loop over the 5-fold models to get predictions
    all_fold_predictions_list <- list()

    for (fold_name in names(model_folds_list)) {
      model_single_fold <- model_folds_list[[fold_name]]

      # 6. Run prediction (using a single fold's model)
      pred_fold_res <- predict_brain_ct_age(
        df_base,
        imputation_data,
        model_single_fold, # Pass the single fold model
        sample_type
      )

      # Store only the prediction vector
      all_fold_predictions_list[[fold_name]] <- pred_fold_res$predictions
    }

    # 4. **NEW**: Average the predictions across the 5 folds
    # Convert list to a data.frame (N_samples x 5_folds)
    predictions_folds_df <- as.data.frame(all_fold_predictions_list)

    # Calculate row-wise means (N_samples x 1)
    final_avg_predictions <- rowMeans(predictions_folds_df, na.rm = TRUE)

    # 5. Build the final results data.frame with the *averaged* predictions
    pred_res_avg <- data.frame(
      predictions = final_avg_predictions,
      ages = df_base$age,
      donors = df_base$donor_id,
      sample_type = sample_type,
      celltype = cell_type
    )

    # Add results to the list (faster than rbind in a loop)
    final_results_list[[cell_type]] <- pred_res_avg
  }

  # After the loop, combine all data.frames at once
  final_pred_res <- dplyr::bind_rows(final_results_list)

  print(paste("--- Finished", processing_message, "Pipeline ---"))
  return(final_pred_res)
}

