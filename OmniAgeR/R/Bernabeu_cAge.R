#' @title Calculate Chronological Age (cAge) using the hybrid model
#'
#' @description This function loads pre-trained cAge model weights and applies them
#'              to a user-provided methylation beta-value matrix.
#'              It implements the hybrid model logic from the paper:
#'              1. Predicts age using the standard linear model.
#'              2. If the predicted age is < 20 years, it re-predicts using the log(age) model.
#'              It also automatically handles missing CpGs and NA values via mean imputation.
#'
#' @param beta.m (matrix) A numeric matrix of beta values.
#'               ASSUMPTION: Rows must be CpGs (named with CpG IDs).
#'               ASSUMPTION: Columns must be samples (named with Sample IDs).
#'
#' @return A named numeric vector of predicted ages. The names of the vector
#'         correspond to the sample IDs (column names) from `beta.m`.
#'
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Bernabeu_cAge.o <- Bernabeu_cAge(hannum_bmiq_m)



Bernabeu_cAge <- function(beta.m) {

  ###### (1) Load Model Data
  data("Bernabeu_cAge_Coef")


  # Separate coefficients for log model (linear and quadratic)
  coef_log_2 <- Bernabeu_cAge_coef_log[grep('_2', rownames(Bernabeu_cAge_coef_log)), , drop = FALSE]
  coef_log_2_simp <- gsub('_2', '', rownames(coef_log_2))
  coef_log <- Bernabeu_cAge_coef_log[!(rownames(Bernabeu_cAge_coef_log) %in% rownames(coef_log_2)), , drop = FALSE]

  # Separate coefficients for non-log model (linear and quadratic)
  coef_2 <- Bernabeu_cAge_coef[grep('_2', rownames(Bernabeu_cAge_coef)), , drop = FALSE]
  coef_2_simp <- gsub('_2', '', rownames(coef_2))
  coef <- Bernabeu_cAge_coef[!(rownames(Bernabeu_cAge_coef) %in% rownames(coef_2)), , drop = FALSE]

  # Get the total list of required CpGs
  cpgs_linear <- union(rownames(coef), rownames(coef_log))
  cpgs_squared <- union(coef_2_simp, coef_log_2_simp)
  all_cpgs <- union(cpgs_linear, cpgs_squared)

  # Use the input matrix directly
  data <- beta.m

  print(paste0(
    "[Bernabeu_cAge] Number of represented Bernabeu_cAge CpGs (max=",
    length(all_cpgs), ")=", sum(all_cpgs %in% rownames(beta.m))))

  ###### (2) QC and Data Prep
  #########################################################################################################
  #message("2. Quality Control and data preparation")

  # --- User provides beta-value matrix with rows=CpGs, cols=Samples ---

  ## 2.1 Subset CpG sites
  #message("2.1 Subsetting CpG sites to those required for predictor calculation")
  coef_data <- data[intersect(rownames(data), all_cpgs), , drop = FALSE]

  ## 2.2 Impute completely missing CpGs
  #message("2.2 Finding CpGs not present in input matrix...")
  if (nrow(coef_data) == length(all_cpgs)) {
    #message("All required sites are present.")
  } else if (nrow(coef_data) == 0) {
    stop("[Bernabeu_cAge] Error: None of the required CpGs are in the dataset. Analysis cannot proceed.")
  } else {
    missing_cpgs <- all_cpgs[!(all_cpgs %in% rownames(coef_data))]
    message(paste("[Bernabeu_cAge]",length(missing_cpgs), "unique sites are missing - adding with mean Beta Value from GS training sample (N = 18,413)"))

    mat <- matrix(nrow = length(missing_cpgs), ncol = ncol(coef_data))
    rownames(mat) <- missing_cpgs
    colnames(mat) <- colnames(coef_data)

    # Fill with training set means
    missing_cpg_means <- Bernabeu_cAge_means[missing_cpgs, "mean"]
    # Use sweep for efficient matrix filling by row
    mat <- sweep(mat, 1, missing_cpg_means, `*`)

    coef_data <- rbind(coef_data, mat)
  }

  ## 2.3 Impute NA Values (within-sample NAs)
  #message("2.3 Converting NA Values to the mean for each probe (using sample means)")

  # Define an internal function to replace NAs with the row mean
  na_to_mean <- function(methyl_row) {
    row_mean <- mean(methyl_row, na.rm = TRUE)

    # Handle rows that are all NA
    if (is.nan(row_mean)) {
      cpg_name <- names(methyl_row)[1] # Get CpG name from the (named) vector
      if (cpg_name %in% rownames(Bernabeu_cAge_means)) {
        # Fallback to training set mean
        row_mean <- Bernabeu_cAge_means[cpg_name, "mean"]
      } else {
        # Should not happen if all_cpgs is correct, but as a failsafe
        warning(paste("[Bernabeu_cAge]: Could not find mean for all-NA probe:", cpg_name, ". Imputing with 0.5"))
        row_mean <- 0.5
      }
    }

    methyl_row[is.na(methyl_row)] <- row_mean
    return(methyl_row)
  }

  # Apply na_to_mean row-wise and transpose back
  # We need to keep track of names for the all-NA fallback
  original_cpg_names <- rownames(coef_data)

  coef_data_fixed_list <- lapply(original_cpg_names, function(cpg_name) {
    row_data <- coef_data[cpg_name, ]
    names(row_data) <- cpg_name # Pass name to na_to_mean
    na_to_mean(row_data)
  })

  coef_data <- do.call(rbind, coef_data_fixed_list)

  rownames(coef_data) <- original_cpg_names
  # Ensure column order is preserved
  colnames(coef_data) <- colnames(beta.m)


  ###### (3) cAge Prediction
  #########################################################################################################
  #message("3. Obtaining cAge predictions")
  #message("3.1. Preparing data for prediction")

  # Ensure coef_data rows are in the exact order required by the model lists
  coef_data <- coef_data[all_cpgs, , drop = FALSE]

  ## Prepare input for the linear (age) model
  scores <- coef_data[rownames(coef), , drop = FALSE]
  scores_quadratic <- (coef_data[coef_2_simp, , drop = FALSE])^2
  rownames(scores_quadratic) <- paste0(rownames(scores_quadratic), "_2")
  scores_linear <- rbind(scores, scores_quadratic)

  ## Prepare input for the log(age) model
  scores_log <- coef_data[rownames(coef_log), , drop = FALSE]
  scores_quadratic_log <- (coef_data[coef_log_2_simp, , drop = FALSE])^2
  rownames(scores_quadratic_log) <- paste0(rownames(scores_quadratic_log), "_2")
  scores_log <- rbind(scores_log, scores_quadratic_log)

  ## 3.2. Calculate cAge using the standard (age) model
  #message("3.2. Calculating cAge using model trained on linear age")

  # Ensure weights and data rows are in the exact same order
  coefficients_ordered <- Bernabeu_cAge_coef[rownames(scores_linear), "Coefficient"]
  pred_linear <- scores_linear * coefficients_ordered
  pred_linear_pp <- colSums(pred_linear, na.rm = TRUE)
  pred_linear_pp <- pred_linear_pp + Bernabeu_cAge_intercept

  ## 3.3. Apply hybrid model logic
  #message("3.3. Applying hybrid logic for individuals < 20 years")

  # Identify individuals < 20 years
  under20s_mask <- pred_linear_pp < 20
  under20s_mask[is.na(under20s_mask)] <- FALSE # Treat NA predictions as >= 20

  under20s_samples <- names(pred_linear_pp)[under20s_mask]

  # Initialize the final predictions vector
  final_predictions <- pred_linear_pp

  # If there are any under 20s, recalculate them
  if (length(under20s_samples) > 0) {
    #message(paste("Re-running prediction for", length(under20s_samples), "sample(s) using log(age) model"))

    # Ensure weights and data rows are in the exact same order
    coefficients_log_ordered <- Bernabeu_cAge_coef_log[rownames(scores_log), "Coefficient"]

    # Subset the log-model data just for the under-20s samples
    scores_log_subset <- scores_log[, under20s_samples, drop = FALSE]

    pred_log <- scores_log_subset * coefficients_log_ordered
    pred_log_pp <- colSums(pred_log, na.rm = TRUE)
    pred_log_pp <- pred_log_pp + Bernabeu_cAge_intercept_log

    # Convert from log(age) back to age and update the final predictions vector
    final_predictions[under20s_samples] <- exp(pred_log_pp)
  } else {
    #message("No samples predicted < 20 years. Log(age) model not used.")
  }

  ###### (4) Finalize and Return Results
  #########################################################################################################
  #message("4. Finalizing predictions.")

  # Ensure the final vector has the same order as the input matrix columns
  final_predictions <- final_predictions[colnames(beta.m)]

  return(final_predictions)
}

