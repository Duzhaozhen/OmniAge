#' @title Calculate DNAmFitAge and related fitness biomarkers
#'
#' @description
#' Calculates all 6 DNAm fitness biomarkers, DNAmFitAge, and FitAgeAcceleration
#' from a DNA methylation matrix and phenotype data.
#'
#' This function serves as a wrapper that:
#' 1. Pre-processes the data and imputes missing CpGs.
#' 2. Calls `internal_DNAmFitnessEstimators` to calculate 6 fitness markers.
#' 3. Calls `internal_FitAgeEstimator` to calculate DNAmFitAge and its acceleration.
#'
#' @param beta_matrix A numeric matrix. **Rows must be CpGs, Columns must be Samples.**
#'   The `colnames` must be the sample IDs.
#' @param age_vector A numeric vector of chronological age for each sample.
#'   The order **must match the column order** of `beta_matrix`.
#' @param sex_vector A character or factor vector of sex for each sample
#'   (e.g., "Male" or "Female"). The order **must match the column order**
#'   of `beta_matrix`.
#' @param grimage_vector A numeric vector of pre-calculated DNAmGrimAge values.
#'   The order **must match the column order** of `beta_matrix`.
#' @param package_name A string with the name of your R package
#'   (e.g., "OmniAgeR"). This is used to locate the internal model data file.
#' @param verbose Logical. If `TRUE` (default), the function will print
#'   progress messages and QC information to the console.
#'
#' @return
#' A `data.frame` with 12 columns:
#' \itemize{
#'   \item `SampleID`: The sample identifiers.
#'   \item `Age`: The input chronological age.
#'   \item `Female`: The input sex, coded as 1 for Female, 0 for Male.
#'   \item `DNAmGait_noAge`, `DNAmGrip_noAge`, `DNAmVO2max`: Fitness biomarkers.
#'   \item `DNAmGait_wAge`, `DNAmGrip_wAge`, `DNAmFEV1_wAge`: Age-adjusted fitness biomarkers.
#'   \item `DNAmGrimAge`: The input DNAmGrimAge.
#'   \item `DNAmFitAge`: The calculated biological fitness age.
#'   \item `FitAgeAccel`: The fitness age acceleration (residual of DNAmFitAge regressed on Age).
#' }
#'
#' @export
#'
#' @references
#' McGreevy KM, Radak Z, Torma F, et al.
#' DNAmFitAge: biological age indicator incorporating physical fitness.
#' \emph{Aging} 2023
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' age <- PhenoTypesHannum_lv$Age
#' sex <- ifelse(PhenoTypesHannum_lv$Sex=="F","Female","Male")
#' GrimAge1.o <- GrimAge1(hannum_bmiq_m,age,sex)
#'
#' DNAmFitAge.o <- DNAmFitAge(hannum_bmiq_m,age,sex,GrimAge1.o$DNAmGrimAge1 )




DNAmFitAge <- function(beta_matrix, age_vector, sex_vector,grimage_vector, verbose = TRUE) {

  if (verbose) {
    print(paste0("[DNAmFitAge] Starting DNAmFitAge calculation..."))
  }

  # --- 0. Load Internal Model Data ---
  data("DNAmFitAgeCoef")

  # --- 1. Validate Inputs & Construct Initial Data Frame ---

  # *** KEY CHANGE HERE ***
  n_samples <- ncol(beta_matrix)

  if (length(age_vector) != n_samples) {
    stop(paste0("[DNAmFitAge] ERROR: Length of 'age_vector' (", length(age_vector),
                ") does not match columns in 'beta_matrix' (", n_samples, ")."))
  }
  if (length(sex_vector) != n_samples) {
    stop(paste0("[DNAmFitAge] ERROR: Length of 'sex_vector' (", length(sex_vector),
                ") does not match columns in 'beta_matrix' (", n_samples, ")."))
  }

  # *** KEY CHANGE HERE ***
  sample_ids <- colnames(beta_matrix)
  if (is.null(sample_ids)) {
    stop("[DNAmFitAge] ERROR: 'beta_matrix' must have colnames specifying Sample IDs.")
  }

  # Convert sex to Female (0/1)
  female_numeric <- ifelse(sex_vector == "Female", 1, 0)
  if (!all(female_numeric %in% c(0, 1))) {
    warning("[DNAmFitAge] 'sex_vector' contains values other than 'Male' or 'Female'. Assuming non-'F' values are Male (0).")
  }

  # Define internal ID column name
  id_col <- "Internal_SampleID"

  # Build base data frame (without CpGs)
  base_data <- data.frame(
    Internal_SampleID = sample_ids,
    Age = age_vector,
    Female = female_numeric,
    stringsAsFactors = FALSE
  )

  # --- 2. Prepare Data (Add CpGs and Pre-process) ---
  if (verbose) {
    print(paste0("[DNAmFitAge] Pre-processing data and imputing missing CpGs..."))
  }

  # *** KEY CHANGE HERE ***
  # Transpose the beta matrix so Rows=Samples, Cols=CpGs, which is
  # what the internal helper functions expect.
  beta_matrix_transposed <- t(beta_matrix)

  # Combine base data with the transposed beta matrix
  data_for_prep <- cbind(base_data, beta_matrix_transposed)

  # Call the helper function (which you must also add to your package)
  # This function expects 'DNAmFitnessModels' to be passed to it.
  data_prep <- internal_data_prep(dataset = data_for_prep,
                                  idvariable = id_col,
                                  DNAmFitnessModels = DNAmFitnessModels,
                                  verbose)


  # --- 3. Estimate Fitness Biomarkers ---
  if (verbose) {
    print(paste0("[DNAmFitAge] Estimating DNAm fitness biomarkers..."))
  }

  data_FitnessEst <- internal_DNAmFitnessEstimators(
    data_prep,
    IDvar = id_col,
    DNAmFitnessModels = DNAmFitnessModels
  )


  # Create a separate GrimAge data frame for merging later
  grimage_df <- data.frame(
    Internal_SampleID =  colnames(beta_matrix),
    DNAmGrimAge = grimage_vector,
    stringsAsFactors = FALSE
  )

  # --- 4. Merge DNAmGrimAge ---
  data_FitAge_prep <- merge(grimage_df,
                            data_FitnessEst,
                            by = id_col)




  # --- 5. Calculate DNAmFitAge ---
  if (verbose) {
    print(paste0("[DNAmFitAge] Estimating DNAmFitAge and FitAgeAcceleration..."))
  }

  FitAge_out <- internal_FitAgeEstimator(data_FitAge_prep, IDvar = id_col)



  #### Combine Results ####
  # Define the 6 fitness biomarker names
  fitness_vars <- c("DNAmGait_noAge", "DNAmGrip_noAge", "DNAmVO2max", "DNAmGait_wAge",
                    "DNAmGrip_wAge", "DNAmFEV1_wAge")

  data_FitnessEst_sub <- data_FitnessEst[, c(id_col, "Age", "Female", fitness_vars)]


  FitAge_out_sub <- FitAge_out[, c(id_col, "DNAmGrimAge", "DNAmFitAge", "FitAgeAccel")]


  all_results_combined <- merge(data_FitnessEst_sub,
                                FitAge_out_sub,
                                by = id_col)


  col_index <- which(colnames(all_results_combined) == id_col)
  colnames(all_results_combined)[col_index] <- "SampleID"

  all_results_combined <- all_results_combined[match(sample_ids, all_results_combined$SampleID), ]

  return(all_results_combined)

}




