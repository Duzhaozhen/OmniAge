#' @title Predict Cortical DNA Methylation Clock Age (2020)
#'
#' @description Predicts DNAm age in cortical samples using the elastic net model by Shireby et al. (2020).
#'
#' @param beta.m A numeric matrix of beta values. Rows should be CpG probes and columns should be individual samples.
#'
#' @return A named numeric vector of predicted cortical DNAm ages for each sample.
#'
#' @seealso The main function \code{\link{EpiAge}} can be used to calculate multiple clocks simultaneously.
#'
#' @references
#' Shireby GL, Davies JP, Francis PT, et al.
#' Recalibrating the epigenetic clock: implications for assessing biological age in the human cortex.
#' \emph{Brain.} 2020
#'
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' CorticalClock.out <- CorticalClock(hannum_bmiq_m)
#'


CorticalClock <- function(beta.m) {

  # --- Step 1: Load and parse coefficients (from package internal data) ---
  data("CorticalClockCoef")

  coef_data <- CorticalClock_coef
  ref_data <- CorticalClock_ref # This is the reference data for imputation
  clock.name <- "CorticalClock"

  # 1a. Parse into the coef.lv list format required by the template
  coef.lv <- list()

  # Intercept
  coef.lv[[1]] <- coef_data[1,2]

  # Coefficients
  coefficients_vec <- coef_data[-1,2]
  names(coefficients_vec) <- coef_data[-1,1]
  coef.lv[[2]] <- coefficients_vec

  # --- Step 2: Define the non-linear transformation function (anti.trafo) ---
  # This function is identical to the one in Horvath 2013/2018
  iageF <- function(ptage.v, adult.age = 20) {
    y.idx <- which(ptage.v <= 0)
    a.idx <- which(ptage.v > 0)
    mage.v <- ptage.v
    mage.v[a.idx] <- ptage.v[a.idx] * (adult.age + 1) + adult.age
    mage.v[y.idx] <- exp(ptage.v[y.idx] + log(adult.age + 1)) - 1
    return(mage.v)
  }


  # --- Step 3: Calculate the linear predictor (Imputation and Weighted Sum) ---

  # 3a. Parse variables from coef.lv
  intercept <- coef.lv[[1]]
  coefficients <- coef.lv[[2]]
  clock_probes <- names(coefficients)

  # 3b. Probe matching and imputation of missing probes
  user_probes <- rownames(beta.m)
  missing_probes <- setdiff(clock_probes, user_probes)

  beta_complete <- beta.m

  print(paste0("[CorticalClock]",
    "Number of represented ", clock.name, " CpGs (max=",
    length(clock_probes), ")=", (length(clock_probes)-length(missing_probes))))


  if (length(missing_probes) > 0) {
    if (is.null(ref_data)) {
      # If no reference data is provided, stop and throw an error
      stop(paste(clock.name, "Error: Required probes are missing and no imputation_ref provided:",
                 paste(missing_probes, collapse=", ")))
    }

    # Check if the reference data contains all missing probes
    if (!all(missing_probes %in% names(ref_data))) {
      stop(paste(clock.name, "Error: imputation_ref is missing data for some required probes."))
    }

    # Get values from the reference data
    ref_values <- ref_data[missing_probes]

    # Create a matrix matching the dimensions of beta.m (cols = N samples)
    ref_matrix <- matrix(
      rep(ref_values, ncol(beta.m)),
      ncol = ncol(beta.m),
      byrow = FALSE
    )
    rownames(ref_matrix) <- names(ref_values)
    colnames(ref_matrix) <- colnames(beta.m)

    # Combine the imputed probes with the original data
    beta_complete <- rbind(beta.m, ref_matrix)


  } #else {
  #print(paste0(clock.name, ": All ", length(clock_probes), " required probes are present."))
  #}

  # 3c. Ensure order and handle NA values

  # Strictly order the beta matrix according to the clock coefficients
  beta_ordered <- beta_complete[clock_probes, , drop = FALSE]

  # Internal function: Impute NA for a single probe (row)
  # Uses the mean of that probe across *all samples*
  imputeNA_row <- function(cpg_vector) {
    if (anyNA(cpg_vector)) {
      cpg_mean <- mean(cpg_vector, na.rm = TRUE)
      # if (is.nan(cpg_mean)) {
      #   # If a probe is NA in all samples, fill with 0.5
      #   print(paste0(clock.name, ": Warning: A probe was NA in all samples. Imputing with 0.5."))
      #   cpg_mean <- 0.5
      # }
      cpg_vector[is.na(cpg_vector)] <- cpg_mean
    }
    return(cpg_vector)
  }

  # Use apply (1 = by row) to impute NAs
  beta_final <- t(apply(beta_ordered, 1, imputeNA_row))
  colnames(beta_final) <- colnames(beta_ordered) # apply drops colnames

  # 3d. Calculate the linear predictor (Weighted Sum)

  # Matrix multiplication: (N_samples x N_probes) %*% (N_probes x 1) -> (N_samples x 1)
  # t(beta_final) has dimensions (N_samples x N_probes)
  # coefficients has dimensions (N_probes x 1)
  pred.v <- (t(beta_final) %*% coefficients) + intercept

  # Clean up and return
  final_predictor <- as.vector(pred.v)
  names(final_predictor) <- rownames(pred.v) # Match sample IDs

  # Assign the result to predTage.v
  predTage.v <- final_predictor

  # --- Step 4: Apply the non-linear transformation ---
  predMage.v <- iageF(predTage.v)

  # --- Step 5: Return the final age vector ---
  return(predMage.v)
}
