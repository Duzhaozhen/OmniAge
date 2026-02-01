#' DNA Methylation Cell-Type Fraction (CTF) Aging Clock
#'
#' @description
#' Predicts biological age based on immune cell type fractions derived from DNA methylation data.
#'
#' @param CTF_m A numeric matrix or data frame where rows are samples and columns are cell types.
#'   Must contain the specific cell types required by the model (e.g., predicted by EpiDISH).
#'
#' @return A named numeric vector of predicted ages.
#' @export
#' @importFrom randomForest randomForest importance
#'
#'
#' @examples
#'
#' download_OmniAgeR_example("TZH_example_CTF")
#' load_OmniAgeR_example("TZH_example_CTF")
#' DNAm_CTF_Clock_o<-DNAm_CTF_Clock(CTF_m = TZH_Frac_m)
#'

DNAm_CTF_Clock <- function(CTF_m){
  # load model
  data("DNAm_CTF_model")

  # 1. Validate input format
  if (!is.matrix(CTF_m) && !is.data.frame(CTF_m)) {
    stop("Error: 'CTF_m' must be a matrix or data frame.")
  }

  # 2. Data validation and preparation
  # Retrieve required features and ensure data integrity (no missing columns or NAs).
  required_features <- rownames(DNAm_CTF_model$importance)

  # 2.1 Check for missing features (columns)
  missing_cols <- setdiff(required_features, colnames(CTF_m))
  if (length(missing_cols) > 0) {
    stop(paste0(
      "[DNAm_CTF_Clock] Input matrix is missing required cell types:\n",
      paste(missing_cols, collapse = ", "),
      "\nPlease ensure you used EpiDISH or similar tools to estimate all 12 immune cell fractions."
    ))
  }

  # 2.2 Subset and align data
  # We subset first to ensure we only check for NAs in the relevant columns.
  data_for_prediction <- as.data.frame(CTF_m)[, required_features, drop = FALSE]

  # 2.3 Check for missing values (NA)
  # Strict validation: The model requires complete data.
  if (any(is.na(data_for_prediction))) {
    stop("[DNAm_CTF_Clock] Input data contains NA values in required cell types. The model requires complete data without missing values.")
  }
  # 3. Perform prediction
  pred_age <- predict(DNAm_CTF_model, newdata = data_for_prediction)

  # 4. Format output
  if (!is.null(rownames(CTF_m))) {
    names(pred_age) <- rownames(CTF_m)
  }

  return(pred_age)

}
