#' @title Estimate Proteomic Organ Ages and Age Gaps
#'
#' @description
#' This function calculates biological organ ages using pre-trained LASSO models
#' based on proteomic data (SomaScan). It performs automated data normalization, 
#' version scaling, and calculates Age Gaps relative to a reference baseline.
#'
#' @param protExp A numeric matrix or data frame of protein expression values.
#' Rows should represent proteins (SeqId) and columns should represent samples.
#' \strong{Column names must exactly match the row names of \code{metadata} in identical order.}
#' Values should be raw Relative Fluorescence Units (RFU).
#' 
#' @param metadata A data frame containing sample metadata. Must include a 
#' 'Sex' column (with 'Female' and 'Male' levels). \strong{An 'Age' column (numeric) 
#' is optional but highly recommended;} without it, Age Gaps cannot be calculated 
#' (will return NA). \strong{Row names must exactly match the column names of 
#' \code{protExp} in identical order.}
#' 
#' @param assayVersion Character string specifying the SomaScan assay version:
#' "v4", "v4.1", or "v5". Default is "v4.1".
#' 
#' @param reference Character string specifying the reference baseline for 
#' Age Gap calculation. Either "original" (training cohort baseline) or 
#' "cohort" (internal normalization based on the input dataset). Ignored if 'Age' is missing.
#' @param verbose Logical, whether to print progress messages to the console.
#'
#' @return A named \code{list} where each element corresponds to a specific organ 
#' or aging model (e.g., "Heart", "Brain", "CognitionAdipose"). Each element is a 
#' \code{data.frame} containing the following columns:
#' \describe{
#'   \item{SampleID}{Character. The unique identifier for each sample, matched from the input.}
#'   \item{Organ}{Character. The name of the organ or aging model.}
#'   \item{ChronologicalAge}{Numeric. The actual chronological age of the sample provided in the metadata.}
#'   \item{PredictedAge}{Numeric. The raw biological age estimated by the LASSO model.}
#'   \item{ExpectedAge}{Numeric. The baseline expected biological age(NA if ChronologicalAge is not provided). If \code{reference = "original"}, this is derived from the pre-trained Knight-ADRC cohort LOWESS curve. If \code{reference = "cohort"}, this is calculated dynamically using a LOESS curve fitted to the provided input dataset.}
#'   \item{AgeGap}{Numeric. The unstandardized difference between PredictedAge and ExpectedAge. (NA if ChronologicalAge is not provided).}
#'   \item{AgeGapZ}{Numeric. The Z-score standardized Age Gap (NA if ChronologicalAge is not provided). Standardized using either the training cohort's standard deviation (\code{reference = "original"}) or the internal standard deviation of the input dataset (\code{reference = "cohort"}).}
#' }
#' 
#' @importFrom stats na.omit approx loess loess.control predict approxfun lowess
#' @references
#' Oh, H.SH., Rutledge, J., Nachun, D. et al.
#' Organ aging signatures in the plasma proteome track health and disease. 
#' \emph{Nature} 2023
#' 
#' @export
#'
#' @examples
#' protExample <- loadOmniAgeRdata(
#'     "omniager_proteomic_somascan_v41_example",
#'     verbose = FALSE
#' )
#' 
#' metadata <- protExample[[1]]
#' protExp <- protExample[[2]]
#'
#' results <- wyssCorayOrganAge(protExp, metadata,assayVersion="v4.1")

wyssCorayOrganAge <- function(protExp, 
                              metadata, 
                              assayVersion = c("v4.1", "v4", "v5"),
                              reference = c("original", "cohort"),
                              verbose = TRUE) {
  
  # ---------------------------------------------------------
  # 1. Parameter Validation and Input Handling
  # ---------------------------------------------------------
  assayVersion <- match.arg(assayVersion)
  reference <- match.arg(reference)
  
  # Sex is mandatory because it is a feature in the LASSO models.
  if (!"Sex" %in% colnames(metadata)) {
    stop("Metadata is missing the required column: 'Sex' (must contain 'Female'/'Male').")
  }
  
  # Check if Age is present
  has_age <- "Age" %in% colnames(metadata)
  if (!has_age && verbose) {
    message("Note: 'Age' column not found in metadata. Only PredictedAge will be calculated. Age Gaps will be NA.")
  }
  
  # Ensure strict sample alignment between expression and metadata
  if (!identical(colnames(protExp), rownames(metadata))) {
    stop(paste0(
      "Sample alignment error: The column names of 'protExp' must exactly match ",
      "the row names of 'metadata' in both content AND order.\n",
      "Please align your datasets before running this function."
    ))
  }
  
  # Encode biological sex (Female = 1, Male = 0)
  metadata$Sex_F <- ifelse(metadata$Sex == "Female", 1, 0)
  
  # ---------------------------------------------------------
  # 2. Load Pre-trained Model Parameters
  # ---------------------------------------------------------
  omniModels <- loadOmniAgeRdata("omniager_proteomic_wysscoray_organages_model")
  
  # ---------------------------------------------------------
  # 3. Data Normalization and Pre-processing
  # ---------------------------------------------------------
  
  # Step A: Version-specific scale factor adjustment
  if (assayVersion != "v4.1") {
    scales <- omniModels$scale_factors[[assayVersion]]
    commonProts <- intersect(rownames(protExp), scales$SeqId)
    
    if (length(commonProts) > 0) {
      sFactors <- scales$ScaleFactor[match(commonProts, scales$SeqId)]
      protExp[commonProts, ] <- sweep(protExp[commonProts, , drop = FALSE], 1, sFactors, `*`)
    }
  }
  
  # Step B: Sanity check for data scale
  global_mean <- mean(as.matrix(protExp), na.rm = TRUE)
  if (global_mean < 500) {
    warning(paste0(
      "Data scale warning: The global mean of the input protein matrix is low (", 
      round(global_mean, 2), ").\n",
      "This model requires raw, un-logged Relative Fluorescence Units (RFU).\n",
      "Verify that input data has not been log-transformed or normalized."
    ), call. = FALSE)
  }
  
  # Step C: Log10 transformation
  protExpLog <- log10(protExp)
  
  # ---------------------------------------------------------
  # 4. Multi-organ Age Prediction Loop
  # ---------------------------------------------------------
  organList <- names(omniModels$organs)
  
  resList <- lapply(organList, function(organ) {
    if (verbose) message(sprintf("Processing organ: %s...", organ))
    
    cfg <- omniModels$organs[[organ]]
    
    targetProts <- cfg$lasso_weights$Feature[cfg$lasso_weights$Feature != "Sex_F"]
    
    if (!all(targetProts %in% rownames(protExpLog))) {
      warning(sprintf("Organ '%s' skipped: insufficient protein features found.", organ))
      return(NULL)
    }
    
    scaler <- cfg$prot_scaler
    idx <- match(targetProts, scaler$SeqId)
    
    z_prot <- scale(t(protExpLog[targetProts, , drop = FALSE]), 
                    center = scaler$mean[idx], 
                    scale = scaler$sd[idx])
    
    inputM <- cbind(Sex_F = metadata[rownames(z_prot), "Sex_F"], z_prot)
    beta <- cfg$lasso_weights$Weight
    
    # Perform linear prediction (This is the ONLY step that happens if Age is missing)
    predAge <- cfg$lasso_intercept + as.numeric(inputM %*% beta)
    
    # ---------------------------------------------------------
    # 5. Age Gap Calculation (Conditional on 'Age' being present)
    # ---------------------------------------------------------
    if (has_age) {
      if (reference == "original") {
        validTable <- na.omit(cfg$lowess_table)
        expAge <- approx(x = validTable$Chronological_Age,
                         y = validTable$Expected_Predicted_Age,
                         xout = metadata$Age,
                         rule = 1)$y
        
        rawGap <- predAge - expAge
        zGap <- rawGap / cfg$agegap_scaler$sd
        
      } else if (reference == "cohort") {
        # Fallback if someone forces cohort with < 30 samples but has Age
        n_samples <- length(metadata$Age)
        if (n_samples < 30) {
          stop(sprintf("Cohort mode for '%s' requires N >= 30 samples with Age.", organ))
        }
        
        lowess_fit <- lowess(x = metadata$Age, y = predAge, f = 2/3, iter = 5)
        lowess_fit_int <- approxfun(x = lowess_fit$x, y = lowess_fit$y, rule = 2)
        
        expAge <- lowess_fit_int(metadata$Age)
        rawGap <- predAge - expAge
        
        if (all(is.na(rawGap))) {
          zGap <- rep(NA_real_, length(rawGap))
        } else {
          mu <- mean(rawGap, na.rm = TRUE)
          sigma <- sqrt(mean((rawGap - mu)^2, na.rm = TRUE)) 
          zGap <- (rawGap - mu) / sigma
        }
      }
    } else {
      # If no Age is provided, set all reference-dependent metrics to NA
      expAge <- rep(NA_real_, length(predAge))
      rawGap <- rep(NA_real_, length(predAge))
      zGap   <- rep(NA_real_, length(predAge))
    }
    
    # Consolidate results
    return(data.frame(
      SampleID = rownames(z_prot),
      Organ = organ,
      ChronologicalAge = if(has_age) metadata$Age else rep(NA_real_, length(predAge)),
      PredictedAge = predAge,
      ExpectedAge = expAge,
      AgeGap = rawGap,
      AgeGapZ = zGap,
      stringsAsFactors = FALSE
    ))
  })
  
  # ---------------------------------------------------------
  # 6. Result Aggregation
  # ---------------------------------------------------------
  resList <- resList[!vapply(resList, is.null, logical(1))]
  names(resList) <- vapply(resList, function(x) x$Organ[1], character(1))
  
  if (verbose) message("Age estimation completed successfully. Returning named list.")
  
  return(resList)
}
