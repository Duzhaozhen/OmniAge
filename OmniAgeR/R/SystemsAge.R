#' @title Calculate Systems Age and 11 System-Specific Scores
#'
#' @description
#' Implements the Systems Age epigenetic clock, a multi-system aging biomarker
#' based on DNA methylation (DNAm). This function calculates a composite 'Systems Age'
#' as well as 11 individual physiological system aging scores from a single
#' blood methylation dataset.
#'
#' @param DNAm A numeric matrix or data.frame of DNA methylation beta values.
#'   **Rows must correspond to CpGs** and **columns to samples**.
#'
#' @param RData A required parameter specifying the location or content of the
#'   pre-trained 'SystemsAge_data' object. This can be either:
#'   \itemize{
#'     \item{1. A **character string** path to the directory containing the
#'       'SystemsAge_data.qs' file (e.g., from `get_OmniAgeR_path()`).}
#'     \item{2. The pre-loaded **list** object itself, loaded via
#'       `load_OmniAgeR_data(object_name = "SystemsAge_data")`.}
#'   }
#'
#'
#' @details
#' The Systems Age framework, developed by Sehgal et al. (2025), quantifies
#' aging heterogeneity across 11 distinct physiological systems. It integrates
#' supervised and unsupervised machine learning with clinical biomarkers and
#' mortality risk to derive system-specific scores.
#'
#' This function applies the pre-trained models to new methylation data. It
#' requires the `SystemsAge_data.qs` object, which contains the necessary
#' PCA models, elastic net coefficients, and scaling parameters from the
#' original study.
#'
#' This data object must first be downloaded using
#' `download_OmniAgeR_data(clocks = "SystemsAge")`.
#'
#' For computational efficiency, especially when processing multiple datasets
#' or running the function in a loop, it is recommended to load the data
#' object once using `load_OmniAgeR_data(object_name = "SystemsAge_data")`
#' and pass the resulting list to the `RData` parameter.
#'
#'
#' @return
#' A data.frame where the first column is `Sample_ID` (derived from the
#' `colnames` of the input `DNAm` matrix) and the subsequent 13 columns
#' contain the calculated scores. These include the 11 system scores
#' (e.g., 'Blood', 'Brain', 'Heart'), the 'Age_prediction' score,
#' and the composite 'SystemsAge' score. All scores are scaled to the
#' unit of years.
#'
#' @references
#' Sehgal, R., Markov, Y., Qin, C. et al.
#' Systems Age: a single blood methylation test to quantify aging heterogeneity across 11 physiological systems.
#' \emph{Nat Aging} (2025).
#'
#'
#' @export
#'
#' @examples
#' # Download the external data
#' download_OmniAgeR_data(clocks = "SystemsAge") #  ZENODO_DOI: "10.5281/zenodo.17162604"
#'
#' # Either path to the data
#' RData <- get_OmniAgeR_path()
#' # OR
#' RData <- load_OmniAgeR_data(object_name = "SystemsAge_data")
#'
#'
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' SystemsAge.o <- SystemsAge(hannum_bmiq_m,RData = RData)
#'
#'



SystemsAge <- function(DNAm, RData = NULL) {

  pheno <- data.frame(Sample_ID = colnames(DNAm))
  DNAm <- t(DNAm)
  # Check RData
  is_valid_path <- checkmate::test_character(RData, len = 1, any.missing = FALSE)
  is_valid_list <- checkmate::test_list(RData, any.missing = FALSE)

  if (!is_valid_path && !is_valid_list) {
    stop(
      paste("[SystemsAge] RData argument must be a valid file path (character) or a pre-loaded data object (list).",
            "It cannot be NULL. Please load the data first.",
            "See `?download_OmniAgeR_data` to download the required files.")
    )
  }


  # handle RData
  if (is.character(RData)) {
    RData <- load_OmniAgeR_data(object_name = "SystemsAge_data", path = RData)
  }



  if (rlang::hash(RData) != "d984914ff6aa17d8a6047fed5f9f6e4d") {
    stop("[SystemsAge] The downloaded SystemsAge data is corrupted or the wrong data (e.g., SystemsAge) was passed. See `?download_methylCIPHER()`.")
  }

  ## Imputation
  DNAm <- PCClocks_impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = RData$imputeMissingCpGs,
    subset = TRUE
  )

  ## Re-align to make sure things lined up with the object
  DNAm <- DNAm[, names(RData$imputeMissingCpGs), drop = F]

  ## Calculate methylation PCs
  DNAmPCs <- predict(RData$DNAmPCA, DNAm)



  ## Calculate DNAm system PCs then system scores
  DNAmSystemPCs <- DNAmPCs[, 1:4017] %*% as.matrix(RData$system_vector_coefficients[1:4017, ])
  system_scores <- matrix(nrow = dim(DNAmSystemPCs)[1], ncol = 11)
  i <- 1
  groups <- c("Blood", "Brain", "Cytokine", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metab", "Lung", "MusculoSkeletal")
  for (group in groups) {
    tf <- grepl(group, colnames(DNAmSystemPCs))
    sub <- DNAmSystemPCs[, tf]
    sub_system_coefficients <- RData$system_scores_coefficients_scale[tf]
    if (length(sub_system_coefficients) == 1) {
      system_scores[, i] <- sub * -1
    } else {
      system_scores[, i] <- sub %*% sub_system_coefficients
    }
    i <- i + 1
  }
  colnames(system_scores) <- groups

  ## Generate predicted chronological age
  Age_prediction <- as.matrix(DNAmPCs) %*% as.matrix(RData$Predicted_age_coefficients[2:4019]) + RData$Predicted_age_coefficients[1]
  Age_prediction <- Age_prediction * RData$Age_prediction_model[2] + (Age_prediction**2) * RData$Age_prediction_model[3] + RData$Age_prediction_model[1]
  Age_prediction <- Age_prediction / 12
  system_scores <- cbind(system_scores, Age_prediction)
  colnames(system_scores)[12] <- "Age_prediction"

  ## Generating overall system index
  colnames(system_scores) <- c("Blood", "Brain", "Inflammation", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metabolic", "Lung", "MusculoSkeletal", "Age_prediction")
  system_PCA <- predict(RData$systems_PCA, system_scores)
  pred <- system_PCA %*% RData$Systems_clock_coefficients
  system_scores <- cbind(system_scores, pred)
  colnames(system_scores)[13] <- "SystemsAge"

  ## Scale system ages
  system_ages <- system_scores
  for (i in c(1:13)) {
    y <- system_ages[, i]
    system_ages[, i] <- (((y - RData$transformation_coefs[i, 1]) / RData$transformation_coefs[i, 2]) * RData$transformation_coefs[i, 4]) + RData$transformation_coefs[i, 3]
    system_ages[, i] <- system_ages[, i] / 12
  }

  if (is.null(pheno)) {
    results <- as.data.frame(system_ages)
  } else {
    results <- cbind(pheno, system_ages)
  }
  row.names(results) <- NULL
  return(results)
}
