# -------------------------------------------------------------------
# 1. Data Preparation (Modified for R package)
#    - Accepts DNAmFitnessModels as an argument
#    - Uses prefixed message()
# -------------------------------------------------------------------

#' Prepare and Impute Data for DNAm Fitness Models
#'
#' Subsets the input data to required columns and imputes any missing
#' CpGs using pre-calculated sex-specific medians.
#'
#' @param dataset A data.frame or matrix containing sample metadata and DNAm
#'   beta values. Must include columns specified by `idvariable`, 'Female',
#'   'Age', and all available CpG beta values.
#' @param idvariable A string representing the name of the unique sample
#'   identifier column in `dataset`.
#' @param DNAmFitnessModels An internal data object (list) containing
#'   model coefficients, the master CpG list (`AllCpGs`), and imputation
#'   medians (`Female_Medians_All`, `Male_Medians_All`).
#' @param verbose Logical. If TRUE, prints messages about the number of
#'   CpGs represented and imputed.
#'
#' @return A data.frame, subsetted to required columns and imputed with
#'   sex-specific medians for any missing CpGs, ready for downstream
#'   model estimation.
#'
#' @noRd


internal_data_prep <- function(dataset, idvariable, DNAmFitnessModels,verbose){

  extract_these <- colnames(dataset)[which(colnames(dataset) %in% DNAmFitnessModels$AllCpGs)]
  output_data <- dataset[, c(idvariable, "Female", "Age", extract_these)]

  if(length(extract_these) != length(DNAmFitnessModels$AllCpGs)){
    data_fem <- output_data[output_data$Female == 1, ]
    data_male <- output_data[output_data$Female == 0, ]

    cpgs_toadd <- !colnames(DNAmFitnessModels$Female_Medians_All) %in% extract_these

    if(nrow(data_fem) != 0){
      data_fem <- data.frame(data_fem, DNAmFitnessModels$Female_Medians_All[, cpgs_toadd])
    }

    if(nrow(data_male) != 0){
      data_male <- data.frame(data_male, DNAmFitnessModels$Male_Medians_All[, cpgs_toadd])
    }


    output_data <- rbind(data_fem, data_male)

    if (verbose) {print(paste0("Number of represented DNAmFitAge",
                               " CpGs (max=",length(DNAmFitnessModels$AllCpGs), ")=", length(extract_these)))
    }


  }
  return(output_data)
}

# -------------------------------------------------------------------
# 2. Single Model Estimator
# -------------------------------------------------------------------

#' Generic DNAm Model Estimator
#'
#' Applies a pre-defined linear model (in tidy format) to a dataset.
#'
#' @details
#' This generic estimator function constructs a design matrix from the
#' `dataset` based on the `term` column of the `TidyModel` object.
#' It then performs matrix multiplication with the `estimate` column
#' (the coefficients) to calculate the DNAm-based score. It automatically
#' includes an intercept term.
#'
#' @param dataset A data.frame containing the necessary predictor
#'   variables (e.g., CpGs) as columns.
#' @param TidyModel A "tidy" model data.frame, typically from the
#'   `DNAmFitnessModels` object. Must contain 'term' (predictor names)
#'   and 'estimate' (coefficient) columns.
#' @param IDvar A string representing the name of the unique sample
#'   identifier column in `dataset`.
#'
#' @return A data.frame with two columns: 'DNAmEstimate' (the calculated
#'   score) and 'ID' (the sample identifier).
#'
#' @noRd
internal_DNAmEstimatorAnyModel <- function(dataset, TidyModel, IDvar){

  int_length <- nrow(dataset)
  mod_length <- length(TidyModel$term)

  Xdat <- dataset[, colnames(dataset) %in% TidyModel$term]
  Xdat <- data.frame("Intercept" = rep(1, int_length), Xdat)
  Xdatnew <- as.matrix(Xdat[, c("Intercept", TidyModel$term[2:mod_length])])
  if(sum(colnames(Xdatnew)[2:mod_length] == TidyModel$term[2:mod_length]) == mod_length-1){
    estimate <- Xdatnew %*% TidyModel$estimate
  }
  if(sum(colnames(Xdatnew)[2:mod_length] == TidyModel$term[2:mod_length]) != mod_length-1){

    message("[DNAmFitAge] ERROR: Not All Columns in New Dataframe for model.")

    estimate <- rep(NA, int_length)
  }

  est_data <- data.frame(DNAmEstimate = estimate, ID = dataset[, {{IDvar}}])

  return(est_data)
}

# -------------------------------------------------------------------
# 3. Fitness Estimator (Modified for R package)
#    - Accepts DNAmFitnessModels as an argument
#    - Calls internal_DNAmEstimatorAnyModel
# -------------------------------------------------------------------

#' Calculate All DNAm-Based Fitness Estimators
#'
#' A dispatcher function that calculates multiple DNAm fitness scores.
#'
#' @details
#' This function splits the input data by sex ('Female' == 1 or 0)
#' and applies the appropriate sex-specific models (e.g., for Gait,
#' Grip, FEV1) by repeatedly calling `internal_DNAmEstimatorAnyModel`.
#' The `VO2maxModel` is applied to both sexes. All individual
#' estimates are then merged by the sample ID into a single data.frame.
#'
#' @param data The prepared data.frame, typically the output from
#'   `internal_data_prep`. Must contain 'Female' and all required
#'   predictor columns.
#' @param IDvar A string representing the name of the unique sample
#'   identifier column.
#' @param DNAmFitnessModels The internal data object (list) containing all
#'   required tidy models (e.g., `Gait_noAge_Females`, `VO2maxModel`, etc.).
#'
#' @return A data.frame containing the original data merged with new
#'   columns for each calculated fitness estimator
#'   (e.g., 'DNAmGait_noAge', 'DNAmGrip_noAge', etc.).
#'
#' @noRd
internal_DNAmFitnessEstimators <- function(data, IDvar, DNAmFitnessModels){

  # Define variables inside the function
  DNAmFitness_Xvars <- c("DNAmGait_noAge", "DNAmGrip_noAge", "DNAmVO2max", "DNAmGait_wAge",
                         "DNAmGrip_wAge", "DNAmFEV1_wAge", "DNAmGrimAge")

  data_fem <- data[data$Female ==1,]
  data_male <- data[data$Female ==0,]

  fem_est1 <- internal_DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_noAge_Females, IDvar)
  fem_est2 <- internal_DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_noAge_Females, IDvar)
  fem_est3 <- internal_DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$VO2maxModel, IDvar)
  fem_est4 <- internal_DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_wAge_Females, IDvar)
  fem_est5 <- internal_DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_wAge_Females, IDvar)
  fem_est6 <- internal_DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$FEV1_wAge_Females, IDvar)

  male_est1 <- internal_DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_noAge_Males, IDvar)
  male_est2 <- internal_DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_noAge_Males, IDvar)
  male_est3 <- internal_DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$VO2maxModel, IDvar)
  male_est4 <- internal_DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_wAge_Males, IDvar)
  male_est5 <- internal_DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_wAge_Males, IDvar)
  male_est6 <- internal_DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$FEV1_wAge_Males, IDvar)

  fem_est1 <- rbind(fem_est1, male_est1)
  fem_est2 <- rbind(fem_est2, male_est2)
  fem_est3 <- rbind(fem_est3, male_est3)
  fem_est4 <- rbind(fem_est4, male_est4)
  fem_est5 <- rbind(fem_est5, male_est5)
  fem_est6 <- rbind(fem_est6, male_est6)

  fem_est1[, DNAmFitness_Xvars[1]] <- fem_est1$DNAmEstimate
  fem_est2[, DNAmFitness_Xvars[2]] <- fem_est2$DNAmEstimate
  fem_est3[, DNAmFitness_Xvars[3]] <- fem_est3$DNAmEstimate
  fem_est4[, DNAmFitness_Xvars[4]] <- fem_est4$DNAmEstimate
  fem_est5[, DNAmFitness_Xvars[5]] <- fem_est5$DNAmEstimate
  fem_est6[, DNAmFitness_Xvars[6]] <- fem_est6$DNAmEstimate

  all_ests <- Reduce(function(x,y) merge(x = x, y = y, by = "ID", all.x = TRUE, all.y = TRUE),
                     list(fem_est1[!is.na(fem_est1$ID), 2:3],
                          fem_est2[!is.na(fem_est2$ID), 2:3],
                          fem_est3[!is.na(fem_est3$ID), 2:3],
                          fem_est4[!is.na(fem_est4$ID), 2:3],
                          fem_est5[!is.na(fem_est5$ID), 2:3],
                          fem_est6[!is.na(fem_est6$ID), 2:3]))

  match1 <- match(data[, IDvar], all_ests$ID)
  data_and_est <- data.frame(data, all_ests[match1, ])

  return(data_and_est)
}

# -------------------------------------------------------------------
# 4. FitAge Estimator (Modified)
#    - Uses prefixed message()
# -------------------------------------------------------------------
#' Calculate DNAmFitAge and FitAgeAccel
#'
#' Calculates the final DNAmFitAge score from its components and the
#' corresponding acceleration value.
#'
#' @details
#' This function calculates `DNAmFitAge` by applying a sex-specific,
#' hard-coded linear combination of the standardized component scores
#' (VO2max, Grip, Gait, GrimAge). It operates only on complete cases
#' for these predictors.
#'
#' It also calculates `FitAgeAccel` by taking the residuals of a
#' linear model regressing `DNAmFitAge` on chronological `Age`.
#'
#' @param data The data.frame containing all required DNAm fitness
#'   estimators, typically the output from `internal_DNAmFitnessEstimators`.
#'   Must include 'Age', 'Female', 'DNAmGait_noAge', 'DNAmGrip_noAge',
#'   'DNAmVO2max', and 'DNAmGrimAge'.
#' @param IDvar A string representing the name of the unique sample
#'   identifier column.
#'
#' @return A data.frame containing the sample `IDvar`, 'Age', all
#'   component scores, the calculated 'DNAmFitAge', and 'FitAgeAccel'.
#'
#' @noRd
internal_FitAgeEstimator <- function(data, IDvar){

  FitAge_Xvars <- c("DNAmGait_noAge", "DNAmGrip_noAge", "DNAmVO2max", "DNAmGrimAge")
  fem <- data[data$Female == 1, c("Age", FitAge_Xvars, IDvar)]
  male <- data[data$Female == 0, c("Age", FitAge_Xvars, IDvar)]

  fem_comcase <- fem[complete.cases(fem), ]
  male_comcase <- male[complete.cases(male), ]

  female_fitest <- 0.1044232 * ((fem_comcase$DNAmVO2max - 46.825091) / (-0.13620215)) +
    0.1742083 * ((fem_comcase$DNAmGrip_noAge - 39.857718) / (-0.22074456)) +
    0.2278776 * ((fem_comcase$DNAmGait_noAge - 2.508547) / (-0.01245682))  +
    0.4934908 * ((fem_comcase$DNAmGrimAge - 7.978487) / (0.80928530))

  male_fitest <- 0.1390346 * ((male_comcase$DNAmVO2max - 49.836389) / (-0.141862925)) +
    0.1787371 * ((male_comcase$DNAmGrip_noAge - 57.514016) / (-0.253179827)) +
    0.1593873 * ((male_comcase$DNAmGait_noAge - 2.349080) / (-0.009380061))  +
    0.5228411 * ((male_comcase$DNAmGrimAge - 9.549733) / (0.835120557))

  fem_X <- data.frame(fem_comcase, DNAmFitAge = female_fitest)
  male_X <- data.frame(male_comcase, DNAmFitAge = male_fitest)

  returned_X <- rbind(fem_X, male_X)
  returned_X$FitAgeAccel <- residuals(lm(DNAmFitAge ~ Age, data = returned_X))

  Female_meanAbsDev <- mean(abs(fem_comcase$Age - female_fitest), na.rm=TRUE)
  #message(paste0("[DNAmFitAge] QC: Female Mean Absolute Deviation (MAD) from Age: ", round(Female_meanAbsDev, 2)))

  Male_meanAbsDev <- mean(abs(male_comcase$Age - male_fitest), na.rm=TRUE)
  #message(paste0("[DNAmFitAge] QC: Male Mean Absolute Deviation (MAD) from Age: ", round(Male_meanAbsDev, 2)))
  return(returned_X)
}


