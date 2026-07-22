#' @title Calculate DNAmFitAge and DNAm fitness biomarkers
#'
#' @description
#' Calculates six DNA methylation-based fitness biomarkers and DNAmFitAge
#' from a DNA methylation beta-value matrix, chronological age, sex and
#' pre-calculated DNAmGrimAge values.
#'
#' @param betaM A numeric DNA methylation beta-value matrix with CpG probe
#'   identifiers as row names and sample identifiers as column names.
#'   \code{colnames(betaM)} must be provided and must be identical to
#'   \code{names(age)}, \code{names(sex)} and
#'   \code{names(grimageVector)}.
#'
#' @param age A numeric vector containing chronological age for each sample.
#'   A named vector is recommended. If names are supplied, they must be
#'   identical to \code{colnames(betaM)}, including sample order. An unnamed
#'   vector is accepted when its length equals \code{ncol(betaM)}; in this
#'   case, values are matched to samples according to the column order of
#'   \code{betaM}, and a warning is issued.
#'
#' @param sex A character or factor vector containing sex for each sample.
#'   Values must be either \code{"Male"} or \code{"Female"}. A named vector
#'   is recommended. If names are supplied, they must be identical to
#'   \code{colnames(betaM)}, including sample order. An unnamed vector is
#'   accepted when its length equals \code{ncol(betaM)}; in this case, values
#'   are matched to samples according to the column order of \code{betaM},
#'   and a warning is issued.
#'
#' @param grimageVector A numeric vector containing pre-calculated
#'   DNAmGrimAge values for each sample. A named vector is recommended.
#'   If names are supplied, they must be identical to
#'   \code{colnames(betaM)}, including sample order. An unnamed vector is
#'   accepted when its length equals \code{ncol(betaM)}; in this case, values
#'   are matched to samples according to the column order of \code{betaM},
#'   and a warning is issued.
#'
#' @param minCoverage A numeric value between 0 and 1 specifying the minimum
#'   proportion of required CpGs that must be present. Default is 0.5.
#'
#' @param verbose Logical. Whether to print progress messages.
#'   Default is \code{TRUE}.
#'   
#' @return A data frame with one row per sample and the following columns:
#' \itemize{
#'   \item \code{SampleID}: Sample identifier.
#'   \item \code{DNAmVO2max}: DNAm-predicted VO2max.
#'   \item \code{DNAmGait_noAge}: Gait-speed biomarker without age.
#'   \item \code{DNAmGrip_noAge}: Grip-strength biomarker without age.
#'   \item \code{DNAmGait_wAge}: Gait-speed biomarker including age.
#'   \item \code{DNAmGrip_wAge}: Grip-strength biomarker including age.
#'   \item \code{DNAmFEV1_wAge}: FEV1 biomarker including age.
#'   \item \code{DNAmGrimAge}: The supplied DNAmGrimAge value.
#'   \item \code{DNAmFitAge}: The calculated DNAmFitAge value.
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
#' hannumExample <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )
#'
#' hannumBmiqM <- hannumExample[[1]]
#' phenoTypesHannum <- hannumExample[[2]]
#' sampleIds <- colnames(hannumBmiqM)
#'
#' age <- setNames(
#'     phenoTypesHannum$Age,
#'     sampleIds
#' )
#'
#' sex <- setNames(
#'     ifelse(
#'         phenoTypesHannum$Sex == "F",
#'         "Female",
#'         "Male"
#'     ),
#'     sampleIds
#' )
#'
#' grimAge1Out <- grimAge1(
#'     betaM = hannumBmiqM,
#'     age = age,
#'     sex = sex
#' )
#'
#' grimageVector <- setNames(
#'     grimAge1Out$DNAmGrimAge1,
#'     grimAge1Out$SampleID
#' )
#'
#' dnamFitAgeOut <- dnamFitAge(
#'     betaM = hannumBmiqM,
#'     age = age,
#'     sex = sex,
#'     grimageVector = grimageVector
#' )
#' 

dnamFitAge <- function(
    betaM,
    age,
    sex,
    grimageVector,
    minCoverage = 0.5,
    verbose = TRUE
) {
  # --- 1. Validate and align input objects ---
  betaM <- .validateBetaMatrix(
    betaM,
    requireColnames = TRUE
  )
  
  sampleIds <- colnames(betaM)
  
  age <- .validateAge(
    age = age,
    sampleNames = sampleIds,
    reorder = FALSE
  )
  
  sex <- .validateSex(
    sex = sex,
    sampleNames = sampleIds,
    allowedValues = c("Male", "Female"),
    reorder = FALSE
  )
  
  grimageVector <- .validateGrimAgeVector(
    grimageVector = grimageVector,
    sampleNames = sampleIds,
    reorder = FALSE
  )
  
  if (!is.numeric(minCoverage) ||
      length(minCoverage) != 1L ||
      is.na(minCoverage) ||
      minCoverage < 0 ||
      minCoverage > 1) {
    stop(
      "`minCoverage` must be a single numeric value between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (!is.logical(verbose) ||
      length(verbose) != 1L ||
      is.na(verbose)) {
    stop(
      "`verbose` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  # --- 2. Load model data ---
  DNAmFitnessModels <- loadOmniAgeRdata(
    "omniager_dnamfitage_coef",
    verbose = verbose
  )
  
  femaleNumeric <- ifelse(
    sex == "Female",
    1,
    0
  )
  
  # --- 3. Prepare input data ---
  betaTrans <- t(betaM)
  
  dataPrep <- .prepareFitAgeData(
    betaM = betaTrans,
    sampleIds = sampleIds,
    femaleVec = unname(femaleNumeric),
    ageVec = unname(age),
    modelData = DNAmFitnessModels,
    minCoverage = minCoverage,
    verbose = verbose
  )
  
  # Return a consistently structured result when coverage is insufficient.
  if (is.null(dataPrep)) {
    return(data.frame(
      SampleID = sampleIds,
      DNAmVO2max = NA_real_,
      DNAmGait_noAge = NA_real_,
      DNAmGrip_noAge = NA_real_,
      DNAmGait_wAge = NA_real_,
      DNAmGrip_wAge = NA_real_,
      DNAmFEV1_wAge = NA_real_,
      DNAmGrimAge = unname(grimageVector),
      DNAmFitAge = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # --- 4. Estimate fitness biomarkers ---
  fitnessEst <- .estimateFitnessMarkers(
    dataPrep,
    DNAmFitnessModels
  )
  
  fitnessEst$DNAmGrimAge <- unname(grimageVector)
  
  # --- 5. Calculate DNAmFitAge ---
  finalResults <- .calculateFinalFitAge(
    fitnessEst
  )
  
  resultOrder <- match(
    sampleIds,
    finalResults$SampleID
  )
  
  if (anyNA(resultOrder)) {
    stop(
      "Internal sample realignment failed during DNAmFitAge calculation.",
      call. = FALSE
    )
  }
  
  finalResults[
    resultOrder,
    ,
    drop = FALSE
  ]
}







#' Prepare and Impute Data for FitAge Calculation
#'
#' @description
#' Filters CpGs, checks coverage, and performs sex-specific median imputation.
#'
#' @param betaM Transposed beta matrix (Samples x CpGs).
#' @param sampleIds Character vector of sample identifiers.
#' @param femaleVec Numeric vector (1 for Female, 0 for Male).
#' @param ageVec Numeric vector of chronological age.
#' @param modelData List containing AllCpGs and sex-specific medians.
#' @param minCoverage Minimum coverage threshold (0-1).
#' @param verbose Logical.
#'
#' @return A data.frame with metadata and complete CpG features, or NULL
#' if coverage is too low.
#' @keywords internal
#' @noRd

.prepareFitAgeData <- function(betaM, sampleIds, femaleVec, ageVec, modelData,
                               minCoverage, verbose) {
    allRequired <- modelData$AllCpGs
    presentCpGs <- intersect(colnames(betaM), allRequired)

    coverageRatio <- length(presentCpGs) / length(allRequired)

    if (verbose) {
        message(sprintf(
            "[DNAmFitAge] Probe Check: Found %d / %d required CpGs (%.1f%%).",
            length(presentCpGs), length(allRequired),
            coverageRatio * 100
        ))
    }

    if (coverageRatio < minCoverage) {
        if (verbose) {
            warning(sprintf(
                "[DNAmFitAge] Aborted: Coverage (%.1f%%) is below your threshold (%.1f%%).",
                coverageRatio * 100, minCoverage * 100
            ))
        }
        return(NULL)
    }

    df <- data.frame(SampleID = sampleIds, Female = femaleVec, Age = ageVec)
    betaSubset <- betaM[, presentCpGs, drop = FALSE]

    missingCpGs <- setdiff(allRequired, presentCpGs)

    if (length(missingCpGs) > 0) {
        if (verbose) {
            message(sprintf(
                "[DNAmFitAge] Imputing %d missing sites using sex-specific medians...",
                length(missingCpGs)
            ))
        }

        #
        imputeMat <- matrix(0, nrow = length(sampleIds), ncol = length(missingCpGs))
        colnames(imputeMat) <- missingCpGs

        isFemale <- (femaleVec == 1)

        if (any(isFemale)) {
            imputeMat[isFemale, ] <- rep(
                as.numeric(modelData$Female_Medians_All[1, missingCpGs]),
                each = sum(isFemale)
            )
        }
        if (any(!isFemale)) {
            imputeMat[!isFemale, ] <- rep(
                as.numeric(modelData$Male_Medians_All[1, missingCpGs]),
                each = sum(!isFemale)
            )
        }
        betaSubset <- cbind(betaSubset, imputeMat)
    }


    return(cbind(df, betaSubset[, allRequired, drop = FALSE]))
}


#' Dispatch and Estimate Individual Fitness Markers
#'
#' @param data Prepared data.frame from .prepareFitAgeData.
#' @param modelData Internal model coefficients list.
#'
#' @return A data.frame with sample IDs and estimated fitness biomarkers.
#' @keywords internal
#' @noRd
.estimateFitnessMarkers <- function(data, modelData) {
    # Define model pairs for dispatch
    clocks <- list(
        DNAmGait_noAge = c("Gait_noAge_Females", "Gait_noAge_Males"),
        DNAmGrip_noAge = c("Grip_noAge_Females", "Grip_noAge_Males"),
        DNAmGait_wAge  = c("Gait_wAge_Females", "Gait_wAge_Males"),
        DNAmGrip_wAge  = c("Grip_wAge_Females", "Grip_wAge_Males"),
        DNAmFEV1_wAge  = c("FEV1_wAge_Females", "FEV1_wAge_Males")
    )

    # Add VO2max (unisex model)
    res <- data[, c("SampleID", "Female", "Age")]

    # Calculate unisex VO2max
    res$DNAmVO2max <- .applyTidyModel(data, modelData$VO2maxModel)

    # Calculate sex-specific markers
    for (marker in names(clocks)) {
        femModel <- modelData[[clocks[[marker]][1]]]
        maleModel <- modelData[[clocks[[marker]][2]]]

        scores <- rep(NA_real_, nrow(data))
        scores[data$Female == 1] <- .applyTidyModel(
            data[data$Female == 1, ],
            femModel
        )
        scores[data$Female == 0] <- .applyTidyModel(
            data[data$Female == 0, ],
            maleModel
        )
        res[[marker]] <- scores
    }

    return(res)
}

#' Apply a Tidy Model for Linear Prediction
#'
#' @param df Feature data.frame.
#' @param tidyMod Data.frame with 'term' and 'estimate' columns.
#'
#' @return A numeric vector of predicted values.
#' @noRd
.applyTidyModel <- function(df, tidyMod) {
    # Extract terms (excluding intercept)
    vars <- tidyMod$term[-1]
    # Matrix multiplication: Intercept + (X %*% weights)
    score <- tidyMod$estimate[1] + (as.matrix(df[, vars]) %*% tidyMod$estimate[-1])
    return(as.vector(score))
}


#' Final Aggregation and FitAge Score Calculation
#'
#' @description
#' Applies sex-specific coefficients to standardize and combine fitness markers.
#'
#' @param data Data.frame containing all component fitness biomarkers.
#'
#' @return A data frame containing the component fitness biomarkers,
#'   supplied DNAmGrimAge values and calculated DNAmFitAge values.
#'   
#' @keywords internal
#' @noRd
.calculateFinalFitAge <- function(data) {
    # Identify complete cases for final aggregation
    compIdx <- stats::complete.cases(data[, c(
        "DNAmGait_noAge", "DNAmGrip_noAge",
        "DNAmVO2max", "DNAmGrimAge"
    )])
    data$DNAmFitAge <- NA_real_

    # Hard-coded coefficients from McGreevy 2023
    #

    # Females
    fIdx <- which(compIdx & data$Female == 1)
    if (length(fIdx) > 0) {
        d <- data[fIdx, ]
        data$DNAmFitAge[fIdx] <- 0.1044232 * ((d$DNAmVO2max - 46.825091) / -0.13620215) +
            0.1742083 * ((d$DNAmGrip_noAge - 39.857718) / -0.22074456) +
            0.2278776 * ((d$DNAmGait_noAge - 2.508547) / -0.01245682) +
            0.4934908 * ((d$DNAmGrimAge - 7.978487) / 0.80928530)
    }

    # Males
    mIdx <- which(compIdx & data$Female == 0)
    if (length(mIdx) > 0) {
        d <- data[mIdx, ]
        data$DNAmFitAge[mIdx] <- 0.1390346 * ((d$DNAmVO2max - 49.836389) / -0.141862925) +
            0.1787371 * ((d$DNAmGrip_noAge - 57.514016) / -0.253179827) +
            0.1593873 * ((d$DNAmGait_noAge - 2.349080) / -0.009380061) +
            0.5228411 * ((d$DNAmGrimAge - 9.549733) / 0.835120557)
    }

    cols_to_remove <- c("Age", "Female")
    data <- data[, !(names(data) %in% cols_to_remove)]
    return(data)
}

#' Validate a pre-calculated DNAmGrimAge vector
#'
#' @param grimageVector A numeric vector containing pre-calculated
#'   DNAmGrimAge values.
#' @param sampleNames Expected sample identifiers.
#' @param reorder Logical. Whether to reorder the vector when the identifiers
#'   match but occur in a different order.
#' @param allowUnnamed Logical. Whether an unnamed vector is accepted and
#'   matched to samples by position. Default is \code{TRUE}.
#'
#' @return A validated numeric vector.
#'
#' @keywords internal
#' @noRd
.validateGrimAgeVector <- function(
    grimageVector,
    sampleNames,
    reorder = FALSE,
    allowUnnamed = TRUE
) {
  if (!is.numeric(grimageVector)) {
    stop(
      "`grimageVector` must be a numeric vector.",
      call. = FALSE
    )
  }
  
  if (length(grimageVector) != length(sampleNames)) {
    stop(
      "`grimageVector` must contain exactly one value for each ",
      "sample in betaM.",
      call. = FALSE
    )
  }
  
  if (anyNA(grimageVector) ||
      any(!is.finite(grimageVector))) {
    stop(
      "`grimageVector` must not contain missing or non-finite values.",
      call. = FALSE
    )
  }
  
  grimageVector <- .validateSampleNames(
    x = grimageVector,
    sampleNames = sampleNames,
    argumentName = "grimageVector",
    reorder = reorder,
    allowUnnamed = allowUnnamed
  )
  
  grimageVector
}
