#' @title Calculate DNAmFitAge and related fitness biomarkers
#'
#' @description
#' Calculates all 6 DNAm fitness biomarkers, DNAmFitAge, and FitAgeAcceleration
#' from a DNA methylation matrix and phenotype data.
#
#'
#' @param betaM A numeric matrix. **Rows must be CpGs, Columns must be Samples.**
#'   The `colnames` must be the sample IDs.
#' @param age A numeric vector of chronological age for each sample.
#'   The order **must match the column order** of `betaM`.
#' @param sex A character or factor vector of sex for each sample
#'   (e.g., "Male" or "Female"). The order **must match the column order**
#'   of `beta_matrix`.
#' @param grimageVector A numeric vector of pre-calculated DNAmGrimAge values.
#'   The order **must match the column order** of `betaM`.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages
#'
#' @return
#' A `data.frame` with 12 columns:
#' \itemize{
#'   \item `SampleID`: The sample identifiers.
#'   \item `DNAmGait_noAge`, `DNAmGrip_noAge`, `DNAmVO2max`: Fitness biomarkers.
#'   \item `DNAmGait_wAge`, `DNAmGrip_wAge`, `DNAmFEV1_wAge`: Age-adjusted fitness biomarkers.
#'   \item `DNAmGrimAge`: The input DNAmGrimAge.
#'   \item `DNAmFitAge`: The calculated biological fitness age.
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
#' hannumBmiqM <- hannumExample[[1]]
#' phenoTypesHannum <- hannumExample[[2]]
#' age <- phenoTypesHannum$Age
#' sex <- ifelse(phenoTypesHannum$Sex == "F", "Female", "Male")
#' GrimAge1O <- grimAge1(hannumBmiqM, age, sex)
#' dnamFitAgeOut <- dnamFitAge(hannumBmiqM, age, sex, GrimAge1O$DNAmGrimAge1)
dnamFitAge <- function(betaM, age, sex, grimageVector, minCoverage = 0,
                       verbose = TRUE) {
    # --- 1. Object conversion and validation ---
    DNAmFitnessModels <- loadOmniAgeRdata(
        "omniager_dnamfitage_coef",
        verbose = verbose
    )

    n_samples <- ncol(betaM)
    sampleIds <- colnames(betaM)
    femaleNumeric <- ifelse(sex == "Female", 1, 0)

    # --- 2. Data Preparation with Coverage Check ---
    # Transpose for sample-wise operations
    betaTrans <- t(betaM)

    # Pass minCoverage to the internal prep function
    dataPrep <- .prepareFitAgeData(
        betaM = betaTrans,
        sampleIds = sampleIds,
        femaleVec = femaleNumeric,
        ageVec = age,
        modelData = DNAmFitnessModels,
        minCoverage = minCoverage,
        verbose = verbose
    )


    if (is.null(dataPrep)) {
        res <- data.frame(SampleID = sampleIds)
        res[, c("DNAmFitAge")] <- NA_real_
        return(res)
    }

    fitnessEst <- .estimateFitnessMarkers(dataPrep, DNAmFitnessModels)
    fitnessEst$DNAmGrimAge <- grimageVector
    finalResults <- .calculateFinalFitAge(fitnessEst)

    return(finalResults[match(sampleIds, finalResults$SampleID), ])
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
#' @return Data.frame including DNAmFitAge and FitAgeAccel.
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
