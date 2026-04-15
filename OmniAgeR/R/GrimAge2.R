#' @title Calculate GrimAge2
#'
#' @description
#' Calculates DNA methylation GrimAge2, a composite biomarker of mortality risk
#' and biological aging. This function implements the updated GrimAge2 model.
#'
#' @details
#' This function calculates DNAm GrimAge2 in a multi-step process. First, it
#' predicts DNAm-based surrogate biomarkers for several plasma proteins from
#' the input beta values. These predicted biomarkers, along with chronological
#' age and sex, are then used to calculate a composite mortality risk score.
#' This score is calibrated to the scale of chronological age to produce the
#' final `DNAmGrimAge2`.
#'
#' @param betaM A numeric matrix of DNA methylation beta values. Rows should
#'   represent CpG sites and columns should represent individual samples.
#' @param age A numeric vector of chronological ages for the samples corresponding
#'   to the columns in `betaM`.
#' @param sex A character vector of sample sexes. Must contain "Male" or "Female"
#'   for each sample.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @return
#' A data.frame containing the following columns:
#' \itemize{
#'   \item `SampleID`: Identifier for each sample.
#'   \item `DNAm...`: Columns for each of the predicted surrogate biomarkers
#'   (e.g., `DNAmADM`, `DNAmGDF15`).
#'   \item `DNAmGrimAge2`: The final calibrated GrimAge2 score.
#' }
#'
#'
#'
#' @references
#' Lu AT, Binder AM, Zhang J, et al.
#' DNA methylation GrimAge version 2.
#' \emph{Aging} 2022
#'
#' @export
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
#' GrimAge2Oout <- grimAge2(betaM = hannumBmiqM, age, sex)
grimAge2 <- function(betaM, age, sex,
                     minCoverage = 0, verbose = TRUE) {
    # 1. Load the model file
    grimage2 <- loadOmniAgeRdata(
        "omniager_grimage2_model",
        verbose = verbose
    )
    # 2. Extract coefficients from model object (grimage2)
    protCoefs <- grimage2[[1]] # CpG weights for proteins/biomarkers
    finalModel <- grimage2[[2]] # Final weights for COX
    calibParams <- grimage2[[3]] # Calibration parameters (gold standard)

    uqCpgs <- unique(protCoefs$var[startsWith(protCoefs$var, "cg")])
    fakeWeights <- setNames(rep(1, length(uqCpgs)), uqCpgs)
    coverage <- .checkCpGCoverage(betaM, fakeWeights, "GrimAge2", minCoverage, verbose)
    # 3. Handle Covariates
    femaleVec <- ifelse(sex == "Female", 1, 0)

    # 4. Phase 1: Predict Surrogate Biomarkers
    # Identify available CpGs (no imputation used)
    availableCpGs <- intersect(protCoefs$var, rownames(betaM))

    if (verbose) {
        nTotalCpGs <- length(setdiff(unique(protCoefs$var), c("Intercept", "Age")))
        message(sprintf(
            "[GrimAge2] Found %d / %d required CpGs (%.1f%%).",
            length(availableCpGs), nTotalCpGs, (length(availableCpGs) / nTotalCpGs) * 100
        ))
    }

    proteinNames <- unique(protCoefs$Y.pred)
    protPredList <- list()

    # Optimized matrix multiplication for protein prediction
    for (pName in proteinNames) {
        pSub <- protCoefs[protCoefs$Y.pred == pName, ]
        presentVars <- intersect(pSub$var, c(availableCpGs, "Intercept", "Age"))
        pSubValid <- pSub[pSub$var %in% presentVars, ]

        # Calculate scores based on original logic: sum of available features
        # Start with Intercept
        score <- if ("Intercept" %in% pSubValid$var) pSubValid$beta[pSubValid$var == "Intercept"] else 0

        # Add Age effect if model requires it
        if ("Age" %in% pSubValid$var) {
            score <- score + (pSubValid$beta[pSubValid$var == "Age"] * age)
        }

        # Add CpG effect (matrix multiplication)
        cpgVars <- intersect(pSubValid$var, availableCpGs)
        if (length(cpgVars) > 0) {
            # Use crossprod or %*% for efficiency. t(betaM) matches samples to rows.
            score <- score + as.vector(t(betaM[cpgVars, , drop = FALSE]) %*%
                pSubValid$beta[match(cpgVars, pSubValid$var)])
        }

        protPredList[[pName]] <- score
    }

    protDf <- as.data.frame(protPredList)

    # 5. Phase 2: Calculate Combined Risk Score (COX)
    # Features: Age, Female, and predicted DNAm Biomarkers
    finalInput <- cbind(Age = age, Female = femaleVec, protDf)
    finalInput$Intercept <- 1

    # Strictly follow original variable matching
    availableFinalVars <- intersect(finalModel$var, colnames(finalInput))
    finalWeightsSub <- finalModel[match(availableFinalVars, finalModel$var), ]

    coxScore <- as.numeric(as.matrix(finalInput[, availableFinalVars]) %*% finalWeightsSub$beta)

    # 6. Phase 3: Calibration to Chronological Age
    #
    coxParams <- calibParams[calibParams$var == "COX", ]
    ageParams <- calibParams[calibParams$var == "Age", ]

    zCox <- (coxScore - coxParams$mean) / coxParams$sd
    grimAge2Score <- (zCox * ageParams$sd) + ageParams$mean

    # 7. Final Formatting and Renaming
    res <- data.frame(
        SampleID = colnames(betaM),
        protDf,
        DNAmGrimAge2 = grimAge2Score,
        stringsAsFactors = FALSE
    )

    # Rename according to original requirements
    renameMap <- c(
        "DNAmadm" = "DNAmADM", "DNAmCystatin_C" = "DNAmCystatinC",
        "DNAmGDF_15" = "DNAmGDF15", "DNAmleptin" = "DNAmLeptin",
        "DNAmpai_1" = "DNAmPAI1", "DNAmTIMP_1" = "DNAmTIMP1",
        "DNAmlog.CRP" = "DNAmlogCRP", "DNAmlog.A1C" = "DNAmlogA1C"
    )

    names(res) <- ifelse(names(res) %in% names(renameMap),
        renameMap[names(res)], names(res)
    )

    return(res)
}
