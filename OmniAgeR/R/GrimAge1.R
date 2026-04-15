#' @title Calculate GrimAge1
#'
#' @description
#' Calculates DNA methylation GrimAge1, a composite biomarker of mortality risk
#' and biological aging.
#'
#' @details
#' This function calculates DNAm GrimAge1 in a multi-step process. First, it
#' predicts DNAm-based surrogate biomarkers for several plasma proteins from
#' the input beta values. These predicted biomarkers, along with chronological
#' age and sex, are then used to calculate a composite mortality risk score.
#' This score is calibrated to the scale of chronological age to produce the
#' final `DNAmGrimAge1`.
#'
#' @param betaM A numeric matrix of DNA methylation beta values. Rows should
#'   represent CpG sites and columns should represent individual samples.
#' @param age A numeric vector of chronological ages for the samples
#'  corresponding to the columns in `betaM`.
#' @param sex A character vector of sample sexes. Must contain "Male" or
#'  "Female" for each sample.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A data.frame containing the following columns:
#' \itemize{
#'   \item `SampleID`: Identifier for each sample.
#'   \item `DNAm...`: Columns for each of the predicted surrogate biomarkers
#'   (e.g., `DNAmADM`, `DNAmGDF15`).
#'   \item `DNAmGrimAge1`: The final calibrated GrimAge1 score.
#' }
#'
#'
#'
#' @references
#' Lu AT, Quach A, Wilson JG, et al.
#' DNA methylation GrimAge strongly predicts lifespan and healthspan
#' \emph{Aging} 2019
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
#' GrimAge1Oout <- grimAge1(betaM = hannumBmiqM, age, sex)
grimAge1 <- function(betaM, age, sex,
                     minCoverage = 0, verbose = TRUE) {
    # 1. Load model weights
    grimage1 <- loadOmniAgeRdata(
        "omniager_grimage1_model",
        verbose = verbose
    )
    # 2. Extract coefficients from model object (grimage1)
    protCoefs <- grimage1[[1]] # CpG weights for proteins
    finalModel <- grimage1[[2]] # Final weights for COX
    calibParams <- grimage1[[3]] # Calibration means/sds
    uqCpgs <- unique(protCoefs$var[startsWith(protCoefs$var, "cg")])
    fakeWeights <- setNames(rep(1, length(uqCpgs)), uqCpgs)
    coverage <- .checkCpGCoverage(betaM, fakeWeights, "GrimAge1", minCoverage, verbose)
    # 3. Handle Covariates (Age and Sex)
    femaleVec <- ifelse(sex == "Female", 1, 0)
    # 4. Phase 1: Predict DNAm Protein Biomarkers
    availableCpGs <- intersect(protCoefs$var, rownames(betaM))
    proteinNames <- unique(protCoefs$Y.pred)
    protPredList <- list()
    for (pName in proteinNames) {
        # Subset coefficients for this specific protein
        pSub <- protCoefs[protCoefs$Y.pred == pName, ]
        # Intersection with available data
        presentVars <- intersect(pSub$var, c(availableCpGs, "Intercept", "Age"))
        pSubValid <- pSub[pSub$var %in% presentVars, ]
        # Calculate score using only present features
        score <- 0
        if ("Intercept" %in% pSubValid$var) {
            score <- pSubValid$beta[pSubValid$var == "Intercept"]
        }
        # Age component (if required)
        if ("Age" %in% pSubValid$var) {
            score <- score + (pSubValid$beta[pSubValid$var == "Age"] * age)
        }
        # CpG component: matrix multiplication of present sites
        cpgVars <- intersect(pSubValid$var, availableCpGs)
        if (length(cpgVars) > 0) {
            # t(betaM) ensures samples are rows for the multiplication
            score <- score + as.vector(t(betaM[cpgVars, , drop = FALSE]) %*%
                pSubValid$beta[match(cpgVars, pSubValid$var)])
        }
        protPredList[[pName]] <- score
    }
    protDf <- as.data.frame(protPredList)
    # 5. Phase 2: Calculate Mortality Risk Score (COX)
    # Feature set: Age, Female, and predicted DNAm Proteins
    finalInput <- cbind(Age = age, Female = femaleVec, protDf)
    finalInput$Intercept <- 1
    # Match variables for the final COX model
    availableFinalVars <- intersect(finalModel$var, colnames(finalInput))
    finalWeightsSub <- finalModel[match(availableFinalVars, finalModel$var), ]
    coxScore <- as.numeric(as.matrix(finalInput[, availableFinalVars]) %*% finalWeightsSub$beta)
    # 6. Phase 3: Calibration to Chronological Age
    coxParams <- calibParams[calibParams$var == "COX", ]
    ageParams <- calibParams[calibParams$var == "Age", ]

    # Formula: DNAmGrimAge = ((COX - MeanCOX)/SdCOX * SdAge) + MeanAge
    zCox <- (coxScore - coxParams$mean) / coxParams$sd
    grimAgeScore <- (zCox * ageParams$sd) + ageParams$mean

    # 7. Final Output
    res <- data.frame(
        SampleID = colnames(betaM),
        protDf,
        DNAmGrimAge1 = grimAgeScore,
        stringsAsFactors = FALSE
    )

    return(res)
}
