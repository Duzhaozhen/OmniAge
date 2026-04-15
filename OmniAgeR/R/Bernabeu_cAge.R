#' @title Calculate Chronological Age (cAge) using the hybrid model
#'
#' @description This function loads pre-trained cAge model weights and applies
#'              them to a user-provided methylation beta-value matrix.
#'              It implements the hybrid model logic from the paper:
#'              1. Predicts age using the standard linear model.
#'              2. If the predicted age is < 20 years, it re-predicts using
#'                 the log(age) model.
#'              It also automatically handles missing CpGs and NA values
#'              via mean imputation.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector of predicted ages. The names of the vector
#'         correspond to the sample IDs (column names) from `betaM`.
#'
#' @references
#' Bernabeu, E., McCartney, D.L., Gadd, D.A. et al.
#' Refining epigenetic prediction of chronological and biological age.
#' \emph{Genome Med.} 2023
#'
#'
#' @export
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' bernabeuCAgeO <- bernabeuCAge(hannumBmiqM)
#' @export

bernabeuCAge <- function(betaM, minCoverage = 0, verbose = TRUE) {
    bernabeuCAgeModel <- loadOmniAgeRdata(
        "omniager_bernabeu_cage_coef",
        verbose = verbose
    )
    bernabeuCAgeCoef <- bernabeuCAgeModel[["Bernabeu_cAge_coef"]]
    bernabeuCAgeCoefLog <- bernabeuCAgeModel[["Bernabeu_cAge_coef_log"]]
    bernabeuCAgeIntercept <- bernabeuCAgeModel[["Bernabeu_cAge_intercept"]]
    bernabeuCAgeInterceptLog <- bernabeuCAgeModel[["Bernabeu_cAge_intercept_log"]]
    bernabeuCAgeMeans <- bernabeuCAgeModel[["Bernabeu_cAge_means"]]

    # 1. Define required features (Linear and Quadratic precursors)
    allRequired <- unique(c(
        gsub("_2", "", rownames(bernabeuCAgeCoef)),
        gsub("_2", "", rownames(bernabeuCAgeCoefLog))
    ))

    # 2. Execute data preprocessing and imputation
    refMeansVec <-
        betaM <- .preprocessEpiClockData(
            betaM = betaM,
            requiredCpGs = allRequired,
            referenceMeans = setNames(
                bernabeuCAgeMeans$mean,
                rownames(bernabeuCAgeMeans)
            ),
            minCoverage = minCoverage,
            filterSamples = FALSE,
            clockName = "bernabeuCAge"
        )


    if (is.null(betaM)) {
        if (verbose) warning("[bernabeuCAge] Prediction aborted due to low coverage.")
        return(setNames(rep(NA_real_, ncol(betaM)), colnames(betaM)))
    }

    # 3. Construct feature matrix (Including quadratic terms)
    cpgsToSquare <- unique(c(
        gsub("_2", "", rownames(bernabeuCAgeCoef)[grep("_2", rownames(bernabeuCAgeCoef))]),
        gsub("_2", "", rownames(bernabeuCAgeCoefLog)[grep("_2", rownames(bernabeuCAgeCoefLog))])
    ))
    betaSq <- betaM[cpgsToSquare, , drop = FALSE]^2
    rownames(betaSq) <- paste0(rownames(betaSq), "_2")
    featureMat <- rbind(betaM, betaSq)

    # 4. Internal predictor closure
    runEngine <- function(fM, coe, icpt, lbl) {
        wts <- setNames(as.numeric(coe[, 1]), rownames(coe))
        .calculateLinearPredictor(fM, list(icpt, wts), lbl, minCoverage = 1, verbose = FALSE)
    }

    # 6. Hybrid model execution
    # Phase A: Primary linear estimation
    predAge <- runEngine(featureMat, bernabeuCAgeCoef, bernabeuCAgeIntercept, "cAgeLinear")
    # Phase B: Log-model correction for samples < 20 years
    under20 <- which(predAge < 20)
    if (length(under20) > 0) {
        predLog <- runEngine(
            featureMat[, under20, drop = FALSE],
            bernabeuCAgeCoefLog, bernabeuCAgeInterceptLog, "cAgeLog"
        )
        predAge[under20] <- exp(predLog)
    }

    return(predAge)
}


#' Universal Preprocessing Engine for Epigenetic Clocks
#'
#' @description
#' A robust, multi-stage preprocessing pipeline that handles:
#' 1. Sample-level NA imputation (Row-mean / Reference-mean).
#' 2. Probe coverage validation via .checkCpGCoverage.
#' 3. Feature-level missing probe imputation using reference baselines.
#' 4. Optional sample-level filtering for low-quality samples.
#'
#' @param betaM Numeric matrix of beta values.
#' @param requiredCpGs Character vector of CpGs required by the model.
#' @param referenceMeans Named numeric vector of training/gold-standard means.
#' @param minCoverage Numeric (0-1). Minimum proportion of probes required.
#' @param filterSamples Logical. If TRUE, removes samples with excessive NAs.
#' @param clockName Character. Name for logging/messages.
#' @param verbose Logical.
#'
#' @return A processed matrix or NULL if the coverage threshold is not met.
#' @keywords internal
#' @noRd
.preprocessEpiClockData <- function(betaM,
                                    requiredCpGs,
                                    referenceMeans,
                                    minCoverage = 0,
                                    filterSamples = FALSE,
                                    clockName = "UnknownClock",
                                    verbose = TRUE) {
    # --- Stage 1: Sample-level NA Imputation (In-sample Row Means) ---
    # We address random missingness first to get an accurate coverage count
    rowMeansV <- rowMeans(betaM, na.rm = TRUE)
    naRows <- which(rowSums(is.na(betaM)) > 0)

    if (length(naRows) > 0) {
        for (i in naRows) {
            # Fallback for all-NA rows: Use reference baseline
            if (is.nan(rowMeansV[i])) {
                cpg <- rownames(betaM)[i]
                rowMeansV[i] <- ifelse(cpg %in% names(referenceMeans),
                    referenceMeans[cpg], 0.5
                )
            }
            betaM[i, is.na(betaM[i, ])] <- rowMeansV[i]
        }
    }

    # --- Stage 2: Coverage Assessment ---
    # Use your standardized helper to check probe overlap
    dummyWeights <- setNames(rep(1, length(requiredCpGs)), requiredCpGs)
    check <- .checkCpGCoverage(betaM, dummyWeights, clockName,
        minCoverage = 0, verbose = verbose
    )

    presentCpGs <- rownames(betaM)[check$betaIdx]
    coverageRatio <- length(presentCpGs) / length(requiredCpGs)

    # --- Stage 3: Quality Threshold & Sample Filtering ---
    if (coverageRatio < minCoverage) {
        if (verbose) {
            warning(sprintf(
                "[%s] Aborted: Coverage (%.1f%%) < threshold.",
                clockName, coverageRatio * 100
            ))
        }
        return(NULL)
    }

    # Optional: Strict sample filtering (Required for DunedinPACE)
    if (filterSamples) {
        # Check missingness proportion per sample across required sites
        # Since we imputed random NAs in Stage 1, this mainly identifies
        # samples with systemically failed probes
        sampleMissingness <- colMeans(is.na(betaM[presentCpGs, , drop = FALSE]))
        keepSamples <- which(sampleMissingness <= (1 - minCoverage))

        if (length(keepSamples) == 0) {
            return(NULL)
        }
        if (length(keepSamples) < ncol(betaM)) {
            if (verbose) {
                message(sprintf(
                    "[%s] Removed %d low-quality samples.",
                    clockName, ncol(betaM) - length(keepSamples)
                ))
            }
            betaM <- betaM[, keepSamples, drop = FALSE]
        }
    }

    # --- Stage 4: Feature-level Imputation (Reference-mean) ---
    # Fill entirely missing probes to complete the linear model
    if (coverageRatio < 1) {
        missingCpGs <- setdiff(requiredCpGs, presentCpGs)
        if (verbose) {
            message(sprintf(
                "[%s] Imputing %d missing probes via reference.",
                clockName, length(missingCpGs)
            ))
        }

        imputeMeans <- referenceMeans[missingCpGs]
        imputeMat <- matrix(rep(imputeMeans, ncol(betaM)),
            ncol = ncol(betaM), byrow = FALSE
        )
        rownames(imputeMat) <- missingCpGs
        colnames(imputeMat) <- colnames(betaM)

        # Merge present data with imputed reference rows
        betaM <- rbind(betaM[presentCpGs, , drop = FALSE], imputeMat)
    } else {
        # Subset exactly to model requirements
        betaM <- betaM[requiredCpGs, , drop = FALSE]
    }

    return(betaM)
}
