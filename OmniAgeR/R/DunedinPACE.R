#' @title Calculate DunedinPACE
#'
#' @description
#' Calculates the DunedinPACE, an epigenetic biomarker that quantifies the pace of biological aging from a single blood sample.
#'
#' @details
#' This function calculates DunedinPACE scores from a DNA methylation beta value matrix using a robust pipeline.
#'
#' The process involves several key steps:
#' * **EPICv2 Array Handling**: Automatically detects and preprocesses Illumina EPICv2 array data, adjusting for its specific characteristics.
#' * **Probe Overlap Check**: Verifies that a sufficient number of required CpG probes are present in your data.
#' * **Missing Data Imputation**: Robustly handles both missing probes and missing values within existing probes.
#' * **Quantile Normalization**: Standardizes the data against a gold-standard reference to reduce batch effects and ensure comparability.
#' * **Score Calculation**: Computes the final DunedinPACE score using the pre-trained model weights.
#'
#' @param betaM A numeric matrix (Rows: CpGs, Cols: Samples)
#' @param minCoverage Numeric (0-1). Minimum required probe coverage for
#'   imputation and calculation. Default is 0.
#' @param verbose Logical. Whether to print progress messages.
#'   Default is \code{TRUE}.
#'
#' @return  A named numeric vector containing the calculated DunedinPACE
#' for each sample.
#'
#' @references
#' Belsky DW, Caspi A, Corcoran DL, et al.
#' DunedinPACE, a DNA methylation biomarker of the pace of aging.
#' \emph{eLife} 2022
#'
#' @export
#' @importFrom preprocessCore normalize.quantiles.use.target
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' dunedinPACEOut <- dunedinPACE(hannumBmiqM)
dunedinPACE <- function(betaM, minCoverage = 0, verbose = TRUE) {
    # Load model data
    modelName <- "DunedinPACE"
    modelSpecs <- loadOmniAgeRdata(
        "omniager_dunedinpace_model",
        verbose = verbose
    )
    # --- 1. Dual-Level Coverage Disclosure ---
    if (verbose) {
        # A. Check direct Model CpGs (the 173 probes used for weighted sum)
        modelProbes <- names(modelSpecs$model_weights[[modelName]])
        nModelFound <- sum(modelProbes %in% rownames(betaM))
        message(sprintf(
            "[%s] Model CpGs: Found %d / %d (%.1f%%)",
            modelName, nModelFound, length(modelProbes),
            (nModelFound / length(modelProbes)) * 100
        ))
    }

    # --- 2. Execute Universal Preprocessing ---
    # We pass 'gold_standard_probes' as the required set because they are
    # necessary for the subsequent Quantile Normalization step.
    processedMat <- .preprocessEpiClockData(
        betaM = betaM,
        requiredCpGs = modelSpecs$gold_standard_probes[[modelName]],
        referenceMeans = modelSpecs$gold_standard_means[[modelName]],
        minCoverage = minCoverage,
        filterSamples = TRUE,
        clockName = modelName,
        verbose = verbose # This will print the Background CpG stats
    )

    if (is.null(processedMat)) {
        if (verbose) warning("[DunedinPACE] Prediction aborted: Insufficient probe coverage.")
        return(setNames(rep(NA_real_, ncol(betaM)), colnames(betaM)))
    }

    # --- 4. Quantile Normalization to Gold Standard ---
    #
    if (verbose) message("[DunedinPACE] Normalizing to gold-standard reference...")
    targetMeans <- modelSpecs$gold_standard_means[[modelName]]


    betasNorm <- preprocessCore::normalize.quantiles.use.target(
        processedMat,
        target = targetMeans
    )
    rownames(betasNorm) <- rownames(processedMat)
    colnames(betasNorm) <- colnames(processedMat)

    # --- 5. Score Calculation ---
    weights <- modelSpecs$model_weights[[modelName]]
    intercept <- modelSpecs$model_intercept[[modelName]]


    finalScores <- .calculateLinearPredictor(
        betaM = betasNorm,
        coefLv = list(intercept, weights),
        clockName = modelName,
        minCoverage = 1,
        verbose = FALSE
    )

    # --- 6. Re-align with original samples ---
    fullResults <- setNames(rep(NA_real_, ncol(betaM)), colnames(betaM))
    fullResults[names(finalScores)] <- finalScores

    return(fullResults)
}
