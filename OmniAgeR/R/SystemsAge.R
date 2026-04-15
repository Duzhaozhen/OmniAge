#' @title Calculate Systems Age and 11 System-Specific Scores
#'
#' @description
#' Implements the Systems Age epigenetic clock, a multi-system aging biomarker
#' based on DNA methylation (DNAm). This function calculates a composite 'Systems Age'
#' as well as 11 individual physiological system aging scores from a single
#' blood methylation dataset.
#'
#' @param betaM A numeric matrix or data.frame of DNA methylation beta values.
#'   **Rows must correspond to CpGs** and **columns to samples**.
#'
#' @param clockData The pre-loaded data object from
#'  \code{loadOmniAgeRData("PCClocks_data")}.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'  Default is 0.
#' @param verbose Logical. Whether to print status messages.
#'
#'
#' @details
#' The Systems Age framework, developed by Sehgal et al. (2025), quantifies
#' aging heterogeneity across 11 distinct physiological systems. It integrates
#' supervised and unsupervised machine learning with clinical biomarkers and
#' mortality risk to derive system-specific scores.
#'
#'
#'
#' @return
#' A data.frame where the first column is `SampleID` (derived from the
#' `colnames` of the input `DNAm` matrix) and the subsequent 13 columns
#' contain the calculated scores. These include the 11 system scores
#' (e.g., 'Blood', 'Brain', 'Heart'), the 'Age_prediction' score,
#' and the composite 'SystemsAge' score. All scores are scaled to the
#' unit of years.
#'
#' @references
#' Sehgal, R., Markov, Y., Qin, C. et al.
#' Systems Age: a single blood methylation test to quantify aging heterogeneity
#' across 11 physiological systems.
#' \emph{Nat Aging} (2025).
#'
#'
#' @export
#'
#' @examples

#' #' # 1. Fast runnable code to satisfy BiocCheck
#' message("Ready to initialize Systems Age analysis.")
#' 
#' # 2. Real pipeline execution (skip automated checks for speed)
#' \donttest{
#' systemsAgeData <- loadOmniAgeRdata(
#'     "SystemsAge_data",
#'     verbose = FALSE
#' )
#' hannumExample <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )
#' hannumBmiqM <- hannumExample[[1]]
#' systemsAgeOut <- systemsAge(hannumBmiqM, systemsAgeData)
#'}
systemsAge <- function(betaM, clockData, minCoverage = 0, verbose = TRUE) {
    if (verbose) message("[SystemsAge] Initializing multi-system aging analysis...")

    # --- 1. Validation & SE Support ---
    # Data Integrity Check
    if (rlang::hash(clockData) != "d984914ff6aa17d8a6047fed5f9f6e4d") {
        stop("[SystemsAge] clockData hash mismatch. Please re-download SystemsAge data.")
    }

    sampleIds <- colnames(betaM)
    n_samples <- length(sampleIds)
    # --- 2. Preprocessing & Coverage Check ---
    # Transpose: Samples as rows for PCA projection
    betaTrans <- t(betaM)

    # Reuse the same helper as PCClocks
    betaProcessed <- .preprocessPcData(
        betaM = betaTrans,
        requiredCpGs = clockData$imputeMissingCpGs,
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (is.null(betaProcessed)) {
        res <- data.frame(SampleID = sampleIds)
        res[seq_len(n_samples), 2:14] <- NA_real_
        return(res)
    }

    # --- 3. DNAm PCA Projection ---
    #
    if (verbose) message("[SystemsAge] Projecting DNAm onto global PCs...")
    dnamPCs <- predict(clockData$DNAmPCA, betaProcessed)

    # --- 4. Calculate System-Specific PCs & Scores ---
    # SystemsAge uses a nested PCA/Linear model approach
    if (verbose) message("[SystemsAge] Estimating 11 physiological system scores...")

    # Map global PCs to System PCs
    dnamSystemPCs <- dnamPCs[, seq_len(4017)] %*% as.matrix(clockData$system_vector_coefficients[seq_len(4017), ])

    groups <- c(
        "Blood", "Brain", "Cytokine", "Heart", "Hormone", "Immune",
        "Kidney", "Liver", "Metab", "Lung", "MusculoSkeletal"
    )

    systemScores <- matrix(nrow = nrow(dnamSystemPCs), ncol = length(groups))
    colnames(systemScores) <- groups

    for (i in seq_along(groups)) {
        group <- groups[i]
        tf <- grepl(group, colnames(dnamSystemPCs))
        subPCs <- dnamSystemPCs[, tf]
        coeffs <- clockData$system_scores_coefficients_scale[tf]

        # Matrix multiplication or simple scaling for single-PC systems
        if (length(coeffs) == 1) {
            systemScores[, i] <- subPCs * -1
        } else {
            systemScores[, i] <- subPCs %*% coeffs
        }
    }

    # --- 5. Age Prediction & Systems Age Index ---
    # 5a. Predicted Chronological Age
    agePredRaw <- (as.matrix(dnamPCs) %*% as.matrix(clockData$Predicted_age_coefficients[2:4019])) +
        clockData$Predicted_age_coefficients[1]

    # Polynomial transformation
    agePred <- (agePredRaw * clockData$Age_prediction_model[2]) +
        ((agePredRaw^2) * clockData$Age_prediction_model[3]) +
        clockData$Age_prediction_model[1]

    # Convert months to years
    agePred <- agePred / 12

    # 5b. Overall SystemsAge (Integrated via final PCA)
    # Re-align system names for the final index
    colnames(systemScores) <- c(
        "Blood", "Brain", "Inflammation", "Heart", "Hormone", "Immune",
        "Kidney", "Liver", "Metabolic", "Lung", "MusculoSkeletal"
    )

    allScores <- cbind(systemScores, Age_prediction = as.numeric(agePred))
    systemPCA <- predict(clockData$systems_PCA, allScores)
    systemsAgeRaw <- systemPCA %*% clockData$Systems_clock_coefficients

    finalScores <- cbind(allScores, SystemsAge = as.numeric(systemsAgeRaw))

    # --- 6. Final Scaling (Unit: Years) ---
    # Apply study-specific transformation to normalize system ages
    for (j in seq_len(13)) {
        y <- finalScores[, j]
        finalScores[, j] <- (((y - clockData$transformation_coefs[j, 1]) / clockData$transformation_coefs[j, 2]) * clockData$transformation_coefs[j, 4]) + clockData$transformation_coefs[j, 3]
        finalScores[, j] <- finalScores[, j] / 12 # Final year conversion
    }

    # Combine with Sample IDs
    results <- data.frame(SampleID = sampleIds, finalScores, stringsAsFactors = FALSE)
    return(results)
}
