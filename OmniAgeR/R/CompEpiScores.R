#' @title Calculate Epigenetic Scores for the Circulating Proteome (EpiScores)
#'
#' @description
#' Computes the 109 validated epigenetic scores (EpiScores) that serve as
#' DNA methylation-based proxies for the levels of circulating plasma proteins,
#' as defined by Gadd et al. (2022).
#'
#' @details
#' This function implements the method from Gadd et al. (2022) to
#' calculate the **109 validated EpiScores** for circulating plasma proteins.
#'
#' **Imputation Handling:**
#' This function features a robust **two-stage imputation process** for
#' handling missing CpG data:
#' 1.  **Sample-level NA Imputation:** Any existing `NA` values within the
#'     provided `betaM` matrix are imputed to the **row-wise mean**
#'     (the mean of that specific CpG across all samples in the input data).
#' 2.  **Missing CpG Imputation:** CpGs required for the scores but
#'     *entirely absent* from `betaM` are added. Their values are imputed
#'     using the **mean beta value from the original training cohort**
#'     (stored in `EpiScoresCoef`).
#'
#'
#' After imputation, the function iterates through each unique protein score,
#' subsetting the appropriate coefficients, and computes a linear predictor
#' using the `calculateLinearPredictor` helper function.
#'
#' @param betaM A numeric matrix of beta values. Rows should be CpG probes and
#' columns should be individual samples.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A list of length 109. Each element of the list corresponds to one
#' calculated EpiScore (one for each protein). The name of each list element is
#' the name of the EpiScore. Each element contains a named numeric vector of
#' the calculated scores for all samples.
#'
#' @references
#' Gadd DA, Hillary RF, McCartney DL, et al.
#' Epigenetic scores for the circulating proteome as tools for
#' disease prediction.
#' \emph{Elife.} 2022
#'
#'
#' @export
#'
#' @examples
#'
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' allEpiscoresOut <- compEpiScores(hannumBmiqM)
compEpiScores <- function(betaM, minCoverage = 0, verbose = TRUE) {
    EpiScoresCoef <- loadOmniAgeRdata(
        "omniager_episcores_coef",
        verbose = verbose
    )

    # 1. Preprocess the NA within the samples (first-stage imputation)
    rowMeansV <- rowMeans(betaM, na.rm = TRUE)
    naRows <- which(rowSums(is.na(betaM)) > 0)
    if (length(naRows) > 0) {
        for (i in naRows) {
            betaM[i, is.na(betaM[i, ])] <- rowMeansV[i]
        }
    }

    resList <- list()
    proteinNames <- unique(EpiScoresCoef$Predictor)

    for (protein in proteinNames) {
        tmpCoef <- EpiScoresCoef[EpiScoresCoef$Predictor == protein, ]


        currentWeights <- setNames(tmpCoef$Coefficient, tmpCoef$CpG_Site)

        #  Perform coverage check
        coverage <- .checkCpGCoverage(
            betaM = betaM,
            allWeights = currentWeights,
            clockName = protein,
            minCoverage = minCoverage,
            verbose = verbose
        )

        if (!coverage$pass) {
            resList[[protein]] <- rep(NA_real_, ncol(betaM))
            next
        }

        # 3. Carry out the second stage interpolation
        requiredCpGs <- names(currentWeights)
        presentCpGs <- names(coverage$weightsSubset)
        missingCpGs <- setdiff(requiredCpGs, presentCpGs)

        if (length(missingCpGs) > 0) {
            trainMeans <- tmpCoef[
                tmpCoef$CpG_Site %in% missingCpGs,
                c("CpG_Site", "Mean_Beta_Value")
            ]

            # Construct the imputation matrix
            imputeMat <- matrix(rep(trainMeans$Mean_Beta_Value, ncol(betaM)),
                ncol = ncol(betaM), byrow = FALSE
            )
            rownames(imputeMat) <- trainMeans$CpG_Site
            colnames(imputeMat) <- colnames(betaM)


            currentBeta <- rbind(betaM[coverage$betaIdx, , drop = FALSE], imputeMat)

            finalWeights <- c(
                coverage$weightsSubset,
                setNames(trainMeans$Mean_Beta_Value, trainMeans$CpG_Site)
            )
        } else {
            currentBeta <- betaM[coverage$betaIdx, , drop = FALSE]
            finalWeights <- coverage$weightsSubset
        }

        # 4. Invoke the computing engine (minCoverage=0)
        resList[[protein]] <- .calculateLinearPredictor(
            betaM = currentBeta,
            coefLv = list(0, finalWeights),
            clockName = protein,
            minCoverage = 0,
            verbose = FALSE
        )
    }

    return(resList)
}
