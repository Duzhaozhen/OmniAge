#' @title Calculate CHIP-related Methylation Scores
#'
#' @description
#' Computes CHIP (Clonal Hematopoiesis of Indeterminate Potential) scores
#' based on predefined CpG signatures and a normalized DNA methylation matrix.
#'
#' @details
#' For each CHIP signature (e.g., DNMT3A, TET2), it identifies common CpGs,
#' calculates Z-scores for each
#' probe across all samples, and then computes the score as:
#' (Mean Z-score of positive probes) - (Mean Z-score of negative probes).
#'
#' Probes with zero variance across samples are excluded.
#'
#' @param betaM A numeric matrix of DNAm beta values (probes as rows). Rows
#' should be Illumina 450k/EPIC CpG identifiers and columns should be samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @return
#' A \code{list} of numeric vectors, one for each CHIP signature.
#'
#' @importFrom stats sd
#'
#' @references
#' Kirmani, S., Huan, T., Van Amburg, J.C. et al.
#' Epigenome-wide DNA methylation association study of CHIP provides insight
#' into perturbed gene regulation.
#' \emph{Nat Commun} 2025
#'
#' @export
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' compchipOut <- compCHIP(hannumBmiqM)
compCHIP <- function(betaM, minCoverage = 0, verbose = TRUE) {
    chipCpGList <- loadOmniAgeRdata(
        "omniager_chip_cpg",
        verbose = verbose
    )

    scoreList <- list()

    for (i in seq_along(chipCpGList)) {
        sigName <- names(chipCpGList)[i]
        sigWeights <- chipCpGList[[i]]

        # 1. Perform coverage check
        coverage <- .checkCpGCoverage(
            betaM = betaM,
            allWeights = sigWeights,
            clockName = sigName,
            minCoverage = minCoverage,
            verbose = verbose
        )

        if (!coverage$pass) {
            scoreList[[sigName]] <- rep(NA_real_, ncol(betaM))
            next
        }

        # 2. Extract the matching data and symbols
        subBeta <- betaM[coverage$betaIdx, , drop = FALSE]
        matchedSigns <- sign(coverage$weightsSubset)

        # 3. Z-score
        rowMeansV <- rowMeans(subBeta, na.rm = TRUE)
        rowSdsV <- apply(subBeta, 1, sd, na.rm = TRUE)

        zeroVarIdx <- which(rowSdsV == 0)
        if (length(zeroVarIdx) > 0) {
            if (verbose) message(sprintf("[%s] Excluding %d constant probes.", sigName, length(zeroVarIdx)))
            subBeta <- subBeta[-zeroVarIdx, , drop = FALSE]
            rowMeansV <- rowMeansV[-zeroVarIdx]
            rowSdsV <- rowSdsV[-zeroVarIdx]
            matchedSigns <- matchedSigns[-zeroVarIdx]
        }

        # Z = (X - Mean) / SD
        zMatrix <- (subBeta - rowMeansV) / rowSdsV

        # 4. Calculate the mean value in groups
        posIdx <- which(matchedSigns == 1)
        negIdx <- which(matchedSigns == -1)

        scoreP <- if (length(posIdx) > 0) colMeans(zMatrix[posIdx, , drop = FALSE]) else 0
        scoreN <- if (length(negIdx) > 0) colMeans(zMatrix[negIdx, , drop = FALSE]) else 0

        # 5. Calculate the final score
        scoreList[[sigName]] <- scoreP - scoreN
    }

    return(scoreList)
}
