#' @title Compute a DNAm-based C-Reactive Protein (CRP) Score
#'
#' @description
#' Computes a DNAm-based surrogate proxy for C-Reactive Protein (CRP) levels,
#' a key biomarker for systemic inflammation (inflammaging).
#'
#' @details
#' This function calculates a DNAm-based CRP proxy score from a normalized
#' beta-valued DNAm data matrix. The score is derived by correlating the
#' methylation values with the signs of pre-defined CpG weights from two
#' internal signatures: (i) crp: a vector containing 1765 CpGs and weights
#' associated with CRP as inferred from a large meta-analysis over 20,000 blood
#' samples where both DNAm and CRP was measured, and where no adjustment for
#' variations between memory and naive T-cells was made. (ii) intCRP: contains
#' a 62 CpG CRP-signature that is adjusted for variations between memory and
#' naive T-cells. The resulting output is a relative score that captures the
#' methylation pattern associated with CRP levels.
#'
#' @param betaM A numeric matrix of DNAm beta values (probes as rows). Rows
#' should be Illumina 450k/EPIC CpG identifiers and columns should be samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @return A \code{list} containing two numeric vectors:
#' \code{CRP} and \code{intCRP}.
#'
#' @references
#' Wielscher, M., Mandaviya, P.R., Kuehnel, B. et al.
#' DNA methylation signature of chronic low-grade inflammation and its role in
#' cardio-respiratory diseases.
#' \emph{Nat Commun} 2022
#'
#' @export
#' @importFrom stats sd cor
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' compcrpOut <- compCRP(hannumBmiqM)
compCRP <- function(betaM, minCoverage = 0, verbose = TRUE) {
    crpCpGList <- loadOmniAgeRdata(
        "omniager_crp_cpg",
        verbose = verbose
    )
    scoreList <- list()
    sigNames <- names(crpCpGList)

    for (i in seq_along(crpCpGList)) {
        currentSigName <- sigNames[i]
        sigWeights <- crpCpGList[[i]]

        # 1. Perform coverage check
        coverage <- .checkCpGCoverage(
            betaM = betaM,
            allWeights = sigWeights,
            clockName = currentSigName,
            minCoverage = minCoverage,
            verbose = verbose
        )

        if (!coverage$pass) {
            scoreList[[currentSigName]] <- rep(NA_real_, ncol(betaM))
            next
        }

        # 2. Z-score
        subBeta <- betaM[coverage$betaIdx, , drop = FALSE]

        rowMeansV <- rowMeans(subBeta, na.rm = TRUE)
        rowSdsV <- apply(subBeta, 1, sd, na.rm = TRUE)

        keepIdx <- which(rowSdsV > 0)
        if (length(keepIdx) < length(rowSdsV)) {
            if (verbose) {
                message(sprintf(
                    "[%s] Excluding %d constant probes.",
                    currentSigName, length(rowSdsV) - length(keepIdx)
                ))
            }
            subBeta <- subBeta[keepIdx, , drop = FALSE]
            rowMeansV <- rowMeansV[keepIdx]
            rowSdsV <- rowSdsV[keepIdx]
            currentWeights <- coverage$weightsSubset[keepIdx]
        } else {
            currentWeights <- coverage$weightsSubset
        }

        zMatrix <- (subBeta - rowMeansV) / rowSdsV

        # 3. Perform coverage check
        scores <- as.vector(cor(zMatrix, sign(currentWeights)))
        names(scores) <- colnames(betaM)

        scoreList[[currentSigName]] <- scores
    }

    return(scoreList)
}
