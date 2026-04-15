#' @title Estimate epiTOC3 scores
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and
#' an age vector (optional) and will return the epiTOC3 scores.
#'
#' @param betaM A numeric matrix of DNAm beta values (probes as rows).
#' @param age Optional numeric vector representing chronological ages
#'   of the samples.
#' @param minCoverage Numeric (0-1). Minimum required probe coverage.
#'   Default is 0.
#' @param verbose Logical. Whether to print coverage statistics.
#'
#' @details
#' Building upon a dynamic model of DNA methylation gain in 170 unmethylated
#' population doubling associated CpGs, epiTOC3 can directly estimate the
#' cumulative number of stem cell divisions in a tissue.
#'
#' @return A list containing the following entries
#'
#' * tnsc: The estimated cumulative number of stem-cell divisions per stem-cell
#'   per year and per sample using the full epiTOC3 model.
#' * tnsc2: The estimated cumulative number of stem-cell divisions per stem-cell
#'   per year and per sample using an approximation of epiTOC3 which assumes all
#'   epiTOC3 CpGs have beta-values exactly 0 in the fetal stage.
#' * irS: This is returned only if the ages are provided, and gives the
#'   estimated average lifetime intrinsic rate of stem-cell division per
#'   sample, as derived from epiTOC3
#' * irS2: As irS, but for the approximation.
#' * irT: The median estimate over all irS values, yielding a median estimate
#'   for the intrinsic rate of stem-cell division for the tissue.
#' * irT2: As irT, but for the approximation.
#' * avETOC3: The simple average over the 170 epiTOC3 sites.
#'
#'
#' @importFrom stats median
#'
#' @examples
#' lungInv <- loadOmniAgeRdata(
#'     "omniager_lung_inv",
#'     verbose = FALSE
#' )
#' lungInvM <- lungInv$bmiq_m
#' phenoDf <- lungInv$PhenoTypes
#' epitoc3Out <- epiTOC3(betaM = lungInvM, age = phenoDf$Age)
#'
#' @export
#'

epiTOC3 <- function(betaM, age = NULL, minCoverage = 0, verbose = TRUE) {
    estParams <- loadOmniAgeRdata(
        "omniager_epitoc3_model",
        verbose = verbose
    )

    dummyWeights <- setNames(seq_len(nrow(estParams)), rownames(estParams))

    # Perform coverage check
    coverageResult <- .checkCpGCoverage(
        betaM = betaM,
        allWeights = dummyWeights,
        clockName = "epiTOC3",
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (!coverageResult$pass) {
        return(list(
            tnsc = rep(NA_real_, ncol(betaM)), tnsc2 = rep(NA_real_, ncol(betaM)),
            irS = if (!is.null(age)) rep(NA_real_, ncol(betaM)) else NULL,
            irS2 = if (!is.null(age)) rep(NA_real_, ncol(betaM)) else NULL,
            irT = NA_real_, irT2 = NA_real_,
            avETOC3 = rep(NA_real_, ncol(betaM))
        ))
    }

    # Extract the matching data and parameters
    matchedParams <- estParams[names(coverageResult$weightsSubset), , drop = FALSE]
    subBeta <- betaM[coverageResult$betaIdx, , drop = FALSE]

    deltaV <- matchedParams[, 1]
    beta0V <- matchedParams[, 2]

    # Core algorithm implementation
    avETOC3 <- colMeans(subBeta, na.rm = TRUE)

    # Full Model: 2 * Mean( (Beta - Beta0) / (delta * (1 - Beta0)) )
    scalingFactor <- 2 / (deltaV * (1 - beta0V))
    tnscV <- colMeans(sweep(subBeta, 1, beta0V, "-") * scalingFactor, na.rm = TRUE)

    # Approximation (beta0 = 0)
    scalingFactor2 <- 2 / deltaV
    tnsc2V <- colMeans(subBeta * scalingFactor2, na.rm = TRUE)

    # Intrinsic Rate
    irS <- NULL
    irS2 <- NULL
    irT <- NULL
    irT2 <- NULL

    if (!is.null(age)) {
        irS <- tnscV / age
        irS2 <- tnsc2V / age
        irT <- median(irS, na.rm = TRUE)
        irT2 <- median(irS2, na.rm = TRUE)
    }

    return(list(
        tnsc = tnscV,
        tnsc2 = tnsc2V,
        irS = irS,
        irS2 = irS2,
        irT = irT,
        irT2 = irT2,
        avETOC3 = avETOC3
    ))
}
