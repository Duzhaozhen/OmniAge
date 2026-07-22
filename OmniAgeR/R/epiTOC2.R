#' @title
#' Estimate epiTOC2 scores
#'
#' @aliases epiTOC2
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and an
#' age vector (optional) and will return the epiTOC2 scores.
#'
#' @param betaM A numeric DNA methylation beta-value matrix with CpG probe
#'   identifiers as row names and samples as columns. CpG identifiers in
#'   \code{rownames(betaM)} are required. Sample identifiers in
#'   \code{colnames(betaM)} are optional when \code{age = NULL}, but are
#'   required when \code{age} is supplied.
#'
#' @param age A numeric vector containing chronological age for each sample.
#'   A named vector is recommended. If names are supplied, they must be
#'   identical to \code{colnames(betaM)}, including sample order. An unnamed
#'   vector is accepted when its length equals \code{ncol(betaM)}; in this
#'   case, values are matched to samples according to the column order of
#'   \code{betaM}, and a warning is issued.
#'
#' @param minCoverage A numeric value between 0 and 1 specifying the minimum
#'   proportion of required CpGs that must be present. Default is 0.5.
#'
#' @param verbose Logical. Whether to print coverage statistics.
#'   Default is \code{TRUE}.
#'
#' @details
#' Building upon a dynamic model of DNA methylation gain in unmethylated
#'  CpG-rich regions, epiTOC2 can directly estimate the cumulative number
#'  of stem cell divisions in a tissue. The details of the algorithm are
#'  described in Teschendorff et al. 2020.
#'
#' @return A list containing the following entries
#'
#' * tnsc: The estimated cumulative number of stem-cell divisions per
#'   stem-cell per year and per sample using the full epiTOC2 model.
#' * tnsc2: The estimated cumulative number of stem-cell divisions per
#'   stem-cell per year and per sample using an approximation of epiTOC2
#'   which assumes all epiTOC2 CpGs have beta-values exactly 0 in the
#'   fetal stage.
#' * irS: This is returned only if the ages are provided, and gives the
#'   estimated average lifetime intrinsic rate of stem-cell division per
#'   sample, as derived from epiTOC2
#' * irS2: As irS, but for the approximation.
#' * irT: The median estimate over all irS values, yielding a median estimate
#'   for the intrinsic rate of stem-cell division for the tissue.
#' * irT2: As irT, but for the approximation.
#' 
#' Sample-level vectors retain \code{colnames(betaM)} when column names are
#' available; otherwise, unnamed numeric vectors are returned.
#'
#' @references
#' Teschendorff AE.
#' A comparison of epigenetic mitotic-like clocks for cancer risk prediction.
#' \emph{Genome Med.} 2020
#'
#' @importFrom stats median
#' 
#' @examples
#' lungInv <- loadOmniAgeRdata(
#'     "omniager_lung_inv",
#'     verbose = FALSE
#' )
#'
#' lungInvM <- lungInv$bmiq_m
#' phenoDf <- lungInv$PhenoTypes
#'
#' age <- setNames(
#'     phenoDf$Age,
#'     colnames(lungInvM)
#' )
#'
#' epitoc2Out <- epiTOC2(
#'     betaM = lungInvM,
#'     age = age
#' )
#'
#' @export
#'

epiTOC2 <- function(betaM, age = NULL, minCoverage = 0.5, verbose = TRUE) {
    
  ageProvided <- !is.null(age)
  
  betaM <- .validateBetaMatrix(
    betaM,
    requireColnames = ageProvided
  )
  
  if (ageProvided) {
    age <- .validateAge(
      age = age,
      sampleNames = colnames(betaM),
      reorder = FALSE
    )
    
    if (any(age <= 0)) {
      stop(
        "`age` must contain values greater than zero because ",
        "epiTOC2 intrinsic division rates are calculated by ",
        "dividing cumulative division estimates by age.",
        call. = FALSE
      )
    }
  }
    
    estParams <- loadOmniAgeRdata(
        "omniager_epitoc2_model",
        verbose = verbose
    )

    dummyWeights <- setNames(seq_len(nrow(estParams)), rownames(estParams))

    # Perform coverage check
    coverageResult <- .checkCpGCoverage(
        betaM = betaM,
        allWeights = dummyWeights,
        clockName = "epiTOC2",
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (!coverageResult$pass) {
      naScores <- setNames(
        rep(NA_real_, ncol(betaM)),
        colnames(betaM)
      )
      
      return(list(
        tnsc = naScores,
        tnsc2 = naScores,
        irS = if (ageProvided) naScores else NULL,
        irS2 = if (ageProvided) naScores else NULL,
        irT = if (ageProvided) NA_real_ else NULL,
        irT2 = if (ageProvided) NA_real_ else NULL
      ))
    }

    # Extract the matching data and parameters
    matchedParams <- estParams[names(coverageResult$weightsSubset), , drop = FALSE]
    subBeta <- betaM[coverageResult$betaIdx, , drop = FALSE]

    deltaV <- as.numeric(matchedParams[, 1])
    beta0V <- as.numeric(matchedParams[, 2])

    # Core algorithm implementation
    # Full Model
    resScale <- 2 / (deltaV * (1 - beta0V))
    tnscV <- colMeans(sweep(subBeta, 1, beta0V, "-") * resScale, na.rm = TRUE)

    # Approximation (Assuming beta0 = 0)
    resScale2 <- 2 / deltaV
    tnsc2V <- colMeans(subBeta * resScale2, na.rm = TRUE)

    names(tnscV) <- colnames(betaM)
    names(tnsc2V) <- colnames(betaM)
    # Intrinsic Rate
    irS <- NULL
    irS2 <- NULL
    irT <- NULL
    irT2 <- NULL

    if (ageProvided) {
      irS <- tnscV / unname(age)
      irS2 <- tnsc2V / unname(age)
      
      names(irS) <- colnames(betaM)
      names(irS2) <- colnames(betaM)
      
      irT <- stats::median(
        irS,
        na.rm = TRUE
      )
      
      irT2 <- stats::median(
        irS2,
        na.rm = TRUE
      )
    }

    return(list(
        tnsc = tnscV,
        tnsc2 = tnsc2V,
        irS = irS,
        irS2 = irS2,
        irT = irT,
        irT2 = irT2
    ))
}
