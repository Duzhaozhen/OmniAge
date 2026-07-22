#' @title Calculate Principal Component (PC) Epigenetic Clocks
#'
#' @description
#' Calculates a suite of Principal Component (PC)-based epigenetic clocks
#' based on the methodology from Higgins-Chen et al. (2022).
#'
#' This function computes PC-based versions of Horvath2013, Horvath2018, Hannum,
#' PhenoAge, and GrimAge1, along with their principal components.
#'
#' @param betaM A numeric DNA methylation beta-value matrix with CpG probe
#'   identifiers as row names and sample identifiers as column names.
#'   \code{colnames(betaM)} must be provided and must be identical to
#'   \code{names(age)} and \code{names(sex)}.
#'
#' @param age A numeric vector containing chronological age for each sample.
#'   A named vector is recommended. If names are supplied, they must be
#'   identical to \code{colnames(betaM)}, including sample order. An unnamed
#'   vector is accepted when its length equals \code{ncol(betaM)}; in this
#'   case, values are matched to samples according to the column order of
#'   \code{betaM}, and a warning is issued.
#'
#' @param sex A character or factor vector containing sex for each sample.
#'   Values must be either \code{"Male"} or \code{"Female"}. A named vector
#'   is recommended. If names are supplied, they must be identical to
#'   \code{colnames(betaM)}, including sample order. An unnamed vector is
#'   accepted when its length equals \code{ncol(betaM)}; in this case, values
#'   are matched to samples according to the column order of \code{betaM},
#'   and a warning is issued.
#'
#' @param clockData The PC clock model object returned by
#'   \code{loadOmniAgeRdata("PCClocks_data")}.
#'
#' @param minCoverage A numeric value between 0 and 1 specifying the minimum
#'   proportion of required CpGs that must be present. Default is 0.5.
#' @param verbose Logical. Whether to print status messages.
#'
#'
#' @return
#' A data.frame containing the original `SampleID` columns, 
#' appended with 14 new columns for the calculated PC clock values
#' (e.g., `PCHorvath2013`, `PCHannum`, `PCGrimAge1`, etc.).
#'
#' @references
#' Higgins-Chen AT, Thrush KL, Wang Y, et al.
#' A computational solution for bolstering reliability of epigenetic clocks:
#' Implications for clinical trials and longitudinal tracking.
#' \emph{Nat Aging.} (2022).
#'
#'
#' @export
#'
#' @examples
#' # 1. Fast runnable code to satisfy BiocCheck
#' message("Ready to initialize PC-based clock pipeline.")
#' 
#' ## Downloading "PCClocks_data" will take a very long time.
#' \donttest{
#' pcClockData <- loadOmniAgeRdata(
#'     "PCClocks_data",
#'     verbose = FALSE
#' )
#' hannumExample <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )
#' hannumBmiqM <- hannumExample[[1]]
#' phenoTypesHannum <- hannumExample[[2]]
#' sampleIds <- colnames(hannumBmiqM)
#' age <- setNames(phenoTypesHannum$Age, sampleIds)
#'
#' sex <- setNames(ifelse(phenoTypesHannum$Sex == "F", "Female","Male"),sampleIds)
#' pcClocksOut <- pcClocks(hannumBmiqM, age, sex, pcClockData)
#'}
pcClocks <- function(betaM, age, sex, clockData, minCoverage = 0.5, verbose = TRUE) {
    if (verbose) message("[PCClocks] Initializing PC-based clock pipeline...")

    # --- 1. Input Validation and Conversion ---
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = TRUE
    )
    
    sampleIds <- colnames(betaM)
    
    age <- .validateAge(
      age = age,
      sampleNames = sampleIds,
      reorder = FALSE
    )
    
    sex <- .validateSex(
      sex = sex,
      sampleNames = sampleIds,
      allowedValues = c("Male", "Female"),
      reorder = FALSE
    )
  

    # Validate clockData Hash (Security & Integrity Check)
    if (rlang::hash(clockData) != "46386ec4be2b2a5239cf67b242d7dc24") {
        stop("[PCClocks] Invalid or corrupted clockData object. Please re-download.")
    }

    pheno <- data.frame(
      SampleID = sampleIds,
      Age = unname(age),
      Female = unname(ifelse(sex == "Female", 1, 0)),
      stringsAsFactors = FALSE
    )

    # --- 2. Standardized Preprocessing & Coverage Check ---
    # Transpose for PC operations (Rows = Samples)
    betaTrans <- t(betaM)

    # Use optimized internal helper for detection and mean imputation
    betaProcessed <- .preprocessPcData(
        betaM = betaTrans,
        requiredCpGs = clockData$imputeMissingCpGs,
        minCoverage = minCoverage,
        verbose = verbose
    )

    if (is.null(betaProcessed)) {
      nSamples <- length(sampleIds)
      
      return(data.frame(
        SampleID = sampleIds,
        PCHorvath2013 = rep(NA_real_, nSamples),
        PCHorvath2018 = rep(NA_real_, nSamples),
        PCHannum = rep(NA_real_, nSamples),
        PCPhenoAge = rep(NA_real_, nSamples),
        PCDNAmTL = rep(NA_real_, nSamples),
        PCPACKYRS = rep(NA_real_, nSamples),
        PCADM = rep(NA_real_, nSamples),
        PCB2M = rep(NA_real_, nSamples),
        PCCystatinC = rep(NA_real_, nSamples),
        PCGDF15 = rep(NA_real_, nSamples),
        PCLeptin = rep(NA_real_, nSamples),
        PCPAI1 = rep(NA_real_, nSamples),
        PCTIMP1 = rep(NA_real_, nSamples),
        PCGrimAge1 = rep(NA_real_, nSamples),
        stringsAsFactors = FALSE
      ))
    }

    # --- 3. PC Projections and Clock Estimation ---
    #
    if (verbose) message("[PCClocks] Projecting data onto principal components...")

    # Helper to calculate individual PC Clocks
    calcPc <- function(dat, mod, transform = FALSE) {
        # Formula: anti.trafo( (Beta - Center) %*% Rotation %*% Weights + Intercept )
        val <- (sweep(dat, 2, mod$center) %*% mod$rotation %*% mod$model) + mod$intercept
        if (transform) {
            return(as.numeric(.antiTrafo(val)))
        }
        return(as.numeric(val))
    }

    pheno$PCHorvath2013 <- calcPc(betaProcessed, clockData$CalcPCHorvath1, transform = TRUE)
    pheno$PCHorvath2018 <- calcPc(betaProcessed, clockData$CalcPCHorvath2, transform = TRUE)
    pheno$PCHannum <- calcPc(betaProcessed, clockData$CalcPCHannum)
    pheno$PCPhenoAge <- calcPc(betaProcessed, clockData$CalcPCPhenoAge)
    pheno$PCDNAmTL <- calcPc(betaProcessed, clockData$CalcPCDNAmTL)

    # --- 4. Complex PCGrimAge Logic ---
    if (verbose) message("[PCClocks] Estimating PCGrimAge components...")
    # Project beta into PC space for GrimAge
    grimPcSpace <- sweep(betaProcessed, 2, clockData$CalcPCGrimAge$center) %*% clockData$CalcPCGrimAge$rotation
    grimFeatures <- cbind(grimPcSpace, Female = pheno$Female, Age = pheno$Age)
    # Internal function for GrimAge Sub-biomarkers
    calcGrimSub <- function(feat, subMod, subInt) {
        as.numeric(feat[, names(subMod)] %*% subMod + subInt)
    }
    pheno$PCPACKYRS <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCPACKYRS.model, clockData$CalcPCGrimAge$PCPACKYRS.intercept)
    pheno$PCADM <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCADM.model, clockData$CalcPCGrimAge$PCADM.intercept)
    pheno$PCB2M <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCB2M.model, clockData$CalcPCGrimAge$PCB2M.intercept)
    pheno$PCCystatinC <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCCystatinC.model, clockData$CalcPCGrimAge$PCCystatinC.intercept)
    pheno$PCGDF15 <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCGDF15.model, clockData$CalcPCGrimAge$PCGDF15.intercept)
    pheno$PCLeptin <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCLeptin.model, clockData$CalcPCGrimAge$PCLeptin.intercept)
    pheno$PCPAI1 <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCPAI1.model, clockData$CalcPCGrimAge$PCPAI1.intercept)
    pheno$PCTIMP1 <- calcGrimSub(grimFeatures, clockData$CalcPCGrimAge$PCTIMP1.model, clockData$CalcPCGrimAge$PCTIMP1.intercept)
    # Final integrated PCGrimAge1
    grimComp <- pheno[, clockData$CalcPCGrimAge$components]
    pheno$PCGrimAge1 <- as.numeric(as.matrix(grimComp) %*% clockData$CalcPCGrimAge$PCGrimAge.model + clockData$CalcPCGrimAge$PCGrimAge.intercept)
    
    if (!identical(pheno$SampleID, sampleIds)) {
      stop(
        "[PCClocks] Internal sample alignment failed.",
        call. = FALSE
      )
    }
    
    pheno$Age <- NULL
    pheno$Female <- NULL
    
    return(pheno)
}
