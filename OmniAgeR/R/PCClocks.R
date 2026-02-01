#' @title Calculate Principal Component (PC) Epigenetic Clocks
#'
#' @description
#' Calculates a suite of Principal Component (PC)-based epigenetic clocks
#' based on the methodology from Higgins-Chen et al. (2022).
#'
#' This function computes PC-based versions of Horvath2013, Horvath2018, Hannum,
#' PhenoAge, and GrimAge1, along with their principal components.
#'
#' @param DNAm A numeric matrix of DNA methylation (beta) values.
#'   **Samples must be in columns** and CpG probes in rows.
#' @param age_v A numeric vector of chronological ages for each sample,
#'   in the same order as the columns of `DNAm`.
#' @param sex_v A character vector of biological sex for each sample, in the
#'   same order as the columns of `DNAm`. Values of "Female" are
#'   encoded as 1; all other values are encoded as 0.
#' @param RData A **required** argument specifying the clock data. This can be:
#'   \itemize{
#'     \item A `character` string: The path to the *directory* containing the
#'       `PCClocks_data.qs2` file.
#'     \item A `list`: The pre-loaded data object returned from
#'       `load_OmniAgeR_data(object_name = "PCClocks_data")`.
#'   }
#'
#' @details
#' PCClocks calculation requires the `PCClocks_data.qs2` object, which can be
#' downloaded using [download_OmniAgeR_data()].
#'
#' The `RData` argument must be provided as either the path to the data
#' directory (e.g., from `get_OmniAgeR_path()`) or the pre-loaded list
#' object. The advantage of passing the pre-loaded object is that the data
#' will not be re-loaded on each function call, improving performance when
#' running the function multiple times (e.g., in a loop).
#'
#' @return
#' A data.frame containing the original `Sample_ID`, `Age`, and `Female`
#' columns, appended with 14 new columns for the calculated PC clock values
#' (e.g., `PCHorvath2013`, `PCHannum`, `PCGrimAge1`, etc.).
#'
#' @references
#' Higgins-Chen AT, Thrush KL, Wang Y, et al.
#' A computational solution for bolstering reliability of epigenetic clocks: Implications for clinical trials and longitudinal tracking.
#' \emph{Nat Aging.} (2022).
#'
#'
#' @export
#'
#' @examples
#'
#' # Download the external data
#' download_OmniAgeR_data(clocks = "PCClocks") #  ZENODO_DOI: "10.5281/zenodo.17162604"
#'
#' # Either path to the data
#' RData <- get_OmniAgeR_path()
#' # OR
#' RData <- load_OmniAgeR_data(object_name = "PCClocks_data")
#'
#'
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' age <- PhenoTypesHannum_lv$Age
#' sex <- ifelse(PhenoTypesHannum_lv$Sex=="F","Female","Male")
#' PCClocks.o <- PCClocks(hannum_bmiq_m,age,sex,RData = RData)
#'


PCClocks <- function(DNAm,age_v,sex_v,RData) {


  pheno <- data.frame(Sample_ID = colnames(DNAm),Age=age_v,
                      Female= ifelse(sex_v == "Female",1,0 ))
  DNAm <- t(DNAm)

  # Check RData
  is_valid_path <- checkmate::test_character(RData, len = 1, any.missing = FALSE)
  is_valid_list <- checkmate::test_list(RData, any.missing = FALSE)

  if (!is_valid_path && !is_valid_list) {
    stop(
      paste("[PCClocks] RData argument must be a valid file path (character) or a pre-loaded data object (list).",
            "It cannot be NULL. Please load the data first.",
            "See `?download_OmniAgeR_data` to download the required files.")
    )
  }



  # handle RData
  # if (is.null(RData)) {
  #   RData <- load_PCClocks_data(object_name = "PCClocks_data")
  # } else if (is.character(RData)) {
  #   RData <- load_PCClocks_data(object_name = "PCClocks_data", path = RData)
  # }
  if (is.character(RData)) {
    RData <- load_OmniAgeR_data(object_name = "PCClocks_data", path = RData)
  }

  if (rlang::hash(RData) != "46386ec4be2b2a5239cf67b242d7dc24") {
    stop("[PCClocks] The downloaded PCClocks data is corrupted or the wrong data (e.g., SystemsAge) was passed. See `?download_methylCIPHER()`.")
  }


  ## Imputation
  DNAm <- PCClocks_impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = RData$imputeMissingCpGs,
    subset = TRUE
  )

  ## Re-align to make sure things lined up with the object
  DNAm <- DNAm[, names(RData$imputeMissingCpGs), drop = F]

  #print("[PCClocks] Calculating PC Clocks now")

  # Calculate PC Clocks
  pheno$PCHorvath2013 <- as.numeric(anti.trafo(sweep(DNAm, 2, RData$CalcPCHorvath1$center) %*% RData$CalcPCHorvath1$rotation %*% RData$CalcPCHorvath1$model + RData$CalcPCHorvath1$intercept))
  pheno$PCHorvath2018 <- as.numeric(anti.trafo(sweep(DNAm, 2, RData$CalcPCHorvath2$center) %*% RData$CalcPCHorvath2$rotation %*% RData$CalcPCHorvath2$model + RData$CalcPCHorvath2$intercept))
  pheno$PCHannum <- as.numeric(sweep(DNAm, 2, RData$CalcPCHannum$center) %*% RData$CalcPCHannum$rotation %*% RData$CalcPCHannum$model + RData$CalcPCHannum$intercept)
  pheno$PCPhenoAge <- as.numeric(sweep(DNAm, 2, RData$CalcPCPhenoAge$center) %*% RData$CalcPCPhenoAge$rotation %*% RData$CalcPCPhenoAge$model + RData$CalcPCPhenoAge$intercept)
  pheno$PCDNAmTL <- as.numeric(sweep(DNAm, 2, RData$CalcPCDNAmTL$center) %*% RData$CalcPCDNAmTL$rotation %*% RData$CalcPCDNAmTL$model + RData$CalcPCDNAmTL$intercept)
  DNAm <- cbind(sweep(DNAm, 2, RData$CalcPCGrimAge$center) %*% RData$CalcPCGrimAge$rotation, Female = pheno$Female, Age = pheno$Age)
  pheno$PCPACKYRS <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCPACKYRS.model)] %*% RData$CalcPCGrimAge$PCPACKYRS.model + RData$CalcPCGrimAge$PCPACKYRS.intercept)
  pheno$PCADM <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCADM.model)] %*% RData$CalcPCGrimAge$PCADM.model + RData$CalcPCGrimAge$PCADM.intercept)
  pheno$PCB2M <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCB2M.model)] %*% RData$CalcPCGrimAge$PCB2M.model + RData$CalcPCGrimAge$PCB2M.intercept)
  pheno$PCCystatinC <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCCystatinC.model)] %*% RData$CalcPCGrimAge$PCCystatinC.model + RData$CalcPCGrimAge$PCCystatinC.intercept)
  pheno$PCGDF15 <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCGDF15.model)] %*% RData$CalcPCGrimAge$PCGDF15.model + RData$CalcPCGrimAge$PCGDF15.intercept)
  pheno$PCLeptin <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCLeptin.model)] %*% RData$CalcPCGrimAge$PCLeptin.model + RData$CalcPCGrimAge$PCLeptin.intercept)
  pheno$PCPAI1 <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCPAI1.model)] %*% RData$CalcPCGrimAge$PCPAI1.model + RData$CalcPCGrimAge$PCPAI1.intercept)
  pheno$PCTIMP1 <- as.numeric(DNAm[, names(RData$CalcPCGrimAge$PCTIMP1.model)] %*% RData$CalcPCGrimAge$PCTIMP1.model + RData$CalcPCGrimAge$PCTIMP1.intercept)
  pheno$PCGrimAge1 <- as.numeric(as.matrix(subset(pheno, select = RData$CalcPCGrimAge$components)) %*% RData$CalcPCGrimAge$PCGrimAge.model + RData$CalcPCGrimAge$PCGrimAge.intercept)


  return(pheno)
}
