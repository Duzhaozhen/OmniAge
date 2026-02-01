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
#' The coefficients for these 109 predictors are loaded from the
#' internal `EpiScoresCoef.rda` data object.
#'
#' **Imputation Handling:**
#' This function features a robust **two-stage imputation process** for
#' handling missing CpG data:
#' 1.  **Sample-level NA Imputation:** Any existing `NA` values within the
#'     provided `beta.m` matrix are imputed to the **row-wise mean**
#'     (the mean of that specific CpG across all samples in the input data).
#' 2.  **Missing CpG Imputation:** CpGs required for the scores but
#'     *entirely absent* from `beta.m` are added. Their values are imputed
#'     using the **mean beta value from the original training cohort**
#'     (stored in `EpiScoresCoef`).
#'
#'
#' After imputation, the function iterates through each unique protein score,
#' subsetting the appropriate coefficients, and computes a linear predictor
#' using the `calculateLinearPredictor` helper function.
#'
#' @param beta.m A numeric matrix of normalized DNAm beta values. Rows must correspond to CpG identifiers and columns to samples. `NA` values are permitted and will be imputed internally
#'
#' @param verbose A logical value (default: `TRUE`) indicating whether to print progress messages and warnings (e.g., about missing CpGs) to the console.
#'
#' @return A list of length 109. Each element of the list corresponds to one calculated EpiScore (one for each protein). The name of each list element is the name of the EpiScore. Each element contains a named numeric vector of the calculated scores for all samples.
#'
#' @references
#' Gadd DA, Hillary RF, McCartney DL, et al.
#' Epigenetic scores for the circulating proteome as tools for disease prediction.
#' \emph{Elife.} 2022
#'
#'
#' @export
#'
#' @examples
#'
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#'
#' all_episcores <- CompEpiScores(hannum_bmiq_m, verbose = TRUE)
#'
#'


CompEpiScores <- function(beta.m,verbose=TRUE) {

  if (verbose) {
    print(paste0("[EpiScores] Starting EpiScores calculation..."))
  }
  # --- Step 1: Load and parse coefficients ---
  data("EpiScoresCoef")
  res_list <- list()
  episcores_predictor <- unique(EpiScoresCoef[[4]])

  beta.m <-  beta.m[intersect(rownames(beta.m), EpiScoresCoef[[1]]),]

  ## Convert NAs to Mean Value for all individuals across each probe
  na_to_mean <-function(methyl){
    methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
    return(methyl)
  }

  beta.m <- t(apply(beta.m,1,function(x) na_to_mean(x)))

  ## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site

  beta.m <- if(nrow(beta.m) == length(unique(EpiScoresCoef[[1]]))) {
    message("[EpiScores] All sites present");
    beta.m
  } else if(nrow(beta.m)==0){
    message("[EpiScores] There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
  } else {
    missing_cpgs = EpiScoresCoef[-which(EpiScoresCoef[[1]] %in% rownames(beta.m)),c("CpG_Site","Mean_Beta_Value")]
    message(paste("[EpiScores]",length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
    mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(beta.m))
    row.names(mat) <- unique(missing_cpgs$CpG_Site)
    colnames(mat) <- colnames(beta.m)
    mat[is.na(mat)] <- 1
    missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) {
      missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
    } else {missing_cpgs
    }
    ids <- unique(row.names(mat))
    missing_cpgs1 <- missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
    mat<-mat*missing_cpgs1$Mean_Beta_Value
    beta.m <- rbind(beta.m,mat)
  }
  # --- Step 2: Calculate the linear predictor ---
  # (Requires the 'calculateLinearPredictor' function)

  for (i in seq_along(episcores_predictor)) {
    tmp_Coef  <- EpiScoresCoef[EpiScoresCoef[[4]] %in% episcores_predictor[i], ]

    Coef_lv <- list()

    # Intercept
    Coef_lv[[1]] <- 0

    # Coefficients
    coefficients <- as.numeric(as.vector(tmp_Coef[,2]))
    names(coefficients) <- as.vector(tmp_Coef[, 1])
    Coef_lv[[2]] <- coefficients



    predage.v <- calculateLinearPredictor(beta.m,
                                          coef.lv = Coef_lv,
                                          clock.name = episcores_predictor[i],
                                          verbose)

    res_list[[episcores_predictor[i]]] <- predage.v
  }

  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(res_list)
}
