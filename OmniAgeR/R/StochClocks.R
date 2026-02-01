#' @title Calculate Stochastic Epigenetic Clocks (Tong et al., 2024)
#'
#' @description
#' Calculates DNAm-Age estimates based on one of three stochastic epigenetic
#' clocks (StocH, StocZ, StocP). These are stochastic analogues of the
#' Horvath, Zhang, and Levine/PhenoAge clocks, respectively.
#'
#' They consist of the same CpGs as the original clocks but were trained on
#' artificial cohorts generated via a stochastic process of age-related DNAm
#' change accrual. Although built using sorted monocyte samples from the
#' MESA study, their predictive ability is largely independent of immune
#' cell type or whole blood, as the DNAm data is standardized prior to
#' clock application.
#'
#' The stochastic clocks may serve as a useful complement to the original clocks,
#' providing a method to assess if specific age-acceleration patterns
#' (or decelerations) could be driven by an underlying stochastic process.
#' This may offer insights into the biological mechanisms of phenotypes
#' linked to epigenetic age acceleration.
#'
#' @details
#' This function relies on the `glmStocALL` data object. It assumes that
#' the command `data("glmStocALL")` will load a list object named
#' `glmStocALL.lo`, which contains the three fitted `glmnet` models.
#'
#' @param data.m A numeric matrix of normalized DNAm beta values.
#'   **Rows must be CpGs and columns must be samples.** Rownames must be
#'   CpG identifiers (e.g., "cg00000029") and colnames must be
#'   sample identifiers.
#' @param ages.v (Optional) A numeric vector of chronological ages for the
#'   samples. Must be the same length as the number of columns in `data.m`.
#'   If provided, this will be used to calculate EAA and IAA. (Default: `NULL`)
#' @param refM.m (Optional) A reference matrix for cell-type deconvolution
#'   (e.g., from the `EpiDISH` package, such as `centEpiFib.m`).
#'   If provided (along with `ages.v`), this will be used to calculate
#'   Intrinsic Age Acceleration (IAA). (Default: `NULL`)
#' @param verbose A boolean (default: `TRUE`). If `TRUE`, the function will
#'   print status messages during calculation, such as the number of
#'   CpGs found for each clock.
#'
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item `mage`: A list of three numeric vectors (StocH, StocP, StocZ)
#'     containing the predicted DNAm ages for each sample.
#'   \item `eaa`: A list of three numeric vectors containing the Extrinsic
#'     Age Acceleration (EAA) residuals, calculated as the residuals from
#'     `lm(mage ~ ages.v)`. This list will be empty if `ages.v` is `NULL`.
#'   \item `iaa`: A list of three numeric vectors containing the Intrinsic
#'     Age Acceleration (IAA) residuals, calculated as the residuals from
#'     `lm(mage ~ ages.v + estimated_cell_fractions)`. This list will be
#'     empty if `ages.v` or `refM.m` is `NULL`.
#'   \item `estF`: A matrix of estimated cell-type proportions as
#'     returned by `EpiDISH::epidish`. This will be `NULL` if `refM.m`
#'     is `NULL`.
#' }
#'
#' @import glmnet
#' @import EpiDISH
#' @export
#'
#' @references
#' Tong, H., Dwaraka, V.B., Chen, Q. et al.
#' Quantifying the stochastic component of epigenetic aging.
#' \emph{Nat Aging} (2024). \doi{10.1038/s43587-024-00636-6}
#'
#' @examples
#' # Load example data
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' # hannum_bmiq_m is a matrix: CpGs (rows) x Samples (columns)
#' # Hannum_example_pheno is a dataframe with an 'Age' column
#'
#' # Calculate clocks with age acceleration
#' StochClocks.out <- StochClocks(hannum_bmiq_m,
#'                                ages.v = PhenoTypesHannum_lv$Age)
#'
#'

StochClocks <- function(data.m,ages.v=NULL,refM.m=NULL, verbose = TRUE){

  # --- Start message ---
  if (verbose) {
    print(paste0("[StochClocks] Starting StochClocks calculation..."))
  }


  data("glmStocALL"); ## load in stochastic clock information
  Stoc_name <- paste0("Stoc",names(glmStocALL.lo) )
  estF.m <- NULL;
  if(!is.null(refM.m)){
    estF.m <- epidish(data.m,ref=refM.m,method="RPC",maxit=500)$est;
  }
  predMage.lv <- list();
  eaa.lv <- list();
  iaa.lv <- list();
  for(c in 1:length(glmStocALL.lo)){
    glm.o <- glmStocALL.lo[[c]];
    coef.m <- coef(glm.o);
    intC <- coef.m[1,ncol(coef.m)]; ### intercept term
    coef.v <- coef.m[-1,ncol(coef.m)]; ### estimated regression coefficients

    total_required_cpgs <- length(coef.v)
    commonCpGs.v <- intersect(names(coef.v), rownames(data.m));
    n_represented <- length(commonCpGs.v)

    if (verbose) {print(paste0("[StochClocks] Number of represented ", Stoc_name[c],
                               " CpGs (max=",total_required_cpgs, ")=", n_represented))
    }

    rep.idx <- match(commonCpGs.v,names(coef.v));
    map.idx <- match(commonCpGs.v,rownames(data.m));
    predMage.lv[[Stoc_name[c]]] <- as.vector(intC + matrix(coef.v[rep.idx],nrow=1) %*% data.m[map.idx,]);
    if(!is.null(ages.v)){
      eaa.lv[[Stoc_name[c]]] <- lm(predMage.lv[[c]] ~ ages.v)$res;
      if(!is.null(refM.m)){
        iaa.lv[[Stoc_name[c]]] <- lm(predMage.lv[[c]] ~ ages.v + estF.m)$res;
      }
    }
  }
  return(list(mage=predMage.lv,eaa=eaa.lv,iaa=iaa.lv,estF=estF.m));
}
