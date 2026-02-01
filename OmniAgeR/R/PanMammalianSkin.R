#' @title Calculate Universal Pan-Mammalian Skin Epigenetic Clocks
#' @description
#' Applies the two universal pan-mammalian BLOOD clocks (Clock 2 and 3)
#' from Lu et al. (2023) to a given dataset of DNA methylation values.
#'
#' @details
#' This function is a specific adaptation for the Skin-only clocks.
#' It requires the 'PanMammalianSkinCoef' data object to be available
#' in the R package's data directory.
#'
#' @param sample_info A data.frame with sample metadata. Must include
#'   'Sample' (unique sample ID) and 'SpeciesLatinName' (Latin species name).
#'   Other columns like 'Age' or 'Tissue' are optional but will be included in the output if present.
#' @param beta.m A matrix or data.frame of methylation beta values in the format
#'   CpGs x Samples. Rownames must be CpG identifiers, and column names
#'   must be sample identifiers ('Sample') matching those in `sample_info`.
#' @param anage_data A data.frame containing the AnAge database information.
#'   Must include 'SpeciesLatinName', 'GestationTimeInYears',
#'   'averagedMaturity.yrs', and 'maxAge'. # MODIFIED (Typo fix)
#' @param verbose A logical value (default: FALSE) indicating whether to print
#'   detailed progress messages to the console. # NEW
#'
#' @return A data.frame containing the 'Sample', 'SpeciesLatinName',
#'   and the calculated ages: 'DNAmAgePanMammalianSkin2',
#'   'DNAmAgePanMammalianSkin3', 'DNAmRelativeAge_Skin', and 'DNAmRelativeAdultAge_Skin'.
#'
#' @export
#'
#' @references
#' Lu, A.T., Fei, Z., Haghani, A. et al.
#' Universal DNA methylation age across mammalian tissues.
#' \emph{Nat Aging.} 2023
#'
#'
#'
#' @examples
#'
#' download_OmniAgeR_example("Tursiops_example")
#' load_OmniAgeR_example("Tursiops_example")
#' data("anage_data")
#'
#' # Run the calculation with progress messages
#'
#' clock_results <- PanMammalianSkin(
#'   sample_info = example_sample_info,
#'   beta.m = example_beta_m,
#'   anage_data = anage_data,
#'   verbose = TRUE # NEW
#' )
#'



PanMammalianSkin <- function(sample_info,
                             beta.m,
                             anage_data,
                             verbose = FALSE) {

  if (verbose) {
    print(paste0("[PanMammalianSkin] Starting Universal Pan-Mammalian Skin Clock calculation..."))
  }

  data("PanMammalianSkinCoef")

  # --- 1. Define Internal Helper Functions ---
  # (Helper functions are unchanged)
  F2_antitrans_clock2 <- function(y, y.maxAge, y.gestation, const = 1) {
    x0 = const * exp(-exp(-1 * y))
    x1 = x0 * (y.maxAge + y.gestation)
    x = x1 - y.gestation
    return(x)
  }
  F1_logli <- function(age1, m1, m2 = m1, c1 = 1) {
    ifelse(age1 >= m1, (age1 - m1) / m2, c1 * log((age1 - m1) / m2 / c1 + 1))
  }
  F2_revtrsf_clock3 <- function(y.pred, m1, m2 = m1, c1 = 1) {
    ifelse(y.pred < 0, (exp(y.pred / c1) - 1) * m2 * c1 + m1, y.pred * m2 + m1)
  }
  F3_loglifn = function(dat1, b1 = 1, max_tage = 4,
                        c1 = 5, c2 = 0.38, c0 = 0) {
    a2 = (dat1$GestationTimeInYears + c0) / (dat1$averagedMaturity.yrs)
    a_Logli = c1 * a2^c2
    dat1$a_Logli = a_Logli
    if ("Age" %in% names(dat1)) {
      age1 = (dat1$maxAge + dat1$GestationTimeInYears) / (dat1$averagedMaturity.yrs + dat1$GestationTimeInYears)
      a1 = age1 / (1 + max_tage)
      dat1$a1_Logli = a1
      x = dat1$Age + dat1$GestationTimeInYears
      t2 = dat1$averagedMaturity.yrs * b1 + dat1$GestationTimeInYears
      x2 = x / t2
      y = F1_logli(x2, a_Logli, a_Logli)
      dat1$LogliAge <- y
    }
    return(dat1)
  }

  # --- 2. Define Internal Variable Names ---

  y.name <- c('Y.pred2', 'Y.pred3')
  age.name <- c('DNAmAgePanMammalianSkin2', 'DNAmAgePanMammalianSkin3')

  # --- 3. Process Input Data ---

  if (!all(c("SpeciesLatinName", "Sample") %in% names(sample_info))) {
    stop("[PanMammalianSkin] `sample_info` must contain 'SpeciesLatinName' and 'Sample' columns.")
  }
  if (is.null(rownames(beta.m)) || is.null(colnames(beta.m))) {
    stop("[PanMammalianSkin] `beta.m` must have rownames (CpG IDs) and colnames (Sample IDs).")
  }
  samples_info <- as.character(sample_info$Sample)
  samples_beta <- colnames(beta.m)
  if (!all(samples_info %in% samples_beta)) {
    missing_in_beta <- samples_info[!samples_info %in% samples_beta]
    stop(paste("[PanMammalianSkin] The following samples are in `sample_info` but not in `beta.m` colnames:",
               paste(missing_in_beta, collapse=", ")))
  }

  # Assign inputs to internal variable names
  info <- sample_info
  anage <- anage_data
  dat.meth0 <- beta.m
  glmnet.list <- PanMammalianSkinCoef # Uses the object from the package's data

  # (Step 1) Merge AnAge data into sample info
  info <- merge(by = 'SpeciesLatinName', info,
                subset(anage, select = c(SpeciesLatinName, GestationTimeInYears,
                                         averagedMaturity.yrs, maxAge)))

  # (Step 2) Generate HighmaxAge variable (for Clock 2)
  MYMAX <- 1.3
  info$HighmaxAge <- MYMAX * info$maxAge
  info$HighmaxAge[info$SpeciesLatinName == 'Homo sapiens'] <- info$maxAge[info$SpeciesLatinName == 'Homo sapiens']
  info$HighmaxAge[info$SpeciesLatinName == 'Mus musculus'] <- info$maxAge[info$SpeciesLatinName == 'Mus musculus']

  # (Step 3) Prepare methylation data

  mycpgs <- unique(c(glmnet.list[[1]]$probe,
                     glmnet.list[[2]]$probe))
  mycpgs <- mycpgs[mycpgs != 'Intercept']

  # Check for missing CpGs
  missing_cpgs <- setdiff(mycpgs, rownames(dat.meth0))

  if (verbose) {
    print(paste0(
      "[PanMammalianSkin] Number of represented PanMammalianSkin CpGs (max=",
      length(mycpgs), ")=", length(mycpgs)-length(missing_cpgs)))
  }



  if (length(missing_cpgs) > 0) {
    # This warning will fire regardless of 'verbose' setting, which is good.
    # warning(paste("[PanMammalianSkin] The following required CpGs are missing from `beta.m` rownames and will be treated as having NA values:",
    #               paste(missing_cpgs, collapse=", ")))
    missing_mat <- matrix(NA_real_,
                          nrow = length(missing_cpgs),
                          ncol = ncol(dat.meth0),
                          dimnames = list(missing_cpgs, colnames(dat.meth0)))
    dat.meth0 <- rbind(dat.meth0, missing_mat)
  }

  dat.meth0_filtered <- dat.meth0[mycpgs, , drop = FALSE]
  dat.meth <- as.data.frame(t(dat.meth0_filtered))
  dat.meth$Sample <- rownames(dat.meth)
  dat.meth$Intercept <- 1

  # (Step 4) Merge methylation data with info
  info <- merge(by = 'Sample', info, dat.meth, all.x = TRUE)


  # --- 4. Run Predictions ---
  for (k in seq_along(glmnet.list)) {
    glmnet_model <- glmnet.list[[k]]
    model_vars <- as.character(glmnet_model$probe)
    model_coefs <- glmnet_model$coef

    info_subset_matrix <- as.matrix(info[, model_vars])
    info_subset_matrix[is.na(info_subset_matrix)] <- 0

    info[, y.name[k]] <- as.numeric(info_subset_matrix %*% model_coefs)
  }

  # --- 5. Inverse Transform Predictions to Age ---
  # (1) Clock 1
  #info[, age.name[1]] <- exp(info[, y.name[1]]) - 2

  # (2) Clock 2
  info$DNAmRelativeAge <- exp(-exp(-1 * info[, y.name[1]]))
  info[, age.name[1]] <- F2_antitrans_clock2(info[, y.name[1]],
                                             info$HighmaxAge,
                                             info$GestationTimeInYears,
                                             const = 1)

  # (3) Clock 3
  info <- F3_loglifn(info)
  info$m1 <- info$a_Logli

  info$DNAmRelativeAdultAge <- F2_revtrsf_clock3(info[, y.name[2]], info$m1)

  info[, age.name[2]] <-
    info$DNAmRelativeAdultAge * (info$averagedMaturity.yrs + info$GestationTimeInYears) - info$GestationTimeInYears

  # --- 6. Format and Return Output ---

  output_cols <- c('Sample', 'SpeciesLatinName', 'MammalNumberHorvath', 'Age', 'Tissue',
                   'DNAmRelativeAge_Skin', 'DNAmRelativeAdultAge_Skin', age.name)

  available_cols <- intersect(output_cols, names(info))

  output <- info[, available_cols]

  return(output)
}
