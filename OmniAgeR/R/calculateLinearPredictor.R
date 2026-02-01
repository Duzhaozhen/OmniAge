#' @title (Internal Helper) Calculate a linear predictor for an epigenetic clock
#'
#' @description
#' A generic helper function to calculate the weighted linear sum
#' (predictor) for any epigenetic clock, given a beta matrix and a
#' coefficient list.
#'
#' @param beta.m A numeric matrix of beta values. (CpGs as rows, Samples as cols).
#' @param coef.lv A list containing the clock coefficients.
#'   - `coef.lv[[1]]` MUST be the numeric intercept.
#'   - `coef.lv[[2]]` MUST be a named numeric vector of CpG weights.
#' @param clock.name A string, the name of the clock for print messages.
#'
#' @return A named numeric vector of predicted ages/values for each sample.
#'   Returns NA for all samples if `min.perc` threshold is not met.
#'
calculateLinearPredictor <- function(beta.m,
                                     coef.lv,
                                     clock.name = "Unknown Clock",
                                     verbose = TRUE) {

  # --- 1. Extract intercept and weights ---
  intercept <- coef.lv[[1]]
  all_weights <- coef.lv[[2]]

  # --- 2. Match required CpGs with those in the input data ---
  required_cpgs <- names(all_weights)

  map.idx <- match(required_cpgs, rownames(beta.m))
  # rep.idx = indices of CpGs that are represented in the data
  rep.idx <- which(!is.na(map.idx))

  n_required <- length(required_cpgs)
  n_represented <- length(rep.idx)

  # --- 3. Print status message ---
  if (verbose) {
    print(paste0("[",clock.name,
      "] Number of represented ", clock.name, " CpGs (max=",
      n_required, ")=", n_represented))
  }

  # --- 3. Subset data and weights ---
  # Get the row indices from beta.m that match the required CpGs
  beta_matrix_indices <- map.idx[rep.idx]

  # Get the weights for those CpGs
  weights_subset <- all_weights[rep.idx]

  # Subset beta matrix
  tmpB.m <- beta.m[beta_matrix_indices, , drop = FALSE]

  # --- 6. Calculate the linear predictor ---
  # t(tmpB.m) -> [samples x CpGs]
  # weights_subset -> [CpGs x 1]
  # Result -> [samples x 1]
  age.v <- as.vector(intercept + t(tmpB.m) %*% weights_subset)

  names(age.v) <- colnames(beta.m)
  return(age.v)
}


