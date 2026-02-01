#' @title Calculate EnsembleAge Epigenetic Clocks
#'
#' @description
#' Calculates epigenetic age using the EnsembleAge multi-clock framework
#' (Haghani et al., 2025). This function computes one of the three available
#' clock versions (Static, Dynamic, or HumanMouse) by loading the
#' pre-calculated coefficients from the `EnsembleAgeCoef` package data object.
#'
#' @param beta.m A numeric matrix of methylation beta values.
#'   **Rows must correspond to CpG probes** (e.g., `cg00000165`) and
#'   **columns to samples** (e.g., `GSM990532`). Rownames (CpG IDs) and
#'   colnames (Sample IDs) are required.
#' @param clock_version A character string specifying which version of the
#'   EnsembleAge clocks to calculate. Valid options are `"Dynamic"`, `"Static"`,
#'   or `"HumanMouse"`. Default is `"HumanMouse"`.
#' @param verbose A logical flag. If `TRUE` (the default), the function will
#'   print status messages during calculation, including messages from the
#'   internal `calculateLinearPredictor` function.
#'
#' @return
#' A `list` where each element is a named numeric vector of predicted values.
#' The names of the list elements correspond to the specific sub-clocks
#' calculated (e.g., `HumanMouse_HumanMouse`, `Static_Static`, etc.).
#'
#' @note
#' **Understanding the `clock_version` output:**
#' * **`"Static"` and `"HumanMouse"`:** The returned values are the final,
#'     calibrated age predictors. The `"HumanMouse"` clock returns age
#'     normalized by maximum species lifespan.
#' * **`"Dynamic"`:** This option returns a list of predictors from all
#'     40 clocks in the 'Dynamic' ensemble library. To complete the full
#'     "EnsembleAge.Dynamic" *methodology* as described in the paper,
#'     the user must perform **additional steps**:
#'     1.  Have a dataset with control and intervention groups.
#'     2.  Use the returned 40 predictors to find "responsive" clocks
#'         (e.g., those with |Z-score| > 2 between groups).
#'     3.  Calculate the final `EnsembleAge.Dynamic` value by taking the
#'         **median** of only those responsive clocks.
#'     This function provides the necessary *basis* for this calculation.
#'
#' @references
#' Haghani, A., Lu, A.T., Yan, Q. et al.
#' EnsembleAge: enhancing epigenetic age assessment with a multi-clock framework.
#' \emph{GeroScience} (2025).
#'
#' @export
#' @importFrom utils data
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' # Ensure it's CpGs=rows, Samples=cols
#' # Calculate the HumanMouse clock version
#' age_predictions <- EnsembleAge(hannum_bmiq_m, clock_version = "HumanMouse")
#'
EnsembleAge <- function(beta.m, clock_version = "HumanMouse", verbose = TRUE) {

  # --- 1. Start message ---
  if (verbose) {
    print(paste0("[EnsembleAge] Starting EnsembleAge_", clock_version, " calculation..."))
  }

  res_list <- list()

  # --- Step 1: Load and parse coefficients ---
  # This lazy-loads the 'EnsembleAgeCoef' data object from the package
  data("EnsembleAgeCoef", envir = environment())

  # Check if clock_version is valid
  if (!clock_version %in% names(EnsembleAgeCoef)) {
    stop(paste("[EnsembleAge] Invalid 'clock_version'. Must be one of:", paste(names(EnsembleAgeCoef), collapse = ", ")))
  }

  temp_coef_list <- EnsembleAgeCoef[[clock_version]]

  # --- Step 2: Loop through each sub-clock ---
  for (i in seq_along(temp_coef_list)) {

    temp_Coef_sub <- temp_coef_list[[i]]
    clock_name <- names(temp_coef_list)[i]

    # Check for valid coefficients
    if (is.null(temp_Coef_sub) || nrow(temp_Coef_sub) < 2) {
      if(verbose) message(paste("[EnsembleAge] Warning: Skipping clock", clock_name, "due to missing/invalid coefficients."))
      next
    }

    Coef_lv <- list()

    # Intercept
    Coef_lv[[1]] <- as.numeric(temp_Coef_sub[1, 2])

    # Coefficients
    coefficients <- as.numeric(as.vector(temp_Coef_sub[2:nrow(temp_Coef_sub), 2]))
    names(coefficients) <- as.vector(temp_Coef_sub[2:nrow(temp_Coef_sub), 1])
    Coef_lv[[2]] <- coefficients

    # --- Step 3: Calculate the linear predictor ---
    # (Requires the 'calculateLinearPredictor' function)
    predage.v <- calculateLinearPredictor(
      beta.m,
      coef.lv = Coef_lv,
      clock.name = paste0(clock_version, "_", clock_name),
      verbose = verbose # Pass the verbose status to the helper function
    )

    res_list[[paste0(clock_version, "_", clock_name)]] <- predage.v
  }

  return(res_list)
}








