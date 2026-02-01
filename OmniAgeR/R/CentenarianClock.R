#' @title Calculate Centenarian Epigenetic Clocks (Eric Dec et al.)
#'
#' @description  Calculates the Centenarian epigenetic clocks (ENCen40 and ENCen100) developed by Eric Dec et al. (2023). This function serves as a wrapper that loads the internal clock coefficients and computes the linear predictors for each clock using the helper function.
#' `calculateLinearPredictor`.
#'
#' @param beta.m a matrix of methylation beta values. Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#'
#' @param verbose A boolean (default: `TRUE`). If `TRUE`, the function will
#'   print status messages, including the number of represented CpGs
#'   for each clock.
#'
#'
#' @return A list containing the predicted scores for each Centenarian clock.
#' \itemize{
#'   \item `ENCen40`: Elastic net clock trained on individuals aged 40+.
#'   \item `ENCen100`: Elastic net clock specifically trained on centenarians (100+).
#' }
#'
#'
#'
#'
#'
#' @export
#'
#' @references
#' Dec, E., Clement, J., Cheng, K. et al.
#' Centenarian clocks: epigenetic clocks for validating claims of exceptional longevity.
#' \emph{GeroScience} 2023
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' CentenarianClock.out <- CentenarianClock(hannum_bmiq_m)


CentenarianClock <- function(beta.m,verbose = TRUE) {
  # --- Start message ---
  if (verbose) {
    print(paste0("[CentenarianClock] Starting CentenarianClock calculation..."))
  }

  res_list <- list()
  # --- Step 1: Load and parse coefficients ---
  data("CentenarianENCoef")
  for (i in seq_along(CentenarianENCoef)) {
    tmp_Coef  <- CentenarianENCoef[[i]]

    Coef_lv <- list()

    # Intercept
    Coef_lv[[1]] <- as.numeric(tmp_Coef[1, 2])

    # Coefficients
    coefficients <- as.numeric(as.vector(tmp_Coef[2:nrow(tmp_Coef), 2]))
    names(coefficients) <- as.vector(tmp_Coef[2:nrow(tmp_Coef), 1])
    Coef_lv[[2]] <- coefficients

    # --- Step 2: Calculate the linear predictor ---
    # (Requires the 'calculateLinearPredictor' function)
    predage.v <- calculateLinearPredictor(beta.m,
                                          coef.lv = Coef_lv,
                                          clock.name = names(CentenarianENCoef)[i],
                                          verbose)

    res_list[[names(CentenarianENCoef)[i]]] <- predage.v
  }


  # --- Step 3: Return final age vector ---
  # (Names are already attached by the helper function)
  return(res_list)
}

