#' @title Predicts DNA methylation age using the Zhang clock(Elastic Net model)
#' #'
#' @description
#' This function takes a matrix of DNA methylation beta values and calculates
#' the epigenetic age for each sample based on the elastic net model developed
#' by Zhang et al. (2019).
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @details
#' This function implements the elastic net epigenetic clock from
#' Zhang et al. (2019). A unique feature of this clock is the
#' pre-processing step. Instead of using raw beta values, it
#' first performs a per-sample standardization. Specifically, it calculates a
#' Z-score for each CpG based on the mean and standard deviation of all
#' measured CpGs within that same sample. The final age is then predicted
#' by taking the weighted sum of these standardized values using the model's
#' pre-defined coefficients.
#'
#' @return A numeric vector of predicted ages, with sample names preserved.
#'
#'
#' @references
#' Zhang Q, Vallerga CL, Walker RM, et al.
#' Improved precision of epigenetic clock estimates across tissues and
#' its implication for biological ageing.
#' \emph{Genome Med.} 2019
#'
#' @export
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' zhangClockOut <- zhangClock(hannumBmiqM)
#'
zhangClock <- function(betaM,
                       minCoverage = 0,
                       verbose = TRUE) {
    zhangClockCoef <- loadOmniAgeRdata(
        "omniager_zhang_clock_coef",
        verbose = verbose
    )

    # --- 1. Per-Sample Standardization ---
    if (verbose) message("[zhangClock] Performing per-sample standardization...")

    sample_means <- colMeans(betaM, na.rm = TRUE)
    sample_sds <- apply(betaM, 2, sd, na.rm = TRUE)

    # Safe handling: Prevent division by zero errors caused by a standard deviation of 0
    sample_sds[sample_sds == 0] <- 1

    # Z-score
    beta_scaled <- sweep(betaM, 2, sample_means, "-")
    beta_scaled <- sweep(beta_scaled, 2, sample_sds, "/")

    # --- 2. Predict Age ---
    predicted_age <- .calLinearClock(
        betaM       = beta_scaled,
        coefData    = zhangClockCoef,
        clockLabel  = "zhangClock",
        minCoverage = minCoverage,
        verbose     = verbose
    )

    return(predicted_age)
}
