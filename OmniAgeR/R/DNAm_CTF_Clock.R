#' DNA Methylation Cell-Type Fraction (CTF) Aging Clock
#'
#' @description
#' Predicts biological age based on immune cell type fractions derived from
#' DNA methylation data.
#'
#' @param ctfM A numeric matrix or data frame where rows are samples and columns are cell types.
#'   Must contain the specific cell types required by the model (e.g., predicted by EpiDISH).
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named numeric vector of predicted ages.
#' @export
#' @importFrom randomForest randomForest importance
#'
#'
#' @examples
#'
#' tzhFracM <- loadOmniAgeRdata(
#'     "omniager_tzh_example_ctf",
#'     verbose = FALSE
#' )[[2]]
#' dnamCTFClockOut <- dnamCTFClock(ctfM = tzhFracM)
#'
dnamCTFClock <- function(ctfM, verbose = TRUE) {
    # --- 1. Load the internal model ---
    dnamCtfModel <- loadOmniAgeRdata(
        "omniager_dnam_ctf_model",
        verbose = verbose
    )
    # --- 2. Verify feature integrity ---
    requiredFeatures <- rownames(dnamCtfModel$importance)

    missingCols <- setdiff(requiredFeatures, colnames(ctfM))
    if (length(missingCols) > 0) {
        stop(
            "[dnamCTFClock] Missing required cell types: ",
            paste(missingCols, collapse = ", "),
            ". Ensure you provide estimated fractions for all required types."
        )
    }

    # --- 3. Data alignment and NA inspection ---
    dataForPred <- ctfM[, requiredFeatures, drop = FALSE]

    if (any(is.na(dataForPred))) {
        stop(
            "[dnamCTFClock] Input contains NA values. ",
            "Random Forest model requires complete data."
        )
    }

    # --- 5. Prediction ---
    predAge <- stats::predict(dnamCtfModel, newdata = dataForPred)

    if (!is.null(rownames(ctfM))) {
        names(predAge) <- rownames(ctfM)
    }

    return(predAge)
}
