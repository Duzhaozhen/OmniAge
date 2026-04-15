#' @title Pipek's Multi-tissue Elastic Net Epigenetic Clock (239 CpGs)
#'
#' @description
#' Implements the "elasticNet (239)" epigenetic clock model proposed by
#' Pipek and Csabai (2023). This model was trained using elastic net
#' regression on a large multi-tissue dataset containing methylation
#' data from Illumina 27K, 450K, and EPIC platforms. It is designed to
#' provide improved accuracy on EPIC array data compared to the original
#' Horvath2013 clock.
#'
#' @details
#' Implements the "elasticNet (239)" model using **239 CpGs** shared
#' across Illumina 27K, 450K, and EPIC arrays.
#'
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed prior to input.
#'   \item **Normalization:** Recommended but not strictly required.
#' }
#'
#' Internally applies Horvath's (2013) log-linear age transformation.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted biological ages. The vector is
#' named using the sample IDs from the \code{rownames} of \code{betaM}.
#'
#' @export
#'
#' @references
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for
#' methylation array data. \emph{J Math Chem} 2023
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' pipekElasticNetOut <- pipekElasticNet(hannumBmiqM)
#'
pipekElasticNet <- function(betaM,
                            minCoverage = 0,
                            verbose = TRUE) {
    pipekElasticNetCoef <- loadOmniAgeRdata(
        "omniager_pipek_elasticnet_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, pipekElasticNetCoef, "pipekElasticNet",
        minCoverage, verbose
    )
    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}


#' @title Pipek's Filtered Horvath Epigenetic Clock (272 CpGs)
#'
#' @description
#' Implements the "filtered H (272)" epigenetic clock model proposed by
#' Pipek and Csabai (2023). This model restricts feature selection to the
#' subset of CpG sites from the original Horvath pan-tissue clock that
#' are also present on the EPIC array.
#'
#' @details
#' Implements the "filtered H (272)" model. Unlike the full ElasticNet
#' model, this clock was trained by limiting candidate features to the
#' intersection of original Horvath probes and the cross-platform
#' (27K/450K/EPIC) probe set.
#'
#' **Key Features:**
#' \itemize{
#'   \item **272 CpGs:** A subset of the original Horvath clock,
#'   re-optimized for better EPIC array compatibility.
#'   \item **Best for:** Datasets pre-filtered to Horvath clock probes
#'   but requiring updated calibration.
#' }
#'
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed.
#'   \item **Normalization:** Recommended.
#' }
#'
#' Internally applies Horvath's (2013) log-linear age transformation.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted biological ages. The vector is
#' named using the sample IDs from the \code{rownames} of \code{betaM}.
#'
#' @export
#'
#' @references
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for
#' methylation array data. \emph{J Math Chem} 2023
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' pipekFilteredhOut <- pipekFilteredh(hannumBmiqM)
#'
pipekFilteredh <- function(betaM,
                           minCoverage = 0,
                           verbose = TRUE) {
    pipekFilteredhCoef <- loadOmniAgeRdata(
        "omniager_pipek_filteredh_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, pipekFilteredhCoef, "pipekFilteredh",
        minCoverage, verbose
    )
    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}


#' @title Pipek's Retrained Horvath Epigenetic Clock (308 CpGs)
#'
#' @description
#' Implements the "retrained H (308)" epigenetic clock model proposed by
#' Pipek and Csabai (2023). This model is a direct recalibration of the
#' original Horvath clock probes using a large multi-platform training set.
#'
#' @details
#' Implements the "retrained H (308)" model. It includes all **308 CpG sites** #' from the original Horvath pan-tissue clock that are present across 27K,
#' 450K, and EPIC platforms.
#'
#' **Key Features:**
#' \itemize{
#'   \item **No Feature Selection:** Coefficients were simply refitted to the
#'   new data without dropping probes.
#'   \item **Robust Update:** Serves as a direct update to the Horvath clock
#'   to correct for accuracy loss on EPIC arrays.
#' }
#'
#' **Input Requirements:**
#' \itemize{
#'   \item **Complete Data:** Missing values must be imputed.
#'   \item **Normalization:** Recommended.
#' }
#'
#' Internally applies Horvath's (2013) log-linear age transformation.
#'
#' @param betaM A numeric matrix of DNA methylation beta values.
#'   `rownames` (CpG probe IDs) and `colnames` (Sample IDs) are required.
#'   The matrix should not contain `NA` values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A numeric vector of predicted biological ages. The vector is
#' named using the sample IDs from the \code{rownames} of \code{betaM}.
#'
#' @export
#'
#' @references
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for
#' methylation array data. \emph{J Math Chem} 2023
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' pipekRetrainedhOut <- pipekRetrainedh(hannumBmiqM)
#'
pipekRetrainedh <- function(betaM,
                            minCoverage = 0,
                            verbose = TRUE) {
    pipekRetrainedhCoef <- loadOmniAgeRdata(
        "omniager_pipek_retrainedh_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, pipekRetrainedhCoef, "pipekRetrainedh",
        minCoverage, verbose
    )
    # transformation
    predAgev <- .antiTrafo(predAgev)
    return(predAgev)
}
