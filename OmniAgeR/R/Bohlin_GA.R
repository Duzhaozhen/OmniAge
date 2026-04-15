#' @title Calculate the Bohlin Gestational Age (Cord Blood)
#'
#' @description
#' Implements the epigenetic clock for predicting gestational age (GA) using
#' newborn cord blood, as described by Bohlin et al. (2016).
#'
#' @param betaM a matrix of methylation beta values.
#' Needs to be rows = samples and columns = CpGs, with rownames and colnames.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#'
#' @return #' A **numeric vector** containing the predicted gestational
#' age (in weeks) for each sample. The vector is named with the sample IDs
#' from the `rownames` of `betaM`.
#'
#' @export
#'
#' @references
#' Bohlin J, Håberg SE, Magnus P, et al.
#' Prediction of gestational age based on genome-wide differentially
#' methylated regions. \emph{Genome Biol.} 2016
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' bohlinGaOut <- bohlinGa(hannumBmiqM)
bohlinGa <- function(betaM,
                     minCoverage = 0,
                     verbose = TRUE) {
    bohlinGACoef <- loadOmniAgeRdata(
        "omniager_bohlin_ga_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, bohlinGACoef, "bohlinGa",
        minCoverage, verbose
    )
    ## Convert the number of days into weeks
    predAgev <- predAgev / 7
    return(predAgev)
}
