#' @title Calculate the Mayne Placental Gestational Age.
#'
#' @description
#' Implements the 62-CpG placental epigenetic clock for estimating gestational
#' age (GA), as described by Mayne et al. (2017).
#'
#' @param betaM A numeric DNA methylation beta-value matrix with CpG probe
#'   identifiers as row names and samples as columns. CpG identifiers in
#'   \code{rownames(betaM)} are required. Sample identifiers in
#'   \code{colnames(betaM)} are optional. The matched CpG values used for
#'   calculation must not contain missing values.
#'   
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return
#' A **numeric vector** containing the predicted gestational age (in weeks) for
#'  each sample.
#'
#' @export
#'
#' @references
#' Mayne BT, Leemaqz SY, Smith AK, Breen J, Roberts CT, Bianco-Miotto T.
#' Accelerated placental aging in early onset preeclampsia pregnancies
#' identified by DNA methylation.
#' \emph{Epigenomics} 2017
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' mayneGaOut <- mayneGa(hannumBmiqM)
mayneGa <- function(betaM,
                    minCoverage = 0.5,
                    verbose = TRUE) {
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
    mayneGACoef <- loadOmniAgeRdata(
        "omniager_mayne_ga_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, mayneGACoef, "mayneGa",
        minCoverage, verbose
    )
    return(predAgev)
}
