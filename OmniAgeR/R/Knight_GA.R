#' @title Calculate the Knight gestational age
#'
#' @description A function to calculate the Knight gestational age
#'
#'
#' @param betaM A numeric DNA methylation beta-value matrix with CpG probe
#'   identifiers as row names and samples as columns. CpG identifiers in
#'   \code{rownames(betaM)} are required. Sample identifiers in
#'   \code{colnames(betaM)} are optional. The matched CpG values used for
#'   calculation must not contain missing values.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.5.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A vector of predicted Gestational Ages (in weeks).
#'
#' @export
#'
#' @references
#' Knight AK, Craig JM, Theda C, et al.
#' An epigenetic clock for gestational age at birth
#' based on blood methylation data.
#' \emph{Genome Biol.} 2016
#'
#' @examples
#' hannumBmiqM <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )[[1]]
#' knightGaOut <- knightGa(hannumBmiqM)
#'
knightGa <- function(betaM,
                     minCoverage = 0.5,
                     verbose = TRUE) {
   
    betaM <- .validateBetaMatrix(
      betaM,
      requireColnames = FALSE
    )
  
    knightCoef <- loadOmniAgeRdata(
        "omniager_knight_ga_coef",
        verbose = verbose
    )

    predAgev <- .calLinearClock(
        betaM, knightCoef, "knightGa",
        minCoverage, verbose
    )
    return(predAgev)
}
