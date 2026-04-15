#' @title Calculate the Knight gestational age
#'
#' @description A function to calculate the Knight gestational age
#'
#' @param betaM A matrix of beta values (CpGs in rows, samples in columns).
#' This matrix must be pre-normalized (e.g., via BMIQ) and imputed.
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#'
#' @return A named vector of predicted Gestational Ages (in weeks).
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
                     minCoverage = 0,
                     verbose = TRUE) {
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
