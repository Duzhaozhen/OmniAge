#' @title Meer Whole Lifespan Multi-Tissue Mouse DNA Methylation Age Predictor
#'
#' @description
#' Predicts the chronological age of mouse samples based on the Whole Lifespan 
#' Multi-tissue (WLMT) DNA methylation clock developed by Meer et al. (2018).
#' 
#' @param betaM A numeric matrix of beta values (ranging from 0 to 1). Rows must 
#' be labeled with mouse genomic coordinates (typically in 'chr:pos' format 
#' based on the mm10/GRCm38 assembly) corresponding to the clock's CpG sites.
#' Columns represent individual samples.
#' 
#' @param minCoverage A numeric value (0-1). The minimum proportion of
#'   required CpGs that must be present. Default is 0.
#' 
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' 
#' @details
#' The clock is a robust age predictor that covers the entire mouse lifespan 
#' and remains unbiased across different age groups and sexes. It utilizes an 
#' elastic net regression model across 435 specific CpG sites. The raw output 
#' is calculated in days and is converted to months (divided by 30.5) by this function.
#' 
#' @return A numeric vector of predicted chronological ages in months
#'
#' @export
#'
#' @references
#' Meer, Margarita V et al.
#' A whole lifespan mouse multi-tissue DNA methylation clock.
#' \emph{Elife} 2018
#' 
#' @examples
#' exampleMat <- loadOmniAgeRdata(
#'     "omniager_mouse_rrbs_example",
#'     verbose = FALSE
#' )
#' predRes <- meerMouse(exampleMat)
#' 

meerMouse <- function(betaM,
                        minCoverage = 0,
                        verbose = TRUE) {
  
  if (!is.numeric(betaM)) {
    stop("Error: 'betaM' must be a numeric matrix or data frame.")
  }
  
  if (any(betaM < 0, na.rm = TRUE) || any(betaM > 1, na.rm = TRUE)) {
    stop("Error: 'betaM' values must be between 0 and 1 (beta values). ",
         "Please check if your input is already in percentage or contains invalid values.")
  }
  
  # 1. Obtain the model data list
  meerMouseCoef <- loadOmniAgeRdata(
    "omniager_meer_mouse_coef",
    verbose = verbose
  )
  
 
  # 2. Calculate the linear DNAm Score using the betaM matrix
  rawAge <- .calLinearClock(
    betaM = betaM * 100, # Convert to percentage
    coefData = meerMouseCoef, 
    clockLabel = "meerMouse",
    minCoverage = minCoverage, 
    verbose = verbose
  )
  #  Transforms age in days to age in months.
  predAgeMonths <- rawAge / 30.5
  
  return(predAgeMonths)
}