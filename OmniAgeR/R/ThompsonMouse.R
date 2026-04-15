#' @title Thompson Mouse Multi-Tissue DNA Methylation Age Predictor
#'
#' @description
#' Predicts the chronological age of mice across multiple tissues and their entire lifespan 
#' using the epigenetic clock developed by Thompson et al. (2018).
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
#' This function applies the epigenetic clock models from Thompson et al. (2018), 
#' which were trained primarily on reduced representation bisulfite sequencing (RRBS) data 
#' mapped to the mm10 mouse genome. The clock is highly robust and applicable across 
#' multiple mouse tissues, including adipose, blood, cerebellum, cortex, heart, kidney, 
#' liver, lung, muscle, and spleen.
#' 
#' @return A numeric vector of predicted chronological ages in months
#'
#' @export
#'
#' @references
#' Thompson, Michael J et al. 
#' A multi-tissue full lifespan epigenetic clock for mice.
#' \emph{Aging} 2018
#'
#' @examples
#' exampleMat <- loadOmniAgeRdata(
#'     "omniager_mouse_rrbs_example",
#'     verbose = FALSE
#' )
#' predRes <- thompsonMouse(exampleMat)

thompsonMouse <- function(betaM,
                      minCoverage = 0,
                      verbose = TRUE) {
  
  # 1. Obtain the model data list
  clockCoef <- loadOmniAgeRdata(
    "omniager_thompson_mouse_coef",
    verbose = verbose
  )
  
  
  # 2. Calculate the linear DNAm Score using the betaM matrix
  predAgeMonths <- .calLinearClock(
    betaM = betaM, # Convert to percentage
    coefData = clockCoef, 
    clockLabel = "thompsonMouse",
    minCoverage = minCoverage, 
    verbose = verbose
  )
 
  
  return(predAgeMonths)
}