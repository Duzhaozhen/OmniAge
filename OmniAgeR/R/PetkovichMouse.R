#' @title Petkovich Blood Mouse DNA Methylation Age Predictor
#'
#' @description
#' Predicts the biological age of mouse based on the blood-based epigenetic clock 
#' developed by Petkovich et al. (2017). This model utilizes the methylation 
#' status of 90 specific CpG sites to estimate age.
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
#' This predictor was constructed using **Elastic Net regression** on partial 
#' DNA methylomes obtained via Reduced Representation Bisulfite Sequencing (RRBS) 
#' from a cohort of 141 C57BL/6 mice.
#' 
#' The estimation process involves:
#' 1. **Linear Summation**: Calculating a raw DNAm score (MSc) by taking the 
#' weighted average of methylation levels across 90 CpG sites.
#' 2. **Non-linear Transformation**: Applying a **Power Law** conversion (i.e.
#' ((MSc - c) / a)^(1/b) ) to account for the decelerating rate of 
#' epigenetic changes observed during  mouse maturation.
#' 3. **Output**: Inverting the power law function to derive the predicted 
#' "Methylation Age".
#' 
#' @return A numeric vector of predictions.
#'
#' @export
#'
#' @references
#' Petkovich, Daniel A et al.
#' Using DNA Methylation Profiling to Evaluate Biological Age and Longevity Interventions.
#' \emph{Cell Metab.} 2017
#'
#' @examples
#' exampleMat <- loadOmniAgeRdata(
#'     "omniager_mouse_rrbs_example",
#'     verbose = FALSE
#' )
#' predRes <- petkovichMouse(exampleMat)
#'


petkovichMouse <- function(betaM,
                           minCoverage = 0,
                           verbose = TRUE) {
  # 1. Obtain the weight coefficients
  petkovichCoef <- loadOmniAgeRdata(
    "omniager_petkovich_mouse_coef",
    verbose = verbose
  )
  
  # 2. Calculate the linear DNAm Score (MSc)
  msc <- .calLinearClock(
    betaM, petkovichCoef, "petkovichMouse",
    minCoverage, verbose
  )
  
  # 3. Perform age conversion
  if (verbose) {
      message("[petkovichMouse] Applying power-law transformation to calculate age in months...")
    }
  a <- 0.1666
  b <- 0.4185
  c <- -1.712
  
  # Apply non-linear power-law conversion to calculate age: ((MSc - c) / a)^(1/b)
  predAgev <- ((msc - c) / a)^(1/b)
  predAgev <- predAgev/30.5  # days to months

  
  names(predAgev) <- colnames(betaM)
  
  return(predAgev)
}
