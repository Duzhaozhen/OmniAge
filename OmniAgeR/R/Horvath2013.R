#' @title Calculate Horvath's DNAm Age (2013)
#'
#' @description A function to calculate the Horvath epigenetic clock age(2013) from a DNA methylation beta value matrix.
#'
#' @param beta.m A numeric matrix of beta values. Rows should be CpG probes and columns should be individual samples.
#'
#' @details Implements the Horvath (2013) pan-tissue clock. The function calculates a weighted linear predictor from 353 CpGs found in the input matrix and then transforms this value using a non-linear function to return the final DNAm age.
#'
#' @return A numeric vector of predicted DNAm ages, with names corresponding to the sample IDs from the input matrix's column names.
#'
#' @seealso The main function \code{\link{EpiAge}} can be used to calculate multiple clocks simultaneously.
#'
#' @export
#'
#' @references
#' Horvath S.
#' DNA methylation age of human tissues and cell types.
#' \emph{Genome Biol.} 2013
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' Horvath2013.out <- Horvath2013(hannum_bmiq_m)
#'


Horvath2013 <- function(beta.m) {

  # Step 1: Load the coefficients and parse them into the required list format.
  # This creates an object called `horv.df` in the function's local environment.
  data('Horvath2013CpG')

  horvageP.lv <- list()
  # Intercept
  horvageP.lv[[1]] <- as.numeric(horv.df[1, 2])
  # Coefficients
  coefficients <- as.numeric(as.vector(horv.df[2:nrow(horv.df), 2]))
  names(coefficients) <- as.vector(horv.df[2:nrow(horv.df), 1])
  horvageP.lv[[2]] <- coefficients

  # Step 2: Define the necessary helper functions inside the main function.
  # This keeps them from cluttering your global workspace.

  # Inverse age transformation function
  iageF <- function(ptage.v, adult.age = 20) {
    y.idx <- which(ptage.v <= 0)
    a.idx <- which(ptage.v > 0)
    mage.v <- ptage.v
    mage.v[a.idx] <- ptage.v[a.idx] * (adult.age + 1) + adult.age
    mage.v[y.idx] <- exp(ptage.v[y.idx] + log(adult.age + 1)) - 1
    return(mage.v)
  }

  # Prediction of transformed age function
  PredTage <- function(beta.m, horvageP.lv) {
    # Match required CpGs with those in the input data
    map.idx <- match(names(horvageP.lv[[2]]), rownames(beta.m))
    rep.idx <- which(!is.na(map.idx))

    print(paste0("[Horvath2013] Number of represented Horvath2013 CpGs (max=",length(horvageP.lv[[2]]),")=",
               length(rep.idx)))

    # Subset beta matrix to only the matched CpGs
    tmpB.m <- beta.m[map.idx[rep.idx], , drop = FALSE]

    # Calculate the weighted sum for each sample
    intercept <- horvageP.lv[[1]]
    weights <- horvageP.lv[[2]][rep.idx]

    # Using matrix multiplication for efficiency instead of a for-loop
    tage.v <- as.vector(intercept + t(tmpB.m) %*% weights)

    return(tage.v)
  }

  # Step 3: Execute the prediction and transformation
  predTage.v <- PredTage(beta.m, horvageP.lv)
  predMage.v <- iageF(predTage.v)

  # Step 4: Add sample names to the output vector and return it
  names(predMage.v) <- colnames(beta.m)
  return(predMage.v)
}
