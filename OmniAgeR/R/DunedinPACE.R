
#' @title Calculate DunedinPACE
#'
#' @description
#' Calculates the DunedinPACE, an epigenetic biomarker that quantifies the pace of biological aging from a single blood sample.
#'
#' @details
#' This function calculates DunedinPACE scores from a DNA methylation beta value matrix using a robust pipeline.
#'
#' The process involves several key steps:
#' * **EPICv2 Array Handling**: Automatically detects and preprocesses Illumina EPICv2 array data, adjusting for its specific characteristics.
#' * **Probe Overlap Check**: Verifies that a sufficient number of required CpG probes are present in your data.
#' * **Missing Data Imputation**: Robustly handles both missing probes and missing values within existing probes.
#' * **Quantile Normalization**: Standardizes the data against a gold-standard reference to reduce batch effects and ensure comparability.
#' * **Score Calculation**: Computes the final DunedinPACE score using the pre-trained model weights.
#'
#' @param betas A numeric matrix of DNA methylation beta values. Rows should represent CpG sites and columns should represent individual samples.
#'
#' @param proportionOfProbesRequired A numeric value between 0 and 1 specifying the minimum proportion of required probes that must be present in the input `betas` matrix to proceed with the calculation. Defaults to `0.8` (80%). This is automatically adjusted to `0.7` for EPICv2 data.
#'
#' @return  A named numeric vector containing the calculated DunedinPACE for each sample.
#'
#' @references
#' Belsky DW, Caspi A, Corcoran DL, et al.
#' DunedinPACE, a DNA methylation biomarker of the pace of aging.
#' \emph{eLife} 2022
#'
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' DunedinPACE.out <- DunedinPACE(hannum_bmiq_m)



DunedinPACE <- function(betas, proportionOfProbesRequired=0.8 ) {
  data("DunedinPACECpG")
  requireNamespace("preprocessCore")

  # loop through models
  model_results <- lapply(mPACE_Models_list$model_names, function(model_name) {
    # make sure it has been converted to a matrix
    if( !is.numeric(as.matrix(betas)) ) { stop("[DunedinPACE] betas matrix/data.frame is not numeric!") }
    probeOverlap <- length(which(rownames(betas) %in% mPACE_Models_list$model_probes[[model_name]])) / length(mPACE_Models_list$model_probes[[model_name]])
    probeOverlap_background <- length(which(rownames(betas) %in% mPACE_Models_list$gold_standard_probes[[model_name]])) / length(mPACE_Models_list$gold_standard_probes[[model_name]])
    # make sure enough of the probes are present in the data file
    print(paste0("[DunedinPACE] Number of represented DunedinPACE CpGs (max=",length(mPACE_Models_list$model_probes[[model_name]]),")=",
               length(which(rownames(betas) %in% mPACE_Models_list$model_probes[[model_name]]))))
    print(paste0("[DunedinPACE] Number of represented Backgroup CpGs of DunedinPACE (max=",length(mPACE_Models_list$gold_standard_probes[[model_name]]),")=",
               length(which(rownames(betas) %in% mPACE_Models_list$gold_standard_probes[[model_name]]))))

    if( probeOverlap < proportionOfProbesRequired | probeOverlap_background < proportionOfProbesRequired ) {
      result <- rep(NA, ncol(betas))
      names(result) <- colnames(betas)
      result
    } else {
      # Work with a numeric matrix of betas
      betas.mat <- as.matrix(betas[which(rownames(betas) %in% mPACE_Models_list$gold_standard_probes[[model_name]]),])
      # If probes don't exist, we'll add them as rows of values based on their mean in the gold standard dataset
      probesNotInMatrix <- mPACE_Models_list$gold_standard_probes[[model_name]][which(mPACE_Models_list$gold_standard_probes[[model_name]] %in% rownames(betas.mat) == F)]
      if( length(probesNotInMatrix) > 0 ) {
        for( probe in probesNotInMatrix ) {
          tmp.mat <- matrix(0, nrow=1, ncol=ncol(betas.mat))
          rownames(tmp.mat) <- probe
          colnames(tmp.mat) <- colnames(betas.mat)
          tmp.mat[probe,] <- rep(mPACE_Models_list$gold_standard_means[[model_name]][probe], ncol(tmp.mat))
          betas.mat <- rbind(betas.mat, tmp.mat)
        }
      }

      # Identify samples with too many missing probes and remove them from the matrix
      samplesToRemove <- colnames(betas.mat)[which(apply(betas.mat, 2, function(x) { 1 - ( length(which(is.na(x))) / length(x) ) < proportionOfProbesRequired}))]
      if( length(samplesToRemove) > 0 ) {
        betas.mat <- betas.mat[,-which(colnames(betas.mat) %in% samplesToRemove)]
      }
      if(ncol(betas.mat) > 0) {
        # Identify missingness on a probe level
        pctValuesPresent <- apply( betas.mat, 1, function(x) { 1 - (length(which(is.na(x))) / length(x)) } )
        # If they're missing values, but less than the proportion required, we impute to the cohort mean
        probesToAdjust <- which(pctValuesPresent < 1 & pctValuesPresent >= proportionOfProbesRequired)
        if( length(probesToAdjust) > 0 ) {
          if( length(probesToAdjust) > 1 ) {
            betas.mat[probesToAdjust,] <- t(apply( betas.mat[probesToAdjust,], 1 , function(x) {
              x[is.na(x)] = mean( x, na.rm = TRUE )
              x
            }))
          } else {
            betas.mat[probesToAdjust,which(is.na(betas.mat[probesToAdjust,]))] <- mean(betas.mat[probesToAdjust,], na.rm=T)
          }
        }
        # If they're missing too many values, everyones value gets replaced with the mean from the Dunedin cohort
        if( length(which(pctValuesPresent < proportionOfProbesRequired)) > 0 ) {
          probesToReplaceWithMean <- rownames(betas.mat)[which(pctValuesPresent < proportionOfProbesRequired)]
          for( probe in probesToReplaceWithMean ) {
            betas.mat[probe,] <- rep(mPACE_Models_list$model_means[[model_name]][probe], ncol(betas.mat))
          }
        }

        # Normalize the matrix to the gold standard dataset
        betas.norm <- preprocessCore::normalize.quantiles.use.target(betas.mat, target=mPACE_Models_list$gold_standard_means[[model_name]])
        rownames(betas.norm) <- rownames(betas.mat)
        colnames(betas.norm) <- colnames(betas.mat)
        # Calculate score:
        score = mPACE_Models_list$model_intercept[[model_name]] + rowSums(t(betas.norm[mPACE_Models_list$model_probes[[model_name]],]) %*% diag(mPACE_Models_list$model_weights[[model_name]]))
        names(score) <- colnames(betas.norm)
        if( length(samplesToRemove) > 0 ) {
          score.tmp <- rep(NA, length(samplesToRemove))
          names(score.tmp) <- samplesToRemove
          score <- c(score, score.tmp)
        }
        score <- score[colnames(betas)]
        score
      } else {
        result <- rep(NA, ncol(betas.mat))
        names(result) <- colnames(betas.mat)
        result
      }
    }
  })
  names(model_results) <- mPACE_Models_list$model_names
  return(model_results$DunedinPACE)
}

