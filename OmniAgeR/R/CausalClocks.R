#' @title Calculate Causal, Damage, and Adaptation Epigenetic Clocks
#'
#' @description Calculates three related epigenetic clocks (Causal, Damage, Adaptation) from a DNA methylation beta value matrix. These clocks were developed to distinguish different aspects of aging.
#'
#' @param beta.m DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @details
#' This function implements the three epigenetic clocks described by Ying et al. (2024) to dissect biological aging into distinct components.
#' The Damage clock (DamAge) is designed to capture the accumulation of molecular damage.
#' The Adaptation clock (AdaptAge) reflects the body's adaptive responses to this damage.
#' The Causal clock (CausalAge) is enriched for CpGs with a causal effect on mortality.
#' Each score is calculated as a weighted linear sum of beta values from its specific set of CpG sites.
#'
#' @return A list containing the predicted scores for the "Causal", "Damage", and "Adaptation" clocks.
#'
#' @seealso The main function \code{\link{EpiAge}} can be used to calculate multiple clocks simultaneously.
#'
#' @references
#' Ying K, Liu H, Tarkhov AE, et al.
#' Causality-enriched epigenetic age uncouples damage and adaptation.
#' \emph{Nat Aging} 2024
#'
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' CausalClock.out <- CausalClock(hannum_bmiq_m)
#'



CausalClock <- function(beta.m){
  # --- Step 1: Load Coefficients ---
  data("CausalCpG");
  # --- Step 2: Calculate Scores for Each Clock ---
  est.lv <- list(); # Initialize a list to store results
  for(i in 1:length(causalClock.l)){

    commonCpG.v <- intersect(rownames(beta.m),causalClock.l[[i]][,1])
    print(paste0("[CausalClock] Number of represented ",c("CausalAge","DamAge","AdaptAge")[i]," CpGs (max=",(nrow(causalClock.l[[i]])-1),")=",
               length(commonCpG.v)))


    tmp.m <- beta.m[match(commonCpG.v,rownames(beta.m)),];
    coefs.v <- as.numeric(causalClock.l[[i]][match(commonCpG.v,causalClock.l[[i]][,1]),2]);
    est.v <- causalClock.l[[i]][1,2] + as.vector(t(tmp.m) %*% matrix(coefs.v,ncol=1));
    est.lv[[i]] <- est.v;

  }
  #names(est.lv) <- names(causalClock.l);
  names(est.lv) <- c("CausalAge","DamAge","AdaptAge")
  return(est.lv);
}
