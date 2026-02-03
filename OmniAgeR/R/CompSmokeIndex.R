#' @title Compute a DNAm-based Smoking Index
#'
#' @description
#' Computes a DNA methylation-based smoking index (S-index) using 1,501 smoking-associated 
#' CpG sites. This index has been shown to correlate with cancer risk across multiple 
#' tissues and provides a robust epigenetic signature of smoking exposure.
#'
#' @details
#' The function calculates the score by first standardizing the methylation data (z-score 
#' transformation) across samples for each CpG site, and then computing a weighted 
#' average based on the provided signature coefficients.
#'
#' @param data.m A numeric matrix of normalized DNA methylation data (e.g.,
#'   beta values from BMIQ). Rows must be CpG probes (with rownames matching
#'   the IDs in `coeffSmkIdx.v`) and columns must be samples.
#'
#' @return A named \code{numeric} vector of smoking index for each sample in \code{data.m}.
#'
#'
#' @importFrom stats sd
#' @export
#'
#' @references
#' Teschendorff, Andrew E et al. 
#' Correlation of Smoking-Associated DNA Methylation Changes in Buccal Cells With DNA Methylation Changes in Epithelial Cancer
#' \emph{JAMA oncology} 2015
#' 
#' 
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' CompSmokeIndex_o <- CompSmokeIndex(hannum_bmiq_m)
#'

CompSmokeIndex <- function(data.m){
  data("coeffSmkIdx");
  common.v <- intersect(rownames(data.m),names(coeffSmkIdx.v));
  #print(paste("Found a fraction of ",length(common.v)/length(coeffSmkIdx.v)," of probes",sep=""));
  cat(paste0("[CompSmokeIndex] Number of represented CompSmokeIndex CpGs (max=",(length(coeffSmkIdx.v)),")=",
             length(common.v),"\n"))
  tmp.m <- data.m[match(common.v,rownames(data.m)),];
  sign.v <- coeffSmkIdx.v[match(common.v,names(coeffSmkIdx.v))];
  av.v <- rowMeans(tmp.m);
  sd.v <- apply(tmp.m,1,sd);
  v.idx <- which(sd.v>0);
  z.m <- (tmp.m[v.idx,]-av.v[v.idx])/sd.v[v.idx];
  si.v <- colMeans(sign.v[v.idx]*z.m);
  return(si.v);
}
