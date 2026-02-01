#' Calculate CHIP-related Methylation Scores
#'
#' @description
#' Computes CHIP (Clonal Hematopoiesis of Indeterminate Potential) scores
#' based on predefined CpG signatures and a normalized DNA methylation matrix.
#'
#' @details
#' For each CHIP signature (e.g., DNMT3A, TET2), it identifies common CpGs, calculates Z-scores for each
#' probe across all samples, and then computes the score as:
#' (Mean Z-score of positive probes) - (Mean Z-score of negative probes).
#'
#' Probes with zero variance across samples are excluded.
#'
#' @param nbeta.m A numeric matrix of normalized DNA methylation data (e.g.,
#'   beta values from BMIQ). Rows must be CpG probes (with rownames matching
#'   the IDs in `chipCpG.lv`) and columns must be samples.
#'
#' @return
#' A list of the same length as `chipCpG.lv`. Each element is a numeric
#' vector containing the calculated CHIP score for each sample (in the same
#' order as the columns of `nbeta.m`).
#'
#' @importFrom stats sd
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' CompCHIP.o <- CompCHIP(hannum_bmiq_m)
#'


CompCHIP <- function(nbeta.m){
  data("chipCpG")
  score.lv <- list();
  for(i in 1:length(chipCpG.lv)){
    common.v <- intersect(names(chipCpG.lv[[i]]),rownames(nbeta.m));
    print(paste0("[CHIP] Number of represented ", names(chipCpG.lv)[i], " CpGs (max=",
                   length(chipCpG.lv[[i]]), ")=", length(common.v)))


    #print(paste("A fraction ",length(common.v)/length(chipCpG.lv[[i]])," of CRP probes have been found. Number found is=",length(common.v),sep=""));
    match(common.v,names(chipCpG.lv[[i]])) -> map1.idx;
    match(common.v,rownames(nbeta.m)) -> map2.idx;
    tmp.m <- nbeta.m[map2.idx,];
    sd.v <- apply(tmp.m,1,sd);
    z.m <- (tmp.m - rowMeans(tmp.m))/sd.v;
    sign.v <- sign(chipCpG.lv[[i]][map1.idx]);
    pos.idx <- which(sign.v==1);
    neg.idx <- which(sign.v==-1);
    scoreP.v <- rep(0,ncol(tmp.m));
    scoreN.v <- rep(0,ncol(tmp.m));
    if(length(pos.idx)==1){
      scoreP.v <- z.m[pos.idx,];
    }
    else if (length(pos.idx)>1){
      scoreP.v <- colMeans(z.m[pos.idx,]);
    }

    if(length(neg.idx)==1){
      scoreN.v <- z.m[neg.idx,];
    }
    else if (length(neg.idx)>1){
      scoreN.v <- colMeans(z.m[neg.idx,]);
    }

    score.lv[[i]] <- scoreP.v - scoreN.v;
  }
  names(score.lv) <- names(chipCpG.lv);
  return(score.lv);
} ### EOF
