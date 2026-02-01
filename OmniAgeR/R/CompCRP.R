#' @title Compute a DNAm-based C-Reactive Protein (CRP) Score
#'
#' @description
#' Computes a DNAm-based surrogate proxy for C-Reactive Protein (CRP) levels, a key biomarker for systemic inflammation (inflammaging).
#'
#' @details
#' This function calculates a DNAm-based CRP proxy score from a normalized beta-valued DNAm data matrix. The score is derived by correlating the methylation values with the signs of pre-defined CpG weights from two internal signatures: (i) crp: a vector containing 1765 CpGs and weights associated with CRP as inferred from a large meta-analysis over 20,000 blood samples where both DNAm and CRP was measured, and where no adjustment for variations between memory and naive T-cells was made. (ii) intCRP: contains a 62 CpG CRP-signature that is adjusted for variations between memory and naive T-cells. The resulting output is a relative score that captures the methylation pattern associated with CRP levels.
#'
#' @param nbeta.m normalized beta-valued DNAm data matrix with rownames the CpG identifiers. Missing values should be imputed.
#'
#' @return A named list containing two numeric vectors of CRP estimates, one for each signature provided in `crpCpG.lv` (`crp` and `intCRP`).
#'
#' @author Andrew Teschendorff (andrew@sinh.ac.cn)
#'
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' CompCRP.out <- CompCRP(hannum_bmiq_m)




CompCRP <- function(nbeta.m){
  data("crpCpG")
  score.lv <- list();
  for(i in 1:length(crpCpG.lv)){
    sigCRP.v <- crpCpG.lv[[i]];
    common.v <- intersect(names(sigCRP.v),rownames(nbeta.m));

    cat(paste0("[CRP] Number of represented ",c("CRP","intCRP")[i]," CpGs (max=",(length(sigCRP.v)),")=",
               length(common.v),"\n"))
    #print(paste("A fraction ",length(common.v)/length(sigCRP.v)," of CRP probes have been found. Number found is=",length(common.v),sep=""));
    match(common.v,names(sigCRP.v)) -> map1.idx;
    match(common.v,rownames(nbeta.m)) -> map2.idx;
    tmp.m <- nbeta.m[map2.idx,];
    z.m <- (tmp.m - rowMeans(tmp.m))/apply(tmp.m,1,sd);
    score.lv[[i]] <- as.vector(cor(z.m,sign(sigCRP.v[map1.idx])));
  }
  names(score.lv) <- names(crpCpG.lv);
  return(score.lv);
}



