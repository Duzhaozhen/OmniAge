#' @title
#' Estimate stemTOC score
#'
#' @aliases stemTOC
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will return the stemTOC score.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @details
#' The function will return the 0.95 upper quantile of the 371 stemTOC CpGs. Compared to stemTOCvitro CpGs, the stemTOC CpGs are filtered for significant DNA hypermethylation with chronological age in large in-vivo datasets
#'
#' @return The stemTOC score of each sample.
#'
#' @references
#' Zhu, T., Tong, H., Du, Z. et al.
#' An improved epigenetic counter to track mitotic age in normal and precancerous tissues.
#' \emph{Nat Commun} 2024
#'
#'
#' @examples
#' download_OmniAgeR_example("LungInv")
#' load_OmniAgeR_example("LungInv")
#' stemTOC.v<-stemTOC(data.m = bmiq.m)
#'
#' @export
#'


stemTOC <- function(data.m){
  data('epiTOCcpgs3')
  common.v <- intersect(rownames(data.m),epiTOCcpgs3.v);
  print(paste("[stemTOC] Number of represented stemTOC CpGs (max=371)=",length(common.v),sep=""))
  stemTOC.v <- apply(data.m[match(epiTOCcpgs3.v,rownames(data.m)),],2,quantile,0.95,na.rm=T)
  return(stemTOC.v);
}

