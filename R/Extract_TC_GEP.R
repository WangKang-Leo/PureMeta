
#' Extracting Tumor Cell-based GEP
#' @param tumor  Bulk GEP of tumor tissues
#' @param normal  Bulk GEP of adjacent normal tissues, or post-treatment tissues with pathologic stage pT0N0
#'
#' @return GEP of tumor cells
#'
#' @importFrom ISOpureR ISOpure.step1.CPE ISOpure.step2.PPE
#'
#' @export
#'
#' @examples
#'
Extract_TC_GEP <- function(tumor,normal) {
  uni_gene=unique(gmt$gene)
  genelist=Reduce(intersect,list(row.names(tumor),uni_gene))
  pro=intersect(row.names(tumor),coding_gene$symbol)
  genelist=Reduce(union,list(genelist,pro))
  tumor=tumor[genelist,]
  normal=normal[genelist,]
  ISOpureS1model=ISOpure.step1.CPE(as.matrix(tumor),as.matrix(normal))
  ISOpureS2model=ISOpure.step2.PPE(as.matrix(tumor),as.matrix(normal),ISOpureS1model)
  purity=ISOpureS1model$alphapurities
  tumor_cell_GEP=as.data.frame(ISOpureS2model$cc_cancerprofiles)
  row.names(tumor_cell_GEP)=row.names(tumor)
  colnames(tumor_cell_GEP)=colnames(tumor)
  return(list(purity=purity,tumor_cell_GEP=tumor_cell_GEP))
}
