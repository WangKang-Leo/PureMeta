
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
  ISOpureS1model=ISOpure.step1.CPE(as.matrix(tumor),as.matrix(normal))
  ISOpureS2model=ISOpure.step2.PPE(as.matrix(tumor),as.matrix(normal),ISOpureS1model)
  purity=ISOpureS1model$alphapurities
  tumor_cell_GEP=as.data.frame(ISOpureS2model$cc_cancerprofiles)
  row.names(tumor_cell_GEP)=row.names(tumor)
  colnames(tumor_cell_GEP)=row.names(normal)
  return(list(purity,tumor_cell_GEP))
}

#' Title
#' @param data GEP of tumor cells
#'
#' @return tumor cell-based metabolic phenotype of each sample
#'
#' @importFrom clusterProfiler GSEA
#'
#' @importFrom gson read.gmt
#'
#' @examples
MetaPhenotype<-function(){
  gmt=read.gmt("E:/Projects/Proteogenomic_PROMIX/dataprocessing/Validation/Korean_RNA_seq/metabolite_genelist.gmt")
  Metabolite_group=data.frame(ID=NA,NES=NA,p.adjust=NA,samplesID=NA,zscore=NA)
  exp=mydata
  for(i in 1:223){
    a=i
    gene=t(exp[i,])
    genelist=gene[,1]
    names(genelist)=row.names(gene)
    genelist=sort(genelist,decreasing=T)
    #GSEA=fgsea(gmt,genelist,nperm=100000,minSize=1,maxSize=10000,gseaParam=1)
    #table=GSEA%>%select("pathway","NES","padj")
    GSEA=GSEA(genelist,TERM2GENE=gmt,minGSSize=1,maxGSSize=1000,pvalueCutoff=1,by="fgsea",nPermSimple = 10000)
    table=GSEA@result%>%select("ID","NES","p.adjust")
    table$samplesID=row.names(exp)[i]
    for (j in names(GSEA@geneSets))
    {id=intersect(as.character(GSEA@geneSets[[j]]),rownames(gene))
    zscore=sum(gene[c(id),])
    table$zscore[table$ID==j]=zscore}
    Metabolite_group=rbind(Metabolite_group,table)
    print(paste0("Finish:",a))}
  return(Metabolite_group[-1,])
}



