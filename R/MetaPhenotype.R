#' Classify tumor into three metabolic phenotype for seven pathways
#'
#' @param TC
#' @param pvalue
#'
#' @return
#' @export
#'
#' @importFrom clusterProfiler GSEA
#'
#' @importFrom dplyr full_join %>% select
#'
#' @import tidyverse
#'
#' @examples
MetaPhenotype<-function(TC,pvalue){
  Metabolite=data.frame(ID=NA,NES=NA,p.adjust=NA,samplesID=NA,zscore=NA)
  exp=t(TC)%>%as.data.frame()
  for(i in 1:ncol(TC)){
    a=i
    gene=t(exp[i,])
    genelist=gene[,1]
    names(genelist)=row.names(gene)
    genelist=sort(genelist,decreasing=T)
    GSEA=GSEA(genelist,TERM2GENE=gmt,minGSSize=1,maxGSSize=1000,pvalueCutoff=1,by="fgsea",nPermSimple = 10000)
    table=GSEA@result%>%select("ID","NES","p.adjust")
    table$samplesID=row.names(exp)[i]
    for (j in names(GSEA@geneSets))
    {id=intersect(as.character(GSEA@geneSets[[j]]),rownames(gene))
    zscore=sum(gene[c(id),])
    table$zscore[table$ID==j]=zscore}
    Metabolite=rbind(Metabolite,table)
    print(paste0("Finish ",a," samples"))}
  Metabolite=Metabolite[-1,]
  GSEA_results=Metabolite
  Metabolite$Amino_acid_TC[Metabolite$ID=="Amino_acid"&Metabolite$p.adjust>=pvalue]="Neutral"
  mean=mean(Metabolite$zscore[Metabolite$ID=="Amino_acid"&Metabolite$p.adjust<pvalue])
  Metabolite$Amino_acid_TC[Metabolite$ID=="Amino_acid"&Metabolite$p.adjust<pvalue&Metabolite$zscore>mean]="Upregulated"
  Metabolite$Amino_acid_TC[Metabolite$ID=="Amino_acid"&Metabolite$p.adjust<pvalue&Metabolite$zscore<=mean]="Downregulated"
  table(Metabolite$Amino_acid_TC[Metabolite$ID=="Amino_acid"])

  Metabolite$Lipid_TC[Metabolite$ID=="Lipid"&Metabolite$p.adjust>=pvalue]="Neutral"
  mean=mean(Metabolite$NES[Metabolite$ID=="Lipid"&Metabolite$p.adjust<pvalue])
  Metabolite$Lipid_TC[Metabolite$ID=="Lipid"&Metabolite$p.adjust<pvalue&Metabolite$NES>mean]="Upregulated"
  Metabolite$Lipid_TC[Metabolite$ID=="Lipid"&Metabolite$p.adjust<pvalue&Metabolite$NES<=mean]="Downregulated"
  table(Metabolite$Lipid_TC[Metabolite$ID=="Lipid"])

  Metabolite$Carbohydrate_TC[Metabolite$ID=="Carbohydrate"&Metabolite$p.adjust>=pvalue]="Neutral"
  mean=mean(Metabolite$zscore[Metabolite$ID=="Carbohydrate"&Metabolite$p.adjust<pvalue])
  Metabolite$Carbohydrate_TC[Metabolite$ID=="Carbohydrate"&Metabolite$p.adjust<pvalue&Metabolite$zscore>mean]="Upregulated"
  Metabolite$Carbohydrate_TC[Metabolite$ID=="Carbohydrate"&Metabolite$p.adjust<pvalue&Metabolite$zscore<=mean]="Downregulated"
  table(Metabolite$Carbohydrate_TC[Metabolite$ID=="Carbohydrate"])

  Metabolite$TCA_cycle_TC[Metabolite$ID=="TCA cycle"&Metabolite$p.adjust>=pvalue]="Neutral"
  mean=mean(Metabolite$zscore[Metabolite$ID=="TCA cycle"&Metabolite$p.adjust<pvalue])
  Metabolite$TCA_cycle_TC[Metabolite$ID=="TCA cycle"&Metabolite$p.adjust<pvalue&Metabolite$zscore>mean]="Upregulated"
  Metabolite$TCA_cycle_TC[Metabolite$ID=="TCA cycle"&Metabolite$p.adjust<pvalue&Metabolite$zscore<=mean]="Downregulated"
  table(Metabolite$TCA_cycle_TC[Metabolite$ID=="TCA cycle"])

  Metabolite$Energy_TC[Metabolite$ID=="Energy"&Metabolite$p.adjust>=pvalue]="Neutral"
  mean=mean(Metabolite$zscore[Metabolite$ID=="Energy"&Metabolite$p.adjust<pvalue])
  Metabolite$Energy_TC[Metabolite$ID=="Energy"&Metabolite$p.adjust<pvalue&Metabolite$zscore>mean]="Upregulated"
  Metabolite$Energy_TC[Metabolite$ID=="Energy"&Metabolite$p.adjust<pvalue&Metabolite$zscore<=mean]="Downregulated"
  table(Metabolite$Energy_TC[Metabolite$ID=="Energy"])

  Metabolite$Nucleotide_TC[Metabolite$ID=="Nucleotide"&Metabolite$p.adjust>=pvalue]="Neutral"
  mean=mean(Metabolite$zscore[Metabolite$ID=="Nucleotide"&Metabolite$p.adjust<pvalue])
  Metabolite$Nucleotide_TC[Metabolite$ID=="Nucleotide"&Metabolite$p.adjust<pvalue&Metabolite$zscore>mean]="Upregulated"
  Metabolite$Nucleotide_TC[Metabolite$ID=="Nucleotide"&Metabolite$p.adjust<pvalue&Metabolite$zscore<=mean]="Downregulated"
  table(Metabolite$Nucleotide_TC[Metabolite$ID=="Nucleotide"])

  Metabolite$Vitamin_TC[Metabolite$ID=="Vitamin_cofactor"&Metabolite$p.adjust>=pvalue]="Neutral"
  mean=mean(Metabolite$zscore[Metabolite$ID=="Vitamin_cofactor"&Metabolite$p.adjust<pvalue])
  Metabolite$Vitamin_TC[Metabolite$ID=="Vitamin_cofactor"&Metabolite$p.adjust<pvalue&Metabolite$zscore>mean]="Upregulated"
  Metabolite$Vitamin_TC[Metabolite$ID=="Vitamin_cofactor"&Metabolite$p.adjust<pvalue&Metabolite$zscore<=mean]="Downregulated"
  table(Metabolite$Vitamin_TC[Metabolite$ID=="Vitamin_cofactor"])

  Amino_acid=Metabolite[Metabolite$ID=="Amino_acid",c("samplesID","Amino_acid_TC")]
  colnames(Amino_acid)=c("samplesID","Amino_acid_TC")
  Lipid=Metabolite[Metabolite$ID=="Lipid",c("samplesID","Lipid_TC")]
  colnames(Lipid)=c("samplesID","Lipid_TC")
  Carbohydrate=Metabolite[Metabolite$ID=="Carbohydrate",c("samplesID","Carbohydrate_TC")]
  colnames(Carbohydrate)=c("samplesID","Carbohydrate_TC")
  TCA_cycle=Metabolite[Metabolite$ID=="TCA cycle",c("samplesID","TCA_cycle_TC")]
  colnames(TCA_cycle)=c("samplesID","TCA_cycle_TC")
  Energy=Metabolite[Metabolite$ID=="Energy",c("samplesID","Energy_TC")]
  colnames(Energy)=c("samplesID","Energy_TC")
  Nucleotide=Metabolite[Metabolite$ID=="Nucleotide",c("samplesID","Nucleotide_TC")]
  colnames(Nucleotide)=c("samplesID","Nucleotide_TC")
  Vitamin=Metabolite[Metabolite$ID=="Vitamin_cofactor",c("samplesID","Vitamin_TC")]
  colnames(Vitamin)=c("samplesID","Vitamin_TC")
  Metabolite_TC=Amino_acid%>%full_join(Lipid,by="samplesID")%>%
    full_join(Carbohydrate,by="samplesID")%>%full_join(TCA_cycle,by="samplesID")%>%
    full_join(Energy,by="samplesID")%>%full_join(Nucleotide,by="samplesID")%>%full_join(Vitamin,by="samplesID")
  sapply(Metabolite_TC, function(x) sum(is.na(x)))
  Metabolite_TC[is.na(Metabolite_TC)]="Neutral"
  return(list(GSEA_results=GSEA_results,Metabolite_phenotype=Metabolite_TC))
}
