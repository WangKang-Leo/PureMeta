library(data.table)
tumor=data.frame(fread("E:/Projects/Proteogenomic_PROMIX/Manuscript/Revision/Analyses/Figure3/GSE62944_RAW/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt"),row.names=1)
normal=data.frame(fread("E:/Projects/Proteogenomic_PROMIX/Manuscript/Revision/Analyses/Figure3/GSE62944_RAW/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_TPM.txt/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_TPM.txt"),row.names=1)
tumor_meta=data.frame(fread("E:/Projects/Proteogenomic_PROMIX/Manuscript/Revision/Analyses/Figure3/GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt/GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt"))
colnames(tumor_meta)=c("sampleID","type")
tumor_meta$sampleID=gsub("-", ".", tumor_meta$sampleID)
normal_meta=data.frame(fread("E:/Projects/Proteogenomic_PROMIX/Manuscript/Revision/Analyses/Figure3/GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt/GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt"))
colnames(normal_meta)=c("sampleID","type")
normal_meta$sampleID=gsub("-", ".", normal_meta$sampleID)
tumor_brca=tumor[,tumor_meta$sampleID[tumor_meta$type=="BRCA"]]
normal_brca=normal[,normal_meta$sampleID[normal_meta$type=="BRCA"]]
saveRDS(tumor_brca[,1:10],file="E:/R_Dev/PureMeta/data_raw/TCGA_BRCA_log2TPM.rds")
saveRDS(normal_brca[,1:5],file="E:/R_Dev/PureMeta/data_raw/TCGA_normal_breast_log2TPM.rds")

library(readxl);library(clusterProfiler);library(data.table)
gmt=read.gmt("E:/Projects/Proteogenomic_PROMIX/dataprocessing/step6/Metabolism_cluster/metabolite.gmt")
uni_gene=unique(gmt$gene)
exp=as.data.frame(fread("E:/Projects/Proteogenomic_PROMIX/dataprocessing/step6/Metabolism_cluster/promix_non_normal.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
cg=names(tail(sort(apply(exp,1,sd)),5000))  #########变化最大得2000gene
genelist=Reduce(intersect,list(row.names(exp),uni_gene))
genelist=Reduce(union,list(genelist,cg))
exp=exp[genelist,]
normal=exp[,c("147_op","315_op","224_op","611_OP","311_op")]
tumor=exp[,!(colnames(exp)%in%colnames(normal))]

library(usethis)
tumor=tumor[,1:30]
use_data(tumor,compress = "gzip",overwrite = TRUE)
normal=normal
use_data(normal,compress = "gzip",overwrite = TRUE)
