###plot boxplot of several genes in one plot

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/onTAD/N4T7_ontad")
options(scipen=200)
library(data.table)
library(dplyr)
require(cowplot)
require(tidyverse)
require(ggplot2)
require(ggsci)
require(ggpubr)

###our own data
expr <- fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments/all_sample_gene_expr.txt")
## calculate common genes
A2B_NT <- c("HOXD1","HOXD-AS2","HOXD3","HOXD4","HOXD8","HOXD9","HOXD10","HOXD11","HOXD12","HOXD13","HAGLR","EVX2","MTX2")
B2A_NT <- c("ITGB1BP1","CPSF3","IAH1","ADAM17","YWHAQ","TAF1B","GRHL1","ASAP2","MBOAT2")
A2B_TM <- c("FOXF1","MTHFSD","FENDRR","IRF8")
A2B_NT_mat <- expr[expr$V1%in%unique(A2B_NT),]
B2A_NT_mat <- expr[expr$V1%in%unique(B2A_NT),]
A2B_TM_mat <- expr[expr$V1%in%unique(A2B_TM),]
##A2B_NT
#gene_expr <- data.frame(scale(t(A2B_NT[,c(2,15,3,16,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,13,26,14,27)]),center=T,scale=F))
gene_expr <- data.frame(t(A2B_NT_mat[,c(2,15,3,16,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,13,26,14,27)]))
colnames(gene_expr) <- A2B_NT_mat$V1
gene_expr <- gene_expr[,c(A2B_NT)]
gene_expr$gene <- rownames(gene_expr)
test <- gather(gene_expr,name,expression,-gene)
test <- separate(test,"gene",c("1","Patient","Group"),"-")[,2:5]
test$Patient <- paste0("CRC-",test$Patient)
colnames(test) <- c("Patient","Tissue","Gene","Expression")
p <- ggboxplot(test, x="Gene",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=45,
               color="Tissue",palette=c("lightblue","orange"),add="jitter",add.params=list(size=0.5))
p# + stat_compare_means(aes(group = Tissue),label="p.signif",method="t.test",paired = T)
##B2A_NT
gene_expr <- data.frame(t(B2A_NT_mat[,c(2,15,3,16,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,13,26,14,27)]))
colnames(gene_expr) <- B2A_NT_mat$V1
gene_expr <- gene_expr[,c(B2A_NT)]
gene_expr$gene <- rownames(gene_expr)
test <- gather(gene_expr,name,expression,-gene)
test <- separate(test,"gene",c("1","Patient","Group"),"-")[,2:5]
test$Patient <- paste0("CRC-",test$Patient)
colnames(test) <- c("Patient","Tissue","Gene","Expression")
p <- ggboxplot(test, x="Gene",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=45,
               color="Tissue",palette=c("lightblue","orange"),add="jitter",add.params=list(size=0.5))
p + stat_compare_means(aes(group = Tissue),label="p.signif",method="t.test",paired = T)
##A2B_TM
gene_expr <- data.frame(t(A2B_TM_mat[,c(16,29,17,30,18,31,19,32,26,33,27,34)]),center=T,scale=F)
colnames(gene_expr) <- A2B_TM_mat$V1
gene_expr <- gene_expr[,c(A2B_TM)]
gene_expr$gene <- rownames(gene_expr)
test <- gather(gene_expr,name,expression,-gene)
test <- separate(test,"gene",c("1","Patient","Group"),"-")[,2:5]
test$Patient <- paste0("CRC-",test$Patient)
colnames(test) <- c("Patient","Tissue","Gene","Expression")
test[test=="M"] <- "Z" 
p <- ggboxplot(test, x="Gene",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=45,
               color="Tissue",palette=c("orange","pink"),add="jitter",add.params=list(size=0.5))
p + stat_compare_means(aes(group = Tissue),label="p.signif",method="t.test",paired = T)

###TCGA data
ref <- as.data.frame(fread("/lustre/user/liclab/ganjb/resource/hg19_refgene/geneID_ref.txt",select=c(10,2)))
colnames(ref)<-c("ENSEMBL","SYMBOL")
expr <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA/tcga.txt"))
for (i in 2:ncol(expr)){
  expr[,i] <- 200000*expr[,i]/sum(expr[,i])
}
apply(expr[,2:10],2,sum)
meta <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA/meta_tcga.txt"))
colnames(expr)[1] <- "ENSEMBL"
expr$ENSEMBL <- separate(expr,"ENSEMBL",".")[[1]]

clu1_pro <- ref[ref$SYMBOL%in%c(A2B_NT),]
clu1_pro <- clu1_pro[is.na(clu1_pro$ENSEMBL)==F,]
mygene <- as.data.frame(merge(clu1_pro,expr,"ENSEMBL"))
tcga_N <- mygene[,paste0(c(meta[meta$tissue=="Normal",1]))]
tcga_T <- mygene[,paste0(c(meta[meta$tissue=="Tumor",1]))]
mygene <- t(cbind(mygene$SYMBOL,tcga_N,tcga_T))
colnames(mygene) <- c(mygene[1,])
mygene <- as.data.frame(mygene[-1,])
mygene <- mygene[,c("FOXF1","MTHFSD","FENDRR","IRF8")]
mygene$gene <- rownames(mygene)
test <- gather(mygene,name,expression,-gene)
test <- separate(test,"gene",c("1","name1","name2","Group"),"-")[,2:6]
test$patient <- paste0(test$name1,"-",test$name2)
test <- test[,-c(1:2)]
#colnames(test) <- c(A2B_NT)
test <- test[test$Group%in%c("11A","01A"),]
test[test=="11A"] <- "N"
test[test=="01A"] <- "T"
test$expression <- as.numeric(test$expression)
p <- ggboxplot(test, x="name",y="expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=45,
               color="Group",palette=c("lightblue","orange"),add="jitter",add.params=list(size=0.2))
p

