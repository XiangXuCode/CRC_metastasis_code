#used to find diff loop between samples; plot heatmap for single loop; find continously altered loops

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/merged_validPairs/hic_results/data/loop/diffloop_apa")
library(data.table)
library(dplyr)
library(VennDiagram)
library(pheatmap)
library(parallel)
options(scipen=200)

###Diff analysis:(loop lost/gained)+(apa diff(1.5 fold) of common loop)
loop_CN <- as.data.frame(fread("../CN/merged_loops.bedpe",select=1:6))
colnames(loop_CN) <- c("chr1","start1","end1","chr2","start2","end2")
loop_CN$len <- loop_CN$end2-loop_CN$end1
loop_CN$chr1 <- paste0("chr",loop_CN$chr1);loop_CN$chr2 <- paste0("chr",loop_CN$chr2)
loop_CN <- loop_CN[loop_CN$len <= 4000000 & loop_CN$len >= 100000,]
#plot(density(loop_CN$len))
loop_CT <- as.data.frame(fread("../CT/merged_loops.bedpe",select=1:6))
colnames(loop_CT) <- c("chr1","start1","end1","chr2","start2","end2")
loop_CT$len <- loop_CT$end2-loop_CT$end1
loop_CT <- loop_CT[loop_CT$len <= 4000000 & loop_CT$len >= 100000,]
loop_CT$chr1 <- paste0("chr",loop_CT$chr1);loop_CT$chr2 <- paste0("chr",loop_CT$chr2)
#plot(density(loop_CT$len))
loop_LT <- as.data.frame(fread("../LT/merged_loops.bedpe",select=1:6))
colnames(loop_LT) <- c("chr1","start1","end1","chr2","start2","end2")
loop_LT$len <- loop_LT$end2-loop_LT$end1
loop_LT <- loop_LT[loop_LT$len <= 4000000 & loop_LT$len >= 100000,]
loop_LT$chr1 <- paste0("chr",loop_LT$chr1);loop_LT$chr2 <- paste0("chr",loop_LT$chr2)
#plot(density(loop_LT$len))
loop_CN <- loop_CN[loop_CN$chr1!="chrY",]
loop_CT <- loop_CT[loop_CT$chr1!="chrY",]
loop_LT <- loop_LT[loop_LT$chr1!="chrY",]
boxplot(list(N=loop_CN$len,T=loop_CT$len,M=loop_LT$len),outline=F,col=c("lightblue","orange","pink"),
        cex.axis=2)
head(loop_CN)
fwrite(loop_CN[,c(1,2,6)],"CN_loop.bed",sep="\t",quote=F,col.names = F,row.names = F)
fwrite(loop_CT[,c(1,2,6)],"CT_loop.bed",sep="\t",quote=F,col.names = F,row.names = F)
fwrite(loop_LT[,c(1,2,6)],"LT_loop.bed",sep="\t",quote=F,col.names = F,row.names = F)

#bedtools intersect -a CN_loop_sorted.bed -b CT_loop_sorted.bed -f 0.90 -r -wa -wb > NT_com_loop.txt
#bedtools intersect -a CT_loop_sorted.bed -b LT_loop_sorted.bed -f 0.90 -r -wa -wb > TM_com_loop.txt
NT_com <- as.data.frame(fread("./NT_com_loop.txt"))
colnames(NT_com) <- c("chr1","start1","end2","chr1","start1","end2")
NT_T <- anti_join(loop_CT[,c(1,2,6)],NT_com[,4:6])
NT_N <- anti_join(loop_CN[,c(1,2,6)],NT_com[,1:3])

loop <- as.data.frame(fread("./NTM_loop_apa.txt",select=c(1:9)))
loop <- loop[,-c(3:5)]
colnames(loop)[1:3] <- c("chr1","start1","end2")
NT_com_N <- merge(loop,NT_com[,1:3],c("chr1","start1","end2"))
NT_com_T <- merge(loop,NT_com[,4:6],c("chr1","start1","end2"))
NT_com <- unique(rbind(NT_com_N,NT_com_T))
#filter loops whose neighboring loop also have same trend 
#NT_T
NT_T_h <- NT_com[NT_com$T>1.5*NT_com$N,1:3]
NT_N_h <- NT_com[NT_com$N>1.5*NT_com$T,1:3]
NT_com <- as.data.frame(fread("./NT_com_loop.txt"))
colnames(NT_com) <- c("chr1_N","start1_N","end2_N","chr1_T","start1_T","end2_T")
colnames(NT_T_h) <- c("chr1_T","start1_T","end2_T")
test_T <- merge(NT_T_h,NT_com,c("chr1_T","start1_T","end2_T"))
colnames(loop)[1:3] <- c("chr1_N","start1_N","end2_N")
test_T <- merge(loop,test_T[,4:6],c("chr1_N","start1_N","end2_N"))
test_T <- test_T[test_T$T>1.5*test_T$N,]
colnames(NT_T_h) <- c("chr1_N","start1_N","end2_N")
test_T2 <- merge(NT_T_h,NT_com,c("chr1_N","start1_N","end2_N"))
colnames(loop)[1:3] <- c("chr1_N","start1_N","end2_N")
test_T2 <- merge(loop,test_T2[,1:3],c("chr1_N","start1_N","end2_N"))
test_T2 <- test_T[test_T2$T>1.5*test_T2$N,]
NT_T_h <- unique(inner_join(test_T,test_T2))
NT_T <- as.data.frame(rbind(as.matrix(NT_T),as.matrix(NT_T_h[,1:3])))
#NT_N
colnames(NT_N_h) <- c("chr1_N","start1_N","end2_N")
test_N <- merge(NT_N_h,NT_com,c("chr1_N","start1_N","end2_N"))
colnames(loop)[1:3] <- c("chr1_T","start1_T","end2_T")
test_N <- merge(loop,test_N[,4:6],c("chr1_T","start1_T","end2_T"))
test_N <- test_N[test_N$N>1.5*test_N$T,]
colnames(NT_N_h) <- c("chr1_T","start1_T","end2_T")
test_N2 <- merge(NT_N_h,NT_com,c("chr1_T","start1_T","end2_T"))
colnames(loop)[1:3] <- c("chr1_T","start1_T","end2_T")
test_N2 <- merge(loop,test_N2[,1:3],c("chr1_T","start1_T","end2_T"))
test_N2 <- test_N2[test_N2$N>1.5*test_N2$T,]
NT_N_h <- unique(inner_join(test_N,test_N2))
NT_N <- as.data.frame(rbind(as.matrix(NT_N),as.matrix(NT_N_h[,1:3])))
fwrite(NT_N,"./test/NT_N_loop.txt",sep="\t",quote=F,col.names=T)
fwrite(NT_T,"./test/NT_T_loop.txt",sep="\t",quote=F,col.names=T)
NT_N <- fread("./test/NT_N_loop.txt")
NT_T <- fread("./test/NT_T_loop.txt")
NT_all <- unique(rbind(loop_CN[,c(1,2,6)],loop_CT[,c(1,2,6)]))
NT_com <- anti_join(anti_join(NT_all,NT_N),NT_T)
NT_N$chr2 <- NT_N$chr1;NT_N$end1 <- as.numeric(NT_N$start1)+10000;NT_N$start2 <- as.numeric(NT_N$end2)-10000;
NT_T$chr2 <- NT_T$chr1;NT_T$end1 <- as.numeric(NT_T$start1)+10000;NT_T$start2 <- as.numeric(NT_T$end2)-10000;
NT_com$chr2 <- NT_com$chr1;NT_com$end1 <- as.numeric(NT_com$start1)+10000;NT_com$start2 <- as.numeric(NT_com$end2)-10000;
NT_N <- rbind(as.matrix(NT_N[,c("chr1","start1","end1")]),as.matrix(NT_N[,c("chr2","start2","end2")]))
NT_T <- rbind(as.matrix(NT_T[,c("chr1","start1","end1")]),as.matrix(NT_T[,c("chr2","start2","end2")]))
NT_com <- rbind(as.matrix(NT_com[,c("chr1","start1","end1")]),as.matrix(NT_com[,c("chr2","start2","end2")]))
fwrite(NT_N,"./test/NT_N_bou.bed",sep="\t",quote=F,col.names=F)
fwrite(NT_T,"./test/NT_T_bou.bed",sep="\t",quote=F,col.names=F)
fwrite(NT_com,"./test/NT_T_com_bou.bed",sep="\t",quote=F,col.names=F)

#TvsM
TM_com <- as.data.frame(fread("./TM_com_loop.txt"))
colnames(TM_com) <- c("chr1","start1","end2","chr1","start1","end2")
TM_M <- anti_join(loop_LT[,c(1,2,6)],TM_com[,4:6])
TM_T <- anti_join(loop_CT[,c(1,2,6)],TM_com[,1:3])

loop <- as.data.frame(fread("./NTM_loop_apa.txt",select=c(1:9)))
loop <- loop[,-c(3:5)]
colnames(loop)[1:3] <- c("chr1","start1","end2")
TM_com_T <- merge(loop,TM_com[,1:3],c("chr1","start1","end2"))
TM_com_M <- merge(loop,TM_com[,4:6],c("chr1","start1","end2"))
TM_com <- unique(rbind(TM_com_T,TM_com_M))
#filter loops whose neighboring loop also have same trend 
#TM_M
TM_M_h <- TM_com[TM_com$M>1.5*TM_com$T,1:3]
TM_T_h <- TM_com[TM_com$T>1.5*TM_com$M,1:3]
TM_com <- as.data.frame(fread("./TM_com_loop.txt"))
colnames(TM_com) <- c("chr1_T","start1_T","end2_T","chr1_M","start1_M","end2_M")
colnames(TM_M_h) <- c("chr1_M","start1_M","end2_M")
test_M <- merge(TM_M_h,TM_com,c("chr1_M","start1_M","end2_M"))
colnames(loop)[1:3] <- c("chr1_T","start1_T","end2_T")
test_M <- merge(loop,test_M[,4:6],c("chr1_T","start1_T","end2_T"))
test_M <- test_M[test_M$M>1.5*test_M$T,]
colnames(TM_M_h) <- c("chr1_T","start1_T","end2_T")
test_M2 <- merge(TM_M_h,TM_com,c("chr1_T","start1_T","end2_T"))
colnames(loop)[1:3] <- c("chr1_T","start1_T","end2_T")
test_M2 <- merge(loop,test_M2[,1:3],c("chr1_T","start1_T","end2_T"))
test_M2 <- test_M[test_M2$M>1.5*test_M2$T,]
TM_M_h <- unique(inner_join(test_M,test_M2))
TM_M <- as.data.frame(rbind(as.matrix(TM_M),as.matrix(TM_M_h[,1:3])))
#TM_T
colnames(TM_T_h) <- c("chr1_T","start1_T","end2_T")
test_T <- merge(TM_T_h,TM_com,c("chr1_T","start1_T","end2_T"))
colnames(loop)[1:3] <- c("chr1_M","start1_M","end2_M")
test_T <- merge(loop,test_T[,4:6],c("chr1_M","start1_M","end2_M"))
test_T <- test_T[test_T$T>1.5*test_T$M,]
colnames(TM_T_h) <- c("chr1_M","start1_M","end2_M")
test_T2 <- merge(TM_T_h,TM_com,c("chr1_M","start1_M","end2_M"))
colnames(loop)[1:3] <- c("chr1_M","start1_M","end2_M")
test_T2 <- merge(loop,test_T2[,1:3],c("chr1_M","start1_M","end2_M"))
test_T2 <- test_T2[test_T2$T>1.5*test_T2$M,]
TM_T_h <- unique(inner_join(test_T,test_T2))
TM_T <- as.data.frame(rbind(as.matrix(TM_T),as.matrix(TM_T_h[,1:3])))
fwrite(TM_M,"./test/TM_M_loop.txt",sep="\t",quote=F,col.names=T)
fwrite(TM_T,"./test/TM_T_loop.txt",sep="\t",quote=F,col.names=T)
TM_M <- fread("./test/TM_M_loop.txt")
TM_T <- fread("./test/TM_T_loop.txt")
TM_all <- unique(rbind(loop_CT[,c(1,2,6)],loop_LT[,c(1,2,6)]))
TM_com <- anti_join(anti_join(TM_all,TM_T),TM_M)
TM_M$chr2 <- TM_M$chr1;TM_M$end1 <- as.numeric(TM_M$start1)+10000;TM_M$start2 <- as.numeric(TM_M$end2)-10000;
TM_T$chr2 <- TM_T$chr1;TM_T$end1 <- as.numeric(TM_T$start1)+10000;TM_T$start2 <- as.numeric(TM_T$end2)-10000;
TM_com$chr2 <- TM_com$chr1;TM_com$end1 <- as.numeric(TM_com$start1)+10000;TM_com$start2 <- as.numeric(TM_com$end2)-10000;
TM_M <- rbind(as.matrix(TM_M[,c("chr1","start1","end1")]),as.matrix(TM_M[,c("chr2","start2","end2")]))
TM_T <- rbind(as.matrix(TM_T[,c("chr1","start1","end1")]),as.matrix(TM_T[,c("chr2","start2","end2")]))
TM_com <- rbind(as.matrix(TM_com[,c("chr1","start1","end1")]),as.matrix(TM_com[,c("chr2","start2","end2")]))
fwrite(TM_M,"./test/TM_M_bou.bed",sep="\t",quote=F,col.names=F)
fwrite(TM_T,"./test/TM_T_bou.bed",sep="\t",quote=F,col.names=F)
fwrite(TM_com,"./test/TM_T_com_bou.bed",sep="\t",quote=F,col.names=F)

#######################################################
###plot heatmap for single loop
#######################################################
i <- 5655 ##select the loop to draw
loop <- as.data.frame(fread("./NTM_loop.txt"))
loop$V2 <- (loop$V2-4*10000)/10000
loop$V3 <- (loop$V3+5*10000)/10000
loop$V5 <- (loop$V5-4*10000)/10000
loop$V6 <- (loop$V6+5*10000)/10000
CN_mat <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/merged_validPairs/hic_results/data/all_CN/matrix/chr4_10Kb_CN.txt",))
CT_mat <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/merged_validPairs/hic_results/data/all_CT/matrix/chr4_10Kb_CT.txt",))
LT_mat <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/merged_validPairs/hic_results/data/all_LT/matrix/chr4_10Kb_LT.txt",))

loop_CN <- CN_mat[loop[i,2]:loop[i,3],loop[i,5]:loop[i,6]]
loop_CT <- CT_mat[loop[i,2]:loop[i,3],loop[i,5]:loop[i,6]]
loop_LT <- LT_mat[loop[i,2]:loop[i,3],loop[i,5]:loop[i,6]]

col <- colorRampPalette(c("white","lightblue"))
max_int <- max(max(loop_CN),max(loop_CT),max(loop_LT))
pheatmap(loop_CN,color = col(max_int[1:max(loop_CN)]),cluster_rows=F,cluster_cols=F,border_color=NA,show_colnames=F,show_rownames=F,
         legend=F)
pheatmap(loop_CT,color = col(max_int[1:max(loop_CT)]),cluster_rows=F,cluster_cols=F,border_color=NA,show_colnames=F,show_rownames=F,
         legend=F)
pheatmap(loop_LT,color = col(max_int[1:max(loop_LT)]),cluster_rows=F,cluster_cols=F,border_color=NA,show_colnames=F,show_rownames=F,
         legend=F)

score_CN <- loop_CN[6,6]/mean(unlist(loop_CN[7:9,7:9]))
print(paste0("score_CN is: ",score_CN))
score_CT <- loop_CT[6,6]/mean(unlist(loop_CT[7:9,7:9]))
print(paste0("score_CT is: ",score_CT))
score_LT <- loop_LT[6,6]/mean(unlist(loop_LT[7:9,7:9]))
print(paste0("score_LT is: ",score_LT))

#######################################################
###find continously altered loops
#######################################################
T_u <- as.data.frame(fread("./NT_T_loop.txt"))
M_u <- as.data.frame(fread("./TM_M_loop.txt"))
T_d <- as.data.frame(fread("./NT_N_loop.txt"))
M_d <- as.data.frame(fread("./TM_T_loop.txt"))
uu <- inner_join(T_u,M_u,c("V1","V2","V3","V4","V5","V6"))
dd <- inner_join(T_d,M_d,c("V1","V2","V3","V4","V5","V6"))
fwrite(uu,"uu_loop.txt",sep="\t",quote=F,col.names=F)
fwrite(dd,"dd_loop.txt",sep="\t",quote=F,col.names=F)
colnames(dd) <- colnames(uu) <- c("c","s","e","c","s","e")
fwrite(rbind(uu[,1:3],uu[,4:6]),"uu_bou.bed",sep="\t",quote=F,col.names=F)
fwrite(rbind(dd[,1:3],dd[,4:6]),"dd_bou.bed",sep="\t",quote=F,col.names=F)
###annotation for diff regions
library(GenomicFeatures)
library("ChIPseeker")
library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
uu <- readPeakFile("./test/TM_M_spe_bou_sorted.bed")
dd <- readPeakFile("./test/TM_T_spe_bou_sorted.bed")
peaks <- list("uu"=uu,"dd"=dd)
# you can define promotor region by yourself
promoter <- getPromoters(TxDb=txdb, upstream=1500, downstream=1500)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
#you can transfer geneID (Entrez，ENSEMBL，SYMBOL，GENENAME) after annotatePea be imported into annoDb
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000), 
                       verbose=FALSE,addFlankGeneInfo=F, flankDistance=5000,annoDb="org.Hs.eg.db")
gene = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
fwrite(as.data.frame(peakAnnoList[["uu"]]),"./test/TM_M_spe_bou_genes_chipseeker.txt",sep="\t",quote=F,col.names=T)
fwrite(as.data.frame(peakAnnoList[["dd"]]),"./test/TM_T_spe_bou_genes_chipseeker.txt",sep="\t",quote=F,col.names=T)

clu1 <- fread("./uu_genes.txt")
clu1_dis <- clu1[clu1$annotation=="Distal Intergenic",]
clu1_gene <- clu1[clu1$annotation!="Distal Intergenic",]
clu1_gene <- separate(clu1_gene,"annotation",c("annotation"),"\\(")
clu1_intron <- clu1_gene[clu1_gene$annotation=="Intron ",]
clu1_gene <- clu1_gene[clu1_gene$annotation!="Intron ",]
clu1_pro <- clu1_gene[clu1_gene$annotation=="Promoter",]
#fwrite(clu1_pro,"./uu_genes_pro.txt",sep="\t",quote=F)

clu1_pro <- clu1_pro[!duplicated(clu1_pro)]
#expr <- fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments/all_sample_gene_expr.txt")
expr <- fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/All_read_count_normalized.txt")
mygene <- match(clu1_pro$SYMBOL,expr$V1)
mygene <- expr[mygene,]
mygene <- mygene[is.na(mygene$V1)!=T,]
#fwrite(mygene,"uu_genes_expr.txt",sep="\t",quote=F,col.names=T)
CN <- mygene[,c(2:14)]
CN$x <- apply(CN,1,sum);CN$x <- CN$x/13
CT <- mygene[,c(15:27)]
CT$x <- apply(CT,1,sum);CT$x <- CT$x/13
wilcox.test(CT$x,CN$x,paired = T)
t.test(CT$x,CN$x,paired = T)
boxplot(list("N"=CN$x,"T"=CT$x))

CT <- mygene[,c(15:21,24,27)]
CT$x <- apply(CT,1,sum);CT$x <- CT$x/9
CN <- mygene[,c(2:8,10,14)]
CN$x <- apply(CN,1,sum); CN$x <- CN$x/9
wilcox.test(CT$x,CN$x,paired = T)
t.test(CT$x,CN$x,paired = T)
boxplot(list("N"=CN$x,"T"=CT$x))

CT <- mygene[,c(15:21,24,27)]
CT$x <- apply(CT,1,sum);CT$x <- CT$x/9
LT <- mygene[,c(29:34)]
LT$x <- apply(LT,1,sum); LT$x <- LT$x/6
wilcox.test(CT$x,LT$x,paired = T)
t.test(CT$x,LT$x,paired = T)
boxplot(list("N"=LT$x,"T"=CT$x))

#venn plot
library(VennDiagram)
library(grDevices)
T_gained <- paste0(T_u$V1,T_u$V2,T_u$V3,T_u$V4,T_u$V5,T_u$V6)
M_gained <- paste0(M_u$V1,M_u$V2,M_u$V3,M_u$V4,M_u$V5,M_u$V6)
venn.plot<-venn.diagram(list("T_gained"=T_gained,"M_gained"=M_gained),fill=c("orange","pink"),filename=NULL,
             col=NA,imagetype="tiff",alpha=0.5,margin=0.05)
pdf(file="uu_venn.pdf")
grid.draw(venn.plot)
dev.off()
T_lost <- paste0(T_d$V1,T_d$V2,T_d$V3,T_d$V4,T_d$V5,T_d$V6)
M_lost <- paste0(M_d$V1,M_d$V2,M_d$V3,M_d$V4,M_d$V5,M_d$V6)
venn.plot<-venn.diagram(list("T_lost"=T_lost,"M_lost"=M_lost),fill=c("orange","pink"),filename=NULL,
                        col=NA,imagetype="tiff",alpha=0.5,margin=0.05)
pdf(file="dd_venn.pdf")
grid.draw(venn.plot)
dev.off()

