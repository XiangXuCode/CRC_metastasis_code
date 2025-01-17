library("DESeq2")
library("ggplot2")

dir.create("Output_DEseq2_PCA")

gene_count = read.csv("Input_prepDE/gene_count_matrix.csv",row.names = 1)
all_sample = colnames(gene_count)
Group = c()
for (i in all_sample){
  Group <-append(Group,unlist(strsplit(i,"[.]"))[3])
}
table(Group)

batch = c()
for (i in all_sample){
  batch <-append(batch,unlist(strsplit(i,"[.]"))[2])
}
table(batch)
sample_group <- data.frame(row.names = all_sample, group = Group, batch = batch)

countMatrix <- as.matrix(gene_count)

dds <- DESeqDataSetFromMatrix(countMatrix,colData = sample_group, design = ~ batch + group)

dds <- DESeq(dds,parallel = T)
vsd <- vst(dds, blind=FALSE)

k = Group %in% c("N","T","M");table(k)

mat <- assay(vsd)
mm <- model.matrix(~group, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat

write.csv(as.data.frame(assay(vsd)), file = "Output_DEseq2_PCA/removeBatchEffect_normalized_count_all.csv")

p4 = plotPCA(vsd, intgroup=c("group"))
ggsave(p4, filename = "Output_DEseq2_PCA/removeBatchEffect_cluster_PCA.pdf")


