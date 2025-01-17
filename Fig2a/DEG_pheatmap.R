library("pheatmap")
library("RColorBrewer")
library("ggplot2")

dir.create("Output_heatmap")

gene_count = read.csv("Input_DEG/combine_significant_normalized_count.csv",row.names = 1)

all_sample = colnames(gene_count)
Group = c()
for (i in all_sample){
  Group <-append(Group,unlist(strsplit(i,"[.]"))[3])
}
table(Group)
sample_group <- data.frame(row.names = all_sample, group = Group)

k = Group %in% c("N","T","M");table(k)
gene_count <- as.data.frame(gene_count)[,k]
gene_count <-subset(gene_count, select = c(1,4,7,9,11,13,15,17,19,22,2,5,8,10,12,14,16,18,20,23,3,6,21))

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

p = pheatmap(gene_count, show_rownames = F, cluster_cols = F, annotation_col = sample_group, scale="row",clustering_callback = callback, cutree_rows = 3, gaps_col = c(10,20))

ggsave(p, filename = "Output_heatmap/DEG_z-score_heatmap.pdf", width = 10, height = 10)


