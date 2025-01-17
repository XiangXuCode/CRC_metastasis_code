library("pheatmap")
library("RColorBrewer")
library("ggplot2")

dir.create("Output_compartment_heatmap")
gene_count = read.csv("Input_compartment/combine_compartment_rm_NA.csv",row.names = 1)

all_sample = colnames(gene_count)
Group = c()
for (i in all_sample){
  Group <-append(Group,unlist(strsplit(i,"[.]"))[3])
}
table(Group)
sample_group <- data.frame(row.names = all_sample, group = Group)

k = Group %in% c("N","T","M");table(k)
gene_count <- as.data.frame(gene_count)[,k]
gene_count <-subset(gene_count, select = c(4,7,10,13,20,1,15,17,5,8,11,14,21,2,16,18,3,6,9,12,19))

gene_count <- t(gene_count)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,7]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"),alpha=0)(paletteLength)

myBreaks <- c(seq(min(gene_count), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(gene_count)/paletteLength, max(gene_count), length.out=floor(paletteLength/2)))

p = pheatmap(gene_count, show_rownames = T, show_colnames = F, cluster_rows = F,cluster_cols = T, gaps_row = c(8,16),cutree_cols = 4,color=myColor, breaks=myBreaks, clustering_callback = callback, fontsize = 15)
ggsave(p, filename = "Output_compartment_heatmap/compartment_heatmap.pdf", width = 15, height = 10)

