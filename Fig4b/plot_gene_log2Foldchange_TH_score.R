library(ggplot2)
library(ggpubr)

T_vs_N_sample = c("CRC-01")
M_vs_T_sample = c("CRC-02")

for(s in T_vs_N_sample){
  dd = read.csv(paste0("Output_gene_log2Foldchange_TH_score/",s,"-HiC/T_vs_N_combine_log2Foldchange.csv"))
  
  p = ggboxplot(dd,x = "group",y = "log2Foldchange",color = "group",palette = c("grey","#C21F30","#15559A"), add = c("mean"),title = paste0(s,"-HiC"),width = 0.8,bxp.errorbar = TRUE,outlier.shape = NA,ylim = c(-0.05,0.08))
  my_comparisons = list(c("stable","raise"),c("stable","reduce"),c("raise","reduce"))

  p2 = p + stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif",
                              method = "t.test")
  ggsave(paste0("Output_gene_log2Foldchange_TH_score/",s,"-HiC/T_vs_N_combine_log2Foldchange_large.pdf"),p2,width=4,height = 7)


  p = ggboxplot(dd,x = "group",y = "log2Foldchange",color = "group",palette = c("grey","#C21F30","#15559A"), add = c("mean"),title = paste0(s,"-HiC"),width = 0.8,bxp.errorbar = TRUE,outlier.shape = NA)
  my_comparisons = list(c("stable","raise"),c("stable","reduce"),c("raise","reduce"))
  p2 = p + stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif",
                              method = "t.test")
  ggsave(paste0("Output_gene_log2Foldchange_TH_score/",s,"-HiC/T_vs_N_combine_log2Foldchange.pdf"),p2,width=4,height = 7)
}


for(s in M_vs_T_sample){
  dd = read.csv(paste0("Output_gene_log2Foldchange_TH_score/",s,"-HiC/M_vs_T_combine_log2Foldchange.csv"))
  
  p = ggboxplot(dd,x = "group",y = "log2Foldchange",color = "group",palette = c("grey","#C21F30","#15559A"),add = "mean",title = paste0(s,"-HiC"),width = 0.8,bxp.errorbar = TRUE,outlier.shape=NA,ylim = c(-0.05,0.08))
  my_comparisons = list(c("stable","raise"),c("stable","reduce"),c("raise","reduce"))
  p2 = p + stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif",
                              method = "t.test")
  ggsave(paste0("Output_gene_log2Foldchange_TH_score/",s,"-HiC/M_vs_T_combine_log2Foldchange_large.pdf"),p2,width=4,height = 7)

  p = ggboxplot(dd,x = "group",y = "log2Foldchange",color = "group",palette = c("grey","#C21F30","#15559A"),add = "mean",title = paste0(s,"-HiC"),width = 0.8,bxp.errorbar = TRUE,outlier.shape=NA)
  my_comparisons = list(c("stable","raise"),c("stable","reduce"),c("raise","reduce"))
  p2 = p + stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif",
                              method = "t.test")
  ggsave(paste0("Output_gene_log2Foldchange_TH_score/",s,"-HiC/M_vs_T_combine_log2Foldchange.pdf"),p2,width=4,height = 7)
}
 
