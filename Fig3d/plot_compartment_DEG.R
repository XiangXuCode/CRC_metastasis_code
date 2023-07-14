library(ggplot2)
library(ggpubr)

T_vs_N_sample = c("CRC-01")
M_vs_T_sample = c("CRC-02")

for(s in T_vs_N_sample){
  dd = read.csv(paste0("Output_compartment_gene_log2Foldchange/",s,"-HiC/T_vs_N_combine_log2Foldchange.csv"))
  
  p = ggboxplot(dd,x = "group",y = "log2Foldchange",color = "group",palette = c("#BACCD9","#15559A","#C21F30","#EEA2A4"),add = "mean",title = paste0(s,"-HiC"),width = 0.8,bxp.errorbar = TRUE)
  my_comparisons = list( c("common_B", "diff_A_to_B"), c("common_B", "diff_B_to_A"), c("common_B", "common_A") , c("diff_A_to_B", "diff_B_to_A") , c("diff_A_to_B", "common_A") , c("diff_B_to_A", "common_A") )
  p2 = p + stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif",
                              method = "t.test")
  
  ggsave(paste0("Output_compartment_gene_log2Foldchange/",s,"-HiC/T_vs_N_combine_log2Foldchange.pdf"),p2,width=5,height = 7)
}

for(s in M_vs_T_sample){
  dd = read.csv(paste0("Output_compartment_gene_log2Foldchange/",s,"-HiC/M_vs_T_combine_log2Foldchange.csv"))
  
  p = ggboxplot(dd,x = "group",y = "log2Foldchange",color = "group",palette = c("#BACCD9","#15559A","#C21F30","#EEA2A4"),add = "mean",title = paste0(s,"-HiC"),width = 0.8,bxp.errorbar = TRUE)
  my_comparisons = list(c("common_B", "diff_A_to_B"), c("common_B", "diff_B_to_A"), c("common_B", "common_A") , c("diff_A_to_B", "diff_B_to_A") , c("diff_A_to_B", "common_A") , c("diff_B_to_A", "common_A") )
  p2 = p + stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif",
                              method = "t.test")
  
  ggsave(paste0("Output_compartment_gene_log2Foldchange/",s,"-HiC/M_vs_T_combine_log2Foldchange.pdf"),p2,width=5,height = 7)
}
