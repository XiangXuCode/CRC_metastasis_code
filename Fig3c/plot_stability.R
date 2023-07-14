library(ggplot2)
library(ggpubr)

dd = read.csv(file = "Output_compartment_stability/compartment_stability.txt", sep = "\t", header = T)

p = ggboxplot(dd,x = "group",y = "SOD",color = "group",palette = c("#BACCD9","#15559A","#C21F30","#EEA2A4"),title = "Compartment_stability",width = 0.6,bxp.errorbar = TRUE)+geom_jitter(width = 0.2, col = "black")

my_comparisons = list( c("without_metastasis_T_vs_N", "metastasis_T_vs_N"), c("without_metastasis_T_vs_N", "metastasis_M_vs_T"), c("without_metastasis_T_vs_N", "metastasis_M_vs_N") , c("metastasis_T_vs_N", "metastasis_M_vs_T") , c("metastasis_T_vs_N", "metastasis_M_vs_N") , c("metastasis_M_vs_T", "metastasis_M_vs_N") )

p2 = p + stat_compare_means(comparisons = my_comparisons,
                            label = "p.signif",
                            method = "wilcox.test")

ggsave(paste0("Output_compartment_stability/compartment_stability.pdf"),p2,width=20,height = 7)
