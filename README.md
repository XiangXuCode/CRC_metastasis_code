# CRC_metastasis_code
Multiomics data analysis code for colorectal cancer.

The code follows the order of the figures:

Fig1:
1c_DEseq2.R is used to apply principal component analysis (PCA) to analyze RNA-seq data. Need to install R module, DEseq2 and required modules.
1d_plot_compartment_track.py is used to plot A/B compartment track. Need to install python module, pyGenomeTracks and required modules.

Fig2:
2a_DEG_pheatmap.R is used to plot differentially expression gene heatmap. Need to install R module, pheatmap, RColorBrewer, ggplot2 and required modules.
2b_Significant_DEG_venn3.py is used to plot venn of differentially expression gene. Need to install python module, numpy, matplotlib and required modules.

Fig3:
3a_compartment_pheatmap.R is used to plot A/B compartment heatmap. Need to install R module, pheatmap, RColorBrewer, ggplot2 and required modules.
3b_compartment_change.py is used to analyze A/B compartment switch ratio. Need to install python module, numpy, matplotlib and required modules.
3c_compartment_stability.py and plot_stability.R are used to analyze and plot A/B compartment stability. Need to install R module, ggplot2, ggpubr and required modules.
3d_compartment_DEG_log2Foldchange.py and plot_compartment_DEG.R are used to analyze and plot the expression in A/B compartment switch region.
3e_Scatter_plot_contour_plot_correlation_signif_DEG.R is used to plot the relation between PCA and gene expression. Need to install R module,ggplot2,MASS,dplyr,viridis,ggpointdensity,ggpubr and required modules.

Fig4:
4a_hicPlotTADs.py is used to plot TAD tracks. Need to install python module, HiCExplorer and required modules. 
4b_combine_log2Foldchange.py and plot_gene_log2Foldchange_TH_score.R are used to analyze and plot gene expression in different TH score type.
4e_Square_domain_recognition.py and Square_domain_expression_barplot.R are used to recognize square domain and analyze square domain anchor gene expression. Need to install R module, ggplot2, ggpubr, rstatix and required modules.

Fig5:
Call_loop.sh is used to call loop. Need to install python module, Juicer and required modules. 
calculate_APA.sh is used to loop APA.
plot&analysis_loop.R is used to analyze and plot loops. Need to install R module, data.table, dplyr, VennDiagram, pheatmap, grDevices and required modules.
multiple_genes_boxplot.R is used to plot multiple genes. Need to install R module, cowplot, tidyverse, ggsci and required modules.

Fig6:
b_1_circos.for_F6a.py is used to map the circle of chromosomal translocations. Need to install python module, pyCircos (https://github.com/ponnhide/pyCircos) and pyCircos required modules. 
b_2_bar_T.M.count.for_F6b.py is used to draw a bar diagram describing the statistics of chromatin translocation regions in multiple samples. Common python modules such as pandas and matplotlib seaborn numpy need to be installed. 
b_3_heatmap.for_F6c.py is used to map structural variations and nearby key genes. Need to install python module, NeoLoopFinder NeoLoopFinder (https://github.com/XiaoTaoWang/) and NeoLoopFinder required modules.
