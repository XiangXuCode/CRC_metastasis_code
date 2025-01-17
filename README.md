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

Fig4:
4a_hicPlotTADs.py is used to plot TAD tracks. Need to install python module, HiCExplorer and required modules. 

Fig5:
Call_loop.sh is used to call loop. Need to install python module, Juicer and required modules. 
calculate_APA.sh is used to loop APA.
plot&analysis_loop.R is used to analyze and plot loops. Need to install R module, data.table, dplyr, VennDiagram, pheatmap, grDevices and required modules.
multiple_genes_boxplot.R is used to plot multiple genes. Need to install R module, cowplot, tidyverse, ggsci and required modules.
