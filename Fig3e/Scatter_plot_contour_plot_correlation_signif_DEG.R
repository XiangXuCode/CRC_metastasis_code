library(ggplot2)
library(MASS)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(ggpubr)

dir.create("Output_scatter_plot_contour_plot_correlation_signif_DEG")

T_vs_N_sample = c("CRC-01")
M_vs_T_sample = c("CRC-02")

for(s in T_vs_N_sample){

  dir.create(paste0("Output_scatter_plot_contour_plot_correlation_signif_DEG/",s,"-HiC"))
  df1 = read.csv(paste0("input_intersect_compartment_score_signif_DEG/",s,"-HiC/T_vs_N"),sep = "\t",row.names = NULL, header = FALSE)
  df1_PCA1_change = df1['V4']
  colnames(df1_PCA1_change) = "Compartment change"

  df1_log2Foldchange = df1['V9']
  colnames(df1_log2Foldchange) = "Expression log2Foldchange"

  dfdata =df1[c('V4','V9')]
  colnames(dfdata)=c('df1_PCA1_change','df1_log2Foldchange')

  p1 <- ggplot(data = dfdata, mapping = aes(x = df1_PCA1_change,
                                         y = df1_log2Foldchange)) + 
   geom_pointdensity() + 
    scale_color_viridis() + 
    geom_smooth(method = lm) +  
    stat_cor(method = "spearman") + 
    xlab("Compartment change") + 
    theme(axis.title.x = element_text(size = 16,
                                      face = "bold", 
                                     vjust = 0.5, 
                                      hjust = 0.5))+
    ylab("Expression log2Foldchange") + 
    theme(axis.title.y = element_text(size = 16,
                                      face = "bold", 
                                      vjust = 0.5, 
                                     hjust = 0.5))+
    theme_bw()+
    theme(panel.grid.major=element_line(colour=NA),
         panel.background = element_rect(fill = "transparent",colour = NA),
         plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank(),
         text=element_text(size=12,  family="serif")) +
         theme(legend.position='none')  
  p2 <- p1 + labs(title = paste0(s,"-HiC"))
  ggsave(paste0("Output_scatter_plot_contour_plot_correlation_signif_DEG/",s,"-HiC/T_vs_N_Scatter_contour.pdf"),p2,width = 8,height = 4)
}

for(s in M_vs_T_sample){
  
  dir.create(paste0("Output_scatter_plot_contour_plot_correlation_signif_DEG/",s,"-HiC"))
  df1 = read.csv(paste0("input_intersect_compartment_score_signif_DEG/",s,"-HiC/M_vs_T"),sep = "\t",row.names = NULL, header = FALSE)
  df1_PCA1_change = df1['V4']
  colnames(df1_PCA1_change) = "Compartment change"
  
  df1_log2Foldchange = df1['V9']
  colnames(df1_log2Foldchange) = "Expression log2Foldchange"
  
  dfdata =df1[c('V4','V9')]
  colnames(dfdata)=c('df1_PCA1_change','df1_log2Foldchange')
  
  p1 <- ggplot(data = dfdata, mapping = aes(x = df1_PCA1_change,
                                            y = df1_log2Foldchange)) + 
    geom_pointdensity() + 
    scale_color_viridis() + 
    geom_smooth(method = lm) + 
    stat_cor(method = "spearman") + 
    xlab("Compartment change") + 
    theme(axis.title.x = element_text(size = 16,
                                      face = "bold", 
                                      vjust = 0.5, 
                                      hjust = 0.5))+
    ylab("Expression log2Foldchange") + 
    theme(axis.title.y = element_text(size = 16,
                                      face = "bold", 
                                      vjust = 0.5, 
                                      hjust = 0.5))+
    theme_bw()+
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank(),
          text=element_text(size=12,  family="serif")) +
    theme(legend.position='none') 
  p2 <- p1 + labs(title = paste0(s,"-HiC"))
  ggsave(paste0("Output_scatter_plot_contour_plot_correlation_signif_DEG/",s,"-HiC/M_vs_T_Scatter_contour.pdf"),p2,width = 8,height = 4)
}
