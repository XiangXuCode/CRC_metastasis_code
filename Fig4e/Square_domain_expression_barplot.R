library(ggplot2)
library(ggpubr)
library(rstatix)

dirs <- list.files(path = "151-combine_process_for_boxplot2")

dd = read.csv(paste0("Input_square_domain_expression/process.csv"))

dir.create("Output_plot_square_domain_expression")

p = ggboxplot(dd,"sample","normalized_count",col= "group",palette = "lancet",add = "mean",,bxp.errorbar = TRUE,shape="group",outlier.shape = NA,add.params =list(size = 0.4, jitter = 0.4),ggtheme = theme_bw(),xlab = "sample",ylab = "normalized count") + rotate_x_text(45)

stat_t_test <- rstatix::t_test(group_by(dd,sample),normalized_count ~ group)
stat_t_test <- rstatix::add_significance(stat_t_test, 'p') 
stat_t_test.test <-  rstatix::add_xy_position(stat_t_test, x = 'sample')
p1 = p + stat_pvalue_manual(stat_t_test.test, 
                        label = 'p.signif',
                        tip.length = 0.005,
                        hide.ns = FALSE)

ggsave(paste0("Output_plot_square_domain_expression/boxplot.pdf"),p1,height=10,width=20)


p = ggboxplot(dd,"sample","normalized_count",col= "group",palette = "lancet",add = "mean",,bxp.errorbar = TRUE,shape="group",ylim = c(3,17),outlier.shape = NA,add.params =list(size = 0.4, jitter = 0.4),ggtheme = theme_bw(),xlab = "sample",ylab = "normalized count") + rotate_x_text(45)

stat_t_test <- rstatix::t_test(group_by(dd,sample),normalized_count ~ group)
stat_t_test <- rstatix::add_significance(stat_t_test, 'p') 
stat_t_test.test <-  rstatix::add_xy_position(stat_t_test, x = 'sample')
p1 = p + stat_pvalue_manual(stat_t_test.test, 
                        label = 'p.signif',
                        tip.length = 0.005,
                        hide.ns = FALSE)
ggsave(paste0("Output_plot_square_domain_expression/boxplot_large.pdf"),p1,height=10,width=20)

