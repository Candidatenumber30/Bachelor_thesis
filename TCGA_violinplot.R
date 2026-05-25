library(ggplot2)
library(ggpubr)
library(ggsci)

# Load merged data file of samples and SRGN expression with cancer types ordered by median

load("TCGA_Pancancer_files.RData")

# Violin plot
lancet <- pal_lancet()(9)
custom_lancet <- colorRampPalette(lancet)
p_TCGA <- ggviolin(merged_ordered, x = "cancer_type", y = "SRGN", 
                   add = c("boxplot"), fill = "cancer_type", 
                   orientation = c("horizontal"))+ 
                   labs(
                    x = "Cancer type", 
                    y = expression("Log (1+"*italic(SRGN)*")"))+
                    theme(
                    axis.text.x = element_text(size=16),
                    axis.text.y = element_text(size = 16),
                    axis.title.x = element_text(size = 20),
                    axis.title.y = element_text(size = 20),
                    legend.position = "none")+ 
                    stat_compare_means(label = "p.signif", size = 10) +
                    labs(fill = "Cancer type") + geom_jitter(size=0.5) +
                    scale_fill_manual(values = custom_lancet(32))
p_TCGA

ggsave("TCGA_Violinplot.pdf", p_TCGA, height = 14, width = 13)
