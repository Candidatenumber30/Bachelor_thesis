library(ggpubr)
library (tidyverse)

# Read Clinical data and LCPM expression data
clin_complete <- qs_read("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCAN_B_Final_Clinical_Data.qs2")
lcpm <- qs_read("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCANB_Final_LCPM.qs2")
lcpm <- as.data.frame(lcpm)

# MAke dataframe of SRGN expression LCPM values
SRGN <- lcpm["SRGN",, drop = FALSE]
SRGN <- t(SRGN)
SRGN <- as.data.frame(SRGN)
SRGN$sample_id <- rownames(SRGN)

# Merge SRGN expression data with clinical data
merged <- merge(clin_complete, SRGN, by.x = "GEX.assay", "sample_id")

# Order Subtypes
merged$SSP.PAM50 <- factor(merged$SSP.PAM50, levels = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like", "Normal-like"))

# Calculate samples in each subtype
counts <- merged %>%
  group_by(SSP.PAM50) %>%
  summarise(n = n())

labels_with_n <- setNames(
  paste0(counts$SSP.PAM50, "\n(n=", counts$n, ")"),
  counts$SSP.PAM50)

# Violin plot
p_violin_PAM50 <- ggviolin(merged, 
  x = "SSP.PAM50", y = "SRGN", add = c("boxplot"), fill = "SSP.PAM50", 
  xlab = "PAM50 subtype",
  axis.text.x = element_text(size=1))+
  stat_compare_means(label = "p.signif", size = 8)+
  labs(fill = "SSP.PAM50", y = expression(Log(1 + italic(SRGN)))) + geom_jitter(size=0.2) +
  scale_fill_manual(values = pal_jco()(5)) + 
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none")+ 
    scale_x_discrete(labels = labels_with_n)
p_violin_PAM50

ggsave("figure_Violin_PAM50.pdf", p_violin_PAM50, width = 9, height = 8)
