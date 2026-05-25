library(qs2)
library(edgeR)
library(data.table)
library(dplyr)
library(tidyverse)

# Read clinical data, CPM and LCPM expression data
lcpm <- qs_read("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCANB_Final_LCPM.qs2")
cpm <- qs_read("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCANB_Final_CPM.qs2")
lcpm <- qs_read("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCANB_Final_LCPM.qs2")
cpm <- as.data.frame(cpm)

# take out CPMs as chunks of 500 patients within each subtype, export as tab delimited txt

# Luminal A
Luma_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Luminal A") %>% pull(GEX.assay)
Luma_cpm <- cpm %>% dplyr::select(all_of(Luma_pt))
chunk_size <- 500
num_patients <- ncol(Luma_cpm)-1
starts <- seq(1, num_patients, by = chunk_size)

for (i in seq_along(starts)) {
  start_col <- starts[i]
  end_col <- min(start_col + chunk_size - 1, num_patients)
  chunk_data <- Luma_cpm[,start_col:end_col, drop = FALSE]
  chunk_data <- chunk_data %>% rownames_to_column(var = "GeneID")
  filename <- paste0("Luma_cpm_chunk_", i, ".txt")
  write.table(chunk_data, filename, sep = "\t", row.names = FALSE)
  cat(paste("Exported:", filename, "\n"))
}
df_Luma3 <- data.table::fread("Luma_cpm_chunk_3.txt")

# Luminal B

sum(clin_complete$SSP.PAM50 == "Luminal A")

Lumb_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Luminal B") %>% pull(GEX.assay)
Lumb_cpm <- cpm %>% dplyr::select(all_of(Lumb_pt))
chunk_size <- 500
num_patients_Lumb <- ncol(Lumb_cpm)-1
starts_Lumb <- seq(2, num_patients_Lumb, by = chunk_size)

for (i in seq_along(starts_Lumb)) {
  start_col_Lumb <- starts_Lumb[i]
  end_col_Lumb <- min(start_col_Lumb + chunk_size - 1, num_patients_Lumb)
  chunk_data <- Lumb_cpm[, start_col_Lumb:end_col_Lumb, drop = FALSE]
  chunk_data <- chunk_data %>% rownames_to_column(var = "GeneID")
  filename <- paste0("Lumb_cpm_chunk_", i, ".txt")
  write.table(chunk_data, filename, sep = "\t", row.names = FALSE)
  cat(paste("Exported:", filename, "\n"))
}


# HER2-enriched

HER2_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "HER2-enriched") %>% pull(GEX.assay)
HER2_cpm <- cpm %>% dplyr::select(all_of(HER2_pt))
#HER2_cpm <- HER2_cpm %>% rownames_to_column(var = "GeneID")
chunk_size <- 500
num_patients_HER2 <- ncol(HER2_cpm)-1
starts_HER2 <- seq(1, num_patients_HER2, by = chunk_size)
for (i in seq_along(starts_HER2)) {
  start_col <- starts_HER2[i]
  end_col <- min(start_col + chunk_size - 1, num_patients_HER2)
  chunk_data <- HER2_cpm[, start_col:end_col, drop = FALSE]
  chunk_data <- chunk_data %>% rownames_to_column(var = "GeneID")
  filename <- paste0("HER2_cpm_chunk_", i, ".txt")
  write.table(chunk_data, filename, sep = "\t", row.names = FALSE)
  cat(paste("Exported:", filename, "\n"))
}

# Basal-like

Basal_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Basal-like") %>% pull(GEX.assay)
Basal_cpm <- cpm %>% dplyr::select(all_of(Basal_pt))
chunk_size <- 500
num_patients_Basal <- ncol(Basal_cpm)-1
starts_Basal <- seq(1, num_patients_Basal, by = chunk_size)
for (i in seq_along(starts_Basal)) {
  start_col <- starts_Basal[i]
  end_col <- min(start_col + chunk_size - 1, num_patients_Basal)
  chunk_data <- Basal_cpm[, start_col:end_col, drop = FALSE]
  chunk_data <- chunk_data %>% rownames_to_column(var = "GeneID")
  filename <- paste0("Basal_cpm_chunk_", i, ".txt")
  write.table(chunk_data, filename, sep = "\t", row.names = FALSE)
  cat(paste("Exported:", filename, "\n"))
}

# Normal-like

Normal_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Normal-like") %>% pull(GEX.assay)
Normal_cpm <- cpm %>% dplyr::select(all_of(Normal_pt))
chunk_size <- 500
num_patients_Normal <- ncol(Normal_cpm)-1
starts_Normal <- seq(1, num_patients_Normal, by = chunk_size)
for (i in seq_along(starts_Normal)) {
  start_col <- starts_Normal[i]
  end_col <- min(start_col + chunk_size - 1, num_patients_Normal)
  chunk_data <- Normal_cpm[, start_col:end_col, drop = FALSE]
  chunk_data <- chunk_data %>% rownames_to_column(var = "GeneID")
  filename <- paste0("Normal_cpm_chunk_", i, ".txt")
  write.table(chunk_data, filename, sep = "\t", row.names = FALSE)
  cat(paste("Exported:", filename, "\n"))
}

files_all <- c("CIBERSORTx_Job56_Results_Lumb_1.csv", "CIBERSORTx_Job57_Results_Lumb_2.csv", "CIBERSORTx_Job58_Results_Lumb_3.csv", 
               "CIBERSORTx_Job6_Results_LumA1.csv", "CIBERSORTx_Job7_Results_LumA2.csv", "CIBERSORTx_Job14_Results_LumA3.csv",
               "CIBERSORTx_Job18_Results_LumA4.csv", "CIBERSORTx_Job10_Results_LumA5.csv", )


# files from cibersortx, all subtypes

files_all <- c("CIBERSORTx_Job56_Results_Lumb_1.csv", "CIBERSORTx_Job57_Results_Lumb_2.csv", "CIBERSORTx_Job58_Results_Lumb_3.csv", 
               "CIBERSORTx_Job6_Results_LumA1.csv", "CIBERSORTx_Job7_Results_LumA2.csv", "CIBERSORTx_Job14_Results_LumA3.csv",
               "CIBERSORTx_Job18_Results_LumA4.csv", "CIBERSORTx_Job10_Results_LumA5.csv", "CIBERSORTx_Job43_Results_HER2_1.csv", 
               "CIBERSORTx_Job44_Results_HER2_2.csv", "CIBERSORTx_Job45_Results_Basal_1.csv", "CIBERSORTx_Job46_Results_Normal_1.csv",
               "CIBERSORTx_Job47_Results_Normal_2.csv")

cibersort_res_all <- files_all %>%
  lapply(readr::read_csv) %>%
  bind_rows()

# merge CIBERSORTx results with clinical data and LCPM SRGN values

all_clinical_merged <- merge(clin_complete, cibersort_res_all, by.x = "GEX.assay", by.y = "Mixture")

SRGN <- c("SRGN")
SRGN_df <- as.data.frame(lcpm)[intersect(rownames(lcpm), SRGN),]
SRGN_df <- as.data.frame(t(SRGN_df))
SRGN_df$Sample_ID <- rownames(SRGN_df)

all_clinical_merged_SRGN <- merge(SRGN_df, all_clinical_merged, by.y = "GEX.assay", by.x = "Sample_ID")

# Correlation plots between SRGN expression and estimated immune cell abundance
library(ggpubr)
library(ggplot2)
names(all_clinical_merged_SRGN)
cell_types <- c(
  "B cells naive", "B cells memory", "Plasma cells",
  "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta", "NK cells resting",
  "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1",
  "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils",
  "Neutrophils"
)

all(cell_types %in% colnames(all_clinical_merged_SRGN))
genes <- c("SRGN")

# clean labels
gene_labels <- c(
  SRGN = "Log(1+SRGN)"
)

SRGN <- all_clinical_merged_SRGN$SRGN
SRGN <- as.numeric(SRGN)

pdf("SRGN_Immune_cell_Correlation_all.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types) {
    
    df <- data.frame(
      x = all_clinical_merged_SRGN[[g]],
      y = all_clinical_merged_SRGN[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "spearman",   
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#4DBBD5") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#E64B35", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()

# files from cibersortx, Basal-like merge with SRGN expression and clinical data

cibersort_res_Basal <- read.csv("CIBERSORTx_Job45_Results_Basal_1.csv")
names(cibersort_res_Basal) <- gsub("\\.", " ", names(cibersort_res_Basal))
Basal_clinical <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Basal-like")
Basal_clinical_merged <- merge(Basal_clinical, cibersort_res_Basal, by.x = "GEX.assay", by.y = "Mixture")
Basal_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Basal-like") %>% pull(GEX.assay)
Basal_gex <- gex %>% dplyr::select(all_of(Basal_pt))
SRGN <- c("SRGN")
df_Basal <- as.data.frame(Basal_gex)[intersect(rownames(Basal_gex), SRGN),]
df_Basal <- as.data.frame(t(df_Basal))
df_Basal$Sample_ID <- rownames(df_Basal)
Basal_clinical_merged <- merge(df_Basal, Basal_clinical_merged, by.y = "GEX.assay", by.x = "Sample_ID")
Basal_clinical_merged <- Basal_clinical_merged %>%
  mutate(SRGN_log2 = log2(1 + SRGN)
  )

# Correlation SRGN with Immune cell populations
names(Basal_clinical_merged)
cell_types <- c(
  "B cells naive", "B cells memory", "Plasma cells",
  "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta", "NK cells resting",
  "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1",
  "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils",
  "Neutrophils"
)

cell_types_Basal <- c(
  "B cells naive", "B cells memory", "Plasma cells",
  "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory  Tregs ", "T cells gamma delta", "NK cells resting",
  "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1",
  "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils",
  "Neutrophils"
)

genes <- c("SRGN_log2")

# clean labels
gene_labels <- c(
  SRGN_log2 = "Log(1+SRGN)"
)

cell_types
cell_types_Basal %in% names(Basal_clinical_merged)

pdf("SRGN_Immune_cell_Correlation_Basal.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types_Basal) {
    
    df <- data.frame(
      x = Basal_clinical_merged[[g]],
      y = Basal_clinical_merged[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "spearman",   # Try with pearson as well
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#4DBBD5") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#E64B35", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()

# files from cibersortx, luma

Luma_1 <- read.csv("CIBERSORTx_Job6_Results_LumA1.csv")
Luma_2 <- read.csv("CIBERSORTx_Job7_Results_LumA2.csv")
Luma_3 <- read.csv("CIBERSORTx_Job10_Results_LumA5.csv")
Luma_4 <- read.csv("CIBERSORTx_Job14_Results_LumA3.csv")
Luma_5 <- read.csv("CIBERSORTx_Job18_Results_LumA4.csv")


files_Luma <- list.files(
  pattern = "CIBERSORTx_Job.*LumA.*\\.csv$",
  full.names = TRUE
)

cibersort_res_LumA <- files_Luma %>%
  lapply(readr::read_csv) %>%
  bind_rows()

luma_clinical <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Luminal A")
luma_clinical_merged <- merge(luma_clinical, cibersort_res_LumA, by.x = "GEX.assay", by.y = "Mixture")
luma_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Luminal A") %>% pull(GEX.assay)
luma_gex <- gex %>% dplyr::select(all_of(luma_pt))

SRGN <- c("SRGN")
Luma_df <- as.data.frame(luma_gex)[intersect(rownames(luma_gex), SRGN),]

Luma_df <- as.data.frame(t(Luma_df))
Luma_df$Sample_ID <- rownames(Luma_df)

luma_clinical_merged <- merge(Luma_df, luma_clinical_merged, by.y = "GEX.assay", by.x = "Sample_ID")

luma_clinical_merged <- luma_clinical_merged %>%
  mutate(SRGN_log2 = log2(1 + SRGN)
  )

library(ggpubr)
library(patchwork)

p1 <- ggdensity(luma_clinical_merged, x = "SRGN_log2") + 
  geom_vline(xintercept = median(luma_clinical_merged$SRGN_log2), linetype = "dashed", color = "red") +
  labs(title = "SRGN Log Values") +
  theme(legend.position = "none")
p1


##3 Doing correlations

library(ggpubr)
library(ggplot2)
names(luma_clinical_merged)
cell_types <- c(
  "B cells naive", "B cells memory", "Plasma cells",
  "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta", "NK cells resting",
  "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1",
  "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils",
  "Neutrophils"
)

genes <- c("SRGN_log2")

# clean labels
gene_labels <- c(
  SRGN_log2 = "SRGN"
)

#dev.off()
pdf("SRGN_Immune_Correlation_Luma.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types) {
    
    df <- data.frame(
      x = luma_clinical_merged[[g]],
      y = luma_clinical_merged[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "spearman",   # Try with pearson as well
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#1F78B4") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#FF0000", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()


# files from cibersortx, lumb

Lumb_1 <- read.csv("CIBERSORTx_Job56_Results_Lumb_1.csv")
Lumb_2 <- read.csv("CIBERSORTx_Job57_Results_Lumb_2.csv")
Lumb_3 <- read.csv("CIBERSORTx_Job58_Results_Lumb_3.csv")

files_Lumb <- c("CIBERSORTx_Job56_Results_Lumb_1.csv", "CIBERSORTx_Job57_Results_Lumb_2.csv", "CIBERSORTx_Job58_Results_Lumb_3.csv")

cibersort_res_Lumb <- files_Lumb %>%
  lapply(readr::read_csv) %>%
  bind_rows()
head(Lumb_clinical$GEX.assay)
class(cibersort_res_Lumb$Mixture)
head(cibersort_res_Lumb$Mixture)
intersect(Lumb_clinical$GEX.assay, cibersort_res_Lumb$Mixture)
setdiff(Lumb_clinical$GEX.assay, cibersort_res_Lumb$Mixture)
setdiff(cibersort_res_Lumb$Mixture, Lumb_clinical$GEX.assay)
any(Lumb_clinical$GEX.assay %in% cibersort_res_Lumb$Mixture)
Lumb_clinical <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Luminal B")
Lumb_clinical_merged <- merge(Lumb_clinical, cibersort_res_Lumb, by.x = "GEX.assay", by.y = "Mixture")
Lumb_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Luminal B") %>% pull(GEX.assay)
Lumb_gex <- gex %>% dplyr::select(all_of(Lumb_pt))
SRGN <- c("SRGN")
df_Lumb <- as.data.frame(Lumb_gex)[intersect(rownames(Lumb_gex), SRGN),]
df_Lumb <- as.data.frame(t(df_Lumb))
df_Lumb$Sample_ID <- rownames(df_Lumb)
Lumb_clinical_merged <- merge(df_Lumb, Lumb_clinical_merged, by.y = "GEX.assay", by.x = "Sample_ID")

Lumb_clinical_merged <- Lumb_clinical_merged %>%
  mutate(SRGN_log2 = log2(1 + SRGN)
  )

pdf("SRGN_Immune_Correlation_Lumb.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types) {
    
    df <- data.frame(
      x = Lumb_clinical_merged[[g]],
      y = Lumb_clinical_merged[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "spearman",   # Try with pearson as well
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#1F78B4") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#FF0000", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()

pdf("SRGN_Immune_Correlation_Lumb_pearson.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types) {
    
    df <- data.frame(
      x = Lumb_clinical_merged[[g]],
      y = Lumb_clinical_merged[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "pearson",   # Try with pearson as well
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#1F78B4") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#FF0000", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()

# files from cibersortx, HER2

files_all <- c("CIBERSORTx_Job56_Results_Lumb_1.csv", "CIBERSORTx_Job57_Results_Lumb_2.csv", "CIBERSORTx_Job58_Results_Lumb_3.csv", 
               "CIBERSORTx_Job6_Results_LumA1.csv", "CIBERSORTx_Job7_Results_LumA2.csv", "CIBERSORTx_Job14_Results_LumA3.csv",
               "CIBERSORTx_Job18_Results_LumA4.csv", "CIBERSORTx_Job10_Results_LumA5.csv", "CIBERSORTx_Job43_Results_HER2_1.csv", 
               "CIBERSORTx_Job44_Results_HER2_2.csv", )


HER2_1 <- read.csv("CIBERSORTx_Job43_Results_HER2_1.csv")
HER2_2 <- read.csv("CIBERSORTx_Job44_Results_HER2_2.csv")
files_HER2 <- c("CIBERSORTx_Job43_Results_HER2_1.csv", "CIBERSORTx_Job44_Results_HER2_2.csv")

cibersort_res_HER2 <- files_HER2 %>%
  lapply(readr::read_csv) %>%
  bind_rows()

HER2_clinical <- clin_complete %>% dplyr::filter(SSP.PAM50 == "HER2-enriched")

HER2_clinical_merged <- merge(HER2_clinical, cibersort_res_HER2, by.x = "GEX.assay", by.y = "Mixture")

HER2_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "HER2-enriched") %>% pull(GEX.assay)
HER2_gex <- gex %>% dplyr::select(all_of(HER2_pt))

SRGN <- c("SRGN")
df_HER2 <- as.data.frame(HER2_gex)[intersect(rownames(HER2_gex), SRGN),]

df_HER2 <- as.data.frame(t(df_HER2))
df_HER2$Sample_ID <- rownames(df_HER2)

HER2_clinical_merged <- merge(df_HER2, HER2_clinical_merged, by.y = "GEX.assay", by.x = "Sample_ID")

HER2_clinical_merged <- HER2_clinical_merged %>%
  mutate(SRGN_log2 = log2(1 + SRGN)
  )


pdf("SRGN_Immune_Correlation_HER2_all.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types) {
    
    df <- data.frame(
      x = HER2_clinical_merged[[g]],
      y = HER2_clinical_merged[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "spearman",   # Try with pearson as well
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#1F78B4") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#FF0000", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()




# files from cibersortx, basal


cibersort_res_Basal <- read.csv("CIBERSORTx_Job45_Results_Basal_1.csv")
names(cibersort_res_Basal) <- gsub("\\.", " ", names(cibersort_res_Basal))
Basal_clinical <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Basal-like")
Basal_clinical_merged <- merge(Basal_clinical, cibersort_res_Basal, by.x = "GEX.assay", by.y = "Mixture")
Basal_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Basal-like") %>% pull(GEX.assay)
Basal_gex <- gex %>% dplyr::select(all_of(Basal_pt))
SRGN <- c("SRGN")
df_Basal <- as.data.frame(Basal_gex)[intersect(rownames(Basal_gex), SRGN),]
df_Basal <- as.data.frame(t(df_Basal))
df_Basal$Sample_ID <- rownames(df_Basal)
Basal_clinical_merged <- merge(df_Basal, Basal_clinical_merged, by.y = "GEX.assay", by.x = "Sample_ID")
Basal_clinical_merged <- Basal_clinical_merged %>%
  mutate(SRGN_log2 = log2(1 + SRGN)
  )

##3 Doing correlations, basal

library(ggpubr)
library(ggplot2)
names(Basal_clinical_merged)
cell_types <- c(
  "B cells naive", "B cells memory", "Plasma cells",
  "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta", "NK cells resting",
  "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1",
  "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils",
  "Neutrophils"
)

cell_types_Basal <- c(
  "B cells naive", "B cells memory", "Plasma cells",
  "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory  Tregs ", "T cells gamma delta", "NK cells resting",
  "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1",
  "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils",
  "Neutrophils"
)

genes <- c("SRGN_log2")

# clean labels
gene_labels <- c(
  SRGN_log2 = "SRGN"
)

cell_types
(cell_types_Basal %in% names(Basal_clinical_merged))
#dev.off()
pdf("SRGN_Immune_Correlation_Basal.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types_Basal) {
    
    df <- data.frame(
      x = Basal_clinical_merged[[g]],
      y = Basal_clinical_merged[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "spearman",   # Try with pearson as well
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#4DBBD5") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#E64B35", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()

nejm <- c("#4DBBD5", "#E64B35")

pdf("SRGN_Immune_Correlation_Basal_pearson.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types_Basal) {
    
    df <- data.frame(
      x = Basal_clinical_merged[[g]],
      y = Basal_clinical_merged[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "pearson",   # Try with pearson as well
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#1F78B4") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#FF0000", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()

# normal-like

Normal_1 <- read.csv("CIBERSORTx_Job46_Results_Normal_1.csv")
Normal_2 <- read.csv("CIBERSORTx_Job47_Results_Normal_2.csv")
files_Normal <- c("CIBERSORTx_Job46_Results_Normal_1.csv", "CIBERSORTx_Job47_Results_Normal_2.csv")

cibersort_res_Normal <- files_Normal %>%
  lapply(readr::read_csv) %>%
  bind_rows()

Normal_clinical <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Normal-like")
Normal_clinical_merged <- merge(Normal_clinical, cibersort_res_Normal, by.x = "GEX.assay", by.y = "Mixture")
Normal_pt <- clin_complete %>% dplyr::filter(SSP.PAM50 == "Normal-like") %>% pull(GEX.assay)
Normal_gex <- gex %>% dplyr::select(all_of(Normal_pt))

SRGN <- c("SRGN")
df_Normal <- as.data.frame(Normal_gex)[intersect(rownames(Normal_gex), SRGN),]

df_Normal <- as.data.frame(t(df_Normal))
df_Normal$Sample_ID <- rownames(df_Normal)

Normal_clinical_merged <- merge(df_Normal, Normal_clinical_merged, by.y = "GEX.assay", by.x = "Sample_ID")
Normal_clinical_merged <- Normal_clinical_merged %>%
  mutate(SRGN_log2 = log2(1 + SRGN)
  )

genes <- c("SRGN_log2")

# clean labels
gene_labels <- c(
  SRGN_log2 = "SRGN"
)

#dev.off()
pdf("SRGN_Immune_Correlation_Normal.pdf", width = 14, height = 10)
for (g in genes) {
  plots <- list()
  plot_count <- 0
  
  for (ct in cell_types) {
    
    df <- data.frame(
      x = Normal_clinical_merged[[g]],
      y = Normal_clinical_merged[[ct]]
    )
    df <- na.omit(df)
    
    if (nrow(df) < 10) next
    if (sd(df$x) == 0 || sd(df$y) == 0) next
    
    plot_count <- plot_count + 1
    
    p <- ggscatter(
      df,
      x = "x",
      y = "y",
      add = "reg.line",
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.method = "spearman",   # Try with pearson as well
      xlab = gene_labels[g],
      ylab = ct,
      title = paste0(gene_labels[g], " vs ", ct)
    ) +
      # make points transparent
      geom_point(size = 1.2, alpha = 0.4, color = "#1F78B4") +
      # make regression line clearly visible
      geom_smooth(method = "lm", se = TRUE, color = "#FF0000", linewidth = 0.8) +
      theme(
        plot.title = element_text(size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)
      )
    
    plots[[plot_count]] <- p
    
    if (plot_count %% 20 == 0) {
      print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
      plots <- list()
      plot_count <- 0
    }
  }
  
  if (length(plots) > 0) {
    print(ggarrange(plotlist = plots, ncol = 5, nrow = 4))
  }
}

dev.off()

rm(b, Basal_clinical, Basal_clinical_merged, Basal_cpm, Basal_df, Basal_gex, c, chunk_data, chunk_dt, cibersort_res, cibersort_res_Basal, cibersort_res_HER2, cibersort_res_LumA, cibersort_res_Lumb, cibersort_res_Normal, combined, d, data, df, df_Basal, df_HER2, df_Lumb, df_Normal, df_sub, dgelist, endpoints, fit_RF_PAM50, h, headers, HER2_clinical, HER2_clinical_merged, HER2_cpm, HER2_gex, j, k, lcpm, luma_clinical, luma_clinical_merged, luma_cpm, luma_gex, Luma_gex, luma2, luma3, Lumb_1, Lumb_2, Lumb_3, LUmb_4, LUmb_5, lumb_clinical, Lumb_clinical, Lumb_clinical_merged, Lumb_clinical_merged_2, lumb_clinical_merged, Lumb_cpm, lumb_gex, Lumb_gex, merged, mergedall, n, n2, Normal_1, Normal_2, Normal_clinical, Normal_clinical_merged, Normal_cpm, Normal_gex, outcomes, p_OS_PAM50, p_OS_PAM50_all, p_OS_SRGN_all, p_RF_PAM50, p_RF_PAM50_all, p_SRGN_OS, palette_map, PAM50_df, plots, q, ref_rows, renamed_treatmet, s, SCAN_B, SSP_PAM50_untreated_df, test, treatment_group, x, age_col, bad_idx, Basal_pt, cell_types, cell_types_Basal, chunk_size, cols_to_keep, ct, density, end_col, end_col_Lumb, ep, event_var, events_used, filename, files, files_Basal, files_HER2, files_Luma)
rm(files_Lumb, files_Normal, form, g, gene_labels, gene_vargenes, gorup_var, group_vars, HER2_pt, i, keep.genes, lines, luma_pt, Luma_pt, lumb_pt, Lumb_pt, n_groups, n_used, Normal_pt, num_patients, num_patients_Basal, num_patients_HER2, num_patients_Lumb, num_patients_Normal, OS_SRGN_cox, outcome_name, p , p_cox_os_all, p_cox_os_untreated, p_density, p_hist_SRGN, p_SSP_PAM50, p_SSP_PAM50_untreated, p_treatgroup, p1, p2, pal, plot_count, plot_RF_cox, plot6_cox, required_cols, RF_SRGN_cox, RF_SRGN_cox_untreated, SRGN, start_col, start_col_Lumb, starts, starts_Basal, starts_HER2, starts_Lumb, starts_Normal, time_var, unique_patient_ids, untreated_patient_ids, vars_needed)
rm(clin_none, Luma_cpm, Lumb_4, Lumb_5, renamed_treatment, gene_var, genes, group_var)
gc()
