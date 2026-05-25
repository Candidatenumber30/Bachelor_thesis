library(data.table)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library("clusterProfiler")

# Read expression data and metadata
dt <- fread("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/GDC-PANCAN.htseq_counts.tsv", sep = "\t")
meta <- fread("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/GDC-PANCAN.basic_phenotype.tsv")
meta_filtered <- meta %>% dplyr::filter(program == "TCGA")
table(meta_filtered$sample_type)

# Keep samples of Primary tumor
meta_filtered <- meta_filtered %>% dplyr::filter(sample_type == "Primary Tumor")

# Check remaining sample types
table(meta_filtered$sample_type)

# Get sample IDs to keep from metadata and columns that match in dt
samples_to_keep <- meta_filtered$sample
x <- intersect(samples_to_keep, names(dt))

# Choose Ensembl_ID of columns that are in filtered metadata
dt_filtered <- dt %>% dplyr::select(x, Ensembl_ID)
head(dt_filtered$Ensembl_ID)

# Add column of Ensemble_IDs where ommitted everything after .
dt_filtered$Ensembl_Clean <- sub("\\..*","", dt_filtered$Ensembl_ID)
head(dt_filtered$Ensembl_Clean)

# Remove duplicates of Ensembl_IDs
ensembl_ids <- unique(dt_filtered$Ensembl_Clean)

# Convert gene IDs
gene_map <- clusterProfiler::bitr(ensembl_ids,
                                  fromType = "ENSEMBL",
                                  toType = "SYMBOL",
                                  OrgDb = org.Hs.eg.db)

gene_map <- as.data.table(gene_map)

# Keep first SYMBOL per ENSEMBL
setorder(gene_map, ENSEMBL, SYMBOL)
gene_map_unique <- gene_map[, .SD[1], by = ENSEMBL]

# Merge Ensembl IDs with expression data
dt_mapped <- merge(dt_filtered, gene_map_unique,
                   by.x = "Ensembl_Clean", by.y = "ENSEMBL")

# Remove duplicate Ensembl IDs
dt_mapped$Ensembl_Clean <- NULL
dt_mapped$Ensembl_ID <- NULL
dt_mapped <- as.data.frame(dt_mapped)
dt_mapped <- dt_mapped[,c(9725,1:9724)]

# Keep distinct sample IDs
dt_final <- distinct(dt_mapped,SYMBOL,.keep_all=TRUE)

save(dt_final, meta_filtered, file = "Pancancer TCGA Files.RData")

load("Pancancer TCGA Files.RData")
counts <- column_to_rownames(dt_final,var = "SYMBOL")
all(names(counts) == meta_filtered$sample)
meta_filtered <- meta_filtered %>% dplyr::filter(sample %in% x)
group <- factor(meta_filtered$project_id)

# TMM normalization
library(edgeR)
d0 <- DGEList(counts, group = group)
keep.exprs <- filterByExpr(d0)
d0 <- d0[keep.exprs,, keep.lib.sizes=FALSE]
d0 <- calcNormFactors(d0, method = "TMM")
lcpm <- cpm(d0, log=TRUE, prior.count = 1)
save(lcpm, meta_filtered, file = "Pancancer_files.RData")

load("Pancancer_files.RData")

# Check if samples in lcpm match samples in metadata
all(names(lcpm) %in% meta_filtered$sample)
lcpm <- as.data.frame(lcpm)

# Take SRGN expression in lcpm, swwap rows and columns
SRGN <- lcpm["SRGN",]
SRGN_long <- t(SRGN)
SRGN_long <- as.data.frame(SRGN_long)

# Add column with sample IDs in SRGN
SRGN_long$sample_id <- rownames(SRGN_long)
head(SRGN_long$sample_id)

# Merge SRGN expression data with metadata by sample
merged <- merge(meta_filtered, SRGN_long, by.y = "sample_id", by.x = "sample")
merged <- as.data.frame(merged)

# Add column with cancer types
merged$cancer_type <- sub(".*\\-","", merged$project_id)

# order cancer types by median
merged_ordered <- merged %>%
  mutate(cancer_type = fct_reorder(cancer_type, SRGN, .fun = median))

save(merged, merged_ordered, file = "TCGA_Pancancer_files.RData")