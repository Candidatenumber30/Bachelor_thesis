library(tidyverse)
library(scCustomize)
library(dplyr)
library(edgeR)
library(genefu)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)

# Load Clinical data, Raw expression data and gene annotation data

load("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Patient Data/7348_Patients/All_patient_data/StringTie prepDE gene count data unadjusted/SCANB.9142.matrixprepDEgenecount.Rdata")
load("C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/Patient Data/7348_Patients/All_patient_data/Gene Annotations/Gene.ID.ann.Rdata")
load("C:/Users/ssa214/UiT Office 365/O365-O365-Thesis-Emma - General/Emma/files/supplementary_data_table_1_20230113.Rdata")


scanb_meta <- supplementary_data_table_1_20230113

raw_expr <- as.data.frame(matrixprepDEgenecount)
raw_expr <- rownames_to_column(raw_expr, var = "ENSEMBL")
raw_expr <- distinct(raw_expr, ENSEMBL, .keep_all = T)


results <- Updated_HGNC_Symbols(input_data = anno$HGNC)

table(duplicated(Gene.ID.ann$HGNC))
anno <- Gene.ID.ann[,c("Gene.ID", "Gene.Name")]
raw_expr <- merge(anno, raw_expr, by.y = "ENSEMBL", by.x = "Gene.ID")
raw_expr <- raw_expr %>% distinct(Gene.Name, .keep_all = T)
raw_expr <- column_to_rownames(raw_expr, var = "Gene.Name")
rownames(raw_expr) <- raw_expr$Gene.Name
raw_expr$Gene.ID <- NULL

# Filtering of Clinical metadata

# 1. Remove LN and Normal
filtered_meta <- scanb_meta %>%
  filter(!(SpecimenType %in% c("LN", "Normal")))

# 2. Keep only samples present in expression matrix
filtered_samples <- intersect(filtered_meta$GEX.assay, colnames(raw_expr))

filtered_expr <- raw_expr %>%
  dplyr::select(all_of(filtered_samples))

filtered_meta <- filtered_meta %>%
  filter(GEX.assay %in% filtered_samples)

# 3. Identify patients with >1 sample
dup_patients <- filtered_meta %>%
  count(Patient) %>%
  filter(n > 1) %>%
  pull(Patient)

length(dup_patients)

# 4. Metadata for duplicated patients only
filtered_meta_dup <- filtered_meta %>%
  filter(Patient %in% dup_patients)

# 5. Keep one sample per duplicated patient using priority
filtered_meta_dup_unique <- filtered_meta_dup %>%
  mutate(
    priority = case_when(
      BiopsyType == "CoreBiopsy" ~ 1,
      BiopsyType == "CoreBiopsy2nd" ~ 2,
      BiopsyType == "OP" ~ 3,
      TRUE ~ 4
    )
  ) %>%
  arrange(Patient, priority) %>%
  group_by(Patient) %>%
  dplyr::slice(1) %>%
  ungroup()

dup_unique_gex <- filtered_meta_dup_unique$GEX.assay

# 6. Keep all samples from patients who were never duplicated
unique_gex <- filtered_meta %>%
  filter(!(Patient %in% dup_patients)) %>%
  pull(GEX.assay)

# 7. Final sample list
gex_to_keep <- c(unique_gex, dup_unique_gex)

length(gex_to_keep)

# 8. Final metadata and expression matrix
final_meta <- filtered_meta %>%
  dplyr::filter(GEX.assay %in% gex_to_keep)

# Filter expression data to include samples in filtered Clinical metadata

final_expr <- filtered_expr %>%
  dplyr::select(all_of(gex_to_keep))

qs_save(final_meta, "C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCAN_B_Clinical_Data.qs2")
qs_save(final_expr, "C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCAN_B_Raw_Expression.qs2")

# TMM normalization
dgelist <- DGEList(as.matrix(final_expr))
keep.genes <- filterByExpr(dgelist)
dgelist <- dgelist[keep.genes, , keep=FALSE]
dgelist <- calcNormFactors(dgelist, method = "TMM")

cpm <- cpm(dgelist, log = F)
lcpm <- cpm(dgelist, log = T, prior.count = 1)

# Save normalized CPM and LCPM expression data
qs_save(cpm, "C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCANB_Linear_CPM.qs2")
qs_save(lcpm, "C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCANB_Log_CPM.qs2")

# Get gene annotations
expr <- as.matrix(t(lcpm))
x <- rownames(lcpm)
anno <- Gene.ID.ann %>% dplyr::filter(Gene.Name %in% colnames(expr)) %>% distinct(Gene.Name, .keep_all = T)
expr <- as.data.frame(expr)

# Make table of gene annotations
all(anno$Gene.Name == colnames(expr))
annot <- data.frame(
  probe         = colnames(expr),
  Gene.symbol = colnames(expr),
  row.names     = colnames(expr),
  EntrezGene.ID = anno$EntrezGene,
  stringsAsFactors = FALSE
)
# Molecular PAM50 subtypying by Genefu
data("pam50.robust")
pam <- genefu::molecular.subtyping(
  sbt.model  = "pam50",
  data       = expr,       # rows = ENTREZID, cols = samples
  annot      = annot,
  #do.mapping = TRUE,        # <-- important
  verbose    = TRUE
)

subtype <- as.data.frame(pam$subtype)
subtype$Patient_ID <- rownames(subtype)
names(subtype) <- c("Genefu_PAM50", "GEX.assay")

# Merge with metadata
meta <- merge(subtype, final_meta, by = "GEX.assay")

qs_save(meta, "C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCAN_B_Clinical_Data.qs2")