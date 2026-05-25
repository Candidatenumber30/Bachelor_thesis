library(tidyverse)
library(survival)
library(survminer)
library(qs2)

pam_colors <- c("LumA" = "#1F78B4",
                "LumB" = "#A6CEE3",
                "Her2" = "sienna1",
                "Basal" = "#E31A1C",
                "Normal-like" = "#B2B2B2")

# Read clinical metadata
clinical <- qs_read("C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCAN_B_Clinical_Data.qs2")

# Variables required for further analysis
vars_needed <- c(
  "GEX.assay",
  "Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)",
  "Genefu_PAM50",
  "pT",
  "LN.spec",
  "OS_days",
  "RFi_days"
)

# Select samples with required variables in the clinical data
clin_complete <- clinical[
  complete.cases(clinical[, vars_needed]),
]

# Recode PAM50 SUbtype names
clin_complete$Genefu_PAM50 <- recode(clin_complete$Genefu_PAM50,
                                     "LumA" = "Luminal A",
                                     "LumB" = "Luminal B",
                                     "Her2" = "HER2-enriched",
                                     "Basal" = "Basal-like",
                                     "Normal" = "Normal-like"
)

clin_complete$SSP.PAM50 <- recode(clin_complete$SSP.PAM50,
                                  "LumA" = "Luminal A",
                                  "LumB" = "Luminal B",
                                  "Her2" = "HER2-enriched",
                                  "Basal" = "Basal-like",
                                  "Normal" = "Normal-like"
)

# Convert OS, RFi and DRFi survival in days into years
clin_complete$OS_years <- clin_complete$OS_days/365
clin_complete$RFi_years <- clin_complete$RFi_days/365
clin_complete$DRFi_years <- clin_complete$DRFi_days/365

# Create column with Lymph node as Negative or Positive
clin_complete$LN_binary <- NA  
clin_complete$LN_binary[!is.na(clin_complete$LN) & clin_complete$LN == 0] <- "Negative"
clin_complete$LN_binary[!is.na(clin_complete$LN) & clin_complete$LN == 1] <- "Positive"
clin_complete$LN_binary <- factor(clin_complete$LN_binary, levels = c("Negative", "Positive"))

# Merge Pathological tumor stage levels into T1, T2, T3_T4 and Tis
clin_complete$pT_group <- dplyr::case_when(
  clin_complete$pT %in% c("T1a","T1b","T1c","T1mi") ~ "T1",
  clin_complete$pT == "T2" ~ "T2",
  clin_complete$pT %in% c("T3","T4") ~ "T3_T4",
  clin_complete$pT %in% "Tis" ~ "Tis",
  TRUE ~ NA_character_
)

table(clin_complete$`Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)`, useNA = "ifany")
table(clin_complete$Genefu_PAM50, useNA = "ifany")
table(clin_complete$pT, useNA = "ifany")
table(clin_complete$LN.spec, useNA = "ifany")
table(clin_complete$LN_binary, useNA = "ifany")
table(is.na(clin_complete$OS_years))
table(is.na(clin_complete$RFi_years))

clin_complete$Sample_ID

# Create age group column with age categories <50, 50-65 and >65
age_mid <- as.numeric(as.character(
  clin_complete$`Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)`
))
clin_complete$Age_group <- cut(
  age_mid,
  breaks = c(-Inf, 50, 65, Inf),
  labels = c("<50", "50-65", ">65"),
  right = TRUE
)

# Order Molecular subtypes
clin_complete$Genefu_PAM50 <- factor(clin_complete$Genefu_PAM50, levels = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like", "Normal-like"))
clin_complete$SSP.PAM50 <- factor(clin_complete$SSP.PAM50, levels = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like", "Normal-like"))

# Save final clinical data
qs_save(clin_complete, file = "C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCAN_B_Final_Clinical_Data.qs2")

# Read raw expression data
gex <- qs_read("C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCAN_B_Raw_Expression.qs2")
gex <- as.data.frame(gex)

# Filter expression data to include only samples that have required clinical variables
samples_to_keep <- clin_complete$GEX.assay
gex <- gex %>% dplyr::select(all_of(samples_to_keep))

# Save final raw expression data
qs_save(gex, "C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCAN_B_Raw_Expression_Final.qs2")

# Read CPM expression data
gex <- qs_read("C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCANB_Linear_CPM.qs2")
gex <- as.data.frame(gex)

# Filter CPM data and save final CPM data
samples_to_keep <- clin_complete$GEX.assay
gex <- gex %>% dplyr::select(all_of(samples_to_keep))
qs_save(gex, "C:/Users/ssa214/UiT Office 365/O365-Knutsen Group - General/Patient Data/SCAN-B/SCANB_Linear_CPM_Final.qs2")
