library(tidyverse)
library(survival)
library(survminer)
library(dplyr)

# Load file with expression data merged with clinical data
load("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Emma/files/SCAN_B_clinical_merged.RData")
merged <- as.data.frame(merged)

# Stratify into high and low SRGN expression groups
merged$SRGN_class <- ifelse(merged$SRGN > median(merged$SRGN), "High", "Low")
table(merged$SRGN_class)

nejm <- c("#4DBBD5", "#E64B35")

# OS whole cohort

# Kaplan-Meier OS
fit_SRGN <- survfit(
  Surv(OS_years, OS_event) ~ SRGN_class, 
  data = merged
)

p_OS_SRGN <- ggsurvplot(fit_SRGN,  linewidth = 1,  
                            linetype = "strata", 
                            break.time.by = 1,
                            palette = nejm,
                            conf.int = TRUE, 
                            risk.table = T, 
                            pval = TRUE, 
                            pval.coord = c(0, 0.82),
                            xlim = c(0,5),
                            ylim= c(0.8, 1.00),
                            ylab = "Overall survival probability",
                            legend.title = expression(italic(SRGN) * " expression"),
                            legend.labs = c("High", "Low")
)

p_OS_SRGN

pdf("KM_OS_SRGN.pdf", width = 11, height = 8)
print(p_OS_SRGN, newpage = FALSE)
dev.off()

# Rename variables
merged <- merged %>%
  rename(
    `pT stage` = pT_group,
    `Age` = Age_group,
    `LN status` = LN_binary,
    `SRGN class` = SRGN_class,
  )

# Cox Hazard forest plot OS
cox_os_SRGN_all <- coxph(
  Surv(OS_years, OS_event) ~
    `Age` +
    `pT stage` +
    `LN status` +
    `SRGN class`,
  data = merged
)
cox_os_SRGN_all

p_cox_os_SRGN_all <- ggforest(cox_os_SRGN_all, 
                          data = merged,
                          cpositions = c(0.02, 0.22, 0.4),
                          fontsize = 1.4,
                          refLabel = "Reference")
p_cox_os_SRGN_all

ggsave("p_cox_os_all.pdf", p_cox_os_SRGN_all, height = 9, width = 10)


# RFi whole cohort

# Kaplan-Meier RFi

fit_RF_SRGN <- survfit(
  Surv(RFi_years, RFi_event) ~ SRGN_class, 
  data = merged
)

p_RF_SRGN <- ggsurvplot(fit_RF_SRGN,  linewidth = 1,  
                        linetype = "strata", 
                        break.time.by = 1,
                        palette = nejm,
                        conf.int = TRUE, 
                        risk.table = T, 
                        pval = TRUE, 
                        pval.coord = c(0, 0.82),
                        xlim = c(0,5),
                        ylim= c(0.8, 1.00),
                        ylab = "Recurrence-free probability",
                        legend.title = expression(italic(SRGN) * " expression"),
                        legend.labs = c("High", "Low")
)
p_RF_SRGN

pdf("KM_RF_SRGN.pdf", width = 11, height = 8)
print(p_RF_SRGN, newpage = FALSE)
dev.off()

# Cox forest plot RFi

cox_RF_all <- coxph(
  Surv(RFi_years, RFi_event) ~
    `Age` +
    `pT stage` +
    `LN status` +
    `SRGN class`,
  data = merged
)

p_cox_RF_all <- ggforest(cox_RF_all, 
                          data = merged,
                          cpositions = c(0.02, 0.22, 0.4),
                          fontsize = 1.4,
                          refLabel = "Reference")
p_cox_RF_all
ggsave("Cox_RF_all.pdf", p_cox_RF_all, height = 9, width = 10)


# Basal-like

# filter Basal-like samples
basal_clinical_merged <- merged %>% dplyr::filter(SSP.PAM50 == "Basal-like")
basal_clinical_merged$SRGN_class <- ifelse(basal_clinical_merged$SRGN > median(basal_clinical_merged$SRGN), "High", "Low")

# Kaplan-Meier OS
fit_basal_SRGN <- survfit(
  Surv(OS_years, OS_event) ~ SRGN_class, 
  data = basal_clinical_merged
)

p_OS_basal_SRGN <- ggsurvplot(fit_basal_SRGN,  linewidth = 1,  
                        linetype = "strata", 
                        break.time.by = 1,
                        palette = nejm,
                        conf.int = TRUE, 
                        risk.table = T, 
                        pval = TRUE, 
                        pval.coord = c(0, 0.6),
                        xlim = c(0,5),
                        ylim= c(0.5, 1.00),
                        ylab = "Overall survival probability",
                        legend.title = expression(italic(SRGN) * " expression"),
                        legend.labs = c("High", "Low")
)

p_OS_basal_SRGN

pdf("KM_OS_basal_SRGN.pdf", width = 11, height = 8)
print(p_OS_basal_SRGN, newpage = FALSE)
dev.off()

# Kaplan-Meier RFi

fit_RF_basal_SRGN <- survfit(
  Surv(RFi_years, RFi_event) ~ SRGN_class, 
  data = basal_clinical_merged
)

p_RF_basal_SRGN <- ggsurvplot(fit_RF_basal_SRGN,  linewidth = 1,  
                        linetype = "strata", 
                        break.time.by = 1,
                        palette = nejm,
                        conf.int = TRUE, 
                        risk.table = T, 
                        pval = TRUE, 
                        pval.coord = c(0, 0.6),
                        xlim = c(0,5),
                        ylim= c(0.5, 1.00),
                        ylab = "Recurrence-free probability",
                        legend.title = expression(italic(SRGN) * " expression"),
                        legend.labs = c("High", "Low")
)
p_RF_basal_SRGN

pdf("Basal_KM_RF_.pdf", width = 11, height = 8)
print(p_RF_basal_SRGN, newpage = FALSE)
dev.off

# Merge T1 and Tis groups

basal_clinical_merged$pT_group <- fct_collapse(
  basal_clinical_merged$pT_group,
  T1 = c("T1", "Tis")
)

# Rename variables
basal_clinical_merged <- basal_clinical_merged %>%
  rename(
    `pT stage` = pT_group,
    `Age` = Age_group,
    `LN status` = LN_binary,
    `SRGN class` = SRGN_class,
  )

# Forest plot of Cox RFi
cox_RF_basal_SRGN <- coxph(
  Surv(RFi_years, RFi_event) ~
    `Age` +
    `pT stage` +
    `LN status` +
    `SRGN class`,
  data = basal_clinical_merged
)

p_cox_RF_basal_SRGN <- ggforest(cox_RF_basal_SRGN, 
                            data = basal_clinical_merged,
                            cpositions = c(0.02, 0.22, 0.4),
                            fontsize = 1.4,
                            refLabel = "Reference")
p_cox_RF_basal_SRGN
ggsave("Cox_RF_basal_SRGN.pdf", p_cox_RF_basal_SRGN, height = 9, width = 10)

# Cox Forest plot OS

cox_OS_basal_SRGN <- coxph(
  Surv(OS_years, OS_event) ~
    `Age` +
    `pT stage` +
    `LN status` +
    `SRGN class`,
  data = basal_clinical_merged
)

p_cox_OS_basal_SRGN <- ggforest(cox_OS_basal_SRGN, 
                                data = basal_clinical_merged,
                                cpositions = c(0.02, 0.22, 0.4),
                                fontsize = 1.4,
                                refLabel = "Reference")
p_cox_OS_basal_SRGN
ggsave("Cox_OS_basal_SRGN.pdf", p_cox_OS_basal_SRGN)



# Luminal A, OS

luma_clinical_merged <- merged %>% dplyr::filter(SSP.PAM50 == "Luminal A")
luma_clinical_merged$SRGN_class <- ifelse(luma_clinical_merged$SRGN > median(luma_clinical_merged$SRGN), "High", "Low")

fit_luma_SRGN <- survfit(
  Surv(OS_years, OS_event) ~ SRGN_class, 
  data = luma_clinical_merged
)

p_OS_luma_SRGN <- ggsurvplot(fit_luma_SRGN,  linewidth = 1,  
                              linetype = "strata", 
                              break.time.by = 1,
                              palette = nejm,
                              conf.int = TRUE, 
                              risk.table = T, 
                              pval = TRUE, 
                              pval.coord = c(0, 0.6),
                              xlim = c(0,5),
                              ylim= c(0.5, 1.00),
                              ylab = "Overall survival probability",
                              legend.title = expression(italic(SRGN) ~ " expression"),
                              legend.labs = c("High", "Low")
)

p_OS_luma_SRGN

pdf("KM_OS_luma_SRGN.pdf", width = 11, height = 8)
print(p_OS_luma_SRGN)
dev.off()


# luminal A, RFi

fit_RF_luma_SRGN <- survfit(
  Surv(RFi_years, RFi_event) ~ SRGN_class, 
  data = luma_clinical_merged
)

p_RF_luma_SRGN <- ggsurvplot(fit_RF_luma_SRGN,  linewidth = 1,  
                              linetype = "strata", 
                              break.time.by = 1,
                              palette = nejm,
                              conf.int = TRUE, 
                              risk.table = T, 
                              pval = TRUE, 
                              pval.coord = c(0, 0.6),
                              xlim = c(0,5),
                              ylim= c(0.5, 1.00),
                              ylab = "Recurrence-free interval survival probability",
                              legend.title = expression(italic(SRGN) ~ " expression"),
                              legend.labs = c("High", "Low")
)
p_RF_luma_SRGN

pdf("KM_RF_luma_SRGN.pdf", width = 11, height = 8)
print(p_RF_luma_SRGN)
dev.off()

luma_clinical_merged$pT_group <- fct_collapse(
  luma_clinical_merged$pT_group,
  T1 = c("T1", "Tis")
)

cox_RF_luma_SRGN <- coxph(
  Surv(RFi_years, RFi_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = luma_clinical_merged
)

p_cox_RF_luma_SRGN <- ggforest(cox_RF_luma_SRGN, 
                                data = luma_clinical_merged,
                                cpositions = c(0.02, 0.22, 0.4),
                                fontsize = 0.8,
                                refLabel = "Reference")
p_cox_RF_luma_SRGN
ggsave("Cox_RF_luma_SRGN.pdf", p_cox_RF_luma_SRGN)

cox_OS_luma_SRGN <- coxph(
  Surv(OS_years, OS_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = luma_clinical_merged
)

p_cox_OS_luma_SRGN <- ggforest(cox_OS_luma_SRGN, 
                                data = luma_clinical_merged,
                                cpositions = c(0.02, 0.22, 0.4),
                                fontsize = 0.8,
                                refLabel = "Reference")
p_cox_OS_basal_SRGN
ggsave("Cox_OS_luma_SRGN.pdf", p_cox_OS_luma_SRGN)


# Luminal B, OS


lumb_clinical_merged <- merged %>% dplyr::filter(SSP.PAM50 == "Luminal B")
lumb_clinical_merged$SRGN_class <- ifelse(lumb_clinical_merged$SRGN > median(lumb_clinical_merged$SRGN), "High", "Low")

fit_lumb_SRGN <- survfit(
  Surv(OS_years, OS_event) ~ SRGN_class, 
  data = lumb_clinical_merged
)

p_OS_lumb_SRGN <- ggsurvplot(fit_lumb_SRGN,  linewidth = 1,  
                             linetype = "strata", 
                             break.time.by = 1,
                             palette = nejm,
                             conf.int = TRUE, 
                             risk.table = T, 
                             pval = TRUE, 
                             pval.coord = c(0, 0.6),
                             xlim = c(0,5),
                             ylim= c(0.5, 1.00),
                             ylab = "Overall survival probability",
                             legend.title = expression(italic(SRGN) ~ " expression"),
                             legend.labs = c("High", "Low")
)

p_OS_lumb_SRGN

pdf("KM_OS_lumb_SRGN.pdf", width = 11, height = 8)
print(p_OS_lumb_SRGN)
dev.off()


# luminal B, RFi

fit_RF_lumb_SRGN <- survfit(
  Surv(RFi_years, RFi_event) ~ SRGN_class, 
  data = lumb_clinical_merged
)

p_RF_lumb_SRGN <- ggsurvplot(fit_RF_lumb_SRGN,  linewidth = 1,  
                             linetype = "strata", 
                             break.time.by = 1,
                             palette = nejm,
                             conf.int = TRUE, 
                             risk.table = T, 
                             pval = TRUE, 
                             pval.coord = c(0, 0.6),
                             xlim = c(0,5),
                             ylim= c(0.5, 1.00),
                             ylab = "Recurrence-free interval survival probability",
                             legend.title = expression(italic(SRGN) ~ " expression"),
                             legend.labs = c("High", "Low")
)
p_RF_lumb_SRGN

pdf("KM_RF_lumb_SRGN.pdf", width = 11, height = 8)
print(p_RF_lumb_SRGN)
dev.off()

lumb_clinical_merged$pT_group <- fct_collapse(
  lumb_clinical_merged$pT_group,
  T1 = c("T1", "Tis")
)

cox_RF_lumb_SRGN <- coxph(
  Surv(RFi_years, RFi_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = lumb_clinical_merged
)

p_cox_RF_lumb_SRGN <- ggforest(cox_RF_lumb_SRGN, 
                               data = lumb_clinical_merged,
                               cpositions = c(0.02, 0.22, 0.4),
                               fontsize = 0.8,
                               refLabel = "Reference")
p_cox_RF_lumb_SRGN
ggsave("Cox_RF_lumb_SRGN.pdf", p_cox_RF_lumb_SRGN)

cox_OS_lumb_SRGN <- coxph(
  Surv(OS_years, OS_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = lumb_clinical_merged
)

p_cox_OS_lumb_SRGN <- ggforest(cox_OS_lumb_SRGN, 
                               data = lumb_clinical_merged,
                               cpositions = c(0.02, 0.22, 0.4),
                               fontsize = 0.8,
                               refLabel = "Reference")
p_cox_OS_basal_SRGN
ggsave("Cox_OS_lumb_SRGN.pdf", p_cox_OS_lumb_SRGN)



# HER2-enriched, OS Kaplan-Meier

her2_clinical_merged <- merged %>% dplyr::filter(SSP.PAM50 == "HER2-enriched")
her2_clinical_merged$SRGN_class <- ifelse(her2_clinical_merged$SRGN > median(her2_clinical_merged$SRGN), "High", "Low")

fit_her2_SRGN <- survfit(
  Surv(OS_years, OS_event) ~ SRGN_class, 
  data = her2_clinical_merged
)

p_OS_her2_SRGN <- ggsurvplot(fit_her2_SRGN,  linewidth = 1,  
                             linetype = "strata", 
                             break.time.by = 1,
                             palette = nejm,
                             conf.int = TRUE, 
                             risk.table = T, 
                             pval = TRUE, 
                             pval.coord = c(0, 0.6),
                             xlim = c(0,5),
                             ylim= c(0.5, 1.00),
                             ylab = "Overall survival probability",
                             legend.title = expression(italic(SRGN) ~ " expression"),
                             legend.labs = c("High", "Low")
)

p_OS_her2_SRGN

pdf("KM_OS_her2_SRGN.pdf", width = 11, height = 8)
print(p_OS_her2_SRGN)
dev.off()

# HER2-enriched, RFi Kaplan-Meier

fit_RF_her2_SRGN <- survfit(
  Surv(RFi_years, RFi_event) ~ SRGN_class, 
  data = her2_clinical_merged
)

p_RF_her2_SRGN <- ggsurvplot(fit_RF_her2_SRGN,  linewidth = 1,  
                             linetype = "strata", 
                             break.time.by = 1,
                             palette = nejm,
                             conf.int = TRUE, 
                             risk.table = T, 
                             pval = TRUE, 
                             pval.coord = c(0, 0.6),
                             xlim = c(0,5),
                             ylim= c(0.5, 1.00),
                             ylab = "Recurrence-free interval survival probability",
                             legend.title = expression(italic(SRGN) ~ " expression"),
                             legend.labs = c("High", "Low")
)
p_RF_her2_SRGN

pdf("KM_RF_her2_SRGN.pdf", width = 11, height = 8)
print(p_RF_her2_SRGN)
dev.off()

her2_clinical_merged$pT_group <- fct_collapse(
  her2_clinical_merged$pT_group,
  T1 = c("T1", "Tis")
)

# HER2-enriched RFi, Cox Hazard forest plot
cox_RF_her2_SRGN <- coxph(
  Surv(RFi_years, RFi_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = her2_clinical_merged
)

p_cox_RF_her2_SRGN <- ggforest(cox_RF_her2_SRGN, 
                               data = her2_clinical_merged,
                               cpositions = c(0.02, 0.22, 0.4),
                               fontsize = 0.8,
                               refLabel = "Reference")
p_cox_RF_her2_SRGN
ggsave("Cox_RF_her2_SRGN.pdf", p_cox_RF_her2_SRGN)

# HER2-enriched OS Cox Hazard ratio forest plot

cox_OS_her2_SRGN <- coxph(
  Surv(OS_years, OS_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = her2_clinical_merged
)

p_cox_OS_her2_SRGN <- ggforest(cox_OS_her2_SRGN, 
                               data = her2_clinical_merged,
                               cpositions = c(0.02, 0.22, 0.4),
                               fontsize = 0.8,
                               refLabel = "Reference")
p_cox_OS_basal_SRGN
ggsave("Cox_OS_her2_SRGN.pdf", p_cox_OS_her2_SRGN)

# Normal-like, OS


normal_clinical_merged <- merged %>% dplyr::filter(SSP.PAM50 == "Normal-like")
normal_clinical_merged$SRGN_class <- ifelse(normal_clinical_merged$SRGN > median(normal_clinical_merged$SRGN), "High", "Low")

fit_normal_SRGN <- survfit(
  Surv(OS_years, OS_event) ~ SRGN_class, 
  data = normal_clinical_merged
)

p_OS_normal_SRGN <- ggsurvplot(fit_normal_SRGN,  linewidth = 1,  
                             linetype = "strata", 
                             break.time.by = 1,
                             palette = nejm,
                             conf.int = TRUE, 
                             risk.table = T, 
                             pval = TRUE, 
                             pval.coord = c(0, 0.6),
                             xlim = c(0,5),
                             ylim= c(0.5, 1.00),
                             ylab = "Overall survival probability",
                             legend.title = expression(italic(SRGN) ~ " expression"),
                             legend.labs = c("High", "Low")
)

p_OS_normal_SRGN

pdf("KM_OS_normal_SRGN.pdf", width = 11, height = 8)
print(p_OS_normal_SRGN)
dev.off()

# Normal-like, RFi

fit_RF_normal_SRGN <- survfit(
  Surv(RFi_years, RFi_event) ~ SRGN_class, 
  data = normal_clinical_merged
)

p_RF_normal_SRGN <- ggsurvplot(fit_RF_normal_SRGN,  linewidth = 1,  
                             linetype = "strata", 
                             break.time.by = 1,
                             palette = nejm,
                             conf.int = TRUE, 
                             risk.table = T, 
                             pval = TRUE, 
                             pval.coord = c(0, 0.6),
                             xlim = c(0,5),
                             ylim= c(0.5, 1.00),
                             ylab = "Recurrence-free interval survival probability",
                             legend.title = expression(italic(SRGN) ~ " expression"),
                             legend.labs = c("High", "Low")
)
p_RF_normal_SRGN

pdf("KM_RF_normal_SRGN.pdf", width = 11, height = 8)
print(p_RF_normal_SRGN)
dev.off()

normal_clinical_merged$pT_group <- fct_collapse(
  normal_clinical_merged$pT_group,
  T1 = c("T1", "Tis")
)

cox_RF_normal_SRGN <- coxph(
  Surv(RFi_years, RFi_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = normal_clinical_merged
)

p_cox_RF_normal_SRGN <- ggforest(cox_RF_normal_SRGN, 
                               data = normal_clinical_merged,
                               cpositions = c(0.02, 0.22, 0.4),
                               fontsize = 0.8,
                               refLabel = "Reference")
p_cox_RF_normal_SRGN
ggsave("Cox_RF_normal_SRGN.pdf", p_cox_RF_normal_SRGN)

cox_OS_normal_SRGN <- coxph(
  Surv(OS_years, OS_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = normal_clinical_merged
)

p_cox_OS_normal_SRGN <- ggforest(cox_OS_normal_SRGN, 
                               data = normal_clinical_merged,
                               cpositions = c(0.02, 0.22, 0.4),
                               fontsize = 0.8,
                               refLabel = "Reference")
p_cox_OS_basal_SRGN
ggsave("Cox_OS_normal_SRGN.pdf", p_cox_OS_normal_SRGN)



# all subtypes, DRFi (Disease-Free survival)

fit_DRFi_SRGN <- survfit(
  Surv(DRFi_years, DRFi_event) ~ SRGN_class, 
  data = merged
)

p_DRFi_SRGN <- ggsurvplot(fit_DRFi_SRGN,  linewidth = 1,  
                        linetype = "strata", 
                        break.time.by = 1,
                        palette = nejm,
                        conf.int = TRUE, 
                        risk.table = T, 
                        pval = TRUE, 
                        pval.coord = c(0, 0.6),
                        xlim = c(0,5),
                        ylim= c(0.5, 1.00),
                        ylab = "Disease free survival probability",
                        legend.title = expression(italic(SRGN) * " expression"),
                        legend.labs = c("High", "Low")
)

p_DRFi_SRGN

pdf("KM_DRFi_SRGN.pdf", width = 11, height = 8)
print(p_DRFi_SRGN)
dev.off()


cox_DRFi_SRGN <- coxph(
  Surv(DRFi_years, DRFi_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = merged
)

p_cox_DRFi_SRGN <- ggforest(cox_DRFi_SRGN, 
                          data = merged,
                          cpositions = c(0.02, 0.22, 0.4),
                          fontsize = 0.8,
                          refLabel = "Reference")
p_cox_DRFi_SRGN
ggsave("Cox_DRFi_SRGN.pdf", p_cox_DRFi_SRGN)

# basal, DRFi

fit_basal_DRFi_SRGN <- survfit(
  Surv(DRFi_years, DRFi_event) ~ SRGN_class, 
  data = basal_clinical_merged
)

p_DRFi_basal_SRGN <- ggsurvplot(fit_basal_DRFi_SRGN,  linewidth = 1,  
                          linetype = "strata", 
                          break.time.by = 1,
                          palette = nejm,
                          conf.int = TRUE, 
                          risk.table = T, 
                          pval = TRUE, 
                          pval.coord = c(0, 0.6),
                          xlim = c(0,5),
                          ylim= c(0.5, 1.00),
                          ylab = "Disease free survival probability",
                          legend.title = expression(italic(SRGN) * " expression"),
                          legend.labs = c("High", "Low")
)

p_DRFi_basal_SRGN

pdf("KM_basal_DRFi_SRGN.pdf", width = 11, height = 8)
print(p_DRFi_basal_SRGN)
dev.off()


cox_basal_DRFi_SRGN <- coxph(
  Surv(DRFi_years, DRFi_event) ~
    Age_group +
    pT_group +
    LN_binary +
    SRGN_class,
  data = basal_clinical_merged
)

summary(cox_basal_DRFi_SRGN)

p_cox_basal_DRFi_SRGN <- ggforest(cox_basal_DRFi_SRGN, 
                            data = basal_clinical_merged,
                            cpositions = c(0.02, 0.22, 0.4),
                            fontsize = 0.8,
                            refLabel = "Reference")
p_cox_basal_DRFi_SRGN
ggsave("Cox_DRFi_SRGN.pdf", p_cox_DRFi_SRGN)