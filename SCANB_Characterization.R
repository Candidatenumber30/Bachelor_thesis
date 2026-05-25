library(qs2)
library(tidyverse)
library(ggsci)
library(patchwork)
# Read clinical data
clinical <- qs_read("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCAN_B_Final_Clinical_Data.qs2")

# Number of patients
unique_patient_ids <- length(unique(clinical$Patient))

# Piechart of PAM50 subtypes distribution
unique(clinical$SSP.PAM50)
PAM50_df <- clinical %>% group_by(SSP.PAM50) %>% count()
PAM50_df$percentage <- round(PAM50_df$n / sum(PAM50_df$n)*100, digits = 0)
PAM50_df <- PAM50_df %>%
  mutate(
    percentage = percentage,
    label_text = paste0(percentage, "%"),
    lab_ypos = cumsum(percentage)-percentage / 2
  )

p_PAM50 <- ggplot(PAM50_df, aes(x = "", y = n, fill = SSP.PAM50)) +
  geom_col() +
  coord_polar(theta = "y") +
  labs(fill = "PAM50 subtype") +
  scale_fill_manual(values = pal_jco()(5), 
  labels = c("LumA" = "Luminal A", "LumB" = "Luminal B", "Her2" = "HER2-enriched", "Basal" = "Basal-like", "Normal" = "Normal-like"))+
  geom_text(aes(x=1.6, label = paste0(
    n,"(",scales::percent(n / sum(n), accuracy = 1),")")), 
    position = position_stack(vjust = 0.5), color = "black", size = 6)+
  theme_void()+
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )
p_PAM50
ggsave("Piechart_PAM50.pdf", p_PAM50,  height = 8, width = 11)

# Piechart of pT stage distribution
pT_df <- clinical %>% group_by(pT_group) %>% count()
p_pT <- ggplot(pT_df, aes(x = "", y = n, fill = pT_group)) +
  geom_col() +
  coord_polar(theta = "y") +
  labs(fill = "pT") +
  scale_fill_manual(values = pal_jco()(8))+
  geom_text(aes(x=1.6, label = paste0(
    n,"(",scales::percent(n / sum(n), accuracy = 1),")")),, 
            position = position_stack(vjust = 0.5), color = "black", size = 6)+
  theme_void()+
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )
p_pT
ggsave("Piechart_pT.pdf", p_pT, height = 8, width = 11)

# Piechart of LN status distribution
LN_df <- clinical %>% group_by(LN) %>% count()
p_LN <- ggplot(renamed_LN, aes(x = "", y = n, fill = LN)) +
  geom_col() +
  coord_polar(theta = "y") +
  labs(fill = "LN") +
  scale_fill_manual(values = pal_jco()(9))+
  geom_text(aes(x=1.6, label = paste0(
    n,"(",scales::percent(n / sum(n), accuracy = 1),")")),
            position = position_stack(vjust = 0.5), color = "black", size = 6)+
  theme_void()+
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
    )
p_LN
ggsave("Piechart_LN.pdf", p_LN, height = 8, width = 11)
