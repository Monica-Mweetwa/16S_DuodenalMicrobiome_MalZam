#title: "Table 1 - Zambia - SAM/ stunting"
#author: Monica N Mweetwa
#date: 01/07/2025
setwd("~/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Stunting_SAM") #Mac
#setwd("C:/Users/Monica/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Stunting_SAM") #Windows
library(readr)
library(data.table) ## to use %like% function
library(tidyverse)

load("Data/RData/Zam_phyloseqObj_AbsASV_Withtree.RData") #ps0_withtree

add <- data.frame(sample_data(ps0_withtree))
#create variable for different degrees of SAM
add2 <- add %>%
  mutate(mal.type = ifelse(study == "EE children", "SAM", "only stunted"),
         Odema_char = as.character(Odema),
         MUAC_char = ifelse(MUAC >= 11.5, ">=11.5cm", "<11.5cm"),
         LAZ_char = ifelse(LAZ_endo > -2, "Not stunted", "Stunted"),
         BFeeding_char2 = ifelse(BFeeding_char == "No", "No_1", "Yes"),
         Sex_num = ifelse(Sex == "Male", 1, 2))
#get sequencing depth
seq_dataframe <- ps0_withtree%>%
  ps_melt() %>%
  dplyr::group_by(PID, study) %>%
  summarise(SeqDepth = sum(Abundance))

add3 <- left_join(add2, seq_dataframe)


#Median IQR - stratified by malnutrition type
#Variables of interest
tab1vars <- c('Sex', 'hiv_char', 'WAZ_endo', 'LAZ_endo', 'WHZ_endo', 'AgeMonths_endo', 
              'CD', 'VH', 'ESA', 'GLP.2_endo','GastricPh_endo', 'I.FABP_new', 'LPS_endo', 
              'sCD14_new', 'Odema_char', 'MUAC_char', 'LAZ_char', 'BFeeding_char2', 'SeqDepth')

tab1 <- tableone::CreateTableOne(data=add3, vars=tab1vars, strata= c("mal.type"))
tab1 <- as.data.frame(print(tab1, nonnormal=tab1vars, contDigits=1))
tab1 <- tibble::rownames_to_column(tab1, "row_names")
write.csv(tab1, "Output/DescriptiveTable_all_2025.csv")

#check for missing values.
pct_missing <-   add3 %>%
  group_by(mal.type)%>%
  select(Sex, hiv_char, WAZ_endo, LAZ_endo, WHZ_endo, AgeMonths_endo, 
         CD, VH, ESA, GLP.2_endo,GastricPh_endo, I.FABP_new, LPS_endo, 
         sCD14_new, Odema_char, MUAC_char, BFeeding_char, SeqDepth) %>% 
  summarise_all(list(name = ~sum(is.na(.)))) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "clinfeat")
write.csv(pct_missing, "Output/NMissingDescriptiveData.csv")


#Correlation matrix of participant features used as covariates in analysis
library("PerformanceAnalytics")

add3 <- add2 %>%
  select(AgeMonths_endo, Sex_num, HIV_status, BFeeding)
chart.Correlation(add3, histogram=TRUE, pch=19)

#Correlation matrix of participant features used as covariates in analysis
library("Hmisc")
#SAM
lm_pred_vars <- add2 %>%
  filter(study == "EE children") %>%
  #mutate(mal.type_num = ifelse(mal.type == "SAM", 1, 2)) %>%
  select(AgeMonths_endo, Sex_num, HIV_status, BFeeding, WAZ_endo, LAZ_endo, WHZ_endo, AgeMonths_endo, 
         CD, VH, ESA, GLP.2_endo,GastricPh_endo, I.FABP_new, LPS_endo, sCD14_new, Odema, MUAC) %>%
  rename("AgeMonths_endo" = "Age","HIV_status" = "HIV status","Sex_num" = "Sex","BFeeding" = "Breast feeding",
         "WAZ_endo" = "WAZ", "LAZ_endo" = "LAZ", "WHZ_endo" = "WLZ", "GLP.2_endo" = "GLP2", "GastricPh_endo" = "Gastric pH",
         "I.FABP_new" = "iFABP", "LPS_endo" = "LPS", "sCD14_new" = "sCD14")
res2 <- rcorr(as.matrix(lm_pred_vars))
res2$r[is.na(res2$r)]=0
res2$p[is.na(res2$p)]=1

col<- colorRampPalette(c("blue", "white", "red"))(20)
corrplot(res2$r, type="upper", order="hclust", col=col,
         p.mat = res2$P, sig.level = 0.05, diag=FALSE,
         title='SAM - Participant Features Correlation Matrix')
#STUNTING
lm_pred_vars <- add2 %>%
  filter(study == "BEECH children") %>%
  #mutate(mal.type_num = ifelse(mal.type == "SAM", 1, 2)) %>%
  select(AgeMonths_endo, Sex_num, HIV_status, BFeeding, WAZ_endo, LAZ_endo, WHZ_endo, AgeMonths_endo, 
         CD, VH, ESA, GLP.2_endo,GastricPh_endo, I.FABP_new, LPS_endo, sCD14_new) %>%
  rename("AgeMonths_endo" = "Age","HIV_status" = "HIV status","Sex_num" = "Sex","BFeeding" = "Breast feeding",
         "WAZ_endo" = "WAZ", "LAZ_endo" = "LAZ", "WHZ_endo" = "WLZ", "GLP.2_endo" = "GLP2", "GastricPh_endo" = "Gastric pH",
         "I.FABP_new" = "iFABP", "LPS_endo" = "LPS", "sCD14_new" = "sCD14")
res2 <- rcorr(as.matrix(lm_pred_vars))
res2$r[is.na(res2$r)]=0
res2$p[is.na(res2$p)]=1

col<- colorRampPalette(c("blue", "white", "red"))(20)
corrplot(res2$r, type="upper", order="hclust", col=col,
         p.mat = res2$P, sig.level = 0.05, diag=FALSE,
         title='STUNTING - Participant Features Correlation Matrix')


#ggsave("Output/ParticipantFeatures_CorrelationPlot_2025Jul.png", plot = last_plot(), dpi = 600, device = "png", 
#       width = 10, height = 10, units = "in")
