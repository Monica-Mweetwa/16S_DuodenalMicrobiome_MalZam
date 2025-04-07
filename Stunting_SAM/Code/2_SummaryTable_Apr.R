#title: "Table 1 - Zambia - SAM/ stunting"
#author: Monica N Mweetwa
#date: 30/12/2024
setwd("~/Dropbox/BEECH_16s_Analysis/Stunting_SAM") #Mac
setwd("C:/Users/Monica/Dropbox/BEECH_16s_Analysis/Stunting_SAM") #Windows
library(readr)
library(data.table) ## to use %like% function
library(tidyverse)

load("Data/RData/phyloseq_dataset_species.RData") #ps0

add <- data.frame(sample_data(ps0))
#create variable for different degrees of SAM
add2 <- add %>%
  mutate(hiv_char = ifelse(study == "BEECH children" & PID != "B116", "Negative", hiv_char),
         mal.type = ifelse(study == "EE children", "SAM", "Stunted"),
         Odema_char = as.character(Odema),
         MUAC_char = ifelse(MUAC >= 11.5, ">=11.5cm", "<11.5cm"),
         LAZ_char = ifelse(LAZ_endo > -2, "Not stunted", "Stunted"))
  

#Median IQR - stratified by malnutrition type
#Variables of interest
tab1vars <- c('Sex', 'hiv_char', 'WAZ_endo', 'LAZ_endo', 'WHZ_endo', 'AgeMonths_endo', 
              'CD', 'VH', 'ESA', 'GLP.2_endo','GastricPh_endo', 'I.FABP_new', 'LPS_endo', 
              'sCD14_new', 'Odema_char', 'MUAC_char', 'LAZ_char')

tab1 <- tableone::CreateTableOne(data=add2, vars=tab1vars, strata= c("mal.type"))
tab1 <- as.data.frame(print(tab1, nonnormal=tab1vars, contDigits=1))
tab1 <- tibble::rownames_to_column(tab1, "row_names")
write.csv(tab1, "Output/DescriptiveTable_all.csv")

#check for missing values.
pct_missing <-   add2 %>%
  group_by(mal.type)%>%
  select(Sex, hiv_char, WAZ_endo, LAZ_endo, WHZ_endo, AgeMonths_endo, 
         CD, VH, ESA, GLP.2_endo,GastricPh_endo, I.FABP_new, LPS_endo, 
         sCD14_new, Odema_char, MUAC_char) %>% 
  summarise_all(list(name = ~sum(is.na(.)))) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "clinfeat")
write.csv(pct_missing, "Output/NMissingDescriptiveData.csv")



