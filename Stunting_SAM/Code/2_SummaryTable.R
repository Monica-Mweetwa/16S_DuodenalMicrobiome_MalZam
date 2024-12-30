#title: "Table 1 - Zambia - SAM/ stunting"
#author: Monica N Mweetwa
#date: 30/12/2024
setwd("/Users/monica_mweetwa/Library/CloudStorage/Dropbox/BEECH_16s_Analysis/Stunting_SAM")
load("Data/RData/phyloseq_dataset_species.RData") #ps0
library(readr)
library(data.table) ## to use %like% function
library(tidyverse)

#Median IQR
#Variables of interest
tab1vars <- c('Sex', 'hiv_char', 'AgeMonths_endo', 'WAZ_endo', 'LAZ_endo', 'WHZ_endo', 'CD', 'VH', 'ESA', 'GLP.2_endo','GastricPh_endo', 'I.FABP_new', 'LPS_endo', 'LBP_endo', 'sCD14_new')

tab1 <- tableone::CreateTableOne(data=dat2, vars=tab1vars, strata= c("mal.type"))
tab1 <- as.data.frame(print(tab1, nonnormal=tab1vars, contDigits=1))
tab1 <- tibble::rownames_to_column(tab1, "row_names")
write.csv(tab1, "Output/DescriptiveTable_all.csv")
