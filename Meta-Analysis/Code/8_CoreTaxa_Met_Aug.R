########### All Malnourished Duodenum Analysis ##################
#Zambia vs Bangladesh vs Pakistan vs Central African Republic vs Madagascar
rm(list=ls())
library(microViz)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(vegan)
library(tidyr)
library(broom)
library(ANCOMBC)
#Import data
setwd("C:/Users/Monica/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis") #Windows
#setwd("~/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis") #Mac
rm(list=ls())
load("Data/RData/BEED_SEEM_BEECH_AbsASV.RData") #ps_SB_norm.ra.scaled

comparisons1 <- list(c("BEECH","BEED"), c("BEECH","SEEM"), c("BEED","SEEM"))
colors_study = c("BEECH" = "#1b9e77", 
                 "BEED" = "#d95f02",
                 "AFRIBIOTA" = "#ba4aa4",
                 "SEEM" = "#66a61e")
############################################################### Data cleaning
#Merge duplicate PIDs
ps_norm_merged <- phyloseq::merge_samples(ps_SB_norm.ra.scaled, "PID")
# replace metadata as NAs get introduced during merging
meta <- sample_data(ps_SB_norm.ra.scaled)
merged_meta <- meta[!duplicated(meta$PID), ]
merged_meta$sampleID2 <- rownames(merged_meta)
rownames(merged_meta) <- NULL
merged_meta2 <- merged_meta %>% column_to_rownames(var = "PID")
sample_data(ps_norm_merged) <- merged_meta2

ps_norm_merged2 = subset_taxa(ps_norm_merged, Class != "Chloroplast")
ps_norm_merged2 = subset_taxa(ps_norm_merged2, Kingdom == "Bacteria")
ps_norm_merged2 = subset_taxa(ps_norm_merged2, Family != "Mitochondria") #remove those assigned as Mitochondria Family
ps_norm_merged2 #3032 taxa

#edit taxonomy table to have ASVTaxa before sequences 
ps_norm_merged3 <- ps_norm_merged2
tax_tab <- as.data.frame(tax_table(ps_norm_merged3)) %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species, ASVTaxa, ASVReadable, ASVSeq)
tax_table(ps_norm_merged3) <- as.matrix(tax_tab)

############################################################### relative abundance plots
# all
ps_norm_merged3 %>%
  ps_filter(study != "ME") %>%
  tax_fix() %>%
  comp_barplot(tax_level = "Genus",
               label = "SAMPLE",
               n_taxa = 15) +
  #coord_flip() + # horizontal bars are often more readable
  facet_grid(
    cols = vars(study),
    scales = "free_x", # this only frees y scale per row in grid faceting
    space = "free_x" # allows bars to be same size by freeing facet heights
  ) +
  theme(axis.text.x = (element_blank()),
        axis.ticks.x = (element_blank()))
ggsave("Output/CoreTaxa/TopGenera_all.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 12, height = 5, units = "in")

#subset of interest
ps_norm_merged3 %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study != "ME") %>%
  tax_fix() %>%
  comp_barplot(tax_level = "Genus",
               label = "SAMPLE",
               n_taxa = 15) +
  #coord_flip() + # horizontal bars are often more readable
  facet_grid(
    cols = vars(study),
    scales = "free_x", # this only frees y scale per row in grid faceting
    space = "free_x" # allows bars to be same size by freeing facet heights
  ) +
  theme(axis.text.x = (element_blank()),
        axis.ticks.x = (element_blank()))
ggsave("Output/CoreTaxa/TopGenera_subset.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 8, height = 5, units = "in")

############################################################### Identify core taxa in each study 
######################################################################## present in at least 80% of samples

###### SEEM
#Filter for those with prevalence >= 80%
seem_phydat_rel <- ps_norm_merged3 %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study == "SEEM") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  dplyr::summarise(avg_abundance_seem = round((mean(Abundance)*100),5),
                   prevalence_seem = round((sum(Abundance > 0) / dplyr::n())*100,5),
                   SD_val_seem = round((sd(Abundance, na.rm = T))*100,5)) %>%
  filter(prevalence_seem >= 80) #16 core genera
write.csv(seem_phydat_rel, "Output/CoreTaxa/coreTax_SEEM_Aug.csv")

###### BEED
#filter for taxa with relative abundance of atleast 0.0001 in 1 or more samples
beed_phydat_rel <- ps_norm_merged3 %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study == "BEED") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  dplyr::summarise(avg_abundance_beed = round((mean(Abundance)*100),5),
                   prevalence_beed = round((sum(Abundance > 0) / dplyr::n())*100,5),
                   SD_val_beed = round((sd(Abundance, na.rm = T))*100,5)) %>%
  filter(prevalence_beed >= 80) #19 core genera
write.csv(beed_phydat_rel, "Output/CoreTaxa/coreTax_BEED_Aug.csv")

###### BEECH
#filter for taxa with relative abundance of atleast 0.0001 in 1 or more samples
beech_phydat_rel <- ps_norm_merged3 %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study == "BEECH") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  dplyr::summarise(avg_abundance_beech = round((mean(Abundance)*100),5),
                   prevalence_beech = round((sum(Abundance > 0) / dplyr::n())*100,5),
                   SD_val_beech = round((sd(Abundance, na.rm = T))*100,5)) %>%
  filter(prevalence_beech >= 80) #17 core genera
write.csv(beech_phydat_rel, "Output/CoreTaxa/coreTax_BEECH_Aug.csv")

########### Use full dataset for BEED and SEEM
###### SEEM
#filter for taxa with relative abundance of atleast 0.0001 in 1 or more samples
seem_phydat_full <- ps_norm_merged3 %>%
  ps_filter(study == "SEEM") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  dplyr::summarise(avg_abundance_beech = round((mean(Abundance)*100),5),
                   prevalence_beech = round((sum(Abundance > 0) / dplyr::n())*100,5),
                   SD_val_beech = round((sd(Abundance, na.rm = T))*100,5)) %>%
  filter(prevalence_beech >= 80) #16 core ASVs
write.csv(seem_phydat_full, "Output/CoreTaxa/coreTax_SEEMFull_Aug.csv")

###### BEED
#filter for taxa with relative abundance of atleast 0.0001 in 1 or more samples
beed_phydat_full <- ps_norm_merged3 %>%
  ps_filter(study == "BEED") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  dplyr::summarise(avg_abundance_beech = round((mean(Abundance)*100),5),
                   prevalence_beech = round((sum(Abundance > 0) / dplyr::n())*100,5),
                   SD_val_beech = round((sd(Abundance, na.rm = T))*100,5)) %>%
  filter(prevalence_beech >= 80) #15 core ASVs
write.csv(beed_phydat_full, "Output/CoreTaxa/coreTax_BEEDFull_Aug.csv")

##### Check which ASVs are common to all groups and are part of the core genera common to all groups

###### SEEM
#Filter for those with prevalence >= 80%
seem_phydat_asv <- ps_norm_merged3 %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study == "SEEM") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "ASVReadable") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  dplyr::summarise(avg_abundance_seem = round((mean(Abundance)*100),5),
                   prevalence_seem = round((sum(Abundance > 0) / dplyr::n())*100,5),
                   SD_val_seem = round((sd(Abundance, na.rm = T))*100,5)) %>%
  filter(prevalence_seem >= 80) #16 core genera
write.csv(seem_phydat_asv, "Output/CoreTaxa/coreASV_SEEM_Aug.csv")

###### BEED
#filter for taxa with relative abundance of atleast 0.0001 in 1 or more samples
beed_phydat_asv <- ps_norm_merged3 %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study == "BEED") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "ASVReadable") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  dplyr::summarise(avg_abundance_beed = round((mean(Abundance)*100),5),
                   prevalence_beed = round((sum(Abundance > 0) / dplyr::n())*100,5),
                   SD_val_beed = round((sd(Abundance, na.rm = T))*100,5)) %>%
  filter(prevalence_beed >= 80) #19 core genera
write.csv(beed_phydat_asv, "Output/CoreTaxa/coreASV_BEED_Aug.csv")

###### BEECH
#filter for taxa with relative abundance of atleast 0.0001 in 1 or more samples
beech_phydat_asv <- ps_norm_merged3 %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study == "BEECH") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "ASVReadable") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  dplyr::summarise(avg_abundance_beech = round((mean(Abundance)*100),5),
                   prevalence_beech = round((sum(Abundance > 0) / dplyr::n())*100,5),
                   SD_val_beech = round((sd(Abundance, na.rm = T))*100,5)) %>%
  filter(prevalence_beech >= 80) #17 core genera
write.csv(beech_phydat_asv, "Output/CoreTaxa/coreASV_BEECH_Aug.csv")

