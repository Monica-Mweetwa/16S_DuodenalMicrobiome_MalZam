########### All Malnourished Duodenum Analysis ##################
################## Differentail abundance analysis ####################
library(microViz)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(vegan)
library(tidyr)
library(broom)
#Import data
#setwd("C:/Users/Monica/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis") #Windows
setwd("~/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis") #Mac
load("Data/RData/BEED_SEEM_BEECH_AbsASV.RData") #ps_SB_norm.ra.scaled

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
ps_norm_merged4 <-  ps_norm_merged3 %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") #remove children that are not stunted

colors_study = c("BEECH" = "#1b9e77", 
                 "BEED" = "#d95f02",
                 "AFRIBIOTA" = "#ba4aa4",
                 "SEEM" = "#66a61e")
colors_country = c("Zambia" = "#1b9e77", 
                   "Bangladesh" = "#d95f02",
                   "CAR" = "#7570b3",
                   "Madagascar" = "#e7298a",
                   "Pakistan" = "#66a61e")

###############################################################
######## Differential abundance analyses - ANCOM - BC2 ##########
#################################################################
library(ANCOMBC)
library(data.table)
############### Associations between core taxa abundance and anthro, Age and Sex #################
##########################################################################################
ancombc_func_cont <- function(dat, vars, main){
  core_tax <- c("Abiotrophia", "Actinomyces","Alloprevotella","Atopobium",
                "Bergeyella","Burkholderia-Caballeronia-Paraburkholderia",
                "Comamonadaceae Family", "Corynebacterium","Dolosigranulum","Fusobacterium",
                "Gemella", "Granulicatella","Haemophilus", "Johnsonella",
                "Lachnoanaerobaculum","Lachnospiraceae Family","Leptotrichia", "Neisseria",
                "Porphyromonas","Prevotella","Prevotella_7","Pseudomonas",
                "Rothia", "Streptobacillus","Streptococcus","Veillonella")
  output = ancombc2(data = dat, tax_level = "Genus",
                    fix_formula = vars, rand_formula = NULL,
                    p_adj_method = "BY", pseudo_sens = TRUE,
                    prv_cut = 0.1, lib_cut = 0, s0_perc = 0.05,
                    group = NULL, struc_zero = F)
  mod <- output$res %>% select(taxon, contains(main))
  mod_core <- mod %>% filter(taxon %in% core_tax)
}

lm_age <- ancombc_func_cont(vars = "Age",main = "Age", dat = ps_norm_merged4)
write.csv(lm_age, "Output/RelAbund_analysis/lm_age_univar_aug.csv")
lm_laz <- ancombc_func_cont(vars = "LAZ",main = "LAZ",dat = ps_norm_merged4)
write.csv(lm_laz, "Output/RelAbund_analysis/lm_laz_univar_aug.csv")
lm_waz <- ancombc_func_cont(vars = "WAZ",main = "WAZ",dat = ps_norm_merged4)
write.csv(lm_waz, "Output/RelAbund_analysis/lm_waz_univar_aug.csv")
lm_wlz <- ancombc_func_cont(vars = "WLZ",main = "WLZ",dat = ps_norm_merged4)
write.csv(lm_wlz, "Output/RelAbund_analysis/lm_wlz_univar_aug.csv")

lm_age <- ancombc_func_cont(vars = "Age + study + Sex",main = "Age", dat = ps_norm_merged4)
write.csv(lm_age, "Output/RelAbund_analysis/lm_age_mult_aug.csv")
lm_laz <- ancombc_func_cont(vars = "LAZ + Sex + study + Age",main = "LAZ",dat = ps_norm_merged4)
write.csv(lm_laz, "Output/RelAbund_analysis/lm_laz_mult_aug.csv")
lm_waz <- ancombc_func_cont(vars = "WAZ + Sex + study + Age",main = "WAZ",dat = ps_norm_merged4)
write.csv(lm_waz, "Output/RelAbund_analysis/lm_waz_mult_aug.csv")
lm_wlz <- ancombc_func_cont(vars = "WLZ + Sex + study + Age",main = "WLZ",dat = ps_norm_merged4)
write.csv(lm_wlz, "Output/RelAbund_analysis/lm_wlz_mult_aug.csv")

##Analysis of each study alone
#SEEM
lm_age <- ancombc_func_cont(vars = "Age",main = "Age", dat = ps_norm_merged4%>% ps_filter(study == "SEEM"))
write.csv(lm_age, "Output/RelAbund_analysis/IndividualStudies/lm_age_univar_seem.csv")
lm_laz <- ancombc_func_cont(vars = "LAZ",main = "LAZ",dat = ps_norm_merged4%>% ps_filter(study == "SEEM"))
write.csv(lm_laz, "Output/RelAbund_analysis/IndividualStudies/lm_laz_univar_seem.csv")
lm_waz <- ancombc_func_cont(vars = "WAZ",main = "WAZ",dat = ps_norm_merged4%>% ps_filter(study == "SEEM"))
write.csv(lm_waz, "Output/RelAbund_analysis/IndividualStudies/lm_waz_univar_seem.csv")
lm_wlz <- ancombc_func_cont(vars = "WLZ",main = "WLZ",dat = ps_norm_merged4%>% ps_filter(study == "SEEM"))
write.csv(lm_wlz, "Output/RelAbund_analysis/IndividualStudies/lm_wlz_univar_seem.csv")

lm_age <- ancombc_func_cont(vars = "Age + Sex",main = "Age", dat = ps_norm_merged4 %>% ps_filter(study == "SEEM"))
write.csv(lm_age, "Output/RelAbund_analysis/IndividualStudies/lm_age_mult_seem.csv")
lm_laz <- ancombc_func_cont(vars = "LAZ + Sex + Age",main = "LAZ",dat = ps_norm_merged4 %>% ps_filter(study == "SEEM"))
write.csv(lm_laz, "Output/RelAbund_analysis/IndividualStudies/lm_laz_mult_seem.csv")
lm_waz <- ancombc_func_cont(vars = "WAZ + Sex + Age",main = "WAZ",dat = ps_norm_merged4 %>% ps_filter(study == "SEEM"))
write.csv(lm_waz, "Output/RelAbund_analysis/IndividualStudies/lm_waz_mult_seem.csv")
lm_wlz <- ancombc_func_cont(vars = "WLZ + Sex + Age",main = "WLZ",dat = ps_norm_merged4 %>% ps_filter(study == "SEEM"))
write.csv(lm_wlz, "Output/RelAbund_analysis/IndividualStudies/lm_wlz_mult_seem.csv")

#BEED
lm_age <- ancombc_func_cont(vars = "Age",main = "Age", dat = ps_norm_merged4%>% ps_filter(study == "BEED"))
write.csv(lm_age, "Output/RelAbund_analysis/IndividualStudies/lm_age_univar_beed.csv")
lm_laz <- ancombc_func_cont(vars = "LAZ",main = "LAZ",dat = ps_norm_merged4%>% ps_filter(study == "BEED"))
write.csv(lm_laz, "Output/RelAbund_analysis/IndividualStudies/lm_laz_univar_beed.csv")
lm_waz <- ancombc_func_cont(vars = "WAZ",main = "WAZ",dat = ps_norm_merged4%>% ps_filter(study == "BEED"))
write.csv(lm_waz, "Output/RelAbund_analysis/IndividualStudies/lm_waz_univar_beed.csv")
lm_wlz <- ancombc_func_cont(vars = "WLZ",main = "WLZ",dat = ps_norm_merged4%>% ps_filter(study == "BEED"))
write.csv(lm_wlz, "Output/RelAbund_analysis/IndividualStudies/lm_wlz_univar_beed.csv")

lm_age <- ancombc_func_cont(vars = "Age + Sex",main = "Age", dat = ps_norm_merged4 %>% ps_filter(study == "BEED"))
write.csv(lm_age, "Output/RelAbund_analysis/IndividualStudies/lm_age_mult_beed.csv")
lm_laz <- ancombc_func_cont(vars = "LAZ + Sex + Age",main = "LAZ",dat = ps_norm_merged4 %>% ps_filter(study == "BEED"))
write.csv(lm_laz, "Output/RelAbund_analysis/IndividualStudies/lm_laz_mult_beed.csv")
lm_waz <- ancombc_func_cont(vars = "WAZ + Sex + Age",main = "WAZ",dat = ps_norm_merged4 %>% ps_filter(study == "BEED"))
write.csv(lm_waz, "Output/RelAbund_analysis/IndividualStudies/lm_waz_mult_beed.csv")
lm_wlz <- ancombc_func_cont(vars = "WLZ + Sex + Age",main = "WLZ",dat = ps_norm_merged4 %>% ps_filter(study == "BEED"))
write.csv(lm_wlz, "Output/RelAbund_analysis/IndividualStudies/lm_wlz_mult_beed.csv")

#BEECH
lm_age <- ancombc_func_cont(vars = "Age",main = "Age", dat = ps_norm_merged4%>% ps_filter(study == "BEECH"))
write.csv(lm_age, "Output/RelAbund_analysis/IndividualStudies/lm_age_univar_beech.csv")
lm_laz <- ancombc_func_cont(vars = "LAZ",main = "LAZ",dat = ps_norm_merged4%>% ps_filter(study == "BEECH"))
write.csv(lm_laz, "Output/RelAbund_analysis/IndividualStudies/lm_laz_univar_beech.csv")
lm_waz <- ancombc_func_cont(vars = "WAZ",main = "WAZ",dat = ps_norm_merged4%>% ps_filter(study == "BEECH"))
write.csv(lm_waz, "Output/RelAbund_analysis/IndividualStudies/lm_waz_univar_beech.csv")
lm_wlz <- ancombc_func_cont(vars = "WLZ",main = "WLZ",dat = ps_norm_merged4%>% ps_filter(study == "BEECH"))
write.csv(lm_wlz, "Output/RelAbund_analysis/IndividualStudies/lm_wlz_univar_beech.csv")

lm_age <- ancombc_func_cont(vars = "Age + Sex",main = "Age", dat = ps_norm_merged4 %>% ps_filter(study == "BEECH"))
write.csv(lm_age, "Output/RelAbund_analysis/IndividualStudies/lm_age_mult_beech.csv")
lm_laz <- ancombc_func_cont(vars = "LAZ + Sex + Age",main = "LAZ",dat = ps_norm_merged4 %>% ps_filter(study == "BEECH"))
write.csv(lm_laz, "Output/RelAbund_analysis/IndividualStudies/lm_laz_mult_beech.csv")
lm_waz <- ancombc_func_cont(vars = "WAZ + Sex + Age",main = "WAZ",dat = ps_norm_merged4 %>% ps_filter(study == "BEECH"))
write.csv(lm_waz, "Output/RelAbund_analysis/IndividualStudies/lm_waz_mult_beech.csv")
lm_wlz <- ancombc_func_cont(vars = "WLZ + Sex + Age",main = "WLZ",dat = ps_norm_merged4 %>% ps_filter(study == "BEECH"))
write.csv(lm_wlz, "Output/RelAbund_analysis/IndividualStudies/lm_wlz_mult_beech.csv")


## BEED - univar fulldataset - pearson
core_tax <- c("Abiotrophia", "Actinomyces","Alloprevotella","Atopobium",
              "Bergeyella","Burkholderia-Caballeronia-Paraburkholderia",
              "Comamonadaceae Family", "Corynebacterium","Dolosigranulum","Fusobacterium",
              "Gemella", "Granulicatella","Haemophilus", "Johnsonella",
              "Lachnoanaerobaculum","Lachnospiraceae Family","Leptotrichia", "Neisseria",
              "Porphyromonas","Prevotella","Prevotella_7","Pseudomonas",
              "Rothia", "Streptobacillus","Streptococcus","Veillonella")

pearson_corr_grouped <- function(var){
  grouped_dat<- dat_for_cor_beed %>%
    group_by(OTU) %>%
    group_split() #creates list of dataframes
  results <- lapply(grouped_dat, function(g_dat) {
    cor_1 <- psych::corr.test(g_dat %>% select(var, LAZ, WLZ, WAZ),
                              method = "pearson", adjust = "fdr")
    cr1 <- cor_1$ci %>%
      rownames_to_column(var = "rel")
    corr_dat <- cor_1$ci2 %>% mutate(var_name = var)
    corr_dat$var_name <- cr1$rel
    corr_dat$group <- unique(g_dat$OTU)
    return(corr_dat)
  })
  combined_results <- bind_rows(results)
  mod_core <- combined_results %>% filter(group %in% core_tax)
  print(mod_core)
}
#Convert abundance to log10 values
dat_for_cor_beed <- ps_norm_merged3 %>% 
  tax_fix() %>%
  tax_transform("identity", rank = "Genus") %>%
  ps_melt() %>%
  filter(study == "BEED") %>%
  mutate(Abundance = ifelse(Abundance == 0, 0.01, Abundance)) %>%
  mutate(Abundance_log10 = log10(Abundance)) #Convert scaled absolute ASV abundance to log10 scale
cor_1_beed <- pearson_corr_grouped(var = "Abundance_log10")
corr_abund_beed <- cor_1_beed %>%
  filter(var_name %like% "Ab_10")
write.csv(corr_abund_beed, "Output/RelAbund_analysis/IndividualStudies/FullDat_Pearson_BEED.csv")

## SEEM - univar fulldataset - pearson
pearson_corr_grouped <- function(var){
  grouped_dat<- dat_for_cor_seem %>%
    group_by(OTU) %>%
    group_split() #creates list of dataframes
  results <- lapply(grouped_dat, function(g_dat) {
    cor_1 <- psych::corr.test(g_dat %>% select(var, LAZ, WLZ, WAZ),
                              method = "pearson", adjust = "fdr")
    cr1 <- cor_1$ci %>%
      rownames_to_column(var = "rel")
    corr_dat <- cor_1$ci2 %>% mutate(var_name = var)
    corr_dat$var_name <- cr1$rel
    corr_dat$group <- unique(g_dat$OTU)
    return(corr_dat)
  })
  combined_results <- bind_rows(results)
  mod_core <- combined_results %>% filter(group %in% core_tax)
  print(mod_core)
}

dat_for_cor_seem <- ps_norm_merged3 %>% 
  tax_fix() %>%
  tax_transform("identity", rank = "Genus") %>%
  ps_melt() %>%
  filter(study == "SEEM") %>%
  mutate(Abundance = ifelse(Abundance == 0, 0.01, Abundance)) %>%
  mutate(Abundance_log10 = log10(Abundance)) #Convert scaled absolute ASV abundance to log10 scale
cor_1_seem <- pearson_corr_grouped(var = "Abundance_log10")
corr_abund_seem <- cor_1_seem %>%
  filter(var_name %like% "Ab_10")
write.csv(corr_abund_seem, "Output/RelAbund_analysis/IndividualStudies/FullDat_Pearson_SEEM.csv")

