########## All Malnourished Duodenum Analysis ##################
#Zambia vs Bangladesh vs Pakistan vs Central African Republic vs Madagascar
rm(list=ls())
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
load("Data/RData/AllDatasets_RelAbundNoSpikeNoDuplicates.RData") #ps0.ra_unique -- AFRIBIOTA, BEED, BEECH, ME, SEEM

############################################
######### Summary Statistics ###############
############################################

#get summary statistics - all
met_dat <- ps0.ra_unique %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  samdat_tbl()
merged_meta <- met_dat[!duplicated(met_dat$PID), ]

tab1vars <- c('WAZ', 'WLZ', 'LAZ', 'Age', 'Sex', 'LAZ_cat', 'WAZ_cat', 'WLZ_cat')
tab1 <- tableone::CreateTableOne(data=merged_meta, vars=tab1vars, strata= c("Country"))
tab1 <- as.data.frame(print(tab1, nonnormal=tab1vars, contDigits=1))
tab1 <- tibble::rownames_to_column(tab1, "row_names")
write.csv(tab1, "~/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis/Output/DescriptiveDat/DescriptiveTable_anthro_all.csv")

#get summary statistics - paired dataset only
met_dat2 <- merged_meta %>%
  filter(study != "AFRIBIOTA")

tab1vars <- c('WAZ', 'WLZ', 'LAZ', 'Age', 'Sex', 'LAZ_cat', 'WAZ_cat', 'WLZ_cat')
tab2 <- tableone::CreateTableOne(data=met_dat2, vars=tab1vars, strata= c("Country"))
tab2 <- as.data.frame(print(tab2, nonnormal=tab1vars, contDigits=1))
tab2 <- tibble::rownames_to_column(tab2, "row_names")
write.csv(tab2, "~/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis/Output/DescriptiveDat/DescriptiveTable_anthro_paired.csv")

colors_study = c("BEECH" = "#1b9e77", 
                 "BEED" = "#d95f02",
                 "AFRIBIOTA" = "#7570b3",
                 "SEEM" = "#66a61e")
comparisons2 <- list(c("BEECH","BEED"), c("BEECH","SEEM"), c("BEED","SEEM"),
                     c("AFRIBIOTA","BEED"), c("AFRIBIOTA","SEEM"), c("AFRIBIOTA","SEEM"))
comparisons1 <- list(c("BEECH","BEED"), c("BEECH","SEEM"), c("BEED","SEEM"))

#Age boxplots
a <- merged_meta %>%
  mutate(Country = ifelse(Country == "Central African Republic", "CAR", Country)) %>%
  mutate(study2 = ifelse(study == "BEECH children", "BEECH", study)) %>%
  ggplot(aes(x = study2, y = Age, colour = study2)) +
  geom_boxplot() +
  scale_color_manual(values = colors_study) +
  labs(y = "Age (Months)",
       x = "STUDY") +
  theme_bw() +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(comparisons = comparisons2, label = "p.value")

#Anthropometry - categorical
b <- merged_meta %>%
  pivot_longer(c(LAZ_cat, WAZ_cat, WLZ_cat), names_to = "anthro", values_to = "Category") %>%
  mutate(anthro2 = ifelse(anthro == "LAZ_cat", "LAZ", anthro)) %>%
  mutate(anthro2 = ifelse(anthro == "WAZ_cat", "WAZ", anthro2)) %>%
  mutate(anthro2 = ifelse(anthro == "WLZ_cat", "WLZ", anthro2)) %>%
  mutate(Cat2 = ifelse((stringr::str_detect(Category, "Severe")), "Severe", Category)) %>%
  mutate(Cat2 = ifelse((stringr::str_detect(Category, "Moderate")), "Moderate", Cat2)) %>%
  mutate(Cat2 = ifelse((stringr::str_detect(Category, "Not")), "No", Cat2)) %>%
  mutate(Cat2 = ifelse(Category == "NA", NA, Cat2)) %>%
  filter(!is.na(Cat2)) %>%
  ggplot(aes(x = study, fill = Cat2)) +
  geom_bar(position = "dodge") +
  labs(x = "STUDY",
       y = "Count",
       fill = "Anthropometry Category") +
  facet_grid(~anthro2, scales = "free_x") +
  theme_bw() +
  theme(legend.position = "bottom")

#Anthropometry - continuous
c <- merged_meta %>%
  pivot_longer(c(LAZ, WAZ, WLZ), names_to = "anthro", values_to = "Category") %>%
  mutate(study2 = ifelse(study == "BEECH children", "BEECH", study)) %>%
  filter(!is.na(Category)) %>%
  ggplot(aes(x = study2, y = Category, colour = study2)) +
  geom_boxplot() +
  scale_color_manual(values = colors_study) +
  labs(y = "Anthropometry Scores",
       x = "STUDY") +
  theme_bw() +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(comparisons = comparisons1, label = "p.value") +
  facet_grid(~anthro, scales = "free_x") 

#Sex
d <- merged_meta %>%
  mutate(Country = ifelse(Country == "Central African Republic", "CAR", Country)) %>%
  ggplot(aes(x = study, fill = Sex)) +
  geom_bar(position = "dodge") +
  labs(y = "Count") +
  theme_bw() +
  theme(legend.position = "bottom") 

ggarrange(a,d,b,c, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"))
ggsave("~/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis/Output/DescriptiveDat/Boxplot_metadata2.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 15, height = 12, units = "in")

