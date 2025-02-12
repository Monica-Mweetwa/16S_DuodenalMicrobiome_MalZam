########### All Malnourished Duodenum Analysis ##################
#Zambia vs Bangladesh vs Pakistan vs Central African Republic vs Madagascar

library(microViz)
library(tidyverse)
library(phyloseq)
#Import data
setwd("C:/Users/Monica/Dropbox/BEECH_16s_Analysis/Meta-Analysis")
load("Data/RData/StuntingCombined_unique.RData")

colors = c("Zambia" = "#1b9e77", 
           "Bangladesh" = "#d95f02",
           "CAR" = "#7570b3",
           "Madagascar" = "#e7298a",
           "Pakistan" = "#66a61e")

#Get log2 linear model taxa vs study 
lm_mod_func = function(var,data,rank){
  require(microViz)
  mod_lm <- data %>%
    tax_fix() %>%
    tax_transform("compositional", rank = rank) %>%
    taxatree_models(type = "lm", 
                    rank = rank,
                    trans = "log2", 
                    trans_args = list(zero_replace = "halfmin"),
                    variables = c(var))
  lm_stats <- taxatree_models2stats(mod_lm)
  lm_stats <- taxatree_stats_p_adjust(data = lm_stats, method = "BY", grouping = "rank")
  lm_stats2 <- lm_stats %>% taxatree_stats_get()
  #Calculate mean abundance and prevalence for each feature
  prop_prev_data <- data %>%
    tax_fix() %>%
    tax_transform("compositional", rank = rank) %>%
    ps_melt() %>%
    group_by(OTU) %>%
    summarise(avg_abundance = mean(Abundance),
              prevalence = sum(Abundance > 0) / dplyr::n()) %>%
    rename("taxon"="OTU")
  lm_stat3 <- left_join(lm_stats2, prop_prev_data, by = "taxon")
  print(lm_stat3)
}

#1. Afribiota vs ALL
rf_dat_laz_unique2 <- rf_dat_laz_unique %>%
  ps_mutate(Age2 = ifelse(study == "Afribiota", "> 2.5 yrs", "< 2.5 yrs"))

lm_models <- lm_mod_func(var = "Age2",
                         data = rf_dat_laz_unique2,
                         rank = "Genus")
write.csv(lm_models, "Output/ST_Combined/lm_sig_mods_catAge.csv")
sig_models <- lm_models %>%
  filter(p.adj.BY.rank < 0.05) %>%
  filter(prevalence > 0.1) %>%
  mutate(study = str_remove(term, "^study"))
sig_taxa <- sig_models %>% select(taxon) %>% distinct()
colors = c("> 2.5 yrs" = "orange", 
           "< 2.5 yrs" = "pink")
#plot relative abundance
plot_data <- rf_dat_laz_unique2 %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt
plot_data2 <- plot_data %>%
  filter(OTU %in% sig_taxa$taxon) %>%
  filter(!is.na(Abundance)) %>%
  filter(Abundance != 0)

ggplot(plot_data2, aes(x=Age2, y=Abundance, colour=Age2)) + 
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = colors, name = "Age") +
  scale_y_log10() +
  theme_classic() + 
  facet_wrap(~OTU, scales = "free_y", ncol = 11) +
  labs(y = "Relative Abundance", x = " ")
ggsave("Output/ST_Combined/lm_sig_mods_age.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 29, height = 8, units = "in")

#2. STUDY
rf_dat_laz_unique3 <- rf_dat_laz_unique %>%
  ps_filter(study != "Afribiota") %>%
  ps_mutate(study2 = ifelse(study == "BEECH children", "BEECH", study))

lm_models <- lm_mod_func(var = "study2",
                         data = rf_dat_laz_unique3,
                         rank = "Genus")
write.csv(lm_models, "Output/ST_Combined/lm_sig_mods_paired.csv")
sig_models <- lm_models %>%
  filter(p.adj.BY.rank < 0.05) %>%
  filter(prevalence > 0.1) %>%
  mutate(study = str_remove(term, "^study"))
sig_taxa <- sig_models %>% select(taxon) %>% distinct()

colors = c("BEECH" = "#1b9e77", 
           "BEED" = "#d95f02",
           "SEEM" = "#66a61e")
#plot relative abundance
plot_data <- rf_dat_laz_unique3 %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt
plot_data2 <- plot_data %>%
  filter(OTU %in% sig_taxa$taxon) %>%
  filter(!is.na(Abundance)) %>%
  filter(Abundance != 0)

ggplot(plot_data2, aes(x=study2, y=Abundance, colour=study2)) + 
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = colors, name = "Study") +
  scale_y_log10() +
  theme_classic() + 
  facet_wrap(~OTU, scales = "free_y", ncol = 7) +
  labs(y = "Relative Abundance", x = " ")
ggsave("Output/ST_Combined/lm_sig_mods_paired.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 26, height = 8, units = "in")
#plot absolute abundance
plot_data <- rf_dat_laz_unique3 %>%
  tax_fix() %>%
  tax_transform("identity", rank = "Genus") %>%
  ps_melt
plot_data2 <- plot_data %>%
  filter(OTU %in% sig_taxa$taxon) %>%
  filter(!is.na(Abundance)) %>%
  filter(Abundance != 0) 

ggplot(plot_data2, aes(x=study2, y=Abundance, colour=study2)) + 
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = colors, name = "Study") +
  scale_y_log10() +
  theme_classic() + 
  facet_wrap(~OTU, scales = "free_y", ncol = 10) +
  labs(y = "Relative Abundance", x = " ")
ggsave("Output/ST_Combined/lm_sig_mods_abs_paired.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 26, height = 8, units = "in")

#PERMAONVA on the paired analysis
rf_dat_laz_unique_comp <- rf_dat_laz_unique %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_filter(study != "Afribiota")

#PERMANOVA 
library(vegan)

#permanova function
permanova_func = function(equation,data, var){
  dist = phyloseq::distance(data, method="bray")
  metadata <- data.frame(sample_data(data))
  require(vegan)
  results <- adonis2(equation, data = metadata, na.action = na.omit)
  coef <- results %>%
    rownames_to_column(var = "Term") %>%
    filter(Term != "Residual") %>%
    mutate(Term = ifelse(Term == "Model", var, paste0(Term))) %>%
    filter(Term != "Total")
  print(coef)
}

perm_age# <- permanova_func(equation = dist ~ Age, var = "Age", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_sex# <- permanova_func(equation = dist ~ Sex, var = "Sex", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_waz <-permanova_func(equation = dist ~ WAZ, var = "WAZ", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_whz# <- permanova_func(equation = dist ~ WLZ, var = "WLZ", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_laz <- permanova_func(equation = dist ~ LAZ, var = "LAZ", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
reg_stats <- Reduce(function(x,y) merge(x, y, all.x = TRUE, all.y = TRUE),
                    list(perm_country, perm_study, perm_age, perm_sex, perm_waz, perm_whz, perm_laz))
#multivar 

perm_waz <-permanova_func(equation = dist ~ WAZ + study, var = "WAZ", data = rf_dat_laz_unique_comp)
perm_laz <- permanova_func(equation = dist ~ LAZ + study, var = "LAZ", data = rf_dat_laz_unique_comp)


#3. Differential abundance with LAZ, WAZ & WLZ
rf_dat_laz_unique3 <- rf_dat_laz_unique %>%
  ps_filter(study != "Afribiota") %>%
  ps_mutate(study2 = ifelse(study == "BEECH children", "BEECH", study))
#LAZ
lm_models <- lm_mod_func(var = "LAZ",
                         data = rf_dat_laz_unique3,
                         rank = "Genus")
write.csv(lm_models, "Output/ST_Combined/lm_sig_mods_laz_paired.csv") # no significant associations
#WAZ
lm_models <- lm_mod_func(var = "WAZ",
                         data = rf_dat_laz_unique3,
                         rank = "Genus")
write.csv(lm_models, "Output/ST_Combined/lm_sig_mods_waz_paired.csv") # no significant associations
#WLZ
lm_models <- lm_mod_func(var = "WLZ",
                         data = rf_dat_laz_unique3,
                         rank = "Genus")
write.csv(lm_models, "Output/ST_Combined/lm_sig_mods_wlz_paired.csv") # no significant associations
