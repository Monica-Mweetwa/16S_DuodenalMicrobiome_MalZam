########### All Malnourished Duodenum Analysis ##################
#Zambia vs Bangladesh vs Pakistan vs Central African Republic vs Madagascar

library(microViz)
library(tidyverse)
library(phyloseq)
#Import data
#setwd("C:/Users/Monica/Dropbox/BEECH_16s_Analysis/Meta-Analysis") #Windows
setwd("~/Dropbox/BEECH_16s_Analysis/Meta-Analysis") #Mac
load("Data/RData/StuntingCombined_unique.RData")


#Make the Zambian group the reference point
rf_dat_laz_unique@sam_data[["Country"]] <- relevel(as.factor(rf_dat_laz_unique@sam_data[["Country"]]), ref="Zambia")
rf_dat_laz_unique <- rf_dat_laz_unique%>%
  ps_mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country)),
            study = ifelse(study == "Afribiota", "AFRIBIOTA", study),
            study = ifelse(study == "BEECH children", "BEECH", study))

############################################
######### Summary Statistics ###############
############################################

#get summary statistics - all
met_dat <- rf_dat_laz_unique %>%
  samdat_tbl()
tab1vars <- c('WAZ', 'WLZ', 'LAZ', 'Age')
tab1 <- tableone::CreateTableOne(data=met_dat, vars=tab1vars, strata= c("Country"))
tab1 <- as.data.frame(print(tab1, nonnormal=tab1vars, contDigits=1))
tab1 <- tibble::rownames_to_column(tab1, "row_names")
write.csv(tab1, "Output/ST_Combined/DescriptiveTable_anthro_all.csv")
#get summary statistics - paired dataset only
met_dat2 <- rf_dat_laz_unique %>%
  samdat_tbl() %>% filter(study != "Afribiota")
tab1vars <- c('WAZ', 'WLZ', 'LAZ', 'Age')
tab2 <- tableone::CreateTableOne(data=met_dat2, vars=tab1vars, strata= c("Country"))
tab2 <- as.data.frame(print(tab2, nonnormal=tab1vars, contDigits=1))
tab2 <- tibble::rownames_to_column(tab2, "row_names")
write.csv(tab2, "Output/ST_Combined/DescriptiveTable_anthro_paired.csv")

colors = c("Zambia" = "#1b9e77", 
           "Bangladesh" = "#d95f02",
           "CAR" = "#7570b3",
           "Madagascar" = "#e7298a",
           "Pakistan" = "#66a61e")
#Age boxplots
met_dat %>%
  ggplot(aes(x = Country, y = Age, colour = Country)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  labs(y = "Age (Months)") +
  theme_bw() +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format")
ggsave("Output/ST_Combined/Boxplot_age.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 6, height = 4, units = "in")

#Anthropometry boxplots
met_dat %>%
  pivot_longer(cols = c(WAZ, WLZ, LAZ), names_to = "anthro", values_to = "vals") %>%
  ggplot(aes(x = Country, y = vals, colour = Country)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  facet_grid(rows = "anthro", scales = "free") +
  labs(y = "Z-Scores") +
  theme_bw() +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format")
ggsave("Output/ST_Combined/Boxplot_anthro.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 6, height = 8, units = "in")

#############################################
############# Relative Abundance ############
#############################################
rel_plot <- rf_dat_laz_unique %>%
  comp_barplot(tax_level = "Genus",
               label = "sampleID2",
               n_taxa = 15) +
  #coord_flip() + # horizontal bars are often more readable
  facet_grid(
    cols = vars(Country),
    scales = "free_x", # this only frees y scale per row in grid faceting
    space = "free_x" # allows bars to be same size by freeing facet heights
  ) +
  theme(axis.text.x = (element_blank()),
        axis.ticks.x = (element_blank()))
ggsave("Output/ST_Combined/RelAbund_gen.png", plot = rel_plot, dpi = 600, device = "png", 
       width = 15, height = 6, units = "in")

#No Afribiota
rel_plot1 <- rf_dat_laz_unique %>%
  ps_filter(study != "AFRIBIOTA") %>%
  comp_barplot(tax_level = "Genus",
               label = "sampleID2",
               n_taxa = 15) +
  #coord_flip() + # horizontal bars are often more readable
  facet_grid(
    cols = vars(Country),
    scales = "free_x", # this only frees y scale per row in grid faceting
    space = "free_x" # allows bars to be same size by freeing facet heights
  ) +
  theme(axis.text.x = (element_blank()),
        axis.ticks.x = (element_blank()))
ggsave("Output/ST_Combined/RelAbund_gen_paired.png", plot = rel_plot1, dpi = 600, device = "png", 
       width = 15, height = 6, units = "in")

############ Identify core taxa in each cohort #################
################################################################

###### get prevalence of each taxa per group
prop_prev_data_1 <- rf_dat_laz_unique %>%
  ps_filter(study == "SEEM") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  summarise(avg_abundance = mean(Abundance),
            prevalence = sum(Abundance > 0) / dplyr::n()) %>%
  rename("taxon"="OTU") %>%
  filter(prevalence > 0.8)

prop_prev_data_2 <- rf_dat_laz_unique %>%
  ps_filter(study == "BEED") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  summarise(avg_abundance = mean(Abundance),
            prevalence = sum(Abundance > 0) / dplyr::n()) %>%
  rename("taxon"="OTU") %>%
  filter(prevalence > 0.8)

prop_prev_data_3 <- rf_dat_laz_unique %>%
  ps_filter(study == "BEECH") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt() %>%
  group_by(OTU) %>%
  summarise(avg_abundance = mean(Abundance),
            prevalence = sum(Abundance > 0) / dplyr::n()) %>%
  rename("taxon"="OTU") %>%
  filter(prevalence > 0.8)

#Venn diagram of core taxa
library(ggvenn)
library(RColorBrewer)
AA <- c("hi","foo", "bar","yep","woo","hoo")
BB <- c("baa","yep", "woo","yes")
CC <- c("yes","foo","hi","woo", "huh")

x <- list(SEEM=prop_prev_data_1$taxon , BEED=prop_prev_data_2$taxon , BEECH=prop_prev_data_3$taxon)

venn_plot <- ggvenn(x, show_elements = T, 
       label_sep = "\n", 
       fill_color = brewer.pal(name="Set2",n=3),
       text_size = 1)
ggsave("Output/ST_Combined/Core_genera_paired.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 12, height = 16, units = "in")

#############################################
############# Alpha diversity ###############
#############################################

ad = estimate_richness(rf_dat_laz_unique, measures = "Shannon")
ad <- rownames_to_column(ad, "PID")
add <- data.frame(sample_data(rf_dat_laz_unique))
add <- rownames_to_column(add, "PID")
add2 = merge(add, ad, by = "PID")

alpha_plot <- ggplot(add2, aes(Country, Shannon, fill = Country)) + 
  geom_boxplot( outlier.shape = NA, size = 0.8, width = 0.8, na.rm = TRUE) +
  theme_bw() + 
  scale_fill_manual(values = colors) +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format") +
  labs(x = " ", 
       y = "Shannon index")  +
  theme(axis.text.x = element_text(vjust = 0.5),
        plot.title = element_text(size = 25, hjust = 0.5),
        legend.position = "none")
ggsave("Output/ST_Combined/AlphaDivComp.png", plot = alpha_plot, dpi = 600, 
       device = "png", width = 7, height = 4, units = "in")
#No Afribiota
alpha_plot1 <- ggplot(add2 %>% filter(study != "AFRIBIOTA"), aes(Country, Shannon, fill = Country)) + 
  geom_boxplot( outlier.shape = NA, size = 0.8, width = 0.8, na.rm = TRUE) +
  theme_bw() + 
  scale_fill_manual(values = colors) +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format") +
  labs(x = " ", 
       y = "Shannon index")  +
  theme(axis.text.x = element_text(vjust = 0.5),
        plot.title = element_text(size = 25, hjust = 0.5),
        legend.position = "none")
ggsave("Output/ST_Combined/AlphaDivComp_paired.png", plot = alpha_plot1, dpi = 600, 
       device = "png", width = 7, height = 4, units = "in")

#Shannon scatterplot with anthropometry
library(ggpubr)
a <- add2 %>% 
  pivot_longer(names_to = "clin_vars", values_to = "vals", cols = c(Age, WAZ,
                                                                    WLZ, LAZ)) %>%
  filter(!is.na(vals)) %>%
  filter(vals != 0) %>%
  ggplot(aes(x = vals, y = Shannon, colour = Country)) +
  scale_color_manual(values = colors) + labs(color = "Country") +
  geom_point() + theme_bw() + theme(legend.position = "none") +
  labs(x = " ",
       y = "Shannon's index") +
  stat_cor(method = "spearman", label.y = 3.5, cor.coef.name = "rho") +
  geom_smooth(method = "lm") +
  facet_grid(Country ~ clin_vars,scales="free",margin="mal.type",drop = TRUE,switch="both")

ggsave("Output/ST_Combined/Alpha_anthro2.png", plot = a, dpi = 600, device = "png", width = 10, height = 8, units = "in")

######################################
########### Beta Diversity ###########
######################################

rf_dat_laz_unique_comp <- rf_dat_laz_unique %>%
  tax_transform("compositional", rank = "Genus") 
  
#PERMANOVA 
library(vegan)
permanova_func = function(equation,data, var){
  dist = phyloseq::distance(data, method="bray")
  metadata <- data.frame(sample_data(data))
  results <- adonis2(equation, data = metadata, na.action = na.omit)
  coef <- results %>%
    rownames_to_column(var = "Term") %>%
    filter(Term != "Residual") %>%
    mutate(Term = ifelse(Term == "Model", var, paste0(Term))) %>%
    filter(Term != "Total")
  print(coef)
}
#all datasets combined
perm_country <- permanova_func(equation = dist ~ Country, var = "Country", data = rf_dat_laz_unique_comp)
write_csv(perm_country, "Output/ST_Combined/perm_country_all.csv")
perm_country_paired <- permanova_func(equation = dist ~ Country, var = "Country", data = rf_dat_laz_unique_comp %>% ps_filter(study != "AFRIBIOTA"))
write_csv(perm_country_paired, "Output/ST_Combined/perm_country_paired.csv")

#PERMANOVA - pairwise
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
dist_all = phyloseq::distance(rf_dat_laz_unique_comp, method="bray")
metadata <- data.frame(sample_data(rf_dat_laz_unique_comp))
pairwise_perm <- pairwise.adonis(dist_all, phyloseq::sample_data(rf_dat_laz_unique_comp)$study)
write_csv(pairwise_perm, "Output/ST_Combined/perm_study_pairwise.csv")
pairwise_perm <- pairwise.adonis(dist_all, phyloseq::sample_data(rf_dat_laz_unique_comp)$Country)
write_csv(pairwise_perm, "Output/ST_Combined/perm_country_pairwise.csv")

#plot PcoA coloured by Country
#all datasets
beta_plot <- rf_dat_laz_unique %>%
  tax_transform("compositional", rank = "Genus") %>%
  dist_calc("bray") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Country", plot_taxa = 5, size = 8) +
  scale_color_manual(values = colors) +
  labs(color = "Country") +
  theme_classic() +
  annotate("text", x = 1, y = 2, size = 3,
           label = paste("PERMANOVA Test:\nR2 = 0.49\np = 0.001"))
ggsave("Output/ST_Combined/BetaPlot_Country.png", plot = beta_plot, dpi = 600, device = "png", width = 8, height = 6, units = "in")

#paired read datasets only
beta_plot2 <- rf_dat_laz_unique %>%
  ps_filter(study != "AFRIBIOTA") %>%
  tax_transform("compositional", rank = "Genus") %>%
  dist_calc("bray") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Country", plot_taxa = 1:5, size = 8) +
  stat_ellipse(aes(linetype = Country, colour = Country)) +
  scale_color_manual(values = colors) +
  theme_bw() +
  ggside::geom_xsidedensity(aes(fill = Country), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = Country), alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = colors) +
  ggside::theme_ggside_void()+
  annotate("text", x = -0.5, y = 1.3, size = 3,
           label = paste("PERMANOVA Test:\nR2 = 0.12\np = 0.001"))
ggsave("Output/ST_Combined/BetaPlot_paired.png", plot = beta_plot2, dpi = 600, device = "png", width = 8, height = 6, units = "in")


##############################################################
######## Log linear modelling ##########
##############################################################
#Get log2 linear model: taxa vs study 
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
    rename("OTU"="taxon")
  lm_stat3 <- left_join(lm_stats2, prop_prev_data, by = "taxon")
  print(lm_stat3)
}

#1. AFRIBIOTA vs ALL
rf_dat_laz_unique2 <- rf_dat_laz_unique %>%
  ps_mutate(Age2 = ifelse(study == "AFRIBIOTA", "> 2.5 yrs", "< 2.5 yrs"))

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

#2. SEEM vs BEED vs BEECH
rf_dat_laz_unique3 <- rf_dat_laz_unique %>%
  ps_filter(study != "AFRIBIOTA") %>%
  ps_mutate(study2 = ifelse(study == "BEECH children", "BEECH", study))

#a. BEECH vs BEED & BEECH vs SEEM
lm_models <- lm_mod_func(var = "study2",
                         data = rf_dat_laz_unique3,
                         rank = "Genus")
write.csv(lm_models, "Output/ST_Combined/lm_sig_mods_paired.csv")

sig_models <- lm_models %>%
  filter(p.adj.BY.rank < 0.05) %>%
  filter(prevalence > 0.1) %>%
  mutate(study = str_remove(term, "^study"))
sig_taxa_beed <- sig_models %>% filter(study == "2BEED") %>% select(taxon) %>% distinct()
sig_taxa_seem <- sig_models %>% filter(study == "2SEEM") %>% select(taxon) %>% distinct()

#b. SEEM vs BEED
lm_models <- lm_mod_func(var = "study2",
                         data = rf_dat_laz_unique3 %>% ps_filter(study2 != "BEECH"),
                         rank = "Genus")
write.csv(lm_models, "Output/ST_Combined/lm_sig_mods_paired_noBEECH.csv")
sig_models <- lm_models %>%
  filter(p.adj.BY.rank < 0.05) %>%
  filter(prevalence > 0.1) %>%
  mutate(study = str_remove(term, "^study"))
sig_taxa_beedseem <- sig_models %>% select(taxon) %>% distinct()

colors = c("BEECH" = "#1b9e77", 
           "BEED" = "#d95f02",
           "SEEM" = "#66a61e")

#plot relative abundance
plot_data <- rf_dat_laz_unique3 %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt
plot_data_beed <- plot_data %>%
  filter(OTU %in% sig_taxa_beed$taxon) %>%
  mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance)) #%>%
plot_data_seem <- plot_data %>%
  filter(OTU %in% sig_taxa_seem$taxon) %>%
  mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance)) 
plot_data_beedseem <- plot_data %>%
  filter(OTU %in% sig_taxa_beedseem$taxon) %>%
  mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance)) 

a <- ggplot(plot_data_beed %>% filter(study != "SEEM"), aes(x=study2, y=Abundance, colour=study2)) + 
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = colors, name = "Study") +
  scale_y_log10() +
  theme_classic() + 
  facet_wrap(~OTU, scales = "free_y", ncol = 8) +
  labs(y = "Relative Abundance", x = " ")

b <- ggplot(plot_data_seem %>% filter(study != "BEED"), aes(x=study2, y=Abundance, colour=study2)) + 
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = colors, name = "Study") +
  scale_y_log10() +
  theme_classic() + 
  facet_wrap(~OTU, scales = "free_y", ncol = 8) +
  labs(y = "Relative Abundance", x = " ")

c <- ggplot(plot_data_beedseem %>% filter(study2 != "BEECH"), aes(x=study2, y=Abundance, colour=study2)) + 
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = colors, name = "Study") +
  scale_y_log10() +
  theme_classic() + 
  facet_wrap(~OTU, scales = "free_y", ncol = 8) +
  labs(y = "Relative Abundance", x = " ")

ggarrange(a, c, b, ncol = 1, labels = c("A.", "B.", "C."), heights = c(2, 1,1))
ggsave("Output/ST_Combined/lm_sig_mods_paired2.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 24, height = 12, units = "in")


################################################################################
######## Combined plots: relative abundance, alpha and beta diversity ##########
################################################################################

#1. paired dataset only
library(ggpubr)
ggarrange(ggarrange(alpha_plot1, beta_plot2, venn_plot, nrow = 1, labels = c("A.", "B.", "C.")), rel_plot1,  
          common.legend = F, ncol = 1, labels = c("", "D.")) 
ggsave("Output/ST_Combined/combined_Figs_paired_Apr.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 18, height = 15, units = "in")

#1. all datasets retrieved
library(ggpubr)
ggarrange(ggarrange(alpha_plot, beta_plot, nrow = 1, labels = c("A.", "B.")), rel_plot,  
          common.legend = F, ncol = 1, labels = c("", "C.")) 
ggsave("Output/ST_Combined/combined_Figs_all.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 14, height = 10, units = "in")


############## Done!
