########### All Malnourished Duodenum Analysis ##################
#Zambia vs Bangladesh vs Pakistan vs Central African Republic vs Madagascar

library(microViz)
library(tidyverse)
library(phyloseq)
#Import data
setwd("C:/Users/Monica/Dropbox/BEECH_16s_Analysis/Meta-Analysis")
load("Data/RData/StuntingCombined_unique.RData")


#Make the Zambian group the reference point
rf_dat_laz_unique@sam_data[["Country"]] <- relevel(as.factor(rf_dat_laz_unique@sam_data[["Country"]]), ref="Zambia")


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
  mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country))) %>%
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
  mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
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

#Relative Abundance
rel_plot <- rf_dat_laz_unique %>%
  ps_mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
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
  ps_mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
  ps_filter(study != "Afribiota") %>%
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


#Alpha diversity 
ad = estimate_richness(rf_dat_laz_unique, measures = "Shannon")
ad <- rownames_to_column(ad, "PID")
add <- data.frame(sample_data(rf_dat_laz_unique))
add <- rownames_to_column(add, "PID")
add2 = merge(add, ad, by = "PID")
add2 <- add2 %>%
  mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))

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
alpha_plot1 <- ggplot(add2 %>% filter(study != "Afribiota"), aes(Country, Shannon, fill = Country)) + 
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
  stat_cor(method = "spearman", label.y = 3.5) +
  geom_smooth(method = "lm") +
  facet_grid(Country ~ clin_vars,scales="free",margin="mal.type",drop = TRUE,switch="both")

ggsave("Output/ST_Combined/Alpha_anthro2.png", plot = a, dpi = 600, device = "png", width = 10, height = 8, units = "in")


#Beta Diversity
rf_dat_laz_unique_comp <- rf_dat_laz_unique %>%
  tax_transform("compositional", rank = "Genus")

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
#all datasets combined
perm_country <- permanova_func(equation = dist ~ Country, var = "Country", data = rf_dat_laz_unique_comp)
perm_study <- permanova_func(equation = dist ~ study, var = "study", data = rf_dat_laz_unique_comp)
perm_age <- permanova_func(equation = dist ~ Age, var = "Age", data = rf_dat_laz_unique_comp)
perm_sex <- permanova_func(equation = dist ~ Sex, var = "Sex", data = rf_dat_laz_unique_comp)
perm_waz <-permanova_func(equation = dist ~ WAZ, var = "WAZ", data = rf_dat_laz_unique_comp)
perm_whz <- permanova_func(equation = dist ~ WLZ, var = "WLZ", data = rf_dat_laz_unique_comp)
perm_laz <- permanova_func(equation = dist ~ LAZ, var = "LAZ", data = rf_dat_laz_unique_comp)
reg_stats <- Reduce(function(x,y) merge(x, y, all.x = TRUE, all.y = TRUE),
                    list(perm_country, perm_study, perm_age, perm_sex, perm_waz, perm_whz, perm_laz))
write.csv(reg_stats, "Output/ST_Combined/permanova_alldatasets.csv")

#paired datasets only
perm_country <- permanova_func(equation = dist ~ Country, var = "Country", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_study <- permanova_func(equation = dist ~ study, var = "study", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_age <- permanova_func(equation = dist ~ Age, var = "Age", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_sex <- permanova_func(equation = dist ~ Sex, var = "Sex", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_waz <-permanova_func(equation = dist ~ WAZ, var = "WAZ", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_whz <- permanova_func(equation = dist ~ WLZ, var = "WLZ", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
perm_laz <- permanova_func(equation = dist ~ LAZ, var = "LAZ", data = rf_dat_laz_unique_comp %>% ps_filter(study != "Afribiota"))
reg_stats <- Reduce(function(x,y) merge(x, y, all.x = TRUE, all.y = TRUE),
                    list(perm_country, perm_study, perm_age, perm_sex, perm_waz, perm_whz, perm_laz))
write.csv(reg_stats, "Output/ST_Combined/permanova_paired.csv")

#plot PcoA coloured by country
#all datasets
beta_plot <- 
  rf_dat_laz_unique %>%
  ps_mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
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

beta_plot2 <- rf_dat_laz_unique_comp %>%
  ps_filter(study != "Afribiota") %>%
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


#Combined plots: relative abundance, alpha and beta diversity

#1. paired dataset only
library(ggpubr)
ggarrange(ggarrange(alpha_plot1, beta_plot2, nrow = 1, labels = c("A.", "B.")), rel_plot1,  
          common.legend = F, ncol = 1, labels = c("", "C.")) 
ggsave("Output/ST_Combined/combined_Figs_paired.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 14, height = 10, units = "in")

#1. all datasets retrieved
library(ggpubr)
ggarrange(ggarrange(alpha_plot, beta_plot, nrow = 1, labels = c("A.", "B.")), rel_plot,  
          common.legend = F, ncol = 1, labels = c("", "C.")) 
ggsave("Output/ST_Combined/combined_Figs_all.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 14, height = 10, units = "in")

