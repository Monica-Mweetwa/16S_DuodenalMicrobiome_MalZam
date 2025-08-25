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
#Import data
#setwd("C:/Users/Monica/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis") #Windows
setwd("~/Dropbox/PhD_Data/BEECH/BEECH_16s_Analysis/Meta-Analysis") #Mac
load("Data/RData/AllDatasets_RelAbundNoSpikeNoDuplicates.RData") #ps0.ra_unique -- AFRIBIOTA, BEED, BEECH, ME, SEEM

##############################################################
############# General comparisons with AFRIBIOTA  ############
##############################################################

#1 . Relative abundance
rel_plot1 <- ps0.ra_unique %>%
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
ggsave("Output/AFR_Comparisons/RelAbund_gen_paired.png", plot = rel_plot1, dpi = 600, device = "png", 
       width = 15, height = 6, units = "in")

#2. Alpha Diversity
ad = estimate_richness(ps0.ra_unique, measures = "Shannon")
ad <- rownames_to_column(ad, "PID2")
add <- data.frame(sample_data(ps0.ra_unique))
add <- rownames_to_column(add, "PID2")
add2 = merge(add, ad, by = "PID2")
comparisons1 <- list(c("BEECH","BEED"), c("BEECH","SEEM"), c("BEED","SEEM"))
comparisons2 <- list(c("BEECH","BEED"), c("BEECH","SEEM"), c("BEED","SEEM"),
                     c("AFRIBIOTA","BEED"), c("AFRIBIOTA","SEEM"), c("AFRIBIOTA","SEEM"))
colors_study = c("BEECH" = "#1b9e77", 
                 "BEED" = "#d95f02",
                 "AFRIBIOTA" = "#ba4aa4",
                 "SEEM" = "#66a61e")

alpha_plot <- add2 %>%
  filter(study != "ME") %>% #remove samples from children with SAM
  filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ggplot(aes(x = study, y = Shannon, fill = study)) +
  scale_fill_manual(values = colors_study) +
  geom_boxplot() + theme_bw() +
  labs(x = " ",
       y = "Shannon's index",
       fill = "STUDY") +
  ggpubr::stat_compare_means(comparisons = comparisons2, label = "p.value")
ggsave("Output/AFR_Comparisons/AlphaDiv.png", plot = alpha_plot, dpi = 600, device = "png", 
       width = 6, height = 6, units = "in")

#3. Beta diversity
#export data to input in PRIMER7 for ERMANOVA analysis
ps0.gen <- ps0.ra_unique %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study != "ME") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus")
out_tab <- as.data.frame(ps0.gen@otu_table)
sam_tab <- as.data.frame(ps0.gen@sam_data)
write.csv(out_tab, "Data/PRIMER7/otu_tab_rel_meta-analysis.csv")
write_csv(sam_tab, "Data/PRIMER7/sam_tab_rel_meta-analysis.csv")

beta_plot <- ps0.ra_unique %>%
  ps_filter(study != "ME") %>% #remove samples from children with SAM
  ps_filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  ps_filter(LAZ_cat != "Not stunted") %>% #remove children that are not stunted
  ps_filter(study != "ME") %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  dist_calc("bray") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "study", plot_taxa = 5, size = 8) +
  theme_classic() +
  stat_ellipse(aes(colour = study)) +
  ggside::geom_xsidedensity(aes(fill = study), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = study), alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = colors_study) +
  scale_color_manual(values = colors_study) +
  labs(color = "STUDY",
       x = "PC1 (53.3%)",
       y = "PC2 (7.3%)") +
  ggside::theme_ggside_void() +
  ggtext::geom_richtext(aes(x = 1, y = 2, label = "PERMANOVA Test: Pseudo F = 33.4, *p* = 0.001 <br> PERMDISP Test: F = 3.78, *p* = 0.042"),
                        size = 3, fill=NA, label.color = NA) 
ggsave("Output/AFR_Comparisons/beta_plot.png", plot = beta_plot, dpi = 600, device = "png", 
       width = 6, height = 6, units = "in")

ggpubr::ggarrange(ggarrange(alpha_plot, beta_plot, nrow = 1, labels = c("(a)", "(b)"), common.legend = T, legend = "left"), 
                  rel_plot1, ncol = 1, labels = c("", "(c)"))

ggsave("Output/AFR_Comparisons/combined_Figs_all2025.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 14, height = 10, units = "in")

#differential abundance
#ps0.ra_unique_age <- ps0.ra_unique %>%
#  ps_mutate(Age2 = ifelse(study == "AFRIBIOTA", "older", "younger"))
#output = ancombc2(data = ps0.ra_unique_age, tax_level = "Genus",
#                  fix_formula = "Age2", rand_formula = NULL,
#                  p_adj_method = "BY", pseudo_sens = TRUE,
#                  prv_cut = 0.1, lib_cut = 0, s0_perc = 0.05,
#                  group = "Age2", struc_zero = TRUE, neg_lb = TRUE)
#res_prim = output$res
#res_prim.maltype <- res_prim %>% select(taxon, contains("mal.type"))

