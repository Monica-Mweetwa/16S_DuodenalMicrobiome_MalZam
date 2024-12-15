########### All Malnourished Duodenum Analysis ##################
#Zambia vs Bangladesh vs Pakistan vs Central African Republic vs Madagascar

#Import data
load("Data/RData/AllMAlStudies.RData")

library(microViz)
library(tidyverse)
library(phyloseq)
rf_dat@sam_data[["Country"]] <- relevel(as.factor(rf_dat@sam_data[["Country"]]), ref="Zambia")

#Look at clinical features of the samples 
met_dat <- rf_dat %>%
  samdat_tbl()
#get summary statistics
tab1vars <- c('WAZ', 'WLZ', 'LAZ', 'Age')
tab1 <- tableone::CreateTableOne(data=met_dat %>% filter(study != "EE children"), vars=tab1vars, strata= c("Country"))
tab1 <- as.data.frame(print(tab1, nonnormal=tab1vars, contDigits=1))
tab1 <- tibble::rownames_to_column(tab1, "row_names")
write.csv(tab1, "Output/ST_Combined/DescriptiveTable_anthro.csv")

met_dat %>%
ggplot(aes(x = Country, y = Age)) +
  geom_boxplot() +
  labs(y = "Age (Months)") +
  theme_bw() +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format")
ggsave("Output/ST_Combined/Boxplot_age.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 6, height = 4, units = "in")
met_dat %>%
  mutate(Country2 = ifelse(study == "EE children", "Zambia (SAM)", paste0(Country))) %>%
  pivot_longer(cols = c(WAZ, WLZ, LAZ), names_to = "anthro", values_to = "vals") %>%
  ggplot(aes(x = Country2, y = vals)) +
  geom_boxplot() +
  facet_grid(rows = "anthro", scales = "free") +
  labs(y = "Scores",
       x = "Country") +
  theme_bw() +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format")
ggsave("Output/ST_Combined/WithSAM/Boxplot_anthro2.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 9, height = 8, units = "in")

met_dat %>%
  mutate(Country2 = ifelse(study == "EE children", "Zambia (SAM)", paste0(Country))) %>%
  mutate(Country2 = ifelse(Country2 == "Central African Republic", "CAR", paste0(Country2)))%>%
  ggplot(aes(x = Country2, y = Age, fill = Sex)) +
  geom_boxplot() +
  labs(y = "Age (Months)",
       x = "Country") +
  theme_bw() +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format")
ggsave("Output/ST_Combined/WithSAM/Boxplot_age2.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 6, height = 4, units = "in")

#Relative Abundance
rf_dat %>%
  ps_mutate(Country2 = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
  ps_mutate(Country2 = ifelse(study == "EE children", "Zambia (SAM)", paste0(Country2)))%>%
  comp_barplot(tax_level = "Genus",
               label = "PID",
               n_taxa = 15) +
  #coord_flip() + # horizontal bars are often more readable
  facet_grid(
    cols = vars(Country2),
    scales = "free_x", # this only frees y scale per row in grid faceting
    space = "free_x" # allows bars to be same size by freeing facet heights
  ) +
  theme(axis.text.x = (element_blank()),
        axis.ticks.x = (element_blank()))
ggsave("Output/ST_Combined//WithSAM/RelAbund_gen.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 15, height = 6, units = "in")

########### Stunted Duodenum Analysis ##################
#Zambia vs Bangladesh vs Pakistan vs Central African Republic vs Madagascar


#Import data
load("Data/RData/StuntingCombined.RData")

library(microViz)
library(tidyverse)
library(phyloseq)
rf_dat_st@sam_data[["Country"]] <- relevel(as.factor(rf_dat_st@sam_data[["Country"]]), ref="Zambia")
library(ggpubr)

#Look at clinical features of the samples 
met_dat <- rf_dat_st %>%
  samdat_tbl()
met_dat %>%
  mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
  ggplot(aes(x = Country, y = Age)) +
  geom_boxplot() +
  labs(y = "Age (Months)",
       x = "Country") +
  theme_bw() +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format")
ggsave("Output/ST_Combined/Boxplot_age.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 6, height = 4, units = "in")
met_dat %>%
  mutate(Country = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
  pivot_longer(cols = c(WAZ, WLZ, LAZ), names_to = "anthro", values_to = "vals") %>%
  ggplot(aes(x = Country, y = vals)) +
  geom_boxplot() +
  facet_grid(rows = "anthro", scales = "free") +
  labs(y = "Age (Months)") +
  theme_bw() +
  ggpubr::stat_compare_means(ref.group = "Zambia", label = "p.format")
ggsave("Output/ST_Combined/Boxplot_anthro.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 6, height = 8, units = "in")

#Relative Abundance
rel.plot <- rf_dat_st %>%
  ps_mutate(Country2 = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
  comp_barplot(tax_level = "Genus",
               label = "PID",
               n_taxa = 15) +
  #coord_flip() + # horizontal bars are often more readable
  facet_grid(
    cols = vars(Country2),
    scales = "free_x", # this only frees y scale per row in grid faceting
    space = "free_x" # allows bars to be same size by freeing facet heights
  ) +
  theme(axis.text.x = (element_blank()),
        axis.ticks.x = (element_blank()))
ggsave("Output/ST_Combined/RelAbund_gen_st.png", plot = rel.plot, dpi = 600, device = "png", 
       width = 12, height = 4, units = "in")

#Average
rf_dat_st %>%
  ps_mutate(study = ifelse(Country == "Bangladesh", "Chen", study)) %>%
  #ps_filter(study == "Afribiota") %>%
  phyloseq::merge_samples(group = "Country") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 15) +
  theme(axis.text = element_text(size = 14))
ggsave("Output/ST_Combined/RelAbund_gen_average.png", plot = last_plot(), dpi = 600, device = "png", 
         width = 12, height = 6, units = "in")
  
rf_dat_st %>%
  ps_mutate(Country2 = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
  comp_barplot(tax_level = "Phylum",
               label = "PID") +
  #coord_flip() + # horizontal bars are often more readable
  facet_grid(
    cols = vars(Country2),
    scales = "free_x", # this only frees y scale per row in grid faceting
    space = "free_x" # allows bars to be same size by freeing facet heights
  ) +
  theme(axis.text.x = (element_blank()),
        axis.ticks.x = (element_blank()))
ggsave("Output/ST_Combined/RelAbund_phy_st.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 12, height = 4, units = "in")


#Get mean values of each genera in the plot
#convert to relative abundance & melt to get values
abund_dat <- rf_dat_st%>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_melt()
abund_dat_val <- abund_dat %>%
  filter(OTU == "Streptococcus" | OTU == "Alloprevotella" | OTU == "Gemella" | OTU == "Veillonella" | 
           OTU == "Granulicatella" | OTU == "Fusobacterium" | OTU == "Haemophilus" | OTU == "Rothia") %>%
  group_by(OTU, Country) %>%
  summarise_at(vars(Abundance),
               list(Mean =  mean,
                    SD = sd,
                    Media = median))
my_comparisons <- list( c("Zambia", "Bangladesh"), 
                        c("Zambia", "Central African Republic"), 
                        c("Zambia", "Madagascar"), 
                        c("Zambia", "Pakistan") )

abund_dat %>%
  filter(OTU == "Streptococcus" | OTU == "Alloprevotella" | OTU == "Gemella" | OTU == "Veillonella" | 
           OTU == "Granulicatella" | OTU == "Fusobacterium" | OTU == "Haemophilus" | OTU == "Rothia") %>%
  ggplot(aes(x = Country, y = Abundance, fill = Country)) +
  geom_boxplot() +
  facet_grid(cols = vars(OTU),
             scales = "free_x", # this only frees y scale per row in grid faceting
  ) +
  labs(y = "Relative Abundance") +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # Add pairwise comparisons p-value
  theme_bw()+
  theme(axis.text.x = (element_blank()),
        axis.ticks.x = (element_blank()))
ggsave("Output/ST_Combined/RelAbund_specificgen_st.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 12, height = 6, units = "in")


#Analyze

#Alpha diversity
ad = estimate_richness(rf_dat_st, measures = "Shannon")
ad <- rownames_to_column(ad, "SampleID2")
add <- data.frame(sample_data(rf_dat_st))
add <- rownames_to_column(add, "SampleID2")
add2 = merge(add, ad, by = "SampleID2")
add2 <- add2 %>%
  mutate(Country2 = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))

  
colors = c("Zambia" = "#1b9e77", 
           "Bangladesh" = "#d95f02",
           "CAR" = "#7570b3",
           "Madagascar" = "#e7298a",
           "Pakistan" = "#66a61e")
alpha_plot <- ggplot(add2, aes(Country2, Shannon, fill = Country2)) + 
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

#Shannon 
library(ggpubr)
a <- ggplot(add2,aes(x = LAZ, y = Shannon, colour = Country2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = colors) + labs(color = "Country") +
  stat_cor(method = "pearson", label.x = -6, label.y = 3) +
  theme_classic() +
  facet_grid(cols = vars(Country2), scales = "free_x")

b <- ggplot(add2,aes(x = WAZ, y = Shannon, colour = Country2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = colors) + labs(color = "Country") +
  stat_cor(method = "pearson", label.x = -4, label.y = 3) +
  theme_classic()  +
  facet_grid(cols = vars(Country2), scales = "free_x")

c <- ggplot(add2,aes(x = WLZ, y = Shannon, colour = Country2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = colors) + labs(color = "Country") +
  stat_cor(method = "pearson", label.x = -3, label.y = 3) +
  theme_classic()  +
  facet_grid(cols = vars(Country2), scales = "free_x")

d <- ggplot(add2,aes(x = Age, y = Shannon, colour = Country2)) +
  geom_point() +
  labs(x = "Age (Months)",
       color = "Country") +
  geom_smooth(method = "lm") +
  scale_color_manual(values = colors) + 
  stat_cor(method = "pearson", label.x = 10, label.y = 3) +
  theme_classic()  +
  facet_grid(cols = vars(Country2), scales = "free_x")

e <- ggplot(add2 %>% filter(!is.na(Sex)),
         aes(x = Sex, y = Shannon, colour = Country2)) +
  geom_boxplot() +
  scale_color_manual(values = colors) + labs(color = "Country") +
  theme_classic()  +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Male", "Female"))) +
  facet_grid(cols = vars(Country2), scales = "free_x")

ggarrange(a, b, c, d, e, 
          common.legend = T, ncol = 1, labels = c("A.", "B.", "C.", "D.", "E.")) 
ggsave("Output/ST_Combined/Alpha_anthro.png", plot = last_plot(), dpi = 600, device = "png", width = 10, height = 16, units = "in")


#Beta diversity
beta_plot <- rf_dat_st %>%
  ps_mutate(Country2 = ifelse(Country == "Central African Republic", "CAR", paste0(Country)))%>%
  tax_transform("compositional", rank = "Genus") %>%
  dist_calc("bray") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Country2", plot_taxa = 5, size = 8) +
  scale_color_manual(values = colors) +
  labs(color = "Country") +
  theme_classic()
ggsave("Output/ST_Combined/BetaPlot_st.png", plot = beta_plot, dpi = 600, device = "png", width = 8, height = 6, units = "in")

#paired only
beta_plot2 <- rf_dat_st %>%
  ps_filter(study != "Afribiota") %>%
  tax_transform("compositional", rank = "Genus") %>%
  dist_calc("bray") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Country", plot_taxa = NULL, size = 8) +
  stat_ellipse(aes(linetype = Country, colour = Country)) +
  scale_color_manual(values = colors) +
  #scale_colour_brewer(palette = "Dark2") +
  theme_classic()
ggsave("Output/ST_Combined/BetaPlot_paired.png", plot = beta_plot2, dpi = 600, device = "png", width = 8, height = 6, units = "in")

#paired only with SAM
load("Data/RData/AllMAlStudies.RData")
rf_dat %>%
  ps_filter(study != "Afribiota") %>%
  tax_transform("compositional", rank = "Genus") %>%
  dist_calc("bray") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "study", plot_taxa = NULL, size = 8) +
  stat_ellipse(aes(linetype = study, colour = study)) +
  scale_colour_brewer(palette = "Dark2") +
  theme_classic()
ggsave("Output/ST_Combined/BetaPlot_paired_withSAM.png", plot = last_plot(), dpi = 600, device = "png", width = 8, height = 6, units = "in")

#plot LAZ
rf_dat_st %>%
  ps_filter(study != "Afribiota") %>%
  ps_filter(!is.na(LAZ)) %>%
  tax_transform("compositional", rank = "Genus") %>%
  dist_calc("bray") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "LAZ", shape = "Country", plot_taxa = 1:5, size = 8) +
  #stat_ellipse(aes(linetype = study, colour = study)) +
  scale_color_continuous(type = "viridis") +
  theme_classic()
ggsave("Output/ST_Combined/BetaPlot_paired_laz.png", plot = last_plot(), dpi = 600, device = "png", width = 8, height = 6, units = "in")

#Check PERMANOVA
permanova_func = function(equation,data){
  dist = phyloseq::distance(data, method="bray")
  metadata <- data.frame(sample_data(data))
  require(vegan)
  results <- adonis2(equation, data = metadata, na.action = na.omit)
  coef <- results %>%
    rownames_to_column(var = "Term") %>%
    filter(Term != "Residual") %>%
    filter(Term != "Total")
  print(coef)
}

#regression analysis
#agglomerate at genus level
rf_paired <- rf_dat_st %>%
  ps_filter(study != "Afribiota")
  
rf_dat_st.gen <- tax_glom(rf_paired, taxrank = "Genus", NArm = TRUE)
#Rename taxa with genus names
genus.names <- make.names(tax_table(rf_dat_st.gen)[,'Genus'], unique = TRUE)
taxa_names(rf_dat_st.gen) <- genus.names
rf_dat_st.gen_comp <- rf_dat_st.gen %>%
  tax_transform("compositional", rank = "Genus")
lm_study <- permanova_func(equation = dist ~ study, data = rf_dat_st.gen_comp)
lm_age <- permanova_func(equation = dist ~ Age, data = rf_dat_st.gen_comp)
lm_sex <- permanova_func(equation = dist ~ Sex, data = rf_dat_st.gen_comp)
lm_waz <-permanova_func(equation = dist ~ WAZ, data = rf_dat_st.gen_comp)
lm_whz <- permanova_func(equation = dist ~ WLZ, data = rf_dat_st.gen_comp)
lm_laz <- permanova_func(equation = dist ~ LAZ, data = rf_dat_st.gen_comp)
reg_stats <- Reduce(function(x,y) merge(x, y, all.x = TRUE, all.y = TRUE),
                    list(lm_study, lm_age, lm_sex, lm_waz, lm_whz, lm_laz))
write.csv(reg_stats, "Output/ST_Combined/Beta_reg_st.csv")


#Create combined plot for alpha, beta and relatiev abundance
ggarrange(ggarrange(alpha_plot, beta_plot, nrow = 1, labels = c("A.", "B.")),
          rel.plot, ncol = 1, labels = c("", "C."))
ggsave("Output/ST_Combined/Combined_Fig2.png", plot = last_plot(), dpi = 600, device = "png", width = 12, height = 10, units = "in")



#Genera associations with Anthro
lm_mod_func = function(var,data){
  require(microViz)
  mod_lm <- data %>%
    tax_fix() %>%
    tax_transform("compositional", rank = "Genus") %>%
    #tax_filter(min_prevalence = 0.1, undetected = 0, use_counts = T) %>%
    taxatree_models(type = "lm", 
                    rank = "Genus",
                    trans = "log2", 
                    trans_args = list(zero_replace = "halfmin"),
                    variables = c(var))
  lm_stats <- taxatree_models2stats(mod_lm)
  lm_stats <- taxatree_stats_p_adjust(data = lm_stats, method = "BY", grouping = "rank")
  lm_stats2 <- lm_stats %>% taxatree_stats_get() %>% filter(p.adj.BY.rank < 0.05) %>% filter(term == var)
  print(lm_stats2)
}

#Pakistan
lm_laz <- lm_mod_func(var = "LAZ", data = rf_dat_st %>% ps_filter(Country == "Pakistan"))#nothing
lm_waz <- lm_mod_func(var = "WAZ", data = rf_dat_st %>% ps_filter(Country == "Pakistan"))#nothing
lm_whz <- lm_mod_func(var = "WLZ", data = rf_dat_st %>% ps_filter(Country == "Pakistan"))#nothing
lm_sex <- lm_mod_func(var = "Sex", data = rf_dat_st %>% ps_filter(Country == "Pakistan"))#nothing
lm_age <- lm_mod_func(var = "Age", data = rf_dat_st %>% ps_filter(Country == "Pakistan"))#nothing
#lm_study <- lm_mod_func(var = "study", data = rf_dat_st)#nothing
write.csv(lm_laz, "Output/ST_Combined/lm_clinical_feats_pk.csv")

#Bangladesh
lm_laz <- lm_mod_func(var = "LAZ", data = rf_dat_st %>% ps_filter(Country == "Bangladesh")) #6 taxa
lm_waz <- lm_mod_func(var = "WAZ", data = rf_dat_st %>% ps_filter(Country == "Bangladesh"))#nothing
lm_whz <- lm_mod_func(var = "WLZ", data = rf_dat_st %>% ps_filter(Country == "Bangladesh"))#nothing
lm_sex <- lm_mod_func(var = "Sex", data = rf_dat_st %>% ps_filter(Country == "Bangladesh"))#nothing
lm_age <- lm_mod_func(var = "Age", data = rf_dat_st %>% ps_filter(Country == "Bangladesh"))#nothing
#lm_study <- lm_mod_func(var = "study", data = rf_dat_st)#nothing
write.csv(lm_laz, "Output/ST_Combined/lm_clinical_feats_bg.csv")

#plot significant taxa
plot_data <- rf_paired %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  ps_get() %>%
  ps_otu2samdat(c("Helicobacter", "Ralstonia", "Abiotrophia",
                  "Lachnoanaerobaculum", "Pseudomonas", "Lactobacillus")) %>%
  samdat_tbl()

plot_data %>% pivot_longer(cols = c(Helicobacter, Ralstonia, Abiotrophia, Lachnoanaerobaculum,
                                    Pseudomonas,Lactobacillus), 
                           names_to = "vars", values_to = "val") %>%
  filter(!is.na(val)) %>%
  filter(val != 0) %>%
  ggplot(aes(x = LAZ, y = val, color = Country)) +
  geom_point() +
  scale_y_log10()  +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_color_manual(values = colors) +
  theme_bw() +
  facet_grid(Country~vars, scales = "free")
ggsave("Output/ST_Combined/laz_genera_lm.png", plot = last_plot(), dpi = 600, device = "png", 
       width = 8, height = 4, units = "in")



