########### All Malnourished Duodenum Analysis ##################
#Zambia vs Bangladesh vs Pakistan vs Central African Republic vs Madagascar

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
load("Data/RData/BEED_SEEM_BEECH_AbsASV.RData") #ps_SB_norm.ra.scaled -- BEED, BEECH, ME, SEEM absolute counts

comparisons1 <- list(c("BEECH","BEED"), c("BEECH","SEEM"), c("BEED","SEEM"))
colors_study = c("BEECH" = "#1b9e77", 
                 "BEED" = "#d95f02",
                 "AFRIBIOTA" = "#ba4aa4",
                 "SEEM" = "#66a61e")

#Merge duplicate PIDs
ps_norm_merged <- phyloseq::merge_samples(ps_SB_norm.ra.scaled, "PID")
# replace metadata as NAs get introduced during merging
meta <- sample_data(ps_SB_norm.ra.scaled)
merged_meta <- meta[!duplicated(meta$PID), ]
merged_meta$sampleID2 <- rownames(merged_meta)
rownames(merged_meta) <- NULL
merged_meta2 <- merged_meta %>% column_to_rownames(var = "PID")
sample_data(ps_norm_merged) <- merged_meta2

#########################################################################
############# Alpha diversity - BEED, BEECH and SEEM only ###############
#########################################################################
ad = estimate_richness(ps_norm_merged, measures = "Shannon")
ad <- rownames_to_column(ad, "PID")
add <- data.frame(sample_data(ps_norm_merged))
add <- rownames_to_column(add, "PID")
add2 = merge(add, ad, by = "PID")
add3 <- add2 %>% 
  filter(study != "ME") %>% #remove samples from children with SAM
  filter(WLZ_cat != "Severely wasted") %>% #remove samples from children with SAM
  filter(LAZ_cat != "Not stunted") #remove children that are not stunted

##linear regression analysis - All##
lm_func = function(equation,data){
  # Fit the linear model and calculate standardized beta coefficients
  lm_mod = glm(formula = equation, data = data)
  # Extract coefficients and confidence intervals
  mod <- summary(lm_mod)
  coef = as.data.frame(mod[["coefficients"]])
  coef <- coef %>%
    mutate(lwr=Estimate-1.96*`Std. Error`,
           upr=Estimate+1.96*`Std. Error`,
           `95% CI` = paste0("(",format(lwr,scientific = T),", ",format(upr,scientific = T),")"),
           Number_of_observations = paste0(nobs(lm_mod))) %>%
    rownames_to_column(var = "Term") %>%
    filter(Term != "(Intercept)")
  print(coef)
}
library(car) ##for VIF function
library(corrplot) ##for plotting the correlation

#All samples together
add4 <- add3 %>%
  mutate(Sex_num = ifelse(Sex == "Male", 1, 2),
         study_num = ifelse(study == "BEECH", 1, 2),
         study_num = ifelse(study == "BEED", 3, study_num)) %>%
  select(Age, Sex_num, study_num, LAZ, WAZ, WLZ)
var_add4 <- cor(add4, use = "na.or.complete") # independent variables correlation matrix 

corrplot(var_add4 ,method='number', is.corr = T,type="upper", 
         title='Meta-Analysis Predictors Correlation Matrix',
         mar=c(0,0,1,0)) # inspect for collinearity
export::graph2ppt(width = 10, height = 6.9, aspectr=16/9 ,center=FALSE,
          offy=(10 *9/16)- 4.9 ,offx=0,vector.graphic = FALSE,
          file="Output/AlphaDiv/RegressionModel_Diagnostics_Met.pptx",
          margins = c(top = 0, right = 0, bottom = 0, left= 0),
          upscale=TRUE, append = TRUE)

#Weak correlations between variables: include all in analysis

#univar
mod1_laz <- lm_func(Shannon~LAZ, data = add3)
write_csv(mod1_laz, "Output/AlphaDiv/Shannon_laz_univar.csv")
mod1_waz <- lm_func(Shannon~WAZ, data = add3)
write_csv(mod1_waz, "Output/AlphaDiv/Shannon_waz_univar.csv")
mod1_wlz <- lm_func(Shannon~WLZ, data = add3)
write_csv(mod1_wlz, "Output/AlphaDiv/Shannon_wlz_univar.csv")
mod1_age <- lm_func(Shannon~Age, data = add3)
write_csv(mod1_age, "Output/AlphaDiv/Shannon_age_univar.csv")
#multivar
mod1_laz <- lm_func(Shannon~LAZ + Age + study + Sex, data = add3)
write_csv(mod1_laz, "Output/AlphaDiv/Shannon_laz.csv")
mod1_waz <- lm_func(Shannon~WAZ + Age + study + Sex, data = add3)
write_csv(mod1_waz, "Output/AlphaDiv/Shannon_waz.csv")
mod1_wlz <- lm_func(Shannon~WLZ + Age + study + Sex, data = add3)
write_csv(mod1_wlz, "Output/AlphaDiv/Shannon_wlz.csv")
mod1_age <- lm_func(Shannon~Age + study + Sex, data = add3)
write_csv(mod1_age, "Output/AlphaDiv/Shannon_age.csv")

#model diagnostics
library(export)
library(DHARMa)
model_dignos_func <- function(equation, label){
  simulationOutput  <- DHARMa::simulateResiduals(glm(equation, data = add3))
  export::graph2ppt(plot(simulationOutput, title = label),
            width = 10,height = 4.9,center=FALSE,aspectr=16/9,
            offy=(10 *9/16)- 4.9,offx=0,vector.graphic = F,
            file="Output/AlphaDiv/RegressionModel_Diagnostics_Met.pptx",   margins=0, upscale=TRUE,append = TRUE)
  }
model_dignos_func(equation = "Shannon~WAZ", label = "glm(Shannon~WAZ)")
model_dignos_func(equation = "Shannon~LAZ", label = "glm(Shannon~LAZ)")
model_dignos_func(equation = "Shannon~WLZ", label = "glm(Shannon~WLZ)")
model_dignos_func(equation = "Shannon~Age", label = "glm(Shannon~Age)")
model_dignos_func <- function(equation, label){
  simulationOutput  <- DHARMa::simulateResiduals(glm(equation, data = add3))
  export::graph2ppt(plot(simulationOutput, title = label),
                    width = 10,height = 4.9,center=FALSE,aspectr=16/9,
                    offy=(10 *9/16)- 4.9,offx=0,vector.graphic = F,
                    file="Output/AlphaDiv/RegressionModel_Diagnostics_Met.pptx",   margins=0, upscale=TRUE,append = TRUE)
  a <- as.data.frame(car::vif(lm(equation, dat = add3)))
  barplot(a$GVIF, main = paste0(label," - VIF Values"),horiz = FALSE,  col = "steelblue", ylim = c(0, 2)) 
  export::graph2ppt(width = 10,height = 4.9,center=FALSE,aspectr=16/9,
                    offy=(10 *9/16)- 4.9,offx=0,vector.graphic = F,
                    file="Output/AlphaDiv/RegressionModel_Diagnostics_Met.pptx",   margins=0, upscale=TRUE,append = TRUE)
}
model_dignos_func(equation = "Shannon~WAZ + Age + study + Sex", label = "glm(Shannon~WAZ + Age + study + Sex)")
model_dignos_func(equation = "Shannon~LAZ + Age + study + Sex", label = "glm(Shannon~LAZ + Age + study + Sex)")
model_dignos_func(equation = "Shannon~WLZ + Age + study + Sex", label = "glm(Shannon~WLZ + Age + study + Sex)")
model_dignos_func(equation = "Shannon~Age + study + Sex", label = "glm(Shannon~Age + study + Sex)")

#Shannon scatterplot with anthropometry
library(ggpubr)
a <- add3 %>% 
  ggplot(aes(x = Age, y = Shannon)) +
  geom_point(aes(colour = study)) + theme_bw() +
  scale_color_manual(values = colors_study) + labs(color = "STUDY") +
  labs(x = "Age (months)",
       y = "Shannon's index",
       subtitle = "β = 0.03 , p = 0.034") +
  geom_smooth(method = "lm") 
b <- add3 %>% 
  ggplot(aes(x = LAZ, y = Shannon)) +
  geom_point(aes(colour = study)) + theme_bw() +
  scale_color_manual(values = colors_study) + labs(color = "STUDY") +
  labs(x = "LAZ",
       y = "",
       subtitle = "β = 0.167 , p = 0.012") +
  geom_smooth(method = "lm") 
c <- add3 %>% 
  ggplot(aes(x = WAZ, y = Shannon)) +
  geom_point(aes(colour = study)) + theme_bw() +
  scale_color_manual(values = colors_study) + labs(color = "STUDY") +
  labs(x = "WAZ",
       y = "",
       subtitle = "β = 0.199 , p = 0.007") +
  geom_smooth(method = "lm") 
ggarrange(a, b, c, nrow = 1, common.legend = T, legend = "right")
ggsave("Output/AlphaDiv/Alpha_anthro_meta.png", plot = last_plot(), dpi = 600, device = "png", width = 8, height = 4, units = "in")


