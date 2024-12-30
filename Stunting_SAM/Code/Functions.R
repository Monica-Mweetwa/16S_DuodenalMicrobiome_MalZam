#title: "16S_DuodenalMicrobiome_MalZam - Functions"
#author: Monica N Mweetwa
#date: 30/12/2024

#regression function for analysis of individual malnutrition classes
lin_uni_func = function(equation,data,class2){
  require(lm.beta)
  data2 <- data %>% filter(mal.type == class2)
  # Fit the linear model and calculate standardized beta coefficients
  lm_mod = lm(formula = equation, data = data2)
  lm_mod_st <- lm.beta(lm_mod)
  # Conduct Shapiro-Wilk test for normality of residuals
  shapiro_test <- shapiro.test(residuals(lm_mod))
  # Extract coefficients and confidence intervals
  mod <- summary(lm_mod_st)
  coef = as.data.frame(mod[["coefficients"]])
  p_sharp = shapiro_test$p.value
  coef <- coef %>%
    mutate(lwr=Estimate-1.96*`Std. Error`,
           upr=Estimate+1.96*`Std. Error`,
           `95% CI` = paste0("(",format(lwr,scientific = T),", ",format(upr,scientific = T),")"),
           Number_of_observations = paste0(nobs(lm_mod)),
           shapiro_test = shapiro_test[["statistic"]][["W"]],
           shapiro_test_p = p_sharp) %>%
    rownames_to_column(var = "Term") %>%
    select(-lwr, -upr)%>%
    filter(Term != "(Intercept)")
  print(coef)
}

#regression function for analysis of comdbined data
lin_uni_func2 = function(equation,data){
  require(lm.beta)
  require(broom)
  # Fit the linear model and calculate standardized beta coefficients
  lm_mod = lm(formula = equation, data = data)
  lm_mod_st <- lm.beta(lm_mod)
  # Conduct Shapiro-Wilk test for normality of residuals
  shapiro_test <- shapiro.test(residuals(lm_mod))
  # Extract coefficients and confidence intervals
  mod <- summary(lm_mod_st)
  coef = as.data.frame(mod[["coefficients"]])
  p_sharp = shapiro_test$p.value
  coef <- coef %>%
    mutate(lwr=Estimate-1.96*`Std. Error`,
           upr=Estimate+1.96*`Std. Error`,
           `95% CI` = paste0("(",format(lwr,scientific = T),", ",format(upr,scientific = T),")"),
           Number_of_observations = paste0(nobs(lm_mod)),
           shapiro_test = shapiro_test[["statistic"]][["W"]],
           shapiro_test_p = p_sharp) %>%
    rownames_to_column(var = "Term") %>%
    select(-lwr, -upr)%>%
    filter(Term != "(Intercept)")
  print(coef)
}

#function to get output of permanova analysis
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

#QQ-plots from linear regression residuals
qq_plot_func <- function (equation,title, data = data) {
  lm_mod = lm(formula = equation, data = data)
  qq_plot <- ggplot(data, aes(sample = residuals(lm_mod))) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(
    title = paste0(title),
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal()
  return(qq_plot)
}

# Create residuals vs. fitted values plot
fit_plot_func <- function (equation, title, data = data) {
  lm_mod = lm(formula = equation, data = data)
  residuals_plot <- ggplot(data, aes(x = fitted(lm_mod), y = residuals(lm_mod))) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Fitted Values",
    y = "Residuals",
    title = paste0(title)
  ) +
  theme_minimal()
}

ggplotRegression <- function (equation, y_label, x_label, data = data) {
  fit = lm(equation, data = data)
  Pval = as.character(ifelse(signif(summary(fit)$coef[2,4], 5)<0.001, "<0.001", round(signif(summary(fit)$coef[2,4], 5), 3)))
  
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    ylab(y_label) +xlab(x_label)+
    labs(title = paste("Adj R2 = ",round(signif(summary(fit)$adj.r.squared, 5), 3),
                       "Intercept =",round(signif(fit$coef[[1]],5 ),3),
                       " Slope =",Pval))+
    theme(plot.title = element_text(size = 6))
  #" P =",round(signif(summary(fit)$coef[2,4], 5),3)))
}



