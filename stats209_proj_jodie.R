###########################################
# Jodie Bhattacharya
# Stats 209 Project Code
# Created: Nov 2021
###########################################
library(haven)
library(tidyverse)
setwd("C:/Users/jodie/OneDrive/Documents/stats 209/project")

set.seed(1)

df <- read_dta("./AER merged.dta")

frt <- function(data, outcome, treatment, test_stat, nreps=1000, one.sided=TRUE, caption="") {
  simulated_test <- list()
  outcome <- deparse(substitute(outcome))
  
  for(i in 1:nreps){
    
    # Create temporary dataframe to permute so we don't modify the original
    reshuffled <- data 
    
    # Permute the Y column with the 'sample()' function. 
    reshuffled[[outcome]] <- sample(reshuffled[[outcome]], size = nrow(reshuffled), 
                                 replace = FALSE)
    
    # Calculate the T stat
    T_sim <- test_stat(reshuffled, outcome, treatment)
    print(T_sim)
    
    # Append simulated T stat to list
    simulated_test[i] <- T_sim
  }    
  
  simulated_test <- unlist(simulated_test)
  
  # plot histogram of simulated values
  t_obs <- test_stat(data, outcome, treatment)
  plot <- ggplot() +
    ylab("Count") + xlab("Simulated T stat") +
    geom_histogram(aes(x = simulated_test), bins = 20) +
    geom_vline(xintercept = t_obs, size = 1, 
               linetype = "dashed", colour = "black") +
    labs(caption=caption)
  print(plot)
  
  # calculate p-value
  if(one.sided) {
    exceed_count <- length(simulated_test[simulated_test >= 
                                            t_obs])
    p_val <- exceed_count / nreps
    return(mean(simulated_test >= t_obs))
  } else {
    abs_simulated_means <- abs(simulated_test)
    abs_diff_means_obs <- abs(t_obs)
    exceed_count <- length(abs_simulated_means[abs_simulated_means >= 
                                                 abs_diff_means_obs])
    p_val <- exceed_count / nreps
    return(p_val)
  }
  
}

# FRT with diff in means ------------------------------------
# comparing treatment and control with `gave` as outcome ----

# treatment is a binary, 1 or 0
diff_in_means <- function(data, outcome, treatment) {
  outcome <- enquo(outcome)
  treatment <- enquo(treatment)

  # Calculate the means for control and treatment
  mean_1_sim <- mean(data %>% filter(!!(treatment) == 1) 
                         %>% pull(!!outcome))
  mean_0_sim <- mean(data %>% filter(!!(treatment) == 0) 
                           %>% pull(!!outcome))
  
  # Calculate diff in means
  mean_diff_sim <- mean_1_sim - mean_0_sim
  return(mean_diff_sim)
}
diff_in_means(df, gave, treatment)
frt(df, gave, treatment, diff_in_means, nreps=5000)
# 1:1
frt(filter(df, ratio!=2, ratio!=3), gave, treatment, diff_in_means, nreps=5000)
# 2:1
frt(filter(df, ratio!=1, ratio!=3), gave, treatment, diff_in_means, nreps=5000)
# 3:1
frt(filter(df, ratio!=1, ratio!=2), gave, treatment, diff_in_means, nreps=5000)
# $25k
frt(filter(df, size==0 || size==1), gave, treatment, diff_in_means, nreps=5000)
# $50k
frt(filter(df, size==0 || size==2), gave, treatment, diff_in_means, nreps=5000)
# $100k
frt(filter(df, size==0 || size==3), gave, treatment, diff_in_means, nreps=5000)
# low
frt(filter(df, ask==0 || ask==1), gave, treatment, diff_in_means, nreps=5000)
# med
frt(filter(df, ask==0 || ask==2), gave, treatment, diff_in_means, nreps=5000)
# high
frt(filter(df, ask==0 || ask==3), gave, treatment, diff_in_means, nreps=5000)
# blue state
frt(filter(df, blue0==1), gave, treatment, diff_in_means, nreps=5000)
# red state
frt(filter(df, red0==1), gave, treatment, diff_in_means, nreps=5000)
# donated in 2005
frt(filter(df, dormant==0), gave, treatment, diff_in_means, nreps=5000)
# didn't donate in 2005
frt(filter(df, dormant==1), gave, treatment, diff_in_means, nreps=5000)


# comparing treatment and control with `amount` as outcome --------
diff_in_means(df, amount, treatment)
frt(df, amount, treatment, diff_in_means, nreps=5000)

# 1:1
frt(filter(df, ratio!=2, ratio!=3), amount, treatment, diff_in_means, nreps=5000)
# 2:1
frt(filter(df, ratio!=1, ratio!=3), amount, treatment, diff_in_means, nreps=5000)
# 3:1
frt(filter(df, ratio!=1, ratio!=2), amount, treatment, diff_in_means, nreps=5000)
# $25k
frt(filter(df, size==0 || size==1), amount, treatment, diff_in_means, nreps=5000)
# $50k
frt(filter(df, size==0 || size==2), amount, treatment, diff_in_means, nreps=5000)
# $100k
frt(filter(df, size==0 || size==3), amount, treatment, diff_in_means, nreps=5000)
# low
frt(filter(df, ask==0 || ask==1), amount, treatment, diff_in_means, nreps=5000)
# med
frt(filter(df, ask==0 || ask==2), amount, treatment, diff_in_means, nreps=5000)
# high
frt(filter(df, ask==0 || ask==3), amount, treatment, diff_in_means, nreps=5000)
# blue state
frt(filter(df, blue0==1), amount, treatment, diff_in_means, nreps=5000)
# red state
frt(filter(df, red0==1), amount, treatment, diff_in_means, nreps=5000)
# donated in 2005
frt(filter(df, dormant==0), amount, treatment, diff_in_means, nreps=5000)
# didn't donate in 2005
frt(filter(df, dormant==1), amount, treatment, diff_in_means, nreps=5000)

# FRT with diff test stat to compare treatment types ---------------------------
coef_match <- function(data, outcome, treatment) {
  outcome <- enquo(outcome)
  treatment <- enquo(treatment)
  
  df3 <- data %>% select(!!outcome, !!treatment)
  model <- lm(df3)
  return(model$coefficients[2])
}
coef_match(df, gave, treatment)
coef_match(reshuffled, gave, treatment)

# `gave` as outcome
frt(df, gave, treatment, coef_match, nreps=5000)

# `amount` as outcome
frt(df, amount, treatment, coef_match, nreps=5000)

# boostrap to compare subgroups ------------------------------------
library(boot)
boot.fn <- function(data, index, group1, group2, outcome) {
  group1 <- enquo(group1)
  group2 <- enquo(group2)
  outcome <- enquo(outcome)
  
  # Calculate the means for 2 groups
  mean_1_sim <- mean(data[index,] %>% filter(!!(group1) == 1) 
                     %>% pull(!!outcome))
  mean_0_sim <- mean(data[index,] %>% filter(!!(group2) == 1) 
                     %>% pull(!!outcome))
  
  # Calculate diff in means
  mean_diff_sim <- mean_1_sim - mean_0_sim
  
  return(mean_diff_sim)
}
boot(boot.fn, data=df, group1=ratio2, group2=ratio3, outcome=gave, 1000)

# counts
boot.fn2 <- function(data, index, group1, group2, outcome) {
  group1 <- enquo(group1)
  group2 <- enquo(group2)
  outcome <- enquo(outcome)
  
  # Calculate the means for 2 groups
  mean_1_sim <- mean(data[index,] %>% filter(!!(group1) == 1) 
                     %>% pull(!!outcome))
  mean_0_sim <- mean(data[index,] %>% filter(!!(group2) == 1) 
                     %>% pull(!!outcome))
  
  # Calculate diff in means
  mean_diff_sim <- mean_1_sim - mean_0_sim
  
  if(mean_diff_sim > 0)
    return(1)
  else
    return(0)
}
boot.mean = boot(boot.fn2, data=df, group1=ratio3, group2=ratio2, outcome=gave, 1000)
with(boot.mean, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)

boot.count = boot(boot.fn2, data=df, group1=ratio3, group2=ratio2, outcome=amount, 1000)
# p-value
with(boot.count, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)

# FRT to compare treatment types ---------------------------
frt_diff <- function(data, outcome, treatment, nreps=1000, one.sided=TRUE, caption="") {
  simulated_test <- list()
  outcome1 <- deparse(substitute(outcome))
  outcome <- enquo(outcome)
  treatment <- enquo(treatment)
  
  for(i in 1:nreps){
    
    # Create temporary dataframe to permute so we don't modify the original
    reshuffled <- data 
    
    # Permute the Y column with the 'sample()' function. 
    reshuffled[[outcome1]] <- sample(reshuffled[[outcome1]], size = nrow(reshuffled), 
                                    replace = FALSE)
    
    # Calculate the T stat
    # Calculate the means for control and treatment
    mean_1_sim <- mean(reshuffled %>% filter(!!(treatment) == 1) 
                       %>% pull(!!outcome))
    mean_0_sim <- mean(reshuffled %>% filter(!!(treatment) == 0) 
                       %>% pull(!!outcome))
    
    # Calculate diff in means
    T_sim <- mean_1_sim - mean_0_sim
    
    # Append simulated T stat to list
    simulated_test[i] <- T_sim
  }    
  
  simulated_test <- unlist(simulated_test)
  
  # plot histogram of simulated values
  mean_1_sim <- mean(data %>% filter(!!(treatment) == 1) 
                     %>% pull(!!outcome))
  mean_0_sim <- mean(data %>% filter(!!(treatment) == 0) 
                     %>% pull(!!outcome))
  t_obs <- mean_1_sim - mean_0_sim
  plot <- ggplot() +
    ylab("Count") + xlab("Simulated T stat") +
    geom_histogram(aes(x = simulated_test), bins = 20) +
    geom_vline(xintercept = t_obs, size = 1, 
               linetype = "dashed", colour = "black") +
    labs(caption=caption)
  print(plot)
  
  # calculate p-value
  if(one.sided) {
    exceed_count <- length(simulated_test[simulated_test >= 
                                            t_obs])
    p_val <- exceed_count / nreps
    return(mean(simulated_test >= t_obs))
  } else {
    abs_simulated_means <- abs(simulated_test)
    abs_diff_means_obs <- abs(t_obs)
    exceed_count <- length(abs_simulated_means[abs_simulated_means >= 
                                                 abs_diff_means_obs])
    p_val <- exceed_count / nreps
    return(p_val)
  }
  
}
# gave as outcome ----
# 3:1 vs 1:1
frt_diff(df %>% filter(ratio==1 | ratio==3), gave, ratio3, nreps=5000)
# 2:1 vs 1:1
frt_diff(df %>% filter(ratio==1 | ratio==2), gave, ratio2, nreps=5000)
# 3:1 vs 2:1
frt_diff(df %>% filter(ratio==2 | ratio==3), gave, ratio3, nreps=5000)
# $100k vs $25k
frt_diff(df %>% filter(size==1 | size==3), gave, size100, nreps=5000)
# $50k vs $25k
frt_diff(df %>% filter(size==1 | size==2), gave, size50, nreps=5000)
# $100k vs $50k
frt_diff(df %>% filter(size==2 | size==3), gave, size100, nreps=5000)
# High vs. Low
frt_diff(df %>% filter(ask==1 | ask==3), gave, askd3, nreps=5000)
# Med vs. Low
frt_diff(df %>% filter(ask==1 | ask==2), gave, askd2, nreps=5000)
# High vs. Med
frt_diff(df %>% filter(ask==3 | ask==2), gave, askd3, nreps=5000)

# amount as outcome ----
# 3:1 vs 1:1
frt_diff(df %>% filter(ratio==1 | ratio==3), amount, ratio3, nreps=5000)
# 2:1 vs 1:1
frt_diff(df %>% filter(ratio==1 | ratio==2), amount, ratio2, nreps=5000)
# 3:1 vs 2:1
frt_diff(df %>% filter(ratio==2 | ratio==3), amount, ratio3, nreps=5000)
# $100k vs $25k
frt_diff(df %>% filter(size==1 | size==3), amount, size100, nreps=5000)
# $50k vs $25k
frt_diff(df %>% filter(size==1 | size==2), amount, size50, nreps=5000)
# $100k vs $50k
frt_diff(df %>% filter(size==2 | size==3), amount, size100, nreps=5000)
# High vs. Low
frt_diff(df %>% filter(ask==1 | ask==3), amount, askd3, nreps=5000)
# Med vs. Low
frt_diff(df %>% filter(ask==1 | ask==2), amount, askd2, nreps=5000)
# High vs. Med
frt_diff(df %>% filter(ask==3 | ask==2), amount, askd3, nreps=5000)

# process results into latex ----------------
library(readr)
library(kableExtra)
treatment_diff <- read_csv("treatment_diff.csv", 
                           col_types = cols(X7 = col_skip(), X8 = col_skip(), 
                                            X9 = col_skip(), X10 = col_skip(), 
                                            X11 = col_skip()))

kable(treatment_diff, "latex", booktabs = TRUE, digits=3, caption="Comparing treatment types: gave as outcome", label="trt_trt_gave",
      col.names=c("","tau", "SE", "p-value", "$R^2$", "p-value")) %>%
  kable_styling(latex_options = "striped", position="center", full_width="F") %>%
  add_header_above(c(" " = 1, "Lin's estimator" = 4, "FRT" = 1))
# footnote("*** Significant at, or below, 1 percent. ** Significant at, or below, 5 percent. * Significant at, or below, 10 percent.") %>%

treatment_diff2 <- read_csv("treatment_diff2.csv")
kable(treatment_diff, "latex", booktabs = TRUE, digits=3, caption="Comparing treatment types: amount as outcome", label="trt_trt_amt",
      col.names=c("","tau", "SE", "p-value", "$R^2$", "p-value")) %>%
  kable_styling(latex_options = "striped", position="center", full_width="F") %>%
  add_header_above(c(" " = 1, "Lin's estimator" = 4, "FRT" = 1))

results1 <- read_csv("results1.csv", col_types = cols(X9 = col_skip(), 
                                                      +     X10 = col_skip(), X11 = col_skip()))

kable(results1, "latex", booktabs = TRUE, digits=3, caption="Treatment effects: gave as outcome", label="results1",
      col.names=c("","Treatment","1:1","2:1","3:1","$25,000","$50,000","$100,000")) %>%
  kable_styling(latex_options = c("striped", "scale_down"), position="center") %>%
  add_header_above(c(" " = 2, "Match ratio" = 3, "Threshold" = 3)) %>%
  footnote("*** Significant at, or below, 1 percent. ** Significant at, or below, 5 percent. * Significant at, or below, 10 percent.")

results2 <- read_csv("results2.csv")
kable(results2, "latex", booktabs = TRUE, digits=3, caption="Treatment effects: gave as outcome (cont.)", label="results2",
      col.names=c("", "Low","Med","High","Blue State","Red State","Donated in 2005","Didn't donate in 2005")) %>%
  kable_styling(latex_options = c("striped", "scale_down"), position="center") %>%
  add_header_above(c(" " = 1, "Example amount" = 3, "Covariates" = 4))

results3 <- read_csv("results3.csv")
kable(results3, "latex", booktabs = TRUE, digits=3, caption="Treatment effects: amount as outcome", label="results3",
      col.names=c("","Treatment","1:1","2:1","3:1","$25,000","$50,000","$100,000")) %>%
  kable_styling(latex_options = c("striped", "scale_down"), position="center") %>%
  add_header_above(c(" " = 2, "Match ratio" = 3, "Threshold" = 3))

results4 <- read_csv("results4.csv")
kable(results4, "latex", booktabs = TRUE, digits=3, caption="Treatment effects: amount as outcome (cont.)", label="results4",
      col.names=c("", "Low","Med","High","Blue State","Red State","Donated in 2005","Didn't donate in 2005")) %>%
  kable_styling(latex_options = c("striped", "scale_down"), position="center") %>%
  add_header_above(c(" " = 1, "Example amount" = 3, "Covariates" = 4)) %>%
  footnote("*** Significant at, or below, 1 percent. ** Significant at, or below, 5 percent. * Significant at, or below, 10 percent.")

ml1 <- read_csv("ml1.csv")
ml1 <- rbind(ml1, c("x", "x", "x"))
kable(ml1, "latex", booktabs = TRUE, digits=3, caption="Lin's estimator classification: `gave` on `treatment`", label="ml1",
      col.names=c("", "Estimated ATE", "SE")) %>%
  kable_styling(latex_options = c("striped"), position="center")

ml2 <- read_csv("ml2.csv")
kable(ml2, "latex", booktabs = TRUE, digits=3, caption="Lin's estimator regression: `amount` on `treatment`", label="ml2",
      col.names=c("", "Estimated ATE", "SE")) %>%
  kable_styling(latex_options = c("striped"), position="center")
