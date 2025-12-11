# R Thornley
# 11/09/2025
# Project: P1_COMPADRE_DRAGNET
# get model diagnostics for the GLMMS to show why we picked the Ordbeta over the Hurdle model

# install the mixed up package directly from Github
# remotes::install_github('m-clark/mixedup')

# load packages
library(tidyverse) 
library(glmmTMB)
library(purrr)
library(easystats)
library(mixedup)

################################################################################
# import the data 

# test the different RE outputs as a function of 
data <- read_csv("results/DRAGNet_T0_T1_all.csv")
data <- data %>% mutate(site_block = paste0(site_name, "_", block))
data <- data %>% mutate(site_block_taxon = paste0(site_name, "_", block, "_", New_taxon))

unique(data$trt)

# create a simple df
DIST <- c("Control", "Disturbance")
DIST <- data %>% filter(trt %in% DIST)

NPK <- c("Control", "NPK")
NPK <- data %>% filter(trt %in% NPK)

###############################################################################

# model options 

# model 1 
mod_1 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = DIST)

# model 2 
mod_2 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = beta_family(link = "logit"),
  ziformula = ~ 1 + trt * year_trt,
  data = DIST)

# model 3
mod_3 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = NPK)

# model 4
mod_4 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = beta_family(link = "logit"),
  ziformula = ~ 1 + trt * year_trt,
  data = NPK)

###############################################################################

#function to extract model R2 and variances
create_model_summary <- function(model){
  
  r2_vals <- r2_nakagawa(model)
  dat <- tibble(
    r2_marginal    = unlist(r2_vals[1]),
    r2_conditional = unlist(r2_vals[2]),
    icc_adjusted   = unlist(icc(model)[[1]]),
    icc_conditional = unlist(icc(model)[[2]]),
    icc_unadjusted = unlist(icc(model)[[3]]),
    var.residual   = unlist(get_variance_residual(model)),
    var.fixed      = unlist(get_variance_fixed(model)),
    var.random     = unlist(get_variance_random(model))
  ) %>%
    mutate(
      var_total     = var.residual + var.fixed + var.random,
      prop_fixed    = 100 * var.fixed / var_total,
      prop_random   = 100 * var.random / var_total,
      prop_residual = 100 * var.residual / var_total
    )
  
  dat
}

mod_1 <- create_model_summary(mod_1)
mod_1$model <- "Ordered_beta"
mod_1$treatment <- "Disturbance"
mod_2 <- create_model_summary(mod_2)    
mod_2$model <- "Hurdle"
mod_2$treatment <- "NPK"
mod_3 <- create_model_summary(mod_3)
mod_3$model <- "Ordered_beta"
mod_3$treatment <- "Disturbance"
mod_4 <- create_model_summary(mod_4)
mod_4$model <- "Hurdle"
mod_4$treatment <- "NPK"


results <- rbind(mod_1, mod_2, mod_3, mod_4)
results$random_effect <- "Species"
results <- results[c(1:3, 10:14)]
names(results)
names(results) <- c("Marginal R2", "Conditional R2", "ICC", "Prop. of var. (fixed)", 
                 "Prop. of var. (random)", "Prop. of var. (residual)", "Model", "Treatment")
results <- results %>% dplyr::select(Model, Treatment, everything())

write_csv(results, "figures/SM_S2c_Model_comparison_table_app.csv")

# this is the result in the SM as of 01/12/2025

