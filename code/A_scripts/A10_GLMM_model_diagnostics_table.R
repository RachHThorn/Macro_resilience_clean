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

# model 1: 
mod_1 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = DIST)

# model 2: 
mod_2 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = beta_family(link = "logit"),
  ziformula = ~ 1 + trt * year_trt,
  data = DIST)

# model 3
mod_3 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = beta_family(link = "logit"),
  ziformula = ~ 1 + trt * year_trt,
  data = DIST)

# model 4
mod_4 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = beta_family(link = "logit"),
  ziformula = ~ 1 + trt * year_trt,
  data = DIST)

# model 5
mod_5 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = NPK)

# model 6
mod_6 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = NPK)

# model 7
mod_7 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = NPK)


###############################################################################

# now these 2 models won't converge and I don't get why -
# leave for now and figure out later

# model 3:
mod_3 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|taxon_site) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = data)

# model 4: 
mod_4 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|taxon_site) + (1|site_block_taxon),
  family = beta_family(link = "logit"),
  ziformula = ~ 1 + trt * year_trt,
  data = data)


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
mod_2 <- create_model_summary(mod_2)    
mod_3 <- create_model_summary(mod_3)
mod_4 <- create_model_summary(mod_4)

results <- rbind(mod_1, mod_2, mod_3, mod_4)

################################################################################

mod_1$model <- "Ordbeta"
mod_2$model <- "Hurdle"

both <- rbind(mod_1, mod_2)
both <- both[c(1:3, 10:13)]
names(both) <- c("Marginal R2", "Conditional R2", "ICC", "Prop. of var. (fixed)", 
                 "Prop. of var. (random)", "Prop. of var. (residual)", "Model")
both <- both %>% dplyr::select(Model, everything())

write_csv(both, "figures/Model_comparison_table_app.csv")

###############################################################################
# RUN THIS CODE FOR ALL MODELS
###############################################################################

data <- read_csv("results/DRAGNet_T0_T1_all.csv")
data$taxon_site
data$New_taxon
experiments <- c("Control", "Disturbance")

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


get_model_taxon_site <- function(data, experiments) {
  
  data <- data %>% filter(trt %in% experiments)
  data <- data %>% mutate(site_block_taxon = paste0(site_name, "_", block, "_", New_taxon))
  
  hurdle_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|taxon_site) + (1|site_block_taxon),
    family = beta_family(link = "logit"),
    ziformula = ~ 1 + trt * year_trt,
    data = data)

  betareg_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|taxon_site) + (1|site_block_taxon),
    family = ordbeta(link = "logit"),
    data = data)

  results <-  bind_rows(
    create_model_summary(hurdle_mod, component = "cond") %>% mutate(model = "Hurdle_cond"),
    create_model_summary(hurdle_mod, component = "zi")   %>% mutate(model = "Hurdle_zi"),
    create_model_summary(betareg_mod)                    %>% mutate(model = "Ordbeta"))
  
  return(results)
}

get_model_taxon_site(data, experiments)


################################################################################
# Wrapper to process each experiment (treatment)
################################################################################

process_time_period <- function(file, time_label) {
  dat <- read_csv(file)
  
  exp_list <- list(
    DIST  = c("Control", "Disturbance"),
    NPK   = c("Control", "NPK"),
    INTER = c("Control", "NPK+Disturbance")
  )
  
  taxon <- map_dfr(names(exp_list), function(e) {
    get_RE_taxon(dat, exp_list[[e]]) %>%
      mutate(time_period = time_label, experiment = e)
  })
  
  taxon_site <- map_dfr(names(exp_list), function(e) {
    get_RE_taxon_site(dat, exp_list[[e]]) %>%
      mutate(time_period = time_label, experiment = e)
  })
  
  list(taxon = taxon, taxon_site = taxon_site)
}
  

  
see <- get_model_taxon_site(data, experiments)
  
# apply function
