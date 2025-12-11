# R Thornley
# 29/08/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S3_GLMM_extract_random_effects_all_DRAG_data

rm(list = ls())

################################################################################
# Instructions
################################################################################

# 1) Create functions that run the GLMMS and extract the random effects at species 
# and species in site level
# 2) Wrapper function to process the data for each experiment (treatment + control)
# 3) Run functions for all time periods and save to file

################################################################################
# 1) Create functions that run the GLMMS and extract the random effects at species 
# and species in site (population) level
################################################################################

get_RE_taxon <- function(data, experiments) {
  
  data <- data %>% filter(trt %in% experiments)
  data <- data %>% mutate(site_block_taxon = paste0(site_name, "_", block, "_", New_taxon))
  
  hurdle_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
    family = beta_family(link = "logit"),
    ziformula = ~ 1 + trt * year_trt,
    data = data
  )
  
  betareg_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
    family = ordbeta(link = "logit"),
    data = data
  )
  
  bind_rows(
    extract_random_coefs(hurdle_mod, component = "cond") %>% mutate(model = "Hurdle_cond"),
    extract_random_coefs(hurdle_mod, component = "zi")   %>% mutate(model = "Hurdle_zi"),
    extract_random_coefs(betareg_mod)                    %>% mutate(model = "Ordbeta")
  )
}

get_RE_taxon_site <- function(data, experiments) {
  
  data <- data %>% filter(trt %in% experiments)
  data <- data %>% mutate(site_block_taxon = paste0(site_name, "_", block, "_", New_taxon))
  
  hurdle_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|taxon_site) + (1|site_block_taxon),
    family = beta_family(link = "logit"),
    ziformula = ~ 1 + trt * year_trt,
    data = data
  )
  
  betareg_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|taxon_site) + (1|site_block_taxon),
    family = ordbeta(link = "logit"),
    data = data
  )
  
  bind_rows(
    extract_random_coefs(hurdle_mod, component = "cond") %>% mutate(model = "Hurdle_cond"),
    extract_random_coefs(hurdle_mod, component = "zi")   %>% mutate(model = "Hurdle_zi"),
    extract_random_coefs(betareg_mod)                    %>% mutate(model = "Ordbeta")
  )
}

################################################################################
# 3) Wrapper function to process the data for each experiment (treatment + control)
################################################################################

process_time_period <- function(file, time_label) {
  dat <- read_csv(file)
  
  exp_list <- list(
    DIST  = c("Control", "Disturbance"),
    NPK   = c("Control", "NPK"),
    INTER = c("Control", "NPK+Disturbance")
  )
  
  taxon <- purrr::map_dfr(names(exp_list), function(e) {
    get_RE_taxon(dat, exp_list[[e]]) %>%
      mutate(time_period = time_label, experiment = e)
  })
  
  taxon_site <- purrr::map_dfr(names(exp_list), function(e) {
    get_RE_taxon_site(dat, exp_list[[e]]) %>%
      mutate(time_period = time_label, experiment = e)
  })
  
  list(taxon = taxon, taxon_site = taxon_site)
}

################################################################################
# 4) Run functions for all time periods and save to file
################################################################################

# use the file paths from 
files <- c(
  "results/DRAGNet_T0_T1_all.csv",
  "results/DRAGNet_T0_T2_all.csv",
  "results/DRAGNet_T0_T3_all.csv")

results <- purrr::imap(files, ~ process_time_period(.x, .y))

all_taxon <- purrr::map(results, "taxon") %>% bind_rows()
all_taxon_site <- purrr::map(results, "taxon_site") %>% bind_rows()

write_csv(all_taxon,      "results/RE_SE_Taxon_all_DRAGNet.csv")
write_csv(all_taxon_site, "results/RE_SE_Taxon_site_all_DRAGNet.csv")
