# R Thornley 
# 09/12/2025
# Project: P1_COMPADRE_DRAGNET
# Re-run regressions with a sub-set of species for the appendix to see if the results change

library(tidyverse)
library(purrr)
library(robustbase) # runs robust OLS models
library(mixedup)

################################################################################
# Instructions
################################################################################

# 1) Load and tidy demographic and species cover change (random effects) data
# 2) Export some useful dfs for summary information relating to final modelling data set
# 3) Create and tidy a nested df for modelling
# 4) Apply robust estimators: use robustbase package
# 5) VISUALISE (ROUGH) results of the models / filter for significant terms

################################################################################
# 1) Load and tidy demographic and species cover change (random effects) data
################################################################################

# Load demography variables and find the mean of each Taxon and trait combination
# for modelling purposes
# log transform demographic variables 

# list of the metrics that benefit from being logged
needs_natural_logging <- c("R0", "percapita_repro", "age_repro", "Lmean", "Lmax", 
                           "T_generation", "Reactivity", "FirstStepAtt")


demo <- read_csv("results/all_COMPADRE_metrics.csv") %>% 
  filter(DRAGNet == TRUE) %>% 
  rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var", "Demo_value" = "value") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  dplyr::select(Taxon, Demo_trait, Demo_value) %>%
  group_by(Taxon, Demo_trait) %>%
  mutate(Sample_size = n(),
         s = sd(Demo_value, na.rm = TRUE),
         se = s/sqrt(Sample_size),
         Demo_value_mean = mean(Demo_value, na.rm = TRUE)) %>%
  mutate(Taxon_label = str_replace_all(Taxon, "_", " ")) %>%
  mutate(Taxon_label = paste0(Taxon_label, " n = ", Sample_size)) %>%
  rename(Demo_se = se, Demo_sd = s)

# create a new variable for the logged values and a log flag column
demo <- demo %>%
  mutate(Demo_value = if_else(Demo_trait %in% needs_natural_logging, log(Demo_value), Demo_value)) %>% 
  mutate(Log_flag = if_else(Demo_trait %in% needs_natural_logging, "logged", "not_logged"))

# n shows us the number of matrices we are using for each data point
demo %>% group_by(Taxon, Demo_trait) %>% count()

# save a vector of taxa names to use in the modelling
comp_taxa <- unique(demo$Taxon) # 42 taxa

# read in effect sizes from the whole of DRAGNet 
taxa <- read_csv("results/RE_SE_Taxon_all_DRAGNet.csv")
unique(taxa$time_period)

# select the random effects for the var of interest only (New_taxon)
# create a new ID variable the makes a note if data are in the compadre data set
# change some labels in the data frame to tidy it up
taxa <- 
  taxa %>% 
  filter(str_detect(group_var, "New_taxon")) %>%
  mutate(ID = if_else(group %in% comp_taxa, "1", "0")) %>%
  mutate(time_period = case_when(time_period == 1 ~ "T0-T1",
                                 time_period == 2 ~ "T0-T2",
                                 time_period == 3 ~ "T0-T3")) %>%
  dplyr::select(group, value, se, model, time_period, experiment, ID) %>%
  rename(Taxon = group, RE = value, RE_se = se) %>%
  filter(ID == "1")
names(taxa)
unique(taxa$Taxon) # 41 unique taxa here

# first simplify the demo data to just the mean values
demo <- 
  demo %>% 
  group_by(Taxon, Demo_trait) %>%
  summarise(Demo_value = mean(Demo_value, na.rm = TRUE))

# check there is only one row per species / demo_trait combo 
demo %>% 
  group_by(Taxon, Demo_trait) %>% 
  count() %>%
  filter(n > 1) 

# join both data sets
both <- taxa %>% left_join(demo, by = "Taxon")
names(both)
unique(both$model)
unique(both$time_period)
unique(both$experiment)
unique(both$Taxon)

both <- both %>% filter(model == "Ordbeta")

################################################################################
# 2) Export some useful dfs for summary information relating to final modelling data set
################################################################################

# export that list of species we are using in the final models
taxa_list <- unique(both$Taxon)
length(taxa_list) # 41 taxa here
# saveRDS(taxa_list, "results/List_taxa_OLS_mods.R")

# save the modelling data
# write_csv(both, "results/OLS_modelling_master_data.csv")

################################################################################
# 3) Create and tidy a nested df for modelling with different levels of data filtering
################################################################################

distances <- read_csv("results/all_distances.csv")

med_dist <- distances %>%
  group_by(Taxon) %>%
  summarise(median_distance = median(distance), .groups = "drop") %>%
  arrange(desc(median_distance))   # highest distance first

med_dist
max_drop <- 9  
max_drop <- 9 # or something else, but must be <= nrow(med_dist) - 1

# default uses the Tukey's bisquare (biweight) function
fit_rlm_safe_lmrob <- 
  purrr::possibly(~ lmrob(RE ~ Demo_value, method = "MM",
                          control = lmrob.control(psi = "bisquare"), data = .x),
                  otherwise = NULL)


# build the augment function manually in case of failed models
augment_lmrob <- function(m) {
  if (is.null(m)) return(tibble())
  df <- model.frame(m)
  tibble::as_tibble(df) |>
    dplyr::mutate(
      .fitted    = as.numeric(fitted(m)),
      .resid     = as.numeric(resid(m)),
      .rweight   = if (!is.null(m$rweights)) as.numeric(m$rweights) else NA_real_,
      .std.resid = .resid / sigma(m)   # sigma(m) is the robust scale for lmrob
    )
}


run_robust_models_for_subset <- function(species_subset, both_data) {
  
  both_sub <- both_data %>% 
    filter(Taxon %in% species_subset)
  
  # nest
  nested_sub <- both_sub %>% 
    group_by(model, time_period, experiment, Demo_trait) %>% 
    nest()
  
  # clean each nested df
  nested_clean <- nested_sub %>%
    mutate(data = purrr::map(
      data,
      ~ .x %>%
        dplyr::select(Taxon, RE, Demo_value) %>%
        filter(is.finite(Demo_value)) %>%
        drop_na(Demo_value)
    ))
  
  # run robust models
  robust_models <- nested_clean %>%
    mutate(mod = purrr::map(data, fit_rlm_safe_lmrob)) %>%
    mutate(results = purrr::map(mod, ~ if (is.null(.x)) tibble() else broom::tidy(.x)))
  
  # extract slope for Demo_value only
  robust_results <- robust_models %>%
    dplyr::select(time_period, experiment, model, Demo_trait, results) %>%
    unnest(results) %>%
    filter(term == "Demo_value")
  
  return(robust_results)
}

n_sp <- nrow(med_dist)
max_drop <- min(9, n_sp - 1)  # just to be safe

sens_grid <- tibble(
  n_dropped = 0:max_drop
) %>%
  mutate(
    species_kept = map(
      n_dropped,
      ~ med_dist$Taxon[(.x + 1):n_sp]  # drop top .x high-distance species
    )
  )

sensitivity_results <- sens_grid %>%
  mutate(
    model_output = map(species_kept, ~ run_robust_models_for_subset(.x, both))
  ) %>%
  unnest(model_output)


sensitivity_results %>%
  filter(model == "Ordbeta") %>%
  arrange(Demo_trait, time_period, experiment, n_dropped) %>%
  head()

#################################################################################
# have a quick look at the results
################################################################################

# 1) first chaneg some names and filter the sensitivity results 
# by the hypotheses we want to test
# create the data set for plotting

sensitivity_results <-
  sensitivity_results %>% mutate(
    New_demo_trait = case_when(
      Demo_trait == "T_generation"    ~ "Generation time",
      Demo_trait == "R0"              ~ "Net reproductive rate",
      Demo_trait == "percapita_repro" ~ "Per capita reproduction",
      Demo_trait == "Lmean"           ~ "Mean life expectancy",
      Demo_trait == "age_repro"       ~ "Age at first reproduction",
      Demo_trait == "Lambda"          ~ "Population growth rate",
      Demo_trait == "Reactivity"      ~ "Reactivity",
      Demo_trait == "FirstStepAtt"    ~ "First step attenuation",
      Demo_trait == "Lmax"            ~ "Maximum longevity", # changed this name
      TRUE                            ~ Demo_trait
    ),
    hypoth = case_when(
      New_demo_trait == "Net reproductive rate"     & experiment == "DIST" ~ "H1",
      New_demo_trait == "Per capita reproduction"   & experiment == "DIST" ~ "H1",
      New_demo_trait == "Generation time"           & experiment == "DIST" ~ "H2",
      New_demo_trait == "Age at first reproduction" & experiment == "DIST" ~ "H2",
      New_demo_trait == "Reactivity"                & experiment == "DIST" ~ "H3",
      New_demo_trait == "First step attenuation"    & experiment == "DIST" ~ "H3",
      New_demo_trait == "Population growth rate"    & experiment == "DIST" ~ "H4",
      New_demo_trait == "Maximum longevity"    & experiment == "NPK"  ~ "H5", # changed this name
      New_demo_trait == "Mean life expectancy"      & experiment == "NPK"  ~ "H5",
      TRUE ~ "none"
    )) %>%
  filter(hypoth != "none") %>%
  mutate(
    sign = if_else(estimate >= 0, "pos", "neg"),
    overlaps0 = (estimate - std.error) <= 0 & (estimate + std.error) >= 0,
    sig = if_else(overlaps0, "ns", "sig"),
    col_cat = factor(paste(sign, sig, sep = "_"),
                     levels = c("neg_ns","neg_sig","pos_ns","pos_sig"))) %>%
  mutate(
    
  )

# 2) The plot out the change in the regression co-efficient as a response
# to dropping species with high distances from the model

plot <- sensitivity_results %>%
  filter(model == "Ordbeta") %>%
  filter(experiment != "INTER") %>%
  left_join(sens_grid, by = "n_dropped") %>%
  ggplot(aes(
    x     = n_dropped,
    y     = estimate,
    group = Demo_trait,
    color = hypoth,
    alpha = sig
  )) +
  geom_line() +
  geom_point() +
  geom_text(
    data        = label_df,
    aes(label = hypoth),
    hjust       = -0.1,
    show.legend = FALSE
    # if label_df doesn't have the same x/y columns you can add:
    # , inherit.aes = FALSE
  ) +
  facet_grid(
    experiment ~ time_period,
    scales   = "free",
    labeller = labeller(
      experiment = c(
        NPK  = "Nutrient Addition",
        DIST = "Disturbance"
      )
    )
  ) +
  theme_bw() +
  labs(
    color = "Hypotheses",
    x     = "Number of high-climatic-distance species dropped",
    y     = "Robust regression coefficient"
  ) +
  scale_alpha_manual(
    name   = "Significance",
    values = c("sig" = 1, "ns" = 0.3),
    labels = c(
      "ns"  = "non-significant at p < 0.05",
      "sig" = "significant at p < 0.05"
    )
  ) +
  scale_x_continuous(breaks = 0:9) +
  coord_cartesian(xlim = c(0, max(sensitivity_results$n_dropped) + 1)) +
  theme(
    strip.background = element_blank(),
    panel.grid       = element_blank(),
    strip.text.x     = element_text(size = 14),
    strip.text.y     = element_text(size = 14)
  )

plot

ggsave("figures/climate_senstivity_coeffcient_plot.pdf", plot, height = 5, width = 8)

