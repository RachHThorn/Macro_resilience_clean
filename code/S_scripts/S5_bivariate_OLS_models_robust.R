# R Thornley
# 08/09/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S5_bivariate_OLS_models_robust

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
comp_taxa <- unique(demo$Taxon) # 41 taxa

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
unique(taxa$Taxon) # 40 unique taxa here

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

################################################################################
# 2) Export some useful dfs for summary information relating to final modelling data set
################################################################################

# export that list of species we are using in the final models
taxa_list <- unique(both$Taxon)
length(taxa_list) # 40 taxa here
saveRDS(taxa_list, "results/List_taxa_OLS_mods.R")

# save the modelling data
write_csv(both, "results/OLS_modelling_master_data.csv")

################################################################################
# 3) Create and tidy a nested df for modelling
################################################################################

nested_both <- both %>% group_by(model, time_period, experiment, Demo_trait) %>% 
  nest()

# clean up data frame by removing rows with NA/ INF etc.
nested_clean <- 
  nested_both %>%
  mutate(data = map(data, ~ .x %>% dplyr::select(Taxon, RE, Demo_value))) %>%
  mutate(data = map(data, ~ .x %>% 
                      filter(if_all(everything(), ~ is.finite(Demo_value))) %>%  # keep only finite values
                      drop_na(Demo_value)))

################################################################################
# 4) APPLY Robust estimators: use robustbase package
################################################################################
# this function / package offers the option of p-values 
# it has also been shown to be more robust on small data sets and where there are outliers in x and y

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

# default uses the Tukey's bisquare (biweight) function
# build a safe version of the lmrob function
fit_rlm_safe_lmrob <- 
  purrr::possibly(~ lmrob(RE ~ Demo_value, method = "MM",
                          control = lmrob.control(psi = "bisquare"), data = .x),
                  otherwise = NULL)

# build the MASS augment function manually
# so it returns something on the lmrob model output
# the default MASS function does not return anything
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

# Now run the models / get results / 
robust_models_3 <- 
  nested_clean %>%
  mutate(mod = map(data, fit_rlm_safe_lmrob)) %>%
  mutate(results = map(mod, ~ if (is.null(.x)) tibble() else broom::tidy(.x))) %>%
  mutate(aug = map(mod, ~ if (is.null(.x)) tibble() else augment_lmrob(.x)))

# save to disk for model checking and plotting
saveRDS(robust_models_3, file = "results/robust_model_results.rds")

################################################################################
# 5) VISUALISE (ROUGH) results of the models / filter for significant terms
################################################################################

# extract the results of these models 
# filter out the estimate and p-value for the New_demo_value term only
robust_results_3 <- robust_models_3 %>%
  dplyr::select(time_period, experiment, model, Demo_trait, results) %>%
  unnest(results) %>%
  filter(term == "Demo_value")

# visualise the effect sizes and their std. error for the Ordbeta model
robust_results_3 %>%
  filter(model == "Ordbeta") %>%
  ggplot(aes(estimate, Demo_trait)) + 
  theme_bw()+
  geom_point() +
  geom_errorbar(aes(xmin = estimate - std.error, xmax = estimate + std.error), size = 0.3) +
  facet_grid(time_period ~ experiment)+
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed")+
  ggtitle("Robust_OLS_effects")
# save if required
# ggsave("figures/Robust_3_effects_logged_vars.jpeg", height = 5, width = 12)

# filter out any significant terms for the ordered beta model
robust_results_3 %>% filter(model == "Ordbeta") %>% filter(p.value < 0.05)
# the p-values are very small with this estimator - is this an issue????

