# R Thornley
# 07/10/2025
# Project: P1_COMPADRE_DRAGNET
# S8_comm_OLS_models_robust_scaled
# Fit the robust interaction models scaled with demo vars not logged

rm(list= ls())

################################################################################
# Instructions
################################################################################

# 1) Load and tidy demographic and species cover change (random effects) data
# 2) Export some useful dfs for summary information relating to final modelling data set
# 3) Load community variables; create and tidy a nested df for modelling
# 4) JOIN COMM DATA WITH DEMO/SPECIES COVER AND FILTER FOR HYPOTHESES
# 5) RUN ROBUST MODELS FOR EACH HYPOTHESIS INCLUDING INTERACTION
# 6) TEST THE SPECIES PROVENANCE DATA INTERACTION

################################################################################
# 1) Load and tidy demographic and species cover change (random effects) data
################################################################################

# Load demography variables and find the mean of each Taxon and trait combination
# for modelling purposes
# log transform demographic variables 

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

# n shows us the number of matrices we are using for each data point
demo %>% group_by(Taxon, Demo_trait) %>% count()

# save a vector of taxa names to use in the next stage of modelling
comp_taxa <- unique(demo$Taxon) # 42 taxa
comp_taxa

# read in effect sizes from the whole of DRAGNet 
taxa <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv")
unique(taxa$time_period)
names(taxa)
unique(taxa$group_var)

# select the random effects for the var of interest only (taxon_site)
taxa <- taxa %>% filter(str_detect(group_var, "taxon_site"))
# tidy up the names of the 
taxa <- taxa %>%
  separate(group, into = c("genus", "species", "site_name"),
           sep = "_", extra = "merge",  # merge all remaining pieces into 'site_name'
           fill  = "right", remove = FALSE) %>%
  mutate(Taxon  = str_c(genus, species, sep = "_"))
taxa$Taxon

# create a new ID variable the makes a note if data are in the compadre data set
# change some labels in the data frame to tidy it up
taxa <- taxa %>% 
  mutate(ID = if_else(Taxon %in% comp_taxa, "1", "0")) %>%
  mutate(time_period = case_when(time_period == 1 ~ "T0-T1",
                                 time_period == 2 ~ "T0-T2",
                                 time_period == 3 ~ "T0-T3")) %>%
  dplyr::select(Taxon, group, value, se, model, time_period, experiment, ID, site_name) %>%
  rename(Taxon_site = group, RE = value, RE_se = se) %>%
  filter(ID == "1")
names(taxa)
unique(taxa$Taxon) # 41 unique taxa here

# now join the data
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
# create a site variable as this is now missing

both %>% group_by(site_name, time_period, experiment, Taxon, Demo_trait, model) %>% tally()

# filter out the ordered beta model results only
both <- both %>% filter(model == "Ordbeta") %>% filter(experiment %in% c("DIST","NPK"))
unique(both$experiment)

# Names tidy and add hypotheses numbers
both <- both %>% 
  mutate(New_demo_trait = case_when(
    Demo_trait == "T_generation"    ~ "Generation time",
    Demo_trait == "R0"              ~ "Net reproductive rate",
    Demo_trait == "percapita_repro" ~ "Per capita reproduction",
    Demo_trait == "Lmean"           ~ "Mean life expectancy",
    Demo_trait == "age_repro"       ~ "Age at first reproduction",
    Demo_trait == "Lambda"          ~ "Population growth rate",
    Demo_trait == "Reactivity"      ~ "Reactivity",
    Demo_trait == "FirstStepAtt"    ~ "First step attenuation",
    Demo_trait == "Lmax"            ~ "Maximum length of life",
    TRUE                            ~ Demo_trait))
both <- both %>% mutate(hypoth = case_when(
  New_demo_trait == "Net reproductive rate"     & experiment == "DIST" ~ "H1",
  New_demo_trait == "Per capita reproduction"   & experiment == "DIST" ~ "H1",
  New_demo_trait == "Generation time"           & experiment == "DIST" ~ "H2",
  New_demo_trait == "Age at first reproduction" & experiment == "DIST" ~ "H2",
  New_demo_trait == "Reactivity"                & experiment == "DIST" ~ "H3",
  New_demo_trait == "First step attenuation"    & experiment == "DIST" ~ "H3",
  New_demo_trait == "Population growth rate"    & experiment == "DIST" ~ "H4",
  New_demo_trait == "Maximum length of life"    & experiment == "NPK"  ~ "H5",
  New_demo_trait == "Mean life expectancy"      & experiment == "NPK"  ~ "H5",
  TRUE ~ "none")) %>%
  filter(hypoth != "none")

unique(both$site_name) # 40 sites
unique(both$Taxon_site) # 139 unique site taxa combinations
unique(both$Taxon) # 39 taxa

################################################################################
# 2) Export some useful dfs for summary information relating to final modelling data set
################################################################################

# export that list of species we are using in the final models
taxa_list <- unique(both$Taxon)
length(taxa_list) # 39 taxa here
saveRDS(taxa_list, "results/List_taxa_OLS_mods.R")

# save the modelling data
write_csv(both, "results/OLS_comm_var_modelling_master_data_T0.csv")

################################################################################
# 3) Load community variables; create and tidy a nested df for modelling
################################################################################

# we hypothesized that:
# 1) whether species were in their native range would impact the ability of the demo vars to predict
# 2) the site level complexity of the community: eg species diversity
# 3) The total cover of the plots - a proxy for competition
# 4) The variability of the CSR strategy in the plots

# we need all these variables to be at the site / time / treatment, level of detail
# there should be 1072 rows for each of these variables
# Load and tidy community variables

cover <- read_csv("results/total_cover_T0.csv") %>% 
  dplyr::select(!experiment) 
cover <- cover %>% group_by(site_name, year_trt, trt) %>% mutate(mean_cover = mean(mean_cover))
unique(cover$site_name) # cover values for 54 sites

# read in the diversity data
div <- read_csv("results/diversity_metrics_dragnet_T0.csv") # plot / site
#read in the strategy data
strategy <- read_csv("results/community_specialisation_results_T0.csv") # plot / site
# join the data
all_comm <- div %>% full_join(strategy, join_by("site_name", "year_trt", "trt"))
all_comm <- all_comm %>% left_join(cover, join_by("site_name", "year_trt", "trt"))

# tidy some names
all_comm <- 
  all_comm %>% 
  mutate(trt = case_when(trt == "Disturbance" ~ "DIST",
                         trt == "NPK+Disturbance"  ~ "INTER", 
                         trt == "NPK"  ~ "NPK", 
                         trt == "Control" ~ "CONTROL")) %>%
  rename("experiment" = "trt") %>%
  dplyr::select(-time_period)

# read in in the species provenance data
native <- read_csv("data/DRAGNet_native_status.csv") %>% 
  dplyr::select(site_name, New_taxon, new_provenance) %>%
  rename(Taxon = New_taxon)

################################################################################
# 4) JOIN COMM DATA WITH DEMO/SPECIES COVER AND FILTER FOR HYPOTHESES
################################################################################

# join the community /demo / species change data
all <- both %>% left_join(all_comm, by = c("site_name", "experiment"))
all <- all %>% left_join(native)

# make it clear that the 'time_period' variable relates to the cover change data not to the community metrics
all <- all %>% rename("time_period_cover_change" = "time_period")
names(all)

# safe_scale function to scale the demo variables without returning NAs
safe_scale <- function(x, tol = 1e-12) {
  x_num <- as.numeric(x)
  mu <- mean(x_num, na.rm = TRUE)
  s  <- sd(x_num,   na.rm = TRUE)
  if (is.na(s) || s < tol) {
    # constant or too few observations → return 0 where not NA
    return(ifelse(is.na(x_num), NA_real_, 0))
  }
  (x_num - mu) / s
}

# scale the data before modelling
all_scaled <- all %>% 
  mutate(across(site_time_rich:mean_cover, ~ as.numeric(safe_scale(.x))))

all_scaled <- all_scaled %>%
  group_by(Demo_trait, experiment, time_period_cover_change) %>%
  mutate(Demo_value = safe_scale(Demo_value)) %>%
  ungroup()

# check scaling look ok
range(all_scaled$variability, na.rm = TRUE)

# save the scaled modelling data set for examining their ranges in a later script
write_csv(all_scaled, "results/scaled_modelling_data_interactions.csv")

################################################################################
# 5) RUN ROBUST MODELS FOR EACH HYPOTHESIS INCLUDING INTERACTION 
# (all site based community variables 
# (not provenance - that is species level - see section 6))
################################################################################

# put the comm values into long format to match the demo vars 
# make nested df
nested_all <- all_scaled %>% 
  dplyr::select(!new_provenance) %>% 
  pivot_longer(cols = site_time_rich:mean_cover, names_to = "Comm_var", values_to = "Comm_value") %>%
  group_by(time_period_cover_change, experiment, Demo_trait, Comm_var, hypoth) %>% 
  nest()

# clean up data frame by removing rows with NA/ INF etc.
nested_clean <- 
  nested_all %>%
  mutate(data = purrr::map(data, ~ .x %>% 
                      filter(if_all(everything(), ~ is.finite(Demo_value))) %>%  # keep only finite values
                      drop_na(Demo_value)))


# default uses the Tukey's bisquare (biweight) function
fit_rlm_safe_lmrob <- 
  purrr::possibly(~ lmrob(RE ~ Demo_value * Comm_value, method = "MM",
                          control = lmrob.control(psi = "bisquare"), data = .x),
                  otherwise = NULL)


# build the augment function manually
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

# run the functions on the nested df
comm_results <- 
  nested_clean %>%
  mutate(mod = purrr::map(data, fit_rlm_safe_lmrob)) %>%
  mutate(results = purrr::map(mod, ~ if (is.null(.x)) tibble() else broom::tidy(.x))) %>%
  mutate(aug = purrr::map(mod, ~ if (is.null(.x)) tibble() else augment_lmrob(.x))) %>%
  mutate(r2 = purrr::map_dbl(mod, ~ if (is.null(.x)) NA_real_ else summary(.x)$r.squared)) %>%
  mutate(confint = purrr::map(mod, ~ if (is.null(.x)) tibble() else as_tibble(confint(.x), rownames = "term")))

# save to disk for model checking and plotting
saveRDS(comm_results, file = "results/robust_model_comm_results_T0_scaled_not_logged.rds")


# extract the results of these models 
# filter out the estimate and p-value for the New_demo_value term only
comm_results <- comm_results %>%
  dplyr::select(time_period_cover_change, experiment, Demo_trait, Comm_var, results, hypoth) %>%
  unnest(results) %>%
  filter(term != "(Intercept)")

# quickly see if there are any significant interaction terms
sig_p <- comm_results %>% filter(p.value < 0.05 & term == "Demo_value:Comm_value")
sig_p # there are quite a few but the estimates are small for a lot of them

# make a variable that refers to each model
comm_results <- comm_results %>% 
  mutate(model = paste0(Demo_trait, "_", Comm_var, "_", time_period_cover_change, "_", 
                        experiment))
comm_results %>%
  group_by(model) %>%
  mutate(flag = p.value < 0.05)

# write_csv(comm_results, "results/Robust_OLS_comm_vars_results_T0_scaled.csv")

# only look at models with significant interaction terms 
# there are some which should be investigated further.
comm_results %>%
  group_by(model) %>%
  filter(term == "Demo_value:Comm_value" & p.value < 0.05) 

################################################################################
# 6) TEST THE SPECIES PROVENANCE DATA INTERACTION
################################################################################

# create a nested df
nested_all <- all %>% 
  dplyr::select(Taxon, Taxon_site, RE, model, Demo_trait, time_period_cover_change,
                experiment, Demo_value, New_demo_trait, hypoth, year_trt, 
                experiment, new_provenance) %>% 
  group_by(time_period_cover_change, experiment, Demo_trait) %>% 
  nest()

# clean up the NA values etc..
nested_clean <- 
  nested_all %>%
  mutate(data = purrr::map(data, ~ .x %>% 
                      filter(if_all(everything(), ~ is.finite(Demo_value))) %>%  # keep only finite values
                      drop_na(Demo_value) %>%
                      filter(if_all(everything(), ~ !is.na(new_provenance)))))

# default uses the Tukey's bisquare (biweight) function
fit_rlm_safe_lmrob <- 
  purrr::possibly(~ lmrob(RE ~ Demo_value * new_provenance, method = "MM",
                          control = lmrob.control(psi = "bisquare"), data = .x),
                  otherwise = NULL)


# build the augment function manually
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

# extra model parameters add
# Now run the models / get results / 
prov_results <- 
  nested_clean %>%
  mutate(mod = purrr::map(data, fit_rlm_safe_lmrob)) %>%
  mutate(results = purrr::map(mod, ~ if (is.null(.x)) tibble() else broom::tidy(.x))) %>%
  mutate(aug = purrr::map(mod, ~ if (is.null(.x)) tibble() else augment_lmrob(.x))) %>%
  mutate(r2 = purrr::map_dbl(mod, ~ if (is.null(.x)) NA_real_ else summary(.x)$r.squared)) %>%
  mutate(confint = purrr::map(mod, ~ if (is.null(.x)) tibble() else as_tibble(confint(.x), rownames = "term")))


# save to disk for model checking and plotting
saveRDS(prov_results, file = "results/robust_model_provenance_results_T0_scaled_not_logged.rds")

# pick out the results
prov_results <- prov_results %>%
  dplyr::select(time_period_cover_change, experiment, Demo_trait, results) %>%
  unnest(results) %>%
  filter(term != "(Intercept)")

write_csv(comm_results, "results/robust_OLS_provenance_results_T0.csv")

# only look at models with significant interaction terms 
# there are some which should be investigated further
comm_results %>%
  mutate(model = paste0(Demo_trait," _", time_period_cover_change, "_", 
                        experiment)) %>%
  group_by(model) %>%
  filter(term == "Demo_value:new_provenanceNAT" & p.value < 0.05) 
# not many here - so probably not worth pursuing

################################################################################

