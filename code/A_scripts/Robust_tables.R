# R Thornley
# 07/11/2025
# Project: P1_COMPADRE_DRAGNET
# Script: Appendices tables for Robust models

################################################################################
# Robust OLS table results
################################################################################

# import the model results 
mods <- readRDS(file = "results/robust_model_results.rds")
names(mods)
# extract the key parameters of the results for a table
estimates <- mods %>%
  filter(model == "Ordbeta") %>% 
  dplyr::select(model, time_period, experiment, Demo_trait, results, r2) %>%
  unnest(c(results)) %>%
  filter(term == "Demo_value")

CIs <- mods %>%
  filter(model == "Ordbeta") %>% 
  dplyr::select(model, time_period, experiment, Demo_trait, confint) %>%
  unnest(c(confint)) %>%
  filter(term == "Demo_value")

both <- estimates %>%
  left_join(CIs) %>%
  rename("CI-2.5%" = '2.5 %', "CI-97.5%" = '97.5 %')

names(both)
both <- 
  both %>% mutate(
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
  ungroup() %>%
  dplyr::select(hypoth, New_demo_trait, time_period, estimate, statistic, p.value, std.error, r2, 11:12) %>%
  rename("Hypothesis" = "hypoth", "Robust regression co-efficient" = "estimate", "Standard error" = "std.error",
         "t-value" = "statistic", "p-value" = "p.value", "Demographic metric" = "New_demo_trait", 
         "Time period" = "time_period") %>%
  arrange(Hypothesis, 'Demographic metric', 'Time period')
names(both)
both <- both %>% mutate(across(4:10, ~ round(.x, 4)))

write_csv(both, "figures/robust_OLS_results.csv")

################################################################################
# Robust interactions table results: comm vars
################################################################################

# import the model results 
comm <- readRDS(file = "results/robust_model_comm_results_T0_scaled_not_logged.rds")
names(comm)
estimates <- comm %>%
  dplyr::select(time_period_cover_change, experiment, Demo_trait, hypoth, Comm_var, results, r2) %>%
  unnest(c(results)) %>%
  filter(term %in% c("Demo_value", "Comm_value", "Demo_value:Comm_value"))

CIs <- comm %>%
  dplyr::select(time_period_cover_change, experiment, Demo_trait, hypoth, Comm_var, confint, r2) %>%
  unnest(c(confint)) %>%
  filter(term %in% c("Demo_value", "Comm_value", "Demo_value:Comm_value"))

both <- estimates %>%
  left_join(CIs) %>%
  rename("CI-2.5%" = '2.5 %', "CI-97.5%" = '97.5 %')

both <- both %>% filter(!Comm_var %in% c("site_time_shan", "site_time_simp", "sigmaonMeans"))
unique(both$Comm_var)

both <- 
  both %>% mutate(
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
    ),
    Comm_var = case_when(
      Comm_var == "variability" ~ "Functional diversity",
      Comm_var == "site_time_rich" ~ "Species richness",
      Comm_var == "site_time_invsimp" ~ "Species diversity",
      Comm_var == "meanSigmas" ~ "Functional specialisation",
      Comm_var == "mean_cover" ~ "Total cover",
      TRUE ~ Comm_var
      ),
    term = case_when(
      term == "Demo_value" ~ "Demographic",
      term == "Comm_value" ~ "Community",
      term == "Demo_value:Comm_value" ~ "Interaction",
      TRUE ~ term
      )
    )%>%
  filter(hypoth != "none") %>%
  ungroup() %>%
  dplyr::select(hypoth, New_demo_trait, Comm_var, time_period_cover_change, term, estimate, statistic, p.value, std.error, r2, 12:13) %>%
  rename("Hypothesis" = "hypoth", "Robust regression coefficient" = "estimate", "Standard error" = "std.error",
         "t-value" = "statistic", "p-value" = "p.value", "Demographic metric" = "New_demo_trait", 
         "Time period" = "time_period_cover_change", "Term" = "term") %>%
  arrange(Hypothesis, 'Demographic metric', 'Time period')

names(both)
both <- both %>% mutate(across(6:12, ~ round(.x, 4)))

write_csv(both, "figures/robust_OLS_interaction_results.csv")

################################################################################
# Robust interactions table results: provenance 
################################################################################

prov <- readRDS(file = "results/robust_model_provenance_results_T0_scaled_not_logged.rds")
names(prov)

estimates <- prov %>%
  dplyr::select(time_period_cover_change, experiment, Demo_trait, results, r2) %>%
  unnest(c(results)) %>%
  filter(term %in% c("Demo_value", "new_provenanceNAT", "Demo_value:new_provenanceNAT"))

CIs <- prov %>%
  dplyr::select(time_period_cover_change, experiment, Demo_trait, confint, r2) %>%
  unnest(c(confint)) %>%
  filter(term %in% c("Demo_value", "new_provenanceNAT", "Demo_value:new_provenanceNAT"))

both <- estimates %>%
  left_join(CIs) %>%
  rename("CI-2.5%" = '2.5 %', "CI-97.5%" = '97.5 %')

both <- 
  both %>% mutate(
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
    ),
    term = case_when(
      term == "Demo_value" ~ "Demographic",
      term == "new_provenanceNAT" ~ "Species provenance",
      term == "Demo_value:new_provenanceNAT" ~ "Interaction",
      TRUE ~ term
    )
  ) %>%
  filter(hypoth != "none") %>%
  ungroup() %>%
  dplyr::select(hypoth, New_demo_trait, time_period_cover_change, term, estimate, statistic, p.value, std.error, r2, 12:13) %>%
  rename("Hypothesis" = "hypoth", "Robust regression coefficient" = "estimate", "Standard error" = "std.error",
         "t-value" = "statistic", "p-value" = "p.value", "Demographic metric" = "New_demo_trait", 
         "Time period" = "time_period_cover_change", "Term" = "term") %>%
  arrange(Hypothesis, 'Demographic metric', 'Time period')

names(both)
both <- both %>% mutate(across(5:9, ~ round(.x, 4)))

write_csv(both, "figures/robust_OLS_interaction_provenance_results.csv")

