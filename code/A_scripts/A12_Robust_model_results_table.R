# R Thornley
# 12/09/2025
# Project: P1_COMPADRE_DRAGNET
# make a table of the Robust OLS results for the appendices
# for the univariate OLS results
# for the community results

################################################################################
# BIVARIATE OLS TABLE
################################################################################

# import the results table
dat <- readRDS("results/Robust_models_Sept/robust_models_3.rds")
head(dat)

# extract the results of these models 
# filter out the estimate and p-value for the New_demo_value term only
dat <- dat %>%
  filter(model == "Ordbeta") %>%
  dplyr::select(time_period, experiment, Demo_trait, results) %>%
  unnest(results) %>%
  filter(model == "Ordbeta") %>%
  filter(experiment %in% c("DIST", "NPK")) %>%
  filter(term == "Demo_value") %>%
  ungroup() %>%
  dplyr::select(-term)

# now we need to filter this table by our hypotheses
dat <- dat %>%
    mutate(
      New_demo_trait = case_when(
        Demo_trait == "T_generation"    ~ "Generation time",
        Demo_trait == "R0"              ~ "Net reproductive rate",
        Demo_trait == "percapita_repro" ~ "Per capita reproduction",
        Demo_trait == "Lmean"           ~ "Mean life expectancy",
        Demo_trait == "age_repro"       ~ "Age at first reproduction",
        Demo_trait == "Lambda"          ~ "Population growth rate",
        Demo_trait == "Reactivity"      ~ "Reactivity",
        Demo_trait == "FirstStepAtt"    ~ "First step attenuation",
        Demo_trait == "Lmax"            ~ "Maximum length of life",
        TRUE                            ~ Demo_trait)) %>%
        mutate(
      hypoth = case_when(
        New_demo_trait == "Net reproductive rate"     & experiment == "DIST" ~ "H1",
        New_demo_trait == "Per capita reproduction"   & experiment == "DIST" ~ "H1",
        New_demo_trait == "Generation time"           & experiment == "DIST" ~ "H2",
        New_demo_trait == "Age at first reproduction" & experiment == "DIST" ~ "H2",
        New_demo_trait == "Reactivity"                & experiment == "DIST" ~ "H3",
        New_demo_trait == "First step attenuation"    & experiment == "DIST" ~ "H3",
        New_demo_trait == "Population growth rate"    & experiment == "DIST" ~ "H4",
        New_demo_trait == "Maximum length of life"    & experiment == "NPK"  ~ "H5",
        New_demo_trait == "Mean life expectancy"      & experiment == "NPK"  ~ "H5",
        TRUE ~ "none"
      ))

dat <- dat %>% filter(!hypoth =="none")

# re-order the table by hypotheses:
dat <- 
  dat %>% 
  ungroup() %>%
  arrange(hypoth) %>% 
  dplyr::select(hypoth, experiment, time_period, New_demo_trait, estimate, std.error, statistic, p.value) %>%
  rename("Hypotheses" = hypoth, "Experiment" = experiment, "Time period" = time_period, 
         "Demographic metric" = New_demo_trait, "Estimate" = estimate,  "Standard Error" = std.error, "Statitsic" = statistic, "P-value" = p.value)
dat
write_csv(dat, "figures/robust_univariate_OLS_table_appendix.csv")   

################################################################################
# INTERACTION TERMS TABLE
################################################################################

# import the results for the community interaction regressions
comm <- read_csv("results/Robust_OLS_comm_vars_results.csv")
comm <- comm %>% dplyr::select(-model)
prov <- read_csv("results/Robust_OLS_provenance_results.csv")
prov <- prov %>% mutate(Comm_var = "Species provenance", hypoth = "H6")

dat <- rbind(comm, prov)

# now we need to filter this table by our hypotheses
dat <- dat %>%
  mutate(
    New_demo_trait = case_when(
      Demo_trait == "T_generation"    ~ "Generation time",
      Demo_trait == "R0"              ~ "Net reproductive rate",
      Demo_trait == "percapita_repro" ~ "Per capita reproduction",
      Demo_trait == "Lmean"           ~ "Mean life expectancy",
      Demo_trait == "age_repro"       ~ "Age at first reproduction",
      Demo_trait == "Lambda"          ~ "Population growth rate",
      Demo_trait == "Reactivity"      ~ "Reactivity",
      Demo_trait == "FirstStepAtt"    ~ "First step attenuation",
      Demo_trait == "Lmax"            ~ "Maximum length of life",
      TRUE                            ~ Demo_trait)) %>%
  mutate(
    hypoth = case_when(
      New_demo_trait == "Net reproductive rate"     & experiment == "DIST" ~ "H1",
      New_demo_trait == "Per capita reproduction"   & experiment == "DIST" ~ "H1",
      New_demo_trait == "Generation time"           & experiment == "DIST" ~ "H2",
      New_demo_trait == "Age at first reproduction" & experiment == "DIST" ~ "H2",
      New_demo_trait == "Reactivity"                & experiment == "DIST" ~ "H3",
      New_demo_trait == "First step attenuation"    & experiment == "DIST" ~ "H3",
      New_demo_trait == "Population growth rate"    & experiment == "DIST" ~ "H4",
      New_demo_trait == "Maximum length of life"    & experiment == "NPK"  ~ "H5",
      New_demo_trait == "Mean life expectancy"      & experiment == "NPK"  ~ "H5",
      TRUE ~ "none"
    ))

unique(dat$Comm_var)
# filter out the comm vars we actually want 
dat <- 
  dat %>% 
  filter(Comm_var %in% c("site_time_rich", "site_time_invsimp", "mean_cover", 
                         "Species provenance", "variability"))

# Tidy up the community var names
dat <- dat %>%
  mutate(
    Comm_var = case_when(
      Comm_var == "site_time_rich"    ~ "Species richness",
      Comm_var ==   "site_time_invsimp"   ~ "Simpson's diversity",
      Comm_var == "mean_cover" ~ "Total plot cover",
      Comm_var == "variability" ~ "Functional diversity",
      TRUE                            ~ Comm_var))

dat <- dat %>%
  mutate(
    term = case_when(
      term == "Demo_value" ~ "Demographic metric",
      term == "Comm_value" ~ "Community variable",
      term == "Demo_value:Comm_value" ~ "Interaction",
      TRUE ~ term
    )
  )


dat <- dat %>%
  ungroup() %>%
  mutate(model = paste0(Demo_trait," _", time_period, "_", 
                        experiment, "_", Demo_trait, "_", Comm_var)) %>%
  group_by(model) %>%
  filter(any(term == "Interaction" & p.value < 0.05, na.rm = TRUE)) %>%
  ungroup()

dat <- dat %>% mutate(hypoth = factor(hypoth, levels = c("H1", "H2", "H3", "H4", "H5")))
dat <- dat %>% dplyr::select(hypoth, time_period, experiment, New_demo_trait, Comm_var, 
                      term, model, estimate, std.error, statistic, p.value) %>%
  arrange(hypoth, time_period, experiment, New_demo_trait, Comm_var, term)

names(dat)
unique(dat$hypoth)
unique(dat$Comm_var)

dat <- dat %>% 
  rename("Hypothesis" = hypoth, "Experiment" = experiment, "Time period" = time_period, 
               "Demographic metric" = New_demo_trait, "Community variable" = Comm_var, "Term" = term,
               "Estimate" = estimate,  "Standard Error" = std.error, 
         "Statitsic" = statistic, "P-value" = p.value) %>%
  dplyr::select(-model)
names(dat)
# round columns to 4 dp
dat <- dat %>% mutate(across(Estimate:`P-value`, ~ if (is.numeric(.x)) round(.x, 4) else .x))

dat <- dat %>% 
  mutate(Experiment = case_when(
  Experiment == "DIST" ~ "Disturbance",
  Experiment == "NPK" ~ "Nutrient Addition",
  TRUE ~ Experiment))

# export results as a table for the appendices
write_csv(dat, "figures/robust_community_OLS_table_appendix.csv")

