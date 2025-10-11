# R Thornley
# 16/07/2025
# Project: P1_COMPADRE_DRAGNET
# Script: Get summary data for appendices tables

library(tidyverse)

################################################################################
# PART 1: Number of taxa in each analysis
################################################################################

# read in the OLS modelling data
dat <- read_csv("results/OLS_modelling_master_data.csv")
names(dat)
# tidy some names and create hypotheses numbers
dat <- dat %>%
    filter(model == "Ordbeta", experiment %in% c("DIST", "NPK")) %>%
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
        New_demo_trait == "Maximum length of life"    & experiment == "NPK"  ~ "H5",
        New_demo_trait == "Mean life expectancy"      & experiment == "NPK"  ~ "H5",
        TRUE ~ "none"
      )
    ) %>%
    filter(hypoth != "none") 

# calculate the number of taxa in each group
Nos_taxa <- dat %>%
  group_by(time_period, experiment, Demo_trait) %>%
  mutate(n_taxa = n_distinct(Taxon)) %>%
  ungroup() %>%
  dplyr::select(time_period, experiment, New_demo_trait, hypoth, n_taxa) %>%
  distinct() %>%
  rename(
    `Time period`        = time_period,
    `Experiment`         = experiment,
    `Demographic metric` = New_demo_trait,
    `Hypothesis`         = hypoth,
    `Number of taxa`     = n_taxa
  ) %>%
  dplyr::select(Hypothesis, `Demographic metric`, Experiment, `Time period`, `Number of taxa`)

# Change the factors so they are in order
Nos_taxa$Hypothesis <- factor(Nos_taxa$Hypothesis,
                    levels = c("H1", "H2", "H3", "H4", "H5"))
Nos_taxa$`Time period` <- factor(Nos_taxa$`Time period`,
                                 levels = c("T0-T1", "T0-T2", "T0-T3"))

Nos_taxa <- Nos_taxa %>% arrange(Hypothesis, `Time period`, Experiment, 'Demographic metric')

# This needs exporting or joining with the output from part 2

###############################################################################
# PART 2: Number of matrices in each analysis
################################################################################

# read in the direct output from the S4 script
demo <- read_csv("results/all_COMPADRE_metrics.csv") %>% 
  filter(DRAGNet == TRUE) %>% 
  rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var", "Demo_value" = "value") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  dplyr::select(Taxon, MatrixID, Demo_trait, Demo_value)

# get the list of species that we use in the final analysis
mod_taxa <- readRDS("results/List_taxa_OLS_mods.R")

# filter the demo data with this list
demo <- demo %>% filter(Taxon %in% mod_taxa)

# number of unique matrices used in this study
demo %>% summarise(n_unique = n_distinct(MatrixID))
# 340 matrices used in the analysis altogether


# read in effect sizes from the whole of DRAGNet 
taxa <- 
  read_csv("results/RE_SE_Taxon_all_DRAGNet.csv") %>%
  filter(str_detect(group_var, "New_taxon")) %>%
  filter(group %in% mod_taxa) %>%
  mutate(time_period = case_when(time_period == 1 ~ "T0-T1",
                                 time_period == 2 ~ "T0-T2",
                                 time_period == 3 ~ "T0-T3")) %>%
  dplyr::select(group, value, se, model, time_period, experiment) %>%
  rename(Taxon = group, RE = value, RE_se = se) %>%
  filter(model == "Ordbeta")
unique(taxa$Taxon) # 40 unique taxa here

# join both data sets
both <- taxa %>% left_join(demo, by = "Taxon")
names(both)
# check this data look reasonable by counting nos of matrices as before
both %>% group_by(experiment, time_period, model, Demo_trait, Taxon) %>% count()

# nos of matrices per analysis
both <- both %>% group_by(experiment, time_period, model, Demo_trait) %>% count()

# we now need to filter this for the rows that are associated with our hypotheses
names(both)
both <- both %>%
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
      New_demo_trait == "Maximum length of life"    & experiment == "NPK"  ~ "H5",
      New_demo_trait == "Mean life expectancy"      & experiment == "NPK"  ~ "H5",
      TRUE ~ "none"
    )
  ) %>%
  filter(hypoth != "none") 
names(both)

# now tidy this table so it can be joined with the data from 
both <- 
  both %>%
  ungroup() %>%
  rename(
  `Time period`        = time_period,
  `Experiment`         = experiment,
  `Demographic metric` = New_demo_trait,
  `Hypothesis`         = hypoth,
  `Number of matrices`     = n
) %>%
  dplyr::select(Hypothesis, `Demographic metric`, Experiment, `Time period`, `Number of matrices`)

# Change the factors so they are in order
both$Hypothesis <- factor(both$Hypothesis,
                              levels = c("H1", "H2", "H3", "H4", "H5"))
both$`Time period` <- factor(both$`Time period`,
                                 levels = c("T0-T1", "T0-T2", "T0-T3"))

Nos_mat <- both %>% arrange(Hypothesis, `Time period`, Experiment, 'Demographic metric')

################################################################################
# join the tables from the two first sections
###############################################################################

Table_1 <- Nos_taxa %>% left_join(Nos_mat)
write_csv(Table_1, "results/A_x_number_species_matrices_per_analysis.csv")

###############################################################################
# PART 3: Number of matrices for each species and demography trait
################################################################################

# read in the direct output from the S4 script
nos_mat_taxa <- read_csv("results/all_COMPADRE_metrics.csv") %>% 
  filter(DRAGNet == TRUE) %>% 
  rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var", "Demo_value" = "value") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  drop_na()

# n shows us the number of matrices we are using for each species
nos_mat_taxa <-
  nos_mat_taxa %>% group_by(Taxon, Demo_trait) %>% count() %>% rename('Number of matrices' = n)

###############################################################################
# PART 4: Number of DRAGNet sites for each species
################################################################################

Drag <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv")
names(Drag)
unique(Drag$group_var)

# read in data; filter by the RE of taxon site, split the name into species and site
Drag <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv") %>%
  filter(group_var == "taxon_site") %>%
  filter(model == "Ordbeta") %>%
  separate(group, into = c("sp1", "sp2", "site_name"), sep = "_") %>%
  mutate(Taxon = paste0(sp1,"_", sp2)) %>%
  dplyr::select(Taxon, site_name, model, time_period, experiment)

# filter by the taxa we are studying 
mod_taxa <- readRDS("results/List_taxa_OLS_mods.R")
Drag <- Drag %>% filter(Taxon %in% mod_taxa)

# Now count the number of sites per taxa
Drag <- Drag %>%
  group_by(Taxon) %>%
  summarise(n_sites = n_distinct(site_name), .groups = "drop") %>%
  rename('Number of DRAGNet sites' = n_sites)

names(Drag)
names(nos_mat_taxa)

Table_2 <- Drag %>% 
  left_join(nos_mat_taxa) %>% 
  rename('Number of COMPADRE matrices' = 'Number of matrices') %>%
  dplyr::select(Taxon, Demo_trait, 'Number of DRAGNet sites', 'Number of COMPADRE matrices')
# Tidy the table names
Table_2 <-
  Table_2 %>%
  mutate(
    'Demographic metric' = case_when(
      Demo_trait == "T_generation"    ~ "Generation time",
      Demo_trait == "R0"              ~ "Net reproductive rate",
      Demo_trait == "percapita_repro" ~ "Per capita reproduction",
      Demo_trait == "Lmean"           ~ "Mean life expectancy",
      Demo_trait == "age_repro"       ~ "Age at first reproduction",
      Demo_trait == "Lambda"          ~ "Population growth rate",
      Demo_trait == "Reactivity"      ~ "Reactivity",
      Demo_trait == "FirstStepAtt"    ~ "First step attenuation",
      Demo_trait == "Lmax"            ~ "Maximum length of life",
      TRUE                            ~ Demo_trait
    )) %>%
  mutate(Taxon = str_replace(Taxon, "_", " ")) %>%
  dplyr::select(Taxon, 'Demographic metric', 'Number of DRAGNet sites', 
                'Number of COMPADRE matrices')
write_csv(Table_2, "results/A_x_summary_table_nos_sites_matrices_per_taxon.csv")

    