library(dplyr)
library(purrr)
library(readr)
library(tidyr)

follow_ups <- 1:3

arms <- tibble(
  experiment = c("DIST", "DIST", "NPK", "NPK"),
  trt = c("Control", "Disturbance", "Control", "NPK")
)

# Demographic metrics tested for each treatment comparison.
metrics <- tibble(
  experiment = rep(c("DIST", "NPK"), c(7, 2)),
  Demo_trait = c(
    "R0", "percapita_repro", "T_generation", "age_repro",
    "Reactivity", "FirstStepAtt", "Lambda", "Lmax", "Lmean"
  )
)

# Pool the species that entered the fitted demographic regressions
# within each treatment comparison and time period.
regression_species <- readRDS("results/robust_model_results.rds") %>%
  ungroup() %>%
  filter(model == "Ordbeta", !map_lgl(mod, is.null)) %>%
  semi_join(metrics, by = c("experiment", "Demo_trait")) %>%
  transmute(experiment, time_period, Taxon = map(data, "Taxon")) %>%
  unnest(Taxon) %>%
  distinct()

# Keep cover records for those species in the corresponding comparison.
# Control records can contribute to both comparisons.
cover <- map_dfr(follow_ups, function(i) {
  read_csv(
    paste0("results/DRAGNet_T0_T", i, "_all.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      time_period = paste0("T0-T", i),
      site_name, block, trt, year_trt, Taxon = New_taxon
    )
}) %>%
  inner_join(arms, by = "trt", relationship = "many-to-many") %>%
  semi_join(regression_species, by = c("experiment", "time_period", "Taxon"))

# Use original plot IDs to check that each arm had at least one
# plot surveyed at both baseline and follow-up.
surveys <- read_csv("data/full-cover-2025-09-11.csv", show_col_types = FALSE) %>%
  filter(trt %in% arms$trt, !is.na(plot)) %>%
  distinct(site_name, block, trt, plot, year_trt)

surveyed_arms <- map_dfr(follow_ups, function(i) {
  surveys %>%
    filter(year_trt %in% c(0, i)) %>%
    group_by(site_name, block, trt, plot) %>%
    filter(n_distinct(year_trt) == 2L) %>%
    ungroup() %>%
    distinct(site_name, block, trt) %>%
    mutate(time_period = paste0("T0-T", i))
})

# A matched block has four arm-by-time combinations in the filtered
# cover data: control and treatment at baseline and follow-up.
matched_blocks <- cover %>%
  distinct(experiment, time_period, site_name, block, trt, year_trt) %>%
  semi_join(surveyed_arms, by = c("time_period", "site_name", "block", "trt")) %>%
  count(experiment, time_period, site_name, block, name = "arm_time_cells") %>%
  filter(arm_time_cells == 4L) %>%
  select(-arm_time_cells)

# The site count includes every site contributing cover data
# for the selected species, including sites without a matched block.
replication_table <- cover %>%
  summarise(
    unique_species = n_distinct(Taxon),
    unique_sites = n_distinct(site_name),
    .by = c(experiment, time_period)
  ) %>%
  left_join(
    count(matched_blocks, experiment, time_period, name = "matched_blocks"),
    by = c("experiment", "time_period")
  ) %>%
  mutate(matched_blocks = coalesce(matched_blocks, 0L)) %>%
  select(experiment, time_period, matched_blocks, unique_species, unique_sites) %>%
  arrange(experiment, time_period)

# Count each matched block once per treatment across all follow-ups.
matched_blocks_across_time <- matched_blocks %>%
  distinct(experiment, site_name, block) %>%
  count(experiment, name = "matched_blocks")

# Count unique species, sites and matched site–blocks across all six rows.
pooled_totals <- cover %>%
  summarise(
    unique_species = n_distinct(Taxon),
    unique_sites = n_distinct(site_name),
    unique_site_blocks = n_distinct(matched_blocks$site_name,
                                    matched_blocks$block)
  )

replication_table
matched_blocks_across_time
pooled_totals



#####################################################

species_counts <- cover %>%
  summarise(
    n_sites = n_distinct(site_name),
    n_blocks = n_distinct(site_name, block),
    .by = c(experiment, time_period, Taxon)
  )

species_matched_counts <- cover %>%
  distinct(experiment, time_period, Taxon, site_name, block, trt, year_trt) %>%
  semi_join(
    surveyed_arms,
    by = c("time_period", "site_name", "block", "trt")
  ) %>%
  count(
    experiment, time_period, Taxon, site_name, block,
    name = "arm_time_cells"
  ) %>%
  filter(arm_time_cells == 4L) %>%
  summarise(
    n_matched_sites = n_distinct(site_name),
    n_matched_blocks = n(),
    .by = c(experiment, time_period, Taxon)
  )

per_species_estimate <- species_counts %>%
  left_join(
    species_matched_counts,
    by = c("experiment", "time_period", "Taxon")
  ) %>%
  mutate(across(c(n_matched_sites, n_matched_blocks), ~ replace_na(.x, 0L))) %>%
  arrange(experiment, time_period, Taxon)

per_species_estimate


supplementary_table <- per_species_estimate %>%
  mutate(
    Species = stringr::str_replace_all(Taxon, "_", " "),
    `Blocks (sites)` = paste0(n_blocks, " (", n_sites, ")")
  ) %>%
  select(Species, experiment, time_period, `Blocks (sites)`) %>%
  pivot_wider(
    names_from = c(experiment, time_period),
    names_glue = "{experiment} {time_period}",
    values_from = `Blocks (sites)`,
    values_fill = "—",
    names_sort = TRUE
  ) %>%
  arrange(Species)

supplementary_table
write_csv(supplementary_table, "results/rep_table_row_1_supp_info.csv")

# Suggested table note. Entries give the number of distinct DRAGNet blocks, 
# with the number of sites in parentheses, contributing cover data to each species estimate. 
# A dash indicates that the species did not enter any fitted demographic regression 
# for that comparison and time period.

# make the demographic metrics

logged_traits <- c(
  "R0", "percapita_repro", "age_repro", "Lmean", "Lmax",
  "T_generation", "Reactivity", "FirstStepAtt"
)

fitted_species_metrics <- readRDS("results/robust_model_results.rds") %>%
  ungroup() %>%
  filter(model == "Ordbeta", !map_lgl(mod, is.null)) %>%
  semi_join(metrics, by = c("experiment", "Demo_trait")) %>%
  transmute(Demo_trait, Taxon = map(data, "Taxon")) %>%
  unnest(Taxon) %>%
  distinct()

matrix_records <- read_csv(
  "results/all_COMPADRE_metrics.csv",
  show_col_types = FALSE
) %>%
  filter(DRAGNet, !is.na(MatrixID), is.finite(value)) %>%
  transmute(
    Taxon = stringr::str_replace_all(SpeciesAccepted, " ", "_"),
    Demo_trait = demo_var,
    MatrixID,
    value
  ) %>%
  filter(!Demo_trait %in% logged_traits | value > 0) %>%
  semi_join(fitted_species_metrics, by = c("Taxon", "Demo_trait"))

matrix_total <- matrix_records %>%
  summarise(
    n_matrices = n_distinct(MatrixID),
    n_species = n_distinct(Taxon)
  )

matrices_by_metric <- matrix_records %>%
  summarise(
    n_matrices = n_distinct(MatrixID),
    n_species = n_distinct(Taxon),
    .by = Demo_trait
  ) %>%
  arrange(Demo_trait)

# Supplementary table with one row per species and one column per metric.
matrix_supplement <- matrix_records %>%
  summarise(n_matrices = n_distinct(MatrixID), .by = c(Taxon, Demo_trait)) %>%
  mutate(Species = stringr::str_replace_all(Taxon, "_", " ")) %>%
  select(Species, Demo_trait, n_matrices) %>%
  pivot_wider(
    names_from = Demo_trait,
    values_from = n_matrices,
    names_sort = TRUE
  ) %>%
  arrange(Species)

matrix_total
matrices_by_metric
matrix_supplement

write_csv(matrix_supplement, "results/rep_table_row_9_supp_info.csv")
