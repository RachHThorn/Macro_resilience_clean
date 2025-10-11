# R Thornley
# 12/09/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S6_get_community_metrics_each_time

################################################################################
# Instructions
################################################################################

# 1) Calculate Species richness and diversity at the DRAGNet site level for three time periods
# 2) Species native or non native status at the DRAGNet site level for three time periods
# 3) The average total cover per site at the DRAGNet site level for three time periods

################################################################################
# 1) Species richness and diversity at a site
################################################################################

get_rich_div <- function(file_path, time_period) {

  # read once
  df <- readr::read_csv(file_path, show_col_types = FALSE) %>%
    mutate(unique = paste(site_name, year_trt, trt, block, sep = "_"))
  
  # ---- diversity (per quadrat) ----
  wide <- df %>%
    dplyr::select(unique, max_cover, New_taxon) %>%
    pivot_wider(
      names_from = New_taxon,
      values_from = max_cover,
      values_fn = ~ mean(.x, na.rm = TRUE)
    ) %>%
    replace(is.na(.), 0)
  
  spec_mat <- wide %>% dplyr::select(-unique)
  
  diversity <- wide %>%
    transmute(
      unique,
      shannon = vegan::diversity(as.matrix(spec_mat), index = "shannon",  MARGIN = 1),
      simpson = vegan::diversity(as.matrix(spec_mat), index = "simpson",  MARGIN = 1),
      invsimp = vegan::diversity(as.matrix(spec_mat), index = "invsimpson", MARGIN = 1)
    )
  
  # ---- richness (per quadrat) ----
  richness <- df %>%
    group_by(site_name, year_trt, trt, unique) %>%
    summarise(richness = n_distinct(New_taxon), .groups = "drop")
  
  # join quadrat-level metrics
  all <- richness %>%
    left_join(diversity, by = "unique")
  
  # ---- collapse to site × time × trt means ----
  all <- all %>%
    group_by(site_name, year_trt, trt) %>%
    summarise(
      site_time_rich    = mean(richness, na.rm = TRUE),
      site_time_simp    = mean(simpson,  na.rm = TRUE),
      site_time_shan    = mean(shannon,  na.rm = TRUE),
      site_time_invsimp = mean(invsimp,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(time_period = time_period) %>%
    dplyr::select(site_name, year_trt, time_period, trt, site_time_rich:site_time_invsimp)
  
  # (optional) count combos if you want to inspect
  # nrow(all)
  
  return(all)
}

# Suppose you have:
file_paths <- c(
  "results/DRAGNet_T0_T1_all.csv",
  "results/DRAGNet_T0_T2_all.csv",
  "results/DRAGNet_T0_T3_all.csv"
)

time_periods <- c("T0-T1", "T0-T2", "T0-T3")

# Apply function to each file/time combination
site_div_list <- map2(file_paths, time_periods, get_rich_div)

# Combine into one tibble
site_div_all <- bind_rows(site_div_list)

# export data
write_csv(site_div_all, "results/diversity_metrics_dragnet_.csv")

#################################################################################
# 2) Species native or non native status at a site
################################################################################

raw_path    <- "data/full-cover-2025-07-16.csv"   # raw DRAGNet with provenance
common_path <- "results/common_species_drag_comp.csv"  # contains column 'x' of focal taxa
out_path    <- "results/native_status.csv"

# Non-taxa labels to drop
none_taxa <- c(
  "Fungi", "Other_litter", "Other_standing_water", "Ground", "Bryophyte",
  "Other_animal_diggings", "Other_woody_overstory", "Lichen",
  "Other_animal_digging", "Other_animal_droppings"
)

dat_raw <- readr::read_csv(raw_path, show_col_types = FALSE)
common  <- readr::read_csv(common_path, show_col_types = FALSE) |> pull(x)

# Title case, underscores, remove non-taxa and ambiguous entries; collapse varieties to species
native <- dat_raw %>%
  mutate(
    New_taxon = str_to_sentence(Taxon),
    New_taxon = str_replace_all(New_taxon, " ", "_")
  ) %>%
  filter(
    !str_detect(New_taxon, "\\.sp"),         # drop genus-level "X sp."
    !str_detect(New_taxon, "_x_"),           # drop hybrids
    !str_detect(New_taxon, "Unknown"),       # drop unknowns
    !New_taxon %in% none_taxa                # drop non-taxa
  ) %>%
  mutate(
    New_taxon = case_when(
      New_taxon == "Helianthemum_nummularium_var._Grandiflorum" ~ "Helianthemum_nummularium",
      New_taxon == "Mimosa_quadrivalvis_var._Platycarpa"        ~ "Mimosa_quadrivalvis",
      New_taxon == "Sebaea_sedoides_var._Schoenlandii"          ~ "Sebaea_sedoides",
      TRUE ~ New_taxon
    )
  ) %>%
  dplyr::select(site_name, New_taxon, local_provenance)

# Keep only focal species (compadre/dragnet set) & one row per site × taxon
all_site_taxa <- native %>%
  filter(New_taxon %in% common) %>%
  distinct(site_name, New_taxon, .keep_all = TRUE)

# Normalise provenance values
# Convert "NULL" / "UNK" to NA_character_ so we can fill from overrides
all_site_taxa <- all_site_taxa %>%
  mutate(local_provenance = na_if(local_provenance, "NULL"),
         local_provenance = na_if(local_provenance, "UNK"))

# Manual overrides for site × taxon with known status
# Use a small lookup tibble instead of a long case_when
override_prov <- tibble::tribble(
  ~site_name,                                      ~New_taxon,               ~override_prov,
  "Bayreuth DRAGNet",                              "Fragaria_vesca",         "NAT",
  "Blandy Experimental Farm",                      "Lactuca_serriola",       "INT",
  "CEREEP - Ecotron IDF",                          "Cerastium_fontanum",     "NAT",
  "CEREEP - Ecotron IDF",                          "Erigeron_canadensis",    "INT",
  "CEREEP - Ecotron IDF",                          "Alliaria_petiolata",     "NAT",
  "CEREEP - Ecotron IDF",                          "Achillea_millefolium",   "NAT",
  "CEREEP - Ecotron IDF",                          "Myosotis_ramosissima",   "NAT",
  "Hazelrigg",                                     "Ranunculus_acris",       "NAT",
  "Piedmont Prairie",                              "Verbascum_thapsus",      "INT",
  "University of Michigan Biological Station",     "Tragopogon_dubius",      "INT",
  "University of Wyoming - SAREC",                 "Tragopogon_dubius",      "INT",
  "University of Wyoming - SAREC",                 "Bromus_tectorum",        "INT",
  "University of Wyoming - SAREC",                 "Lactuca_serriola",       "INT",
  "Wytham Woods",                                  "Ranunculus_bulbosus",    "NAT"
)

# Join overrides; use override when present, else original value
new <- all_site_taxa %>%
  left_join(override_prov, by = c("site_name", "New_taxon")) %>%
  mutate(new_provenance = dplyr::coalesce(override_prov, local_provenance)) %>%
  dplyr::select(site_name, New_taxon, local_provenance, new_provenance)

# Quick checks
message("Distinct taxa: ", dplyr::n_distinct(new$New_taxon))   # e.g., 74
message("Distinct sites: ", dplyr::n_distinct(new$site_name))  # e.g., 43

# Distribution before/after fixes (optional plots)
ggplot(all_site_taxa, aes(local_provenance)) + geom_bar() +
  labs(title = "Provenance (raw)", x = "local_provenance", y = "count")

ggplot(new, aes(new_provenance)) + geom_bar() +
  labs(title = "Provenance (after overrides)", x = "new_provenance", y = "count")

write_csv(new, out_path)

################################################################################
# 3) The average total cover per site per treatment and time point 
################################################################################
# get the sum of cover values for total cover per quadrat
# then take a site, time, experiment mean

get_total_cover <- function(file_path, time_period) {
  data <- read.csv(file_path)
  
  data <- data %>% 
    mutate(unique = paste(site_name, year_trt, trt, block, sep = "_")) %>% 
    group_by(unique) %>%
    mutate(total_q_cover = sum(max_cover, na.rm = TRUE)) %>%  
    group_by(site_name, year_trt, trt) %>%
    summarise(mean_cover = mean(total_q_cover, na.rm = TRUE), .groups = "drop") %>%
    mutate(experiment = case_when(
      trt == "Disturbance"       ~ "DIST",
      trt == "NPK+Disturbance"   ~ "INTER", 
      trt == "NPK"               ~ "NPK", 
      trt == "Control"           ~ "CONTROL",
      TRUE                       ~ trt
    ),
    time_period = time_period)
  
  return(data)
}

# file paths and time periods
file_paths <- c(
  "results/DRAGNet_T0_T1_all.csv", 
  "results/DRAGNet_T0_T2_all.csv", 
  "results/DRAGNet_T0_T3_all.csv"
)

time_periods <- c("T0-T1", "T0-T2", "T0-T3")

# apply function to all inputs with purrr::map2
cover_results <- purrr::map2(file_paths, time_periods, get_total_cover)

# bind all results together
cover_all <- bind_rows(cover_results)

# save combined file if needed
write_csv(cover_all, "results/total_cover.csv")



