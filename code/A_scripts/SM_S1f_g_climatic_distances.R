# R Thornley 
# 09/12/2025
# Project: P1_COMPADRE_DRAGNET
# Appendices data overlap: distance metrics for climate data and re-run OLS models
# to perform a sensitivity analysis

rm(list= ls())

library(tidyverse)
library(purrr)
library(robustbase) # runs robust OLS models
library(mixedup)
library(Rcompadre) # package to interact with COMPADRE database
library(popbio) # basic MPM analyses
library(popdemo) # advanced MPM transient measures
library(Rage) # life history metrics from MPMs

################################################################################
# Instructions
################################################################################

# 1) Get locations for final modelling data
# 2) READ IN THE DATA - extract the climate variables for DRAGNet and COMPADRE
# 3) CALCULATE DISTANCE METRICS
# 4) VISUALISE DISTANCES per species
# 5) Examine mean and median distances
# 6) Load and tidy demographic and species cover change (random effects) data
# 7) Create and tidy a nested df for modelling and apply robust regressions
# 8) VISUALISE (ROUGH) results of the models / filter for significant terms

################################################################################
# 1) Get locations for final modelling data
###############################################################################

# compadre metrics
demo <- read_csv("results/all_COMPADRE_metrics.csv")
# filter for only the overlap species found in the initial text search
demo <- demo %>% filter(DRAGNet == TRUE)
# get a vector of species names in the correct format
demo_species <- demo %>% 
  rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var", "Demo_value" = "value") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  pull(Taxon) %>%
  unique()
demo_species # 42 species 

# two stages of GLMMs modelling outputs
# 1) species level estimates
GLMM_sp <- read_csv("results/RE_SE_Taxon_all_DRAGNet.csv")
names(GLMM_sp)
# filter this for the compadre names
GLMM_sp <- GLMM_sp %>% filter(group %in% demo_species) 
# what are the unique species names in this data set
GLMM_sp_species <- GLMM_sp %>% pull(group) %>% unique() # 41 species

# 2) population level estimates
GLMM_pop <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv")
names(GLMM_pop)
# select the random effects for the var of interest only (taxon_site)
# and tidy the names so we ahve taxona dn site separately 
# (the output from the GLMM gives us these as one variable)
GLMM_pop <- GLMM_pop %>% filter(str_detect(group_var, "taxon_site"))
# tidy up the names of the 
GLMM_pop <- GLMM_pop %>%
  separate(group, into = c("genus", "species", "site_name"),
           sep = "_", extra = "merge",  # merge all remaining pieces into 'site_name'
           fill  = "right", remove = FALSE) %>%
  mutate(Taxon  = str_c(genus, species, sep = "_")) %>%
  dplyr::select(Taxon, site_name) %>%
  filter(Taxon %in% demo_species)
GLMM_pop

GLMM_pop_species <- GLMM_pop %>% pull(Taxon) %>% unique() # 41 species
length(GLMM_pop_species)
GLMM_pop_sites <- GLMM_pop %>% pull(site_name) %>% unique() # 40 sites
length(GLMM_pop_sites)

# NOW get the coordinates for the 1) DRAGNet and 2) COMPADRE data sets

# 1) DRAGNet

# load the DRAGNet meta data to get the Lat / long of the sites 
meta_drag <- read_csv("data/site-info-drag-2024-07-29.csv")
meta_drag <- meta_drag %>%
  filter(site_name %in% GLMM_pop_sites) %>% 
  dplyr::select(site_name, site_code, latitude, longitude)
# check we now have 40 sites
n_distinct(meta_drag$site_name) # 40 sites

# we now need to merge this data with the GLMM data file
# that we have used in our second pop. level analysis
# remember this data is for all times and experiments
# we want unique site and taxa combinations
GLMM_pop_unique <- GLMM_pop %>% group_by(site_name, Taxon) %>% summarise() # 148 site Taxon combinations
both <- GLMM_pop_unique %>% left_join(meta_drag)
n_distinct(both$Taxon) # 41 taxa
n_distinct(both$site_name) # 40 site name
write_csv(both, "results/DRAGNet_species_locations_final.csv")

# 2) COMPADRE
# load the compadre meta data
meta_comp <- read_csv("results/meta_data_compadre.csv")
n_distinct(meta_comp$SpeciesAccepted) # 460 species 

meta_comp <- 
  meta_comp %>% 
  dplyr::select(MatrixID, SpeciesAccepted, Lat, Lon, Country, Continent) %>%
  rename(Taxon = SpeciesAccepted) %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) 
n_distinct(meta_comp$Taxon) # 460 species 

small <- meta_comp %>%
  filter(Taxon %in% GLMM_pop_species) %>%
  dplyr::select(Taxon, Lat, Lon) %>%
  count(Taxon, Lat, Lon, name = "number_matrices")

# unique number of species in data set
n_distinct(small$Taxon) # 41 species 
# unique number of locations in data set
n_distinct(small$Lat)

small %>% summarise(sum(number_matrices))

# the problem is that a lot the matrices don't have a location associated with them

# check that is is not a coding issue but a data base issue with the missing compadre 

# Download the whole compadre data base using the cdb_fetch command from Rcompadre
Com_dat <- cdb_fetch("compadre")

# compadre metrics
demo <- read_csv("results/all_COMPADRE_metrics.csv")
# filter for only the overlap species found in the initial text search
demo <- demo %>% filter(DRAGNet == TRUE)
n_distinct(demo$MatrixID) # 373 matrices

# get a vector of Matrix IDs that we used in our analysis
demo_mat <- demo %>% 
  pull(MatrixID) %>%
  unique()

# now we can filter the COMPADRE data using this vector
Com_subset <- Com_dat %>%
  filter(MatrixID %in% demo_mat)
names(Com_subset)
# extract the meta data from the compdare object
Com_meta <- Com_subset %>% 
  as.tibble() %>% 
  dplyr::select(SpeciesAccepted, Lat, Lon, Continent, MatrixID, Country, Ecoregion, DOI_ISBN)

see <- Com_meta %>% arrange(DOI_ISBN)

Com_meta <- Com_meta %>%
  rename("Taxon" = "SpeciesAccepted") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_"))
n_distinct(Com_meta$Taxon) # 42 species here

# if I drop NA from the lat / lon cols
Com_meta_clean <- Com_meta %>% drop_na(Lat, Lon)
n_distinct(Com_meta_clean$Taxon) # 34 species which we have data for

setdiff(union(Com_meta$Taxon, both$Taxon), intersect(Com_meta$Taxon, both$Taxon))
# the species that is in compadre but that we don't have a GLMM estimate for is "Cirsium_undulatum"
# get rid of that before we save the df

Com_meta <- Com_meta %>% filter(!Taxon == "Cirsium_undulatum")
n_distinct(Com_meta$Taxon) # now we have 41 taxa
# check that there is now no dif between the Taxon names
setdiff(union(Com_meta$Taxon, both$Taxon), intersect(Com_meta$Taxon, both$Taxon))

write_csv(Com_meta, "results/COMPADRE_species_locations_final.csv")

################################################################################
# 2) READ IN THE DATA - extract the climate variables for DRAGNet and COMPADRE
###############################################################################

# Read in location files
dragnet <- read_csv("results/DRAGNet_species_locations_final.csv")
names(dragnet)
dragnet <- dragnet %>% rename("Lat" = "latitude", "Lon" = "longitude")

compadre <- read_csv("results/COMPADRE_species_locations_final.csv")
names(compadre) # 
see_compadre <- compadre[complete.cases(compadre[, c("Lon", "Lat", "Taxon")]), ]
n_distinct(see_compadre$Taxon) # so there are 8 taxa that don't have any location data in COMPADRE

# Read in World Clim data (need the geodata library to do this)
message("Downloading WorldClim bioclimatic data (~10 arcmin)...")
bio <- geodata::worldclim_global(var = "bio", res = 10, path = tempdir())
#Ensure naming is tractable
names(bio) <- paste0("bio", 1:19)

#Clean coordinates and species NA rows
dragnet$Lon <- as.numeric(dragnet$Lon); dragnet$Lat <- as.numeric(dragnet$Lat)
dragnet <- dragnet[complete.cases(dragnet[, c("Lon", "Lat", "Taxon")]), ]

compadre$Lon <- as.numeric(compadre$Lon); compadre$Lat <- as.numeric(compadre$Lat)
compadre <- compadre[complete.cases(compadre[, c("Lon", "Lat", "Taxon")]), ]

#Crop to extent of points (with small buffer)
ext <- terra::ext(min(c(dragnet$Lon, compadre$Lon)) - 2,
                  max(c(dragnet$Lon, compadre$Lon)) + 2,
                  min(c(dragnet$Lat, compadre$Lat)) - 2,
                  max(c(dragnet$Lat, compadre$Lat)) + 2)
bio <- terra::crop(bio, ext)

#Function to extract climate variables for each of the rasters
extract_clim <- function(df, bio_rast){
  pts <- cbind(df$Lon, df$Lat)
  vals <- terra::extract(bio_rast, pts)[, -1, drop = FALSE]  # drop ID
  out <- cbind(df, vals)
  out <- out[complete.cases(out[, grep("^bio", names(out))]), , drop = FALSE]
  return(out)
}

dr_clim_all <- extract_clim(dragnet, bio)
n_distinct(dr_clim_all$Taxon) # we have 41 taxa here
cp_clim_all <- extract_clim(compadre, bio)
n_distinct(cp_clim_all$Taxon) # but only 33 taxa here
# columns of this data set are for the 16 bioclim variables

################################################################################
# 2) CALCULATE DISTANCE METRICS
###############################################################################

# FUNCTION to calculate the distances for 2 different matrices from DRAGNet and COMPADRE
pca_distance_one_species <- function(sp, dr_clim_all, cp_clim_all, clim_vars,
                                     var_threshold = 0.95) {
  
  # Subset to one species
  dr_sp <- dr_clim_all %>% filter(as.character(Taxon) == sp)
  cp_sp <- cp_clim_all %>% filter(as.character(Taxon) == sp)
  
  # Skip species with no data in one of the datasets
  if (nrow(dr_sp) == 0 || nrow(cp_sp) == 0) return(NULL)
  
  # Climate matrices (unique rows only)
  A <- dr_sp %>%
    dplyr::select(all_of(clim_vars)) %>%
    distinct() %>%
    as.matrix()
  
  B <- cp_sp %>%
    dplyr::select(all_of(clim_vars)) %>%
    distinct() %>%
    as.matrix()
  
  # Skip species where distinct() wipes everything
  if (nrow(A) == 0 || nrow(B) == 0) return(NULL)
  
  # Combine A and B to learn a common PCA basis
  AB <- rbind(A, B)
  
  # 🔹 Remove zero-variance (constant) columns per species
  sds  <- apply(AB, 2, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  
  # If after filtering there isn't enough variation, skip this species
  if (sum(keep) < 2) {
    message("Skipping species ", sp, " (not enough climatic variation after filtering).")
    return(NULL)
  }
  
  AB_red <- AB[, keep, drop = FALSE]
  
  # Reconstruct A and B from reduced matrix
  A_red <- AB_red[1:nrow(A), , drop = FALSE]
  B_red <- AB_red[(nrow(A) + 1):nrow(AB_red), , drop = FALSE]
  
  # PCA on reduced climate space
  pca <- prcomp(AB_red, center = TRUE, scale. = TRUE)
  
  eigvals <- pca$sdev^2
  cumvar  <- cumsum(eigvals) / sum(eigvals)
  K <- which(cumvar >= var_threshold)[1]
  if (is.na(K)) K <- length(eigvals)  # fallback (just in case)
  
  scores <- pca$x[, 1:K, drop = FALSE]
  
  A_pc <- scores[1:nrow(A_red), , drop = FALSE]
  B_pc <- scores[(nrow(A_red) + 1):nrow(scores), , drop = FALSE]
  
  # Pairwise distances A × B (n × m)
  D_pc <- proxy::dist(A_pc, B_pc, method = "euclidean")
  D_pc <- as.matrix(D_pc)
  
  # Long format: one row per DR–CP pair
  dist_long <- as.data.frame(as.table(D_pc))
  colnames(dist_long) <- c("dr_idx", "cp_idx", "distance")
  
  dist_long <- dist_long %>%
    mutate(
      dr_idx = as.integer(dr_idx),
      cp_idx = as.integer(cp_idx),
      Taxon  = sp
    ) %>%
    relocate(Taxon)
  
  return(dist_long)
}

# list the columns you want
clim_vars <- c("bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",
               "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
               "bio17", "bio18", "bio19")

# Species present in both datasets
species_vec <- intersect(
  unique(dr_clim_all$Taxon),
  unique(cp_clim_all$Taxon)
)

all_distances_long <- map_dfr(
  species_vec,
  ~ pca_distance_one_species(.x, dr_clim_all, cp_clim_all, clim_vars)
)

all_distances_long
names(all_distances_long)

see_pairs <- all_distances_long %>%
  group_by(Taxon) %>%
  summarise(n_unique_pairs = nrow(distinct(across(c(dr_idx, cp_idx)))))
max(see_pairs$n_unique_pairs) # 24
min(see_pairs$n_unique_pairs) # 1
mean(see_pairs$n_unique_pairs) # 4.5



################################################################################
# 3) VISUALISE DISTANCES per species
###############################################################################

plot <- all_distances_long %>%
  group_by(Taxon) %>%
  mutate(mean_dist = mean(distance, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = reorder(Taxon, mean_dist), y = distance)) +
  geom_boxplot() +
  labs(x = "Species (ordered by mean climatic distance)",
       y = "Distance",
       title = "DRAGNet and COMPADRE climatic distances per species") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
plot

ggsave("figures/SM_S1f_boxplot_climate_distances.pdf", plot, width = 8, height = 5)

###############################################################################
# 4) Examine mean and median distances
###############################################################################
# MEAN DISTANCES per Taxa
# visualise the mean distances
all_distances_long %>% group_by(Taxon) %>% 
  summarise(mean_distance = mean(distance)) %>%
  ggplot(aes(x= reorder(Taxon, mean_distance), mean_distance))+
  geom_col()

# save df as mean distances
mean_distances <- all_distances_long %>% group_by(Taxon) %>% 
  summarise(mean_distance = mean(distance))
n_distinct(mean_distances$Taxon) 
range(mean_distances$mean_distance)


# MEDIAN DISTANCES per Taxa
# visualise the median distances
all_distances_long %>% group_by(Taxon) %>% 
  summarise(median_distance = median(distance)) %>%
  ggplot(aes(x= reorder(Taxon, median_distance), median_distance))+
  geom_col()

# save df as median distances
median_distances <- all_distances_long %>% group_by(Taxon) %>% 
  summarise(median_distance = median(distance))
range(median_distances$median_distance)

# save the master distance file for plotting and use in applying the extra modelling
write_csv(all_distances_long, "results/all_distances.csv")

################################################################################
# 5) Load and tidy demographic and species cover change (random effects) data
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
# 6) Create and tidy a nested df for modelling with different levels of data filtering
################################################################################

all_distances_long <- read_csv("results/all_distances.csv")

med_dist <- all_distances_long %>%
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
    species_kept = purrr::map(
      n_dropped,
      ~ med_dist$Taxon[(.x + 1):n_sp]  # drop top .x high-distance species
    )
  )

sensitivity_results <- sens_grid %>%
  mutate(
    model_output = purrr::map(species_kept, ~ run_robust_models_for_subset(.x, both))
  ) %>%
  unnest(model_output)


sensitivity_results %>%
  filter(model == "Ordbeta") %>%
  arrange(Demo_trait, time_period, experiment, n_dropped) %>%
  head()

#################################################################################
# 7) Visualise results for the figures in the SM
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

# create a summary labelling data set
label_df <- sensitivity_results %>%
  filter(model == "Ordbeta") %>%
  filter(experiment != "INTER") %>%
  left_join(sens_grid, by = "n_dropped") %>%
  group_by(experiment, time_period, Demo_trait, hypoth) %>%  # one label per line
  slice_max(order_by = n_dropped, n = 1, with_ties = FALSE) %>%
  ungroup()


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
    nudge_x = 0.3,
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

ggsave("figures/SM_S1g_climate_senstivity_coeffcient_plot.pdf", plot, height = 5, width = 8)


