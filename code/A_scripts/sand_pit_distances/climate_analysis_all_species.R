# R Thornley 
# 09/12/2025
# Project: P1_COMPADRE_DRAGNET
# Appendices data overlap: distance metrics for climate data 

library(tidyverse)
library(terra)
library(sp)
library(geodata)
library(proxy) # for the distance metrics
library(purrr)

################################################################################
# 1) READ IN THE DATA - extract the climate variables for DRAGNet and COMPADRE
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

ggsave("figures/boxplot_climate_distances.pdf", plot, width = 8, height = 5)

###############################################################################
# 4) find mean distances and create cut offs
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
