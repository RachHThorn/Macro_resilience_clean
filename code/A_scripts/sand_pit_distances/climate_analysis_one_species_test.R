# R Thornley 
# 09/12/2025
# Project: P1_COMPADRE_DRAGNET
# Appendices data overlap: distance metrics for climate data 

library(tidyverse)
library(terra)
library(sp)
library(geodata)
library(proxy) # for the distance metrics

################################################################################
# 1) READ IN THE DATA - extract the climate variables for DRAGNet and COMPADRE
###############################################################################

# Read in location files
dragnet <- read_csv("results/DRAGNet_species_locations_final.csv")
names(dragnet)
dragnet <- dragnet %>% rename("Lat" = "latitude", "Lon" = "longitude")

compadre <- read_csv("results/COMPADRE_species_locations_final.csv")
names(compadre)

out_dir <- "results/PCA_clim_plots_DRAGNet_COMPADRE"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
cp_clim_all <- extract_clim(compadre, bio)
# columns of this data set are for the 16 bioclim variables

################################################################################
# 2) prep the data for distance analysis: ONE SPECIES TEST
#################################################################################

# we would like to do a PCA on the species 
dr_clim_all %>% group_by(Taxon) %>% tally() %>% filter(n>3) # only 14 species have enough clim data to create the convex hull....
# therefore it would better to do another method - distance metrics could be an option
unique(dr_clim_all$Taxon)

# try this for one species only
sp <- "Centaurea_jacea"

# DRAGNet and COMPADRE subsets
dr_sp <- dr_clim_all %>%
  filter(as.character(Taxon) == sp)

cp_sp <- cp_clim_all %>%
  filter(as.character(Taxon) == sp)

## Optional quick check of the number of rows for each species
nrow(dr_sp)
nrow(cp_sp)

# list the columns you want
clim_vars <- c("bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",
               "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
               "bio17", "bio18", "bio19")

clim_combined <- bind_rows(
  dr_sp %>% mutate(dataset = "DR"),
  cp_sp %>% mutate(dataset = "CP")
)

clim_combined <- clim_combined %>% 
  dplyr::select(dataset, tidyselect::all_of(clim_vars))

################################################################################
# 3) create multi-variate distances: ONE SPECIES TEST
################################################################################

# create two separate matrices for the different data sets with unique rows only
# scale before transforming
A <- dr_sp %>% dplyr::select(all_of(clim_vars)) %>% distinct() %>% scale() %>% as.matrix()  # dataset 1 (n × p)
B <- cp_sp %>% dplyr::select(all_of(clim_vars)) %>% distinct() %>% 
  scale(center = attr(A, "scaled:center"), scale  = attr(A, "scaled:scale")) %>% as.matrix()
# use same scaling as A dataset 2 (m × p)

# do a PCA first - then find distances

# What this is doing conceptually
# Join A and B, to learn the common PCA rotation and scaling.
# Split them back into A_pc and B_pc (now K-dimensional, with K < n),
# Compute pairwise distances between matrices A and B, exactly as you wanted.
# Because K is much smaller than p and ≤ (n_total − 1), you can avoid the singularity problem 
# when we try to use that affects Mahalanobis in the original p-dimensional space.


# Combine A and B just to fit a *common* PCA
AB <- rbind(A, B)
pca <- prcomp(AB, center = TRUE, scale. = TRUE)

# Choose number of PCs (e.g. enough to explain 95% variance)
eigvals <- pca$sdev^2
cumvar  <- cumsum(eigvals) / sum(eigvals)
K <- which(cumvar >= 0.95)[1]   # first K PCs reaching 95%

# Project A and B into PC space (same rotation)
scores <- pca$x[, 1:K, drop = FALSE]

A_pc <- scores[1:nrow(A), , drop = FALSE]
B_pc <- scores[(nrow(A) + 1):nrow(scores), , drop = FALSE]

## Find distances between A and B in PC space (n × m)
D_pc <- proxy::dist(A_pc, B_pc, method = "euclidean")
D_pc <- as.matrix(D_pc)

dim(D_pc)  # n rows (A) × m cols (B)

################################################################################
# 3) Apply to all species
#################################################################################



