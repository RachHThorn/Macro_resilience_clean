# R Thornley
# 28/10/2025
# Project: P1_COMPADRE_DRAGNET
# Script: Ax_climatic_overlap_species
# find in multi-variate climate space whether the species locations
# in dragnet overlap with those in compadre for our focal species

library(tidyverse)

################################################################################
# 1) IMPORT DATA ON SPECIES LOCATIONS and EXPORT FOR ROB
################################################################################

# DRAGNet site and taxa data load (all taxa)
# using the data for the second stage of modelling as input
Drag_dat <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv") %>%
  filter(group_var == "taxon_site") %>%
  filter(model == "Ordbeta") %>%
  separate(group, into = c("sp1", "sp2", "site_name"), sep = "_") %>%
  mutate(Taxon = paste0(sp1,"_", sp2)) %>%
  dplyr::select(Taxon, site_name) %>%
  distinct(Taxon, site_name)

# Load the list of taxa used in the modelling and filter the DRAGNet data for these taxa
mod_taxa <- readRDS("results/List_taxa_OLS_mods.R")
Drag_dat <- 
  Drag_dat %>% filter(Taxon %in% mod_taxa) 

unique(Drag_dat$site_name) # 40 sites
unique(Drag_dat$Taxon) # 39 taxa

# load the DRAGNet meta data to get the Lat / long of these taxa/site combinations
# we have used in our analysis
meta <- read_csv("data/site-info-drag-2024-07-29.csv") %>%
  filter(site_name %in% list_sites) %>% 
  dplyr::select(site_name, site_code, latitude, longitude)

# now join the meta data with the dragnet data 
cover_locations <- Drag %>% left_join(meta, by = "site_name")
# check we still have the same number of sites and taxa
n_distinct(cover_locations$Taxon) # 39 taxa
n_distinct(cover_locations$site_name) # 40 sites

# save this file for the spatial overlap analysis
write_csv(cover_locations, "results/DRAGNet_species_locations_final.csv")

# now load the compadre meta data to get the species locations where the MPMs were created
# filter for the species we have modelled and summarise to unique location, species combinations
# listing the number of matrices
matrix_locations <- read_csv("results/meta_data_compadre.csv") %>%
  dplyr::select(MatrixID, SpeciesAccepted, Lat, Lon, Country, Continent) %>%
  rename(Taxon = SpeciesAccepted) %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  filter(Taxon %in% mod_taxa) %>%
  dplyr::select(Taxon, Lat, Lon) %>%
  drop_na() %>%
  count(Taxon, Lat, Lon, name = "number_matrices")
# export this df as a file for the spatial overlap analysis
write_csv(matrix_locations, "results/COMPADRE_species_locations_final.csv")


# save the modelling data
dat <- read_csv("results/OLS_modelling_master_data.csv")
names(dat)
