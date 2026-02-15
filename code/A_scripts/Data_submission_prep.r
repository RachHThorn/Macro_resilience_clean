# Author: R Thornley 
# Date Final Version: 15/01/2026
# Github Repro: Macro_resilience_clean
# Script Name: Data_submission_prep.r

################################################################################
# Description
################################################################################

# script to create the dataset submitted with the paper and the data author list
# it is important to have this data filtered to the site level as this is what dictates
# data ownership in the DRAGNet master data set
# we cannot publish the whole DRAGNet data set with the paper as this is against network regulations

################################################################################
# Instructions
################################################################################

# 1) Filter the DRAGNet data set for only the sites used in the final analysis
# 2) Tidy and Export the final data set and data authorship list

################################################################################
# 1) Filter the data exported from S1 further for only the sites used in the final analysis
################################################################################

# read in dragnet data from script 1
new_drag <- read_csv("results/tidy_drag_from_S1.csv")

# read in the final taxa lists from the final analyses (outputs of scripts 5 and 8)
OLS_taxa <- readRDS("results/List_taxa_OLS_models.R")
see_1 <- new_drag %>% filter(New_taxon %in% OLS_taxa)
species_1 <- n_distinct(see_1$New_taxon) # 39 taxa
sites_1 <- n_distinct(see_1$site_name) # 40 sites

# what happens if I filter the master data using the taxa list from the second stage of modelling
site_taxa <- readRDS("results/List_taxa_OLS_sites_mods.R")
see_2 <- new_drag %>% filter(New_taxon %in% site_taxa)
species_2 <- n_distinct(see_2$New_taxon) # 39 taxa
sites_2 <- n_distinct(see_2$site_name) # 40 sites

# check if the results from these two data sets are the same
identical(sites_1, sites_2) # ok the two site lists are the same
identical(species_1, species_2) # ok the two species lists are the same

# get output of these taxa as a csv for use in the DOI
taxa_list <- OLS_taxa %>% as.data.frame()
taxa_list <- taxa_list %>%
  rename(taxon_scientific_name = ".") %>%
  mutate(taxon_scientific_name = str_replace(taxon_scientific_name, "_", " "))
rownames(taxa_list) <- NULL
taxa_list$general_taxonomic_coverage <- "This is a plant"
taxa_list$taxon_rank <- "species"
write_csv(taxa_list, "results/List_taxa_data_storage.csv")

################################################################################
# 2) Tidy and Export the final data set and data authorship for saving
################################################################################

# for DRAGNet submission I need latitude and longitude for the sites too
# get this from the site info file provided by the network
site_info <- read_csv("data/site-info-2025-12-08.csv")

# join this to the final dragnet df I am using / tidy names
final_dat <- site_info %>% 
  dplyr::select(site_name, site_code, latitude, longitude) %>%
  right_join(see_1) %>%
  rename(taxon = New_taxon, treatment = trt)
# change taxon to remove the underscore
final_dat <- final_dat %>% mutate(taxon = str_replace(taxon, "_", " "))

# export file
write_csv(final_dat, "results/final_DOI_data_29_01_2026.csv")

# I also need a geographical coverage input for the EDI submission
rows_required <- c("geographicDescription","northBoundingCoordinate", "southBoundingCoordinate",
              "eastBoundingCoordinate", "westBoundingCoordinate", "minimumAltitude",
              "maximumAltitude,altitudeUnits")
small <- final_dat %>% dplyr::select(site_code, latitude, longitude)
bbox <- with(small, c(
  xmin = min(longitude, na.rm = TRUE),
  ymin = min(latitude, na.rm = TRUE),
  xmax = max(longitude, na.rm = TRUE),
  ymax = max(latitude, na.rm = TRUE)
))

bbox
# xmin = west bounding box
# ymin = south bounding box
# xmax = east bonding box
# ymax = north bounding box

# Get the list of data authors of the data used in the final analysis

# get a final list of sites from the DOI data set
final_sites_list <- unique(final_dat$site_code)

# filter the DRAGNet coordinator list for the data authors
data_authors <- read_csv("data/site-and-network-coordinators-2025-08-21.csv")
network_coord <- data_authors %>% filter(network_coordinator == "yes")
site_coord <- data_authors %>% filter(site_code %in% final_sites_list)
n_distinct(site_coord$site_code) # Millville is missing off this file

identical(final_sites_list, site_coord$site_code)
setdiff(final_sites_list, site_coord$site_code)

# mill_dn.us
see_1 %>% filter(site_code == "mill_dn.us")
site_coord %>% filter(site_code == "mill_dn.us")
data_authors %>% filter(site_code == "mill_dn.us")
site_info %>% filter(site_code == "mill_dn.us")

# the site and network coordinators file is missing information for the mill_dn.us site
# I need to add this manually

# authors list combine
final_authors <- network_coord %>% rbind(site_coord) %>%
  distinct() %>% arrange(lastname)
unique(final_authors$site_code) 
unique(final_authors$site_name) 
names(final_authors)
# Add missing row
final_authors <- final_authors %>%
  add_row(lastname = "Kulmatiski", firstname = "Andrew", 
          site_name = "Millville", site_code = "mill_dn.us") %>%
  arrange(lastname)
# there are 76 data authors

# save this df
write_csv(final_authors, "results/data_authorship_list.csv")



