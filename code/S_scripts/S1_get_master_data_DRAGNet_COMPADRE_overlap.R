# Author: R Thornley 
# Date Final Version: 15/01/2026
# Github Repro: Macro_resilience_clean
# Script Name: S1_get_master_data_DRAGNet_COMPADRE_overlap.R

################################################################################
# Description
################################################################################

# This script processes the full DRAGNet data set to the subset of data presented in 
# the main body of the manuscript
# by filtering for the treatments, time frames analysed (years 0-3)
# and by filtering for the initial overlap in species between DRAGNet and COMPADRE
# note: the final analysis may not use all these data

################################################################################
# Instructions
################################################################################
# 1) Remove all non taxa entries from DRAGNet master data set / tidy years
# 2) Download COMPADRE to see what the initial overlap is
# 3) Use the overlap vector to filter the DRAGNet data to create the master data that gets passed onto S2


rm(list = ls())

################################################################################
# 1) Remove all non taxa entries from DRAGNet master data set / tidy years
################################################################################

# create vector of none taxa entries in cover data
none_taxa <- c("Fungi", "Other_litter", "Other_standing_water", "Ground", "Bryophyte", 
               "Other_animal_diggings", "Other_woody_overstory", "Lichen",
               "Other_animal_droppings", "Other_rock",
               "Other_soil_biocrust", "Other_standing_dead")

# read in master raw data file and tidy taxon entries and time points
drag <- 
  read.csv("data/full-cover-2025-09-11.csv") %>%
  mutate(New_taxon = str_to_sentence(Taxon)) %>%
  mutate(New_taxon = str_replace_all(New_taxon, " ", "_")) %>%
  filter(!str_detect(New_taxon, ".sp")) %>% # get rid of entries not to taxon level
  filter(!str_detect(New_taxon, "_x_")) %>% # get rid of any hybrids
  filter(!str_detect(New_taxon, "Unknown")) %>% # get rid of unknown species
  filter(!New_taxon %in% none_taxa) %>% # get rid of non taxa entries 
  mutate(New_taxon = case_when(New_taxon == "Helianthemum_nummularium_var._Grandiflorum" ~ "Helianthemum_nummularium",
                               New_taxon == "Mimosa_quadrivalvis_var._Platycarpa" ~ "Mimosa_quadrivalvis",
                               New_taxon == "Sebaea_sedoides_var._Schoenlandii" ~ "Sebaea_sedoides", 
                               TRUE ~ New_taxon)) %>%
  arrange(New_taxon) %>%
  mutate(year_trt = case_when(year_trt == -1 ~ "-T1",
                              year_trt == 0 ~ "T0",
                              year_trt == 1 ~ "T1",
                              year_trt == 2 ~ "T2",
                              year_trt == 3 ~ "T3",
                              year_trt == 4 ~ "T4",
                              year_trt == 5 ~ "T5",
                              year_trt == 6 ~ "T6"))

# filter for the first 3 years and the baseline year only
# and the experimental treatments that we are examining
# drop "NPK+Disturbance" and "NPK_Cessation"
# select the columns that are needed for the analysis

drag <- drag %>% 
  filter(year_trt %in% c("T0", "T1", "T1", "T2", "T3")) %>%
  filter(!trt %in% c("NPK+Disturbance", "NPK_Cessation")) %>% 
  droplevels() %>%
  dplyr::select(site_name, site_code, New_taxon, block, trt, year_trt, max_cover, year) 
unique(drag$trt)
unique(drag$year_trt)
drag_taxa <- drag %>% distinct(New_taxon) %>% pull()

# which years does the data run from and to
# min is 2018
# max is 2025
min(drag$year)
max(drag$year)
unique(drag$trt)
unique(drag$year_trt)


################################################################################
# 2) Download COMPADRE - and filter DRAGNet for initial taxonomic overlap
################################################################################

# Download the compadre data base using the cdb_fetch command from the Rcompadre package
Com_dat <- cdb_fetch("compadre")

# the SpeciesAccepted variable is the official taxon name
# get the list of taxa in COMPADRE and tidy the names
Com_taxa <- Com_dat$SpeciesAccepted
Com_taxa <- unique(Com_taxa)
Com_taxa <- str_replace_all(Com_taxa, " ", "_")

# Find shared taxa between the two data sets
common <- Com_taxa[Com_taxa %in% drag_taxa] 
common
length(common) # there are 76 shared taxa

# output the list of shared taxa as it is needed to filter the modelling data frames
write.csv(common, "results/common_species_drag_comp.csv")

################################################################################
# 3) Use the overlap vector to filter the DRAGNet data to create the master
# data that gets passed onto S2
################################################################################

new_drag <- drag %>% filter(New_taxon %in% common)
n_distinct(new_drag$site_name) # we have 48 sites
n_distinct(new_drag$New_taxon) # we have 76 taxon

write_csv(new_drag, "results/tidy_drag_from_S1.csv")

# NOTE: 
# The final analysis uses a smaller data set from 40 sites and from 39 species
# It is important that only this smaller data set is shared publicly
# I'm required by the network to give data authorship to the site co-ordinators 
# whose data is included in the main analysis only
# therefore I am going to filter this further - for only the sites that are included in the analysis
# although in reality some of the code will break using the DOI data as we don't get down to this
# level of data resolution until much later in the processing pipeline

