# R Thornley 
# 29/10/2025
# Project: P1_COMPADRE_DRAGNET
# Script: Work out taxonomic overlap for climate plotting (Rob: spatial analysis)

library(tidyverse) # general data wrangling
library(Rcompadre) # package to interact with COMPADRE database
library(popbio) # basic MPM analyses
library(popdemo) # advanced MPM transient measures
library(Rage) # life history metrics from MPMs

###############################################################################
# LOAD compadre metrics and the results of the GLMMs (species responses at species and population levels)
################################################################################

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

################################################################################
# NOW get the coordinates for the 1) DRAGNet and 2) COMPADRE data sets
################################################################################

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

################################################################################
# check that is is not a coding issue but a data base issue with the missing compadre 
################################################################################

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

# just read these files in a check we are happy
comp <-read_csv("results/COMPADRE_species_locations_final.csv")

n_distinct(comp$Taxon) # 41
n_distinct(comp$Lat) # 55
n_distinct(comp$Lon) # 55


drag <- read_csv("results/DRAGNet_species_locations_final.csv")
names(drag)
n_distinct(drag$Taxon) # 41
n_distinct(drag$latitude) # 40
n_distinct(drag$longitude) # 40





