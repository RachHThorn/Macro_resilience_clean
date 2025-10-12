# R Thornley 
# 05/11/2024
# Project: P1_COMPADRE_DRAGNET
# Script: S1_get_summary_information_DRAGNET_COMPADRE

################################################################################
# Instructions
################################################################################

# 1) Get list of DRAGNet taxa
# 2) Get overlap species list DRAGNet and COMPADRE
# 3) Get a filtered COMPADRE data base with the 82 shared taxa locations for map plotting
# 4) Get list of blocks in DRAGNet surveyed per site and per time point

################################################################################
# 1) Get list of DRAGNet taxa
################################################################################

# Read in latest DRAGNet data and tidy up vector of unique taxa to match format of compadre

# create vector of none taxa entries in cover data to use in filtering
none_taxa <- c("Fungi", "Other_litter", "Other_standing_water", "Ground", "Bryophyte", 
               "Other_animal_diggings", "Other_woody_overstory", "Lichen",
               "Other_animal_droppings", "Other_rock",
               "Other_soil_biocrust", "Other_standing_dead")

# create a tidy character vector of unique taxa
Drag_taxa <- read_csv("data/full-cover-2025-09-11.csv") %>%
  dplyr::select(Taxon) %>%
  distinct(Taxon) %>%
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
  dplyr::select(New_taxon) %>%
  pull(New_taxon)

length(Drag_taxa) # there are 1301 taxa in the DRAGNet data set
Drag_taxa # all these names now look ok = all taxa and no subspecies present
write.csv(Drag_taxa, "results/all_drag_taxa.csv") # export a list of taxa found in DRAGNet 

################################################################################
# 2) Get overlap species list DRAGNet and COMPADRE
################################################################################

# Download the compadre data base using the cdb_fetch command from the Rcompadre package
Com_dat <- cdb_fetch("compadre")

# the SpeciesAccepted variable is the official taxon name
# get the list of taxa in COMPADRE and tidy the names
Com_taxa <- Com_dat$SpeciesAccepted
Com_taxa <- unique(Com_taxa)
Com_taxa <- str_replace_all(Com_taxa, " ", "_")

# Find shared taxa between the two data sets
common <- Com_taxa[Com_taxa %in% Drag_taxa] 
length(common) # there are 84 shared taxa
write.csv(common, "results/common_species_drag_comp.csv")

################################################################################
# 3) Get a filtered compadre data base with the 82 shared taxa locations for plotting
################################################################################

Com_dat <- cdb_fetch("compadre")
names(Com_dat)
Com_dat <- Com_dat %>% mutate(SpeciesAccepted = str_replace_all(SpeciesAccepted, " ", "_"))
Com_dat$SpeciesAccepted
Com_dat <- Com_dat %>% filter(SpeciesAccepted %in% common)
names(Com_dat)
Com_dat <- Com_dat %>% as.tibble() %>% dplyr::select(MatrixID, SpeciesAccepted, Lat, Lon, Country, Continent, DOI_ISBN)
Com_dat <- Com_dat %>% drop_na(Lat) %>% distinct(Lat, Lon, SpeciesAccepted, .keep_all=TRUE)

write_csv(Com_dat, "results/matrix_locations_dragnet.csv")

################################################################################
# 4) Get list of blocks in DRAGNet surveyed per site and per time point
###############################################################################

# dragnet master data load
drag <- read.csv("data/full-cover-2025-09-11.csv")

blocks_surveyed <- drag %>% 
  group_by(site_name, year_trt) %>% 
  mutate(block = as.character(block)) %>% 
  summarise(n_block = n_distinct(block)) %>%
  pivot_wider(names_from = year_trt, values_from = n_block) %>% 
  replace(is.na(.),0) %>%
  dplyr::select(site_name, "-1", "0", "1", "2", "3", "4", "5") %>%
  print(n = 56)

colnames(blocks_surveyed)
# output table the number of blocks surveyed in each time point
write.csv(blocks_surveyed, "results/dragnet_blocks_surveyed.csv") 

# interpretation of col names in the table created
# blocks represented in -1 - means they have participated in the seed study and is a pre-treatment year
# 0 - year pre-treatment
# 1, 2, 3, etc years post treatment


