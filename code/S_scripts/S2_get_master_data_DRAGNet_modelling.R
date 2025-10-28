# R Thornley 
# 27/08/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S2_get_master_data_DRAGNET_modelling

################################################################################
# Instructions
################################################################################

# 1) Remove all non taxa entries from DRAGNet master data set / tidy years
# 2) Create function that filters for the list of sites that we want and tidy names / cover values
# 3) Apply function to all the time points and experiments and export the two data sets (all drag and filtered drag)

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
         

################################################################################
# 2) Filter for the list of sites that we want and tidy names / cover values
###############################################################################

# file path to the df of sites and sampling years produced in S1
wanted_file_path <- "results/dragnet_blocks_surveyed.csv"

# function that processes the dagnet data filters it for sites needed for each 
# time frame we are interested in and preps the data for modelling

process_dragnet <- function(drag, wanted_file_path, years) {
  # Load site list
  wanted <- read_csv(wanted_file_path) %>%
    dplyr::select(2, 4:9) %>%
    setNames(c("site_name", "T0", "T1", "T2", "T3", "T4", "T5")) %>%
    dplyr::select(site_name, T0, !!sym(paste0(years[length(years)]))) %>%
    pull(site_name)

  # Filter data
  dat <- drag %>%
    dplyr::filter(site_name %in% wanted) %>%
    dplyr::filter(year_trt %in% years) %>%
    dplyr::select(site_name, New_taxon, block, trt, year_trt, max_cover) %>%
    group_by(site_name) %>%
    complete(New_taxon, year_trt, trt, block) %>%
    replace(is.na(.), 0) %>%
    ungroup() %>%
    mutate(
      trt = factor(trt, levels = c("Control","Disturbance","NPK+Disturbance","NPK","NPK_Cessation")),
      year_trt = factor(year_trt, levels = years),
      site_name = factor(site_name),
      New_taxon = factor(New_taxon),
      taxon_site = factor(paste0(New_taxon, "_", site_name)),
      new_max_cover = case_when(max_cover == 0 ~ 0,
                                max_cover > 0 ~ pmin(max_cover/100, 0.99),
                                TRUE ~ max_cover)
    ) %>%
    filter(trt != "NPK_Cessation") %>%
    group_by(site_name, New_taxon, taxon_site, trt, block) %>%
    filter(any(new_max_cover != 0)) %>%
    ungroup()
  
  return(dat)
}

# load in the species overlap file produced in S1
species <- 
  read_csv("results/common_species_drag_comp.csv") %>%
  pull(x)
length(species) # 85 species that overlap


################################################################################
# 3) Run for each time step and export
################################################################################

# T0 and T1
dat_1 <- process_dragnet(drag, wanted_file_path,
                         years = c("T0","T1"))
write_csv(dat_1, "results/DRAGNet_T0_T1_all.csv")

# Filter for COMPADRE overlap and save
dat_1_overlap <- dat_1 %>% filter(New_taxon %in% species)
write_csv(dat_1_overlap, "results/DRAGNet_T0_T1_overlap.csv")

# T0 and T2       
dat_2 <- process_dragnet(drag, wanted_file_path,
                         years = c("T0","T2"))
write_csv(dat_2, "results/DRAGNet_T0_T2_all.csv")
          
# Filter for COMPADRE overlap and save
dat_2_overlap <- dat_2 %>% filter(New_taxon %in% species)
write_csv(dat_2_overlap, "results/DRAGNet_T0_T2_overlap.csv")

# T0 and T3
dat_3 <- process_dragnet(drag, wanted_file_path,
                         years = c("T0","T3"))
write_csv(dat_3, "results/DRAGNet_T0_T3_all.csv")

# Filter for COMPADRE overlap and save
dat_3_overlap <- dat_3 %>% filter(New_taxon %in% species)
write_csv(dat_3_overlap, "results/DRAGNet_T0_T3_overlap.csv")

################################################################################