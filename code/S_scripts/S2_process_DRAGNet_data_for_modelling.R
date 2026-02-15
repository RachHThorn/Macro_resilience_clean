# Author: R Thornley 
# Date Final Version: 15/01/2026
# Github Repro: Macro_resilience_clean
# Script Name: S2_process_DRAGNet_data_for_modelling.R

################################################################################
# Description
################################################################################

# This script prepares the DRAGNet data for the GLMM modelling 

################################################################################
# Instructions
################################################################################
# 1) Find the unique combinations of blocks we require
# 2) Create paired values by adding zero values only where required
# 3) Further data wrangling: convert percentages to proportions
# 4) Output modelling data files that correspond to the 3 time periods we are analysing

rm(list = ls())

################################################################################
# 1) Find the unique combinations of blocks we require
################################################################################

# create paired data throughout the dataset by working out the unique combinations of quadrats
# alongside the years when sites have submitted data
# there are some years when the sites have not submitted data yet and these are not real zeros
# interpretation of col names in the table created
# blocks represented in -1 - means they have participated in the seed study and is a pre-treatment year
# 0 - year pre-treatment
# 1, 2, 3, etc years post treatment

# read in the full DRAGNet dataset and find the unique blocks listed
blocks_surveyed <- 
  read.csv("data/full-cover-2025-09-11.csv") %>% 
  group_by(site_name, year_trt) %>% 
  mutate(block = as.character(block)) %>% 
  summarise(n_block = n_distinct(block)) %>%
  pivot_wider(names_from = year_trt, values_from = n_block) %>% 
  replace(is.na(.),0) %>%
  dplyr::select(site_name, "-1", "0", "1", "2", "3", "4", "5")

colnames(blocks_surveyed)

################################################################################
# 2) FUNCTION: create paired values and other data wrangling
################################################################################

# read in the dragnet data from S1
drag <- read_csv("results/tidy_drag_from_S1.csv")
drag$year_trt

# Function: 
# first we filter for only the sites we want where the zeros are real
# not for unsurveyed sites / sites with missing data in the master data set
# we do this for the years specified for each master data set
# 

# Function that processes each dataset in turn
process_dragnet <- function(drag, blocks_surveyed, years) {
  # use the blocks_surveyed object to filter the dataframe for the unique combinations
  wanted <- 
    blocks_surveyed %>%
    setNames(c("site_name", "-T1", "T0", "T1", "T2", "T3", "T4", "T5")) %>%
    dplyr::select(site_name, T0, !!sym(paste0(years[length(years)]))) %>%
    pull(site_name)
  # filter the master dragnet data for the list of sites and years we want
  # get complete cases for the unique taxon/plot combinations
  # replace the NAs created with zeros
  # tidy factor levels
  # convert cover values to proportions / converting 100% to 0.99
  # filter for any plots that have no all zeros 
  dat <- drag %>%
    dplyr::filter(site_name %in% wanted) %>%
    dplyr::filter(year_trt %in% years) %>%
    dplyr::select(site_name, New_taxon, block, trt, year_trt, max_cover) %>%
    group_by(site_name) %>%
    complete(New_taxon, year_trt, trt, block) %>%
    replace(is.na(.), 0) %>%
    ungroup() %>%
    mutate(
      trt = factor(trt, levels = c("Control","Disturbance","NPK")),
      year_trt = factor(year_trt, levels = years),
      site_name = factor(site_name),
      New_taxon = factor(New_taxon),
      taxon_site = factor(paste0(New_taxon, "_", site_name)),
      new_max_cover = case_when(max_cover == 0 ~ 0,
                                max_cover > 0 ~ pmin(max_cover/100, 0.99),
                                TRUE ~ max_cover)
    ) %>%
    group_by(site_name, New_taxon, taxon_site, trt, block) %>%
    filter(any(new_max_cover != 0)) %>%
    ungroup() %>%
  
  return(dat)
}

################################################################################
# 3) Run FUNCTION for each time step and export files
# no overlap here between compadre and dragnet
################################################################################

# T0 and T1
dat_1 <- process_dragnet(drag, blocks_surveyed,
                         years = c("T0","T1"))
write_csv(dat_1, "results/DRAGNet_T0_T1_all.csv")

# T0 and T2       
dat_2 <- process_dragnet(drag, blocks_surveyed,
                         years = c("T0","T2"))
write_csv(dat_2, "results/DRAGNet_T0_T2_all.csv")

# T0 and T3
dat_3 <- process_dragnet(drag, blocks_surveyed,
                         years = c("T0","T3"))
write_csv(dat_3, "results/DRAGNet_T0_T3_all.csv")

################################################################################
# # 4) if you want to look at the initial overlap between compadre and DRAGNet
################################################################################

# # Filter for COMPADRE overlap and save
# dat_1_overlap <- dat_1 %>% filter(New_taxon %in% species)
# write_csv(dat_1_overlap, "results/DRAGNet_T0_T1_overlap.csv")
# 
# # Filter for COMPADRE overlap and save
# dat_2_overlap <- dat_2 %>% filter(New_taxon %in% species)
# write_csv(dat_2_overlap, "results/DRAGNet_T0_T2_overlap.csv")
# 
# # Filter for COMPADRE overlap and save
# dat_3_overlap <- dat_3 %>% filter(New_taxon %in% species)
# write_csv(dat_3_overlap, "results/DRAGNet_T0_T3_overlap.csv")
# 
# #################################################################################
