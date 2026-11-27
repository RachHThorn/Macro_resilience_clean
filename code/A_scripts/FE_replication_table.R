# Author: R Thornley 
# Date Final Version: 26/11/2026
# Github Repro: Macro_resilience_clean
# Script Name: FE_replication_table.R

################################################################################
# Description
################################################################################

# script to calculate the replication values required for submission of the 
# manuscript to Functional Ecology

################################################################################
# Instructions
################################################################################

# 1) Row 1 of table: disturbance/control pairs for T0-T1; + site and species total count for this analysis
# 2) Row 2 of table: disturbance/control pairs for T0-T2; + site and species total count for this analysis
# 3) Row 3 of table: disturbance/control pairs for T0-T3; + site and species total count for this analysis

# 4) Row 4 of table: nutrient/control pairs for T0-T1; + site and species total count for this analysis
# 5) Row 5 of table: nutrient/control pairs for T0-T2; + site and species total count for this analysis
# 6) Row 6 of table: nutrient/control pairs for T0-T3; + site and species total count for this analysis

# Row 1-6 - species and time-point specific replicates for SM table

# 7) Row 7 of table: number of control plots used for community complexity metrics; + site counts
# 8) Row 8 of table: number of native and non-native records, counting each species once at each site

# 9) Row 9 of table: number of distinct matrices across COMPADRE for the 39 species; 

# Row 9 - species and demographic metric-specific replicates for SM table

################################################################################
# R1-R6
################################################################################

# get list of species used in the models
taxa <- read_csv("results/list_taxa_OLS_models.csv") %>%
  mutate(taxa = stringr::str_replace_all(taxa, " ", "_")) %>%
  pull(taxa)

R1 <- read_csv("results/DRAGNet_T0_T1_all.csv") %>% 
  filter(New_taxon %in% taxa) %>%
  group_by(trt, year_trt) %>%
  summarise(n_site = n_distinct(site_name),
            n_taxon = n_distinct(New_taxon))


R1 %>%
  filter(year_trt == "T1") %>%
  group_by(trt, year_trt) %>%
  summarise(n_site = n_distinct(site_name),
            n_taxon = n_distinct(New_taxon))

R2 <- read_csv("results/DRAGNet_T0_T2_all.csv") %>% 
  filter(New_taxon %in% taxa) %>%
  filter(year_trt == "T2") %>%
  group_by(trt, year_trt) %>%
  summarise(n_site = n_distinct(site_name),
            n_taxon = n_distinct(New_taxon))
            
R3 <- read_csv("results/DRAGNet_T0_T3_all.csv") %>% 
  filter(New_taxon %in% taxa) %>%
  filter(year_trt == "T3") %>%
  group_by(trt, year_trt) %>%
  summarise(n_site = n_distinct(site_name),
            n_taxon = n_distinct(New_taxon))      


################################################################################
# SM table of species and time-point specific counts
################################################################################

R1 <- read_csv("results/DRAGNet_T0_T1_all.csv") %>% 
  filter(trt == "Disturbance" & year_trt == "T1") %>% 
  group_by(New_taxon) %>% tally() %>% rename("T0_T1" = n)

R2 <- read_csv("results/DRAGNet_T0_T2_all.csv") %>% 
  filter(trt == "Disturbance" & year_trt == "T2") %>% 
  group_by(New_taxon) %>% tally() %>% rename("T0_T2" = n)

R3 <- read_csv("results/DRAGNet_T0_T3_all.csv") %>% 
  filter(trt == "Disturbance" & year_trt == "T3") %>% 
  group_by(New_taxon) %>% tally() %>% rename("T0_T3" = n)

DIST_time <- R1 %>% left_join(R2) %>% left_join(R3) %>% mutate(across(everything(), ~ tidyr::replace_na(.x, 0)))
DIST_time$Treatment <- "Disturbance"


R1 <- read_csv("results/DRAGNet_T0_T1_all.csv") %>% 
  filter(trt == "NPK" & year_trt == "T1") %>% 
  group_by(New_taxon) %>% tally() %>% rename("T0_T1" = n)

R2 <- read_csv("results/DRAGNet_T0_T2_all.csv") %>% 
  filter(trt == "NPK" & year_trt == "T2") %>% 
  group_by(New_taxon) %>% tally() %>% rename("T0_T2" = n)

R3 <- read_csv("results/DRAGNet_T0_T3_all.csv") %>% 
  filter(trt == "NPK" & year_trt == "T3") %>% 
  group_by(New_taxon) %>% tally() %>% rename("T0_T3" = n)

NPK_time <- R1 %>% left_join(R2) %>% left_join(R3) %>% mutate(across(everything(), ~ tidyr::replace_na(.x, 0)))
NPK_time$Treatment <- "Nutrient"

both <- DIST_time %>% rbind(NPK_time)

