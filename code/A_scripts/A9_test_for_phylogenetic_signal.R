# R Thornley
# 05/03/2025
# Test for phylo signal within our demographic vars using the compadre tree

library(tidyverse)
library(phytools)

rm(list= ls())

# Load phylogenetic tree (as an object "tree") 
com_tree <- read.tree("data/COMPADRE-COMADRE_Phylo_June_16_2019.tre")
plot(com_tree)
com_tree$tip.label

# Load demography variables
demo <- read_csv("results/March_2025/all_COMPADRE_metrics.csv")

# pull out the taxa that is in DRAGNet only
demo <- demo %>% filter(DRAGNet == TRUE)
# rename some variables
demo <- demo %>% rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var")
# reformat the Taxon variable
demo <- demo %>% mutate(Taxon = str_replace_all(Taxon, " ", "_")) 

################################################################################

# pivot the data longer and create a nested df for applying the function
dat_nested <- demo %>%
  group_by(Demo_trait) %>% 
  nest()

# test the function for one of the nested elements of the data frame
#data <- dat_nested %>% pull(data) %>% pluck(8)

# function that tests for phylo dependency for each trait and exports results
get_phylo_signal <-function(data) {
  
  species <- data %>% dplyr::select(Taxon, value) %>% drop_na() %>% distinct(Taxon) %>% pull(Taxon)
  pruned_tree <- drop.tip(com_tree, setdiff(com_tree$tip.label, species))
  plot(pruned_tree)
  
  data <- data %>% dplyr::select(Taxon, value) %>% drop_na() %>% as.data.frame()
  data <- setNames(data$value, data$Taxon)
  
  result_k <- phylosig(pruned_tree, data, method = "K", test = TRUE, nsim = 999)
  test_value <- result_k$K
  p_value <- result_k$P
  result_k <- as.data.frame(cbind(test_value, p_value))
  result_k$test_name <- "Blomberg_k"
  
  result_lambda <- phylosig(pruned_tree, data, method = "lambda", test = TRUE, nsim = 999)
  test_value <- result_lambda$lambda
  p_value <- result_lambda$P
  result_lambda <- as.data.frame(cbind(test_value, p_value))
  result_lambda$test_name <- "Pagel_lambda"
  
  result <- result_k %>% rbind(result_lambda)
  return(result)
  
} 

# now apply function to the nested data
results <- dat_nested %>% dplyr::mutate(phylo_depend = purrr::map(data, ~ get_phylo_signal(.x))) %>%
  dplyr::select(Demo_trait, phylo_depend) %>% unnest(cols = c(phylo_depend))
results

# get this in a nice table format
library(kableExtra)

table_results <- results %>% dplyr::select(Demo_trait, test_name, test_value, p_value) %>% 
  mutate(Demo_trait = case_when(Demo_trait == "T_generation" ~ "Generation time",
                                      Demo_trait == "R0" ~ "Net reproductive rate",
                                      Demo_trait == "percapita_repro" ~ "Per Capita Reproduction",
                                      Demo_trait == "Lmean" ~ "Mean Life Expectancy",
                                      Demo_trait == "age_repro" ~ "Age At First Reproduction",
                                      Demo_trait == "Lambda" ~ "Population Growth Rate",
                                      Demo_trait == "Reactivity" ~ "Reactivity", 
                                      Demo_trait == "FirstStepAtt" ~ "First step attenuation",
                                      Demo_trait == "Lmax" ~ "Maximum Life Expectancy",
                                      TRUE ~ Demo_trait)) %>%
  rename("Demographic variable" = Demo_trait, "Phylogenetic test" = test_name, 
       "Phylogenetic test value" = test_value, "p value" = p_value)

write_csv(table_results, "results/March_2025/phylo_dependency.csv")

###############################################################################

# Interpret results:
# - If the p-value is low (e.g., < 0.05), there is significant phylogenetic signal
# - The K value itself indicates the strength of the signal (values closer to 1 indicate stronger signal)               