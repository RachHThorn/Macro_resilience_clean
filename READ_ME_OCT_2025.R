# R Thornley
# Macro-ecology DRAGNet and COMPADRE overlap paper
# 2024-2025

################################################################################
# Required Packages for project #
################################################################################

# the mixedup package is not available fro latest versions of R on CRAN so needs 
# installing directly from the github repro
install.packages("remotes")
remotes::install_github('m-clark/mixedup', force = TRUE)
install.packages("phytools")
install.packages("kableExtra")

# load packages
packages <- c("rstan", "assertthat", "remotes", "mixedup", "tidyverse", "Rcompadre",
              "popbio", "popdemo", 'Rage', "glmmTMB", "lme4", 'mixedup', 'DHARMa', 'broom', 
              "MASS", 'robustbase', 'vegan', 'rtry', 'ggtern', 'patchwork', 'gt', 
              "phytools", "ape", "kableExtra", "gt", "sf", "rnaturalearth", "rnaturalearthdata",
              "scales", "magick", "ggh4x")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

################################################################################

library(tidyverse) # general data wrangling
library(Rcompadre) # package to interact with COMPADRE database
library(popbio) # basic MPM analyses
library(popdemo) # advanced MPM transient measures
library(Rage) # life history metrics from MPMs

library(glmmTMB) # GLMMs package - runs the ordered beta-reg and hurdle models
library(lme4) # provides functions to interrogate GLMMs
library(mixedup) #  provides functions to interrogate GLMMs
library(DHARMa) # residual simulation and model checking for GLMMs
library(broom) # tidy models and extract coefficients
library(MASS) # check if I need this!!!!
# note you need dplyr::select in this code as the MASS package creates a conflict
# or run the below to ensure that stand alone select is the dplyr version
select <- dplyr::select
library(robustbase) # runs robust OLS models
library(mixedup)

library(vegan) # create community diversity metrics
library(rtry) # package that loads data downloaded from the TRY plant data base

library(phytools) # package for looking at the phylo dependency of the species
library(ape) # package for reading in phylo trees

library(ggtern) # plotting CSR strategies in ternary plot (this won't load now in my version of R!)
# this package has not been updated to work with the new ggplot and may not load
library(patchwork) # helps in creating plot layouts within ggplot
library(gt) # creates tables
library(kableExtra)
library(ggh4x)

library(rnaturalearth) # plotting maps 
library(rnaturalearthdata) # plotting maps
library(scales) # for displaying maps 
library(magick) # for trimming and displaying maps 

library(gt) # for modifying tables

library(sf) # plotting vector data for maps

################################################################################
# Data Processing and Analysis Scripts # (folder: S_scripts)
################################################################################

# script 1: Get a list of species that overlap in DRAGNet and COMPADRE etc.
# Make sure we are working with the latest core data
source("code/S_scripts/S1_get_summary_information_DRAGNet_COMPADRE.R")

# script 2: Get_master_data_DRAGNet_modelling
source("code/S_scripts/S2_get_master_data_DRAGNet_modelling.R")

# script 3: Run GLMMs and extract Random Effect estimates at species and species within site level
# WARNING - this code can hang a bit due to model calculations
source("code/S_scripts/S3_GLMM_extract_random_effects_all_DRAGNet_cover_data.R")

# script 4 : Get demographic variables from MPMs found in COMPADRE 
source("code/S_scripts/S4_get_demographic_metrics_COMPADRE.R")

# script 5 : Model DRAGNet Taxon effects from GLMMs with COMPADRE demographic variables
source("code/S_scripts/S5_bivariate_OLS_models_robust.R")

# script 6 : Get community complexity variables (richness, diversity, plot level cover)
source("code/S_scripts/S6_get_community_metrics_T0.R")

# script 7a : Get community context variables from functional databases
# this script can't be run entirely from source as the files from TRY are large and 
# are not stored in Github - I have commented out the section that can't be run
source("code/S_scripts/S7_species_strategies_calculate_T0.R")

# script 8 : Test for phylogenetic signal in demographic traits
source("code/S_scripts/S8_test_for_phylogenetic_signal.R")

# script 9 : Model DRAGNet Taxon effects from GLMMs with COMPADRE demographic variables and community complexity variables
source("code/S_scripts/S9_comm_OLS_models_robust_scaled.R")

################################################################################
# Scripts to make figures for main body of manuscript # (folder: F_scripts)
################################################################################

# Figure 1: Map figure with hypotheses
source("code/F_scripts/F1_map.R")

# Figure 2 is a schematic of the analysis work flow and is produced outside R

# Figure 3: Robust OLS 
source("code/F_scripts/F3_robust_OLS_effects.R")

# Figure 4: Visualise Community variable interaction results
source("code/F_scripts/F4_plot_robust_interactions.R")

# Table 1: Variable description table
source("code/T1_variables.R")

################################################################################
# Appendices / supporting information figures and tables # (folder: A_scripts)
################################################################################

# A1: A summary description of the DRAGnet protocol

# A2: get sample sizes with zeros (DRAGNet)
source("code/A_scripts/A2_get_sample_sizes_with_zeros.R")

# A3 - equations of the beta reg models

# A4: spread of the demographic metrics
source("code/A_scripts/A4_spread_compadre_data.R")

# A5_A6: data summary tables
source("code/A_scripts/A5_A6_data_summary_table.R")

# A7: PCA Compadre
source("code/A_scripts/PCA_compadre.R")

# A8: effect size contextualisation from GLMM models
source("code/A_scripts/A8_random_effects_visualise.R")

# A9: effect size contextualisation from GLMM models
source("code/A_scripts/A9_test_for_phylogenetic_signal.R")

# A10: GLMM diagnostics
source("code/A_scripts/A10_GLMM_model_diagnostics_signal.R")

# A11: Example residuals of OLS models
source("code/A_scripts/A11_residual_patterns_interactions.R")

# A12: Model results tables (bivariate OLS and interaction models)
source("code/A_scripts/A12_robust_model_results_table.R")

# A13: Predictions of interactions by hypothesis
source("code/A_scripts/A13_results_interactions_by_hypoth.R")

# 
