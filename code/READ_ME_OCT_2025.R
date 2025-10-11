# R Thornley
# Macro-ecology DRAGNet and COMPADRE overlap paper
# 2024-2025

################################################################################
# Required Packages for project #
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

library(vegan) # create community diversity metrics
library(rtry) # package that loads data downloaded from the TRY plant data base

library(ggtern) # plotting CSR strategies in ternary plot
library(patchwork) # helps in creating plot layouts within ggplot
library(gt) # creates tables

################################################################################
# Data Processing and Analysis Scripts # (folder: S_scripts)
################################################################################

# script 1: Get a list of species that overlap in DRAGNet and COMPADRE etc.
# Make sure we are working with the latest core data
source("code/S1_summary_information_DRAGNet_COMPADRE.R")

# script 2: Get_master_data_DRAGNet_modelling
source("code/S2_get_master_data_DRAGNet_modelling.R")

# script 3: Run GLMMs and extract Random Effect estimates at species and species within site level
source("code/S3_GLMM_extract_random_effects_all_DRAGNet_cover_data.R")

# script 4 : Get demographic variables from MPMs found in COMPADRE 
source("code/S4_get_demographic_metrics_COMPADRE.R")

# script 5 : Model DRAGNet Taxon effects from GLMMs with COMPADRE demographic variables
source("code/S5_bivariate_OLS_models_robust_final.R")

# script 6a : Get community complexity variables (richness, diversity, plot level cover)
source("code/S6_get_community_metrics_each_time.R")

# script 6b : Get community complexity variables (richness, diversity, plot level cover)
source("code/S6_get_community_metrics_T0.R")

# script 7a : Get community context variables from functional databases
source("code/S_scripts/S7_species_strategies_calculate_each_time.R")

# script 7b : Get community context variables from functional databases
source("code/S_scripts/S7_species_strategies_calculate_T0.R")

# script 8 : Test for phylogenetic signal in demographic traits
source("code/S8_test_for_phylogenetic_signal.R")

# script 9 : Model DRAGNet Taxon effects from GLMMs with COMPADRE demographic variables and community complexity variables
source("code/S9_comm_OLS_models_robust_scaled.R")

################################################################################
# Scripts to make figures for main body of manuscript # (folder: F_scripts)
################################################################################

# Figure 1: Map figure
source("code/F1_map.R")

# Figure 3: Robust OLS 
source("code/F3_robust_OLS_effects.R")

# Figure 4: Visualise Community variable interaction results
source("code/F4_plot_robust_interactions.R")

# Table 1: Variable description table
source("code/T1_variables.R")

################################################################################
# Appendices / supporting information figures and tables # (folder: A_scripts)
################################################################################

# script A1: get sample sizes with zeros (DRAGNet)
source("A1_get_sample_sizes_with_zeros.R")

# script A2: effect size contextualisation from GLMM models
source("A2_random_effects_visualise.R")

# script A4: spread of the demographic metrics
source("A4_spread_compadre_data.R")

# script A5: phylogenetic signal
source("A5_test_for_phylogenetic_signal.R")

# script Ax: results tables for the robust models
source("Ax_Robust_models_results_table.R")

# script Ax: PCA compadre
source("Ax_PCA_compadre.R")



