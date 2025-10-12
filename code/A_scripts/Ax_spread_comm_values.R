# R Thornley
# 06/10/2025
# Spread community values - for establishing the ranges of prediction values 
# to work with in interpretation of the interaction effects in the community models

# this is not yet included in the appendix but someone may ask for it
# we may need to show that the community variables are not correlated 
# especially funcitonal specialisation and diversity

################################################################################
# 1) Load community variables
################################################################################

# Load and tidy community variables

cover <- read_csv("results/total_cover_T0.csv") %>% 
  dplyr::select(!experiment) 
cover <- cover %>% group_by(site_name, year_trt, trt) %>% mutate(mean_cover = mean(mean_cover))
unique(cover$site_name) # cover values for 54 sites

# read in the diveristy data
div <- read_csv("results/diversity_metrics_dragnet_T0.csv") # plot / site
#read in the strategy data
strategy <- read_csv("results/community_specialisation_results_T0.csv") # plot / site
# join the data
all_comm <- div %>% full_join(strategy, join_by("site_name", "year_trt", "trt"))
all_comm <- all_comm %>% left_join(cover, join_by("site_name", "year_trt", "trt"))

# tidy some names
all_comm <- 
  all_comm %>% 
  mutate(trt = case_when(trt == "Disturbance" ~ "DIST",
                         trt == "NPK+Disturbance"  ~ "INTER", 
                         trt == "NPK"  ~ "NPK", 
                         trt == "Control" ~ "CONTROL")) %>%
  rename("experiment" = "trt") %>%
  dplyr::select(-time_period)

# read in in the species provenance data
native <- read_csv("results/native_status.csv") %>% 
  dplyr::select(site_name, New_taxon, new_provenance) %>%
  rename(Taxon = New_taxon)

################################################################################
# 2) Visualise community variables: spread values
################################################################################

# check the ranges of values for the interaction modelling
see <- all_comm %>%
  mutate(across(site_time_rich:sigmaonMeans,
                ~ ifelse(is.infinite(.x), NA_real_, .x))) %>%
  group_by(year_trt, experiment) %>%
  summarise(across(site_time_rich:sigmaonMeans, ~paste(range(.x, na.rm = TRUE), collapse = " - ")))


all_comm %>% ggplot(aes(site_time_shan)) +
  geom_density() +
  ggtitle("shannon diversity spread")
# shannon 0.5, 1, 2

all_comm %>% ggplot(aes(meanSigmas)) +
  geom_density()+
  ggtitle("meanSigmas: site functional specialisation")
# meanSigmas 0.3, 0.45, 0.55

all_comm %>% ggplot(aes(variability)) +
  geom_density()+
  ggtitle("variability: site functional variability")
# variability 0.1, 0.25, 0.4


################################################################################
# 3) Correlation of some variables for interpretation of interaction effects
################################################################################

# import the scaled demo and comm vars we used for the interaction modelling
new_comm <- read_csv("results/scaled_modelling_data_interactions.csv")
names(new_comm)
unique(new_comm$model)
unique(new_comm$Taxon)
unique(new_comm$Taxon_site)
unique(new_comm$site_name)
unique(new_comm$experiment)
new_comm <-
  new_comm %>% 
  dplyr::select(site_name, year_trt, experiment, meanSigmas, variability, site_time_invsimp, mean_cover) 

# these two variables across the data set are not well correlated
ggplot(new_comm, aes(meanSigmas, variability, colour = experiment)) +
  geom_point()


