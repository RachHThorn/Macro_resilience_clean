# R Thornley
# 15/07/2025
# A2: random effects visualise 
# random effects plots to compare the GLMM results in the context of all the data
# NOTE: this has been simplified for the paper to show only T0-T1 for the Ordered Beta models

library(tidyverse)
library(ggpubr)
library(ggrepel)

################################################################################

# read in effect sizes from the whole of DRAGNet to contextualise the effects from
taxa <- read_csv("results/RE_SE_Taxon_all_DRAGNet.csv")

# read in the taxon RE from the compadre and dragnet overlap
comp_taxa <- read_csv("results/common_species_drag_comp.csv") %>% pull(x)

# select the random effects for the var of interest only (New_taxon)
# create a new ID variable
# change some labels in the data frame to tidy it up
taxa <- 
  taxa %>% 
  filter(str_detect(group_var, "New_taxon")) %>%
  mutate(ID = if_else(group %in% comp_taxa, "1", "0")) %>%
  mutate(ID = factor(ID, levels = c("0", "1"))) %>%
  mutate(time_period = case_when(time_period == 1 ~ "T0-T1",
                                        time_period == 2 ~ "T0-T2",
                                        time_period == 3 ~ "T0-T3")) %>%
  dplyr::select(group, value, se, model, time_period, experiment, ID) %>%
  rename(taxon = group) %>%
  mutate(label = if_else(ID == "1", taxon, NA_character_))

################################################################################

# make sure the labels for the experiments are ordered correctly
custom_labels = c(DIST = "Disturbance", NPK = "NPK")

# put the labels on the y axis instead
# this is not perfect - trying to get the y axis non overlapping
taxa_labels <- taxa %>% filter(model == "Ordbeta", time_period == "T0-T1", ID == "1")

plot <- 
  taxa %>% 
  filter(model == "Ordbeta") %>%
  filter(time_period == "T0-T1") %>%
  filter(experiment %in% c("DIST", "NPK")) %>% 
  mutate(experiment = factor(experiment, levels = c("DIST", "NPK"))) %>%
  ggplot(aes(x = value, y = reorder(taxon, value), colour = ID)) +
  theme_classic() +
  geom_point(aes(colour = ID), size = 0.3) +
  geom_errorbar(aes(xmin = value - se, xmax = value + se, colour = ID), size = 0.3) +
  # Only show y-axis labels for taxa with label == 1
  scale_y_discrete(breaks = taxa_labels$taxon) +
  ylab("Taxon") +
  xlab("Standardised effects with SE") +
  facet_wrap(~experiment, scales = "free_y", labeller = labeller(experiment = custom_labels)) +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  scale_colour_manual(values = c("grey", "red"), labels = c("all DRAGNet", "COMPADRE overlap")) +
  theme(legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(colour = "red", size = 5.5))
  
plot
ggsave("figures/S1d_GLMM_effect_size_plot.jpeg", height = 8, width = 10)

###############################################################################





plot <- 
  taxa %>% filter(model == "Ordbeta") %>%
  filter(time_period == "T0-T1") %>%
  filter(experiment %in% c("DIST", "NPK")) %>% 
  mutate(experiment = factor(experiment, levels = c("DIST", "NPK"))) %>%
  ggplot(aes(value, reorder(taxon, value)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = value - se, xmax = value + se, colour = ID), size = 0.3)+
  geom_text_repel(
    aes(label = label), size = 3, 
    nudge_x = -0.1,
    direction = "y",
    hjust = 1,
    na.rm = TRUE) +  # labels only where ID==1
  ylab("Taxon") +
  xlab("Standardised effects with SE")+
  facet_wrap(~experiment, scales = "free_y", labeller = labeller(experiment = custom_labels))+
  geom_vline(aes(xintercept = 0), colour ="red", linetype = "dashed")+
  theme(axis.text.y = element_blank())+
  scale_colour_manual(values = c("black", "red"), labels = c("all DRAGNet", "COMPADRE overlap"))+
  theme(legend.title = element_blank(), axis.ticks.y = element_blank(), legend.position = "top")
plot
