# R Thornley
# 30/08/2025
# Visualisation of the demographic metrics

library(tidyverse)
library(ggpubr)

################################################################################
# LOAD DATA AND TIDY
################################################################################

# Load demography variables and find the mean and SE of each Taxon and trait
demo <- read_csv("results/all_COMPADRE_metrics.csv") %>% 
  filter(DRAGNet == TRUE) %>% 
  rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var", "Demo_value" = "value") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  select(Taxon, Demo_trait, Demo_value) %>%
  group_by(Taxon, Demo_trait) %>%
  mutate(Sample_size = n(),
         s = sd(Demo_value, na.rm = TRUE),
         se = s/sqrt(Sample_size),
         Demo_value_mean = mean(Demo_value, na.rm = TRUE),
         Demo_value_log = log(Demo_value),
         Taxon_label = str_replace_all(Taxon, "_", " "),
         Taxon_label = paste0(Taxon_label, " n = ", Sample_size))

# list the demography traits we are dealng with
print(unique(demo$Demo_trait))
# 1) "T_generation"    
# 2) "R0"              
# 3) "percapita_repro" 
# 4) "Lmean"           
# 5) "age_repro"       
# 6) "Lambda"         
# 7) "Reactivity"      
# 8) "FirstStepAtt"    
# 9) "Lmax"  

# see if there are zero values in this data set
demo %>%
  group_by(Demo_trait, Taxon) %>%
  summarise(n_total = n(), n_zeros = sum(Demo_value == 0, na.rm = TRUE),
            prop_zeros = mean(Demo_value == 0, na.rm = TRUE)) %>%
  filter(n_zeros > 0)
# there are a couple
demo %>% filter(Demo_value == 0)

################################################################################
# BOXPLOTS SHOWING SPREAD OF DATA FOR EACH OF THE DEMO VARS
################################################################################

# "T_generation"   
p1 <-  
  demo %>%
  filter(Demo_trait == "T_generation") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Generation Time")
p1
# ggsave("figures/T_generation_boxplot.jpeg", p1, height = 5, width = 6)

# "R0" 
p2 <-  
  demo %>%
  filter(Demo_trait == "R0") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reproductive output")
p2
# ggsave("figures/Reproductive_rate_boxplot.jpeg", p2, height = 5, width = 6)

# "percapita_repro" 
p3 <-  
  demo %>%
  filter(Demo_trait == "percapita_repro") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Per capita reproduction")
p3
# ggsave("figures/Percapita_reproduction_boxplot.jpeg", p3, height = 5, width = 6)

# Lmean
p4 <-  
  demo %>%
  filter(Demo_trait == "Lmean") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Mean Life Expectancy")
p4
# ggsave("figures/Mean_life_expectancy_boxplot.jpeg", p4, height = 5, width = 6)

# "age_repro"
demo$Demo_trait
p5 <-  
  demo %>%
  filter(Demo_trait == "age_repro") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Age at first reproduction")
p5
# ggsave("figures/Reproductive_output_boxplot.jpeg", p5, height = 5, width = 6)

# "Lambda"   
p6 <- 
  demo %>%
  filter(Demo_trait == "Lambda") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Lambda")
p6
# ggsave("figures/Lambda_boxplot.jpeg", p6, height = 5, width = 6)

# "Reactivity"    
p7 <- 
  demo %>%
  filter(Demo_trait == "Reactivity") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reactivity")
p7
# ggsave("figures/Reactivity_boxplot.jpeg", p7, height = 5, width = 6)

#  "FirstStepAtt"   
p8 <- 
  demo %>%
  filter(Demo_trait == "FirstStepAtt") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("FirstStepAtt")
p8
# ggsave("figures/FirstStepAtt_boxplot.jpeg", p8, height = 5, width = 6)

# "Lmax" 
p9 <- 
  demo %>%
  filter(Demo_trait == "Lmax") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Maximum Life Expectancy")
p9
# ggsave("figures/Maximum_life_expectancy_boxplot.jpeg", p9, height = 5, width = 6)

combined <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,
                      ncol = 3, nrow = 3,
                      common.legend = FALSE)

combined <- annotate_figure(
  combined,
  bottom = text_grob("Taxon", color = "black", face = "bold", size = 14),
  left = text_grob("Demographic Metric", color = "black", face = "bold", size = 14, rot = 90)
)
combined
# ggsave("figures/all_demo_metrics_boxplots.jpeg", combined, height = 18, width = 12)

################################################################################
# BOXPLOTS: VISUALISE THE VALUES LOGGED
################################################################################

# "T_generation"   
p1_again <-  
  demo %>%
  filter(Demo_trait == "T_generation") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Generation Time")
p1_again
# ggsave("figures/T_generation_boxplot.jpeg", p1, height = 5, width = 6)

# "R0" 
p2_again <-  
  demo %>%
  filter(Demo_trait == "R0") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reproductive output (log)")
p2_again

# "percapita repro"
p3_again <-  
  demo %>%
  filter(Demo_trait == "percapita_repro") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Per capita reproduction")
p3_again

# "Lmean"
p4_again <-  
  demo %>%
  filter(Demo_trait == "Lmean") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Mean Life Expectancy (log)")
p4_again

# "age_repro"
p5_again <-  
  demo %>%
  filter(Demo_trait == "age_repro") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Age at first reproduction")
p5_again
# ggsave("figures/Reproductive_output_boxplot.jpeg", p5, height = 5, width = 6)

# "Lambda"   
p6_again <- 
  demo %>%
  filter(Demo_trait == "Lambda") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Lambda")
p6_again
# ggsave("figures/Lambda_boxplot.jpeg", p6, height = 5, width = 6)

# "Reactivity"    
p7_again <- 
  demo %>%
  filter(Demo_trait == "Reactivity") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reactivity (log)")
p7_again

#  "FirstStepAtt"   
p8_again <- 
  demo %>%
  filter(Demo_trait == "FirstStepAtt") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("FirstStepAtt")
p8_again

# "Lmax"
p9_again <- 
  demo %>%
  filter(Demo_trait == "Lmax") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Maximum Life Expectancy (log)")
p9_again


# Combine the logged and no-logged plots 
# NOTE (may need to re-order and rename these)
combined <- ggarrange(p1_again, p2_again, p3_again, p4_again, 
                      p5_again, p6_again, p7_again, p8_again, p9_again,
                      ncol = 3, nrow = 3,
                      common.legend = FALSE)

combined <- annotate_figure(
  combined,
  bottom = text_grob("Taxon", color = "black", face = "bold", size = 14),
  left = text_grob("Demographic Metric", color = "black", face = "bold", size = 14, rot = 90)
)
combined
ggsave("figures/all_demo_metrics_logged_boxplots.jpeg", combined, height = 18, width = 12)


################################################################################
# HISTOGRAMS OF MEAN DEMO VARIABLES
################################################################################

# take the mean of the non-logged values
demo$Demo_value_mean
small_demo <- 
  demo %>% 
  select(Taxon, Demo_trait, Demo_value_mean) %>%
  filter(!is.infinite(Demo_value_mean)) %>%
  group_by(Taxon, Demo_trait) %>% 
  summarise(Demo_value_mean = mean(Demo_value_mean, na.rm = TRUE), .groups = "drop") %>%
  drop_na()

# Non-transformed mean values
Means_density <-  
  small_demo %>%
  ggplot(aes(Demo_value_mean))+
  geom_density()+
  theme_classic()+
  facet_wrap(~Demo_trait, scales = "free")+
  ggtitle("Density of Demo var Means")
Means_density
# ggsave("figures/Non_transformed_means_demo_vars_density.jpeg", Means_density, height = 5, width = 6)

# take the mean of the logged values
small_demo <- 
  demo %>% 
  select(Taxon, Demo_trait, Demo_value_log) %>%
  filter(!is.infinite(Demo_value_log)) %>%
  group_by(Taxon, Demo_trait) %>% 
  summarise(Demo_value_log_mean = mean(Demo_value_log, na.rm = TRUE), .groups = "drop") %>%
  drop_na()

# Logged mean values
Logged_means_density <-  
  small_demo %>%
  ggplot(aes(Demo_value_log_mean))+
  geom_density()+
  theme_classic()+
  facet_wrap(~Demo_trait, scales = "free")+
  ggtitle("Density of Demo var Logged Means")
Logged_means_density
# ggsave("figures/Logged_means_demo_vars_density.jpeg", Logged_means_density, height = 5, width = 6)

################################################################################
# NUMBER OF ZEROS IN DATA FOR MODEL SUITABILITY
################################################################################

# Do some of the demography variables have a lot of zeros
# There are no zeros in the means
names(demo)
demo %>% group_by(Demo_trait) %>% 
  summarise(n_total = n(),
  n_zeros = sum(Demo_value_mean == 0, na.rm = TRUE),
  prop_zeros = mean(Demo_value_mean == 0, na.rm = TRUE))
  