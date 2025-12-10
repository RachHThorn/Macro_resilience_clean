# R Thornley (adapted from Chrissy Hernandez code)
# 03/12/2024
# PCA plot and extract loadings

library(tidyverse)
library(missMDA)
library(ggfortify)

dat <- read_csv("results/LHT_PCA_raw_data.csv")
names(dat)

meta <- dat %>% dplyr::select(1:28)
LHT <- dat %>% dplyr::select(29:37)

################################################################################

# method 1: a traditional PCA on complete cases only

################################################################################

# Calculate the PCA:
pca <- prcomp(LHT, center = TRUE, scale. = TRUE)

# the amount of total variance in your data that these loadings represent
variancePCA <- summary(pca)$importance[2,]
variancePCA
variancePCA*100 

# visualise this as a scree plot
variancePCA %>% 
  t() %>% 
  as.data.frame() %>% 
  pivot_longer(cols= PC1:PC9) %>%
  mutate(value = value*100) %>%
  mutate(group = "1") %>%
  ggplot(aes(name, value, group = group))+
  geom_point()+
  geom_path()+
  ggtitle("scree plot") +
  xlab("")+
  ylab("percent variance captured")


# LHT loadings for PC1 and PC2:
pca$rotation[,1:2] # get the loadings of PC 1 and 2
# get all loadings and save them as df
loads <- as.data.frame(pca$x)

dim(meta)
dim(loads)

# output the data we need for modelling for just the DRAGnet species
# we only have 22 DRAGNet species here

meta <- meta %>%
  dplyr::select(SpeciesAccepted, DRAGnet) %>%
  cbind(loads) %>%
  dplyr::rename("Taxon" = "SpeciesAccepted") %>%
  mutate(Taxon = str_replace(Taxon, " ", "_")) %>%
  filter(DRAGnet == TRUE) %>%
  group_by(Taxon) %>% 
  summarise(across(c(PC1:PC9), mean)) 

write_csv(meta, "results/LHT_PCA_loadings_complete_cases.csv")

################################################################################
# Visualise
################################################################################

# visualise the PCA using ggfortify
# Prepare the data
names(dat)
LHT <- dat[, 29:37]

# Principal component analysis
# Calculate the PCA:
pca <- prcomp(LHT, center = TRUE, scale. = TRUE)
# Plot
p <- autoplot(pca, x = 2, y = 1, loadings = TRUE, loadings.label = TRUE,
         data = dat, colour = 'DRAGnet',
         loadings.colour = "black",
         loadings.label.colour = "black",
         loadings.label.size = 5,
         label.vjust = 5)+ 
  theme_classic()+
  theme(text = element_text(size = 20))+
  theme(legend.title = element_blank())+
  theme(legend.position = "top")+
  scale_colour_manual(values = c("green2","palevioletred4"), labels = c("Non-DRAGNet", "DRAGNet"))
p
ggsave("figures/pca_context_plot_S1c.jpeg", p, height = 5, width = 5)

################################################################################
# Calculate PCA from incomplete cases
################################################################################

dat <- read_csv("results/LHT_PCA_loadings_imputed.csv")
names(dat)
unique(dat$SpeciesAccepted) # 600 species

dat %>% filter(DRAGnet == TRUE) %>% summarise(count = n_distinct(SpeciesAccepted)) 
# why do we only have 38 species???

dat %>% filter(DRAGnet == TRUE) %>% group_by(SpeciesAccepted) %>% tally() %>% arrange(-n)
# but we have a lot of matrices for these species

# just tidy this up a bit
# we need to change any inf values to NA
# dat <- dat %>% mutate(across(29:37, ~na_if(abs(.), Inf)))

output %>% filter(SpeciesAccepted %in% Drag_taxa_new)

# we need to get rid of ridiculous data here and put NA values instead
check <-
  output %>% 
  dplyr::select(SpeciesAccepted, DRAGNet, T_generation:Lmax) %>%
  mutate(ID = row_number()) %>%
  pivot_longer(cols = T_generation:Lmax, names_to = "demo_var", values_to = "value") %>%
  mutate(value = case_when(demo_var == "T_generation" & value > 100 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "R0" & value > 2000 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "percapita_repro" & value > 2000 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "Lmean" & value > 200 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "age_repro" & value > 20 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "Lambda" & value > 2 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "Reactivity" & value > 5000 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "FirstStepAtt" & value > 5000 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "L_max" & value < 80 ~ NA, TRUE ~ value))



# %>% group_by(ID) %>% pivot_wider(names_from = demo_var, values_from = value)


# we can then get the mean values per species
# probably best to get the mean values per species and then look to impute the missing cases
dat <- 
  dat %>% 
  filter(DRAGnet == TRUE) %>%
  group_by(SpeciesAccepted) %>% 
  summarise(across(T_generation:Lmax, ~ mean(., na.rm = TRUE)))

# split the data into the numeric bit we pass to the pca and the labels we need for later
meta <- dat %>% dplyr::select(1)
LHT <- dat %>% dplyr::select(2:9)

###############################################################################

# Method 2: non complete cases PCA
#https://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html#:~:text=PCA%20with%20missing%20values,-Then%2C%20before%20modeling&text=The%20package%20missMDA%20allows%20the,on%20the%20imputed%20data%20set
# method taken from this blog

# not sure what this bit is doing here
nb <- estim_ncpPCA(LHT, method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
# (available methods include GCV to approximate CV)
nb$ncp 
plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")

res.comp <- imputePCA(LHT) # iterativePCA algorithm
res.comp
res.comp$completeObs[1:3,]
imp <- cbind.data.frame(res.comp$completeObs)

# now run the pca on the imputed data
pca <- prcomp(imp, center = TRUE, scale. = TRUE)
variancePCA <- summary(pca)$importance[2,]

# visualise variance per pca as a scree plot
variancePCA %>% 
  t() %>% 
  as.data.frame() %>% 
  pivot_longer(cols= PC1:PC8) %>%
  mutate(value = value*100) %>%
  mutate(group = "1") %>%
  ggplot(aes(name, value, group = group))+
  theme_bw()+
  geom_point()+
  geom_path()+
  ggtitle("scree plot") +
  xlab("")+
  ylab("percent variance captured")

################################################################################
# visualise the PCA using ggfortify
# Prepare the data

autoplot(pca, x = 2, y = 1, loadings = TRUE, loadings.label = TRUE,
         data = dat, 
         loadings.colour = "black",
         loadings.label.colour = "black")+ 
  theme_classic()+
  theme(legend.title = element_blank())+
  theme(legend.position = "top")
# this needs sorting more but is ok for now

################################################################################

# get all loadings and save them as df
loads <- as.data.frame(pca$x)

meta %>%
  dplyr::select(SpeciesAccepted) %>%
  cbind(loads) %>%
  dplyr::rename("Taxon" = "SpeciesAccepted") %>%
  mutate(Taxon = str_replace(Taxon, " ", "_")) %>%
  write_csv("results/DEC_2024/LHT_PCA_loadings_imputed.csv")

