# R Thornley (adapted from Chrissy Hernandez code)
# 11/12/2024
# PCA plot and extract loadings

rm(list = ls())

library(tidyverse)
library(missMDA)
library(ggfortify)

dat <- read_csv("results/all_COMPADRE_metrics.csv")
dat <- dat %>% pivot_wider(values_from = "value", names_from = "demo_var")

# split the pcs data from the meta data
meta <- dat %>% dplyr::select(1:4)
LHT <- dat %>% dplyr::select(5:13)

# Calculate the PCA:
LHT <- LHT %>% drop_na()
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

################################################################################

# visualise the PCA using ggfortify
small_dat <- dat %>% drop_na()
names(small_dat)
small_dat$DRAGNet


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
ggsave("figures/pca_context_plot_BES.jpeg", p, height = 5, width = 5)
# this needs sorting more but is ok for now
