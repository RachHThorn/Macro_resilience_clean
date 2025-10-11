# R Thornley (adapted from Chrissy Hernandez code)
# 03/12/2024
# PCA plot and extract loadings

library(tidyverse)
library(missMDA)
library(ggfortify)

dat <- read_csv("results/DEC_2024/LHT_PCA_raw_data.csv")
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

meta %>%
  dplyr::select(SpeciesAccepted, DRAGnet) %>%
  cbind(loads) %>%
  dplyr::rename("Taxon" = "SpeciesAccepted") %>%
  mutate(Taxon = str_replace(Taxon, " ", "_")) %>%
  filter(DRAGnet == TRUE) %>%
  group_by(Taxon) %>% 
  summarise(across(c(PC1:PC9), mean)) %>%
  write_csv("results/DEC_2024/LHT_PCA_loadings_complete_cases.csv")

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
ggsave("figures/pca_context_plot_BES.jpeg", p, height = 5, width = 5)
# this needs sortomg more but is ok for now

###############################################################################

# can we do this using ggplot
loads <- as.data.frame(pca$x)
plot_dat <- 
  meta %>% # the original data set
  dplyr::select(SpeciesAccepted, DRAGnet) %>%
  cbind(loads) %>% # join to the loadings values
  dplyr::rename("Taxon" = "SpeciesAccepted") %>%
  mutate(Taxon = str_replace(Taxon, " ", "_"))

plot_dat %>%
  ggplot(aes(PC1, PC2, colour = DRAGnet))+
  geom_point()
# we get our coloured dots ok here but not our arrows

###############################################################################

#https://clauswilke.com/blog/2020/09/07/pca-tidyverse-style/
# copied from the above blog - but it doesn't do arrows and dots on the same plot at the 
# moment - so I will stick with the ggfortify package

library(broom)

arrow_style <- arrow(
  angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
)
pca %>% 
  tidy(matrix = "rotation") %>%
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(xend = 0, yend = 0, arrow = arrow_style) +
  geom_text(
    aes(label = column),
    hjust = 1, nudge_x = -0.02, 
    color = "#904C2F"
  ) +
  xlim(-1.25, .5) + ylim(-.5, 1) +
  coord_fixed()+ # fix aspect ratio to 1:1
  theme_bw() 


