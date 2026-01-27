# R Thornley (adapted from Chrissy Hernandez code)
# 11/12/2024
# PCA plot and extract loadings

rm(list = ls())

# Load libraries
library(tidyverse)
library(ggfortify)


# read in COMPADRE metrics
dat <- read_csv("results/all_COMPADRE_metrics.csv") %>%
  mutate(value = as.numeric(value)) %>%     # ensure numeric
  pivot_wider(
    names_from = demo_var,
    values_from = value,
    values_fn = mean                          # if duplicates exist
  )

# split meta + PCA columns
meta_cols <- names(dat)[1:4]
lht_cols  <- names(dat)[5:13]

# make ONE dataset used for both PCA + plotting
small_dat <- dat %>%
  drop_na(all_of(lht_cols))                  # only drop NA in PCA vars

# run PCA on the same rows you will plot
pca <- prcomp(small_dat[, lht_cols], center = TRUE, scale. = TRUE)

# optional: make sure DRAGNet is a factor with 2 levels
# (adjust column name here to match your file exactly)
small_dat <- small_dat %>%
  mutate(DRAGNet = factor(DRAGNet, levels = c(0, 1), labels = c("Non-DRAGNet","DRAGNet")))
# If your column is actually DRAGnet (lowercase n), change DRAGNet -> DRAGnet above and below.

# plot
p1 <- autoplot(
  pca,
  x = 1, y = 2,
  data = small_dat,
  colour = "DRAGNet",
  loadings = TRUE,
  loadings.label = TRUE,
  loadings.colour = "black",
  loadings.label.colour = "black",
  loadings.label.size = 5
) +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "top") +
  scale_colour_manual(values = c("green2","palevioletred4"))

p1

# this plot has the default 
ggsave("figures/pca_context_plot_rough_labels.jpeg", p1, height = 5, width = 5)

p2 <- autoplot(
  pca,
  x = 1, y = 2,
  data = small_dat,
  colour = "DRAGNet",
  loadings = TRUE,          # keeps the arrows
  loadings.label = FALSE,  # removes text labels
  loadings.colour = "black"
) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  scale_colour_manual(values = c("green2","palevioletred4"))

p2

# this plot has empty labels for putting in powerpoint manually (as symbols are annoying)
ggsave("figures/pca_context_plot_no_labels.jpeg", p2, height = 5, width = 5)
