# R Thornley
# 18/11/2026
# Project: P1_COMPADRE_DRAGNET
# Script: F2_demographic_variation
# Variation in demographic metrics across species and PCA LH trait space

rm(list = ls())

################################################################################
# INSTRUCTIONS
################################################################################
# 1) CREATE FIGURE WIDE LABELS / COLOURS
# 2) LOAD DEMOGRAPHIC METRICS AND THE GROWTH FORM DATA
# 3) CREATE DENSITY PLOTS WITH ALL DEMOGRAPHIC VALUES
# 4) LH TRAIT SPACE PCA 
# 5) BRING PANELS TOGETHER IN FINAL PLOT

################################################################################
# 1) CREATE FIGURE WIDE LABELS / COLOURS
################################################################################

trait_labels <- c(
  T_generation = "Generation time",
  R0 = "Net reproductive rate",
  percapita_repro = "Per-capita reproduction",
  Lmean = "Mean life expectancy",
  age_repro = "Age at first reproduction",
  Lmax = "Maximum longevity",
  Lambda = "Population growth rate",
  Reactivity = "Reactivity",
  FirstStepAtt = "First-step attenuation"
)
# Match the trait transformations in the figure-preparation script.
log_traits <- c("T_generation", "R0", "percapita_repro", "Lmean",
                "age_repro", "Lmax", "Reactivity")
growth_colours <- c(
  forb_annual = "#D55E00", forb_biennial = "#CC79A7",
  forb_perennial = "#0072B2", graminoid_perennial = "#009E73",
  woody_perennial = "#555555"
)
growth_shapes <- setNames(c(21, 24, 22, 23, 25), names(growth_colours))
growth_labels <- c("Annual forb", "Biennial forb", "Perennial forb",
                   "Perennial graminoid", "Woody perennial")

################################################################################
# 2) LOAD DEMOGRAPHIC METRICS AND THE GROWTH FORM DATA
################################################################################

demo <- read.csv(file.path("results/all_COMPADRE_metrics.csv"))
gift <- read.csv(file.path("data/GIFT_data_updated_for_PCA.csv"))
# Amend the spelling differences for one species
gift$Species[gift$Species == "Andropogon gerardi"] <- "Andropogon gerardii"
stopifnot(!anyDuplicated(gift$Species))

# Use all available finite filtered values for each trait.
demo <- subset(demo, DRAGNet & demo_var %in% names(trait_labels) & is.finite(value))
stopifnot(!anyDuplicated(demo[c("SpeciesAccepted", "MatrixID", "demo_var")]))
species_means <- aggregate(value ~ SpeciesAccepted + demo_var, demo, mean)
species_means$plot_value <- species_means$value
use_log <- species_means$demo_var %in% log_traits
species_means$plot_value[use_log] <- log1p(species_means$value[use_log])
stopifnot(all(is.finite(species_means$plot_value)))
species_means$Final_growth_form <- gift$Final_growth_form[
  match(species_means$SpeciesAccepted, gift$Species)
]
stopifnot(all(species_means$Final_growth_form %in% names(growth_colours)))

# The number of species can differ by trait.
counts <- table(factor(species_means$demo_var, levels = names(trait_labels)))
species_means$trait <- factor(
  species_means$demo_var, levels = names(trait_labels),
  labels = sprintf("%s   ·   %s   ·   n = %d", trait_labels,
                   ifelse(names(trait_labels) %in% log_traits, "log(1 + x)", "original scale"),
                   counts)
)
# Retain every finite filtered matrix value for the density curves.
density_values <- demo
density_values$plot_value <- density_values$value
density_use_log <- density_values$demo_var %in% log_traits
density_values$plot_value[density_use_log] <-
  log1p(density_values$value[density_use_log])

# Use exactly the same facet labels as the species-mean dots.
# The n in these labels continues to count species.
density_values$trait <- species_means$trait[
  match(density_values$demo_var, species_means$demo_var)
]
stopifnot(all(is.finite(density_values$plot_value)),
          !anyNA(density_values$trait))

# Offset points vertically to show overlaps without moving their trait values.
species_means$dot_y <- -0.10
for (panel in levels(species_means$trait)) {
  rows <- which(species_means$trait == panel)
  rows <- rows[order(species_means$plot_value[rows])]
  species_means$dot_y[rows] <- rep(c(-0.06, -0.14, -0.22), length.out = length(rows))
}

################################################################################
# 3) CREATE DENSITY PLOTS WITH ALL DEMOGRAPHIC VALUES
################################################################################

filtered_distributions <- ggplot(species_means, aes(x = plot_value)) +
  # CHANGED: density curves use matrix values; geom_point uses species_means.
  geom_density(data = density_values,
               aes(y = after_stat(scaled)), fill = "grey88",
               colour = "grey45", linewidth = 0.45, trim = TRUE) +
  geom_point(aes(y = dot_y, fill = Final_growth_form, shape = Final_growth_form),
             size = 1.65, colour = "grey20", stroke = 0.35) +
  facet_wrap(vars(trait), ncol = 1, scales = "free_x") +
  scale_fill_manual(values = growth_colours, breaks = names(growth_colours),
                    labels = growth_labels, name = "Growth form") +
  scale_shape_manual(values = growth_shapes, breaks = names(growth_shapes),
                     labels = growth_labels, name = "Growth form") +
  guides(fill = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  coord_cartesian(ylim = c(-0.34, 1.08)) +
  labs(title = "Demographic metric variation among species",
       subtitle = "Filtered COMPADRE data",
       # CHANGED: the x-axis now covers matrix values and species means.
       x = "Metric value (scale shown in each panel)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.line.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), strip.background = element_blank(),
    strip.text = element_text(hjust = 0, face = "bold"),
    legend.position = "bottom", legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

# preview part A of the figure - denisty curves
filtered_distributions

################################################################################
# 4) LH TRAIT SPACE PCA 
################################################################################

# Fit the six-trait PCA as in the comparison script, retaining its orientation.
trait_labels <- c(
  T_generation = "Generation time",
  R0 = "Net reproductive rate",
  percapita_repro = "Per capita reproduction",
  Lmean = "Mean life expectancy",
  age_repro = "Age at first reproduction",
  Lmax = "Maximum longevity"
)
traits <- names(trait_labels)
# As in the companion figure, transform after taking arithmetic species means.
log_traits <- c("age_repro", "Lmax", "Lmean", "percapita_repro",
                "R0", "T_generation")

filtered <- read.csv(file.path("results/all_COMPADRE_metrics.csv"))
unfiltered <- read.csv(file.path("results/unfiltered_COMPADRE_metrics.csv"))
filtered <- subset(filtered, DRAGNet & demo_var %in% traits)
unfiltered <- unfiltered[unfiltered$DRAGNet, c("MatrixID", "SpeciesAccepted", traits)]

# Compare the same source matrices. The filtered file stores one row per metric.
stopifnot(!anyDuplicated(filtered[c("MatrixID", "demo_var")]),
          !anyDuplicated(unfiltered$MatrixID),
          setequal(filtered$MatrixID, unfiltered$MatrixID))
matrix_rows <- match(filtered$MatrixID, unfiltered$MatrixID)
stopifnot(all(filtered$SpeciesAccepted == unfiltered$SpeciesAccepted[matrix_rows]))
filtered_values <- matrix(NA_real_, nrow(unfiltered), length(traits),
                          dimnames = list(NULL, traits))
filtered_values[cbind(matrix_rows, match(filtered$demo_var, traits))] <- filtered$value
filtered_wide <- data.frame(unfiltered[c("MatrixID", "SpeciesAccepted")], filtered_values)

fit_pca <- function(data) {
  species_means <- aggregate(
    data[traits], list(species = data$SpeciesAccepted),
    function(values) mean(values[is.finite(values)])
  )
  input <- as.matrix(species_means[traits])
  rownames(input) <- species_means$species
  transformed <- intersect(log_traits, traits)
  input[, transformed] <- log1p(input[, transformed])
  keep <- rowSums(!is.finite(input)) == 0
  prcomp(input[keep, , drop = FALSE], center = TRUE, scale. = TRUE)
}

fits <- list(Filtered = fit_pca(filtered_wide), Unfiltered = fit_pca(unfiltered))
# PCA axis signs are arbitrary. Orient each unfiltered axis by the corresponding
# loading vector so a sign flip alone doesn't look like a change.
signs <- ifelse(colSums(fits$Filtered$rotation[, 1:2] *
                        fits$Unfiltered$rotation[, 1:2]) < 0, -1, 1)
fits$Unfiltered$x[, 1:2] <- sweep(fits$Unfiltered$x[, 1:2], 2, signs, "*")
fits$Unfiltered$rotation[, 1:2] <- sweep(fits$Unfiltered$rotation[, 1:2], 2, signs, "*")

pca_unfiltered <- fits$Unfiltered
pca_scores <- data.frame(pca_unfiltered$x[, 1:2],
                         Species = rownames(pca_unfiltered$x), dataset = "Unfiltered")
pca_scores$Final_growth_form <- gift$Final_growth_form[
  match(pca_scores$Species, gift$Species)
]
stopifnot(!anyNA(pca_scores$Final_growth_form))

fits <- list(Unfiltered = pca_unfiltered)
scores <- pca_scores
scores$dataset <- factor(scores$dataset, levels = names(fits))
growth_colours <- c(forb_annual = "#D55E00", forb_biennial = "#CC79A7",
                    forb_perennial = "#0072B2", graminoid_perennial = "#009E73",
                    woody_perennial = "#555555")
growth_labels <- c("Annual forb", "Biennial forb", "Perennial forb",
                   "Perennial graminoid", "Woody perennial")

trait_labels <- c(T_generation = "Generation time", R0 = "Net reproductive rate",
                  percapita_repro = "Per-capita reproduction",
                  Lmean = "Mean life expectancy", age_repro = "Age at first reproduction",
                  Lmax = "Maximum longevity")

loadings <- do.call(rbind, lapply(names(fits), function(dataset) {
  rotation <- fits[[dataset]]$rotation[names(trait_labels), 1:2]
  data.frame(rotation, dataset, metric = seq_along(trait_labels))
}))
loadings$dataset <- factor(loadings$dataset, levels = names(fits))
arrow_scale <- 0.65 * min(max(abs(scores$PC1)) / max(abs(loadings$PC1)),
                          max(abs(scores$PC2)) / max(abs(loadings$PC2)))
loadings$x <- loadings$PC1 * arrow_scale
loadings$y <- loadings$PC2 * arrow_scale
loadings$label_y <- loadings$y * 1.12
pca_growth_forms <- ggplot(scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.35) +
  geom_segment(
    data = loadings, aes(x = 0, y = 0, xend = x, yend = y),
    inherit.aes = FALSE, colour = "grey30", linewidth = 0.45,
    arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")
  ) +
  geom_point(aes(fill = Final_growth_form, shape = Final_growth_form), size = 2,
             colour = "grey25", stroke = 0.35) +
  geom_label(
    data = loadings, aes(x * 1.12, label_y, label = metric),
    inherit.aes = FALSE, size = 8.5, size.unit = "pt", colour = "grey20",
    label.padding = grid::unit(0.07, "lines"), linewidth = 0
  ) +
  scale_fill_manual(values = growth_colours, breaks = names(growth_colours),
                    labels = growth_labels, na.value = "grey70") +
  scale_shape_manual(values = growth_shapes, breaks = names(growth_shapes),
                     labels = growth_labels, name = "Growth form") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 2.9)),
         shape = guide_legend(nrow = 3, byrow = TRUE)) +
  coord_equal() +
  labs(title = "B   Multivariate life-history space",
       subtitle = "Unfiltered COMPADRE data (n = 39 species)",
       x = sprintf("PC1 (%.1f%%)", 100 * pca_unfiltered$sdev[1]^2 / sum(pca_unfiltered$sdev^2)),
       y = sprintf("PC2 (%.1f%%)", 100 * pca_unfiltered$sdev[2]^2 / sum(pca_unfiltered$sdev^2)),
       fill = "Growth form",
       caption = paste(vapply(
         split(sprintf("%d  %s", seq_along(trait_labels), trait_labels),
               rep(1:3, each = 2)),
         paste, collapse = "    ", FUN.VALUE = character(1)
       ), collapse = "\n")) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top", legend.title.position = "top",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 8.5),
        legend.key.width = grid::unit(4, "mm"),
        legend.spacing.x = grid::unit(1, "mm"),
        plot.title = element_text(face = "bold"),
        panel.spacing = grid::unit(7, "mm"),
        plot.caption = element_text(hjust = 0, size = 8.5))

# Keep both headings outside the plots so the PCA aspect ratio cannot move them.
panel_heading <- function(title, subtitle) {
  grid::grobTree(
    grid::textGrob(title, x = grid::unit(2, "mm"), y = grid::unit(8, "mm"),
                   just = "left", gp = grid::gpar(fontsize = 12, fontface = "bold")),
    grid::textGrob(subtitle, x = grid::unit(2, "mm"), y = grid::unit(3, "mm"),
                   just = "left", gp = grid::gpar(fontsize = 9))
  )
}

# preview the PCS plot
pca_growth_forms

################################################################################
# 5) BRING PANELS TOGETHER IN FINAL PLOT
################################################################################

filtered_distributions <- filtered_distributions +
  labs(title = NULL, subtitle = NULL) +
  theme(legend.position = "none",
        strip.text = element_text(size = 8.5, hjust = 0, face = "bold"),
        axis.text = element_text(size = 8.5),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5))

pca_growth_forms <- pca_growth_forms +
  labs(title = NULL, subtitle = NULL) +
  theme(axis.text = element_text(size = 8.5),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5))


panel_heading <- function(title, x_mm = 2) {
  grid::textGrob(
    title,
    x = grid::unit(x_mm, "mm"),
    y = grid::unit(1, "npc") - grid::unit(1, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontface = "bold", lineheight = 1)
  )
}

heading_a <- panel_heading(
  "A   Demographic metric variation\n      among species"
)

heading_b <- panel_heading(
  "B   Multivariate life-history space",
  x_mm = 12
)


# Reserve the lower right for the journal caption.
right <- gridExtra::arrangeGrob(
  pca_growth_forms, grid::nullGrob(), ncol = 1, heights = c(3, 1)
)
combined <- gridExtra::arrangeGrob(
  heading_a, heading_b, filtered_distributions, right,
  ncol = 2, widths = c(1, 1.12),
  heights = grid::unit.c(grid::unit(11, "mm"), grid::unit(1, "null"))
)

# Preview in the RStudio Plots pane when running the script interactively.
if (interactive()) {
  grid::grid.newpage()
  grid::grid.draw(combined)
}


# Export directly at the intended print size, with a separate raster proof.
ggsave("figures/F2_180mm.pdf", combined, width = 180, height = 170, units = "mm")
ggsave("figures/F2_180mm.tiff", combined, width = 180, height = 170,
       units = "mm", dpi = 600, bg = "white", compression = "lzw")
ggsave("figures/F2_180mm.jpg", combined, width = 180, height = 170,
       units = "mm", dpi = 300, bg = "white", quality = 95)
