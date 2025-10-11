# R Thornley
# 09/09/2025
# Project: P1_COMPADRE_DRAGNET
# Script: F3_robust_OLS_effects

################################################################################
# Instructions
################################################################################

# 1) LOAD and TIDY uni-variate robust model results data / tidy names
# 2) BUILD FIGURE
# 3) SAVE FIGURE

###############################################################################
# 1) LOAD and TIDY uni-variate robust model results data / tidy names
###############################################################################

dat <- readRDS("results/robust_model_results.rds") %>%
  dplyr::select(time_period, experiment, model, Demo_trait, results) %>%
  tidyr::unnest(results) %>%
  dplyr::filter(term == "Demo_value") %>%
  filter(model == "Ordbeta", experiment %in% c("DIST", "NPK")) %>%
  mutate(
    New_demo_trait = case_when(
      Demo_trait == "T_generation"    ~ "Generation time",
      Demo_trait == "R0"              ~ "Net reproductive rate",
      Demo_trait == "percapita_repro" ~ "Per capita reproduction",
      Demo_trait == "Lmean"           ~ "Mean life expectancy",
      Demo_trait == "age_repro"       ~ "Age at first reproduction",
      Demo_trait == "Lambda"          ~ "Population growth rate",
      Demo_trait == "Reactivity"      ~ "Reactivity",
      Demo_trait == "FirstStepAtt"    ~ "First step attenuation",
      Demo_trait == "Lmax"            ~ "Maximum length of life",
      TRUE                            ~ Demo_trait
    ),
    hypoth = case_when(
      New_demo_trait == "Net reproductive rate"     & experiment == "DIST" ~ "H1",
      New_demo_trait == "Per capita reproduction"   & experiment == "DIST" ~ "H1",
      New_demo_trait == "Generation time"           & experiment == "DIST" ~ "H2",
      New_demo_trait == "Age at first reproduction" & experiment == "DIST" ~ "H2",
      New_demo_trait == "Reactivity"                & experiment == "DIST" ~ "H3",
      New_demo_trait == "First step attenuation"    & experiment == "DIST" ~ "H3",
      New_demo_trait == "Population growth rate"    & experiment == "DIST" ~ "H4",
      New_demo_trait == "Maximum length of life"    & experiment == "NPK"  ~ "H5",
      New_demo_trait == "Mean life expectancy"      & experiment == "NPK"  ~ "H5",
      TRUE ~ "none"
    )
  ) %>%
  filter(hypoth != "none") %>%
  mutate(
    sign = if_else(estimate >= 0, "pos", "neg"),
    overlaps0 = (estimate - std.error) <= 0 & (estimate + std.error) >= 0,
    sig = if_else(overlaps0, "ns", "sig"),
    col_cat = factor(paste(sign, sig, sep = "_"),
                     levels = c("neg_ns","neg_sig","pos_ns","pos_sig")))


###############################################################################
# 2) FUNCTION TO BUILD FIGURE
################################################################################

# Specify colour palette and use ggplot to make the plots
make_h_plot <- function(dat, Hx, show_y = TRUE) {
  cols <- c(
    "neg_ns" = "#9ecae1",  # pastel blue
    "neg_sig"= "#08519c",  # bold blue
    "pos_ns" = "#fcae91",  # pastel red
    "pos_sig"= "#cb181d"   # bold red
  )
  
  dat %>%
    filter(hypoth == Hx) %>%
    ggplot(aes(x = estimate, y = time_period, color = col_cat)) +
    theme_bw() +
    geom_point(size = 2) +
    # ggplot2 >= 3.4; for older ggplot2, use ggstance::geom_errorbarh(...)
    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
    facet_wrap(~ New_demo_trait, ncol = 1, scales = "free_x") +  # free x per facet
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = cols, guide = "none") +
    labs(title = Hx, x = NULL, y = NULL) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold", size = 11),
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      # Keep x-axis ticks/labels visible for ALL panels
      axis.text.x = element_text(),
      axis.ticks.x = element_line()
    ) +
    { if (!show_y) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) else theme() }
}

###############################################################################
# 3) FUNCTION TO CREATE MULTIPLE PLOT LAYOUT
################################################################################

build_rm_figure_2x3 <- function(RM_df,
                                figure_title = "RM Effects: Ordbeta (DIST/NPK)") {
  
  H_avail <- unique(dat$hypoth)
  get_or_spacer <- function(Hx, show_y) {
    if (Hx %in% H_avail) make_h_plot(dat, Hx, show_y = show_y) else plot_spacer()
  }
  
  # Top row
  p_A <- get_or_spacer("H1", show_y = TRUE)    # row1 col1
  p_B <- get_or_spacer("H2", show_y = FALSE)   # row1 col2
  p_C <- get_or_spacer("H3", show_y = FALSE)   # row1 col3
  
  # Bottom row
  p_D <- get_or_spacer("H4", show_y = TRUE)    # row2 col1
  p_E <- get_or_spacer("H5", show_y = TRUE)    # row2 col2 <-- keep y axis visible here
  p_F <- plot_spacer()                         # row2 col3 (blank)
  
  design <- "
  ABC
  DEF
  "
  
  (p_A + p_B + p_C + p_D + p_E + p_F) +
    plot_layout(design = design, guides = "collect") +
    plot_annotation(
      title = figure_title,
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
}

################################################################################
# 3) SAVE FIGURE
################################################################################

# make the figure the size and layout you want
fig_RM3 <- build_rm_figure_2x3(dat)

# Print it
fig_RM3

# Save as TIFF (good for journals)
ggsave("figures/fig_RM3_effects.tiff", plot = fig_RM3,
       width = 12, height = 8, dpi = 600, compression = "lzw")

# final plot without a title for the manuscript
# Save as TIFF (good for journals)
ggsave("figures/fig_3_effects_final.tiff", plot = fig_RM3,
       width = 8, height = 6, dpi = 600, compression = "lzw")
