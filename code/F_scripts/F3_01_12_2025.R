# R Thornley
# 09/09/2025
# Project: P1_COMPADRE_DRAGNET
# Script: F3_01_12_2025

################################################################################
# Instructions
################################################################################

# 1) LOAD and TIDY uni-variate robust model results data / tidy names
# 2) BUILD FIGURE
# 3) CREATE PLOT LAYOUT
# 4) SAVE FIGURE

###############################################################################
# 1) LOAD and TIDY uni-variate robust model results data / tidy names
###############################################################################

# import the model results 
mods <- readRDS(file = "results/robust_model_results.rds")
names(mods)
# extract the key parameters of the results for a table
estimates <- mods %>%
  filter(model == "Ordbeta") %>% 
  dplyr::select(model, time_period, experiment, Demo_trait, results, r2) %>%
  unnest(c(results)) %>%
  filter(term == "Demo_value")

CIs <- mods %>%
  filter(model == "Ordbeta") %>% 
  dplyr::select(model, time_period, experiment, Demo_trait, confint) %>%
  unnest(c(confint)) %>%
  filter(term == "Demo_value")

both <- estimates %>%
  left_join(CIs) %>%
  rename("CI-2.5%" = '2.5 %', "CI-97.5%" = '97.5 %')

names(both)

# create the data set for plotting
dat <- 
  both %>% mutate(
    New_demo_trait = case_when(
      Demo_trait == "T_generation"    ~ "Generation time",
      Demo_trait == "R0"              ~ "Net reproductive rate",
      Demo_trait == "percapita_repro" ~ "Per capita reproduction",
      Demo_trait == "Lmean"           ~ "Mean life expectancy",
      Demo_trait == "age_repro"       ~ "Age at first reproduction",
      Demo_trait == "Lambda"          ~ "Population growth rate",
      Demo_trait == "Reactivity"      ~ "Reactivity",
      Demo_trait == "FirstStepAtt"    ~ "First step attenuation",
      Demo_trait == "Lmax"            ~ "Maximum longevity", # changed this name
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
      New_demo_trait == "Maximum longevity"    & experiment == "NPK"  ~ "H5", # changed this name
      New_demo_trait == "Mean life expectancy"      & experiment == "NPK"  ~ "H5",
      TRUE ~ "none"
    )) %>%
  filter(hypoth != "none") %>%
  mutate(
    sign = if_else(estimate >= 0, "pos", "neg"),
    overlaps0 = (estimate - std.error) <= 0 & (estimate + std.error) >= 0,
    sig = if_else(overlaps0, "ns", "sig"),
    col_cat = factor(paste(sign, sig, sep = "_"),
                     levels = c("neg_ns","neg_sig","pos_ns","pos_sig"))) %>%
  mutate(
    
  )

###############################################################################
# 2) FUNCTION TO BUILD FIGURE
################################################################################

make_h_plot <- function(dat, Hx, show_y = TRUE,
                        trait_filter = NULL,
                        title_override = NULL) {
  cols <- c(
    "neg_ns" = "#9ecae1",  # pastel blue
    "neg_sig"= "#08519c",  # bold blue
    "pos_ns" = "#fcae91",  # pastel red
    "pos_sig"= "#cb181d"   # bold red
  )
  
  plot_dat <- dat %>%
    filter(hypoth == Hx)
  
  # Optionally restrict to one (or more) traits
  if (!is.null(trait_filter)) {
    plot_dat <- plot_dat %>%
      filter(New_demo_trait %in% trait_filter)
  }
  
  panel_title <- if (!is.null(title_override)) title_override else Hx
  
  plot_dat %>%
    ggplot(aes(x = estimate, y = time_period, color = col_cat)) +
    theme_bw() +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = estimate - std.error,
                       xmax = estimate + std.error),
                   height = 0.2) +
    facet_wrap(~ New_demo_trait, ncol = 1, scales = "free_x") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = cols, guide = "none") +
    labs(title = panel_title, x = NULL, y = NULL) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold", size = 11),
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(),
      axis.ticks.x = element_line()
    ) +
    { if (!show_y) theme(axis.text.y = element_blank(),
                         axis.ticks.y = element_blank()) else theme() }
}

################################################################################
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
  
  # Split H5 into two panels: H5a (Maximum longevity) and H5b (Mean life expectancy)
  p_E <- make_h_plot(dat,
                     Hx             = "H5",
                     show_y         = TRUE,
                     trait_filter   = "Maximum longevity",
                     title_override = "H5")
  
  p_F <- make_h_plot(dat,
                     Hx             = "H5",
                     show_y         = FALSE,
                     trait_filter   = "Mean life expectancy",
                     title_override = "H5")
  
  design <- "
  ABC
  DEF
  "
  
  (p_A + p_B + p_C + p_D + p_E + p_F) +
    plot_layout(design = design, guides = "collect")
  
}


################################################################################
# 4) SAVE FIGURE
################################################################################

# make the figure the size and layout you want
fig_RM3 <- build_rm_figure_2x3(dat)

# Print it
fig_RM3

# final plot without a title for the manuscript
# Save as TIFF (good for journals)
ggsave("figures/Fig_3_01_12_2025.tiff", plot = fig_RM3,
       width = 8, height = 6, dpi = 600, compression = "lzw")
