# R Thornley
# 09/09/2025
# Table of Variables for manuscript

library(tidyverse)
library(gt)

################################################################################
#  Build text data as a plain data.frame
################################################################################

vars <- data.frame(
  type = c(
    rep("Demographic: Life History Trait", 6),
    "Demographic: Emergent property",
    rep("Demographic: Transient Metric", 2),
    rep("Community Context", 6)
  ),
  variable = c(
    "Net Reproductive Rate",
    "Per Capita Reproduction",
    "Generation Time",
    "Age At First Reproduction",
    "Mean Life Expectancy",
    "Maximum Length of Life",
    "Population Growth Rate",
    "Reactivity",
    "First Step Attenuation",
    "Species richness",
    "Inverse Simpson’s diversity index",
    "Total cover",
    "Species provenance",
    "Functional diversity",
    "Functional specialisaton"
  ),
  definition = c(
    "The mean number of recruits produced during the mean life expectancy of an individual in the population.",
    "The number of offspring produced per individual.",
    "Number of years necessary for the individuals of a population to be fully replaced by new ones.",
    "Age at which an individual experiences the first reproductive event.",
    "The average life span of an individual in a population.",
    "The longest expected lifespan of an individual in a population.",
    "The asymptotic growth rate of the population.",
    "The maximum population growth achieved by a given population in a single timestep, relative to stable growth rate. A measure of the initial population recovery following a disturbance event when the population is still within a transient period.",
    "The population decline achieved by a given population structure in a single timestep, relative to stable growth rate. A measure of the initial population decline following a disturbance event when the population is still within a transient period.",
    "The number of species recorded in a DRAGNet plot.",
    "A measure of species diversity in a DRAGNet plot which incorporates both species richness and evenness.",
    "The cumulative (total) cover value of each of the surveyed plots.",
    "Whether a species is found in it's native range",
    "The level of functional strategy diversity of the DRAGNet site.",
    "The level of functional strategy specialisation of the DRAGNet site"
  ),
  source = c(
    rep("COMPADRE", 9),
    rep("DRAGNet", 4),
    rep("DRAGNet / TRY", 2)),
  stringsAsFactors = FALSE
)

################################################################################
# ORDER VARS
################################################################################
# Order and precompute even-row indices for manual striping
vars_ordered <- vars %>%
  mutate(type = factor(
    type,
    levels = c(
      "Demographic: Life History Trait",
      "Demographic: Emergent property",
      "Demographic: Transient Metric",
      "Community Context"
    )
  )) %>%
  arrange(type)

even_idx <- seq(2, nrow(vars_ordered), 2)  # 2,4,6,...

################################################################################
# MAKE TABLE BLACK AND WHITE: NO TAGS
################################################################################

tab_bw <- vars_ordered |>
  gt(groupname_col = "type") |>
  cols_label(
    variable   = "Variable name",
    definition = "Definition",
    source     = "Data source"
  ) |>
  tab_options(
    table.font.names          = c("Arial","Helvetica","sans-serif"),
    table.width               = pct(100),
    column_labels.font.weight = "bold",
    data_row.padding          = px(6),
    row_group.padding         = px(8)
  ) |>
  cols_width(
    variable   ~ px(230),
    definition ~ px(540),
    source     ~ px(110)
  ) |>
  # Strong structure
  opt_table_lines() |>
  tab_style(
    style = list(
      cell_text(weight = "bold", size = px(13)),
      cell_borders(sides = "bottom", color = "#666666", weight = px(2))
    ),
    locations = cells_row_groups()
  ) |>
  # Manual zebra striping (version-proof)
  tab_style(
    style = cell_fill(color = "#F2F2F2"),
    locations = cells_body(rows = even_idx)
  ) |>
  # Narrow left guide rule
  tab_style(
    style = cell_borders(sides = "left", color = "#888888", weight = px(3)),
    locations = cells_body(rows = TRUE)
  ) |>
  # Firm rule under column labels
  tab_style(
    style = cell_borders(sides = "bottom", color = "#444444", weight = px(2)),
    locations = cells_column_labels(everything())
  ) 
tab_bw



# export the table (size is to cover a half width of a column in a journal)
gtsave(tab_bw, "figures/T1_variables_table_bw_half_width.pdf", vwidth = 1200)

gtsave(
  data = tab_bw,
  "figures/T1_variables_table_bw_half_width.pdf",
  vwidth = 700)
gtsave(
  tab_bw,
  filename = "figures/T1_final.png",
  vwidth = 1600,   # controls width (pixels)
  vheight = 3000,  # increase if the table is long
  expand = 10      # add margin
)





