# R Thornley
# 08/09/2025
# Map of DRAGNet and COMPADRE sites

# at the moment this code produces an ok figure
# I want some key changes to be made
# The amount of data keys for both the maps - how are they decided?
# For compadre - I would like to start with the smallest dot representing 1, 
# the medium sized dot representing the median value and the largest dot representing the 
# sites wiht the most data - at the moment the first two are 4.
# the amout of data key on the DRAGNet map looks more logical.
# I would also like to re-position the three panels.
# to fit a full width half page on a jounral; 
# the top row would be the two maps side by side
# the second row would be the hypothesis figure centered and larger to cover as much as the
# full width as possible


# Load Packages
# install.packages(c("sf","ggplot2","rnaturalearth","rnaturalearthdata","biscale","sfdep","classInt","patchwork","ggrepel","ggspatial"))
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
install.packages("biscale")
library(biscale)
library(sfdep)
library(classInt)
library(patchwork)
library(ggrepel)
library(ggspatial)
library(tidyverse)
install.packages("ggnewscale")
library(ggnewscale)

# we need the location of each of the matrices used in the whole analyses
# we also need the location of each of the sites used in the analyses

################################################################################
# GET THE DRAGNet sites and locations
################################################################################

Drag <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv")
names(Drag)
unique(Drag$group_var)

# read in data; filter by the RE of taxon site, split the name into species and site
Drag <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv") %>%
  filter(group_var == "taxon_site") %>%
  filter(model == "Ordbeta") %>%
  separate(group, into = c("sp1", "sp2", "site_name"), sep = "_") %>%
  mutate(Taxon = paste0(sp1,"_", sp2)) %>%
  dplyr::select(Taxon, site_name, model, time_period, experiment)

# filter by the taxa we are studying 
mod_taxa <- readRDS("results/List_taxa_OLS_mods.R")
Drag <- Drag %>% filter(Taxon %in% mod_taxa)

# now filter by our hypotheses - drop INTER
Drag <- Drag %>% filter(experiment %in% c("DIST", "NPK"))

# this gives us the number of species across
Drag <- Drag %>% 
  count(site_name, name = "number_species_occ_in_site") %>%
  arrange(-number_species_occ_in_site)

# get list of sites we have used
list_sites <- unique(Drag$site_name)

# read in the DRAGNet meta data for the sites - with the site locations
meta <- read_csv("data/site-info-drag-2024-07-29.csv")
names(meta)
meta <- meta %>% filter(site_name %in% list_sites) %>% 
  dplyr::select(site_name, site_code, latitude, longitude)

# join the Drag df with the meta df
cover_locations <- Drag %>% left_join(meta)

################################################################################
# GET THE COMPADRE sites and locations
################################################################################

# get the list of species that we use in the final analysis
taxa_wanted <- readRDS("results/List_taxa_OLS_mods.R")

# read in the compadre meta data
# and select just the matrix ID and the location columns that we need
matrix_locations <- read_csv("results/meta_data_compadre.csv") %>%
  dplyr::select(MatrixID, SpeciesAccepted, Lat, Lon, Country, Continent) %>%
  rename(Taxon = SpeciesAccepted) %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  filter(Taxon %in% taxa_wanted)

# this data frame can be simplified further
# we can also create a column with 
matrix_locations <- 
  matrix_locations %>%
  dplyr::select(Taxon, Lat, Lon) %>%
  drop_na() %>%
  count(Taxon, Lat, Lon, name = "number_matrices")

################################################################################
# Rough plots of the global locations of the data sets
################################################################################

world <- ne_countries(scale = "medium", returnclass = "sf")
theme_set(theme_classic())

comp_fig <- 
  ggplot(data = world) +
  geom_sf()+
  geom_point(data = matrix_locations, aes(x = Lon, y = Lat, size = number_matrices), 
             colour = "green")
comp_fig

names(cover_locations)
drag_fig <- 
  ggplot(data = world) +
  geom_sf()+
  geom_point(data = cover_locations, aes(x = longitude, y = latitude, size = number_species_occ_in_site), 
             colour = "orange")
drag_fig

################################################################################
# Separate maps (no overlap) with hypothesis figure
###############################################################################
get_breaks <- function(z) {
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  br <- unique(as.numeric(quantile(z, c(0.5, 0.85, 0.97), na.rm = TRUE)))
  br <- br[br > 0 & is.finite(br)]
  if (length(br) < 2) { br <- sort(unique(pretty(z, n = 3))); br <- br[br > 0] }
  if (length(br) < 2) return(NULL)
  br
}

# Map function (one-line title, compact legend on top, taller panel, no axes) 
make_point_map <- function(df, lon, lat, size_col, fill_color, title_text,
                           max_pt_mm = 3.8, base_pt = 8.5, title_pt = 8) {
  
  pts  <- st_as_sf(df, coords = c(lon, lat), crs = 4326) |> st_transform(crs_robin)
  brks <- get_breaks(pts[[size_col]])
  lab_fmt <- label_number(accuracy = 1, big.mark = ",")
  
  p <- ggplot() +
    geom_sf(data = world_p, fill = "grey96", color = "grey80", linewidth = 0.2) +
    geom_sf(
      data = pts, aes(size = .data[[size_col]]),
      shape = 21, fill = fill_color, color = "black",
      stroke = 0.20, alpha = 0.85
    ) +
    coord_sf(
      crs = crs_robin,
      xlim = c(crop_bbox["xmin"], crop_bbox["xmax"]),
      ylim = c(crop_bbox["ymin"], crop_bbox["ymax"]),
      expand = FALSE
    ) +
    { if (is.null(brks))
      scale_size_area(name = "Amount of data", max_size = max_pt_mm)
      else
        scale_size_area(name = "Amount of data", max_size = max_pt_mm,
                        breaks = brks, labels = lab_fmt)
    } +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_classic(base_size = base_pt) +
    theme(
      plot.title = element_text(size = title_pt, face = "bold",
                                lineheight = 1.0, margin = margin(b = 1)),
      axis.text  = element_blank(),                        # remove axes for tighter/cleaner look
      axis.ticks = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size = base_pt),
      legend.text  = element_text(size = base_pt - 0.5),
      legend.key.width  = unit(6, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.spacing.x  = unit(1, "mm"),
      legend.spacing.y  = unit(0.5, "mm"),
      legend.margin     = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0),
      plot.margin       = margin(1, 2, 1, 2),              # small & symmetric
      aspect.ratio      = 0.55                             # make the map panel taller
    )
  p
}

# Panels A & B (one-line titles, smaller font) 
p_matrices <- make_point_map(
  matrix_locations, "Lon", "Lat", "number_matrices",
  col_green,  "A) COMPADRE — points sized by number of matrices"
)

p_cover <- make_point_map(
  cover_locations, "longitude", "latitude", "number_species_occ_in_site",
  col_orange, "B) DRAGNet — points sized by number of species records"
)

# Panel C: TIFF graphic, trimmed, aligned to left with same inner margins
graphic_path <- "figures/Hypothesis_figure.tiff"
img <- image_read(graphic_path) |> image_trim()

# give panel C a left/right margin similar to A/B so the left edges align visually
panel_left  <- unit(2, "pt")
panel_right <- unit(2, "pt")

p_graphic <- ggdraw() +
  draw_image(img, x = 0, y = 1, width = 1, height = 1, hjust = 0, vjust = 1, interpolate = TRUE) +
  draw_label("C) Overview of hypotheses and study variables",
             x = 0.01, y = 0.985, hjust = 0, vjust = 1,
             fontface = "bold", size = 8) +
  theme(plot.margin = margin(1, convertUnit(panel_right, "pt", valueOnly=TRUE),
                             1, convertUnit(panel_left,  "pt", valueOnly=TRUE)))

# Stack (tight gaps)
p_stack3 <- (
  p_matrices / p_cover / p_graphic
) + plot_layout(ncol = 1, heights = c(1, 1, 1)) &
  theme(plot.margin = margin(1, 1, 1, 1))
p_stack3
# Export
ragg::agg_tiff("figures/map_hypoth_sept_final.tiff",
               width = 90, height = 190, units = "mm", res = 600, compression = "lzw")
print(p_stack3); dev.off()

################################################################################
# Separate maps (no overlap) with hypothesis figure - example 2
###############################################################################

# ── Extra packages used later ──────────────────────────────────────────────────
library(magick)
library(cowplot)
library(ragg)
library(patchwork)

# ── Figure controls (easy to tweak) ────────────────────────────────────────────
FIG_WIDTH_MM  <- 180          # target full-width half-page (adjust per journal)
FIG_HEIGHT_MM <- 130          # overall height; tweak after you preview
TOP_ROW_REL   <- 0.9          # relative height of the top row (maps)
BOT_ROW_REL   <- 1.1          # relative height of the bottom row (hypothesis)
MAX_DOT_MM    <- 3.8          # consistent max dot size across maps

col_green  <- "#37a354"
col_orange <- "#f39c12"

world <- ne_countries(scale = "medium", returnclass = "sf")
theme_set(theme_classic())

# If not already defined earlier in your script:
crs_robin <- "+proj=robin"
world_p   <- sf::st_transform(world, crs_robin)

# Crop box (fallback if not defined upstream) — tweak if you want tighter framing
if (!exists("crop_bbox")) {
  wb <- st_bbox(world_p)
  crop_bbox <- c(xmin = wb["xmin"], xmax = wb["xmax"], ymin = wb["ymin"], ymax = wb["ymax"])
}

# ── Helpers for breaks ─────────────────────────────────────────────────────────
get_quantile_breaks <- function(z) {
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  br <- unique(as.numeric(stats::quantile(z, c(0.5, 0.85, 0.97), na.rm = TRUE)))
  br <- br[br > 0 & is.finite(br)]
  if (length(br) < 2) {
    br <- sort(unique(pretty(z, n = 3))); br <- br[br > 0]
  }
  if (length(br) < 2) return(NULL)
  br
}

get_compadre_breaks <- function(z) {
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  med <- as.integer(round(stats::median(z)))
  mx  <- max(z)
  br  <- unique(sort(c(1L, med, mx)))
  # guard against duplicates (e.g., if median == 1 or 1 == max in tiny datasets)
  if (length(br) < 3) {
    add <- setdiff(sort(unique(pretty(z, n = 3))), br)
    br  <- sort(unique(c(br, head(add, 3 - length(br)))))
  }
  br
}

# ── Generic map maker with optional explicit breaks ────────────────────────────
make_point_map <- function(df, lon, lat, size_col, fill_color, title_text,
                           legend_title = "Amount of data",
                           breaks = NULL,
                           max_pt_mm = MAX_DOT_MM, base_pt = 8.5, title_pt = 8) {
  
  pts  <- st_as_sf(df, coords = c(lon, lat), crs = 4326) |> st_transform(crs_robin)
  lab_fmt <- scales::label_number(accuracy = 1, big.mark = ",")
  
  ggplot() +
    geom_sf(data = world_p, fill = "grey96", color = "grey80", linewidth = 0.2) +
    geom_sf(
      data = pts, aes(size = .data[[size_col]]),
      shape = 21, fill = fill_color, color = "black",
      stroke = 0.20, alpha = 0.85
    ) +
    coord_sf(
      crs = crs_robin,
      xlim = c(crop_bbox["xmin"], crop_bbox["xmax"]),
      ylim = c(crop_bbox["ymin"], crop_bbox["ymax"]),
      expand = FALSE
    ) +
    {
      if (is.null(breaks)) {
        scale_size_area(name = legend_title, max_size = max_pt_mm)
      } else {
        scale_size_area(name = legend_title, max_size = max_pt_mm,
                        breaks = breaks, labels = lab_fmt)
      }
    } +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_classic(base_size = base_pt) +
    theme(
      plot.title = element_text(size = title_pt, face = "bold",
                                lineheight = 1.0, margin = margin(b = 1)),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size = base_pt),
      legend.text  = element_text(size = base_pt - 0.5),
      legend.key.width  = unit(6, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.spacing.x  = unit(1, "mm"),
      legend.spacing.y  = unit(0.5, "mm"),
      legend.margin     = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0),
      plot.margin       = margin(1, 2, 1, 2),
      aspect.ratio      = 0.55
    )
}

# ── Build the two maps with the new key logic ──────────────────────────────────

# COMPADRE: breaks = {1, median, max} of number_matrices
comp_breaks <- get_compadre_breaks(matrix_locations$number_matrices)

p_matrices <- make_point_map(
  matrix_locations, "Lon", "Lat", "number_matrices",
  col_green,
  "A) COMPADRE — points sized by number of matrices",
  legend_title = "Matrices per site",
  breaks = comp_breaks
)
p_matrices

# DRAGNet: keep quantile-style legend
drag_breaks <- get_quantile_breaks(cover_locations$number_species_occ_in_site)

p_cover <- make_point_map(
  cover_locations, "longitude", "latitude", "number_species_occ_in_site",
  col_orange,
  "B) DRAGNet — points sized by number of species records",
  legend_title = "Species records per site",
  breaks = drag_breaks
)
p_cover

# ── Panel C: hypothesis graphic, trimmed and labeled ───────────────────────────
graphic_path <- "figures/Hypothesis_figure.tiff"
img <- image_read(graphic_path) |> image_trim()

p_graphic <- ggdraw() +
  draw_image(img, x = 0, y = 1, width = 1, height = 1, hjust = 0, vjust = 1, interpolate = TRUE) +
  draw_label("C) Overview of hypotheses and study variables",
             x = 0.01, y = 0.985, hjust = 0, vjust = 1,
             fontface = "bold", size = 9) +
  theme(plot.margin = margin(2, 2, 2, 2))
p_graphic

# ── Patchwork: maps (top row) + hypothesis (full-width bottom row) ─────────────
design <- "
AB
CC
"

p_full <- (p_matrices + p_cover + p_graphic) +
  plot_layout(design = design,
              heights = c(TOP_ROW_REL, BOT_ROW_REL))

# Preview
p_full

# ── Export at journal-friendly size ────────────────────────────────────────────
ragg::agg_tiff("figures/map_hypoth_halfpage_fullwidth.tiff",
               width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM,
               units = "mm", res = 600, compression = "lzw")
print(p_full); dev.off()

################################################################################
# 
################################################################################

# ── Packages ───────────────────────────────────────────────────────────────────
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggspatial)
library(magick)
library(cowplot)
library(ragg)

# ── Data prep (your existing wrangling) ────────────────────────────────────────
# DRAGNet
Drag <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv") %>%
  filter(group_var == "taxon_site") %>%
  filter(model == "Ordbeta") %>%
  separate(group, into = c("sp1", "sp2", "site_name"), sep = "_") %>%
  mutate(Taxon = paste0(sp1,"_", sp2)) %>%
  dplyr::select(Taxon, site_name, model, time_period, experiment)

mod_taxa <- readRDS("results/List_taxa_OLS_mods.R")
Drag <- Drag %>% filter(Taxon %in% mod_taxa)
Drag <- Drag %>% filter(experiment %in% c("DIST", "NPK"))

Drag <- Drag %>% 
  count(site_name, name = "number_species_occ_in_site") %>%
  arrange(-number_species_occ_in_site)

list_sites <- unique(Drag$site_name)

meta <- read_csv("data/site-info-drag-2024-07-29.csv") %>%
  filter(site_name %in% list_sites) %>% 
  dplyr::select(site_name, site_code, latitude, longitude)

cover_locations <- Drag %>% left_join(meta, by = "site_name")

# COMPADRE
taxa_wanted <- readRDS("results/List_taxa_OLS_mods.R")

matrix_locations <- read_csv("results/meta_data_compadre.csv") %>%
  dplyr::select(MatrixID, SpeciesAccepted, Lat, Lon, Country, Continent) %>%
  rename(Taxon = SpeciesAccepted) %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  filter(Taxon %in% taxa_wanted) %>%
  dplyr::select(Taxon, Lat, Lon) %>%
  drop_na() %>%
  count(Taxon, Lat, Lon, name = "number_matrices")

# Basemap + projection #
# ── Basemap + projection (with tighter crop) ───────────────────────────────────
world     <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
crs_robin <- "+proj=robin"
world_p   <- sf::st_transform(world, crs_robin)

# Build a projected bbox by transforming a WGS84 lat–lon window to Robin
make_robin_bbox <- function(lat_min = -60, lat_max = 80,
                            lon_min = -180, lon_max = 180,
                            crs_target = crs_robin) {
  poly_ll <- sf::st_polygon(list(matrix(
    c(lon_min, lat_min,
      lon_max, lat_min,
      lon_max, lat_max,
      lon_min, lat_max,
      lon_min, lat_min),
    ncol = 2, byrow = TRUE))) |>
    sf::st_sfc(crs = 4326)
  bb_proj <- sf::st_transform(poly_ll, crs_target) |> sf::st_bbox()
  c(xmin = bb_proj["xmin"], xmax = bb_proj["xmax"],
    ymin = bb_proj["ymin"], ymax = bb_proj["ymax"])
}

# Set a sensible crop to remove empty polar space; tweak if you want tighter
crop_bbox <- make_robin_bbox(lat_min = -60, lat_max = 80)

# ── Hypothesis panel (add spacer above to avoid any overlap) ───────────────────
graphic_path <- "figures/Hypothesis_figure.tiff"
img <- magick::image_read(graphic_path) |> magick::image_trim()

p_graphic_inner <- cowplot::ggdraw() +
  cowplot::draw_image(img, x = 0, y = 1, width = 1, height = 1,
                      hjust = 0, vjust = 1, interpolate = TRUE) +
  cowplot::draw_label("C) Overview of hypotheses and study variables",
                      x = 0.01, y = 0.96,   # lower than before
                      hjust = 0, vjust = 1,
                      fontface = "bold", size = 9) +
  theme(plot.margin = margin(6, 2, 2, 2))  # modest top padding

# Wrap with a spacer row that pushes the whole content down slightly
p_graphic <- patchwork::wrap_plots(
  patchwork::plot_spacer(),   # top spacer
  p_graphic_inner,
  ncol = 1,
  heights = c(0.10, 0.90)     # increase first value to push content further down
)

# ── Figure controls ────────────────────────────────────────────────────────────
FIG_WIDTH_MM  <- 180   # target full-width (adjust per journal)
FIG_HEIGHT_MM <- 130   # overall height (adjust after preview)
TOP_ROW_REL   <- 0.9   # relative height of top row (maps)
BOT_ROW_REL   <- 1.1   # relative height of bottom row (hypothesis)
MAX_DOT_MM    <- 3.8   # max symbol size (consistent across maps)

col_green  <- "#37a354"
col_orange <- "#f39c12"

# ── Legend break helpers ───────────────────────────────────────────────────────
get_quantile_breaks <- function(z) {
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  br <- unique(as.numeric(stats::quantile(z, c(0.5, 0.85, 0.97), na.rm = TRUE)))
  br <- br[br > 0 & is.finite(br)]
  if (length(br) < 2) {
    br <- sort(unique(pretty(z, n = 3))); br <- br[br > 0]
  }
  if (length(br) < 2) return(NULL)
  br
}

get_compadre_breaks <- function(z) {
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  med <- as.integer(round(stats::median(z)))
  mx  <- max(z)
  br  <- unique(sort(c(1L, med, mx)))
  if (length(br) < 3) {
    add <- setdiff(sort(unique(pretty(z, n = 3))), br)
    br  <- sort(unique(c(br, head(add, 3 - length(br)))))
  }
  br
}

# ── Generic point map with 4-side border and no axis lines ────────────────────
make_point_map <- function(df, lon, lat, size_col, fill_color, title_text,
                           legend_title = "Amount of data",
                           breaks = NULL,
                           max_pt_mm = MAX_DOT_MM, base_pt = 8.5, title_pt = 8) {
  
  pts  <- st_as_sf(df, coords = c(lon, lat), crs = 4326) |> st_transform(crs_robin)
  lab_fmt <- scales::label_number(accuracy = 1, big.mark = ",")
  
  ggplot() +
    geom_sf(data = world_p, fill = "grey96", color = "grey80", linewidth = 0.2) +
    geom_sf(
      data = pts, aes(size = .data[[size_col]]),
      shape = 21, fill = fill_color, color = "black",
      stroke = 0.20, alpha = 0.85
    ) +
    coord_sf(
      crs = crs_robin,
      xlim = c(crop_bbox["xmin"], crop_bbox["xmax"]),
      ylim = c(crop_bbox["ymin"], crop_bbox["ymax"]),
      expand = FALSE
    ) +
    {
      if (is.null(breaks)) {
        scale_size_area(name = legend_title, max_size = max_pt_mm)
      } else {
        scale_size_area(name = legend_title, max_size = max_pt_mm,
                        breaks = breaks, labels = lab_fmt)
      }
    } +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_classic(base_size = base_pt) +
    theme(
      plot.title       = element_text(size = title_pt, face = "bold",
                                      lineheight = 1.0, margin = margin(b = 1)),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),                                  # remove partial axis lines
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.3),  # add 4-side frame
      legend.position  = "top",
      legend.direction = "horizontal",
      legend.title     = element_text(size = base_pt),
      legend.text      = element_text(size = base_pt - 0.5),
      legend.key.width  = unit(6, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.spacing.x  = unit(1, "mm"),
      legend.spacing.y  = unit(0.5, "mm"),
      legend.margin     = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0),
      plot.margin       = margin(1, 2, 1, 2),
      aspect.ratio      = 0.55
    )
}

# ── Build maps with required legend logic ──────────────────────────────────────
comp_breaks <- get_compadre_breaks(matrix_locations$number_matrices)
p_matrices <- make_point_map(
  matrix_locations, "Lon", "Lat", "number_matrices",
  col_green,
  "A) COMPADRE — points sized by number of matrices",
  legend_title = "Matrices per site",
  breaks = comp_breaks
)

drag_breaks <- get_quantile_breaks(cover_locations$number_species_occ_in_site)
p_cover <- make_point_map(
  cover_locations, "longitude", "latitude", "number_species_occ_in_site",
  col_orange,
  "B) DRAGNet — points sized by number of species records",
  legend_title = "Species records per site",
  breaks = drag_breaks
)

# ── Hypothesis panel (nudged down to avoid title overlap) ─────────────────────
graphic_path <- "figures/Hypothesis_figure.tiff"
img <- image_read(graphic_path) |> image_trim()

p_graphic <- ggdraw() +
  draw_image(img, x = 0, y = 1, width = 1, height = 1,
             hjust = 0, vjust = 1, interpolate = TRUE) +
  draw_label("C) Overview of hypotheses and study variables",
             x = 0.01, y = 0.94,  # moved down from 0.985 → 0.94
             hjust = 0, vjust = 1,
             fontface = "bold", size = 9) +
  theme(plot.margin = margin(10, 2, 2, 2))  # extra top padding

# ── Layout: AB on top, CC full-width bottom ───────────────────────────────────
design <- "
AB
CC
"
p_full <- (p_matrices + p_cover + p_graphic) +
  plot_layout(design = design,
              heights = c(TOP_ROW_REL, BOT_ROW_REL))
p_full

# ── Export ─────────────────────────────────────────────────────────────────────
ragg::agg_tiff("figures/map_hypoth_halfpage_fullwidth.tiff",
               width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM,
               units = "mm", res = 600, compression = "lzw")
print(p_full); dev.off()

################################################################################
# INDEPENDENT TITLES
################################################################################

# ── Packages ───────────────────────────────────────────────────────────────────
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggspatial)
library(magick)
library(cowplot)
library(ragg)

# ── Data prep (your existing wrangling) ────────────────────────────────────────
# DRAGNet
Drag <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv") %>%
  filter(group_var == "taxon_site") %>%
  filter(model == "Ordbeta") %>%
  separate(group, into = c("sp1", "sp2", "site_name"), sep = "_") %>%
  mutate(Taxon = paste0(sp1,"_", sp2)) %>%
  dplyr::select(Taxon, site_name, model, time_period, experiment)

mod_taxa <- readRDS("results/List_taxa_OLS_mods.R")
Drag <- Drag %>% filter(Taxon %in% mod_taxa)
Drag <- Drag %>% filter(experiment %in% c("DIST", "NPK"))

Drag <- Drag %>% 
  count(site_name, name = "number_species_occ_in_site") %>%
  arrange(-number_species_occ_in_site)

list_sites <- unique(Drag$site_name)

meta <- read_csv("data/site-info-drag-2024-07-29.csv") %>%
  filter(site_name %in% list_sites) %>% 
  dplyr::select(site_name, site_code, latitude, longitude)

cover_locations <- Drag %>% left_join(meta, by = "site_name")

# COMPADRE
taxa_wanted <- readRDS("results/List_taxa_OLS_mods.R")

matrix_locations <- read_csv("results/meta_data_compadre.csv") %>%
  dplyr::select(MatrixID, SpeciesAccepted, Lat, Lon, Country, Continent) %>%
  rename(Taxon = SpeciesAccepted) %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  filter(Taxon %in% taxa_wanted) %>%
  dplyr::select(Taxon, Lat, Lon) %>%
  drop_na() %>%
  count(Taxon, Lat, Lon, name = "number_matrices")

# ── Basemap + projection (with tighter crop) ───────────────────────────────────
world     <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
crs_robin <- "+proj=robin"
world_p   <- sf::st_transform(world, crs_robin)

# Build a projected bbox by transforming a WGS84 lat–lon window to Robin
make_robin_bbox <- function(lat_min = -60, lat_max = 80,
                            lon_min = -180, lon_max = 180,
                            crs_target = crs_robin) {
  poly_ll <- sf::st_polygon(list(matrix(
    c(lon_min, lat_min,
      lon_max, lat_min,
      lon_max, lat_max,
      lon_min, lat_max,
      lon_min, lat_min),
    ncol = 2, byrow = TRUE))) |>
    sf::st_sfc(crs = 4326)
  bb_proj <- sf::st_transform(poly_ll, crs_target) |> sf::st_bbox()
  c(xmin = bb_proj["xmin"], xmax = bb_proj["xmax"],
    ymin = bb_proj["ymin"], ymax = bb_proj["ymax"])
}

# Crop to remove empty polar space (tweak lat range if desired)
crop_bbox <- make_robin_bbox(lat_min = -60, lat_max = 80)

# ── Figure controls ────────────────────────────────────────────────────────────
FIG_WIDTH_MM  <- 180   # target full-width (adjust per journal)
FIG_HEIGHT_MM <- 130   # overall height (adjust after preview)
TOP_ROW_REL   <- 0.9   # relative height of top row (maps)
BOT_ROW_REL   <- 1.1   # relative height of bottom row (hypothesis)
MAX_DOT_MM    <- 3.8   # max symbol size (consistent across maps)

col_green  <- "#37a354"
col_orange <- "#f39c12"

# ── Legend break helpers ───────────────────────────────────────────────────────
get_quantile_breaks <- function(z) {
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  br <- unique(as.numeric(stats::quantile(z, c(0.5, 0.85, 0.97), na.rm = TRUE)))
  br <- br[br > 0 & is.finite(br)]
  if (length(br) < 2) {
    br <- sort(unique(pretty(z, n = 3)))
    br <- br[br > 0]
  }
  if (length(br) < 2) return(NULL)
  br
}


get_compadre_breaks <- function(z) {
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  med <- as.integer(round(stats::median(z)))
  mx  <- max(z)
  br  <- unique(sort(c(1L, med, mx)))
  if (length(br) < 3) {
    add <- setdiff(sort(unique(pretty(z, n = 3))), br)
    br  <- sort(unique(c(br, head(add, 3 - length(br)))))
  }
  br
}

# ── Generic point map with 4-side border and no axis lines ────────────────────
make_point_map <- function(df, lon, lat, size_col, fill_color, title_text,
                           legend_title = "Amount of data",
                           breaks = NULL,
                           max_pt_mm = MAX_DOT_MM, base_pt = 8.5, title_pt = 8) {
  
  pts  <- st_as_sf(df, coords = c(lon, lat), crs = 4326) |> st_transform(crs_robin)
  lab_fmt <- scales::label_number(accuracy = 1, big.mark = ",")
  
  ggplot() +
    geom_sf(data = world_p, fill = "grey96", color = "grey80", linewidth = 0.2) +
    geom_sf(
      data = pts, aes(size = .data[[size_col]]),
      shape = 21, fill = fill_color, color = "black",
      stroke = 0.20, alpha = 0.85
    ) +
    coord_sf(
      crs = crs_robin,
      xlim = c(crop_bbox["xmin"], crop_bbox["xmax"]),
      ylim = c(crop_bbox["ymin"], crop_bbox["ymax"]),
      expand = FALSE
    ) +
    {
      if (is.null(breaks)) {
        scale_size_area(name = legend_title, max_size = max_pt_mm)
      } else {
        scale_size_area(name = legend_title, max_size = max_pt_mm,
                        breaks = breaks, labels = lab_fmt)
      }
    } +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_classic(base_size = base_pt) +
    theme(
      plot.title       = element_text(size = title_pt, face = "bold",
                                      lineheight = 1.0, margin = margin(b = 1)),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),                                  # remove partial axis lines
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.3),  # add 4-side frame
      legend.position  = "top",
      legend.direction = "horizontal",
      legend.title     = element_text(size = base_pt),
      legend.text      = element_text(size = base_pt - 0.5),
      legend.key.width  = unit(6, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.spacing.x  = unit(1, "mm"),
      legend.spacing.y  = unit(0.5, "mm"),
      legend.margin     = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0),
      plot.margin       = margin(1, 2, 1, 2),
      aspect.ratio      = 0.55
    )
}

# ── Build maps with required legend logic ──────────────────────────────────────
comp_breaks <- get_compadre_breaks(matrix_locations$number_matrices)
p_matrices <- make_point_map(
  matrix_locations, "Lon", "Lat", "number_matrices",
  col_green,
  "A) COMPADRE — points sized by number of matrices",
  legend_title = "Matrices per site",
  breaks = comp_breaks
)

drag_breaks <- get_quantile_breaks(cover_locations$number_species_occ_in_site)
p_cover <- make_point_map(
  cover_locations, "longitude", "latitude", "number_species_occ_in_site",
  col_orange,
  "B) DRAGNet — points sized by number of species records",
  legend_title = "Species records per site",
  breaks = drag_breaks
)

# ── Panel C: title as a separate row (no overlap) ─────────────────────────────
graphic_path <- "figures/Hypothesis_figure.tiff"
img <- magick::image_read(graphic_path) |> magick::image_trim()

# Image-only panel
p_graphic_img <- cowplot::ggdraw() +
  cowplot::draw_image(img, x = 0, y = 1, width = 1, height = 1,
                      hjust = 0, vjust = 1, interpolate = TRUE) +
  theme(plot.margin = margin(2, 2, 2, 2))

# Title-only panel (appears above the image)
p_graphic_title <- ggplot() +
  theme_void() +
  labs(title = "C) Overview of hypotheses and study variables") +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 0,
                              face = "bold", size = 9,
                              margin = margin(b = 2))
  )

# Stack title above image inside Panel C
p_graphic <- patchwork::wrap_plots(
  p_graphic_title,
  p_graphic_img,
  ncol = 1,
  heights = c(0.08, 0.92)  # adjust to move title higher/lower
)

# ── Layout: AB on top, CC full-width bottom ───────────────────────────────────
design <- "
AB
CC
"
p_full <- (p_matrices + p_cover + p_graphic) +
  plot_layout(design = design,
              heights = c(TOP_ROW_REL, BOT_ROW_REL))
p_full

# ── Export ─────────────────────────────────────────────────────────────────────
ragg::agg_tiff("figures/map_hypoth_halfpage_fullwidth.tiff",
               width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM,
               units = "mm", res = 600, compression = "lzw")
print(p_full); dev.off()

