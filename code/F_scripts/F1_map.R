# R Thornley
# 08/09/2025
# Project: P1_COMPADRE_DRAGNET
# Script: F1_map
# Map of DRAGNet and COMPADRE sites

################################################################################
# INSTRUCTIONS
################################################################################
# 1) IMPORT DATA ON SPECIES LOCATIONS and NOS OF MATRICES
# 2) LOAD MAP BACKGROUNDS and GET CRS / Bounding boxes
# 3) CREATE CONTROLS ON FIGURE LAYOUT
# 4) PUT TOGETHER FINAL FIGURE 

################################################################################
# 1) IMPORT DATA ON SPECIES LOCATIONS and NOS OF MATRICES
################################################################################

# DRAGNet data load
Drag_dat <- read_csv("results/RE_SE_Taxon_site_all_DRAGNet.csv") %>%
  filter(group_var == "taxon_site") %>%
  filter(model == "Ordbeta") %>%
  separate(group, into = c("sp1", "sp2", "site_name"), sep = "_") %>%
  mutate(Taxon = paste0(sp1,"_", sp2)) %>%
  dplyr::select(Taxon, site_name, model, time_period, experiment)

# Load the list of taxa used in the modelling and filter the DRAGNet data for these taxa
mod_taxa <- readRDS("results/List_taxa_OLS_mods.R")
Drag <- 
  Drag_dat %>% filter(Taxon %in% mod_taxa) %>%
  filter(experiment %in% c("DIST", "NPK")) %>%
  count(site_name, name = "number_species_occ_in_site") %>%
  arrange(-number_species_occ_in_site)
# get a unique list of sites we are using
list_sites <- unique(Drag$site_name)

# load the DRAGNet meta data
meta <- read_csv("data/site-info-drag-2024-07-29.csv") %>%
  filter(site_name %in% list_sites) %>% 
  dplyr::select(site_name, site_code, latitude, longitude)

cover_locations <- Drag %>% left_join(meta, by = "site_name")

# COMPADRE data load and filter for the taxa we used in the final modelling
taxa_wanted <- readRDS("results/List_taxa_OLS_mods.R")

matrix_locations <- read_csv("results/meta_data_compadre.csv") %>%
  dplyr::select(MatrixID, SpeciesAccepted, Lat, Lon, Country, Continent) %>%
  rename(Taxon = SpeciesAccepted) %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  filter(Taxon %in% taxa_wanted) %>%
  dplyr::select(Taxon, Lat, Lon) %>%
  drop_na() %>%
  count(Taxon, Lat, Lon, name = "number_matrices")

################################################################################
# 2) LOAD MAP BACKGROUNDS and GET CRS / Bounding boxes
################################################################################

# Load the Basemap from rnaturalearth and specify the global projection
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

# Crop to remove empty polar space (tweak lat/long range if desired)
crop_bbox <- make_robin_bbox(lat_min = -60, lat_max = 80)

################################################################################
# 3) CREATE CONTROLS ON FIGURE LAYOUT
################################################################################

# ── Figure controls ────────────────────────────────────────────────────────────
FIG_WIDTH_MM  <- 180   # target full-width (adjust per journal)
FIG_HEIGHT_MM <- 130   # overall height (adjust after preview)
TOP_ROW_REL   <- 0.9   # relative height of top row (maps)
BOT_ROW_REL   <- 1.1   # relative height of bottom row (hypothesis)
MAX_DOT_MM    <- 3.8   # max symbol size (consistent across maps)

# colour palette specify
col_green  <- "#37a354"
col_orange <- "#f39c12"

################################################################################
# 4) HELPER FUNCTIONS 
################################################################################

# Breaks for the DRAGNet map categories
# function that specifies the number of 'breaks' required for the nos of sites included
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

# Breaks for the COMPADRE map categories
# function that specifies the number of 'breaks' required for the nos of matrices included
# quantile breaks
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

# Breaks specified manually
get_compadre_breaks_manual <- function(z, breaks = c(1, 10, 36)) {
  # keep only finite positive values
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  
  # ensure min and max coverage
  min_z <- min(z)
  max_z <- max(z)
  
  br <- sort(unique(c(min_z, breaks, max_z)))
  
  br
}

matrix_locations$number_matrices
cover_locations$number_species_occ_in_site

################################################################################
# 5a) MAP CREATION FUNCTION WITH BREAKS DISCRETELY SPECIFIED as quantile functions
################################################################################

# function to make the global maps for both data sets with discrete breaks
make_point_map_discrete <- function(df, lon, lat, size_col, fill_color, title_text,
                           legend_title = "Amount of data",
                           breaks = NULL,
                           max_pt_mm = MAX_DOT_MM, base_pt = 8.5, title_pt = 8) {
  
  pts  <- st_as_sf(df, coords = c(lon, lat), crs = 4326) |> sf::st_transform(crs_robin)
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

################################################################################
# 5b) MAP FUNCTION WITH BREAKS CONTIUNOUSLY SPECIFIED
################################################################################

# function to make the global maps for both data sets with continuous breaks
make_point_map_cont <- function(df, lon, lat, size_col, fill_color, title_text,
                           legend_title = "Amount of data",
                           max_pt_mm = MAX_DOT_MM, min_pt_mm = 1.25,
                           base_pt = 8.5, title_pt = 8,
                           trans = "sqrt") {
  # build sf points in Robinson projection
  pts <- sf::st_as_sf(df, coords = c(lon, lat), crs = 4326) |>
    sf::st_transform(crs_robin)
  
  # label formatter for legend
  lab_fmt <- scales::label_number(accuracy = 1, big.mark = ",")
  
  # keep only finite positive values for sizing; others plotted but not sized
  vals <- pts[[size_col]]
  vals_ok <- is.finite(vals) & vals > 0
  pts[[size_col]][!vals_ok] <- NA_real_
  
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
    # continuous size scaling
    scale_size_continuous(
      name   = legend_title,
      range  = c(min_pt_mm, max_pt_mm),   # size in mm for smallest/largest dots
      trans  = trans,                     # e.g., "sqrt" (default), "log10", or "identity"
      labels = lab_fmt
    ) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_classic(base_size = base_pt) +
    theme(
      plot.title       = element_text(size = title_pt, face = "bold",
                                      lineheight = 1.0, margin = margin(b = 1)),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.3),
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


#################################################################################
# 5c) BUILD MAPS WITH BREAKS SPECIFIED MANUALLY
################################################################################

# use the helper functions to get the breaks required for the nos of matrices
comp_breaks <- get_compadre_breaks_manual(matrix_locations$number_matrices)
# now build the map from the function
p_matrices <- make_point_map_discrete(
  matrix_locations, "Lon", "Lat", "number_matrices",
  col_green,
  "A) COMPADRE — points sized by number of matrices",
  legend_title = "Matrices per site",
  breaks = comp_breaks)
p_matrices

# use the helper functions to get the breaks required for the nos of species in each site
drag_breaks <- get_quantile_breaks(cover_locations$number_species_occ_in_site)
# now build the map from the function
p_cover <- make_point_map(
  cover_locations, "longitude", "latitude", "number_species_occ_in_site",
  col_orange,
  "B) DRAGNet — points sized by number of species records",
  legend_title = "Species records per site",
  breaks = drag_breaks
)
p_cover

#################################################################################
# BUILD MAPS WITH BREAKS CONTINUOUS
################################################################################

# use the helper functions to get the breaks required for the nos of matrices
comp_breaks <- get_compadre_breaks(matrix_locations$number_matrices)
# now build the map from the function
p_matrices <- make_point_map_cont(
  matrix_locations, "Lon", "Lat", "number_matrices",
  col_green,
  "A) COMPADRE — points sized by number of matrices",
  legend_title = "Matrices per site")
p_matrices
matrix_locations


################################################################################
# 4) PUT TOGETHER FINAL FIGURE 
# load the hypothesis panel (created in another software)
# and make then accompanying maps
################################################################################

# ── Panel A: title as a separate row (no overlap) ─────────────────────────────
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

################################################################################
# ── Layout: AB on top, CC full-width bottom ───────────────────────────────────
################################################################################

design <- "
AB
CC
"
p_full <- (p_graphic + p_matrices + p_cover) +
  patchwork::plot_layout(design = design,
              heights = c(TOP_ROW_REL, BOT_ROW_REL))
p_full

# ── Export ─────────────────────────────────────────────────────────────────────
ragg::agg_tiff("figures/Fig_1_map.tiff",
               width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM,
               units = "mm", res = 600, compression = "lzw")
print(p_full); dev.off()

