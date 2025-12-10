# R Thornley
# 05/11/2025
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
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
crs_robin <- "+proj=robin"
world_p <- sf::st_transform(world, crs_robin)

# Build a projected bbox by transforming a WGS84 lat–lon window to Robin
make_robin_bbox <- function(lat_min = -30, lat_max = 75,
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
crop_bbox <- make_robin_bbox(lat_min = -10, lat_max = 75)

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
get_dragnet_breaks_manual <- function(z, breaks = c(1, 10, 45)) {
  # keep only finite positive values
  z <- z[is.finite(z) & z > 0]
  if (!length(z)) return(NULL)
  
  # ensure min and max coverage
  min_z <- min(z)
  max_z <- max(z)
  
  br <- sort(unique(c(min_z, breaks, max_z)))
  
  br
}

# Breaks for the COMPADRE map categories
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

################################################################################
# 5) MAP CREATION FUNCTION WITH BREAKS DISCRETELY/MANUALLY SPECIFIED
################################################################################

# function to make the global maps for both data sets with discrete breaks
make_point_map_discrete <- function(df, lon, lat, size_col, fill_color, title_text,
                                    legend_title = "",
                                    breaks = NULL,
                                    max_pt_mm = MAX_DOT_MM * 1.5, 
                                    min_pt_mm = 1.5,
                                    base_pt = 8.5, 
                                    title_pt = 16) {
  
  pts  <- sf::st_as_sf(df, coords = c(lon, lat), crs = 4326) |> sf::st_transform(crs_robin)
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
      legend.title     = element_text(size = base_pt + 1), # change title size here
      legend.text      = element_text(size = base_pt + 8), # change label size here
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
# 6) BUILD MAPS WITH BREAKS SPECIFIED MANUALLY
################################################################################

# use the helper functions to get the breaks required for the nos of matrices
comp_breaks <- get_compadre_breaks_manual(matrix_locations$number_matrices)
# now build the map from the function
p_matrices <- make_point_map_discrete(
  matrix_locations, "Lon", "Lat", "number_matrices",
  col_green,
  "A) COMPADRE — Matrix population models per site",
  legend_title = "",
  breaks = comp_breaks)
p_matrices


# now build the map from the function
p_matrices <- make_point_map_discrete(
  matrix_locations, "Lon", "Lat", "number_matrices",
  col_green,
  legend_title = "",
  title_text = "",
  breaks = comp_breaks)
p_matrices


ggsave(filename = "figures/matrices_map_01_12.png", plot =  p_matrices,
       width = 180, height = 130, units = "mm", dpi = 300)


# use the helper functions to get the breaks required for the nos of species in each site
drag_breaks <- get_dragnet_breaks_manual(cover_locations$number_species_occ_in_site)
# now build the map from the function
p_cover <- make_point_map_discrete(
  cover_locations, "longitude", "latitude", "number_species_occ_in_site",
  col_orange,
  "B) DRAGNet — Number of species records per site",
  breaks = drag_breaks
)
p_cover

p_cover <- make_point_map_discrete(
  cover_locations, "longitude", "latitude", "number_species_occ_in_site",
  col_orange,
  title_text = "",
  breaks = drag_breaks
)
p_cover

ggsave(filename = "figures/cover_map_01_12.png", plot =  p_cover,
       width = 180, height = 130, units = "mm", dpi = 300)



#################################################################################
# 7) IMPORT HYPOTHESIS FIGURE and CREATE TITLE FOR IT
################################################################################

# ── Panel A: title as a separate row (no overlap) ─────────────────────────────
graphic_path <- "figures/Fig_1A_06_11_2025.tiff"
img <- magick::image_read(graphic_path) |> magick::image_trim()
img

# Image-only panel
p_graphic_img <- cowplot::ggdraw() +
  cowplot::draw_image(img, x = 0, y = 1, width = 1, height = 1,
                      hjust = 0, vjust = 1, interpolate = TRUE) +
  theme(plot.margin = margin(2, 2, 2, 2))
p_graphic_img


# Title-only panel (appears above the image)
p_graphic_title <- ggplot() +
  theme_void() +
  labs(title = "A) Overview of hypotheses and study variables") +
  theme(
    plot.title = element_text(hjust = 0.55, vjust = 0,
                              face = "bold", size = 20,
                              family = "Arial",
                              margin = margin(b = 2))
  )
p_graphic_title

# Stack title above image inside Panel C
p_graphic <- patchwork::wrap_plots(
  p_graphic_title,
  p_graphic_img,
  ncol = 1,
  heights = c(0.005, 0.7)  # adjust to move title higher/lower
)
p_graphic

################################################################################
# 8) PUT TOGETHER FINAL FIGURE 
# load the hypothesis panel (created in another software)
# and make then accompanying maps
################################################################################

# Layout: A (hypotheses) on top - full width , BC on bottom 
p_graphic
p_cover
p_matrices

###############################################################################
# this doesn't work
##############################################################################

library(patchwork)

design <- "AABC"
p_full <- (p_graphic + p_matrices + p_cover) + patchwork::plot_layout(design = design)
p_full

p_full <- (p_graphic + (p_matrices + p_cover)) +
  patchwork::plot_layout(design = design,
                         heights = c(TOP_ROW_REL, BOT_ROW_REL))
p_full

p_full <- p_graphic / (p_matrices | p_cover) +
  patchwork::plot_layout(heights = c(TOP_ROW_REL, BOT_ROW_REL))
p_full


# Export
ragg::agg_tiff("figures/Fig_1_map.tiff",
               width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM,
               units = "mm", res = 600, compression = "lzw")
print(p_full); dev.off()


plot_element_1 <- (p_cover | p_matrices)
plot_element_2 <- p_graphic
plot_element_3 <- plot_element_2 / plot_element_1 +
  plot_layout (nrow = 2,
               heights = c(1, 1),
               widths = c(1,1))
plot_element_3

################################################################################
# using ggarrange
################################################################################

library(ggpubr)

# bottom row: p_matrices | p_cover
bottom_row <- ggarrange(
  p_matrices, p_cover,
  ncol = 2, nrow = 1,
  widths = c(1, 1),
  font.label = list(size = 12, face = "bold")
)
bottom_row

# full figure: p_graphic on top, bottom_row below
p_full <- ggarrange(
  p_graphic_img, bottom_row,
  ncol = 1, nrow = 2,
  heights = c(TOP_ROW_REL, BOT_ROW_REL),  # your row height ratios
  align = "v",
  common.legend = FALSE
)

p_full


# full figure: p_graphic on top, bottom_row below
p_full <- ggarrange(
  p_graphic, bottom_row,
  ncol = 1, nrow = 2,
  align = "v",
  common.legend = FALSE,        # set TRUE if you want a shared legend
  legend = "top"                # where to place the shared legend if TRUE
)

p_full

