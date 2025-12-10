# Climatic PCA to compadre species in DRAGNet vs COMPADRE for representativity
# Created by Rob Salguero-Gomez, rob.salguero@biology.ox.ac.uk
# 29 Oct 2025

#Install / load dependencies
required_pkgs <- c("terra", "sp", "geodata")
to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if(length(to_install)) install.packages(to_install)
library(terra)
library(sp)
library(geodata)

#File paths
dragnet_csv <- "results/DRAGNet_species_locations_final.csv"
compadre_csv <- "results/COMPADRE_species_locations_final.csv"
out_dir <- "results/PCA_clim_plots_DRAGNet_COMPADRE"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#Helper function to standardise coord column names (lon/lat)
standardise_coords <- function(df){
  names_lower <- tolower(names(df))
  lon_idx <- which(names_lower %in% c("lon","longitude","long","x"))[1]
  lat_idx <- which(names_lower %in% c("lat","latitude","y"))[1]
  if(is.na(lon_idx) || is.na(lat_idx)) stop("Could not detect lon/lat columns. Found columns: ", paste(names(df), collapse = ", "))
  names(df)[lon_idx] <- "lon"
  names(df)[lat_idx] <- "lat"
  df
}

#Helper function to detect species column robustly
detect_species_column <- function(df){
  names_lower <- tolower(names(df))
  #Common candidate names in order of preference
  candidates <- c("species", "scientificname", "scientific_name", "species_name", "sp", "taxon", "taxon_name")
  for (cand in candidates) {
    idx <- which(names_lower == cand)
    if(length(idx) == 1) return(names(df)[idx])
  }
  #Fallback: try to find a column with many repeated values but not all unique and with non-numeric values
  non_numeric_cols <- sapply(df, function(x) !(is.numeric(x) || all(is.na(x))))
  if(any(non_numeric_cols)){
    #Choose the non-numeric column with the largest number of repeated values (likely species)
    reps <- sapply(names(df)[non_numeric_cols], function(nm){
      v <- df[[nm]]
      #Measure distinctness: lower distinct count relative to rows = candidate
      (length(unique(na.omit(v))) / max(1, sum(!is.na(v))))
    })
    best <- names(reps)[which.min(reps)]
    #Only accept if distinctness < 0.9 (i.e., many repeats) and not lon/lat
    if(reps[best] < 0.9 && !(tolower(best) %in% c("lon","latitude","lat","longitude","long","x","y"))){
      return(best)
    }
  }
  # unable to detect
  return(NA_character_)
}

#Read and prepare data
dr <- read.csv(dragnet_csv, stringsAsFactors = FALSE, check.names = FALSE)
cp <- read.csv(compadre_csv, stringsAsFactors = FALSE, check.names = FALSE)

#Standardise lon/lat (will error if not found)
dr <- standardise_coords(dr)
cp <- standardise_coords(cp)

#Detect species column for each dataframe
spcol_dr <- detect_species_column(dr)
spcol_cp <- detect_species_column(cp)

if(is.na(spcol_dr)) stop("Could not detect species column in DRAGNet file. Columns found: ", paste(names(dr), collapse = ", "))
if(is.na(spcol_cp)) stop("Could not detect species column in COMPADRE file. Columns found: ", paste(names(cp), collapse = ", "))

#Normalise the species column name to 'species' but only after we've detected the correct column
names(dr)[names(dr) == spcol_dr] <- "species"
names(cp)[names(cp) == spcol_cp] <- "species"

#Clean coordinates and species NA rows
dr$lon <- as.numeric(dr$lon); dr$lat <- as.numeric(dr$lat)
cp$lon <- as.numeric(cp$lon); cp$lat <- as.numeric(cp$lat)
dr <- dr[complete.cases(dr[, c("lon", "lat", "species")]), ]
cp <- cp[complete.cases(cp[, c("lon", "lat", "species")]), ]

#Download climate data (WorldClim via geodata)
message("Downloading WorldClim bioclimatic data (~10 arcmin)...")
bio <- geodata::worldclim_global(var = "bio", res = 10, path = tempdir())
#Ensure naming is tractable
names(bio) <- paste0("bio", 1:19)

#Crop to extent of points (with small buffer)
ext <- terra::ext(min(c(dr$lon, cp$lon)) - 2,
                  max(c(dr$lon, cp$lon)) + 2,
                  min(c(dr$lat, cp$lat)) - 2,
                  max(c(dr$lat, cp$lat)) + 2)
bio <- terra::crop(bio, ext)

#Function to extract climate variables for each dataset
extract_clim <- function(df, bio_rast){
  pts <- cbind(df$lon, df$lat)
  vals <- terra::extract(bio_rast, pts)[, -1, drop = FALSE]  # drop ID
  out <- cbind(df, vals)
  out <- out[complete.cases(out[, grep("^bio", names(out))]), , drop = FALSE]
  return(out)
}

dr_clim_all <- extract_clim(dr, bio)
cp_clim_all <- extract_clim(cp, bio)



#Build species list from actual species column (not site names)
species_list <- unique(as.character(dr_clim_all$species))
message(sprintf("Found %d unique species in DRAGNet file.", length(species_list)))

if(length(species_list) == 0) stop("No species found in DRAGNet after extraction. Check species column and coordinate overlap with climate grid.")

#Loop over species and compute PCA/hull/plots
for(sp in species_list){
  message("Processing species: ", sp)
  dr_sp <- dr_clim_all[as.character(dr_clim_all$species) == sp, , drop = FALSE]
  if(nrow(dr_sp) < 3){
    warning(sprintf("Skipping %s: fewer than 3 DRAGNet points (%d).", sp, nrow(dr_sp)))
    next
  }
  cp_sp <- cp_clim_all[as.character(cp_clim_all$species) == sp, , drop = FALSE]
  
  clim_dr <- dr_sp[, grep("^bio", names(dr_sp)), drop = FALSE]
  #Remove zero-variance columns
  varvals <- apply(clim_dr, 2, var, na.rm = TRUE)
  clim_dr <- clim_dr[, varvals > 0, drop = FALSE]
  if(ncol(clim_dr) < 2){
    warning(sprintf("Skipping %s: not enough climatic variables with variance after filtering.", sp))
    next
  }
  
  pca <- prcomp(clim_dr, center = TRUE, scale. = TRUE)
  dr_scores <- predict(pca, newdata = clim_dr)
  dr_xy <- dr_scores[, 1:2, drop = FALSE]
  
  #Project COMPADRE points if any
  cp_xy <- matrix(nrow = 0, ncol = 2)
  inside_flag <- logical(0)
  if(nrow(cp_sp) > 0){
    clim_cp <- cp_sp[, colnames(clim_dr), drop = FALSE]
    valid_cp <- complete.cases(clim_cp)
    if(any(valid_cp)){
      cp_scores <- predict(pca, newdata = clim_cp[valid_cp, , drop = FALSE])
      cp_xy <- cp_scores[, 1:2, drop = FALSE]
      #Point-in-polygon test
      hull_idx <- chull(dr_xy[,1], dr_xy[,2])
      hull_coords <- dr_xy[hull_idx, , drop = FALSE]
      pip <- sp::point.in.polygon(cp_xy[,1], cp_xy[,2], hull_coords[,1], hull_coords[,2])
      inside_flag <- pip >= 1  # 1/2/3 are inside/boundary/vertex
      #Keep only the cp_sp rows that had valid climate (if needed to save summary)
      cp_sp <- cp_sp[valid_cp, , drop = FALSE]
    } else {
      cp_xy <- matrix(nrow = 0, ncol = 2)
    }
  } else {
    hull_idx <- chull(dr_xy[,1], dr_xy[,2])
    hull_coords <- dr_xy[hull_idx, , drop = FALSE]
  }
  
  #If hull_idx not defined (should be) ensure it is computed
  if(!exists("hull_idx")) hull_idx <- chull(dr_xy[,1], dr_xy[,2])
  hull_coords <- dr_xy[hull_idx, , drop = FALSE]
  
  #Plot
  safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", sp)
  pngfile <- file.path(out_dir, paste0("PCA_clim_DRAGNet_COMPADRE_", safe_name, ".png"))
  png(pngfile, width = 1200, height = 1000, res = 150)
  par(mar = c(5,5,4,2) + 0.1)
  
  #Determine limits
  all_x <- dr_xy[,1]
  all_y <- dr_xy[,2]
  if(nrow(cp_xy) > 0){ all_x <- c(all_x, cp_xy[,1]); all_y <- c(all_y, cp_xy[,2]) }
  xlim <- range(all_x, na.rm = TRUE); ylim <- range(all_y, na.rm = TRUE)
  xpad <- max(0.1, diff(xlim) * 0.08); ypad <- max(0.1, diff(ylim) * 0.08)
  xlim <- c(xlim[1]-xpad, xlim[2]+xpad); ylim <- c(ylim[1]-ypad, ylim[2]+ypad)
  
  plot(NA, xlim = xlim, ylim = ylim, xlab = "PC1", ylab = "PC2",
       main = paste0("Climatic PCA: ", sp))
  points(dr_xy[,1], dr_xy[,2], pch = 20, col = "darkgrey")
  polygon(c(hull_coords[,1], hull_coords[1,1]),
          c(hull_coords[,2], hull_coords[1,2]),
          border = "blue", lwd = 2, col = rgb(0,0,1, alpha = 0.15))
  if(nrow(cp_xy) > 0){
    cols <- ifelse(inside_flag, "darkgreen", "red")
    points(cp_xy[,1], cp_xy[,2], pch = 21, bg = cols, col = "black", cex = 1.2)
    legend("topright", legend = c("DRAGNet (sites)", "DRAGNet convex hull", "COMPADRE inside", "COMPADRE outside"),
           pch = c(20, NA, 21, 21), pt.bg = c(NA, NA, "darkgreen", "red"),
           lty = c(NA, 1, NA, NA), col = c("darkgrey", "blue", "black", "black"), bty = "n")
  } else {
    legend("topright", legend = c("DRAGNet (sites)", "DRAGNet convex hull"),
           pch = c(20, NA), lty = c(NA,1), col = c("darkgrey","blue"), bty = "n")
  }
  var_pc1 <- round(100 * (pca$sdev[1]^2 / sum(pca$sdev^2)), 1)
  var_pc2 <- round(100 * (pca$sdev[2]^2 / sum(pca$sdev^2)), 1)
  mtext(paste0("PC1: ", var_pc1, "%; PC2: ", var_pc2, "%"), side = 3, line = 0.5, cex = 0.9)
  dev.off()
  
  message("Saved plot to: ", pngfile)
  if(nrow(cp_xy) > 0){
    message(sprintf("Species %s: DRAGNet pts = %d, COMPADRE pts = %d; inside = %d; outside = %d",
                    sp, nrow(dr_xy), nrow(cp_xy), sum(inside_flag), sum(!inside_flag)))
  } else {
    message(sprintf("Species %s: DRAGNet pts = %d; no COMPADRE points (with climate) for this species.", sp, nrow(dr_xy)))
  }
}

message("Done. Plots are in: ", normalizePath(out_dir))
