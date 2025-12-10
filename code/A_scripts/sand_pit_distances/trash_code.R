
################################################################################
# PCA analysis 
################################################################################

pca <- prcomp(clim_combined %>% dplyr::select(all_of(clim_vars)), scale. = TRUE)

scores <- bind_cols(
  clim_combined %>% dplyr::select(dataset),
  as.data.frame(pca$x[, 1:2]) %>% rename(PC1 = PC1, PC2 = PC2)
)

get_hull <- function(df){
  df[chull(df$PC1, df$PC2), ]
}

dr_hull <- scores %>% filter(dataset == "DR") %>% get_hull()
cp_hull <- scores %>% filter(dataset == "CP") %>% get_hull()

library(sf)

poly_dr <- st_polygon(list(as.matrix(dr_hull[, c("PC1", "PC2")]))) |> st_sfc()
poly_cp <- st_polygon(list(as.matrix(cp_hull[, c("PC1", "PC2")]))) |> st_sfc()

area_dr <- st_area(poly_dr)
area_cp <- st_area(poly_cp)
area_cp
area_overlap <- st_area(st_intersection(poly_dr, poly_cp))

overlap_ratio <- as.numeric(area_overlap) / as.numeric(st_area(st_union(poly_dr, poly_cp)))
overlap_ratio


################################################################################
# 3) create multi-variate distances: ONE SPECIES TEST - all distance metrics
################################################################################

# create two separate matrices for the different data sets with unique rows only
# scale before transforming
A <- dr_sp %>% dplyr::select(all_of(clim_vars)) %>% distinct() %>% scale() %>% as.matrix()  # dataset 1 (n × p)
B <- cp_sp %>% dplyr::select(all_of(clim_vars)) %>% distinct() %>% 
  scale(center = attr(A, "scaled:center"), scale  = attr(A, "scaled:scale")) %>% as.matrix()
# use same scaling as A dataset 2 (m × p)


# different ways of finding distances # 

# 1) Euclidean distance
D_euc <- proxy::dist(A, B, method = "euclidean")
D_euc$Species <- 
  D <- as.matrix(D)
D


# 2) Mahalanobis distances
# we need to rejoin the two matrices for this as both must undergo the same scaling / 
# PCA rotation and covariance matrix structure
AB <- rbind(A, B)
S  <- cov(AB)
Sinv <- solve(S)

mah_fun <- function(x, y) {
  diff <- x - y
  sqrt(t(diff) %*% Sinv %*% diff)
}

D_mah <- proxy::dist(A, B, method = mah_fun)
# this method fails as there are too many dimensions


# 3) do a PCA first - then find distances
# Combine A and B just to fit a *common* PCA

# What this is doing conceptually
# Join A and B, to learn the common PCA rotation and scaling.
# Split them back into A_pc and B_pc (now K-dimensional, with K < n),
# Compute pairwise distances between matrices A and B, exactly as you wanted.
# Because K is much smaller than p and ≤ (n_total − 1), you can avoid the singularity problem 
# when we try to use that affects Mahalanobis in the original p-dimensional space.

AB <- rbind(A, B)
pca <- prcomp(AB, center = TRUE, scale. = TRUE)

# Choose number of PCs (e.g. enough to explain 95% variance)
eigvals <- pca$sdev^2
cumvar  <- cumsum(eigvals) / sum(eigvals)
K <- which(cumvar >= 0.95)[1]   # first K PCs reaching 95%

# Project A and B into PC space (same rotation)
scores <- pca$x[, 1:K, drop = FALSE]

A_pc <- scores[1:nrow(A), , drop = FALSE]
B_pc <- scores[(nrow(A) + 1):nrow(scores), , drop = FALSE]

## Find distances between A and B in PC space (n × m)
D_pc <- proxy::dist(A_pc, B_pc, method = "euclidean")
D_pc <- as.matrix(D_pc)

dim(D_pc)  # n rows (A) × m cols (B)

# ok lets stick with method 3 for this analysis


################################################################################
# robust way to calculate cut-offs

# compute robust cutoffs
med_med <- median(med_dist$median_distance, na.rm = TRUE)
mad_med <- mad(med_dist$median_distance, na.rm = TRUE)

upper_cutoff_species <- med_med + 3 * mad_med

median_distances_flagged <- med_dist%>%
  mutate(
    is_outlier = median_distance > upper_cutoff_species
  )

species_to_drop <- median_distances_flagged %>%
  filter(is_outlier) %>%
  pull(Taxon)

distances <- all_distances_long %>%
  filter(!Taxon %in% species_to_drop)

