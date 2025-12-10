# R Thornley
# 12/12/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S7_species_strategies_calculate_T0

rm(list = ls())

################################################################################
# Instructions
################################################################################

# 1) get list of DRAGNet species and desirable traits
# 2) source the necessary leaf traits from TRY
# 3) put this data into the Pierce et al. 2017 spreadsheet (StrateFy) to get CSR per species
# 4) use these values to get a measure of community CSR distance using the function 'specialization' from Ricotta et al. 2023

###############################################################################
# TRY DATA LOAD AND FILTER FOR REQUIRED TRAITS - EXPORT FOR SPREADSHEET #
################################################################################
# 
# # Load in the data ordered through TRY for all the species in DRAGNet
# # use the rtry package to load and query data
# # warning: data sets are quite large (GBs)
# 
# dat_1 <- rtry::rtry_import("data/37407.txt")
# dat_2 <- rtry_import("data/38952.txt")
# dat <- rbind(dat_1, dat_2)
# dat <- dat %>% arrange(TraitID)
# 
# n_distinct(dat$AccSpeciesID) # 928 species - DRAGNet / TRY overlap
# unique(dat$AccSpeciesName)
# unique(dat$TraitID) # this is a trait number
# unique(dat$TraitName) # this is a trait descriptor
# 
# IDs <- dat %>% dplyr::select(TraitID) %>% drop_na() %>% arrange(TraitID) %>% pull(TraitID)
# unique(IDs) # we have selected 179 traits
# 
# # for the CSR strategies we need 3 leaf traits
# # LDMC
# # Leaf Area
# # SLA
# 
# # look at a list of the trait names
# unique(dat$TraitName)
# 
# # Leaf Dry Mass per Leaf fresh mass (LDMC)
# LDMC <- dat %>% 
#   dplyr::select(AccSpeciesName, AccSpeciesID, ObservationID, TraitID, TraitName, StdValue, UnitName) %>%
#   filter(TraitName == "Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)") %>%
#   mutate(StdValue = as.numeric(StdValue))%>%
#   group_by(AccSpeciesName) %>% 
#   mutate(LDMC_mean = mean(StdValue)) %>%
#   drop_na()
# LDMC$UnitName # g/g - ratio or a percentage
# 
# LDMC_mean <- LDMC %>%
#   group_by(AccSpeciesName) %>%
#   summarise(LDMC_mean = mean(LDMC_mean))
# # 546 obs
# 
# # this needs to be converted to a precentage not a ratio
# LDMC_mean$LDMC_mean <- LDMC_mean$LDMC_mean*100
# 
# # Specific Leaf Area (SLA)
# # this trait is captured by a few different categories
# wanted_SLA <- c("Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded",
#                 "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included",
#                 "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded") 
# SLA <- dat %>% 
#   dplyr::select(AccSpeciesName, AccSpeciesID, ObservationID, TraitID, TraitName, StdValue, UnitName) %>%
#   filter(TraitName %in% wanted_SLA) %>%
#   mutate(StdValue = as.numeric(StdValue))%>%
#   group_by(AccSpeciesName) %>% 
#   mutate(SLA_mean = mean(StdValue)) %>%
#   drop_na()
# SLA$UnitName # mm2 mg-1
# 
# SLA_mean <- SLA %>%
#   group_by(AccSpeciesName) %>%
#   summarise(SLA_mean = mean(SLA_mean))
# # 622 obs
# 
# # Leaf Area codes (1, 3108, 3109, 3117, 3110, 3111, 3112, 3113, 3114)
# wanted_LA <- c("Leaf area (not yet refined if leaf or leaflet)",
#                "Leaf area (in case of compound leaves: leaf, petiole excluded)",
#                "Leaf area (in case of compound leaves: leaflet, petiole and rachis excluded)",
#                "Leaf area (in case of compound leaves: leaf, petiole included)",
#                "Leaf area (in case of compound leaves: leaflet, petiole and rachis included)", 
#                "Leaf area (in case of compound leaves: leaf, undefined if petiole in- or excluded)", 
#                "Leaf area (in case of compound leaves: leaflet, undefined if petiole and rachis are in- or excluded)",
#                "Leaf area (in case of compound leaves undefined: leaf or leaflet; petiole and rachis in- excluded)")
# LA <- dat %>% 
#   dplyr::select(AccSpeciesName, AccSpeciesID, ObservationID, TraitID, TraitName, StdValue, UnitName) %>%
#   dplyr::filter(TraitName %in% wanted_LA) %>%
#   mutate(StdValue = as.numeric(StdValue))%>%
#   group_by(AccSpeciesName) %>% 
#   mutate(LA_mean = mean(StdValue)) %>%
#   drop_na()
# LA$UnitName # mm2
# 
# LA_mean <- LA %>%
#   group_by(AccSpeciesName) %>%
#   summarise(LA_mean = mean(LA_mean))
# # 514 obs
# 
# all <- LDMC_mean %>% left_join(SLA_mean) %>% 
#   left_join(LA_mean) %>% drop_na() %>% # some of the species don't have data for all 3 categories
#   dplyr::select(AccSpeciesName, LA_mean, LDMC_mean, SLA_mean)
# # export this data for use with the excel spreadsheet StrateFy
# all$AccSpeciesName
# write_csv(all, "results/TRY_DRAGNet_StateFy_data_more.csv")
# 
# ################################################################################
# # VIEW CREATED CSR STRATEGIES FROM StateFy spreadsheet #
# ################################################################################
# 
# # the file imported is from the StateFy Excel spreadsheet
# # strategies were calculated for all the 
# # Look at strategies and overlap with compadre
# 
# # read in the CSR strategies as calculated by the StateFy spreadsheet
# csr <- read_csv("results/CSR_results_StateFy.csv")
# n_distinct(csr$AccSpeciesName) # we have the leaf trait values for 428 species
# names(csr)
# 
# # visualise these strategies
# ggtern(data = csr, aes(x = percent_r, y = percent_c, z = percent_s)) + 
#   geom_point() +
#   ggtitle("CSR for 428 species in DRAGNet")
# # colour the points by csr categorisation
# ggtern(data = csr, aes(x = percent_r, y = percent_c, z = percent_s, colour = `Strategy class`)) + 
#   geom_point() +
#   theme_bw() +
#   labs(title("CSR for 428 species in DRAGNet"))
# 
# # tidy up this figure 
# csr %>% rename('Percent C' = percent_c, 'Percent R' = percent_r, 'Percent S' = percent_s )
# 
# # we need a list of the 40 species we are using in our models
# # we have this as an output from thee robust OLS modelling file
# focal <- readRDS("results/List_taxa_OLS_mods.R")
# length(focal) # 40 species 
# 
# # the csr file needs the names changing so they have an underscore
# csr <- csr %>% mutate(AccSpeciesName = str_replace(AccSpeciesName, " ", "_"))
# focal_csr <- csr %>% filter(AccSpeciesName %in% focal)
# dim(focal_csr)
# # for 30 of my focal species I have csr strategies - not that many
# 
# # get a plot of those 55 strategies
# ggtern(data = focal_csr, aes(x = percent_r, y = percent_c, z = percent_s, colour = `Strategy class`)) + 
#   geom_point() +
#   theme_bw()+
#   ggtitle("CSR for 30 focal species with compadre/dragnet overlap")
# 
# # NOTE: for an appendix we could also visualise this triangle as we did the other data overlap figures
# # make this plot so that DRAGNet species specialisation is highlighted


################################################################################
# DREGREE OF CSR SPECIALIASATION #
################################################################################

# to use the 'specialization' function from Ricotta et al. we need two matrices
# 1) CSR table (a table of the csr values) with taxon as rownames; cols c, s, r values from Pierce et al. 
# 2) Community table (a table of species abundances per quadrat); taxon are colnames / plots are rownames

# 1) CSR table (rows = taxa; cols = C,S,R)
csr <- read_csv("data/CSR_results_StateFy.csv", show_col_types = FALSE) %>%
  mutate(AccSpeciesName = str_replace_all(AccSpeciesName, " ", "_"))

CSRtable <- csr %>%
  dplyr::select(AccSpeciesName, percent_c, percent_s, percent_r) %>%
  rename(Taxon = AccSpeciesName, C = percent_c, S = percent_s, R = percent_r) %>%
  as.data.frame()

unique(CSRtable$Taxon)

# robust name parsing
CSRtable <- CSRtable %>%
  mutate(across(c(C, S, R), ~ readr::parse_number(as.character(.x))))
rownames(CSRtable) <- CSRtable$Taxon
CSRtable$Taxon <- NULL
CSRtable <- as.matrix(CSRtable)
stopifnot(is.numeric(CSRtable), !is.null(rownames(CSRtable)), !is.null(colnames(CSRtable)))
csr_taxa <- rownames(CSRtable)

# 2) DRAGNet community long table, nested by site

# Helper to normalise year_trt once, but keep your pipe layout
norm_year <- function(x) paste0("T", readr::parse_number(as.character(x)))

# T0 only

comm_by_site_0 <- read_csv("results/DRAGNet_T0_T1_all.csv", show_col_types = FALSE) %>%
  filter(year_trt == "T0") %>%
  mutate(
    New_taxon = str_replace_all(New_taxon, " ", "_"),
    unique    = paste(site_name, year_trt, trt, block, sep = "_")
  ) %>%
  dplyr::select(site_name, New_taxon, unique, max_cover) %>%
  filter(New_taxon %in% csr_taxa) %>%
  group_by(site_name) %>% tidyr::nest() %>% ungroup()


# specialization() (Ricotta et al.)
specialization <- function(comm, CSRtable, tol = 1e-8){
  if(!inherits(comm, "data.frame") && !inherits(comm, "matrix")) stop("Incorrect definition of argument comm")
  if(!inherits(CSRtable, "data.frame") && !inherits(CSRtable, "matrix")) stop("Incorrect definition of argument CSRtable")
  if(!is.numeric(comm)) stop("Incorrect definition of argument comm")
  if(any(comm < -tol))  stop("Incorrect definition of argument comm")
  comm[comm < -tol] <- 0
  if(!is.numeric(CSRtable)) stop("Incorrect definition of argument CSRtable")
  if(any(CSRtable < -tol))  stop("Incorrect definition of argument CSRtable")
  
  if(is.null(rownames(comm)) | is.null(colnames(comm)) |
     is.null(rownames(CSRtable)) | is.null(colnames(CSRtable)))
    stop("Tables comm and CSRtable must have rownames and colnames")
  
  if(!all(colnames(comm) %in% rownames(CSRtable))) stop("Some species names in comm are missing in CSRtable")
  if(!all(colnames(CSRtable) %in% c("C","S","R")) | ncol(CSRtable) != 3) stop("Incorrect definition of argument CSRtable")
  
  CSRtable <- CSRtable[colnames(comm), ]
  CSRtable[CSRtable < -tol] <- 0
  
  CSR   <- sweep(CSRtable, 1, rowSums(CSRtable), "/")
  Pcomm <- sweep(comm, 1, rowSums(comm), "/")
  
  FUN <- function(vx) {
    BIND <- rbind.data.frame(rep(1/3, 3), vx)
    den  <- as.vector(as.dist(dist(matrix(c(1,0,0, 1/3,1/3,1/3), 2, 3, byrow=TRUE))))
    as.vector(dist(BIND) / den)
  }
  
  sigmaSpe     <- apply(CSR, 1, FUN); names(sigmaSpe) <- rownames(CSR)
  CSRMeans     <- Pcomm %*% CSR
  sigmaonMeans <- apply(CSRMeans, 1, FUN)
  meanSigmas   <- Pcomm %*% sigmaSpe
  variability  <- (meanSigmas - sigmaonMeans) / meanSigmas
  
  sigmaCom <- cbind.data.frame(
    meanSigmas   = as.numeric(meanSigmas),
    sigmaonMeans = as.numeric(sigmaonMeans),
    variability  = as.numeric(variability)
  )
  rownames(sigmaCom) <- rownames(comm)
  
  list(sigmaSpe = sigmaSpe, sigmaCom = sigmaCom)
}

# 3) Helper: compute specialization per site (kept your shape)
get_special <- function(comm_df, CSRtable) {
  wide <- comm_df %>%
    pivot_wider(
      id_cols     = New_taxon,      # rows = taxa
      names_from  = unique,         # cols = quadrats
      values_from = max_cover,
      values_fn   = sum,
      values_fill = 0
    ) %>%
    as.data.frame()
  
  rownames(wide) <- wide$New_taxon
  wide$New_taxon <- NULL
  
  comm_mat <- t(as.matrix(wide))
  storage.mode(comm_mat) <- "numeric"
  
  if (nrow(comm_mat) == 0 || ncol(comm_mat) == 0)
    return(tibble(ID = character(), meanSigmas = numeric(),
                  sigmaonMeans = numeric(), variability = numeric()))
  
  # drop empty rows/cols
  comm_mat <- comm_mat[rowSums(comm_mat) > 0, , drop = FALSE]
  comm_mat <- comm_mat[, colSums(comm_mat) > 0, drop = FALSE]
  if (nrow(comm_mat) == 0 || ncol(comm_mat) == 0)
    return(tibble(ID = character(), meanSigmas = numeric(),
                  sigmaonMeans = numeric(), variability = numeric()))
  
  # align species
  spp <- intersect(colnames(comm_mat), rownames(CSRtable))
  if (length(spp) == 0)
    return(tibble(ID = rownames(comm_mat),
                  meanSigmas = NA_real_, sigmaonMeans = NA_real_, variability = NA_real_))
  
  comm_mat <- comm_mat[, spp, drop = FALSE]
  CSR_sub  <- CSRtable[spp, , drop = FALSE]
  
  spec <- specialization(comm_mat, CSR_sub, tol = 1e-8)
  as_tibble(spec$sigmaCom, rownames = "ID")
}

safe_get_special <- purrr::possibly(get_special,
                                    otherwise = tibble::tibble(
                                      ID = character(),
                                      meanSigmas = numeric(),
                                      sigmaonMeans = numeric(),
                                      variability = numeric()
                                    ))

# 4) Apply to each site and unnest (unchanged layout)

special_0 <- comm_by_site_0 %>%
  mutate(comm_dat = purrr::map(data, safe_get_special, CSRtable = CSRtable)) %>%
  dplyr::select(site_name, comm_dat) %>%
  tidyr::unnest(comm_dat)
special_0$time_period <- "T0"

# take the mean site, time, treatment strategy values
# in line with the other community variables
# first create the appropriate columns using separate
all <- special_0 %>% separate(ID, into = c("new_site", "year_trt", "trt","block"),
                        sep = "_", extra = "merge",  # merge all remaining pieces into 'site_name'
                        fill  = "right", remove = FALSE) %>%
  dplyr::select(!new_site) %>% 
  group_by(site_name, year_trt, trt) %>%
  mutate(meanSigmas = mean(meanSigmas), sigmaonMeans = mean(sigmaonMeans), variability = mean(variability)) %>%
  dplyr::select(!block) %>%
  dplyr::select(!ID) %>%
  distinct()

write_csv(all, "results/community_specialisation_results_T0.csv")

################################################################################



