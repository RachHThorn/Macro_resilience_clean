# R Thornley
# 19/11/2024
# Project: P1_COMPADRE_DRAGNET
# Script: S4_get_demographic_metrics_COMPADRE.R

rm(list = ls())

################################################################################
# Instructions
################################################################################

# 1) Download and filter the compadre data base
# 2) Calculate the Life History Traits and resilience metrics
# 3) Check data and filter out unreasonable values / save results

################################################################################
# 1) Download and filter the compadre data base
################################################################################ 

# Download the compadre data base using the cdb_fetch command from Rcompadre
Com_dat <- cdb_fetch("compadre")

# not all the matrices in compadre are good enough quality to give us the metrics we need 
# these can be filtered
# function cdb_flag gives us info on potential issues with the matrices
# it adds a column to the database so we can filter it
Com_flags <- cdb_flag(Com_dat)
Com_flags$check_NA_A # missing values in matA
Com_flags$check_NA_U # missing values in matU
# similar check_NA vars do the same
Com_flags$check_irreducible #  matA is irreducible - usually biologically unreasonable - cause some calculations to fail
Com_flags$check_ergodic # Non-Ergodic matrices are usually biologically unreasonable - cause some calculations to fail
Com_flags$check_surv_gte_1 # contains matrix with survival rates greater than 1 - biologically unreasonable

names(Com_flags)

# We need to drop the matrices that have particular issues
# we need matrices that are defensible
# we have to have matrices that can be split into transitions related to reproduction vs. survival in
# order to run the calculations
Com_clean <- subset(Com_flags,
                    check_ergodic == T &
                      check_irreducible == T &
                      check_primitive == T &
                      MatrixTreatment == "Unmanipulated" &
                      MatrixSplit == "Divided" &
                      MatrixDimension > 1 &
                      check_NA_A == FALSE & 
                      check_NA_U == FALSE & 
                      check_NA_F == FALSE &
                      SurvivalIssue < 1.00)

# we can also query the data for species from the DRAGNet list
# read in the shared list of species found in both compadre and dragnet from Script 1
Drag_taxa <- read_csv("results/common_species_drag_comp.csv")
Drag_taxa <- Drag_taxa %>% pull(x)
Drag_taxa_new <- str_replace_all(Drag_taxa, "_", " ") # replace underscore with space
Drag_taxa_new
# now we can filter the COMPADRE data for the DRAGNet species
shared_data_more <- cdb_check_species(Com_clean, Drag_taxa_new, return_db = TRUE)
unique(shared_data_more$SpeciesAccepted) # 73 species
dim(shared_data_more) # 373 matrices 

################################################################################
# 2) Calculate the Life History Traits and resilience metrics
################################################################################

# Calculate demographic metrics for as many MPMs found in compadre as possible
# then filter for the species that overlap the two data bases

#Converting projection interval to numeric and assuming an annual basis when not stated
Com_clean$ProjectionInterval <- as.numeric(Com_clean$ProjectionInterval)
Com_clean$ProjectionInterval[which(is.na(Com_clean$ProjectionInterval))] <- 1

# save length of data object for later
long <- dim(Com_clean)[1] 

# Pull out the metadata we want:
output <- cdb_flatten(Com_clean)[,c("MatrixID", "SpeciesAuthor", "SpeciesAccepted", "Family", 
                             "Class", "Kingdom", "Authors", "DOI_ISBN", 
                             "YearPublication", "Country", "Continent",
                             "Ecoregion", "MatrixPopulation", "MatrixTreatment",
                             "MatrixStartYear", "MatrixStartMonth", "MatrixStartSeason", 
                             "MatrixEndYear", "MatrixEndMonth", "MatrixEndSeason",
                             "MatrixComposite", "Lat", "Lon", "MatrixSplit", 
                             "MatrixFec", "MatrixDimension", "OrganismType")]

# output the meta data from the data base so the matrix locations can be plotted out
write_csv(output, "results/meta_data_compadre.csv")

# Add a column for whether the spp is in dragnet:
output$DRAGNet<- ifelse(output$SpeciesAccepted %in% Drag_taxa_new, TRUE, FALSE)

# Add columns that we're going to fill in with the loop:
output$T_generation<- NA
output$R0<- NA
output$percapita_repro<- NA
output$Lmean<- NA
output$age_repro<- NA
output$Lambda<- NA
output$Reactivity<- NA
output$FirstStepAtt<- NA

dat <- Com_clean

for (i in 1:long){
  
  #The calculations here employed define the beginning of life when an individual become established. Thus, we do not consider transitions from the "prop" stages
  lifeStages <- dat$mat[[i]]@matrixClass$MatrixClassOrganized
  prop <- (which(lifeStages == "prop"))
  active <- (which(lifeStages == "active" | lifeStages=="Active"))
  dorm <- (which(lifeStages == "dorm"))
  
  # This block of code converts everything to a 1-year timestep:
  Umat <- as.matrix(dat$mat[[i]]@matU)^(1/dat$ProjectionInterval[i])
  Umat[is.na(Umat)] <- 0
  Fmat <- as.matrix(dat$mat[[i]]@matF)^(1/dat$ProjectionInterval[i])
  Fmat[is.na(Fmat)] <- 0
  Cmat <- as.matrix(dat$mat[[i]]@matC)^(1/dat$ProjectionInterval[i])
  Cmat[is.na(Cmat)] <- 0
  Amat <- Umat+Fmat+Cmat
  
  # This block of code re-arranges the matrices so that the stages go in the order prop, active, dorm
  try(reArrangeVector <- c(prop,active,dorm))
  try(Umat <- as.matrix(Umat[reArrangeVector,reArrangeVector]))
  rownames(Umat) <- colnames(Umat)
  try(Fmat <- as.matrix(Fmat[reArrangeVector,reArrangeVector]))
  rownames(Fmat) <- colnames(Fmat)
  try(Cmat <- as.matrix(Cmat[reArrangeVector,reArrangeVector]))
  rownames(Cmat) <- colnames(Cmat)
  try(Amat <- as.matrix(Amat[reArrangeVector,reArrangeVector]))
  rownames(Amat) <- colnames(Amat)
  
  #First non-propagule stage -- THIS PART DOES NOT WORK
  lifeStages<- lifeStages[reArrangeVector]
  notProp <- min(which(lifeStages != "prop"))
  if (is.infinite(notProp)) {notProp <- 1}
  
  #Various life history traits from Rage
  output[i,"Lmean"] <- life_expect_mean(Umat, start=notProp) #Mean life expectancy obtained by markov chain approaches. 
  output[i,"Lmax"] <- longevity(Umat, start=notProp) #The age a which survivorship falls to some critical proportion.
  output$age_repro[i] <- mature_age(Umat,Fmat,notProp) # Mean age at maturity (in the same units as the matrix population model sampling periodicity)
  
  #Net reproductive rate (R0; Caswell 2001, p 126)"The mean number of offspring by which an individual will be replaced by the end of its live.
  output$R0[i] <- net_repro_rate(Umat,Fmat,notProp)
  
  #lambda
  output$Lambda[i] <- max(Re(eigen(Amat, only.values=TRUE)$values))
  
  #Generation time (T; Caswell 2001, p 129)
  output$T_generation[i]<- Rage::gen_time(Umat, Fmat, method = "age_diff") 
  # Note that this generation time measure does not take into account the way
  # that we're starting life at the first non-propagule stage. I think we could use method = 'cohort' and then 
  
  # Average per-capita reproduction:
  wmean<- Re(eigen(Amat)$vectors[,1]) # stable stage distribution
  # Rescale to sum to 1 (proportions):
  wmean<- wmean/sum(wmean)
  # Multiply the stable distribution by Fmat to get a cohort of offspring:
  offspring<- Fmat%*%wmean
  # Calculate the proportion of the population that is sexually reproducing:
  repro_active<- which(colSums(Fmat)>0)
  # The per-capita reproduction is the sum of the offspring cohort, divided by
  # the reproductively active individuals:
  output$percapita_repro[i]<- sum(offspring)/sum(wmean[repro_active])
  
  #Reactivity
  output$Reactivity[i] <- reac(Amat,bound="upper")
  
  #First step attenuation
  output$FirstStepAtt[i] <- reac(Amat,bound="lower")
}

output

################################################################################
# 3) Check data and filter out unreasonable values / save results
################################################################################

# we can now filter these matrices for the species in DRAGNet
names(output)
output %>% group_by(SpeciesAccepted) %>% tally() # we have metrics for 460 species
nos_matrix <- output %>% filter(DRAGNet == TRUE) %>% group_by(SpeciesAccepted) %>% tally()
nos_matrix # we have 42 species
nos_matrix %>% filter(n > 1) # we have 25 species for which we have more than one matrix

# some of the metrics are not reasonable values though
# just tidy this up a bit
# we need to change any inf values to NA
# dat <- output %>% mutate(across(29:37, ~na_if(abs(.), Inf)))
names(output)

# we need to get rid of ridiculous data here and put NA values instead
check <-
  output %>% 
  dplyr::select(SpeciesAccepted, MatrixID, DRAGNet, T_generation:Lmax) %>%
  mutate(ID = row_number()) %>%
  pivot_longer(cols = T_generation:Lmax, names_to = "demo_var", values_to = "value") %>%
  mutate(value = case_when(demo_var == "T_generation" & value > 100 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "R0" & value > 2000 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "percapita_repro" & value > 2000 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "Lmean" & value > 200 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "age_repro" & value > 20 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "Lambda" & value > 2 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "Reactivity" & value > 5000 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "FirstStepAtt" & value > 5000 ~ NA, TRUE ~ value)) %>%
  mutate(value = case_when(demo_var == "L_max" & value < 80 ~ NA, TRUE ~ value))

# just check the number of species we have here now
check %>% filter(DRAGNet == TRUE) %>% summarise(n_species = n_distinct(SpeciesAccepted))
# there are just 42 taxa 
# now save the dataframe (this is all metrics for all species in compadre but with a DRAGNet flag)
write_csv(check, "results/all_COMPADRE_metrics.csv")




