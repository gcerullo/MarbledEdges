# Now for this script we are gonna use the real covariates to predict for all of the sample points that were actually sampled?  

library(terra)

model <- readRDS("Models/pc1_interaction_model.rds")
source('scripts/02_OrganiseMurreletData.R')

#==================================
#Information on Murrelet survey sites: (Valente 2022)
#Most surveys were conducted around proposed timber harvest sites with 1 survey station per 8–10 ha.
#Site-level sampling effort varied with number of stations, but station size
#(200-m radius circle) was standardized (Evans Mack et al., 2003)
# so we conducted all analyses at the station level. 
#Stations were surveyed iteratively until murrelet breeding activity was recorded at the site 
#or the minimum number of required surveys was conducted#
#(5–9 surveys in each of 2 years) (Evans Mack et al., 2003).

#now for each point between 1990 and 2020 across PNW we extract the following: 
#dist coast
#hab amount 100 and 2000 (from murrelet SDM) (FROM SDMs)
#edge density 100 and 2000 (need canopy cover CANCOV_CON from https://lemma.forestry.oregonstate.edu/data/structure-maps)
#ownershipodf       
#ownershipor       
#ownershipwa  
#canopy density 
#canopy cover (to determine hard edge proportion). We then used annual GNN map products (Ohmann & Gregory, 2002 [https://lemmadownload.forestry.oregonstate.edu/]) to identify pixels with open canopies, including urban areas and forests with conifer cover ≤40%. Murrelet habitat pixels adjacent to an open canopy cell were further classified as a hard edge, and we used the distribution of hard edges to represent regional fragmentation (Figure 2).

#Get Mammu SDMS  ####
# List all subfolders (production targets) under the 'Rasters' directory
sdm_folder <- list.files("Rasters/MAMU_SDMs", full.names = TRUE, recursive = FALSE)
sdm_spat <- rast(sdm_folder)
plot(sdm_spat[[36]])

#get the ownership 
ownership <- vect("Rasters/land_ownership/All_merge.shx")
plot(ownership)



#for stations that were consistently zeros 
sum(analysisSurveys$detected) #747 detections 
unique(analysisSurveys$id) 
zero_stations <- analysisSurveys %>% group_by(id) %>% summarise(tot_detections = sum(detected)) %>%  
  filter(tot_detections ==0)

#Find the sites (and respective stations) that were allowed to be harvested
#PROBLEM - stations are nested in sites; but I have no way of knowing which stations belong 
#to which site. Yet site-level harvest is determined by cumulative station-level information. 

#1. Get actual site-level informartion 
#2. Filter sites that were likely harvested (ie meet rules for enough murrelet absence)
#3. Predict occupancy for these sites from models built both with without this data (to avoid circularity): naive model, interaction model, 

#plot how edge proportion and ocean condition have changed over time for the marbled murrlet. 

#for each point make two predictions - the naive occupancy and the fragmentation*ocean condition prediction

#During bad ocean years, is the magnified effect of fragmentation leading to underestimating occupancy; and is this especially the case with distance to 
#coast? 






