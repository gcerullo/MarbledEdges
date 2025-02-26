#organise murrelet data for analysis 

# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)
library(ggcorrplot)


#Information on predictor variables ####
{
  
  #ownership:  federal (e.g., national parks and forests), state (e.g., state parks and forests), private industrial (predominantly forest industry lands), private nonindustrial (e.g., small family-owned land holdings), and other (e.g., lands owned by local governments, Indigenous peoples, or under private conservation easements)
  
  #MODEL COVARIATES 
  # meanCanopy100: Mean canopy cover within a 100-meter radius
  # meanConDens100: Mean coniferous tree density within 100 meters
  # meanCoastDist: Mean distance to coast from each site,
  # meanhabAmountDich2000: Mean dichotomized habitat amount within 2000 meters (landscape-level)
  # meanHabAmountDich100: Mean dichotomized habitat amount within 100 meters.
  # meanYoungPlusNonFor2000: Mean proportion of young forest plus non-forest areas within 2000 meters, indicating early-successional habitat.
  # meanYoungPlusNonFor100: Mean proportion of young forest plus non-forest areas within 100 meters.
  # meanEdgeRook2000: Mean edge density (rook length) within 2000 meters
  # meanEdgeRook100: Mean edge density within 100 meters.
  # meanDoy: Mean day of the year (DOY) for surveys, accounting for seasonal changes.
  # meanEdgeArea100: Mean edge area within 100 meters
  # meanEdgeArea2000: Mean edge area within 2000 meters
  
  # sdCanopy100: Standard deviation of canopy cover within 100 meters.
  # sdConDens100: Standard deviation of coniferous tree density within 100 meters.
  # sdCoastDist: Standard deviation of distance to coast.
  # sdhabAmountDich2000: Standard deviation of dichotomized habitat amount within 2000 meters.
  # sdHabAmountDich100: Standard deviation of dichotomized habitat amount within 100 meters.
  # sdYoungPlusNonFor2000: Standard deviation of young forest plus non-forest areas within 2000 meters.
  # sdYoungPlusNonFor100: Standard deviation of young forest plus non-forest areas within 100 meters.
  # sdEdgeRook2000: Standard deviation of edge density within 2000 meters.
  # sdEdgeRook100: Standard deviation of edge density within 100 meters.
  # sdDoy: Standard deviation of day of the year (DOY) for surveys.
  # sdEdgeArea100: Standard deviation of edge area within 100 meters.
  # sdEdgeArea2000: Standard deviation of edge area within 2000 meters.
}

#read in inputs 
#save outputs as csvs 
landscapeVars  <- read.csv("Inputs/LandscapeVariables.csv")
siteData  <- read.csv("Inputs/siteData.csv")
surveys  <- read.csv("Inputs/MurreletSurveys.csv") %>% dplyr::select(-X) #multipe repeats a some stations
pointsInRoiAndSamplingWindow  <- read.csv("Inputs/pointsInRoiAndSamplingWindow.csv")

#read in PC1 data
PC1_all <- read.csv("Outputs/PC1_scaled_inverted.csv") %>%  
  dplyr::select(Value_scaled, Year) %>%  #dplyr::select PC1 which has already been inverted and scaled
  rename(PC1 = Value_scaled, 
         year = Year) %>%  
  mutate(year_t1 = year +1, 
#give PC1 value for the previous year; joining on year_t1 thus actually gives PC1 vals for the year before
                  PC1_t1 = PC1)

#PC1 in prev year
PC1_t1 <- PC1_all %>% dplyr::select(year_t1, PC1_t1)
earliest_PC1_year <- PC1_t1 %>%
  filter(year_t1 == min(year_t1)) %>% 
  pull(year_t1)

## Filter `siteData` to include only sites with a specific edge condition (e.g., gt2kmFakeEdge==1) 
# and create new variables `edgeArea100` and `edgeArea2000` based on specific edge and habitat values
analysisSites = siteData %>% 
  filter(gt2kmFakeEdge == 1) %>% 
  mutate(edgeArea100 = (edgeRook_100_40 * (pi * 100^2)) / (habAmountDich_100 * (pi * 100^2) + 1)) %>% 
  mutate(edgeArea2000 = (edgeRook_2000_40 * (pi * 2000^2)) / (habAmountDich_2000 * (pi * 2000^2) + 1)) %>%  
  # join by PC1 in prev year
  left_join(PC1_t1, by = c("year"= "year_t1"))%>%
  #only keep data for which we have data on PC1 value
  filter(year >= earliest_PC1_year)

#join site data with murrelet surveys and PC1 
analysisSurveys = analysisSites %>% 
  dplyr::select(id, year) %>% 
  left_join(surveys, by = c('id', 'year')) %>% 
  
  #only keep data for which we have data on PC1 value
  filter(year >= earliest_PC1_year)

#quick summary
detections <- analysisSurveys %>% filter(detected == 1)
detections %>% count() #915 individual detections from 31,879 surveys


# Calculate means and standard deviations for key variables in `analysisSites` and `analysisSurveys`
# This list will store these statistics for later use or reference
{meansAndSds = 
  list(meanCanopy100 = mean(analysisSites$canopy100, na.rm=T),
                   meanConDens100 = mean(analysisSites$conDens100, na.rm=T),
                   meanCoastDist = mean(analysisSites$coastDist, na.rm=T),
                   meanhabAmountDich2000 = mean(analysisSites$habAmountDich_2000, na.rm=T),
                   meanHabAmountDich100 = mean(analysisSites$habAmountDich_100, na.rm=T),
                   meanYoungPlusNonFor2000 = mean(analysisSites$youngPlusNonFor_2000_40, na.rm=T),
                   meanYoungPlusNonFor100 = mean(analysisSites$youngPlusNonFor_100_40, na.rm=T),
                   meanEdgeRook2000 = mean(analysisSites$edgeRook_2000_40, na.rm=T),
                   meanEdgeRook100 = mean(analysisSites$edgeRook_100_40, na.rm=T),
                   meanDoy = mean(analysisSurveys$doy, na.rm=T),
                   meanEdgeArea100 = mean(analysisSites$edgeArea100, na.rm=T),
                   meanEdgeArea2000 = mean(analysisSites$edgeArea2000, na.rm=T),
                   sdCanopy100 = sd(analysisSites$canopy100, na.rm=T),
                   sdConDens100 = sd(analysisSites$conDens100, na.rm=T),
                   sdCoastDist = sd(analysisSites$coastDist, na.rm=T),
                   sdhabAmountDich2000 = sd(analysisSites$habAmountDich_2000, na.rm=T),
                   sdHabAmountDich100 = sd(analysisSites$habAmountDich_100, na.rm=T),
                   sdYoungPlusNonFor2000 = sd(analysisSites$youngPlusNonFor_2000_40, na.rm=T),
                   sdYoungPlusNonFor100 = sd(analysisSites$youngPlusNonFor_100_40, na.rm=T),
                   sdEdgeRook2000 = sd(analysisSites$edgeRook_2000_40, na.rm=T),
                   sdEdgeRook100 = sd(analysisSites$edgeRook_100_40, na.rm=T),
                   sdDoy = sd(analysisSurveys$doy, na.rm=T),
                   sdEdgeArea100 = sd(analysisSites$edgeArea100, na.rm=T),
                   sdEdgeArea2000 = sd(analysisSites$edgeArea2000, na.rm=T))
  }

# Standardize dplyr::selected variables in `analysisSites` (centered and scaled) for later modeling
#PC1 has already been scaled in prev script (01)
analysisSites = analysisSites %>% 
  mutate(scaleCanopy100 = scale(canopy100, center=T, scale=T),
         scaleConDens100 = scale(conDens100, center=T, scale=T),
         scaleCoastDist = scale(coastDist, center=T, scale=T),
         scaleHabAmount2000 = scale(habAmountDich_2000, center=T, scale=T),
         scaleHabAmount100 = scale(habAmountDich_100, center=T, scale=T),
         scaleYoungPlusNonFor2000 = scale(youngPlusNonFor_2000_40, center=T, scale=T),
         scaleYoungPlusNonFor100 = scale(youngPlusNonFor_100_40, center=T, scale=T),
         scaleEdgeDens2000 = scale(edgeRook_2000_40, center=T, scale=T),
         scaleEdgeDens100 = scale(edgeRook_100_40, center=T, scale=T),
         scaleEdgeArea100 = scale(edgeArea100, center=T, scale=T),
         scaleEdgeArea2000 = scale(edgeArea2000, center=T, scale=T),
         scaleUtmN = scale(utmN, center=T, scale=T)) %>% 
  arrange(id)%>% 
  dplyr::select(-X)

hist(analysisSites$habAmountDich_2000, breaks = 100)
hist(analysisSites$edgeRook_2000_40, breaks =100)
# Standardize day-of-year variable `doy` in `analysisSurveys`
# Creates a squared version `scaleDoy2` for non-linear effects in later modeling
analysisSurveys = analysisSurveys %>% 
  mutate(scaleDoy = scale(doy)) %>% 
  mutate(scaleDoy2 = scaleDoy^2)

# Reshape detection data (`y`) and arrange it for occupancy modeling
# Each row represents a site with columns for different surveys (detections)
y = analysisSurveys %>% 
  pivot_wider(id_cols = id, names_from = survey, values_from = detected) %>% 
  arrange(id) %>% 
  dplyr::select(-id)

# Prepare survey-level covariates for the occupancy model, creating a list of matrices
surveyCovs = list(
  'scaleDoy' = analysisSurveys %>% 
    pivot_wider(id_cols = id, names_from = survey, values_from = scaleDoy) %>% 
    arrange(id) %>% 
    dplyr::select(-id),
  'scaleDoy2' = analysisSurveys %>% 
    pivot_wider(id_cols = id, names_from = survey, values_from = scaleDoy2) %>% 
    arrange(id) %>% 
    dplyr::select(-id)
)

# Compute correlation matrix for standardized variables in `analysisSites`
cor_matrix <- analysisSites %>% 
  dplyr::select(scaleCanopy100, scaleConDens100, scaleEdgeDens100, scaleCoastDist, scaleHabAmount100,
         scaleHabAmount2000, scaleEdgeDens2000, scaleEdgeArea100, scaleEdgeArea2000, PC1_t1) %>% 
  cor()

# Create an `unmarkedFrameOccu` object `analysisData` for occupancy modeling
# This includes detection data (`y`), site-level covariates (`siteCovs`), and observation-level covariates (`obsCovs`)
analysisData = unmarkedFrameOccu(y = y, siteCovs = analysisSites, obsCovs = surveyCovs)

#save  unmarked dataframe ####
saveRDS(analysisData, file = "Outputs/analysisDataUnmarked.rds")
