#organise murrelet data for analysis 

# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)
library(ggcorrplot)
 options(scipen = 999)

 #inputs 
can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif") #canopy cover of conifers
PC1 <- read.csv("Outputs/PC1_scaled_inverted.csv") 

#murrelet survey detection data
landscapeVars  <- read.csv("Inputs/LandscapeVariables.csv")
siteData  <- read.csv("Inputs/siteData.csv")
surveys  <- read.csv("Inputs/MurreletSurveys.csv") %>% dplyr::select(-X) #multipe repeats a some stations
pointsInRoiAndSamplingWindow  <- read.csv("Inputs/pointsInRoiAndSamplingWindow.csv") %>%  
  mutate(crs = "EPSG:32610")

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


#take a look at the stratification of survey points 

# Create SpatVector in UTM Zone 10N (EPSG:32610)
points_utm <- vect(pointsInRoiAndSamplingWindow, geom = c("utmE", "utmN"), crs = "EPSG:32610")
# Transform to EPSG:5070 (NAD83 / Conus Albers)
points_5070 <- project(points_utm, "EPSG:5070")
plot(can_cov, main = "Sampled points on can cov")
plot(points_5070, add = TRUE, col = "blue", pch = 16, cex = 0.5)

#read in PC1 data calculated in script 1 #
PC1_all <- PC1 %>%  
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


#  create new variables `edgeArea100` and `edgeArea2000` based on specific edge and habitat values
analysisSites = siteData %>% 
  mutate(edgeArea100 = (edgeRook_100_40 * (pi * 100^2)) / (habAmountDich_100 * (pi * 100^2) + 1)) %>% 
  mutate(edgeArea2000 = (edgeRook_2000_40 * (pi * 2000^2)) / (habAmountDich_2000 * (pi * 2000^2) + 1)) %>%  
  # join by PC1 in prev year
  left_join(PC1_t1, by = c("year"= "year_t1"))%>%
  #only keep data for which we have data on PC1 value for 1 yr bfore
  filter(year >= earliest_PC1_year) 


#join site data with murrelet surveys and PC1 
analysisSurveys = analysisSites %>% 
  dplyr::select(id, year,utmE,utmN) %>% 
  left_join(surveys, by = c('id', 'year')) %>% 
  #only keep data for which we have data on PC1 value
  filter(year >= earliest_PC1_year) %>% 
  mutate(year = as.factor(year))


#quick summary
detections <- analysisSurveys %>% filter(detected == 1)
detections %>% count() #761 detections 
 analysisSites %>% dplyr::select(id) %>% unique() %>% count()  #20149
analysisSurveys %>% count() #29323 surveys


# Calculate means and standard deviations for key variables in `analysisSites` and `analysisSurveys`
# This list will store these statistics for later use or reference
meansAndSds <- data.frame(
  meanCanopy100 = mean(analysisSites$canopy100, na.rm = TRUE),
  meanCoastDist = mean(analysisSites$coastDist, na.rm = TRUE),
  meanHabAmountDich2000 = mean(analysisSites$habAmountDich_2000, na.rm = TRUE),
  meanHabAmountDich100 = mean(analysisSites$habAmountDich_100, na.rm = TRUE),
  meanEdgeDens2000 = mean(analysisSites$edgeRook_2000_40, na.rm = TRUE),
  meanEdgeDens100 = mean(analysisSites$edgeRook_100_40, na.rm = TRUE),
  meanDoy = mean(analysisSurveys$doy, na.rm = TRUE),
  meanConDens100 = mean(analysisSites$conDens100, na.rm = TRUE),
  sdCanopy100 = sd(analysisSites$canopy100, na.rm = TRUE),
  sdConDens100 = sd(analysisSites$conDens100, na.rm = TRUE),
  sdCoastDist = sd(analysisSites$coastDist, na.rm = TRUE),
  sdHabAmountDich2000 = sd(analysisSites$habAmountDich_2000, na.rm = TRUE),
  sdHabAmountDich100 = sd(analysisSites$habAmountDich_100, na.rm = TRUE),
  sdYoungPlusNonFor2000 = sd(analysisSites$youngPlusNonFor_2000_40, na.rm = TRUE),
  sdYoungPlusNonFor100 = sd(analysisSites$youngPlusNonFor_100_40, na.rm = TRUE),
  sdEdgeRook2000 = sd(analysisSites$edgeRook_2000_40, na.rm = TRUE),
  sdEdgeRook100 = sd(analysisSites$edgeRook_100_40, na.rm = TRUE),
  sdDoy = sd(analysisSurveys$doy, na.rm = TRUE)
)

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
hist(analysisSites$coastDist, breaks =100)


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

#---------------------------------
#Get a summary of the elevational istribution of points points: 
#---------------------------------
#Get an elevation DEM in the correct projection (takes 20mins or so)
# can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif") #canopy cover of conifers
# elevation <- rast("Rasters/DEMs/dem30m.tif")
# 
# #Step 1: Reproject elevation to EPSG:5070 (same as can_cov)
# elev_proj <- project(elevation, can_cov, method = "bilinear")
# 
# # Step 2: Crop to the extent of can_cov
# elev_crop <- crop(elev_proj, can_cov)
#
#elevational_range <- terra::extract(elev_crop, points_5070, bind = TRUE)
#elevational_range_df <-   as.data.frame(elevational_range) #%>% 
  # dplyr::select(id, cancov_2020) %>%  
  # rename(point_id =  id) 
#write.csv(elevational_range_df, "Inputs/pointsInRoiAndSamplingWindow_withDEM.csv")

elevational_range_df <- read.csv("Inputs/pointsInRoiAndSamplingWindow_withDEM.csv")
hist(elevational_range_df$dem30m) 
(max)

#plot covariates of model 
analysisSurveys

detected_sites <- analysisSurveys %>% dplyr::select(id, detected)

allsites <- analysisSites %>% 
dplyr::select(coastDist,habAmountDich_100,habAmountDich_2000,edgeRook_100_40,edgeRook_2000_40) %>%  
  rename(
    `Distance to coast`     = coastDist,
    `Hab amount 100 m`      = habAmountDich_100,
    `Hab amount 2000 m`     = habAmountDich_2000,
    `Edge amount 100 m`     = edgeRook_100_40,
    `Edge amount 2000 m`    = edgeRook_2000_40
  ) %>% 
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>% 
  ggplot( aes(x = Value)) +
  geom_histogram(fill = "#999999", color = "white", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal(base_size = 18) +  # Minimal theme with larger base font size
  theme(
    legend.position = "top",  # Position the legend at the top
    legend.title = element_blank(),  # Increase legend title size for clarity
    legend.text = element_text(size = 14),  # Increase legend text size for better readability
    axis.title = element_text(size = 16),  # Increase axis title font size
    axis.text = element_text(size = 14),  # Increase axis label font size
    panel.grid.major = element_blank(),  # No major gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a thin border around the plot
  )

detectedsites <- analysisSites %>% 
  

#if we only want to consider sites where there was a detection
left_join(detected_sites, by = "id") %>%  
  filter(detected == 1) %>% 
dplyr::select(coastDist,habAmountDich_100,habAmountDich_2000,edgeRook_100_40,edgeRook_2000_40) %>%  
  rename(
    `Distance to coast`     = coastDist,
    `Hab amount 100 m`      = habAmountDich_100,
    `Hab amount 2000 m`     = habAmountDich_2000,
    `Edge amount 100 m`     = edgeRook_100_40,
    `Edge amount 2000 m`    = edgeRook_2000_40
  ) %>% 
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>% 
  ggplot( aes(x = Value)) +
  geom_histogram(fill = "#999999", color = "white", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal(base_size = 18) +  # Minimal theme with larger base font size
  theme(
    legend.position = "top",  # Position the legend at the top
    legend.title = element_blank(),  # Increase legend title size for clarity
    legend.text = element_text(size = 14),  # Increase legend text size for better readability
    axis.title = element_text(size = 16),  # Increase axis title font size
    axis.text = element_text(size = 14),  # Increase axis label font size
    panel.grid.major = element_blank(),  # No major gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a thin border around the plot
  )


#save ouptuts #####
#save  unmarked dataframe ####
saveRDS(analysisData, file = "Outputs/analysisDataUnmarked.rds")

#Export clearance survey site covariates figure ####s
ggsave(
  filename = "Figures/clearance_survey_all_site_covars.png",               # File path and name
  plot = allsites ,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

ggsave(
  filename = "Figures/clearance_survey_mamu_detected_site_covars.png",               # File path and name
  plot = detectedsites,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

