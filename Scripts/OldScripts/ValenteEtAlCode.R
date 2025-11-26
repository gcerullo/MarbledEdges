
# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)

#read in murrelet data 
load("Inputs/ValenteEtAlEnvironment.RData")

#save outputs as csvs 
write.csv(landscapeVars, "Inputs/LandscapeVariables.csv")
write.csv(siteData, "Inputs/SiteData.csv")
write.csv(surveys, "Inputs/MurreletSurveys.csv")
write.csv(pointsInRoiAndSamplingWindow, "Inputs/pointsInRoiAndSamplingWindow.csv")

#export site id and location to extract age locations for each point: 
spatialSites <- analysisSites %>%  select(id, year, utmN, utmE)
plot(spatialSites)
write.csv(spatialSites,"Outputs/murrelet_samplingSites.csv")

#Information on predictor variables ####
        

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


# Analysis setup
# Filter `siteData` to include only sites with a specific edge condition (e.g., gt2kmFakeEdge==1) 
# and create new variables `edgeArea100` and `edgeArea2000` based on specific edge and habitat values
analysisSites = siteData %>% 
  filter(gt2kmFakeEdge == 1) %>% 
  mutate(edgeArea100 = (edgeRook_100_40 * (pi * 100^2)) / (habAmountDich_100 * (pi * 100^2) + 1)) %>% 
  mutate(edgeArea2000 = (edgeRook_2000_40 * (pi * 2000^2)) / (habAmountDich_2000 * (pi * 2000^2) + 1))

# Join the `analysisSites` data with `surveys` on `id` and `year`
# Creates `analysisSurveys` with only relevant survey information
analysisSurveys = analysisSites %>% 
  select(id, year) %>% 
  left_join(surveys, by = c('id', 'year'))

# Calculate means and standard deviations for key variables in `analysisSites` and `analysisSurveys`
# This list will store these statistics for later use or reference
meansAndSds = list(meanCanopy100 = mean(analysisSites$canopy100, na.rm=T),
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

# Standardize selected variables in `analysisSites` (centered and scaled) for later modeling
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
  arrange(id)

# Standardize day-of-year variable `doy` in `analysisSurveys`
# Creates a squared version `scaleDoy2` for non-linear effects in later modeling
analysisSurveys = analysisSurveys %>% 
  mutate(scaleDoy = scale(doy, center = TRUE, scale = TRUE)) %>% 
  mutate(scaleDoy2 = scaleDoy^2)

# Reshape detection data (`y`) and arrange it for occupancy modeling
# Each row represents a site with columns for different surveys (detections)
y = analysisSurveys %>% 
  pivot_wider(id_cols = id, names_from = survey, values_from = detected) %>% 
  arrange(id) %>% 
  select(-id)

# Prepare survey-level covariates for the occupancy model, creating a list of matrices
surveyCovs = list(
  'scaleDoy' = analysisSurveys %>% 
    pivot_wider(id_cols = id, names_from = survey, values_from = scaleDoy) %>% 
    arrange(id) %>% 
    select(-id),
  'scaleDoy2' = analysisSurveys %>% 
    pivot_wider(id_cols = id, names_from = survey, values_from = scaleDoy2) %>% 
    arrange(id) %>% 
    select(-id)
)

# Compute correlation matrix for standardized variables in `analysisSites`
analysisSites %>% 
  select(scaleCanopy100, scaleConDens100, scaleEdgeDens100, scaleCoastDist, scaleHabAmount100,
         scaleHabAmount2000, scaleEdgeDens2000, scaleEdgeArea100, scaleEdgeArea2000) %>% 
  cor()

# Create an `unmarkedFrameOccu` object `analysisData` for occupancy modeling
# This includes detection data (`y`), site-level covariates (`siteCovs`), and observation-level covariates (`obsCovs`)
analysisData = unmarkedFrameOccu(y = y, siteCovs = analysisSites, obsCovs = surveyCovs)

#check resurveys

# Extract the detection history matrix
y_mat <- analysisData@y  # Rows = sites, Columns = visits

# Number of sites
n_sites <- nrow(y_mat)

# Number of visits per site (i.e., non-NA columns)
visits_per_site <- apply(y_mat, 1, function(x) sum(!is.na(x)))

# Summarize
summary_stats <- list(
  total_sites = n_sites,
  mean_visits = mean(visits_per_site),
  median_visits = median(visits_per_site),
  min_visits = min(visits_per_site),
  max_visits = max(visits_per_site),
  visits_distribution = table(visits_per_site)
)

summary_stats
sum(visits_per_site >= 2)

#####
# Occupancy Modeling
# Check the presence/absence pattern by calculating the maximum detection across surveys per site
table(apply(analysisData@y, 1, FUN = 'max', na.rm = TRUE))

# Various occupancy models with different covariate structures
# The `occu` function is used to estimate occupancy probability and detection probability

# 1. Model with only detection covariates
a = Sys.time()
p = occu(~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 1, data = analysisData)
Sys.time() - a

# 2. Model including detection and `scaleCoastDist` as a site covariate
a = Sys.time()
dist = occu(~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ scaleCoastDist, data = analysisData,
            starts = c(coef(p)[1], 0, coef(p)[2:11]))
Sys.time() - a

# 3. Model with year as an additional covariate
a = Sys.time()
year = occu(~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ scaleCoastDist + as.factor(year), data = analysisData,
            starts = c(coef(dist)[1:2], rep(0, 28), coef(dist)[3:12]))
Sys.time() - a

# Save intermediate results to file
save(list = ls(), file = 'Models/ManuscriptResults.RData')

# 4. Model with additional habitat amount and edge density covariates at different scales
a = Sys.time()
model1 = occu(~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ scaleCoastDist + as.factor(year) + scaleHabAmount100 + scaleEdgeDens100 + scaleHabAmount2000 + scaleEdgeDens2000, data = analysisData,
              starts = c(coef(year)[1:30], rep(0, 4), coef(year)[31:40]))
Sys.time() - a

# Save results
save(list = ls(), file = 'Models/ManuscriptResults.RData')

# 5. Model including interaction terms between habitat amount and edge density at multiple scales
a = Sys.time()
model2 = occu(~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ scaleCoastDist + as.factor(year) + scaleHabAmount100 + scaleEdgeDens100 + scaleHabAmount100 * scaleEdgeDens100 + scaleHabAmount2000 + scaleEdgeDens2000 + scaleHabAmount2000 * scaleEdgeDens2000, data = analysisData,
              starts = c(coef(model1)[1:32], 0, coef(model1)[33:34], 0, coef(model1)[35:44]))
Sys.time() - a

# Save final model results
save(list = ls(), file = 'Models/ManuscriptResults.RData')

# 6. Additional model exploring coastal distance and fragmentation interactions
a = Sys.time()
model3 = occu(~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ scaleCoastDist + as.factor(year) + scaleHabAmount100 + scaleEdgeDens100 + scaleCoastDist * scaleEdgeDens100 + scaleHabAmount2000 + scaleEdgeDens2000 + scaleCoastDist * scaleEdgeDens2000, data = analysisData,
              starts = coef(model2))
Sys.time() - a

# Save results for coastal distance interaction model
save(list = ls(), file = 'Models/ManuscriptResults.RData')
