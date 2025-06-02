library(tidyverse)
library(MatchIt)

#read in inputs predicted from 9 ####
frag_linear_loss <- readRDS("hab_points_fragmentation_linear_hab_loss_cumalative_final_model_8thMay2025.rds")
final2020 <- readRDS("Outputs/PNW_2020_extracted_covars.rds") %>%   #read in starting occupancy for 2020 from scrippt 8
  as.data.frame() #%>%  

long_lat <- final2020 %>% select(point_id, x, y)
#add in information on elevation at each point 
pt_elevation <- read.csv("Outputs/elevation_at_each_point.csv") 

#---------------------------------------------------------

#what was the elevational range of murrelets in dataset from pre-clearance surveys? 
elevational_range_df <- read.csv("Inputs/pointsInRoiAndSamplingWindow_withDEM.csv")
hist(elevational_range_df$dem30m) 
max(elevational_range_df$dem30m, na.rm = TRUE) #no surveys above 1636.42 metres! 

#get quantiles of elevation  
Allelev = as.numeric(quantile(elevational_range_df$dem30m, 1, na.rm = TRUE))
q95elev = as.numeric(quantile(elevational_range_df$dem30m, 0.95, na.rm = TRUE))
q90elev = as.numeric(quantile(elevational_range_df$dem30m, 0.90, na.rm = TRUE))
q80elev = as.numeric(quantile(elevational_range_df$dem30m, 0.80, na.rm = TRUE))

### add key information 
frag_linear_loss <- frag_linear_loss %>% left_join(pt_elevation) %>%  left_join(long_lat) %>%  
#carry out initial filter to remove extremes and mask urban/snow etc 
  filter(cancov_2020 > -1) %>% #remove non-forest habitat points
  filter(OceanYear == "Good Ocean Years") %>% 
  filter( distance_to_coastline < 100000) %>% #remove points v far from coast 
  filter(dem30m < elev_thresh) %>% 
  filter(hab_loss_amount == 0)

#scale back to original values
frag_linear_loss <- frag_linear_loss %>%
  mutate(habAmountDich_2000 = (scaleHabAmount2000 * meansAndSds$sdHabAmountDich2000) + meansAndSds$meanHabAmountDich2000) %>%  
  mutate(edgeRook_2000_40 = (scaleEdgeDens2000 * meansAndSds$sdEdgeRook2000) + meansAndSds$meanEdgeDens2000) %>%  
  mutate(habAmountDich_100 = (scaleHabAmount100 * meansAndSds$sdHabAmountDich100) + meansAndSds$meanHabAmountDich100) %>%  
  mutate(edgeRook_100_40 = (scaleEdgeDens100 * meansAndSds$sdEdgeRook100) + meansAndSds$meanEdgeDens100)


#===============================
#carry out matching 
#============================
binary_df <- frag_linear_loss %>%  
 filter(ownership %in% c("Private Industrial", "Federal")) %>%
  mutate(ownership_binary = case_when(
    ownership == "Federal" ~ 1,
    ownership == "Private Industrial" ~ 0
  ))  %>%  select(point_id, distance_to_coastline, dem30m, ownership_binary)

names(frag_linear_loss)

# No matching; constructing a pre-match matchit object
m.out0 <- matchit(ownership_binary ~ distance_to_coastline + dem30m,
                  data = binary_df,
                  method = NULL,
                  distance = "glm")

# Checking balance prior to matching
summary(m.out0)

#We can see severe imbalances as measured by the standardized mean differences (Std. Mean Diff.), variance ratios (Var. Ratio), and empirical cumulative distribution function (eCDF) statistics. Values of standardized mean differences and eCDF statistics close to zero and values of variance ratios close to one indicate good balance, and here many of them are far from their ideal values.
#Values of standardized mean differences and eCDF statistics close to zero #
#and values of variance ratios close to one indicate good balance, and here many of them are far from their ideal values.

# 1:1 NN PS matching w/o replacement
m.out1 <- matchit(ownership_binary ~ distance_to_coastline + dem30m,
                  data = binary_df,
                  method = "nearest",
                  distance = "glm")

# Checking balance after NN matching
summary(m.out1, un = FALSE)

plot(m.out1, type = "density", interactive = FALSE,
     which.xs = ~distance_to_coastline + dem30m)
#black line = treated/federal 
#grey line = untreated/private industrial 

# Full matching on a probit PS
m.out2 <- matchit(ownership_binary ~ distance_to_coastline + dem30m,
                  data = binary_df,
                  method = "full",
                  distance = "glm",
                  link = "probit")
m.out2