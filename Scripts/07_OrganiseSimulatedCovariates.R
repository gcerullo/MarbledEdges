# Organise Simulated Landscape Covariates

#To fit with our models we need: 
  #1. Distance to coast ; scaleCoastDist          
  #2. Habitat amount 100 and 2km;  scaleHabAmount100 and scaleHabAmount2000                    
  #3. Edge density 100 and 2km ; scaleEdgeDens100 and scaleEdgeDens2000    
  
library(tidyverse)
library(data.table)

#read in inputs 
hab <- readRDS("Outputs/HabAmount10000pts_Simulated_055_LandscapeParams.rds")
edge <- readRDS("Outputs/EdgeDensity10000pts_Simulated_055_LandscapeParams.rdsedge_density_list.rds")
model <- readRDS("Models/pc1_interaction_model.rds")
covariates <- readRDS("Outputs/ScaledCovariates.rds") %>%   #covariates from script 04
  select(PC1_t1,scaleCoastDist,scaleDoy,scaleDoy2,scaleDoy2,OceanYear) %>%  unique()
real_murrelet_site_data <- read.csv("Inputs/siteData.csv") #to make sure we match this structure


# List all TIFF files in the folder (recursive = FALSE to avoid subfolders)
tif_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)
# Read multiple rasters at once into a single SpatRaster object
landscapes<- rast(tif_files)
plot(landscapes)

#process habitat amount ####
#process data. NB 1 is forest and 0 = plantation.  
hab_df <- hab %>% rbindlist(idcol = "landscape_name")

#num survey points per landscape and buffer_size
hab_df %>% group_by(buffer_size,landscape_name) %>% summarize(count = n(), .groups = 'drop')

hab_df<- hab_df %>% 
  pivot_wider(
    names_from = buffer_size, 
    values_from = c(frac_0, frac_1)
  ) %>%    
  mutate(habAmountDich_100 = frac_1_buffer_1, 
         habAmountDich_2000 = frac_1_buffer_20,
         point_id = id) %>% 
  select(landscape_name, point_id,habAmountDich_100,habAmountDich_2000)
  
#process edge density ####

edge_df <- edge %>%  rbindlist(idcol = "landscape_name")

#num survey points per landscape and buffer_size
edge_df %>% group_by(buffer_size,landscape_name) %>% summarize(count = n(), .groups = 'drop')

#process edge data and also combine datasets
edge_df <- edge_df %>% 
  pivot_wider(
    names_from = buffer_size, 
    values_from = c(frac_0, frac_1)
  ) %>%    
  mutate(edgeRook_100_40 = frac_1_buffer_1, 
         edgeRook_2000_40 = frac_1_buffer_20,
         point_id = id) %>% 
  select(landscape_name, point_id,edgeRook_100_40,edgeRook_2000_40)

all_df <- hab_df %>% left_join(edge_df)

#Compare real and simulated covariates ####
#nuderstand real murrelet site strucuture 
hist(real_murrelet_site_data$habAmountDich_100)
hist(real_murrelet_site_data$habAmountDich_2000)
hist(real_murrelet_site_data$edgeRook_100_40)
hist(real_murrelet_site_data$edgeRook_2000_40)


#and for our simulated landscapes 
hist(hab_df$habAmountDich_100)
hist(hab_df$habAmountDich_2000)
hist(all_df$edgeRook_100_40)
hist(all_df$edgeRook_2000_40)

#oraganise data for prediction with correct names and scaling 
prediction_df <- all_df %>% mutate(
  scaleHabAmount100 = scale(habAmountDich_100), 
  scaleHabAmount2000 = scale(habAmountDich_2000), 
  scaleEdgeDens100 = scale(edgeRook_100_40), 
  scaleEdgeDens2000 = scale(edgeRook_2000_40)) %>%  
  select(landscape_name, point_id,scaleHabAmount100,scaleHabAmount2000,scaleEdgeDens100,scaleEdgeDens2000) %>% 
  cross_join(covariates)


# Predict Occupancy with standard errors (for error ribbon)
predictions <- predict(
  model,  
  newdata = prediction_df,
  type = "state", 
  se.fit = TRUE  # Obtain standard errors for predictions
) %>% 
  rename(Occupancy = Predicted)

# Add predicted occupancy and standard errors to the data
prediction_df <- prediction_df %>% cbind(predictions) %>%  
  mutate(lower_CI = Occupancy - 1.96 * SE, 
         upr_CI = Occupancy + 1.96 * SE )

plot_data <- prediction_df %>% 
  select(landscape_name, point_id,Occupancy,OceanYear,lower_CI,upr_CI) %>% 
  mutate(landscape_numeric = as.numeric(gsub("patches_", "", landscape_name))) %>% 
  mutate(landscape_name = fct_reorder(landscape_name, landscape_numeric))

#rapid plots 
# Create the boxplot
plot_data %>%  
  ggplot( aes(x = landscape_name, y = Occupancy)) +
  geom_jitter(aes(color = landscape_name), size = 2, width = 0.1, alpha = 0.2) +  # Optional: adds points
  geom_boxplot(fill = "#56B4E9", color = "black", width = 0.5, outlier.shape = 21, outlier.size = 2) + 
  theme_classic(base_size = 14) +
  facet_wrap(~OceanYear)+# Clean theme for Nature/Science style
  labs(x = "Increasing landscape fragmentation ->", y = "Occupancy", title = "Occupancy Across Landscapes") +
  theme(
    text = element_text(size = 16, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

unique(plot_data$landscape_name)
