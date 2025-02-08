# Organise Simulated Landscape Covariates

#To fit with our models we need: 
  #1. Distance to coast ; scaleCoastDist          
  #2. Habitat amount 100 and 2km;  scaleHabAmount100 and scaleHabAmount2000                    
  #3. Edge density 100 and 2km ; scaleEdgeDens100 and scaleEdgeDens2000    
  
library(tidyverse)
library(data.table)
library(terra)
library(raster)
library(RColorBrewer)
library(stringr)

#read in inputs 
#----------------------------------------------------------------------------------
#read in simulated inputs for correct production target!!

#prodcution 0.24
hab <- readRDS("Outputs/production_0.24_HabAmount10000pts_SimulatedParams.rds") #from script 06
edge <- readRDS("Outputs/production_0.24_EdgeDensity10000pts_SimulatedParams.rds") #from script 06
pt_value <- readRDS("Outputs/production_0.24_point_forest_or_plantation.rds") %>% as.data.frame() %>%  
  dplyr::select(landscape_name, point_id,raster_value.lyr.1) %>%  
  rename(forest_plantation = raster_value.lyr.1) %>%  
  mutate(point_id = paste0("p",point_id))
# Define the folder where the landscape TIFF files are stored
input_folder <- "Rasters/production_0.24"

#production 0.58
# hab <- readRDS("Outputs/production_0.58_HabAmount10000pts_SimulatedParams.rds") #from script 06
# edge <- readRDS("Outputs/production_0.58_EdgeDensity10000pts_SimulatedParams.rdsedge_density_list.rds") #from script 06
# pt_value <- readRDS("Outputs/production_0.58_point_forest_or_plantation.rds") %>% as.data.frame() %>%  
#   dplyr::select(landscape_name, point_id,raster_value.lyr.1) %>%  
#   rename(forest_plantation = raster_value.lyr.1) %>%  
#   mutate(point_id = paste0("p",point_id))
# Define the folder where the landscape TIFF files are stored
#input_folder <- "Rasters/production_0.58"


tif_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)
# Read multiple rasters at once into a single SpatRaster object
landscapes<- rast(tif_files)
plot(landscapes)

#model, covariate and site inputs
model <- readRDS("Models/pc1_interaction_model.rds") #from script 03
covariates <- readRDS("Outputs/ScaledCovariates.rds") %>%   #covariates from script 04
  dplyr::select(PC1_t1,scaleCoastDist,scaleDoy,scaleDoy2,scaleDoy2,OceanYear) %>%  unique()
real_murrelet_site_data <- read.csv("Inputs/siteData.csv") # from script 02

#----------------------------------------------------------------------------------
#process habitat amount ####
#process data. NB 0 is forest and 1 = plantation.  
hab_df <- hab %>% rbindlist(idcol = "landscape_name")

#num survey points per landscape and buffer_size
hab_df %>% group_by(buffer_size,landscape_name) %>% summarize(count = n(), .groups = 'drop')

hab_df<- hab_df %>% 
  pivot_wider(
    names_from = buffer_size, 
    values_from = c(frac_0, frac_1)
  ) %>%    
  mutate(habAmountDich_100 = frac_0_buffer_1, 
         habAmountDich_2000 = frac_0_buffer_20,
         point_id = id) %>% 
  dplyr::select(landscape_name, point_id,habAmountDich_100,habAmountDich_2000)
  
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
  dplyr::select(landscape_name, point_id,edgeRook_100_40,edgeRook_2000_40)

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
  dplyr::select(landscape_name, point_id,scaleHabAmount100,scaleHabAmount2000,scaleEdgeDens100,scaleEdgeDens2000) %>% 
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
  #add information about whether the point location is forest or plantation 
  #left_join(pt_value) %>%  
  #filter only points that are actually forest 
 #  filter(forest_plantation ==0) %>% 
  dplyr::select(landscape_name, point_id,Occupancy,OceanYear,lower_CI,upr_CI) %>% 
  mutate(landscape_numeric = as.numeric(gsub("patches_", "", landscape_name))) %>% 
  mutate(landscape_name = fct_reorder(landscape_name, landscape_numeric)) 

#rapid plots 
# Create the boxplot
p024 <- plot_data %>%  
  ggplot( aes(x = landscape_name, y = Occupancy)) +
  geom_jitter(color = 'lightgrey', size = 2, width = 0.1, alpha = 0.05) +  
  geom_boxplot(fill = "#56B4E9", color = "black", width = 0.5, outlier.shape = NA, outlier.size = NA) + 
  theme_classic(base_size = 14) +
  facet_wrap(~OceanYear)+# Clean theme for Nature/Science style
  labs(x = "Increasing landscape fragmentation ->", y = "Occupancy", title = "Occupancy in forest points across landscape (P = 0.24)") +
  theme(
    text = element_text(size = 16, family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

#plot terra raster plots #####
#reorder raster to be in increasing stages of fragmentation for plotting purposes 

file_names <- basename(sources(landscapes))
numeric_values <- as.numeric(gsub(".*_(\\d+)\\.tif", "\\1", file_names))
sorted_indices <- order(numeric_values)
landscapes<- landscapes[[sorted_indices]]


# Create "Figures" folder if it doesn't exist
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

# Extract and clean landscape names
landscape_names <- sources(landscapes)  # Extract file names
landscape_names <- basename(landscape_names)  # Remove file paths
landscape_names <- gsub("\\.tif$", "", landscape_names)  # Remove ".tif"

# Define number of rows and columns
num_cols <- 6
num_rows <- 1

# Define output file path
output_path <- file.path("Figures", "All_Landscapes_024.png")

# Open PNG graphics device
png(output_path, width = 2500, height = 1000, res = 200)  # Adjusted for 5x2 layout

# Set up multi-panel layout
par(mfrow = c(num_rows,num_cols),  # 2 rows, 5 columns
    mar = c(1, 1, 2, 1),  # Small margins for space
    oma = c(4, 4, 4, 4),  # Outer margins for caption
    bg = "white")  # White background

# Loop through the first 10 layers (to fit the 5x2 grid)
for (i in 1:min( nlyr(landscapes))) {
  plot(landscapes[[i]], 
       col = c( "darkgreen","tan"),  # Custom colors
       axes = FALSE, 
       box = FALSE, 
       legend = FALSE,
       main = landscape_names[i])  # Add title
}

# Close the graphics device
dev.off()

message("All landscapes saved as a single figure in 'Figures' folder.")

#Export figs #####

ggsave("Figures/Occupancy_Boxplot_p024.png", plot = p024, width = 10, height = 6, dpi = 300)
