# Predict occupancy for hypoethetical landscapes with different configurations 

# Load necessary libraries
library(tidyverse)
library(data.table)
library(unmarked)
library(lme4)
library(terra)
library(raster)
library(RColorBrewer)
library(stringr)
library(Matrix)

source("scripts/01_OceanPredictors.R")
source("scripts/02_OrganiseMurreletData.R")

# Define the production targets
production_targets <- c("0.24", "0.58")
model <- readRDS("Models/final_model_5thMay2025.rds")

# Model, covariate, and site data inputs
covariates <- readRDS("Outputs/ScaledCovariates.rds") %>%   
  dplyr::select(PC1_t1, scaleCoastDist, scaleDoy, scaleDoy2, scaleDoy2, OceanYear) %>%  unique()
real_murrelet_site_data <- read.csv("Inputs/siteData.csv")

#read in means and SDs of covariates from original model to enable coorrect scalings
meansAndSds

# #determine which distance to coast you want to make the simlations for #####
coast_mean <- mean(siteCovs(analysisData)$coastDist, na.rm = TRUE)
coast_sd <- sd(siteCovs(analysisData)$coastDist, na.rm = TRUE)
current_dist <- covariates$scaleCoastDist %>% unique() 
# Reverse scaling 
original_value <- (current_dist * coast_sd) + coast_mean
print(original_value)

# #add different distances to the coast
# covariates <- covariates %>%  crossing(coast_levels)

# Loop through production targets
for (target in production_targets) {
  
  # Read in simulated inputs for the correct production target
  hab <- readRDS(paste0("Outputs/production_", target, "_HabAmount10000pts_SimulatedParams.rds"))
  edge <- readRDS(paste0("Outputs/production_", target, "_EdgeDensity10000pts_SimulatedParams.rds"))
  pt_value <- readRDS(paste0("Outputs/production_", target, "_point_forest_or_plantation.rds")) %>% 
    as.data.frame() %>%  
    dplyr::select(landscape_name, point_id, raster_value.lyr.1) %>%  
    rename(forest_plantation = raster_value.lyr.1) %>%  
    mutate(point_id = paste0("p", point_id))
  
  # Define the folder where the landscape TIFF files are stored
  input_folder <- paste0("Rasters/production_", target)
  
  # Get a list of .tif files
  tif_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)
  
  # Read the rasters into a SpatRaster object
  landscapes <- rast(tif_files)
  plot(landscapes)
  
  
  # Process habitat amount data
  hab_df <- hab %>% rbindlist(idcol = "landscape_name")
  hab_df <- hab_df %>% 
    pivot_wider(
      names_from = buffer_size, 
      values_from = c(frac_0, frac_1)
    ) %>%    
    mutate(habAmountDich_100 = frac_0_buffer_1, 
           habAmountDich_2000 = frac_0_buffer_20,
           point_id = id) %>% 
    dplyr::select(landscape_name, point_id, habAmountDich_100, habAmountDich_2000)
  
  # Process edge density data
  edge_df <- edge %>%  rbindlist(idcol = "landscape_name")
  edge_df <- edge_df %>% 
    pivot_wider(
      names_from = buffer_size, 
      values_from = c(frac_0, frac_1)
    ) %>%    
    mutate(edgeRook_100_40 = frac_1_buffer_1, 
           edgeRook_2000_40 = frac_1_buffer_20,
           point_id = id) %>% 
    dplyr::select(landscape_name, point_id, edgeRook_100_40, edgeRook_2000_40)
  
  # Combine habitat and edge data
  all_df <- hab_df %>% left_join(edge_df)
  
  
  # Prepare data for prediction with scaled covariates
  prediction_df <- all_df %>% mutate(
    scaleHabAmount100 = (habAmountDich_100  - meansAndSds$meanHabAmountDich100) /meansAndSds$sdHabAmountDich100,
    scaleHabAmount2000 = (habAmountDich_2000 - meansAndSds$meanHabAmountDich2000) /meansAndSds$sdHabAmountDich2000,
    scaleEdgeDens100 =  (edgeRook_100_40 - meansAndSds$meanEdgeDens100) /meansAndSds$sdEdgeRook100, 
    scaleEdgeDens2000 = (edgeRook_2000_40 - meansAndSds$meanEdgeDens2000)/meansAndSds$sdEdgeRook2000) %>% 
   dplyr::select(landscape_name, point_id, scaleHabAmount100, scaleHabAmount2000, scaleEdgeDens100, scaleEdgeDens2000) %>% 
    cross_join(covariates)
  
  # Predict occupancy with standard errors (for error ribbon)
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
    left_join(pt_value) %>%  
    #filter only points that are actually forest 
    filter(forest_plantation ==0) %>% 
    dplyr::select(landscape_name, point_id, Occupancy, OceanYear, lower_CI, upr_CI) %>% 
    mutate(landscape_numeric = as.numeric(gsub("patches_", "", landscape_name))) %>% 
    arrange(landscape_numeric) %>%  # Sort from smallest to largest  
    mutate(landscape_numeric= as.factor(landscape_numeric)) %>% 
    #mutate(landscape_name = fct_reorder(landscape_name)) 
    mutate(OceanYear = factor(OceanYear, 
                            levels = c("Bad Ocean Years", "Good Ocean Years")))
  
  # # Rapid plots
  # p_plot <- plot_data %>%  
  #   ggplot( aes(x = landscape_name, y = Occupancy)) +
  #   geom_jitter(color = 'lightgrey', size = 2, width = 0.1, alpha = 0.05) +  
  #   geom_boxplot(fill = "#56B4E9", color = "black", width = 0.5, outlier.shape = NA, outlier.size = NA) + 
  #   theme_classic(base_size = 14) +
  #   facet_wrap(~OceanYear) +
  #   labs(x = "Increasing landscape fragmentation ->", y = "Occupancy", title = paste("Occupancy in forest points across landscape (P =", target, ")")) +
  #   theme(
  #     text = element_text(size = 16, family = "serif"),
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     axis.title = element_text(size = 14),
  #     plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
  #     legend.position = "none"
  #   )
  
  p_plot <- plot_data %>%  
    ggplot(aes(x = landscape_numeric, y = Occupancy)) +
    # Plot the points with jitterdodge (so they stay behind the boxplot)
    #geom_point(aes(color = OceanYear), size = 2, alpha = 0.05, position = position_jitterdodge(jitter.width = 0.1)) +  # Points with jitter
    # Plot the boxplots
    geom_boxplot(aes(fill = OceanYear), color = "black", width = 0.5, outlier.shape = NA) + 
    # Custom fill colors for OceanYear
    scale_fill_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + 
    scale_color_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + # Color points
 #   facet_wrap(~CoastDist, ncol = 4) +  # Facet by Coastal Distance
    theme_minimal(base_size = 18) +  # Minimal theme with larger base font size
    ylim(0,1)+
    labs(x = "Number of patches (degree of fragmentation per se)", y = "Occupancy") +
    theme(
        legend.position = "top",  # Position the legend at the top
        legend.title = element_blank(),  # Increase legend title size for clarity
        legend.text = element_text(size = 14),  # Increase legend text size for better readability
        axis.title = element_text(size = 16),  # Increase axis title font size
        axis.text = element_text(size = 14),  # Increase axis label font size
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # No major gridlines
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a thin border around the plot
      )
    
  
  
  # # Save occupancy boxplot
  # ggsave(paste0("Figures/ForestPointOccupancyInSimulatedLandscapes_P", target, ".png"),
  #        plot = p_plot, 
  #        width = 10,                            # Width in inches (publication size)
  #        height = 8,                           # Height in inches (publication size)
  #        dpi = 300,  
  #        bg = "white")  # White background)
  # 
  # Save occupancy boxplot as PDF
  ggsave(
    filename = paste0("Figures/ForestPointOccupancyInSimulatedLandscapes_P", target, ".pdf"),
    plot = p_plot, 
    width = 10,           # Width in inches
    height = 8,           # Height in inches
    bg = "white"          # Optional: white background (default for PDF is transparent)
  )
  
  # Export terra raster plots
  file_names <- basename(sources(landscapes))
  numeric_values <- as.numeric(gsub(".*_(\\d+)\\.tif", "\\1", file_names))
  sorted_indices <- order(numeric_values)
  landscapes <- landscapes[[sorted_indices]]
  
  # Create "Figures" folder if it doesn't exist
  if (!dir.exists("Figures")) {
    dir.create("Figures")
  }
  
  # Extract and clean landscape names
  landscape_names <- sources(landscapes)  # Extract file names
  landscape_names <- basename(landscape_names)  # Remove file paths
  landscape_names <- gsub("\\.tif$", "", landscape_names)  # Remove ".tif"
  
  # Define number of rows and columns
  num_cols <- 4
  num_rows <- 4
  #PNG
  # # Define output file path
  # output_path <- file.path("Figures", paste0("All_Landscapes_", target, ".png"))
  # 
  # # Open PNG graphics device
  # png(output_path, width = 2500, height = 1000, res = 200)
  
  # Define output file path
  output_path <- file.path("Figures", paste0("All_Landscapes_", target, ".pdf"))
  
  # Open PDF graphics device
  pdf(output_path, width = 12.5, height = 5)  # widths in inches
  
  
  # Set up multi-panel layout
  par(mfrow = c(num_rows, num_cols), 
      mar = c(1, 1, 2, 1), 
      oma = c(4, 4, 4, 4), 
      bg = "white")  # White background
  
  # Loop through the layers in the landscape
  for (i in 1:nlyr(landscapes)) {
    plot(landscapes[[i]], 
         col = c("darkgreen", "tan"), 
         axes = FALSE, 
         box = FALSE, 
         legend = FALSE, 
         main = landscape_names[i])
  }
  
  # Close the graphics device
  dev.off()
  
  message(paste("All landscapes for production target", target, "saved as a single figure in 'Figures' folder."))
}
