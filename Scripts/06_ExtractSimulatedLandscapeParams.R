#Extract covariates from simulated rasters. To fit with our models we need: 
#1. Distance to coast ; scaleCoastDist          
#2. Habitat amount 100 and 2km;  scaleHabAmount100 and scaleHabAmount2000                    
#3. Edge density 100 and 2km ; scaleEdgeDens100 and scaleEdgeDens2000    

#mutate(edgeArea100 = (edgeRook_100_40 * (pi * 100^2)) / (habAmountDich_100 * (pi * 100^2) + 1)) %>% 
#mutate(edgeArea2000 = (edgeRook_2000_40 * (pi * 2000^2)) / (habAmountDich_2000 * (pi * 2000^2) + 1)) %>%  

#NOTE; when we read in the simulated rasters,rasters are based on number of cells rather than 
#area per se. So I can treat each cell as a 100X100 cell 
library(terra)
library(sf)
library(exactextractr)
library(data.table)
library(landscapemetrics)
# library(future.apply)  # For parallel processing
library(unmarked)
library(tidyverse)
# library(purrr)
library(ggplot2)


#read in final model
model <- readRDS("Models/pc1_interaction_model.rds")

# Define the folder where the TIFF files are stored
input_folder <- "Rasters/production_0.55"

# List all TIFF files in the folder (recursive = FALSE to avoid subfolders)
tif_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)

# Read multiple rasters at once into a single SpatRaster object
landscapes<- rast(tif_files)
plot(landscapes)
plot(landscapes[[1]])

#define key paramas####
buffer_sizes <- c(1, 20)  # Buffer sizes in meters
buffer_distances <- c(1, 20)  # Buffer sizes in meters

#--------------------------------------
#Develop sampling points across landscape #####
#--------------------------------------

test <- landscapes[[7]]
plot(test)

#extract the centre of each 100m cell 


#subset stratified points across the rasters

# Define the number of total points to stratify
n_points <- 100

# Calculate the spacing for rows and columns (1000 rows, 1000 columns)
rows <- seq(1, nrow(test), length.out = sqrt(n_points))
cols <- seq(1, ncol(test), length.out = sqrt(n_points))
# Create grid of row and column indices
grid_indices <- expand.grid(row = rows, col = cols)
# Convert to spatial coordinates
cell_centers <- xyFromCell(test, (grid_indices$row - 1) * ncol(test) + grid_indices$col)
# Create SpatVector of points
points <- vect(cell_centers, type = "points", crs = crs(test))
# Plot the raster and stratified points
plot(test)  # Plot the raster
plot(points, add = TRUE, col = "red", pch = 16)  # Overlay the stratified points


# remove points on the edge of the raster so that buffers are inside the boundary
# Set the edge buffer distance
edge_buffer <- 20  # 20 cells from the edges


# Extract coordinates from the SpatVector
coords <- crds(points)  # Extract coordinates (returns a matrix with x and y)

# Apply the filtering conditions
valid_coords <- coords[
  coords[, 1] > edge_buffer & coords[, 1] <= (1000 - edge_buffer) &  # x filter
    coords[, 2] > edge_buffer & coords[, 2] <= (1000 - edge_buffer),   # y filter
]

# Convert the filtered coordinates back into a SpatVector
points <- vect(valid_coords, type = "points", crs = crs(points))

# Add IDs to the filtered points
points$id <- paste0(1:nrow(points))

# Plot the raster and the filtered valid points
plot(test)  # Plot the raster
plot(points, add = TRUE, col = "red", pch = 16)  # Overlay the filtered points






#--------------------------------------------------
#Calculate habitat amount functions ####
#--------------------------------------------------

#1000*1000 = 1,000,000 
#so each square is 1ha (100x100) to make 1Mha
#EXTRACTION TO BUFFERS ####
#vector based extraction of buffer data
# Function to process a single raster and calculate habitat amount
process_extraction <- function(raster, buffer_distances, points) {
  habAmount <- list()
  
  # Ensure points have unique IDs
  points$id <- paste0("p", 1:nrow(points))
  
  for (buffer_distance in buffer_distances) {
    # Create buffers around each point and retain IDs
    buffers <- buffer(points, width = buffer_distance)
    buffers$id <- points$id  # Retain point IDs in the buffer object
    
    # Convert buffers to sf format
    buffers_sf <- st_as_sf(buffers)
    
    # Create an extent polygon from the raster
    raster_ext <- st_as_sf(st_as_sfc(st_bbox(raster)))
    
    # Clip the buffers to the raster extent
    buffers_clipped <- st_intersection(buffers_sf, raster_ext)
    
    # Perform extraction - 'frac' calculates the fraction of each value
    extraction <- exact_extract(raster, buffers_clipped, 'frac', progress = TRUE)
    extraction$id <- buffers$id  # Add ID back in 
    
    # Store results in the list
    habAmount[[paste0("buffer_", buffer_distance)]] <- extraction
  }
  
  # Collapse the nested list into a long data frame
 habAmountlong <-  rbindlist(habAmount, idcol = "buffer_size", fill = TRUE) %>% 
   select(buffer_size, frac_0, frac_1, id)
  
  return(habAmountlong)
}

#all_habAmount <- process_extraction(raster = test, buffer_distances = buffer_distances, points = points) 


# Process habitat amount for multiple rasters and store outputs in a list
process_all_landscapes <- function(rasters, buffer_distances, points) {
  results_list <- list()
  
  # Loop through each layer in the SpatRaster object
  for (i in 1:nlyr(rasters)) {
    cat("Processing raster", i, "of", nlyr(rasters), "\n")
    
    raster <- rasters[[i]]  # Extract individual raster layer
    
    #get raster name based on file name 
    raster_name <- sources(raster) %>% basename() %>% tools::file_path_sans_ext()  
    
    # Apply the function and store results
    habAmount <- process_extraction(raster = raster, buffer_distances = buffer_distances, points)
    
    results_list[[raster_name]] <- habAmount
  }
  
  return(results_list)
}

# Example Usage
landscapes # Assuming `landscape` is a SpatRaster object with multiple layers

# Process all rasters
all_habAmount <- process_all_landscapes(rasters = landscapes, buffer_distances = buffer_distances, points = points)

#Save habitat amounts outputs #####
saveRDS(all_habAmount, "Outputs/HabAmount10000pts_Simulated_055_LandscapeParams.rds")

#--------------------------------------------------
#CALCULATE proportional EDGE DENSITY ####
#NB - this method calculates the number of cells that are forest 
#and that have a hard edge as a proportion of a cells in a buffer  
#--------------------------------------------------


# Define a function to compute forest edge cells for a single layer
forest_edge_cells <- function(raster_layer) {
  
  # Create a mask for forest cells (value == 1)
  forest_mask <- raster_layer == 1
  
  # Identify cells that share a border with any different class in NSEW direction  
  class_difference <- boundaries(raster_layer, inner = TRUE, directions = 4, classes = TRUE)
  
  # Compute forest edge cells (forest cells sharing a hard edge with another class)
  forest_edges <- class_difference * forest_mask
  
  return(forest_edges)
}

# Create an empty list to store forest edge rasters for each layer
forest_edge_list <- list()

# Loop through each layer in the SpatRaster
for (i in 1:nlyr(landscapes)) {
  # Apply the function to each layer
  forest_edge_layer <- forest_edge_cells(landscapes[[i]])
  # Add the result to the list
  forest_edge_list[[i]] <- forest_edge_layer
}

# Combine the list of layers back into a SpatRaster
landscape_forest_edges <- rast(forest_edge_list)

# Copy layer names from the original landscapes SpatRaster to the new SpatRaster
names(landscape_forest_edges) <- names(landscapes)
plot(landscape_forest_edges[[4]])

#now extract the fractional coverage of edge and non-edge cells for each buffer

process_extraction_edge <- function(raster, buffer_distances, points) {
  edgeAmount <- list()
  
  # Ensure points have unique IDs
  points$id <- paste0("p", 1:nrow(points))
  
  for (buffer_distance in buffer_distances) {
    # Create buffers around each point and retain IDs
    buffers <- buffer(points, width = buffer_distance)
    buffers$id <- points$id  # Retain point IDs in the buffer object
    
    # Convert buffers to sf format
    buffers_sf <- st_as_sf(buffers)
    
    # Create an extent polygon from the raster
    raster_ext <- st_as_sf(st_as_sfc(st_bbox(raster)))
    
    # Clip the buffers to the raster extent
    buffers_clipped <- st_intersection(buffers_sf, raster_ext)
    
    # Perform extraction - 'frac' calculates the fraction of each value
    extraction <- exact_extract(raster, buffers_clipped, 'frac', progress = TRUE)
    extraction$id <- buffers$id  # Add ID back in 
    
    # Store results in the list
    habAmount[[paste0("buffer_", buffer_distance)]] <- extraction
  }
  
  # Collapse the nested list into a long data frame
  edgeAmount <-  rbindlist(habAmount, idcol = "buffer_size", fill = TRUE) %>% 
    select(buffer_size, frac_0, frac_1, id)
  
  return(edgeAmount)
}


#  Usage
landscape_forest_edges #

# Process all rasters
all_edges <- process_all_landscapes(rasters = landscape_forest_edges, buffer_distances = buffer_distances, points = points)

#---------------------------------------------------


#Quick plot of distribution of edge density values for different landscapes & buffers ####


# Group the data by landscape_name and buffer_size, and create a histogram of 'value'
edge_density_df %>%
  group_by(landscape_name, buffer_size) %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 500, fill = "blue", color = "black", alpha = 0.7) +  # Adjust binwidth as needed
  facet_wrap(~landscape_name + buffer_size, scales = "free_y") +  # Facet by landscape_name and buffer_size
  labs(title = "Histogram of Values Grouped by Landscape and Buffer Size",
       x = "Value",
       y = "Frequency") +
  theme_minimal()


#--------------------------------------------------
#CALCULATE LINEAR EDGE DENSITY ####
#NB - this method calculates the linear length of hard edges using landscape metrics 
#--------------------------------------------------

#Description of edge density in landscapemetrics package:
#The formula for Edge Density (ED) is:
# (𝐸/𝐷) *10,000
#where𝐸 = total edge length in the landscape (measured in meters)
#𝐴 = total landscape area (measured in square meters)
#The factor 10,000 is used to convert the result into meters per hectare (m/ha), as 1 hectare = 10,000 m².

# 
# # Assuming 'landscapes' is a single SpatRaster with multiple layers
# # and 'points' is a SpatVector of points
# 
# # Function to calculate edge density for a single point with a specific buffer size
# calculate_edge_density <- function(point, raster, buffer_size) {
#   # Step 1: Create buffer around the point
#   buffer <- buffer(point, width = buffer_size)  # Create buffer around the point
#   
#   # Step 2: Mask the raster with the buffer
#   masked_raster <- mask(raster, buffer)  # Mask raster by the buffer
#   
#   # Step 3: Calculate edge density within the masked raster
#   edge_density <- lsm_l_ed(masked_raster, count_boundary = FALSE, directions = 4)  # rook's edge density
#   
#   # Extract point ID from the point object
#   point_id <- point$id
#   
#   # Add buffer size and point information to the result
#   edge_density$buffer_size <- buffer_size
#   edge_density$point_id <- point_id  # Assign point ID to the result
#   
#   return(edge_density)
# }
# 
# # Preallocate the list to store results (number of points * number of buffer sizes * number of raster layers)
# num_points <- nrow(points)
# num_buffer_sizes <- length(buffer_sizes)
# total_results <- num_points * num_buffer_sizes * nlyr(landscapes)
# edge_density_results <- vector("list", total_results)
# 
# # Initialize counter for results list
# result_index <- 1
# 
# start_time <- Sys.time()
# # Initialize a list to store dataframes for each raster layer
# edge_density_list <- list()
# 
# # Start the clock to track processing time
# start_time <- Sys.time()
# 
# # Loop through each layer of the landscapes (SpatRaster with multiple layers)
# for (layer_index in 1:nlyr(landscapes)) {
#   # Extract the current raster layer
#   raster_layer <- landscapes[[layer_index]]  # Extract current layer
#   
#   # Get the name of the raster layer from its source 
#   raster_name <- sources(raster_layer) %>% basename() %>% tools::file_path_sans_ext()    
#  
#   message("Processing layer ", layer_index, " (", raster_name, ") of ", nlyr(landscapes))
#   
#   # Initialize a list to store edge density results for this layer
#   edge_density_results <- list()
#   
#   # Loop through each point and buffer size to calculate edge density
#   for (point in 1:num_points) {
#     point_sf <- points[point, ]  # Extract individual point
#     cat("Processing point", point, "of", num_points, "\n")  # Print progress for points
#     for (buffer_size in buffer_sizes) {
#       # Calculate edge density for the point with the given buffer size and current layer
#       result <- calculate_edge_density(point_sf, raster_layer, buffer_size)  # No point_id argument here
#       
#       # Store the result in the list for this raster layer
#       edge_density_results[[length(edge_density_results) + 1]] <- result
#     }
#   }
#   
#   # Convert the list of results into a data frame
#   edge_density_df <- do.call(rbind, edge_density_results)
#   
#   # Assign the dataframe to the list using the raster name as the key
#   edge_density_list[[raster_name]] <- edge_density_df
#   
# }
# 
# end_time <- Sys.time()
# message("Edge density calculation completed in ", end_time - start_time, " seconds.")
# 
# edge_density_df <- rbindlist(edge_density_list, idcol = "landscape_name" )
# # Save the list of dataframes to an RDS file
# saveRDS(edge_density_list, "Outputs/EdgeDensity10000pts_Simulated_055_LandscapeParams.rdsedge_density_list.rds")
