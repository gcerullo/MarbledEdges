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
library(future.apply)  # For parallel processing
library(unmarked)
library(tidyverse)
library(purrr)

#read in final model
model <- readRDS("Models/pc1_interaction_model.rds")

# Define the folder where the TIFF files are stored
input_folder <- "Rasters/production_0.55"

# List all TIFF files in the folder (recursive = FALSE to avoid subfolders)
tif_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)

# Read multiple rasters at once into a single SpatRaster object
landscape <- rast(tif_files)
plot(landscape[[1]])

#----------------------------
#CALCULATE HABITAT AMOUNT 
#----------------------------
test <- landscape[[10]]
plot(test)
#extract the centre of each 100m cell 

cell_centers <- xyFromCell(test, 1:ncell(test))
points <- vect(cell_centers, type="points")

#test points 
#points <- points[500050:500150, ] #test function with few points
plot(points, add=TRUE, col="red",  cex=0.0000000000000000001)
plot(points, add=TRUE, col="red")

# Assign a unique ID to each point
points$id <- paste0("p", 1:nrow(points))


#1000*1000 = 1,000,000 
#so each square is 1ha (100x100) to make 1Mha
#EXTRACTION TO BUFFERS ####
#vector based extraction of buffer data

#1Mha landscape; so each square is 1ha
# Initialize lists to store results
buffer_distances <- c(1,20)

raster <- test
process_extraction <- function(raster, buffer_distances) {
  results <- list()
  
  for (buffer_distance in buffer_distances) {
    # Create buffers around each point
    buffers <- buffer(points, width = buffer_distance)
    
    # Convert spatVector buffers to sf format
    buffers_sf <- st_as_sf(buffers)
    
    # Create an extent polygon from the raster
    raster_ext <- st_as_sf(st_as_sfc(st_bbox(raster)))
    
    # Clip the buffers_sf to the raster extent
    buffers_clipped <- st_intersection(buffers_sf, raster_ext)
    
    # Perform extraction
    extraction <- exact_extract(raster, buffers_clipped, progress = TRUE)
    
    # Loop through each extraction and assign corresponding ID
    extraction_with_id <- lapply(1:length(extraction), function(i) {
      # Add the ID to each individual extraction dataframe
      extraction_df <- extraction[[i]]
      extraction_df$id <- buffers_clipped$id[i]  # Add the id for this buffer 
     
    #   # #apply the summary of buffer habitat coverage
    #   extraction_df <- extraction_df %>%
    #     #get proportional coverage of each intersecting cell of value 1 or 0 (forest/plantation)
    #     group_by(id, value) %>%
    #     summarise(
    #       #how many cells of 0/1 does the buffer intersect?
    #       sum_class_cells = sum(value),
    #       #what's the sum of fractional coverage of cells (e.g. does the buffer cover all of the cell (=1), summed
    #       total_class_coverage_fraction = sum(coverage_fraction),
    #       .groups = "drop"  # Ensures grouping does not persist
    #     )  
    #    return(extraction_df)  # Return modified extraction dataframe with id
    # })
      
      # **Use data.table for fast summarization**
      setDT(extraction_df)[, .(
        sum_class_cells = sum(value),
        total_class_coverage_fraction = sum(coverage_fraction)
      ), by = .(id, value)]
    })
    
    # Store results in a list where each point keeps its own dataframe
    results[[paste0("buffer_", buffer_distance)]] <- extraction_with_id
  }
  
  
  return(results)
}


# Example Usage
results_test <- process_extraction(raster = test, buffer_distances = c(1, 20))

#--------------------------------------------------
#CALCULATE EDGE DENSITY ####
#--------------------------------------------------
library(terra)  # For raster and vector manipulation
library(landscapemetrics)  # For landscape metrics (including edge density)

points
raster
# Assuming points (SpatVector) and raster (SpatRaster) are already loaded
# Example buffer sizes (adjust based on your analysis)
buffer_sizes <- c(1, 20)  # Buffer sizes in meters

# Function to calculate edge density for a single point with a specific buffer size
calculate_edge_density <- function(point, raster, buffer_size) {
  # Step 1: Create buffer around the point
  buffer <- buffer(point, width = buffer_size)  # Create buffer around the point
  
  # Step 2: Mask the raster with the buffer
  masked_raster <- mask(raster, buffer)  # Mask raster by the buffer
  
  # Step 3: Calculate edge density within the masked raster
  edge_density <- lsm_l_ed(masked_raster, count_boundary = FALSE, directions = 4)  # rook's edge density
  
  # Add buffer size and point information to the result
  edge_density$buffer_size <- buffer_size
  edge_density$point_id <- point$id
  
  return(edge_density)
}

# Preallocate the list to store results (number of points * number of buffer sizes)
num_points <- nrow(points)
num_buffer_sizes <- length(buffer_sizes)
total_results <- num_points * num_buffer_sizes
edge_density_results <- vector("list", total_results)

# Initialize counter for results list
result_index <- 1

# Loop through each point and buffer size to calculate edge density
for (point in 1:num_points) {
  point_sf <- points[point, ]  # Extract individual point
  
  for (buffer_size in buffer_sizes) {
    # Calculate edge density for the point with the given buffer size
    result <- calculate_edge_density(point_sf, raster, buffer_size)
    
    # Store the result in the preallocated list
    edge_density_results[[result_index]] <- result
    
    # Increment result index
    result_index <- result_index + 1
    
    # Print progress for each point and buffer size
    message("Processed point ", point, " with buffer size ", buffer_size, " meters.")
  }
}

# Optionally, you can print the final progress message when everything is done
message("Edge density calculation completed for all points and buffer sizes.")


# Combine the results into a single data frame or tibble
edge_density_df <- do.call(rbind, edge_density_results)

# View the results
print(edge_density_df)
# 
# #--------------------------------------------------
# #CALCULATE EDGE DENSITY PARRALELISED 
# #--------------------------------------------------
# 
# # Set up parallel processing
# plan(multisession, workers = parallel::detectCores() - 1)  # Use all but one core
# 
# # Example buffer sizes (adjust based on your analysis)
# buffer_sizes <- c(1, 20)  # Buffer sizes in meters
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
#   edge_density <- lsm_l_ed(masked_raster, count_boundary = FALSE, directions = 4)  # DIRECTIONS = rook's edge density
#   
#   # Add buffer size and point information to the result
#   edge_density$buffer_size <- buffer_size
#   edge_density$point_id <- point$id
#   
#   return(edge_density)
# }
# 
# # Preallocate a data.table for storing results (more efficient for large datasets)
# edge_density_dt <- data.table(point_id = integer(0), buffer_size = integer(0), edge_density = numeric(0))
# 
# # Function to calculate edge density for all buffer sizes for a given point
# calculate_for_point <- function(point, raster, buffer_sizes) {
#   results <- list()
#   for (buffer_size in buffer_sizes) {
#     result <- calculate_edge_density(point, raster, buffer_size)
#     # Add the result to the list
#     results[[length(results) + 1]] <- list(
#       point_id = point$id,
#       buffer_size = buffer_size,
#       edge_density = result$value  # Assuming 'value' holds the edge density metric
#     )
#   }
#   return(results)
# }
# 
# # Use future.apply to process points in parallel
# result_list <- future_lapply(1:nrow(points), function(point_index) {
#   point <- points[point_index, ]  # Extract individual point
#   calculate_for_point(point, raster, buffer_sizes)
# })
# 
# # Flatten the result list and convert it to a data.table
# edge_density_dt <- rbindlist(unlist(result_list, recursive = TRUE), use.names = TRUE)
# 
# # Optionally, you can print the final progress message when everything is done
# message("Edge density calculation completed for all points and buffer sizes.")
# 
# # Print first few rows of results
# print(head(edge_density_dt))
# 
# #--------------------------------------------------
# #VISUALISE WHAT'S HAPPENING IN THE CALCULATION of edge density 
# #--------------------------------------------------
# # Select one point from your points SpatVector and one buffer size
# points
# test
# 
# point <- points[40000, ]  # Just take the first point for simplicity
# buffer_size <- 20 # Buffer size of 100 meters
# 
# # Step 1: Create the buffer around the point
# buffer <- buffer(point, width = buffer_size)  # Create buffer using terra::buffer
# #add 
# 
# # Step 2: Visualize the original point and the buffer
# # Convert the buffer to an sf object for ggplot visualization
# buffer_sf <- st_as_sf(buffer)
# 
# # Convert the point to an sf object for visualization
# point_sf <- st_as_sf(point)
# 
# # Plot the point and buffer
# ggplot() +
#   geom_sf(data = buffer_sf, fill = "red", color = "blue", alpha = 0.5) +  # Buffer area
#   geom_sf(data = point_sf, color = "red", size = 3) +  # Point location
#   theme_minimal() +
#   labs(title = "Point and Buffer Visualization", subtitle = paste("Buffer size:", buffer_size, "meters")) +
#   coord_sf()
# 
# # Step 3: Mask the raster by the buffer
# masked_raster <- mask(raster, buffer)
# 
# # Plot the original raster and the masked area
# plot(raster, main = "Original Raster", col = terrain.colors(10))  # Plot original raster
# plot(masked_raster, main = "Masked Raster (with Buffer)", col = terrain.colors(10))  # Plot masked raster
# 
# # Step 4: Calculate edge density within the masked area
# edge_density <- lsm_l_ed(masked_raster, count_boundary = FALSE, directions = 4)
# 
# # Print edge density result
# print(edge_density)
# 
# 
