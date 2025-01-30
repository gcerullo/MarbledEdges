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
plot(points, add=TRUE, col="red",  cex=0.0000000000000000001)
# Assign a unique ID to each point
points$id <- 1:nrow(points)

#1000*1000 = 1,000,000 
#so each square is 1ha (100x100) to make 1Mha
#EXTRACTION TO BUFFERS ####
#vector based extraction of buffer data

#1Mha landscape; so each square is 1ha
# Initialize lists to store results
results_habAmount <- list()

buffer_distances <- c(1,20)

#TO CORRECT!!!!!!!!!!!!!!!
#The 1million cell raster is currently at 1m2 resolution; actually each cell should correspond to a ha!!!
#Come back and make sure tihs is correct before extracting the buffer!!!!!!!! 

raster <- test
#crs(raster) <- "+proj=utm +zone=11 +datum=WGS84"
# Function to process extraction
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
    
    # Add the ID of the sample point to the extraction
    extraction_with_id <- cbind(extraction, id = buffers_clipped$id)
    
    # Store results with buffer distance information
    results[[paste0("buffer_", buffer_distance)]] <- extraction_with_id
  }
  return(results)
}

# # Process results
results_test <- process_extraction(test, buffer_distances)

#--------------------------------------------------
#TO DO; SUMMARISE HABITAT AMOUNT PER CELL
#--------------------------------------------------

# results_test[[100]]
# 
# #Capture cells within each  ####
# 
# # Function to summarize LIDAR extractions
# summarize_extractions <- function(results_LiDAR) {
#   
#   # Use map to apply a function to each element of the results_LiDAR list
#   map(results_LiDAR, ~{
#     
#     # Access the extractions and station info within each buffer list
#     results_test <- .x$extraction
#     station_info <- .x$station_info
#     
#     # Use map2 to iterate over extractions and station info simultaneously
#     map2(extractions, station_info, ~{
#       # Calculate summary statistics for each extraction DataFrame
#       df <- as_tibble(.x)
#       summarise(df, 
#                 Mean = mean(value, na.rm = TRUE), 
#                 SD = sd(value, na.rm = TRUE), 
#                 Median = median(value, na.rm = TRUE)) %>%
#         # Add station info as a new column for reference
#         mutate(StationID_Year = .y)
#     }) %>%
#       # Optional: Combine all mini summaries into a single tibble
#       bind_rows()
#   })
# }
# 
# # Process the results_LiDAR data
# summarized_lidar_data <- summarize_lidar_extractions(results_LiDAR) %>%  
#   #add in information on buffer size
#   rbindlist(idcol = "buffer_size") %>% 
#   #make buffer size column numeric
#   mutate(buffer_size = str_extract(buffer_size, "\\d+") %>% as.numeric()) %>% 
#   #bring buffer size back to metres 
#   mutate(buffer_size = buffer_size/metre_to_ft)

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

#--------------------------------------------------
#CALCULATE EDGE DENSITY PARRALELISED 
#--------------------------------------------------

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)  # Use all but one core

# Example buffer sizes (adjust based on your analysis)
buffer_sizes <- c(1, 20)  # Buffer sizes in meters

# Function to calculate edge density for a single point with a specific buffer size
calculate_edge_density <- function(point, raster, buffer_size) {
  # Step 1: Create buffer around the point
  buffer <- buffer(point, width = buffer_size)  # Create buffer around the point
  
  # Step 2: Mask the raster with the buffer
  masked_raster <- mask(raster, buffer)  # Mask raster by the buffer
  
  # Step 3: Calculate edge density within the masked raster
  edge_density <- lsm_l_ed(masked_raster, count_boundary = FALSE, directions = 4)  # DIRECTIONS = rook's edge density
  
  # Add buffer size and point information to the result
  edge_density$buffer_size <- buffer_size
  edge_density$point_id <- point$id
  
  return(edge_density)
}

# Preallocate a data.table for storing results (more efficient for large datasets)
edge_density_dt <- data.table(point_id = integer(0), buffer_size = integer(0), edge_density = numeric(0))

# Function to calculate edge density for all buffer sizes for a given point
calculate_for_point <- function(point, raster, buffer_sizes) {
  results <- list()
  for (buffer_size in buffer_sizes) {
    result <- calculate_edge_density(point, raster, buffer_size)
    # Add the result to the list
    results[[length(results) + 1]] <- list(
      point_id = point$id,
      buffer_size = buffer_size,
      edge_density = result$value  # Assuming 'value' holds the edge density metric
    )
  }
  return(results)
}

# Use future.apply to process points in parallel
result_list <- future_lapply(1:nrow(points), function(point_index) {
  point <- points[point_index, ]  # Extract individual point
  calculate_for_point(point, raster, buffer_sizes)
})

# Flatten the result list and convert it to a data.table
edge_density_dt <- rbindlist(unlist(result_list, recursive = TRUE), use.names = TRUE)

# Optionally, you can print the final progress message when everything is done
message("Edge density calculation completed for all points and buffer sizes.")

# Print first few rows of results
print(head(edge_density_dt))

#--------------------------------------------------
#VISUALISE WHAT'S HAPPENING IN THE CALCULATION of edge density 
#--------------------------------------------------
# Select one point from your points SpatVector and one buffer size
points
test

point <- points[40000, ]  # Just take the first point for simplicity
buffer_size <- 20 # Buffer size of 100 meters

# Step 1: Create the buffer around the point
buffer <- buffer(point, width = buffer_size)  # Create buffer using terra::buffer
#add 

# Step 2: Visualize the original point and the buffer
# Convert the buffer to an sf object for ggplot visualization
buffer_sf <- st_as_sf(buffer)

# Convert the point to an sf object for visualization
point_sf <- st_as_sf(point)

# Plot the point and buffer
ggplot() +
  geom_sf(data = buffer_sf, fill = "red", color = "blue", alpha = 0.5) +  # Buffer area
  geom_sf(data = point_sf, color = "red", size = 3) +  # Point location
  theme_minimal() +
  labs(title = "Point and Buffer Visualization", subtitle = paste("Buffer size:", buffer_size, "meters")) +
  coord_sf()

# Step 3: Mask the raster by the buffer
masked_raster <- mask(raster, buffer)

# Plot the original raster and the masked area
plot(raster, main = "Original Raster", col = terrain.colors(10))  # Plot original raster
plot(masked_raster, main = "Masked Raster (with Buffer)", col = terrain.colors(10))  # Plot masked raster

# Step 4: Calculate edge density within the masked area
edge_density <- lsm_l_ed(masked_raster, count_boundary = FALSE, directions = 4)

# Print edge density result
print(edge_density)


