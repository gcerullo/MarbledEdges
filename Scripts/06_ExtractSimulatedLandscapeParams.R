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
landscapes<- rast(tif_files)
plot(landscapes[[1]])

#----------------------------
#CALCULATE HABITAT AMOUNT 
#----------------------------
test <- landscapes[[1]]
plot(test)
#extract the centre of each 100m cell 

cell_centers <- xyFromCell(test, 1:ncell(test))
points <- vect(cell_centers, type="points")

#test points 
points <- points[500050:500950, ] #test function with few points
plot(points, add=TRUE, col="red",  cex=0.0000000000000000001)
plot(points, add=TRUE, col="red")

# Assign a unique ID to each point
points$id <- paste0("p", 1:nrow(points))

#1000*1000 = 1,000,000 
#so each square is 1ha (100x100) to make 1Mha
#EXTRACTION TO BUFFERS ####
#vector based extraction of buffer data

#--------------------------------------------------
#Calculate habitat amount ####
#--------------------------------------------------
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

all_habAmount <- process_extraction(raster = test, buffer_distances = buffer_distances, points = points) 


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
buffer_distances <- c(1, 20)

# Process all rasters
all_habAmount <- process_all_landscapes(rasters = landscapes, buffer_distances = buffer_distances, points = points)

#Save habitat amounts outputs #####
saveRDS(all_habAmount, "Outputs/Simulated_055_LandscapeParams.rds")
#--------------------------------------------------
#CALCULATE EDGE DENSITY ####
#--------------------------------------------------

#Description of edge density in landscapemetrics package:
#The formula for Edge Density (ED) is:
# (𝐸/𝐷) *10,000
#where𝐸 = total edge length in the landscape (measured in meters)
#𝐴 = total landscape area (measured in square meters)
#The factor 10,000 is used to convert the result into meters per hectare (m/ha), as 1 hectare = 10,000 m².

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

start_time <- Sys.time()

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

end_time <- Sys.time()
message("Edge density calculation completed in ", end_time - start_time, " seconds.")

# Combine the results into a single data frame or tibble
edge_density_df <- do.call(rbind, edge_density_results)

# View the results
print(edge_density_df)
