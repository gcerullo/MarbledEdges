# Load necessary libraries
library(terra)
library(sf)
library(exactextractr)
library(data.table)
library(landscapemetrics)
library(future.apply)  # For parallel processing
library(purrr)  # For functional programming

# Set up parallel processing (use all but one core for efficiency)
plan(multisession, workers = parallel::detectCores() - 1)  # Adjust based on available cores

# Load model and raster files
model <- readRDS("Models/pc1_interaction_model.rds")
input_folder <- "Rasters/production_0.55"
tif_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)

# Read rasters into a SpatRaster object
landscape <- rast(tif_files)
plot(landscape[[1]])

# Function to create points with IDs from raster
create_points_with_ids <- function(raster) {
  cell_centers <- xyFromCell(raster, 1:ncell(raster))
  points <- vect(cell_centers, type = "points")
  points$id <- 1:nrow(points)
  return(points)
}

# Generate points based on raster
test_raster <- landscape[[10]]  # Use the 10th raster as an example
points <- create_points_with_ids(test_raster)

# Function for buffer extraction (to extract data within a buffer)
process_extraction <- function(raster, points, buffer_distances) {
  map(buffer_distances, ~{
    buffer <- buffer(points, width = .x)  # Create buffer
    buffers_sf <- st_as_sf(buffer)
    raster_ext <- st_as_sf(st_as_sfc(st_bbox(raster)))  # Get raster extent
    buffers_clipped <- st_intersection(buffers_sf, raster_ext)
    
    # Extract values within the clipped buffer area
    extraction <- exact_extract(raster, buffers_clipped, progress = TRUE)
    extraction_with_id <- cbind(extraction, id = buffers_clipped$id)
    
    # Return the result with the buffer distance
    list(buffer_distance = .x, extraction = extraction_with_id)
  })
}


# Process extraction results with multiple buffer sizes
buffer_distances <- c(1, 20)
results_test <- process_extraction(test_raster, points, buffer_distances)

#THIS SUMMARY IS CURRENTLY WRONG AND NEEDS CHANGING!!!!!!!!
# Summarize extractions (mean, SD, and median) 
summarize_extractions <- function(extraction_results) {
  map_dfr(extraction_results, ~{
    df <- as_tibble(.x$extraction)
    summarise(df, 
              Mean = mean(value, na.rm = TRUE), 
              SD = sd(value, na.rm = TRUE), 
              Median = median(value, na.rm = TRUE)) %>%
      mutate(buffer_size = .x$buffer_distance)
  })
}

# Summarize extraction results into a data.table
summarized_data <- summarize_extractions(results_test)

#----------------------------------------------------------------------------------
# Function to calculate edge density for a single point with a specific buffer size
#----------------------------------------------------------------------------------

# Function to calculate edge density for a point and buffer size
calculate_edge_density <- function(point, raster, buffer_size) {
  # Create buffer around the point
  buffer <- buffer(point, width = buffer_size)
  
  # Mask the raster with the buffer
  masked_raster <- mask(raster, buffer)
  
  # Calculate edge density within the masked raster
  edge_density <- lsm_l_ed(masked_raster, count_boundary = FALSE, directions = 4)
  
  # Add buffer size and point information to the result
  edge_density$buffer_size <- buffer_size
  edge_density$point_id <- point$id
  
  return(edge_density)
}

# Function to calculate edge density for all buffer sizes for a given point
calculate_for_point <- function(point, raster, buffer_sizes) {
  # Efficiently process each buffer size and combine results
  map_dfr(buffer_sizes, function(buffer_size) {
    edge_density <- calculate_edge_density(point, raster, buffer_size)
    data.table(point_id = point$id, buffer_size = buffer_size, edge_density = edge_density$value)
  })
}

# Use future_lapply to process points in parallel
result_list <- future_lapply(1:nrow(points), function(point_index) {
  point <- points[point_index, ]  # Extract individual point
  calculate_for_point(point, test_raster, buffer_sizes)
}, future.seed = TRUE)  # Ensure reproducibility in parallel processes

# Flatten the result list and convert it to a data.table
edge_density_dt <- rbindlist(result_list, use.names = TRUE)

# Print progress and the first few rows of results
message("Edge density calculation completed for all points and buffer sizes.")
print(head(edge_density_dt))
