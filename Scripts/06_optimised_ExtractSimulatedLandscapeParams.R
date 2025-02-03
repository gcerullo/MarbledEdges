library(terra)
library(sf)
library(exactextractr)
library(furrr)
library(data.table)
library(purrr)


library(terra)
library(sf)
library(data.table)
library(purrr)
library(dplyr)
library(exactextractr)

raster <- test  # Assign raster once outside the function

points <- points[1:10, ] #test function with few points

# Compute raster extent once (avoiding repeated computation)
raster_ext <- st_as_sf(st_as_sfc(st_bbox(raster)))

process_extraction <- function(raster, buffer_distances) {
  results <- vector("list", length(buffer_distances))  # Preallocate results list
  names(results) <- paste0("buffer_", buffer_distances)  # Assign buffer names
  
  # Convert points once before looping
  points_sf <- st_as_sf(points)
  
  for (buffer_distance in buffer_distances) {
    # Create buffers around each point and ensure valid geometries
    buffers_sf <- st_buffer(points_sf, dist = buffer_distance)
    buffers_sf <- st_make_valid(buffers_sf)  # Fix invalid geometries
    
    # Clip buffers to raster extent
    buffers_clipped <- st_intersection(buffers_sf, raster_ext)
    
    # Perform extraction (fail-safe)
    extraction <- tryCatch(
      exact_extract(raster, buffers_clipped, progress = TRUE),
      error = function(e) {
        message("Error in exact_extract(): ", e)
        return(NULL)
      }
    )
    
    if (is.null(extraction)) next  # Skip if extraction failed
    
    # Attach IDs directly
    extraction_with_id <- Map(function(df, id) {
      df$id <- id
      df
    }, extraction, buffers_clipped$id)
    
    results[[paste0("buffer_", buffer_distance)]] <- extraction_with_id
  }
  
  return(results)
}

# Process results
process_results <- function(results_test) {
  extracted_data <- flatten(results_test)  # Flatten nested list
  
  # Convert to data.table efficiently
  final_dt <- rbindlist(lapply(extracted_data, function(df) {
    setDT(df)[, .(total_coverage = sum(value * coverage_fraction)), by = .(value, id)]
  }), use.names = TRUE, fill = TRUE)
  
  return(final_dt)
}

# Run processing
results_test <- process_extraction(raster, buffer_distances)
final_dt <- process_results(results_test)
