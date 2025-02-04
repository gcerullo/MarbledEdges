# Load necessary libraries
library(terra)
library(sf)
library(exactextractr)
library(data.table)
library(landscapemetrics)
library(future.apply)  
library(unmarked)
library(tidyverse)
library(purrr)

# Define input/output folders
input_folder <- "Rasters/production_0.55"
output_folder <- "Results/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# List all raster files
tif_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)

# Load multiple rasters
landscape <- rast(tif_files)

# Define buffer distances
buffer_distances <- c(1, 20)

# Extract points from raster
get_points <- function(raster) {
  cell_centers <- xyFromCell(raster, 1:ncell(raster))
  points <- vect(cell_centers, type = "points")
  points$id <- paste0("p", 1:nrow(points))
  return(points)
}

# Function to process habitat amount
process_habitat_amount <- function(raster, points, buffer_distances) {
  results <- list()
  
  for (buffer_distance in buffer_distances) {
    buffers <- buffer(points, width = buffer_distance)
    buffers_sf <- st_as_sf(buffers)
    raster_ext <- st_as_sf(st_as_sfc(st_bbox(raster)))
    buffers_clipped <- st_intersection(buffers_sf, raster_ext)
    
    extraction <- exact_extract(raster, buffers_clipped, progress = TRUE)
    
    extraction_with_id <- lapply(1:length(extraction), function(i) {
      extraction_df <- extraction[[i]]
      extraction_df$id <- buffers_clipped$id[i]
      setDT(extraction_df)[, .(
        sum_class_cells = sum(value),
        total_class_coverage_fraction = sum(coverage_fraction)
      ), by = .(id, value)]
    })
    
    results[[paste0("buffer_", buffer_distance)]] <- extraction_with_id
  }
  
  return(results)
}

# Function to calculate edge density
calculate_edge_density <- function(point, raster, buffer_size) {
  buffer <- buffer(point, width = buffer_size)
  masked_raster <- mask(raster, buffer)
  edge_density <- lsm_l_ed(masked_raster, count_boundary = FALSE, directions = 4)
  edge_density$buffer_size <- buffer_size
  edge_density$point_id <- point$id
  return(edge_density)
}

# Loop through all rasters
for (i in 1:nlyr(landscape)) {
  raster <- landscape[[i]]
  raster_name <- names(landscape)[i]
  
  message("Processing raster: ", raster_name)
  
  points <- get_points(raster)
  
  # Process habitat amount
  habAmount <- process_habitat_amount(raster, points, buffer_distances)
  saveRDS(habAmount, file = paste0(output_folder, "habAmount_", raster_name, ".rds"))
  
  # Process edge density
  edge_density_results <- list()
  for (point in 1:nrow(points)) {
    for (buffer_size in buffer_distances) {
      edge_density_results[[paste0(point, "_", buffer_size)]] <- calculate_edge_density(points[point, ], raster, buffer_size)
    }
  }
  
  edge_density_df <- do.call(rbind, edge_density_results)
  saveRDS(edge_density_df, file = paste0(output_folder, "edgeDensity_", raster_name, ".rds"))
  
  message("Completed raster: ", raster_name)
}
