#test my approach to calculating cells that share a boundary in NSEW direction with a different class value


# Create a small simulated raster for testing purposes showing closed canopy habitat (1) and open non-habitat (2)
raster_data <- matrix(c(
  0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,
  0, 0, 0, 0, 0,0,0,1,0,0,0,0,1,
  0, 0, 1, 0, 0,0,1,1,1,0,1,0,0,
  0, 0, 0, 0, 0,1,1,1,0,0,0,0,0,
  0, 0, 0, 0, 0,0,1,1,0,0,0,0,0
), nrow = 5, byrow = TRUE)

# Convert the matrix to a SpatRaster
raster <- rast(raster_data)
plot(raster)

# Set up a global counter for progress tracking
processed_cells <- 0
total_cells <- ncell(raster)


boundary_raster <- focal(raster,  
                         w = c(3, 3),  
                         fun = function(x, ...) {
                           
                           # Update processed cell counter
                           processed_cells <<- processed_cells + 1
                           if (processed_cells %% 1000 == 0) {
                             progress <- round((processed_cells / total_cells) * 100, 2)
                             message("Processing: ", progress, "%")
                           }
                           
                           center_value <- x[5]  # Explicitly extract center, even if NA
                           
                           # # Print for debugging
                           # message("Focal window: ", paste(x, collapse = ", "),
                           #         " | Center: ", center_value)
                           
                           if (is.na(center_value)) {
                             return(NA)  # Keep NA if the focal cell itself is NA
                           }
                           
                           # Extract NSEW neighbors
                           neighbors <- x[c(2, 4, 6, 8)]  # N, S, E, W
                           
                           # Remove NA neighbors before comparison
                           valid_neighbors <- neighbors[!is.na(neighbors)]
                           
                           # Check if any non-NA neighbor is different
                           if (any(valid_neighbors != center_value)) {
                             return(1)  # Boundary cell
                           } else {
                             return(0)  # Not a boundary
                           }
                         }, pad = TRUE)  # Keep padding enabled



ncell(raster)
# Plot the results
par(mfrow = c(1, 2))  # Arrange plots side by side
plot(raster, main = "Original Raster")
plot(boundary_raster, main = "Boundary Raster with Updated Values")


forest_with_edge <- raster*boundary_raster
plot(forest_with_edge, main = "Forest cells
     meet hard edge")
