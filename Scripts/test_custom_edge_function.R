#test my approach to calculating cells that share a boundary in NSEW direction with a different class value

library(terra)

# Example: Simulated raster data
set.seed(123)
habitat <- rast(nrows = 10, ncols = 10, vals = sample(c(0, 1), 100, replace = TRUE))
canopy <- rast(nrows = 10, ncols = 10, vals = sample(c(0, 1), 100, replace = TRUE, prob = c(0.8, 0.2)))

# Create binary masks
habitat_mask <- habitat == 1
canopy_mask <- canopy == 1

# Check adjacency: Identify cells with at least one open canopy neighbor
adjacent_open <- focal(canopy_mask, w = matrix(c(0, 1, 0,
                                                 1, 0, 1,
                                                 0, 1, 0), nrow = 3), fun = max, fill = NA, na.rm = TRUE)

# Exclude central cells that are open canopy
habitat_only_mask <- habitat_mask & !canopy_mask

# Final result: Habitat cells adjacent to open canopy
result <- habitat_only_mask & adjacent_open

# Plot results
par(mfrow = c(1, 3))
plot(habitat, main = "Habitat Raster (1 = Habitat, 0 = Non-Habitat)", col = c("white", "green"))
plot(canopy, main = "Canopy Raster (1 = Open, 0 = Closed)", col = c("white", "blue"))
plot(result, main = "Result: Habitat Adjacent to Open Canopy", col = c("white", "red"))

#-----------------------------
#With real data:: 
habitat_binary2 <- ifel(SDM2020 >= 45, 8, 4) #8 is habitat 
open_canopy <- ifel(can_cov <= 4000, 5, 2) #5 is open canopy, 2 is closed canopy
habitat_mask <- habitat_binary2 == 8
canopy_mask <- open_canopy == 5

plot(habitat_mask)
plot(canopy_mask)


# Check adjacency: Identify cells with at least one open canopy neighbor
adjacent_open <- focal(canopy_mask, w = matrix(c(0, 1, 0,
                                                 1, 0, 1,
                                                 0, 1, 0), nrow = 3), fun = max, fill = NA, na.rm = TRUE)

# Exclude central cells that are open canopy
habitat_only_mask <- habitat_mask & !canopy_mask

# Final result: Habitat cells adjacent to open canopy
result <- habitat_only_mask & adjacent_open
plot(result)
writeRaster(result, "Rasters/test_edge_approach.tif",overwrite=TRUE)

