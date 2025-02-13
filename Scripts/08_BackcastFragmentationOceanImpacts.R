# Now for this script we are gonna use the real covariates to predict for all of the sample points that were actually sampled?  

library(terra)
library(sf)
library(exactextractr)

model <- readRDS("Models/pc1_interaction_model.rds")
source('scripts/02_OrganiseMurreletData.R')

#==================================
#Information on Murrelet survey sites: (Valente 2022)
#Most surveys were conducted around proposed timber harvest sites with 1 survey station per 8–10 ha.
#Site-level sampling effort varied with number of stations, but station size
#(200-m radius circle) was standardized (Evans Mack et al., 2003)
# so we conducted all analyses at the station level. 
#Stations were surveyed iteratively until murrelet breeding activity was recorded at the site 
#or the minimum number of required surveys was conducted#
#(5–9 surveys in each of 2 years) (Evans Mack et al., 2003).

#now for each point between 1990 and 2020 across PNW we extract the following: 
#dist coast
#hab amount 100 and 2000 (from murrelet SDM)
#edge density 100 and 2000
#ownershipodf      
#ownershipor       
#ownershipwa  
#canopy density 

#GNN 
#-TPHC_GE_3: Density of live conifers >2.5 dbh
#- CANCOV: Canopy cover of all live trees 
#- CANCOV_HDW: Canopy cover of all conifers. 

#======================
#read in ownership shapefile
ownership <- vect("Rasters/land_ownership/All_merge.shx")
#plot(ownership)

#organise raster data ####

# Step 1: Load the GNN variables (canopy cover) raster
can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_2021.tif")
can_cov_conifer <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_2021.tif")

# Align the vector to the extent of the raster using crop
ownership <- crop(ownership, ext_raster)
#plot(ownership_cropped)

# Step 2: Visualize the GNN variable raster
plot(can_cov)

# Step 3: Get the extent and CRS of the can_cov raster
ext_raster <- ext(can_cov)  # Extent of can_cov
crs_raster <- crs(can_cov)  # CRS of can_cov

# Step 4: Load the SDM (Species Distribution Model) 2019 raster
tif_files <- list.files(path = "Rasters/MAMU_SDMs", pattern = "\\.tif$", full.names = TRUE)
tif2019 <- tif_files[[32]]  
SDM2019 <- rast(tif2019)

# Step 5: Visualize the SDM raster
plot(SDM2019)

# Step 6: Ensure the SDM raster has the same CRS and extent as the GNN raster (can_cov)
# Crop SDM to the extent of the can_cov raster
SDM_crop <- crop(SDM2019, ext_raster)

# Step 7: Optionally, check if the CRS of the SDM raster matches the can_cov CRS
# If not, you can reproject SDM to match the can_cov CRS
if (crs(SDM_crop) != crs(can_cov)) {
  SDM_crop <- project(SDM_crop, crs_raster)
}

# Step 8: Ensure the resolution of the SDM raster matches the GNN raster (can_cov)
# Compare both x and y resolutions separately
if (any(res(SDM_crop) != res(can_cov))) {
  SDM_crop <- resample(SDM_crop, can_cov, method = "bilinear")  # Use resampling if needed
}

# Step 9: Now, both SDM_crop and can_cov are aligned
# Plot the cropped SDM to confirm the alignment visually
plot(SDM_crop)
plot(can_cov, add = TRUE)  # Overlay the can_cov raster to check alignment

#-------------------------------------------
#Create stratified points across landscapes
#-------------------------------------------

# Set the CRS of ownership and grid_points to match canopy_cover's CRS
ownership <- project(ownership, crs_raster)
#generate 500m grid 

# Define the resolution of the grid (500m x 500m)
grid_res_x <- 10000  # Resolution in x-direction (500m)
grid_res_y <- 10000# Resolution in y-direction (500m)

# Create a raster grid of 500m resolution over the extent of the original raster
grid_raster <- rast(nrows = ceiling((ext_raster[4] - ext_raster[3]) / grid_res_y),
                    ncols = ceiling((ext_raster[2] - ext_raster[1]) / grid_res_x),
                    extent = ext_raster,
                    res = c(grid_res_x, grid_res_y))
#ensure grid raster is in correct crs 
# Set the CRS of the grid raster to match the canopy cover raster
crs(grid_raster) <- crs_raster

# Convert the raster grid into polygons (vector cells)
grid_vect <- as.polygons(grid_raster)

# Extract the centroids of each polygon (grid cell)
grid_points <- centroids(grid_vect)

# Plot the original raster and overlay the grid points
plot(can_cov, main = "Raster with 500m Grid Points")
plot(grid_points, add = TRUE, col = "blue", pch = 16, cex = 0.5)
plot(grid_points)

points_within_non_na <- terra::extract(can_cov, grid_points, na.rm = TRUE, xy = TRUE) %>%  
  na.omit() 
points <- vect(points_within_non_na,geom = c("x", "y"), crs = crs(can_cov))


plot(can_cov, main = "Raster with 500m Grid Points")
plot(points, add = TRUE, col = "blue", pch = 16, cex = 0.5)

#====================================
#Hab amount 
#======================================
#define key paramas####
buffer_sizes <- c(100, 2000)  # Buffer sizes in meters
buffer_distances <- c(100, 2000)  # Buffer sizes in meters

# Add IDs to the filtered points
points$id <- paste0(1:nrow(points))

# Plot the raster and the filtered valid points
plot(SDM_crop)  # Plot the raster
plot(points, add = TRUE, col = "red", pch = 16)  # Overlay the filtered points

#convert SDM to binary habitat or non habitat (>.45 threshold)
habitat_binary <- ifel(SDM_crop >= 45, 1, 0)

#extract habitat amount for different buffers of each point: 
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
    
    # Identify buffer ids that remain (after removing those buffers that didn't intersect)
    buffer_ids_remaining <- buffers$id[(buffers$id %in% buffers_clipped$id)]
    
    # Perform extraction - 'frac' calculates the fraction of each value
    extraction <- exact_extract(raster, buffers_clipped, 'frac', progress = TRUE)
    extraction$id <- buffer_ids_remaining  # Add ID back in 
    
    # Store results in the list
    habAmount[[paste0("buffer_", buffer_distance)]] <- extraction
  }
  
  # Collapse the nested list into a long data frame
  habAmountlong <-  rbindlist(habAmount, idcol = "buffer_size", fill = TRUE) %>% 
    dplyr::select(buffer_size, frac_0, frac_1, id)
  
  return(habAmountlong)
}
habAmount <- process_extraction(raster = habitat_binary, buffer_distances = buffer_distances, points)


