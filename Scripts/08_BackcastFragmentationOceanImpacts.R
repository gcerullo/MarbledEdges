# Now for this script we are gonna use the real covariates to predict for all of the sample points that were actually sampled?  

library(terra)
library(sf)
library(exactextractr)
library(data.table)
library(rnaturalearth)

#read in inputs ####
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

#now for each point in 2020 across PNW we extract the following: 
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


# Step 2: Visualize the GNN variable raster
plot(can_cov)

# Step 3: Get the extent and CRS of the can_cov raster
ext_raster <- ext(can_cov)  # Extent of can_cov
crs_raster <- crs(can_cov)  # CRS of can_cov

#plot(ownership_cropped)

# Step 4: Load the SDM (Species Distribution Model) 2019 raster
tif_files <- list.files(path = "Rasters/MAMU_SDMs", pattern = "\\.tif$", full.names = TRUE)
tif2019 <- tif_files[[32]]  
SDM2019 <- rast(tif2019)

# Step 5: Visualize the SDM raster
plot(SDM2019)

# Step 6: Ensure the SDM raster has the same CRS and extent and alignment as the GNN raster (can_cov)
# Crop SDM to the extent of the can_cov raster
# Align habitat_binary2 to pen_canopy
SDM2019 <- resample(SDM2019, can_cov, method = "bilinear")


# # Step 7: Optionally, check if the CRS of the SDM raster matches the can_cov CRS
# # If not, you can reproject SDM to match the can_cov CRS
# if (crs(SDM2019) != crs(can_cov)) {
#   SDM2019 <- project(SDM2019, crs_raster)
# }
# 
# # Step 8: Ensure the resolution of the SDM raster matches the GNN raster (can_cov)
# # Compare both x and y resolutions separately
# if (any(res(SDM2019) != res(can_cov))) {
#   SDM2019 <- resample(SDM2019, can_cov, method = "bilinear")  # Use resampling if needed
# }
# Step 9: Now, both SDM2019 and can_cov are aligned
# Plot the cropped SDM to confirm the alignment visually
plot(SDM2019)
plot(can_cov, add = TRUE)  # Overlay the can_cov raster to check alignment

#-------------------------------------------
#Create stratified points across landscapes
#-------------------------------------------


#generate 500m grid 

# Define the resolution of the grid (1000m x 1000m)
grid_res_x <- 500  # Resolution in x-direction (1000m)
grid_res_y <- 500# Resolution in y-direction (1000m)

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
plot(can_cov, main = "Raster with 500 Grid Points")
#plot(grid_points, add = TRUE, col = "blue", pch = 16, cex = 0.5)
#plot(grid_points)


points_within_non_na <- terra::extract(can_cov, grid_points, na.rm = TRUE, xy = TRUE) %>%  
  na.omit() 
points <- vect(points_within_non_na,geom = c("x", "y"), crs = crs(can_cov))

#get a binary raster of murrelet habitat


# Add IDs to the filtered points
points$id <- paste0(1:nrow(points))

plot(can_cov, main = "Raster with 500m Grid Points")
plot(points, add = TRUE, col = "blue", pch = 16, cex = 0.5)

#====================================
#Hab amount 
#======================================
#define key paramas####
buffer_sizes <- c(100, 2000)  # Buffer sizes in meters
buffer_distances <- c(100, 2000)  # Buffer sizes in meters


# Plot the raster and the filtered valid points
plot(SDM2019)  # Plot the raster
plot(points, add = TRUE, col = "red", pch = 16)  # Overlay the filtered points

#convert SDM to binary habitat or non habitat (>.45 threshold)
habitat_binary <- ifel(SDM2019 >= 45, 1, 0)

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
saveRDS(habAmount, "Outputs/PNW_2020_habamount.rds")

#====================================
#Edge amount 
#======================================

#Confusion matrix
#habitat      #open canopy         meaning 
#   8              5                40 =   habitat that is also open canopy (should be v rare/nonexistent)
#   8              2                16 =   habitat is not in open canopy (ie closed canopy habitat)
#   4              5                20 =   non-habitat open canopy 
#   4              2                8 =    non-habitat closed canopy

# identify forest edge cells - which are layers with conifer canopy coverage <= 40%
open_canopy <- ifel(can_cov_conifer <= 4500, 5, 2) #5 is open canopy, 2 is closed canopy
#plot(open_canopy)
habitat_binary2 <- ifel(SDM2019 >= 45, 8, 4)     #8 is habitat, 4 is non-habitat 
#plot(habitat_binary2)
intermediate_raster <- open_canopy * habitat_binary2
#plot(intermediate_raster)
#take only closed canopy habitat (16) and potential hardedges (20
# Filter raster: keep only values 16 or 20, set other values to NA
habitat_and_edge <- intermediate_raster
habitat_and_edge[!(habitat_and_edge == 16 | habitat_and_edge == 20)] <- NA
plot(habitat_and_edge)
# Count the number of non-NA values in the raster
freq_table <- terra::freq(habitat_and_edge)


# #find closed canopy habitat that bordered hard edges (or hard edges that border closed caopy)  
# class_difference <- boundaries(habitat_and_edge, inner = TRUE, directions = 4, classes = FALSE)
#The above is too slow so I have written a custom approach below (only need to run once): 

# Initialize a counter for the processed cells
processed_cells <- 0

#----------------------------------------------------------------

# Set up a global counter for progress tracking
processed_cells <- 0
total_cells <- ncell(habitat_and_edge)


boundary_raster <- focal(habitat_and_edge,  
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
# writeRaster(boundary_raster, "Rasters/v2_intermediate_boundaries_PNW_habitat_nonhab_edged_murrelet.tif",overwrite=TRUE)

# OLD VERSION: 
#boundary_raster <- focal(habitat_and_edge, w = matrix(c(0, 1, 0, 1, 0, 0, 1, 0, 0), 3, 3), 
#                          fun = function(x, ...) {
#                            # Update processed cell counter
#                            processed_cells <<- processed_cells + 1
#                            
#                            # Print progress every 1000 cells processed
#                            if (processed_cells %% 1000 == 0) {
#                              progress <- round((processed_cells / total_cells) * 100, 2)
#                              message("Processing: ", progress, "%")
#                            }
#                            
#                            # Check if the focal cell (center) is NA
#                            center_value <- x[5]
#                            if (is.na(center_value)) {
#                              return(NA)  # Skip if the focal cell is NA
#                            }
#                            
#                            # Check if any of the 4-connected neighbors (N, S, E, W) are different from the focal cell
#                            # Only consider neighbors that are NOT NA
#                            neighbors <- x[c(1, 3, 7, 9)]  # N, S, E, W
#                            valid_neighbors <- neighbors[!is.na(neighbors)]  # Remove NA neighbors
#                            
#                            if (any(valid_neighbors != center_value)) {
#                              return(1)  # Boundary cell if any non-NA neighbor is different
#                            } else {
#                              return(0)  # Not a boundary cell
#                            }
#                          }, pad = TRUE)
#writeRaster(boundary_raster, "Rasters/intermediate_boundaries_PNW_habitat_nonhab_edged_murrelet.tif")
#----------------------------------------------------------------

#can start here: nb, this shows all pixels where murrelet habitat meets non-habitat 
boundary_raster <- rast("Rasters/v2_intermediate_boundaries_PNW_habitat_nonhab_edged_murrelet.tif")
plot(boundary_raster)
terra::freq(boundary_raster)

#filters to only murrelet habitat that meets an edge
edgeforesthabitat <- boundary_raster*habitat_binary

plot(habitat_binary)
plot(edgeforesthabitat) 
terra::freq(habitat_binary)  #35,735,491 cells of habitat 

terra::freq(edgeforesthabitat) # 1085544 cells of habitat that meet an edge  
1085544/35735491 
#now take only boundary cells that are in murrelet habitat 
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
    
    # Identify buffer ids that remain (after removing those buffers that didn't intersect)
    buffer_ids_remaining <- buffers$id[(buffers$id %in% buffers_clipped$id)]
    
    
    # Perform extraction - 'frac' calculates the fraction of each value
    extraction <- exact_extract(raster, buffers_clipped, 'frac', progress = TRUE)
    extraction$id <- buffer_ids_remaining  # Add ID back in 
    
    # Store results in the list
    edgeAmount[[paste0("buffer_", buffer_distance)]] <- extraction
  }
  
  # Collapse the nested list into a long data frame
  edgeAmount <-  rbindlist(edgeAmount, idcol = "buffer_size", fill = TRUE) %>% 
    dplyr::select(buffer_size, frac_0, frac_1, id)
  
  return(edgeAmount)
}

#edgeforesthabitat then 1 = forest in an edge, and ignore 0 
edgeAmount <- process_extraction_edge(raster = edgeforesthabitat, buffer_distances = buffer_distances, points)
saveRDS(edgeAmount, "Outputs/v2_PNW_2020_edgeamount.rds")


####################################
######################################

#        ORGANISE COVARIATE DATA 

#####################################
#####################################
#====================================
#organise edge and habitat data ####
#====================================

#read in inputs 
habAmount <- readRDS("Outputs/PNW_2020_habAmount.rds")
edgeAmount <- readRDS("Outputs/v2_PNW_2020_edgeamount.rds")


# Process habitat amount data
hab_df <- habAmount %>% 
  pivot_wider(
    names_from = buffer_size, 
    values_from = c(frac_0, frac_1)
  ) %>%    
  mutate(habAmountDich_100 = frac_1_buffer_100, 
         habAmountDich_2000 = frac_1_buffer_2000,
         point_id = id) %>% 
  dplyr::select( point_id, habAmountDich_100, habAmountDich_2000)

hist(hab_df$habAmountDich_100)
# Process edge density data

edge_df <- edgeAmount  %>% 
  pivot_wider(
    names_from = buffer_size, 
    values_from = c(frac_0, frac_1)
  ) %>%    
  mutate(edgeRook_100_40 = frac_1_buffer_100, 
         edgeRook_2000_40 = frac_1_buffer_2000,
         point_id = id) %>% 
  dplyr::select(point_id, edgeRook_100_40, edgeRook_2000_40)
hist(edge_df$edgeRook_2000_40)

# Combine habitat and edge data
all_df <- hab_df %>% left_join(edge_df)   
 

#Get which points are actually in murrelet habitat 
point_level_habitat <- terra::extract(SDM2019, points) %>%  
  rename(point_leve_hab_probability = probability, 
         point_id = ID) %>%  
        mutate(point_id = as.character(point_id)) %>% 
  mutate(point_id =  paste0("p",point_id))

all_df <- all_df %>%left_join(point_level_habitat)

#average edge in good habitat
all_df %>% filter(point_leve_hab_probability >=45) %>%  
  summarise(meanEdge = mean(edgeRook_2000_40))

#====================================
# get distance to coastline ####
#====================================
coastline_data <- ne_coastline(scale = 110)  # scale = 110 is for 1:110m resolution
plot(coastline_data)


# Convert the coastline sf object to a terra SpatVector
coastline_vector <- vect(coastline_data)

# Get the extent of the coastline data
coastline_extent <- ext(coastline_vector)

# Define the bounding box for North America (Longitude: -170 to -50, Latitude: 24 to 85)
north_america_bbox <- ext(-170, -50, 24, 85)
# Define the bounding box for the Pacific Northwest (Longitude: -130 to -116, Latitude: 40 to 50)
pacific_nw_bbox <- ext(-170, -100, 24, 60)

# Crop the coastline vector to North America
coastline_vector <- crop(coastline_vector, pacific_nw_bbox)
plot(coastline_vector)

#reproject coastline to EPSG:5070
coastline_vector <- project(coastline_vector, "EPSG:5070")

# Plot the coastline vector (in blue)
plot(coastline_vector, col = "blue", main = "Points Over Coastline")

# Add the points on top of the coastline plot (in red)
plot(points, add = TRUE, col = "red", pch = 20)

# Calculate the distance from each point to the nearest coastline geometry
distance_to_coastline <- distance(coastline_vector, points)

# Extract the minimum distance (closest distance) for each point from the distance matrix
closest_distances <- apply(distance_to_coastline, 2, min, na.rm = TRUE)

# Create a data frame with Point ID and Distance to Coastline
distance_df <- data.frame(
  point_id = points$id,  
  distance_to_coastline = closest_distances) %>%  
  mutate(point_id =  paste0("p",point_id))


#add distance (m) to all_df 
all_df <- all_df %>% left_join(distance_df)

#====================================
# get ownership ####
#====================================
#In the valente paper, ownership is a data source covariate to account for potential detection heterogeneity introduced by local land management practices. 
#In the model; #ownership = factor: blm, odf. oregon. washington.
ownership
summary(ownership)

#If we want, we can extract the actual ownership information for each point, tho we don't need this for our models 
unique(ownership$Own_simple)


#.........................


# Set the raster resolution (to rasterise shapefile)
resolution <- 30  # Set a suitable resolution for your analysis
# Create an empty raster template with the same extent and CRS as your SpatVector
unique(ownership$Own_simple)
ownership_raster_template <- terra::rast(ownership,
                                         resolution = resolution)
# Rasterize the polygon layer
ownership_raster <- terra::rasterize(ownership, ownership_raster_template, field = "Own_simple")
ownership_raster_pj <- terra::project(ownership_raster, crs_raster) #match crs
# Resample ownership_raster to match the resolution of can_cov
ownership_raster_resampled <- terra::resample(ownership_raster_pj, can_cov, method = "near")
plot(ownership_raster_resampled)
#now extract the points 
extracted_ownership_values <- terra::extract(ownership_raster_resampled, points, bind = TRUE) %>% 
  as.data.frame() %>% 
  select(id, Own_simple) %>% 
  rename(ownership = Own_simple)
saveRDS(extracted_ownership_values, "Outputs/ownership_random_points_pnw.rds")

## The slow v accurate way:
## #match pt crs to ownership (faster than reproject the massive shapefile)
## points_pj <- project(points, crs(ownership))
## ownership_data <- terra::extract(ownership, points_pj)
## saveRDS(ownership_data, "Outputs/ownership_random_points_pnw.rds")
## .........................


#====================================
#finally assembly and scaling of data ####
#====================================
ownership_data <- readRDS("Outputs/ownership_random_points_pnw.rds") %>%  
  mutate(point_id =  paste0("p",id)) 

all_df <- all_df %>% left_join(ownership_data)


#what proportion of murrelet habitat is edgey? 
all_df %>% 
  filter(point_leve_hab_probability >= 45) %>%  
  summarise(median_edge2000 = mean(edgeRook_2000_40))

final2020 <- all_df

saveRDS(final2020, "PNW_2020_extracted_covars.rds")
