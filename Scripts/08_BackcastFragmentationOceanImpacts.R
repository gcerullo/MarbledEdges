# Now for this script we are gonna use the real covariates to predict for all of the sample points that were actually sampled?  

library(terra)
library(sf)
library(exactextractr)
library(data.table)
library(rnaturalearth)
library(tidyverse)
library(summarytools)

#read in inputs ####
#source('scripts/02_OrganiseMurreletData.R')
#rm(pointsInRoiAndSamplingWindow)
#rm(analysisData)
#==================================
#for each point in 2020 across PNW we extract the following: 
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
can_cov_all_trees <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_2020.tif")
can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif") #canopy cover of conifers

# Step 2: Visualize the GNN variable raster
plot(can_cov)

# Step 3: Get the extent and CRS of the can_cov raster
ext_raster <- ext(can_cov)  # Extent of can_cov
crs_raster <- crs(can_cov)  # CRS of can_cov

# Step 4: Load the SDM (Species Distribution Model) 2019 raster
tif_files <- list.files(path = "Rasters/MAMU_SDMs", pattern = "\\.tif$", full.names = TRUE)
tif2020 <- tif_files[[36]]  
print(tif2020)
SDM2020 <- rast(tif2020)

# Step 5: Visualize the SDM raster
plot(SDM2020)

# Step 6: Ensure the SDM raster has the same CRS and extent and alignment as the GNN raster (can_cov)
# Crop SDM to the extent of the can_cov raster
# Align habitat_binary2 to pen_canopy
SDM2020 <- resample(SDM2020, can_cov, method = "bilinear")


# # Step 7: Optionally, check if the CRS of the SDM raster matches the can_cov CRS
# # If not, you can reproject SDM to match the can_cov CRS
# if (crs(SDM2020) != crs(can_cov)) {
#   SDM2020 <- project(SDM2020, crs_raster)
# }
# 
# # Step 8: Ensure the resolution of the SDM raster matches the GNN raster (can_cov)
# # Compare both x and y resolutions separately
# if (any(res(SDM2020) != res(can_cov))) {
#   SDM2020 <- resample(SDM2020, can_cov, method = "bilinear")  # Use resampling if needed
# }
# Step 9: Now, both SDM2020 and can_cov are aligned
# Plot the cropped SDM to confirm the alignment visually
plot(SDM2020)
plot(can_cov, add = TRUE)  # Overlay the can_cov raster to check alignment

#check write rasters and check in qgis 

 #writeRaster(SDM2020, "Rasters/sdm_mamu2020_pj.tif", overwrite=TRUE)
 # writeRaster(can_cov, "Rasters/can_cov2019.tif", overwrite=TRUE)

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


#only keep points that fall within the murrelets range 
points_within_non_na <- terra::extract(SDM2020, grid_points, na.rm = TRUE, xy = TRUE) %>%  
  na.omit() %>%  
  filter(probability > 0)
rm()
  
points <- vect(points_within_non_na,geom = c("x", "y"), crs = crs(SDM2020))

#get a binary raster of murrelet habitat


# Add IDs to the filtered points
points$id <- paste0("p", 1:nrow(points))
points$ID <- NULL


plot(SDM2020, main = "Raster with 500m Grid Points")
plot(points, add = TRUE, col = "blue", pch = 16, cex = 0.5)

writeVector(grid_points, "Vectors/grid_points.shp", overwrite=TRUE)

#====================================
#Hab amount 
#======================================
#define key paramas####
buffer_sizes <- c(100, 2000)  # Buffer radius sizes in meters
buffer_distances <- c(100, 2000)  # Buffer radius sizes in meters


# Plot the raster and the filtered valid points
plot(SDM2020)  # Plot the raster
plot(points, add = TRUE, col = "red", pch = 16)  # Overlay the filtered points

#convert SDM to binary habitat or non habitat (>.45 threshold)
SDM2020

habitat_binary <- ifel(SDM2020 >= 45, 1, 0)

#extract habitat amount for different buffers of each point: 
process_extraction <- function(raster, buffer_distances, points) {
  habAmount <- list()
  
  # Ensure points have unique IDs
 # points$id <- paste0("p", 1:nrow(points))
  
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
#CORRECT EDGE AMOUUNT

# identify forest edge cells - which are layers with conifer canopy coverage <= 40%
#>-1 excludes non-forest eg urban areas etc 
open_canopy <- ifel(can_cov <= 4000 & can_cov >= 0, 1, 0)
# writeRaster(open_canopy, "Rasters/open_canopy_5_is_open.tif",overwrite=TRUE)

#plot(open_canopy)
habitat_binary <- ifel(SDM2020 >= 45, 1, 0)    

habitat_mask <- habitat_binary == 1
canopy_mask <- open_canopy == 1
plot(habitat_mask)
# Check adjacency: Identify cells with at least one open canopy neighbor
adjacent_open <- focal(canopy_mask, w = matrix(c(0, 1, 0,
                                                 1, 0, 1,
                                                 0, 1, 0), nrow = 3), fun = max, fill = NA, na.rm = TRUE)

# Exclude central cells that are open canopy
habitat_only_mask <- habitat_mask & !canopy_mask

# Final result: Habitat cells adjacent to open canopy
result <- habitat_only_mask & adjacent_open
plot(result)
#writeRaster(result, "Rasters/murrelet_hard_edges.tif",overwrite=TRUE)


#extract edge amount per buffer size 

process_extraction_edge <- function(raster, buffer_distances, points) {
  edgeAmount <- list()
  
  # Ensure points have unique IDs
#  points$id <- paste0("p", 1:nrow(points))
  
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
murrelet_hard_edges <- rast("Rasters/murrelet_hard_edges.tif")

edgeAmount <- process_extraction_edge(raster = murrelet_hard_edges, buffer_distances = buffer_distances, points)
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

#check we have same num of rows 
print(habAmount)
print(edgeAmount)

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
hist(hab_df$habAmountDich_2000)

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
hist(edge_df$edgeRook_100_40)
hist(edge_df$edgeRook_2000_40)

# Combine habitat and edge data
all_df <- hab_df %>% left_join(edge_df)   


#get points with ID and xy coords as dataframe
point_df <- as.data.frame(points, geom = "xy") %>%  
  rename(point_id = id, 
         point_leve_habitat = probability)



all_df <- all_df %>%left_join(point_df, by = "point_id")


#average edge in good habitat
all_df %>% filter(point_leve_habitat >=45) %>%  
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
  distance_to_coastline = closest_distances) 

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
resolution <- 30  # Set a suitable resolution for  analysis
# Create an empty raster template with the same extent and CRS as  SpatVector
unique(ownership$Own_simple)
ownership_raster_template <- terra::rast(ownership,
                                         resolution = resolution)
# Rasterize the polygon layer
ownership_raster <- terra::rasterize(ownership, ownership_raster_template, field = "Own_simple")
ownership_raster_pj <- terra::project(ownership_raster, crs_raster) #match crs

# Resample ownership_raster to match the resolution of can_cov
ownership_raster_resampled <- terra::resample(ownership_raster_pj, can_cov, method = "near")
plot(ownership_raster_resampled)

names(ownership_raster_resampled) <- "Ownership_Data"

ownership_df <- ownership_raster_resampled %>% as.data.frame() 
plot

#writeRaster(ownership_raster_resampled, "Rasters/test_reprojected_ownership.tif")
ownership_raster_resampled <- rast("Rasters/test_reprojected_ownership.tif")

checkPoints <- vect(all_df, geom = c("x", "y"), crs = "EPSG:5070")  # Specify a CRS, e.g., WGS84
plot(ownership_raster_resampled, main = "Raster with 500m Grid Points")
plot(checkPoints, add = TRUE, col = "blue", pch = 16, cex = 0.5)

extracted_ownership_values <- terra::extract(ownership_raster_resampled, points, bind = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::select(id, Ownership_Data) %>% 
  rename(ownership = Ownership_Data)
saveRDS(extracted_ownership_values, "Outputs/ownership_random_points_pnw.rds")

#why are so many NAs being created?!! 
extracted_ownership_values %>% group_by(ownership) %>% count()


## The slow v accurate way:
## #match pt crs to ownership (faster than reproject the massive shapefile)
## points_pj <- project(points, crs(ownership))
## ownership_data <- terra::extract(ownership, points_pj)
## saveRDS(ownership_data, "Outputs/ownership_random_points_pnw.rds")
## .........................



#--------------------------------
# add information on canopy cover so that we can use a forest/non-forest mask to filter points 
#now extract the points 

can_cov_all_trees
points
extracted_cancov_values <- terra::extract(can_cov_all_trees, points, bind = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::select(id, cancov_2020) %>%  
  rename(point_id =  id) 

# writeRaster(can_cov_all_trees, "Rasters/test_cancov_alltrees.tif")
plot(can_cov, main = "Canopy Cover")
plot(points, add = TRUE, col = "blue", pch = 16, cex = 0.5)

#------------------------------------------
# add informaton on elevation 
# #quickly get elevation for all of the points 
# final2020_vect = vect(final2020, geom = c("x", "y"), crs = "EPSG:5070")
# 
# #Get an elevation DEM in the correct projection (takes 20mins or so)
# elevation <- rast("Rasters/DEMs/dem30m.tif")
# 
# #Step 1: Reproject elevation to EPSG:5070 (same as can_cov)
# elev_proj <- project(elevation, can_cov, method = "bilinear")
# 
# # Step 2: Crop to the extent of can_cov
# elev_crop <- crop(elev_proj, can_cov)
# pt_elevations <- terra::extract(elev_crop, final2020_vect, bind = TRUE) %>% as.data.frame() %>%  
#   dplyr::select(point_id, dem30m)
# write.csv(pt_elevations, "Outputs/elevation_at_each_point.csv")


#====================================
#final assembly and scaling of data ####
#====================================
ownership_data <- readRDS("Outputs/ownership_random_points_pnw.rds") %>%  
  rename(point_id =  id) 

all_df <- all_df %>% left_join(ownership_data) %>%  left_join(extracted_cancov_values)
dfSummary(all_df) #ownership (32% points miss data;  cancov_2020   20% of points miss data 

ownership_NAs <- all_df %>% filter(!is.na(ownership)) 

#check why we have NAs 
write.csv(ownership_NAs, "points_withNA_ownership.csv")
# Convert to SpatVector
ownership_NAs_vect <- vect(ownership_NAs, geom = c("x", "y"), crs = "EPSG:5070")

# Set up a 1-row, 2-column plotting layout
par(mfrow = c(1, 2))
# Plot the raster
plot(ownership_raster_resampled, main = "Ownership Raster")
# Plot the vector
plot(ownership_NAs_vect, main = "Ownership Vector", add = TRUE, col = "red")

# Reset plotting layout (optional)
par(mfrow = c(1, 1))

#remove points we don't have ownership or canopy cover data for 
dfSummary(all_df)
final2020 <- all_df %>% 
filter(!is.na(ownership)) %>%  
  filter(!is.na(cancov_2020)) 

#add elevation data 
pt_elevation <- read.csv("Outputs/elevation_at_each_point.csv") 
final2020<- final2020 %>% left_join(pt_elevation)
final2020 <- final2020 %>% select()
dfSummary(final2020)

checkPoints <- vect(final2020, geom = c("x", "y"), crs = "EPSG:5070")  # Specify a CRS, e.g., WGS84
plot(SDM2020, main = "Raster with 500m Grid Points")
plot(checkPoints, add = TRUE, col = "blue", pch = 16, cex = 0.5)

#

# #---------------------------------------------------------------------------
# #quickly summarise data #### (can read in final2020 at bottom)
# #---------------------------------------------------------------------------
# #what was the elevational range of murrelets in dataset from pre-clearance surveys? 
# elevational_range_df <- read.csv("Inputs/pointsInRoiAndSamplingWindow_withDEM.csv")
# hist(elevational_range_df$dem30m) 
# max(elevational_range_df$dem30m, na.rm = TRUE) #no surveys above 1636.42 metres! 
# q95elev = as.numeric(quantile(elevational_range_df$dem30m, 0.95, na.rm = TRUE))
# q90elev = as.numeric(quantile(elevational_range_df$dem30m, 0.90, na.rm = TRUE))
# q80elev = as.numeric(quantile(elevational_range_df$dem30m, 0.80, na.rm = TRUE))
# 
# quickplotSummaryMedianSD <- function(x, plot_title = "Median and Standard Deviation for Each Variable") {
#   
#   quickplotSummary <- x %>%
#     filter(!ownership == "Unknown") %>% 
#     filter(cancov_2020 >0 ) %>%  
#     filter(distance_to_coastline > 100000) %>%  
#     filter(dem30m < q90elev) %>% 
#     
#     
#     group_by(ownership) %>%
#     summarise(
#       across(
#         c(habAmountDich_100, habAmountDich_2000, edgeRook_100_40, edgeRook_2000_40, point_leve_habitat),
#         list(
#           median = median,
#           sd = sd,
#           n = ~n()
#         ),
#         .names = "{.col}_{.fn}"
#       )
#     ) %>%
#     pivot_longer(
#       cols = -ownership,
#       names_to = c("variable", "stat"),
#       names_pattern = "(.*)_(median|sd|n)",
#       values_to = "value"
#     ) %>%
#     pivot_wider(
#       names_from = stat,
#       values_from = value
#     ) %>%
#     mutate(
#       lower_ci = median - 1.96 * (sd / sqrt(n)),
#       upper_ci = median + 1.96 * (sd / sqrt(n))
#     )
#   
#   # plot summaries
#   ggplot(quickplotSummary, aes(x = ownership, y = median, color = ownership)) +
#     geom_point(size = 3) +
#     geom_errorbar(aes(ymin = median - lower_ci, ymax = median + upper_ci), width = 0.2) +
#     facet_wrap(~ variable, scales = "free_y") +
#     labs(
#       title = plot_title,
#       x = "Ownership",
#       y = "Median ± SD"
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       legend.position = "none"
#     )
# }
# 
# allpoints_median <- quickplotSummaryMedianSD(final2020, plot_title = "All Points: Median and SD")
# 
# allforest_median <- quickplotSummaryMedianSD(
#      final2020 %>% filter(cancov_2020 > 0),
# 
#     plot_title = "Forest Points : Median and SD")
# 
# allhabpoints_median <- quickplotSummaryMedianSD(
#   final2020 %>% filter(point_leve_habitat > 45),
#   plot_title = "Habitat > 45: Median and SD"
# )
# allhab100kmpoints_median <- quickplotSummaryMedianSD(
#   final2020 %>% filter(point_leve_habitat > 45, distance_to_coastline < 100000),
#   plot_title = "Habitat > 45 & <100km from Coast: Median and SD"
# )
# 
# 
# #With 95% CI 
# quickplotSummaryMean95 <- function(x, plot_title = "Mean and 95% CI for Each Variable") {
#   
#   quickplotSummary <- x %>%
#     filter(!ownership == "Unknown") %>% 
#     filter(cancov_2020 >0 ) %>%  
#     filter(distance_to_coastline > 100000) %>%  
#     filter(dem30m < q90elev) %>% 
#    
#     
#     group_by(ownership) %>%
#     summarise(
#       across(
#         c(habAmountDich_100, habAmountDich_2000, edgeRook_100_40, edgeRook_2000_40, point_leve_habitat),
#         list(
#           mean = mean,
#           sd = sd,
#           n = ~n()
#         ),
#         .names = "{.col}_{.fn}"
#       )
#     ) %>%
#     pivot_longer(
#       cols = -ownership,
#       names_to = c("variable", "stat"),
#       names_pattern = "(.*)_(mean|sd|n)",
#       values_to = "value"
#     ) %>%
#     pivot_wider(
#       names_from = stat,
#       values_from = value
#     ) %>%
#     mutate(
#       lower_ci = mean - 1.96 * (sd / sqrt(n)),
#       upper_ci = mean + 1.96 * (sd / sqrt(n))
#     )
#   
#   # plot summaries with 95% CI
#   ggplot(quickplotSummary, aes(x = ownership, y = mean, color = ownership)) +
#     geom_point(size = 3) +
#     geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
#     facet_wrap(~ variable, scales = "free_y") +
#     labs(
#       title = plot_title,
#       x = "Ownership",
#       y = "Mean ± 95% CI"
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       legend.position = "none"
#     )
# }
# 
# allpoints_mean <- quickplotSummaryMean95(final2020, plot_title = "All Points: Mean and 95% CI")
# 
# allhabpoints_mean <- quickplotSummaryMean95(
#   final2020, #%>% filter(point_leve_habitat > 45),
#   plot_title = "AllPoints: Mean and 95% CI"
# )
# 
# allforest_meanDistanceOnly <- quickplotSummaryMean95(
#   final2020 %>% filter( distance_to_coastline < 100000),
#   #final2020 %>% filter(cancov_2020 > 4000),
#   plot_title = "AllPoints <84km : Mean and 95% CI")
# 
# allforest_meanDistanceOnlyCanCov <- quickplotSummaryMean95(
#   final2020 %>% filter( distance_to_coastline < 84000, cancov_2020 >0),
#   #final2020 %>% filter(cancov_2020 > 4000),
#   plot_title = "Forest Points (>0% canopy) <84km : Mean and 95% CI")
# 
# cowplot::plot_grid(allforest_meanDistanceOnly,allforest_meanDistanceOnlyCanCov)
# 
# allhab100kmpoints_mean <- quickplotSummaryMean95(
#   final2020 %>% filter(point_leve_habitat > 45, distance_to_coastline < 100000),
#   plot_title = "Habitat > 45 & <100km from Coast: Mean and 95% CI"
# )
# 
# allpoints_mean
# allhabpoints_mean
# allforest_mean
# allhab100kmpoints_mean
#   #  filter(distance_to_coastline < 100000) %>% 
#   #filter(point_leve_habitat > 0) %>%


#EXPORT OUTPUT #####
saveRDS(final2020, "Outputs/PNW_2020_extracted_covars.rds")
write.csv(final2020, "Outputs/PNW_2020_extracted_covars.csv")
#read in for quick plotting 
pt_elevations = read.csv("Outputs/elevation_at_each_point.csv")
final2020 <- read.csv("Outputs/PNW_2020_extracted_covars.csv") %>%  left_join(pt_elevations)


#Export for visual checks in Qgis 
# checkPoints <- vect(final2020, geom = c("x", "y"), crs = "EPSG:5070")  # Specify a CRS, e.g., WGS84
# plot(SDM2020, main = "Raster with 500m Grid Points")
# plot(checkPoints, add = TRUE, col = "blue", pch = 16, cex = 0.5)
# 
# #add in information on elevation at each point 
# pt_elevation <- read.csv("Outputs/elevation_at_each_point.csv") 
# 
# lowhabitat_points <- final2020 %>% 
#   filter(habAmountDich_2000 <0.05)
# 
# lowhabitat_points_federal <- final2020 %>% 
#   filter(habAmountDich_2000 <0.05) %>% 
#   filter(ownership == "Federal") 
# lowhabitat_points_state <- final2020 %>% 
#   filter(habAmountDich_2000 <0.05) %>% 
#   filter(ownership == "State")
# lowhabitat_points_private_industrial <- final2020 %>% 
#   filter(habAmountDich_2000 <0.05) %>% 
#   filter(ownership == "Private Industrial")
# lowhabitat_points_privateNonindustrial <- final2020 %>% 
#   filter(habAmountDich_2000 <0.05) %>% 
#   filter(ownership == "Private Non-industrial")
# 
# 
# lowhabitat_points_federal_100000m_elev794m <- final2020 %>% 
#   filter(habAmountDich_2000 <0.05) %>% 
#   filter(ownership == "Federal") %>%  
#   filter(distance_to_coastline <100000) %>%  
#   filter(cancov_2020 >-1) %>%  
#   left_join(pt_elevations) %>% 
#   filter(dem30m < q90elev)
# 
# lowhabitat_points_private_industrial_100000m_elev794m <- final2020 %>% 
#   filter(habAmountDich_2000 <0.05) %>% 
#   filter(ownership == "Private Industrial") %>%  
#   filter(distance_to_coastline <100000) %>%  
#   filter(cancov_2020 >-1) %>%  
#   left_join(pt_elevations) %>%  
#   filter(dem30m < q90elev)
# 
# checkPoints <- vect(lowhabitat_points_federal_84000m, geom = c("x", "y"), crs = "EPSG:5070")  # Specify a CRS, e.g., WGS84
# plot(SDM2020, main = "Raster with 500m Grid Points")
# plot(checkPoints, add = TRUE, col = "blue", pch = 16, cex = 0.5)
# 
# #test what's happening in Qgis
# write.csv(final2020, "PNW_2020_extracted_covars.csv")
# #writeRaster(ownership_raster_resampled, "Rasters/ownership_raster_resampled.tif", overwrite=TRUE)
# write.csv(lowhabitat_points, "Rasters/points_less005_2km_hab_amount.csv")
# write.csv(lowhabitat_points_federal, "Rasters/federal_points_less005_2km_hab_amount.csv")
# write.csv(lowhabitat_points_state, "Rasters/state_points_less005_2km_hab_amount.csv")
# write.csv(lowhabitat_points_private_industrial, "Rasters/privateIndustrial_points_less005_2km_hab_amount.csv")
# write.csv(lowhabitat_points_privateNonindustrial, "Rasters/privateNonIndustrial_points_less005_2km_hab_amount.csv")
# write.csv(lowhabitat_points_federal_100000m_elev794m, "Rasters/lowhabitat_points_federal_100000m_elev794m.csv")
# write.csv(lowhabitat_points_private_industrial_100000m_elev794m, "Rasters/lowhabitat_points_private_industrial_100000m_elev794m.csv")

