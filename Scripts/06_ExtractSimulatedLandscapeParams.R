#Extract covariates from simulated rasters. To fit with our models we need: 
#1. Distance to coast ; scaleCoastDist          
#2. Habitat amount 100 and 2km;  scaleHabAmount100 and scaleHabAmount2000                    
#3. Edge density 100 and 2km ; scaleEdgeDens100 and scaleEdgeDens2000         
library(terra)
library(sf)
library(exactextractr)
library(data.table)

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

#100*100 = 10,000 squares
#so each square is 1km2 to make 1Mha
#EXTRACTION TO BUFFERS ####
#vector based extraction of buffer data

#1Mha landscape; so each square is 1ha
# Initialize lists to store results
results_habAmount <- list()

buffer_distances <- c(1,20)

#TO CORRECT!!!!!!!!!!!!!!!
#The 1million cell raster is currently at 1m2 resolution; actually each cell should correspond to a ha!!!
#Come back and make sure tihs is correct before extracting the buffer!!!!!!!! 



buffer_distance <- 100
raster <- test
crs(raster) <- "+proj=utm +zone=11 +datum=WGS84"
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
    
    # Store results with buffer distance information
    
    # Name each list element by StationID_Year
    # named_extraction <- setNames(extraction, buffers_clipped$StationID_Year)
    
    
    # # Store results with buffer distance information
    #  results[[paste0("buffer_", buffer_distance)]] <- list(
    #    extraction = extraction,
    #   
    #   #keep statation information
    #   # station_info = buffers_clipped$StationID_Year
    #  )
  }
  return(extraction)
}

# Process results
results_test <- process_extraction(test, buffer_distances)
results_test[[100]]

#Capture cells within each  ####

# Function to summarize LIDAR extractions
summarize_extractions <- function(results_LiDAR) {
  
  # Use map to apply a function to each element of the results_LiDAR list
  map(results_LiDAR, ~{
    
    # Access the extractions and station info within each buffer list
    results_test <- .x$extraction
    station_info <- .x$station_info
    
    # Use map2 to iterate over extractions and station info simultaneously
    map2(extractions, station_info, ~{
      # Calculate summary statistics for each extraction DataFrame
      df <- as_tibble(.x)
      summarise(df, 
                Mean = mean(value, na.rm = TRUE), 
                SD = sd(value, na.rm = TRUE), 
                Median = median(value, na.rm = TRUE)) %>%
        # Add station info as a new column for reference
        mutate(StationID_Year = .y)
    }) %>%
      # Optional: Combine all mini summaries into a single tibble
      bind_rows()
  })
}

# Process the results_LiDAR data
summarized_lidar_data <- summarize_lidar_extractions(results_LiDAR) %>%  
  #add in information on buffer size
  rbindlist(idcol = "buffer_size") %>% 
  #make buffer size column numeric
  mutate(buffer_size = str_extract(buffer_size, "\\d+") %>% as.numeric()) %>% 
  #bring buffer size back to metres 
  mutate(buffer_size = buffer_size/metre_to_ft)