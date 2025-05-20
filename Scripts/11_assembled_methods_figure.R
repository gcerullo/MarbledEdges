#11. Build nice summary figure for the manuscript (components assembled in inkscape)

#read in packages 
library(tidyverse)
library(terra)
source('scripts/02_OrganiseMurreletData.R')

#====================================

#READ IN PARAMS 
#Survey sites
analysisSurveys
#SDM Map 
#organise raster data ####

# MAMU SDM
can_cov_all_trees <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_2020.tif")
can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif") #canopy cover of conifers
plot(can_cov)
ext_raster <- ext(can_cov)  # Extent of can_cov
crs_raster <- crs(can_cov)  # CRS of can_cov
tif_files <- list.files(path = "Rasters/MAMU_SDMs", pattern = "\\.tif$", full.names = TRUE)
tif2020 <- tif_files[[36]]  
SDM2020 <- rast(tif2020)
SDM2020 <- resample(SDM2020, can_cov, method = "bilinear")
plot(SDM2020)

#edge amount
hard_edges <- rast( "Rasters/murrelet_hard_edges.tif")
plot(hard_edges)

#habitat amount





########################################################

#Nice figure 

# Convert raster to a data frame for ggplot2
SDM_2020_df <- as.data.frame(SDM2020, xy = TRUE, na.rm = TRUE)
colnames(SDM_2020_df) <- c("x", "y", "value") %>%  
  filter(value >0)

# Optional: Load a shapefile for country boundaries or study area
boundaries <- st_read("path_to_your_shapefile.shp")

# Create the map
ggplot() +
  geom_raster(data = SDM_2020_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = boundaries, fill = NA, color = "black", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "SDM 2020") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  coord_sf()

