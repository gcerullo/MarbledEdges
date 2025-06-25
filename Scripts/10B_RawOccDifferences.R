#plots both binary 90th percentile top edge-reduction areas, and with raw occ increases and SE 
library(tidyverse)
library(terra)

#read in params ####
results_all_edges <- readRDS("Outputs/ReduceLandscapeAndLocalEdges.rds") %>%  
  dplyr::select(point_id, Occupancy, OceanYear, decrease_factor)

final2020 <- readRDS("Outputs/PNW_2020_extracted_covars.rds") %>%   #read in starting occupancy for 2020 from scrippt 8
  as.data.frame() #%>%  
raw_diff_se <- readRDS("Outputs/raw_diff_se_05edge_reduction.rds") %>%  
  rename( raw_diff = estim, 
          SE = se)# %>% 
raw_diff_se <- raw_diff_se %>% left_join(results_all_edges) %>% 
  filter(decrease_factor %in% c(0, 0.5))

long_lat <- final2020 %>% dplyr::select(point_id, x, y)

raw_diff_se <- raw_diff_se %>% left_join(long_lat)
can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif") #canopy cover of conifers
can_cov_binary <- ifel(can_cov > 0, 1, NA)
plot(can_cov_binary)

#--------------------------------
#carry out initial filter
raw_diff_se <- raw_diff_se %>% 
  filter(cancov_2020 > -1) %>% #remove non-forest habitat points
  filter( distance_to_coastline < 60000) %>% #remove points v far from coast 
  filter(dem30m < 600)  #remove points above a given altitude

#get the occupancy values at each point and edge reduction amount
decreaseEdgegoodYears <- raw_diff_se %>%
  filter(OceanYear == "Good Ocean Years") %>% 
  mutate(OccChange50quantile = as.numeric(ntile(raw_diff, 10))) %>% 
  filter(decrease_factor == 0.5)


decreaseEdgebadYears <- raw_diff_se %>%
  filter(OceanYear == "Bad Ocean Years") %>% 
  mutate(OccChange50quantile = as.numeric(ntile(raw_diff, 10)))%>% 
  filter(decrease_factor == 0.5)

#Plot histogram of change in occupancy by ownership from reducing edge 50%
#calculate difference in occupancy if we removed 50% of the edge 
plot_data_histogram <- function(data, x_var){
  filtered_data <- data %>%
    filter(!ownership %in% c("Unknown", NA)) %>%
    filter(ownership %in% c("Federal", "Private Industrial", "Private Conservation", 
                            "Indian", "Private Non-industrial", "State")) %>%
    mutate(ownership = factor(ownership, levels = c("Private Conservation", "Federal", 
                                                    "State", "Indian", 
                                                    "Private Non-industrial", 
                                                    "Private Industrial")))
  
  # Compute medians and counts
  summary_stats <- filtered_data %>%
    group_by(ownership) %>%
    summarise(
      median_value = median({{ x_var }}, na.rm = TRUE),
      count = n_distinct(point_id),  # Adjust if point_id exists
      .groups = "drop"
    )
  
  # Create styled plot
  ggplot(filtered_data, aes(x = {{ x_var }})) +
    geom_histogram(binwidth = 1, fill = "#4C9A2A", color = "black", alpha = 0.7) +
    geom_vline(data = summary_stats, aes(xintercept = median_value),
               linetype = "dashed", color = "red", linewidth = 1) +
    geom_text(data = summary_stats,
              aes(x = Inf, y = Inf, label = paste0("n = ", count)),
              hjust = 1.2, vjust = 1.5, size = 5, inherit.aes = FALSE) +
    labs(
      x = "Percentiles of raw occupancy increase from reducing edge amount",
      y = "Frequency"
    ) +
    facet_wrap(~ ownership, scales = "free_y", ncol = 3, nrow = 2) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "top",  # Position the legend at the top
      legend.title = element_blank(),  # Increase legend title size for clarity
      legend.text = element_text(size = 14),  # Increase legend text size for better readability
      axis.title = element_text(size = 16),  # Increase axis title font size
      axis.text = element_text(size = 14),  # Increase axis label font size
      panel.grid.major = element_blank(),  # No major gridlines
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
}

quantilesOccChange50quantileGoodYears <- plot_data_histogram(decreaseEdgegoodYears, OccChange50quantile)
quantilesOccChange50quantileBadYears <- plot_data_histogram(decreaseEdgebadYears, OccChange50quantile)

# quantilesOccChange100quantileGoodYears <- plot_data_histogram(decreaseEdgegoodYears, OccChange100quantile)
# quantilesOccChange100quantileBadYears <- plot_data_histogram(decreaseEdgegoodYears, OccChange100quantile)
plot_data_histogram(decreaseEdgegoodYears, OccChange50quantile)

#-------------------------------------------


########VIRIDIS
# Template raster
good_sf <- vect(decreaseEdgegoodYears, geom = c("x", "y"), crs = "EPSG:5070")
template <- rast(good_sf, resolution = 1000)
crs(template) <- "EPSG:5070"
background <- terra::resample(can_cov_binary, template, method = "bilinear")  


map_edge_priority <- function(
    good_data,   
    bad_data,
    quant_col = "OccChange50quantile", 
    occchange_col = "OccChange50",
    se_col = "SE_0.5",
    quantile_threshold = 9,
    background_raster,
    output_prefix = "Figures/priority_rsw_occ"
) {
  library(terra)
  library(viridis)  # For viridis color scales
  
  # Filter top priority areas
  good_top <- good_data[good_data[[quant_col]] > quantile_threshold, ]
  bad_top  <- bad_data[bad_data[[quant_col]] > quantile_threshold, ]
  
  # Convert to spatial objects
  good_sf <- vect(good_top, geom = c("x", "y"), crs = "EPSG:5070")
  bad_sf  <- vect(bad_top, geom = c("x", "y"), crs = "EPSG:5070")
  
  # Template raster
  template <- rast(good_sf, resolution = 1000)
  crs(template) <- "EPSG:5070"

  
  # # add your occupancy layers on top  
  # plot(good_occ, col = occ_pal, breaks = good_occ_breaks, legend = TRUE, add = TRUE)  
  
  # Rasterize for each input
  good_occ <- rasterize(good_sf, template, field = occchange_col, fun = "first", background = NA)
  good_se  <- rasterize(good_sf, template, field = se_col, fun = "first", background = NA)
  
  bad_occ  <- rasterize(bad_sf, template, field = occchange_col, fun = "first", background = NA)
  bad_se   <- rasterize(bad_sf, template, field = se_col, fun = "first", background = NA)
  
  # Compute decile (10-bin quantile) breaks from raster values
  deciles <- seq(0, 1, length.out = 11)
  
  good_occ_breaks <- quantile(values(good_occ), probs = deciles, na.rm = TRUE)
  good_se_breaks  <- quantile(values(good_se),  probs = deciles, na.rm = TRUE)
  
  bad_occ_breaks  <- quantile(values(bad_occ),  probs = deciles, na.rm = TRUE)
  bad_se_breaks   <- quantile(values(bad_se),   probs = deciles, na.rm = TRUE)
  
  # Define 10-color viridis palettes
  occ_pal <- viridis(10, option = "viridis")   # or "plasma", "inferno", etc.
  se_pal  <- viridis(10, option = "magma")     # different palette for SE for contrast
  
  # Export maps
  png(paste0(output_prefix, "occ_se_map.png"), width = 3000, height = 3000, res = 300)
  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1), oma = c(4, 2, 2, 2))
  
  # Export maps
  pdf(paste0(output_prefix, "occ_se_map.pdf"), width = 10, height = 10) # Use PDF export
  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1), oma = c(4, 2, 2, 2))
  
  # Plot good year occupancy change
  plot(background_raster, col = "#E5E5E5",
       legend = FALSE, axes = FALSE, box = FALSE)
  plot(good_occ, col = occ_pal, breaks = good_occ_breaks, legend = TRUE, add = TRUE)
  title("Good Year Occupancy Change", line = 0.2, cex.main = 1)
  
  # Plot good year SE
  plot(background_raster, col = "#E5E5E5",
       legend = FALSE, axes = FALSE, box = FALSE)
  plot(good_se, col = se_pal, breaks = good_se_breaks, legend = TRUE, add = TRUE)
  title("Good Year Standard Error", line = 0.2, cex.main = 1)
  
  # Plot bad year occupancy change
  plot(background_raster, col = "#E5E5E5",
       legend = FALSE, axes = FALSE, box = FALSE)
  plot(bad_occ, col = occ_pal, breaks = bad_occ_breaks, legend = TRUE, add = TRUE)
  title("Bad Year Occupancy Change", line = 0.2, cex.main = 1)
  
  # Plot bad year SE
  plot(background_raster, col = "#E5E5E5",
       legend = FALSE, axes = FALSE, box = FALSE)
  plot(bad_se, col = se_pal, breaks = bad_se_breaks, legend = TRUE, add = TRUE)
  title("Bad Year Standard Error", line = 0.2, cex.main = 1)
  
  mtext("Change in Occupancy 
         from reducing edginess 50% in priority areas",
        side = 1, line = 2, outer = TRUE, cex = 1.5, font = 2)
  dev.off()
  
  message("Maps saved to: ", output_prefix, "*")
}
map_edge_priority(
  good_data = decreaseEdgegoodYears,
  bad_data = decreaseEdgebadYears,
  quant_col = "OccChange50quantile",
  occchange_col = "raw_diff",
  se_col = "SE",
  background_raster = background,
  output_prefix = "Figures/top90_lands_edge_50pc_raw_occ"
)


#--------------------------------------------------------
#plot BINARY map of top 90 areas
map_edge_priority_binary <- function(
    good_data,   #decide to plot for good or bad years
    bad_data,
    quant_col = "OccChange50quantile",  #e.g. look at Occ
    occchange_col = "raw_diff",
    quantile_threshold = 9,
    background_raster,
    output_prefix = "Figures/TESTpriority_"
) {
  library(terra)
  
  # Filter top priority areas
  good_top <- good_data[good_data[[quant_col]] > quantile_threshold, ]
  bad_top <- bad_data[bad_data[[quant_col]] > quantile_threshold, ]
  
  # Convert to spatial objects
  good_sf <- vect(good_top, geom = c("x", "y"), crs = "EPSG:5070")
  bad_sf  <- vect(bad_top, geom = c("x", "y"), crs = "EPSG:5070")
  
  # Template raster
  template <- rast(good_sf, resolution = 1000)
  crs(template) <- "EPSG:5070"
  
  # Rasterize for each input
  good_raster_bin <- !is.na(rasterize(good_sf, template, field = quant_col)) * 1
  bad_raster_bin  <- !is.na(rasterize(bad_sf, template, field = quant_col)) * 1
  good_raster_bin[good_raster_bin == 0] <- NA
  bad_raster_bin[bad_raster_bin == 0] <- NA
  
  both_raster_bin <- mask(bad_raster_bin, good_raster_bin)
  both_raster_bin[both_raster_bin == 0] <- NA
  
  good_occ <- rasterize(good_sf, template, field = occchange_col, fun = "first", background = NA)
  bad_occ  <- rasterize(bad_sf, template, field = occchange_col, fun = "first", background = NA)


  
  
  # Create output directory
  dir.create(dirname(output_prefix), showWarnings = FALSE)
  
  # Export map: binary priorities
  png(paste0(output_prefix, "binary_map.png"), width = 1500, height = 1200, res = 300)
  par(mfrow = c(1, 3), mar = c(0.5, 0.5, 3, 0.5))
  
  # Export map: binary priorities
  pdf(paste0(output_prefix, "binary_map.pdf"), width = 10, height = 8) # Use PDF export
  par(mfrow = c(1, 3), mar = c(0.5, 0.5, 3, 0.5))
  
  # Named list of rasters and their colors
  priority_rasters <- list(
    "Good Ocean Year Priorities" = list(r = good_raster_bin, col = "#56B4E9"),
    "Bad Ocean Year Priorities" = list(r = bad_raster_bin, col =   "#D55E00"),
    "Both" = list(r = both_raster_bin, col = "darkgreen")
  )
  
  
  
  # Plot loop
  for (plot_name in names(priority_rasters)) {
    plot(background_raster,
         col = "#E5E5E5",
         legend = FALSE, axes = FALSE, box = TRUE)
    
    plot(priority_rasters[[plot_name]]$r,
         col = priority_rasters[[plot_name]]$col,
         legend = FALSE, axes = FALSE, box = FALSE, add = TRUE)
    
    title(plot_name, line = 0.5, cex.main = 1)
  }
  
  mtext("Priority areas for reducing edge amount", 
        side = 1, line = 2, outer = TRUE, cex = 1.5, font = 2)
  dev.off()


  message("Maps saved to: ", output_prefix, "*")
}

map_edge_priority_binary(
  good_data = decreaseEdgegoodYears,
  bad_data = decreaseEdgebadYears,
  quant_col = "OccChange50quantile", 
  occchange_col = "raw_diff",
  background_raster = background,
  output_prefix = "Figures/BINARY_edge_50pc_"
)


# #-----------------------------------------------
# Make tifs for exporting priority areas
# # #-----------------------------------------------
# # #priority lands
  TopLandsToReduceEdgeGoodYears <- decreaseEdgegoodYears %>% filter(OccChange50quantile >9)  
# 
  TopLandsToReduceEdgeBadYears <- decreaseEdgebadYears %>% filter(OccChange50quantile >9) 
# # 
# # # Convert the data to a spatial object
  good_sf <- vect(TopLandsToReduceEdgeGoodYears, geom = c("x", "y"), crs = "EPSG:5070")
  bad_sf <- vect(TopLandsToReduceEdgeBadYears, geom = c("x", "y"), crs = "EPSG:5070")
# # 
# # #Define the extent and resolution
# # # Create a raster with the desired extent and resolution
  raster_template <- rast(good_sf, resolution = 1000) 
  crs(raster_template) <- "EPSG:5070"  # Set CRS to EPSG 5070
# # 
# # #reducing edge here would be most beneficial in good ocean years 
# # #................................................................
# # #plot quantiles
  good_raster <- rasterize(good_sf, raster_template, field = "OccChange50quantile", fun = "first", background = NA)
 good_raster <- !is.na(good_raster) * 1 #get 1,0
 good_raster[good_raster == 0] <- NA

# #...............................................................
 bad_raster <- rasterize(bad_sf, raster_template, field = "OccChange50quantile", fun = "first", background = NA)

 bad_raster <- !is.na(bad_raster) * 1 #get 1,0
 bad_raster[bad_raster == 0] <- NA
#
# 
 both_have_values <- bad_raster %>% mask(good_raster)
 both_have_values <- as.numeric(both_have_values)
 both_have_values[both_have_values == 0] <- NA
# 
 plot(good_raster)
 plot(bad_raster)
 plot(both_have_values)
# 

writeRaster(good_raster, "Rasters/GoodYear_odds_top90thAbsoluteOccupancyIncreaseFromReducingedge50pc.tif",overwrite=TRUE)
writeRaster(bad_raster, "Rasters/BadYear_odds_top90thAbsoluteOccupancyIncreaseFromReducingedge50pc.tif",overwrite=TRUE)
writeRaster(both_have_values, "Rasters/GoodAndBad_odds_top90thAbsoluteOccupancyIncreaseFromReducingedge50pc.tif",overwrite=TRUE)

#area of restoration priority
good_raster %>% as.data.frame() %>% count()
bad_raster %>% as.data.frame() %>% count()
both_have_values %>%  as.data.frame() %>% count()


#EXPORTS ---------------------------------------------------
#export quantiles by actor bar charts 
ggsave(
  filename = "Figures/reducingEdge50pc_RAWoccIncrease_percentiles_Good_years.tif",               # File path and name
  plot = quantilesOccChange50quantileGoodYears,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

ggsave(
  filename = "Figures/reducingEdge50pc_RAWoccIncrease_percentiles_Bad_years_.tif",               # File path and name
  plot = quantilesOccChange50quantileBadYears,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)
