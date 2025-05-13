#make a map of where reductions in edge density would most increase occupancy

library(tidyverse)
library(terra)

#read in params ####
results_all_edges <- readRDS("Outputs/ReduceLandscapeAndLocalEdges.rds")

final2020 <- readRDS("Outputs/PNW_2020_extracted_covars.rds") %>%   #read in starting occupancy for 2020 from scrippt 8
  as.data.frame() #%>%  

long_lat <- final2020 %>% select(point_id, x, y)

results_all_edges <- results_all_edges %>% left_join(long_lat)
can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif") #canopy cover of conifers


#--------------------------------
#carry out initial filter
results_all_edges <- results_all_edges %>% 
  filter(cancov_2020 > -1) %>% #remove non-forest habitat points
  filter( distance_to_coastline < 100000) %>% #remove points v far from coast 
 filter(dem30m < 600)  #remove points above a given altitude

#get the occupancy values at each point and edge reduction amount
decreaseEdgegoodYears <- results_all_edges %>% select(point_id, ownership, Occupancy, decrease_factor, SE,OceanYear,x,y) %>% 
  filter(OceanYear == "Good Ocean Years") %>% 
  pivot_wider(
    names_from = decrease_factor,
    values_from = c(Occupancy, SE),
    names_sep = "_"  ) %>%  
  #calculate change in occupancy and SE
  mutate(OccChange50 = Occupancy_0.5 - Occupancy_0, 
         OccChange50_SE = sqrt(SE_0.5^2 + SE_0^2), 
         OccChange100 = Occupancy_1 - Occupancy_0, 
         OccChange100_SE = sqrt(SE_1^2 + SE_0^2)         ) %>%  
  mutate(OccChange50quantile = as.numeric(ntile(OccChange50, 10)), 
         OccChange100quantile = as.numeric(ntile(OccChange100, 10))) #ten quantiles of change



decreaseEdgebadYears <- results_all_edges %>% select(point_id, ownership, Occupancy, decrease_factor, SE,OceanYear,x,y) %>% 
  filter(OceanYear == "Bad Ocean Years") %>% 
  pivot_wider(
    names_from = decrease_factor,
    values_from = c(Occupancy, SE),
    names_sep = "_") %>%  
  #calculate change in occupancy and SE
  mutate(OccChange50 = Occupancy_0.5 - Occupancy_0, 
         OccChange50_SE = sqrt(SE_0.5^2 + SE_0^2), 
         OccChange100 = Occupancy_1 - Occupancy_0, 
         OccChange100_SE = sqrt(SE_1^2 + SE_0^2)         ) %>%  
  mutate(OccChange50quantile = as.numeric(ntile(OccChange50, 10)), 
         OccChange100quantile = as.numeric(ntile(OccChange100, 10))) #ten quantiles of change


#PLot histogram of change in occupancy by ownership from reducing edge 50%
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
    title = "How much would reducing edge 50% increase occupancy by ownership?",
    x = "Occupancy Change Quantile",
    y = "Frequency"
  ) +
  facet_wrap(~ ownership, scales = "free_y", ncol = 3, nrow = 2) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "ivory", color = NA)
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = 1:10,
    labels = scales::label_number(accuracy = 1)
  )
}

quantilesOccChange50quantile <- plot_data_histogram(decreaseEdgegoodYears, OccChange50quantile)
plot_data_histogram(decreaseEdgegoodYears, OccChange50quantile)

#-----------------------------------------------
#Map priorities for reducing edges #####

TopLandsToReduceEdgeGoodYears <- decreaseEdgegoodYears %>% filter(OccChange50quantile >9) 
TopLandsToReduceEdgeBadYears <- decreaseEdgebadYears %>% filter(OccChange50quantile >9) 


# Convert the data to a spatial object
good_sf <- vect(TopLandsToReduceEdgeGoodYears, geom = c("x", "y"), crs = "EPSG:5070")
bad_sf <- vect(TopLandsToReduceEdgeBadYears, geom = c("x", "y"), crs = "EPSG:5070")

#Define the extent and resolution
# Create a raster with the desired extent and resolution
raster_template <- rast(good_sf, resolution = 1000) 
crs(raster_template) <- "EPSG:5070"  # Set CRS to EPSG 5070

#reducing edge here would be most beneficial in good ocean years 
#................................................................
#plot quantiles
good_raster <- rasterize(good_sf, raster_template, field = "OccChange50quantile", fun = "first", background = NA)
#uncertainty and raw occupancy change
good_raster_OccChange <- rasterize(good_sf, raster_template, field = "OccChange50", fun = "first", background = NA)
good_raster_SE <- rasterize(good_sf, raster_template, field = "SE_0.5", fun = "first", background = NA)
plot(good_raster_OccChange)
plot(good_raster_SE)

good_raster <- !is.na(good_raster) * 1 #get 1,0
# Replace all 0s with NA
good_raster[good_raster == 0] <- NA
plot(good_raster)

#reducing edge here would be most beneficial in bad ocean years 
#...............................................................
bad_raster <- rasterize(bad_sf, raster_template, field = "OccChange50quantile", fun = "first", background = NA)

#uncertainty and raw occupancy change
bad_raster_OccChange <- rasterize(bad_sf, raster_template, field = "OccChange50", fun = "first", background = NA)
bad_raster_SE <- rasterize(bad_sf, raster_template, field = "SE_0.5", fun = "first", background = NA)
plot(bad_raster_OccChange)
plot(bad_raster_SE)


bad_raster <- !is.na(bad_raster) * 1 #get 1,0
bad_raster[bad_raster == 0] <- NA
plot(bad_raster)



#reducing edge here would be most beneficial in both good and bad ocean years
both_have_values <- bad_raster %>% mask(good_raster)
plot(both_have_values)

both_have_values <- as.numeric(both_have_values)
both_have_values[both_have_values == 0] <- NA

plot(good_raster)
plot(bad_raster)
plot(both_have_values)

writeRaster(good_raster, "Rasters/GoodYear_top90thpercentilOccupancyIncreaseFromReducingedge50pc.tif",overwrite=TRUE)
writeRaster(bad_raster, "Rasters/BadYear_top90thpercentilOccupancyIncreaseFromReducingedge50pc.tif",overwrite=TRUE)
writeRaster(both_have_values, "Rasters/GoodAndBad_top90thpercentilOccupancyIncreaseFromReducingedge50pc.tif",overwrite=TRUE)

#-------------------------------------------
#Final plots ####
graphics.off()

# Make sure the 'Figures' directory exists or create it
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

# Open PNG device for exporting the plot
png("Figures/priority_areas_figure.png", width = 3000, height = 1200, res = 300)


#background <- resample(can_cov, good_raster, method = "bilinear")

# Set up 1-row, 3-column layout
par(mfrow = c(1, 3), mar = c(1, 1, 1, 1))

# Plot 1: Good
plot(background, 
     col = gray.colors(100, start = 1, end = 0), 
     main = "Good Ocean Year Priorities", 
     legend = FALSE, 
     axes = FALSE, 
     box = TRUE)

plot(good_raster,
     col = "blue",
     legend = FALSE,
     axes = FALSE,
     box = FALSE,
     add = TRUE)

# Plot 2: Bad
plot(background, 
     col = gray.colors(100, start = 1, end = 0), 
     main = "Bad Ocean Year Priorities", 
     legend = FALSE, 
     axes = FALSE, 
     box = TRUE)

plot(bad_raster,
     col = "orange",
     legend = FALSE,
     axes = FALSE,
     box = FALSE,
     add = TRUE)

# Plot 3: Both
plot(background, 
     col = gray.colors(100, start = 1, end = 0), 
     main = "Both", 
     legend = FALSE, 
     axes = FALSE, 
     box = TRUE)

plot(bad_raster,
     col = "darkgreen",
     legend = FALSE,
     axes = FALSE,
     box = TRUE,
     add = TRUE)

# Add a shared title
mtext("Priority areas for reducing edge amount", 
      side = 1, line = 2, outer = TRUE, cex = 1.5, font = 2)


dev.off()


#Plots of uncertainty and raw occ change estimates 

# Export to PNG
png("Figures/occ_change_se_panels.png", width = 3000, height = 3000, res = 300)

# Set layout: 2 rows, 2 columns
par(mfrow = c(2, 2), mar = c(1, 1, 2, 1), oma = c(4, 2, 2, 2))  # Adjust for minimal space

# Define color scales
occ_pal <- colorRampPalette(c("blue", "white", "red"))(100)   # For occupancy change
se_pal  <- colorRampPalette(c("white", "darkorange"))(100)    # For standard error

# Plot 1: Good Year Occupancy Change
plot(background, col = gray.colors(100, start = 1, end = 0),
     legend = FALSE, axes = FALSE, box = FALSE)
plot(good_raster_OccChange, col = occ_pal, legend = TRUE, add = TRUE)
title("Good Year Occupancy Change", line = 0.5, cex.main = 1)

# Plot 2: Good Year SE
plot(background, col = gray.colors(100, start = 1, end = 0),
     legend = FALSE, axes = FALSE, box = FALSE)
plot(good_raster_SE, col = se_pal, legend = TRUE, add = TRUE)
title("Good Year Standard Error", line = 0.5, cex.main = 1)

# Plot 3: Bad Year Occupancy Change
plot(background, col = gray.colors(100, start = 1, end = 0),
     legend = FALSE, axes = FALSE, box = FALSE)
plot(bad_raster_OccChange, col = occ_pal, legend = TRUE, add = TRUE)
title("Bad Year Occupancy Change", line = 0.5, cex.main = 1)

# Plot 4: Bad Year SE
plot(background, col = gray.colors(100, start = 1, end = 0),
     legend = FALSE, axes = FALSE, box = FALSE)
plot(bad_raster_SE, col = se_pal, legend = TRUE, add = TRUE)
title("Bad Year Standard Error", line = 0.5, cex.main = 1)

# Shared bottom title
mtext("Occupancy Change and Uncertainty Across Ocean Years",
      side = 1, line = 2, outer = TRUE, cex = 1.5, font = 2)

# Finish
dev.off()


#EXPORTS ---------------------------------------------------
#export quantiles:
ggsave(
  filename = "Figures/reducingEdge50pc_occIncrease_percentiles.tif",               # File path and name
  plot = quantilesOccChange50quantile,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)
