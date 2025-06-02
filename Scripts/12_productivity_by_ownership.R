#plotsite productive by distance to coast 
library(terra)
library(tidyverse)

#read in ownership shapefile
ownership <- vect("Rasters/land_ownership/All_merge.shx")
final2020 <- readRDS("Outputs/PNW_2020_extracted_covars.rds") %>%   #read in starting occupancy for 2020 from scrippt 8
  as.data.frame()
point_id <- final2020 %>% dplyr::select(point_id, x,y)
long_lat <- vect(point_id, geom = c("x", "y"), crs = "EPSG:5070") # Assuming coordinates are in WGS84
#read in the forest productivity 
#organise raster data ####

# Step 1: Load the GNN variables (canopy cover) raster
can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif") #canopy cover of conifers
# Step 3: Get the extent and CRS of the can_cov raster
ext_raster <- ext(can_cov)  # Extent of can_cov
crs_raster <- crs(can_cov)  # CRS of can_cov

# #reample productivity
#productivity <- rast("Rasters/ForestSiteProductivity/Forest_Site_Productivity.tif") 
# productivity_proj <- project(productivity, crs(can_cov))
# # Resample the reprojected raster to align with can_cov
# productivity_pj <- resample(productivity_proj, can_cov, method = "bilinear")
# productivity_pj <- resample(productivity, can_cov, method = "median")
# plot(productivity_pj)
#writeRaster(productivity_pj, "Rasters/ForestSiteProductivity/Forest_Site_Productivity_pj.tif")

productivity <- rast("Rasters/ForestSiteProductivity/Forest_Site_Productivity_pj.tif") 
plot(productivity)
# productivity Code Description (i invert to so that big = more wood)
# 1 = 225+ cubic feet/acre/year.
# 2 = 165-224 cubic feet/acre/year.
# 3 = 120-164 cubic feet/acre/year.
# 4 = 85-119 cubic feet/acre/year.
# 5 = 50-84 cubic feet/acre/year.
# 6 = 20-49 cubic feet/acre/year.
# 7 = 0-19 cubic feet/acre/year.

#I invert this so that: 
# 7 = 225+ cubic feet/acre/year.
# 6 = 165-224 cubic feet/acre/year.
# 5 = 120-164 cubic feet/acre/year.
# 4 = 85-119 cubic feet/acre/year.
# 3 = 50-84 cubic feet/acre/year.
# 2 = 20-49 cubic feet/acre/year.
# 1 = 0-19 cubic feet/acre/year.

#Extract productivity for each point 
points_within_non_na <- terra::extract(productivity, long_lat, na.rm = TRUE, xy = TRUE) %>%  
  na.omit() %>%  
  as.data.frame()

#add information on productivty for each point 
final2020 <- final2020 %>% left_join(points_within_non_na) %>%  
  rename(productivity=Band_1)

invert_values <- function(values) {
  values_no_na <- na.omit(values)  # Remove NA values
  max_val <- max(values_no_na)    # Calculate max without NA
  min_val <- min(values_no_na)    # Calculate min without NA
  inverted <- max_val + min_val - values
  return(inverted)
}

# Apply the function to the SpatRaster
#productivity_inverted <- app(productivity, fun = invert_values)


filtered_data <- final2020 %>%
  filter(distance_to_coastline <100000) %>% 
  filter(ownership %in% c("Federal", "Private Industrial", "Private Conservation", 
                          "Private Non-industrial", "Indian", "State")) %>%
  mutate(ownership = factor(ownership, levels = c("Private Conservation", "Federal", 
                                                  "State", "Indian", 
                                                  "Private Non-industrial", 
                                                  "Private Industrial"))) %>%  
  mutate(inverted_productivity = invert_values(productivity)) %>% 
  #add categorical productivity
  mutate(cat_productivity = case_when(
    inverted_productivity == 7 ~ "225+",
    inverted_productivity == 6 ~ "165–224",
    inverted_productivity == 5 ~ "120–164",
    inverted_productivity == 4 ~ "85–119",
    inverted_productivity == 3 ~ "50–84",
    inverted_productivity == 2 ~ "20–49",
    inverted_productivity == 1 ~ "0–19",
    TRUE ~ NA_character_
  )) %>%
  mutate(cat_productivity = factor(cat_productivity, levels = c(
    "0–19", "20–49", "50–84", "85–119", "120–164", "165–224", "225+"
  )))

# Plot

gam <- ggplot(filtered_data, aes(x = distance_to_coastline / 1000, y = inverted_productivity, color = ownership)) +
  geom_smooth(method = "gam", se = TRUE, size = 1.2, alpha = 0.2) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c(
    "Private Conservation" = "#1b9e77",
    "Federal" = "#66a61e",
    "State" = "#7570b3",
    "Indian" = "#e7298a",
    "Private Non-industrial" = "#d95f02",
    "Private Industrial" = "#e6ab02"
  )) +
    labs(
      x = "Distance to Coastline (km)",
      y = expression("Potential productivity (ft"^3*"/acre/year)"),  # Customize y-axis label
      color = "Ownership" )+ 
  facet_wrap(~ ownership, ncol = 3) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c(
      "0–19", "20–49", "50–84", "85–119",
      "120–164", "165–224", "225+"
    ),
    limits = c(1, 8)
  ) +
  theme(
    axis.title.x = element_blank(),  # Removes x-axis title

    axis.text.x = element_blank(),  # Removes x-axis values
    axis.ticks.x = element_blank(), # Removes x-axis ticks
    legend.position = "none",       # Legend at the top
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),  # Increase legend text size
    axis.title = element_text(size = 16),   # Larger axis titles
    axis.text = element_text(size = 14),    # Larger axis labels
    strip.text = element_text(size = 14, face = "bold"),  # Bold facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)  # Add border
  )


# Plot
lm <- ggplot(filtered_data, aes(x = distance_to_coastline/1000, y = inverted_productivity, color = ownership)) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, alpha = 0.2) +  # Smooth line with confidence interval
  #geom_point(alpha = 0.3, size = 1) +  # Optional: Add scatter points for raw data
  theme_classic(base_size = 14) +  # Nature-style theme
  scale_color_manual(values = c(
    "Private Conservation" = "#1b9e77",
    "Federal"= "#66a61e",
    "State" = "#7570b3",
    "Indian" = "#e7298a",
    "Private Non-industrial" =  "#d95f02",
    "Private Industrial" = "#e6ab02"
  ))  +
  labs(
    x = "Distance to Coastline (km)",
    y = expression("Potential productivity (ft"^3*"/acre/year)"),  # Customize y-axis label
    color = "Ownership" )+ 
  facet_wrap(~ ownership, ncol = 3) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c(
      "0–19", "20–49", "50–84", "85–119",
      "120–164", "165–224", "225+"
    ),
    limits = c(1, 8)
  ) +
  theme(
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis values
    axis.ticks.x = element_blank(), # Removes x-axis ticks
    legend.position = "none",       # Legend at the top
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),  # Increase legend text size
    axis.title = element_text(size = 16),   # Larger axis titles
    axis.text = element_text(size = 14),    # Larger axis labels
    strip.text = element_text(size = 14, face = "bold"),  # Bold facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)  # Add border
  )





#==========================
data_summary <- filtered_data %>%
  group_by(ownership, distance_to_coastline = round(distance_to_coastline / 1000)) %>%
  summarize(
    mean_productivity = mean(inverted_productivity, na.rm = TRUE),
    sd_productivity = sd(inverted_productivity, na.rm = TRUE),
    .groups = "drop"
  )

# Add lower and upper bounds for SD
data_summary <- data_summary %>%
  mutate(
    lower_bound = mean_productivity - sd_productivity,
    upper_bound = mean_productivity + sd_productivity
  )

# Plot with standard deviation ribbons
meanSD <- ggplot(data_summary, aes(x = distance_to_coastline, y = mean_productivity, color = ownership)) +
  geom_line(size = 1.2) +  # Mean trend line
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = ownership), alpha = 0.2, color = NA) +  # SD ribbon
  theme_classic(base_size = 14) +
  scale_color_manual(values = c(
    "Private Conservation" = "#1b9e77",
    "Federal"= "#66a61e",
    "State" = "#7570b3",
    "Indian" = "#e7298a",
    "Private Non-industrial" =  "#d95f02",
    "Private Industrial" = "#e6ab02"
  )) +
  scale_fill_manual(values = c(
    "Private Conservation" = "#1b9e77",
    "Federal"= "#66a61e",
    "State" = "#7570b3",
    "Indian" = "#e7298a",
    "Private Non-industrial" =  "#d95f02",
    "Private Industrial" = "#e6ab02"
  )) +
  labs(
    x = "Distance to Coastline (km)",
    y = expression("Potential productivity (ft"^3*"/acre/year)"),  # Customize y-axis label
    color = "Ownership" )+ 
  facet_wrap(~ ownership, ncol = 3) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c(
      "0–19", "20–49", "50–84", "85–119",
      "120–164", "165–224", "225+"
    ),
    limits = c(1, 8)
  ) +
  
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 16),  # Larger axis titles
    axis.text = element_text(size = 14),  # Larger axis labels
    strip.text = element_text(size = 14, face = "bold"),  # Bold facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)  # Add border
  )


# Combine plots into a grid with titles
yields_distCoast_diffMethods <- cowplot::plot_grid(
  lm, gam, meanSD,
  nrow = 3,
  labels = c("Linear Model", "GAM", "Mean with SD"),
  label_size = 14,
  label_fontface = "bold"
)


#Export figure 
ggsave(
  filename = "Figures/yields_distCoast_diffMethods.png",  # File path and name
  plot = yields_distCoast_diffMethods,                                         # The plot object
  width = 8.27,                                                  # A4 width in inches
  height = 11.69,                                                # A4 height in inches
  dpi = 300,                                                     # Print-quality resolution
  units = "in",                                                  # Units in inches
  device = "png",                                                # Save as PNG
  bg = "white"                                                   # White background
)
