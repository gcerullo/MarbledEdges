#This code computes PC1 metrics for murrelet analysis 


#read in packages 
library(tidyverse)
library(terra)

#read in inputs 

#ocean conditions 
#Raw data on the ocean conditions based on 16 PCAd variables; downloaded from: 
#https://www.fisheries.noaa.gov/west-coast/science-data/ocean-conditions-indicators-trends
ocean_cond <- read.csv("Inputs/2023-Stoplight-RAWDATA.csv") %>%  
  #remove additional data not used in PCA construction
  slice(1:21)


#--------------------------------------------------------------------
#OCEAN CONDITIONS 
#--------------------------------------------------------------------

#forest_age information 

#Clean data #### 

ocean_cond <- ocean_cond %>%  
  # Remove rows where all columns are NA
  filter(if_all(everything(), ~ !is.na(.))) %>%
  # Convert all columns starting with "X" to character
  mutate(across(starts_with("X"), as.character)) %>% 
  # Pivot to long format
  pivot_longer(
    cols = starts_with("X"),  # Select columns that start with 'X' (year columns)
    names_to = "Year",        # New column for year
    values_to = "Value"       # New column for values
  ) %>%
  # Clean up Year column to remove 'X' prefix
  mutate(Year = str_remove(Year, "^X"))


#filter PC1 data and scale and invert values (so that high values indicate good ocean years) 
PC1 <- ocean_cond %>% filter( Ecosystem.Indicators == "Principal Component scores (PC1)") %>%  
  mutate(Values = as.numeric(Value)) %>%  
  #invert values
  mutate(Value_inverted = max(Values) - Values) %>%  # Subtract each value from the maximum so that big values are good years, low values are bad years
  #scale values
  mutate(Value_scaled = scale(Value_inverted))     # Standardize values; Centers the data to a mean of 0 and scales it to a standard deviation of 1.
  
#Is there a significant change through time of  
PC1 %>%  ggplot( aes(x = Year, y = Value_scaled)) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear regression line
  labs(title = "Linear Regression: Year vs PC1_scaled", x = "Year", y = "PC1_scaled") +
  theme_minimal()

min(PC1$Value_scaled)

PC1_quantiles <-  PC1 %>%
  summarise(
    q05 = quantile(`Value_scaled`, 0.05),
    q10 = quantile(`Value_scaled`, 0.1),
    q20 = quantile(`Value_scaled`, 0.2),
    q50 = quantile(`Value_scaled`, 0.5),
    q80 = quantile(`Value_scaled`, 0.5),
    q90 = quantile(`Value_scaled`, 0.9),
    q95 = quantile(`Value_scaled`, 0.95) 
  )


PC1_plot <- PC1 %>%
  ggplot(aes(x = Year, y = Value_scaled, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  
  # Plot horizontal lines from legend_lines with mapped colors
  geom_hline(data = legend_lines, aes(yintercept = y, color = type), linetype = "dotted", size = 2) +
  
  labs(y = "Ocean condition (scaled PC1)", color = NULL) +
  
  scale_color_manual(values = c("Good ocean year" = "#2c7fb8", "Bad ocean year" = "#d95f0e")) +
  
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

#Get an elevation DEM in the correct projection
# can_cov <- rast("Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif") #canopy cover of conifers
# elevation <- rast("Rasters/DEMs/dem30m.tif")
# 
# #Step 1: Reproject elevation to EPSG:5070 (same as can_cov)
# elev_proj <- project(elevation, can_cov, method = "bilinear")
# 
# # Step 2: Crop to the extent of can_cov
# elev_crop <- crop(elev_proj, can_cov)
# plot(elev_crop)

##########################################################################
#save outputs #### 

#write.csv(PC1, "Outputs/PC1_scaled_inverted.csv")

ggsave(
  filename = "Figures/PC1_good_bad_ocean_10th_90th_percentile.png",               # File path and name
  plot = PC1_plot,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

