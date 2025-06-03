#Functions for building scenarios 
library(tidyverse)


# Define a function to read in Harris yield curve data and organise and process each  badly organised dataset
process_data <- function(data, colnms) {
  
  #give colnanes 
  colnms <- c("Age", "mbfac", "rd20ac", "thin", "MgCac", "pfac", "adv3", "sht", "tpa", "ddi", "HEWA", "WIFL", "WIWA",  "groupID")
  
  #conversion factors
  bf_to_m3 = 0.0023597372
  acre_to_ha = 0.404686
  
  data %>%
    mutate(all_NA = ifelse(rowSums(is.na(.)) == ncol(.), TRUE, FALSE)) %>%  # Identify rows where all columns are NA
    mutate(groupID = cumsum(all_NA) + 1) %>%  # Increment the groupID after each NA row
    fill(groupID, .direction = "down") %>% 
    drop_na() %>% 
    select(-all_NA) %>%
    filter(!str_detect(.[[1]], "AGE")) %>%  # Filter rows where the first column does not contain "AGE"
    setNames(colnms) %>%   # Assign column names
    mutate(across(everything(), as.numeric)) %>%   # Convert all columns to numeric
    mutate(Age = (Age *5)-5) %>%  # Age is currently given in 5 yr increments - make true age 
    
    #convert thousands of board feet per acre to m3/ha
    mutate(
      bfac = mbfac *1000, 
      m3ac = bfac * bf_to_m3, 
      m3ha = m3ac/acre_to_ha
    ) %>% select(-c(bfac,m3ac )) %>%  
    
    #get thin volumes
    mutate(
      thin_m3ha = thin * 1000 * bf_to_m3 / acre_to_ha
    )
  
  
}

#.................................................................................
# Define a function to create the yield curve plots for the Elliot, extensive and intensive forestry plot_mbfac_curve <- function(data, dataset_name, ribbon_color = "blue") {
# 
plot_mbfac_curve <- function(data, dataset_name, ribbon_color = "blue") {
  
  # Prepare the data for the median curve and standard deviation
  summary_data <- data %>%
    group_by(Age) %>%
    summarize(median_mbfac = median(mbfac, na.rm = TRUE),
              sd_mbfac = sd(mbfac, na.rm = TRUE),
              median_m3ha = median(m3ha, na.rm = TRUE),
              sd_m3ha = sd(m3ha, na.rm = TRUE))  # Calculate the standard deviation

  # Base ggplot object
  p <- ggplot() +
    
    # Group-specific raw lines for each stand (no extrapolation)
    geom_line(data = data, aes(x = Age, y = m3ha, group = groupID, color = factor(groupID)), 
              size = 0.8, alpha = 0.2, show.legend = FALSE) +  # Slight transparency, hide legend
    
    # Add ribbon for median ± standard deviation
    geom_ribbon(data = summary_data, aes(x = Age, ymin = median_m3ha - sd_m3ha, ymax = median_m3ha + sd_m3ha),
                fill = ribbon_color, alpha = 0.4) +  # Ribbon for standard deviation range
    
    # Median curve across all groups
    geom_line(data = summary_data, aes(x = Age, y = median_m3ha),
              color = "black", size = 1.2, linetype = "solid") +  # Thicker solid black line for median
    
    # Add linear regression line through the origin (intercept = 0) in red
    geom_smooth(data = data, aes(x = Age, y = m3ha), 
                method = "lm", formula = y ~ x - 1,  # Formula for forcing intercept through origin
                color = "red", se = FALSE, size = 1, linetype = "dashed") +  # Red dashed linear regression line
    
    ylim(0, 1000) +
    xlim(0, 500) +
    
    # Customize the appearance
    theme_minimal(base_size = 16) +  # Clean theme with larger base font size
    theme(panel.grid = element_blank(),  # Remove gridlines for a cleaner look
          axis.line = element_line(color = "black"),  # Add axis lines
          axis.title = element_text(face = "bold"),  # Bold axis titles
          legend.position = "none") +  # Remove legend for simplicity
    
    # Labels
    labs(x = "Age", y = "m3/ha", 
         title = paste("Age vs m3/ha Curves for", dataset_name)) +
    
    # Customize color scale
    scale_color_manual(values = scales::hue_pal()(length(unique(data$groupID))))  # Different color for each group
  
  return(p)
}

#............................
#define a function for visualising starting landscape composition in the scenarioParams.R script.

# Function to create and plot raster for each starting landscape
plot_landscape <- function(df, landscape_name) {
  
  # Filter data for the specific landscape
  landscape_data <- df %>% filter(startingL == landscape_name)
  
  # Calculate the number of cells proportional to parcel size
  num_cells <- round(landscape_data$area / landscape_size_ha * 100)  # Scale to 100 cells total
  
  # Check if there are clearcut area, and create the values vector accordingly
  if ("clearcut" %in% landscape_data$habitat) {
    # If there is both clearcut and forest
    values <- c(rep(1, num_cells[landscape_data$habitat == "clearcut"]),
                rep(2, num_cells[landscape_data$habitat == "forest"]))
  } else {
    # Only forest, no clearcut
    values <- rep(2, num_cells[landscape_data$habitat == "forest"])
  }
  
  # Create a raster object with the calculated number of cells
  r <- rast(nrow = 1, ncol = length(values), vals = values)
  
  # Determine which colors to use, depending on the values present
  if (all(values == 2)) {
    # Only forest (value 2), use only green
    plot(r, col = "green", legend = FALSE, main = landscape_name)
  } else {
    # Mix of clearcut and forest, use grey for clearcut and green for forest
    plot(r, col = c("grey", "green"), legend = FALSE, main = landscape_name)
  }
  
  # Add forest age as text labels (centered within the forest area)
  forest_positions <- which(landscape_data$habitat == "forest")
  if (length(forest_positions) > 0) {
    text(x = (sum(num_cells[landscape_data$habitat == "clearcut"]) + num_cells[forest_positions] / 2), 
         y = 1.2, labels = landscape_data$age[forest_positions], cex = 1.2, col = "black")
  }
}

#.........................................................

# function to determine the relative (to max-even-flow) production of scenarios
determine_relative_production_targets <- function(scenario_dt, production_targets_dt, startingL_prefix) {
  
  # Ensure scenario and production targets are data.tables
  setDT(scenario_dt)
  setDT(production_targets_dt)
  
  # Filter production targets based on startingL prefix ("High", "Medium" or "Low")
  production_targets_filtered <- production_targets_dt %>%
    filter(startsWith(startingL, startingL_prefix)) %>%
    as.data.table()
  
  # Join simulated scenarios that production targets, filtering totalScenario_m3 which fall between 
  # production_target_m3_min and production_target_m3_max 
  results <- scenario_dt[production_targets_filtered, 
                         on = .(totalScenario_m3 >= production_target_m3_min, 
                                totalScenario_m3 <= production_target_m3_max), 
                         nomatch = 0, 
                         allow.cartesian = TRUE]
  
  # Create a new column that combines R, E, I, startingL, and age into one
  results[, composition := paste("R", R, "E", E, "I", I, startingL, age, sep = "_")]
  
  # Select unique rows based on the composition
  unique(results, by = "composition")
}

#...........................
# Create a function that produces a wide range of production targets, once we have calculated P=1
#for a given starting landscape
update_production_target <- function(df, target) {
  df %>%
    mutate(production_target = target, 
           production_volume = production_volume * production_target)
}

#...........................
## Function to plot how many scenarios we've created as histograms grouped by age
plot_scenario_number_histograms <- function(results_dt, title_prefix) {
  # Get unique ages in the dataset
  unique_ages <- unique(results_dt$age)
  
  # Loop through each age and plot a histogram
  for (age_val in unique_ages) {
    # Filter by age
    age_group <- results_dt[age == age_val]
    
    # Plot histogram for production_target, grouped by age
    hist(age_group$production_target, 
         main = paste(title_prefix, "Age:", age_val), 
         xlab = "Production Target", 
         ylab = "Frequency",
         col = "lightblue",
         border = "black")
  }
}




# #.............................................................
# #Create a function that estimates yields for Extensive harvesting
# #According to a range of starting Params
# 
# 
# calculate_harvested_volume <- function(A0 =  250,  # Initial age of the forest in the starting landscape
#                                        A_harvest = 90,  # Minimum age for harvesting forest extensive
#                                        V0 = 200,  # Initial standing volume of timber at A0
#                                        g = 1,  # Annual growth rate in m³/ha/year
#                                        H = 0.6,  # Percentage of volume to be harvested from a given stand
#                                        T = 150,  # Total scenario simulation time in years
#                                        delta_t = 5,  # Time interval between harvests
#                                        total_area = 100000 # Total landscape area (in ha)
# ) 
# {
#   # Calculate the year in the scenario when harvesting starts (when age reaches A_harvest).
#   #If forest is old enoough to begin harvesting immediately (leading to -ve t_harvest start), then set t_harvest_start as 0
#   t_harvest_start <- max(A_harvest - A0, 0)
#   
#   # Initialize variables to store results
#   scenario_years <- seq(t_harvest_start, T, by = delta_t)  # Years of the scenario when harvesting occurs
#   
#   # Filter out any negative years (before scenario starts)
#   scenario_years <- scenario_years[scenario_years >= 0]
#   
#   # Calculate the area dedicated to each harvest stagger
#   area_per_stagger <- total_area / length(scenario_years)
#   
#   # Initialize harvested volumes
#   harvested_volumes <- numeric(length(scenario_years))
#   
#   # Loop through each harvest year and calculate harvested volume
#   for (i in seq_along(scenario_years)) {
#     t_n <- scenario_years[i]
#     forest_age_at_tn <- A0 + t_n  # Forest age at the time of harvest
#     timber_volume_at_tn <- V0 + g * (forest_age_at_tn - A0)  # Timber volume at the time of harvest
#     harvested_volumes[i] <- H * timber_volume_at_tn * area_per_stagger  # Volume extracted per harvest stagger
#   }
#   
#   # Create a data frame for clarity
#   results <- data.frame(
#     Scenario_Year = scenario_years,
#     Forest_Age = A0 + scenario_years,
#     Timber_Volume_per_ha = V0 + g * (scenario_years - t_harvest_start),
#     Area_Harvested = rep(area_per_stagger, length(scenario_years)),
#     Harvested_Volume = harvested_volumes
#   )
#   
#   # Print the results
#   print(results)
#   
#   # Return total harvested volume
#   total_harvested_volume <- sum(harvested_volumes)
#   cat("Total Harvested Volume over the scenario:", total_harvested_volume, "m³\n")
#   
#   return(list(results = results, total_harvested_volume = total_harvested_volume))
# }
# 
# # Test the function with default parameters (A0 = 50, A_harvest = 90, V0 = 200, etc.)
# calculate_harvested_volume()
# 
# # Example 1: Initial forest age is 60, growth rate is 1.5 m³/ha/year, harvest 70%, total area 2000 ha
# calculate_harvested_volume(A0 = 60, V0 = 500, g = 1.5, H = 0.70, total_area = 2000, delta_t = 5)
# 
# # Example 2: Harvesting every 10 years, initial forest age 70, 80% of timber harvested, total area 1500 ha
# calculate_harvested_volume(A0 = 70, V0 = 250, g = 2, H = 0.80, total_area = 1500, delta_t = 10)
# 
# 
# # Example 1: Initial forest age is 60, growth rate is 1.5 m³/ha/year, harvest 70%
# x <- calculate_harvested_volume(A0 = 150, V0 = 500, g = 1.5, H = 0.80, T = 150, delta_t = 5)
# z2 <- as.data.frame(x[1])
# # Example 2: Harvesting every 10 years, initial forest age 70, 80% of timber harvested
# calculate_harvested_volume(A0 = 70, V0 = 250, g = 2, H = 0.8, T = 150, delta_t = 10)
# 
