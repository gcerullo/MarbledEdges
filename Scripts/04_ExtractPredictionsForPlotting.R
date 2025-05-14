# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)
library(ggcorrplot)

source("scripts/01_DeterminePredictors.R")
source("scripts/02_OrganiseMurreletData.R")
# Read in unmarked object
analysisData <- readRDS("Outputs/analysisDataUnmarked.rds")
model <- readRDS("Models/final_model_5thMay2025.rds")



# Step 1: Extract site-level covariates from the unmarked object
# -------------------------------------------------------------
# Step 1: Extract Mean Values for `scaleDoy` and `scaleDoy2`
doy_means <- analysisSurveys %>%
  summarize(
    mean_scaleDoy = mean(scaleDoy, na.rm = TRUE),
    mean_scaleDoy2 = mean(scaleDoy^2, na.rm = TRUE)
  )

# Step 2: Define Prediction Grid
edge_values <- seq(
  min(siteCovs(analysisData)$scaleEdgeDens2000, na.rm = TRUE),
  max(siteCovs(analysisData)$scaleEdgeDens2000, na.rm = TRUE),
  length.out = 100
)

q05 <- PC1_quantiles %>% pull(q10)
q10 <- PC1_quantiles %>% pull(q10)
q20 <- PC1_quantiles %>% pull(q20)
q50 <- PC1_quantiles %>% pull(q50)
q80 <- PC1_quantiles %>% pull(q80)
q90 <- PC1_quantiles %>% pull(q90)
q95 <- PC1_quantiles %>% pull(q95)


pc1_levels <- c(q10,
                #q50,
                q90)
                #takes 10th percentile
               

predict_data <- expand.grid(
  scaleEdgeDens2000 = edge_values,
  PC1_t1 = pc1_levels
)

# Add other covariates
predict_data <- predict_data %>%
  mutate(
    scaleCoastDist = mean(siteCovs(analysisData)$scaleCoastDist, na.rm = TRUE),
    scaleHabAmount100 = mean(siteCovs(analysisData)$scaleHabAmount100, na.rm = TRUE),
    scaleEdgeDens100 = mean(siteCovs(analysisData)$scaleEdgeDens100, na.rm = TRUE),
    scaleHabAmount2000 = mean(siteCovs(analysisData)$scaleHabAmount2000, na.rm = TRUE),
    scaleDoy = doy_means$mean_scaleDoy,
    scaleDoy2 = doy_means$mean_scaleDoy2
  )



# Step 3: Predict Occupancy with standard errors (for error ribbon)
predictions <- predict(
  model,  
  newdata = predict_data,
  type = "state", 
  se.fit = TRUE  # Obtain standard errors for predictions
)

# Add predicted occupancy and standard errors to the data
predict_data$Occupancy <- predictions$Predicted
predict_data$SE <- predictions$SE

# Calculate 95% confidence intervals for occupancy (using 1.96 * SE for a 95% CI)
predict_data$LowerCI <- predict_data$Occupancy - 1.96 * predict_data$SE
predict_data$UpperCI <- predict_data$Occupancy + 1.96 * predict_data$SE

# Get mean and SD of the original (unscaled) Edge Density variable
edge_mean <- mean(siteCovs(analysisData)$edgeRook_2000_40, na.rm = TRUE)
edge_sd <- sd(siteCovs(analysisData)$edgeRook_2000_40, na.rm = TRUE)

# Add unscaled values to the dataset
predict_data$UnscaledEdgeDens2000 <- (predict_data$scaleEdgeDens2000 * edge_sd) + edge_mean

# Add labels for Ocean Year conditions
predict_data$OceanYear <- factor(
  predict_data$PC1_t1,
  levels = pc1_levels,
  labels = c("Bad Ocean Years", 
           # "Neutral Ocean Years", 
             "Good Ocean Years")
)

# Plot figure with error ribbons ####

# Define a  color palette with muted tones
nature_palette <- c("Bad Ocean Years" = "#D55E00",   # Grey for "Bad" condition
                   # "Neutral Ocean Years" = "#A0A0A0", # Lighter grey for "Neutral"
                    "Good Ocean Years" = "#56B4E9")   # Muted blue for "Good" condition

ocean_murrelet_by_ocean_and_edge <- 
  ggplot(predict_data, aes(x = UnscaledEdgeDens2000, y = Occupancy, color = OceanYear)) +
  geom_line(size = 1.5) +  # Thicker lines for clearer visibility
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = OceanYear), 
              alpha = 0.3, color = NA) +  # Add error ribbons
  scale_color_manual(
    values = nature_palette,  # Apply the Nature-like muted color palette
    name = "Ocean Conditions"
  ) +
  scale_fill_manual(
    values = nature_palette,  # Apply same palette to ribbon fill
    name = "Ocean Conditions"
  ) +
  labs(
    x = "Proportion of Edge (in 2000m buffer)",  # Clear x-axis label
    y = "Predicted Occupancy Probability",  # Clear y-axis label
  ) +
  theme_minimal(base_size = 18) +  # Minimal theme with larger base font size
  theme(
    legend.position = "top",  # Position the legend at the top
    legend.title = element_blank(),  # Increase legend title size for clarity
    legend.text = element_text(size = 14),  # Increase legend text size for better readability
    axis.title = element_text(size = 16),  # Increase axis title font size
    axis.text = element_text(size = 14),  # Increase axis label font size
    panel.grid.major = element_blank(),  # No major gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a thin border around the plot
  )

ocean_murrelet_by_ocean_and_edge

#------------------------------------------------
#GET ESTIMATES OF INTERACTING PC1 & FRAGMENTATION FOR DIFFERENT DISTANCES FROM THE COAST #####
#------------------------------------------------
coast_mean <- mean(siteCovs(analysisData)$coastDist, na.rm = TRUE)
coast_sd <- sd(siteCovs(analysisData)$coastDist, na.rm = TRUE)

coast_levels <- c(
  (1000 - coast_mean) / coast_sd,
  (20000 - coast_mean) / coast_sd,
  (40000 - coast_mean) / coast_sd, 
  (88000 - coast_mean) / coast_sd
  )

predict_data2 <- expand.grid(
  scaleEdgeDens2000 = edge_values,
  PC1_t1 = pc1_levels,
  scaleCoastDist =coast_levels  # Now using scaled coastal distances
)

# Add other covariates
predict_data2 <- predict_data2 %>%
  mutate(
    scaleHabAmount100 = mean(siteCovs(analysisData)$scaleHabAmount100, na.rm = TRUE),
    scaleEdgeDens100 = mean(siteCovs(analysisData)$scaleEdgeDens100, na.rm = TRUE),
    scaleHabAmount2000 = mean(siteCovs(analysisData)$scaleHabAmount2000, na.rm = TRUE),
    scaleDoy = doy_means$mean_scaleDoy,
    scaleDoy2 = doy_means$mean_scaleDoy2
  )

predictions2 <- predict(
  model,  
  newdata = predict_data2,
  type = "state", 
  se.fit = TRUE
)

# Add predicted values to `predict_data`
predict_data2$Occupancy <- predictions2$Predicted
predict_data2$SE <- predictions$SE

# Calculate 95% confidence intervals
predict_data2$LowerCI <- predict_data2$Occupancy - 1.96 * predict_data2$SE
predict_data2$UpperCI <- predict_data2$Occupancy + 1.96 * predict_data2$SE

predict_data2$CoastDist <- (predict_data2$scaleCoastDist * coast_sd) + coast_mean

# Convert CoastDist to a factor for faceting
predict_data2$CoastDist <- factor(predict_data2$CoastDist, 
                                  levels = c(1000, 20000, 40000,88000), 
                                  labels = c("On Coast",
                                             "20 km inland", 
                                             "40 km inland", 
                                             "80 km inland"))

#add labels 
# Add labels for Ocean Year conditions
predict_data2$OceanYear <- factor(
  predict_data2$PC1_t1,
  levels = pc1_levels,
  labels = c("Bad Ocean Years", 
             # "Neutral Ocean Years", 
             "Good Ocean Years" ))


# Add unscaled values to the dataset
predict_data2$UnscaledEdgeDens2000 <- (predict_data2$scaleEdgeDens2000 * edge_sd) + edge_mean


ocean_murrelet_hab_occ_varrying_dist <- 
  ggplot(predict_data2, aes(x = UnscaledEdgeDens2000, y = Occupancy, color = OceanYear)) +
  geom_line(size = 1.5) +  
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = OceanYear), 
              alpha = 0.3, color = NA) +  
  scale_color_manual(
    values = nature_palette,  
    name = "Ocean Conditions"
  ) +
  scale_fill_manual(
    values = nature_palette,  
    name = "Ocean Conditions"
  ) +
  labs(
    x = "Proportion of Edge (in 2000m buffer)",
    y = "Predicted Occupancy Probability",
  ) +
  facet_wrap(~CoastDist, ncol = 4) +  # Facet by Coastal Distance
  theme_minimal(base_size = 18) +  # Minimal theme with larger base font size
  theme(
    legend.position = "top",  # Position the legend at the top
    legend.title = element_blank(),  # Increase legend title size for clarity
    legend.text = element_text(size = 14),  # Increase legend text size for better readability
    axis.title = element_text(size = 16),  # Increase axis title font size
    axis.text = element_text(size = 14),  # Increase axis label font size
    panel.grid.major = element_blank(),  # No major gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a thin border around the plot
  )

ocean_murrelet_hab_occ_varrying_dist

#------------------------------------------------
#Check occupancy for different amounts of habitat and edge
# -------------------------------------------------------------
# Step 1: Extract Mean Values for `scaleDoy` and `scaleDoy2`
doy_means <- analysisSurveys %>%
  summarize(
    mean_scaleDoy = mean(scaleDoy, na.rm = TRUE),
    mean_scaleDoy2 = mean(scaleDoy^2, na.rm = TRUE)
  )

# Step 2: Define Prediction Grid
edge_values <- seq(
  min(siteCovs(analysisData)$scaleEdgeDens2000, na.rm = TRUE),
  max(siteCovs(analysisData)$scaleEdgeDens2000, na.rm = TRUE),
  length.out = 100
)

# Step 2: Define hab values Grid
hab_values <- seq(
  min(siteCovs(analysisData)$scaleHabAmount2000, na.rm = TRUE),
  max(siteCovs(analysisData)$scaleHabAmount2000, na.rm = TRUE),
  length.out = 100
)

q90 <- as.numeric(quantile(hab_values, 0.9, na.rm = TRUE)) 
q10 <- as.numeric(quantile(hab_values, 0.1, na.rm = TRUE)) 
q50 <-as.numeric(quantile(hab_values, 0.5, na.rm = TRUE)) 


hab_levels <- c(q10,
                q50,
                q90)
#takes 10th percentile


predict_data <- expand.grid(
  scaleEdgeDens2000 = edge_values,
  scaleHabAmount2000 = hab_levels
)

# Add other covariates
predict_data <- predict_data %>%  mutate(
  scaleCoastDist = mean(siteCovs(analysisData)$scaleCoastDist, na.rm = TRUE),
  scaleHabAmount100 = mean(siteCovs(analysisData)$scaleHabAmount100, na.rm = TRUE),
  scaleEdgeDens100 = mean(siteCovs(analysisData)$scaleEdgeDens100, na.rm = TRUE),
  PC1_t1 = mean(siteCovs(analysisData)$PC1_t1, na.rm = TRUE),
  scaleDoy = doy_means$mean_scaleDoy,
  scaleDoy2 = doy_means$mean_scaleDoy2
)



# Step 3: Predict Occupancy with standard errors (for error ribbon)
predictions <- predict(
  model,  
  newdata = predict_data,
  type = "state", 
  se.fit = TRUE  # Obtain standard errors for predictions
)

# Add predicted occupancy and standard errors to the data
predict_data$Occupancy <- predictions$Predicted
predict_data$SE <- predictions$SE

# Calculate 95% confidence intervals for occupancy (using 1.96 * SE for a 95% CI)
predict_data$LowerCI <- predict_data$Occupancy - 1.96 * predict_data$SE
predict_data$UpperCI <- predict_data$Occupancy + 1.96 * predict_data$SE

# Get mean and SD of the original (unscaled) Edge Density variable
edge_mean <- mean(siteCovs(analysisData)$edgeRook_2000_40, na.rm = TRUE)
edge_sd <- sd(siteCovs(analysisData)$edgeRook_2000_40, na.rm = TRUE)

# Add unscaled values to the dataset
predict_data$UnscaledEdgeDens2000 <- (predict_data$scaleEdgeDens2000 * edge_sd) + edge_mean



predict_data$habAmount <- factor(
  predict_data$scaleHabAmount2000,
  levels = hab_levels,
  labels = c("Low habitat", 
             "Medium habitat", 
             "High habitat")
)

# Plot figure with error ribbons ####

# Define a  color palette with muted tones
nature_palette <- c("Low habitat" = "#D55E00",   # Grey for "Bad" condition
                    "Medium habitat " = "#A0A0A0", # Lighter grey for "Neutral"
                    "High habitat" = "#56B4E9")   # Muted blue for "Good" condition

ocean_murrelet_hab_occ <- 
  ggplot(predict_data, aes(x = UnscaledEdgeDens2000, y = Occupancy, color = habAmount)) +
  geom_line(size = 1.5) +  # Thicker lines for clearer visibility
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = habAmount), 
              alpha = 0.3, color = NA) +  # Add error ribbons
  scale_color_manual(
    values = nature_palette,  # Apply the Nature-like muted color palette
    name = "habAmount"
  ) +
  scale_fill_manual(
    values = nature_palette,  # Apply same palette to ribbon fill
    name = "habAmount"
  ) +
  labs(
    x = "Proportion of Edge (in 2000m buffer)",  # Clear x-axis label
    y = "Predicted Occupancy Probability",  # Clear y-axis label
  ) +
  theme_minimal(base_size = 18) +  # Minimal theme with larger base font size
  theme_minimal(base_size = 18) +  # Minimal theme with larger base font size
  theme(
    legend.position = "top",  # Position the legend at the top
    legend.title = element_blank(),  # Increase legend title size for clarity
    legend.text = element_text(size = 14),  # Increase legend text size for better readability
    axis.title = element_text(size = 16),  # Increase axis title font size
    axis.text = element_text(size = 14),  # Increase axis label font size
    panel.grid.major = element_blank(),  # No major gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a thin border around the plot
 )
ocean_murrelet_hab_occ


#Save outputs #### 

# Save the plot using ggsave
ggsave(
  filename = "Figures/murrelet_by_ocean_and_edge.png",               # File path and name
  plot = ocean_murrelet_by_ocean_and_edge,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

# Save the plot using ggsave
ggsave(
  filename = "Figures/murrelet_occ_by_varrying_dist_and_ocean_condition.png",               # File path and name
  plot = ocean_murrelet_hab_occ_varrying_dist,          
  width = 16,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)


# Save the plot using ggsave
ggsave(
  filename = "Figures/occupancy_by_edge_and_habitat_amount.png",               # File path and name
  plot = ocean_murrelet_hab_occ,          
  width = 16,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)
