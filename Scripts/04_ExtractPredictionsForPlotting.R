# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)
library(ggcorrplot)

#read in unmarked object
analysisData <- readRDS("Outputs/analysisDataUnmarked.rds")

#read in interaction model
model <- readRDS("Models/pc1_interaction_model.rds")

# Step 1: Extract site-level covariates from the unmarked object
# -------------------------------------------------------------
# Use siteCovs() to extract site-level covariates from the unmarkedFrameOccu object
site_covs <- siteCovs(analysisData)

# Check the structure of the site-level covariates
str(site_covs)

# Define the range of edge density values for predictions
edge_values <- seq(
  min(site_covs$scaleEdgeDens2000, na.rm = TRUE),
  max(site_covs$scaleEdgeDens2000, na.rm = TRUE),
  length.out = 100
)

# Define levels of ocean year conditions (PC1_t1) for good, neutral, and bad years
min(site_covs$PC1_t1) # worst year -1.822797 
max(site_covs$PC1_t1) # best year 1.647497

pc1_levels <- c(-1.82, 0, 1.64)  # Example: Low, Neutral, High

# Step 2: Create a data frame for predictions
# -------------------------------------------
# Use expand.grid() to create all combinations of edge density and PC1 levels
predict_data <- expand.grid(
  scaleEdgeDens2000 = edge_values,
  PC1_t1 = pc1_levels
)

# Add other covariates at their mean values (from the site-level covariates)
predict_data <- predict_data %>%
  mutate(
    scaleCoastDist = mean(site_covs$scaleCoastDist, na.rm = TRUE),
    scaleHabAmount100 = mean(site_covs$scaleHabAmount100, na.rm = TRUE),
    scaleEdgeDens100 = mean(site_covs$scaleEdgeDens100, na.rm = TRUE),
    scaleHabAmount2000 = mean(site_covs$scaleHabAmount2000, na.rm = TRUE),
    scaleDoy = mean(site_covs$scaleDoy, na.rm = TRUE),
    scaleDoy2 = mean(site_covs$scaleDoy^2, na.rm = TRUE)
  )

# Step 3: Predict occupancy probabilities using the unmarked model
# ----------------------------------------------------------------
# Use the predict() function to calculate predicted occupancy
predict_data$Occupancy <- predict(
  your_model,  # Replace 'your_model' with your fitted unmarked occupancy model
  newdata = predict_data,
  type = "state"  # 'state' indicates we're predicting occupancy probabilities
)$Predicted

# Add labels for ocean year conditions
predict_data$OceanYear <- factor(
  predict_data$PC1_t1,
  levels = pc1_levels,
  labels = c("Bad Ocean Years", "Neutral Ocean Years", "Good Ocean Years")
)

# Step 4: Plot the interaction between edge density and ocean year conditions
# ---------------------------------------------------------------------------
# Use ggplot2 to visualize the interaction effect
ggplot(predict_data, aes(x = scaleEdgeDens2000, y = Occupancy, color = OceanYear)) +
  geom_line(size = 1.2) +  # Plot occupancy predictions as lines
  scale_color_manual(
    values = c("red", "blue", "green"),  # Colors for bad, neutral, and good years
    name = "Ocean Conditions"  # Legend title
  ) +
  labs(
    x = "Edge Density (2000m, scaled)",  # X-axis label
    y = "Predicted Occupancy Probability",  # Y-axis label
    title = "Interaction of Edge Density and Ocean Conditions on Murrelet Occupancy",
    subtitle = "Predictions for Good, Neutral, and Bad Ocean Years"
  ) +
  theme_minimal(base_size = 15) +  # Minimalistic theme
  theme(
    legend.position = "top",  # Place legend on top
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
