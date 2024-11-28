# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)
library(ggcorrplot)

source("scripts/02_OrganiseMurreletData.R")
# Read in unmarked object
analysisData <- readRDS("Outputs/analysisDataUnmarked.rds")
# Read in interaction model
model <- readRDS("Models/pc1_interaction_model.rds")

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

minPC1 <- min(siteCovs(analysisData)$PC1_t1)
maxPC1 <- max(siteCovs(analysisData)$PC1_t1)

pc1_levels <- c(minPC1,
               # 0,
                maxPC1)  # Example levels bad, (medium), and good ocean year conditions
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
nature_palette <- c("Bad Ocean Years" = "#4C4C4C",   # Grey for "Bad" condition
                    #"Neutral Ocean Years" = "#A0A0A0", # Lighter grey for "Neutral"
                    "Good Ocean Years" = "#1F77B4")   # Muted blue for "Good" condition

ocean_murrelet_hab_occ <- 
  ggplot(predict_data, aes(x = scaleEdgeDens2000, y = Occupancy, color = OceanYear)) +
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
    x = "Edge Density (2000m, scaled)",  # Clear x-axis label
    y = "Predicted Occupancy Probability",  # Clear y-axis label
    title = "Effect of Edge Density and Ocean Conditions on Murrelet Occupancy"  # Concise title
  ) +
  theme_minimal(base_size = 18) +  # Minimal theme with larger base font size
  theme(
    legend.position = "top",  # Position the legend at the top
    legend.title = element_text(size = 16),  # Increase legend title size for clarity
    legend.text = element_text(size = 14),  # Increase legend text size for better readability
    axis.title = element_text(size = 16),  # Increase axis title font size
    axis.text = element_text(size = 14),  # Increase axis label font size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Bold title and centered
    plot.subtitle = element_text(size = 14, hjust = 0.5),  # Subtitle size
    panel.grid.major = element_blank(),  # No major gridlines
    panel.grid.minor = element_blank(),  # No minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a thin border around the plot
  )

#Save outputs #### 
# Define the output folder and file name
output_folder <- "Figures"
output_file <- file.path(output_folder, "murrelet_occupancy_plot.png")

# Create the folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Save the plot using ggsave
ggsave(
  filename = output_file,               # File path and name
  plot = ocean_murrelet_hab_occ,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

