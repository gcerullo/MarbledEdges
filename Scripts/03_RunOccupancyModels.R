#Fit occupancy models 

# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)
library(ggcorrplot)

#read in unmarked object
analysisData <- readRDS("Outputs/analysisDataUnmarked.rds")

# --- OCCUPANCY MODELING WITH INTERACTIONS ---

# Step 1: Presence/Absence Data Preparation
# ------------------------------------------
# Summarize detection data by calculating the maximum detection (1 = present, 0 = absent)
# across surveys per site. Use the `apply` function to check detection for each site.
presence_absence_summary <- table(apply(analysisData@y, 1, FUN = 'max', na.rm = TRUE))
print(presence_absence_summary)

# --- BEGIN OCCUPANCY MODELING ---
# The `occu` function estimates two probabilities:
# 1. Detection probability (p): Probability of detecting a species if present.
# 2. Occupancy probability (ψ): Probability that a site is occupied.
      #Structure is as follows: ~ detection predictors ~ occupancy predictors
      #The first ~: Defines covariates for detection probability (p).
      #The second ~: Defines covariates for occupancy probability (ψ).

# Step 2: Simple Detection Model
# --------------------------------
# This model only includes detection covariates. Site occupancy is assumed constant (~1).
start_time <- Sys.time()
simple_detection_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 1,
  data = analysisData
)
print(Sys.time() - start_time)  # Print runtime

plogis(-1.6) # baseline occupancy is about 16% (need to transform back to native scale; as outcome is on logit scale)

# Step 3: Adding a Site Covariate - Distance to Coast
# ---------------------------------------------------
# Introduce `PC1` as a site covariate to explain occupancy.
dist_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2  ~ scaleCoastDist, # detection
                                                                        # occupancy
  data = analysisData,
  starts <- c(
    coef(simple_detection_model)[1],    # Occupancy intercept - we use baseline occ from prev model
    0,                                 # Placeholder for new occupancy covariate (`scaleCoastDist`)
    coef(simple_detection_model)[2:10] # Detection coefficients for the 9 
  ))
  
  
 # starts = c(coef(simple_detection_model)[1], 0, coef(simple_detection_model)[2:11]))

  starts = c(-1.6, 0, coef(simple_detection_model)[1:9]))  # Provide 11 starting values

#-1.6: Initial value for the occupancy intercept (from simple_detection_model).
#0: Initial value for scaleCoastDist (new parameter for occupancy).
#coef(simple_detection_model)[1:9]: Coefficients for detection (9 parameters from simple_detection_model).

print(Sys.time() - start_time)


# Step 4: Adding Year as a Factor
# --------------------------------
# Include year as a categorical covariate to account for temporal effects on occupancy.
start_time <- Sys.time()
year_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2
  ~ scaleCoastDist + as.factor(year),
  data = analysisData,
  starts = c(
    coef(dist_model)[1:2],  # Existing occupancy coefficients
    rep(0, 27),            # Starting values for year dummy variables
    coef(dist_model)[3:11]  # Detection coefficients
  )
)
print(Sys.time() - start_time)

# Save intermediate results to a file for reproducibility
save(list = ls(), file = 'Models/ManuscriptResults.RData')

# Step 5: Habitat Amount and Edge Metrics at Multiple Scales
# ------------------------------------------------------------
# Add site-level habitat amount and edge density metrics at both 100m and 2000m scales.
start_time <- Sys.time()
model_with_habitat <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    scaleCoastDist + as.factor(year) + scaleHabAmount100 + scaleEdgeDens100 + scaleHabAmount2000 + scaleEdgeDens2000,
  data = analysisData,
  starts = c(coef(year_model)[1:30], rep(0, 4), coef(year_model)[31:40])
)
print(Sys.time() - start_time)

# Save updated results
save(list = ls(), file = 'Models/ManuscriptResults.RData')

# Step 6: Adding Habitat-Edge Interactions
# -----------------------------------------
# Include interaction terms between habitat amount and edge density metrics at both scales.
start_time <- Sys.time()
model_with_interactions <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    scaleCoastDist + as.factor(year) + scaleHabAmount100 + scaleEdgeDens100 + scaleHabAmount100 * scaleEdgeDens100 + 
    scaleHabAmount2000 + scaleEdgeDens2000 + scaleHabAmount2000 * scaleEdgeDens2000,
  data = analysisData,
  starts = c(coef(model_with_habitat)[1:32], 0, coef(model_with_habitat)[33:34], 0, coef(model_with_habitat)[35:44])
)
print(Sys.time() - start_time)

# Save results for interaction model
save(list = ls(), file = 'Models/ManuscriptResults.RData')

# Step 7: Adding Coastal Distance and Edge Interactions
# ------------------------------------------------------
# Explore interactions between `scaleCoastDist` and edge density metrics.
start_time <- Sys.time()
coastal_interaction_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    scaleCoastDist + as.factor(year) + scaleHabAmount100 + scaleEdgeDens100 + scaleCoastDist * scaleEdgeDens100 + 
    scaleHabAmount2000 + scaleEdgeDens2000 + scaleCoastDist * scaleEdgeDens2000,
  data = analysisData,
  starts = coef(model_with_interactions)  # Use previous model's coefficients
)
print(Sys.time() - start_time)

# Save results
save(list = ls(), file = 'Models/ManuscriptResults.RData')

# Step 8: Adding Interactions with PCA1
# ---------------------------------------
# Incorporate `PC1` as a site covariate and include interaction terms with edge density.
start_time <- Sys.time()
pca_interaction_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    scaleCoastDist + as.factor(year) + scaleHabAmount100 + scaleEdgeDens100 + PC1 + 
    scaleEdgeDens100 * PC1 + scaleHabAmount2000 + scaleEdgeDens2000 + PC1 * scaleEdgeDens2000,
  data = analysisData,
  starts = c(coef(coastal_interaction_model)[1:32], 0, coef(coastal_interaction_model)[33:44])
)
print(Sys.time() - start_time)

# Save final results
save(list = ls(), file = 'Models/ManuscriptResults.RData')
