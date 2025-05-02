#Fit occupancy models 

# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)
library(ggcorrplot)
library(AICcmodavg)

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

plogis(-1.7) # baseline occupancy is about 15% (need to transform back to native scale; as outcome is on logit scale)

# Step 3: Adding a Site Covariate - Distance to Coast
# ---------------------------------------------------
# Introduce `scaleCoastDist` and PC1 in year t-1 as a site covariate to explain occupancy.
dist_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 # detection
  ~ PC1_t1 + scaleCoastDist,  # occupancy
                                                                       
  data = analysisData,
  starts = c(
    coef(simple_detection_model)[1],  # Occupancy intercept - we use baseline occ from prev model
    0,0,                                 # Placeholder for new occupancy covariate (`scaleCoastDist`)
    coef(simple_detection_model)[2:10] # Detection coefficients for the 9 detection coefficients 
  ))
  

# Step 5: introduce Habitat Amount and Edge Metrics at Multiple Scales
# ------------------------------------------------------------
# Add site-level habitat amount and edge density metrics at both 100m and 2000m scales.
start_time <- Sys.time()
model_with_habitat <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    PC1_t1 + scaleCoastDist  + scaleHabAmount100 + scaleEdgeDens100 + scaleHabAmount2000 + scaleEdgeDens2000,
  data = analysisData,
  starts = c(
    coef(dist_model)[1:3],
    rep(0, 4),
    coef(dist_model)[4:12])
)
print(Sys.time() - start_time)

# Save updated results
#save(list = ls(), file = 'Models/ManuscriptResults.RData')

# Step 6: Adding Habitat-Edge Interactions
# -----------------------------------------
# Include interaction terms between habitat amount and edge density metrics at both scales.
start_time <- Sys.time()
model_with_interactions <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    PC1_t1 + scaleCoastDist + scaleHabAmount100 + scaleEdgeDens100 + scaleHabAmount100 * scaleEdgeDens100 + 
    scaleHabAmount2000 + scaleEdgeDens2000 + scaleHabAmount2000 * scaleEdgeDens2000,
  data = analysisData,
  starts = c(coef(model_with_habitat)[1:7],
          rep(0,2),
          coef(model_with_habitat)[8:16]))

print(Sys.time() - start_time)

# Save results for interaction model
#save(list = ls(), file = 'Models/ManuscriptResults.RData')
# 
# # Step 7: Adding PC1 and Edge Interactions
# # ------------------------------------------------------
# Explore interactions between `scaleCoastDist` and edge density metrics - how important is distance in modulating
#species' responses to edges .
start_time <- Sys.time()
coastal_interaction_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~
    PC1_t1 + scaleCoastDist  + scaleHabAmount100 + scaleEdgeDens100 + scaleCoastDist * scaleEdgeDens100 +
    scaleHabAmount2000 + scaleEdgeDens2000 + scaleCoastDist * scaleEdgeDens2000,
  data = analysisData,
  starts = coef(model_with_interactions)  # Use previous model's coefficients
)
print(Sys.time() - start_time)

# Save results
#save(list = ls(), file = 'Models/ManuscriptResults.RData')
print(coastal_interaction_model)
# Step 8: Adding Interactions with PCA1
# ---------------------------------------
# Incorporate `PC1` as a site covariate and include interaction terms with edge density.
start_time <- Sys.time()
pc1_interaction_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    PC1_t1 + scaleCoastDist  + 
    scaleHabAmount100 + scaleEdgeDens100 +
    scaleHabAmount2000 + scaleEdgeDens2000 + 
    scaleEdgeDens100 * PC1_t1+
    PC1_t1 * scaleEdgeDens2000,
  data = analysisData,
  starts = c(coef(coastal_interaction_model))
)
print(Sys.time() - start_time)


#Step 9: Adding a 3 way interaction (distance*PC1*edge2000) 
length(coef(pc1_interaction_model))

start_time <- Sys.time()
multiple_interaction_model <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    PC1_t1 + scaleCoastDist  + 
    scaleHabAmount100 + scaleHabAmount2000 + 
    scaleEdgeDens100 + scaleEdgeDens2000 +
    scaleCoastDist * scaleEdgeDens2000 * PC1_t1,
  data = analysisData,
  starts = c(coef(pc1_interaction_model),0,0)
)
print(Sys.time() - start_time)

# Save results
#save(list = ls(), file = 'Models/ManuscriptResults.RData')

#test for overfitting 


# which is the best model?
model_list <- fitList(model_with_interactions, multiple_interaction_model, pc1_interaction_model)

model_list_all <- fitList(simple_detection_model,dist_model,model_with_habitat, model_with_interactions,
                          coastal_interaction_model,
                          pc1_interaction_model,multiple_interaction_model)

model_list_subs <- fitList(pc1_interaction_model, coastal_interaction_model, multiple_interaction_model) 

modSel(model_list)
modSel(model_list_all)
modSel(model_list_subs)

multiple_interaction_model
WHY_multiple_interaction_model <- readRDS("Models/pc1_3way_v2_interaction_model.rds")
WHY_pc1_interaction_model <- readRDS("Models/pc1_interaction_model.rds")
WHY_pc1_interaction_model
pc1_interaction_model

save(list = ls(), file = 'Models/TESTManuscriptResults.RData')


#compare ALL models

#test for overfitting: 

#check I am not overfitting 
#define the kfold cross-validation for 1-way interaction model 

#kfolds coastal 
k_fold_results_coastal <- crossVal(
  object = coastal_interaction_model     ,    # Your fitted model
  method = "Kfold",                       # Specify k-fold validation
  folds = 10,                             # Number of folds (can adjust as needed)
  statistic = unmarked:::RMSE_MAE ,                   # Use default RMSE and MAE statistics
  parallel = FALSE)


# Define the k-fold cross-validation
k_fold_results_pc1 <- crossVal(
  object = pc1_interaction_model,    # Your fitted model
  method = "Kfold",                       # Specify k-fold validation
  folds = 10,                             # Number of folds (can adjust as needed)
  statistic = unmarked:::RMSE_MAE ,                   # Use default RMSE and MAE statistics
  parallel = FALSE)


# Define the k-fold cross-validation for 2-way model 
k_fold_results <- crossVal(
  object = multiple_interaction_model,    # Your fitted model
  method = "Kfold",                       # Specify k-fold validation
  folds = 10,                             # Number of folds (can adjust as needed)
  statistic = unmarked:::RMSE_MAE ,                   # Use default RMSE and MAE statistics
  parallel = FALSE)

# View results - seems like there is not much difference in model performance in terms of RMSE and MAE 
print(k_fold_results_coastal)
print(k_fold_results_pc1)
print(k_fold_results)


# Save final results
save(list = ls(), file = 'Models/ManuscriptResults.RData')
saveRDS(pc1_interaction_model,"Models/pc1_interaction_model.rds")
saveRDS(multiple_interaction_model,"Models/pc1_3way_v2_interaction_model.rds")


