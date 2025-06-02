#Fit occupancy models 

# Load required libraries for data manipulation and occupancy modeling
library(tidyverse)
library(unmarked)
library(terra)
library(ggcorrplot)
library(AICcmodavg)
library(flextable)
library(officer)


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
#Fit occupancy models 


generate_starts <- function(previous_model, new_terms, new_terms_starts) {
  # Extract coefficients from the previous model
  prev_coefs <- coef(previous_model)
  
  # Ensure new terms and their starts are the same length
  if (length(new_terms) != length(new_terms_starts)) {
    stop("Length of new_terms must match the length of new_terms_starts.")
  }
  
  # Create named vector for new terms
  new_coefs <- setNames(new_terms_starts, new_terms)
  
  # Combine previous coefficients and new initialized values
  c(prev_coefs, new_coefs)
}

# Model 1: Simple Detection Model
# --------------------------------
simple_detection_model <- occu(
  formula = ~ ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + 
    scaleDoy + scaleDoy2 ~ 1,
  data = analysisData
)

plogis(-1.7) # baseline occupancy is about 15% (need to transform back to native scale; as outcome is on logit scale)


# Model 2: Adding Distance to Coast
# -----------------------------------
new_terms_dist <- c("PC1_t1", "scaleCoastDist")
#set vaguely informative starts based on hypothesis
new_terms_starts <- c(0, 0) #hypothesis; positive effect of good ocean year; negative effect of dist to coast

starts_dist_model <- generate_starts(
  previous_model = simple_detection_model, 
  new_terms = new_terms_dist,
  new_terms_starts = new_terms_starts 
)

dist_model <- occu(
  formula = ~ ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + 
    scaleDoy + scaleDoy2 ~ PC1_t1 + scaleCoastDist, 
  data = analysisData, 
  starts = starts_dist_model
)
coef(dist_model)


# Model 3: Adding Habitat and Edge Metrics
# ----------------------------------------
new_terms_habitat <- c("scaleHabAmount100", "scaleEdgeDens100", 
                       "scaleHabAmount2000", "scaleEdgeDens2000")
new_terms_starts <- c(0, 0, 0, 0) #hypothesis; negative effect of landscape-level edge; all else-positive 

starts_habitat <- generate_starts(
  previous_model = dist_model, 
  new_terms = new_terms_habitat, 
  new_terms_starts = new_terms_starts
)

model_with_habitat <- occu(
  formula = ~ ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + 
    scaleDoy + scaleDoy2 ~ PC1_t1 + scaleCoastDist + scaleHabAmount100 + 
    scaleEdgeDens100 + scaleHabAmount2000 + scaleEdgeDens2000, 
  data = analysisData, 
  starts = starts_habitat
)

# Model 4: interaction between edge and amount
# ----------------------------------------
new_terms_edge_amount <- c("scaleHabAmount100 * scaleEdgeDens100",
                           "scaleHabAmount2000 * scaleEdgeDens2000")
new_terms_starts <- c(0,0)
starts_edge_amount <- generate_starts(
  previous_model = model_with_habitat, 
  new_terms = new_terms_edge_amount, 
  new_terms_starts = new_terms_starts
)

##TEST: 
model_edge_amount_interactions <- occu(
  formula = ~ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + scaleDoy + scaleDoy2 ~ 
    PC1_t1 + scaleCoastDist + scaleHabAmount100 + scaleEdgeDens100 + scaleHabAmount100 * scaleEdgeDens100 + 
    scaleHabAmount2000 + scaleEdgeDens2000 + scaleHabAmount2000 * scaleEdgeDens2000,
  data = analysisData,
  starts = starts_edge_amount)


# Model 5: Adding interaction between PC1 and 2km edge density 
# -----------------------------------------
new_terms_pc1_interactions <- c("scaleEdgeDens2000:PC1_t1")
new_terms_starts <- c(0) 
starts_pc1_interactions <- generate_starts(
  previous_model = model_with_habitat, 
  new_terms = new_terms_pc1_interactions,
  new_terms_starts = new_terms_starts
)

length(starts_pc1_interactions)
pc1_interaction_model <- occu(
  formula = ~ ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + 
    scaleDoy + scaleDoy2 ~ PC1_t1 + scaleCoastDist + scaleHabAmount100 + 
    scaleEdgeDens100  + scaleHabAmount2000 + scaleEdgeDens2000 + 
    scaleEdgeDens2000*PC1_t1, 
  data = analysisData, 
  starts = starts_pc1_interactions
)

# Model 6: Adding interaction between PC1 and 2km edge density & PC1_t1:scaleCoastDist

new_terms_PC1_two_way <- c("scaleEdgeDens2000:PC1_t1", "PC1_t1:scaleCoastDist")
new_terms_starts <- c(0, 0) 
starts_PC1_two_way <- generate_starts(
  previous_model = model_with_habitat, 
  new_terms = new_terms_PC1_two_way,
  new_terms_starts = new_terms_starts
)

PC1_two_way <-  occu(
  formula = ~ ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + 
    scaleDoy + scaleDoy2 ~ PC1_t1 + scaleCoastDist + scaleHabAmount100 + 
    scaleEdgeDens100  + scaleHabAmount2000 + scaleEdgeDens2000 + 
    scaleEdgeDens2000*PC1_t1 + PC1_t1*scaleCoastDist,
  data = analysisData, 
  starts = starts_PC1_two_way
)


# Model 7: Adding interaction between PC1 and 2km edge, PC1_t1:scaleCoastDist & scaleCoastDist:scaleEdgeDens2000
new_terms_multiple_two_way <- c("PC1_t1:scaleCoastDist", "scaleCoastDist:scaleEdgeDens2000","PC1_t1:scaleEdgeDens2000")
new_terms_starts <- c(0,0,0) 
starts__multiple_two_way <- generate_starts(
  previous_model = model_with_habitat, 
  new_terms = new_terms_multiple_two_way,
  new_terms_starts = new_terms_starts
)

multiple_two_way <-  occu(
  formula = ~ ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + 
    scaleDoy + scaleDoy2 ~ PC1_t1 + scaleCoastDist + scaleHabAmount100 + 
    scaleEdgeDens100  + scaleHabAmount2000 + scaleEdgeDens2000 + 
    scaleEdgeDens2000*PC1_t1 +scaleCoastDist*scaleEdgeDens2000 + PC1_t1*scaleCoastDist,
  data = analysisData, 
  starts = starts__multiple_two_way
)

# Model 8: Adding a three-way interaction between scaleCoastDist * scaleEdgeDens2000 * PC1_t1, 

# -----------------------------------
new_terms_threeway <- c("scaleEdgeDens2000:PC1_t1",
                        "scaleCoastDist:scaleEdgeDens2000",
                        "PC1_t1:scaleCoastDist",
                        "scaleEdgeDens2000:PC1_t1:scaleCoastDist")

new_terms_starts <- c(0, 0, 0, 0)
starts_threeway <- generate_starts(
  previous_model = model_with_habitat, 
  new_terms = new_terms_threeway,
  new_terms_starts = new_terms_starts
)

threeway_interaction_model <- occu(
  formula = ~ ownership + scaleCanopy100 + scaleConDens100 + scaleEdgeDens100 + 
    scaleDoy + scaleDoy2 ~
    PC1_t1 + scaleCoastDist +
    scaleHabAmount100 +   scaleEdgeDens100 + 
    scaleHabAmount100 + scaleHabAmount2000+
    scaleCoastDist * scaleEdgeDens2000 * PC1_t1, 
  data = analysisData, 
  starts = starts_threeway
)

#Compare model performance 

model_list_all <- fitList(simple_detection_model,dist_model,model_with_habitat, model_edge_amount_interactions, 
                          pc1_interaction_model,PC1_two_way, multiple_two_way,threeway_interaction_model)

mod_performance <- modSel(model_list_all)

#------------------------------------------------------------------
#Compare by K-fold cross validation for the 4 best-performing models
#------------------------------------------------------------------


# Define the k-fold cross-validation 
k_fold_mod1 <- crossVal(
  object = simple_detection_model,        #  fitted model
  method = "Kfold",                       # Specify k-fold validation
  folds = 10,                             # Number of folds (can adjust as needed)
  statistic = unmarked:::RMSE_MAE ,       # Use default RMSE and MAE statistics
  parallel = FALSE)

k_fold_mod2 <- crossVal(
  object = dist_model,    
  method = "Kfold",                      
  folds = 10,                             
  statistic = unmarked:::RMSE_MAE ,      
  parallel = FALSE)

k_fold_mod3 <- crossVal(
  object = model_with_habitat,    
  method = "Kfold",                      
  folds = 10,                             
  statistic = unmarked:::RMSE_MAE ,      
  parallel = FALSE)

k_fold_mod4 <- crossVal(
  object = model_edge_amount_interactions,   
  method = "Kfold",                       
  folds = 10,                            
  statistic = unmarked:::RMSE_MAE ,       
  parallel = FALSE)

k_fold_mod5 <- crossVal(
  object = pc1_interaction_model,    
  method = "Kfold",                      
  folds = 10,                             
  statistic = unmarked:::RMSE_MAE ,       
  parallel = FALSE)

k_fold_mod6 <- crossVal(
  object = PC1_two_way,    
  method = "Kfold",                       
  folds = 10,                          
  statistic = unmarked:::RMSE_MAE ,       
  parallel = FALSE)

k_fold_mod7 <- crossVal(
  object = multiple_two_way,    
  method = "Kfold",                       
  folds = 10,                            
  statistic = unmarked:::RMSE_MAE ,       
  parallel = FALSE)

k_fold_mod8 <- crossVal(
  object = threeway_interaction_model,    
  method = "Kfold",                       
  folds = 10,                             
  statistic = unmarked:::RMSE_MAE ,      
  parallel = FALSE)

 kfold_list <- list(k_fold_mod1,k_fold_mod2, k_fold_mod3, k_fold_mod4, k_fold_mod5,k_fold_mod6, k_fold_mod7, k_fold_mod8)
 saveRDS(kfold_list, "Models/Kfold_all_model_performance.rds")


#define the kfold cross-validation for 1-way interaction model

# 
# 
# # Define the k-fold cross-validation Model 3 
# k_fold_results_modelwithhabitat <- crossVal(
#   object = model_with_habitat,    # Your fitted model
#   method = "Kfold",                       # Specify k-fold validation
#   folds = 10,                             # Number of folds (can adjust as needed)
#   statistic = unmarked:::RMSE_MAE ,       # Use default RMSE and MAE statistics
#   parallel = FALSE)
# 
# # Define the k-fold cross-validation Model 5
# 
# k_fold_results_pc1 <- crossVal(
#   object = pc1_interaction_model,    
#   method = "Kfold",                 
#   folds = 10,                      
#   statistic = unmarked:::RMSE_MAE ,                  
#   parallel = FALSE)
# 
# # Define the k-fold cross-validation Model 7 
# 
# # Define the k-fold cross-validation
# k_fold_results_multiple_two_way <- crossVal(
#   object = multiple_two_way,    
#   method = "Kfold",                       
#   folds = 10,                             
#   statistic = unmarked:::RMSE_MAE ,                   
#   parallel = FALSE)
# 
# 
# # Define the k-fold cross-validation model 8 
# k_fold_results_threeway_interaction <- crossVal(
#   object = threeway_interaction_model,    
#   method = "Kfold",                      
#   folds = 10,                             
#   statistic = unmarked:::RMSE_MAE ,                 
#   parallel = FALSE)
# 
# # View results - seems like there is not much difference in model performance in terms of RMSE and MAE
# print(k_fold_results_modelwithhabitat)
# print(k_fold_results_pc1)
# print(k_fold_results_multiple_two_way)
# print(k_fold_results_threeway_interaction)
# 
# kfold_list <- list(k_fold_results_modelwithhabitat,k_fold_results_pc1,k_fold_results_multiple_two_way,k_fold_results_threeway_interaction)

# #Save outputs
# #save K-fold cross-validation performance 
# saveRDS(kfold_list, "Models/Kfold_model_performance.rds")
# kfold_list <- readRDS("Models/Kfold_model_performance.rds")
# print(kfold_list)


#Save best model ####
saveRDS(multiple_two_way, "Models/final_model_5thMay2025.rds" )

#---------------------------------------------------------
#Make tables of model performance 
#---------------------------------------------------------
#NB:
#simple_detection = model 1 
# dist_model = model 2 
#model_with_habitat = model 3
#model_edge_amount_interactions = model 4 
#pc1_interaction_model = model 5
#PC1_two_way = model 6
#multiple_two_way = model 7 
#threeway_interaction_model = model 8

# make a table of model coefficients####
models
# Extract and clean up the model comparison table
model_table <- mod_performance@Full %>%
  select(model, nPars, AIC, delta, AICwt, cumltvWt) %>%
  mutate(
    AIC = round(AIC, 2),
    delta = round(delta, 2),
    AICwt = ifelse(AICwt < 0.001, "<0.001", round(AICwt, 3)),
    cumltvWt = round(cumltvWt, 2)
  ) %>%
  rename(
    Model = model,
    `No. Parameters` = nPars,
    AIC = AIC,
    `ΔAIC` = delta,
    `AIC Weight` = AICwt,
    `Cumulative Weight` = cumltvWt
  ) %>%  
  mutate(Model = recode(Model,
                        "simple_detection_model" = "Model 1",
                        "dist_model" = "Model 2",
                        "model_with_habitat" = "Model 3",
                        "model_edge_amount_interactions" = "Model 4",
                        "pc1_interaction_model" = "Model 5",
                        "PC1_two_way" = "Model 6",
                        "multiple_two_way" = "Model 7",
                        "threeway_interaction_model" = "Model 8"))


# Create a formatted flextable
ft <- flextable(model_table)
ft <- autofit(ft)
ft <- set_caption(ft, "Table 1. Model comparison table showing AIC-based support for candidate occupancy models.")

# Export to Word a table of model comparison
doc <- read_docx() %>%
  body_add_par("Model Selection Table", style = "heading 1") %>%
  body_add_flextable(ft)

#export model comparison
#print(doc, target = "Models/Tables/model_comparison_table.docx")

#------------------------------------------
# Print parametres for the best performing model 
#------------------------------------------

model <- readRDS("Models/final_model_5thMay2025.rds")

# Extract coefficients

# Helper function to extract, rename, tag, and compute 95% CI
extract_component <- function(summary_component, component_name) {
  as.data.frame(summary_component) %>%
    rownames_to_column(var = "Term") %>%
    rename(
      Estimate = Estimate,
      StdError = SE,
      z = z,
      PValue = 'P(>|z|)'
    ) %>%
    mutate(
      Component = component_name,
      CI_lower = Estimate - 1.96 * StdError,
      CI_upper = Estimate + 1.96 * StdError
    ) %>%
    select(Term, Estimate, StdError, CI_lower, CI_upper, z, PValue,Component)
}
# Extract, combine, and format
model_summary <- summary(model)
final_table <- bind_rows(
  extract_component(model_summary$state, "Occupancy"),
  extract_component(model_summary$det, "Detection")
) %>%
  mutate(
    Estimate = round(Estimate, 3),
    StdError = round(StdError, 3),
    z = round(z, 3),
    CI_lower = round(CI_lower, 3),
    CI_upper = round(CI_upper, 3),
    PValue = case_when(
      PValue < 0.001 ~ "< 0.001",
      TRUE ~ format(round(PValue, 3), nsmall = 3)
  ))
final_table

# Create a formatted flextable
ft2 <- flextable(final_table) %>%
  set_table_properties(width = .95, layout = "autofit") %>%
  fontsize(size = 9, part = "all") %>%
  autofit() %>%
  set_caption("Table 2. Parameter estimates from the top-ranked occupancy model, including 95% confidence intervals and p-values.")

doc2 <- read_docx() %>%
  body_add_par("Best-performing Model Coefficients Table", style = "heading 1") %>%
  body_add_flextable(ft2)

#Export best performing table
#print(doc2, target = "Models/Tables/best_model_coefficients_table.docx")



