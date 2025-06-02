library(tidyverse)
library(unmarked)
#Read in inputs 

#read in 10 points 
run_dataset <- read.csv("10_points_dusty.csv") %>%  
  arrange(point_id, decrease_factor) 
run_dataset <- run_dataset %>%
  mutate(across(where(is.character), as.factor))

str(run_dataset)
  
#read in model and function
source("Functions/se_diff_occu.R")#custom function for calculating % change and SE on logis scale
source("Functions/se_odds_ratio_occu.R")#custom function for calculating % change and SE on logis scale

model <- readRDS("Models/final_model_5thMay2025.rds")

#make predictions of fitted occupancies
occupancy_predictions <- 
  unmarked::predict(
    model,
    newdata = run_dataset,
    type = "state",
    se.fit = FALSE
  ) %>% rename(Occupancy = Predicted) %>%  
  cbind(run_dataset) %>%  
  arrange(point_id, decrease_factor) %>%  
  dplyr::select(point_id, decrease_factor, Occupancy,SE)

results <- list()

# Loop through unique point IDs
for (i in unique(run_dataset$point_id)) {
  # Filter the dataset for the current point ID
  current_data <- run_dataset %>% filter(point_id == i)
  
  # Apply the se_diff_occ function
  #estim <- se_diff_occ(
  estim <- se_odds_ratio_occ(
  
      newdat = current_data,
    fit = model,
    type = "state",
    order = c(2, 1)
  )
  
  # Store the result
  results[[as.character(i)]] <- estim
}

# Optionally, combine results into a single data frame if needed
final_results <- bind_rows(results, .id = "point_id")

#add percentage change to fitted occupancy estimates 
occupancy_predictions_pc <- occupancy_predictions %>%  
  left_join(final_results)


