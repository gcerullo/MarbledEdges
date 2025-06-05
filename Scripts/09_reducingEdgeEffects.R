
library(terra)
library(sf)
library(exactextractr)
library(data.table)
library(tidyverse)
library(mgcv)
library(furrr)
library(parallel)
library(skimr)
library(summarytools)
library(ggh4x) 
source("scripts/02_OrganiseMurreletData.R")
source("Functions/se_diff_occu.R")#custom function for calculating % change and SE on logis scale
source("Functions/se_odds_ratio_occu.R")#calculate log odds ration 
source("Functions/raw_diff_occu.R")#calculate log odds ration 

#----------------------------------------------------------------
# read in inputs #####
covariates <- readRDS("Outputs/ScaledCovariates.rds") %>%   
  dplyr::select(PC1_t1, scaleDoy, scaleDoy2, scaleDoy2, OceanYear) %>%  unique() %>%  
  #Lets assume the worst PC1 (2016) conditions on record; #-1.880589
  mutate(PC1_t1 = case_when(
    OceanYear == "Bad Ocean Years" ~ -1.880589,
    TRUE ~ PC1_t1  # keeps the original value of PC1_t1 when OceanYear is not "Bad Ocean Years"
  ))
model <- readRDS("Models/final_model_5thMay2025.rds")

final2020 <- readRDS("Outputs/PNW_2020_extracted_covars.rds") %>%   #read in starting occupancy for 2020 from scrippt 8
  as.data.frame() #%>%  

long_lat <- final2020 %>% dplyr::select(point_id, x, y)
#add in information on elevation at each point 
pt_elevation <- read.csv("Outputs/elevation_at_each_point.csv") 

#read in means and SDs of covariates from original model to enable coorrect scalings
meansAndSds

# #remind ourselves of what model has seen unscaled
# detected_sites <- analysisSurveys %>% dplyr::select(id, detected)
# analysisSites %>% 
#   
#   #-------------------------------------------
# #if we only want to consider sites where there was a detection
# left_join(detected_sites, by = "id") %>%  
#   filter(detected == 1) %>% 
#   #-------------------------------------------
# dplyr::select(coastDist,habAmountDich_100,habAmountDich_2000,edgeRook_100_40,edgeRook_2000_40) %>%  
#   pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>% 
#   ggplot( aes(x = Value)) +
#   geom_histogram(fill = "blue", color = "white", alpha = 0.7) +
#   facet_wrap(~ Variable, scales = "free") +
#   theme_minimal() +
#   labs(title = "Histograms of Dataframe Columns", x = "Value", y = "Frequency")

#-------------------------------------------------------------
#read in functions ####

#----------------------------------
#functions ####
#----------------------------------
#function for scaling data and making model predictions dataset ####

prep_and_predict_parallelised <- function(x) {
  # Prepare data for prediction with scaled covariates
  prediction_df <- x %>% mutate(
    scaleHabAmount100 = (habAmountDich_100  - meansAndSds$meanHabAmountDich100) /meansAndSds$sdHabAmountDich100,
    scaleHabAmount2000 = (habAmountDich_2000 - meansAndSds$meanHabAmountDich2000) /meansAndSds$sdHabAmountDich2000,
    scaleEdgeDens100 =  (edgeRook_100_40 - meansAndSds$meanEdgeDens100) /meansAndSds$sdEdgeRook100, 
    scaleEdgeDens2000 = (edgeRook_2000_40 - meansAndSds$meanEdgeDens2000)/meansAndSds$sdEdgeRook2000,
    scaleCoastDist = (distance_to_coastline - meansAndSds$meanCoastDist) /meansAndSds$sdCoastDist
  ) %>%  
    dplyr::select(
      point_id,
      scaleHabAmount100,
      scaleHabAmount2000,
      scaleEdgeDens100,
      scaleEdgeDens2000,
      scaleCoastDist,
      decrease_factor,
      cancov_2020
    ) %>%
    cross_join(covariates)
  
  # Function for predicting a chunk of data
  predict_chunk <- function(data_chunk) {
    predict(
      model,
      newdata = data_chunk,
      type = "state",
      se.fit = FALSE
    ) %>%
      rename(Occupancy = Predicted)
  }
  
  # Split data into manageable chunks
  chunk_size <- 350000  # Adjust based on memory and cores
  chunks <- split(prediction_df, (seq_len(nrow(prediction_df)) - 1) %/% chunk_size)
  
  # Set up parallel processing
  plan(multisession, workers = parallel::detectCores() - 1)  # Use all but one core
  
  # Predict occupancy in parallel
  predictions <- future_map_dfr(chunks, predict_chunk)
  
  # Combine predictions back with the original data
  prediction_df <- prediction_df %>%
    bind_cols(predictions) %>%
    mutate(
      lower_CI = Occupancy - 1.96 * SE,
      upr_CI = Occupancy + 1.96 * SE
    )
  
  # Add back in other information
  add_data_back <- x %>%
    dplyr::select(
      point_id,
      point_leve_habitat,
      ownership,
      distance_to_coastline
    ) %>%
    unique()
  
  prediction_df <- prediction_df %>%
    left_join(add_data_back, by = "point_id")
  
  return(prediction_df)
}

ownership_only_variation_plot <- function(x){
  #add categorical distance to coast info 
  frag_df <- x 
  
  summary_stats <- frag_df %>%
    group_by(ownership) %>%
    summarise(
      count = n_distinct(point_id),
      .groups = "drop") %>%
    filter(ownership %in% c("Federal", "Private Industrial","Private Conservation", "Private Non-industrial","Indian", "State"))
  

  #BOXPLOTS
  frag_df %>%
    mutate(
      decrease_factor = case_when(
        decrease_factor == 0~ 0,
        decrease_factor == 0.25 ~ 25,
        decrease_factor == 0.5 ~ 50,
        decrease_factor == 0.75 ~ 75,
        decrease_factor == 1 ~ 100,
        TRUE ~ NA_real_    )) %>%  # Default case for unexpected values
    filter(ownership %in% c("Federal", "Private Industrial","Private Conservation",  "Private Non-industrial","Indian", "State")) %>% 
    mutate(ownership = factor(ownership, levels = c("Private Conservation", "Federal", "State", "Indian","Private Non-industrial","Private Industrial"))) %>% 
    
    group_by(decrease_factor) %>% 
    ggplot(aes(x = as.factor(decrease_factor), y = Occupancy, fill = as.factor(OceanYear))) + 
    #geom_point(aes(color = OceanYear), size = 2, alpha = 0.02, position = position_jitterdodge(jitter.width = .25)) +  # Points with jitter
    geom_boxplot(color = "black", outlier.shape = NA) +  # Boxplot with black outline
    #geom_violin(color = "black", alpha = 0.8, draw_quantiles = c(0.25, 0.5,0.75)) + 
    # Custom fill colors for OceanYear
    scale_fill_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + 
    scale_color_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + # Color points
    
    
    geom_text(
      data = summary_stats,
      aes(x = Inf, y = Inf, label = paste0("n = ", count)),
      hjust = 1.2, vjust = 1.5, size = 4, inherit.aes = FALSE
    )+
    
    theme_classic(base_size = 14) +  # Nature-style theme
    ylim(0,0.7)+
    labs(
      x = "% Reduction in Edge Density",  # Label for x-axis
      y = "Occupancy"   # Label for y-axis
    ) +
    facet_wrap(~ ownership, ncol = 3) +  # Customize facet labels
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
}

plot_ownership_distribution <- function(data, x_var, x_axis_title, xlims =NULL, elev_thresh = NULL) {
  
  # Filter and categorize data
  filtered_data <- data %>%
    filter(!ownership %in% c("Unknown", NA)) %>%  
    filter(ownership %in% c("Federal", "Private Industrial","Private Conservation", "Indian", "Private Non-industrial", "State")) %>% 
    filter(cancov_2020 >-1) %>%  
    filter(dem30m < elev_thresh)  %>%  
    mutate(ownership = factor(ownership, levels = c("Private Conservation", "Federal", "State", "Indian","Private Non-industrial","Private Industrial")))  
  
  # Compute median for each facet group
  summary_stats <- filtered_data %>%
    group_by(ownership) %>%
    summarise(
      median_value = median(.data[[x_var]], na.rm = TRUE),
      count = n_distinct(point_id),
      .groups = "drop"
    )
  
  # Create histogram plot with facet and median lines
  ggplot(filtered_data, aes(x = .data[[x_var]])) +  
    geom_histogram(fill = "#4C9A2A", color = "black", alpha = 0.7) +  
    geom_vline(data = summary_stats, aes(xintercept = median_value), 
               linetype = "dashed", color = "red", linewidth = 1) +  
    
    geom_text(
      data = summary_stats,
      aes(x = Inf, y = Inf, label = paste0("n = ", count)),
      hjust = 1.2, vjust = 1.5, size = 5, inherit.aes = FALSE
    ) +  
    
    labs(
      x = x_axis_title,  
      y = "Frequency"
    ) +
    
    facet_wrap(~ ownership, scales = "free_y", ncol = 3, nrow = 2) +  
    theme_minimal(base_size = 16) +  
    theme(
      legend.position = "top",  # Position the legend at the top
      legend.title = element_blank(),  # Increase legend title size for clarity
      legend.text = element_text(size = 14),  # Increase legend text size for better readability
      axis.title = element_text(size = 16),  # Increase axis title font size
      axis.text = element_text(size = 14),  # Increase axis label font size
      panel.grid.major = element_blank(),  # No major gridlines
      panel.border = element_rect(color = "black", fill = NA, size = 1)
      )+
    # Add a thin border around the plot
    scale_x_continuous(expand = c(0, 0), 
                       limits = xlims,
                       labels = scales::label_number(accuracy = 0.1)  # Format to 1 decimal place
    )
}
#----------------------------------
#Change edge amount 
#----------------------------------

#Consider the effect of 
# Define the decrease percentages for edge amount for all ownerships 
decrease_edge <- c(0, 0.25, 0.5, 0.75, 1)


#FOR ALL EDGES 
# Apply the decreases conditionally for ALL edges
results_all_edges <- lapply(decrease_edge, function(decrease) {
  final2020 %>%
    mutate(
      edgeRook_100_40 = edgeRook_100_40 * (1 - decrease),
      edgeRook_2000_40 = edgeRook_2000_40 * (1 - decrease)
    ) %>%
    mutate(decrease_factor = decrease) # Add a column to track the reduction factor
})

# Combine all the results into one data frame for comparison
results_all_edges <- bind_rows(results_all_edges)
results_all_edges_processed <- results_all_edges %>%  prep_and_predict_parallelised() %>%  
                left_join(pt_elevation)

#saveRDS(results_all_edges_processed, "Outputs/ReduceLandscapeAndLocalEdges.rds" )


#FOR LANDSCAPE EDGES ONLY 

# Apply the decreases conditionally for LANDSCAPE EDGED ONLY 
results_landscape_edges <- lapply(decrease_edge, function(decrease) {
  final2020 %>%
    mutate(
      edgeRook_2000_40 = edgeRook_2000_40 * (1 - decrease)
    ) %>%
    mutate(decrease_factor = decrease) # Add a column to track the reduction factor
})

results_landscape_edges <- bind_rows(results_landscape_edges)
results_landscape_edges_processed <- results_landscape_edges %>%  prep_and_predict_parallelised() %>%  
  left_join(pt_elevation)
#saveRDS(results_landscape_edges_processed, "Outputs/ReduceLandscapeEdges.rds" )

#===========================================================
# #calculate uncertainty around percentage change in edge amount
 

final2020 <- readRDS("Outputs/PNW_2020_extracted_covars.rds")    #read in starting occupancy for 2020 from scrippt 8
  
results_all_edges <- readRDS("Outputs/ReduceLandscapeAndLocalEdges.rds") 
results_all_edges <- results_all_edges %>%  select(point_id, Occupancy, OceanYear, decrease_factor) %>%  
  filter(decrease_factor %in% c(0,0.5))
results_all_edges_good <- results_all_edges %>% filter(OceanYear == "Good Ocean Years") 
results_all_edges_bad <- results_all_edges %>% filter(OceanYear == "Bad Ocean Years") 


# # construct model matrix for the predicted occupancy
predict_df <- lapply(decrease_edge, function(decrease) {
  final2020 %>%
    mutate(
      edgeRook_100_40 = edgeRook_100_40 * (1 - decrease),
      edgeRook_2000_40 = edgeRook_2000_40 * (1 - decrease)
    ) %>%
    mutate(decrease_factor = decrease)  # Add a column to track the reduction factor
 
    })

predict_df <- bind_rows(predict_df)
predict_df <- predict_df %>% 
  filter(cancov_2020 > -1) %>% #remove non-forest habitat points
  filter( distance_to_coastline < 60000) %>% #remove points v far from coast 
  filter(dem30m < 600) %>%  #remove points above a given altitude
   mutate(
    scaleHabAmount100 = (habAmountDich_100  - meansAndSds$meanHabAmountDich100) /meansAndSds$sdHabAmountDich100,
    scaleHabAmount2000 = (habAmountDich_2000 - meansAndSds$meanHabAmountDich2000) /meansAndSds$sdHabAmountDich2000,
    scaleEdgeDens100 =  (edgeRook_100_40 - meansAndSds$meanEdgeDens100) /meansAndSds$sdEdgeRook100, 
    scaleEdgeDens2000 = (edgeRook_2000_40 - meansAndSds$meanEdgeDens2000)/meansAndSds$sdEdgeRook2000,
    scaleCoastDist = (distance_to_coastline - meansAndSds$meanCoastDist) /meansAndSds$sdCoastDist
  ) %>%  
  # dplyr::select(
  #   point_id,
  #   scaleHabAmount100,
  #   scaleHabAmount2000,
  #   scaleEdgeDens100,
  #   scaleEdgeDens2000,
  #   scaleCoastDist,
  #   decrease_factor,
  #   cancov_2020
  # ) %>%
  cross_join(covariates) %>%  
  mutate(ownership = sample(c("or", "wa"), size = n(), replace = TRUE)) %>% 
  mutate(scaleCanopy100 = meansAndSds$meanCanopy100) %>%  
  mutate(scaleConDens100 = meansAndSds$meanCanopy100) %>%  
  mutate(ownership = as.factor(ownership))


# custom se_diff_occ params; calculates percentage change for a 50% decrease in edge 
# custom se_log_odds_occ params; calculates log odds change for a 50% decrease in edge 
#' @param newdat A data frame containing the covariate values for prediction. Must contain the two rows we want to contrast.
#' @param fit A fitted model object of class `unmarkedFitOccu` (from the `unmarked` package).
#' @param type A character string indicating which component to analyze: `"state"` for occupancy or `"det"` for detection. Defaults to `"state"`.
#' @param order A numeric vector of length 2 specifying the row indices in `newdat` to compare. 
predict_df50 <- predict_df %>%  filter(decrease_factor %in% c(0,0.5))

predict_df50_good <- predict_df50 %>% filter(OceanYear == "Good Ocean Years") 
predict_df50_bad <- predict_df50 %>% filter(OceanYear == "Bad Ocean Years") 


#select one #HARD CODED DECISION 
dataset <- predict_df50_good
dataset <- predict_df50_bad

run_dataset <- dataset

# Get unique point IDs
point_ids <- unique(run_dataset$point_id)
results <- vector("list", length(point_ids))

results<- list()

# Loop through unique point IDs
for (point_id in unique(run_dataset$point_id)) {
  # Filter the dataset for the current point ID
  current_data <- run_dataset %>% filter(point_id == !!point_id)
  
  # Apply the se_log_odds_occ function
  estim <- se_odds_ratio_occ(
    newdat = current_data,
    fit = model,
    type = "state",
    order = c(2, 1)
  )
  
  # Store the result
  results[[as.character(point_id)]] <- estim
}

# Optionally, combine results into a single data frame if needed
final_results <- bind_rows(results, .id = "point_id")
final_results_good <- final_results %>% left_join(results_all_edges_good)
#saveRDS(final_results_good, "Outputs/good_years_50pcEdgeReduction_logOdds_SE.rds")

final_results_bad <- bind_rows(results, .id = "point_id")
final_results_bad <- final_results_bad %>% left_join(results_all_edges_bad)
#saveRDS(final_results_bad, "Outputs/bad_years_50pcEdgeReduction_logOdds_SE.rds")


#combine results ---------------------
results_good <- readRDS("Outputs/good_years_50pcEdgeReduction_logOdds_SE.rds")
results_bad <- readRDS("Outputs/bad_years_50pcEdgeReduction_logOdds_SE.rds")


ggplot(results_good, aes(x = estim)) +
  geom_histogram(fill = "steelblue", color = "black") +
  xlim(-2, 25) +
  theme_classic() +
  labs(
    x = "log odds change from baseline",
    y = "Frequency",
    title = "Histogram of Estimates in Good Years")

ggplot(results_good, aes(x = estim)) +
  geom_histogram(fill = "steelblue", color = "black") +
  xlim(-2, 25) +
  theme_classic() +
  labs(
    x = "log odds change from baseline",
    y = "Frequency",
    title = "Histogram of Estimates in Good Years")

results_bad <- results_bad %>%  mutate(OceanYear = "Bad Ocean Years") 
results_good <- results_good %>%  mutate(OceanYear = "Good Ocean Years") 
log_odds_se <- rbind(results_bad, results_good)
log_odds_se <- log_odds_se %>%  left_join(final2020)

#saveRDS(log_odds_se, "Outputs/log_odds_se_05edge_reduction.rds")


####
#For raw occupancy difference only 
# 
# #select one #HARD CODED DECISION 
# dataset <- predict_df50_good
# dataset <- predict_df50_bad
# 
# run_dataset <- dataset
# 
# # Get unique point IDs
# point_ids <- unique(run_dataset$point_id)
# results <- vector("list", length(point_ids))
# 
# results<- list()
# 
# # Loop through unique point IDs
# for (point_id in unique(run_dataset$point_id)) {
#   # Filter the dataset for the current point ID
#   current_data <- run_dataset %>% filter(point_id == !!point_id)
#   
#   # Apply the se_log_odds_occ function
#   estim <- se_raw_diff_occ(
#     newdat = current_data,
#     fit = model,
#     type = "state",
#     order = c(2, 1)
#   )
#   
#   # Store the result
#   results[[as.character(point_id)]] <- estim
# }
# 
# # Optionally, combine results into a single data frame if needed
# final_results <- bind_rows(results, .id = "point_id")
# final_results_good <- final_results %>% left_join(results_all_edges_good)
# #saveRDS(final_results_good, "Outputs/good_years_50pcEdgeReduction_raw_occ_diff_SE.rds")
# 
# final_results_bad <- bind_rows(results, .id = "point_id")
# final_results_bad <- final_results_bad %>% left_join(results_all_edges_bad)
# #saveRDS(final_results_bad, "Outputs/bad_years_50pcEdgeReduction_raw_occ_diff_SE.rds")
# results_bad <- final_results_bad %>%  mutate(OceanYear = "Bad Ocean Years") 
# results_good <- final_results_good %>%  mutate(OceanYear = "Good Ocean Years") 
# raw_diff_se <- rbind(results_bad, results_good)
# raw_diff_se <- raw_diff_se %>%  left_join(final2020)
# saveRDS(raw_diff_se, "Outputs/raw_diff_se_05edge_reduction.rds")


##


#=============================================================
#Plot the distribution of edge and habitat amount for different ownership classifications

results_all_edges <- readRDS("Outputs/ReduceLandscapeAndLocalEdges.rds")
results_landscape <- readRDS("Outputs/ReduceLandscapeEdges.rds")
hist(results_all_edges$Occupancy)

#scale back to original values
results_all_edges <- results_all_edges %>%
  mutate(habAmountDich_2000 = (scaleHabAmount2000 * meansAndSds$sdHabAmountDich2000) + meansAndSds$meanHabAmountDich2000) %>%  
  mutate(edgeRook_2000_40 = (scaleEdgeDens2000 * meansAndSds$sdEdgeRook2000) + meansAndSds$meanEdgeDens2000) %>%  
  mutate(habAmountDich_100 = (scaleHabAmount100 * meansAndSds$sdHabAmountDich100) + meansAndSds$meanHabAmountDich100) %>%  
  mutate(edgeRook_100_40 = (scaleEdgeDens100 * meansAndSds$sdEdgeRook100) + meansAndSds$meanEdgeDens100)


#decide initial filtering params
final_filt <- function(df, distance, height) {
  df %>%  filter(cancov_2020 > -1) %>% #remove non-forest habitat points
    filter( distance_to_coastline < distance) %>% #remove points v far from coast 
    filter(dem30m < height)  #remove points above a given altitude
} 

#plot for all edges Figure 3 
plot_all_edges <- results_all_edges %>%  final_filt(distance = 60000, #maximum detection in our dataframe
                                  height = 600) %>%  
  ownership_only_variation_plot()

test_all <- results_all_edges %>%  dplyr::select(point_id)
#plot for landscape edges 

plot_landscape_edges  <- results_landscape %>%  final_filt(distance = 60000, #maximum detection in our dataframe
                                  height = 600) %>%  
  ownership_only_variation_plot()


#now plot distribution of key edge and habitat covariates ownership #####
hab100Data <-  
  results_all_edges %>%  filter(decrease_factor == 0) %>%  
  filter(distance_to_coastline < 60000) %>%  
  filter(cancov_2020 >-1) %>%  
  plot_ownership_distribution(x_var = "habAmountDich_100",
                              x_axis_title = "Habitat amount in 100 m radius",
                              xlims =NULL,
                              elev_thresh = 600)

hab2000Data <-  
  results_all_edges %>%  filter(decrease_factor == 0) %>%  
  filter(distance_to_coastline < 60000) %>%  
  filter(cancov_2020 >-1) %>%  
  plot_ownership_distribution(x_var = "habAmountDich_2000",
                              x_axis_title = "Habitat amount in 2000 m radius",
                              xlims =NULL,
                              elev_thresh = 600)


edge100Data <-  
  results_all_edges %>%  filter(decrease_factor == 0) %>%  
  filter(distance_to_coastline < 60000) %>%  
  filter(cancov_2020 >-1) %>%  
  plot_ownership_distribution(x_var = "edgeRook_100_40",
                              x_axis_title = "Edge amount in 100 m radius",
                              xlims =NULL,
                              elev_thresh = 600)

edge2000Data <-  
  results_all_edges %>%  filter(decrease_factor == 0) %>%  
  filter(distance_to_coastline < 60000) %>%  
  filter(cancov_2020 >-1) %>%  
  plot_ownership_distribution(x_var = "edgeRook_2000_40",
                              x_axis_title = "Edge amount in 2000 m radius",
                              xlims =NULL,
                              elev_thresh = 600)


#EXPORT PLOTS 

#Reduce all edges 
ggsave(
  filename = "Figures/reduceALLedges_60000dist_600elev.pdf",               # File path and name
  plot = plot_all_edges,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "pdf",                       # Output format
  bg = "white"                          # Set background to white
)

#reduce just landscape edges 
ggsave(
  filename = "Figures/reduceLANDSCAPEedges_60000dist_600elev.pdf",               # File path and name
  plot = plot_landscape_edges,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "pdf",                       # Output format
  bg = "white"                          # Set background to white
)

#Export covariate structure by actor #####
ggsave(
  filename = "Figures/byOwnership_edge100.png",               # File path and name
  plot = edge100Data,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)


ggsave(
  filename = "Figures/byOwnership_edge2000.png",               # File path and name
  plot = edge2000Data,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)


#at 100m 
ggsave(
  filename = "Figures/byOwnership_hab100.png",               # File path and name
  plot = hab100Data,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

ggsave(
  filename = "Figures/byOwnership_hab2000.png",               # File path and name
  plot = hab2000Data ,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

# #get quantiles of elevation  
# Allelev = as.numeric(quantile(elevational_range_df$dem30m, 1, na.rm = TRUE))
# q95elev = as.numeric(quantile(elevational_range_df$dem30m, 0.95, na.rm = TRUE))
# q90elev = as.numeric(quantile(elevational_range_df$dem30m, 0.90, na.rm = TRUE))
# q80elev = as.numeric(quantile(elevational_range_df$dem30m, 0.80, na.rm = TRUE))


