#What if we reduce edge amount? 
#And what if we use a model without edge at all? 


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

source("scripts/02_OrganiseMurreletData.R")

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

long_lat <- final2020 %>% select(point_id, x, y)
#add in information on elevation at each point 
pt_elevation <- read.csv("Outputs/elevation_at_each_point.csv") 

#read in means and SDs of covariates from original model to enable coorrect scalings
meansAndSds

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

saveRDS



#Plot data 


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
        factor_values == 0.25 ~ 25,
        factor_values == 0.5 ~ 50,
        factor_values == 0.75 ~ 75,
        factor_values == 1 ~ 100,
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
    ylim(0,0.3)+
    labs(
      x = "% Reduction in Edge Density",  # Label for x-axis
      y = "Occupancy"   # Label for y-axis
    ) +
    facet_wrap(~ ownership, ncol = 3) +  # Customize facet labels
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",  # Place the legend at the bottom
      legend.title = element_blank()  # Remove the legend title
    ) +
    guides(
      fill = guide_legend(title = NULL),  # Remove the fill legend title
      color = guide_legend(title = NULL)  # Remove the color legend title
    )
}


result_df %>%  filter(cancov_2020 > -1) %>% #remove non-forest habitat points
    filter( distance_to_coastline < 50000) %>% #remove points v far from coast 
    filter(dem30m < 900) %>%  
  ownership_only_variation_plot

model





