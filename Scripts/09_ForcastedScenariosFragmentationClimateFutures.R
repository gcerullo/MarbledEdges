#calculate fututure murrelet occuupancy under different ocean and fragmentation futures 

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

#-----------------------------------------
#define params ####
minimum_habitat_amount <- 0.05 # to be included, needs to have more than 
#this amount of habitat (removes site with loads of non-habitat)
#-----------------------------------------

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

#quickly check data locations
checkPoints <- vect(final2020, geom = c("x", "y"), crs = "EPSG:5070")  # Specify a CRS, e.g., WGS84
plot(SDM2020, main = "Raster with 500m Grid Points")
plot(checkPoints, add = TRUE, col = "blue", pch = 16, cex = 0.5)

#read in means and SDs of covariates from original model to enable coorrect scalings
meansAndSds


#_______________________________________________________________________________________________
#remind ourselves of range of covariates values that the original occupancy the model has seen: 
modelled_covs <- analysisSites %>% dplyr::select(scaleCoastDist,scaleEdgeDens100,scaleEdgeDens2000, scaleHabAmount100, scaleHabAmount2000, 
                                          PC1_t1)
happy_estimate_range <- modelled_covs%>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>% 
  ggplot( aes(x = Value)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "white", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms of Dataframe Columns", x = "Value", y = "Frequency")

#----------------------------------
#functions ####
#----------------------------------
#function for scaling data and making model predictions dataset ####
prep_and_predict <- function(x){
  
    # Prepare data for prediction with scaled covariates
  
 # scaling = original - mean / std 
  prediction_df <- x %>% mutate(
    scaleHabAmount100 = (habAmountDich_100  - meansAndSds$meanHabAmountDich100) /meansAndSds$sdHabAmountDich100,
    scaleHabAmount2000 = (habAmountDich_2000 - meansAndSds$meanHabAmountDich2000) /meansAndSds$sdHabAmountDich2000,
    scaleEdgeDens100 =  (edgeRook_100_40 - meansAndSds$meanEdgeDens100) /meansAndSds$sdEdgeRook100, 
    scaleEdgeDens2000 = (edgeRook_2000_40 - meansAndSds$meanEdgeDens2000)/meansAndSds$sdEdgeRook2000,
    scaleCoastDist = (distance_to_coastline - meansAndSds$meanCoastDist) /meansAndSds$sdCoastDist
    ) %>%  
    dplyr::select(point_id, scaleHabAmount100, scaleHabAmount2000, scaleEdgeDens100, scaleEdgeDens2000,scaleCoastDist,hab_loss_amount,cancov_2020) %>% 
    cross_join(covariates)
  
  # Predict occupancy with standard errors (for error ribbon)
  predictions <- predict(
    model,  
    newdata = prediction_df,
    type = "state", 
    se.fit = FALSE  # Obtain standard errors for predictions
  ) %>% 
    rename(Occupancy = Predicted)
  
  # Add predicted occupancy and standard errors to the data
  prediction_df <- prediction_df %>% cbind(predictions) %>%  
    mutate(lower_CI = Occupancy - 1.96 * SE, 
           upr_CI = Occupancy + 1.96 * SE )
  
  #add back in other info on ownership and SDM model 
  add_data_back <- x %>%  
    dplyr::select(point_id, point_leve_habitat, ownership,distance_to_coastline) %>% unique()
  
  prediction_df <- prediction_df %>%  left_join(add_data_back)
  
  return(prediction_df)
}

#to use if we use all points (and not just the poins that >.45 hab suitabili)
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
      hab_loss_amount,
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



#=================================================
#explore relationship between habitat amount and edge density

hab_points <- final2020  #include all points

 #check mean distance to coast 
 hab_points %>% summarise(mean_dist = mean(distance_to_coastline)) # mean suitable habitat is > 94km Km from coast
                                                                   # This is because murrelet SDMs dont account for coast dist
 
 max_edge <- max(final2020$edgeRook_2000_40)
# 
# hab_points_35 <- hab_points %>% filter(distance_to_coastline <= 35000)
# 
# #Viualise currently across PNW what the relationship at 2km2 is between habitat amount and edge density
# hab_points2020 <- hab_points  %>%
#   ggplot( aes(x = habAmountDich_2000, y = edgeRook_2000_40)) +
#   geom_point(size = 3, alpha = 0.1, shape = 1) +
#   geom_smooth(method = "loess", span = 0.75, se = FALSE, color = "red") +
#   labs(x = "Total Habitat (2km buffer)",
#        y = "Proportion of Habitat as Edge (2km buffer)") +
#   theme_minimal(base_size = 16) +  
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),  
#     axis.text = element_text(color = "black"),
#     axis.title = element_text(face = "bold")
#   ) +
#   scale_x_continuous(expand = c(0, 0)) +  
#   scale_y_continuous(expand = c(0, 0)) +  
#   theme(
#     panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "ivory")  
#   )
# 
# 
# # Fit the loess model
# loess_model <- loess(edgeRook_2000_40 ~ habAmountDich_2000, data = hab_points, span = 0.75)
# 
# # Create a dataframe with smoothed values
# loess_edge_amount <- data.frame(
#   habAmountDich_2000 = hab_points$habAmountDich_2000,
#   loess_edge_amount = predict(loess_model)
# )

#=================================================
#predict scenarios of future fragmentation amount
#=================================================
# a circle with a 2km2 radius is 12.57² = 1257 ha. 
buff2km_ha_area = 1257
#add different amounts of percentage increase in fragmentation from current (with no change in habitat loss)
hab_loss_amount <- data.frame(hab_loss_amount = c(0, 0.1,0.2,0.3,0.4,0.5)) %>% 
                                 mutate(hab_loss_amount_ha= hab_loss_amount * buff2km_ha_area)

#hab_loss_amount <- data.frame(hab_loss_amount = c(0, 30,60,90,120,150)) / 314 #GIVE IN hectares in instead

# hab_points_fragmentation <- hab_points %>% cross_join(hab_loss_amount) %>%  
#   mutate(edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40*hab_loss_amount))
         
# #asssume fragmentation increases linearly with decline habitat loss 
# hab_points_fragmentation_linear_hab_loss <- hab_points %>% cross_join(hab_loss_amount) %>%  
#   mutate(edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40*hab_loss_amount), 
# # #hab amount declines linearly with increasing fragmentation
# habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*hab_loss_amount))

#asssume fragmentation increases linearly with decline habitat loss
# - AND CONSIDER AMOUNT INSTEAD OF PERCENTAGE  
# Combine and calculate fragmentation and edge amount
hab_points_fragmentation_linear_hab_loss_cumalative <- hab_points %>%
  cross_join(hab_loss_amount) %>%
  mutate(
    
    #if habitat amount is < 0.0001 then make 0; functionally the same but low vals lead to odd proportional changes
    habAmountDich_2000 = if_else(habAmountDich_2000 < 0.0001, 0, habAmountDich_2000),  
    
    
    # Habitat amount declines linearly with increasing fragmentation
    habAmountDich_2000_temp = (habAmountDich_2000) - hab_loss_amount,  # Remove habitat in fixed amounts
    
    # Ensure total habitat loss is not negative
    habAmountDich_2000_temp = if_else(habAmountDich_2000_temp < 0, 0.0001, habAmountDich_2000_temp),  
    
    # Calculate proportional habitat change
    prop_hab_change = (habAmountDich_2000_temp - habAmountDich_2000) / habAmountDich_2000, #leads to Inf values if habAmountDich_2000 = 0 
     
    #nsure Inf values are dealt with correctly
    prop_hab_change = ifelse(is.infinite(prop_hab_change), 0, prop_hab_change),     # Apply inverse proportional change to edge area
    prop_edge_change = -prop_hab_change,
    
    # Update edge area based on habitat change
    edgeRook_2000_40_temp = edgeRook_2000_40 + (edgeRook_2000_40 * prop_edge_change),
    
    # Ensure edge area does not exceed the max value observed in the landscape
    edgeRook_2000_40 = if_else(edgeRook_2000_40 > max_edge, max_edge, edgeRook_2000_40_temp) , 
    
    #make sure hab amount is named correctly
    habAmountDich_2000 = habAmountDich_2000_temp
  )

# xtest <- hab_points_fragmentation_linear_hab_loss_cumalative %>% filter(point_id == "p642021")
#dfSummary(hab_points_fragmentation_linear_hab_loss_cumalative)
#overlyhigh <- hab_points_fragmentation_linear_hab_loss_cumalative %>%  filter(prop_hab_change >10)
# ##############
# ################
# 
# #Visualise different fragmentation and hab amount relationships
# 
# #non-linear quadratic effect
# hab_loss_amount_test <- seq(0, 1, length.out = 100)
# habitat_amount <- 1 - hab_loss_amount_test
# edge_area <- hab_loss_amount_test^2
# ggplot(data.frame(hab_loss_amount_test, habitat_amount, edge_area), aes(x = hab_loss_amount_test)) +
#   geom_line(aes(y = habitat_amount, color = "Habitat Amount"), size = 1) +
#   geom_line(aes(y = edge_area, color = "Edge Area"), size = 1, linetype = "dashed") +
#   labs(x = "Percent Change", y = "Amount", title = "Quadratic Relationship Between Habitat Amount and Edge Area") +
#   scale_color_manual(values = c("Habitat Amount" = "blue", "Edge Area" = "red")) +
#   theme_minimal() +
#   theme(legend.title = element_blank(), legend.position = "top")
# 
# # quadratic inverted U-shaped relationship
# habitat_amount <- 1 - hab_loss_amount_test  # Linear habitat loss
# edge_area <- hab_loss_amount_test - hab_loss_amount_test^2  # Inverted U-shape for edge area
# 
# # Create a dataframe
# data_plot <- data.frame(hab_loss_amount = hab_loss_amount_test, 
#                         habitat_amount = habitat_amount, 
#                         edge_area = edge_area)
# 
# # Plot the relationships
# ggplot(data_plot, aes(x = hab_loss_amount)) +
#   geom_line(aes(y = habitat_amount, color = "Habitat Amount"), size = 1) +
#   geom_line(aes(y = edge_area, color = "Edge Area"), size = 1, linetype = "dashed") +
#   labs(x = "Percent Change", y = "Amount", title = "Inverted U-Shaped Relationship Between Habitat Amount and Edge Area") +
#   scale_color_manual(values = c("Habitat Amount" = "blue", "Edge Area" = "red")) +
#   theme_minimal() +
#   theme(legend.title = element_blank(), legend.position = "top")
# 

#==========================================================
#plot figures for different assumptions of habitat loss, fragmentation and ocean, by ownership 
#------------------------------------------------------------

interacting_plots <- function(x){
#add categorical distance to coast info 
frag_df <- x %>%  
mutate(coast_categorical = case_when(
  
  distance_to_coastline < 20000 ~ "<20km",
  distance_to_coastline >= 20000 & distance_to_coastline < 40000 ~ "20-40 km",
  distance_to_coastline >= 40000 & distance_to_coastline < 60000 ~ "40-60 km",
  distance_to_coastline >= 60000 ~ ">60 km"
  
  # distance_to_coastline < 10000 ~ "coast (<10km)",
  # distance_to_coastline >= 10000 & distance_to_coastline < 23000 ~ "10-23 km",
  # distance_to_coastline >= 23000 ~ ">23 km"
))%>%
#  mutate(coast_categorical = factor(coast_categorical, levels = c("coast (<10km)", "10-23 km", ">23 km"))) %>% 
  mutate(coast_categorical = factor(coast_categorical, levels = c("<20km", "20-40 km", "40-60 km", ">60 km"))) %>% 
  
   #if plotting area mutate into hectares (a buffer is 314 hectate)
   mutate(hab_loss_amount = hab_loss_amount*314) 

#BOXPLOTS
frag_df %>%
  filter(ownership %in% c("Federal", "Private Industrial", "Private Non-industrial")) %>% 
  group_by(hab_loss_amount) %>% 
  ggplot(aes(x = as.factor(hab_loss_amount), y = Occupancy, fill = as.factor(OceanYear))) + 
  geom_point(aes(color = OceanYear), size = 2, alpha = 0.02, position = position_jitterdodge(jitter.width = .25)) +  # Points with jitter
  geom_boxplot(color = "black") +  # Boxplot with black outline
  # Custom fill colors for OceanYear
  scale_fill_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + 
  scale_color_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + # Color points
  theme_classic(base_size = 14) +  # Nature-style theme
  labs(
    x = "Projected increase in future habitat loss (ha ) per 2km buffer",  # Label for x-axis
    y = "Occupancy"   # Label for y-axis
  ) +
  facet_wrap(coast_categorical~ ownership, ncol = 3) +  # Customize facet labels
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",  # Place the legend at the bottom
    legend.title = element_blank()  # Remove the legend title
  ) +
  guides(
    fill = guide_legend(title = NULL),  # Remove the fill legend title
    color = guide_legend(title = NULL)  # Remove the color legend title
  )
}


ownership_only_variation_plot <- function(x){
  #add categorical distance to coast info 
  frag_df <- x %>%  
    #if plotting area mutate into hectares (a buffer is 314 hectate)
    mutate(hab_loss_amount_ha = round(buff2km_ha_area*hab_loss_amount))
  
  #BOXPLOTS
  frag_df %>%
    filter(ownership %in% c("Federal", "Private Industrial", "Private Non-industrial", "State")) %>% 
    mutate(ownership = factor(ownership, levels = c("State", "Federal", "Private Industrial", "Private Non-industrial"))) %>% 
    mutate(ownership = factor(ownership, levels = c("State", "Federal", "Private Non-industrial", "Private Industrial"))) %>%  
  
    group_by(hab_loss_amount) %>% 
    ggplot(aes(x = as.factor(hab_loss_amount_ha), y = Occupancy, fill = as.factor(OceanYear))) + 
    #geom_point(aes(color = OceanYear), size = 2, alpha = 0.02, position = position_jitterdodge(jitter.width = .25)) +  # Points with jitter
    geom_boxplot(color = "black", outlier.shape = NA) +  # Boxplot with black outline
    # Custom fill colors for OceanYear
    scale_fill_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + 
    scale_color_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + # Color points
    theme_classic(base_size = 14) +  # Nature-style theme
    ylim(0,0.5)+
    labs(
      x = "Projected increase in future habitat loss (ha) per 2km buffer",  # Label for x-axis
      y = "Occupancy"   # Label for y-axis
    ) +
    facet_wrap(~ ownership, ncol = 4) +  # Customize facet labels
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",  # Place the legend at the bottom
      legend.title = element_blank()  # Remove the legend title
    ) +
    guides(
      fill = guide_legend(title = NULL),  # Remove the fill legend title
      color = guide_legend(title = NULL)  # Remove the color legend title
    )
}

dfSummary(frag_linear_loss)
#prep-and-predict (use for smaller number of points - ie if filter hab >0.45)
frag_linear_loss <-  hab_points_fragmentation_linear_hab_loss_cumalative %>%  
 # filter( distance_to_coastline < 84000, cancov_2020 >0) %>% 
  prep_and_predict_parallelised()  #fragmentation and linear loss of 2km habitat 

saveRDS(frag_linear_loss, "hab_points_fragmentation_linear_hab_loss_cumalative_final_model_5thMay2025.rds")

head(frag_linear_loss)

#-----------------------------------------------------------------------------------
#CAN START HERE ####
frag_linear_loss <- readRDS("hab_points_fragmentation_linear_hab_loss_cumalative_final_model_5thMay2025.rds")
dfSummary(frag_linear_loss)

frag_linear_loss <- frag_linear_loss %>%  filter( distance_to_coastline < 100000) 

#check covariate structure
frag_linear_loss %>%  
  dplyr::select(scaleCoastDist,scaleEdgeDens100,scaleEdgeDens2000, scaleHabAmount100, scaleHabAmount2000, 
                PC1_t1) %>% 
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>% 
  ggplot( aes(x = Value)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "white", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms of Dataframe Columns", x = "Value", y = "Frequency")

hist(frag_linear_loss$scaleCoastDist)
highOcc <- frag_linear_loss %>%  
  filter(Occupancy >0.7) %>%  
  filter(hab_loss_amount == 0) %>%  
  filter(ownership == "State")

#compare to what we can confidently predict over in our model
happy_estimate_range 

#scale back to original hab amount 
frag_linear_loss <- frag_linear_loss %>% mutate(habAmountDich_2000 = (scaleHabAmount2000 * meansAndSds$sdHabAmountDich2000) + meansAndSds$meanHabAmountDich2000)

#quickly check the consequences of applying thresholds for removal of data with almost 
# no surrounding habitat on federal and state lands (as these are often snowy/ icy )
plot <- frag_linear_loss %>%
  filter(cancov_2020 > 0) %>% 

  #-------------------------------
  #THIS SECTION IS DEALING WITH POINTS THAT FALL IN FED.STATES LAND BUT WITH V LOW SURROUNDING CANOPY COVER 
  #filter out species points
  mutate(to_exclude = habAmountDich_2000 < minimum_habitat_amount & 
           hab_loss_amount == 0 & 
           ownership %in% c("State", "Federal")) %>%  
  # if any point_id has a to_exclude value of TRUE, all rows for that point_id also have to_exclude as TRUE
  group_by(point_id) %>% 
  mutate(to_exclude = any(to_exclude)) %>% 
  ungroup() %>%  
# Filter out excluded points
  filter(!to_exclude) %>%  
  #-------------------------------
  ownership_only_variation_plot()
plot
names(frag_linear_loss)



#----------------------------------------------------------------------
#plots of edge amount and habitat amount by landscape ownership ####
#-------------------------------------------------

  
plot_ownership_distribution <- function(data, x_var, x_axis_title, title) {
  
  # Filter and categorize data
  filtered_data <- data %>%
    filter(!ownership %in% c("Unknown", NA)) %>%  
    filter(ownership %in% c("Federal", "State", "Private Industrial", "Private Non-industrial")) %>% 
    filter(distance_to_coastline <84000) 

    # Conditionally filter based on x_var
    if (x_var == "habAmountDich_2000") {
      filtered_data <- filtered_data %>% filter(.data[[x_var]] > minimum_habitat_amount)
    }
  
    
    # mutate(coast_categorical = case_when(
    #   distance_to_coastline < 20000 ~ "<20km",
    #   distance_to_coastline >= 20000 & distance_to_coastline < 40000 ~ "20-40 km",
    #   distance_to_coastline >= 40000 & distance_to_coastline < 60000 ~ "40-60 km",
    #   distance_to_coastline >= 60000 ~ ">60 km"
    # )) %>%
    # mutate(coast_categorical = factor(coast_categorical, levels = c("<20km", "20-40 km","40-60 km", ">60 km")))
  
  # Compute median for each facet group
  summary_stats <- filtered_data %>%
    group_by(ownership) %>%
    summarise(
      median_value = median(.data[[x_var]], na.rm = TRUE),
      count = n(),
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
      title = title,
      x = x_axis_title,  
      y = "Frequency"
    ) +
    
    facet_wrap(~ ownership, scales = "free_y", ncol = 4, nrow = 1) +  
    theme_minimal(base_size = 16) +  
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),  
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    ) +
    scale_x_continuous(expand = c(0, 0)) +  
   # scale_y_continuous(expand = c(0, 0)) +  
  #  scale_y_continuous(expand = c(0, 0), limits = c(0, 600)) +  
    
    theme(
      panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "ivory")  
    )
}



edge_ownership <- plot_ownership_distribution(hab_points, 
                                              x_var = "edgeRook_2000_40",
                                              x_axis_title = "Proportion Edge Amount (2km)",
                                              title = "Distribution of Edge Amount in 2km Buffers by Ownership")

habitat_ownership <- plot_ownership_distribution(hab_points, 
                                                 x_var ="habAmountDich_2000", 
                                                 x_axis_title = "Proportion Habitat Amount (2km)",
                                                 title = "Distribution of Habitat Amount in 2km Buffers by Ownership")
edge_ownership
habitat_ownership

#test for 100 m buffers
edge_ownership100 <- plot_ownership_distribution(hab_points, "edgeRook_100_40",
                                              x_axis_title = "Proportion Edge Amount (100 m)",
                                              "Distribution of Edge Amount in 100 m  Buffers by Ownership")

habitat_ownership100 <- plot_ownership_distribution(hab_points, "habAmountDich_100", 
                                                 x_axis_title = "Proportion Habitat Amount (100 m)",
                                                 "Distribution of Habitat Amount in 100 m Buffers by Ownership")
edge_ownership100
habitat_ownership100

# #=================================================
# #What does ownership look like close to the coast? 
# #=================================================
final2020 %>% 
  filter(distance_to_coastline < 40000) %>%  
  group_by(ownership) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(x = ownership, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Number of Observations by Ownership",
       x = "Ownership Type",
       y = "Count") +
  theme_minimal()

#save figures

# Save the plot using ggsave
# ggsave(
#   filename = "Figures/points_habAmount_by_edge.png",               # File path and name
#   plot = hab_points2020,          
#   width = 10,                            # Width in inches (publication size)
#   height = 8,                           # Height in inches (publication size)
#   dpi = 300,                            # Resolution for publication (300 DPI)
#   units = "in",                         # Units for width and height
#   device = "png",                       # Output format
#   bg = "white"                          # Set background to white
# )


# ggsave(
#   filename = "Figures/forescasted_frag_only_plot.png",               # File path and name
#   plot = frag_only_plot,          
#   width = 10,                            # Width in inches (publication size)
#   height = 8,                           # Height in inches (publication size)
#   dpi = 300,                            # Resolution for publication (300 DPI)
#   units = "in",                         # Units for width and height
#   device = "png",                       # Output format
#   bg = "white"                          # Set background to white
# )

ggsave(
  filename = "Figures/forescasted_frag_with_linear_2km_hab_loss_plot.png",               # File path and name
  plot = frag_linear_loss_plot,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

#ownersip without distance stratification: 
ggsave(
  filename = "Figures/ownership_allpoints_forescasted_frag_with_linear_2km_hab_loss_plot.png",               # File path and name
  plot = owner_frag_linear_loss_plot,          
  width = 16,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                          # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

#ownersip without distance stratification: 
ggsave(
  filename = "Figures/ownership__less100km_forescasted_frag_with_linear_2km_hab_loss_plot.png",               # File path and name
  plot = owner_frag_linear_loss_plot_100km,          
  width = 16,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                          # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)
# 
# 
# ggsave(
#   filename = "Figures/forescasted_CUMALATIVE_frag_with_linear_2km_hab_loss_plot.png",               # File path and name
#   plot = frag_linear_loss_cumalative_plot,          
#   width = 10,                            # Width in inches (publication size)
#   height = 8,                           # Height in inches (publication size)
#   dpi = 300,                            # Resolution for publication (300 DPI)
#   units = "in",                         # Units for width and height
#   device = "png",                       # Output format
#   bg = "white"                          # Set background to white
# )
# 
# ggsave(
#   filename = "Figures/forescasted_CUMALATIVE2_frag_with_linear_2km_hab_loss_plot.png",               # File path and name
#   plot = frag_linear_loss_cumalative_plot2,          
#   width = 10,                            # Width in inches (publication size)
#   height = 8,                           # Height in inches (publication size)
#   dpi = 300,                            # Resolution for publication (300 DPI)
#   units = "in",                         # Units for width and height
#   device = "png",                       # Output format
#   bg = "white"                          # Set background to white
# )
# 
# ggsave(
#   filename = "Figures/forescasted_quadratic_frag_with__2km_hab_loss_plot.png",               # File path and name
#   plot = frag_quadratic_loss_plot,          
#   width = 10,                            # Width in inches (publication size)
#   height = 8,                           # Height in inches (publication size)
#   dpi = 300,                            # Resolution for publication (300 DPI)
#   units = "in",                         # Units for width and height
#   device = "png",                       # Output format
#   bg = "white"                          # Set background to white
# )
# 
# ggsave(
#   filename = "Figures/forescasted_nShaped_frag_with__2km_hab_loss_plot.png",               # File path and name
#   plot = frag_nShaped_loss_plot,          
#   width = 10,                            # Width in inches (publication size)
#   height = 8,                           # Height in inches (publication size)
#   dpi = 300,                            # Resolution for publication (300 DPI)
#   units = "in",                         # Units for width and height
#   device = "png",                       # Output format
#   bg = "white"                          # Set background to white
# )
# 
# 
# ggsave(
#   filename = "Figures/forescasted_nShapedLoess2_frag_with__2km_hab_loss_plot.png",               # File path and name
#   plot = frag_nShapedLoess_loss_plot,          
#   width = 10,                            # Width in inches (publication size)
#   height = 8,                           # Height in inches (publication size)
#   dpi = 300,                            # Resolution for publication (300 DPI)
#   units = "in",                         # Units for width and height
#   device = "png",                       # Output format
#   bg = "white"                          # Set background to white
# )

#Export fragmentation and habitat amount for points in hab_points (117,290 pts)
ggsave(
  filename = "Figures/edge_amount_by_actor.png",               # File path and name
  plot = edge_ownership,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

ggsave(
  filename = "Figures/hab_amount_by_actor.png",               # File path and name
  plot = habitat_ownership,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

ggsave(
  filename = "Figures/edge_01_increase.png",               # File path and name
  plot = frag_01edge_plot,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

# #Viualise currently across PNW what the relationship at 2km2 is between habitat amount and edge density
#  hab_points  %>%
#    ggplot( aes(x = habAmountDich_2000, y = edgeRook_2000_40)) +
#    geom_point(size = 3) +
#    geom_smooth(method = "loess", span = 0.75, se = FALSE, color = "blue") +
#    labs(title = "Relationship Between Total Habitat and Edge Density",
#         x = "Total Habitat (2km buffer)",
#         y = "Proportion of Habitat as Edge (2km buffer") +
#  theme_minimal()
# 
# # 
# gam_model <- gam(edgeRook_2000_40 ~ s(habAmountDich_2000), data = hab_points)
# summary(gam_model) #habitat amount explains about 11.9% of edge density; so there's a lot more that's important
# plot(gam_model, shade = TRUE, rug = TRUE, main = "GAM Fit: Habitat Amount vs. Edge Density")
# # 
# # # Create prediction data
# pred_data <- data.frame(habAmountDich_2000 = seq(min(hab_points$habAmountDich_2000),
#                                                max(hab_points$habAmountDich_2000),
#                                                length.out = 200))
# # Predict with confidence intervals
# pred_data$edgeRook_2000_40 <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
# 
# # Calculate the upper and lower bounds of the confidence interval
# pred_data$fit_upper <- pred_data$edgeRook_2000_40$fit + 1.96 * pred_data$edgeRook_2000_40$se.fit
# pred_data$fit_lower <- pred_data$edgeRook_2000_40$fit - 1.96 * pred_data$edgeRook_2000_40$se.fit
# 
# # Convert variables to numeric
# final2020$habAmountDich_2000 <- as.numeric(final2020$habAmountDich_2000)
# pred_data$edgeRook_2000_40 <- as.numeric(pred_data$edgeRook_2000_40$fit)
# pred_data$fit_upper <- as.numeric(pred_data$fit_upper)
# pred_data$fit_lower <- as.numeric(pred_data$fit_lower)
# 
# 
# # Scatterplot + GAM fitted curve with confidence intervals
# 
# # Scatterplot + GAM fitted curve with confidence intervals
# ggplot(hab_points, aes(x = habAmountDich_2000, y = edgeRook_2000_40)) +
#   geom_point(alpha = 0.2) +  # Raw data points
#   geom_line(data = pred_data, aes(x = habAmountDich_2000, y = edgeRook_2000_40),
#             color = "blue", size = 1.2) +  # GAM fit
#   geom_ribbon(data = pred_data, aes(x = habAmountDich_2000,
#                                     ymin = fit_lower, ymax = fit_upper),
#               fill = "blue", alpha = 0.2) +  # Confidence intervals
#   labs(title = "GAM Fit with Confidence Intervals: Habitat Amount vs. Edge Density",
#        x = "Habitat Amount (habAmountDich_2000)",
#        y = "Edge Density (edgeRook_2000_40)") +
#   theme_minimal()

# 
#  
#  # FROM VALENTE
# #understand rough relationship of Valente et al habitat amount and edge density.
#  # Create data for Total Habitat and Habitat Edge
# rough_relationship <- data.frame(
#    Year = c(1990, 1995, 2000, 2005, 2010, 2015, 2020),
#    Total_Habitat = c(2.40, 2.30, 2.20, 2.15, 2.00, 1.85, 1.90),  # Million hectares
#    Habitat_Edge = c(0.14, 0.147, 0.150, 0.152, 0.154, 0.16, 0.170)  # Proportion of total habitat
#  )
#  
# # Plot Total Habitat over Time
# ggplot(rough_relationship, aes(x = Year)) +
#   geom_line(aes(y = Total_Habitat), color = "blue") +
#   geom_point(aes(y = Total_Habitat), color = "blue") +
#   labs(
#     title = "Total Murrelet Habitat Over Time",
#     x = "Year",
#     y = "Total Habitat (Million Hectares)"
#   ) +
#   theme_minimal()
# 
# # Plot Habitat Edge over Time
# ggplot(rough_relationship, aes(x = Year)) +
#   geom_line(aes(y = Habitat_Edge), color = "red") +
#   geom_point(aes(y = Habitat_Edge), color = "red") +
#   labs(
#     title = "Habitat Edge Over Time",
#     x = "Year",
#     y = "Habitat Edge (Proportion of Total Habitat)"
#   ) +
#   theme_minimal() 
# 
# 
# 
# ggplot(rough_relationship, aes(x = Total_Habitat, y = Habitat_Edge)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "loess", span = 0.75, se = FALSE, color = "blue") +
#   labs(title = "Relationship Between Total Habitat and Edge Density",
#        x = "Total Habitat (Million ha)",
#        y = "Proportion of Habitat as Edge") +
#   theme_minimal()
# 


#---------------------------
#EXPLORATORY FIGURES:
#---------------------------
# #GEOM POINT
#  # Summarize data by OceanYear
#  summary_data2020 <- predict2020 %>%
#    filter(point_leve_habitat >= 45) %>%  
#    group_by(OceanYear) %>%
#    summarize(
#      mean_occupancy = mean(Occupancy, na.rm = TRUE),   # Mean Occupancy
#      mean_SE = sqrt(sum(SE^2, na.rm = TRUE) / n()),   # Pooled SE
#      lower_CI = mean_occupancy - qnorm(0.975) * mean_SE, # Lower 95% CI
#      upper_CI = mean_occupancy + qnorm(0.975) * mean_SE  # Upper 95% CI
#    )
#  
# 
#  # Plot summarized data
#  ggplot(summary_data2020, aes(x = OceanYear, y = mean_occupancy)) +
#    geom_point(size = 3, color = "blue") +  # Points for mean occupancy
#    geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, color = "darkgray") +  # Error bars for CI
#    labs(
#      title = "Mean Occupancy by Ocean Year with Confidence Intervals",
#      x = "Ocean Year",
#      y = "Mean Occupancy"
#    ) +
#    theme_minimal() +
#    theme(
#      text = element_text(size = 12),
#      plot.title = element_text(hjust = 0.5, size = 16)
#    )
#  

 # Plot Occupancy vs. Distance to Coastline
 # ggplot(predict2020, aes(x = scaleCoastDist, y = Occupancy)) +
 #   geom_point(size = 2, color = "blue", alpha = 0.6) +  # Scatter plot of points
 #   geom_smooth(method = "lm", color = "red", se = TRUE) +  # Linear trend line with confidence interval
 #   labs(
 #     title = "Occupancy vs. Distance to Coastline",
 #     x = "Scaled Distance to Coastline",
 #     y = "Occupancy"
 #   ) +
 #   theme_minimal() +
 #   theme(
 #     text = element_text(size = 12),
 #     plot.title = element_text(hjust = 0.5, size = 16)
 #   )
 #
