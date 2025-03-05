#calculate fututure murrelet occuupancy under different ocean and fragmentation futures 

library(terra)
library(sf)
library(exactextractr)
library(data.table)
library(tidyverse)
library(mgcv)

# read in inputs #####
covariates <- readRDS("Outputs/ScaledCovariates.rds") %>%   
  dplyr::select(PC1_t1, scaleDoy, scaleDoy2, scaleDoy2, OceanYear) %>%  unique() %>%  
  #Lets assume the worst PC1 (2016) conditions on record; #-1.880589
  mutate(PC1_t1 = case_when(
    OceanYear == "Bad Ocean Years" ~ -1.880589,
    TRUE ~ PC1_t1  # keeps the original value of PC1_t1 when OceanYear is not "Bad Ocean Years"
  ))

model <- readRDS("Models/pc1_interaction_model.rds")
model <- readRDS("Models/pc1_3wayinteraction_model.rds")
final2020 <- readRDS("PNW_2020_extracted_covars.rds") %>%   #read in starting occupancy for 2020 from scrippt 8
as.data.frame() %>%  
  #only keep points with an SDM value
filter(!is.na(point_leve_hab_probability))

#----------------------------------
#functions ####
#----------------------------------

#function for scaling data and making model predictions dataset ####
prep_and_predict <- function(x){
  
  # Prepare data for prediction with scaled covariates
  prediction_df <- x %>% mutate(
    scaleHabAmount100 = scale(habAmountDich_100), 
    scaleHabAmount2000 = scale(habAmountDich_2000), 
    scaleEdgeDens100 = scale(edgeRook_100_40), 
    scaleEdgeDens2000 = scale(edgeRook_2000_40), 
    scaleCoastDist = scale(distance_to_coastline)) %>%  
    dplyr::select(point_id, scaleHabAmount100, scaleHabAmount2000, scaleEdgeDens100, scaleEdgeDens2000,scaleCoastDist,percent_change) %>% 
    cross_join(covariates)
  
  # Predict occupancy with standard errors (for error ribbon)
  predictions <- predict(
    model,  
    newdata = prediction_df,
    type = "state", 
    se.fit = TRUE  # Obtain standard errors for predictions
  ) %>% 
    rename(Occupancy = Predicted)
  
  # Add predicted occupancy and standard errors to the data
  prediction_df <- prediction_df %>% cbind(predictions) %>%  
    mutate(lower_CI = Occupancy - 1.96 * SE, 
           upr_CI = Occupancy + 1.96 * SE )
  
  #add back in other info on ownership and SDM model 
  add_data_back <- x %>%  
    select(point_id, point_leve_hab_probability, ownership,distance_to_coastline) %>% unique()
  
  prediction_df <- prediction_df %>%  left_join(add_data_back)
  
  return(prediction_df)
}

#=================================================
#explore relationship between habitat amount and edge density
#=================================================
hab_points <- final2020 %>% filter(point_leve_hab_probability >=45) 
#check mean distance to coast 
hab_points %>% summarise(mean_dist = mean(distance_to_coastline)) # mean suitable habitat is > 48.6Km from coast
                                                                  # This is because murrelet SDMs dont account for coast dist

max_edge <- max(final2020$edgeRook_2000_40)

hab_points_35 <- hab_points %>% filter(distance_to_coastline <= 35000)

#Viualise currently across PNW what the relationship at 2km2 is between habitat amount and edge density
hab_points2020 <- hab_points  %>%
  ggplot( aes(x = habAmountDich_2000, y = edgeRook_2000_40)) +
  geom_point(size = 3, alpha = 0.1, shape = 1) +
  geom_smooth(method = "loess", span = 0.75, se = FALSE, color = "red") +
  labs(x = "Total Habitat (2km buffer)",
       y = "Proportion of Habitat as Edge (2km buffer)") +
  theme_minimal(base_size = 16) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) +  
  theme(
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "ivory")  
  )


# Fit the loess model
loess_model <- loess(edgeRook_2000_40 ~ habAmountDich_2000, data = hab_points, span = 0.75)

# Create a dataframe with smoothed values
loess_edge_amount <- data.frame(
  habAmountDich_2000 = hab_points$habAmountDich_2000,
  loess_edge_amount = predict(loess_model)
)

#=================================================
#predict scenarios of future fragmentation amount
#=================================================
# a circle with a 2km2 radius is 3.14km² = 314 ha. So 0.1 = 31 ha, 62 ha 
#add different amounts of percentage increase in fragmentation from current (with no change in habitat loss)
percent_change <- data.frame(percent_change = c(0, 0.1,0.2,0.3,0.4,0.5))

percent_change <- data.frame(percent_change = c(0, 30,60,90,120,150)) / 314 #GIVE IN hectares in instead
hab_points_fragmentation <- hab_points %>% cross_join(percent_change) %>%  
  mutate(edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40*percent_change))
         
# #asssume fragmentation increases linearly with decline habitat loss 
# hab_points_fragmentation_linear_hab_loss <- hab_points %>% cross_join(percent_change) %>%  
#   mutate(edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40*percent_change), 
# # #hab amount declines linearly with increasing fragmentation
# habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*percent_change))

#asssume fragmentation increases linearly with decline habitat loss
# - AND CONSIDER AMOUNT INSTEAD OF PERCENTAGE  
# Combine and calculate fragmentation and edge amount
hab_points_fragmentation_linear_hab_loss_cumalative <- hab_points %>%
  cross_join(percent_change) %>%
  mutate(
    # Habitat amount declines linearly with increasing fragmentation
    habAmountDich_2000_temp = habAmountDich_2000 - percent_change,  # Remove habitat in fixed amounts
    
    # Ensure total habitat loss is not negative
    habAmountDich_2000_temp = if_else(habAmountDich_2000_temp < 0, 0, habAmountDich_2000_temp),
    
    # Calculate proportional habitat change
    prop_hab_change = (habAmountDich_2000_temp - habAmountDich_2000) / habAmountDich_2000,
    
    # Apply inverse proportional change to edge area
    prop_edge_change = -prop_hab_change,
    
    # Update edge area based on habitat change
    edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40 * prop_edge_change),
    
    # Ensure edge area does not exceed 0.15
    edgeRook_2000_40 = if_else(edgeRook_2000_40 > max_edge, max_edge, edgeRook_2000_40), 
    
    #make sure hab amount is named correctly
    habAmountDich_2000 = habAmountDich_2000_temp
  )

hab_points35_fragmentation_linear_hab_loss_cumalative <- hab_points_fragmentation_linear_hab_loss_cumalative %>%  
  filter(distance_to_coastline <= 35000)
# 
# #asssume fragmentation increases in steps of 0.01 and that habitat reduces linearly with increasing edge habitat loss
# # - AND CONSIDER AMOUNT INSTEAD OF PERCENTAGE  
# edge_increase <- data.frame(edge_increase = c(0.01, 0.02,0.03,0.04,0.05,0.06))
# 
# hab_points_fragmentation_increase001 <- hab_points %>% cross_join(edge_increase) %>%  
#   mutate(
#     #ie if you lose 25% of habitat, you gain 25% of edge
#     edgeRook_2000_40 = edgeRook_2000_40 + edge_increase, 
#     # hab amount declines linearly with increasing fragmentation
#     habAmountDich_2000 = habAmountDich_2000 - (edge_increase *habAmountDich_2000)
#   ) %>% mutate(percent_change = edge_increase)
# 
# hab_points45_fragmentation_increase001 <- hab_points_45 %>% cross_join(edge_increase) %>%  
#   mutate(
#     #ie if you lose 25% of habitat, you gain 25% of edge
#     edgeRook_2000_40 = edgeRook_2000_40 + edge_increase, 
#     # hab amount declines linearly with increasing fragmentation
#     habAmountDich_2000 = habAmountDich_2000 - (edge_increase *habAmountDich_2000)
#   ) %>% mutate(percent_change = edge_increase)
# 
# 
# #above if you have less habitat than can be lost, no deforestation happens. 
# #here, if habAmountDich_2000 <  than percentage loss, you lose all remaining habitat, and edge and hab amount become 0
# hab_points2_fragmentation_linear_hab_loss_cumalative <- hab_points %>% 
#   cross_join(percent_change) %>%  
#   mutate(
#     # Only mutate if habAmountDich_2000 >= percent_change
#     edgeRook_2000_40 = if_else(
#       habAmountDich_2000 >= percent_change, 
#       edgeRook_2000_40 + (edgeRook_2000_40 * (percent_change / habAmountDich_2000)), 
#       0  # Set to 0 if condition is not met
#     ),
#     
#     habAmountDich_2000 = if_else(
#       habAmountDich_2000 >= percent_change, 
#       habAmountDich_2000 - percent_change, 
#       0  # Set to 0 if condition is not met
#     )
#   )
# 
# 
# edge = 0.0249112 
# hab = 0.400001
# edge + (edge*(0.4/hab))
# 0.5/hab
# #asssume fragmentation increases quadratically with decline habitat loss (see example figure below)
# hab_points_fragmentation_quadratic_hab_loss <- hab_points %>% cross_join(percent_change) %>%  
#   mutate(edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40*percent_change^2), 
#          # #hab amount declines linearly 
#          habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*percent_change))
# 
# #assume fragmentation increases in roughly n-shape with declining habitat area 
# #Example of what a quadratic change in edge area looks like with declining habitat 
# hab_points_fragmentation_nShape_hab_loss <- hab_points %>% cross_join(percent_change) %>%
#   mutate(
#     edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40 * (percent_change - percent_change^2)),
#     # #hab amount declines linearly 
#     habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*percent_change)
#     )
# hist(hab_points_fragmentation_nShape_loess_hab_loss$edgeRook_2000_40, breaks = 100)
# 
# #assume fragmentation increases in loess-derived n-shape with declining habitat area
# #here instead we use the loess n-shape of hab amount to edge amount from actually murrelet habitat to modify 
# #edge amount directly. We: 
# #1 Moidfy habitat amount using percentage 
# #2. For the specific habitat amount left after, we calculate edge amount based on the loess curve shape 
# #. We do this by joining loess_derived edge amount to habitat_amount (based ona rolling join, where we join to the closest habitat amount we have a loess value for)
# hab_points_percent <- hab_points %>% cross_join(percent_change) %>%  
#   mutate(habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*percent_change))
# hist(hab_points_percent$habAmountDich_2000, breaks = 100)  
# # Convert data frames to data.table and perform rolling join
# hab_points_fragmentation_nShape_loess_hab_loss <- setDT(loess_edge_amount)[setDT(hab_points_percent),
#                                                                            on = .(habAmountDich_2000), roll = "nearest", mult = "first"] %>%  
#   as.data.frame() %>%  
#   select(-edgeRook_2000_40) %>% # remove actual edgeRook_2000_40 - we replace with loess derived estimates 
#   rename( edgeRook_2000_40 = loess_edge_amount)
# 
# 
# ##############
################
#COME BACK TO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Logic test

#WRONG - if you lose 100% of your habitat (percentage change =1), the edge amount stays the same! 

#example run through - is this working how I expect...Or would a decline in habit from 100% to 75% currently 
# make the same edge area as a decline from 50% to 25% (ie not considering starting context appropriately)
# Create the final dataframe directly
test_df <- expand_grid(
  edgeRook_2000_40 = seq(0, 1, by = 0.01),
  habAmountDich_2000 = seq(0, 1, by = 0.1),
  percent_change = seq(0, 1, by = 0.2)
)

test_df2 <- test_df %>%  
  mutate( # #hab amount declines linearly 
    v2_habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*percent_change),
    v2_edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40 * (percent_change - percent_change^2))
  )

##############
################

#Visualise different fragmentation and hab amount relationships

#non-linear quadratic effect
percent_change_test <- seq(0, 1, length.out = 100)
habitat_amount <- 1 - percent_change_test
edge_area <- percent_change_test^2
ggplot(data.frame(percent_change_test, habitat_amount, edge_area), aes(x = percent_change_test)) +
  geom_line(aes(y = habitat_amount, color = "Habitat Amount"), size = 1) +
  geom_line(aes(y = edge_area, color = "Edge Area"), size = 1, linetype = "dashed") +
  labs(x = "Percent Change", y = "Amount", title = "Quadratic Relationship Between Habitat Amount and Edge Area") +
  scale_color_manual(values = c("Habitat Amount" = "blue", "Edge Area" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "top")

# quadratic inverted U-shaped relationship
habitat_amount <- 1 - percent_change_test  # Linear habitat loss
edge_area <- percent_change_test - percent_change_test^2  # Inverted U-shape for edge area

# Create a dataframe
data_plot <- data.frame(percent_change = percent_change_test, 
                        habitat_amount = habitat_amount, 
                        edge_area = edge_area)

# Plot the relationships
ggplot(data_plot, aes(x = percent_change)) +
  geom_line(aes(y = habitat_amount, color = "Habitat Amount"), size = 1) +
  geom_line(aes(y = edge_area, color = "Edge Area"), size = 1, linetype = "dashed") +
  labs(x = "Percent Change", y = "Amount", title = "Inverted U-Shaped Relationship Between Habitat Amount and Edge Area") +
  scale_color_manual(values = c("Habitat Amount" = "blue", "Edge Area" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "top")


#==========================================================
#plot figures for different assumptions of habitat loss, fragmentation and ocean, by ownership 

interacting_plots <- function(x){
#add categorical distance to coast info 
frag_df <- x %>%  
mutate(coast_categorical = case_when(
  distance_to_coastline < 10000 ~ "coast (<10km)",
  distance_to_coastline >= 10000 & distance_to_coastline < 23000 ~ "10-23 km",
  distance_to_coastline >= 23000 ~ ">23 km"
))%>%
  mutate(coast_categorical = factor(coast_categorical, levels = c("coast (<10km)", "10-23 km", ">23 km"))) %>% 
   #if plotting area mutate into hectares (a buffer is 314 hectate)
   mutate(percent_change = percent_change*314) 

#BOXPLOTS
frag_df %>%
  filter(ownership %in% c("Federal", "Private Industrial", "Private Non-industrial")) %>% 
  group_by(percent_change) %>% 
  ggplot(aes(x = as.factor(percent_change), y = Occupancy, fill = as.factor(OceanYear))) + 
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
  facet_wrap(coast_categorical~ ownership) +  # Customize facet labels
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

#run different scenarios
frag_only <- prep_and_predict(hab_points_fragmentation) #just fragmentation
frag_linear_loss <-  prep_and_predict(hab_points_fragmentation_linear_hab_loss) #fragmentation and linear loss of 2km habitat 
frag_linear_loss_cumalative <-  prep_and_predict(hab_points_fragmentation_linear_hab_loss_cumalative) #hab loss refers to actual AMOUNT of habitat loss. so 0.1 of 2k (or200 ha) = 20hectare, 0.5 = 100ha 
frag_linear_loss_cumalative2 <-  prep_and_predict(hab_points35_fragmentation_linear_hab_loss_cumalative) #hab and edge amount can go down to 0; hab loss refers to actual AMOUNT of habitat loss. so 0.1 of 2k (or200 ha) = 20hectare, 0.5 = 100ha 
frag_quadratic_loss <-  prep_and_predict(hab_points_fragmentation_quadratic_hab_loss) #2km hab loss and quadratic incease in edges
frag_nShaped_loss <-  prep_and_predict(hab_points_fragmentation_nShape_hab_loss) #2km hab loss and nshaped incease in edges
frag_nShapedLoess_loss <-  prep_and_predict(hab_points_fragmentation_nShape_loess_hab_loss) #2km hab loss and nshaped incease in edges
frag_01edge <- prep_and_predict(hab_points_fragmentation_increase001)
frag45_01edge <- prep_and_predict(hab_points45_fragmentation_increase001)

#build plots 
frag_only_plot <- interacting_plots(frag_only)
frag_linear_loss_plot <-interacting_plots(frag_linear_loss)
frag_linear_loss_cumalative_plot <- interacting_plots(frag_linear_loss_cumalative)
frag_linear_loss_cumalative_plot2 <-  interacting_plots(frag_linear_loss_cumalative2) #hab loss refers to actual AMOUNT of habitat loss. so 0.1 of 2k (or200 ha) = 20hectare, 0.5 = 100ha 
frag_quadratic_loss_plot <- interacting_plots(frag_quadratic_loss)
frag_nShaped_loss_plot <- interacting_plots(frag_nShaped_loss)
frag_nShapedLoess_loss_plot <- interacting_plots(frag_nShapedLoess_loss)
frag_01edge_plot <- interacting_plots(frag_01edge)
frag45_01edge_plot <- interacting_plots(frag45_01edge)


#----------------------------------------------------------------------
#plots of edge amount and habitat amount by landscape ownership 
model <- readRDS("Models/pc1_interaction_model.rds")


plot_ownership_distribution <- function(data, x_var, x_axis_title, title) {
  
  # Filter and categorize data
  filtered_data <- data %>%
    filter(!ownership %in% c("Unknown", NA)) %>%  
    filter(ownership %in% c("Federal", "State", "Private Industrial", "Private Non-industrial")) %>% 
    mutate(coast_categorical = case_when(
      distance_to_coastline < 10000 ~ "<10km",
      distance_to_coastline >= 10000 & distance_to_coastline < 23000 ~ "10-23 km",
      distance_to_coastline >= 23000 ~ ">23 km"
    )) %>%
    mutate(coast_categorical = factor(coast_categorical, levels = c("<10km", "10-23 km", ">23 km")))
  
  # Compute median for each facet group
  medians <- filtered_data %>%
    group_by(coast_categorical, ownership) %>%
    summarise(median_value = median(.data[[x_var]], na.rm = TRUE), .groups = "drop")
  
  # Create histogram plot with facet and median lines
  ggplot(filtered_data, aes(x = .data[[x_var]])) +  
    geom_histogram(fill = "#4C9A2A", color = "black", alpha = 0.7, bins = 100) +  
    geom_vline(data = medians, aes(xintercept = median_value), 
               linetype = "dashed", color = "red", linewidth = 1) +  
    labs(
      title = title,
      x = x_axis_title,  
      y = "Frequency"
    ) +
    facet_wrap(~interaction(coast_categorical, ownership), scales = "free_y", ncol = 3, nrow = 4) +  
    theme_minimal(base_size = 16) +  
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),  
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    ) +
    scale_x_continuous(expand = c(0, 0)) +  
    scale_y_continuous(expand = c(0, 0)) +  
    theme(
      panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "ivory")  
    )
}

#make ownership plots 
edge_ownership <- plot_ownership_distribution(hab_points, "edgeRook_2000_40",
                                              x_axis_title = "Proportion Edge Amount (2km)",
                                              "Distribution of Edge Amount in 2km Buffers by Ownership")

habitat_ownership <- plot_ownership_distribution(hab_points, "habAmountDich_2000", 
                                                 x_axis_title = "Proportion Habitat Amount (2km)",
                                                 "Distribution of Habitat Amount in 2km Buffers by Ownership")

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
ggsave(
  filename = "Figures/points_habAmount_by_edge.png",               # File path and name
  plot = hab_points2020,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)


ggsave(
  filename = "Figures/forescasted_frag_only_plot.png",               # File path and name
  plot = frag_only_plot,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

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

ggsave(
  filename = "Figures/forescasted_CUMALATIVE_frag_with_linear_2km_hab_loss_plot.png",               # File path and name
  plot = frag_linear_loss_cumalative_plot,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

ggsave(
  filename = "Figures/forescasted_CUMALATIVE2_frag_with_linear_2km_hab_loss_plot.png",               # File path and name
  plot = frag_linear_loss_cumalative_plot2,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

ggsave(
  filename = "Figures/forescasted_quadratic_frag_with__2km_hab_loss_plot.png",               # File path and name
  plot = frag_quadratic_loss_plot,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

ggsave(
  filename = "Figures/forescasted_nShaped_frag_with__2km_hab_loss_plot.png",               # File path and name
  plot = frag_nShaped_loss_plot,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)


ggsave(
  filename = "Figures/forescasted_nShapedLoess2_frag_with__2km_hab_loss_plot.png",               # File path and name
  plot = frag_nShapedLoess_loss_plot,          
  width = 10,                            # Width in inches (publication size)
  height = 8,                           # Height in inches (publication size)
  dpi = 300,                            # Resolution for publication (300 DPI)
  units = "in",                         # Units for width and height
  device = "png",                       # Output format
  bg = "white"                          # Set background to white
)

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
#    filter(point_leve_hab_probability >= 45) %>%  
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
