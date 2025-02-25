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
#Viualise currently across PNW what the relationship at 2km2 is between habitat amount and edge density
hab_points  %>%
  ggplot( aes(x = habAmountDich_2000, y = edgeRook_2000_40)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", span = 0.75, se = FALSE, color = "blue") +
  labs(title = "Relationship Between Total Habitat and Edge Density",
       y = "Total Habitat (2km buffer)",
       x = "Proportion of Habitat as Edge (2km buffer") +
  theme_minimal()

#=================================================
#predict scenarios of future fragmentation amount
#=================================================
#add different amounts of percentage increase in fragmentation from current (with no change in habitat loss)

percent_change <- data.frame(percent_change = c(0, 0.25, 0.5,0.75,1))
hab_points_fragmentation <- hab_points %>% cross_join(percent_change) %>%  
  mutate(edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40*percent_change))
         
#asssume fragmentation increases linearly with decline habitat loss 
hab_points_fragmentation_linear_hab_loss <- hab_points %>% cross_join(percent_change) %>%  
  mutate(edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40*percent_change), 
# #hab amount declines linearly with increasing fragmentation
habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*percent_change))

#asssume fragmentation increases quadratically with decline habitat loss (see example figure below)
hab_points_fragmentation_quadratic_hab_loss <- hab_points %>% cross_join(percent_change) %>%  
  mutate(edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40*percent_change^2), 
         # #hab amount declines linearly 
         habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*percent_change))

#assume fragmentation increases in roughly n-shape with declining habitat area 
#Example of what a quadratic change in edge area looks like with declining habitat 
hab_points_fragmentation_nShape_hab_loss <- hab_points %>% cross_join(percent_change) %>%
  mutate(
    edgeRook_2000_40 = edgeRook_2000_40 + (edgeRook_2000_40 * (percent_change - percent_change^2)),
    # #hab amount declines linearly 
    habAmountDich_2000 = habAmountDich_2000 - (habAmountDich_2000*percent_change)
    )

#COME BACK TO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#example run through - is this working how I expect...Or would a decline in habit from 100% to 75% currently 
# make the same edge area as a decline from 50% to 25% (ie not considering starting context appropriately)
x = 1 # select a given row of data
print(percent_change)
y =2 #select a percentage change row
hab_points$edgeRook_2000_40[x]#edge amount
hab_points$habAmountDich_2000[x]#hab amount

edgeRook_2000_40_test = hab_points$edgeRook_2000_40[x] + (hab_points$edgeRook_2000_40[x] * (percent_change$percent_change[y] - percent_change$percent_change[y]^2))
print(edgeRook_2000_40_test)
# #hab amount declines linearly 
habAmountDich_2000_test = hab_points$habAmountDich_2000[x] - (hab_points$habAmountDich_2000[x]*percent_change$percent_change[y])
print(habAmountDich_2000_test)

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
  mutate(coast_categorical = factor(coast_categorical, levels = c("coast (<10km)", "10-23 km", ">23 km"))) 

#BOXPLOTS
frag_df %>%
  filter(ownership %in% c("Federal", "Private Industrial", "Private Non-industrial")) %>% 
  group_by(percent_change) %>% 
  ggplot(aes(x = as.factor(percent_change), y = Occupancy, fill = as.factor(OceanYear))) + 
  geom_boxplot(color = "black") +  # Boxplot with black outline
  geom_point(aes(color = OceanYear), size = 2, alpha = 0.05, position = position_jitterdodge(jitter.width = 0.1)) +  # Points with jitter
  # Custom fill colors for OceanYear
  scale_fill_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + 
  scale_color_manual(values = c("Good Ocean Years" = "#56B4E9", "Bad Ocean Years" = "#D55E00")) + # Color points
  theme_classic(base_size = 14) +  # Nature-style theme
  labs(
    x = "Proportion increase in future fragmentation per 2km",  # Label for x-axis
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
frag_quadratic_loss <-  prep_and_predict(hab_points_fragmentation_quadratic_hab_loss) #2km hab loss and quadratic incease in edges
frag_nShaped_loss <-  prep_and_predict(hab_points_fragmentation_nShape_hab_loss) #2km hab loss and nshaped incease in edges

#build plots 
frag_only_plot <- interacting_plots(frag_only)
frag_linear_loss_plot <-interacting_plots(frag_linear_loss)
frag_quadratic_loss_plot <- interacting_plots(frag_quadratic_loss)
frag_nShaped_loss_plot <- interacting_plots(frag_nShaped_loss)


#save figures

# Save the plot using ggsave
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

# 
#


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
