#Compute additional variables for marbled murrelet models

#read in packages 
library(tidyverse)

#read in inputs 

#ocean conditions 
#Raw data on the ocean conditions based on 16 PCAd variables; downloaded from: 
#https://www.fisheries.noaa.gov/west-coast/science-data/ocean-conditions-indicators-trends
ocean_cond <- read.csv("Inputs/2023-Stoplight-RAWDATA.csv") %>%  
  #remove additional data not used in PCA construction
  slice(1:21)


#--------------------------------------------------------------------
#OCEAN CODITIONS 
#--------------------------------------------------------------------




#forest_age information 

#Clean data #### 

ocean_cond <- ocean_cond %>%  
  # Remove rows where all columns are NA
  filter(if_all(everything(), ~ !is.na(.))) %>%
  # Convert all columns starting with "X" to character
  mutate(across(starts_with("X"), as.character)) %>% 
  # Pivot to long format
  pivot_longer(
    cols = starts_with("X"),  # Select columns that start with 'X' (year columns)
    names_to = "Year",        # New column for year
    values_to = "Value"       # New column for values
  ) %>%
  # Clean up Year column to remove 'X' prefix
  mutate(Year = str_remove(Year, "^X"))


#filter PC1 data and scale and invert values (so that high values indicate good ocean years) 
PC1 <- ocean_cond %>% filter( Ecosystem.Indicators == "Principal Component scores (PC1)") %>%  
  mutate(Values = as.numeric(Value)) %>%  
  #invert values
  mutate(Value_inverted = max(Values) - Values) %>%  # Subtract each value from the maximum so that big values are good years, low values are bad years
  #scale values
  mutate(Value_scaled = scale(Value_inverted))     # Standardize values; Centers the data to a mean of 0 and scales it to a standard deviation of 1.
  

PC1_quantiles <-  PC1 %>%
  summarise(
    q10 = quantile(`Value_scaled`, 0.1),
    q90 = quantile(`Value_scaled`, 0.9)
  )

#visualise data #### 
PC1 %>%  
  ggplot(aes(x = Year, y = Value_scaled, group = 1)) +  # Set group = 1 for all data points
  geom_line(color = "#2c7fb8", linewidth = 1) +  # Line plot with custom color
  geom_point(color = "#d95f0e", size = 3) +      # Points with custom color
  theme_minimal(base_size = 14) +               # Clean minimal theme
  labs(
    title = "Scaled and inverted PCA1 Scores Over Time (high = good year)",
    x = "Year",
    y = "Value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

#Oregon nest and murrelet data


##########################################################################
#save outputs #### 

write.csv(PC1, "Outputs/PC1_scaled_inverted.csv")
