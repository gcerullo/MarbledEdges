# This code is for detailing key paramatres for building scenarios for different 
# starting landscapes for the Elliott 

library(terra)
library(data.table)
#library(rflsgen) #for generating landscapes with different fragmentation patterns
library(tidyverse)
library(landscapeR)

#read in custom functions
source("Inputs/Scenarios/Functions.R")
source("Inputs/Scenarios/CalculateYields.R")

#Read in inputs ####
# predicted yield curves for Elliot stands
#Elliott <- read.csv("Inputs/Scenarios/yields_exst.csv", na = c("", "NA")) %>% process_data()

options(scipen = 999)

#DEFINE PARAMS
#-------------------------------------------------------------------------------------
#Set landscape size ####
#-------------------------------------------------------------------------------------
landscape_size_ha <- 1000000 # 1 million hectares

#-------------------------------------------------------------------------------------
#Set scenario length ####
#-------------------------------------------------------------------------------------
scenarioLength <- 150 # years 


#-------------------------------------------------------------------------------------
#Define MAI of forest and yield error allowed in scenario construction  ####
#-------------------------------------------------------------------------------------
# #this is needed to calculate the even flow vol from starting landscapes
# #median of MAI for Elliot stands at 150 yrs
# MAI <- Elliott %>% filter(Age == 150) %>% 
#   mutate(MAI = m3ha/150) %>%  
#   summarise(medianMAI = median (MAI)) %>% pull()

#Define error about yields which is permitted 
error <- 0.0002 # +/- tiny error 


#--------------------------------------------------------------------------------
#Build starting landscapes ####
#-------------------------------------------------------------------------------------

  NoCl_old <- tibble(
    age = c(300),
    habitat = "forest",
    area = landscape_size_ha, 
    startingL = "NoCl_old")
  
  # Combine all tibbles into a single dataframe called startingL
  startingL <- bind_rows(
                         NoCl_old)


#-------------------------------------------------------------------------------------
#Visualize the composition of the starting landscape
#-------------------------------------------------------------------------------------

# Get the unique starting landscapes
landscapes <- unique(startingL$startingL)

# Set up plot layout (3 rows, 3 columns for 9 landscapes)
par(mfrow = c(3, 3))

# Loop through each landscape and plot
for (landscape in landscapes) {
  plot_landscape(startingL, landscape)
}

#-------------------------------------------------------------------------------------
#Calculate maximum production potential of each starting landscape (ie P = 1)
#-------------------------------------------------------------------------------------
#RETURN TO; THE MAX EVEN FLOW IS GOING TO VARY ACCORDING TO THE AGE OF THE FOREST
#AT SCENARIO START I THINK 

#NB to do this we calculate the max even flow that could be achieved for remaining forested areas in the landscape;
#i.e, the volume of wood that would be derivable each year over 150 yrs in a sustained way 
#max_SL_production_volume gives the 150 yr vol

# maxEvenFlow <-  startingL %>% 
#   filter(habitat == "forest") %>% 
#   mutate(
#          
#          #max even flow is given as the area of forest in the landscape * MAI
#          max_SL_production_volume = area * MAI *scenarioLength) %>% 
#   
#   dplyr::select(startingL, max_SL_production_volume)  

#Volume 

Elliot_summary %>% print(n=100)

evenFlow_params <- tibble(
  age = c(90,150, 300),
  growing_stock_ha = c(391.32059,549.19579,716.67762))


# max even flow if mangaged continuously with 60yr thin = MAI + growthStock/SimulationLength 
maxEvenFlow <-  startingL %>%
  filter(habitat == "forest") %>%
  left_join(evenFlow_params) %>% 
  mutate(
    forestGrowingStock = growing_stock_ha*area, 
    forestMAI = forestGrowingStock/scenarioLength, 
    forestEvenFlow_m3yr = forestMAI + (forestGrowingStock/scenarioLength),
    
    #max even flow (m3_forest_over150 years
    max_SL_production_volume = forestEvenFlow_m3yr *scenarioLength) %>%
  dplyr::select(startingL,age, max_SL_production_volume)  



# #-------------------------------------------------------------------
# #calculate range of production targets ####
# #-------------------------------------------------------------------
# Define the sequence of production targets
production_targets <- tibble(production_target = seq(0.01, 2, by = 0.01))
production_targets <- production_targets %>% merge(maxEvenFlow) %>% 
  mutate(production_target_m3 = production_target * max_SL_production_volume, 
         #add error permitted about yields 
         production_target_m3_max = production_target_m3+ (production_target_m3*error), 
         production_target_m3_min = production_target_m3- (production_target_m3*error)
  )


# #-------------------------------------------------------------------
#Generate many different triad scenarios via brute-force approach####
# #-------------------------------------------------------------------

# This approach generates many scenarios based on the proportional 
# coverage of reserve, intensive, and extensive treatments
options(scipen=999)

#Check that what you are doing is computationally feasible 

#how many scenarios would this generate: 
dec_places <- 0.1 # enter how many decimal places you want
num_comp <- 7 #enter how many unique compartments
(1/dec_places)^num_comp # this tells you how many different combinations there are


#Triad approach (3 compartments)
triad <- expand.grid(R = seq(0, 1, by = 0.005), 
                     E = seq(0, 1, by = 0.005), 
                     I = seq(0, 1, by = 0.005)) %>% 
  # Filter combinations where the row sum equals 1 
  filter(rowSums(across(everything())) == 1)

# Create a generalized function for adding 0.001 increments and 'filling' 
# between 0.005 steps 

add_increments <- function(df, inc_col1, dec_col1) {
  df %>%
    mutate(!!inc_col1 := .data[[inc_col1]] + 0.001, 
           !!dec_col1 := .data[[dec_col1]] - 0.001) %>%
    #ensure each row falls between 0 and 1
    filter(if_all(everything(), ~ . >= 0 & . <= 1))
}

# Apply transformations using map
triad_all <- map_dfr(1:6, ~ {
  case_when(
    .x == 1 ~ accumulate(1:4, ~ add_increments(.x, "R", "E"), .init = triad),
    .x == 2 ~ accumulate(1:4, ~ add_increments(.x, "R", "I"), .init = triad),
    .x == 3 ~ accumulate(1:4, ~ add_increments(.x, "E", "I"), .init = triad),
    .x == 4 ~ accumulate(1:4, ~ add_increments(.x, "E", "R"), .init = triad),
    .x == 5 ~ accumulate(1:4, ~ add_increments(.x, "I", "R"), .init = triad),
    .x == 6 ~ accumulate(1:4, ~ add_increments(.x, "I", "E"), .init = triad)
  ) %>% bind_rows()
})

# Remove duplicates
triad_all <- unique(triad_all)  # triad_all shows the proportion of the landscape 
# that is made up of each regime 

#NB - when you look at what a 0.001 step means;
#for a 1Mha landscape, regimes can vary in 1000 ha (10km2) increments

#RULES: Management regimes can go anywhere 
triad_all

#----------------------------------------------------------
#Create S.M.P (Starting landscape, management composition, 
#production target scenarios)
#----------------------------------------------------------
#DETERMINE AREA OF EACH MAGMENT REGIME ACROSS THE 1MHA SCENARIO 

#extract only the basic information on startingL to avoid duplication of scenarios below
startingL_basic <- startingL %>% 
  group_by(startingL) %>%
  mutate(area = sum(area)) %>% 
  filter(habitat =="forest") %>%  
  dplyr::select(-habitat)

# #no historical clearance 
 SL_noCL <- startingL_basic %>% 
  filter(startsWith(startingL, "NoCl"))

noCL_scenarios <- triad_all %>% crossing(SL_noCL) %>%  
  mutate(R_area = R*area, 
         E_area = E*area,
         I_area = I*area,
         ScenarioID = row_number())

#---------------------------------------------------------------------------------
# add yields attained in each scenario (which is a row) ####
#---------------------------------------------------------------------------------
yield_df

#RULE - all plantation comes from forest conversion 
noCL_scenarios <- noCL_scenarios %>% left_join(yield_df, by = "age") %>%  
  #                     convert to 1000ha (10km2) block into      #all plantations from forest lands
  mutate(totalScenario_I_m3 = (I_area/1000) *                   IForest_m3_1000ha, 
         totalScenario_E_m3 = (E_area/1000) * Ext_m3_1000ha, 
         totalScenario_R_m3 = 0,
         totalScenario_m3 = totalScenario_I_m3 + totalScenario_E_m3 + totalScenario_R_m3
  )
names(noCL_scenarios)

#-------------------------------------------------------------------------------------------
# For each scenario, determine relative (to max-even-flow) production target 
#-------------------------------------------------------------------------------------------

noCL_results   <- determine_relative_production_targets(noCL_scenarios, production_targets, "NoCl")

#only keep scenarios that are fully from intensive + R (ie sparing only for different production target)
sparing <- noCL_results %>% ungroup %>%filter(E == 0) %>%  dplyr::select(R, I, age, startingL, totalScenario_m3, production_target, composition)


#-------------------------------------------------------------------------------------------
#generate spatially explicit raster #### 
#-------------------------------------------------------------------------------------------
#For a given production target see what combination of Reserve and Plantation you would need to meet the target
print(sparing) # we select P = 0.24 an P=0.58 
#------------------------------------------------------
#loop through different patch sizes

# Create an empty raster with the desired dimensions
#r <- rast(nrows = nrows, ncols = ncols, resolution = 1)
#ext(r) <- c(0, ncols, 0, nrows)  # Define the extent

# Define the number of plantation cells (70% of total area)
nrows = 1000
ncols = 1000
total_area <-nrows*ncols #landscaoe is 1Mha and made up of 100m2 (1ha) pixels 

#Note - if i wanted 1mha made up of 30m^2 cells that would be 3333.33 x 3333.33 cells. 
#Might want to do this as we are using 100m edge density as a predictor so min cell size may need to be < 100x100
nrows = 3333
ncols = 3333
generate_landscape <- function(prop_plantation_cells, total_area, num_steps = 4000, step_size = 250) {
  # Create an empty list to store the rasters
  raster_list <- list()
  
  # Loop through different `num` of patches values from 1 to num_steps, in steps of step_size
  for (i in seq(1, num_steps, by = step_size)) {
    # Set up the landscape resolution and dimensions (same as before)
    m <- matrix(0, nrows, ncols)  # 1MHA
    r <- rast(m)
    
    # Number of patches
    num <- i 
    
    # Calculate the size for each patch
    size <- prop_plantation_cells / i  # Each patch should take an equal share of the plantation area
    
    # Create the landscape with the specified number of patches and size
    rr <- makeClass(r, num, size)  # Generate the class (patches)
    
    # Save the generated landscape to the list
    raster_list[[paste("patches_", num, sep = "")]] <- rr  # Using dynamic names for the list
    
    # Plot the generated landscape
    plot(rr)
    title(paste("Landscape with", num, "patches"))
  }
  
  # View all of the rasters in a grid
  par(mfrow = c(3, 4))  # Set the plotting area to 3 rows x 4 columns (adjust as needed)
  
  for (i in seq_along(raster_list)) {
    plot(raster_list[[i]])  # Plot the current raster in the list
    title(names(raster_list)[i])  # Add the title with the name of the raster
  }
  
  return(raster_list)
}

# Example usage:
# Define the proportion of plantation area (which determines the production target)
prop_plantation_cells024 <- 0.319 * total_area
prop_plantation_cells058 <- 0.771 * total_area

# Call the function with different plantation areas
p024 <- generate_landscape(prop_plantation_cells024, total_area)
p058<- generate_landscape(prop_plantation_cells058, total_area)

for (z in p024) {
  plot(z)
}

for (z in p058) {
  plot(z)
}


#Output all rasters tifs in a folder with the production target name

# Create the 'production_0.24' folder if it doesn't exist
output_folder24 <- "Rasters/production_0.24"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Create the 'production_0.58' folder if it doesn't exist
output_folder58 <- "Rasters/production_0.58"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Loop through the list `x` and save each raster as a .tif file
lapply(names(p024), function(name) {
  # Create a filename for each raster (e.g., "patches_1.tif")
  output_file <- file.path(output_folder24, paste0(name, ".tif"))
  
  # Save the raster as a TIFF file
  writeRaster(p024[[name]], filename = output_file, overwrite = TRUE)
})

lapply(names(p058), function(name) {
  # Create a filename for each raster (e.g., "patches_1.tif")
  output_file <- file.path(output_folder58, paste0(name, ".tif"))
  
  # Save the raster as a TIFF file
  writeRaster(p058[[name]], filename = output_file, overwrite = TRUE)
})

