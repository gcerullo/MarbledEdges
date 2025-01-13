# This code is for detailing key paramatres for building scenarios for different 
# starting landscapes for the Elliott 

library(terra)
library(data.table)
#library(rflsgen) #for generating landscapes with different fragmentation patterns
library(tidyverse)
library(landscapeR)
library(terra)

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
#NB; these vary in the amounts of historical clearance, 
#and in the age of remaining forest at the scenario start 

{# {
#   #High historical clearance landscapes 
#   HighCl_young <-  tibble(
#     age = c(0, 90),
#     habitat = c("clearcut", "forest"),
#     area = c(landscape_size_ha*0.5, landscape_size_ha*0.5), 
#     startingL = "HighCl_young") 
#   
#   HighCl_medium <- tibble(
#     age = c(0, 150),
#     habitat = c("clearcut", "forest"),
#     area = c(landscape_size_ha*0.5, landscape_size_ha*0.5),
#     startingL = "HighCl_medium")
#   
#   HighCl_old <- tibble(
#     age = c(0, 300),
#     habitat = c("clearcut", "forest"),
#     area = c(landscape_size_ha*0.5, landscape_size_ha*0.5),
#     startingL = "HighCl_old")
#   
#   #Medium historical clearance landscapes
#   MedCl_young <-  tibble(
#     age = c(0, 90),
#     habitat = c("clearcut", "forest"),
#     area = c(landscape_size_ha*0.25, landscape_size_ha*0.75), 
#     startingL = "MedCl_young") 
#   
#   MedCl_medium <- tibble(
#     age = c(0, 150),
#     habitat = c("clearcut", "forest"),
#     area = c(landscape_size_ha*0.25, landscape_size_ha*0.75), 
#     startingL = "MedCl_medium")
#   
#   MedCl_old <- tibble(
#     age = c(0, 300),
#     habitat = c("clearcut", "forest"),
#     area = c(landscape_size_ha*0.25, landscape_size_ha*0.75), 
#     startingL = "MedCl_old")
#   
#   #No historical clearance landscapes 
#   NoCl_young <-  tibble(
#     age = c(90),
#     habitat = "forest",
#     area = landscape_size_ha, 
#     startingL = "NoCl_young") 
# 
# NoCl_medium <- tibble(
#   age = c(150),
#   habitat = "forest",
#   area = landscape_size_ha, 
#   startingL = "NoCl_medium")
}  
  NoCl_old <- tibble(
    age = c(300),
    habitat = "forest",
    area = landscape_size_ha, 
    startingL = "NoCl_old")
  
  # Combine all tibbles into a single dataframe called startingL
  startingL <- bind_rows(#HighCl_young, HighCl_medium, HighCl_old, 
                         #MedCl_young, MedCl_medium, MedCl_old, 
                         # NoCl_young, NoCl_medium,
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
#   select(startingL, max_SL_production_volume)  

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
  select(startingL,age, max_SL_production_volume)  



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
  select(-habitat)

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
sparing <- noCL_results %>% ungroup %>%filter(E == 0) %>%  select(R, I, age, startingL, totalScenario_m3, production_target, composition)
sparing 

#-------------------------------------------------------------------------------------------
#generate spatially explicit raster #### 
#-------------------------------------------------------------------------------------------
#lets start with a medium production target of 0.55 (26.9% forest,77.1% plantation) - view(sparing)

#------------------------------------------------------
#loop through different patch sizes

# Create an empty raster with the desired dimensions
#r <- rast(nrows = nrows, ncols = ncols, resolution = 1)
#ext(r) <- c(0, ncols, 0, nrows)  # Define the extent

# Define the number of plantation cells (70% of total area)
n_plantation_cells <- round(0.7 * total_area_km2)


# Create an empty list to store the rasters
raster_list <- list()

# Loop through different `num` of patches values from 1 to 1000, in steps of 100
for (i in seq(1, 500, by = 50)) {
  
  # Set up the landscape resolution and dimensions (same as before)
  nrows <- 1000
  ncols <- 1000
  m <- matrix(0,  nrows, ncols) #1MHA 
  total_area_km2 <- nrows*ncols
  r <- rast(m)
  
  # Number of patches
  num <- i 
  
  # Calculate the size for each patch
  size <- n_plantation_cells / i  # Each patch should take an equal share of the plantation area
  
  # Create the landscape with the specified number of patches and size
  rr <- makeClass(r, num, size)  # Generate the class (patches)
  
  # Save the generated landscape to the list
  raster_list[[paste("patches_", num, sep = "")]] <- rr  # Using dynamic names for the list
  
  # Plot the generated landscape
  plot(rr)
  title(paste("Landscape with", num, "patches"))
}

# View all of the rasters
par(mfrow = c(3, 4))  # Set the plotting area to 3 rows x 4 columns (adjust as needed)

for (i in seq_along(raster_list)) {
  plot(raster_list[[i]])  # Plot the current raster in the list
  #title(names(raster_list)[i])  # Add the title with the name of the raster
}


#for a single landscape


m <- matrix(0,  1000, 1000) #1MHA
r <- rast(m)
#ext(r) <- c(0, 10, 0, 10)

#rr <- makePatch(r, size=500, rast=TRUE)
plot(rr)

num <- 1000
size <- (1000000*0.7)/num
num*size
rr <- makeClass(r, num, size)
plot(rr)

?makeClass()




# 
# # Define landscape dimensions (approximately 1Mha with 30m pixels)
# nb_rows <- 3334 # Number of rows
# nb_cols <- 3333 # Number of columns
# 
# ##TEST
# cls_a <- flsgen_create_class_targets(
#   "Plantation",
#   NP = c(1, 10),  #number of patches
#   AREA = c(1, 1000000), #patch area min cells, max cells 
#   CA = c(1, 40000), #total class area min and max 
#   PLAND = c(77.1,77.1)
#  # MESH = c(225, 225) #effective mesh size
# )
# 
# cls_b <- flsgen_create_class_targets(
#   "Forest",
#   NP = c(1, 1000),
#   AREA = c(1, 100000),
# )
# 
# 
# 200*200
# ls_targets <- flsgen_create_landscape_targets(
#   nb_rows = 200, 
#   nb_cols = 200,
#   classes = list(cls_a, cls_b)
# )
# 
# structure <- flsgen_structure(ls_targets) 
# landscape <- flsgen_generate(structure_str = structure) %>%
#   #replace all -1s with 0s so that we only have two classes (third class is needed to optimse speedy creation of landscapes)
#   ifel(. == -1, 0, .)
# 
# plot(landscape)
# #
# plot(landscape)
# #--------------------------------------------------------------
# library(terra)
# library(tidyverse)
# 
# # Parameters
# n_landscapes <- 100
# landscape_area <- 10000 * 100 # 10,000 ha with 30x30m pixels
# pixel_size <- 30 # 30m x 30m resolution
# plantation_pct <- 0.7
# forest_pct <- 0.3
# 
# # Function to create a landscape with fixed coverage and increasing fragmentation
# generate_landscape <- function(fragmentation_degree) {
#   # Create an empty raster with specified resolution
#   n_pixels <- landscape_area / (pixel_size^2)
#   r <- rast(ncol = sqrt(n_pixels), nrow = sqrt(n_pixels), res = pixel_size)
#   
#   # Total number of plantation and forest pixels
#   n_plantation <- round(n_pixels * plantation_pct)
#   n_forest <- round(n_pixels * forest_pct)
#   
#   # Assign plantation (0) and forest (1) to the raster
#   r[] <- c(rep(0, n_plantation), rep(1, n_forest))
#   
#   # Shuffle the values to randomize their locations
#   r[] <- sample(r[], length(r[]), replace = FALSE)
#   
#   # Fragmentation logic: Creating fragmented forest blocks
#   forest_cells <- which(r[] == 1)  # Find forest cells
#   n_forest_cells <- length(forest_cells)
#   
#   # Define number of forest patches based on fragmentation degree
#   n_patches <- max(1, round(n_forest_cells * (1 - fragmentation_degree)))
#   
#   # Generate patch labels for the forest cells
#   patch_labels <- sample(1:n_patches, n_forest_cells, replace = TRUE)
#   
#   # Assign the patch labels to the forest cells to simulate fragmentation
#   r[forest_cells] <- patch_labels
#   
#   return(r)
# }
# 
# # Generate landscapes with increasing fragmentation
# landscapes <- map(1:n_landscapes, ~generate_landscape(. / n_landscapes))
# 
# # Inspect one landscape
# plot(landscapes[[10]])
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # # Define landscape dimensions (approximately 1Mha with 30m pixels)
# # nb_rows <- 3334 # Number of rows
# # nb_cols <- 3333 # Number of columns
# # 
# # ##TEST
# # cls_a <- flsgen_create_class_targets(
# #   "Plantation",
# #   NP = c(1, 10),  #number of patches
# #   AREA = c(1, 1000000), #patch area min cells, max cells 
# #   CA = c(1, 40000), #total class area min and max 
# #   PLAND = c(77.1,77.1)
# #  # MESH = c(225, 225) #effective mesh size
# # )
# # 
# # cls_b <- flsgen_create_class_targets(
# #   "Forest",
# #   NP = c(1, 1000),
# #   AREA = c(1, 100000),
# # )
# # 
# # 
# # 200*200
# # ls_targets <- flsgen_create_landscape_targets(
# #   nb_rows = 200, 
# #   nb_cols = 200,
# #   classes = list(cls_a, cls_b)
# # )
# # 
# # structure <- flsgen_structure(ls_targets) 
# # landscape <- flsgen_generate(structure_str = structure) %>%
# #   #replace all -1s with 0s so that we only have two classes (third class is needed to optimse speedy creation of landscapes)
# #   ifel(. == -1, 0, .)
# # 
# # plot(landscape)
# # #
# # plot(landscape)
# # #--------------------------------------------------------------
# # library(terra)
# # library(tidyverse)
# # 
# # # Parameters
# # n_landscapes <- 100
# # landscape_area <- 10000 * 100 # 10,000 ha with 30x30m pixels
# # pixel_size <- 30 # 30m x 30m resolution
# # plantation_pct <- 0.7
# # forest_pct <- 0.3
# # 
# # # Function to create a landscape with fixed coverage and increasing fragmentation
# # generate_landscape <- function(fragmentation_degree) {
# #   # Create an empty raster with specified resolution
# #   n_pixels <- landscape_area / (pixel_size^2)
# #   r <- rast(ncol = sqrt(n_pixels), nrow = sqrt(n_pixels), res = pixel_size)
# #   
# #   # Total number of plantation and forest pixels
# #   n_plantation <- round(n_pixels * plantation_pct)
# #   n_forest <- round(n_pixels * forest_pct)
# #   
# #   # Assign plantation (0) and forest (1) to the raster
# #   r[] <- c(rep(0, n_plantation), rep(1, n_forest))
# #   
# #   # Shuffle the values to randomize their locations
# #   r[] <- sample(r[], length(r[]), replace = FALSE)
# #   
# #   # Fragmentation logic: Creating fragmented forest blocks
# #   forest_cells <- which(r[] == 1)  # Find forest cells
# #   n_forest_cells <- length(forest_cells)
# #   
# #   # Define number of forest patches based on fragmentation degree
# #   n_patches <- max(1, round(n_forest_cells * (1 - fragmentation_degree)))
# #   
# #   # Generate patch labels for the forest cells
# #   patch_labels <- sample(1:n_patches, n_forest_cells, replace = TRUE)
# #   
# #   # Assign the patch labels to the forest cells to simulate fragmentation
# #   r[forest_cells] <- patch_labels
# #   
# #   return(r)
# # }
# # 
# # # Generate landscapes with increasing fragmentation
# # landscapes <- map(1:n_landscapes, ~generate_landscape(. / n_landscapes))
# # 
# # # Inspect one landscape
# # plot(landscapes[[10]])
