# MarbledEdges - Gianluca Cerullo, 2nd June 2025

These scripts are to replicate in Rstudio an analysis exploring how ocean conditions and edge amount influence marbled murrelet occupancy. The scripts should be run in consecutive order. After downloading the MarbledEdges.Rproj project, you will need to reproduce the following file structure for scripts to run smoothly:

-MarbledEdges.Rproj

-Scripts

-Inputs

-Outputs

-Models

-Rasters

-Figures

#Scripts: 

#-------------------------
#01_OceanPredictors
#-------------------------

#script description 
Downloads, subsets and scales the PC1 variable denoting ocean conditions in 
the Pacific Northwest based on climate and biological variables. 

#inputs: 
#16 PCA variables; downloaded from: 
#https://www.fisheries.noaa.gov/west-coast/science-data/ocean-conditions-indicators-trends

#outputs
Outputs/PC1_scaled_inverted.csv

#-------------------------
#02_OrganiseMurreletData
#-------------------------

#script description 
Organises murrelet survey data and covariates ready for occupancy modelling

#Inputs
Inputs/LandscapeVariables.csv
Inputs/siteData.csv"
Inputs/MurreletSurveys.csv")
Inputs/pointsInRoiAndSamplingWindow.csv
Outputs/PC1_scaled_inverted.csv" 

#Outputs
Outputs/analysisDataUnmarked.rds #unmarked-ready rds file

#-------------------------
#03_RunOccupancyModels
#-------------------------

#script description: 
Fit occupancy models, carry out model selection and comparison, carry out k-fold validation, 
produce model performance tables

#Inputs
Outputs/analysisDataUnmarked.rds

#Outputs
Models/final_model_5thMay2025.rds #best performing model
Models/Tables/model_performance_table.csv #table of model performance

#-------------------------
#04_ExtractPredictionsForPlotting
#-------------------------

#Script Description
Uses best-fitting model to plot impacts edge amount and ocean condition on murrlet occupancy

#Inputs 
source("scripts/01_DeterminePredictors.R")
source("scripts/02_OrganiseMurreletData.R")
Outputs/analysisDataUnmarked.rds
Models/final_model_5thMay2025.rds

#Outputs
Fig1

#-------------------------
#05_BuildLandscapes
#-------------------------

#Script Description 
Builds hypothetical 1Mha simulated logging landscapes (outputted as rasters) 
matched by timber production target but varying in fragmentation per se

#Inputs
source("Inputs/Scenarios/Functions.R")
source("Inputs/Scenarios/CalculateYields.R")

#Outputs 
NB:rasters for different landscapes are stored in folder denoting the proportional
prodcution target. Make sure you have a folder called "Rasters"

Rasters/production_0.24
Rasters/production_0.55
Rasters/production_0.58

#-------------------------
#06_ExtractSimulatedLandscapeParams
#-------------------------

#Script Description 
Quantifies edge amount and habitat amount around forest points for simulated 
1Mha landscapes
NB: There is a hard-coded decision to make regarding which production target landscape to explore. 

#Outputs
Outputs/EdgeDensity10000pts_SimulatedParams.rds
Outputs/HabAmount10000pts_SimulatedParams.rds

#-------------------------
#07_OrganiseSimulatedCovariates
#-------------------------

#Script Description 
Use best-fitting model to quantify the occupancy of forest points for simulated 1Mha landscapes

#Inputs
Models/pc1_interaction_model.rds

#Outputs
Figure 2, showing occupancy predictions for simulated landscapes

#-------------------------
#08_BackcastFragmentationOceanImpacts
#-------------------------

#Script Description 
Stratify points across the PNW for the year 2020 and extract model covariates 

#Inputs
vect("Rasters/land_ownership/All_merge.shx") #landownership shapefile

#GNN raster data must be requested at: https://lemmadownload.forestry.oregonstate.edu/
#canopy cover of all trees
Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_2020.tif
#canopy cover of conifers
Rasters/GNN_2021/2025_02_11_cerullo/rasters/cancov_con_2020.tif
#SDM of murrelets for year 2020
Rasters/MAMU_SDMs/MAMU_maxent_cloglog_2020_20220208.tif

#Outputs (there are temporary files produced which are not listed here)
Outputs/elevation_at_each_point.csv #elevation at all points
Outputs/PNW_2020_extracted_covars.csv #covariates for all points


#-------------------------
#09_BackcastFragmentationOceanImpacts
#-------------------------

#Script Description 
Scales covariates for filtered points across the PNW, and then uses best-fitting model
to predict occupancy for each point under current 2020 configurations and with different
amounts of edge reduction. Also, plots covariates (e.g. edge amount and habitat amount) for 2020 lands.

#Outputs
Outputs/log_odds_se_05edge_reduction.rds #log odds of occupancy from reducing edge 50%
Outputs/ReduceLandscapeAndLocalEdges.rds #what happens if we reduce edge density

#-------------------------
#10_MapEdgeReductionPrioritiesLogOdds
#-------------------------
#Script Description 
Plot priority opportunities by actor and map priorities with standard errors 

#Outputs
Rasters/GoodYear_odds_top90thpercentilOccupancyIncreaseFromReducingedge50pc.tif
Rasters/BadYear_odds_top90thpercentilOccupancyIncreaseFromReducingedge50pc.tif"
Rasters/GoodAndBad_odds_top90thpercentilOccupancyIncreaseFromReducingedge50pc.tif"

Figures showing how occupancy is distributed by actor  


#-------------------------
#10B_RawOccDifferences
#-------------------------
#Script Description 
Plots binary and raw occupancy differences with SE

#Output

Figures and maps showing how absolute occupancy increase is distributed by actor if we reduce edges 50%  




#-------------------------
#11_productivity_by_ownership
#-------------------------
#Script Description 
Shows how potential timber yields vary with increasing coastal distance  

#Inputs
#potential yield data are downloaded from: https://www.arcgis.com/home/item.html?id=b1b10e7890ea4116b863ae790d9b718c

#Ouptus
Figure showing timber yield potential with coastal distance


