#organise nest and rangewide murrelet data ####

library(tidyverse)
library(lme4)

#read in inputs 
nest_oregon <- read_excel('Inputs/nest_data/NestandInlandSiteSummary_2016-2022_2023-06-01.xlsx') #oregon 
nest_params <- read_csv("Inputs/nest_data/LandscapeData.csv") #ethans nest characteristics
all_captures <- read_excel("Inputs/nest_data/CaptureDatabase_2017-2022_2022-12-13.xlsx") %>%  
  rename(usgs_band = "USGS Band #")
rangewide_nests <- read_excel("Inputs/nest_data/AllMAMUNestDatabase_2023-11-03.xlsx")

#get tag number by site
site_by_tag <- nest_oregon %>% select(SiteName, UtmE, UtmN, BandID_TaggedBird1) %>% 
  rename(usgs_band = BandID_TaggedBird1)

#add individual bird data to Ethan's summary
nest_params <- nest_params 
#site info 
#NB: prob best not to trust Ethan's habitat calculations; they are for a fixed year (2019) SDM of quality, and the 100m ones are calculated wrong. 

# site_hab_features<- nest_params %>% select(Year, SiteName, 
#                                habAmount500m, habAmount1km,
#                                edgeAmount500m,edgeAmount1km)
# all_captures %>% left_join(site_hab_features, by = "")



unique(nest_oregon$BandID_TaggedBird1)
unique(all_captures$`USGS Band #`)
names()
names(nest_params)
