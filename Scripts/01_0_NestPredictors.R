#organise nest and rangewide murrelet data ####

library(tidyverse)
library(lme4)

#read in inputs 
nest_oregon <- read_excel('Inputs/nest_data/NestandInlandSiteSummary_2016-2022_2023-06-01.xlsx') #oregon 
nest_params <- read_csv("Inputs/nest_data/LandscapeData.csv") #ethans nest characteristics
all_captures <- read_excel("Inputs/nest_data/CaptureDatabase_2017-2022_2022-12-13.xlsx") %>%  
  rename(usgs_band = "USGS Band #")
rangewide_nests <- read_excel("Inputs/nest_data/AllMAMUNestDatabase_2023-11-03.xlsx")

#site info !!!!!!!!!
#NB: prob best not to trust Ethan's habitat calculations; they are for a fixed year (2019) SDM of quality, and the 100m ones are calculated wrong. 

#get tag number by site
site_by_tag <- nest_oregon %>% select(SiteName, UtmE, UtmN, BandID_TaggedBird1) %>% 
  rename(usgs_band = BandID_TaggedBird1)

#add individual bird data to Ethan's summary
nest_params <- nest_params %>% left_join(site_by_tag) 

#add previous ocean condition 
nest_params <- nest_params %>% left_join(PC1_t1, by = c("Year"= "year_t1")) %>% 
  as.data.frame()

#take the nests where we know the fates 
nests37 <- nest_params %>% filter(!is.na(Fate))

#how many nests surveyed per year? 

nests37 %>% group_by(Year) %>% count()

#quick model 
model1 <- glmer(survival ~ PC1_t1 + (1 | SiteName), data = nests37, family = "binomial")
summary(model1)
glm(survival ~ scale(edgeAmount5km), data = nests37, family = "binomial")
glm(survival ~ scale(edgeAmount5km) + PC1_t1, data = nests37, family = "binomial")
glm(survival ~ scale(edgeAmount5km) * PC1_t1, data = nests37, family = "binomial")

# site_hab_features<- nest_params %>% select(Year, SiteName, 
#                                habAmount500m, habAmount1km,
#                                edgeAmount500m,edgeAmount1km)
# all_captures %>% left_join(site_hab_features, by = "")

unique(nest_oregon$BandID_TaggedBird1)
unique(all_captures$`USGS Band #`)
names()
names(nest_params)


