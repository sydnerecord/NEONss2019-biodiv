# title: "Tick Pathogen ecocomDP"
# author: Melissa Chen, Wynne Moss, Brendan Hobart, Matt Bitters
# date: 9 July 2020

library("tidyverse")

#### Load ####
tick_path <- read.csv("data/tick_pathogen.csv")
### Required columns #### 
req_col_observation <- c("observation_id"
                         , "event_id", "package_id", "location_id","taxon_id"
                         , "observation_datetime"
                         , "variable_name","value","unit")
req_col_location <- c("location_id"
                      , "location_name"
                      , "latitude", "longitude", "elevation")

req_col_taxon <- c("taxon_id"
                   , "taxon_rank","taxon_name")

req_col_location_ancillary <- c("location_ancillary_id", "location_id"
                                , "variable_name","value") # for nlcd class

req_col_observation_ancillary <- c("observation_ancillary_id", "event_id"
                                   , "variable_name", "value") # for host species

#### Split by host type ####
tick_path2 <- tick_path %>% separate(subsampleID, into=c("plotID2","date2","host","lifestage"), remove=FALSE, sep = "[.]")

#### Observation table ####
dpID <- "DP1.10092.001" # package ID

observation_raw <- data.frame(event_id = tick_path2$subsampleID)
observation_raw$package_id <- dpID
observation_raw$location_id <- tick_path2$namedLocation
observation_raw$taxon_id <- tick_path2$testPathogenName
observation_raw$observation_datetime <- tick_path2$collectDate # should be UTC, or "Z"
observation_raw$positive <- ifelse(tick_path2$testResult=="Positive",1,0)
observation_raw$host <- tick_path2$host
observation_raw$lifeStage <- tick_path2$lifestage

observation_collapsed <- observation_raw %>% group_by(event_id, package_id, location_id, taxon_id, host, observation_datetime, lifeStage) %>%
  summarize(numberTested = n(), value = sum(positive)/numberTested, variable_name="InfectionRate",unit="proportion") %>%
  ungroup() %>%mutate(observation_id = paste0("obs",row_number())) 

observation <- observation_collapsed %>% select(one_of(req_col_observation))

#### Observation Ancillary table ####
observation_ancillary <- observation_collapsed %>% select(event_id, host, numberTested, lifeStage) %>% 
  mutate(numberTested=as.character(numberTested)) %>%
  pivot_longer(cols=c(host,numberTested,lifeStage), names_to = "variable_name", values_to="value") %>%
  mutate(observation_ancillary_id = paste0("oban",row_number())) %>%
  select(one_of(req_col_observation_ancillary))

#### Location table ####
location_raw <- data.frame(location_name = tick_path2$namedLocation)
location_raw$latitude <- tick_path2$decimalLatitude
location_raw$longitude <- tick_path2$decimalLongitude
location_raw$elevation <- tick_path2$elevation
location_raw$latlongUncertainty <- tick_path2$coordinateUncertainty
location_raw$elevationUncertainty <- tick_path2$elevationUncertainty
location_raw$domainID <- tick_path2$domainID
location_raw$siteID <- tick_path2$siteID
location_raw$nlcdClass <- tick_path2$nlcdClass

location_collapsed <- distinct(location_raw) %>% mutate(location_id = paste0("loc",row_number()))

location <- location_collapsed %>% select(one_of(req_col_location))

#### Location Ancillary table ####

location_ancillary <- location_collapsed %>% select(location_id, latlongUncertainty, elevationUncertainty, domainID, siteID, nlcdClass) %>%
  mutate(latlongUncertainty=as.character(latlongUncertainty), elevationUncertainty=as.character(elevationUncertainty)) %>%
  pivot_longer(cols=c(latlongUncertainty, elevationUncertainty, domainID, siteID, nlcdClass), names_to = "variable_name", values_to="value") %>%
  mutate(location_ancillary_id=paste0("loan",row_number())) %>% 
  select(one_of(req_col_location_ancillary))

#### Taxon table ####
req_col_taxon

taxon_raw <- data.frame(taxon_name=tick_path2$testPathogenName)

# Taxonomies are irregular; look for keywords to find rank
genus <- grep(" sp.$", taxon_raw$taxon_name)
taxon_raw$taxon_rank <- "species"
taxon_raw$taxon_rank[genus] <- "genus"

taxon_collapsed <- distinct(taxon_raw) %>% mutate(taxon_id=paste0("tax",row_number()))

taxon <- taxon_collapsed %>% select(one_of(req_col_taxon))

## Should I include taxon ancillary with which hosts each tends to infect?


write.csv(observation, file="data/observation_tickpathogen.csv", row.names = FALSE, quote=FALSE)
write.csv(observation_ancillary, file="data/observation_ancillary_tickpathogen.csv", row.names = FALSE, quote=FALSE)
write.csv(location, file="data/location_tickpathogen.csv", row.names = FALSE, quote=FALSE)
write.csv(location_ancillary, file="data/location_ancillary_tickpathogen.csv", row.names = FALSE, quote=FALSE)
write.csv(taxon, file="data/taxon_tickpathogen.csv", row.names = FALSE, quote=FALSE)

