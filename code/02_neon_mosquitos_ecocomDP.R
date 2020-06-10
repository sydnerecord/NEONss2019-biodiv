#title: format mosquito data into ecocomDP tables
#author: Natalie Robinson
#date: 6/10/2020

#packages used
library(dplyr)

# run 01_data_mosquitos to get source dataframe

#################################################################################
# Format data for ecocomDP

# location table
table_location <- mos_dat %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  )

# taxon table 
table_taxon <- mos_dat %>%
  select(taxonID, taxonRank, scientificName) %>%
  rename(taxon_id = taxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName) %>%
  distinct() 

# observations table
table_observations <- mos_dat %>% 
  rename(observation_id = uid,
         event_id = eventID,
         location_id = namedLocation,
         latitude = decimalLatitude,
         longitude = decimalLongitude,
         observation_datetime = startCollectDate,
         taxon_id = taxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName) %>%
  mutate(package_id = NA) #%>%


###################################
# write out in ecocomDP format
###
readr::write_csv(
  table_location,
  'table_location.csv')

readr::write_csv(
  table_taxon,
  'table_taxon.csv')

readr::write_csv(
  table_observations,
  'table_observations.csv')
