#title: ecocomDP-zooplankton
#author: Stephanie Parker
#date: 06/10/2020

#packages used
library(tidyverse)

# Install and load devtools
# install.packages("devtools")
library(devtools)

# Install and load dev version of ecocomDP
# install_github("EDIorg/ecocomDP", ref = 'development')
library(ecocomDP)

# Install and load neonUtilities
# install_github("NEONScience/NEON-utilities/neonUtilities", dependencies=TRUE)
library(neonUtilities)

#get zooplankton data from 01_data_zooplankton

# REQUIRED TABLES -- format for 

# Location table
table_location <- zoo_fielddata %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  )

# Taxon table
table_taxon <- zoo_taxonomyProcessed %>%
  select(taxonID, taxonRank, scientificName) %>%
  distinct() %>%
  rename(taxon_id = taxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName)

# Observation table
table_observation <- zoo_taxonomyProcessed %>% 
  select(uid,
         sampleID,
         namedLocation, 
         collectDate,
         zooVolumePerBottle,
         zooSubsampleVolume,
         individualCount,
         adjCountPerBottle,
         taxonID) %>%
  left_join(zoo_fielddata %>% select(sampleID, towsTrapsVolume)) %>%
  mutate(variable_name = 'density',
         value = adjCountPerBottle / towsTrapsVolume,
         unit = 'count per liter') %>% rename(observation_id = uid,
                                              event_id = sampleID,
                                              # package_id = NA,
                                              location_id = namedLocation,
                                              observation_datetime = collectDate,
                                              taxon_id = taxonID) %>%
  mutate(package_id = NA) %>%
  select(observation_id, event_id, package_id,
         location_id, observation_datetime,
         taxon_id, variable_name, value, unit)

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
  table_observation,
  'table_observation.csv')
