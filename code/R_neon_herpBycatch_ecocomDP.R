# R_EXAMPLE_neon_macroinverts_ecocomDP_rev20191015.R
# modified to get herp bycatch from the pitfall traps
# author Matthew Helmus 

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

include_zeros <- FALSE # Still need to figure out the zeros... is there a NULL taxonID?

#################################################################################
# beetle dpid to get herp bycatch
my_dpid <- 'DP1.10022.001'
my_site_list <- c('BART')

all_tabs <- neonUtilities::loadByProduct(
  dpID = my_dpid,
  site = my_site_list,
  check.size = FALSE,
  package = "expanded")

# get location and lat long
bet_fielddata <- all_tabs$bet_fielddata

# get actual sampling data
herp_sorting<- all_tabs$bet_sorting %>% filter(sampleType == "vert bycatch herp")

# index <- intersect(colnames(herp_sorting), colnames(bet_fielddata))
if(include_zeros) {
  index <- unique(bet_fielddata$namedLocation)
} else { 
  index <- intersect(herp_sorting$namedLocation, bet_fielddata$namedLocation)
}

#inv_taxonomyProcessed <- all_tabs$bet_expertTaxonomistIDProcessed
herp_taxonomyProcessed <- herp_sorting %>% 
  dplyr::select(taxonID, scientificName, taxonRank) %>% 
  unique()

# REQUIRED TABLES -- format for 

# location
table_location <- bet_fielddata %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  ) %>% filter(location_id %in% index)


# taxon
table_taxon <- herp_taxonomyProcessed %>%
  select(taxonID, taxonRank, scientificName) %>%
  distinct() %>%
  rename(taxon_id = taxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName)

# observation
table_observation <- herp_sorting %>% 
  select(uid,
         sampleID,
         namedLocation, 
         collectDate,
         individualCount,
         taxonID) %>%
  rename(value = individualCount,
         observation_id = uid,
         event_id = sampleID,
         # package_id = NA,
         location_id = namedLocation,
         observation_datetime = collectDate,
         taxon_id = taxonID) %>%
  mutate(variable_name = 'abundance',
         unit = 'count per pitfall trap',
         package_id = NA) %>%
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
