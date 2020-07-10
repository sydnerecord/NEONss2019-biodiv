# Reformat NEON Tick data to EcocomDP formatting
# Authors: Wynne Moss, Melissa Chen, Brendan Hobart, Matt Bitters
# Date: 7/10/2020

### Script to reformat NEON tick abundances
# uses the data from script 01_data_tick.R (downloads, cleans, links NEON tick data)
# follows same basic protocol as Mosquito data (Natalie Robinson) 


### To do:
# correct for tick drag length
# separate tables for life stages?


### Load packages
# devtools::install_github("https://github.com/EDIorg/ecocomDP")
library(ecocomDP)


### Read in data
# see script 01 for details
tick_long <- readRDS("data/tck_longform.Rdata")
tick_site_env <- readRDS("data/tck_sitexspecies_env.Rdata")

### Create location table

tick_location <- tick_site_env %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  )

# 284 tick plots (matches the env data; each row corresponds to a tick drag plot)
nrow(tick_location)
length(unique(tick_location$location_id))
length(unique(tick_site_env$namedLocation))
length(unique(tick_site_env$plotID))


### Taxon table
tick_taxon <- tick_long %>%
  select(acceptedTaxonID, taxonRank, scientificName) %>%
  rename(taxon_id = acceptedTaxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName) %>%
  distinct() 

tick_taxon
