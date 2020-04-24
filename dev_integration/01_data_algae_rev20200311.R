# contact: Lara Jansen (ljansen@pdx.edu)/ Eric Sokol (esokol@battelleecology.org)

# Daijiang here: I don't understand why left join parentSampleID with sampleID
# all sampleID matched but not sampleID and parentsampleID.
# source("code/00_pkg_functions.R")

# the algae data product ID:
# algae DP1.20166.001

# # all data -- could be big
# alg_allTabs <- loadByProduct(dpID = "DP1.20166.001", 
#                              site = "all", package = "expanded", 
#                              check.size = TRUE)

library(neonUtilities)
library(tidyverse)


# Try this:
# only 2 sites, restricted in time
alg_allTabs <- neonUtilities::loadByProduct(
  dpID = "DP1.20166.001", 
  site = c("MAYF", "PRIN"),
  startdate = "2016-1", 
  enddate = "2018-11",  
  package = "expanded", check.size = TRUE)

# pull out the field data into a data.frame
alg_field_data <- alg_allTabs$alg_fieldData
alg_tax_long <- alg_allTabs$alg_taxonomyProcessed
alg_biomass <-alg_allTabs$alg_biomass


## create table_observation
table_observation <- alg_tax_long %>%
  left_join(alg_biomass) %>% 
  left_join(alg_field_data) %>%
  select(
    uid,
    sampleID,
    siteID,
    collectDate,
    algalParameterValue,
    algalParameterUnit,
    algalParameter,
    perBottleSampleVolume,
    fieldSampleVolume,
    algalSampleType,
    benthicArea,
    acceptedTaxonID,
    scientificName
  ) %>%
  filter(algalParameterUnit == 'cellsPerBottle') %>% # filter to cells per bottle
  mutate(
    density = case_when(
        algalSampleType %in% c('seston') ~ algalParameterValue / perBottleSampleVolume,
        TRUE ~ (algalParameterValue / perBottleSampleVolume) * (fieldSampleVolume / benthicArea) #add phytoplankton back in when applicable
      ),
    cell_density_standardized_unit = case_when(
      algalSampleType == 'phytoplankton' ~ 'cells/mL',
      TRUE ~ 'cells/m2'
    )
  )
