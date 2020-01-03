#title: ecocomDP-algae
#author: Lara Jansen
#date: 1/02/2020

#packages used
library(dplyr)
library(neonUtilities)


#get algae data from 01_data_algae

#Location table
table_location <- alg_field_data %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  )
#Taxon table
table_taxon <- alg_tax_long %>%
  select(acceptedTaxonID, taxonRank, scientificName) %>%
  distinct() %>%
  rename(taxon_id = acceptedTaxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName)
View(table_taxon)

#Observation table

###change NA's to 0's in alg_tax_long data for calculations only
alg_tax_long$perBottleSampleVolume[is.na(alg_tax_long$perBottleSampleVolume)] <- 0
#join algae biomass and taxonomy data 
alg_tax_biomass <-alg_biomass1t %>%
  left_join(alg_tax_long, by=c("parentSampleID" = "sampleID")) %>%
  mutate(perBSVol=
           case_when(
             perBottleSampleVolume ==0 ~ estBSVolume,
             perBottleSampleVolume >0 ~ perBottleSampleVolume))
#create the table
table_observation1 <- alg_tax_biomass2 %>% 
  left_join(alg_field_data, by = c("sampleID" = "parentSampleID")) %>%
  select(uid.x,
         sampleID,
         siteID.x, 
         collectDate.x,
         algalParameterValue,
         algalParameterUnit,
         algalParameter,
         perBSVol,
         fieldSampleVolume.y,
         algalSampleType,
         benthicArea,
         acceptedTaxonID) %>%
  filter(algalParameterUnit=='cellsPerBottle')%>%
  mutate(density=
           case_when(
             algalSampleType %in% c('seston') ~ algalParameterValue / perBSVol,
             TRUE ~ (algalParameterValue / perBSVol) * (fieldSampleVolume.y / benthicArea) #add phytoplankton back in when applicable
           ),cell_density_standardized_unit = case_when(
             algalSampleType == 'phytoplankton' ~ 'cells/mL',
             TRUE ~ 'cells/m2'))

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
  table_observation1,
  'table_observation.csv')