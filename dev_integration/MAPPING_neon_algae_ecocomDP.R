#title: ecocomDP-algae
#author: Lara Jansen
#date: 1/02/2020

#packages used
library(dplyr)
library(neonUtilities)

# # in user provided
# neon_site_id = c("MAYF", "PRIN"),
# data_startdate = "2016-1", 
# data_enddate = "2018-11",  
# data_package_type = "expanded", 
# check.size = FALSE


# in hard coded by data product
dpID <- "DP1.20120.001"

metadata_all <- neonUtilities::getProductInfo(my_dpid)

landing_page_url <- paste0("https://data.neonscience.org/data-products/", my_dpid)

#################

# Try this:
# only 2 sites, restricted in time
all_tabs_in <- neonUtilities::loadByProduct(
  dpID = "DP1.20166.001", 
  site = c("MAYF", "PRIN"),
  startdate = "2016-1", 
  enddate = "2018-11",  
  package = "expanded", 
  check.size = FALSE)

# note, scientific name encoding is messed up

# variables_in <- all_tabs_in[[grep('variables',names(all_tabs_in))]]
# variables_in$table %>% unique()

field_data_in <- all_tabs_in$alg_fieldData
tax_long_in <- all_tabs_in$alg_taxonomyProcessed
biomass_in <-all_tabs_in$alg_biomass


#get algae data from 01_data_algae

#Location table
table_location <- field_data_in %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  )
#Taxon table
table_taxon <- tax_long_in %>%
  select(acceptedTaxonID, taxonRank, scientificName) %>%
  distinct() %>%
  rename(taxon_id = acceptedTaxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName)
# View(table_taxon)

#Observation table

###change NA's to 0's in tax_long_in data for calculations only
# tax_long_in$perBottleSampleVolume[is.na(tax_long_in$perBottleSampleVolume)] <- 0
#join algae biomass and taxonomy data 
alg_tax_biomass <- biomass_in %>%
  left_join(tax_long_in, by=c("parentSampleID" = "sampleID")) 

# alg_tax_biomass <-biomass_in1t %>%
#   left_join(tax_long_in, by=c("parentSampleID" = "sampleID")) %>%
#   mutate(perBSVol=
#            case_when(
#              perBottleSampleVolume ==0 ~ estBSVolume,
#              perBottleSampleVolume >0 ~ perBottleSampleVolume))



#create the table
table_observation <- alg_tax_biomass %>% 
  left_join(field_data_in, by = c("sampleID" = "parentSampleID")) %>%
  select(uid.x,
         sampleID,
         siteID.x, 
         collectDate.x,
         algalParameterValue,
         algalParameterUnit,
         algalParameter,
         perBottleSampleVolume,
         fieldSampleVolume.y,
         algalSampleType,
         benthicArea,
         acceptedTaxonID) %>%
  filter(algalParameterUnit=='cellsPerBottle')%>%
  mutate(density=
           case_when(
             algalSampleType %in% c('seston') ~ algalParameterValue / perBottleSampleVolume,
             TRUE ~ (algalParameterValue / perBottleSampleVolume) * (fieldSampleVolume.y / benthicArea) #add phytoplankton back in when applicable
           ),cell_density_standardized_unit = case_when(
             algalSampleType == 'phytoplankton' ~ 'cells/mL',
             TRUE ~ 'cells/m2'))

out_list <- list(
  metadata_neon = metadata_all,
  table_location = table_location,
  table_taxon = table_taxon,
  table_observation = table_observation)

# ##################################
# # write out in ecocomDP format
# ###
# readr::write_csv(
#   table_location,
#   'table_location.csv')
# 
# readr::write_csv(
#   table_taxon,
#   'table_taxon.csv')
# 
# readr::write_csv(
#   table_observation,
#   'table_observation.csv')


# # troubleshooting encoding issues for sci names
# table_taxon %>% filter(taxon_id == 'NEONDREX6001')
# restR::find.scientific.name(type = 'ALGAE',taxonID = 'NEONDREX6001')
# tax_long_in$scientificName %>% unique()
# all_tabs_in$alg_taxonomyRaw %>% filter(taxonID == 'NEONDREX6001') %>% slice(1) %>% select(scientificName)
