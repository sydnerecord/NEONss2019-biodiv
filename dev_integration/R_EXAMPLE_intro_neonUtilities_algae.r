# contact: Eric Sokol (esokol@battelleecology.org)

# the NEON package from CRAN
library(neonUtilities)

# the algae data product ID:
# algae DP1.20166.001

# # all data -- could be big
# alg_allTabs <- loadByProduct(dpID = "DP1.20166.001", 
#                              site = "all", package = "expanded", 
#                              check.size = TRUE)

# Try this:
# only 2 sites, restricted in time
alg_allTabs <- loadByProduct(dpID = "DP1.20166.001", 
                             site = c("MAYF", "PRIN"),
                             startdate = "2016-1", enddate = "2018-11",  
                             package = "expanded", check.size = TRUE)

# pull out the field data into a data.frame
alg_field_data <- alg_allTabs$alg_fieldData

# want to use acceptedTaxonID
# consider taxon resolution
# consider algalParameterUnit
# algalParameterValue is the "count" or "density
# pull out the taxonomy data
alg_tax_long <- alg_allTabs$alg_taxonomyProcessed

# or pull by data table...
alg_fieldData_v2 <- neonUtilities::getDatatable(
  dpid = 'DP1.20166.001',
  data_table_name = 'alg_fieldData',
  sample_location_list = c("MAYF", "PRIN"))
##look at data structure
str(alg_tax_long)
head(alg_tax_long)
summary(alg_tax_long)
str(alg_field_data)
head(alg_field_data)

####make workable tables ###
# location
table_location <- alg_field_data %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  )

View(table_location)
# taxon
table_taxon <- alg_tax_long %>%
  select(acceptedTaxonID, taxonRank, scientificName) %>%
  distinct() %>%
  rename(taxon_id = acceptedTaxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName)
View(table_taxon)
# observation
str(alg_tax_long)
str(alg_field_data)
summary(alg_tax_long$algalParameterUnit)


#break up by algal sample type-benthic vs. phyto

str(alg_field_data)
alg_tax_long %>%
  mutate(
type=case_when(
  algalSampleType %in% c('seston', 'phytoplankton') ~ algalParameterValue / perBottleSampleVolume,
  TRUE ~ (algalParameterValue / perBottleSampleVolume) * (fieldSampleVolume_taxonomy / benthicArea)
),cell_density_standardized_unit = case_when(
  algalSampleType == 'phytoplankton' ~ 'cells/mL',
  TRUE ~ 'cells/m2'))

summary(alg_field_data$algalSampleType)
##
str(alg_tax_long)
str(alg_field_data)

#parentSampleID vs sampleID
alg_field_data$sampleID<-alg_field_data$parentSampleID
View(alg_field_data)
#issue of child records?
head(alg_field_data,3)
head(alg_tax_long,3)
summary(alg_tax_long)
#
names(table_observation)
head(table_observation0)
labels(table_observation1)
table_observation1 <- alg_tax_long %>% 
left_join(alg_field_data, by = c("sampleID" = "parentSampleID")) %>%
  select(uid.x,
         sampleID,
         siteID.x, 
         collectDate.x,
         algalParameterValue,
         algalParameterUnit,
         algalParameter,
         perBottleSampleVolume,
         fieldSampleVolume,
         algalSampleType,
         benthicArea,
         acceptedTaxonID) %>%
      mutate(density=
             case_when(
           algalSampleType %in% c('seston') ~ algalParameterValue / perBottleSampleVolume,
           TRUE ~ (algalParameterValue / perBottleSampleVolume) * (fieldSampleVolume / benthicArea) #add phytoplankton back in when applicable
         ),cell_density_standardized_unit = case_when(
           algalSampleType == 'phytoplankton' ~ 'cells/mL',
           TRUE ~ 'cells/m2'))
View(table_observation1)

names(table_observation)
levels(table_observation1$algalSampleType)
###Aggregate by bout ####
head(table_observation1)
summary(table_observation0$density)
summary(table_observation0$algalSampleType)
#test mutate to create density
test<-table_observation %>%
  mutate(density=algalParameterValue/perBottleSampleVolume)
head(test)
  #mutate(package_id = NA) %>%
  #select(observation_id, event_id, package_id,
        # location_id, observation_datetime,
         #taxon_id, variable_name, value, unit)

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
