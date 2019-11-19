# contact: Lara Jansen (ljansen@pdx.edu)/ Eric Sokol (esokol@battelleecology.org)

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
head(alg_field_data,2)
head(alg_tax_long,2)

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

##add biomass table for table_observation as it has additional estimated 'perBottleSampleVolumes'
alg_biomass <-alg_allTabs$alg_biomass
head(alg_biomass,2)

#add up labSampleVolume + preservativeVolume
alg_biomass1t<- alg_biomass %>%
  mutate(estBSVolume=preservativeVolume+labSampleVolume) %>%
  filter(analysisType=='taxonomy')
head(alg_biomass1t,2)
dim(alg_biomass1t)

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
#### observation ###

#ADD biomass data as additional 'perBottleSampleVolume'
###change NA's to 0's in alg_field data
head(alg_field_data,3)
head(alg_tax_long,2)
summary(alg_tax_long)
alg_tax_long$perBottleSampleVolume[is.na(alg_tax_long$perBottleSampleVolume)] <- 0
##add biomass estBSVolume 
head(alg_biomass1,2)

readr::write_csv(
  alg_biomass1,
  'algbiomass.csv')
alg_tax_biomass <-alg_tax_long %>%
  left_join(alg_biomass1t, by=c("sampleID" = "parentSampleID")) %>%
  mutate(perBSVol=
           case_when(
             perBottleSampleVolume ==0 ~ estBSVolume,
             perBottleSampleVolume >0 ~ perBottleSampleVolume))

##create table_observation

table_observation1 <- alg_tax_biomass %>% 
left_join(alg_field_data, by = c("sampleID" = "parentSampleID")) %>%
  select(uid.x,
         sampleID,
         siteID.x, 
         collectDate.x,
         algalParameterValue,
         algalParameterUnit,
         algalParameter,
         perBSVol,
         fieldSampleVolume.x,
         algalSampleType,
         benthicArea,
         acceptedTaxonID) %>%
  filter(algalParameterUnit=='cellsPerBottle')%>%
      mutate(density=
             case_when(
           algalSampleType %in% c('seston') ~ algalParameterValue / perBSVol,
           TRUE ~ (algalParameterValue / perBSVol) * (fieldSampleVolume.x / benthicArea) #add phytoplankton back in when applicable
         ),cell_density_standardized_unit = case_when(
           algalSampleType == 'phytoplankton' ~ 'cells/mL',
           TRUE ~ 'cells/m2'))
#check
head(table_observation1,2)



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
###Summary/Overall density by bout, site and algal sample type####
library(tidyr)
df <- na.omit(df)
summary_tableobs<-table_observation1 %>%
  group_by(siteID.x,collectDate.x,algalSampleType) %>%
  filter(!is.na(algalSampleType)) %>%
  summarise(overall_density=sum(density,NA,na.rm = TRUE))

levels(table_observation1$siteID.x)
levels(summary_tableobs$siteID.x)

readr::write_csv(
  summary_tableobs,
  'table_sumtotal.csv')