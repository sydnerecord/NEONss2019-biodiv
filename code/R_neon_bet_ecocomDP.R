## Author - Natalie Robinson

library(tidyverse)
library(devtools)
library(githubinstall)
library(neonUtilities)
library(dplyr)
library(tidyr)
options(stringsAsFactors = F)

# install_github("DeveloperName/PackageName") # not available for R3.6.1
#githubinstall("EDIorg/ecocomDP", ref = 'development')

# Install and load dev version of ecocomDP
# install_github("EDIorg/ecocomDP", ref = 'development') # not available for R3.6.1
# library(ecocomDP)

#################################################################################
#Get data
# bet dpid
my_dpid <- 'DP1.10022.001'
my_site_list <- c('HARV', 'ABBY')

all_tabs <- neonUtilities::loadByProduct(
  dpID = my_dpid,
  site = my_site_list,
  check.size = TRUE)

# download field data for all dates for two neon sites -- much more manageable 
bet_field <- all_tabs$bet_fielddata 
bet_sort <- all_tabs$bet_sorting 
bet_pool <- all_tabs$bet_archivepooling
bet_pin <- all_tabs$bet_parataxonomistID

#Taxonomy - required separate call - takes about few mins
bet_taxonomy_expert <- all_tabs$bet_expertTaxonomistIDProcessed

#################################################################################
# Join tables: 1) bet_pin:subsampleID == bet_sort:subsampleID, bet_sort:sampleID == sampleID
#              2) bet_pool:subsampleIDList contains bet_sort:subsampleID, bet_sort:sampleID == sampleID
bet_fieldSort <- left_join(select(bet_field,-uid),
                           select(bet_sort,etOHChangeDate,identificationQualifier,identificationReferences,identifiedBy,identifiedDate,
                                  individualCount,morphospeciesID,processingDate,recordedBy,remarks,sampleID,
                                  sampleType,subsampleCode,subsampleID,taxonID),
                           by='sampleID') %>%
  rename(recordedBy_field=recordedBy.x,
         recordedBy_sort=recordedBy.y,
         remarks_field=remarks.x,
         remarks_sort=remarks.y)




# Copied from plant div code, as template
# Remove duplicate taxa between nested subplots (each taxon should be represented once for the bout/plotID/year). 
# 1) If a taxon/date/bout/plot combo is present in 1m2 data, remove from 10/100
div_1m2_pla <- mutate(div_1m2_pla, primaryKey = paste(plotID,boutNumber,substr(endDate,1,4),taxonID,subplotID,sep='_'))
div_10_100_m2 <- mutate(div_10_100_m2, primaryKey = paste(plotID,boutNumber,substr(endDate,1,4),taxonID,subplotID,sep='_'))
for (key in unique(sub('_[^_]+$','',div_1m2_pla$primaryKey))){
  div_10_100_m2 <- filter(div_10_100_m2,sub('_[^_]+$','',primaryKey) != key)
} 

# 2) If a taxon/date/bout/plot combo is present in 10m2, remove from 100
for (key in unique(sub('_[^_]+$','',div_10_100_m2$primaryKey[nchar(div_10_100_m2$subplotID) > 2]))){
  div_10_100_m2 <- filter(div_10_100_m2,!(nchar(div_10_100_m2$subplotID) == 2 & sub('_[^_]+$','',primaryKey) == key))
} 

##Double check - in10_notIn100 + in100_notIn10 + inBoth should equal all_unique
#in10_notIn100 <- length(which(!unique(sub('_[^_]+$','',div_10_100_m2$primaryKey[nchar(div_10_100_m2$subplotID) > 2])) %in% unique(sub('_[^_]+$','',div_10_100_m2$primaryKey[nchar(div_10_100_m2$subplotID) == 2]))))
#in100_notIn10 <- length(which(!unique(sub('_[^_]+$','',div_10_100_m2$primaryKey[nchar(div_10_100_m2$subplotID) == 2])) %in% unique(sub('_[^_]+$','',div_10_100_m2$primaryKey[nchar(div_10_100_m2$subplotID) > 2]))))
#inBoth <- length(intersect(unique(sub('_[^_]+$','',div_10_100_m2$primaryKey[nchar(div_10_100_m2$subplotID) > 2])),unique(sub('_[^_]+$','',div_10_100_m2$primaryKey[nchar(div_10_100_m2$subplotID) == 2]))))
#all_unique <- length(unique(sub('_[^_]+$','',div_10_100_m2$primaryKey)))

#################################################################################
# Format data
# location
table_location_1m2 <- div_1m2_pla %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  )

table_location_10_100m2 <- div_10_100_m2 %>%
  select(namedLocation, decimalLatitude, decimalLongitude, elevation) %>%
  distinct() %>%
  rename(
    location_id = namedLocation,
    latitude = decimalLatitude,
    longitude = decimalLongitude
  )

# taxon
table_taxon <- div_taxonomyProcessed %>%
  select(taxonID, acceptedTaxonID, taxonRank, scientificName) %>%
  distinct() 

# Join and format
table_1m2observations <- div_1m2_pla %>% 
  select(uid,
         namedLocation, 
         subplotID,
         boutNumber,
         endDate,
         taxonID,
         identificationQualifier,
         morphospeciesID,
         percentCover,
         heightPlantOver300cm
         ) %>%
  left_join(table_taxon,by='taxonID') %>%
  mutate(taxonID = acceptedTaxonID) %>%
  select(-acceptedTaxonID) %>%
  rename(observation_id = uid,
         event_id = boutNumber,
         location_id = namedLocation,
         observation_datetime = endDate,
         taxon_id = taxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName) %>%
  mutate(package_id = NA) #%>%


# All taxonID and acceptedTaxonID values matched prior to removing acceptedTaxonID - is that accurate?
# 1) Create table_1m2observations through the left join, then test:
# which(table_1m2observations$taxonID %in% table_taxon$taxonID[table_taxon$taxonID != table_taxon$acceptedTaxonID]) #Returns 0

table_10_100m2observations <- div_10_100_m2 %>% 
  select(uid,
         namedLocation, 
         subplotID,
         boutNumber,
         endDate,
         taxonID,
         identificationQualifier
  ) %>%
  left_join(table_taxon,by='taxonID') %>%
  mutate(taxonID = acceptedTaxonID) %>%
  select(-acceptedTaxonID) %>%
  rename(observation_id = uid,
         event_id = boutNumber,
         location_id = namedLocation,
         observation_datetime = endDate,
         taxon_id = taxonID,
         taxon_rank = taxonRank,
         taxon_name = scientificName) %>%
  mutate(package_id = NA) #%>%

###################################
# write out in ecocomDP format
###
readr::write_csv(
  table_location_1m2,
  'table_location_1m2.csv')

readr::write_csv(
  table_location_10_100m2,
  'table_location_10_100m2.csv')

readr::write_csv(
  table_taxon,
  'table_taxon.csv')

readr::write_csv(
  table_1m2observations,
  'table_1m2observations.csv')

readr::write_csv(
  table_10_100m2observations,
  'table_10_100m2observations.csv')