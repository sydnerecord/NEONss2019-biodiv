#title: download mosquito data and clean
#author: Natalie Robinson
#date: 6/08/2020

#packages used
library(dplyr)
library(tidyr)
library(neonUtilities)
options(stringsAsFactors = F)

#get mosquito data -------------
# mos dpid
mos_code <-'DP1.10043.001'

#only two sites - 
# site_ids = c("NIWO","DSNY")-EX
mos_allTabs <- loadByProduct(dpID = mos_code, 
                             #site = site_ids,
                             startdate = "2016-1", enddate = "2018-11",  
                             package = "expanded", check.size = TRUE)


# Identify names of desired tables and pull data into corresponding data frames
tables <- names(mos_allTabs)[grepl('mos',names(mos_allTabs))]

for (t in tables){
  assign(t,mos_allTabs[[which(names(mos_allTabs)==t)]])
}

#################################################################################
# Clean data

# Remove rows with missing critical data ########################################
# Note - a missing mos_trapping:eventID was usually a collection cup with no target taxa present.  
#          In 23 cases, there was no eventID and targetTaxaPresent was 'U' (unknown)

# Define important data colums
cols_oi_trap <- c('collectDate','eventID','namedLocation','sampleID')
cols_oi_sort <- c('collectDate','sampleID','namedLocation','subsampleID')
cols_oi_arch <- c('startCollectDate','archiveID','namedLocation')
cols_oi_exp_proc <- c('collectDate','subsampleID','namedLocation')

# Remove rows with missing important data, replacing blank with NA and factors with characters
mos_trapping <- mutate_all(mos_trapping, list(~na_if(.,""))) %>% 
  drop_na(all_of(cols_oi_trap)) %>%
  mutate_if(is.factor, as.character)
mos_sorting <- mutate_all(mos_sorting, list(~na_if(.,""))) %>% 
  drop_na(all_of(cols_oi_sort)) %>%
  mutate_if(is.factor, as.character)
mos_archivepooling <- mutate_all(mos_archivepooling, list(~na_if(.,""))) %>% 
  drop_na(all_of(cols_oi_arch)) %>%
  mutate_if(is.factor, as.character)
mos_expertTaxonomistIDProcessed <- mutate_all(mos_expertTaxonomistIDProcessed, list(~na_if(.,""))) %>% 
  drop_na(all_of(cols_oi_exp_proc)) %>%
  mutate_if(is.factor, as.character)


# Look for duplicate sampleIDs #################################################
length(which(duplicated(mos_trapping$sampleID)))  # 0
length(which(duplicated(mos_sorting$subsampleID))) # 0
length(which(duplicated(mos_archivepooling$archiveID))) # 0

# Make sure no "upstream" records are missing primary "downstream" reference
length(which(!mos_sorting$sampleID %in% mos_trapping$sampleID))  # 0
length(which(!mos_expertTaxonomistIDProcessed$subsampleID %in% mos_sorting$subsampleID))  # 0

# Clear expertTaxonomist:individualCount if it is 0 and there is no taxonID. These aren't ID'ed samples
mos_expertTaxonomistIDProcessed$individualCount[mos_expertTaxonomistIDProcessed$individualCount == 0 & is.na(mos_expertTaxonomistIDProcessed$taxonID)] <- NA


#################################################################################
# Join data

# Add trapping info to sorting table -------------------------------------------- 
# Note - 59 trapping records have no associated sorting record and don't come through the left_join (even though targetTaxaPresent was set to Y or U)
mos_dat <- mos_sorting %>%
  select(-c(collectDate, domainID, namedLocation, plotID, setDate, siteID)) %>%
  left_join(select(mos_trapping,-uid),by = 'sampleID') %>%
  rename(sampCondition_sorting = sampleCondition.x,
         sampCondition_trapping = sampleCondition.y,
         remarks_sorting = remarks,
         dataQF_trapping = dataQF)

# Verify sample barcode consistency before removing one column
which(mos_dat$sampleCode.x != mos_dat$sampleCode.y)

# Consistency is fine with barcodes
mos_dat <- rename(mos_dat,sampleCode = sampleCode.x) %>% 
  select(-sampleCode.y)
  

# Join expert ID data ------------------------------------------------------------
mos_dat <- mos_dat %>%
  left_join(select(mos_expertTaxonomistIDProcessed,
                   -c(uid,collectDate,domainID,namedLocation,plotID,setDate,siteID,targetTaxaPresent)),
            by='subsampleID') %>%
  rename(remarks_expertID = remarks)

# Verify sample barcode and labName consistency before removing one column
which(mos_dat$subsampleCode.x != mos_dat$subsampleCode.y)
which(mos_dat$laboratoryName.x != mos_dat$laboratoryName.y)

# Rename columns and add estimated total individuals for each subsample/species/sex with identification, where applicable
#  Estimated total individuals = # individuals iD'ed * (total subsample weight/ subsample weight)
mos_dat <- rename(mos_dat,subsampleCode = subsampleCode.x, laboratoryName = laboratoryName.x) %>% 
  select(-c(subsampleCode.y, laboratoryName.y)) %>%
  mutate(estimated_totIndividuals = ifelse(!is.na(individualCount),round(individualCount * (totalWeight/subsampleWeight)), NA))


# Add archive data ---------------------------------------------------------------- 
mos_dat <- mos_dat %>%
  left_join(select(mos_archivepooling,-c(domainID,uid,namedLocation,siteID)),by = 'archiveID')

# Verify sample barcode consistency before removing one column
which(mos_dat$archiveIDCode.x != mos_dat$archiveIDCode.y)

# Consistency is fine with barcodes
mos_dat <- rename(mos_dat,archiveIDCode = archiveIDCode.x) %>% 
  select(-archiveIDCode.y)


# Write data if desired
# write.csv(mos_dat,'full_mos_dataset.csv',row.names=F,na='')
