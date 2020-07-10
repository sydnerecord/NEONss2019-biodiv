# Download and clean tick pathogen data 
# Authors:  Melissa Chen, Wynne Moss, Brendan Hobart, Matt Bitters
# Date: 7/1/20  (Happy Canada Day!)

## Note: I do NOT cross-reference subsampleIDs or testing numbers with the tick field data. ##
## OUTPUT: filtered data file; log file (optional; default FALSE)

library(tidyverse)
library(neonUtilities)
library(lubridate)

#### Load files ####
keepLOG <- TRUE

if(!file.exists("./data/tick_pathogen.csv")){
  tickpathdat <- loadByProduct(dpID="DP1.10092.001", site="all", package = "expanded", check.size = F)
  
  tick_path <- tickpathdat$tck_pathogen
  tick_pathqa <- tickpathdat$NEON.University_of_Massachusetts_Laboratory_of_Medical_Zoology.tck_pathogenqa
}


# tick_path <- read.csv("data_raw/filesToStack10092/stackedFiles/tck_pathogen.csv", na.strings = c("","NA","na","n/a","N/A"))
# tick_pathqa <- read.csv("data_raw/filesToStack10092/stackedFiles/NEON.University_of_Massachusetts_Laboratory_of_Medical_Zoology.tck_pathogenqa.csv", na.strings = c("","NA","na","n/a","N/A"))

keepLOG=FALSE # Change if you want to record log of tick filtering 
# OPTIONAL: Set up sink to capture data filtering
if ( keepLOG ) {
  sink(file="data/LOG_tickpath.txt")
}

#### Raw dataset inspection ####
## Record original number of tests
print(paste0("## Original number of rows: "))
nrow(tick_path)

## Record original number of test results that are not NA
print(paste0("## Original number of rows (NA test results removed): "))
nrow(tick_path %>% filter(!is.na(testResult)))

#### Quality check file ####
## Remove any samples that have flagged quality checks
remove_qa_uid_batchID <- tick_pathqa[tick_pathqa$criteriaMet != "Y",c("uid","batchID")]
print(paste0("## Samples that did not pass quality checks: "))
remove_qa_uid_batchID

## Record which batchIDs in tick_path are not in quality file
print(paste0("## Batches NOT found in quality file: "))
unique(tick_path[!(tick_path$batchID %in% tick_pathqa$batchID),"batchID"])
# Note: in the qa file it says NEON_2017831, but in the tick_path file it's NEON_20170831. I think these are the same thing.

#### Tick Pathogen Quality assessment ####
## Checking uid is unique
Is_all_ID_uniqud <- nrow(tick_path) == length(unique(tick_path$uid))
if ( !Is_all_ID_uniqud ) {
  warningCondition(message = "uid is not unique. Check file manually.")
}
print(paste0("## uids are unique?: ", Is_all_ID_uniqud))


## Check there isn't any important missing data 
# List of important variables
checkVar <- c("domainID","siteID","plotID", "decimalLatitude", "decimalLongitude"
              , "plotType", "nlcdClass", "elevation", "collectDate"
              , "subsampleID", "batchID", "testingID", "testPathogenName")
# Check which variables are missing
which_missing_var <- sapply(checkVar, FUN=function(v){any(is.na(tick_path[,v]))} )
print(paste0("## Which variables are missing data?: "))
which_missing_var
# Identify samples with missing values
removeSamples_missingvalue <- tick_path[is.na(tick_path[,names(which_missing_var)[which_missing_var]]),c("uid","namedLocation","collectDate","sampleCondition","remarks")]
# Will need to filter these out.
removeSamples_missingvalue

## Look at sample conditions; filter out all non-OKAY
print(paste0("## Sample conditions: "))
table(tick_path$sampleCondition) # We'll get rid of all "non-okay" samples
removeSamples_sampleCondition <- tick_path[tick_path$sampleCondition!="OK",c("uid","namedLocation","collectDate","sampleCondition","remarks")]
removeSamples_sampleCondition

## How many test results are NA?
print(paste0("## Number of test results that are NA: "))
sum(is.na(tick_path$testResult)) # a lot. I think some are 2019 data, and others are just unknown/untested
removeSamples_testResult <- tick_path[is.na(tick_path$testResult),]

## Check DNA quality of each testing batch
print(paste0("## HardTick DNA Quality: "))
tick_path %>% filter(testPathogenName=="HardTick DNA Quality") %>% select(testResult,testPathogenName) %>% table(useNA="ifany")
# Only positive DNA test results-- filter out others
removeSamples_DNAQual <- tick_path %>% filter(testPathogenName=="HardTick DNA Quality", testResult!="Positive") %>% select(uid, namedLocation, collectDate, sampleCondition, remarks)
removeSamples_DNAQual

## Check number of ticks tested per site per year
print(paste0("## Site x year for all tested ticks: "))
tick_path %>% filter(!is.na(testResult), !(testPathogenName %in% c("HardTick DNA Quality","Ixodes pacificus"))) %>%mutate(year = year(collectDate)) %>% 
  pivot_wider(id_cols=c(siteID, year, testingID), names_from=testPathogenName, values_from=testResult)%>% 
  select(siteID, year) %>% table()
# It is NOT true that only a subset of 130 ticks were tested from each site and year... do they mean PLOT?
# or maybe they mean "by species"?
# By species-- try the two main ones
print(paste0("## Site x year for Ambame only: "))
tick_path %>%
  filter(!is.na(testResult), !(testPathogenName %in% c("HardTick DNA Quality","Ixodes pacificus"))) %>% 
  separate(subsampleID, into=c("plotID","date","species","lifestage"), sep="[.]", remove=FALSE) %>%
  mutate(year = year(collectDate)) %>%
  pivot_wider(id_cols=c(siteID, year, species, testingID), names_from=testPathogenName, values_from=testResult)%>% 
  filter(species=="AMBAME")  %>% select(siteID, year) %>% table()

print(paste0("## Site x year for Ixosca only: "))
tick_path %>%
  filter(!is.na(testResult), !(testPathogenName %in% c("HardTick DNA Quality","Ixodes pacificus"))) %>% 
  separate(subsampleID, into=c("plotID","date","species","lifestage"), sep="[.]", remove=FALSE) %>%
  mutate(year = year(collectDate)) %>%
  pivot_wider(id_cols=c(siteID, year, species, testingID), names_from=testPathogenName, values_from=testResult)%>% 
  filter(species=="IXOSCA")  %>% select(siteID, year) %>% table()

# By plot
print(paste0("## Plot x year: "))
tick_path %>%
  filter(!is.na(testResult), !(testPathogenName %in% c("HardTick DNA Quality","Ixodes pacificus"))) %>%
  separate(subsampleID, into=c("plotID","date","species","lifestage"), sep="[.]", remove=FALSE) %>%
  mutate(year = year(collectDate)) %>%
  pivot_wider(id_cols=c(plotID, year, species, testingID), names_from=testPathogenName, values_from=testResult)%>% 
  select(plotID, year) %>% table()


### Which species were tested for what pathogens?
# Create table
print(paste0("## Which species was tested for which pathogens: "))
tick_path %>% filter(!is.na(testResult)) %>%separate(subsampleID, into=c("plotID","date","species","lifestage"), sep="[.]", remove=FALSE) %>%
  select(testPathogenName, species) %>% table()
## Get rid of hardtick DNA and Ixodes pacificus tests
remove_testPathogenName <- tick_path %>% filter(testPathogenName %in% c("HardTick DNA Quality","Ixodes pacificus"))

if (keepLOG) {
  sink()
}

# Note that this does NOT follow their description from online, which says certain species were tested for certain pathogens.
# In reality, it seems almost all individuals were tested for all pathogens (with a few exceptions)

#### Filtering ####
tick_path_filt <- tick_path %>% filter( !(uid %in% c(remove_qa_uid_batchID$uid
                                 , removeSamples_missingvalue$uid
                                 , removeSamples_sampleCondition$uid
                                 , removeSamples_testResult$uid
                                 , removeSamples_DNAQual$uid
                                 ,remove_testPathogenName$uid)))
# Re-formatting to species x sample
tick_path_samplexspecies <- tick_path_filt %>% select(uid, testPathogenName, testResult) %>%
  mutate(testResult=ifelse(testResult=="Positive",1,0)) %>%
  pivot_wider(id_cols=uid, names_from=testPathogenName, values_from=testResult)

# Re-formatting to sample x env
tick_path_samplexenv <- tick_path_filt %>% 
  rename(latitude = decimalLatitude,
         longidtude = decimalLongitude) %>%  select(-c(dataQF, remarks, individualCount, testPathogenName, sampleCondition, testResult))

write.csv(tick_path_filt, "data/tick_pathogen.csv", row.names = FALSE, quote = FALSE)

# saveRDS(tick_path_filt, "data/tick_path_filt_long.RData") # Optional output
# saveRDS(tick_path_samplexspecies, "data/tick_path_filt_samplexspecies.RData") # Optional output
# saveRDS(tick_path_samplexenv, "data/tick_path_filt_samplexenv.RData") # Optional output



