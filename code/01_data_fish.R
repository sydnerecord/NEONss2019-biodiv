
##############################################################################################
# Stephanie Parker, Thilina Surasinghe (sparker@battelleecology.org, tsurasinghe@bridgew.edu)
# last updated 2019-11-04
# uses esokol's 'get_and_format_NEON_inv_data.R'
##############################################################################################

# setwd("C:/Users/tsurasinghe/OneDrive - Bridgewater State University/Research2019/NEON")

#Clear all existing data
rm(list=ls())

#Close graphics devices
graphics.off()

# sk R to return memory to the operating system, might free up memory
gc()

# set options
options(stringsAsFactors = FALSE)


## getting all required packages

if(!require(installr)) {
  install.packages("installr"); require(installr)} #load / install+load installr
library(installr)

# need for rtools
install.packages("pkgbuild")

# install Rtools- can also install rtools manually. 
install.Rtools(choose_version = TRUE, check = TRUE, GUI = TRUE,
               page_with_download_url = "https://cran.r-project.org/bin/windows/Rtools/")

install.packages("devtools", type = "win.binary")
install.packages("neonUtilities")
install.packages("devtools")
install.packages("raster")
install.packages("backports")
install.packages("digest")
install.packages("xlsx")
devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
install.packages("BiocManager")
BiocManager::install("rhdf5")
# Install and load ecocomDP
install_github("EDIorg/ecocomDP")
install.packages("httr")
install.packages("jsonlite")
install.packages("downloader")
install.packages("dplyr")
install.packages("magrittr")

# # Packages used, load packages
library(neonUtilities)
library(tidyverse)
library(googledrive)
library(devtools)
library(ecocomDP)
# install_github("EDIorg/ecocomDP", ref = 'development')
library(lubridate)
library(iNEXT)
library(xlsx)
library(geoNEON)
library(pkgbuild)
library(rhdf5)
library(httr)
library(jsonlite)
library(dplyr, quietly=T)
library(downloader)
library(magrittr)

# update, without prompts for permission/clarification
update.packages(ask = FALSE)

#################################################################################
# fish dpid
my_dpid_fish <- 'DP1.20107.001'

# some aquatic sites-- for this need all sites
#my_site_list <- c('POSE', 'ARIK')

# get taxon table from API, may take a few minutes to load
# Fish electrofishing, gill netting, and fyke netting counts 
# http://data.neonscience.org/data-product-view?dpCode=DP1.20107.001
full_taxon_fish <- neonUtilities::getTaxonTable(taxonType = 'FISH', recordReturnLimit = NA, stream = "true") 

# -- make ordered taxon_rank list for a reference (subspecies is smallest rank, kingdom is largest)

# a much simple table with useful levels of taxonomic resolution
taxon_rank_fish <- c('superclass', 'class', 'subclass', 'infraclass', 'superorder',
                     'order', 'suborder', 'infraorder', 'section', 'subsection',
                     'superfamily', 'family', 'subfamily', 'tribe', 'subtribe', 'genus',
                     'subgenus','speciesGroup','species','subspecies') %>% rev() 


# get data FISH via api -- will take a while --  package = "basic" is also possible 
all_fish <- neonUtilities::loadByProduct(dpID = my_dpid_fish,  site = "all", startdate = NA, enddate = NA,
                                         package = "expanded", avg = "all", check.size = FALSE)

# The object returned by loadByProduct() is a named list of data frames. 
#To work with each of them, select them from the list using the $ operator.

# tibble::view() # to view tables in a seperate window

# join field data table with the pass tables from all_fish


fsh_dat1 <- dplyr::left_join(x = all_fish$fsh_perPass, y = all_fish$fsh_fieldData, by = c('reachID')) %>% 
 dplyr::filter(is.na(samplingImpractical) | samplingImpractical == "") #remove records where fish couldn't be collected, both na's and blank data

# get rid of dupe col names and .x suffix
fsh_dat1 <- fsh_dat1[,!grepl('\\.y',names(fsh_dat1))]
names(fsh_dat1) <- gsub('\\.x','',names(fsh_dat1))

# add individual fish counts
fsh_dat_indiv <- dplyr::left_join(all_fish$fsh_perFish, fsh_dat1, by = "eventID") 

# get rid of dupe col names and .x suffix
fsh_dat_indiv <- fsh_dat_indiv[,!grepl('\\.y',names(fsh_dat_indiv))]
names(fsh_dat_indiv) <- gsub('\\.x','',names(fsh_dat_indiv))

# fill in missing reachID with event ID
fsh_dat_indiv$reachID <- ifelse(is.na(fsh_dat_indiv$reachID), 
                                substr(fsh_dat_indiv$eventID, 1, 16), fsh_dat_indiv$reachID) 


# add bulk fish counts
fsh_dat_bulk <- dplyr::left_join(all_fish$fsh_bulkCount, fsh_dat1, by = "eventID")

# get rid of dupe col names and .x suffix
fsh_dat_bulk <- fsh_dat_bulk[,!grepl('\\.y',names(fsh_dat_bulk))]
names(fsh_dat_bulk) <- gsub('\\.x','',names(fsh_dat_bulk))

#fill in missing reachID
fsh_dat_bulk$reachID <- ifelse(is.na(fsh_dat_bulk$reachID), 
                               substr(fsh_dat_bulk$eventID, 1, 16), fsh_dat_bulk$reachID) 


# combine indiv and bulk counts
fsh_dat <-dplyr::bind_rows(fsh_dat_indiv, fsh_dat_bulk)

# add count = 1 for indiv data
fsh_dat$count <- ifelse(is.na(fsh_dat$bulkFishCount), 1, fsh_dat$bulkFishCount)

# need to coovert POSIXt format into as.character and then bacl to date-time format
fsh_dat$startDate <- dplyr::if_else(is.na(fsh_dat$startDate),
                  lubridate::as_datetime(substr(as.character(fsh_dat$passStartTime), 1, 10)), fsh_dat$startDate) # 1-10: number of characters on date
fsh_dat$siteID <- dplyr::if_else(is.na(fsh_dat$siteID), 
                         substr(fsh_dat$eventID, 1, 4), fsh_dat$siteID) # four characters on site

# get species and finer res
fsh_dat$taxonRank_ordered <- factor(
  fsh_dat$taxonRank,
  levels = taxon_rank_fish,
  ordered = TRUE) 

# get all records that have rank <= species
fsh_dat_fine <- fsh_dat %>%
  dplyr::filter(taxonRank_ordered <= 'species')


# grouping vars for aggregating density measurements 
my_grouping_vars <- c('domainID','siteID','aquaticSiteType','namedLocation',
                      'startDate','reachID','eventID','samplerType', 'aquaticSiteType', 'netSetTime', 'netEndTime',
                      'netDeploymentTime', 'netLength', 'netDepth', 'fixedRandomReach','measuredReachLength','efTime', 
                      'efTime2', "passStartTime", "passEndTime", 'netDeploymentTime', 'scientificName', 'passNumber') 
# added a few metrics to quantify catch per unit effort such as 'passNumber', efish time, net deployment time 

# aggregate densities for each species group, pull out year and month from StartDate

require(dplyr)

fsh_dat_aggregate <- fsh_dat_fine %>%
  dplyr::select(!!c(my_grouping_vars, 'count')) %>%
  dplyr::group_by_at(dplyr::vars(my_grouping_vars)) %>%  
  dplyr::summarize(
    number_of_fish = sum(count),
    n_obs = dplyr::n()) %>%
    dplyr::mutate(
      year = startDate %>% lubridate::year(),
      month = startDate %>% lubridate::month()
  ) %>% dplyr::ungroup()

# some aquaticSiteType are NA, replace NAs if-based wildcarding namedLOcation 
# grepl is for wildcarding 
fsh_dat_aggregate$aquaticSiteType <- dplyr::if_else(is.na(fsh_dat_aggregate$aquaticSiteType),
                                             if_else(grepl("fish.point", fsh_dat_aggregate$namedLocation), 'stream','lake'), fsh_dat_aggregate$aquaticSiteType)

# sampler type is also missing in a few cases, reaplace from eventID, with a wildcard
fsh_dat_aggregate$samplerType <- dplyr::if_else(is.na(fsh_dat_aggregate$samplerType), 
                                    dplyr::if_else(grepl("e-fisher", fsh_dat_aggregate$eventID), 'electrofisher', 
                                       dplyr::if_else(grepl("gill", fsh_dat_aggregate$eventID), 'gill net', 'mini-fyke net')), fsh_dat_aggregate$samplerType) 
                                              

# make sure that e-fish samples have "" or NA for netset and netend times, this will leave actual missing data as na
# change netset/netend times to a datetime format if needed: lubridate::as_datetime(fsh_xx$netSetTime, format="%Y-%m-%d T %H:%M", tz="GMT")
## then calculate the netdeploymenttime (in hours) from netset and netend time
# to calculate pass duration (in mins) and also caculate average efish time (secs); also if efishtime is zero, --> na's 

fsh_aggregate_mod <- fsh_dat_aggregate %>% dplyr::mutate(netSetTime = dplyr::if_else(condition = samplerType ==  "electrofisher" | samplerType == "two electrofishers",
                    true = lubridate::as_datetime(NA), false = netSetTime), 
                    netEndTime = dplyr::if_else(condition = samplerType ==  "electrofisher" | samplerType == "two electrofishers",
                    true = lubridate::as_datetime(NA), false = netEndTime),
                    netDeploymentTime = dplyr::case_when(samplerType == grepl("electrofish", samplerType) ~ netDeploymentTime, 
                    is.na(netDeploymentTime) ~ as.numeric(difftime(netEndTime, netSetTime, tz = "GMT", units = "hours")), 
                    TRUE ~ netDeploymentTime),
                    mean_efishtime = base::rowMeans(dplyr::select(., c("efTime", "efTime2")), na.rm = T),
                    mean_efishtime = dplyr::case_when(mean_efishtime == 0 ~ NA_real_,TRUE ~ mean_efishtime))

# with the above changes, we have efish time and and net duration to calculate catch per unit effort before moving to wide format
# CUPE with efish time, calculated by = (total number of fish/average e-fish time insecs * 3600) as fish captured per 1-hr of e-fishing
# CUPE with gill nets and fykenets calculated by = (total number of fish/netDeployment time in hours * 24) as fish captured per 1-day-mong net deployment of e-fishing
fsh_aggregate_mod2 <- fsh_aggregate_mod %>% dplyr::mutate(CUPE = dplyr::if_else(condition = samplerType ==  "electrofisher" | samplerType == "two electrofishers",
                            true = number_of_fish/mean_efishtime * 3600, 
                            false = dplyr::if_else(condition = samplerType ==  "mini-fyke net" | samplerType == "mini-fyke net",
                            true = number_of_fish/netDeploymentTime * 24, false = as.numeric(NA))))

# make wide for total observations without catch per unit efforts
fsh_dat_wide_total.obs <- fsh_aggregate_mod %>% 
  dplyr::group_by(year, month, siteID, namedLocation, reachID, fixedRandomReach, aquaticSiteType, samplerType) %>%
  tidyr::spread(key = scientificName, value = number_of_fish, fill = 0, convert = FALSE, drop = TRUE, sep = NULL) 

# make wide for catch per unit efforts
fsh_dat_wide_CUPE <- fsh_aggregate_mod2 %>% 
  dplyr::group_by(year, month, siteID, namedLocation, reachID, fixedRandomReach, aquaticSiteType, samplerType) %>%
  tidyr::spread(key = scientificName, value = CUPE, fill = 0, convert = FALSE, drop = TRUE, sep = NULL) 


## writing the data into csv and txt formats

install.packages("readr")
library(readr)

save(fsh_aggregate_mod, file = "fsh_aggregate_mod.RData")

save(fsh_dat_wide_total.obs, file = "fsh_dat_wide_totobs.RData")

save(fsh_dat_wide_CUPE, file = "fsh_dat_wide_CUPE.RData")

save(fsh_aggregate_mod2, file = "fsh_dat_aggregate_mod2.RData")

# ## writing the data into csv and txt formats

readr::write_csv(fsh_dat_wide_CUPE, "fsh_data_wide_CUPE.csv")

readr::write_csv(fsh_dat_wide_total.obs, "fsh_data_wide_total_obs.csv")

readr::write_csv(fsh_aggregate_mod2, "fsh_data_long.csv")


