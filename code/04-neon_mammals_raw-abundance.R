### This code downloads, unzips and restructures small mammal data (raw abundances) from NEON 
### for each NEON site, we provide mammal raw abundance aggregated per day (across all years), per month (bout) (across all years), and per year
### The final tidy data structure has locations and times as rows and species as columns, with raw abundance provided for each cell.

###########################################################################
###  Load functions
###########################################################################
#install.packages("neonUtilities")
#devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
#install.packages("BiocManager")
#BiocManager::install("rhdf5")
require(neonUtilities)
require(geoNEON)
require(raster)
require(rhdf5)
require(here)
require(tidyverse)
require(lubridate)
select <- dplyr::select


###########################################################################
## Download and unzip mammal data from NEON
###########################################################################
zipsByProduct(dpID="DP1.10072.001", site = "all", startdate = NA, enddate = NA,
             package = "basic", avg = "all", check.size = TRUE, savepath = NA,
             load = F)
stackByTable("/filesToStack10072/", folder=T)


###########################################################################
## Read in data
###########################################################################
dat.mam <- read.delim(file="/filesToStack10072/stackedFiles/mam_pertrapnight.csv", sep=",")
head(dat.mam)

### Create columns for day, month, and year of capture
t1 <- t(as.data.frame(str_split(dat.mam$collectDate, "-")))
dat.mam$year <- as.numeric(t1[,1])
dat.mam$month <- as.numeric(t1[,2])
dat.mam$day <- as.numeric(t1[,3])

### Since data collection usually takes place on several consecutive days within a given month, we consider each month to be a separate bout 
dat.mam$bout <- apply(dat.mam[,62:63], 1, paste, collapse="_") #there are 60 bouts (unique across all years), but the number will differ for each year
bouts <- unique(dat.mam.c1$bout)

### Here we provide the code that summarizes raw abundances per day, month (bout), and year
### We can remove NA (no captures) altogether at this point
dat.mam <- filter(dat.mam.c1, scientificName != "")
spnames <- unique(dat.mam$scientificName) 
spnames <- as.data.frame(spnames)

### Remove species that are bycatch (non-target), dead, or escapted while processing
### remove all where the fate of the individual, unless marked and released, is
### 'dead' = dead, 'escaped' = escaped while handling, 'nontarget' = released, non-target species, 
### should 'released' (= target or opportunistic species released without full processing) be also removed?
dat.mam <- filter(dat.mam, fate != "dead")
dat.mam <- filter(dat.mam, fate != "escaped")
dat.mam <- filter(dat.mam, fate != "nontarget")
#dat.mam <- filter(dat.mam, fate != "released") 

### We can also remove records no id'ed to species
# there are a couple of species that are not id'ed to species and yet denoted as "/" (either or)
# let's remove those first
spnames.sep <- spnames  %>%
  separate(col=spnames, into = c("sciname", "potsp"), sep="/", remove = FALSE) %>%
  distinct()
spnames.sep <- spnames.sep %>%
  filter(is.na(potsp) == TRUE)

# let's remove those denoted as sp.
spnames.sep <- spnames.sep  %>%
  separate(col=sciname, into = c("genus", "species", "subspecies"), sep=" ", remove = FALSE) %>%
  select(spnames,genus,species,subspecies)

spnames.sep <- spnames.sep %>%
  filter(species != "sp.")
  
dat.mam <- dat.mam %>%
  filter(scientificName %in% spnames.sep$spnames)

### Remove recaptures -- Y and U (unknown); only retain N
dat.mam <- dat.mam %>%
  filter(recapture == "N")


### Get raw abundances per day
m1 <-  
  dat.mam %>%
  select(siteID, collectDate, year, month, day, collectDate, scientificName) %>%
  group_by(siteID, year, month, day, collectDate, scientificName) %>%
  summarise(
    count=n()
  ) 

m1_piv <- m1  %>%
  pivot_wider(id_cols=c(siteID, year, month, day, collectDate), names_from = scientificName, values_from = count) %>%
  arrange(siteID, year, month, day) 


### Get raw abundances per month (bout)
m2 <-  
  dat.mam %>%
  select(siteID, collectDate, year, month, bout, scientificName) %>%
  group_by(siteID, year, month, bout, scientificName) %>%
  summarise(
    count=n()
  ) 

m2_piv <- m2  %>%
  pivot_wider(id_cols=c(siteID, year, month, bout), names_from = scientificName, values_from = count) %>%
  arrange(siteID, year, month) 


### Get raw abundances per year
m3 <-  
  dat.mam %>%
  select(siteID, collectDate, year, scientificName) %>%
  group_by(siteID, year, scientificName) %>%
  summarise(
    count=n()
  ) 

m3_piv <- m3  %>%
  pivot_wider(id_cols=c(siteID, year), names_from = scientificName, values_from = count) %>%
  arrange(siteID, year) 


### Replace NA with 0--indication of species not found (captured) at the site
m1_piv[is.na(m1_piv)] <- 0
m2_piv[is.na(m2_piv)] <- 0
m3_piv[is.na(m3_piv)] <- 0

saveRDS(m1_piv, file="mammals_abs-abund-day.rds")
saveRDS(m2_piv, file="mammals_abs-abund-month.rds")
saveRDS(m3_piv, file="mammals_abs-abund-year.rds")
