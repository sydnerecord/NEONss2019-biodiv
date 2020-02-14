### This code downloads, unzips and restructures small mammal data (raw abundances) from NEON 
### for each NEON site, we provide mammal raw abundance aggregated per day (across all years), per month (bout) (across all years), and per year
### The final tidy data structure has locations and times as rows and species as columns, with raw abundance provided for each cell.

library(neonUtilities)
library(tidyverse)
library(lubridate)

if(!file.exists("data/mammal_basic.rds")){
  d = loadByProduct(dpID="DP1.10072.001", site = "all", package = "basic", avg = "all")
  saveRDS(d, file = "data/mammal_basic.rds")
} else {
  d = readRDS("data/mammal_basic.rds")
}

dat.mam = mutate(d$mam_pertrapnight, 
                 collectDate = ymd(collectDate),
                 year = year(collectDate),
                 month = month(collectDate),
                 day = day(collectDate),
                 # Since data collection usually takes place on several consecutive 
                 # days within a given month, we consider each month to be a separate bout 
                 bout = paste(year, month, sep = "_")) %>% 
  as_tibble()

### Here we provide the code that summarizes raw abundances per day, month (bout), and year
### We can remove NA (no captures) altogether at this point
dat.mam <- filter(dat.mam, scientificName != "", !is.na(scientificName))

### Remove species that are bycatch (non-target), dead, or escapted while processing
### remove all where the fate of the individual, unless marked and released, is
### 'dead' = dead, 'escaped' = escaped while handling, 'nontarget' = released, non-target species, 
### should 'released' (= target or opportunistic species released without full processing) be also removed?
dat.mam <- filter(dat.mam, !fate %in% c("dead", "escaped", "nontarget"))
#dat.mam <- filter(dat.mam, fate != "released") 

### We can also remove records no id'ed to species
# there are a couple of species that are not id'ed to species and yet denoted as "/" (either or)
# let's remove those first
dat.mam <- filter(dat.mam, grepl(pattern = "/|sp[.]|", scientificName))

### Remove recaptures -- Y and U (unknown); only retain N
dat.mam <- filter(dat.mam, recapture == "N")

### Get raw abundances per day
m1 <- dat.mam %>%
  select(siteID, year, month, day, scientificName) %>%
  group_by(siteID, year, month, day, scientificName) %>%
  tally(name = "count") %>% 
  ungroup()

m1_piv <- m1 %>%
  pivot_wider(id_cols = c(siteID, year, month, day), 
              names_from = scientificName, values_from = count,
              values_fill = list(count = 0)) %>%
  arrange(siteID, year, month, day) 

saveRDS(m1_piv, file="mammals_abs-abund-day.rds")

# Comments from DL: by month, and by year data can easily be derived from day data
# i.e.
group_by(m1, siteID, year, month, scientificName) %>% 
  summarise(count = sum(count, na.rm = TRUE)) %>% 
  ungroup()

