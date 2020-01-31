### This code downloads, unzips and restructures small mammal data from NEON for 
### the use in spatial capture-recapture (SCR) models SCR models call for a number 
### of decisions (period of closure, startification, etc.) contingnet upon the 
### specific research question. In this scenario, we prep data for a stratified 
### SCR (we stratify by NEON site); we assume closure within a month (bout) and 
### that population is open across months. Each year will be fitted separately 
### Other scenarios are possible (e.g., population closed within a year, open 
### across years, etc.)

###########################################################################
###  Load functions
###########################################################################
install.packages("neonUtilities")
devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
install.packages("BiocManager")
BiocManager::install("rhdf5")
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
zipsByProduct(dpID = "DP1.10072.001", site = "all", startdate = NA, enddate = NA,
              package = "basic", avg = "all", check.size = TRUE, savepath = NA,
              load = F)
stackByTable("filesToStack10072/", folder=T)
# just loadByProduct()?

###########################################################################
## Read in data
###########################################################################
dat.mam <- read_csv(file="filesToStack10072/stackedFiles/mam_pertrapnight.csv")
head(dat.mam)

### 
dat.mam = mutate(dat.mam, year = lubridate::year(collectDate),
                 month = lubridate::month(collectDate),
                 day = lubridate::day(collectDate))

###########################################################################
## Use collectDate to create columns that designate a bout (here, month) and repeated visit within a bout for Scenario 1
## Scenario 1: Stratified spatial capture-recapture model, for each year separately, closed within a bout (month), open across months
## Repeated visits are all days that fall within a bout (month)
###########################################################################
dat.mam.c1 <- dat.mam
dat.mam.c1$bout <- apply(dat.mam.c1[,62:63], 1, paste, collapse="_") 
bouts <- unique(dat.mam.c1$bout)

### Loop over
dat.mam.c1$rep <- 0

for (k in 1:length(bouts)){
b1 <- dat.mam.c1 %>% filter(bout == bouts[k])
sts <- unique(b1$siteID)
  for (j in 1:length(sts)){
  b1sts <- b1 %>% filter(siteID == sts[j])
  reps <- cbind(unique(b1sts$day),seq(1,length(unique(b1sts$day)),1))
    for (i in 1:nrow(reps)){
    dat.mam.c1$rep[dat.mam.c1$bout == bouts[k] & dat.mam.c1$siteID == sts[j] & dat.mam.c1$day == reps[i,1]] <- reps[i,2]
}}}

dat.mam.c1 <- dat.mam.c1 %>%
  unite("bout_rep", bout:rep, remove=FALSE)

saveRDS(dat.mam.c1, file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/filesToStack10072/stackedFiles/mam_pertrapnight_collectDate_Bouts_C1.rds")
dat.mam.c1 <- readRDS(file="/Users/jarzyna.1/Documents/RESEARCH_NEON/Data/filesToStack10072/stackedFiles/mam_pertrapnight_collectDate_Bouts_C1.rds")

### Organize the data to call on later
datorg <- dat.mam.c1 %>%
  select(endDate, year, month, bout, rep, bout_rep) %>%
  unique() %>%
  arrange(endDate, rep)

###########################################################################
## Do for each year separately to cut down on storage/computational costs and because SCR models will be run separately for each year
###########################################################################
years <- unique(dat.mam.c1$year)
years <- years[order(years)]
  
for (j in 1: length(years)){
dat.mam.c1y <- dat.mam.c1 %>% 
  filter(year == years[j])

year_j <- years[j]

###########################################################################
## Cut to only the relevant parts
###########################################################################
### Cut to a few columns to make more manageable
yup <- dat.mam.c1y %>%
	select(uid, siteID, plotID, trapCoordinate, trapStatus, bout_rep, tagID) %>%
  arrange(siteID, plotID, trapCoordinate, bout_rep)

### Add a trap ID column that combines Site, plot, trap in a single term
yup <- yup %>%
	mutate(trap_id = paste0(plotID,"_",trapCoordinate))

###########################################################################
## Extract all the unique dates and sites
###########################################################################
### get all possible dates
all_dates <- sort(unique(yup$bout_rep))
### get all possible sites
sites_unique <- unique(yup$siteID)

###########################################################################
## Checks on trap status assumptions
###########################################################################

### Be careful, there are rows with trap not set
sort(unique(dat.mam$trapStatus))
### There are also cases where the trapStatus is blank
dat.mam %>%
	filter(trapStatus == "")

### Assumptions we make are as follows:
### We assume the trap is not set if status is "1-trap not set"
### We assume trap is not set if trap status is blank
### We assume trap is not set if there is no data.
### We assume trap was set for option 2 (malfunction)
### We assume trap was set for options 3, 4, 5, 6

###########################################################################
## Extract at a site level
###########################################################################
for (i in seq(1, length(sites_unique))) {
site_i <- sites_unique[i]
cat(paste0("Loop ", i, "  Site  ", site_i, "\n"))

### Cut to only one site at a time
### Then get rid of site, plot, trap, leaving just the combined id
### Run distinct command because there are some duplicates
site_subset <- yup %>%
	filter(siteID == site_i) %>%
	select(trap_id, bout_rep, trapStatus, tagID) %>%
	distinct() #Retain only unique/distinct rows from an input tbl.

### Check results
head(site_subset)

###########################################################################
## Extract all the unique traps and tags
###########################################################################
### get all possible traps
trap_unique <- unique(site_subset$trap_id)

### Get all possible tags
tag_unique <- unique(site_subset$tagID)
tag_unique <- tag_unique[tag_unique != ""]

###########################################################################
## Process for trap being out
###########################################################################
### Expand out to all dates and group by trap id
trap_subset <- site_subset %>%
	select(trap_id, bout_rep, trapStatus) %>%
	group_by(trap_id) %>% 
	complete(trap_id, bout_rep = all_dates) #this ensures that for each trap we have all possible bouts_rep

### Create a column for trap set
### We set it to FALSE if "1-trap not set" or if there is an NA induced by not having a row for that day
### We say the trap was set even for option 2 (malfunction)
trap_subset <- trap_subset %>%
	mutate(trap_set = case_when( is.na(trapStatus) ~ 0,
		trapStatus == "1 - trap not set" ~ 0,
		trapStatus == "" ~ 0,
		TRUE ~ 1 
	)) 

### Check the results
head(trap_subset)
dim(trap_subset)

### There might be some duplicates (particularly if more than 1 capture in trap), so let's get rid of them
trap_subset <- trap_subset %>%
	distinct()

### Check the results
head(trap_subset)
dim(trap_subset)

### We now have a dataframe with each trap on each day and whether it was set
trap_subset <- trap_subset %>%
	select(-trapStatus)

###########################################################################
## Process for each tag
###########################################################################
for (k in seq(1, length(tag_unique))) {

tag_k <- tag_unique[k]
cat(paste0(tag_unique[k], "\n"))
cat(paste0("k = ", k, "\n"))
cat(paste0(format(Sys.time(), "%a %b %d %X %Y"), "\n"))

### Separate out the trapstatus column into number and text
capture_subset <- site_subset  %>%
	filter(tagID == tag_k) %>%	
	separate(col=trapStatus, into = c("trapstatus_num", "trapstatus_text"), sep="-", remove = FALSE) %>%
	mutate(trapstatus_num = as.numeric(trapstatus_num)) %>%
	distinct()

### Create a capture column which is 1 for status 4 (>1 capture in the trap), 5 (capture), or wherever there is a tagID
capture_subset <- capture_subset %>%
	mutate(capture = case_when( trapstatus_num == 4 | trapstatus_num == 5 ~ 1,
		tagID != "" ~ 1,
		TRUE ~ 0
	)) 

### Join with trap information
results_i <- trap_subset %>%
	left_join(capture_subset, by = c("trap_id", "bout_rep"))%>%
	mutate(tagID = tag_k) #this bascially joins capture_subset (the one tag) with trap_subset, showing where and when was this individual caught; if not caught==NA

### Create a new column for presence, which is NA if trap not set, 0 if trap set but no capture and 1 for capture
### Add a -9999 for errors
results_i <- results_i %>%
	mutate(presence = case_when( trap_set == 1 & capture == 1 ~ 1,
		trap_set == 1 & capture == 0 ~ 0,
		trap_set == 1 & is.na(capture) ~ 0,
		trap_set == 0 ~ NA_real_,	
		TRUE ~ -9999 
	)) 

#results_i %>% filter(trap_set == 1) %>% head()
#results_i %>% filter(presence == -9999) %>% head()

### Some initial reorganization
results_i_wide <- results_i %>%
	select(trap_id, bout_rep, tagID, presence)  %>%
	#arrange(bout_rep) %>%
	distinct() %>%
  mutate(bout_rep = factor(bout_rep, levels=unique(datorg$bout_rep))) %>%
  arrange(bout_rep)

### Convert this long table into a wide table and sort based on site, plot, and trap
results_i_wide <- results_i_wide  %>%
	pivot_wider(id_cols=c(tagID, trap_id), names_from = bout_rep, values_from = presence) %>%
	arrange(tagID, trap_id) 

head(results_i_wide)
dim(results_i_wide)

### Separate out the trap id column
results_i_wide <- results_i_wide %>%
	separate(col=trap_id, into = c("siteID", "plotID", "trapCoordinate"), sep="_", remove = FALSE)

#results_i_wide[1:10,1:8]

### After first loop over tagID, add to bottom
if (k == 1) {
	site_wide <- results_i_wide
} else {
	site_wide <- site_wide %>%
		bind_rows(results_i_wide)
}

}

### After the firs loop over site add to bottom
if (i == 1) {
  allsites_wide <- site_wide
} else {
  allsites_wide <- allsites_wide %>%
    bind_rows(site_wide)
}

}

### Rearrange into a 3D array (3rd dimension is a bout, 2nd dimension is a repeated survey)
names <- colnames(allsites_wide)
names <- as.data.frame(names[6:length(names)])
colnames(names) <- "names"
names <- names %>%
  separate(col=names, into = c("year", "month", "rep"), sep="_", remove = FALSE)

months <- unique(names$month)
months2 <- as.data.frame(months)
months2$rep <- 0
for (t in 1:length(months)){
  months2[t,2] <- nrow(names[names$month == months[t],])
}

maxs <- max(months2$rep) + 5
s1_app_all <- array(NA,c(nrow(allsites_wide), maxs,nrow(months2)))
s1_app_all <- list()

for (s in 1:nrow(months2)){
  no <- months2[s,1]
  nms <- names %>%
    filter(month == no)
  s1 <- allsites_wide %>%
    select(tagID, trap_id, siteID, plotID, trapCoordinate, c(nms[,1]))
  
  if (ncol(s1) == maxs) {
    s1_app <- s1
  } else {
    f1 <- matrix(NA, nrow(s1), (maxs - ncol(s1)))
    f1 <- as.data.frame(f1)
    
    s1_app <- s1 %>%
      bind_cols(f1)
    
    colnames(s1_app)[6:(5+max(months2$rep))] <- paste0("rep_", 1:max(months2$rep))
    
  }
      s1_app_all[[s]] <- s1_app
  }

### Save a single object to a file
saveRDS(s1_app_all, paste0(year_j, "_allsites_table_C1.rds"))
}
