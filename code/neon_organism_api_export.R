## TITLE:         NEON Organismal Data: read in via API, export data
## AUTHOR:        Phoebe Zarnetske, Quentin Read 
## COLLABORATORS: Sydne Record (Bryn Mawr), Ben Baiser (UFL), Angela Strecker (PSU), 
##                John M. Grady (MSU/Bryn Mawr), Jonathan Belmaker (Tel Aviv U), Mao-Ning Tuanmu (Academia Sinica),
##                Lydia Beaudrot (Rice U), Kate Thibault 
## DATA:          NEON organismal data: all species, all years, all sites
## PROJECT:       "NEON's continental-scale biodiversity"
## DATE:          initiated: June 18, 2018; last run:

## This script reads in NEON's organismal data across all available sites, 
# computes diversity measures per site and year, and cumulatively,
# and exports those data. The API portion of the script is based on
# QDR's neon_api_grad_lab_rawcode.R available at: https://github.com/NEON-biodiversity/teaching/tree/master/grad_lab

#Clear all existing data
rm(list=ls())

#Close graphics devices
graphics.off()

#set working directory
setwd("/Volumes/GoogleDrive/My Drive/Research/ScalingUp/NEON_EAGER/Manuscript4_NEON_Organisms") # GD location
#setwd("/Volumes/neon/final_data") # HPCC location

#Install/load packages
for (package in c("ggplot2", "lme4", "dplyr")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
# This code for ggplot2 sets the theme to mostly black and white 
# (Arial font, and large font, base size=24)
theme_set(theme_bw(12))
theme_update(axis.text.x = element_text(size = 10, angle = 90),
             axis.text.y = element_text(size = 10))

## Code below from https://github.com/NEON-biodiversity/teaching/tree/master/grad_lab/neon_api_grad_lab_rawcode.R


#### R functions for pulling NEON data from the server ####
##*******************************************************##

## Function to display what data are available

# This function takes a NEON product code `productCode` as an argument, gets a list of the files that are available, and displays a representative set of file names from one site-month combination.

display_neon_filenames <- function(productCode) {
  require(httr)
  require(jsonlite)
  req <- GET(paste0("http://data.neonscience.org/api/v0/products/", productCode))
  avail <- fromJSON(content(req, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)
  urls <- unlist(avail$data$siteCodes$availableDataUrls)
  get_filenames <- function(x) fromJSON(content(GET(x), as = 'text'))$data$files$name
  files_test <- sapply(urls, get_filenames, simplify = FALSE)
  files_test[[which.max(sapply(files_test,length))]]
}

## Function to pull all data for a given data product

# The first argument, `productCode` is a NEON product code, and the second, `nametag`, is an identifying string that tells the function which CSV to download for each site-month combination. There are usually a lot of metadata files that we aren't interested in for now that go along with the main data file for each site-month combination, and the `nametag` argument tells the function which file is the one that really has the data we want (details below). The `pkg` argument defaults to download the "basic" data package which is usually all we would want. Finally, the `bind` argument defaults to `TRUE` which means return a single data frame, not a list of data frames.

# There are two steps to what the function does: first it queries the API to get a list of URLs of the CSV files available for all site-month combinations for the desired data product. Second it loops through the subset of those URLs that match the `nametag` argument and tries to download them all. If one gives an error, there is a `try()` function built in so that the function will just skip that file instead of quitting.

pull_all_neon_data <- function(productCode, nametag, pkg = 'basic', bind = TRUE) {
  require(httr)
  require(jsonlite)
  require(dplyr)
  
  # Get list of URLs for all site - month combinations for that data product.
  req <- GET(paste0("http://data.neonscience.org/api/v0/products/", productCode))
  avail <- fromJSON(content(req, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)
  urls <- unlist(avail$data$siteCodes$availableDataUrls)
  
  # Loop through and get the data from each URL.
  res <- list()
  
  pb <- txtProgressBar(min=0, max=length(urls), style=3)
  count <- 0
  
  for (i in urls) {
    count <- count + 1
    setTxtProgressBar(pb, count)
    # Get the URLs for the site-month combination
    req_i <- GET(i)
    files_i <- fromJSON(content(req_i, as = 'text'))
    urls_i <- files_i$data$files$url
    # Read data from the URLs given by the API, skipping URLs that return an error.
    data_i <- try(read.delim(
      grep(paste0('(.*',nametag,'.*', pkg, '.*)'), urls_i, value = TRUE), 
      sep = ',', stringsAsFactors = FALSE), TRUE)
    if (!inherits(data_i, 'try-error')) res[[length(res) + 1]] <- data_i
  }
  
  close(pb)
  
  # Return as a single data frame or as a list of data frames, 
  # depending on what option was selected.
  if (bind) {
    do.call(rbind, res)
  } else {
    res
  }
}

## Function to get spatial information (coordinates) for a site or plot

# The spatial locations and metadata for sites and plots are stored in a different location on NEON's API from the actual data. This function should be called for a single site at a time (`siteID` argument). The second argument, `what`, is a string. The default is `"site"` which will return a single row of a data frame with spatial location for the entire site as a single point. If `what` is set to another string such as `"bird"` it will go through the spatial location data URLs, find all that have `"bird"` in the name, pull the spatial information from them, and return a data frame with one row per bird plot. In either case the data frame has 19 columns (the number of location attributes NEON has listed for each site or plot).

get_site_locations <- function(siteID, what = 'site') {
  require(httr)
  require(jsonlite)
  require(purrr)
  # URLs of all spatial information about the site
  req <- GET(paste0("http://data.neonscience.org/api/v0/locations/", siteID))
  site_loc <- fromJSON(content(req, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)
  
  if (what == 'site') {
    # If only coordinates of the entire site are needed, return them
    return(data.frame(site_loc$data[1:19]))
  } else {
    # If "what" is some type of plot, find all URLs for that plot type
    urls <- grep(what, site_loc$data$locationChildrenUrls, value = TRUE)
    # Get the coordinates for each of those plots from each URL and return them
    loc_info <- map_dfr(urls, function(url) {
      req <- GET(url)
      loc <- fromJSON(content(req, as = 'text'), simplifyDataFrame = TRUE, flatten = TRUE)
      loc[[1]][1:19]
    })
    return(loc_info)
  }
}

#### Download organism data ####
##****************************##

# You can look in the data product catalog (http://data.neonscience.org/data-product-catalog) and manually figure out what the product codes are for small mammal trap data and for bird point count data, but I've provided them here. The `DP1` in the code indicates that this is Level 1 data. For Level 1 data, quality controls were run (Level 0 would be `DP0` meaning completely raw data) but the actual values are still raw values measured in the field, not some kind of calculated quantity (Level 2 and higher would be derived values).

# Breeding landbird point counts
# http://data.neonscience.org/data-product-view?dpCode=DP1.10003.001
bird_code <- 'DP1.10003.001'
# Fish electrofishing, gill netting, and fyke netting counts 
# http://data.neonscience.org/data-product-view?dpCode=DP1.20107.001
fish_code <- 'DP1.20107.001'
# Aquatic plant, bryophyte, lichen, and macroalgae point counts in wadeable streams
# http://data.neonscience.org/data-product-view?dpCode=DP1.20072.001
aquaplant <- 'DP1.20072.001'
# Ground beetles sampled from pitfall traps
# http://data.neonscience.org/data-product-view?dpCode=DP1.10022.001
beetle_code <- 'DP1.10022.001'
# Macroinvertebrate collection
# http://data.neonscience.org/data-product-view?dpCode=DP1.20120.001
macroinv_code <- 'DP1.20120.001'
# Mosquitoes sampled from CO2 traps 
# http://data.neonscience.org/data-product-view?dpCode=DP1.10043.001
mosquito_code <- 'DP1.10043.001'
# Periphyton, seston, and phytoplankton collection 
# http://data.neonscience.org/data-product-view?dpCode=DP1.20166.001
periphyton_code <- 'DP1.20166.001'
# Riparian composition and structure
# http://data.neonscience.org/data-product-view?dpCode=DP1.20275.001
riparian_code<-'DP1.20275.001'
# Plant presence and percent cover 
# http://data.neonscience.org/data-product-view?dpCode=DP1.10058.001
plant_code<-'DP1.10058.001'
# Small mammal box trapping
# http://data.neonscience.org/data-product-view?dpCode=DP1.10072.001
mammal_code <- 'DP1.10072.001'
# Soil microbe community composition
# http://data.neonscience.org/data-product-view?dpCode=DP1.10081.001
microbe_code<-'DP1.10081.001'
# Ticks sampled using drag cloths 
# http://data.neonscience.org/data-product-view?dpCode=DP1.10093.001
tick_code <-'DP1.10093.001'
# Woody plant vegetation structure
# http://data.neonscience.org/data-product-view?dpCode=DP1.10098.001
woody_code<-'DP1.10098.001'
# Zooplankton collection 
# http://data.neonscience.org/data-product-view?dpCode=DP1.20219.001
zoop_code<-'DP1.20219.001'

# Let's take a look at what files are available for NEON small mammal trapping data for a given site-month combination. Running this takes a minute or two and requires an internet connection because we are querying the API.

display_neon_filenames(mammal_code)

# You can see that there are a lot of files available for one site. However the one we are interested in is the file containing the mammals caught per trap per night in the basic data package (expanded data package contains other variables that might be needed for quality control but that we are not interested in here). Let's pull that CSV file for all site-month combinations and combine it into one huge data frame that we can run analysis on. We specify we want everything belonging to the mammal code that contains the string `pertrapnight` in the file name, and by default only get the basic data package. Running this code on your own machine will take quite a few minutes since it has to download a lot of data, but you should get a progress bar showing how much time is remaining.

mammal_data <- pull_all_neon_data(productCode = mammal_code, 
                                  nametag = 'pertrapnight')

mammal_data_plot <- pull_all_neon_data(productCode = mammal_code, 
                                       nametag = 'perplotnight')

# Now let's take a look at what is in that data frame . . . 

str(mammal_data)

# The mammal data frame is huge, with over 600K rows. Most of the rows record trap-nights where no mammal was captured. Let's get rid of those.

nrow(mammal_data)
table(mammal_data$trapStatus)

# You can see that only status 4 and 5 correspond to one or more mammals caught in the trap. Filter the data frame to only keep those rows. We use the function `grepl()` which matches a regular expression to a vector of strings and returns `TRUE` if they match. The regular expression `"4|5"` means any string with the numerals 4 or 5 in it.

mammal_data <- mammal_data %>%
  filter(grepl('4|5', trapStatus))

nrow(mammal_data)

# Make a year column so we can subset by year.
mammal_data$year<-strptime(mammal_data$collectDate, "%Y-%m-%d")$year+1900

# How many datasets are available?
table(bird_data$siteID, lubridate::year(bird_data$startDate))

table(mammal_data$siteID, lubridate::year(mammal_data$collectDate))

# export file to HPCC as final data
write.csv(mammal_data,file="/Volumes/neon/raw_data/organismal_data_june2018/mammal_data.csv",row.names=F)

## Bird data

# Next, let's look at what data are available for birds.

display_neon_filenames(bird_code)

# Since the string `count` is in the name of the data file that we want for each site-month combination (the raw point count data for birds), we use that to pull point count data for each month and stick it all into one big data frame.

bird_data <- pull_all_neon_data(productCode = bird_code, 
                                nametag = 'count')

# Let's see what is in that data frame . . . 

str(bird_data)

# Find the species richness at each site by counting the number of unique taxa.

bird_richness <- bird_data %>%
  group_by(siteID) %>%
  summarize(richness = length(unique(taxonID)))

write.csv(bird_data,file="/Volumes/neon/raw_data/organismal_data_june2018/mammal_data.csv",row.names=F)
