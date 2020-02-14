# install.packages("neonUtilities")
library(neonUtilities)
# ?loadByProduct
# ?getDatatable

#### Download organism data ####
##****************************##

# You can look in the data product catalog (https://data.neonscience.org/apps/browse) and manually figure out what the product codes are for small mammal trap data and for bird point count data, but I've provided them here. The `DP1` in the code indicates that this is Level 1 data. For Level 1 data, quality controls were run (Level 0 would be `DP0` meaning completely raw data) but the actual values are still raw values measured in the field, not some kind of calculated quantity (Level 2 and higher would be derived values).

# Breeding landbird point counts
# http://data.neonscience.org/data-product-view?dpCode=DP1.10003.001
bird_code <- 'DP1.10003.001'
# Fish electrofishing, gill netting, and fyke netting counts 
# http://data.neonscience.org/data-product-view?dpCode=DP1.20107.001
fish_code <- 'DP1.20107.001'
# Aquatic plant, bryophyte, lichen, and macroalgae point counts in wadeable streams
# http://data.neonscience.org/data-product-view?dpCode=DP1.20072.001
aqua_plant <- 'DP1.20072.001'
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
riparian_code<- 'DP1.20275.001'
# Plant presence and percent cover 
# http://data.neonscience.org/data-product-view?dpCode=DP1.10058.001
plant_code <- 'DP1.10058.001'
# Small mammal box trapping
# http://data.neonscience.org/data-product-view?dpCode=DP1.10072.001
mammal_code <- 'DP1.10072.001'
# Soil microbe community composition
# http://data.neonscience.org/data-product-view?dpCode=DP1.10081.001
microbe_code <-'DP1.10081.001'
# Ticks sampled using drag cloths 
# http://data.neonscience.org/data-product-view?dpCode=DP1.10093.001
tick_code <- 'DP1.10093.001'
# Woody plant vegetation structure
# http://data.neonscience.org/data-product-view?dpCode=DP1.10098.001
woody_code <- 'DP1.10098.001'
# Zooplankton collection 
# http://data.neonscience.org/data-product-view?dpCode=DP1.20219.001
zoop_code <- 'DP1.20219.001'

# Let's take a look at what files are available for NEON small mammal trapping data for a given site-month combination. Running this takes a minute or two and requires an internet connection because we are querying the API.


my_site_list <- c('OSBS')

d <- neonUtilities::loadByProduct(
  dpID = plant_code,
  site = my_site_list,
  check.size = TRUE)
str(d)
head(d$variables)
head(d$div_1m2Data)
