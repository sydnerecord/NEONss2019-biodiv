#title: download algae data and clean
#author: Lara Jansen
#date: 1/02/2020

#packages used
library(dplyr)
library(neonUtilities)

#get algae data -------------
# algaeDiv dpid
periphyton_code <-'DPI.20166.001'

#only two sites - site_id = c("MAYF","PRIN")-EX
alg_allTabs <- loadByProduct(dpID = "DP1.20166.001", 
                             site = c("MAYF", "PRIN"),
                             startdate = "2016-1", enddate = "2018-11",  
                             package = "expanded", check.size = TRUE)
# pull out the field data into a data.frame
alg_field_data <- alg_allTabs$alg_fieldData

# pull out the taxonomy data
alg_tax_long <- alg_allTabs$alg_taxonomyProcessed

###Cleaning for ecocomDP ##
#pull biomass table for (ecocomDP) table_observation as it has additional estimated 'perBottleSampleVolume's
alg_biomass <-alg_allTabs$alg_biomass
#add up labSampleVolume + preservativeVolume = estimated 'perBottleSampleVolume'
alg_biomass1t<- alg_biomass %>%
  mutate(estBSVolume=preservativeVolume+labSampleVolume) %>%
  filter(analysisType=='taxonomy')
