#title: download zooplankton data 
#author: Stephanie Parker
#date: 06/10/2020

#packages used
library(tidyverse)

# Install and load devtools
# install.packages("devtools")
library(devtools)


# Install and load neonUtilities
# install_github("NEONScience/NEON-utilities/neonUtilities", dependencies=TRUE)
library(neonUtilities)

#################################################################################
#get zooplankton data#
# zooplankton dpid
my_dpid <- 'DP1.20219.001'
my_site_list <- c('PRLA', 'CRAM') #data for 2 sites is more manageable!

all_tabs <- neonUtilities::loadByProduct(
  dpID = my_dpid,
  site = my_site_list,
  check.size = TRUE)


# download field data for all dates for two neon sites 
zoo_fielddata <- all_tabs$zoo_fieldData

# download zooplankton counts for two sites 
zoo_taxonomyProcessed <- all_tabs$zoo_taxonomyProcessed