# clean out workspace
rm(list = ls())
gc()

# options
options(stringsAsFactors = FALSE)

# load packages
library(tidyverse)

devtools::install_local(path = 'C:/Users/esokol/Documents/Git/FORKS/NEONss2019-biodiv/dev_integration/neonBiodivTools')
library(neonBiodivTools)

##############################
# NEON biodiversity WG data folder on google
# google_id = '1WLWIW9K5o4tHbBBdCjb1skaWDs4Jn7vj'

my_google_id <- '1WLWIW9K5o4tHbBBdCjb1skaWDs4Jn7vj' %>% googledrive::as_id()
my_google_id %>% googledrive::drive_ls()

#fsh_data_long.csv id = '1EU9Xmu4TzuH7zlWZgMjpbkIiMV-OYZvJ'

dat_in <- neonBiodivTools::read_from_google_drive(
  file_name_string = 'fsh_data_long.csv',
  my_path_to_googledirve_directory = my_google_id,
  keep_local_copy_of_file = FALSE)


##############################
