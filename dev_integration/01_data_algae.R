# contact: Lara Jansen (ljansen@pdx.edu)/ Eric Sokol (esokol@battelleecology.org)

# Daijiang here: I don't understand why left join parentSampleID with sampleID
# all sampleID matched but not sampleID and parentsampleID.
source("code/00_pkg_functions.R")

# the algae data product ID:
# algae DP1.20166.001

# # all data -- could be big
# alg_allTabs <- loadByProduct(dpID = "DP1.20166.001", 
#                              site = "all", package = "expanded", 
#                              check.size = TRUE)

# Try this:
# only 2 sites, restricted in time
alg_allTabs <- loadByProduct(dpID = "DP1.20166.001", 
                             site = c("MAYF", "PRIN"),
                             startdate = "2016-1", enddate = "2018-11",  
                             package = "expanded", check.size = TRUE)

# pull out the field data into a data.frame
alg_field_data <- alg_allTabs$alg_fieldData
alg_tax_long <- alg_allTabs$alg_taxonomyProcessed
alg_biomass <-alg_allTabs$alg_biomass

alg_biomass1t <- alg_biomass %>%
  mutate(estBSVolume = preservativeVolume + labSampleVolume) %>%
  filter(analysisType == 'taxonomy') # filter to taxonomy type

alg_tax_long$perBottleSampleVolume[is.na(alg_tax_long$perBottleSampleVolume)] <- 0
# not sure whether the above line is a good idea

alg_tax_biomass2 <- left_join(alg_tax_long, alg_biomass1t) %>%
  mutate(perBSVol=
           case_when(
             perBottleSampleVolume == 0 ~ estBSVolume,
             perBottleSampleVolume > 0 ~ perBottleSampleVolume))

## create table_observation
table_observation1 <- left_join(alg_tax_biomass2, alg_field_data) %>%
  select(
    uid,
    sampleID,
    siteID,
    collectDate,
    algalParameterValue,
    algalParameterUnit,
    algalParameter,
    perBSVol,
    fieldSampleVolume,
    algalSampleType,
    benthicArea,
    acceptedTaxonID,
    scientificName
  ) %>%
  filter(algalParameterUnit == 'cellsPerBottle') %>% # filter to cells per bottle
  mutate(
    density = case_when(
        algalSampleType %in% c('seston') ~ algalParameterValue / perBSVol,
        TRUE ~ (algalParameterValue / perBSVol) * (fieldSampleVolume / benthicArea) #add phytoplankton back in when applicable
      ),
    cell_density_standardized_unit = case_when(
      algalSampleType == 'phytoplankton' ~ 'cells/mL',
      TRUE ~ 'cells/m2'
    )
  )
