# NEONss2019-biodiv

Biodiversity group (WG#17) at NEON Science Summit 2019 led by PL Zarnetske, S Record, D Li

## Links

* Google Drive (Phoebe's): https://drive.google.com/drive/folders/1kQuCMKh7lmeTaeBHQbEY44ivGZUif3KS
  * Includes data and manuscript(s)
* Google Drive (Science Summit Group 17 subdir): https://drive.google.com/drive/folders/1xtnmuqH2wXoV4sKoiRqyCpPwgLoqUv-3
* NEON Science Summit 2019 Wiki https://github.com/earthlab/neon-science-summit/wiki
* Slack channel for WG#14 & WG#17  https://neonworkinggroup.slack.com
* The related group WG#14 (biodiversity change group) github https://github.com/karinorman/temporalNEON

## Organizing code

In the `code` folder, let's try to put all the packages and custom functions into the `00_pkg_functions.R`.

Each taxon group can creat a file named as `01_data_taxa.R`. For example, for plants, we can put all the code to download, clean, and save plant sampling data in `01_data_plant.R`. Within each file, we can group code by sections such as "Download data", "Clean data", "Save data", etc.

After we have all uploaded code to clean different taxanomic groups, we can then create a file named as something like `02_neon_ecocomDP.R` to convert all cleaned data to the ecocommDP format.

```r
neonUtilities:::table_types

# You can look in the data product catalog (http://data.neonscience.org/data-product-catalog) and manually figure out what the product codes are for small mammal trap data and for bird point count data, but I've provided them here. The `DP1` in the code indicates that this is Level 1 data. For Level 1 data, quality controls were run (Level 0 would be `DP0` meaning completely raw data) but the actual values are still raw values measured in the field, not some kind of calculated quantity (Level 2 and higher would be derived values).

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
```
