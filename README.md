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
