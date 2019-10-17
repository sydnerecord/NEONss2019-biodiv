loc = read_csv("data/table_location.csv") %>% rename(loc = location_id)
taxa = read_csv("data/table_taxon.csv") %>% rename(sp = taxon_id)
comm = read_csv("data/table_observation.csv")
names(comm)
comm = select(comm, loc = location_id, datetime = observation_datetime, 
              sp = taxon_id, var = variable_name, value, unit)
table(comm$var)
unique(comm$var)
comm = select(comm, -var, -unit) %>% rename(ct_per_sq_m = value)

setdiff(unique(comm$sp), taxa$sp)
setdiff(taxa$sp, unique(comm$sp))
n_distinct(comm$sp)

comm = left_join(comm, taxa, by = "sp") %>% 
  select(-sp) %>% 
  rename(sp = taxon_name)

comm = mutate(comm, yr = lubridate::year(datetime))

# ----------------------------------------
# Load all organismal data and do QC on it
# ----------------------------------------

# QC for all taxa
# ---------------------
# we skipped this step for now. Oct/16/2019

# ----------------------------  
# Richness estimator functions
# ----------------------------

# Added by QDR, 19 June 2018
# Modified by QDR, 20 June 2018: Create functions for cumulative richness by year and for reshaping data to plot.

xtest = estimator_chao1(comm$sp)
xtest2 = estimator_asymp(comm$sp, q = 0, nboot = 99)

richness_cumulative(comm, sp_name = "sp", site_name = "loc", value_name = "ct_per_sq_m", year_name = "yr")
richness_cumulative(comm, sp_name = "sp", site_name = "loc", value_name = "ct_per_sq_m", year_name = "yr", cumul_by_year = F)

comm2 = bind_rows(mutate(comm, plot = "plot1"), mutate(comm, plot = "plot2"))
richness_cumulative(comm2, sp_name = "sp", site_name = "loc", value_name = "ct_per_sq_m", year_name = "yr", plot_name = "plot", grp_by = "site_plot")
richness_cumulative(comm2, sp_name = "sp", site_name = "loc", value_name = "ct_per_sq_m", year_name = "yr", plot_name = "plot", cumul_by_year = F, grp_by = "site_plot")



# Simplified function for phylogenetic and functional diversity, with null models (added 29 June)
# Modified from nasa code, and required for function below
# Only does it based on incidence (no abundance weighting)
pd_fd <- function(sp_list, pddist, fddist, nnull = 99, phylo_spp = NULL, func_spp = NULL) {
  
  require(picante)
  
  # convert species list to a table
  m <- t(as.matrix(table(sp_list)))
  
  # Get rid of species that aren't in phylogenetic and functional diversity.
  if (!is.null(phylo_spp)) mphy <- m[, dimnames(m)[[2]] %in% phylo_spp, drop = FALSE] else mphy <- m
  if (!is.null(func_spp)) mfunc <- m[, dimnames(m)[[2]] %in% func_spp, drop = FALSE] else mfunc <- m
  
  if (dim(mphy)[2] > 1) {
    # Calculate PD and do null models
    MPD <- ses.mpd(mphy, pddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
    MNTD <- ses.mntd(mphy, pddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
    # Get all z-scores
    MPD_z <- (MPD$mpd.obs[1] - mean(MPD$mpd.rand.mean, na.rm=TRUE))/sd(MPD$mpd.rand.mean, na.rm=TRUE)
    MNTD_z <- (MNTD$mntd.obs[1] - mean(MNTD$mntd.rand.mean, na.rm=TRUE))/sd(MNTD$mntd.rand.mean, na.rm=TRUE)
    MPD_obs <- MPD$mpd.obs[1]
    MNTD_obs <- MNTD$mntd.obs[1]
  } else {
    MPD_z <- MNTD_z <- MPD_obs <- MNTD_obs <- NA
  }
  if (dim(mfunc)[2] > 1) {
    # Calculate FD and do null models
    MPDfunc <- ses.mpd(mfunc, fddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
    MNTDfunc <- ses.mntd(mfunc, fddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
    # Get all z-scores
    MPDfunc_z <- (MPDfunc$mpd.obs[1] - mean(MPDfunc$mpd.rand.mean, na.rm=TRUE))/sd(MPDfunc$mpd.rand.mean, na.rm=TRUE)
    MNTDfunc_z <- (MNTDfunc$mntd.obs[1] - mean(MNTDfunc$mntd.rand.mean, na.rm=TRUE))/sd(MNTDfunc$mntd.rand.mean, na.rm=TRUE)
    MPDfunc_obs <- MPDfunc$mpd.obs[1]
    MNTDfunc_obs <- MNTDfunc$mntd.obs[1]
  } else {
    MPDfunc_z <- MNTDfunc_z <- MPDfunc_obs <- MNTDfunc_obs <- NA
  }
  # Concatenate observed values and z-scores into a vector and return.
  return(data.frame(MPD_phy = MPD_obs, MPD_phy_z = MPD_z, 
                    MNTD_phy = MNTD_obs, MNTD_phy_z = MNTD_z,
                    MPD_func = MPDfunc_obs, MPD_func_z = MPDfunc_z, 
                    MNTD_func = MNTDfunc_obs, MNTD_func_z = MNTDfunc_z))
  
}

# Function to do functional and phylogenetic cumulatively, instead of taxonomic (added 29 June)
# Requires p_dist and f_dist which are matrices with row and column names corresponding to the taxonID, and entries being the pairwise distances
# Also set number of null iterations
PDFD_cumulative <- function(dat, p_dist, f_dist, column = 'taxonID', group = 'site', plot_name = 'plotID', by_year = TRUE, n_iter = 99) {
  
  physpp <- dimnames(p_dist)[[1]]
  funspp <- dimnames(f_dist)[[1]]
  
  dat <- dat %>%
    rename(sp = !!column, plotID = !!plot_name) %>%
    select(siteID, plotID, year, sp)
  if (group == 'site') dat <- dat %>% group_by(siteID)
  if (group == 'plot') dat <- dat %>% group_by(siteID, plotID)
  out <- dat %>%
    do(cbind(year = min(.$year):max(.$year),
             map_dfr(min(.$year):max(.$year), function(yr) pd_fd(.$sp[.$year <= yr], p_dist, f_dist, n_iter, physpp, funspp))
    ))
  if (!by_year) {
    out %>% filter(year == max(year)) %>% select(-year)
  } else {
    out
  }
}

# 

# ----------------------------------------------------------
# Run richness estimators for all taxa by site and year
# ----------------------------------------------------------

# Added by QDR, 19 June 2018
# Modified by QDR, 20 June 2018: Add other taxa

# Site-level
mammal_richness_cumulative <- richness_cumulative(mammal_data)
bird_richness_cumulative <- richness_cumulative(bird_data)
aq_plant_richness_cumulative <- richness_cumulative(aq_plant_data, column = 'scientificName', plot_name = 'namedLocation')
macroinvert_richness_cumulative <- richness_cumulative(macroinvert_data, column = 'acceptedTaxonID', plot_name = 'namedLocation')
fish_richness_cumulative <- richness_cumulative(fish_data, plot_name = 'namedLocation')
mosq_richness_cumulative <- richness_cumulative(mosq_data, column = 'scientificName')
plant_richness_cumulative <- richness_cumulative(plant_data)
tree_richness_cumulative <- richness_cumulative(tree_data)
beetle_richness_cumulative <- richness_cumulative(beetle_data)
zoop_richness_cumulative <- richness_cumulative(zoop_data, plot_name = 'namedLocation')
phytop_richness_cumulative <- richness_cumulative(phytop_data, column = 'scientificName', plot_name = 'namedLocation')
littletree_richness_cumulative <- richness_cumulative(little_tree_data, column = 'binomial')

# Plot-level
mammal_richness_cumulative_plot <- richness_cumulative(mammal_data, group='plot')
bird_richness_cumulative_plot <- richness_cumulative(bird_data, group='plot')
aq_plant_richness_cumulative_plot <- richness_cumulative(aq_plant_data, group='plot', column = 'scientificName', plot_name = 'namedLocation')
macroinvert_richness_cumulative_plot <- richness_cumulative(macroinvert_data, group='plot', column = 'acceptedTaxonID', plot_name = 'namedLocation')
fish_richness_cumulative_plot <- richness_cumulative(fish_data,group='plot', plot_name = 'namedLocation')
mosq_richness_cumulative_plot <- richness_cumulative(mosq_data, group='plot', column = 'scientificName')
plant_richness_cumulative_plot <- richness_cumulative(plant_data,group='plot')
tree_richness_cumulative_plot <- richness_cumulative(tree_data,group='plot')
beetle_richness_cumulative_plot <- richness_cumulative(beetle_data,group='plot')
zoop_richness_cumulative_plot <- richness_cumulative(zoop_data,group='plot', plot_name = 'namedLocation')
phytop_richness_cumulative_plot <- richness_cumulative(phytop_data, group='plot', column = 'scientificName', plot_name = 'namedLocation')
littletree_richness_cumulative_plot <- richness_cumulative(little_tree_data, group='plot', column = 'binomial')

# Site level PD and FD: birds and trees
# Also do mammals, note that not all species are represented in phylogeny
bird_pdfd_cumulative <- PDFD_cumulative(bird_data, p_dist = bird_p_dist, f_dist = bird_f_dist, n_iter = 999)
littletree_pdfd_cumulative <- little_tree_data %>%
  filter(!siteID %in% c('GUAN','LAJA','NOGP','JORN')) %>% # Sites in PR or with too few "Little" trees
  PDFD_cumulative(p_dist = tree_p_dist, f_dist = tree_f_dist, n_iter = 999)
mam_pdfd_cumulative <- PDFD_cumulative(mammal_data, p_dist = mam_p_dist, f_dist = mam_f_dist, n_iter = 999)

# Site level PD/FD: single value across all years
bird_pdfd_total <- PDFD_cumulative(bird_data, p_dist = bird_p_dist, f_dist = bird_f_dist, n_iter = 999, by_year = FALSE)
littletree_pdfd_total <- little_tree_data %>%
  filter(!siteID %in% c('GUAN','LAJA','NOGP','JORN')) %>% # Sites in PR or with too few "Little" trees
  PDFD_cumulative(p_dist = tree_p_dist, f_dist = tree_f_dist, n_iter = 999, by_year = FALSE)
mam_pdfd_total <- PDFD_cumulative(mammal_data, p_dist = mam_p_dist, f_dist = mam_f_dist, n_iter = 999, by_year = FALSE)


# Export richness data by taxonomic group: site-level
# --------------------------------------
write.csv(mammal_richness_cumulative,file.path(export_path,"mammal_richness_cumulative.csv"),row.names=F)
write.csv(bird_richness_cumulative,file.path(export_path,"bird_richness_cumulative.csv"),row.names=F)
write.csv(beetle_richness_cumulative,file.path(export_path,"beetle_richness_cumulative.csv"),row.names=F)
write.csv(macroinvert_richness_cumulative,file.path(export_path,"macroinvert_richness_cumulative.csv"),row.names=F)
write.csv(mosq_richness_cumulative,file.path(export_path,"mosq_richness_cumulative.csv"),row.names=F)
write.csv(phytop_richness_cumulative,file.path(export_path,"phytop_richness_cumulative.csv"),row.names=F)
write.csv(plant_richness_cumulative,file.path(export_path,"plant_richness_cumulative.csv"),row.names=F)
write.csv(tree_richness_cumulative,file.path(export_path,"tree_richness_cumulative.csv"),row.names=F)
write.csv(zoop_richness_cumulative,file.path(export_path,"zoop_richness_cumulative.csv"),row.names=F)
# Export functional and phylogenetic
write.csv(bird_pdfd_cumulative, file.path(export_path, 'bird_pdfd_cumulative.csv'), row.names = FALSE)
write.csv(littletree_pdfd_cumulative, file.path(export_path, 'tree_pdfd_cumulative.csv'), row.names = FALSE)
write.csv(mam_pdfd_cumulative, file.path(export_path, 'mammal_pdfd_cumulative.csv'), row.names = FALSE)
write.csv(bird_pdfd_total, file.path(export_path, 'bird_pdfd_total.csv'), row.names = FALSE)
write.csv(littletree_pdfd_total, file.path(export_path, 'tree_pdfd_total.csv'), row.names = FALSE)
write.csv(mam_pdfd_total, file.path(export_path, 'mammal_pdfd_total.csv'), row.names = FALSE)

# Also combine everything into one data frame and write to a single CSV for easier loading.
alltaxa_richness_cumulative <- rbind(
  data.frame(taxon = 'mammal', mammal_richness_cumulative),
  data.frame(taxon = 'bird', bird_richness_cumulative),
  data.frame(taxon = 'beetle', beetle_richness_cumulative),
  data.frame(taxon = 'mosquito', mosq_richness_cumulative),
  data.frame(taxon = 'all plants', plant_richness_cumulative),
  data.frame(taxon = 'all woody plants', tree_richness_cumulative),
  data.frame(taxon = 'trees', littletree_richness_cumulative),
  data.frame(taxon = 'macroinvertebrate', macroinvert_richness_cumulative),
  data.frame(taxon = 'fish', fish_richness_cumulative),
  data.frame(taxon = 'zooplankton', zoop_richness_cumulative),
  data.frame(taxon = 'phytoplankton', phytop_richness_cumulative),
  data.frame(taxon = 'aquatic plants', aq_plant_richness_cumulative)
)

write.csv(alltaxa_richness_cumulative, file.path(export_path, 'alltaxa_richness_cumulative.csv'), row.names = FALSE)

## Edit 26 June 2018 PLZ
# Export the plot-level richness (when including plotID above in richness calculation)
# Export richness data by taxonomic group: plot-level
# --------------------------------------
write.csv(mammal_richness_cumulative_plot,file.path(export_path,"mammal_richness_cumulative_plot.csv"),row.names=F)
write.csv(bird_richness_cumulative_plot,file.path(export_path,"bird_richness_cumulative_plot.csv"),row.names=F)
write.csv(beetle_richness_cumulative_plot,file.path(export_path,"beetle_richness_cumulative_plot.csv"),row.names=F)
write.csv(macroinvert_richness_cumulative_plot,file.path(export_path,"macroinvert_richness_cumulative_plot.csv"),row.names=F)
write.csv(mosq_richness_cumulative_plot,file.path(export_path,"mosq_richness_cumulative_plot.csv"),row.names=F)
write.csv(phytop_richness_cumulative_plot,file.path(export_path,"phytop_richness_cumulative_plot.csv"),row.names=F)
write.csv(plant_richness_cumulative_plot,file.path(export_path,"plant_richness_cumulative_plot.csv"),row.names=F)
write.csv(tree_richness_cumulative_plot,file.path(export_path,"tree_richness_cumulative_plot.csv"),row.names=F)
write.csv(zoop_richness_cumulative_plot,file.path(export_path,"zoop_richness_cumulative_plot.csv"),row.names=F)


# -------------------------------------------------------------------------------------------------------------------------
## End of Code work 21 June 2018 
# Start with Small Mammals
## Work through diversity measure code next: use simple way or Chao
## Simple Richness: sum up by year, site

# Richness by Site and Year
mammal_richness_yr <- mammal_data %>%
  group_by(siteID, year) %>%
  summarize(richness = length(unique(taxonID)))

bird_richness_yr <- bird_data %>%
  group_by(siteID, year) %>%
  summarize(richness = length(unique(taxonID)))

# Overall Richness by Site
mammal_richness_site <- mammal_data %>%
  group_by(siteID) %>%
  summarize(richness = length(unique(taxonID)))

bird_richness_site <- bird_data %>%
  group_by(siteID) %>%
  summarize(richness = length(unique(taxonID))) 

## Repeat for Chao... evenness... others?

## Repeat, but demonstrate improving estimates (or at least different estimates)
## based on repeat sampling with some sites having multiple sampling bouts in a yr, 
## and multiple years. Demonstrate with rarefaction curves? 
## Should we do something with mark-recap data to get at abundance? there are
## several R packages 

## Beta diversity... ideas? 

## Compare with IUCN range maps

## 
