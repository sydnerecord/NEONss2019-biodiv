# Plots of cumulative richness by taxon and year
# QDR / NEON biodiversity / 24 June 2018

# Modified 26 June 2018: Add Alaska maps, Puerto Rico maps, and scatter plots.

# Data to make these plots were downloaded to google drive so script can be run locally

library(dplyr)
library(ggplot2)
library(reshape2)

# Set file paths
google_drive <- 'C:/Users/Q/google_drive/NEON_EAGER' # if you are Q
google_drive <- '/Volumes/GoogleDrive/My Drive/Research/ScalingUp/NEON_EAGER' # if you are Phoebe
setwd(file.path(google_drive, "Manuscript4_NEON_Organisms")) # GD location
data_path <- file.path(google_drive, 'Manuscript4_NEON_Organisms/data')
fig_path <- file.path(google_drive, 'Manuscript4_NEON_Organisms/figs')

# This code for ggplot2 sets the theme to mostly black and white 
# (Arial font, and large font, base size=24)
theme_set(theme_bw(12))
theme_update(axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             strip.background = element_blank())

# ---------------------------------------------
# Define cumulative richness plotting functions
# ---------------------------------------------

# Function to reshape output of richness_cumulative for plotting with ggplot2

richness_shapeplotdat <- function(rich) {
  rich %>%
    melt(id.vars = c('siteID', 'year')) %>%
    mutate(type = case_when(grepl('chao1', variable) ~ 'chao1',
                            grepl('asymp', variable) ~ 'asymp',
                            TRUE ~ 'observed'),
           stat = case_when(grepl('var|stderr', variable) ~ 'var',
                            grepl('min', variable) ~ 'CImin',
                            grepl('max', variable) ~ 'CImax',
                            TRUE ~ 'estimate')) %>%
    dcast(siteID + year + type ~ stat) %>%
    group_by(siteID) %>%
    mutate(n_year = length(unique(year))) %>%
    ungroup
}

# Function to make a plot of cumulative richness by site and year
# Modified 02 Aug: add option to plot only some of the estimators

richness_plot <- function(dat, min_n_years, title, legend_pos, y_max, which_lines = c('observed', 'asymp', 'chao1')) {
  pd <- position_dodge(width = 0.2)
  
  # Modified 02 Aug: here add something to reorder the sites by the maximum observed richness
  rich_order <- dat %>% filter(type == 'observed') %>% arrange(estimate) %>% select(siteID) %>% unique %>% t %>% c
  dat <- dat %>% mutate(siteID = factor(siteID, levels = rich_order))
  
  dat %>%
    filter(n_year >= min_n_years, type %in% which_lines) %>%
    ggplot(aes(x = year, color = type)) +
    facet_wrap(~ siteID) +
    geom_errorbar(aes(ymin = CImin, ymax = CImax), width = 0.5, position = pd) +
    geom_line(aes(y = estimate), position = pd) +
    geom_point(aes(y = estimate), position = pd) +
    ggtitle(title) +
    scale_y_continuous(expand = c(0,0), name = 'richness') + 
    scale_x_continuous(breaks=min(dat$year):max(dat$year)) +
    theme(legend.position = legend_pos, 
          axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    coord_cartesian(ylim = c(0, y_max))
}

# -----------------------------------------------------------------
# Make plots of the cumulative richness estimates through the years
# -----------------------------------------------------------------

# Load data
richness_cumulative <- read.csv(file.path(data_path, 'alltaxa_richness_cumulative.csv'), stringsAsFactors = FALSE)

# Modification 02 Aug: Load NEON spatial data so that we can add the state abbreviation to the site IDs
site_states <- read.csv(file.path(data_path, 'spatial_data_csvs/All_Neon_TOS_Centroid_V4.csv'), stringsAsFactors = FALSE) %>%
  select(siteID, state) %>%
  unique %>%
  mutate(label = paste0(siteID, ' (', trimws(state), ')'))

all_labels <- site_states$label[match(richness_cumulative$siteID, site_states$siteID)]
richness_cumulative <- richness_cumulative %>%
  mutate(siteID = if_else(!is.na(all_labels), all_labels, siteID))

# Load background richness by NEON site for mammals, birds, and trees included in Little's range maps
mammal_site_background <- read.csv(file.path(data_path, 'IUCN_mammal_by_NEON_site.csv'), stringsAsFactors = FALSE)
bird_site_background <- read.csv(file.path(data_path, 'BOTW_bird_by_NEON_site.csv'), stringsAsFactors = FALSE)
tree_site_background <- read.csv(file.path(data_path, 'Little_tree_by_NEON_site.csv'), stringsAsFactors = FALSE)

# Calculate background richness by NEON site from the background species by site matrices
mammal_background_totals <- data.frame(siteID = names(mammal_site_background)[-(1:2)],
                                       background_richness = colSums(mammal_site_background[,-(1:2)])) %>%
  left_join(site_states) %>%
  mutate(siteID = if_else(!is.na(label), label, siteID))
bird_background_totals <- data.frame(siteID = names(bird_site_background)[-(1:3)],
                                       background_richness = colSums(bird_site_background[,-(1:3)])) %>%
  left_join(site_states) %>%
  mutate(siteID = if_else(!is.na(label), label, siteID))
tree_background_totals <- data.frame(siteID = names(tree_site_background)[-(1:2)],
                                       background_richness = colSums(tree_site_background[,-(1:2)])) %>%
  left_join(site_states) %>%
  mutate(siteID = if_else(!is.na(label), label, siteID))


# Reshape data for better plotting
# Only plot sites that have at least 3 years of samples
# For some taxa, might need to do fewer years

# Add a horizontal dotted line to show background richness for the taxa where we have it (mammal, bird, tree)

# Set option of which ones to plot
# Add 'asymp' as well if you want to add it back in
lines_to_plot <- c('observed','chao1')

p_mam <- richness_cumulative %>%
  filter(taxon == 'mammal') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  left_join(mammal_background_totals) %>%
  richness_plot(min_n_years = 3, title = 'Cumulative small mammal richness', legend_pos = c(0.9, 0.05), y_max = 27, which_lines = lines_to_plot) +
    geom_hline(aes(yintercept = background_richness), linetype = 'dashed', size = 1)
p_bird <- richness_cumulative %>%
  filter(taxon == 'bird') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  left_join(bird_background_totals) %>%
  richness_plot(min_n_years = 3, title = 'Cumulative bird richness', legend_pos = 'bottom', y_max = 150, which_lines = lines_to_plot) +
    geom_hline(aes(yintercept = background_richness), linetype = 'dashed', size = 1)
p_tree_little <- richness_cumulative %>%
  filter(taxon == 'trees') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  left_join(tree_background_totals) %>%
  richness_plot(min_n_years = 3, title = 'Cumulative tree richness (species included in Little)', legend_pos = 'bottom', y_max = 115, which_lines = lines_to_plot) +
    geom_hline(aes(yintercept = background_richness), linetype = 'dashed', size = 1)

p_aqplant <- richness_cumulative %>%
  filter(taxon == 'aquatic plants') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 3, title = 'Cumulative aquatic plant richness', legend_pos = c(0.9, 0.05), y_max = 50, which_lines = lines_to_plot)
p_macroinv <- richness_cumulative %>%
  filter(taxon == 'macroinvertebrate') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 3, title = 'Cumulative benthic macroinvertebrate richness', legend_pos = c(0.9, 0.05), y_max = 300, which_lines = lines_to_plot)
p_fish <- richness_cumulative %>%
  filter(taxon == 'fish') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 2, title = 'Cumulative fish richness', legend_pos = c(0.9, 0.05), y_max = 40, which_lines = lines_to_plot)
p_mosq <- richness_cumulative %>%
  filter(taxon == 'mosquito') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 3, title = 'Cumulative mosquito richness', legend_pos = c(0.9, 0.05), y_max = 75, which_lines = lines_to_plot)
p_plant <- richness_cumulative %>%
  filter(taxon == 'all plants') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 3, title = 'Cumulative terrestrial plant richness', legend_pos = c(0.9, 0.05), y_max = 1000, which_lines = lines_to_plot)
p_tree <- richness_cumulative %>%
  filter(taxon == 'all woody plants') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 3, title = 'Cumulative tree & woody plant richness', legend_pos = c(0.9, 0.05), y_max = 100, which_lines = lines_to_plot)
p_beetle <- richness_cumulative %>%
  filter(taxon == 'beetle') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 3, title = 'Cumulative ground beetle richness', legend_pos = c(0.9, 0.05), y_max = 200, which_lines = lines_to_plot)
p_zoop <- richness_cumulative %>%
  filter(taxon == 'zooplankton') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 3, title = 'Cumulative zooplankton richness', legend_pos = c(0.9, 0.05), y_max = 100, which_lines = lines_to_plot)
p_phytop <- richness_cumulative %>%
  filter(taxon == 'phytoplankton') %>% select(-taxon) %>%
  richness_shapeplotdat() %>%
  richness_plot(min_n_years = 3, title = 'Cumulative phytoplankton richness', legend_pos = c(0.9, 0.05), y_max = 250, which_lines = lines_to_plot)


pdf(file.path(fig_path, 'cumulative_richness_alltaxa.pdf'), height = 7, width = 9)
  p_mam; p_bird; p_tree_little
  p_aqplant; p_macroinv; p_fish; p_mosq; p_plant; p_tree; p_beetle; p_zoop; p_phytop
dev.off()
pdf(file.path(fig_path, 'cumulative_richness_alltaxa_withoutasymp.pdf'), height = 7, width = 9)
  p_mam; p_bird; p_tree_little
  p_aqplant; p_macroinv; p_fish; p_mosq; p_plant; p_tree; p_beetle; p_zoop; p_phytop
dev.off()

# -------------------------------------------------------------------------
# Functional/phylo diversity plots for birds and (Little) trees and mammals
# -------------------------------------------------------------------------

# Modified 02 Aug: add the state names to the labels and order by maximum functional diversity

# Load data
bird_pdfd <- read.csv(file.path(data_path, 'bird_pdfd_cumulative.csv'), stringsAsFactors = FALSE) %>%
  left_join(site_states) %>%
  mutate(siteID = if_else(!is.na(label), label, siteID)) %>%
  select(-state, -label)
tree_pdfd <- read.csv(file.path(data_path, 'tree_pdfd_cumulative.csv'), stringsAsFactors = FALSE) %>%
  left_join(site_states) %>%
  mutate(siteID = if_else(!is.na(label), label, siteID)) %>%
  select(-state, -label)
mammal_pdfd <- read.csv(file.path(data_path, 'mammal_pdfd_cumulative.csv'), stringsAsFactors = FALSE) %>%
  left_join(site_states) %>%
  mutate(siteID = if_else(!is.na(label), label, siteID)) %>%
  select(-state, -label)

# Pick 2 colors for plot
pd_cols <- c('goldenrod','forestgreen')

# Bird pd and fd

pdat_phyfun_bird <- bird_pdfd %>%
  melt(id.vars = c('siteID', 'year')) %>%
  group_by(siteID) %>%
  mutate(n_year = length(unique(year))) %>%
  ungroup %>%
  filter(n_year >= 3, variable %in% c('MPD_phy_z', 'MPD_func_z'))

func_order <- pdat_phyfun_bird %>% filter(variable == 'MPD_func_z') %>% arrange(value) %>% select(siteID) %>% unique %>% t %>% c
pdat_phyfun_bird <- pdat_phyfun_bird %>% mutate(siteID = factor(siteID, levels = func_order))

p_phyfun_bird <- ggplot(pdat_phyfun_bird, aes(x = year, y = value, color = variable)) +
  facet_wrap(~ siteID) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  scale_color_manual(name = 'diversity type', labels = c('functional', 'phylogenetic'), values = pd_cols) +
  ggtitle('Bird cumulative functional and phylogenetic diversity') +
  scale_y_continuous(expand = c(0,0), name = 'PD and FD as effect size', limits = c(-10, 0)) + 
  scale_x_continuous(breaks=min(bird_pdfd$year):max(bird_pdfd$year)) +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 90, vjust = 0.5))

# Tree pd and fd

pdat_phyfun_tree <- tree_pdfd %>%
  melt(id.vars = c('siteID', 'year')) %>%
  group_by(siteID) %>%
  mutate(n_year = length(unique(year))) %>%
  ungroup %>%
  filter(n_year >= 3, variable %in% c('MPD_phy_z', 'MPD_func_z'), !siteID %in% c('ONAQ')) # ONAQ has no valid diversity values

func_order <- pdat_phyfun_tree %>% filter(variable == 'MPD_func_z') %>% arrange(value) %>% select(siteID) %>% unique %>% t %>% c
pdat_phyfun_tree <- pdat_phyfun_tree %>% mutate(siteID = factor(siteID, levels = func_order))

p_phyfun_tree <- ggplot(pdat_phyfun_tree, aes(x = year, y = value, color = variable)) +
  facet_wrap(~ siteID) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  scale_color_manual(name = 'diversity type', labels = c('functional', 'phylogenetic'), values = pd_cols) +
  ggtitle('Tree cumulative functional and phylogenetic diversity', 'only species on Little\'s maps') +
  scale_y_continuous(expand = c(0,0), name = 'PD and FD as effect size', limits = c(-5, 5)) + 
  scale_x_continuous(breaks=min(tree_pdfd$year):max(tree_pdfd$year)) +
  theme(legend.position = c(0.9, 0.05), axis.text.x = element_text(angle = 90, vjust = 0.5))

# Mammal pd and fd

pdat_phyfun_mammal <- mammal_pdfd %>%
  melt(id.vars = c('siteID', 'year')) %>%
  group_by(siteID) %>%
  mutate(n_year = length(unique(year))) %>%
  ungroup %>%
  filter(n_year >= 3, variable %in% c('MPD_phy_z', 'MPD_func_z'))

func_order <- pdat_phyfun_mammal %>% filter(variable == 'MPD_func_z') %>% arrange(value) %>% select(siteID) %>% unique %>% t %>% c
pdat_phyfun_mammal <- pdat_phyfun_mammal %>% mutate(siteID = factor(siteID, levels = func_order))

p_phyfun_mammal <- ggplot(pdat_phyfun_mammal, aes(x = year, y = value, color = variable)) +
  facet_wrap(~ siteID) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  scale_color_manual(name = 'diversity type', labels = c('functional', 'phylogenetic'), values = pd_cols) +
  ggtitle('Mammal cumulative functional and phylogenetic diversity') +
  scale_y_continuous(expand = c(0,0), name = 'PD and FD as effect size', limits = c(-5, 5)) + 
  scale_x_continuous(breaks=min(tree_pdfd$year):max(tree_pdfd$year)) +
  theme(legend.position = c(0.92, 0.05), axis.text.x = element_text(angle = 90, vjust = 0.5))

pdf(file.path(fig_path, 'cumulative_PDFD_mammal_bird_tree.pdf'), height = 7, width = 9)
p_phyfun_mammal
p_phyfun_bird
p_phyfun_tree
dev.off()


# -------------------------------------------------------------------------  
# Make maps of local richness at NEON sites superimposed on background grid
# -------------------------------------------------------------------------

# Get site locations from NEON API

library(raster)
library(httr)
library(jsonlite)

site_locs <- GET('http://data.neonscience.org/api/v0/locations/sites') %>%
  content(as = 'text') %>%
  fromJSON(simplifyDataFrame = TRUE, flatten = TRUE) 
site_locs <- site_locs$data[,1:19] %>%
  filter(!is.na(locationDecimalLatitude) & !is.na(locationDecimalLongitude))

site_coords <- site_locs %>%
  dplyr::select(locationName, locationDecimalLongitude, locationDecimalLatitude, locationElevation) %>%
  setNames(c('siteID', 'lon', 'lat', 'elevation'))

# Remove the state labels from the siteIDs again so that things will match.
richness_cumulative <- richness_cumulative %>%
  mutate(siteID = substr(siteID, 1, 4))
mammal_background_totals <- mammal_background_totals %>%
  mutate(siteID = substr(siteID, 1, 4))
bird_background_totals <- bird_background_totals %>%
  mutate(siteID = substr(siteID, 1, 4))
tree_background_totals <- tree_background_totals %>%
  mutate(siteID = substr(siteID, 1, 4))
mammal_pdfd <- mammal_pdfd %>%
  mutate(siteID = substr(siteID, 1, 4))
bird_pdfd <- bird_pdfd %>%
  mutate(siteID = substr(siteID, 1, 4))
tree_pdfd <- tree_pdfd %>%
  mutate(siteID = substr(siteID, 1, 4))


# Function for drawing a map
# Set color scale manually so that the same colors can be used for US, AK, and PR
richness_map <- function(tax_group, background, legend_name = NA, legend_breaks = NA, lats, lons, scale_max, legend_pos = 'bottom', show_grid = TRUE) {
  
  # Get observed and asymptotic richness from final year for mammals
  rich_mapdat <- richness_cumulative %>%
    filter(taxon == tax_group) %>% dplyr::select(-taxon) %>%
    richness_shapeplotdat() %>%
    group_by(siteID) %>%
    filter(year == max(year)) %>%
    ungroup %>%
    left_join(site_coords)
  
  # Convert mammal raster to format for ggplot2
  if (show_grid) {
    background_df <- background %>%
      as('SpatialPixelsDataFrame') %>%
      as.data.frame %>%
      setNames(c('richness', 'lon', 'lat'))
  }
  
  # Scale the observed richness by the background richness
  rich_mapdat$scaled_estimate <- rich_mapdat$estimate * scale_max/max(rich_mapdat$estimate, na.rm = TRUE)
  
  rybcolors <- rev(RColorBrewer::brewer.pal(9, 'RdYlBu'))
  
  if (show_grid) {
    p <- ggplot(background_df %>% filter(richness > 0), aes(x=lon, y=lat)) + geom_tile(aes(fill = richness))
  } else {
    p <- ggplot(rich_mapdat %>% filter(type == 'observed'), aes(x=lon, y=lat))
  }
  
  # To do 26 June: replace borders() with manually downloaded high-resolution country borders.
  p +
    borders('state') + borders('world', xlim = lons, ylim = lats) +
    geom_point(data = rich_mapdat %>% filter(type == 'observed'), aes(fill = scaled_estimate), size = 5, pch = 21) +
    geom_text(data = rich_mapdat %>% filter(type == 'observed'), aes(label = estimate), size = 3) +
    scale_fill_gradientn(colors = rybcolors, name = legend_name, breaks = legend_breaks, limits = c(0, scale_max)) +
    coord_map(xlim = lons, ylim = lats) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = legend_pos)
  
}

# Load rasters
# Lower 48
mam_raster <- raster(file.path(data_path, 'IUCN_mammal_grid.tif'))
bird_raster <- raster(file.path(data_path, 'BOTW_bird_grid.tif'))
tree_raster <- raster(file.path(data_path, 'Little_tree_grid.tif'))
# Alaska
akmam_raster <- raster(file.path(data_path, 'IUCN_mammal_grid_AK.tif'))
akbird_raster <- raster(file.path(data_path, 'BOTW_bird_grid_AK.tif'))
aktree_raster <- raster(file.path(data_path, 'Little_tree_grid_AK.tif'))
# Puerto Rico (birds only): Little trees are not included for PR, and there are no native terrestrial mammals in PR
# The only thing you would get there are mice and rats and mongeese, so they did not trap mammals there
prbird_raster <- raster(file.path(data_path, 'BOTW_bird_grid_PR.tif'))

# Latitude and longitude limits of maps
latus <- c(25, 50)
lonus <- c(-125, -67)
latpr <- c(17.8, 18.6)
lonpr <- c(-68, -65)
latak <- c(51, 72)
lonak <- c(-179.9, -129)

# Maximum values of richness to use for scale
mam_max <- cellStats(mam_raster, stat = 'max')
bird_max <- cellStats(bird_raster, stat = 'max')
tree_max <- cellStats(tree_raster, stat = 'max')

# Draw maps
# Lower 48
map_mam <- richness_map(tax_group = 'mammal', 
                        background = mam_raster, 
                        legend_name = 'Small mammal richness', 
                        legend_breaks = c(0,5,10,15,20,25), 
                        lats = latus, 
                        lons = lonus,
                        scale_max = mam_max)
map_bird <- richness_map(tax_group = 'bird', 
                         background = bird_raster, 
                         legend_name = 'Bird richness', 
                         legend_breaks = c(0,100,200,300), 
                         lats = latus, 
                         lons = lonus,
                         scale_max = bird_max)
map_tree <- richness_map(tax_group = 'trees', 
                         background = tree_raster, 
                         legend_name = 'Tree richness', 
                         legend_breaks = c(0,50,100,150,200), 
                         lats = latus, 
                         lons = lonus,
                         scale_max = tree_max)
# Alaska
akmap_mam <- richness_map(tax_group = 'mammal', 
                        background = akmam_raster, 
                        lats = latak, 
                        lons = lonak,
                        scale_max = mam_max,
                        legend_pos = 'none')
akmap_bird <- richness_map(tax_group = 'bird', 
                         background = akbird_raster, 
                         lats = latak, 
                         lons = lonak,
                         scale_max = bird_max,
                         legend_pos = 'none')
akmap_tree <- richness_map(tax_group = 'trees', 
                         background = aktree_raster, 
                         lats = latak, 
                         lons = lonak,
                         scale_max = tree_max,
                         legend_pos = 'none')
# Puerto Rico (only bird map is used)
prmap_bird <- richness_map(tax_group = 'bird', 
                         background = prbird_raster, 
                         lats = latpr, 
                         lons = lonpr,
                         scale_max = bird_max,
                         legend_pos = 'none')
prmap_tree <- richness_map(tax_group = 'trees', 
                         background = NA, 
                         lats = latpr, 
                         lons = lonpr,
                         scale_max = tree_max,
                         legend_pos = 'none',
                         show_grid = FALSE) # No raster grid for PR trees 

# Observed vs regional pool plots (all regions)
# Get final observed richness for each site
rich_plotdat <- richness_cumulative %>%
  filter(taxon %in% c('mammal','bird','trees')) %>%
  group_by(siteID) %>%
  filter(year == max(year)) %>%
  ungroup %>%
  dplyr::select(taxon, siteID, richness)

scatter_mam <- rich_plotdat %>%
  filter(taxon=='mammal') %>%
  left_join(mammal_background_totals) %>%
  ggplot(aes(x = background_richness, y = richness)) +
    geom_point() +
    labs(x = 'regional species pool', y = 'observed richness')
scatter_bird <- rich_plotdat %>%
  filter(taxon=='bird') %>%
  left_join(bird_background_totals) %>%
  ggplot(aes(x = background_richness, y = richness)) +
    geom_point() +
    labs(x = 'regional species pool', y = 'observed richness')
scatter_tree <- rich_plotdat %>%
  filter(taxon=='trees', !siteID %in% c('GUAN')) %>% # Get rid of PR for tree plot
  left_join(tree_background_totals) %>%
  ggplot(aes(x = background_richness, y = richness)) +
    geom_point() +
    labs(x = 'regional species pool', y = 'observed richness')

# Save maps
library(cowplot)

# Arrange plots with US map on top, AK in bottom left, and scatter plot on top of PR in bottom right
# Note: the gridding makes each of these calls take a long time to run (maybe >30min total)
# The relative heights and widths can easily be tweaked to change layout of grid

# Bottom right: scatter plot and PR inset map
# Placeholder of a blank ggplot() is put where PR map would go
mam_bottomright <- plot_grid(scatter_mam, ggplot() + panel_border(remove = TRUE), nrow = 2, rel_heights = c(1, 0.5)) # PR has no mammal map
bird_bottomright <- plot_grid(scatter_bird, prmap_bird + panel_border(colour = 'black'), nrow = 2, rel_heights = c(1, 0.5))
tree_bottomright <- plot_grid(scatter_tree, ggplot() + panel_border(remove = TRUE), nrow = 2, rel_heights = c(1, 0.5)) # Don't plot PR tree map either

# Add bottom left: AK inset map
mam_bottomrow <- plot_grid(akmap_mam + panel_border(colour = 'black'), mam_bottomright, nrow = 1, rel_widths = c(1.1, 1)) 
bird_bottomrow <- plot_grid(akmap_bird + panel_border(colour = 'black'), bird_bottomright, nrow = 1, rel_widths = c(1.1, 1))
tree_bottomrow <- plot_grid(akmap_tree + panel_border(colour = 'black'), tree_bottomright, nrow = 1, rel_widths = c(1.1, 1))

# Entire plot: US map on top of bottom row
mam_plot <- plot_grid(map_mam + panel_border(colour = 'black') + theme(legend.key.width = unit(0.3, 'inches'), legend.justification = 'center') , mam_bottomrow, nrow = 2, rel_heights = c(1, 0.6))
bird_plot <- plot_grid(map_bird + panel_border(colour = 'black') + theme(legend.justification = 'center'), bird_bottomrow, nrow = 2, rel_heights = c(1, 0.6))
tree_plot <- plot_grid(map_tree + panel_border(colour = 'black') + theme(legend.justification = 'center'), tree_bottomrow, nrow = 2, rel_heights = c(1, 0.6))


ggsave(file.path(fig_path, 'paneled_mammal_background_map.png'), mam_plot, height = 9, width = 9, dpi = 400)
ggsave(file.path(fig_path, 'paneled_bird_background_map.png'), bird_plot, height = 9, width = 9, dpi = 400)
ggsave(file.path(fig_path, 'paneled_tree_background_map.png'), tree_plot, height = 9, width = 9, dpi = 400)


# -------------------------------------------------------
# Added 1 July: maps for taxa without background richness
# -------------------------------------------------------

taxa_to_plot <- c("beetle", "mosquito", "all plants", "macroinvertebrate", "fish", "zooplankton", "phytoplankton", "aquatic plants")

map_plots <- list()

for (grp in taxa_to_plot) {
  
  # Maximum value for scales
  grp_max <- max(richness_cumulative$richness[richness_cumulative$taxon == grp], na.rm = TRUE) 
  # Breaks for scales
  scale_breaks <- seq(0, grp_max, 10^(floor(log10(grp_max))))
  
  #USA
  
  usmap <- richness_map(tax_group = grp, 
                        background = NA, 
                        legend_name = paste(grp, 'richness'), 
                        legend_breaks = scale_breaks,
                        lats = latus, 
                        lons = lonus,
                        scale_max = grp_max,
                        show_grid = FALSE)
  
  # Alaska
  akmap <- richness_map(tax_group = grp, 
                        background = NA, 
                        lats = latak, 
                        lons = lonak,
                        scale_max = grp_max,
                        legend_pos = 'none',
                        show_grid = FALSE)
  
  # Puerto Rico 
  prmap <- richness_map(tax_group = grp, 
                        background = NA, 
                        lats = latpr, 
                        lons = lonpr,
                        scale_max = grp_max,
                        legend_pos = 'none',
                        show_grid = FALSE)
  
  map_bottomrow <- plot_grid(akmap + panel_border(colour = 'black'), 
                             prmap + panel_border(colour = 'black'), 
                             nrow = 1, rel_widths = c(1.1, 0.9))
  map_plots[[length(map_plots) + 1]] <- plot_grid(usmap + 
                                                    panel_border(colour = 'black') + 
                                                    theme(legend.justification = 'center',
                                                          legend.key.width = unit(0.4, 'inches')), 
                                                  map_bottomrow, 
                                                  nrow = 2, rel_heights = c(1, 0.6))
}

pdf(file.path(fig_path, 'maps_richness_othertaxa.pdf'), height = 9, width = 9)
  sapply(1:8, function(i) print(map_plots[[i]]))
dev.off()

# --------------------------------------------
# Map of phylogenetic and functional diversity
# --------------------------------------------

# This section added 2 July

# Load background PD and FD by grid
mam_pdfd_background_lower48 <- read.csv(file.path(data_path, 'mammal_background_PDFD_lower48.csv'), stringsAsFactors = FALSE)
mam_pdfd_background_AK <- read.csv(file.path(data_path, 'mammal_background_PDFD_AK.csv'), stringsAsFactors = FALSE)
tree_pdfd_background_lower48 <- read.csv(file.path(data_path, 'littletree_background_PDFD_lower48.csv'), stringsAsFactors = FALSE)
tree_pdfd_background_AK <- read.csv(file.path(data_path, 'littletree_background_PDFD_AK.csv'), stringsAsFactors = FALSE)
bird_pdfd_background_lower48 <- read.csv(file.path(data_path, 'bird_background_PDFD_lower48.csv'), stringsAsFactors = FALSE)
bird_pdfd_background_AK <- read.csv(file.path(data_path, 'bird_background_PDFD_AK.csv'), stringsAsFactors = FALSE)

# Load background PD and FD by site (for scatter plot)
mammal_pdfd_site_background <- read.csv(file.path(data_path, 'IUCN_mammal_PDFD_by_NEON_site.csv'), stringsAsFactors = FALSE)
bird_pdfd_site_background <- read.csv(file.path(data_path, 'BOTW_bird_PDFD_by_NEON_site.csv'), stringsAsFactors = FALSE)
tree_pdfd_site_background <- read.csv(file.path(data_path, 'Little_tree_PDFD_by_NEON_site.csv'), stringsAsFactors = FALSE)

# Load observed PD and FD across all years (non-cumulative) by site
mammal_pdfd_total <- read.csv(file.path(data_path, 'mammal_pdfd_total.csv'), stringsAsFactors = FALSE)
bird_pdfd_total <- read.csv(file.path(data_path, 'bird_pdfd_total.csv'), stringsAsFactors = FALSE)
tree_pdfd_total <- read.csv(file.path(data_path, 'tree_pdfd_total.csv'), stringsAsFactors = FALSE)

# (Load everything else above)

# Function for making PD or FD map
# --------------------------------

pdfd_map <- function(site_values, background, column_name, legend_name = NA, legend_breaks = NA, lats, lons, legend_pos = 'bottom', show_grid = TRUE) {
  
  # Get PD or FD for the sites
  site_mapdat <- site_values %>%
    rename(diversity = !!column_name) %>%
    filter(!is.na(diversity)) %>%
    group_by(siteID) %>%
    dplyr::select(siteID, diversity) %>%
    ungroup %>%
    left_join(site_coords)
  
  background <- background %>%
    rename(diversity = !!column_name) %>%
    dplyr::select(lon, lat, diversity)
  
  # Do not do scaling here because of the negative values.

  rybcolors <- rev(RColorBrewer::brewer.pal(9, 'RdYlBu'))
  
  p <- ggplot(background, aes(x = lon, y = lat))
  
  if (show_grid) p <- p + geom_tile(aes(fill = diversity))

  p +
    borders('state') + borders('world', xlim = lons, ylim = lats) +
    geom_point(data = site_mapdat, aes(fill = diversity), size = 5, pch = 21) +
    scale_fill_gradientn(colors = rybcolors, name = legend_name, breaks = legend_breaks, na.value = 'transparent') +
    coord_map(xlim = lons, ylim = lats) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = legend_pos)
  
}

# Create map objects
# ------------------

# Phylogenetic diversity maps
# Lower 48
map_mam_pd_48 <- pdfd_map(site_values = mammal_pdfd_total, 
                          background = mam_pdfd_background_lower48,
                          column_name = 'MPD_phy_z', 
                          legend_name = 'Small mammal\nphylogenetic diversity', 
                          legend_breaks = (-3):2, 
                          lats = latus, 
                          lons = lonus, 
                          legend_pos = 'bottom')
map_bird_pd_48 <- pdfd_map(site_values = bird_pdfd_total, 
                          background = bird_pdfd_background_lower48,
                          column_name = 'MPD_phy_z', 
                          legend_name = 'Bird phylogenetic\ndiversity', 
                          legend_breaks = seq(-20, 0, by = 4), 
                          lats = latus, 
                          lons = lonus, 
                          legend_pos = 'bottom')
map_tree_pd_48 <- pdfd_map(site_values = tree_pdfd, 
                          background = tree_pdfd_background_lower48,
                          column_name = 'MPD_phy_z', 
                          legend_name = 'Tree phylogenetic\ndiversity', 
                          legend_breaks = seq(-6, 2, by = 2), 
                          lats = latus, 
                          lons = lonus, 
                          legend_pos = 'bottom')

# Alaska
map_mam_pd_AK <- pdfd_map(site_values = mammal_pdfd, 
                          background = mam_pdfd_background_AK,
                          column_name = 'MPD_phy_z', 
                          lats = latak, 
                          lons = lonak, 
                          legend_pos = 'none')
map_bird_pd_AK <- pdfd_map(site_values = bird_pdfd, 
                          background = bird_pdfd_background_AK,
                          column_name = 'MPD_phy_z', 
                          lats = latak, 
                          lons = lonak, 
                          legend_pos = 'none')
map_tree_pd_AK <- pdfd_map(site_values = tree_pdfd, 
                           background = tree_pdfd_background_AK,
                           column_name = 'MPD_phy_z', 
                           lats = latak, 
                           lons = lonak, 
                           legend_pos = 'none')

# Scatter plots
scatter_mam_pd <- mammal_pdfd_total %>%
  left_join(mammal_pdfd_site_background %>% setNames(c('siteID', paste0('background_', names(.)[-1])))) %>%
  ggplot(aes(x = background_MPD_phy_z, y = MPD_phy_z)) +
  geom_point() +
  labs(x = 'regional species pool PD', y = 'observed PD')
scatter_bird_pd <- bird_pdfd_total %>%
  left_join(bird_pdfd_site_background %>% setNames(c('siteID', paste0('background_', names(.)[-1])))) %>%
  ggplot(aes(x = background_MPD_phy_z, y = MPD_phy_z)) +
  geom_point() +
  labs(x = 'regional species pool PD', y = 'observed PD')
scatter_tree_pd <- tree_pdfd_total %>%
  left_join(tree_pdfd_site_background %>% setNames(c('siteID', paste0('background_', names(.)[-1])))) %>%
  ggplot(aes(x = background_MPD_phy_z, y = MPD_phy_z)) +
  geom_point() +
  labs(x = 'regional species pool PD', y = 'observed PD')

# Arrange into panels and save

library(cowplot)

# Bottom row: AK map and scatterplot
mam_bottomrow <- plot_grid(map_mam_pd_AK + panel_border(colour = 'black'), scatter_mam_pd, nrow = 1, rel_widths = c(1.1, 1)) 
bird_bottomrow <- plot_grid(map_bird_pd_AK + panel_border(colour = 'black'), scatter_bird_pd, nrow = 1, rel_widths = c(1.1, 1))
tree_bottomrow <- plot_grid(map_tree_pd_AK + panel_border(colour = 'black'), scatter_tree_pd, nrow = 1, rel_widths = c(1.1, 1))

# Entire plot: US map on top of bottom row
mam_plot <- plot_grid(map_mam_pd_48 + panel_border(colour = 'black') + theme(legend.key.width = unit(0.4, 'inches'), legend.justification = 'center') , mam_bottomrow, nrow = 2, rel_heights = c(1, 0.6))
bird_plot <- plot_grid(map_bird_pd_48 + panel_border(colour = 'black') + theme(legend.key.width = unit(0.4, 'inches'), legend.justification = 'center'), bird_bottomrow, nrow = 2, rel_heights = c(1, 0.6))
tree_plot <- plot_grid(map_tree_pd_48 + panel_border(colour = 'black') + theme(legend.key.width = unit(0.4, 'inches'), legend.justification = 'center'), tree_bottomrow, nrow = 2, rel_heights = c(1, 0.6))


ggsave(file.path(fig_path, 'paneled_mammal_background_map_PD.png'), mam_plot, height = 9, width = 9, dpi = 400)
ggsave(file.path(fig_path, 'paneled_bird_background_map_PD.png'), bird_plot, height = 9, width = 9, dpi = 400)
ggsave(file.path(fig_path, 'paneled_tree_background_map_PD.png'), tree_plot, height = 9, width = 9, dpi = 400)

# Functional diversity maps
# Lower 48
map_mam_fd_48 <- pdfd_map(site_values = mammal_pdfd_total, 
                          background = mam_pdfd_background_lower48,
                          column_name = 'MPD_func_z', 
                          legend_name = 'Small mammal\nfunctional diversity', 
                          legend_breaks = (-2):2, 
                          lats = latus, 
                          lons = lonus, 
                          legend_pos = 'bottom')
map_bird_fd_48 <- pdfd_map(site_values = bird_pdfd_total, 
                           background = bird_pdfd_background_lower48,
                           column_name = 'MPD_func_z', 
                           legend_name = 'Bird functional\ndiversity', 
                           legend_breaks = seq(-10,0,by=2), 
                           lats = latus, 
                           lons = lonus, 
                           legend_pos = 'bottom')
map_tree_fd_48 <- pdfd_map(site_values = tree_pdfd, 
                           background = tree_pdfd_background_lower48,
                           column_name = 'MPD_func_z', 
                           legend_name = 'Tree functional\ndiversity', 
                           legend_breaks = (-3):3, 
                           lats = latus, 
                           lons = lonus, 
                           legend_pos = 'bottom')

# Alaska
map_mam_fd_AK <- pdfd_map(site_values = mammal_pdfd, 
                          background = mam_pdfd_background_AK,
                          column_name = 'MPD_func_z', 
                          lats = latak, 
                          lons = lonak, 
                          legend_pos = 'none')
map_bird_fd_AK <- pdfd_map(site_values = bird_pdfd, 
                           background = bird_pdfd_background_AK,
                           column_name = 'MPD_func_z', 
                           lats = latak, 
                           lons = lonak, 
                           legend_pos = 'none')
map_tree_fd_AK <- pdfd_map(site_values = tree_pdfd, 
                           background = tree_pdfd_background_AK,
                           column_name = 'MPD_func_z', 
                           lats = latak, 
                           lons = lonak, 
                           legend_pos = 'none')

# Scatter plots
scatter_mam_fd <- mammal_pdfd_total %>%
  left_join(mammal_pdfd_site_background %>% setNames(c('siteID', paste0('background_', names(.)[-1])))) %>%
  ggplot(aes(x = background_MPD_func_z, y = MPD_func_z)) +
  geom_point() +
  labs(x = 'regional species pool FD', y = 'observed FD')
scatter_bird_fd <- bird_pdfd_total %>%
  left_join(bird_pdfd_site_background %>% setNames(c('siteID', paste0('background_', names(.)[-1])))) %>%
  ggplot(aes(x = background_MPD_func_z, y = MPD_func_z)) +
  geom_point() +
  labs(x = 'regional species pool FD', y = 'observed FD')
scatter_tree_fd <- tree_pdfd_total %>%
  left_join(tree_pdfd_site_background %>% setNames(c('siteID', paste0('background_', names(.)[-1])))) %>%
  ggplot(aes(x = background_MPD_func_z, y = MPD_func_z)) +
  geom_point() +
  labs(x = 'regional species pool FD', y = 'observed FD')

# Arrange into panels and save

library(cowplot)

# Bottom row: AK map and scatterplot
mam_bottomrow <- plot_grid(map_mam_fd_AK + panel_border(colour = 'black'), scatter_mam_fd, nrow = 1, rel_widths = c(1.1, 1)) 
bird_bottomrow <- plot_grid(map_bird_fd_AK + panel_border(colour = 'black'), scatter_bird_fd, nrow = 1, rel_widths = c(1.1, 1))
tree_bottomrow <- plot_grid(map_tree_fd_AK + panel_border(colour = 'black'), scatter_tree_fd, nrow = 1, rel_widths = c(1.1, 1))

# Entire plot: US map on top of bottom row
mam_plot <- plot_grid(map_mam_fd_48 + panel_border(colour = 'black') + theme(legend.key.width = unit(0.4, 'inches'), legend.justification = 'center') , mam_bottomrow, nrow = 2, rel_heights = c(1, 0.6))
bird_plot <- plot_grid(map_bird_fd_48 + panel_border(colour = 'black') + theme(legend.key.width = unit(0.4, 'inches'), legend.justification = 'center'), bird_bottomrow, nrow = 2, rel_heights = c(1, 0.6))
tree_plot <- plot_grid(map_tree_fd_48 + panel_border(colour = 'black') + theme(legend.key.width = unit(0.4, 'inches'), legend.justification = 'center'), tree_bottomrow, nrow = 2, rel_heights = c(1, 0.6))


ggsave(file.path(fig_path, 'paneled_mammal_background_map_FD.png'), mam_plot, height = 9, width = 9, dpi = 400)
ggsave(file.path(fig_path, 'paneled_bird_background_map_FD.png'), bird_plot, height = 9, width = 9, dpi = 400)
ggsave(file.path(fig_path, 'paneled_tree_background_map_FD.png'), tree_plot, height = 9, width = 9, dpi = 400)

