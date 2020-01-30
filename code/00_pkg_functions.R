# Packages used
library(tidyverse)
library(googledrive)
# devtools::install_github("EDIorg/ecocomDP")
library(ecocomDP)
library(neonUtilities)
library(lubridate)
library(iNEXT)

#' Write file to Google Drive
#' 
#' @param data_to_write Data frame you want to write out.
#' @param write_filename File name to be saved.
#' @param my_path_to_googledirve_directory Path in Google Drive directory.
#' @param keep_file_in_working_dir Do you want to keep file locally?
#' 
write_to_google_drive <- function(
  data_to_write,
  write_filename,
  my_path_to_googledirve_directory, 
  keep_file_in_working_dir = FALSE 
){
  # get list of files
  my_list_of_files <- googledrive::drive_ls(my_path_to_googledirve_directory)
  
  # make a new output filename
  # temp write local
  readr::write_csv(data_to_write, write_filename)
  
  # conditional depending on if we need to overwrite or create new
  if(!write_filename %in% my_list_of_files$name){
    drive_upload(write_filename,
                 path = my_path_to_googledirve_directory,
                 name = write_filename,
                 type = NULL,
                 verbose = TRUE)
    message(paste0('Created ',write_filename, ' in ', my_path_to_googledirve_directory))
  }else{
    google_id <- my_list_of_files %>% filter(name == write_filename) %>% select(id) %>% unlist()
    drive_update(file = as_id(google_id),
                 media = write_filename)
    message(paste0('Updated ',write_filename, ' in ', my_path_to_googledirve_directory))
  }
  
  #remove local file
  if(!keep_file_in_working_dir) file.remove(write_filename)
}

#' Chao1 richness estimator
#' 
#' @param x A vector of species IDs.
#' @return A data frame of Chao 1 estiamtion and variance and 95% CI.
#' 
estimator_chao1 <- function(x) {
  xcomm <- table(x)
  S_obs <- length(xcomm) # Number of species observed
  f1 <- sum(xcomm == 1) # Number of singletons
  f2 <- sum(xcomm == 2) # Number of doubletons
  chao1 <- S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1)) # Calculate chao1 estimator
  var_chao1 <- f2 * ( ((f1/f2)/4)^4 + (f1/f2)^3 + ((f1/f2)/2)^2 ) # Variance of estimator
  if (!is.finite(var_chao1)) var_chao1 <- 0 # If no doubletons, variance is zero
  data.frame(chao1 = chao1, 
             chao1_var = var_chao1,
             chao1_CImin = max(S_obs, chao1 - 1.96 * sqrt(var_chao1)),
             chao1_CImax = chao1 + 1.96 * sqrt(var_chao1)
             )
}


#' Get the asymptotic estimator for richness, and the bounds of its 95% conf int
#' 
#' @param x A vector of species IDs.
#' @param ... Other options for iNEXT::iNEXT().
#' @return A data frame with estimated asymptote of species richness, SE, and 95% CI.
#' 
estimator_asymp <- function(x, ...){
  xcomm <- table(x)
  inext_out <- iNEXT(x = list(as.numeric(xcomm)), datatype = "abundance", ...) # run iNEXT on the community
  filter(inext_out$AsyEst, Diversity == 'Species richness') %>% 
    select(asymp_est = Estimator, 
           asymp_est_stderr = s.e., 
           asymp_est_CImin = LCL, 
           asymp_est_CImax = UCL)
}

#' Function to get cumulative richness estimators by site and year
#' 
#' Sequentially add years and see what happens to the cumulative observed richness and richness estimators.
#' Each year's result represents all data up to that point in time.

# Modified 26 June: add option to do by site or by plot and to keep either the final observed values or all years
# (group argument and by_year = TRUE or FALSE)
# Modified 29 June: add option for different plot name for the aquatic taxa (uses namedLocation)

#' Calculate cumulative richness
#' 
#' @param dat The name of the data frame.
#' @param sp_name Column name of species/taxon, a single string.
#' @param site_name Column name of site, a single string.
#' @param plot_name Column name of plot (if any), a single string.
#' @param value_name Column name of abundance/cover/density/presence/absence, a single string.
#' @param year_name Column name of year, a single string.
#' @param grp_by Whether to group by site or group by both site and group.
#' @param cumul_by_year To calculate richness culumatively through sample years? Default is `TRUE`.
#' 
#' @return A data frame of richness, Chao1 etimates, asymptote estimates and their uncertainties.
#' 
richness_cumulative <- function(dat, sp_name = "taxonID", site_name = "siteID", 
                                plot_name = "plotID", value_name = "value",
                                year_name = "year",
                                grp_by = c("site", "site_plot"), 
                                cumul_by_year = TRUE) {
  grp_by = match.arg(grp_by)
  # select relevant columns and rename to a standardized set of names
  if(plot_name %in% names(dat)){
    dat = dat[, c(site_name, plot_name, year_name, sp_name, value_name)] %>% 
      setNames(c("site", "plot", "year", "sp", "abund"))
  } else {
    dat = dat[, c(site_name, year_name, sp_name, value_name)] %>% 
      setNames(c("site", "year", "sp", "abund"))
  }
  
  # remove abund == 0
  dat = filter(dat, abund > 0)
  
  if (grp_by == 'site') dat <- group_by(dat, site)
  if (grp_by == 'site_plot') dat <- group_by(dat, site, plot)
  
  div_f = function(x){ # x is a data frame
    bind_cols(
      tibble(richness = n_distinct(x$sp)),
      estimator_chao1(x$sp),
      estimator_asymp(x$sp)
    )
  }
  
  if(!cumul_by_year){ # use all year's data
    out = dat %>% do(div_f(.)) %>% ungroup()
  } else {# cumulative by year
    yrs = sort(unique(dat$year))
    out = map_dfr(yrs, function(yr){
      # cat("yr = ", yr, "\t")
      filter(dat, year <= yr) %>% do(div_f(.)) %>% ungroup() %>% 
        add_column(year = yr)
    })
  }
  
  out 
}
