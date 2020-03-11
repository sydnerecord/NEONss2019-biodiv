# Modified 26 June: add option to do by site or by plot and to keep either the final observed values or all years
# (group argument and by_year = TRUE or FALSE)
# Modified 29 June: add option for different plot name for the aquatic taxa (uses namedLocation)

#' @title cumulative richness estimators by site and year
#'
#' @description Function to get cumulative richness estimators by site and year. Sequentially add
#' years and see what happens to the cumulative observed richness and richness estimators. Each
#' year's result represents all data up to that point in time.
#'
#'
#' @import iNEXT dplyr tibble
#'
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
#' @export
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
