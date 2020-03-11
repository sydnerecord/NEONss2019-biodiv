#' @title asymptotic estimator for richness
#'
#'
#' @description Get the asymptotic estimator for richness, and the bounds of its 95% conf int
#'
#'
#' @import iNEXT dplyr tibble
#'
#'
#' @param x A vector of species IDs.
#' @param ... Other options for iNEXT::iNEXT().
#'
#'
#' @return A data frame with estimated asymptote of species richness, SE, and 95% CI.
#'
#'
#' @export
estimator_asymp <- function(x, ...){
  xcomm <- table(x)
  inext_out <- iNEXT(x = list(as.numeric(xcomm)), datatype = "abundance", ...) # run iNEXT on the community
  filter(inext_out$AsyEst, Diversity == 'Species richness') %>%
    select(asymp_est = Estimator,
           asymp_est_stderr = s.e.,
           asymp_est_CImin = LCL,
           asymp_est_CImax = UCL)
}
