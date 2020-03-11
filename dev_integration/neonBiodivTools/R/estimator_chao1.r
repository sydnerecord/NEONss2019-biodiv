#' @title Chao1 richness estimator
#'
#' @import iNEXT dplyr tibble
#'
#' @param x A vector of species IDs.
#'
#' @return A data frame of Chao 1 estiamtion and variance and 95% CI.
#'
#' @export
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
