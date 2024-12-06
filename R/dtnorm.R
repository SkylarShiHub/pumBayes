#' @title Density Function for Truncated Normal Distribution
#' @description
#' This function calculates the density of a truncated normal distribution at specified points.
#' @param x A numeric vector of quantiles at which to evaluate the density.
#' @param mean A numeric value specifying the mean of the normal distribution (default is 0).
#' @param sd A numeric value specifying the standard deviation of the normal distribution (default is 1, must be positive).
#' @param lower A numeric value specifying the lower bound of truncation (default is -Inf).
#' @param upper A numeric value specifying the upper bound of truncation (default is Inf).
#' @return A numeric vector of density values corresponding to the input `x`.
#' The values are normalized to ensure the total probability within the truncation bounds equals 1.
#' Values outside the truncation bounds are set to 0.
#' @examples
#' dtnorm(c(-1, 0, 1), mean = 0, sd = 1, lower = -1, upper = 1)
#' @export
dtnorm <- function(x, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  # Validate standard deviation
  if (sd <= 0) stop("Standard deviation must be positive.")
  
  # Calculate the untruncated density
  density <- dnorm(x, mean = mean, sd = sd)
  
  # Calculate the normalization constant
  Z <- pnorm(upper, mean = mean, sd = sd) - pnorm(lower, mean = mean, sd = sd)
  
  # Calculate the truncated density
  density_truncated <- ifelse(x >= lower & x <= upper, density / Z, 0)
  
  return(density_truncated)
}
