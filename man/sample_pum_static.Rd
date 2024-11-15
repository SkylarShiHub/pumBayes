% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_pum_static.R
\name{sample_pum_static}
\alias{sample_pum_static}
\title{Generate posterior samples from the static probit unfolding model}
\usage{
sample_pum_static(
  vote_info,
  hyperparams,
  iter_config,
  pos_leg = 0,
  verbose = FALSE
)
}
\arguments{
\item{vote_info}{A vote matrix (or a rollcall object) in which rows represent members and columns represent issues.}

\item{hyperparams}{A list of hyperparameter values:
- `beta_mean`: Prior mean for beta.
- `beta_var`: Variance of beta.
- `alpha_mean`: A vector of two components representing the prior means of `alpha1` and `alpha2`.
- `alpha_scale`: Scale parameter for `alpha1` and `alpha2`.
- `delta_mean`: A vector of two components representing the prior means of `delta1` and `delta2`.
- `delta_scale`: Scale parameter for `delta1` and `delta2`.}

\item{iter_config}{A list of iteration configurations:
- `num_iter`: Total number of iterations. It is recommended to set this to at least 30,000 to ensure reliable results.
- `start_iter`: Iteration number to start retaining data after burn-in.
- `keep_iter`: Interval at which iterations are kept for posterior samples.
- `flip_rate`: Probability of directly flipping signs in the M-H step, rather than resampling from the priors.}

\item{pos_leg}{Name of the legislator whose position is kept positive.}

\item{verbose}{Logical. If `TRUE`, prints progress and additional information during the sampling process.}
}
\value{
A list primarily containing:
  - `beta`: A matrix of posterior samples for `beta`.
  - `alpha1`: A matrix of posterior samples for `alpha1`.
  - `alpha2`: A matrix of posterior samples for `alpha2`.
  - `delta1`: A matrix of posterior samples for `delta1`.
  - `delta2`: A matrix of posterior samples for `delta2`.
  - `vote_info`: The vote object after data preprocessing.
}
\description{
This function generates posterior samples of all parameters based on the static probit unfolding model.
}
\examples{
hyperparams <- list(beta_mean = 0, beta_var = 1, alpha_mean = c(0, 0),
                    alpha_scale = 5, delta_mean = c(-2, 10), delta_scale = sqrt(10))
iter_config <- list(num_iter = 10, start_iter = 0, keep_iter = 1, flip_rate = 0.1)
post_samples <- sample_pum_static(h116.c, hyperparams,
                                  iter_config, pos_leg = grep("SCALISE", rownames(h116.c$votes)),
                                  verbose = FALSE)
}
