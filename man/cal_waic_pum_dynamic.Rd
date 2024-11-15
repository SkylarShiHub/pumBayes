% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal_waic_pum_dynamic.R
\name{cal_waic_pum_dynamic}
\alias{cal_waic_pum_dynamic}
\title{Calculate WAIC of the dynamic probit unfolding model}
\usage{
cal_waic_pum_dynamic(vote_m, years_v, chain_run)
}
\arguments{
\item{vote_m}{A vote matrix in which rows represent members and columns represent issues.
The entries can only be 0 (indicating ”No”), 1 (indicating ”Yes”), or NA (indicating missing data).}

\item{years_v}{A vector representing the time period for each vote in the model}

\item{chain_run}{The posterior samples obtained from `sample_pum_dynamic` in `pumBayes`}
}
\description{
This function is used to get the WAIC value of the dynamic probit unfolding model
}
\examples{
waic_pum_dy = cal_waic_pum_dynamic(mqVotes, mqTime, post_samples_pum_dy)
}
