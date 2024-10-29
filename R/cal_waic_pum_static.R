#' @title Calculate WAIC of the static probit unfolding model
#' @description This function is used to get the WAIC value of the static probit unfolding model
#' @param chain_run The object output from 'sample_pum_static' in 'pumBayes'.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return The WAIC value of the static probit unfoldinng model
#' @examples
#' waic_pum_static = cal_waic_pum_static(chain_run = post_samples_pum)
#' @export

cal_waic_pum_static <- function(chain_run){

  # Check and process input vote object
  vote_info = chain_run$vote_info

  if (is.matrix(vote_info)) {
    vote_m <- vote_info

  } else if (is.list(vote_info) && "votes" %in% names(vote_info)) {
    vote_m <- vote_info$votes

  } else {

    stop("Vote matrix not found.")
  }

  leg_draws <- as.matrix(chain_run$beta)
  alpha_draws0 = as.matrix(cbind(chain_run$alpha1, chain_run$alpha2))
  alpha_draws <- t(apply(alpha_draws0, 1, function(row) {
    as.vector(t(matrix(row, ncol = 2)))
  }))
  rm(alpha_draws0)
  delta_draws0 = as.matrix(cbind(chain_run$delta1, chain_run$delta2))
  delta_draws <- t(apply(delta_draws0, 1, function(row) {
    as.vector(t(matrix(row, ncol = 2)))
  }))
  rm(delta_draws0)
  block_m <- cbind(1:nrow(vote_m) - 1, 1)

  pum_waic <- calc_waic_probit_bggum_three_utility_block_rcpp(
    leg_draws, alpha_draws, delta_draws, vote_m,
    rep(1,ncol(vote_m)), block_m)
  rm(leg_draws)
  rm(alpha_draws)
  rm(delta_draws)
  return(sum(pum_waic))
}
