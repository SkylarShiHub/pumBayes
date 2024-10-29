#' @title Calculate WAIC of the dynamic probit unfolding model
#' @description This function is used to get the WAIC value of the dynamic probit unfolding model
#' @param vote_m A vote matrix in which rows represent members and columns represent issues.
#' The entries can only be 0 (indicating ”No”), 1 (indicating ”Yes”), or NA (indicating missing data).
#' @param years_v A vector representing the time period for each vote in the model
#' @param chain_run The posterior samples obtained from `sample_pum_dynamic` in `pumBayes`
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @examples
#' waic_dy = cal_waic_pum_dynamic(mqVotes, mqTime, post_samples_dy)
#' @export
cal_waic_pum_dynamic = function(vote_m, years_v, chain_run){
  block_m <- do.call(rbind, lapply(1:nrow(vote_m), function(i) {
    interested_inds <- which(!is.na(vote_m[i,]))
    cbind(i - 1, unique(sort(years_v[interested_inds])))
  }))
  leg_pos = as.matrix(chain_run$beta)
  alpha_draws0 = as.matrix(cbind(chain_run$alpha1, chain_run$alpha2))
  alpha_draws = t(apply(alpha_draws0, 1, function(row) {
    as.vector(t(matrix(row, ncol = 2)))
  }))
  rm(alpha_draws0)
  delta_draws0 = as.matrix(cbind(chain_run$delta1, chain_run$delta2))
  delta_draws = t(apply(delta_draws0, 1, function(row) {
    as.vector(t(matrix(row, ncol = 2)))
  }))
  rm(delta_draws0)
  dynamic_pum_waic_corrseed <- calc_waic_probit_bggum_three_utility_block_rcpp(
    leg_pos[, -grep("RHJackson_beta_9", colnames(leg_pos))], alpha_draws,
    delta_draws, mqVotes, mqTime, block_m)
  return(sum(dynamic_pum_waic_corrseed))
}
