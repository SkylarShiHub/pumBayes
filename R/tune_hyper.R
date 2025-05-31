#' @title Generate Probability Samples for Voting "Yes"
#' @description This function generates probability samples for Voting "Yes".
#' It uses predefined hyperparameters and simulates data based on the specified number of members (`n_leg`) and issues (`n_issue`).
#' @param hyperparams A list of hyperparameter values:
#'   - `beta_mean`: The prior mean of the `beta` parameter, representing legislator positions.
#'   - `beta_var`: The prior variance of `beta`.
#'   - `alpha_mean`: A vector of length two, specifying the prior means of the item discrimination parameters, `alpha1` and `alpha2`.
#'   - `alpha_scale`: The scale parameter for `alpha1` and `alpha2`.
#'   - `delta_mean`: A vector of length two, indicating the prior means of the item difficulty parameters, `delta1` and `delta2`.
#'   - `delta_scale`: The scale parameter for `delta1` and `delta2`.
#' @param n_leg Integer, representing the number of legislators (members) to be simulated.
#' @param n_issue Integer, indicating the number of issues to be simulated.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return A numeric vector containing the simulated probabilities of voting "Yes" for legislators across issues.
#' @examples
#' hyperparams = list(beta_mean = 0, beta_var = 1, alpha_mean = c(0, 0),
#'                    alpha_scale = 5, delta_mean = c(-2, 10),
#'                    delta_scale = sqrt(10))
#' theta = tune_hyper(hyperparams, n_leg = 10, n_issue = 10)
#' @export
tune_hyper <- function(hyperparams = hyperparams, n_leg, n_issue) {
  samples <- matrix(0, nrow = n_issue, ncol = 4)
  for (i in 1:n_issue) {
    if (runif(1) < 0.5) {
      alpha_j1 <- truncated_t_sample(10000000, hyperparams$alpha_mean[1], hyperparams$alpha_scale, TRUE)
      alpha_j2 <- truncated_t_sample(10000000, hyperparams$alpha_mean[2], hyperparams$alpha_scale, FALSE)
      delta_j1 <- sample_t(10000000, hyperparams$delta_mean[1], hyperparams$delta_scale)
      delta_j2 <- sample_t(10000000, hyperparams$delta_mean[2], hyperparams$delta_scale)
    } else {
      alpha_j1 <- truncated_t_sample(10000000, hyperparams$alpha_mean[1], hyperparams$alpha_scale, FALSE)
      alpha_j2 <- truncated_t_sample(10000000, hyperparams$alpha_mean[2], hyperparams$alpha_scale, TRUE)
      delta_j1 <- sample_t(10000000, -hyperparams$delta_mean[1], hyperparams$delta_scale)
      delta_j2 <- sample_t(10000000, -hyperparams$delta_mean[2], hyperparams$delta_scale)
    }

    samples[i, ] <- c(alpha_j1, alpha_j2, delta_j1, delta_j2)
  }
  beta = rnorm(n_leg,hyperparams$beta_mean, sqrt(hyperparams$beta_var))
  samples = list(beta = beta, alpha1 = samples[,1], alpha2 = samples[,2],
                 delta1 = samples[,3], delta2 = samples[,4])
  mat <- matrix(1, nrow = n_leg, ncol = n_issue)
  probability = get_prob_mat(mat, samples$beta, samples$alpha1, samples$alpha2,
                         samples$delta1, samples$delta2)
  mu_string <- paste("(", paste(hyperparams$delta_mean, collapse = ", "), ")", sep = "")
  return(as.vector(probability))
}


get_prob_mat <- function(vote, beta, alpha1, alpha2, delta1, delta2) {

  prob <- matrix(NA, nrow = nrow(vote), ncol = ncol(vote))

  for (j in 1:ncol(vote)) {

    term1 <- -alpha1[j] * (beta - delta1[j]) / sqrt(2)
    term2 <- -alpha2[j] * (beta - delta2[j]) / sqrt(2)
    bvnd_vals <- bvndvec(term1, term2, rep(0.5,length(beta)))
    prob[, j] <- ifelse(is.na(vote[, j]), NA,
                        ifelse(vote[, j] == 1, bvnd_vals, 1 - bvnd_vals))
  }
  rownames(prob) <- rownames(vote)
  colnames(prob) <- colnames(vote)
  return(prob)
}
