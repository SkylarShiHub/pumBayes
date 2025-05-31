#' @title Generate posterior samples from the static probit unfolding model
#' @description This function generates posterior samples of all parameters based on the static probit unfolding model.
#' @param vote_info A logical vote matrix (or a rollcall object) in which rows represent members and columns represent issues.
#' @param hyperparams A list of hyperparameter values:
#'   - `beta_mean`: Prior mean for beta.
#'   - `beta_var`: Variance of beta.
#'   - `alpha_mean`: A vector of two components representing the prior means of `alpha1` and `alpha2`.
#'   - `alpha_scale`: Scale parameter for `alpha1` and `alpha2`.
#'   - `delta_mean`: A vector of two components representing the prior means of `delta1` and `delta2`.
#'   - `delta_scale`: Scale parameter for `delta1` and `delta2`.
#' @param control A list of MCMC configurations:
#'   - `num_iter`: Total number of iterations. It is recommended to set this to at least 30,000 to ensure reliable results.
#'   - `burn_in`: The number of initial iterations to discard as part of the burn-in period before retaining samples.
#'   - `keep_iter`: Interval at which iterations are kept for posterior samples.
#'   - `flip_rate`: Probability of directly flipping signs in the M-H step, rather than resampling from the priors.
#' @param pos_leg Name of the legislator whose position is kept positive.
#' @param verbose Logical. If `TRUE`, prints progress and additional information during the sampling process.
#' @param pre_run A list containing the output from a previous run of the function. If provided, the last iteration of the previous run will be used as the initial point of the new run. Defaults to `NULL`.
#' @param appended Logical. If `TRUE`, the new samples will be appended to the samples from the previous run. Defaults to `FALSE`.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return A list primarily containing:
#'   - `beta`: A matrix of posterior samples for `beta`.
#'   - `alpha1`: A matrix of posterior samples for `alpha1`.
#'   - `alpha2`: A matrix of posterior samples for `alpha2`.
#'   - `delta1`: A matrix of posterior samples for `delta1`.
#'   - `delta2`: A matrix of posterior samples for `delta2`.
#'   - `vote_info`: The input vote object.
#' @examples
#' \donttest{
#' # Long-running example
#' data(h116)
#' h116.c = preprocess_rollcall(h116)
#' hyperparams <- list(beta_mean = 0, beta_var = 1, alpha_mean = c(0, 0),
#'                     alpha_scale = 5, delta_mean = c(-2, 10), delta_scale = sqrt(10))
#' control <- list(num_iter = 2, burn_in = 0, keep_iter = 1, flip_rate = 0.1)
#' h116.c.pum <- sample_pum_static(h116.c, hyperparams,
#'                                   control, pos_leg = grep("SCALISE", rownames(h116.c$votes)),
#'                                   verbose = FALSE, pre_run = NULL, appended = FALSE)
#' }
#' @export
sample_pum_static <- function(vote_info, hyperparams, control,
                              pos_leg = 0, verbose = FALSE, pre_run = NULL, appended = FALSE) {

  # 1. Check and process input vote object
  if (is.matrix(vote_info)) {

    if (all(is.na(vote_info) | (vote_info %in% c(0, 1) & is.numeric(vote_info)))) {
      vote_m <- vote_info
    } else if (all(is.logical(vote_info))) {
      vote_m <- vote_info
      vote_m[vote_m == TRUE] <- 1
      vote_m[vote_m == FALSE] <- 0
    } else if (all(vote_info %in% c("T", "F", "NA"))) {
      vote_m <- vote_info
      vote_m[vote_m == "T"] <- 1
      vote_m[vote_m == "F"] <- 0
      vote_m[vote_m == "NA"] <- NA
    } else {
      invalid_values <- vote_info[!(is.na(vote_info) | vote_info %in% c(0, 1, TRUE, FALSE, "T", "F", "NA"))]
      if (length(invalid_values) > 0) {
        stop(paste("Invalid value found in your vote matrix:", paste(invalid_values, collapse = ", ")))
      }
    }


  } else if (is.list(vote_info) && "votes" %in% names(vote_info)) {
    vote_m <- vote_info$votes
    rownames_old <- rownames(vote_m)
    colnames_old <- colnames(vote_m)
    yea_values <- vote_info$codes$yea
    nay_values <- vote_info$codes$nay
    vote_m <- matrix(ifelse(vote_m %in% yea_values, 1,
                            ifelse(vote_m %in% nay_values, 0, NA)),
                     nrow = nrow(vote_m), ncol = ncol(vote_m))
    rownames(vote_m) <- rownames_old
    colnames(vote_m) <- colnames_old

  } else {
    stop("vote_info is not in a valid type. It has to be a matrix or a rollcall object.")
  }

# 2. Positive index

  if (length(pos_leg) != 1) {
    stop("`pos_leg` should contain exactly one legislator.")
  } else {
    pos_ind = pos_leg
  }

  # # output
  vote_out = vote_m

  total_iter = (control$num_iter - control$burn_in) %/% control$keep_iter

  # Initialize parameters from previous run or default
  if (!is.null(pre_run)) {
    init_info <- init_data_rcpp(
      vote_m,
      leg_pos_init = as.numeric(tail(pre_run$beta, 1)),
      alpha_pos_init = as.numeric(cbind(tail(pre_run$alpha1, 1), tail(pre_run$alpha2, 1))),
      delta_pos_init = as.numeric(cbind(tail(pre_run$delta1, 1), tail(pre_run$delta2, 1))),
      y_star_m_1_init = NULL,
      y_star_m_2_init = NULL,
      y_star_m_3_init = NULL,
      total_iter)
  } else {
    init_info <- init_data_rcpp(
      vote_m, leg_pos_init = NULL, alpha_pos_init = NULL, delta_pos_init = NULL,
      y_star_m_1_init = NULL, y_star_m_2_init = NULL, y_star_m_3_init = NULL,
      total_iter)
  }


  # c++
  draw_info <- sample_probit_static_rcpp(
    vote_m, init_info[[1]], init_info[[2]], init_info[[3]], init_info[[4]],
    init_info[[5]], init_info[[6]], init_info[[7]],
    init_info[[8]], init_info[[9]], hyperparams$beta_mean, sqrt(hyperparams$beta_var),
    hyperparams$alpha_mean, diag(2) * (hyperparams$alpha_scale^2),
    hyperparams$delta_mean, diag(2) * (hyperparams$delta_scale^2), 10000000,
    control$num_iter, control$burn_in, control$keep_iter, control$flip_rate,
    pos_ind - 1, verbose)

  all_param_draw = draw_info[[1]]
  if (!is.null(rownames(vote_m))){
    leg_names <- rownames(vote_m)
    # leg_names <- sapply(rownames(vote_m), function(name) {paste(name, "beta", sep = "_")})
  } else if (!is.null(rownames(vote_info$legis.data))){
    leg_names <- rownames(vote_info$legis.data)
  } else {
    leg_names <- rep("", nrow(vote_m))
  }

  if (is.null(colnames(vote_m))) {
    colnames(vote_m) <- sapply(1:ncol(vote_m), function(i) {
      paste("vote", i, sep = "_")
    })
  } # no operation performed
  alpha_vote_names_1 <- sapply(colnames(vote_m), function(name) {
    paste(name, "alpha", "1", sep = "_")
  })
  alpha_vote_names_2 <- sapply(colnames(vote_m), function(name) {
    paste(name, "alpha", "2", sep = "_")
  })
  delta_vote_names_1 <- sapply(colnames(vote_m), function(name) {
    paste(name, "delta", "1", sep = "_")
  })
  delta_vote_names_2 <- sapply(colnames(vote_m), function(name) {
    paste(name, "delta", "2", sep = "_")
  })
  colnames(all_param_draw) <-
    c(leg_names, alpha_vote_names_1, alpha_vote_names_2, delta_vote_names_1, delta_vote_names_2)

  beta_list <- as.data.frame(all_param_draw[, leg_names])
  alpha1_list <- as.data.frame(all_param_draw[, alpha_vote_names_1])
  alpha2_list <- as.data.frame(all_param_draw[, alpha_vote_names_2])
  delta1_list <- as.data.frame(all_param_draw[, delta_vote_names_1])
  delta2_list <- as.data.frame(all_param_draw[, delta_vote_names_2])

  # Append samples if required
  if (!is.null(pre_run) && appended) {
    beta_list <- rbind(pre_run$beta, beta_list)
    alpha1_list <- rbind(pre_run$alpha1, alpha1_list)
    alpha2_list <- rbind(pre_run$alpha2, alpha2_list)
    delta1_list <- rbind(pre_run$delta1, delta1_list)
    delta2_list <- rbind(pre_run$delta2, delta2_list)
  }

  return(list(
    beta = beta_list,
    alpha1 = alpha1_list,
    alpha2 = alpha2_list,
    delta1 = delta1_list,
    delta2 = delta2_list,
    votes = vote_info
  ))
}


init_y_star_m <- function(vote_m) {
  y_star_m_1 <- vote_m
  y_star_m_2 <- vote_m
  y_star_m_3 <- vote_m

  y_star_m_1[which(vote_m == 1)] = 0
  y_star_m_2[which(vote_m == 1)] = 1
  y_star_m_3[which(vote_m == 1)] = 0

  no_votes <- which(vote_m == 0)
  sample_upper <- rbinom(length(no_votes), size = 1, prob = 0.5)
  y_star_m_1[no_votes[which(sample_upper == 1)]] = 0
  y_star_m_2[no_votes[which(sample_upper == 1)]] = 0
  y_star_m_3[no_votes[which(sample_upper == 1)]] = 1

  y_star_m_1[no_votes[which(sample_upper == 0)]] = 1
  y_star_m_2[no_votes[which(sample_upper == 0)]] = 0
  y_star_m_3[no_votes[which(sample_upper == 0)]] = 0

  return(list(y_star_m_1, y_star_m_2, y_star_m_3))
}


init_data_rcpp <- function(vote_m, leg_pos_init, alpha_pos_init, delta_pos_init,
                           y_star_m_1_init, y_star_m_2_init, y_star_m_3_init, total_iter) {

  if (!is.null(leg_pos_init)) {
    leg_pos_m <-
      matrix(leg_pos_init, nrow = total_iter, ncol = nrow(vote_m), byrow = T)
  } else {
    leg_pos_m <- matrix(0, nrow = total_iter, ncol = nrow(vote_m)) # if it is not imputed
  }

  if (!is.null(alpha_pos_init)) {
    alpha_pos_m <-
      matrix(t(alpha_pos_init), nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T) # vector of 2
  } else {
    alpha_pos_m <-
      matrix(c(rep(1, ncol(vote_m)), rep(-1, ncol(vote_m))),
             nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  }

  if (!is.null(delta_pos_init)) {
    delta_pos_m <-
      matrix(t(delta_pos_init), nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T) # vector of 2
  } else {
    delta_pos_m <-
      matrix(0, nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  }

  if (!is.null(y_star_m_1_init)) {
    y_star_m_1 <- y_star_m_1_init
    y_star_m_2 <- y_star_m_2_init
    y_star_m_3 <- y_star_m_3_init
  } else {
    y_star_info <- init_y_star_m(vote_m)
    y_star_m_1 <- y_star_info[[1]]
    y_star_m_2 <- y_star_info[[2]]
    y_star_m_3 <- y_star_info[[3]]
  }

  all_params_draw <- cbind(leg_pos_m, alpha_pos_m, delta_pos_m) # initial values of beta, alpha, delta
  beta_start_ind = 0; # used in C++
  alpha_start_ind = nrow(vote_m);
  alpha_2_start_ind = alpha_start_ind + ncol(vote_m);
  delta_start_ind = alpha_2_start_ind + ncol(vote_m);
  delta_2_start_ind = delta_start_ind + ncol(vote_m);

  return(list(all_params_draw, y_star_m_1, y_star_m_2,
              y_star_m_3, beta_start_ind,
              alpha_start_ind, alpha_2_start_ind,
              delta_start_ind, delta_2_start_ind))
}
