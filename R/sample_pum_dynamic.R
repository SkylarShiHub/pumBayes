#' @title Generate posterior samples from the dynamic probit unfolding model
#' @description This function generates posterior samples for all parameters based on the dynamic probit unfolding model.
#' @param vote_info A logical vote matrix where rows represent members and columns represent issues.
#' The entries should be FALSE ("No"), TRUE ("Yes"), or NA (missing data).
#' @param years_v A vector representing the time period for each vote in the model.
#' @param hyperparams A list of hyperparameter values including:
#'   - `beta_mean`: Prior mean of beta.
#'   - `beta_var`: Prior variance of beta.
#'   - `alpha_mean`: A vector of 2 values for the prior means of alpha1 and alpha2.
#'   - `alpha_scale`: Scale parameter for alpha1 and alpha2.
#'   - `delta_mean`: A vector of 2 values for the prior means of delta1 and delta2.
#'   - `delta_scale`: Scale parameter for delta1 and delta2.
#'   - `rho_mean`: Prior mean of the autocorrelation parameter `rho`.
#'   - `rho_sigma`: Standard deviation of the prior for `rho`.
#' @param control A list specifying the MCMC configurations, including:
#'   - `num_iter`: Total number of iterations.
#'   - `burn_in`: The number of initial iterations to discard as part of the burn-in period before retaining samples.
#'   - `keep_iter`: Interval at which samples are retained.
#'   - `flip_rate`: Probability of directly flipping signs in the M-H step, rather than resampling from the priors.
#'   - `sd_prop_rho`: Proposal standard deviation for `rho` in the Metropolis-Hastings step.
#' @param sign_refs A list containing sign constraints, including:
#'   - `pos_inds`: Indices of members constrained to have positive values.
#'   - `neg_inds`: Indices of members constrained to have negative values.
#'   - `pos_year_inds`: List of years corresponding to each `pos_ind`.
#'   - `neg_year_inds`: List of years corresponding to each `neg_ind`.
#' @param verbose Logical. If `TRUE`, prints progress and additional information during the sampling process.
#' @param pre_run A list containing the output from a previous run of the function. If provided, the last iteration of the previous run will be used as the initial point of the new run. Defaults to `NULL`.
#' @param appended Logical. If `TRUE`, the new samples will be appended to the samples from the previous run. Defaults to `FALSE`.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return A list containing:
#'   - `beta`: A data frame of posterior samples for beta.
#'   - `alpha1`: A data frame of posterior samples for alpha1.
#'   - `alpha2`: A data frame of posterior samples for alpha2.
#'   - `delta1`: A data frame of posterior samples for delta1.
#'   - `delta2`: A data frame of posterior samples for delta2.
#'   - `rho`: A data frame of posterior samples for rho.
#' @examples
#' \donttest{
#' # Long-running example
#' data(scotus.1937.2021)
#' hyperparams = list(alpha_mean = c(0, 0), alpha_scale = 5,
#'                    delta_mean = c(-2, 10), delta_scale = sqrt(10),
#'                    rho_mean = 0.9, rho_sigma = 0.04)
#' control = list(num_iter = 2, burn_in = 0, keep_iter = 1, flip_rate = 0.1, sd_prop_rho = 0.1)
#' sign_refs = list(pos_inds = c(39, 5), neg_inds = c(12, 29),
#'                  pos_year_inds = list(1:31, 1), neg_year_inds = list(1:29, 1:24))
#' scotus.pum = sample_pum_dynamic(mqVotes, mqTime, hyperparams, control, sign_refs,
#'                                 verbose = FALSE, pre_run = NULL, appended = FALSE)
#' }
#' @export
sample_pum_dynamic <- function(
    vote_info, years_v, hyperparams, control, sign_refs, verbose = FALSE, pre_run = NULL, appended = FALSE) {

  # check inputs
  if(!(ncol(vote_info) == length(years_v))){
    stop("The number of columns in `vote_info` does not match the length of `years_v`.")
  }

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

  } else {
    stop("`vote_info` should be a matrix.")
  }


  total_iter = (control$num_iter - control$burn_in) %/% control$keep_iter

  # Initialize parameters from previous run or default
  if (!is.null(pre_run)) {
    init_info <- init_data_gp_rcpp(
      vote_m, years_v,
      leg_pos_init = as.numeric(tail(pre_run$beta, 1)),
      alpha_pos_init = as.numeric(cbind(tail(pre_run$alpha1, 1), tail(pre_run$alpha2, 1))),
      delta_pos_init = as.numeric(cbind(tail(pre_run$delta1, 1), tail(pre_run$delta2, 1))),
      rho_init = unlist(pre_run$rho)[length(unlist(pre_run$rho))],
      y_star_m_1_init = NULL,
      y_star_m_2_init = NULL,
      y_star_m_3_init = NULL,
      total_iter,
      sign_refs$pos_inds, sign_refs$neg_inds,
      sign_refs$pos_year_inds, sign_refs$neg_year_inds)
  } else {
    init_info <- init_data_gp_rcpp(
      vote_m, years_v, leg_pos_init = NULL, alpha_pos_init = NULL, delta_pos_init = NULL,
      rho_init = NULL, y_star_m_1_init = NULL, y_star_m_2_init = NULL,
      y_star_m_3_init = NULL, total_iter,
      sign_refs$pos_inds, sign_refs$neg_inds,
      sign_refs$pos_year_inds, sign_refs$neg_year_inds)
  }


  #Init info
  draw_info <- sample_probit_dynamic_rcpp(
    vote_m, init_info[[1]], init_info[[2]], init_info[[3]], init_info[[4]],
    init_info[[11]], init_info[[12]], years_v,
    init_info[[13]], unlist(init_info[[14]]), init_info[[6]], init_info[[7]],
    init_info[[8]], init_info[[9]], init_info[[10]],
    hyperparams$alpha_mean, diag(2) * (hyperparams$alpha_scale^2),
    hyperparams$delta_mean, diag(2) * (hyperparams$delta_scale^2),
    hyperparams$rho_mean, hyperparams$rho_sigma, control$sd_prop_rho, 10000000,
    control$num_iter, control$burn_in, control$keep_iter,
    control$flip_rate, init_info[[15]], init_info[[16]],
    init_info[[17]], init_info[[18]], verbose)

  all_param_draw = draw_info[[1]]
  leg_names <- unlist(sapply(1:nrow(vote_m), function(i) {
    sapply(init_info[[14]][[i]], function(year) {
      paste(rownames(vote_m)[i], "beta", year, sep = "_")
    })
  }))

  if (is.null(colnames(vote_m))) {
    colnames(vote_m) <- sapply(1:ncol(vote_m), function(i) {
      paste("vote", i, sep = "_")
    })
  }
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
    c(leg_names, alpha_vote_names_1, alpha_vote_names_2,
      delta_vote_names_1, delta_vote_names_2, "rho")

  beta_list <- as.data.frame(all_param_draw[, leg_names])
  alpha1_list <- as.data.frame(all_param_draw[, alpha_vote_names_1])
  alpha2_list <- as.data.frame(all_param_draw[, alpha_vote_names_2])
  delta1_list <- as.data.frame(all_param_draw[, delta_vote_names_1])
  delta2_list <- as.data.frame(all_param_draw[, delta_vote_names_2])
  rho_list <- as.data.frame(all_param_draw[, "rho"])
  colnames(rho_list) <- "rho"

  # Append samples if required
  if (!is.null(pre_run) && appended) {
    beta_list <- rbind(pre_run$beta, beta_list)
    alpha1_list <- rbind(pre_run$alpha1, alpha1_list)
    alpha2_list <- rbind(pre_run$alpha2, alpha2_list)
    delta1_list <- rbind(pre_run$delta1, delta1_list)
    delta2_list <- rbind(pre_run$delta2, delta2_list)
    rho_list <- rbind(pre_run$rho, rho_list)
  }

  return(list(
    beta = beta_list,
    alpha1 = alpha1_list,
    alpha2 = alpha2_list,
    delta1 = delta1_list,
    delta2 = delta2_list,
    rho = rho_list))
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


init_data_gp_rcpp <- function(
    vote_m, years_v, leg_pos_init, alpha_pos_init, delta_pos_init, rho_init,
    y_star_m_1_init, y_star_m_2_init, y_star_m_3_init, total_iter,
    pos_ind_list, neg_ind_list, pos_ind_years_list, neg_ind_years_list) {

  judge_start_inds <- rep(0, nrow(vote_m))
  judge_end_inds <- rep(0, nrow(vote_m))

  case_judge_years_ind_m <- vote_m
  years_served <- rep(0, nrow(vote_m))
  judge_year_v <- c()
  judge_start_inds <- c()
  judge_end_inds <- c()
  start_ind = 0
  end_ind = 0
  for (i in 1:nrow(vote_m)) {
    interested_inds <- which(!is.na(vote_m[i,]))
    case_judge_years_ind_m[i,interested_inds] <-
      years_v[interested_inds] - min(years_v[interested_inds])
    judge_year_v <-
      c(judge_year_v, list(min(years_v[interested_inds]):max(years_v[interested_inds])))
    end_ind = start_ind +
      max(years_v[interested_inds]) - min(years_v[interested_inds])
    judge_start_inds <- c(judge_start_inds, start_ind)
    judge_end_inds <- c(judge_end_inds, end_ind)
    start_ind = end_ind + 1
  }

  if (!is.null(leg_pos_init)) {
    leg_pos_m <-
      matrix(leg_pos_init, nrow = total_iter, ncol = length(leg_pos_init), byrow = T)
  } else {
    leg_pos_m <- matrix(0, nrow = total_iter, ncol = max(judge_end_inds) + 1)
  }

  if (!is.null(alpha_pos_init)) {
    alpha_pos_m <-
      matrix(t(alpha_pos_init), nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  } else {
    alpha_pos_m <-
      matrix(c(rep(1, ncol(vote_m)),rep(-1, ncol(vote_m))),
             nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  }

  if (!is.null(delta_pos_init)) {
    delta_pos_m <-
      matrix(t(delta_pos_init), nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  } else {
    delta_pos_m <-
      matrix(0, nrow = total_iter, ncol = 2 * ncol(vote_m), byrow = T)
  }

  if (!is.null(rho_init)) {
    rho_v <- rep(rho_init, nrow = total_iter)
  } else {
    rho_v <- rep(0.9, nrow = total_iter)
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

  all_params_draw <- cbind(leg_pos_m, alpha_pos_m, delta_pos_m, rho_v)
  beta_start_ind = 0;
  alpha_start_ind = max(judge_end_inds) + 1;
  alpha_2_start_ind = alpha_start_ind + ncol(vote_m);
  delta_start_ind = alpha_2_start_ind + ncol(vote_m);
  delta_2_start_ind = delta_start_ind + ncol(vote_m);
  rho_start_ind = delta_2_start_ind + ncol(vote_m);

  pos_ind_judge_list <- vector(mode = "integer")
  pos_ind_judge_year_list <- vector(mode = "integer")
  if (length(pos_ind_list) > 0) {
    for (i in 1:length(pos_ind_list)) {
      tmp_judge_list <-
        judge_start_inds[pos_ind_list[i]]:
        judge_end_inds[pos_ind_list[i]]
      tmp_judge_year_list <- judge_year_v[[pos_ind_list[i]]]
      if (length(pos_ind_years_list) > 0) {
        tmp_judge_list <- tmp_judge_list[pos_ind_years_list[[i]]]
        tmp_judge_year_list <-
          tmp_judge_year_list[pos_ind_years_list[[i]]]
      }
      pos_ind_judge_list <- c(pos_ind_judge_list, tmp_judge_list)
      pos_ind_judge_year_list <-
        c(pos_ind_judge_year_list, tmp_judge_year_list)
    }
  }
  neg_ind_judge_list <- vector(mode = "integer")
  neg_ind_judge_year_list <- vector(mode = "integer")
  if (length(neg_ind_list) > 0) {
    for (i in 1:length(neg_ind_list)) {
      tmp_judge_list <-
        judge_start_inds[neg_ind_list[i]]:
        judge_end_inds[neg_ind_list[i]]
      tmp_judge_year_list <- judge_year_v[[neg_ind_list[i]]]
      if (length(neg_ind_years_list) > 0) {
        tmp_judge_list <- tmp_judge_list[neg_ind_years_list[[i]]]
        tmp_judge_year_list <-
          tmp_judge_year_list[neg_ind_years_list[[i]]]
      }
      neg_ind_judge_list <- c(neg_ind_judge_list, tmp_judge_list)
      neg_ind_judge_year_list <-
        c(neg_ind_judge_year_list, tmp_judge_year_list)
    }
  }

  return(list(all_params_draw, y_star_m_1, y_star_m_2,
              y_star_m_3, beta_start_ind,
              alpha_start_ind, alpha_2_start_ind,
              delta_start_ind, delta_2_start_ind, rho_start_ind,
              judge_start_inds, judge_end_inds,
              case_judge_years_ind_m, judge_year_v,
              pos_ind_judge_list, neg_ind_judge_list,
              pos_ind_judge_year_list, neg_ind_judge_year_list))
}
