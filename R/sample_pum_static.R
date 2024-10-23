#' @title Generate posterior samples from the static probit unfolding model
#' @description This function is used to get samples of all the parameters from their posterior distributions
#' based on the static probit unfolding model
#' @param vote_info A vote matrix (or a rollcall object including votes) in which rows represent members and columns represent issues.
#' The entries can only be 0 (indicating ”No”), 1 (indicating ”Yes”), or NA (indicating missing data).
#' @param hyperparams A list of hyperparameter values. beta_mean is a value for the prior mean of beta.
#' beta_var is a value for the variance of beta. alpha_mean is a vector of 2 components representing the prior
#' means of alpha1 and alpha2. alpha_scale is a value for the scale parameters of alpha1 and alpha2.
#' delta_mean is a vector of 2 components representing the prior means of delta1 and delta2.
#' delta_scale is a value for the scale parameters of delta1 and delta2. nu is a fixed value making phi follow
#' Gamma(nu/2,nu/2).
#' @param iter_config A list corresponding to iteration configurations, num_iter is the total number of iterations.
#' start_iter is the iteration number to start retaining data after the burn-in period. keep_iter is the interval
#' at which iterations are kept for posterior samples. flip_rate is the probability of fliping signs directly in the M-H step.
#' @param leg_pos_init A vector of initial values for all the positions of legislators. The
#' length of this vector should be equal to the number of iterations times the number of members
#' @param alpha_pos_init A vector of initial values for all the alpha parameters. This vector
#' should have a length of 2 times the number of iterations times the number of issues,
#' as each issue has two alpha parameters.
#' @param delta_pos_init A vector of initial values for all the delta parameters. The length
#' should be 2 times the number of iterations times the number of issues
#' @param y_star_m_1_init Matrices of initial values for the three auxiliary variables, with
#' dimensions equal to the number of members by the number of issues. Same for the other two parameters.
#' @param pos_ind The index of the reference member whose position is kept positive.
#' @param start_val A vector of all initial values for each parameter.
#' @param sample_beta A Boolean value indicating whether to sample beta in the model.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return A list mainly containing:
#' - `beta`: A matrix of parameter draws for beta.
#' - `alpha1`: A matrix of parameter draws of alpha1.
#' - `alpha2`: A matrix of parameter draws of alpha2.
#' - `delta1`: A matrix of parameter draws of delta1.
#' - `delta2`: A matrix of parameter draws of delta2.
#' @examples
#' hyperparams = list(beta_mean = 0, beta_var = 1, alpha_mean = c(0, 0),
#'                    alpha_scale = 5, delta_mean = c(-2, 10), delta_scale = sqrt(10))
#' iter_config = list(num_iter = 10, start_iter = 0, keep_iter = 1, flip_rate = 0.1)
#' data_process = list(leg_rm = NULL, issue_rm = NULL, combine_ind = NULL, combine_party = NULL)
#' post_samples = sample_pum_static(house_votes_m, data_process, hyperparams, iter_config, pos_ind = 82)
#' @export
sample_pum_static <- function(vote_info, data_process=list(leg_rm = NULL, issue_rm = NULL,
                                                              combine_ind = NULL, combine_party = NULL),
                                 hyperparams, iter_config, pos_ind = 0) {

  if (is.matrix(vote_info)) {
    vote_m <- vote_info
  } else if (is.list(vote_info) && "votes" %in% names(vote_info)) {
    vote_m <- vote_info$votes
    rownames_old <- rownames(vote_m)
    colnames_old <- colnames(vote_m)
    yea_values <- c(1, 2, 3)
    nay_values <- c(4, 5, 6)
    vote_m <- matrix(ifelse(vote_m %in% yea_values, 1,
                            ifelse(vote_m %in% nay_values, 0, NA)),
                     nrow = nrow(vote_m), ncol = ncol(vote_m))
    rownames(vote_m) <- rownames_old
    colnames(vote_m) <- colnames_old

  } else {
    stop("vote_info is not in a valid type. It has to be a matrix or a rollcall object.")
  }

  # Delete all unanimously stance columns
  unanimous_votes <- which(apply(vote_m, 2, function(col) {length(table(col)) == 1}))
  if (length(unanimous_votes) > 0) {
    vote_m <- vote_m[, -unanimous_votes]
  }

  # Delete rows where missing values is greater than or equal to 60%
  absent_members <- which(rowMeans(is.na(vote_m)) >= 0.6)
  if (length(absent_members) > 0) {
    vote_m <- vote_m[-absent_members,]
  }

  leg_rm = data_process$leg_rm
  issue_rm = data_process$issue_rm
  combine_ind = data_process$combine_ind
  combine_party = data_process$combine_party

  if (!is.null(leg_rm)) {
    vote_m <- vote_m[-leg_rm, , drop = FALSE]
  }

  if (!is.null(issue_rm)) {
    vote_m <- vote_m[, -issue_rm, drop = FALSE]
  }

  if (!is.null(combine_ind)){
    for (i in seq_along(combine_ind)) {

      indices <- combine_ind[[i]]
      row1 <- vote_m[indices[1], ]
      row2 <- vote_m[indices[2], ]

      vote_m[indices[1], ] <- ifelse(is.na(row1), row2, row1)

      if (!is.null(combine_party)){
        old_name <- rownames(vote_m)[indices[1]]
        open_bracket_pos <- regexpr("\\(", old_name)[1]
        new_name <- paste0(
          substr(old_name, 1, open_bracket_pos),
          combine_party[i],
          substr(old_name, open_bracket_pos + 2, nchar(old_name)))
        rownames(vote_m)[indices[1]] <- new_name
      } else {
        namelist = c(rownames(vote_m)[indices[1]], rownames(vote_m)[indices[2]])
        rownames(vote_m)[indices[1]] <- namelist[grep("\\(I ", namelist)]
      }

      vote_m <- vote_m[-indices[2], , drop = FALSE]
    }
  }

  total_iter = (iter_config$num_iter - iter_config$start_iter) %/% iter_config$keep_iter
  init_info <- init_data_rcpp(
    vote_m, leg_pos_init = NULL, alpha_pos_init = NULL,
    delta_pos_init = NULL, y_star_m_1_init = NULL, y_star_m_2_init = NULL,
    y_star_m_3_init = NULL, total_iter)

  # c++
  draw_info <- sample_probit_static_rcpp(
    vote_m, init_info[[1]], init_info[[2]], init_info[[3]], init_info[[4]],
    init_info[[5]], init_info[[6]], init_info[[7]],
    init_info[[8]], init_info[[9]], hyperparams$beta_mean, sqrt(hyperparams$beta_var),
    hyperparams$alpha_mean, diag(2) * (hyperparams$alpha_scale^2),
    hyperparams$delta_mean, diag(2) * (hyperparams$delta_scale^2), 10000000,
    iter_config$num_iter, iter_config$start_iter, iter_config$keep_iter, iter_config$flip_rate,
    pos_ind - 1)

  all_param_draw = draw_info[[1]]
  leg_names <- sapply(rownames(vote_m), function(name) {paste(name, "beta", sep = "_")})
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

  return(list(
    beta = beta_list,
    alpha1 = alpha1_list,
    alpha2 = alpha2_list,
    delta1 = delta1_list,
    delta2 = delta2_list,
    other_info = draw_info[-1],
    vote_mat = vote_m,
    leg_names = sub(" .*", "", rownames(vote_m)),
    party = sub(".*\\((.).*", "\\1", rownames(vote_m))
  ))
}

#' @title initialize three auxiliary parameters y
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

#' @title initialize member and position parameters and get starting points
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
