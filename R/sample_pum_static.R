#' @title Generate posterior samples from the static probit unfolding model
#' @description This function generates posterior samples of all parameters based on the static probit unfolding model.
#' @param vote_info A vote matrix (or a rollcall object) in which rows represent members and columns represent issues.
#' @param data_process A list of parameters for data preprocessing:
#'   - `leg_rm`: A list of legislators to be removed.
#'   - `combine_leg_name` and `combine_leg_party`: Vectors for merging legislators' votes by name and party.
#'   - `lop_leg`: The proportion of missing votes that warrants removing a legislator.
#'   - `lop_issue`: The threshold of minority votes below which an issue will be excluded.
#' @param hyperparams A list of hyperparameter values:
#'   - `beta_mean`: Prior mean for beta.
#'   - `beta_var`: Variance of beta.
#'   - `alpha_mean`: A vector of two components representing the prior means of `alpha1` and `alpha2`.
#'   - `alpha_scale`: Scale parameter for `alpha1` and `alpha2`.
#'   - `delta_mean`: A vector of two components representing the prior means of `delta1` and `delta2`.
#'   - `delta_scale`: Scale parameter for `delta1` and `delta2`.
#' @param iter_config A list of iteration configurations:
#'   - `num_iter`: Total number of iterations. It is recommended to set this to at least 30,000 to ensure reliable results.
#'   - `start_iter`: Iteration number to start retaining data after burn-in.
#'   - `keep_iter`: Interval at which iterations are kept for posterior samples.
#'   - `flip_rate`: Probability of directly flipping signs in the M-H step, rather than resampling from the priors.
#' @param pos_leg Name of the legislator whose position is kept positive.
#' @param verbose Logical. If `TRUE`, prints progress and additional information during the sampling process.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return A list primarily containing:
#'   - `beta`: A matrix of posterior samples for `beta`.
#'   - `alpha1`: A matrix of posterior samples for `alpha1`.
#'   - `alpha2`: A matrix of posterior samples for `alpha2`.
#'   - `delta1`: A matrix of posterior samples for `delta1`.
#'   - `delta2`: A matrix of posterior samples for `delta2`.
#'   - `vote_info`: The vote object after data preprocessing.
#' @examples
#' hyperparams <- list(beta_mean = 0, beta_var = 1, alpha_mean = c(0, 0),
#'                     alpha_scale = 5, delta_mean = c(-2, 10), delta_scale = sqrt(10))
#' iter_config <- list(num_iter = 10, start_iter = 0, keep_iter = 1, flip_rate = 0.1)
#' data_process <- list(leg_rm = NULL, combine_leg_name = NULL, combine_leg_party = NULL,
#'                      lop_leg = 0.6, lop_issue = 0)
#' post_samples <- sample_pum_static(house_votes_m, data_process, hyperparams,
#'                                   iter_config, pos_leg = "SCALISE", verbose = FALSE)
#' @export
sample_pum_static <- function(vote_info,
                              data_process=list(leg_rm = NULL, combine_leg_name = NULL, combine_leg_party = NULL,
                                                lop_leg = 0.6, lop_issue = 0),
                              hyperparams, iter_config, pos_leg = "SCALISE", verbose = FALSE) {
  leg_rm = data_process$leg_rm
  combine_leg_name = data_process$combine_leg_name
  combine_leg_party = data_process$combine_leg_party
  lop_leg = data_process$lop_leg
  lop_issue = data_process$lop_issue

  ## checks before data preprocess
  # 1. Check and process input vote object
  if (is.matrix(vote_info)) {

    if (all(is.na(vote_info) | vote_info %in% c(0, 1))) {
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

  # 2. Check for name conflicts
  common_names <- intersect(leg_rm, combine_leg_name)
  if (length(common_names) > 0) {
    stop(paste("The following names appear in both leg_rm and combine_leg_name:", paste(common_names, collapse = ", ")))
  }

  # 3. Check for input combined legislators
  if (!(all(is.na(combine_leg_name)) && all(is.na(combine_leg_party))) &&
      length(combine_leg_name) != length(combine_leg_party)) {
    stop("Legislator names and parties to be combined must have the same length.")
  }
  unmatched_names <- sapply(combine_leg_name, function(name) {
    if (length(grep(name, rownames(vote_m))) == 0) {
      return(name)
    } else {
      return(NULL)
    }
  })
  unmatched_names <- unlist(unmatched_names)
  if (length(unmatched_names) > 0) {
    stop(paste("The following names are not found in the given name list:", paste(unmatched_names, collapse = ", ")))
  }
  invalid_parties <- setdiff(combine_leg_party, c("D", "R", "I"))
  if (length(invalid_parties) > 0) {
    stop(paste("The following parties are invalid:", paste(invalid_parties, collapse = ", ")))
  }

  ## data preprocess
  # 1. Remove legislators
  if (!is.null(leg_rm)) {
    vote_m <- vote_m[-grep(leg_rm, rownames(vote_m)), , drop = FALSE]
  }
  # 2. combine legislators
  if (!is.null(combine_leg_name)){
    for (i in seq_along(combine_leg_name)) {

      indices <- grep(combine_leg_name[i], rownames(vote_m))
      row1 <- vote_m[indices[1], ]
      row2 <- vote_m[indices[2], ]

      vote_m[indices[1], ] <- ifelse(is.na(row1), row2, row1)

      if (!is.null(combine_leg_party)){
        old_name <- rownames(vote_m)[indices[1]]
        open_bracket_pos <- regexpr("\\(", old_name)[1]
        new_name <- paste0(
          substr(old_name, 1, open_bracket_pos),
          combine_leg_party[i],
          substr(old_name, open_bracket_pos + 2, nchar(old_name)))
        rownames(vote_m)[indices[1]] <- new_name
      } else {
        namelist = c(rownames(vote_m)[indices[1]], rownames(vote_m)[indices[2]])
        rownames(vote_m)[indices[1]] <- namelist[grep("\\(I ", namelist)]
      }

      vote_m <- vote_m[-indices[2], , drop = FALSE]
    }
  }

  # 3. Exclude lopsided votes issues
  # lop_issue expressed as the proportion of non-missing votes on the minority side.
  to_remove <- which(apply(vote_m, 2, function(col) {

    col_no_na <- col[!is.na(col)]
    vote_count <- table(col_no_na)

    if (length(vote_count) == 1) {
      return(TRUE)
    }

    minority_count <- min(vote_count)
    minority_proportion <- minority_count / length(col_no_na)

    return(minority_proportion < lop_issue)
  }))
  if (length(to_remove) > 0) {
    vote_m <- vote_m[, -to_remove]
  }

  # 4. Remove legislators where missing values is greater than or equal to 'lop_leg'
  absent_members <- which(rowMeans(is.na(vote_m)) >= lop_leg)
  if (length(absent_members) > 0) {
    vote_m <- vote_m[-absent_members,]
  }

  # 5. Positive index
  matched_rows <- grep(pos_leg, rownames(vote_m), value = TRUE)
  if (length(pos_leg) != 1) {
    stop("`pos_leg` should contain exactly one legislator.")
  } else if (length(matched_rows) == 0){
    stop(paste("The `pos_leg`:", pos_leg, "is not found in the final vote matrix row names."))
  } else if (length(matched_rows) == 1){
    pos_ind = grep(pos_leg, rownames(vote_m))
  }

  # output
  if (is.matrix(vote_info)) {
    vote_out = vote_m
  } else if (is.list(vote_info) && "votes" %in% names(vote_info)) {
    vote_out = vote_info
    # votes
    vote_out$votes = vote_m
    # codes
    vote_out$codes$yea <- 1
    vote_out$codes$nay <- 0
    if (!is.null(vote_out$codes$notInLegis)) {
      vote_out$codes$notInLegis <- NA
    }
    if (!is.null(vote_out$codes$missing)) {
      vote_out$codes$missing <- NA
    }
    # n and m
    vote_out$n = nrow(vote_out$votes)
    vote_out$m = ncol(vote_out$votes)
    # legis.data
    if (!is.null(vote_out$legis.data)) {
      vote_out$legis.data <- vote_out$legis.data[, !(names(vote_out$legis.data) %in% c("icpsrLegis", "partyCode"))]

      vote_out$legis.data <- vote_out$legis.data[-grep(leg_rm, rownames(vote_out$legis.data)),]
      for (i in seq_along(combine_leg_name)) {
        name <- combine_leg_name[i]
        party <- combine_leg_party[i]
        matching_rows <- grep(name, rownames(vote_out$legis.data))
        matching_party_rows <- grep(paste0("\\(", party), rownames(vote_out$legis.data)[matching_rows])
        rows_to_delete <- matching_rows[-matching_party_rows]
        vote_out$legis.data <- vote_out$legis.data[-rows_to_delete, ]
      }
    }

  }

  # All preparations are done!

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
    pos_ind - 1, verbose)

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
    vote_info = vote_out
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
