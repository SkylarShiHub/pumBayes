#' @title Calculate Probabilities for Probit Unfolding Models
#' @description This function computes the probability matrix for both static and dynamic Probit Unfolding Models.
#' Specifically, it calculates the probabilities of voting "Yea" for each legislator (member),
#' issue, (and time period) based on the posterior samples of model parameters.
#' @param vote_info A logical vote matrix (or a rollcall object) in which rows represent members and columns represent issues.
#' The entries should be FALSE ("No"), TRUE ("Yes"), or NA (missing data).
#' @param years_v A vector representing the time period for each vote in the model. This is defultly set as `NULL` for a static model.
#' @param post_samples A list of posterior samples of parameters obtained from MCMC.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return An array of probabilities with three dimensions. The first one represents to members, the second one refers to issues,
#' and the third one refers to MCMC iterations.
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
#' h116.c.pum.predprob = predict_pum(h116.c, years_v = NULL, h116.c.pum)
#' }
#' @export
predict_pum <- function(vote_info, years_v = NULL, post_samples) {
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

  if (is.null(years_v)){
    block_m = cbind(1:nrow(vote_m), 1)
    years_v = rep(1,ncol(vote_m))
  } else{
    block_m <- do.call(rbind, lapply(1:nrow(vote_m), function(i) {
      interested_inds <- which(!is.na(vote_m[i,]))
      cbind(i, unique(sort(years_v[interested_inds])))
    }))
  }

  leg_pos <- as.matrix(post_samples$beta)
  match_index <- grep("RHJackson_beta_9", colnames(leg_pos))
  if (length(match_index) > 0) {
    leg_pos <- leg_pos[, -match_index]
  }

  alpha_draws0 <- as.matrix(cbind(post_samples$alpha1, post_samples$alpha2))
  alpha_draws <- t(apply(alpha_draws0, 1, function(row) {
    as.vector(t(matrix(row, ncol = 2)))
  }))
  rm(alpha_draws0)
  delta_draws0 <- as.matrix(cbind(post_samples$delta1, post_samples$delta2))
  delta_draws <- t(apply(delta_draws0, 1, function(row) {
    as.vector(t(matrix(row, ncol = 2)))
  }))
  rm(delta_draws0)

  n_samples <- nrow(leg_pos)
  n_rows <- nrow(vote_m)
  n_cols <- ncol(vote_m)
  prob_array <- array(NA, dim = c(n_rows, n_cols, n_samples),
                      dimnames = list(rownames(vote_m), colnames(vote_m), NULL))

  for (iter in 1:n_samples) {

    for (ind in 1:nrow(block_m)) {
      i <- block_m[ind, 1]
      year <- block_m[ind, 2]

      interested_cases <- which(years_v == year)
      for (j in interested_cases) {
        if (is.na(vote_m[i, j])) {
          next
        }

        mean_1 <- alpha_draws[iter, 2 * j - 1] *
          (leg_pos[iter, ind] - delta_draws[iter, 2 * j - 1])
        mean_2 <- alpha_draws[iter, 2 * j] *
          (leg_pos[iter, ind] - delta_draws[iter, 2 * j])

        yea_prob <- bvnd(-mean_1 / sqrt(2), -mean_2 / sqrt(2), 0.5)
        yea_prob <- pmin(pmax(yea_prob, 1e-9), 1 - 1e-9)

        prob_array[i, j, iter] <- yea_prob
      }
    }
  }
  return(prob_array)
}
