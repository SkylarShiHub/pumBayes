#' @title Calculate Probabilities for the IDEAL Model
#' @description This function computes the probability matrix for the IDEAL Model.
#' Specifically, it calculates the probabilities of voting "Yea" for each legislator (member),
#' issue, (and time period) based on the posterior samples of model parameters.
#' @param vote_info A logical vote matrix (or a rollcall object) in which rows represent members and columns represent issues.
#' The entries should be FALSE ("No"), TRUE ("Yes"), or NA (missing data).
#' @param post_samples Posterior samples obtained from function 'ideal' in 'pscl' package.
#' @return An array of probabilities with three dimensions. The first one represents to members, the second one refers to issues,
#' and the third one refers to MCMC iterations.
#' @examples
#' prob_ideal = predict_ideal(post_samples_ideal, h116.c)
#' @export
predict_ideal <- function(vote_info, post_samples) {
  # Check and process input vote object
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

  beta1_cols <- colnames(post_samples$beta[,,1])
  vote_cols <- colnames(vote_m)
  common_cols <- intersect(vote_cols, beta1_cols)
  case_vote_m <- vote_m[, common_cols]
  leg_cols = colnames(post_samples$x[,,1])
  vote_rows <- rownames(case_vote_m)
  common_leg <- intersect(vote_rows, leg_cols)
  case_vote_m <- case_vote_m[common_leg,]

  result = cbind(post_samples$x[,,1], post_samples$beta[,,1], post_samples$beta[,,2])
  n_samples <- nrow(result)
  n_rows <- nrow(case_vote_m)
  n_cols <- ncol(case_vote_m)

  prob_array <- array(NA, dim = c(n_rows, n_cols, n_samples),
                      dimnames = list(rownames(case_vote_m), colnames(case_vote_m), NULL))

  theta_inds <- seq(1, n_rows)
  beta_inds <- seq(length(theta_inds) + 1, length(theta_inds) + n_cols)
  alpha_inds <- seq(length(theta_inds) + n_cols + 1, length(theta_inds) + 2 * n_cols)

  for (k in 1:n_samples) {
    row <- result[k, ]
    case_alpha_m <- matrix(row[alpha_inds], nrow = n_rows, ncol = n_cols, byrow = TRUE)
    case_beta_m <- matrix(row[beta_inds], nrow = n_rows, ncol = n_cols, byrow = TRUE)
    case_theta_m <- matrix(row[theta_inds], nrow = n_rows, ncol = n_cols)

    justice_probs <- pnorm(-case_alpha_m + case_beta_m * case_theta_m)
    justice_probs[justice_probs < 1e-7] <- 1e-7
    justice_probs[justice_probs > (1 - 1e-7)] <- 1 - 1e-7

    prob_array[, , k] <- justice_probs
  }

  return(prob_array)
}
