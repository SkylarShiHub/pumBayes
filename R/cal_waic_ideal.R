#' @title Calculate WAIC of the IDEAL model
#' @description This function is used to get the WAIC value of the IDEAL model
#' @param chain_results Posterior samples obtained from function 'ideal' in 'pscl' package.
#' @param vote_info Rollcall object used as input of the 'ideal' function.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return The WAIC value of IDEAL model
#' @examples
#' waic_ideal = cal_waic_ideal(post_samples_ideal, h116)
#' @export

cal_waic_ideal <- function(
    chain_results, vote_info, block = F, block_vote = F) {

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

  beta1_cols <- colnames(chain_results$beta[,,1])
  vote_cols <- colnames(vote_m)
  common_cols <- intersect(vote_cols, beta1_cols)
  case_vote_m <- vote_m[, common_cols]
  leg_cols = colnames(chain_results$x[,,1])
  vote_rows <- rownames(case_vote_m)
  common_leg <- intersect(vote_rows, leg_cols)
  case_vote_m <- case_vote_m[common_leg,]

  if (block) {
    num_votes = nrow(case_vote_m)
  } else if (block_vote) {
    num_votes = ncol(case_vote_m)
  } else {
    num_votes = sum(!is.na(as.vector(case_vote_m)))
  }
  mean_prob <- rep(-Inf, num_votes)
  mean_log_prob <- rep(0, num_votes)
  log_prob_var <- rep(0, num_votes)
  num_iter = 0
  result = cbind(chain_results$x[,,1], chain_results$beta[,,1], chain_results$beta[,,2])
  theta_inds <- seq(1, nrow(case_vote_m))
  beta_inds <- seq(length(theta_inds)+1,length(theta_inds) + ncol(case_vote_m))
  alpha_inds <- seq(length(theta_inds)+ncol(case_vote_m)+1,length(theta_inds)+2*ncol(case_vote_m))
  for (j in 1:nrow(result)) {
      row <- result[j,]
      case_alpha_m <-
        matrix(row[alpha_inds], nrow = nrow(case_vote_m),
               ncol = ncol(case_vote_m), byrow = T)
      case_beta_m <-
        matrix(row[beta_inds], nrow = nrow(case_vote_m),
               ncol = ncol(case_vote_m), byrow = T)
      case_theta_m <-
        matrix(row[theta_inds], nrow = nrow(case_vote_m),
               ncol = ncol(case_vote_m))
      justice_probs <-
        pnorm(-case_alpha_m  + case_beta_m * case_theta_m)
      justice_probs[justice_probs < 1e-9] =
        1e-9
      justice_probs[justice_probs > (1 - 1e-9)] =
        1 - 1e-9
      log_prob <- case_vote_m * log(justice_probs) +
        (1 - case_vote_m) * log(1 - justice_probs)
      if (block) {
        log_prob <- rowSums(log_prob, na.rm = T)
      } else if (block_vote) {
        log_prob <- colSums(log_prob, na.rm = T)
      } else {
        log_prob <- log_prob[!is.na(log_prob)]
      }
      mean_prob <- pmax(mean_prob, log_prob) + log(1 + exp(pmin(mean_prob, log_prob)
                                                           - pmax(mean_prob, log_prob)))
      next_mean_log_prob = (num_iter * mean_log_prob + log_prob) / (num_iter + 1)
      log_prob_var = log_prob_var +
        (log_prob - mean_log_prob) * (log_prob - next_mean_log_prob)
      mean_log_prob = next_mean_log_prob
      num_iter = num_iter + 1
      delta <- log_prob - mean_log_prob
      mean_log_prob <- mean_log_prob + delta / num_iter
      log_prob_var <- log_prob_var + delta * (log_prob - mean_log_prob)
    }
  return(sum(mean_prob - log(num_iter) - (log_prob_var) / (num_iter - 1)))
}
