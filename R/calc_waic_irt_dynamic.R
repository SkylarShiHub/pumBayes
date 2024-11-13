#' @title Calculate WAIC of the dynamic ideal point model
#' @description This function is used to get the WAIC value of the dynamic ideal point model prorosed by Martin and Quinn
#' @param vote_m A vote matrix in which rows represent members and columns represent issues.
#' The entries can only be 0 (indicating ”No”), 1 (indicating ”Yes”), or NA (indicating missing data).
#' @param years_v A vector representing the time period for each vote in the model
#' @param chain_run The posterior samples obtained from `MCMCdynamicIRT1d` in `MCMCpack`
#' @useDynLib pumBayes
#' @examples
#' waic_irt_dy = calc_waic_irt_dynamic(mqVotes, mqTime, post_samples_irt_dy)
#' @export
calc_waic_irt_dynamic <- function(
    case_vote_m, case_year_v, chain_results) {

  num_points <- sum(apply(case_vote_m, 1, function(row) {
    interested_inds <- which(!is.na(row))
    length(table(case_year_v[interested_inds]))
  }))
  mean_prob <- rep(0, num_points)
  mean_log_prob <- rep(0, num_points)
  log_prob_var <- rep(0, num_points)
  num_iter = 0
  result = chain_results
  alpha_inds <- grep("alpha", colnames(result))
  beta_inds <- grep("beta", colnames(result))
  for (j in 1:nrow(result)) {
    row <- result[j,]
    case_alpha_m <-
      matrix(row[alpha_inds], nrow = nrow(case_vote_m),
             ncol = ncol(case_vote_m), byrow = T)
    case_beta_m <-
      matrix(row[beta_inds], nrow = nrow(case_vote_m),
             ncol = ncol(case_vote_m), byrow = T)
    ideology_m <- matrix(NA,  nrow = nrow(case_vote_m),
                         ncol = ncol(case_vote_m))
    for (k in 1:nrow(case_vote_m)) {
      judge <- rownames(case_vote_m)[k]
      judge_ind_v <- grep(paste("theta", judge, "", sep = "."),
                          colnames(result))
      time_ind <- sapply(
        strsplit(colnames(result)[judge_ind_v], "\\."), function(j_ind) {
          as.numeric(gsub("t", "", j_ind[3]))
        })
      ideology_m[k, which(case_year_v %in% time_ind)] <-
        rep(row[judge_ind_v], table(case_year_v[case_year_v %in% time_ind]))
    }
    justice_probs <-
      pnorm(-case_alpha_m  + case_beta_m * ideology_m)

    justice_probs[justice_probs < 1e-9] =
      1e-9
    justice_probs[justice_probs > (1 - 1e-9)] =
      1 - 1e-9
    log_prob <- case_vote_m * log(justice_probs) +
      (1 - case_vote_m) * log(1 - justice_probs)

    data_prob_m <- as.data.frame(cbind(case_year_v, t(log_prob)))
    split_data <- split(data_prob_m[-1], data_prob_m$case_year_v)
    log_prob <- lapply(split_data, function(df) {
      apply(df, 2, function(x) if (all(is.na(x))) NA else sum(x, na.rm = TRUE))
    })

    log_prob <- do.call(rbind, log_prob)
    non_na_rows <- which(!is.na(log_prob), arr.ind = TRUE)[, "row"]
    years_vector <- as.numeric(rownames(log_prob)[non_na_rows])
    log_prob <- log_prob[!is.na(log_prob)]
    mean_prob <- mean_prob + exp(log_prob)
    next_mean_log_prob = (num_iter * mean_log_prob + log_prob) / (num_iter + 1)
    log_prob_var = log_prob_var +
      (log_prob - mean_log_prob) * (log_prob - next_mean_log_prob)
    mean_log_prob = next_mean_log_prob
    num_iter = num_iter + 1
  }

  waic_result <- -2 * (log(mean_prob / num_iter) - (log_prob_var) / (num_iter - 1))

  # Sum WAIC by year
  waic_by_year <- tapply(waic_result, years_vector, sum, na.rm = TRUE)

  return(waic_by_year)
}
