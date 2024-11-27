#' @title Calculate a block version of Watanabe-Akaike Information Criterion (WAIC)
#' @description This function is used to get the WAIC value when blocking members
#' @param vote_info A logical vote matrix (or a rollcall object) in which rows represent members and columns represent issues.
#' The entries should be FALSE ("No"), TRUE ("Yes"), or NA (missing data).
#' @param years_v A vector representing the time period for each vote in the model. This is defultly set as `NULL` for a static model.
#' @param prob_array An array of probabilities with three dimensions.
#' @return The block WAIC value for a static PUM or a vector of WAIC by time for a dynamic PUM.
#' @examples
#' waic_pum_dy = calc_waic(prob_pum_dy, years_v = NULL, h116.c)
#' @export
calc_waic <- function(vote_info, years_v = NULL, prob_array) {
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

  n_rows <- dim(prob_array)[1]
  n_cols <- dim(prob_array)[2]
  n_samples <- dim(prob_array)[3]

  issue_cols <- colnames(prob_array[,,1])
  vote_cols <- colnames(vote_m)
  common_cols <- intersect(vote_cols, issue_cols)
  case_vote_m <- vote_m[, common_cols]
  leg_cols = rownames(prob_array[,,1])
  vote_rows <- rownames(case_vote_m)
  common_leg <- intersect(vote_rows, leg_cols)
  case_vote_m <- case_vote_m[common_leg,]

  block_leg = TRUE
  block_issue = FALSE

  if (is.null(years_v)){
    # static
    if (block_leg) {
      num_votes = nrow(case_vote_m)
    } else if (block_issue) {
      num_votes = ncol(case_vote_m)
    } else {
      num_votes = sum(!is.na(as.vector(case_vote_m)))
    }

    mean_prob <- rep(-Inf, num_votes)
    mean_log_prob <- rep(0, num_votes)
    log_prob_var <- rep(0, num_votes)
    num_iter = 0

    for (j in 1:n_samples){
      log_prob <- case_vote_m * log(prob_array[,,j]) +
        (1 - case_vote_m) * log(1 - prob_array[,,j])
      if (block_leg) {
        log_prob <- rowSums(log_prob, na.rm = T)
      } else if (block_issue) {
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
    return(-2*sum(mean_prob - log(num_iter) - (log_prob_var) / (num_iter - 1)))

  } else {
    # dynamic
    # Initialize WAIC components
    num_points <- sum(apply(case_vote_m, 1, function(row) {
      interested_inds <- which(!is.na(row))
      length(table(years_v[interested_inds]))
    }))
    # mean_prob <- rep(0, num_points)
    mean_prob = rep(-Inf, num_points)
    mean_log_prob <- rep(0, num_points)
    log_prob_var <- rep(0, num_points)
    num_iter <- 0

    for (j in 1:n_samples) {
      log_prob <- case_vote_m * log(prob_array[, , j]) +
        (1 - case_vote_m) * log(1 - prob_array[, , j])

      # Group log probabilities by year
      data_prob_m <- as.data.frame(cbind(years_v, t(log_prob)))
      split_data <- split(data_prob_m[-1], data_prob_m$years_v)
      log_prob <- lapply(split_data, function(df) {
        apply(df, 2, function(x) if (all(is.na(x))) NA else sum(x, na.rm = TRUE))
      })
      log_prob <- do.call(rbind, log_prob)

      # Filter non-NA rows
      non_na_rows <- which(!is.na(log_prob), arr.ind = TRUE)[, "row"]
      years_vector <- as.numeric(rownames(log_prob)[non_na_rows])
      log_prob <- log_prob[!is.na(log_prob)]

      # Update WAIC components
      if (num_iter == 0){
        mean_prob = log_prob
      }
      # mean_prob <- mean_prob + exp(log_prob)
      mean_prob = pmax(mean_prob, log_prob) +
        log(1 + exp(pmin(mean_prob, log_prob) - pmax(mean_prob, log_prob)))
      next_mean_log_prob <- (num_iter * mean_log_prob + log_prob) / (num_iter + 1)
      log_prob_var <- log_prob_var + (log_prob - mean_log_prob) * (log_prob - next_mean_log_prob)
      mean_log_prob <- next_mean_log_prob
      num_iter <- num_iter + 1
    }

    # Calculate WAIC
    waic_result <- -2 * (mean_prob -log(num_iter) - (log_prob_var) / (num_iter - 1))
    # Aggregate WAIC by year
    waic_by_year <- tapply(waic_result, years_vector, sum, na.rm = TRUE)

    return(waic_by_year)

  }
}

calc_waic_dynamic <- function(prob_array, vote_m, years_v) {
  n_rows <- dim(prob_array)[1]
  n_cols <- dim(prob_array)[2]
  n_samples <- dim(prob_array)[3]

  # 初始化 WAIC 组件
  num_points <- sum(apply(vote_m, 1, function(row) {
    interested_inds <- which(!is.na(row))
    length(table(years_v[interested_inds]))
  }))
  mean_prob <- rep(0, num_points)
  mean_log_prob <- rep(0, num_points)
  log_prob_var <- rep(0, num_points)
  num_iter <- 0

  # 遍历样本，更新 WAIC 组件
  for (j in 1:n_samples) {
    log_prob <- vote_m * log(prob_array[, , j]) +
      (1 - vote_m) * log(1 - prob_array[, , j])

    # 按年份分组
    data_prob_m <- as.data.frame(cbind(years_v, t(log_prob)))
    split_data <- split(data_prob_m[-1], data_prob_m$years_v)
    log_prob <- lapply(split_data, function(df) {
      apply(df, 2, function(x) if (all(is.na(x))) NA else sum(x, na.rm = TRUE))
    })
    log_prob <- do.call(rbind, log_prob)

    # 过滤非 NA 行
    non_na_rows <- which(!is.na(log_prob), arr.ind = TRUE)[, "row"]
    years_vector <- as.numeric(rownames(log_prob)[non_na_rows])
    log_prob <- log_prob[!is.na(log_prob)]

    # 更新 WAIC 组件
    # mean_prob <- mean_prob + exp(log_prob)
    mean_prob <- pmax(mean_prob, log_prob) +
      log(1 + exp(pmin(mean_prob, log_prob) - pmax(mean_prob, log_prob)))
    next_mean_log_prob <- (num_iter * mean_log_prob + log_prob) / (num_iter + 1)
    log_prob_var <- log_prob_var + (log_prob - mean_log_prob) * (log_prob - next_mean_log_prob)
    mean_log_prob <- next_mean_log_prob
    num_iter <- num_iter + 1
  }

  # 计算 WAIC
  # waic_result <- -2 * (log(mean_prob / num_iter) - (log_prob_var) / (num_iter - 1))
  waic_result <- -2 * (mean_prob -log(num_iter) - (log_prob_var) / (num_iter - 1))
  # 按年份汇总 WAIC
  waic_by_year <- tapply(waic_result, years_vector, sum, na.rm = TRUE)

  return(waic_by_year)
}
