#' @title Calculate Probabilities for Dynamic Item Response Theory Model
#' @description This function computes the probability matrix for a dynamic item response theory (IRT) model.
#' Specifically, it calculates the probabilities of voting "Yea" for each legislator (member),
#' issue, and time period based on the posterior samples of model parameters.
#' @param post_samples MCMC results obtained from `wnominate' function in `wnominate' package.
#' @param vote_info A logical vote matrix where rows represent members and columns represent issues.
#' The entries should be FALSE ("No"), TRUE ("Yes"), or NA (missing data).
#' @param years_v A vector representing the time period for each vote in the model.
#' @examples
#' prob_irt_dy = predict_irt(post_samples_irt_dy, mqVotes, mqTime)
#' @export
predict_irt <- function(vote_info, years_v, post_samples) {

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

  n_rows <- nrow(vote_m) # number of legislators
  n_cols <- ncol(vote_m) # number of issues
  n_samples <- nrow(post_samples) # number of samples

  # Identify alpha and beta indices
  alpha_inds <- grep("alpha", colnames(post_samples))
  beta_inds <- grep("beta", colnames(post_samples))

  # Initialize 3D array for probabilities
  prob_array <- array(NA, dim = c(n_rows, n_cols, n_samples),
                      dimnames = list(rownames(vote_m), colnames(vote_m), paste0("Sample_", 1:n_samples)))

  for (j in 1:n_samples) {
    row <- post_samples[j,]

    # Generate alpha and beta matrices
    case_alpha_m <- matrix(row[alpha_inds], nrow = n_rows, ncol = n_cols, byrow = TRUE)
    case_beta_m <- matrix(row[beta_inds], nrow = n_rows, ncol = n_cols, byrow = TRUE)

    # Create ideology matrix, combining members and time
    ideology_m <- matrix(NA, nrow = n_rows, ncol = n_cols)
    for (k in 1:n_rows) {
      judge <- rownames(vote_m)[k]
      judge_ind_v <- grep(paste("theta", judge, "", sep = "."), colnames(post_samples))

      time_ind <- sapply(
        strsplit(colnames(post_samples)[judge_ind_v], "\\."), function(j_ind) {
          as.numeric(gsub("t", "", j_ind[3]))
        })

      ideology_m[k, which(years_v %in% time_ind)] <-
        rep(row[judge_ind_v], table(years_v[years_v %in% time_ind]))
    }

    # Calculate justice probabilities
    justice_probs <- pnorm(-case_alpha_m + case_beta_m * ideology_m)
    justice_probs[justice_probs < 1e-9] <- 1e-9
    justice_probs[justice_probs > (1 - 1e-9)] <- 1 - 1e-9

    # Store probabilities in array
    prob_array[, , j] <- justice_probs
  }

  return(prob_array)
}



