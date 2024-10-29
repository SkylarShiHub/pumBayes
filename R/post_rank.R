#' @title Generate Quantile Ranks for Legislators
#' @description This function calculates quantile ranks for each legislator based on posterior samples of beta parameters from MCMC.
#' The function can handle any specified quantiles, such as median (0.5), and is flexible to support other quantiles provided as input.
#' @param beta A matrix of posterior samples of beta obtained from MCMC, with columns representing legislators.
#' @param quantiles A numeric vector specifying the quantiles to be calculated for the ranks (default is `c(0.5)` for median rank).
#' @return A data frame containing the legislators' names, party affiliations, states, and their ranks at each specified quantile.
#' If the median is included, it will be named `median` in the output. The output data frame is sorted in ascending order based on the values in the median column.
#' @examples
#' rank_results = post_rank(beta = post_samples$beta, quantiles = c(0.5))
#' @export

post_rank = function(beta, quantiles = c(0.5)){

  rank_matrix <- apply(beta, 1, function(x) rank(-x)) # Descending ranking
  rank_matrix <- t(rank_matrix)

  rank_quantiles <- apply(rank_matrix, 2, function(x) {
    quantile(x, probs = quantiles, names = FALSE)
  })
  if (length(quantiles) == 1) {
    rank_quantiles <- t(matrix(rank_quantiles))
  }
  all_members <- gsub("_beta", "", colnames(beta))
  party_vector <- ifelse(grepl("\\(D", all_members), "D",
                         ifelse(grepl("\\(R", all_members), "R",
                                ifelse(grepl("\\(I|Indep", all_members), "I", NA)))
  name_vector <- sub(" \\(.*\\)", "", all_members)
  state_vector <- sapply(all_members, function(x) {
    state_info <- sub(".*(.{4})\\)", "\\1", x)
    return(state_info)
  })
  rank_data <- data.frame(name = name_vector, party = party_vector, state = state_vector, t(rank_quantiles))
  colnames(rank_data)[-(1:3)] <- quantiles
  colnames(rank_data) <- gsub("^0.5$", "median", colnames(rank_data))

  sorted_rank_data <- rank_data[order(rank_data$median), ]
  sorted_rank_data$rank <- rank(sorted_rank_data$median, ties.method = "first")
  rownames(sorted_rank_data) <- sorted_rank_data$rank
  return(sorted_rank_data)
}
