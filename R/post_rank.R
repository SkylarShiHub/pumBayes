#' @title Generate Quantile Ranks for Legislators
#' @description This function calculates quantile ranks for each legislator based on posterior samples of beta parameters from MCMC.
#' The function can handle any specified quantiles, such as median (0.5), and is flexible to support other quantiles provided as input.
#' @param vote_info The same vote object used for getting posterior samples.
#' @param beta A matrix of posterior samples of beta obtained from MCMC, with columns representing legislators.
#' @param quantiles A numeric vector specifying the quantiles to be calculated for the ranks (default is `c(0.5)` for median rank).
#' @return A data frame containing the legislators' names, party affiliations, states, and their ranks at each specified quantile.
#' If the median is included, it will be named `median` in the output. The output data frame is sorted in ascending order based on the values in the median column.
#' @examples
#' rank_results = post_rank(vote_info = h116.c, beta = post_samples_pum$beta, quantiles = c(0.5))
#' @export

post_rank = function(vote_info, beta, quantiles = c(0.5)){

  rank_matrix <- apply(beta, 1, function(x) rank(-x)) # Descending ranking
  rank_matrix <- t(rank_matrix)

  rank_quantiles <- apply(rank_matrix, 2, function(x) {
    quantiles_result <- quantile(x, probs = quantiles, names = FALSE)
    rounded_result <- round(quantiles_result)
    positive_result <- ifelse(rounded_result > 0, rounded_result, NA)
    return(positive_result)
  })
  if (length(quantiles) == 1) {
    rank_quantiles <- t(matrix(rank_quantiles))
  }

  if (is.matrix(vote_info)) {
    name_vector = rownames(vote_info)
  } else if (is.list(vote_info) && "legis.data" %in% names(vote_info)) {
    name_vector = rownames(vote_info$legis.data)
  }

  rank_data <- data.frame(name = name_vector, t(rank_quantiles))
  colnames(rank_data)[-1] <- quantiles

  rank_data$`0.5` <- rank(rank_data$`0.5`, ties.method = "first")
  rownames(rank_data) = NULL
  return(rank_data)
}
