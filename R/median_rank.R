#' @title Get median rank plots
#' @description Generates the median rank plot for all members based on posterior samples
#' @param beta A matrix of posterior samples of beta obtained from MCMC
#' @param col_blue A character string indicating whether the "D" party should be colored blue or "R".
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @return A list containing a ggplot object of the rank plot and a data frame of sorted rank data.
#' @examples
#' rank_results = median_rank(beta = post_samples$beta, col_blue = "D")
#' @export
median_rank = function(beta, col_blue = "D"){
  # Create a matrix to store ranks
  rank_matrix <- apply(beta, 1, function(x) rank(-x)) # Descending ranking
  rank_matrix <- t(rank_matrix)

  # Find the median rank of each member
  rank_median <- apply(rank_matrix, 2, median)

  median_value <- median(rank_median)
  posterior_prob_median <- apply(rank_matrix, 2, function(x) {
    mean(x == median_value)
  })
  # 5 or 6 legislators with the highest posterior probabilities of having the median rank
  top_prob_legislators <- order(posterior_prob_median, decreasing = TRUE)[1:5]
  top_legislators_names <- sub("_.*", "", colnames(beta)[top_prob_legislators])

  all_members <- gsub("_beta", "", colnames(beta))
  party_vector <- ifelse(grepl("\\(D", all_members), "D",
                         ifelse(grepl("\\(R", all_members), "R",
                                ifelse(grepl("\\(I", all_members), "I", NA)))
  name_vector <- sub(" \\(.*\\)", "", all_members)
  rank_data <- data.frame(names = name_vector, party = party_vector, median = rank_median)

  # Sort the rank data
  sorted_rank_data <- rank_data[order(-rank_data$median), ]
  sorted_rank_data$party <- factor(sorted_rank_data$party)
  sorted_rank_data$rank <- rank(-sorted_rank_data$median, ties.method = "first")
  rownames(sorted_rank_data) <- gsub("_beta", "", rownames(sorted_rank_data))

  if (col_blue == "D") {
    party_colors <- c("D" = "#A6CEE3", "R" = "#FB9A99", "I" = "black") # Set colors matching "Paired"
    party_shapes <- c("D" = 19, "R" = 17, "I" = 18)
  } else if (col_blue == "R") {
    party_colors <- c("D" = "#FB9A99", "R" = "#A6CEE3", "I" = "black") # Reverse for Republicans as blue
    party_shapes <- c("D" = 17, "R" = 19, "I" = 18)
  }

  top_labels <- sorted_rank_data[1:5, ]
  bottom_labels <- sorted_rank_data[(nrow(sorted_rank_data) - 4):nrow(sorted_rank_data), ]
  median_index <- which(sorted_rank_data$rank == floor(median(sorted_rank_data$rank)) | sorted_rank_data$rank == ceiling(median(sorted_rank_data$rank)))
  median_label <- sorted_rank_data[median_index, , drop = FALSE]

  rank_plot = ggplot(sorted_rank_data, aes(x = median, y = rank, color = party, shape = party)) +
    geom_point(size = 2, alpha = 0.6) +
    geom_text_repel(data = top_labels, aes(label = names),
                    size = 3, col = "black") +
    geom_text_repel(data = bottom_labels, aes(label = names),
                    size = 3, col = "black") +
    geom_text_repel(data = median_label, aes(label = names),
                    size = 3, col = "black") +
    scale_color_manual(values = party_colors) +
    scale_shape_manual(values = party_shapes) +
    # scale_color_manual(values = party_colors) +
    labs(x = "Median", y = "Rank") +
    scale_y_reverse() +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  return(list(rank_plot = rank_plot, rank_data = sorted_rank_data, top5_prob_median = top_legislators_names))
}
