#' @title Generate Data for Item Characteristic Curves
#' @description This function calculates the data needed to plot the item characteristic
#' curve for a specific issue based on posterior samples.
#' @param vote_num The vote number of the issue to be reviewed. This refers to numbers
#' in the column names of the input vote matrix, not the clerk session vote number.
#' @param post_samples A list of posterior samples of parameters obtained from `sample_pum_static` in `pumBayes`.
#' @importFrom Rcpp sourceCpp
#' @useDynLib pumBayes
#' @return A data frame containing `beta_samples`, mean probabilities (`means`),
#' and confidence intervals (`ci_lower` and `ci_upper`) for the input issue,
#' which can be used to plot the item characteristic curve.
#' @examples
#' item_data <- item_char(vote_num = 5, post_samples = post_samples)
#' @export
item_char = function(vote_num, post_samples){
  beta = as.matrix(post_samples$beta)
  alpha1 = as.matrix(post_samples$alpha1)
  alpha2 = as.matrix(post_samples$alpha2)
  delta1 = as.matrix(post_samples$delta1)
  delta2 = as.matrix(post_samples$delta2)

  col_index <- grep(paste0("(^|\\s)", vote_num , "_"), colnames(alpha1))
  beta_samples = seq(min(beta), max(beta), length.out = 500)
  prob_mat = matrix(nrow = nrow(alpha1), ncol = length(beta_samples))

  for (i in (1:length(beta_samples))){
    term1 <- -alpha1[,col_index] * (beta_samples[i] - delta1[,col_index]) / sqrt(2)
    term2 <- -alpha2[,col_index] * (beta_samples[i] - delta2[,col_index]) / sqrt(2)
    prob_mat[,i] <- bvndvec(term1, term2, rep(0.5,length(alpha1[,col_index])))
  }
  means <- colMeans(prob_mat)
  ci_lower <- apply(prob_mat, 2, quantile, probs = 0.05)
  ci_upper <- apply(prob_mat, 2, quantile, probs = 0.95)

  plot_data <- data.frame(
    beta_samples = beta_samples,
    means = means,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
  return(plot_data)
}
