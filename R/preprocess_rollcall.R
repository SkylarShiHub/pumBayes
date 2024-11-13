#' @title Preprocess Roll Call Data
#' @description This function is used to preprocess roll call data for analysis.
#' It allows users to remove legislators, combine legislators with specified indices, exclude lopsided votes based on minority voting proportions, and filter out legislators with excessive missing votes.
#' @param x A rollcall object
#' @param data_preprocess A list of parameters for preprocessing:
#'   \itemize{
#'     \item `leg_rm` (default = NULL): A vector of indices specifying legislators to be removed.
#'     \item `combine_leg_index` (default = NULL): A list of vectors where each vector specifies the indices of legislators to be combined.
#'     \item `combine_leg_party` (default = NULL): A vector specifying the party affiliations for combined legislators.
#'     \item `lop_leg` (default = 0.6): A threshold indicating the maximum allowable proportion of missing votes for each legislator. Legislators with a proportion of missing votes greater than this value are removed.
#'     \item `lop_issue` (default = 0): A threshold for the proportion of non-missing votes on the minority side. Voting issues with a minority proportion lower than this value are excluded.
#'   }
#' @return A rollcall object that has been processed
#' @examples
#' load("data/h116.rda")
#' h116.c = preprocess_rollcall(h116)
#' @export
preprocess_rollcall <- function(x, data_preprocess=list(leg_rm = NULL, combine_leg_index = NULL,
                                                        combine_leg_party = NULL,
                                                        lop_leg = 0.6, lop_issue = 0)){
  leg_rm = data_preprocess$leg_rm
  combine_leg = data_preprocess$combine_leg_index
  combine_leg_party = data_preprocess$combine_leg_party
  lop_leg = data_preprocess$lop_leg
  lop_issue = data_preprocess$lop_issue

  ## checks before data preprocess
  if (!(is.list(x))){
    stop("x is not a list.")
  }

  if (!("codes" %in% names(x))) {
    stop("x dose not include a 'codes' component.")
  } else if (!("yea" %in% names(x$codes))){
    stop("'yea' dose not exist in the 'codes' component.")
  } else if (!("nay" %in% names(x$codes))){
    stop("'nay' dose not exist in the 'codes' component.")
  }


  # 1. Check input vote object
  if (!("votes" %in% names(x))) {
    stop("x does not include 'votes' component.")
  } else if (all(is.na(x$votes))){
    stop("The 'votes' component in x is entirely NA.")
  }

  vote_m <- x$votes

  # 2. Check for name conflicts
  common_index <- intersect(leg_rm, unlist(combine_leg))
  if (length(common_index) > 0) {
    stop(paste("The following index appear in both leg_rm and combine_leg:",
               common_index))
  }

  # 3. Check for input combined legislators
  if (!(all(is.na(combine_leg)) && all(is.na(combine_leg_party))) &&
      length(combine_leg) != length(combine_leg_party)) {
    stop("Legislator names and parties to be combined must have the same length.")
  }

  unmatched_index <- sapply(combine_leg, function(index) {
    if (all(index <= nrow(vote_m) & index >= 1)) {
      return(NULL)
    } else {
      return(index)
    }
  })
  unmatched_index <- unlist(unmatched_index)
  if (length(unmatched_index) > 0) {
    stop(paste("The following index are not found :", paste(unmatched_index)))
  }

  ## data preprocess
  # 1. Remove legislators
  if (!is.null(leg_rm)) {
    vote_m[leg_rm, ] <- NA
  }

  legis.data.out = x$legis.data
  rm.index = c()
  # 2. combine legislators
  # if there is legis object, do the same operations
  # check the party object : warning party information will not be updated
  if (!is.null(combine_leg)){
    for (i in seq_along(combine_leg)) {

      indices <- combine_leg[[i]]
      rows = matrix(nrow = 0, ncol = ncol(vote_m))

      for (j in 1: length(indices)){
        rows <- rbind(rows, vote_m[indices[j], ])
      }

      combined_row <- apply(rows, 2, function(x) {
        first_non_na <- x[!is.na(x)]
        if (length(first_non_na) > 0) {
          return(first_non_na[1])
        } else {
          return(NA)
        }
      })

      vote_m[indices[2],] = combined_row
      # vote_m = vote_m[-indices[1],]

      # output
      if ("legis.data" %in% names(x)){
        if ("party" %in% names(x$legis.data) && !is.null(x$legis.data$party)){
          rows_to_check <- legis.data.out[unlist(combine_leg[i]), ]
          party_value <- combine_leg_party[i]
          matching_indices <- grep(party_value, rows_to_check$party)
          rm.index[i] = unlist(combine_leg[i])[-matching_indices]
        } else {
          warning("Party information not found and will not be updated.")
        }
      }

    }
    # remove
    vote_m = vote_m[-sapply(combine_leg, function(x) x[1]),]
    if (!is.null(rm.index)){
      legis.data.out = legis.data.out[-rm.index,]
    } else {
      legis.data.out = legis.data.out[-sapply(combine_leg, function(x) x[1]),]
    }
  }

  # 3. Exclude lopsided votes issues
  # lop_issue expressed as the proportion of non-missing votes on the minority side.

  # h116$vote.data exist and not NULL
  to_remove <- which(apply(vote_m, 2, function(col) {

    col_no_na <- col[!is.na(col)]
    vote_count <- table(col_no_na)
    yea_sum <- sum(vote_count[names(vote_count) %in% x$codes$yea])
    nay_sum <- sum(vote_count[names(vote_count) %in% x$codes$nay])
    minority_proportion <- min(yea_sum, nay_sum) / sum(yea_sum, nay_sum)

    return(minority_proportion < lop_issue)
  }))

  if (length(to_remove) > 0) {
    vote_m <- vote_m[, -to_remove]
  }

  # 4. Remove legislators where missing values is greater than or equal to 'lop_leg'

  # legislator part
  vote_m_copy <- vote_m
  vote_m_copy[!(vote_m_copy %in% c(x$codes$yea, x$codes$nay))] <- NA
  absent_members <- which(rowMeans(is.na(vote_m_copy)) > lop_leg)
  if (length(absent_members) > 0) {
    vote_m <- vote_m[-absent_members,]
    legis.data.out = legis.data.out[-absent_members,]
  }


  xout = list(votes = vote_m,
              codes = x$codes,
              n = dim(vote_m)[1],
              m = dim(vote_m)[2])

  if ("legis.data" %in% names(x)) {
    xout$legis.data = legis.data.out
  }

  if ("vote.data" %in% names(x) && !is.null(x$vote.data)){
    xout$vote.data = x$vote.data
  }

  if ("desc" %in% names(x) && !is.null(x$desc)){
    xout$desc = x$desc
  }

  if ("source" %in% names(x) && !is.null(x$source)){
    xout$source = x$source
  }
  class(xout) <- class(x)
  return(xout)
}

