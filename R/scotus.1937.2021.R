#' U.S. Supreme Court Voting Data (1937-2021)
#'
#' This dataset contains voting records from the U.S. Supreme Court between 1937 and 2021.
#' Loading `data(scotus.1937.2021)` will load the following two independent objects into the environment:
#'
#' @docType data
#' @name scotus.1937.2021
#'
#' @format The dataset consists of the following two objects:
#' \describe{
#'   \item{mqVotes}{A `48 × 6108` matrix, where each row represents a judge and each column represents a case. Entries are:
#'     \itemize{
#'       \item `1` (`TRUE`): The judge voted to reverse the lower court decision.
#'       \item `0` (`FALSE`): The judge voted to uphold the lower court decision.
#'       \item `NA`: The judge did not vote on the case.
#'     }
#'   }
#'   \item{mqTime}{A numeric vector of length `6108`, indicating the time period associated with each case.}
#' }
#'
#' @usage data(scotus.1937.2021)
#'
#' @source The data were obtained from the Martin-Quinn Scores Database, maintained by Washington University in St. Louis.
#'   The dataset can be accessed and downloaded from \url{http://mqscores.wustl.edu/replication.php}.
#'
#' @examples
#' data(scotus.1937.2021)
#' str(mqVotes)
#' str(mqTime)
#'
#' @keywords datasets
NULL

