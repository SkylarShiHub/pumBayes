#' 116th U.S. House of Representatives Roll Call Votes
#'
#' This dataset contains roll call voting records from the 116th U.S. House of Representatives.
#' The data was obtained using the `readKH()` function from Voteview, which reads roll call
#' vote data from the Voteview database.
#'
#' @docType data
#' @name h116
#'
#' @format A list with 8 elements:
#' \describe{
#'   \item{votes}{A `452 × 952` matrix of roll call votes, where each row represents a legislator and each column represents a vote.}
#'   \item{codes}{A list containing vote codes:}
#'   \describe{
#'     \item{yea}{Codes representing 'Yes' votes.}
#'     \item{nay}{Codes representing 'No' votes.}
#'     \item{notInLegis}{Codes for members not in the legislature.}
#'     \item{missing}{Codes for missing votes.}
#'   }
#'   \item{n}{Integer, number of legislators (452).}
#'   \item{m}{Integer, number of votes (952).}
#'   \item{legis.data}{A data frame (`452 × 6`) containing legislator information:}
#'   \describe{
#'     \item{state}{State abbreviation of each legislator.}
#'     \item{icpsrState}{ICPSR state code.}
#'     \item{cd}{Congressional district.}
#'     \item{icpsrLegis}{ICPSR legislator ID.}
#'     \item{party}{Party affiliation (`"D"`, `"R"`, `"I"`.).}
#'     \item{partyCode}{Numerical party code (`200` for Democrats, `100` for Republicans, `328` for Independents.).}
#'   }
#'   \item{vote.data}{Currently NULL (reserved for additional vote metadata).}
#'   \item{desc}{Description: `"116th U.S. House of Representatives"`.}
#'   \item{source}{URL for the original data source.}
#' }
#'
#' @usage data(h116)
#'
#' @source Jeffrey B. Lewis, Keith Poole, Howard Rosenthal, Adam Boche, Aaron Rudkin, and Luke Sonnet.
#'   Voteview: Congressional roll-call votes database. \url{https://voteview.com/}, 2024. Accessed: 2024-07-15.
#'
#' @examples
#' data(h116)
#' str(h116)
#'
#' @keywords datasets
NULL
