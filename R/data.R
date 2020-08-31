#' @title US House elections dataset from Lee (2008)
#'
#' @description Data on elections to the U.S. House of Representatives
#'
#' @details Lee (2008) uses the data on elections to the U.S. House of
#'   Representatives  (1946-1998) to analyze the incumbency advantage. For a
#'   party's candidate who barely won the election (and became incumbent) or
#'   barely lost the previous election, their districts usually share many
#'   common characteristics. The electoral success for the two groups of common
#'   characteristics. The electoral success for the two groups of candidates in
#'   the next election can be used to identify the causal incumbency advantage.
#'   We subtract the number 50, equivalent to fifty percent, from the vote share
#'   data so that the cutoff is 0 for treatment assignment. Each row of the data
#'   frame represents a party's candidate.
#'
#' @format A data frame with 6,558 rows and 2 variables. \describe{
#'   \item{\code{voteshare}}{the treatment outcome variable \code{y}, vote
#'   share in the next election} \item{\code{margin}}{the covariate \code{x},
#'   democratic margin of victory at the previous election} } Each row
#'
#' @source Mostly Harmless Econometrics data archive at
#'   \url{https://economics.mit.edu/faculty/angrist/data1/mhe}
#'
#' @references{
#'
#' \cite{Lee, D. S. (2008) "Randomized experiments from non-random selection in
#' U.S. House elections," Journal of Econometrics, 142 (2), 675-697.}
#'
#' }
"lee"

#' @title The Head Start dataset from Ludwig and Miller (2007)
#'
#' @description The Head Start policy dataset on the mortality rates for
#'   children
#'
#' @details Ludwig and Miller (2007) studies the impact of Head Start funding on
#'   children's schooling and health. This subset of their data focuses on the
#'   change in mortality rates for children due to the Head Start funding cutoff
#'   set by the Office of Economic Opportunities (OEF). Changes in morality
#'   rates near the funding cutoff can be used to identify the causal effect of
#'   such policy. We subtract the average poverty rate (59.198) from the
#'   poverty rate data so that the cutoff is 0 for treatment assignment. Each
#'   row in the data frame represents a county.
#'
#' @format A data frame with 3,103 rows and 2 variables: \describe{
#'   \item{\code{mortality}}{the treatment outcome variable \code{y}, morality
#'   rates per 10,000 children between 5 and 9 years old between 1973-1983.}
#'   \item{\code{poverty}}{the covariate \code{x}, a county's poverty rate in
#'   1960 relative to the 300th poorest county (59.198)} }
#'
#' @references{
#'
#' \cite{Ludwig, J. and D. L. MIller (2007) "Does Head Start improve children's
#' life chances? Evidence from a regression discontinuity design," Quarterly
#' Journal of Economics, 122 (1), 159-208.}
#'
#' }
"headstart"

