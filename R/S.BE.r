#' @export
#'
#' @title
#' Bernoulli Sampling
#' @description
#' Draws a Bernoulli sample from a finite population of size \code{N}.
#' Each unit is independently selected with the same inclusion probability
#' \code{prob}.
#' @return
#' A vector of length \code{N} where selected units contain their population
#' index and non-selected units contain \code{0}.
#' @details
#' The sample size under Bernoulli sampling is random, following a
#' Binomial(\code{N}, \code{prob}) distribution. To extract the selected
#' indices, use \code{sam[sam != 0]}.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param prob Scalar. Inclusion probability, must satisfy \code{0 < prob <= 1}.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{E.BE}}, \code{\link{S.PO}}, \code{\link{S.SI}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' # Vector U contains the label of a population of size N=5
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' # Draws a Bernoulli sample without replacement of expected size n=3
#' # The inlusion probability is 0.6 for each unit in the population
#' sam <- S.BE(5,0.6)
#' sam
#' # The selected sample is
#' U[sam]
#'
#' ############
#' ## Example 2
#' ############
#' # Uses the Lucy data to draw a Bernoulli sample
#'
#' data(Lucy)
#' attach(Lucy)
#' N <- dim(Lucy)[1]
#' # The population size is 2396. If the expected sample size is 400
#' # then, the inclusion probability must be 400/2396=0.1669
#' sam <- S.BE(N,0.01669)
#' # The information about the units in the sample is stored in an object called data
#' data <- Lucy[sam,]
#' data
#' dim(data)

S.BE <- function(N, prob) {
  sam <- matrix(0, N, 1)
  U   <- runif(N)
  for (k in 1:N) {
    if (U[k] <= prob)
      sam[k] <- k
  }
  return(sam)
}