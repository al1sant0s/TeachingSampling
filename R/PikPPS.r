#' @export
#'
#' @title
#' Inclusion Probabilities Proportional to Size
#' @description
#' Computes first-order inclusion probabilities proportional to an auxiliary
#' size variable \code{x} for a without-replacement sample of size \code{n}.
#' A sequential truncation algorithm ensures all probabilities are at most 1.
#' @return
#' A numeric vector of length \code{N} with the first-order inclusion
#' probability for each unit in the population. Values are in \code{(0, 1]}.
#' @details
#' The initial probabilities \eqn{\pi_k = n x_k / \sum x} may exceed 1 for
#' large units. The algorithm iteratively sets those probabilities to 1 and
#' redistributes the remaining sample size among the other units until all
#' probabilities are valid. The result satisfies \eqn{\sum \pi_k = n}.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param n Desired sample size.
#' @param x Vector of length \code{N} with positive auxiliary size values
#'   for each unit in the population.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.piPS}}, \code{\link{PikSTPPS}}, \code{\link{PikHol}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' x <- c(30,41,50,170,43,200)
#' n <- 3
#' # Two elements yields values bigger than one
#' n*x/sum(x)
#' # With this functions, all of the values are between zero and one
#' PikPPS(n,x)
#' # The sum is equal to the sample size
#' sum(PikPPS(n,x))
#'
#' ############
#' ## Example 2
#' ############
#' # Vector U contains the label of a population of size N=5
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' # The auxiliary information
#' x <- c(52, 60, 75, 100, 50)
#' # Gives the inclusion probabilities for the population accordin to a 
#' # proportional to size design without replacement of size n=4
#' pik <- PikPPS(4,x)
#' pik
#' # The selected sample is
#' sum(pik)
#'
#' ############
#' ## Example 3
#' ############
#' # Uses the Lucy data to compute teh vector of inclusion probabilities 
#' # accordind to a piPS without replacement design
#' data(Lucy)
#' attach(Lucy)
#' # The sample size
#' n=400
#' # The selection probability of each unit is proportional to the variable Income
#' pik <- PikPPS(n,Income)
#' # The inclusion probabilities of the units in the sample
#' pik
#' # The sum of the values in pik is equal to the sample size
#' sum(pik)
#' # According to the design some elements must be selected
#' # They are called forced inclusion units
#' which(pik==1)

PikPPS <- function(n, x) {
  pik <- n * x/sum(x)
  while ((sum(pik > 1)) != 0) {
    s      <- which(pik >= 1)
    new    <- (1:length(pik))[-s]
    pik[s] <- 1
    txnew  <- sum(x[s])
    for (k in new) {
      pik[k] <- (n - length(s)) * x[k]/(sum(x) - txnew)
    }
  }
  pik
}