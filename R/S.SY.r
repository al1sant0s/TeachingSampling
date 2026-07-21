#' @export
#'
#' @title
#' Systematic Sampling
#' @description
#' Draws a systematic sample from a finite population of size \code{N} using
#' a fixed sampling interval \code{a}. A random start \code{r} is chosen
#' uniformly from \code{1} to \code{a}, and every \code{a}-th unit thereafter
#' is selected.
#' @return
#' A vector containing the population indices of the selected units.
#' @details
#' The random start \code{r} is drawn from \code{sample(a, 1)}, and then
#' units \eqn{r, r+a, r+2a, \ldots} are selected. If \code{N} is not a
#' multiple of \code{a}, the sample size varies by one unit depending on the
#' random start. Use \code{\link{E.SY}} to estimate population totals.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param a Sampling interval (skip). The expected sample size is
#'   approximately \code{N/a}.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{E.SY}}, \code{\link{S.SI}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' # Vector U contains the label of a population of size N=5
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' # The population of size N=5 is divided in a=2 groups
#' # Draws a Systematic sample. 
#' sam <- S.SY(5,2)
#' sam
#' # The selected sample is
#' U[sam]
#' # There are only two possible samples
#'
#' ############
#' ## Example 2
#' ############
#' # Uses the Lucy data to draw a Systematic sample
#' data(Lucy)
#' attach(Lucy)
#'
#' N <- dim(Lucy)[1]
#' # The population is divided in 6 groups
#' # The selected sample
#' sam <- S.SY(N,6)
#' # The information about the units in the sample is stored in an object called data
#' data <- Lucy[sam,]
#' data
#' dim(data)

S.SY <- function(N, a) {
  r <- sample(a, 1)
  c <- N - a * floor(N/a)
  if (r <= c)
    n <- floor((N/a)) + 1
  else
    n <- floor(N/a)
  sam <- matrix(0, n, 1)
  for (k in 1:n) {
    sam[k] <- r + (a * (k - 1))
  }
  sam
}