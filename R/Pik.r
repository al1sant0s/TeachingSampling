#' @export
#'
#' @title
#' First-Order Inclusion Probabilities from a Sampling Design
#' @description
#' Computes the first-order inclusion probabilities for each unit in a finite
#' population, given the probability of each possible sample and the indicator
#' matrix of the sampling support.
#' @return
#' A row vector (1 x N matrix) of first-order inclusion probabilities
#' \eqn{\pi_k = P(k \in s)} for each unit \eqn{k} in the population.
#' @details
#' The inclusion probability of unit \eqn{k} is computed as the sum of the
#' probabilities of all samples that contain unit \eqn{k}:
#' \deqn{\pi_k = \sum_{s \ni k} p(s)}
#' The indicator matrix \code{Ind} (output of \code{\link{Ik}}) has one row
#' per possible sample and one column per population unit, with entry 1 if
#' unit \eqn{k} is in sample \eqn{s} and 0 otherwise.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param p Vector of probabilities for each possible sample in the support.
#'   Must sum to 1.
#' @param Ind Indicator matrix of the sampling support, as returned by
#'   \code{\link{Ik}}. Rows are samples, columns are population units.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{Ik}}, \code{\link{Pikl}}, \code{\link{PikPPS}}
#'
#' @examples
#' # Population of size N = 5, sample size n = 2
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' N <- length(U)
#' n <- 2
#' # Sample probabilities (one per possible sample)
#' p <- c(0.13, 0.2, 0.15, 0.1, 0.15, 0.04, 0.02, 0.06, 0.07, 0.08)
#' Ind <- Ik(N, n)
#' pik <- Pik(p, Ind)
#' pik
#' # Check: inclusion probabilities sum to n
#' sum(pik)

Pik <- function(p, Ind) {
  multip <- p * Ind
  pik    <- colSums(multip)
  t(pik)
}