#' @export
#'
#' @title
#' Variance-Covariance Matrix of the Sample Membership Indicators for Fixed Size Without Replacement Sampling Designs
#' @description
#' Computes the matrix \eqn{\Delta_{kl} = \pi_{kl} - \pi_k \pi_l} for all
#' pairs of units in a finite population. This matrix appears in the exact
#' Horvitz-Thompson variance formula.
#' @return
#' An \code{N x N} matrix where entry \eqn{(k, l)} equals
#' \eqn{\pi_{kl} - \pi_k \pi_l}. Diagonal entries equal
#' \eqn{\pi_k(1 - \pi_k)}.
#' @details
#' The matrix \eqn{\Delta} is central to the Horvitz-Thompson variance
#' estimator:
#' \deqn{V(\hat{t}_{y,\pi}) = \sum_k \sum_l \Delta_{kl} \frac{y_k}{\pi_k}
#' \frac{y_l}{\pi_l}}
#' It requires computing both first-order (\code{\link{Pik}}) and
#' second-order (\code{\link{Pikl}}) inclusion probabilities, so it is only
#' feasible for small populations.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size. Recommended \code{N <= 15}.
#' @param n Sample size.
#' @param p Vector of probabilities for each possible sample in the support.
#'   Must sum to 1.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{Pik}}, \code{\link{Pikl}}, \code{\link{VarHT}}
#'
#' @examples
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' N <- length(U)
#' n <- 2
#' p <- c(0.13, 0.2, 0.15, 0.1, 0.15, 0.04, 0.02, 0.06, 0.07, 0.08)
#' sum(p)
#' # Variance-Covariance matrix of the sample membership indicators
#' Deltakl(N, n, p)

Deltakl <- function(N, n, p) {
  Ind   <- Ik(N, n)
  P1    <- as.matrix(Pik(p, Ind))
  Delta <- Pikl(N, n, p) - (t(P1) %*% P1)
  return(Delta)
}