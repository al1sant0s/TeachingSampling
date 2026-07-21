#' @export
#'
#' @title
#' Exact Variance of the Horvitz-Thompson Estimator
#' @description
#' Computes the exact variance of the Horvitz-Thompson estimator of the
#' population total for a given fixed-size without-replacement sampling design,
#' using the full sampling support.
#' @return
#' A scalar: the exact variance of the Horvitz-Thompson estimator
#' \eqn{V(\hat{t}_{y,\pi})}.
#' @details
#' The exact Horvitz-Thompson variance is:
#' \deqn{V(\hat{t}_{y,\pi}) = \sum_{k=1}^N \sum_{l=1}^N \Delta_{kl}
#' \frac{y_k}{\pi_k} \frac{y_l}{\pi_l}}
#' where \eqn{\Delta_{kl} = \pi_{kl} - \pi_k \pi_l}. This requires
#' enumerating the full support and is only feasible for small populations
#' (\code{N <= 15}).
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param y Vector of length \code{N} with the population values of the
#'   variable of interest.
#' @param N Population size. Recommended \code{N <= 15}.
#' @param n Sample size.
#' @param p Vector of probabilities for each possible sample in the support.
#'   Must sum to 1.
#'
#' @references
#' Horvitz, D.G. and Thompson, D.J. (1952). A generalization of sampling
#' without replacement from a finite universe.
#' \emph{Journal of the American Statistical Association}, 47, 663-685.\cr
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.
#'
#' @seealso \code{\link{Deltakl}}, \code{\link{VarSYGHT}}, \code{\link{HT}}
#'
#' @examples
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' y1 <- c(32, 34, 46, 89, 35)
#' y2 <- c(1, 1, 1, 0, 0)
#' N <- length(U)
#' n <- 2
#' p <- c(0.13, 0.2, 0.15, 0.1, 0.15, 0.04, 0.02, 0.06, 0.07, 0.08)
#' # Theoretical variance of the HT estimator
#' VarHT(y1, N, n, p)
#' VarHT(y2, N, n, p)

VarHT <- function(y, N, n, p) {
  Ind      <- Ik(N, n)
  pi1      <- as.matrix(Pik(p, Ind))
  pi2      <- Pikl(N, n, p)
  Delta    <- Deltakl(N, n, p)
  y        <- t(as.matrix(y))
  ykylexp  <- t(y/pi1) %*% (y/pi1)
  A        <- (Delta) * (ykylexp)
  Var      <- sum(A)
  return(Var)
}