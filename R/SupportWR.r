#' @export
#'
#' @title
#' Sampling Support for With-Replacement Designs
#' @description
#' Enumerates all distinct unordered outcomes (multisets) of size \code{m}
#' drawn with replacement from a population of size \code{N}.
#' @return
#' A matrix with \code{choose(N+m-1, m)} rows and \code{m} columns. Each
#' row contains the (sorted) indices of one possible unordered outcome.
#' If \code{ID} is provided, population labels replace indices.
#' @details
#' The number of distinct unordered with-replacement outcomes of size \code{m}
#' from \code{N} units is \eqn{\binom{N+m-1}{m}}. This is much smaller than
#' the \eqn{N^m} ordered outcomes. The algorithm uses a nested loop to
#' generate all non-decreasing sequences of length \code{m} from
#' \eqn{\{1, \ldots, N\}}.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param m Number of draws (sample size with replacement).
#' @param ID Optional vector of population labels of length \code{N}.
#'   If \code{FALSE} (default), integer indices are returned.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{IkWR}}, \code{\link{nk}}, \code{\link{p.WR}}
#'
#' @examples
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' N <- length(U)
#' m <- 2
#' # With-replacement support
#' SupportWR(N, m)
#' SupportWR(N, m, ID = U)
#' y <- c(32, 34, 46, 89, 35)
#' SupportWR(N, m, ID = y)

SupportWR <- function(N, m, ID = FALSE) {
  S  <- 0
  a  <- rep(1, m)
  P1 <- a
  S  <- S + 1
  k  <- m
  while (k > 0) {
    while (a[k] < N) {
      a[k] <- a[k] + 1
      P1   <- rbind(P1, a)
      S    <- S + 1
    }
    if (k > 1) k <- k - 1
    if (a[k] < N) {
      a[k]     <- a[k] + 1
      k1       <- k + 1
      a[k1:m]  <- a[k]
      P1       <- rbind(P1, a)
      S        <- S + 1
      k        <- m
    } else {
      if (k == 1) k <- 0
    }
  }
  nr  <- choose(N + m - 1, m)
  P1  <- matrix(P1, nrow = nr)
  sam <- matrix(ID[P1], nrow = nr)
  if (is.logical(ID) == TRUE) return(P1)
  else return(sam)
}