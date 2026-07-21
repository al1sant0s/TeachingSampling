#' @export
#' @import stats
#'
#' @title
#' Sampling Support for Fixed-Size Without-Replacement Designs
#' @description
#' Enumerates all possible samples of size \code{n} from a population of
#' size \code{N}, returning the complete sampling support as a matrix.
#' @return
#' A matrix with \code{choose(N, n)} rows and \code{n} columns. Each row
#' contains the indices (or labels if \code{ID} is provided) of the units
#' in one possible sample. Samples are listed in lexicographic order.
#' @details
#' This function uses a combinatorial algorithm to enumerate all
#' \code{choose(N, n)} subsets of size \code{n} from \eqn{\{1, \ldots, N\}}.
#' It is intended for small populations only. For \code{N > 15} it becomes
#' very slow.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size. Recommended \code{N <= 15}.
#' @param n Sample size.
#' @param ID Optional vector of population labels of length \code{N}.
#'   If provided, labels replace integer indices in the output.
#'   If \code{FALSE} (default), integer indices are returned.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{Ik}}, \code{\link{SupportWR}}, \code{\link{SupportRS}}
#'
#' @examples
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' N <- length(U)
#' n <- 2
#' # Ten possible samples of size n=2
#' Support(N, n)
#' # Labeled support
#' Support(N, n, ID = U)
#' # Support showing values of y
#' y <- c(32, 34, 46, 89, 35)
#' Support(N, n, ID = y)

Support <- function(N, n, ID = FALSE) {
  m   <- matrix(0, choose(N, n), n)
  sam <- matrix(0, choose(N, n), n)
  for (i in 1:n) {
    a <- 0
    t <- i
    for (r in 1:choose(N, n)) {
      a <- a + 1
      B <- choose(N - t, n - i)
      if (a > B) {
        a <- 1
        t <- t + 1
      }
      if (t > N - n + i) {
        t <- m[r, i - 1] + 1
      }
      m[r, i]   <- t
      sam[r, i] <- ID[t]
    }
  }
  if (is.logical(ID) == TRUE) return(m)
  else return(sam)
}