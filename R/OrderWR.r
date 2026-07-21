#' @export
#'
#' @title
#' Ordered With-Replacement Sampling Support
#' @description
#' Enumerates all ordered sequences of \code{m} draws from a population of
#' size \code{N} with replacement. Unlike \code{\link{SupportWR}}, this
#' function considers order, so sequences that differ only in draw order are
#' treated as distinct outcomes.
#' @return
#' A matrix with \code{N^m} rows and \code{m} columns, where each row is one
#' ordered sequence of draws. If \code{ID} is provided, population labels are
#' substituted for indices.
#' @details
#' The total number of ordered with-replacement sequences of size \code{m}
#' from \code{N} units is \eqn{N^m}. This grows rapidly and the function
#' should only be used for small \code{N} and \code{m}.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param m Number of draws.
#' @param ID Optional vector of population labels of length \code{N}.
#'   If provided, labels are substituted for integer indices in the output.
#'   If \code{FALSE} (default), integer indices are returned.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{SupportWR}}, \code{\link{IkWR}}
#'
#' @examples
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' N <- length(U)
#' # Five possible ordered samples of size m=1
#' OrderWR(N, 1)
#' OrderWR(N, 1, ID = U)
#' # 25 possible ordered samples of size m=2
#' OrderWR(N, 2)
#' OrderWR(N, 2, ID = U)
#' # Note: ordered samples differ from unordered (SupportWR)
#' OrderWR(N, 2)
#' SupportWR(N, 2)

OrderWR <- function(N, m, ID = FALSE) {
  b <- c(1:N)
  grilla <- function(a) {
    A    <- seq(1:length(a))
    unoA <- rep(1, length(A))
    B    <- seq(1:length(a))
    unoB <- rep(1, length(B))
    P1   <- kronecker(A, unoB)
    P2   <- kronecker(unoA, B)
    grid <- matrix(cbind(P1, P2), ncol = 2)
    return(grid)
  }
  if (m == 1) sam <- as.matrix(b)
  if (m == 2) sam <- grilla(b)
  if (m > 2) {
    sam <- grilla(b)
    for (l in 3:m) {
      Sam1 <- rep(0, l)
      for (j in 1:dim(sam)[1]) {
        for (k in 1:length(b)) {
          Sam1 <- rbind(Sam1, c(sam[j, ], b[k]))
        }
      }
      sam <- Sam1[-1, ]
    }
  }
  if (is.logical(ID) == TRUE) return(sam)
  else {
    a   <- dim(sam)
    val <- matrix(NA, a[1], a[2])
    for (ii in 1:(dim(val)[1])) {
      for (jj in 1:(dim(val)[2])) {
        val[ii, jj] <- ID[sam[ii, jj]]
      }
    }
    return(val)
  }
}