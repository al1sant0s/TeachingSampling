#' @export
#'
#' @title
#' Probability Proportional to Size Without-Replacement Sampling (piPS)
#' @description
#' Draws a without-replacement sample of size \code{n} using a sequential
#' algorithm that produces inclusion probabilities proportional to an
#' auxiliary size variable \code{x}.
#' @return
#' A matrix with \code{n} rows and two columns:
#' \itemize{
#'   \item Column 1: population indices of the selected units.
#'   \item Column 2: first-order inclusion probabilities of the selected units.
#' }
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param n Sample size.
#' @param x Vector of length \code{N} with positive auxiliary size values.
#' @param e Optional vector of \code{N} uniform random variates in \code{(0,1)}.
#'   If omitted, \code{runif(N)} is used.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{E.piPS}}, \code{\link{PikPPS}}, \code{\link{S.STPPS}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' x <- c(52, 60, 75, 100, 50)
#' # Draws a piPS sample without replacement of size n=3
#' res <- S.piPS(3, x)
#' res
#' sam <- res[, 1]
#' U[sam]
#' ############
#' ## Example 2
#' ############
#' # Uses the Lucy data
#' data(Lucy)
#' attach(Lucy)
#' res <- S.piPS(400, Income)
#' sam <- res[, 1]
#' Pik.s <- res[, 2]
#' data <- Lucy[sam, ]
#' dim(data)

S.piPS <- function(n, x, e = runif(length(x))) {
  if (length(x) != 1) {
    N   <- length(x)
    x1  <- sort(x, decreasing = TRUE)
    Pik <- PikPPS(n, x1)
    V   <- cumsum(Pik)
    nk  <- matrix(0, N, 1)
    d   <- matrix(0, N, 1)
    I   <- matrix(0, N, 1)
    sam <- matrix(0, N, 1)
    if (e[1] < Pik[1]) {
      I[1]   <- 1
      sam[1] <- 1
    }
    for (k in 2:N) {
      nk[k] <- nk[k - 1] + I[k - 1]
      d[k]  <- Pik[k] * (n - nk[k])/(n - V[k - 1])
      if (e[k] <= d[k]) {
        I[k]   <- 1
        sam[k] <- cumsum(I[1:(k - 1)])[(k - 1)] + I[k]
      }
    }
    samp  <- rev(order(x))[which(sam != 0)]
    Pik1  <- PikPPS(n, x)
    Pik.s <- Pik1[samp]
    return(cbind(samp, Pik.s))
  }
  if (length(x) == 1) {
    Pik.s <- 1
    samp  <- 1
    return(cbind(samp, Pik.s))
  }
}