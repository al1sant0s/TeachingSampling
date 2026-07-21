#' @export
#'
#' @title
#' Estimation of a Population Quantile
#' @description
#' Computes a weighted quantile estimator for finite populations. When
#' inclusion probabilities are provided, the estimator uses the
#' Horvitz-Thompson weights \eqn{d_k = 1/\pi_k}; otherwise, equal weights
#' are assumed (simple random sampling).
#' @return
#' A numeric vector of length equal to the number of variables in \code{y},
#' containing the estimated quantile for each variable.
#' @details
#' The estimator is based on the weighted empirical cumulative distribution
#' function. For each variable, units are sorted by their observed value,
#' cumulative weights are computed, and the quantile is located by
#' interpolation.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every unit in the selected sample.
#' @param Qn Scalar in \eqn{(0, 1)}. The desired quantile level
#'   (e.g. \code{0.5} for the median, \code{0.25} for the first quartile).
#' @param Pik Optional vector of first-order inclusion probabilities. If
#'   omitted, equal probabilities are assumed.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{E.SI}}, \code{\link{E.piPS}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' y <- c(32, 34, 46, 89, 35)
#' x <- c(52, 60, 75, 100, 50)
#' z <- cbind(y, x)
#' Pik <- c(0.58, 0.34, 0.48, 0.33, 0.27)
#' E.Quantile(y, 0.5)
#' E.Quantile(x, 0.25)
#' E.Quantile(z, 0.75)
#' E.Quantile(z, 0.5, Pik)
#' ############
#' ## Example 2
#' ############
#' data(Lucy)
#' attach(Lucy)
#' m <- 400
#' res <- S.PPS(m, Income)
#' sam <- res[, 1]
#' pk.s <- res[, 2]
#' Pik.s <- 1 - (1 - pk.s)^m
#' data <- Lucy[sam, ]
#' attach(data)
#' estima <- data.frame(Income, Employees, Taxes)
#' E.Quantile(estima, 0.5, Pik.s)

E.Quantile <- function(y, Qn, Pik) {
  y <- as.data.frame(y)
  Total <- rep(NA, dim(y)[2])
  if (missing(Pik))
    Pik <- rep(1, dim(y)[1])
  if (any(Pik < 0))
    stop("Probabilities must be positive.")
  w <- 1/Pik
  n <- length(w)
  for (i in 1:dim(y)[2]) {
    ord  <- order(y[, i])
    x    <- y[ord, i]
    w    <- (1/Pik)[ord]
    wcum <- cumsum(w)
    wsum <- wcum[n]
    wper <- wsum * Qn
    lows <- (wcum <= wper)
    k    <- sum(lows)
    if (k != 0 && k != n) {
      wlow  <- wcum[k]
      whigh <- wsum - wlow
      if (whigh > wper)
        Total[i] <- x[k + 1]
      else
        Total[i] <- (wlow * x[k] + whigh * x[k + 1])/wsum
    }
    if (k == 0) Total[i] <- x[1]
    if (k == n) Total[i] <- x[n]
  }
  return(Total)
}