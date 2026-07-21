#' @export
#'
#' @title
#' Estimation of the Population Total under Probability-Proportional-to-Size Sampling With Replacement
#' @description
#' Computes the Hansen-Hurwitz estimator of the population total under a
#' probability proportional to size with-replacement (PPS-WR) sampling design.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error of the total.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling.
#' }
#' @details
#' The Hansen-Hurwitz estimator is \eqn{\hat{t} = (1/m)\sum_{i=1}^m y_i/p_i},
#' where \eqn{p_i} is the selection probability of the \eqn{i}-th draw and
#' \eqn{m} is the number of draws.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every selected unit (with possible repetitions).
#' @param pk Vector of selection probabilities for each draw in the sample.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.PPS}}, \code{\link{HH}}, \code{\link{E.piPS}}
#'
#' @examples
#' data('Lucy')
#' attach(Lucy)
#' m   <- 400
#' res <- S.PPS(m, Employees)
#' sam <- res[, 1]
#' pk  <- res[, 2]
#' y   <- data.frame(Income = Income[sam], Taxes = Taxes[sam])
#' E.PPS(y, pk)

E.PPS <- function(y, pk) {
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  m <- length(pk)
  for (k in 1:dim(y)[2]) {
    ty   <- sum(y[, k]/pk)/m
    Vty  <- (1/m) * (1/(m - 1)) * sum((y[, k]/pk - ty)^2)
    CVe  <- 100 * sqrt(Vty)/ty
    N    <- (1/m) * sum(1/pk)
    VMAS <- (N^2) * (1 - (m/N)) * var(y[, k])/(m)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}