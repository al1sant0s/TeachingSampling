#' @export
#'
#' @title
#' Estimation of the Population Total under Simple Random Sampling With
#' Replacement
#' @description
#' Computes the Hansen-Hurwitz estimator of the population total under a
#' simple random sampling with replacement (WR) design.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error of the total.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling
#'     without replacement.
#' }
#' @details
#' Under simple random sampling with replacement with \code{m} draws, the
#' Hansen-Hurwitz estimator is \eqn{\hat{t} = (N/m)\sum_{i=1}^m y_i}.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param m Number of draws (sample size with replacement).
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every draw in the sample (repetitions allowed).
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.WR}}, \code{\link{E.SI}}
#'
#' @examples
#' data('Lucy')
#' attach(Lucy)
#' N   <- nrow(Lucy)
#' m   <- 400
#' sam <- S.WR(N, m)
#' y   <- data.frame(Income = Income[sam], Taxes = Taxes[sam])
#' E.WR(N, m, y)

E.WR <- function(N, m, y) {
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  for (k in 1:dim(y)[2]) {
    ty   <- (N/m) * sum(y[, k])
    Vty  <- (N^2/m) * var(y[, k])
    CVe  <- 100 * sqrt(Vty)/ty
    VMAS <- (N^2) * (1 - (m/N)) * var(y[, k])/(m)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}