#' @export
#'
#' @title
#' Estimation of the Population Total under Systematic Sampling Without Replacement
#' @description
#' Computes the Horvitz-Thompson estimator of the population total under a
#' systematic sampling design with sampling interval \code{a}.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error of the total.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling.
#' }
#' @details
#' Under systematic sampling the sample size is \eqn{n = N/a}. Because only
#' one systematic sample is observed, the variance cannot be estimated without
#' assumptions. Here, the variance is approximated by treating the systematic
#' sample as a simple random sample of the same size, which is a common
#' conservative approximation.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param a Sampling interval (skip). The expected sample size is \code{N/a}.
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every unit in the selected sample.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.SY}}, \code{\link{E.SI}}
#'
#' @examples
#' data('Lucy')
#' attach(Lucy)
#' N   <- nrow(Lucy)
#' a   <- 10
#' sam <- S.SY(N, a)
#' y   <- data.frame(Income = Income[sam], Taxes = Taxes[sam])
#' E.SY(N, a, y)

E.SY <- function(N, a, y) {
  n <- N/a
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  for (k in 1:dim(y)[2]) {
    ty   <- a * sum(y[, k])
    Vty  <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    CVe  <- 100 * sqrt(Vty)/ty
    VMAS <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}