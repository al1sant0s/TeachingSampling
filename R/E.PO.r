#' @export
#'
#' @title
#' Estimation of the Population Total under Poisson Sampling Without Replacement
#' @description
#' Computes the Horvitz-Thompson estimator of the population total under a
#' Poisson sampling design, where each unit is independently selected with
#' its own inclusion probability.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error of the total.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling.
#' }
#' @details
#' Under Poisson sampling, units are selected independently, so the exact
#' variance of the Horvitz-Thompson estimator has a simple closed form:
#' \eqn{V(\hat{t}) = \sum_k (1 - \pi_k)(y_k/\pi_k)^2}.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every unit in the selected sample.
#' @param Pik Vector of first-order inclusion probabilities for each unit
#'   in the sample.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.PO}}, \code{\link{E.piPS}}
#'
#' @examples
#' data('Lucy')
#' attach(Lucy)
#' N <- nrow(Lucy)
#' n <- 400
#' Pik <- PikPPS(n, Employees)
#' sam <- S.PO(N, Pik)
#' sam <- sam[sam != 0]
#' y   <- data.frame(Income = Income[sam], Taxes = Taxes[sam])
#' E.PO(y, Pik[sam])

E.PO <- function(y, Pik) {
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  n <- length(Pik)
  for (k in 1:dim(y)[2]) {
    ty   <- sum(y[, k]/Pik)
    Vty  <- sum((1 - Pik) * ((y[, k]/Pik)^2))
    CVe  <- 100 * sqrt(Vty)/ty
    N    <- sum(1/Pik)
    VMAS <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}