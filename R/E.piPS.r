#' @export
#'
#' @title
#' Estimation of the Population Total under Probability-Proportional-to-Size Sampling Without Replacement
#' @description
#' Computes the Horvitz-Thompson estimator of the population total under a
#' without-replacement probability proportional to size (piPS) sampling design.
#' The variance is estimated using the Horvitz-Thompson variance approximation
#' based on first-order inclusion probabilities.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error of the total.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling.
#' }
#' @details
#' When all inclusion probabilities are equal (i.e. \code{sum(Pik) == n}),
#' the variance is set to zero, reflecting an equal-probability design.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every unit in the selected sample.
#' @param Pik Vector of first-order inclusion probabilities for each
#'   unit in the sample.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.piPS}}, \code{\link{PikPPS}}, \code{\link{E.PO}}
#'
#' @examples
#' data('Lucy')
#' attach(Lucy)
#' N <- nrow(Lucy)
#' n <- 400
#' x <- Employees
#' res <- S.piPS(n, x)
#' sam <- res[, 1]
#' Pik <- res[, 2]
#' y   <- data.frame(Income = Income[sam], Taxes = Taxes[sam])
#' E.piPS(y, Pik)

E.piPS <- function(y, Pik) {
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  n <- length(Pik)
  for (k in 1:dim(y)[2]) {
    ty <- sum(y[, k]/Pik)
    ck    <- (1 - Pik) * (n/(n - 1))
    P1    <- sum(ck * y[, k]/Pik)
    P2    <- sum(ck)
    ystar <- Pik * P1/P2
    P3    <- ck/(Pik^2)
    if (sum(Pik) == n) {
      Vty <- 0
    } else {
      Vty <- sum(P3 * ((y[, k] - ystar)^2))
    }
    CVe  <- 100 * sqrt(Vty)/ty
    N    <- sum(1/Pik)
    VMAS <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}