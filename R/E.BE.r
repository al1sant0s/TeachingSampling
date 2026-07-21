#' @export
#'
#' @title
#' Estimation of the Population Total under Bernoulli Sampling Without Replacement
#' @description
#' Computes the Horvitz-Thompson estimator of the population total under a
#' Bernoulli sampling design, where each unit in the population is independently
#' selected with the same inclusion probability.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error of the total.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling.
#' }
#' @details
#' Under Bernoulli sampling, the sample size is random. The inclusion
#' probability is constant and equal to \code{prob} for all units. The
#' variance estimator accounts for the randomness of the sample size.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every unit in the selected sample.
#' @param prob Scalar. The (constant) inclusion probability used in the
#'   Bernoulli sampling design. Must satisfy \code{0 < prob <= 1}.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.BE}}, \code{\link{E.SI}}
#'
#' @examples
#' data('Lucy')
#' attach(Lucy)
#' N <- nrow(Lucy)
#' prob <- 0.1
#' sam  <- S.BE(N, prob)
#' sam  <- sam[sam != 0]
#' y    <- data.frame(Income = Income[sam], Employees = Employees[sam])
#' E.BE(y, prob)

E.BE <- function(y, prob) {
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  
  for (k in 1:dim(y)[2]) {
    ty   <- sum(y[, k])/prob
    Vty  <- (1/prob) * ((1/prob) - 1) * sum(y[, k]^2)
    CVe  <- 100 * sqrt(Vty)/ty
    n    <- length(y[, k])
    N    <- n/prob
    VMAS <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}