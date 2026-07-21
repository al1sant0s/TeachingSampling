#' @export
#'
#' @title
#' Estimation of the Population Total under Stratified Simple Random Sampling Without Replacement
#' @description
#' Computes the Horvitz-Thompson estimator of the population total under a
#' stratified simple random sampling without replacement (STSI) design.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling.
#' }
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param S Vector identifying the stratum membership of each unit in the sample.
#' @param Nh Integer vector with the population size of each stratum.
#' @param nh Integer vector with the sample size of each stratum.
#' @param y Vector, matrix or data frame of variables of interest.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.STSI}}, \code{\link{E.SI}}, \code{\link{E.STpiPS}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' data(Lucy)
#' attach(Lucy)
#' N1 <- summary(Level)[[1]]
#' N2 <- summary(Level)[[2]]
#' N3 <- summary(Level)[[3]]
#' Nh <- c(N1, N2, N3)
#' n1 <- N1; n2 <- 100; n3 <- 200
#' nh <- c(n1, n2, n3)
#' sam <- S.STSI(Level, Nh, nh)
#' data <- Lucy[sam, ]
#' attach(data)
#' estima <- data.frame(Income, Employees, Taxes)
#' E.STSI(Level, Nh, nh, estima)
#' ############
#' ## Example 2
#' ############
#' # The variable SPAM is a domain of interest
#' Doma <- Domains(SPAM)
#' SPAM.no  <- estima * Doma[, 1]
#' SPAM.yes <- estima * Doma[, 2]
#' E.STSI(Level, Nh, nh, Doma)
#' E.STSI(Level, Nh, nh, SPAM.no)
#' E.STSI(Level, Nh, nh, SPAM.yes)

E.STSI <- function(S, Nh, nh, y) {
  S <- as.factor(S)
  S <- as.factor(as.integer(S))
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  H <- length(Nh)
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  N <- sum(Nh)
  n <- sum(nh)
  fh <- nh/Nh
  wh <- Nh/N
  for (k in 1:dim(y)[2]) {
    ty   <- 0
    Vty  <- 0
    for (h in 1:H) {
      yh   <- y[which(S == h), k]
      ty   <- ty + Nh[h] * mean(yh)
      Vty  <- Vty + Nh[h]^2 * (1 - fh[h]) * var(yh)/nh[h]
    }
    CVe  <- 100 * sqrt(Vty)/ty
    VMAS <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}