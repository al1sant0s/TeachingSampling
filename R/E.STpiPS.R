#' @export
#'
#' @title
#' Estimation of the Population Total under Stratified Probability Proportional to Size Sampling Without Replacement
#' @description
#' Computes the Horvitz-Thompson estimator of the population total under a
#' stratified without-replacement probability proportional to size (piPS)
#' sampling design.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling.
#' }
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param y Vector, matrix or data frame of variables of interest.
#' @param Pik Vector of first-order inclusion probabilities for each unit
#'   in the sample.
#' @param S Vector identifying the stratum membership of each unit in the sample.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.STpiPS}}, \code{\link{E.piPS}}, \code{\link{E.STSI}}
#'
#' @examples
#' # Uses the Lucy data to draw a stratified random sample
#' # according to a piPS design in each stratum
#' data(Lucy)
#' attach(Lucy)
#' N1 <- summary(Level)[[1]]
#' N2 <- summary(Level)[[2]]
#' N3 <- summary(Level)[[3]]
#' nh <- c(N1, 100, 200)
#' S  <- Level
#' x  <- Employees
#' res <- S.STpiPS(S, x, nh)
#' sam <- res[, 1]
#' pik <- res[, 2]
#' data <- Lucy[sam, ]
#' attach(data)
#' estima <- data.frame(Income, Employees, Taxes)
#' E.STpiPS(estima, pik, Level)

E.STpiPS <- function(y, Pik, S) {
  S <- as.factor(S)
  S <- as.factor(as.integer(S))
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  H <- length(levels(S))
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  n <- length(Pik)
  N <- sum(1/Pik)
  for (k in 1:dim(y)[2]) {
    ty  <- 0
    Vty <- 0
    for (h in 1:H) {
      yh   <- y[which(S == h), k]
      pikh <- Pik[which(S == h)]
      nh   <- length(pikh)
      ck   <- (1 - pikh) * (nh/(nh - 1))
      P1   <- sum(ck * yh/pikh)
      P2   <- sum(ck)
      ystar <- pikh * P1/P2
      P3   <- ck/(pikh^2)
      Vty  <- Vty + sum(P3 * ((yh - ystar)^2))
      ty   <- ty + sum(yh/pikh)
    }
    CVe  <- 100 * sqrt(Vty)/ty
    VMAS <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}