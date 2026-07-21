#' @export
#'
#' @title
#' Estimation of the Population Total under Stratified Probability Proportional to Size Sampling With Replacement
#' @description
#' Computes the Hansen-Hurwitz estimator of the population total under a
#' stratified PPS with-replacement (STPPS) sampling design.
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
#' @param pk Vector of selection probabilities for each draw in the sample.
#' @param mh Integer vector with the number of draws within each stratum.
#' @param S Vector identifying the stratum membership of each unit in the sample.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.STPPS}}, \code{\link{E.PPS}}, \code{\link{E.STpiPS}}
#'
#' @examples
#' # Uses the Lucy data to draw a stratified random sample
#' # according to a PPS design in each stratum
#' data(Lucy)
#' attach(Lucy)
#' m1 <- 83; m2 <- 100; m3 <- 200
#' mh <- c(m1, m2, m3)
#' res <- S.STPPS(Level, Income, mh)
#' sam <- res[, 1]
#' pk  <- res[, 2]
#' data <- Lucy[sam, ]
#' attach(data)
#' estima <- data.frame(Income, Employees, Taxes)
#' E.STPPS(estima, pk, mh, Level)

E.STPPS <- function(y, pk, mh, S) {
  S <- as.factor(S)
  S <- as.factor(as.integer(S))
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  H <- length(mh)
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  N <- sum(1/pk)
  n <- length(pk)
  for (k in 1:dim(y)[2]) {
    ty  <- 0
    Vty <- 0
    for (h in 1:H) {
      yh  <- y[which(S == h), k]
      pkh <- pk[which(S == h)]
      HHh <- sum(yh/pkh)/mh[h]
      ty  <- ty + HHh
      Vty <- Vty + (1/mh[h]) * (1/(mh[h] - 1)) * sum((yh/pkh - HHh)^2)
    }
    CVe  <- 100 * sqrt(Vty)/ty
    VMAS <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}