#' @export
#'
#' @title
#' Estimation of Regression Coefficients under Simple Random Sampling Without Replacement
#' @description
#' Computes the weighted least squares estimator of regression coefficients
#' for a finite population under simple random sampling without replacement.
#' Both the estimated coefficients and their estimated standard errors are
#' returned.
#' @return
#' A three-dimensional array with dimensions \code{[3, P, Q]}, where
#' \code{P} is the number of auxiliary variables and \code{Q} is the number
#' of variables of interest. The three rows correspond to:
#' \itemize{
#'   \item \code{Beta estimation}: Estimated regression coefficient.
#'   \item \code{Standard Error}: Estimated standard error.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#' }
#' @details
#' The estimator uses a working model with weights \eqn{V = 1/(\pi_k c_k)},
#' where \eqn{\pi_k = n/N} under simple random sampling and \eqn{c_k} is an
#' optional variance-stabilising constant. The variance is estimated using
#' the residual-based sandwich approach of Sarndal et al. (1992).
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param n Sample size.
#' @param y Vector, matrix or data frame of variables of interest (response).
#' @param x Vector, matrix or data frame of auxiliary variables (predictors).
#' @param ck Optional variance-stabilising constant. Default is \code{1}
#'   (homoscedastic model).
#' @param b0 Logical. If \code{TRUE}, an intercept column of ones is
#'   prepended to \code{x}. Default is \code{FALSE}.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{GREG.SI}}, \code{\link{E.SI}}
#'
#' @examples
#' ######################################################################
#' ## Example 1: Linear models involving continuous auxiliary information
#' ######################################################################
#'
#' # Draws a simple random sample without replacement
#' data(Lucy)
#' attach(Lucy)
#'
#' N <- dim(Lucy)[1]
#' n <- 400
#' sam <- S.SI(N, n)
#' # The information about the units in the sample 
#' # is stored in an object called data
#' data <- Lucy[sam,]
#' attach(data)
#' names(data)
#'
#' ########### common mean model 
#'
#' estima<-data.frame(Income, Employees, Taxes)
#' x <- rep(1,n)
#' E.Beta(N, n, estima,x,ck=1,b0=FALSE)
#'
#'
#' ########### common ratio model 
#'
#' estima<-data.frame(Income)
#' x <- data.frame(Employees)
#' E.Beta(N, n, estima,x,ck=x,b0=FALSE)
#'
#' ########### Simple regression model without intercept
#'
#' estima<-data.frame(Income, Employees)
#' x <- data.frame(Taxes)
#' E.Beta(N, n, estima,x,ck=1,b0=FALSE)
#'
#' ########### Multiple regression model without intercept
#'
#' estima<-data.frame(Income)
#' x <- data.frame(Employees, Taxes)
#' E.Beta(N, n, estima,x,ck=1,b0=FALSE)
#'
#' ########### Simple regression model with intercept
#'
#' estima<-data.frame(Income, Employees)
#' x <- data.frame(Taxes)
#' E.Beta(N, n, estima,x,ck=1,b0=TRUE)
#'
#' ########### Multiple regression model with intercept
#'
#' estima<-data.frame(Income)
#' x <- data.frame(Employees, Taxes)
#' E.Beta(N, n, estima,x,ck=1,b0=TRUE)
#'
#' ###############################################################
#' ## Example 2: Linear models with discrete auxiliary information
#' ###############################################################
#'
#' # Draws a simple random sample without replacement
#' data(Lucy)
#' attach(Lucy)
#'
#' N <- dim(Lucy)[1]
#' n <- 400
#' sam <- S.SI(N,n)
#' # The information about the sample units is stored in an object called data
#' data <- Lucy[sam,]
#' attach(data)
#' names(data)
#' # The auxiliary information
#' Doma<-Domains(Level)
#'
#' ########### Poststratified common mean model
#'
#' estima<-data.frame(Income, Employees, Taxes)
#' E.Beta(N, n, estima,Doma,ck=1,b0=FALSE)
#'
#' ########### Poststratified common ratio model
#'
#' estima<-data.frame(Income, Employees)
#' x<-Doma*Taxes
#' E.Beta(N, n, estima,x,ck=1,b0=FALSE)

E.Beta <- function(N, n, y, x, ck = 1, b0 = FALSE) {
  if (b0 == TRUE) {
    x <- as.data.frame(cbind(1, x))
  }
  Q <- dim(as.matrix(y))[2]
  P <- dim(as.matrix(x))[2]
  Total <- array(NA, c(3, P, Q))
  rownames(Total) <- c("Beta estimation", "Standard Error", "CVE")
  colnames(Total) <- names(x)
  dimnames(Total)[[3]] <- names(y)
  Pik <- rep(n/N, n)
  for (q in 1:Q) {
    yq <- as.matrix(y[, q])
    x  <- as.matrix(x)
    ck <- as.numeric(unlist(ck))
    V  <- 1/(Pik * ck)
    bq <- solve(t(V * x) %*% x) %*% (t(V * x) %*% yq)
    ek <- yq - x %*% bq
    uk <- c(ek) * x
    Varuk <- (N^2/n) * (1 - (n/N)) * var(uk)
    P1    <- solve(t(V * x) %*% x)
    Vbeta <- as.matrix(P1) %*% as.matrix(Varuk) %*% as.matrix(P1)
    Vbeta <- diag(Vbeta)
    CVe   <- 100 * sqrt(Vbeta)/bq
    if (Q == 1) {
      Total[1, , ] <- bq
      Total[2, , ] <- sqrt(Vbeta)
      Total[3, , ] <- CVe
    }
    if (P == 1 & Q > 1) {
      Total[1, , ][q] <- bq
      Total[2, , ][q] <- sqrt(Vbeta)
      Total[3, , ][q] <- CVe
    }
    if (Q > 1 & P > 1) {
      Total[1, , ][, q] <- bq
      Total[2, , ][, q] <- sqrt(Vbeta)
      Total[3, , ][, q] <- CVe
    }
  }
  return(Total)
}