#' @export
#'
#' @title
#' Generalised Regression Estimator under Simple Random Sampling
#' @description
#' Computes the Generalised Regression (GREG) estimator of the population
#' total under simple random sampling without replacement, using known
#' population totals of auxiliary variables to improve efficiency.
#' @return
#' A matrix with three rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: GREG estimated population total.
#'   \item \code{Standard Error}: Estimated standard error.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#' }
#' @details
#' The GREG estimator is:
#' \deqn{\hat{t}_{GREG} = \hat{t}_{HT} + (\mathbf{t}_x -
#' \hat{\mathbf{t}}_{x,HT})^T \hat{\boldsymbol{\beta}}}
#' where \eqn{\hat{\boldsymbol{\beta}}} are the regression coefficients
#' estimated from the sample, \eqn{\mathbf{t}_x} are the known population
#' totals, and variance is estimated from the residuals.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param n Sample size.
#' @param y Vector, matrix or data frame of variables of interest.
#' @param x Vector, matrix or data frame of auxiliary variables observed
#'   in the sample.
#' @param tx Vector of known population totals for the auxiliary variables.
#' @param b Matrix of regression coefficients (e.g. from \code{\link{E.Beta}}).
#' @param b0 Logical. If \code{TRUE}, an intercept column is prepended to
#'   \code{x}. Default is \code{FALSE}.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{E.Beta}}, \code{\link{E.SI}}, \code{\link{Wk}}
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
#' sam <- S.SI(N,n)
#' # The information about the units in the sample is stored in an object called data
#' data <- Lucy[sam,]
#' attach(data)
#' names(data)
#'
#' ########### common mean model
#'
#' estima<-data.frame(Income, Employees, Taxes)
#' x <- rep(1,n)
#' model <- E.Beta(N, n, estima, x, ck=1,b0=FALSE)
#' b <- t(as.matrix(model[1,,]))
#' tx <- c(N)
#' GREG.SI(N,n,estima,x,tx, b, b0=FALSE)
#'
#' ########### common ratio model
#'
#' estima<-data.frame(Income)
#' x <- data.frame(Employees)
#' model <- E.Beta(N, n, estima, x, ck=x,b0=FALSE)
#' b <- t(as.matrix(model[1,,]))
#' tx <- sum(Lucy$Employees)
#' GREG.SI(N,n,estima,x,tx, b, b0=FALSE)
#'
#' ########### Simple regression model without intercept
#'
#' estima<-data.frame(Income, Employees)
#' x <- data.frame(Taxes)
#' model <- E.Beta(N, n, estima, x, ck=1,b0=FALSE)
#' b <- t(as.matrix(model[1,,]))
#' tx <- sum(Lucy$Taxes)
#' GREG.SI(N,n,estima,x,tx, b, b0=FALSE)
#'
#' ########### Multiple regression model without intercept
#'
#' estima<-data.frame(Income)
#' x <- data.frame(Employees, Taxes)
#' model <- E.Beta(N, n, estima, x, ck=1, b0=FALSE)
#' b <- as.matrix(model[1,,])
#' tx <- c(sum(Lucy$Employees), sum(Lucy$Taxes))
#' GREG.SI(N,n,estima,x,tx, b, b0=FALSE) 
#'
#' ########### Simple regression model with intercept
#'
#' estima<-data.frame(Income, Employees)
#' x <- data.frame(Taxes)
#' model <- E.Beta(N, n, estima, x, ck=1,b0=TRUE) 
#' b <- as.matrix(model[1,,])
#' tx <- c(N, sum(Lucy$Taxes))
#' GREG.SI(N,n,estima,x,tx, b, b0=TRUE) 
#'
#' ########### Multiple regression model with intercept
#'
#' estima<-data.frame(Income)                               
#' x <- data.frame(Employees, Taxes)
#' model <- E.Beta(N, n, estima, x, ck=1,b0=TRUE)
#' b <- as.matrix(model[1,,])
#' tx <- c(N, sum(Lucy$Employees), sum(Lucy$Taxes))            
#' GREG.SI(N,n,estima,x,tx, b, b0=TRUE) 
#'
#' ####################################################################
#' ## Example 2: Linear models with discrete auxiliary information
#' ####################################################################
#'
#' # Draws a simple random sample without replacement
#' data(Lucy)
#'
#' N <- dim(Lucy)[1]
#' n <- 400
#' sam <- S.SI(N,n)
#' # The information about the units in the sample is stored in an object called data
#' data <- Lucy[sam,]
#' attach(data)
#' names(data)
#'
#' # The auxiliary information is discrete type
#' Doma<-Domains(Level)
#'
#' ########### Poststratified common mean model
#'
#' estima<-data.frame(Income, Employees, Taxes)
#' model <- E.Beta(N, n, estima, Doma, ck=1,b0=FALSE)
#' b <- t(as.matrix(model[1,,]))
#' tx <- colSums(Domains(Lucy$Level))
#' GREG.SI(N,n,estima,Doma,tx, b, b0=FALSE) 
#'
#' ########### Poststratified common ratio model 
#'
#' estima<-data.frame(Income, Employees)
#' x <- Doma*Taxes
#' model <- E.Beta(N, n, estima, x ,ck=1,b0=FALSE)
#' b <- as.matrix(model[1,,])
#' tx <- colSums(Domains(Lucy$Level)*Lucy$Taxes)
#' GREG.SI(N,n,estima,x,tx, b, b0=FALSE) 
#'
#' ######################################################################
#' ## Example 3: Domains estimation trough the postestratified estimator
#' ######################################################################
#'
#' # Draws a simple random sample without replacement
#' data(Lucy)
#'
#' N <- dim(Lucy)[1]
#' n <- 400
#' sam <- S.SI(N,n)
#' # The information about the units in the sample is stored in an object called data
#' data <- Lucy[sam,]
#' attach(data)
#' names(data)
#'
#' # The auxiliary information is discrete type
#' Doma<-Domains(Level)
#'
#' ########### Poststratified common mean model for the 
#'   # Income total in each poststratum ###################
#'
#' estima<-Doma*Income
#' model <- E.Beta(N, n, estima, Doma, ck=1, b0=FALSE)
#' b <- t(as.matrix(model[1,,]))
#' tx <- colSums(Domains(Lucy$Level))
#' GREG.SI(N,n,estima,Doma,tx, b, b0=FALSE) 
#'
#' ########### Poststratified common mean model for the 
#'   # Employees total in each poststratum ###################
#'
#' estima<-Doma*Employees
#' model <- E.Beta(N, n, estima, Doma, ck=1,b0=FALSE)
#' b <- t(as.matrix(model[1,,]))
#' tx <- colSums(Domains(Lucy$Level))
#' GREG.SI(N,n,estima,Doma,tx, b, b0=FALSE) 
#'
#' ########### Poststratified common mean model for the 
#'   # Taxes total in each poststratum ###################
#'
#' estima<-Doma*Taxes
#' model <- E.Beta(N, n, estima, Doma, ck=1, b0=FALSE)
#' b <- t(as.matrix(model[1,,]))
#' tx <- colSums(Domains(Lucy$Level))
#' GREG.SI(N,n,estima,Doma,tx, b, b0=FALSE)

GREG.SI <- function(N, n, y, x, tx, b, b0 = FALSE) {
  y   <- as.data.frame(y)
  x   <- as.matrix(x)
  pik <- rep(n/N, n)
  dk  <- 1/pik
  if (b0 == TRUE) {
    x <- as.matrix(cbind(1, x))
  }
  Total <- matrix(NA, nrow = 3, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE")
  colnames(Total) <- names(y)
  for (k in 1:dim(y)[2]) {
    xHT <- t(x) %*% dk
    yHT <- sum(y[, k] * dk)
    ty  <- yHT + (tx - t(xHT)) %*% as.matrix(b[, k])
    e   <- y[, k] - (x %*% as.matrix(b[, k]))
    Vty <- (N^2) * (1 - (n/N)) * var(e)/(n)
    CVe <- 100 * sqrt(Vty)/ty
    Total[, k] <- c(ty, sqrt(Vty), CVe)
  }
  return(Total)
}