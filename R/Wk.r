#' @export
#'
#' @title
#' GREG Generalised Weights
#' @description
#' Computes the generalised regression (GREG) weights for each unit in the
#' sample. These weights incorporate both the sampling design weights and a
#' calibration adjustment based on known population totals of auxiliary
#' variables.
#' @return
#' A numeric vector of length \code{n} with the GREG weight for each unit
#' in the sample.
#' @details
#' The GREG weight for unit \eqn{k} is:
#' \deqn{w_k = \frac{1}{\pi_k} + \mathbf{x}_k^T
#' \left(\sum_s \frac{v_k \mathbf{x}_k \mathbf{x}_k^T}{\pi_k}\right)^{-1}
#' (\mathbf{t}_x - \hat{\mathbf{t}}_{x,\pi})}
#' where \eqn{v_k = 1/(\pi_k c_k)} and \eqn{c_k} is a variance-stabilising
#' constant. The GREG estimator is then \eqn{\hat{t}_{GREG} = \sum_s w_k y_k}.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param x Vector or matrix of auxiliary variables observed in the sample.
#' @param tx Vector of known population totals of the auxiliary variables.
#' @param Pik Vector of first-order inclusion probabilities for each unit
#'   in the sample.
#' @param ck Vector of variance-stabilising constants. Typically \code{ck = 1}
#'   (homoscedastic) or \code{ck = x} (heteroscedastic).
#' @param b0 Logical. If \code{TRUE}, an intercept column is prepended to
#'   \code{x}. Default is \code{FALSE}.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{GREG.SI}}, \code{\link{E.Beta}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' # Without replacement sampling
#' # Vector U contains the label of a population of size N=5
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' # Vector x is the auxiliary information and y is the variables of interest
#' x<-c(32, 34, 46, 89, 35)
#' y<-c(52, 60, 75, 100, 50)
#' # pik is some vector of inclusion probabilities in the sample
#' # In this case the sample size is equal to the population size
#' pik<-rep(1,5)
#' w1<-Wk(x,tx=236,pik,ck=1,b0=FALSE)
#' sum(x*w1)
#' # Draws a sample size without replacement
#' sam <- sample(5,2)
#' pik <- c (0.8,0.2,0.2,0.5,0.3)
#' # The auxiliary information an variable of interest in the selected smaple
#' x.s<-x[sam]
#' y.s<-y[sam]
#' # The vector of inclusion probabilities in the selected smaple
#' pik.s<-pik[sam]
#' # Calibration weights under some specifics model
#' w2<-Wk(x.s,tx=236,pik.s,ck=1,b0=FALSE)
#' sum(x.s*w2)
#'
#' w3<-Wk(x.s,tx=c(5,236),pik.s,ck=1,b0=TRUE)
#' sum(w3)
#' sum(x.s*w3)
#'
#' w4<-Wk(x.s,tx=c(5,236),pik.s,ck=x.s,b0=TRUE)
#' sum(w4)
#' sum(x.s*w4)
#'
#' w5<-Wk(x.s,tx=236,pik.s,ck=x.s,b0=FALSE)
#' sum(x.s*w5)
#'
#' ######################################################################
#' ## Example 2: Linear models involving continuous auxiliary information
#' ######################################################################
#'
#' # Draws a simple random sample without replacement
#' data(Lucy)
#' attach(Lucy)
#'
#' N <- dim(Lucy)[1]
#' n <- 400
#' Pik <- rep(n/N, n)
#' sam <- S.SI(N,n)
#' # The information about the units in the sample is stored in an object called data
#' data <- Lucy[sam,]
#' attach(data)
#' names(data)
#'
#' ########### common ratio model ###################
#'
#' estima<-data.frame(Income)
#' x <- Employees
#' tx <- sum(Lucy$Employees)
#' w <- Wk(x, tx, Pik, ck=1, b0=FALSE)
#' sum(x*w)
#' tx
#' # The calibration estimation
#' colSums(estima*w)
#'
#' ########### Simple regression model without intercept ###################
#'
#' estima<-data.frame(Income, Employees)
#' x <- Taxes
#' tx <- sum(Lucy$Taxes)
#' w<-Wk(x,tx,Pik,ck=x,b0=FALSE)
#' sum(x*w)
#' tx
#' # The calibration estimation
#' colSums(estima*w)
#'
#' ########### Multiple regression model without intercept ###################
#'
#' estima<-data.frame(Income)
#' x <- cbind(Employees, Taxes)
#' tx <- c(sum(Lucy$Employees), sum(Lucy$Taxes))
#' w <- Wk(x,tx,Pik,ck=1,b0=FALSE)
#' sum(x[,1]*w)
#' sum(x[,2]*w)
#' tx
#' # The calibration estimation
#' colSums(estima*w)
#'
#' ########### Simple regression model with intercept ###################
#'
#' estima<-data.frame(Income, Employees)
#' x <- Taxes
#' tx <- c(N,sum(Lucy$Taxes))
#' w <- Wk(x,tx,Pik,ck=1,b0=TRUE)
#' sum(1*w)
#' sum(x*w)
#' tx
#' # The calibration estimation
#' colSums(estima*w)
#'
#' ########### Multiple regression model with intercept ###################
#'
#' estima<-data.frame(Income)
#' x <- cbind(Employees, Taxes)
#' tx <- c(N, sum(Lucy$Employees), sum(Lucy$Taxes))
#' w <- Wk(x,tx,Pik,ck=1,b0=TRUE)
#' sum(1*w)
#' sum(x[,1]*w)
#' sum(x[,2]*w)
#' tx
#' # The calibration estimation
#' colSums(estima*w)
#'
#' ####################################################################
#' ## Example 3: Linear models involving discrete auxiliary information
#' ####################################################################
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
#' # Vector of inclusion probabilities for units in the selected sample
#' Pik<-rep(n/N,n)
#' # The auxiliary information is discrete type
#' Doma<-Domains(Level)
#'
#' ########### Poststratified common mean model ###################
#'
#' estima<-data.frame(Income, Employees, Taxes)
#' tx <- colSums(Domains(Lucy$Level))
#' w <- Wk(Doma,tx,Pik,ck=1,b0=FALSE)
#' sum(Doma[,1]*w)
#' sum(Doma[,2]*w)
#' sum(Doma[,3]*w)
#' tx
#' # The calibration estimation
#' colSums(estima*w)
#'
#' ########### Poststratified common ratio model ###################
#'
#' estima<-data.frame(Income, Employees)
#' x<-Doma*Taxes
#' tx <- colSums(Domains(Lucy$Level))
#' w <- Wk(x,tx,Pik,ck=1,b0=FALSE)
#' sum(x[,1]*w)
#' sum(x[,2]*w)
#' sum(x[,3]*w)
#' tx
#' # The calibration estimation
#' colSums(estima*w)

Wk <- function(x, tx, Pik, ck, b0 = FALSE) {
  if (b0 == TRUE)  x <- as.matrix(cbind(1, x))
  if (b0 == FALSE) x <- as.matrix(x)
  tx    <- as.matrix(tx)
  txpi  <- as.matrix(t(x) %*% (1/Pik))
  V     <- 1/(Pik * ck)
  result <- (1/Pik) + ((V * x) %*% solve(t(V * x) %*% x) %*% (tx - txpi))
  return(result)
}