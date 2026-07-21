#' @export
#'
#' @title
#' Hansen-Hurwitz Estimator of the Population Total
#' @description
#' Computes the Hansen-Hurwitz (HH) estimator of the population total under
#' a with-replacement sampling design, given the sample observations and
#' their selection probabilities.
#' @return
#' A numeric vector or matrix with the estimated total for each variable
#' of interest.
#' @details
#' The Hansen-Hurwitz estimator is:
#' \deqn{\hat{t}_{HH} = \frac{1}{m}\sum_{i=1}^m \frac{y_i}{p_i}}
#' where \eqn{p_i} is the selection probability of the \eqn{i}-th draw
#' and \eqn{m} is the number of draws. This estimator is design-unbiased
#' under any with-replacement sampling design.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param y Vector or matrix of values of the variable(s) of interest for
#'   units in the sample (with possible repetitions).
#' @param pk Vector of selection probabilities for each draw in the sample.
#'
#' @references
#' Hansen, M.H. and Hurwitz, W.N. (1943). On the theory of sampling from
#' finite populations. \emph{Annals of Mathematical Statistics}, 14, 333-362.\cr
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.
#'
#' @seealso \code{\link{E.PPS}}, \code{\link{HT}}, \code{\link{S.PPS}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' # Vector U contains the label of a population of size N=5
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' # Vectors y1 and y2 give the values of the variables of interest
#' y1<-c(32, 34, 46, 89, 35)
#' y2<-c(1,1,1,0,0)
#' y3<-cbind(y1,y2)
#' # The population size is N=5
#' N <- length(U)
#' # The sample size is m=2
#' m <- 2
#' # pk is the probability of selection of every single unit
#' pk <- c(0.35, 0.225, 0.175, 0.125, 0.125)
#' # Selection of a random sample with replacement
#' sam <- sample(5,2, replace=TRUE, prob=pk)
#' # The selected sample is
#' U[sam]
#' # The values of the variables of interest for the units in the sample
#' y1[sam]
#' y2[sam]
#' y3[sam,]
#' # The Hansen-Hurwitz estimator
#' HH(y1[sam],pk[sam])
#' HH(y2[sam],pk[sam])
#' HH(y3[sam,],pk[sam])
#'
#'
#' ############
#' ## Example 2
#' ############
#' # Uses the Lucy data to draw a simple random sample with replacement
#' data(Lucy)
#' attach(Lucy)
#'
#' N <- dim(Lucy)[1]
#' m <- 400
#' sam <- sample(N,m,replace=TRUE)
#' # The vector of selection probabilities of units in the sample
#' pk <- rep(1/N,m)
#' # The information about the units in the sample is stored in an object called data
#' data <- Lucy[sam,]
#' attach(data)
#' names(data)
#' # The variables of interest are: Income, Employees and Taxes
#' # This information is stored in a data frame called estima
#' estima <- data.frame(Income, Employees, Taxes)
#' HH(estima, pk)
#'
#' ################################################################
#' ## Example 3 HH is unbiased for with replacement sampling designs
#' ################################################################
#'
#' # Vector U contains the label of a population of size N=5
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' # Vector y1 and y2 are the values of the variables of interest
#' y<-c(32, 34, 46, 89, 35)
#' # The population size is N=5
#' N <- length(U)
#' # The sample size is m=2
#' m <- 2
#' # pk is the probability of selection of every single unit
#' pk <- c(0.35, 0.225, 0.175, 0.125, 0.125)
#' # p is the probability of selection of every possible sample
#' p <- p.WR(N,m,pk)
#' p
#' sum(p)
#' # The sample membership matrix for random size without replacement sampling designs
#' Ind <- nk(N,m)
#' Ind
#' # The support with the values of the elements
#' Qy <- SupportWR(N,m, ID=y)                 
#' Qy
#' # The support with the values of the elements
#' Qp <- SupportWR(N,m, ID=pk)                 
#' Qp
#' # The HT estimates for every single sample in the support
#' HH1 <- HH(Qy[1,], Qp[1,])[1,]
#' HH2 <- HH(Qy[2,], Qp[2,])[1,]
#' HH3 <- HH(Qy[3,], Qp[3,])[1,] 
#' HH4 <- HH(Qy[4,], Qp[4,])[1,] 
#' HH5 <- HH(Qy[5,], Qp[5,])[1,] 
#' HH6 <- HH(Qy[6,], Qp[6,])[1,] 
#' HH7 <- HH(Qy[7,], Qp[7,])[1,]
#' HH8 <- HH(Qy[8,], Qp[8,])[1,]
#' HH9 <- HH(Qy[9,], Qp[9,])[1,]
#' HH10 <- HH(Qy[10,], Qp[10,])[1,]
#' HH11 <- HH(Qy[11,], Qp[11,])[1,]
#' HH12 <- HH(Qy[12,], Qp[12,])[1,]
#' HH13 <- HH(Qy[13,], Qp[13,])[1,]
#' HH14 <- HH(Qy[14,], Qp[14,])[1,]
#' HH15 <- HH(Qy[15,], Qp[15,])[1,]
#' # The HT estimates arranged in a vector
#' Est <- c(HH1, HH2, HH3, HH4, HH5, HH6, HH7, HH8, HH9, HH10, HH11, HH12, HH13,
#' HH14, HH15)
#' Est
#' # The HT is actually desgn-unbiased
#' data.frame(Ind, Est, p)
#' sum(Est*p)
#' sum(y)

HH <- function(y, pk) {
  y  <- t(as.matrix(y))
  pk <- as.matrix(pk)
  m  <- length(pk)
  result <- (1/m) * (y %*% (1/pk))
  result
}