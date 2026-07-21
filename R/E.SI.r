#' @export
#'
#' @title
#' Estimation of the Population Total under Simple Random Sampling Without
#' Replacement
#' @description
#' Computes the Horvitz-Thompson estimator of the population total under a
#' simple random sampling without replacement (SI) design.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error of the total.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect (always 1 under SI, included for
#'     consistency with other estimators).
#' }
#' @details
#' Under simple random sampling without replacement, the Horvitz-Thompson
#' estimator reduces to \eqn{\hat{t}_y = N \bar{y}_s}, the expansion
#' estimator. The design effect is always 1 because SI is the reference design.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size.
#' @param n Sample size.
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every unit in the selected sample.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{S.SI}}, \code{\link{E.STSI}}, \code{\link{GREG.SI}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' # Uses the Lucy data to draw a random sample of units according to a SI design
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
#' # The variables of interest are: Income, Employees and Taxes
#' # This information is stored in a data frame called estima
#' estima <- data.frame(Income, Employees, Taxes)
#' E.SI(N,n,estima)
#'
#' ############
#' ## Example 2
#' ############
#' # Following with Example 1. The variable SPAM is a domain of interest
#' Doma <- Domains(SPAM)
#' # This function allows to estimate the size of each domain in SPAM
#' estima <- data.frame(Doma)
#' E.SI(N,n,Doma)
#'
#' ############
#' ## Example 3
#' ############
#' # Following with Example 1. The variable SPAM is a domain of interest
#' Doma <- Domains(SPAM)
#' # This function allows to estimate the parameters of the variables of interest 
#' # for every category in the domain SPAM
#' estima <- data.frame(Income, Employees, Taxes)
#' SPAM.no <-  cbind(Doma[,1], estima*Doma[,1])
#' SPAM.yes <- cbind(Doma[,1], estima*Doma[,2])
#' # Before running the following lines, notice that:
#' # The first column always indicates the population size
#' # The second column is an estimate of the size of the category in the domain SPAM
#' # The remaining columns estimates the parameters of interest 
#' # within the corresponding category in the domain SPAM
#' E.SI(N,n,SPAM.no)
#' E.SI(N,n,SPAM.yes)
#'
#' ############
#' ## Example 4
#' ############
#' # Following with Example 1. The variable SPAM is a domain of interest 
#' # and the variable ISO is a populational subgroup of interest
#' Doma <- Domains(SPAM)
#' estima <- Domains(Zone)
#' # Before running the following lines, notice that:
#' # The first column indicates wheter the unit 
#' # belongs to the first category of SPAM or not
#' # The remaining columns indicates wheter the unit 
#' # belogns to the categories of Zone
#' SPAM.no <-  data.frame(SpamNO=Doma[,1], Zones=estima*Doma[,1])
#' # Before running the following lines, notice that:
#' # The first column indicates wheter the unit 
#' # belongs to the second category of SPAM or not
#' # The remaining columns indicates wheter the unit 
#' # belogns to the categories of Zone
#' SPAM.yes <- data.frame(SpamYES=Doma[,2], Zones=estima*Doma[,2])
#' # Before running the following lines, notice that:
#' # The first column always indicates the population size
#' # The second column is an estimate of the size of the 
#' # first category in the domain SPAM
#' # The remaining columns estimates the size of the categories 
#' # of Zone within the corresponding category of SPAM
#' # Finnaly, note that the sum of the point estimates of the last 
#' # two columns gives exactly the point estimate in the second column
#' E.SI(N,n,SPAM.no)
#' # Before running the following lines, notice that:
#' # The first column always indicates the population size
#' # The second column is an estimate of the size of the 
#' # second category in the domain SPAM
#' # The remaining columns estimates the size of the categories 
#' # of Zone within the corresponding category of SPAM
#' # Finnaly, note that the sum of the point estimates of the last two 
#' # columns gives exactly the point estimate in the second column
#' E.SI(N,n,SPAM.yes)

E.SI <- function(N, n, y) {
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  pik <- matrix(n/N, nrow = n, ncol = 1)
  dk  <- 1/pik
  for (k in 1:dim(y)[2]) {
    ty   <- sum(y[, k] * dk)
    Vty  <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    CVe  <- 100 * sqrt(Vty)/ty
    VMAS <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}