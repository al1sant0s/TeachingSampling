#' @export
#'
#' @title
#' Estimation of the Population Total under Two Stage Simple Random Sampling Without Replacement
#' @description
#' Computes the Horvitz-Thompson estimator of the population total under a
#' two-stage simple random sampling without replacement design, where both
#' Primary Sampling Units (PSUs) and Secondary Sampling Units (SSUs) are
#' selected by simple random sampling without replacement.
#' @return
#' A matrix with four rows and one column per variable of interest:
#' \itemize{
#'   \item \code{Estimation}: Estimated population total.
#'   \item \code{Standard Error}: Estimated standard error of the total.
#'   \item \code{CVE}: Estimated coefficient of variation (in percentage).
#'   \item \code{DEFF}: Design effect with respect to simple random sampling.
#' }
#' @details
#' The variance estimator decomposes into two components: the between-PSU
#' component and the within-PSU component, following the classical two-stage
#' variance decomposition of Sarndal et al. (1992).
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param NI Population size of Primary Sampling Units (PSUs).
#' @param nI Sample size of Primary Sampling Units (PSUs).
#' @param Ni Vector of population sizes of Secondary Sampling Units within
#'   each selected PSU.
#' @param ni Vector of sample sizes of Secondary Sampling Units within
#'   each selected PSU.
#' @param y Vector, matrix or data frame containing the values of the
#'   variables of interest for every unit in the selected sample.
#' @param PSU Vector identifying the PSU membership of each unit in the sample.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{E.1SI}}, \code{\link{E.UC}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' # Uses Lucy data to draw a twostage simple random sample 
#' # accordind to a 2SI design. Zone is the clustering variable
#' data(Lucy)
#' attach(Lucy)
#' summary(Zone)
#' # The population of clusters or Primary Sampling Units
#' UI<-c("A","B","C","D","E")
#' NI <- length(UI)
#' # The sample size is nI=3
#' nI <- 3
#' # Selects the sample of PSUs
#' samI<-S.SI(NI,nI)
#' dataI<-UI[samI]
#' dataI   
#' # The sampling frame of Secondary Sampling Unit is saved in Lucy1 ... Lucy3
#' Lucy1<-Lucy[which(Zone==dataI[1]),]
#' Lucy2<-Lucy[which(Zone==dataI[2]),]
#' Lucy3<-Lucy[which(Zone==dataI[3]),]
#' # The size of every single PSU
#' N1<-dim(Lucy1)[1]
#' N2<-dim(Lucy2)[1]
#' N3<-dim(Lucy3)[1]
#' Ni<-c(N1,N2,N3)
#' # The sample size in every PSI is 135 Secondary Sampling Units
#' n1<-135
#' n2<-135
#' n3<-135
#' ni<-c(n1,n2,n3)
#' # Selects a sample of Secondary Sampling Units inside the PSUs
#' sam1<-S.SI(N1,n1)
#' sam2<-S.SI(N2,n2)
#' sam3<-S.SI(N3,n3)
#' # The information about each Secondary Sampling Unit in the PSUs
#' # is saved in data1 ... data3
#' data1<-Lucy1[sam1,]
#' data2<-Lucy2[sam2,]
#' data3<-Lucy3[sam3,]
#' # The information about each unit in the final selected sample is saved in data
#' data<-rbind(data1, data2, data3)
#' attach(data)
#' # The clustering variable is Zone
#' Cluster <- as.factor(as.integer(Zone))
#' # The variables of interest are: Income, Employees and Taxes
#' # This information is stored in a data frame called estima
#' estima <- data.frame(Income, Employees, Taxes)
#' # Estimation of the Population total
#' E.2SI(NI,nI,Ni,ni,estima,Cluster)
#'
#' ########################################################
#' ## Example 2 Total Census to the entire population
#' ########################################################
#' # Uses Lucy data to draw a cluster random sample
#' # accordind to a SI design ...
#' # Zone is the clustering variable
#' data(Lucy)
#' attach(Lucy)
#' summary(Zone)
#' # The population of clusters
#' UI<-c("A","B","C","D","E")
#' NI <- length(UI)
#' # The sample size equals to the population size of PSU
#' nI <- NI
#' # Selects every single PSU
#' samI<-S.SI(NI,nI)
#' dataI<-UI[samI]
#' dataI   
#' # The sampling frame of Secondary Sampling Unit is saved in Lucy1 ... Lucy5
#' Lucy1<-Lucy[which(Zone==dataI[1]),]
#' Lucy2<-Lucy[which(Zone==dataI[2]),]
#' Lucy3<-Lucy[which(Zone==dataI[3]),]
#' Lucy4<-Lucy[which(Zone==dataI[4]),]
#' Lucy5<-Lucy[which(Zone==dataI[5]),]
#' # The size of every single PSU
#' N1<-dim(Lucy1)[1]
#' N2<-dim(Lucy2)[1]
#' N3<-dim(Lucy3)[1]
#' N4<-dim(Lucy4)[1]
#' N5<-dim(Lucy5)[1]
#' Ni<-c(N1,N2,N3,N4,N5)
#' # The sample size of Secondary Sampling Units equals to the size of each PSU
#' n1<-N1
#' n2<-N2
#' n3<-N3
#' n4<-N4
#' n5<-N5
#' ni<-c(n1,n2,n3,n4,n5)
#' # Selects every single Secondary Sampling Unit inside the PSU
#' sam1<-S.SI(N1,n1)
#' sam2<-S.SI(N2,n2)
#' sam3<-S.SI(N3,n3)
#' sam4<-S.SI(N4,n4)
#' sam5<-S.SI(N5,n5)
#' # The information about each unit in the cluster is saved in Lucy1 ... Lucy5
#' data1<-Lucy1[sam1,]
#' data2<-Lucy2[sam2,]
#' data3<-Lucy3[sam3,]
#' data4<-Lucy4[sam4,]
#' data5<-Lucy5[sam5,]
#' # The information about each Secondary Sampling Unit
#' # in the sample (census) is saved in data
#' data<-rbind(data1, data2, data3, data4, data5)
#' attach(data)
#' # The clustering variable is Zone
#' Cluster <- as.factor(as.integer(Zone))
#' # The variables of interest are: Income, Employees and Taxes
#' # This information is stored in a data frame called estima
#' estima <- data.frame(Income, Employees, Taxes)
#' # Estimation of the Population total
#' E.2SI(NI,nI,Ni,ni,estima,Cluster)
#' # Sampling error is null

E.2SI <- function(NI, nI, Ni, ni, y, PSU) {
  y <- cbind(1, y)
  y <- as.data.frame(y)
  names(y)[1] <- "N"
  PSU <- as.factor(PSU)
  
  Total <- matrix(NA, nrow = 4, ncol = dim(y)[2])
  rownames(Total) <- c("Estimation", "Standard Error", "CVE", "DEFF")
  colnames(Total) <- names(y)
  
  f <- ni/Ni
  F <- nI/NI
  
  for (k in 1:dim(y)[2]) {
    ysum <- tapply(y[, k], PSU, sum)
    s2i  <- tapply(y[, k], PSU, var)
    ti   <- (1/f) * ysum
    ty   <- (1/F) * sum(ti)
    part.1 <- NI^2/nI * (1 - F) * var(ti)
    part.2 <- NI/nI * sum(Ni^2/ni * (1 - f) * s2i)
    Vty    <- part.1 + part.2
    CVe    <- 100 * sqrt(Vty)/ty
    n      <- length(y[, k])
    N      <- (NI/nI) * sum(Ni)
    VMAS   <- (N^2) * (1 - (n/N)) * var(y[, k])/(n)
    DEFF   <- Vty/VMAS
    Total[, k] <- c(ty, sqrt(Vty), CVe, DEFF)
  }
  return(Total)
}