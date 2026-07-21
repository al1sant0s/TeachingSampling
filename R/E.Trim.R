#' @export
#' 
#' @title
#' Weight Trimming and Redistribution
#' @description 
#' This function performs a method of trimming sampling weights based on the 
#' evenly redistribution of the net ammount of weight loss among units whose
#' weights were not trimmed. This way, the sum of the timmed sampling weights 
#' remains the same as the original weights.
#' @return
#' A numeric vector of trimmed sampling weights with the same length as
#' \code{dk} and the same sum as the original weights.
#' @details
#' Weights below \code{L} or above \code{U} are replaced by the respective
#' bound. The net weight lost by trimming is then redistributed evenly among
#' all units whose weights were not trimmed. The process iterates until no
#' weight falls outside \code{[L, U]}, ensuring that the sum of the trimmed
#' weights equals the sum of the original weights.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com> with
#'   contributions from Javier Nunez <javier_nunez at inec.gob.ec>
#' @param dk Vector of original sampling weights.
#' @param L Lower bound for weights.
#' @param U Upper bound for weights.
#' 
#' @references 
#' Valliant, R. et. al. (2013), \emph{Practical Tools for Designing and
#' Weigthing Survey Samples}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas 
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#' 
#' @examples
#' 
#' # Example 1
#' dk <- c(1, 1, 1, 10) 
#' summary(dk)
#' L <- 1
#' U <- 3.5 * median(dk)
#' dkTrim <- E.Trim(dk, L, U)
#' sum(dk)
#' sum(dkTrim)
#' 
#' # Example 2
#' dk <- rnorm(1000, 10, 10)
#' L <- 1
#' U <- 3.5 * median(dk)
#' dkTrim <- E.Trim(dk, L, U)
#' sum(dk)
#' sum(dkTrim)
#' summary(dk)
#' summary(dkTrim)
#' hist(dk)
#' hist(dkTrim)

E.Trim <- function(dk, L, U){
  i <- dk > U | dk < L
  dkL <- ifelse(dk < L, L, dk)
  dkLU <- ifelse(dkL > U, U, dkL)
  s <- sum(dk - dkLU)
  
  while (sum(i) != 0) {
    dk <- dkLU + (s/(length(dk) - sum(i))) * (1 - i)
    i <- dk > U | dk < L
    dkL <- ifelse(dk < L, L, dk)
    dkLU <- ifelse(dkL > U, U, dkL)
    s <- sum(dk - dkLU)
  }
  return(dkLU)
}