#' @export
#'
#' @title
#' Complete Sampling Support for All Sample Sizes
#' @description
#' Enumerates all possible non-empty subsets of a population of size \code{N},
#' covering all sample sizes from 1 to \code{N}. The result includes the
#' empty set as the first row.
#' @return
#' A matrix with \eqn{2^N} rows and \code{N} columns. Each row is one subset,
#' with \code{NA} used as padding for subsets smaller than \code{N}. The first
#' row represents the empty set (all zeros).
#' @details
#' This function stacks the outputs of \code{\link{Support}} for all sample
#' sizes \eqn{n = 1, \ldots, N}. It is only feasible for small populations
#' (\code{N <= 10}) due to exponential growth.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param N Population size. Recommended \code{N <= 10}.
#' @param ID Optional vector of population labels of length \code{N}.
#'   If provided, labels replace integer indices in the output.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{Support}}, \code{\link{IkRS}}
#'
#' @examples
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' N <- length(U)
#' # Complete support for all sample sizes
#' SupportRS(N)
#' # Labeled support
#' SupportRS(N, ID = U)

SupportRS <- function(N, ID = FALSE) {
  sam <- matrix(NA, ncol = N, nrow = 1)
  for (k in 1:N) {
    sam <- rbind(sam,
                 cbind(Support(N, k),
                       matrix(NA, ncol = N - k, nrow = choose(N, k))))
  }
  if (is.logical(ID) == TRUE) return(sam)
  else {
    sam <- matrix(ID[SupportRS(N)], nrow = 2^N)
    return(sam)
  }
}