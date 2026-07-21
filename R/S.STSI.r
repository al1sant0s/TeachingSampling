#' @export
#'
#' @title
#' Stratified Simple Random Sampling Without Replacement
#' @description
#' Draws a stratified simple random sample without replacement from a finite
#' population. Within each stratum, units are selected by simple random
#' sampling without replacement.
#' @return
#' A sorted vector of population indices of the selected units, of length
#' \code{sum(nh)}.
#' @details
#' The function selects \code{nh[h]} units from stratum \eqn{h} using
#' \code{base::sample}, and returns all selected indices sorted in ascending
#' order. Use \code{\link{E.STSI}} to estimate population totals from this
#' sample.
#' @author Hugo Andres Gutierrez Rojas <hagutierrezro at gmail.com>
#' @param S Vector of length \code{N} identifying the stratum membership of
#'   each unit in the population.
#' @param Nh Integer vector of length \code{H} with the population size of
#'   each stratum.
#' @param nh Integer vector of length \code{H} with the sample size of each
#'   stratum. Must satisfy \code{nh[h] <= Nh[h]} for all \code{h}.
#'
#' @references
#' Sarndal, C-E. and Swensson, B. and Wretman, J. (1992),
#' \emph{Model Assisted Survey Sampling}. Springer.\cr
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#'
#' @seealso \code{\link{E.STSI}}, \code{\link{S.SI}}, \code{\link{S.STpiPS}}
#'
#' @examples
#' ############
#' ## Example 1
#' ############
#' U <- c("Yves", "Ken", "Erik", "Sharon", "Leslie")
#' Strata <- c("A", "A", "A", "B", "B")
#' Nh <- c(3, 2)
#' nh <- c(2, 1)
#' sam <- S.STSI(Strata, Nh, nh)
#' sam
#' U[sam]
#' ############
#' ## Example 2
#' ############
#' data(Lucy)
#' attach(Lucy)
#' N1 <- summary(Level)[[1]]
#' N2 <- summary(Level)[[2]]
#' N3 <- summary(Level)[[3]]
#' Nh <- c(N1, N2, N3)
#' nh <- c(70, 100, 200)
#' sam <- S.STSI(Level, Nh, nh)
#' data <- Lucy[sam, ]
#' dim(data)

S.STSI <- function(S, Nh, nh) {
  S   <- as.factor(S)
  S   <- as.factor(as.integer(S))
  cum <- cumsum(nh)
  sam <- matrix(0, sum(nh))
  for (k in 1:length(nh)) {
    h     <- which(S == k)
    sam.h <- sample(Nh[k], nh[k])
    if (k == 1) sam[1:nh[k]] <- h[sam.h]
    if (k > 1)  sam[(cum[k-1]+1):(cum[k])] <- h[sam.h]
  }
  sort(sam)
}