#' @export
#'
#' @title
#' Kish Allocation for Stratified Sampling
#' @description
#' Computes the optimal sample size allocation across strata using the
#' Kish (1992) compromise allocation method, which interpolates between
#' uniform and proportional allocation through a design effect parameter \code{I}.
#'
#' @param n     Integer. Total desired sample size.
#' @param N_h   Named numeric vector. Population sizes for each stratum
#'   \eqn{h = 1, \ldots, H}.
#' @param I     Non-negative numeric. Intraclass correlation coefficient (ICC)
#'   or design effect parameter controlling the allocation:
#'   \itemize{
#'     \item \code{I = 0}      → Uniform allocation (equal sample per stratum).
#'     \item \code{I = Inf}    → Proportional allocation (proportional to \eqn{N_h}).
#'     \item \code{0 < I < Inf} → Compromise between uniform and proportional.
#'     \item Recommended value: \code{I = 0.5} (Kish, 1992).
#'   }
#'
#' @return A named integer vector of length \eqn{H} with the allocated sample
#'   sizes per stratum. The values sum to approximately \code{n} (rounding may
#'   cause a difference of ±1).
#'
#' @details
#' The Kish compromise allocation assigns sample sizes as:
#' \deqn{
#'   n_h = n \cdot \frac{\sqrt{I \, W_h^2 + H^{-2}}}
#'   {\sum_{h=1}^{H} \sqrt{I \, W_h^2 + H^{-2}}}
#' }
#' where \eqn{W_h = N_h / N} is the stratum weight and \eqn{H} is the number
#' of strata. This formulation nests two classical allocations as limiting
#' cases: when \eqn{I = 0} the numerator reduces to \eqn{1/H} (uniform),
#' and as \eqn{I \to \infty} it is dominated by \eqn{W_h} (proportional).
#'
#' @references
#' Kish, L. (1992). Weighting for unequal \eqn{P_i}.
#' \emph{Journal of Official Statistics}, 8(2), 183–200.
#'
#' @author Yury Vanessa Ochoa Montes <yury.ochoa@urosario.edu.co>
#'
#' @seealso
#' \code{\link{E.STSI}} for estimation under stratified sampling,
#' \code{\link{S.STSI}} for stratified simple random sampling.
#'
#' @examples
#' N_h <- c(
#'   Corozal     = 41847,
#'   Orange_Walk = 48175,
#'   Belize      = 57658,
#'   Cayo        = 78473,
#'   Stann_Creek = 31347,
#'   Toledo      = 31711
#' )
#'
#' # Uniform allocation (I = 0)
#' kish_allocation(n = 3096, N_h = N_h, I = 0)
#'
#' # Proportional allocation (I -> Inf)
#' kish_allocation(n = 3096, N_h = N_h, I = 1e6)
#'
#' # Kish recommended compromise (I = 0.5)
#' kish_allocation(n = 3096, N_h = N_h, I = 0.5)

kish_allocation <- function(n, N_h, I = 0.5) {
  
  if (!is.numeric(n) || length(n) != 1L || n <= 0 || n != round(n))
    stop("`n` must be a single positive integer.", call. = FALSE)
  
  if (!is.numeric(N_h) || any(N_h <= 0))
    stop("`N_h` must be a numeric vector of positive stratum sizes.", call. = FALSE)
  
  if (!is.numeric(I) || length(I) != 1L || I < 0)
    stop("`I` must be a single non-negative number.", call. = FALSE)
  
  if (n > sum(N_h))
    stop("`n` cannot exceed the total population size sum(N_h).", call. = FALSE)
  
  H   <- length(N_h)
  W_h <- N_h / sum(N_h)
  
  num <- sqrt(I * W_h^2 + 1 / H^2)
  n_h <- round(n * num / sum(num))
  
  if (!is.null(names(N_h))) names(n_h) <- names(N_h)
  
  n_h
}