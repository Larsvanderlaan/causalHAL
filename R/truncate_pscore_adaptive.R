

#' @export
#'
truncate_pscore_adaptive <- function(A, pi, min_trunc_level = 1e-8) {
  cutoffs <- seq(0.1, min_trunc_level, length = 500)
  risks <- sapply(cutoffs, function(cutoff) {
    pi <- pmin(pi, 1 - cutoff)
    pi <- pmax(pi, cutoff)
    alpha <- A/pi - (1-A)/(1-pi) #Riesz-representor
    alpha1 <- 1/pi
    alpha0 <- - 1/(1-pi)
    mean(alpha^2 - 2*(alpha1 - alpha0))
  })
  cutoff <- cutoffs[which.min(risks)[1]]
  pi <- pmin(pi, 1 - cutoff)
  pi <- pmax(pi, cutoff)
}
