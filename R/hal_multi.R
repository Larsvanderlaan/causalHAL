




n <- 2000
W <- replicate(5,runif(n, -1, 1))
pi <-  plogis(2*rowMeans(sin(5*W)))
A <- rbinom(n, 1, pi)
Y <- rnorm(n, rowMeans(sin(5*W))+ (A - pi)* rowMeans(sin(5*W)), 0.5 )



hal_fit <- fit_hal(as.matrix(W), cbind(A,Y), formula = formula, family = "mgaussian", max_degree = max_degree, smoothness_orders = smoothness_orders, num_knots = num_knots, ...)





hal_multitask <- function(W, A, Y, formula = NULL, max_degree =3, smoothness_orders =1, num_knots = c(), ...) {
  n <- 2000
  W <- replicate(5,runif(n, -1, 1))
  pi <-  plogis(2*rowMeans(sin(5*W)))
  A <- rbinom(n, 1, pi)
  Y <- rnorm(n, rowMeans(sin(5*W))+ (A - pi)* rowMeans(sin(5*W)), 0.5 )
  hal_fit <- fit_hal(as.matrix(W), cbind(A,Y), formula = formula, family = "mgaussian", max_degree = max_degree, smoothness_orders = smoothness_orders, num_knots = num_knots, ...)

}
