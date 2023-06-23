

############################################
# Adaptive propensity score truncation
############################################

truncate_pscore_adaptive <- function(A, pi, min_trunc_level = 1e-8) {
  risk_function <- function(cutoff, level) {
    pi <- pmax(pi, cutoff)
    pi <- pmin(pi, 1 - cutoff)
    alpha <- A/pi - (1-A)/(1-pi) #Riesz-representor
    alpha1 <- 1/pi
    alpha0 <- - 1/(1-pi)
    mean(alpha^2 - 2*(alpha1 - alpha0))
  }
  cutoff <- optim(1e-5, fn = risk_function, method = "Brent", lower = min_trunc_level, upper = 0.5, level = 1)$par
  pi <- pmin(pi, 1 - cutoff)
  pi <- pmax(pi, cutoff)
}



############################################
# Custom Riesz loss random forests for inverse propensity score
############################################



compute_inverse_prop <- function(X, A, params) {


  myobjective <- function(preds, dtrain) {
    A <- getinfo(dtrain, "label")
    grad <- 2 * (A* preds - 1)
    hess <-  2 * A
    return(list(grad = grad, hess = hess))
  }

  # Custom Metric
  evalerror <- function(preds, dtrain) {
    A <- getinfo(dtrain, "label")
    err <- A * preds^2  - 2 * preds
    print(mean(err))
    return(list(metric = "MyError", value = mean(err)))
  }




  # Train two models for alpha_1 := 1/P(A=1 | W)

  dtrain_alpha <- xgb.DMatrix(data=as.matrix(X),label=as.matrix(A))

  xgb_alpha <- xgb.train(params = params
                           , data = dtrain_alpha
                           , nrounds = 1, # random forests has 1 round, no boosting.
                           verbose = 0)

  alpha <- predict(xgb_alpha, X) # estimate of alpha_1 := 1/P(A=1 | W)
  return(list(alpha = alpha))
}


library(xgboost)
# L(alpha) = A alpha(W)^2 - 2 alpha(W)
# 2 * alpha * (A - 1)

# Custom objective function (squared error)


n <- 500
X <- as.matrix(runif(n, -1, 1))
m <- ncol(X)
pi <- plogis(3*X)
A <- rbinom(n, 1, pi)

# Set parameters so xgboost is similar to random forest implementation of ranger R package.
params <- list(
  objective = myobjective,
  eval_metric = evalerror,
  learning_rate = 1,
  num_parallel_tree = 500,
  subsample = 0.63,
  colsample_bynode = floor(sqrt(m)) / m,
  reg_lambda = 0,
  max_depth = 10, # max depth should be tuned using cross-validation
  min_child_weight = 2
)

out <- compute_inverse_prop(X, A, params)



############################################
# Calibration for inverse propensity score
############################################





n <- 2000
X <- as.matrix(runif(n, -1, 1))
m <- ncol(X)
pi <- plogis(3*X)
A <- rbinom(n, 1, pi)



# Given an initial estimator of the propensity score, use isotonic regression with Riesz-loss to
# obtain calibrated predictions of inverse propensity score.
calibrate_invpscore <- function(pi, A) {

  params <- list(
    objective = myobjective,
    eval_metric = evalerror,
    learning_rate = 1,
    num_parallel_tree = 1,
    subsample = 1,
    colsample_bynode = floor(sqrt(m)) / m,
    reg_lambda = 0,
    max_depth = 100,
    min_child_weight = 2,
    monotone_constraints = -1
  )
  # Use xgboost to perform isotonic regression with custom loss
  invpscore_1 <- compute_inverse_prop(pi, A, params)$alpha
  invpscore_0 <- compute_inverse_prop(1-pi, 1-A, params)$alpha
  return(list(invpscore_1 = invpscore_1, invpscore_0 = invpscore_0))
}

plot(1/pi, calibrate_invpscore(pi,A)$invpscore_1)
plot(1/(1-pi), calibrate_invpscore(pi,A)$invpscore_0)

