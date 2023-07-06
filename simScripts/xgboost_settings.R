
library(xgboost)
objective_invprop <- function(preds, dtrain) {
  A <- xgboost::getinfo(dtrain, "label")
  grad <- 2 * (A* preds - 1)
  hess <-  2 * A
  return(list(grad = grad, hess = hess))
}

evalerror_invprop <- function(preds, dtrain) {
  A <- xgboost::getinfo(dtrain, "label")
  err <- A * preds^2  - 2 * preds
  return(list(metric = "MyError", value = mean(err)))
}


make_glmnet_inv_prop <- function(nrounds,  alpha  =1, type ) {
  if(type == "pi") {
    objective <- objective_invprop
    eval_metric <- evalerror_invprop
  } else {
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
  }
  return(Lrnr_xgboost$new(nrounds=nrounds,
                          booster="gblinear",
                          updater = "coord_descent",
                          feature_selector = "cyclic",
                          objective = objective,
                          eval_metric = eval_metric,
                          lambda = 0,
                          alpha = alpha
  ))

}







make_ranger_inv_prop <- function(max_depth, type) {
  if(type == "pi") {
    objective <- objective_invprop
    eval_metric <- evalerror_invprop
  } else {
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
  }

  return(Lrnr_xgboost$new(
    nrounds=1,
    objective = objective,
    eval_metric = eval_metric,
    eta = 1,
    num_parallel_tree = 500,
    subsample = 0.63,
    colsample_bynode = 1,
    lambda = 0,
    max_depth = max_depth,
    min_child_weight = 2
  ))

}




make_xgboost_inv_prop <- function(nrounds, max_depth, type, eta = 0.2 ) {
  if(type == "pi") {
    objective <- objective_invprop
    eval_metric <- evalerror_invprop
  } else {
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
  }
  return(Lrnr_xgboost$new(nrounds=nrounds,
                   max_depth = max_depth,
                   objective = objective,
                   eval_metric = eval_metric,
                   eta = eta
  ))

}







make_misspecified <- function( type) {

  if(type == "pi") {
    objective <- objective_invprop
    eval_metric <- evalerror_invprop
  } else {
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
  }
  return(Lrnr_xgboost$new(nrounds=1,
                   max_depth = 2,
                   eta = 1,
                   objective = objective,
                   eval_metric = eval_metric
                   ))

}

get_learners <- function() {


  lrnr_invpi_xg <- Stack$new(list(
    make_ranger_inv_prop(max_depth = 4, type = "pi"),
    make_ranger_inv_prop(max_depth = 5,type = "pi"),
    make_ranger_inv_prop(max_depth = 6,  type = "pi"),
    make_ranger_inv_prop(max_depth = 7, type = "pi")))
   # make_xgboost_inv_prop(nrounds = 20, max_depth = 3, type = "pi"),
  #  make_xgboost_inv_prop(nrounds = 20, max_depth = 4, type = "pi"),
  #  make_xgboost_inv_prop(nrounds = 20, max_depth = 5, type = "pi"),
  #  make_xgboost_inv_prop(nrounds = 20, max_depth = 6, type = "pi")))
  lrnr_invpi_xg <- Pipeline$new(Lrnr_cv$new(lrnr_invpi_xg), Lrnr_cv_selector$new(loss_squared_error))




  lrnr_mu_xg <- Stack$new(list(
    make_ranger_inv_prop(max_depth = 4, type = "mu"),
    make_ranger_inv_prop(max_depth = 5, type = "mu"),
    make_ranger_inv_prop(max_depth = 6, type = "mu"),
    make_ranger_inv_prop(max_depth = 7, type = "mu")))
    #make_xgboost_inv_prop(nrounds = 20, max_depth = 3, type = "mu"),
    #make_xgboost_inv_prop(nrounds = 20, max_depth = 4, type = "mu"),
    #make_xgboost_inv_prop(nrounds = 20, max_depth = 5, type = "mu"),
    #make_xgboost_inv_prop(nrounds = 20, max_depth = 6, type = "mu")))
  lrnr_mu_xg <- Pipeline$new(Lrnr_cv$new(lrnr_mu_xg), Lrnr_cv_selector$new(loss_squared_error))
  #lrnr_mu_xg <- Lrnr_cv$new(Lrnr_ranger$new())
  list(pi = lrnr_invpi_xg, mu = lrnr_mu_xg, misp_pi = Lrnr_cv$new(make_misspecified("pi")),
       misp_mu = Lrnr_cv$new(make_misspecified("mu")))
}



Stack$new(Lrnr_gam$new(), Lrnr_earth$new(), Lrnr_hal9001$new())
# crossfit estimator
lrnr_list <- get_learners()
lrnr_mu <- lrnr_list$mu
lrnr_pi <- lrnr_list$pi
lrnr_misp_pi <- lrnr_list$misp_pi
lrnr_misp_mu <- lrnr_list$misp_mu
