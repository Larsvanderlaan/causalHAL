# causalHAL: Adaptive Debiased Machine Learning with HAL

This (in-development) package implements adaptive debiased machine learning estimators for the ATE in data-driven linear and partially linear regression models using the highly adaptive lasso. The theory for these methods is provided in the working paper: `https://arxiv.org/abs/2307.12544`.

`vignette.Rmd` contains example code for running the partially linear ADMLE of the ATE using the highly adaptive lasso (HAL) or lasso (via glmnet).

The R, sh, and sbatch scripts used to run the simulations in the paper can be found in the folder `simulationScripts`. 
Note, at this point in time, the code documentation is fairly poor.

# Motivation and framework 

Debiased machine learning estimators for nonparametric inference of smooth functionals of the data-generating distribution can suffer from excessive variability and instability. For this reason, practitioners may resort to simpler models based on parametric or semiparametric assumptions. However, such simplifying assumptions may fail to hold, and estimates may then be biased due to model misspecification. To address this problem, we propose Adaptive Debiased Machine Learning (ADML), a nonparametric framework that combines data-driven model selection and debiased machine learning techniques to construct asymptotically linear, adaptive, and superefficient estimators for pathwise differentiable functionals. 

By learning model structure directly from data, ADML avoids the bias introduced by model misspecification and remains free from the restrictions of parametric and semiparametric models. While they may exhibit irregular behavior for the target parameter in a nonparametric statistical model, we demonstrate that ADML estimators provides regular and locally uniformly valid inference for a projection-based oracle parameter. Importantly, this oracle parameter agrees with the original target parameter for distributions within an unknown but correctly specified oracle statistical submodel that is learned from the data. This finding implies that there is no penalty, in a local asymptotic sense, for conducting data-driven model selection compared to having prior knowledge of the oracle submodel and oracle parameter.  

 


## Install

```{r}
devtools::install_github("tlverse/hal9001")
devtools::install_github("Larsvanderlaan/causalHAL")
```

## Run ADMLE using HAL and glmnet

```{r}
library(causalHAL)
library(hal9001)
seed <- rnorm(1)
n <- 1000
d <- 4
pos_const <- 1
W <- replicate(d, runif(n, -1, 1))
colnames(W) <- paste0("W", 1:d)
pi0 <- plogis(pos_const * ( W[,1] + sin(4*W[,1]) +   W[,2] + cos(4*W[,2]) + W[,3] + sin(4*W[,3]) + W[,4] + cos(4*W[,4]) ))
A <- rbinom(n, 1, pi0)
mu0 <-  sin(4*W[,1]) + sin(4*W[,2]) + sin(4*W[,3])+  sin(4*W[,4]) + cos(4*W[,2])
tau <- 1 + W[,1] + abs(W[,2]) + cos(4*W[,3]) + W[,4]
Y <- rnorm(n,  mu0 + A * tau, 0.5)


# User-supplied estimate of propensity score pi = P(A=1|W)
pi.hat <- pi0
# User-supplied estimate of treatment-marginalized outcome regression m = E(Y|W)
m.hat <- mu0 * pi.hat + (mu0 + tau) * (1-pi.hat)

# ADMLE for ATE using partially linear model with HAL.
# Fits additive piece-wise linear spline model for CATE with 50 knot points per covariate using highly adaptive lasso (see tlverse/hal9001 github R package)
set.seed(seed)
ADMLE_fit <- fit_cate_hal_partially_linear(W, A, Y, 
                                           m.hat = m.hat,
                                           pi.hat = pi.hat,
                                           smoothness_orders_cate = 1, num_knots_cate = c(50), max_degree_cate = 1)
# Provides estimates and CI for ATE
inference_ate(ADMLE_fit)

# Same analysis but using glmnet implementation with hal9001-basis design matrix.
# May not reproduce estimates exactly but should be close.
# For those not familiar with hal9001 package, the below code may be easier to play around with.
basis_list <- hal9001::enumerate_basis(W, smoothness_orders = 1, num_knots = 50, max_degree = 1)
tau_basis <- hal9001::make_design_matrix(W, basis_list)
set.seed(seed)
ADMLE_fit <- fit_cate_lasso_partially_linear(tau_basis, A, Y, 
                                             m.hat = m.hat,
                                             pi.hat = pi.hat, standardize = FALSE)

# Provides estimates and CI for ATE
inference_ate(ADMLE_fit)
 ```

