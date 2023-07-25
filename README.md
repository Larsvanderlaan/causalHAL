# causalHAL: Adaptive Debiased Machine Learning with HAL

Implements adaptive debiased machine learning estimators for the ATE in adaptive linear and partially linear regression models using the highly adaptive lasso.


# Motivation and framework 

Debiased machine learning estimators for nonparametric inference of smooth functionals of the data-generating distribution can suffer from excessive variability and instability. For this reason, practitioners may resort to simpler models based on parametric or semiparametric assumptions. However, such simplifying assumptions may fail to hold, and estimates may then be biased due to model misspecification. To address this problem, we propose Adaptive Debiased Machine Learning (ADML), a nonparametric framework that combines data-driven model selection and debiased machine learning techniques to construct asymptotically linear, adaptive, and superefficient estimators for pathwise differentiable functionals. 

Working paper: https://arxiv.org/abs/2307.12544
