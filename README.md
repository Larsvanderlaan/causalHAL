# causalHAL: Adaptive Debiased Machine Learning with HAL

This (in-development) package implements adaptive debiased machine learning estimators for the ATE in data-driven linear and partially linear regression models using the highly adaptive lasso. The theory for these methods is provided in the working paper: `https://arxiv.org/abs/2307.12544`.

`vignette.Rmd` contains example code for running the partially linear ADMLE of the ATE using the highly adaptive lasso (HAL) or lasso (via glmnet).

The R, sh, and sbatch scripts used to run the simulations in the paper can be found in the folder `simulationScripts`. 
Note, at this point in time, the code documentation is fairly poor.

# Motivation and framework 

Debiased machine learning estimators for nonparametric inference of smooth functionals of the data-generating distribution can suffer from excessive variability and instability. For this reason, practitioners may resort to simpler models based on parametric or semiparametric assumptions. However, such simplifying assumptions may fail to hold, and estimates may then be biased due to model misspecification. To address this problem, we propose Adaptive Debiased Machine Learning (ADML), a nonparametric framework that combines data-driven model selection and debiased machine learning techniques to construct asymptotically linear, adaptive, and superefficient estimators for pathwise differentiable functionals. 

By learning model structure directly from data, ADML avoids the bias introduced by model misspecification and remains free from the restrictions of parametric and semiparametric models. While they may exhibit irregular behavior for the target parameter in a nonparametric statistical model, we demonstrate that ADML estimators provides regular and locally uniformly valid inference for a projection-based oracle parameter. Importantly, this oracle parameter agrees with the original target parameter for distributions within an unknown but correctly specified oracle statistical submodel that is learned from the data. This finding implies that there is no penalty, in a local asymptotic sense, for conducting data-driven model selection compared to having prior knowledge of the oracle submodel and oracle parameter.  
