# Forster_Warmuth_counterfactual_regression

This repo contains the code and data for the simulation and real data application results in [this paper](https://arxiv.org/abs/2307.16798).

The file to source is `toSource.R`, this contains all the main functions used in the experiment. In particular, the function `series_df()` can be used to implement either the Forster--Warmuth estimator, or the least squares estimator. It support basis function `bs`, `ns`, which stand for B-splines for polynomial or natural splines or `poly` (orthogonal polynomial basis). `series_cv_new()` is used to automatically choose the number of basis functions.


## Synthetic data.
1. save_plotdata_uniform.R generates the data where the covariate X follows a uniform distribution, and uses series_df function from `toSourse.R` to generate estimates using various methods including dry, Forster¡ª=Warmth and least squares estimators with different basis. It saves the output into `plotdata_uniform_2000.RData`.
2. `save_plotdata_mixture.R` generates the data where the covariate X follows a heavy-tailed distribution, and uses `series_df` function from `toSourse.R `to generate estimates and gives output `plotdata_mixture_2000.RData`, which is for sample size 2000. `plotdata_mixture_400.RData` is likewise generated for sample size n=400.
3. Finally, `ggplot_uniform.R` takes in `plotdata_uniform_2000.RData` and produces output `ggplot_uniform.pdf`; And ggplot_mixture.R takes in both `plotdata_mixture_2000.RData` and `plotdata_mixture_400.RData`, and gives output `ggplot_mixture.pdf`.

Real data.
`data_prep.R` cleans up variables in `rhc.csv`, which is publicly available at [https://hbiostat.org/data/](https://hbiostat.org/data/).
1. `save_data_pseudo_unconfound.R` produces the data including the pseudo outcomes in the unconfoundedness setting. All 71 variables are used to form the covariate X.  The output is in `pseudo_unconfound.RData`. 
1.1 `save_plotdata_unconfound_cv_bs_poly.R` takes in `pseudo_unconfound.RData`,
And give mean and std for the bs and polynomial basis using cross_validation, the main function is series_cv_new, which is in `toSource.R`. It outputs the plotting data `unconfound_plotdata.RData`.
2. `save_data_pseudo_proxy.R` produces the data including the pseudo outcomes in the proximal setting. The pseudo-outcomes are generated from CF_cate function in `toSourse.R`. The output is in `pseudo_proxy.RData`.
2.1. `save_plotdata_proxy_cv_bs_poly.R` takes in `pseudo_proxy.RData`,
And give mean and std for the bs and polynomial basis using cross_validation, the main function is `series_cv_new`, which is in `toSource.R`. It outputs the plotting data `proxy_plotdata.RData`.
3. Finally, `ggplot_unconfound_proxy.R` takes in `unconfound_plotdata.RData` and `proxy_plotdata.RData`. This produces `ggplot_poly_bs.pdf`.

