# bsitar 0.2.0


### New feature

- Added ```add_model_criterion()``` function to compute fit criteria such as "loo", "waic", "kfold", 
 "loo_subsample", "bayes_R2" (Bayesian R-squared), "loo_R2" (LOO-adjusted R-squared), and "marglik" 
 (log marginal likelihood). The computed fit criteria are added to the model object for later use including 
 comparison of models.  The ```add_model_criterion()``` is a wrapper around the ```add_criterion()``` 
 function available from the  **brms** package


### Bugfixes

- bsitar(): The 'sigma_cov_init_beta = random' argument was setting wrong initial values for the covariates  
  included in the 'sigma' formula. The initials for Intercept ('sigma_init_beta') were used for covariates too.



### Miscellaneous
- Replaced Bayesian SITAR model fit shown as example from 'berkeley_mfit' applied to 20 randomly selected males   between 6 and 20 years of age ('berkeley_mdata') with 'berkeley_exfit' that is fit to 70 females between 8 and 20 years of age ('berkeley_exdata'). This is done to use the example model ('berkeley_exfit') in vignette that provide a detailed comparison between non Bayesian SITAR model fit (using the 'sitar' package) and Bayesian SITAR model fit (using the 'bsitar' package). The vignette included in the 'sitar' package analysed the exact same data (70 females between 8 and 20 years of age). 
- The ``bsitar::bsitar()```received options ```file```, ``file_refit```, and ``file_compress``` to save and retreive fitted objects. See  ```brms::brm``` help file for details. 
- Minor corrections/changes to make R code more efficient.


# bsitar 0.1.1


---

### New feature

- Added 'optimize_model' function to perform model optimization by fitting 
  model with varying degree of freedom and by transforming 'x' (predictor) 
  and 'y' (outcome) varibales. The allowed transformations for 'x' and 'y' 
  variables are 'log' (logarithmic) and 'sqrt' (square root) transformation. 
  The 'optimize_model' performs comparison of resulting model fits based on 
  user specified  criteria such as the Watanabe–Akaike information criterion
  ('waic'), leave-one-out cross-validation ('loo') and the Bayesian R square 
  ('bayes_R2'). Please see help file of 'optimize_model' to see documentation.

### Minor changes

- Updated documentation
- Added vignette about the Bayesian SITAR model fit to height data (univariate model)
- Added issue tracker url: https://github.com/Sandhu-SS/bsitar/issues


### Bugfixes

- plot_curves(): Fixed bug to remove warning "Duplicated aesthetics after name 
standardisation: group in when plotting together unadjusted and adjusted curves." 
- plot_curves(): Fixed bug to remove palettes error when tried plotting all 
four curves together (distance, velocity, adjusted and unadjusted)
- plot_curves(): and growthparameters() fixed issues relating to not using 
'dplyr::all_of()' when within the 'dplyr::select()'


### Miscellaneous
- Added more utility functions for internal use (for consistency and efficiency) 
- Minor corrections/changes to make R code more efficient and consistent across sub modules.



