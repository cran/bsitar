
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}


###############################################################################
# Test bsitar with rcs settings
###############################################################################

test_that("bsitar works with rcs settings and vf functions", {
  skip_on_cran()
  
  testthat::skip_if_not(!test_univariate_fit_cov)
  
  test_fit <- univariate_fit_cov
  
  name_loop <- c('ar', 'arma', 'unstr')
  
  autocor_formula_list <- list()
  
  autocor_formula_list[['ar']] <- "~ar(p = 1)"
  
  autocor_formula_list[['arma']] <- "~arma(p = 1, q = 1)"
  
  autocor_formula_list[['unstr']] <-  "~unstr(age, id)"
  

  update_model_args <- as.list(test_fit$model_info$call.bgmfit)[-1]
  update_model_args[['model']] <- test_fit
  
  # Imp - data must be new
  update_model_args[['newdata']] <- test_fit$model_info$bgmfit.data
  
  update_model_args[['save_pars']] <- save_pars(group = TRUE, latent = TRUE, 
                                                all = TRUE, manual = NULL)
  
  # Fit each form and collect models
  test_fit_list <- list()
  for (sfml_i in name_loop) {
    
    print(sfml_i)
    
    update_model_args[['refresh']] <- 0
    update_model_args[['silent']]  <- 2

    update_model_args[['autocor_formula']] <- 
      autocor_formula_list[[sfml_i]]
    
    suppressWarnings(suppressMessages({
      test_fit_list[[sfml_i]] <- do.call(update_model, update_model_args)
    }))
    
  }
  
  
  
  # Loop over each model to test
  for (sfml_i in name_loop) {
    
   # print(sfml_i)
    
    test_fit <- test_fit_list[[sfml_i]]
    
  #   round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean() %>% print()
  #   get_growthparameters(test_fit, re_formula = NA, by = 'sex') %>% 
  #     dplyr::pull(Estimate) %>% mean() %>% print()
  # 
  # }
    
    if(sfml_i == 'ar') {
      true_sbetas <- 14.659
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(10.4725, 10.4725)
      test_gparms <- get_growthparameters(test_fit, 
                                          re_formula = NA, by = 'sex') %>% 
        dplyr::pull(Estimate) %>% mean() 
    }
    if(sfml_i == 'arma') {
      true_sbetas <- 14.156
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(9.0875, 9.0875)
      test_gparms <- get_growthparameters(test_fit, 
                                          re_formula = NA, by = 'sex') %>% 
        dplyr::pull(Estimate) %>% mean()
    }
    if(sfml_i == 'unstr') {
      true_sbetas <- 14.26
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(9.05 , 9.05)
      test_gparms <- get_growthparameters(test_fit, 
                                          re_formula = NA, by = 'sex') %>% 
        dplyr::pull(Estimate) %>% mean()
    }
    
    expect_equal(true_sbetas, test_sbetas, 
                 tolerance = 0.01)
    expect_equal(test_gparms, true_gparms[1], tolerance = 0.01)

    
  }
 
 
  
})



