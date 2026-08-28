
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
  
  sigma_formula_manual_loop <- c('ba', 'fp', 'me', 're', 
                                 'cp', 've', 'vp', 'ls')
  
  sigma_formula_manual_list <- list()
  
  sigma_formula_manual_list[['ba']] <- "list(nlf(sigma ~ z, method = 'ba') + 
         lf(z ~ 1 + splines2::nsk(age_vf) + (1 | gr(id))))"
  
  sigma_formula_manual_list[['fp']] <- "nlf(sigma ~ vf(param1, 
                                                             param2, identity()), 
                                                  method = 'fp') +
                         lf(param1 + param2 ~ 1)"
  
  sigma_formula_manual_list[['me']] <-  "nlf(sigma ~ vf(param1, param2, 
                                                              identity()), 
                                                   method = 'me') +
                         lf(param1 + param2 ~ 1)"
  
  sigma_formula_manual_list[['re']] <-  "nlf(sigma ~ vf(param1, 
                                                              param2, 
                                                              identity(), 
                                                              outcome), 
                                                   method = 're') +
                         lf(param1 + param2 ~ 1)"
  
  sigma_formula_manual_list[['cp']] <-  "nlf(sigma ~ vf(param1, param2, 
                                                             param3, age_vf),
                                                  method = 'cp') +
                         lf(param1 + param2 + param3 ~ 1)"
  
  sigma_formula_manual_list[['ve']] <-  "nlf(sigma ~ vf(param1, param2, 
                                                             age_vf), 
                                                  method = 've') +
                         lf(param1 + param2 ~ 1)"
  
  sigma_formula_manual_list[['vp']] <-  "nlf(sigma ~ vf(param1, 
                                                             param2, age_vf),
                                                  method = 'vp') +
                         lf(param1 + param2 ~ 1)"
  
  
  sigma_formula_manual_list[['ls']] <-  "nlf(sigma ~ ls(x, sigmaa, sigmab, sigmac, 
                                        sigmas1, sigmas2, sigmas3, sigmas4),
                             loop = FALSE) +
                         lf(sigmaa ~ 1+(1 |110| gr(id, by = NULL))) + 
                         lf(sigmab ~ 1+(1 |110| gr(id, by = NULL))) + 
                         lf(sigmac ~ 1+(1 |110| gr(id, by = NULL))) + 
                         lf(sigmas1 + sigmas2 + sigmas3 + sigmas4 ~ 1)"
 
 
  update_model_args <- as.list(test_fit$model_info$call.bgmfit)[-1]
  update_model_args[['model']] <- test_fit
  update_model_args[['newdata']] <- test_fit$model_info$bgmfit.data %>% 
    dplyr::mutate(age_vf = age)
  
  
  
  # Fit each form and collect models
  test_fit_list <- list()
  for (sfml_i in sigma_formula_manual_loop) {
    
   # print(sfml_i)
    
    update_model_args[['refresh']] <- 0
    update_model_args[['silent']]  <- 2
    if(sfml_i == 'ls') {
      update_model_args[['sigmax']]     <- 'age_vf'
      update_model_args[['sigmafixed']] <- 'a+b+c'
    }
    # update_model_args[['sigma_formula_manual']] <- 
    #   gsub_space(sigma_formula_manual_list[[sfml_i]])
    
    update_model_args[['sigma_formula_manual']] <- 
      sigma_formula_manual_list[[sfml_i]]
    
    suppressWarnings(suppressMessages({
      test_fit_list[[sfml_i]] <- do.call(update_model, update_model_args)
    }))
    
  }
  
  
  
  # Loop over each model to test
  for (sfml_i in sigma_formula_manual_loop) {
  #  print(sfml_i)
    test_fit <- test_fit_list[[sfml_i]]
    
  #   round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean() %>% print()
  #   get_growthparameters(test_fit, re_formula = NA) %>% print()
  #   
  # }
    
    if(sfml_i == 'ba') {
      true_sbetas <- 12.0125
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(11.89, 6.62)
      test_gparms <- get_growthparameters(test_fit, re_formula = NA)
    }
    if(sfml_i == 'fp') {
      true_sbetas <- 11.91667
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(11.57, 6.42)
      test_gparms <- get_growthparameters(test_fit, re_formula = NA)
    }
    if(sfml_i == 'me') {
      true_sbetas <- 11.88083
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(11.55 , 6.42  )
      test_gparms <- get_growthparameters(test_fit, re_formula = NA)
    }
    if(sfml_i == 're') {
      true_sbetas <- 11.93167
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(11.57, 6.41)
      test_gparms <- get_growthparameters(test_fit, re_formula = NA)
    }
    if(sfml_i == 'cp') {
      true_sbetas <- 10.98308
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(11.92, 6.72)
      test_gparms <- get_growthparameters(test_fit, re_formula = NA)
    }
    if(sfml_i == 've') {
      true_sbetas <- 11.87833
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(11.59, 6.42 )
      test_gparms <- get_growthparameters(test_fit, re_formula = NA)
    }
    if(sfml_i == 'vp') {
      true_sbetas <- 11.98
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(11.66, 6.37)
      test_gparms <- get_growthparameters(test_fit, re_formula = NA)
    }
    if(sfml_i == 'ls') {
      true_sbetas <- 8.104118
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(11.11, 6.78)
      test_gparms <- get_growthparameters(test_fit, re_formula = NA)
    }
 
    expect_equal(true_sbetas, test_sbetas, tolerance = 0.01)
    expect_equal(round(test_gparms$Estimate[1], 2), true_gparms[1], tolerance = 0.01)
    expect_equal(round(test_gparms$Estimate[2], 2), true_gparms[2],  tolerance = 0.01)
    
  }
 
 
  
})



