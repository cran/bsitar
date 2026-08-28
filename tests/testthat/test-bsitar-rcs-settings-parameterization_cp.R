
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}


###############################################################################
# Test bsitar with rcs settings  parameterization cp NULL
###############################################################################

test_that("bsitar works with rcs settings and parameterization cp NULL", {
  skip_on_cran()
  
  testthat::skip_if_not(!test_univariate_fit_cov)
  
  
  
  test_fit <- univariate_fit_cov
  
  name_loop <- c('parameterization')
  
  parameterization_list <- list()
  
  parameterization_list [['parameterization']] <- 'cp'
  
  
  update_model_args <- as.list(test_fit$model_info$call.bgmfit)[-1]
  update_model_args[['model']] <- test_fit
  
  # Imp - data must be new
  update_model_args[['newdata']] <- test_fit$model_info$bgmfit.data %>% 
    dplyr::mutate(sexid = sex) %>% 
    dplyr::mutate(logage = log(age))
  
  update_model_args[['save_pars']] <- save_pars(group = TRUE, latent = TRUE, 
                                                all = TRUE, manual = NULL)
  
  
  # Fit each form and collect models
  test_fit_list <- list()
  for (sfml_i in name_loop) {
    
    # print(sfml_i)
    
    update_model_args[['refresh']] <- 0
    update_model_args[['silent']]  <- 2
    
    update_model_args[['parameterization']] <- 
      parameterization_list[[sfml_i]]
    
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
    
    if(sfml_i == 'parameterization') {
      true_sbetas <- 14.11
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(8.9775, 8.9775)
      test_gparms <- get_growthparameters(test_fit, 
                                          re_formula = NA, by = 'sex') %>% 
        dplyr::pull(Estimate) %>% mean() 
    }
    
    
    expect_equal(true_sbetas, test_sbetas, 
                 tolerance = 0.01)
    expect_equal(test_gparms, true_gparms[1], tolerance = 0.01)
    
    
  }
  

})






###############################################################################
# Test bsitar with rcs settings  parameterization cp group
###############################################################################

test_that("bsitar works with rcs settings and parameterization cp group", {
  skip_on_cran()
  
  testthat::skip_if_not(!test_univariate_fit_cov)
  
  test_fit <- univariate_fit_cov
  
  name_loop <- c('parameterization')
  
  parameterization_list <- list()
  
  parameterization_list [['parameterization']] <- 'cp'
  
  
  update_model_args <- as.list(test_fit$model_info$call.bgmfit)[-1]
  update_model_args[['model']] <- test_fit
  
  # Imp - data must be new
  update_model_args[['newdata']] <- test_fit$model_info$bgmfit.data %>% 
    dplyr::mutate(sexid = sex) %>% 
    dplyr::mutate(logage = log(age))
  
  update_model_args[['save_pars']] <- save_pars(group = TRUE, latent = TRUE, 
                                                all = TRUE, manual = NULL)
  
  update_model_args[['a_formula_gr_str']] <- 
    "list(
      (1 | 11 | gr(id, by = sexid))
    )"
  update_model_args[['b_formula_gr_str']] <- 
    "list(
      (1 | 11 | gr(id, by = sexid))
    )"
  update_model_args[['c_formula_gr_str']] <- 
    "list(
      (1 | 11 | gr(id, by = sexid))
    )"
  update_model_args[['d_formula_gr_str']] <- 
    "list(
      (1 | 11 | gr(id, by = sexid))
    )"
  
  
  update_model_args[['sigma_formula']] <- 
    "list(~ 1 + sex * logage)"
  
  update_model_args[['sigma_formula_gr_str']] <- 
    "list(
      (1 + logage| 22 | gr(id, by = sexid))
    )"
  
  # Fit each form and collect models
  test_fit_list <- list()
  for (sfml_i in name_loop) {
    
    # print(sfml_i)
    
    update_model_args[['refresh']] <- 0
    update_model_args[['silent']]  <- 2
    
    update_model_args[['parameterization']] <- 
      parameterization_list[[sfml_i]]
    
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
    
    if(sfml_i == 'parameterization') {
      true_sbetas <- 9.653571
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(9.195, 9.195)
      test_gparms <- get_growthparameters(test_fit, 
                                          re_formula = NA, by = 'sex') %>% 
        dplyr::pull(Estimate) %>% mean() 
    }
    
    
    expect_equal(true_sbetas, test_sbetas, 
                 tolerance = 0.01)
    expect_equal(test_gparms, true_gparms[1], tolerance = 0.01)
    
    
  }
  
  
})



###############################################################################
# Test bsitar with rcs settings  parameterization cp group multi_normal
###############################################################################

test_that("bsitar works with rcs settings and parameterization cp group multi_normal", {
  skip_on_cran()
  
  testthat::skip_if_not(!test_univariate_fit_cov)
  
  test_fit <- univariate_fit_cov
  
  name_loop <- c('parameterization')
  
  parameterization_list <- list()
  
  parameterization_list [['parameterization']] <- 'multi_normal'
  
  
  update_model_args <- as.list(test_fit$model_info$call.bgmfit)[-1]
  update_model_args[['model']] <- test_fit
  
  # Imp - data must be new
  update_model_args[['newdata']] <- test_fit$model_info$bgmfit.data %>% 
    dplyr::mutate(sexid = sex) %>% 
    dplyr::mutate(logage = log(age))
  
  update_model_args[['save_pars']] <- save_pars(group = TRUE, latent = TRUE, 
                                                all = TRUE, manual = NULL)
  
  update_model_args[['a_formula_gr_str']] <- 
    "list(
      (1 | 11 | gr(id, by = sexid))
    )"
  update_model_args[['b_formula_gr_str']] <- 
    "list(
      (1 | 11 | gr(id, by = sexid))
    )"
  update_model_args[['c_formula_gr_str']] <- 
    "list(
      (1 | 11 | gr(id, by = sexid))
    )"
  update_model_args[['d_formula_gr_str']] <- 
    "list(
      (1 | 11 | gr(id, by = sexid))
    )"
  
  
  update_model_args[['sigma_formula']] <- 
    "list(~ 1 + sex * logage)"
  
  update_model_args[['sigma_formula_gr_str']] <- 
    "list(
      (1 + logage| 22 | gr(id, by = sexid))
    )"
  
  # Fit each form and collect models
  test_fit_list <- list()
  for (sfml_i in name_loop) {
    
    # print(sfml_i)
    
    update_model_args[['refresh']] <- 0
    update_model_args[['silent']]  <- 2
    
    update_model_args[['parameterization']] <- 
      parameterization_list[[sfml_i]]
    
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
    
    if(sfml_i == 'parameterization') {
      true_sbetas <- 9.653571
      test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2) %>% mean()
      true_gparms <- c(9.195, 9.195)
      test_gparms <- get_growthparameters(test_fit, 
                                          re_formula = NA, by = 'sex') %>% 
        dplyr::pull(Estimate) %>% mean() 
    }
    
    
    expect_equal(true_sbetas, test_sbetas, 
                 tolerance = 0.01)
    expect_equal(test_gparms, true_gparms[1], tolerance = 0.01)
    
    
  }
  
  
})


