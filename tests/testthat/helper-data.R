

  ##############################################################################
  # Create test data
  # The data will allow testing univariate and multivariate bsitar models
  #   - univariate model
  #   - multivariate model
  ##############################################################################


  ##############################################################################
  # test_data
  ##############################################################################

  test_data <- bsitar::berkeley %>% 
    dplyr::select(id, age, height, weight, sex) %>% 
    dplyr::relocate(height, .before = weight) %>% 
    dplyr::relocate(sex, .before = id) %>% 
    dplyr::filter(age <= 20)
  
  if(nlevels(test_data$sex) != 2) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  if(nlevels(test_data$id) != 136) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  if(round(mean(test_data$age), 2) != 13.36) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  if(round(mean(test_data$height), 2) != 157.56) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  if(round(mean(test_data$weight), 2) != 50.16) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  
  
  ##############################################################################
  # test_data_male & test_data_female
  ##############################################################################
  
  test_data_male   <- test_data %>% dplyr::filter(sex == 'Male') %>% 
    droplevels()

  test_data_female <- test_data %>% dplyr::filter(sex == 'Female') %>% 
    droplevels()
  
  
  ##############################################################################
  # Set global settings
  ##############################################################################
  
  test_univariate_fit_cov    <- TRUE
  test_multivariate_fit_cov  <- FALSE

  save_and_use_models        <- FALSE
  
  # Set skip_test_local_rcmd_check FALSE to run tests on local and GitHub 
  skip_test_local_rcmd_check <- TRUE
 
  
  set.seed(113)
  draw_ids   <- 1:5
  mvar_resp <- 'height'
  uvar_resp <- NULL
  
  set_future_globals_maxSize <- 5.0 * 1e9  # 5.0 GB
  
  ##############################################################################
  # check if model files exists
  ##############################################################################

  testfile_exists <- function(filename) {
    files <- list.files(
      path       = testthat::test_path(),   # root of tests/testthat
      recursive  = TRUE,
      full.names = FALSE
    )
    any(basename(files) == filename)
  }
  
  
  ##############################################################################
  # Wraper for testthat::expect_equal()
  ##############################################################################
  
  informative_expect_equal <- function(object, expected, ...) {
    obj_quo   <- rlang::enquo(object)
    exp_label <- quasi_label(rlang::enquo(expected))
    dots <- list(...)
    if(!is.null(dots)) {
      for (i in names(dots)) {
        assign(i, dots[[i]])
      }
    }
    message("Testing: ", info)
    testthat::expect_equal(object, expected, ...)
    message("✓ ", "Success")
    return(invisible(NULL))
  }
  
  
  ##############################################################################
  # custom skip helper that runs tests only on CI and not on local R CMD Checks
  ##############################################################################
  
  skip_local_run_ci <- function() {
    ci <- Sys.getenv("CI")
    if (identical(ci, "") || identical(ci, "false") || identical(ci, "0")) {
      testthat::skip("Skipping: not on CI")
    }
  }
  
  
  ##############################################################################
  # custom skip helper that runs tests only on CI and not on local R CMD Checks
  ##############################################################################
  
  skip_run_ci <- function() {
    # ci <- Sys.getenv("CI")
    # if (identical(ci, "") || identical(ci, "false") || identical(ci, "0")) {
    #   testthat::skip("Skipping: not on CI")
    # }
    testthat::skip("Skipping: model does not exits")
  }
  
  
  ##############################################################################
  ################################ Model files #################################
  ################### Skip Test Model Fit during load_all = TRUE ###############
  ##############################################################################
  
  # Note that now all tests will be run when skip_test_local_rcmd_check = FALSE
  
  
  ##############################################################################
  ################################ Model files #################################
  ########################## save_and_use_models = TRUE ########################
  ##############################################################################
  
  if(save_and_use_models & !skip_test_local_rcmd_check) {
    ############################################################################
    # test_fit_nsk - with sample_prior = "only" and covariate
    ############################################################################
    
    if(test_univariate_fit_cov & 
       !testfile_exists("univariate_fit_cov.rds")) {
      suppressWarnings(suppressMessages({
        univariate_fit_cov <-
          bsitar::bsitar(x = age, y = height, id = id, data = test_data,
                         df = 4,
                         a_formula = ~ 1 + sex,
                         b_formula = ~ 1 + sex,
                         c_formula = ~ 1 + sex,
                         multivariate = list(mvar = FALSE, rescore = TRUE),
                         univariate_by = FALSE,
                         # stype = list(type = 'nsk', normalize = TRUE),
                         # xoffset = mean,
                         # random = a+b+c,
                         # a_prior_beta = normal(lm, ysd, autoscale = FALSE),
                         # b_prior_beta = normal(0, 1.0, autoscale = FALSE),
                         # c_prior_beta = normal(0, 0.25, autoscale = FALSE),
                         # d_prior_beta = normal(0, 1, autoscale = FALSE),
                         # s_prior_beta = normal(lm, 1, autoscale = 1),
                         # a_cov_prior_beta = normal(0, 10, autoscale = FALSE),
                         # b_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                         # c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                         # d_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                         # s_cov_prior_beta = normal(lm, lm, autoscale = FALSE),
                         # a_prior_sd = normal(0, ysd, autoscale = FALSE),
                         # b_prior_sd = normal(0, 1.0, autoscale = FALSE),
                         # c_prior_sd = normal(0, 0.1, autoscale = FALSE),
                         # d_prior_sd = normal(0, 1, autoscale = FALSE),
                         # rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                         # control = list(adapt_delta = 0.8, max_treedepth = 15),
                         # warmup = NULL,
                         # thin = 5,
                         chains = 1,
                         cores = 1,
                         iter = 2000,
                         threads = NULL,
                         backend = "rstan",
                         sample_prior = "only",
                         init = NULL,
                         vcov_init_0 = TRUE,
                         refres = 0,
                         silent = 2,
                         seed = 123)
        
        saveRDS(univariate_fit_cov, 
                testthat::test_path("models", 
                                    "univariate_fit_cov.rds"), compress = 'xz')
        saveRDS(univariate_fit_cov, 
                testthat::test_path("models", 
                                    "univariate_fit_cov_uncompressed.rds"))
        
      })) # suppressWarnings(suppressMessages({
    } # if(test_univariate_fit_cov) {
    
    
    
    ############################################################################
    # test_fit_nsk - with sample_prior = "only" and covariate
    ############################################################################
    
    if(test_multivariate_fit_cov & 
       !testfile_exists("multivariate_fit_cov.rds")) {
      suppressWarnings(suppressMessages({
        multivariate_fit_cov <-
          bsitar::bsitar(x = age, y = list(height, weight), id = id, 
                         data = test_data,
                         df = 4,
                         a_formula = ~ 1 + sex,
                         b_formula = ~ 1 + sex,
                         c_formula = ~ 1 + sex,
                         multivariate = list(mvar =TRUE, rescore = TRUE),
                         univariate_by = FALSE,
                         # stype = list(type = 'nsk', normalize = TRUE),
                         # xoffset = mean,
                         # random = a+b+c,
                         # a_prior_beta = normal(lm, ysd, autoscale = FALSE),
                         # b_prior_beta = normal(0, 1.0, autoscale = FALSE),
                         # c_prior_beta = normal(0, 0.25, autoscale = FALSE),
                         # d_prior_beta = normal(0, 1, autoscale = FALSE),
                         # s_prior_beta = normal(lm, 1, autoscale = 1),
                         # a_cov_prior_beta = normal(0, 10, autoscale = FALSE),
                         # b_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                         # c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                         # d_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                         # s_cov_prior_beta = normal(lm, lm, autoscale = FALSE),
                         # a_prior_sd = normal(0, ysd, autoscale = FALSE),
                         # b_prior_sd = normal(0, 1.0, autoscale = FALSE),
                         # c_prior_sd = normal(0, 0.1, autoscale = FALSE),
                         # d_prior_sd = normal(0, 1, autoscale = FALSE),
                         # rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                         # control = list(adapt_delta = 0.8, max_treedepth = 15),
                         # warmup = NULL,
                         # thin = 5,
                         chains = 1,
                         cores = 1,
                         iter = 2000,
                         threads = NULL,
                         backend = "rstan",
                         sample_prior = "only",
                         init = NULL,
                         vcov_init_0 = TRUE,
                         refres = 0,
                         silent = 2,
                         seed = 123)
        
        saveRDS(multivariate_fit_cov, 
                testthat::test_path("models", 
                                    "multivariate_fit_cov.rds"), compress = 'xz')
        saveRDS(multivariate_fit_cov, 
                testthat::test_path("models", 
                                    "multivariate_fit_cov_uncompressed.rds"))
        
      })) # # suppressWarnings(suppressMessages({
    } # if(test_multivariate_fit_cov) {
    
  } # if(save_and_use_models & !skip_test_local_rcmd_check) {
  
  


  ##############################################################################
  ################################ Model files #################################
  ########################## save_and_use_models = FALSE ########################
  ##############################################################################
  
  if(!save_and_use_models & !skip_test_local_rcmd_check) {
    ############################################################################
    # test_fit_nsk - with sample_prior = "only" and covariate
    ############################################################################
    
    if(test_univariate_fit_cov) {
      suppressWarnings(suppressMessages({
        univariate_fit_cov <-
          bsitar::bsitar(x = age, y = height, id = id, data = test_data,
                         df = 4,
                         a_formula = ~ 1 + sex,
                         b_formula = ~ 1 + sex,
                         c_formula = ~ 1 + sex,
                         multivariate = list(mvar = FALSE, rescore = TRUE),
                         univariate_by = FALSE,
                         # stype = list(type = 'nsk', normalize = TRUE),
                         # xoffset = mean,
                         # random = a+b+c,
                         # a_prior_beta = normal(lm, ysd, autoscale = FALSE),
                         # b_prior_beta = normal(0, 1.0, autoscale = FALSE),
                         # c_prior_beta = normal(0, 0.25, autoscale = FALSE),
                         # d_prior_beta = normal(0, 1, autoscale = FALSE),
                         # s_prior_beta = normal(lm, 1, autoscale = 1),
                         # a_cov_prior_beta = normal(0, 10, autoscale = FALSE),
                         # b_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                         # c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                         # d_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                         # s_cov_prior_beta = normal(lm, lm, autoscale = FALSE),
                         # a_prior_sd = normal(0, ysd, autoscale = FALSE),
                         # b_prior_sd = normal(0, 1.0, autoscale = FALSE),
                         # c_prior_sd = normal(0, 0.1, autoscale = FALSE),
                         # d_prior_sd = normal(0, 1, autoscale = FALSE),
                         # rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                         # control = list(adapt_delta = 0.8, max_treedepth = 15),
                         # warmup = NULL,
                         # thin = 5,
                         chains = 1,
                         cores = 1,
                         iter = 2000,
                         threads = NULL,
                         backend = "rstan",
                         sample_prior = "only",
                         init = NULL,
                         vcov_init_0 = TRUE,
                         refres = 0,
                         silent = 2,
                         seed = 123)
      
      })) # suppressWarnings(suppressMessages({
    } # if(test_univariate_fit_cov) {
    
    
    
    ############################################################################
    # test_fit_nsk - with sample_prior = "only" and covariate
    ############################################################################
    
    if(test_multivariate_fit_cov) {
      suppressWarnings(suppressMessages({
        multivariate_fit_cov <-
          bsitar::bsitar(x = age, y = list(height, weight), id = id, 
                         data = test_data,
                         df = 4,
                         a_formula = ~ 1 + sex,
                         b_formula = ~ 1 + sex,
                         c_formula = ~ 1 + sex,
                         multivariate = list(mvar =TRUE, rescore = TRUE),
                         univariate_by = FALSE,
                         # stype = list(type = 'nsk', normalize = TRUE),
                         # xoffset = mean,
                         # random = a+b+c,
                         # a_prior_beta = normal(lm, ysd, autoscale = FALSE),
                         # b_prior_beta = normal(0, 1.0, autoscale = FALSE),
                         # c_prior_beta = normal(0, 0.25, autoscale = FALSE),
                         # d_prior_beta = normal(0, 1, autoscale = FALSE),
                         # s_prior_beta = normal(lm, 1, autoscale = 1),
                         # a_cov_prior_beta = normal(0, 10, autoscale = FALSE),
                         # b_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                         # c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                         # d_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                         # s_cov_prior_beta = normal(lm, lm, autoscale = FALSE),
                         # a_prior_sd = normal(0, ysd, autoscale = FALSE),
                         # b_prior_sd = normal(0, 1.0, autoscale = FALSE),
                         # c_prior_sd = normal(0, 0.1, autoscale = FALSE),
                         # d_prior_sd = normal(0, 1, autoscale = FALSE),
                         # rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                         # control = list(adapt_delta = 0.8, max_treedepth = 15),
                         # warmup = NULL,
                         # thin = 5,
                         chains = 1,
                         cores = 1,
                         iter = 2000,
                         threads = NULL,
                         backend = "rstan",
                         sample_prior = "only",
                         init = NULL,
                         vcov_init_0 = TRUE,
                         refres = 0,
                         silent = 2,
                         seed = 123)
        
        
      })) # # suppressWarnings(suppressMessages({
    } # if(test_multivariate_fit_cov) {
    
  } # if(!save_and_use_models & !skip_test_local_rcmd_check) {
  
  
  
  
  set_skip_run_ci <- FALSE
  if(!save_and_use_models) {
    if(test_univariate_fit_cov) {
     if(!exists('univariate_fit_cov')) set_skip_run_ci <- TRUE
    }
    if(test_multivariate_fit_cov) {
      if(!exists('multivariate_fit_cov')) set_skip_run_ci <- TRUE
    }
  }
  
  
  