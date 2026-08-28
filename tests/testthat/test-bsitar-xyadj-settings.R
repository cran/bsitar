

# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}
 

###############################################################################
# Test bsitar with xyadj settings
###############################################################################

test_that("bsitar works with genquant_xyadj settings", {
  skip_on_cran()
  
  test_scode <- bsitar(x = age, y = height, id = id, 
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
                       get_stancode = TRUE,
                       get_standata = FALSE, 
                       genquant_xyadj = TRUE,
                       chains = 1, cores = 1, iter = 10, 
                       backend = "rstan",  
                       sample_prior = "no",
                       threads = threading(NULL),
                       # init = '0',
                       init = NULL, # Don't use default random with init_r = 0.5
                       vcov_init_0 = TRUE,
                       refresh = 0, silent = 2,
                       seed = 123)
  
  expect_type(test_scode, "character")
  
  
  test_sdata <- bsitar(x = age, y = height, id = id, 
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
                       get_stancode = FALSE,
                       get_standata = TRUE, 
                       genquant_xyadj = TRUE,
                       chains = 1, cores = 1, iter = 10, 
                       backend = "rstan",  
                       sample_prior = "no",
                       threads = threading(NULL),
                       # init = '0',
                       init = NULL, # Don't use default random with init_r = 0.5
                       vcov_init_0 = TRUE,
                       refresh = 0, silent = 2,
                       seed = 123)
  
  expect_type(test_sdata, "list")
  
  
  suppressWarnings(suppressMessages({
    test_fit <- bsitar(x = age, y = height, id = id,
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
                       get_stancode = FALSE,
                       get_standata = FALSE,
                       genquant_xyadj = TRUE,
                       chains = 1, cores = 1, iter = 10,
                       backend = "rstan",
                       sample_prior = "no",
                       threads = threading(NULL),
                       # init = '0',
                       init = NULL, # Don't use default random with init_r = 0.5
                       vcov_init_0 = TRUE,
                       refresh = 0, silent = 2,
                       seed = 123)
  }))
  
  # Rstan developmental and CRAN versions give different result, so GitHub error
  
  # tomeanx_true <- mean(brms::posterior_summary(test_fit, variable = 'tomeanx_true')[,1])
  # tomeanx_false <-mean(brms::posterior_summary(test_fit, variable = 'tomeanx_false')[,1])
  # tomeany_true <-mean(brms::posterior_summary(test_fit, variable = 'tomeany_true')[,1])
  # tomeany_false <-mean(brms::posterior_summary(test_fit, variable = 'tomeany_false')[,1])
  # 
  # true_tomeanx_true  <- 13.53201
  # true_tomeanx_false <- 13.5299
  # true_tomeany_true  <- 160.0596
  # true_tomeany_false <- 160.0585
  # 
  # 
  # expect_equal(tomeanx_true, true_tomeanx_true, tolerance = 0.01)
  # 
  # expect_equal(tomeanx_false, true_tomeanx_false, tolerance = 0.01)
  # 
  # expect_equal(tomeany_true, true_tomeany_true, tolerance = 0.01)
  # 
  # expect_equal(tomeany_false, true_tomeany_false, tolerance = 0.01)

})
