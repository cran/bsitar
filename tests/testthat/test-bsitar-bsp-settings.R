

# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}
 

###############################################################################
# Test bsitar with bsp settings
###############################################################################

test_that("bsitar works fully with bsp settings", {
  skip_on_cran()
  
  test_scode <- bsitar(x = age, y = height, id = id, 
                       data = test_data_male,  df = 4,
                       stype = list(type = 'bsp', normalize = TRUE),
                       get_stancode = TRUE,
                       get_standata = FALSE, 
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
                       stype = list(type = 'bsp', normalize = TRUE),
                       get_stancode = FALSE,
                       get_standata = TRUE, 
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
                       stype = list(type = 'bsp', normalize = TRUE),
                       get_stancode = FALSE,
                       get_standata = FALSE,
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
  

  true_sbetas <- c(129.12, 0.00, -0.01, 6.45, 38.75, 57.64, 45.17)
  
  test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2)
  
  expect_equal(true_sbetas, test_sbetas, tolerance = 0.01)
  
  test_gparms <- get_growthparameters(test_fit, re_formula = NA)
  
  expect_equal(round(test_gparms$Estimate[1], 2), 11.98, tolerance = 0.01)
  expect_equal(round(test_gparms$Estimate[2], 2), 6.35, tolerance = 0.01)

})
