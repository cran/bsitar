
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}


###############################################################################
# Test bsitar with rcs QR settings
###############################################################################

test_that("bsitar works with rcs settings and decomp QR", {
  skip_on_cran()
 
  test_scode <- bsitar(x = age, y = height, id = id, 
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
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
                       # parameterization = "cp",
                       decomp = "QR",
                       seed = 123)
  
  expect_type(test_scode, "character")
  
  
  test_sdata <- bsitar(x = age, y = height, id = id, 
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
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
                       # parameterization = "cp",
                       decomp = "QR",
                       seed = 123)
  
  expect_type(test_sdata, "list")
  
  
  suppressWarnings(suppressMessages({
    test_fit <- bsitar(x = age, y = height, id = id,
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
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
                       # parameterization = "cp",
                       decomp = "QR",
                       seed = 123)
  }))
  
  
  # Rstan developmental and CRAN versions give different result, so GitHub error
  
  # true_sbetas <- c(158.13,  -0.01,   0.02,  15.06,   0.65,  -3.38,  -0.21)
  # 
  # test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2)
  # 
  # expect_equal(true_sbetas, test_sbetas, tolerance = 0.01)
  # 
  # test_gparms <- get_growthparameters(test_fit, re_formula = NA)
  # 
  # expect_equal(round(test_gparms$Estimate[1], 2), 12.72, tolerance = 0.01)
  # expect_equal(round(test_gparms$Estimate[2], 2), 6.55,  tolerance = 0.01)

})



