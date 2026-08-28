
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
   skip_local_run_ci()
}

if(set_skip_run_ci) {
  skip_run_ci()
}


###############################################################################
# Test prior_sensitivity
###############################################################################

test_that("test-hypothesis_test", {
  skip_on_cran()

  ##############################################################################
  # set model -> save_and_use_models = TRUE
  ##############################################################################
  
  if(save_and_use_models) {
    if(test_univariate_fit_cov) {
      fit               = readRDS(testthat::test_path("models", 
                                                      "univariate_fit_cov.rds")) 
      resp              = uvar_resp
    } else if(test_multivariate_fit_cov) {
      fit               = readRDS(testthat::test_path("models", 
                                                      "multivariate_fit_cov.rds"))  
      resp              = mvar_resp
    } else {
      skip(message = 
             "Both test_univariate_fit_cov and test_multivariate_fit_cov FALSE")
    }
  } # if(save_and_use_models) {
  
  
  ##############################################################################
  # set model -> save_and_use_models = FALSE
  ##############################################################################
  
  if(!save_and_use_models) {
    if(test_univariate_fit_cov) {
      fit               = univariate_fit_cov
      resp              = uvar_resp
    } else if(test_multivariate_fit_cov) {
      fit               = multivariate_fit_cov
      resp              = mvar_resp
    } else {
      skip(message = 
             "Both test_univariate_fit_cov and test_multivariate_fit_cov FALSE")
    }
  } # if(!save_and_use_models) {
  
  
  
  # Use berkeley_exfit
  fit <- berkeley_exfit
  
  testthat::test_that("prior_sensitivity validates plot argument", {
    fake_fit <- structure(list(), class = "bgmfit")
    
    expect_error(
      prior_sensitivity(fake_fit, plot = "bad-plot-type"),
      "argument is of length zero"
    )
  })
  
  
  
  
  testthat::test_that("prior_conflict returns expected names", {
    x <- list(
      sensitivity = data.frame(
        variable = c("a", "a", "b", "b"),
        component = c("prior", "likelihood", "prior", "likelihood"),
        prior = c(0.10, 0.10, 0.10, 0.01),
        likelihood = c(0.10, 0.10, 0.10, 0.01),
        sensitivity = c(0.10, 0.10, 0.10, 0.01)
      )
    )
    class(x) <- c("prior_sensitivity", "bgmfit")
    
    out <- prior_conflict(x, return_table = FALSE)
    
    testthat::expect_named(
      out,
      c(
        "prior_flagged",
        "likelihood_flagged",
        "prior_likelihood_flagged"
      )
    )
    testthat::expect_true("a" %in% out$prior_flagged)
    testthat::expect_true("b" %in% out$likelihood_flagged)
  })
  
  
  
  
  testthat::test_that("prior_sensitivity validates class list", {
    out <- prior_sensitivity(model = fit, variable = "b_a_Intercept", 
                             return_table = FALSE)

    expect_equal(TRUE, is.list(out)
    )
  })
  
  testthat::test_that("prior_sensitivity validates class data.frame", {
    out <- prior_sensitivity(model = fit, variable = "b_a_Intercept", 
                             return_table = TRUE)
    
    expect_equal(TRUE, inherits(out, "data.frame"))
  })

  
  testthat::test_that("prior_sensitivity validates class flextable", {
    out <- prior_sensitivity(model = fit, variable = "b_a_Intercept", 
                             return_table = TRUE, flex_table = TRUE)
    
    expect_equal(TRUE, inherits(out, "flextable"))
  })
  
  
  
  testthat::test_that("prior_sensitivity validates class list", {
    out <- prior_sensitivity(model = fit, variable = "b_a_Intercept", 
                             return_table = FALSE, return_conflict = TRUE)
    
    expect_equal(TRUE, is.list(out)
    )
  })
  
  testthat::test_that("prior_sensitivity validates class data.frame", {
    out <- prior_sensitivity(model = fit, variable = "b_a_Intercept", 
                             return_table = TRUE, flex_table = FALSE, return_conflict = TRUE)
    
    expect_equal(TRUE, inherits(out, "data.frame"))
  })
  
  
  testthat::test_that("prior_sensitivity validates class flextable", {
    out <- prior_sensitivity(model = fit, variable = "b_a_Intercept", 
                             return_table = TRUE, flex_table = TRUE, return_conflict = TRUE)
    
    expect_equal(TRUE, inherits(out, "flextable"))
  })
  
  
}) 



