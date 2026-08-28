
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
   skip_local_run_ci()
}

if(set_skip_run_ci) {
  skip_run_ci()
}


###############################################################################
# Test prior_conflict
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
  
  ps <- prior_sensitivity(fit, variable = "b_a_Intercept", 
                          return_table = FALSE, return_conflict = FALSE)
  
  testthat::test_that("prior_conflict validates class data.frame", {
    out <- prior_conflict(ps, 
                             return_table = TRUE)
    
    expect_equal(TRUE, inherits(out, "data.frame"))
  })

  
  testthat::test_that("prior_conflict validates class flextable", {
    out <- prior_conflict(ps, 
                             return_table = TRUE, flex_table = TRUE)
    
    expect_equal(TRUE, inherits(out, "flextable"))
  })
  
  
  
  testthat::test_that("prior_conflict validates class list", {
    out <- prior_conflict(ps, 
                             return_table = FALSE, return_conflict = TRUE)
    
    expect_equal(TRUE, is.list(out)
    )
  })
  
  testthat::test_that("prior_conflict validates class data.frame", {
    out <- prior_conflict(ps, 
                             return_table = TRUE, flex_table = FALSE, return_conflict = TRUE)
    
    expect_equal(TRUE, inherits(out, "data.frame"))
  })
  
  
  testthat::test_that("prior_conflict validates class flextable", {
    out <- prior_conflict(ps, 
                             return_table = TRUE, flex_table = TRUE, return_conflict = TRUE)
    
    expect_equal(TRUE, inherits(out, "flextable"))
  })
  
  
}) 



