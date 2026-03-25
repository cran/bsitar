
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
    skip_local_run_ci()
}

if(set_skip_run_ci) {
  skip_run_ci()
}


###############################################################################
# Test marginals vs marginaleffects
###############################################################################

test_that("test-hypothesis_test", {
  skip_on_cran()
  
  testthat::skip("Skipping: ubuntu fails on CI. object 'drawid' not found, at bsitar/R/utils-helper-28.R:108:3")
  
  ##############################################################################
  # set options
  ##############################################################################
  
  # For current settings (model and draw_ids) min size should be at least 980MB
  oopts <- options(future.globals.maxSize = set_future_globals_maxSize) 
  on.exit(options(oopts), add = TRUE)
  
  estimate_center   <- "mean"
  estimate_interval <- "eti"
  ec_ <- getOption("marginaleffects_posterior_center")
  options("marginaleffects_posterior_center" = estimate_center)
  on.exit(options("marginaleffects_posterior_center" = ec_), add = TRUE)
  
  ei_ <- getOption("marginaleffects_posterior_interval")
  options("marginaleffects_posterior_interval" = estimate_interval)
  on.exit(options("marginaleffects_posterior_interval" = ei_), add = TRUE)
  
  ec_agg <- getOption("marginaleffects_posterior_center")
  ei_agg <- getOption("marginaleffects_posterior_interval")
  
  if(is.null(ec_agg)) ec_agg <- "mean"
  if(is.null(ei_agg)) ei_agg <- "eti"
  
  
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
  
  
  
  ##############################################################################
  # set options
  ##############################################################################
  
  test_tolerance <- 0.01
  
  # Need to re-assign functions to this local test environment
  fit <- bsitar::expose_model_functions(fit, expose = F)
  
  # brms 
  draw_ids          = draw_ids
  ndraws            = NULL
  re_formula        = NA
  dpar              = NULL
  nlpar             = NULL
  incl_autocor      = TRUE
  allow_new_levels  = FALSE
  sample_new_levels = "uncertainty"
  
  ##############################################################################
  # set options
  ##############################################################################
  
  # Key arguments -> by, cov, deriv, variables, condition, comparison, use_d1
  xvar              = 'age'
  cov               = NULL
  deriv             = 1 
  variables         = NULL
  by                = FALSE
  byfun             = NULL
  comparison        = NULL # deriv = 0 = 'difference'; deriv > 0 = 'dydx'
  
  # Keep it NULL / FALSE to match marginaleffects
  use_d1            = NULL # TRUE use _d1, FALSE USES d0 - dydx, NULL == FALSE
  
  ##############################################################################
  # set options
  ##############################################################################
  
  model             <- fit
  newdata           <- bsitar:::get.newdata(model)
  
  # marginal future
  set_future_method  <- 'future' 
  set_future_session <- 'sequential' 
  set_future_splits  <- TRUE
  set_future_cores   <- 2
  
  brms_args_names <- c('draw_ids', 'ndraws', 're_formula', 
                       'resp', 'dpar', 'nlpar', 'incl_autocor',
                       'allow_new_levels', 'sample_new_levels')
  
  
  print_console     <- FALSE # set it TRUE outside testthat 
  
  test_str_cat <- "hypothesis_test vs package: "
  
  ##############################################################################
  # Clear results
  ##############################################################################
  
  
  remove_out_names <- 
    c('brms_model_out', 'hypothesis_test_brms_model_out',
      'marginal_model_out', 'hypothesis_test_marginal_model_out',
      'marginal_model_eqpd_out', 'hypothesis_test_marginal_model_eqpd_out',
      'hypothesis_test_marginal_model_viadraws_out', 
      'hypothesis_test_marginal_model_direct_out',
      'hypothesis_test_marginal_model_viadraws_eqpd_out', 
      'hypothesis_test_marginal_model_direct_eqpd_out',
      'bayestestR_model_eq_out', 'hypothesis_test_bayestestR_model_eq_out',
      'bayestestR_model_pd_out', 'hypothesis_test_bayestestR_model_pd_out')
  
  for (i in remove_out_names) {
    if(exists(i))      rm(i)
  }
  
  
  ##############################################################################
  # hypothesis_test - brms:::hypothesis.brmsfit
  ##############################################################################
  
  ## Example 1: Hypothesis testing (brms-style)
  
  set_hypothesis <- c("a_Intercept = 0", "b_Intercept > 1")
  
  # brms reference
  brms_model_out <- 
    brms::hypothesis(model, hypothesis = set_hypothesis, draw_ids = draw_ids)
  
  # hypothesis_test (auto-detects brms engine)
  hypothesis_test_brms_model_out <- 
    hypothesis_test(model, draw_ids = draw_ids, hypothesis_str = set_hypothesis)
  
  all.equal(brms_model_out$hypothesis,hypothesis_test_brms_model_out$hypothesis)
  
  
  ##############################################################################
  # hypothesis_test - marginaleffects::hypotheses
  ##############################################################################
  
  set_parameter <- c("apgv", "pgv")
  set_range <- list(apgv = 1, pgv = c(0.5, 0.5))
  set_null <- list(apgv = 0, pgv = 0)
  
  # Directly using get_growthparameters()
  marginal_model_out <- 
    get_growthparameters(model, parameter = set_parameter,
                              re_formula = NA,
                              draw_ids = draw_ids)
  # Uncomment for ROPE/p-direction tests:
  # equivalence_test = list(range = set_range),
  # p_direction = list(null = set_null)
  
  hypothesis_test_marginal_model_out <- 
    hypothesis_test(model, parameter = set_parameter, draw_ids = draw_ids)
  
  all.equal(marginal_model_out, hypothesis_test_marginal_model_out)
  
  
  
  ##############################################################################
  # hypothesis_test - marginaleffects::hypotheses with eq pd
  ##############################################################################
  
  set_parameter <- c("apgv", "pgv")
  set_range <- list(apgv = 1, pgv = c(0.5, 0.5))
  set_null <- list(apgv = 0, pgv = 0)
  
  # Directly using get_growthparameters()
  marginal_model_eqpd_out <- 
    get_growthparameters(model, parameter = set_parameter,
                              re_formula = NA,
                              equivalence_test = list(range = set_range),
                              p_direction = list(null = set_null),
                              draw_ids = draw_ids)
  # Uncomment for ROPE/p-direction tests:
  # equivalence_test = list(range = set_range),
  # p_direction = list(null = set_null)
  
  hypothesis_test_marginal_model_eqpd_out <- 
    hypothesis_test(model, parameter = set_parameter,
                    range = set_range,
                    null = set_null,
                    # equivalence_test = list(range = set_range),
                    # p_direction = list(null = set_null),
                    draw_ids = draw_ids)
  
  all.equal(marginal_model_eqpd_out, hypothesis_test_marginal_model_eqpd_out)
  
  
  ##############################################################################
  # hypothesis_test - marginaleffects::hypotheses - via draws
  ##############################################################################
  
  marginal_model_viadraws_input <- 
    get_growthparameters(model, parameter = set_parameter,
                              pdrawsp = TRUE,
                              equivalence_test = list(range = set_range),
                              p_direction = list(null = set_null),
                              draw_ids = draw_ids)
  
  
  hypothesis_test_marginal_model_viadraws_out <- 
    hypothesis_test(marginal_model_viadraws_input, parameter = set_parameter,
                    # equivalence_test = list(range = set_range),
                    # p_direction = list(null = set_null),
                    reformat = TRUE,
                    draw_ids = draw_ids)
  
  hypothesis_test_marginal_model_direct_out <- 
    hypothesis_test(model, parameter = set_parameter, draw_ids = draw_ids)
  
  all.equal(hypothesis_test_marginal_model_viadraws_out, 
            hypothesis_test_marginal_model_direct_out)
  
  
  ##############################################################################
  # hypothesis_test - marginaleffects::hypotheses - via draws eqpd
  ##############################################################################
  
  marginal_model_viadraws_eqpd_input <- 
    get_growthparameters(model, parameter = set_parameter,
                              pdrawsp = TRUE,
                              equivalence_test = list(range = set_range),
                              p_direction = list(null = set_null),
                              draw_ids = draw_ids)
  
  
  hypothesis_test_marginal_model_viadraws_eqpd_out <- 
    hypothesis_test(marginal_model_viadraws_input, parameter = set_parameter,
                    equivalence_test = list(range = set_range),
                    p_direction = list(null = set_null),
                    reformat = TRUE,
                    draw_ids = draw_ids)
  
  hypothesis_test_marginal_model_direct_eqpd_out <- 
    hypothesis_test(model, parameter = set_parameter,
                    range = set_range,
                    null = set_null,
                    # equivalence_test = list(range = set_range),
                    # p_direction = list(null = set_null),
                    reformat = TRUE,
                    draw_ids = draw_ids)
  
  
  all.equal(hypothesis_test_marginal_model_viadraws_eqpd_out, 
            hypothesis_test_marginal_model_direct_eqpd_out)
  
  
  ##############################################################################
  # hypothesis_test - bayestestR:::equivalence_test.brmsfit - eq
  ##############################################################################
  
  set_parameters <- c("b_a_Intercept", "b_b_Intercept")
  set_effects <- 'fixed'
  set_range <- list(b_a_Intercept = c(100, 150), b_b_Intercept = c(-2, 2))
  
  # bayestestR reference - If you installed bayestestR package
  bayestestR_model_eq_out <- 
    bayestestR::equivalence_test(model, parameters = set_parameters,
                                 effects = set_effects,
                                 range = set_range,
                                 draw_ids = draw_ids)
  
  bayestestR_model_eq_out <- bayestestR_model_eq_out %>% 
    dplyr::filter(Parameter %in% set_parameters) %>% 
    dplyr::filter(Effects %in% set_effects)
  
  # hypothesis_test (auto-detects bayestestR engine)
  # If see package is installed (install.packages("see")), then you can plot 
  # the results by setting plot = TRUE
  hypothesis_test_bayestestR_model_eq_out <- 
    hypothesis_test(model, parameters = set_parameters,
                    range = set_range,
                    plot = FALSE,
                    draw_ids = draw_ids)
  
  all.equal(bayestestR_model_eq_out, hypothesis_test_bayestestR_model_eq_out)
  

  ##############################################################################
  # hypothesis_test - bayestestR:::equivalence_test.brmsfit - pd
  ##############################################################################
  
  set_parameters <- c("b_a_Intercept", "b_b_Intercept")
  set_effects <- 'fixed'
  set_range <- list(b_a_Intercept = c(100, 150), b_b_Intercept = c(-2, 2))
  set_null <- list(b_a_Intercept = 0, b_b_Intercept = 0)
  
  # bayestestR reference - If you installed bayestestR package
  bayestestR_model_pd_out <- 
    bayestestR::p_direction(model, parameters = set_parameters,
                            effects = set_effects,
                            range = set_range,
                            draw_ids = draw_ids)
  
  bayestestR_model_pd_out <- bayestestR_model_pd_out %>% 
    dplyr::filter(Parameter %in% set_parameters) %>% 
    dplyr::filter(Effects %in% set_effects)
  
  # hypothesis_test (auto-detects bayestestR engine)
  # If see package is installed (install.packages("see")), then you can plot 
  # the results by setting plot = TRUE
  hypothesis_test_bayestestR_model_pd_out <- 
    hypothesis_test(model, parameters = set_parameters,
                    null = set_null,
                    equivalence_test = FALSE,
                    p_direction = TRUE,
                    plot = FALSE,
                    draw_ids = draw_ids)
  
  all.equal(bayestestR_model_pd_out, hypothesis_test_bayestestR_model_pd_out)
  
  ##############################################################################
  # run all tests
  ##############################################################################
  
  paste0(test_str_cat, " ", " brms")
  expect_equal(brms_model_out,   
               hypothesis_test_brms_model_out,   tolerance = 0.01)
  
  paste0(test_str_cat, " ", " marginal")
  expect_equal(marginal_model_out,   
               hypothesis_test_marginal_model_out,   tolerance = 0.01)
  
  paste0(test_str_cat, " ", " eqpd")
  expect_equal(marginal_model_eqpd_out,   
               hypothesis_test_marginal_model_eqpd_out,   tolerance = 0.01)
  
  paste0(test_str_cat, " ", " viadraws")
  expect_equal(hypothesis_test_marginal_model_viadraws_out,  
               hypothesis_test_marginal_model_direct_out,   tolerance = 0.01)
  
  paste0(test_str_cat, " ", " viadraws_eqpd")
  expect_equal(hypothesis_test_marginal_model_viadraws_eqpd_out,   
               hypothesis_test_marginal_model_direct_eqpd_out, tolerance = 0.01)
  
  paste0(test_str_cat, " ", " bayestestR_eq")
  expect_equal(bayestestR_model_eq_out,   
               hypothesis_test_bayestestR_model_eq_out,   tolerance = 0.01)
  
  paste0(test_str_cat, " ", " bayestestR_pd")
  expect_equal(bayestestR_model_pd_out,   
               hypothesis_test_bayestestR_model_pd_out,   tolerance = 0.01)
  
  
  
}) # test_that("test-marginals-slopes-bycov-byvariable", {






