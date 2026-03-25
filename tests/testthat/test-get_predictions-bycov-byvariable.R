
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

test_that("test-get_predictions-bycov-byvariable", {
  skip_on_cran()
  
  # skip_if(!exists('fit_cov'), message = "Model 'fit_cov' not found")
  
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
  
  model = fit
  draw_ids = draw_ids
  ndraws = NULL
  re_formula = NA
  xvar <- 'age'
  cov  <- 'sex'
  newdata <- bsitar:::get.newdata(model)
  
  ###############################################################################
  # marginaleffects::predictions vs get_predictions - average = FALSE
  ###############################################################################
  
  variables = cov
  vcov = TRUE
  conf_level = 0.95
  type = NULL
  by = c(xvar, cov)
  byfun = NULL
  wts = FALSE
  transform = NULL
  hypothesis = NULL
  equivalence = NULL
  df = Inf
  numderiv = "fdforward"
  
  # marginaleffects::predictions
  marginaleffects_predictions_args <- list()
  marginaleffects_predictions_args[['model']] <- model
  marginaleffects_predictions_args[['newdata']] <- newdata
  marginaleffects_predictions_args[['variables']] <- variables
  marginaleffects_predictions_args[['vcov']] <- vcov
  marginaleffects_predictions_args[['conf_level']] <- conf_level
  marginaleffects_predictions_args[['type']] <- type
  marginaleffects_predictions_args[['by']] <- by
  marginaleffects_predictions_args[['byfun']] <- byfun
  marginaleffects_predictions_args[['wts']] <- wts
  marginaleffects_predictions_args[['transform']] <- transform
  marginaleffects_predictions_args[['hypothesis']] <- hypothesis
  marginaleffects_predictions_args[['equivalence']] <- equivalence
  marginaleffects_predictions_args[['df']] <- df
  marginaleffects_predictions_args[['numderiv']] <- numderiv
  #
  marginaleffects_predictions_args[['draw_ids']] <- draw_ids
  marginaleffects_predictions_args[['ndraws']] <- ndraws
  marginaleffects_predictions_args[['re_formula']] <- re_formula
  
  
  marginaleffects_out <- do.call(marginaleffects::predictions, 
                                 marginaleffects_predictions_args) %>%
    data.frame()
  
  
  get_predictions_args <- marginaleffects_predictions_args
  
  get_predictions_args[['method']]        <- 'pkg'
  get_predictions_args[['reformat']]      <- FALSE
  get_predictions_args[['average']]       <- FALSE
  get_predictions_args[['newdata_fixed']] <- 0
  
  
  marginal_out <- do.call(get_predictions, 
                          get_predictions_args) %>% data.frame()
  
  
  get_predictions_args_custom                  <- get_predictions_args
  get_predictions_args_custom[['method']]      <- 'custom'
  get_predictions_args_custom[['model_deriv']] <- T
  marginal_out_custom_mdT <- do.call(get_predictions, 
                                     get_predictions_args_custom) %>% data.frame()
  
  
  get_predictions_args_custom[['model_deriv']] <- F
  marginal_out_custom_mdF <- do.call(get_predictions, 
                                     get_predictions_args_custom) %>% data.frame()
  
  out_mdT <- round(marginal_out_custom_mdT$estimate, 2)
  out_mdF <- round(marginal_out_custom_mdF$estimate, 2)
  
  
  out_1 <- round(marginaleffects_out$estimate, 2)
  out_2 <- round(marginal_out$estimate, 2)
  out_3 <- out_mdF
  
  
  expect_equal(out_1,   out_2,   tolerance = 0.01)
  expect_equal(out_mdT, out_mdF, tolerance = 0.01)
  expect_equal(out_2,   out_3,   tolerance = 0.01)
  
  
  # if(!identical(out_1, out_2)) {
  #   stop()
  # }
  # 
  # if(!identical(out_mdT, out_mdF)) {
  #   stop()
  # }
  # 
  # if(!identical(out_2, out_3)) {
  #   stop()
  # }
  
  
  ###############################################################################
  # marginaleffects::avg_predictions vs get_predictions - average = TRUE
  ###############################################################################
  # check if any new argument needed for avg_predictions than predictions
  setdiff(methods::formalArgs(marginaleffects::predictions),
          methods::formalArgs(marginaleffects::avg_predictions))
  
  
  marginaleffects_avg_predictions_args <- marginaleffects_predictions_args
  
  marginaleffects_out <- do.call(marginaleffects::predictions, 
                                 marginaleffects_avg_predictions_args) %>% 
    data.frame()
  
  
  get_predictions_args <- marginaleffects_avg_predictions_args
  
  get_predictions_args[['method']]        <- 'pkg'
  get_predictions_args[['reformat']]      <- FALSE
  get_predictions_args[['average']]       <- TRUE
  get_predictions_args[['newdata_fixed']] <- 0
  
  marginal_out <- do.call(get_predictions, 
                          get_predictions_args) %>% data.frame()
  
  
  get_predictions_args_custom                  <- get_predictions_args
  get_predictions_args_custom[['method']]      <- 'custom'
  get_predictions_args_custom[['model_deriv']] <- T
  marginal_out_custom_mdT <- do.call(get_predictions, 
                                     get_predictions_args_custom) %>% data.frame()
  
  
  get_predictions_args_custom[['model_deriv']] <- F
  marginal_out_custom_mdF <- do.call(get_predictions, 
                                     get_predictions_args_custom) %>% data.frame()
  
  out_mdT <- round(marginal_out_custom_mdT$estimate, 2)
  out_mdF <- round(marginal_out_custom_mdF$estimate, 2)
  
  out_1 <- round(marginaleffects_out$estimate, 2)
  out_2 <- round(marginal_out$estimate, 2)
  out_3 <- out_mdF
  
  
  expect_equal(out_1,   out_2,   tolerance = 0.01)
  expect_equal(out_mdT, out_mdF, tolerance = 0.01)
  expect_equal(out_2,   out_3,   tolerance = 0.01)
  
  
  # if(!identical(out_1, out_2)) {
  #   stop()
  # }
  # 
  # if(!identical(out_mdT, out_mdF)) {
  #   stop()
  # }
  # 
  # if(!identical(out_2, out_3)) {
  #   stop()
  # }
  
  
  ###############################################################################
  # marginaleffects::plot_predictions vs get_predictions - average = FALSE
  # with by 
  ###############################################################################
  
  # check if any new argument needed for avg_predictions than predictions
  args_remove <-
    setdiff(methods::formalArgs(marginaleffects::predictions),
            methods::formalArgs(marginaleffects::plot_predictions))
  
  args_add <-
    setdiff(methods::formalArgs(marginaleffects::plot_predictions),
            methods::formalArgs(marginaleffects::predictions))
  
  
  marginaleffects_plot_predictions_args <- marginaleffects_predictions_args
  for (i in args_remove) {
    marginaleffects_plot_predictions_args[[i]] <- NULL
  }
  
  # args_add
  marginaleffects_plot_predictions_args[['condition']] <- NULL
  marginaleffects_plot_predictions_args[['points']] <- 0
  marginaleffects_plot_predictions_args[['rug']] <- FALSE
  marginaleffects_plot_predictions_args[['gray']] <- FALSE
  marginaleffects_plot_predictions_args[['draw']] <- TRUE
  # need to change by
  marginaleffects_plot_predictions_args[['by']] <- by
  marginaleffects_plot_predictions_args[['condition']] <- NULL
  
  
  marginaleffects_plot <- do.call(marginaleffects::plot_predictions, 
                                  marginaleffects_plot_predictions_args)
  
  marginaleffects_plot_predictions_args[['draw']] <- FALSE
  marginaleffects_out <- do.call(marginaleffects::plot_predictions, 
                                 marginaleffects_plot_predictions_args) %>% 
    data.frame()
  
  
  get_predictions_args <- marginaleffects_plot_predictions_args
  
  get_predictions_args[['method']]        <- 'pkg'
  get_predictions_args[['reformat']]      <- FALSE
  get_predictions_args[['average']]       <- TRUE
  get_predictions_args[['newdata_fixed']] <- 0
  get_predictions_args[['plot']]          <- TRUE
  
  marginal_plot <- do.call(get_predictions, 
                           get_predictions_args)
  
  
  get_predictions_args[['plot']]          <- FALSE
  marginal_out <- do.call(get_predictions, 
                          get_predictions_args) %>% data.frame()
  
  
  get_predictions_args_custom                  <- get_predictions_args
  get_predictions_args_custom[['method']]      <- 'custom'
  get_predictions_args_custom[['model_deriv']] <- T
  marginal_out_custom_mdT <- do.call(get_predictions, 
                                     get_predictions_args_custom) %>% data.frame()
  
  
  get_predictions_args_custom[['model_deriv']] <- F
  marginal_out_custom_mdF <- do.call(get_predictions, 
                                     get_predictions_args_custom) %>% data.frame()
  
  
  # Note that for plot, length differs between marginaleffects and marginals
  # But rounded mean should be same
  out_1 <- round(mean(marginaleffects_out$estimate), 1)
  out_2 <- round(mean(marginal_out$estimate), 1)
  out_3 <- round(mean(marginal_out_custom_mdF$estimate), 1)
  
  
  expect_equal(out_1,   out_2,   tolerance = 0.01)
  expect_equal(out_mdT, out_mdF, tolerance = 0.01)
  expect_equal(out_2,   out_3,   tolerance = 0.01)
  
  
  # if(!identical(out_1, out_2)) {
  #   stop()
  # }
  # 
  # if(!identical(out_mdT, out_mdF)) {
  #   stop()
  # }
  # 
  # if(!identical(out_2, out_3)) {
  #   stop()
  # }
  
  
  ###############################################################################
  # marginaleffects::plot_predictions vs get_predictions - average = FALSE
  # with condition 
  ###############################################################################
  
  # check if any new argument needed for avg_predictions than predictions
  args_remove <-
    setdiff(methods::formalArgs(marginaleffects::predictions),
            methods::formalArgs(marginaleffects::plot_predictions))
  
  args_add <-
    setdiff(methods::formalArgs(marginaleffects::plot_predictions),
            methods::formalArgs(marginaleffects::predictions))
  
  
  marginaleffects_plot_predictions_args <- marginaleffects_predictions_args
  for (i in args_remove) {
    marginaleffects_plot_predictions_args[[i]] <- NULL
  }
  
  # args_add
  marginaleffects_plot_predictions_args[['condition']] <- NULL
  marginaleffects_plot_predictions_args[['points']] <- 0
  marginaleffects_plot_predictions_args[['rug']] <- FALSE
  marginaleffects_plot_predictions_args[['gray']] <- FALSE
  marginaleffects_plot_predictions_args[['draw']] <- TRUE
  # need to change by
  marginaleffects_plot_predictions_args[['by']] <- NULL
  marginaleffects_plot_predictions_args[['condition']] <- by
  
  marginaleffects_plot <- do.call(marginaleffects::plot_predictions, 
                                  marginaleffects_plot_predictions_args)
  
  marginaleffects_plot_predictions_args[['draw']] <- FALSE
  marginaleffects_out <- do.call(marginaleffects::plot_predictions, 
                                 marginaleffects_plot_predictions_args) %>% 
    data.frame()
  
  get_predictions_args <- marginaleffects_plot_predictions_args
  
  get_predictions_args[['method']]        <- 'pkg'
  get_predictions_args[['reformat']]      <- FALSE
  get_predictions_args[['average']]       <- TRUE
  get_predictions_args[['newdata_fixed']] <- 0
  get_predictions_args[['plot']]          <- TRUE
  
  
  marginal_plot <- do.call(get_predictions, 
                           get_predictions_args)
  
  get_predictions_args[['plot']]          <- FALSE
  
  
  marginal_out <- do.call(get_predictions, 
                          get_predictions_args) %>% data.frame()
  
  
  get_predictions_args_custom                  <- get_predictions_args
  get_predictions_args_custom[['method']]      <- 'custom'
  get_predictions_args_custom[['model_deriv']] <- T
  marginal_out_custom_mdT <- do.call(get_predictions, 
                                     get_predictions_args_custom) %>% data.frame()
  
  
  get_predictions_args_custom[['model_deriv']] <- F
  marginal_out_custom_mdF <- do.call(get_predictions, 
                                     get_predictions_args_custom) %>% data.frame()
  
  
  # Note that for plot, length differs between marginaleffects and marginals
  # But rounded mean should be same
  out_1 <- round(mean(marginaleffects_out$estimate), 1)
  out_2 <- round(mean(marginal_out$estimate), 1)
  out_3 <- round(mean(marginal_out_custom_mdF$estimate), 1)
  
  
  expect_equal(out_1,   out_2,   tolerance = 0.01)
  expect_equal(out_mdT, out_mdF, tolerance = 0.01)
  expect_equal(out_2,   out_3,   tolerance = 0.01)
  
  
  # if(!identical(out_1, out_2)) {
  #   stop()
  # }
  # 
  # if(!identical(out_mdT, out_mdF)) {
  #   stop()
  # }
  # 
  # if(!identical(out_2, out_3)) {
  #   stop()
  # }
  
  
  
  
  
})
