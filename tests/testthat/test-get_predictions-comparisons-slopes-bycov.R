
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

test_that("test-get_predictions-comparisons-slopes-bycov", {
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
  
  marginal_deriv_args <- 1
  
  ###############################################################################
  # marginaleffects::slopes vs get_predictions - average = FALSE
  ###############################################################################
  
  variables = xvar # note variables set as xvar for slope
  vcov = TRUE
  conf_level = 0.95
  type = NULL
  by = cov
  byfun = NULL
  wts = FALSE
  transform = NULL
  hypothesis = NULL
  equivalence = NULL
  df = Inf
  numderiv = "fdforward"
  
  # marginaleffects::slopes
  marginaleffects_slopes_args <- list()
  marginaleffects_slopes_args[['model']] <- model
  marginaleffects_slopes_args[['newdata']] <- newdata
  marginaleffects_slopes_args[['variables']] <- variables
  marginaleffects_slopes_args[['vcov']] <- vcov
  marginaleffects_slopes_args[['conf_level']] <- conf_level
  marginaleffects_slopes_args[['type']] <- type
  marginaleffects_slopes_args[['by']] <- by
  marginaleffects_slopes_args[['byfun']] <- byfun
  marginaleffects_slopes_args[['wts']] <- wts
  marginaleffects_slopes_args[['transform']] <- transform
  marginaleffects_slopes_args[['hypothesis']] <- hypothesis
  marginaleffects_slopes_args[['equivalence']] <- equivalence
  marginaleffects_slopes_args[['df']] <- df
  marginaleffects_slopes_args[['numderiv']] <- numderiv
  #
  marginaleffects_slopes_args[['draw_ids']] <- draw_ids
  marginaleffects_slopes_args[['ndraws']] <- ndraws
  marginaleffects_slopes_args[['re_formula']] <- re_formula
  
  
  marginaleffects_out <- do.call(marginaleffects::slopes, 
                                 marginaleffects_slopes_args) %>%
    data.frame()
  
  
  get_predictions_args <- marginaleffects_slopes_args
  
  get_predictions_args[['method']]        <- 'pkg'
  get_predictions_args[['reformat']]      <- FALSE
  get_predictions_args[['average']]       <- FALSE
  get_predictions_args[['newdata_fixed']] <- 0
  get_predictions_args[['deriv']]         <- marginal_deriv_args
  
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
  # marginaleffects::avg_slopes vs get_predictions - average = TRUE
  ###############################################################################
  # check if any new argument needed for avg_slopes than slopes
  setdiff(methods::formalArgs(marginaleffects::slopes),
          methods::formalArgs(marginaleffects::avg_slopes))
  
  
  marginaleffects_avg_slopes_args <- marginaleffects_slopes_args
  
  marginaleffects_out <- do.call(marginaleffects::slopes, 
                                 marginaleffects_avg_slopes_args) %>% 
    data.frame()
  
  
  get_predictions_args <- marginaleffects_avg_slopes_args
  
  get_predictions_args[['method']]        <- 'pkg'
  get_predictions_args[['reformat']]      <- FALSE
  get_predictions_args[['average']]       <- TRUE
  get_predictions_args[['newdata_fixed']] <- 0
  get_predictions_args[['deriv']]         <- marginal_deriv_args
  
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
  # marginaleffects::plot_slopes vs get_predictions - average = FALSE
  # with by 
  ###############################################################################
  
  # check if any new argument needed for avg_slopes than slopes
  args_remove <-
    setdiff(methods::formalArgs(marginaleffects::slopes),
            methods::formalArgs(marginaleffects::plot_slopes))
  
  args_add <-
    setdiff(methods::formalArgs(marginaleffects::plot_slopes),
            methods::formalArgs(marginaleffects::slopes))
  
  
  marginaleffects_plot_slopes_args <- marginaleffects_slopes_args
  for (i in args_remove) {
    marginaleffects_plot_slopes_args[[i]] <- NULL
  }
  
  # args_add
  marginaleffects_plot_slopes_args[['condition']] <- NULL
  marginaleffects_plot_slopes_args[['points']] <- NULL # 0
  marginaleffects_plot_slopes_args[['rug']] <- FALSE
  marginaleffects_plot_slopes_args[['gray']] <- FALSE
  marginaleffects_plot_slopes_args[['draw']] <- TRUE
  # need to change by
  marginaleffects_plot_slopes_args[['by']] <- c(xvar, cov)
  marginaleffects_plot_slopes_args[['condition']] <- NULL
  
  
  marginaleffects_plot <- do.call(marginaleffects::plot_slopes, 
                                  marginaleffects_plot_slopes_args)
  
  marginaleffects_plot_slopes_args[['draw']] <- FALSE
  marginaleffects_out <- do.call(marginaleffects::plot_slopes, 
                                 marginaleffects_plot_slopes_args) %>% 
    data.frame()
  
  
  get_predictions_args <- marginaleffects_plot_slopes_args
  
  get_predictions_args[['method']]        <- 'pkg'
  get_predictions_args[['reformat']]      <- FALSE
  get_predictions_args[['average']]       <- TRUE
  get_predictions_args[['newdata_fixed']] <- 0
  get_predictions_args[['plot']]          <- TRUE
  get_predictions_args[['deriv']]         <- marginal_deriv_args
  
  
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
  
  out_mdT <- round(marginal_out_custom_mdT$estimate, 2)
  out_mdF <- round(marginal_out_custom_mdF$estimate, 2)
  
  # out_1 <- round(marginaleffects_out$estimate, 2)
  # out_2 <- round(marginal_out$estimate, 2)
  # out_3 <- out_mdT
  
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
  # marginaleffects::plot_slopes vs get_predictions - average = FALSE
  # with condition 
  ###############################################################################
  
  # check if any new argument needed for avg_slopes than slopes
  args_remove <-
    setdiff(methods::formalArgs(marginaleffects::slopes),
            methods::formalArgs(marginaleffects::plot_slopes))
  
  args_add <-
    setdiff(methods::formalArgs(marginaleffects::plot_slopes),
            methods::formalArgs(marginaleffects::slopes))
  
  
  marginaleffects_plot_slopes_args <- marginaleffects_slopes_args
  for (i in args_remove) {
    marginaleffects_plot_slopes_args[[i]] <- NULL
  }
  
  # args_add
  marginaleffects_plot_slopes_args[['condition']] <- NULL
  marginaleffects_plot_slopes_args[['points']] <- NULL # 0
  marginaleffects_plot_slopes_args[['rug']] <- FALSE
  marginaleffects_plot_slopes_args[['gray']] <- FALSE
  marginaleffects_plot_slopes_args[['draw']] <- TRUE
  # need to change by
  marginaleffects_plot_slopes_args[['by']] <- NULL
  marginaleffects_plot_slopes_args[['condition']] <- c(xvar, cov)
  
  marginaleffects_plot <- do.call(marginaleffects::plot_slopes, 
                                  marginaleffects_plot_slopes_args)
  
  marginaleffects_plot_slopes_args[['draw']] <- FALSE
  marginaleffects_out <- do.call(marginaleffects::plot_slopes, 
                                 marginaleffects_plot_slopes_args) %>% 
    data.frame()
  
  get_predictions_args <- marginaleffects_plot_slopes_args
  
  get_predictions_args[['method']]        <- 'pkg'
  get_predictions_args[['reformat']]      <- FALSE
  get_predictions_args[['average']]       <- TRUE
  get_predictions_args[['newdata_fixed']] <- 0
  get_predictions_args[['plot']]          <- TRUE
  get_predictions_args[['deriv']]         <- marginal_deriv_args
  
  
  
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
  
  get_predictions_args_custom_plot <- get_predictions_args_custom
  get_predictions_args_custom_plot[['plot']] <- TRUE
  do.call(get_predictions, 
          get_predictions_args_custom_plot)
  
  out_mdT <- round(marginal_out_custom_mdT$estimate, 2)
  out_mdF <- round(marginal_out_custom_mdF$estimate, 2)
  
  # out_1 <- round(marginaleffects_out$estimate, 2)
  # out_2 <- round(marginal_out$estimate, 2)
  # out_3 <- out_mdT
  
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
