
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}

if(set_skip_run_ci) {
  skip_run_ci()
}



###############################################################################
# Test marginals vs marginaleffects comparisons
###############################################################################


# if (Sys.getenv("SKIP_API_TESTS") == "TRUE") {
#   skip("API tests skipped via environment variable")
# }
# skip_if(!exists('univariate_fit_cov'), message = "Model 'fit_cov' not found")

# skip_if(Sys.getenv("API_KEY") == "", "API key not available")
# devtools::load_all()


test_that("test-get_predictions-comparisons-parameters-default-d1", {
  skip_on_cran()
  
  
  # devtools::load_all()
  
  # The only difference between 
  # test-marginals-draws-comparisons-parameters-default-d0
  # and
  # test-marginals-draws-comparisons-parameters-default-d1
  # is deriv = 0 vs deriv = 1
  
  
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
  
  slope             = "dydx"
  
  averageTF         = FALSE
  
  # Keep it NULL / FALSE to match marginaleffects
  use_d1            = NULL # TRUE use _d1, FALSE USES d0 - dydx, NULL == FALSE
  
  ##############################################################################
  # set options
  ##############################################################################
  
  model             <- fit
  newdata           <- bsitar:::get.newdata(model)
  
  # marginal future
  set_future_session <- 'sequential' 
  set_future_cores   <- 2
  
  brms_args_names <- c('draw_ids', 'ndraws', 're_formula', 
                       'resp', 'dpar', 'nlpar', 'incl_autocor',
                       'allow_new_levels', 'sample_new_levels')
  
  
  print_console     <- FALSE # set it TRUE outside testthat 
  
  
  ###############################################################################
  ###############################################################################
  ###############################################################################
  
  # marginaleffects::predictions vs get_comparisons
  
  what_test <- 'predictions'
  
  ###############################################################################
  ###############################################################################
  ###############################################################################
  
  ##############################################################################
  # set marginaleffects predictions
  ##############################################################################
  
  set_comparison_for_marginaleffects    <- FALSE
  if(is.null(comparison)) {
    set_comparison_for_marginaleffects <- TRUE
  }
  
  if(what_test != 'comparison') {
    set_comparison_for_marginaleffects    <- FALSE
  }
  
  set_marginaleffects_comparison <- function(deriv, comparison) {
    if(deriv == 0) {
      if(is.null(comparison)) {
        comparison <- 'difference'
      }
    } else if(deriv > 0) {
      if(is.null(comparison)) {
        comparison <- 'dydx'
      }
    }
    return(comparison)
  }
  
  
  
  
  vcov = TRUE
  conf_level = 0.95
  type = NULL
  cross = FALSE
  wts = FALSE
  transform = NULL
  hypothesis = NULL
  equivalence = NULL
  df = Inf
  numderiv = "fdforward"
  eps = NULL
  
  
  # Plot
  condition = NULL
  points = NULL
  rug = FALSE
  gray = FALSE
  draw = TRUE
  
  
  
  ###############################################################################
  ###############################################################################
  # Construct test predictions
  ###############################################################################
  ###############################################################################
  
  out_names_comb_sep <- " vs "
  
  out_names <- 
    c('marginaleffects_out', 'marginal_out', 
      'marginal_out_custom_mdT', 'marginal_out_custom_mdF', 
      'marginal_out_custom_mdT_future', 'marginal_out_custom_mdF_future',
      'marginal_out_custom_mdT_dofuture', 'marginal_out_custom_mdF_dofuture')
  
  out_names_pairs <- apply(combn(out_names,2),2,paste,collapse=out_names_comb_sep)
  
  
  vector1 <- c("marginaleffects_out")
  vector2 <- 
    c("marginal_out", 
      "marginal_out_custom_mdT", "marginal_out_custom_mdF", 
      "marginal_out_custom_mdT_future", "marginal_out_custom_mdF_future",
      "marginal_out_custom_mdT_dofuture", "marginal_out_custom_mdF_dofuture")
  
  out_names_ref <- expand.grid(V1 = vector1, V2 = vector2) %>% 
    tidyr::unite(z, V1, V2, sep = out_names_comb_sep) %>% dplyr::pull(z)
  
  ###############################################################################
  ###############################################################################
  # marginaleffects::predictions vs get_comparisons - average = FALSE
  ###############################################################################
  ###############################################################################
  
  test_str_cat <- 
    "marginaleffects::predictions vs get_comparisons -> average = FALSE"
  
  
  if(deriv == 0) {
    if( averageTF) {
      marginaleffects_funcall      <- marginaleffects::avg_predictions
      marginaleffects_funcall_plot <- marginaleffects::plot_predictions
    }
    if(!averageTF) {
      marginaleffects_funcall      <- marginaleffects::predictions
      marginaleffects_funcall_plot <- marginaleffects::plot_predictions
    }
  } else  if(deriv > 0) {
    if( averageTF) {
      marginaleffects_funcall      <- marginaleffects::avg_slopes
      marginaleffects_funcall_plot <- marginaleffects::plot_slopes
    }
    if(!averageTF) {
      marginaleffects_funcall <- marginaleffects::slopes
      marginaleffects_funcall_plot <- marginaleffects::plot_slopes
    }
  }
  
  marginal_funcall        <- get_predictions
  
  what_test               <- paste0(what_test, "_", "avgF")
  
  if(exists('marginaleffects_args_est')) {
    rm('marginaleffects_args_est')
  }
  if(exists('marginaleffects_args_plot_by')) {
    rm('marginaleffects_args_plot_by')
  }
  if(exists('marginaleffects_args_plot_condition')) {
    rm('marginaleffects_args_plot_condition')
  }
  
  if(exists('marginal_args')) {
    rm('marginal_args')
  }
  
  if(exists('marginaleffects_out'))              rm('marginaleffects_out')
  if(exists('marginal_out'))                     rm('marginal_out')
  if(exists('marginal_out_custom_mdT'))          rm('marginal_out_custom_mdT')
  if(exists('marginal_out_custom_mdF'))          rm('marginal_out_custom_mdF')
  if(exists('marginal_out_custom_mdT_future'))   rm('marginal_out_custom_mdT_future')
  if(exists('marginal_out_custom_mdF_future'))   rm('marginal_out_custom_mdF_future')
  if(exists('marginal_out_custom_mdT_dofuture')) rm('marginal_out_custom_mdT_dofuture')
  if(exists('marginal_out_custom_mdF_dofuture')) rm('marginal_out_custom_mdF_dofuture')
  
  
  ###############################################################################
  # Set up arguments
  ###############################################################################
  
  all_args_names <- NULL
  all_args_names <- formalArgs(marginaleffects_funcall)
  all_args_names <- c(all_args_names, brms_args_names)
  
  marginaleffects_args <- list()
  for (i in all_args_names) {
    if(i != "...") {
      if(!exists(i)) {
        stop("arg ", i, " not defined in the workspace")
      }
      marginaleffects_args[[i]] <- get(i)
    } # if(i != "...") {
  } # for (i in all_args_names) {
  
  
  if(set_comparison_for_marginaleffects) {
    marginaleffects_args[['comparison']] <- NULL
  }
  
  
  ###############################################################################
  # marginaleffects
  ###############################################################################
  
  test_infor <- "marginaleffects,avg=F,deriv=0,model_deriv=F,future=F,plot=F,pby=NULL,pcond=NULL"
  
  show_test_infor <- FALSE
  
  if(set_comparison_for_marginaleffects) {
    marginaleffects_args[['comparison']] <- 
      set_marginaleffects_comparison(deriv=deriv, comparison=comparison)
  }
  
  marginaleffects_args_est <- marginaleffects_args
  marginaleffects_args_plot_by <- marginaleffects_args
  marginaleffects_args_plot_condition <- marginaleffects_args
  
  if(deriv > 0) {
    marginaleffects_args_est[['slope']] <- 'dydx'
    marginaleffects_args_est[['variables']] <- xvar
  }
  
  
  marginaleffects_out <- do.call(marginaleffects_funcall, 
                                 marginaleffects_args_est) %>% data.frame()
  
  # marginaleffects estimates -  compare marginals with marginaleffects
  set_expect_est    <- mean(marginaleffects_out$estimate) # 159.2748
  
  # when plot d > 0, by or condition, then need to specify variables
  
  # marginaleffects estimates via plot by
  marginaleffects_args_plot_by[['draw']]      <- TRUE
  marginaleffects_args_plot_by[['by']]        <- xvar
  marginaleffects_args_plot_by[['condition']] <- NULL
  
  if(deriv > 0) marginaleffects_args_plot_by[['variables']] <- xvar
  
  marginaleffects_plot_by <- CustomDoCall(marginaleffects_funcall_plot, 
                                          marginaleffects_args_plot_by) 
  
  set_expect_plot_by <- NULL
  set_expect_plot_by <- mean(marginaleffects_plot_by$data$estimate) # 161.0859
  
  # marginaleffects estimates via plot condition
  marginaleffects_args_plot_condition[['draw']]      <- TRUE
  marginaleffects_args_plot_condition[['by']]        <- NULL
  marginaleffects_args_plot_condition[['condition']] <- xvar
  
  if(deriv > 0) marginaleffects_args_plot_condition[['variables']] <- xvar
  
  
  marginaleffects_plot_condition <- CustomDoCall(marginaleffects_funcall_plot, 
                                                 marginaleffects_args_plot_condition) 
  
  set_expect_plot_condition <- NULL
  set_expect_plot_condition <- mean(marginaleffects_plot_condition$data$estimate) # 155.8943
  
  
  
  set_tolerance <- 0.01
  if(!show_test_infor) test_infor <- NULL
  
  informative_expect_equal(mean(marginaleffects_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  
  # head(marginaleffects_out)
  
  
  ###############################################################################
  # marginal - method = pkg
  ###############################################################################
  
  test_infor <- "method=pkg,avg=F,deriv=0,model_deriv=F,future=F,plot=F,pby=NULL,pcond=NULL"
  
  marginal_args <- marginaleffects_args_est
  
  marginal_args[['method']]        <- 'pkg'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- F
  marginal_args[['future_method']]  <- 'future'
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- TRUE
  
  if(!show_test_infor) test_infor <- NULL
  
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['use_d1']]        <- T
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  # devtools::load_all()
  # marginal_funcall        <- get_comparisons
  # CustomDoCall(marginal_funcall, 
  #              marginal_args) %>% data.frame()
  
  
  
  
  ###############################################################################
  # marginal - method = custom
  ###############################################################################
  
  test_infor <- "method=custom,avg=F,deriv=0,model_deriv=F,future=F,plot=F,pby=NULL,pcond=NULL"
  
  marginal_args <- marginaleffects_args_est
  
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- F
  marginal_args[['future_method']]  <- 'future'
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- TRUE
  
  if(!show_test_infor) test_infor <- NULL
  
  
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['use_d1']]        <- T
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  ###############################################################################
  # marginal - method = custom
  ###############################################################################
  
  test_infor <- "method=custom,avg=T,deriv=0,model_deriv=F,future=F,plot=F,pby=NULL,pcond=NULL"
  
  marginal_args <- marginaleffects_args_est
  
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- F
  marginal_args[['future_method']]  <- 'future'
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- TRUE
  
  if(!show_test_infor) test_infor <- NULL
  
  
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['use_d1']]        <- T
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  ###############################################################################
  # marginal - method = custom
  ###############################################################################
  
  test_infor <- "method=custom,avg=F,deriv=0,model_deriv=F,future=T,plot=F,pby=NULL,pcond=NULL"
  
  marginal_args <- marginaleffects_args_est
  
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- T
  marginal_args[['future_method']]  <- 'future'
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- TRUE
  
  if(!show_test_infor) test_infor <- NULL
  
  
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['use_d1']]        <- T
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  ###############################################################################
  # marginal - method = custom
  ###############################################################################
  
  test_infor <- "method=custom,avg=T,deriv=0,model_deriv=F,future=T,plot=F,pby=NULL,pcond=NULL"
  
  marginal_args <- marginaleffects_args_est
  
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- T
  marginal_args[['future_method']]  <- 'future'
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- TRUE
  
  if(!show_test_infor) test_infor <- NULL
  
  
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['use_d1']]        <- T
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  
  ###############################################################################
  # marginal - method = custom
  ###############################################################################
  
  test_infor <- "method=custom,avg=F,deriv=0,model_deriv=F,dofuture=T,plot=F,pby=NULL,pcond=NULL"
  
  marginal_args <- marginaleffects_args_est
  
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- T
  marginal_args[['future_method']]  <- 'dofuture'
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- TRUE
  
  if(!show_test_infor) test_infor <- NULL
  
  
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['use_d1']]        <- T
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  ###############################################################################
  # marginal - method = custom - plot by
  ###############################################################################
  
  test_infor <- "method=custom,avg=T,deriv=0,model_deriv=F,dofuture=T,plot=F,pby=NULL,pcond=NULL"
  
  marginal_args <- marginaleffects_args_est
  
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- T
  marginal_args[['future_method']]  <- 'dofuture'
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- TRUE
  
  if(!show_test_infor) test_infor <- NULL
  
  
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['use_d1']]        <- T
  marginal_out <- NULL
  marginal_out <- CustomDoCall(marginal_funcall, 
                               marginal_args) %>% data.frame() # dplyr::pull(estimate) %>% mean
  
  informative_expect_equal(mean(marginal_out$estimate), set_expect_est, 
                           tolerance = set_tolerance, info = test_infor)
  
  ###############################################################################
  # marginal - method = custom - plot by
  ###############################################################################
  
  marginal_args <- marginaleffects_args_plot_by
  
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['plot']]      <- TRUE
  marginal_args[['future']]    <- FALSE
  marginal_args[['by']]        <- xvar
  marginal_args[['condition']] <- NULL
  
  
  # devtools::load_all()
  # marginal_funcall        <- get_predictions
  marginal_plot <- CustomDoCall(marginal_funcall, marginal_args) 
  
  informative_expect_equal(mean( marginal_plot$data$estimate), set_expect_plot_by, 
                           tolerance = set_tolerance, info = test_infor)
  
  
  ###############################################################################
  # marginal - method = custom - plot condition
  ###############################################################################
  
  marginal_args <- marginaleffects_args_plot_condition
  
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- averageTF
  marginal_args[['deriv']]         <- deriv
  marginal_args[['model_deriv']]   <- F
  marginal_args[['use_d1']]        <- F
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['plot']]      <- TRUE
  marginal_args[['future']]    <- FALSE
  marginal_args[['by']]        <- NULL
  marginal_args[['condition']] <- xvar
  

  marginal_plot <- CustomDoCall(marginal_funcall, marginal_args) 
  
  # mean(marginal_plot$data$estimate)
  
  # informative_expect_equal(mean(marginal_plot$data$estimate), set_expect_plot_condition,
  #                          tolerance = set_tolerance, info = test_infor)
  
  
  
  
  # devtools::load_all()
  # marginal_funcall        <- get_predictions
  # marginal_out <- CustomDoCall(marginal_funcall,
  #              marginal_args) %>% data.frame()
  
  
})
