
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


test_that("test-get_growthparameters", {
  skip_on_cran()
  
  
  
  # devtools::load_all()
  
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
      fit               = berkeley_exfit # univariate_fit_cov
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
  
  ###############################################################################
  # set marginaleffects comparisons
  ###############################################################################
  
  set_comparison_for_marginaleffects    <- FALSE
  if(is.null(comparison)) {
    set_comparison_for_marginaleffects <- TRUE
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
  
  
  ###############################################################################
  ###############################################################################
  ###############################################################################
  
  # marginaleffects::comparisons vs get_comparisons
  
  what_test <- 'comparison'
  
  ###############################################################################
  ###############################################################################
  ###############################################################################
  
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
  # Construct test comparisons
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
  # marginaleffects::comparisons vs get_comparisons - average = FALSE
  ###############################################################################
  ###############################################################################
  
  test_str_cat <- 
    "marginaleffects::comparisons vs get_comparisons -> average = FALSE"
  
  
  marginaleffects_funcall <- marginaleffects::comparisons
  marginal_funcall        <- get_growthparameters
  marginal_args_average   <- FALSE
  
  what_test               <- paste0(what_test, "_", "avgF")
  
  if(exists('marginaleffects_args'))      rm('marginaleffects_args')
  if(exists('marginal_args'))             rm('marginal_args')
  
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
  
  
  
  ###############################################################################
  # marginaleffects
  ###############################################################################
  
  if(set_comparison_for_marginaleffects) {
    marginaleffects_args[['comparison']] <- 
      set_marginaleffects_comparison(deriv=deriv, comparison=comparison)
  }
  
  # marginaleffects_out <- do.call(marginaleffects_funcall, 
  #                                marginaleffects_args) %>% data.frame()
  
  if(set_comparison_for_marginaleffects) {
    marginaleffects_args[['comparison']] <- NULL
  }
  
  
  ###############################################################################
  # set parameters
  ###############################################################################
  
  marginaleffects_args[['parameter']] <- c('apgv', 'pgv')
  
  
  # marginaleffects_args[['variables']]  <- 'sex'
  # marginaleffects_args[['by']]         <- 'sex'
  
  
  ###############################################################################
  # marginal - method = pkg future = F
  ###############################################################################
  
  marginal_args <- marginaleffects_args
  
  marginal_args[['method']]        <- 'pkg'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- marginal_args_average
  marginal_args[['deriv']]         <- deriv
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- F
  marginal_args[['future_method']]  <- set_future_method
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- set_future_splits
  
  
  # marginal_args[['variables']]  <- 'sex'
  # marginal_args[['by']]         <- 'sex'
  marginal_args[['hypothesis']]  <- NULL # ~ pairwise | 'parameter'
  
  # marginal_args[['parameter']]     <- c('all')
  
  # devtools::load_all()
  # marginal_funcall        <- get_growthparameters
  marginal_out_FF <- do.call(marginal_funcall, 
                             marginal_args) # %>% data.frame()
  
  
  
  ###############################################################################
  # marginal - method = pkg future = T
  ###############################################################################
  
  
  marginal_args[['future']]         <- T
  marginal_out_FT <- do.call(marginal_funcall, 
                             marginal_args) # %>% data.frame()
  
  
  ###############################################################################
  # marginal - method = custom, model_deriv = T, future = F
  ###############################################################################
  
  
  marginal_args <- marginaleffects_args
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- marginal_args_average
  marginal_args[['deriv']]         <- deriv
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- F
  marginal_args[['future_method']]  <- set_future_method
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- set_future_splits
  
  
  # marginal_args[['variables']]  <- 'sex'
  # marginal_args[['by']]         <- 'sex'
  marginal_args[['hypothesis']]  <- NULL # ~ pairwise | 'parameter'
  
  # marginal_args[['parameter']]     <- c('all')
  
  
  marginal_out_custom_mdT_FF <- do.call(marginal_funcall, 
                                        marginal_args) %>% data.frame()
  
  
  
  ###############################################################################
  # marginal - method = custom, model_deriv = T, future = T
  ###############################################################################
  
  
  marginal_args <- marginaleffects_args
  
  marginal_args[['model_deriv']]   <- T
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- marginal_args_average
  marginal_args[['deriv']]         <- deriv
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- T
  marginal_args[['future_method']]  <- set_future_method
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- set_future_splits
  
  
  # marginal_args[['variables']]  <- 'sex'
  # marginal_args[['by']]         <- 'sex'
  marginal_args[['hypothesis']]  <- NULL # ~ pairwise | 'parameter'
  
  # marginal_args[['parameter']]     <- c('all')
  
  
  marginal_out_custom_mdT_FT <- do.call(marginal_funcall, 
                                        marginal_args) %>% data.frame()
  
  
  
  ###############################################################################
  # marginal - method = custom, model_deriv = T, future = F
  ###############################################################################
  
  
  marginal_args <- marginaleffects_args
  
  marginal_args[['model_deriv']]   <- F
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- marginal_args_average
  marginal_args[['deriv']]         <- deriv
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- F
  marginal_args[['future_method']]  <- set_future_method
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- set_future_splits
  
  
  # marginal_args[['variables']]  <- 'sex'
  # marginal_args[['by']]         <- 'sex'
  marginal_args[['hypothesis']]  <- NULL # ~ pairwise | 'parameter'
  
  # marginal_args[['parameter']]     <- c('all')
  
  
  marginal_out_custom_mdF_FF <- do.call(marginal_funcall, 
                                        marginal_args) %>% data.frame()
  
  
  
  ###############################################################################
  # marginal - method = custom, model_deriv = T, future = T
  ###############################################################################
  
  
  marginal_args <- marginaleffects_args
  
  marginal_args[['model_deriv']]   <- F
  marginal_args[['method']]        <- 'custom'
  marginal_args[['reformat']]      <- FALSE
  marginal_args[['average']]       <- marginal_args_average
  marginal_args[['deriv']]         <- deriv
  marginal_args[['newdata_fixed']] <- 0
  
  marginal_args[['future']]         <- T
  marginal_args[['future_method']]  <- set_future_method
  marginal_args[['future_session']] <- set_future_session
  marginal_args[['cores']]          <- set_future_cores
  marginal_args[['future_splits']]  <- set_future_splits
  
  
  # marginal_args[['variables']]  <- 'sex'
  # marginal_args[['by']]         <- 'sex'
  marginal_args[['hypothesis']]  <- NULL # ~ pairwise | 'parameter'
  
  # marginal_args[['parameter']]     <- c('all')
  
  
  marginal_out_custom_mdF_FT <- do.call(marginal_funcall, 
                                        marginal_args) %>% data.frame()
  
  
  
  
  expect_equal(marginal_out_FF, marginal_out_FT, 
               tolerance = test_tolerance)
  expect_equal(marginal_out_FF, marginal_out_custom_mdT_FF,
               tolerance = test_tolerance)
  expect_equal(marginal_out_FF, marginal_out_custom_mdT_FT, 
               tolerance = test_tolerance)
  expect_equal(marginal_out_FF, marginal_out_custom_mdF_FF,
               tolerance = test_tolerance)
  expect_equal(marginal_out_FF, marginal_out_custom_mdF_FT, 
               tolerance = test_tolerance)
  
  
  
})
