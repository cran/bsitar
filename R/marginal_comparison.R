

#' Compare growth curves
#' 
#'@description The \strong{marginal_comparison()} function estimates and
#'  compare growth curves such as distance and velocity. This function is a wrapper around the
#'  [marginaleffects::comparisons()] and [marginaleffects::avg_comparisons()].
#'  The [marginaleffects::comparisons()] computes unit-level (conditional)
#'  estimates whereas [marginaleffects::avg_comparisons()] return average
#'  (marginal) estimates. A detailed explanation is available
#'  [here](https://marginaleffects.com). 
#' 
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @param datagrid Generate a grid of user-specified values for use in the
#'   \code{newdata} argument in various functions of the \pkg{marginaleffects}
#'   package. This is useful to define where in the predictor space we want to
#'   evaluate the quantities of interest. See [marginaleffects::datagrid()] for
#'   details. The default value for the \code{datagrid} is \code{NULL} implying
#'   that no custom grid is constructed. To set a data grid, the argument should
#'   be a data.frame constructed by using the [marginaleffects::datagrid()]
#'   function, or else a named list which are internally used for setting up the
#'   grid. For the user convenience, we also allow setting an empty list
#'   \code{datagrid = list()} in which case essential arguments such as
#'   \code{model}, \code{newdata} are taken up from the respective arguments
#'   specified elsewhere. Further, the level 1 predictor (such as age) and any
#'   covariate included in the model fit (e.g., gender) are also automatically
#'   inferred from the \code{model} object.
#' 
#' @param digits An integer (default \code{2}) to set the decimal places for the
#'   estimates. The \code{digits} is passed on to the
#'   [base::round()] function.
#'   
#' @param average A logical to indicate whether to internally call the
#'    [marginaleffects::comparisons()] or the
#'    [marginaleffects::avg_comparisons()] function. If \code{FALSE} (default),
#'    [marginaleffects::comparisons()] is called otherwise
#'    [marginaleffects::avg_comparisons()] when \code{average = TRUE}.
#'
#' @param plot A logical to specify whether to plot comparisons by calling the
#'   [marginaleffects::plot_comparisons()] function (\code{FALSE}) or not
#'   (\code{FALSE}). If \code{FALSE} (default), then
#'   [marginaleffects::comparisons()] or [marginaleffects::avg_comparisons()]
#'   are called to compute predictions (see \code{average} for details)
#'   
#' @param deriv A numeric to specify whether to estimate parameters based on the
#'   differentiation of the distance curve or the model based first derivative.
#'   Please see argument \code{variables} for more details.  
#'   
#' @param reformat A logical (default \code{TRUE}) to reformat the  output
#'   returned by the \code{marginaleffects} as a data.frame with column names
#'   re-defined as follows: \code{conf.low} as \code{Q2.5}, and \code{conf.high}
#'   as \code{Q97.5} (assuming that \code{conf_int = 0.95}). Also, following
#'   columns are dropped from the data frame: \code{term}, \code{contrast},
#'   \code{tmp_idx}, \code{predicted_lo}, \code{predicted_hi}, \code{predicted}.
#' 
#' @inheritParams  growthparameters.bgmfit
#' @inheritParams  growthparameters_comparison.bgmfit
#' @inheritParams  marginal_draws.bgmfit
#' @inheritParams  marginaleffects::comparisons
#' @inheritParams  marginaleffects::avg_comparisons
#' @inheritParams  marginaleffects::plot_comparisons
#' @inheritParams  marginaleffects::datagrid
#' @inheritParams  brms::fitted.brmsfit
#'
#' @return A data frame objects with estimates and CIs for computed parameter(s)
#' 
#' @export marginal_comparison.bgmfit
#' @export
#' 
#' @seealso [marginaleffects::comparisons()]
#'   [marginaleffects::avg_comparisons()]
#'   [marginaleffects::plot_comparisons()]
#' 
#' @references
#' \insertAllCited{}
#' 
#' @inherit berkeley author
#' 
#' @examples
#' 
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' marginal_comparison(model, parameter = 'apgv', draw_ids = 1)
#' }
#' 
marginal_comparison.bgmfit <- function(model,
                                   resp = NULL,
                                   ndraws = NULL,
                                   draw_ids = NULL,
                                   newdata = NULL,
                                   datagrid = NULL,
                                   re_formula = NA,
                                   allow_new_levels = FALSE,
                                   sample_new_levels = "gaussian",
                                   xrange = 1,
                                   digits = 2,
                                   numeric_cov_at = NULL,
                                   aux_variables = NULL,
                                   levels_id = NULL,
                                   avg_reffects = NULL,
                                   idata_method = NULL,
                                   ipts = NULL,
                                   seed = 123,
                                   future = FALSE,
                                   future_session = 'multisession',
                                   cores = NULL,
                                   average = FALSE, 
                                   plot = FALSE, 
                                   showlegends = NULL, 
                                   variables = NULL,
                                   deriv = NULL,
                                   deriv_model = NULL,
                                   comparison = "difference",
                                   type = NULL,
                                   by = FALSE,
                                   conf_level = 0.95,
                                   transform = NULL,
                                   cross = FALSE,
                                   wts = NULL,
                                   hypothesis = NULL,
                                   equivalence = NULL,
                                   eps = NULL,
                                   reformat = NULL,
                                   estimate_center = NULL,
                                   estimate_interval = NULL,
                                   dummy_to_factor = NULL, 
                                   verbose = FALSE,
                                   expose_function = FALSE,
                                   usesavedfuns = NULL,
                                   clearenvfuns = NULL,
                                   envir = NULL,
                                   ...) {
  
  
  if(!is.null(estimate_center)) {
    ec_ <- getOption("marginaleffects_posterior_center")
    options("marginaleffects_posterior_center" = estimate_center)
    on.exit(options("marginaleffects_posterior_center" = ec_), add = TRUE)
  }
  if(!is.null(estimate_interval)) {
    ei_ <- getOption("marginaleffects_posterior_interval")
    options("marginaleffects_posterior_interval" = estimate_interval)
    on.exit(options("marginaleffects_posterior_interval" = ei_), add = TRUE)
  }
  
  try(zz <- insight::check_if_installed(c("marginaleffects"), 
                                        minversion = 
                                          get_package_minversion(
                                            'marginaleffects'
                                          ), 
                                        prompt = FALSE,
                                        stop = FALSE))
  
  if(!isTRUE(zz)) {
    message("Please install the latest version of the 'marginaleffects' package",
            "\n ",
            "remotes::install_github('vincentarelbundock/marginaleffects')")
    return(invisible(NULL))
  }
  
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- parent.frame()
  }
  
  
  if(is.null(usesavedfuns)) {
    if(!is.null(model$model_info$exefuns[[1]])) {
      usesavedfuns <- TRUE
    } else if(is.null(model$model_info$exefuns[[1]])) {
      if(expose_function) {
        model <- expose_model_functions(model, envir = envir)
        usesavedfuns <- TRUE
      } else if(!expose_function) {
        usesavedfuns <- FALSE
      }
    }
  } else { # if(!is.null(usesavedfuns)) {
    if(!usesavedfuns) {
      if(expose_function) {
        model <- expose_model_functions(model, envir = envir)
        usesavedfuns <- TRUE
      }
    } else if(usesavedfuns) {
      check_if_functions_exists(model, checks = TRUE, 
                                usesavedfuns = usesavedfuns)
    }
  }
  
  
  
  if(is.null(ndraws)) {
    ndraws <- brms::ndraws(model)
  }
  
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  xvar_  <- paste0('xvar', resp_rev_)
  xvar   <- model$model_info[[xvar_]]
  cov_   <- paste0('cov', resp_rev_)
  cov    <- model$model_info[[cov_]]
  uvarby <- model$model_info$univariate_by
  
  # Note here, newdata is not model$data but rather model$model_info$bgmfit.data
  # This was must for univariate_by
  if(is.null(newdata)) {
    #newdata <- model$model_info$bgmfit.data
  }
  
  if(!is.na(uvarby)) {
    uvarby_ind <- paste0(uvarby, resp)
    varne <- paste0(uvarby, resp)
   # newdata <- newdata %>% dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
  }
  
  
  if(is.null(deriv) & is.null(deriv_model)) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(deriv == 0 & is.null(deriv_model)) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(deriv == 1 & is.null(deriv_model)) {
    deriv <- 1
    deriv_model <- TRUE
  } else if(is.null(deriv) & !deriv_model) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(is.null(deriv) & deriv_model) {
    deriv <- 1
    deriv_model <- TRUE
  }
   
  # The deriv_model is a placeholder in marginaleffects
  # if(is.null(deriv_model)) {
  #   deriv_model <- TRUE
  # }
  
  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  
  if(idata_method == 'm1') {
    stop("For marginaleffects based functions, the " ,
         " \n",
         " 'idata_method' argument must be either NULL or 'm2'" )
  }
  
  
  if (is.null(eps)) eps <- 1e-6
  
  
  # Initiate non formalArgs()
  term <- NULL;
  contrast <- NULL;
  tmp_idx <- NULL;
  predicted_lo <- NULL;
  predicted_hi <- NULL;
  predicted <- NULL;
  conf.high <- NULL;
  conf.low <- NULL;
  estimate <- NULL;
  `:=` <- NULL;
  `.` <- NULL;
  
  
 
  
  conf <- conf_level
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  # set_names_  <- c('Estimate', 'Est.Error', probtitles)
  set_names_  <- c('Estimate', probtitles)
  
  if(!is.null(model$model_info$decomp)) {
    if(model$model_info$decomp == "QR") deriv_model<- FALSE
  }
  
  expose_method_set <- model$model_info[['expose_method']]
  
  model$model_info[['expose_method']] <- 'NA' # Over ride method 'R'
  
  o <- post_processing_checks(model = model,
                              xcall = match.call(),
                              resp = resp,
                              envir = envir,
                              deriv = deriv, 
                              all = FALSE,
                              verbose = verbose)
  
  oall <- post_processing_checks(model = model,
                                 xcall = match.call(),
                                 resp = resp,
                                 envir = envir,
                                 deriv = deriv, 
                                 all = TRUE,
                                 verbose = FALSE)
  
  
  test <- setupfuns(model = model, resp = resp,
                    o = o, oall = oall, 
                    usesavedfuns = usesavedfuns, 
                    deriv = deriv, envir = envir, 
                    deriv_model = deriv_model, 
                    ...)
  
  if(is.null(test)) return(invisible(NULL))
  
  if(!isTRUE(
    check_pkg_version_exists('brms', 
                             minversion = get_package_minversion('brms'),
                             prompt = FALSE,
                             stop = FALSE,
                             verbose = FALSE))) {
    if(is.null(check_if_functions_exists(model, o, model$xcall,
                                         usesavedfuns = usesavedfuns))) {
      return(invisible(NULL))
    }
  }
  
  
  xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
  scall <- sys.calls()
  
  get_xcall <- function(xcall, scall) {
    scall <- scall[[length(scall)]]
    if(any(grepl("marginal_comparison", scall, fixed = T)) |
       any(grepl("marginal_comparison.bgmfit", scall, fixed = T))) {
      xcall <- "marginal_comparison"
    } else {
      xcall <- xcall
    } 
  }
  
  if(!is.null(model$xcall)) {
    if(model$xcall == "marginal_comparison") {
      xcall <- "marginal_comparison"
    }
  } else {
    scall <- sys.calls()
    xcall <- get_xcall(xcall, scall)
  }
  
  
  xcall <- xcall
  
  check_if_package_installed(model, xcall = xcall)
  
  model$xcall <- xcall
  
  
  arguments <- get_args_(as.list(match.call())[-1], xcall)
  arguments$model <- model
  arguments$usesavedfuns <- usesavedfuns
  
  
  get.cores_ <- get.cores(arguments$cores)
  arguments$cores <- setincores <-  get.cores_[['max.cores']]
  .cores_ps <- get.cores_[['.cores_ps']]
  
  if (future) {
    if (future_session == 'multisession') {
      future::plan('multisession', workers = setincores)
    } else if (future_session == 'multicore') {
      future::plan('multicore', workers = setincores)
    }
  }
  
  
  
  full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                  # fargs = formals(), 
                                  fargs = arguments, 
                                  dargs = list(...), 
                                  verbose = verbose)
  
  full.args$model <- model
  full.args$deriv_model <- deriv_model
  
  full.args$newdata <- newdata
  newdata           <- do.call(get.newdata, full.args)
  
  if(!is.na(uvarby)) {
    uvarby_ind <- paste0(uvarby, resp)
    varne <- paste0(uvarby, resp)
     newdata <- newdata %>% dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
  }
  full.args$newdata <- newdata
  
  # keeping ... cause marginaleffects:: argument is missing, with no default
  full.args[["..."]] <- NULL
  
  
  comparisons_arguments <- full.args
  
  # Drop that not required for marginaleffects::
  exclude_args <- as.character(quote(
    c(
      xrange,
      digits,
      numeric_cov_at,
      acg_asymptote,
      levels_id,
      avg_reffects,
      deriv,
      deriv_model,
      ipts,
      seed,
      future,
      future_session,
      dummy_to_factor,
      verbose,
      expose_function,
      usesavedfuns,
      clearenvfuns,
      envir,
      plot,
      showlegends,
      average,
      parameter,
      estimate_center,
      estimate_interval
    )
  ))[-1]
  
  if(plot) {
    exclude_args <- c(exclude_args, "cross")
  }
  
  for (exclude_argsi in exclude_args) {
    comparisons_arguments[[exclude_argsi]] <- NULL
  }
  
  
  if (!is.null(variables)) {
    if (!is.list(variables)) {
      set_variables <- variables
    } else if (is.list(variables)) {
      set_variables <- variables
      if(is.null(set_variables[[xvar]])) {
        if(deriv == 0) set_variables[[xvar]] <- eps
        if(deriv > 0)  set_variables[[xvar]] <- 0
      } else if(!is.null(set_variables[[xvar]])) {
        if(eval(set_variables[[xvar]]) !=0) {
          if(verbose) {
            message("The value of ", xvar, " is not same as used in the ",
                    " \n", 
                    " model fit. Please check if this is intended")
          }
        }
      }
    }
  } else if (is.null(variables)) {
    if(deriv == 0) set_variables <- list(eps)
    if(deriv > 0)  set_variables <- list(0)
    names(set_variables) <- xvar
  } 
  
  
  
  if(is.null(by)) {
    if(is.null(cov)) {
      set_group <- FALSE
    } else if(!is.null(cov)) {
      set_group <- cov
      if (!set_group %in% cov) {
        stop('by must be one of the ', cov)
      } 
    }
  } else if(!is.null(by)) {
    if (!isFALSE(by)) {
      set_group <- by
    } else if (isFALSE(by)) {
      set_group <- FALSE
    }
  }
  
 
  
  
  
  
    
    comparisons_arguments$variables  <- set_variables
    comparisons_arguments$by         <- set_group

    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    
    
    if(is.null(showlegends)) {
      if(is.null(comparisons_arguments$re_formula)) {
        showlegends <- FALSE
      } else {
        showlegends <- TRUE
      }
    }
    
    
    # Set up datagrid
    
    if(!is.null(datagrid)) {
      if(is.data.frame(datagrid)) {
        set_datagrid <- datagrid
        comparisons_arguments$newdata <- set_datagrid
      } else if(is.list(datagrid)) {
        if(is.null(datagrid[['model']])) setmodel <- model else setmodel <- datagrid$model
        if(is.null(datagrid[['newdata']])) setnewdata <- newdata else setnewdata <- datagrid$newdata
        if(is.null(datagrid[['grid_type']])) setgrid_type <- "mean_or_mode" else setgrid_type <- datagrid$grid_type
        if(is.null(datagrid[[xvar]])) setxvar <- newdata[[xvar]] else setxvar <- datagrid$newdata[[xvar]]
        datagrid_arguments <- list(model = setmodel,
                                   newdata = setnewdata,
                                   grid_type = setgrid_type)
        datagrid_arguments[[xvar]] <- setxvar
        if(setgrid_type == "mean_or_mode") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- set_group
        } else if(setgrid_type == "balanced") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- NULL
          # correctly set comparisons_arguments[['by']]  too 
          comparisons_arguments[['by']] <- NULL
        }
        set_datagrid <- do.call(marginaleffects::datagrid, datagrid_arguments)
        comparisons_arguments$newdata <- set_datagrid
      } else {
        stop("datagrid should be a data frame or named list")
      }
    } else if(is.null(datagrid)) {
      comparisons_arguments$newdata <- comparisons_arguments$newdata
    }
    
    # The datagrid argument is not allowed. It served its purpose by defining 
    # the newdata. So remove it from the arguments
    
    comparisons_arguments[['datagrid']] <- NULL
    
    # Somehow draw_ids not passed correctly if not specified explicitly as arg
    get_draw_ids <- comparisons_arguments[['draw_ids']]
    if(is.null(eval(get_draw_ids))) {
      set_draw_ids <- NULL
    } else if(is.numeric(eval(get_draw_ids))) {
      set_draw_ids <- get_draw_ids
    } else if(!eval(get_draw_ids)) {
      set_draw_ids <- NULL
    }
    comparisons_arguments[['draw_ids']] <- set_draw_ids
    
    # comparisons_arguments$average <- NULL
    # comparisons_arguments$parameter <- NULL
   # comparisons_arguments$... <- NULL
    
    
    suppressWarnings({
      if(!plot) {
        if(!average) {
          out <- do.call(marginaleffects::comparisons, comparisons_arguments)
        } else if(average) {
          out <- do.call(marginaleffects::avg_comparisons, comparisons_arguments)
        }
      }
      if(plot) {
        if(isFALSE(set_group)) comparisons_arguments$by <- NULL
        out <- do.call(marginaleffects::plot_comparisons, comparisons_arguments)
        outp <- out
        if(!showlegends) outp <- outp + ggplot2::theme(legend.position = 'none')
        return(outp)
      }
    })
    
    out_sf <- out
  
  
  

  
  
  out_sf <- out_sf %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                         ~ round(., digits = digits))) %>%
    data.frame()
  
 
  if(is.null(reformat)) {
    if(is.null(hypothesis) && is.null(equivalence)) {
      reformat <- TRUE
    } else {
      reformat <- FALSE
    }
  }
  
  if (reformat) {
    out_sf <- out_sf %>% 
      dplyr::rename(!!as.symbol(set_names_[1]) := estimate) %>% 
      dplyr::rename(!!as.symbol(set_names_[2]) := conf.low) %>% 
      dplyr::rename(!!as.symbol(set_names_[3]) := conf.high) 
      data.frame()
    
    remove_cols_ <- c('term', 'contrast', 'tmp_idx', 'predicted_lo', 
                      'predicted_hi', 'predicted', 'rowid')
    
    out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]
    row.names(out_sf) <- NULL
  }
  
  return(out_sf)
}




#' @rdname marginal_comparison.bgmfit
#' @export
marginal_comparison <- function(model, ...) {
  UseMethod("marginal_comparison")
}


