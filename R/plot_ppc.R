


#' @title Perform posterior predictive distribution checks
#' 
#' @details The \strong{plot_ppc()} function is a wrapper around the
#'   [brms::pp_check()] function, which allows for the visualization of
#'   posterior predictive checks.
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @inheritParams growthparameters.bgmfit
#' @inheritParams brms::pp_check.brmsfit
#' @inheritParams fitted_draws.bgmfit
#' @inheritParams bayesplot::ppc_dens_overlay
#' @inheritParams bayesplot::ppc_ecdf_overlay
#' @inheritParams bayesplot::ppc_freqpoly
#' @inheritParams bayesplot::ppc_violin_grouped
#' @inheritParams bayesplot::ppc_hist
#' @inherit brms::pp_check.brmsfit description
#' 
#' @param ... Additional arguments passed to the [brms::pp_check.brmsfit()] 
#'   function. Please refer to [brms::pp_check.brmsfit()] for details.
#' 
#' @return A \code{ggplot} object that can be further customized using the 
#'   \pkg{ggplot2} package.
#' 
#' @export
#' 
#' @inherit berkeley author
#'
#' @examples
#' 
#' \donttest{
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation, which takes time, the Bayesian SITAR model is fit to 
#' # the 'berkeley_exdata' and saved as an example fit ('berkeley_exfit').
#' # See the 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether the model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' plot_ppc(model, ndraws = 100)
#' }
#' 
plot_ppc.bgmfit <-
  function(model,
           type,
           ndraws = NULL,
           dpar = NULL,
           draw_ids = NULL,
           prefix = c("ppc", "ppd"),
           group = NULL,
           x = NULL,
           newdata = NULL,
           resp = NULL,
           size = 0.25,
           alpha = 0.7,
           trim = FALSE,
           bw = "nrd0",
           adjust = 1,
           kernel = "gaussian",
           n_dens = 1024,
           pad = TRUE,
           discrete = FALSE,
           binwidth = NULL,
           bins = NULL,
           breaks = NULL,
           freq = TRUE,
           y_draw = c("violin", "points", "both"),
           y_size = 1,
           y_alpha = 1,
           y_jitter = 0.1,
           verbose = FALSE,
           deriv_model = NULL,
           dummy_to_factor = NULL, 
           expose_function = FALSE,
           usesavedfuns = NULL,
           clearenvfuns = NULL,
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- parent.frame()
    }
    
    # Depending on dpar 'mu' or 'sigma', subset model_info
    model <- getmodel_info(model = model, dpar = dpar)
    

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
    } else { 
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
    
    
    check_if_package_installed(model, xcall = NULL)
    
    
    if(is.null(ndraws)) {
      ndraws <- brms::ndraws(model)
    }
    
    if(is.null(deriv_model)) {
      deriv_model <- TRUE
    }
    
  
    
    full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                           fargs = formals(), 
                                           dargs = list(...), 
                                           verbose = verbose)
    
    full.args$model <- model
    full.args$deriv <- deriv <- 0
    
    if(!is.null(model$xcall)) {
      arguments <- get_args_(as.list(match.call())[-1], model$xcall)
      newdata <- newdata
    } else {
      newdata <- do.call(get.newdata, full.args)
    }
    
    
    if(!is.null(model$model_info$decomp)) {
      if(model$model_info$decomp == "QR") deriv_model<- FALSE
    }
    
    expose_method_set <- model$model_info[['expose_method']]
    
    model$model_info[['expose_method']] <- 'NA' # Over ride method 'R'
    
    o <- post_processing_checks(model = model,
                                xcall = match.call(),
                                resp = resp,
                                envir = envir,
                                deriv = 'deriv', 
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
                               minimum_version = get_package_minversion('brms'), 
                               prompt = FALSE,
                               stop = FALSE,
                               verbose = FALSE))) {
      if(is.null(check_if_functions_exists(model, o, model$xcall,
                                           usesavedfuns = usesavedfuns))) {
        return(invisible(NULL))
      }
    }
    
    
    misc <- c("verbose", "usesavedfuns", "clearenvfuns", 
              "envir", "fullframe", "dummy_to_factor")
    calling.args <- post_processing_args_sanitize(model = model,
                                          xcall = match.call(),
                                          resp = resp,
                                          envir = envir,
                                          deriv = deriv, 
                                          dots = list(...),
                                          misc = misc,
                                          verbose = verbose)
   
    
    calling.args$object <- full.args$model
    if(is.null(calling.args$newdata)) {
      if(!is.null(newdata)) calling.args$newdata <- newdata
    }
   
    
    . <- do.call(brms::pp_check, calling.args)
    
   
    # Restore function(s)
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    
    if(!is.null(eval(full.args$clearenvfuns))) {
      if(!is.logical(eval(full.args$clearenvfuns))) {
        stop('clearenvfuns must be NULL or a logical')
      } else {
        setcleanup <- eval(full.args$clearenvfuns)
      }
    }
    
    if(is.null(eval(full.args$clearenvfuns))) {
      if(is.null(eval(full.args$usesavedfuns))) {
        full.args$usesavedfuns <- usesavedfuns
      }
      if(eval(full.args$usesavedfuns)) {
        setcleanup <- TRUE 
      } else {
        setcleanup <- FALSE
      }
    }
    
    # Cleanup environment if requested
    if(setcleanup) {
      suppressWarnings({
        tempgenv <- envir
        for (oalli in names(oall)) {
          if(exists(oalli, envir = tempgenv )) {
            remove(list=oalli, envir = tempgenv)
          }
        }
        tempgenv <- test
        for (oalli in names(oall)) {
          if(exists(oalli, envir = tempgenv )) {
            remove(list=oalli, envir = tempgenv)
          }
        }
      })
    } # if(setcleanup) {
    
    .
  }


#' @rdname plot_ppc.bgmfit
#' @export
plot_ppc <- function(model, ...) {
  UseMethod("plot_ppc")
}


