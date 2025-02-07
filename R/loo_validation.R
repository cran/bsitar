

#' @title Perform leave-one-out (LOO) cross-validation
#' 
#' @description The \strong{loo_validation()} function is a wrapper around the
#'   [brms::loo()] function to perform approximate leave-one-out cross-validation
#'   based on the posterior likelihood. See [brms::loo()] for more details.
#'   
#' @details The function supports model comparisons using [loo::loo_compare()]
#'   for comparing information criteria across models. For \code{bgmfit}
#'   objects, \code{LOO} is simply an alias for \code{loo}. Additionally, you
#'   can use [brms::add_criterion()] to store information criteria in the fitted
#'   model object for later use.
#' 
#' @param compare A logical flag indicating if the information criteria of the
#'   models should be compared using [loo::loo_compare()].
#'   
#' @param moment_match A logical flag to indicate whether
#'   [loo::loo_moment_match()] should be applied to problematic observations.
#'   Defaults to \code{FALSE}. For most models, moment matching will only work
#'   if \code{save_pars = save_pars(all = TRUE)} was set when fitting the model
#'   with [brms::brm()]. See [brms::loo_moment_match()] for more details.
#' 
#' @param reloo A logical flag indicating whether [brms::reloo()] should be
#'   applied to problematic observations. Defaults to \code{FALSE}.
#'
#' @param moment_match_args An optional \code{list} of additional arguments
#'   passed to [loo::loo_moment_match()].
#'
#' @param reloo_args An optional \code{list} of additional arguments passed to
#'   [brms::reloo()].
#' 
#' @inherit brms::loo params 
#' @inheritParams growthparameters.bgmfit
#' @inheritParams plot_ppc.bgmfit
#' 
#' @param ... Additional arguments passed to the [brms::loo()] function. 
#'   Please see \code{brms::loo} for details on various options available.
#' 
#' @return If only one model object is provided, an object of class \code{loo}
#'   is returned. If multiple objects are provided, an object of class
#'   \code{loolist} is returned.
#' 
#' @export loo_validation.bgmfit
#' @export
#' 
#' @seealso [brms::loo()] 
#' 
#' @inherit berkeley author
#'
#' @examples
#' 
#' \donttest{
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether model fit object 'berkeley_exfit' exists
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Perform leave-one-out cross-validation
#' loo_validation(model, cores = 1)
#' }
#' 
loo_validation.bgmfit <-
  function(model,
           compare = TRUE,
           resp = NULL,
           dpar = NULL,
           pointwise = FALSE,
           moment_match = FALSE,
           reloo = FALSE,
           k_threshold = 0.7,
           save_psis = FALSE,
           moment_match_args = list(),
           reloo_args = list(),
           model_names = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           cores = 1,
           deriv_model = NULL,
           verbose = FALSE,
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
    
    if(!is.null(ndraws)) {
      if(ndraws == 1) stop("ndraws must be greater than 1")
    }
    
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
    
    

    calling.args$x <- full.args$model
    calling.args$object <- NULL
    calling.args$model <- NULL
  
    if(is.null(calling.args$newdata)) {
      if(!is.null(newdata)) calling.args$newdata <- newdata
    }
    
    
    . <- brms::loo(model,
                   compare = compare,
                   resp = resp,
                   pointwise = pointwise,
                   moment_match = moment_match,
                   reloo = reloo,
                   k_threshold = k_threshold,
                   save_psis = save_psis,
                   moment_match_args = moment_match_args,
                   reloo_args = reloo_args,
                   model_names = model_names,
                   ndraws = ndraws,
                   cores = cores, 
                    ...)
    

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



#' @rdname loo_validation.bgmfit
#' @export
loo_validation <- function(model, ...) {
  UseMethod("loo_validation")
}

