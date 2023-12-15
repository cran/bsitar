

#' Perform leave-one-out (loo) cross-validation
#' 
#' @description The \strong{loo_validation()} is a wrapper around the
#'   [brms::loo()] function to perform approximate leave-one-out
#'   cross-validation based on the posterior likelihood. See [brms::loo()] for
#'   details.
#' 
#' @inherit brms::loo details 
#' 
#' @inheritParams growthparameters.bgmfit
#' 
#' @inheritParams plot_ppc.bgmfit
#' 
#' @param deriv Must be \code{NULL}. 
#' 
#' @param ... Additional arguments passed to the [brms::loo()] function. 
#' Please see \code{brms::loo} for details on various options available.
#' 
#' @return If only one model object is provided, then an object of class
#'   \code{loo} is returned. If multiple objects are provided, an object of
#'   class \code{loolist}.
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
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid fitting the model which takes time, the model  
#' # fit has already been saved as 'berkeley_mfit.rda' file.
#' # See examples section of the main function for details on the model fit.
#' 
#' model <- berkeley_mfit
#' 
#' \donttest{
#' loo_validation(model, cores = 1)
#' }
#' 
#' 
loo_validation.bgmfit <-
  function(model,
           resp = NULL,
           cores = 1,
           deriv = 0,
           usesavedfuns = FALSE,
           clearenvfuns = FALSE,
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      envir <- parent.frame()
    }
    
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             envir = envir,
                             deriv = deriv,
                             all = FALSE)
    
   
    if(deriv == 0) {
      getfunx <- model$model_info[['exefuns']][[o[[2]]]]
      assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    }
    
    
    
    
    if(!usesavedfuns) {
      if(is.null(check_if_functions_exists(model, o, model$xcall))) {
        return(invisible(NULL))
      }
    }
    
    
    if(usesavedfuns) {
      if(is.null(check_if_functions_exists(model, o, model$xcall))) {
        oall <-
          post_processing_checks(model = model,
                                 xcall = match.call(),
                                 resp = resp,
                                 envir = envir,
                                 deriv = deriv,
                                 all = TRUE)
        tempgenv <- envir
        oalli_c <- c()
        oalli_c <- c(oalli_c, paste0(o[[1]], "0"))
        for (oalli in names(oall)) {
          if(!grepl(o[[1]], oalli)) {
            oalli_c <- c(oalli_c, oalli)
          }
        }
        for (oalli in oalli_c) {
          assign(oalli, oall[[oalli]], envir = tempgenv)
        }
        assign(o[[1]], getfunx, envir = tempgenv)
      }
    }
    
    
    . <- brms::loo(model, resp = resp, cores = cores ,...)
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    
    if(!is.null(clearenvfuns)) {
      if(!is.logical(clearenvfuns)) {
        stop('clearenvfuns must be NULL or a logical')
      } else {
        setcleanup <- clearenvfuns
      }
    }
    
    if(setcleanup) {
      tempgenv <- envir
      for (oalli in names(oall)) {
        if(exists(oalli, envir = tempgenv )) {
          remove(list=oalli, envir = tempgenv)
        }
      }
    }
    
    .
  }


#' @rdname loo_validation.bgmfit
#' @export
loo_validation <- function(model, ...) {
  UseMethod("loo_validation")
}

