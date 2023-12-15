

#' Perform posterior predictive distribution checks
#' 
#' @details The \strong{plot_ppc()} is a wrapper around the [brms::pp_check()].
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @inherit growthparameters.bgmfit params
#' @inherit brms::pp_check.brmsfit description 
#' @inherit fitted_draws.bgmfit params
#' 
#' @param ... Additional arguments passed to the [brms::pp_check.brmsfit()] 
#' function. Please see [brms::pp_check.brmsfit()] for details.
#' 
#' @return A ggplot object that can be further customized using the
#' ggplot2 package.
#' 
#' @export
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
#' plot_ppc(model, ndraws = 100)
#' 
plot_ppc.bgmfit <-
  function(model,
           resp = NULL,
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
    

    . <- brms::pp_check(model, resp = resp, ...)
    
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


#' @rdname plot_ppc.bgmfit
#' @export
plot_ppc <- function(model, ...) {
  UseMethod("plot_ppc")
}


