

#' Visualize conditional effects of predictor
#'
#' @details The \strong{plot_conditional_effects()} is a wrapper around the
#'   [brms::conditional_effects()]. The [brms::conditional_effects()] function
#'   from the \pkg{brms} package can used to plot the fitted (distance) curve
#'   when response (e.g., height) is not transformed. However, when the outcome
#'   is log or square root transformed, the [brms::conditional_effects()] will
#'   return the fitted curve on the log or square root scale whereas the
#'   \strong{plot_conditional_effects()} will return the fitted curve on the
#'   original scale. Furthermore, the \strong{plot_conditional_effects()} also
#'   plots the velocity curve on the original scale after making required
#'   back-transformation. Apart from these differences, both these functions
#'   ([brms::conditional_effects] and \strong{plot_conditional_effects()} work
#'   in the same manner. In other words, user can specify all the arguments
#'   which are available in the [brms::conditional_effects()].
#'   
#' @param model An object of class \code{bgmfit}. function.
#'
#' @param resp A character string to specify response variable when processing
#'   posterior draws for the univariate-by-subgroup and multivariate models (see
#'   [bsitar::bsitar()] for details on fitting univariate-by-subgroup and
#'   multivariate models). For univariate model, \code{resp = NULL} (default).
#'   Note that argument \code{resp} must be specified for the
#'   univariate-by-subgroup and multivariate models otherwise it will result in
#'   an error. On the other hand, argument \code{resp} must be \code{NULL} for
#'   the univariate model. The default setting is \code{resp = NULL} assuming a
#'   univariate model.
#'
#' @param deriv An integer to indicate whether to estimate distance curve or
#'   derivatives (velocity and acceleration curves). Default \code{deriv = 0} is
#'   for the distance curve, \code{deriv = 1} for velocity curve, and
#'   \code{deriv = 2} for the acceleration curve.
#'
#' @param deriv_model A logical (default \code{TRUE}) to indicate whether to
#'   estimate model based derivatives or from the differentiation of the
#'   distance curve. When model is fit with \code{decomp = 'QR'}, the only
#'   approach available to estimate derivatives by the  differentiation of the
#'   distance curve.
#' 
#' @inherit brms::conditional_effects description
#' @inherit growthparameters.bgmfit params
#' @inherit fitted_draws.bgmfit params
#'
#' @param ... Additional arguments passed to the [brms::conditional_effects()]
#'   function. Please see [brms::conditional_effects()] for details.
#'
#' @return An object of class 'brms_conditional_effects' which is a named list
#'   with one data.frame per effect containing all information required to
#'   generate conditional effects plots. See brms::conditional_effects for
#'   details.
#'
#' @export plot_conditional_effects.bgmfit
#' @export
#' 
#' @seealso [brms::conditional_effects()]
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
#' # Population average distance curve
#' plot_conditional_effects(model, deriv = 0, re_formula = NA)
#' 
#' \donttest{
#' # Individual-specific distance curves
#' plot_conditional_effects(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' plot_conditional_effects(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' plot_conditional_effects(model, deriv = 1, re_formula = NULL)
#' }
#' 
plot_conditional_effects.bgmfit <-
  function(model,
           resp = NULL,
           deriv = 0,
           deriv_model = TRUE,
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
    
    if(!deriv_model) {
      if(deriv == 1 | deriv == 2) {
        summary <- FALSE
        getfunx <- model$model_info[['exefuns']][[o[[1]]]]
        assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
      }
    }
    
    if(deriv_model) {
      if(deriv == 1 | deriv == 2) {
        getfunx <- model$model_info[['exefuns']][[o[[2]]]]
        assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
      }
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
    
    
    
    if(!deriv_model) {
      xvar  <- model$model_info$xvar
      idvar <- model$model_info$groupvar
      if(length(idvar) > 1) idvar <- idvar[1]
      yvar  <- 'yvar'
      out_    <- brms::conditional_effects(model, resp = resp, ...)
      datace <- out_[[1]] %>% dplyr::select(dplyr::all_of(names(model$data)))
      datace[[idvar]] <- unique(levels(model$data[[idvar]]))[1]
      outx <- fitted_draws(model, resp = resp, newdata = datace,
                           deriv = deriv, ...)
      out_[[1]][['estimate__']] <- outx[, 1]
      out_[[1]][['se__']] <- outx[, 2]
      out_[[1]][['lower__']] <- outx[, 3]
      out_[[1]][['upper__']] <- outx[, 4]
      . <- out_
    }
    
    if(deriv_model) {
     . <- brms::conditional_effects(model, resp = resp, ...)
    }
    
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


#' @rdname plot_conditional_effects.bgmfit
#' @export
plot_conditional_effects <- function(model, ...) {
  UseMethod("plot_conditional_effects")
}





