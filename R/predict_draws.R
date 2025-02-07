



#' @title Predicted values from the posterior predictive distribution
#'
#' @description The \strong{predict_draws()} function is a wrapper around the
#'   [brms::predict.brmsfit()] function, which obtains predicted values (and
#'   their summary) from the posterior distribution. See
#'   [brms::predict.brmsfit()] for details.
#'
#' @details The \strong{predict_draws()} function computes the fitted values
#'   from the posterior distribution. The [brms::predict.brmsfit()] function
#'   from the \pkg{brms} package can be used to obtain predicted (distance)
#'   values when the outcome (e.g., height) is untransformed. However, when the
#'   outcome is log or square root transformed, the [brms::predict.brmsfit()]
#'   function will return the fitted curve on the log or square root scale. In
#'   contrast, the \strong{predict_draws()} function returns the fitted values
#'   on the original scale. Furthermore, \strong{predict_draws()} also computes
#'   the first derivative (velocity), again on the original scale, after making
#'   the necessary back-transformation. Aside from these differences, both
#'   functions ([brms::predict.brmsfit()] and \strong{predict_draws()}) work
#'   similarly. In other words, the user can specify all the options available
#'   in [brms::predict.brmsfit()].
#'
#' @param ... Additional arguments passed to the [brms::predict.brmsfit()]
#'   function. Please see [brms::predict.brmsfit()] for details on the various
#'   options available.
#'
#' @return An array of predicted response values. See [brms::predict.brmsfit()]
#'   for details.
#'
#' @inherit growthparameters.bgmfit params
#' @inherit fitted_draws.bgmfit params
#' @inherit brms::predict.brmsfit params
#'
#' @export predict_draws.bgmfit
#' @export
#'
#' @seealso [brms::predict.brmsfit()]
#'
#' @inherit berkeley author
#'
#' @examples
#' 
#' \donttest{
#'
#' # Fit Bayesian SITAR model
#'
#' # To avoid mode estimation, which takes time, the Bayesian SITAR model is fit 
#' # to the 'berkeley_exdata' and saved as an example fit ('berkeley_exfit').
#' # See the 'bsitar' function for details on 'berkeley_exdata' and 
#' # berkeley_exfit'.
#'
#' # Check and confirm whether the model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#'
#' model <- berkeley_exfit
#' 
#' # Population average distance curve
#' predict_draws(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' predict_draws(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' predict_draws(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' predict_draws(model, deriv = 1, re_formula = NULL)
#'  }
#'  
predict_draws.bgmfit <-
  function(model,
           newdata = NULL,
           resp = NULL,
           dpar = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           re_formula = NA,
           allow_new_levels = FALSE,
           sample_new_levels = "uncertainty",
           incl_autocor = TRUE,
           numeric_cov_at = NULL,
           levels_id = NULL,
           avg_reffects = NULL,
           aux_variables = NULL,
           ipts = 10,
           deriv = 0,
           deriv_model = TRUE,
           summary = TRUE,
           robust = FALSE,
           transform = NULL,
           probs = c(0.025, 0.975),
           xrange = NULL,
           xrange_search = NULL,
           parms_eval = FALSE,
           parms_method = 'getPeak',
           idata_method = NULL,
           verbose = FALSE,
           fullframe = NULL,
           dummy_to_factor = NULL, 
           expose_function = FALSE,
           usesavedfuns = NULL,
           clearenvfuns = NULL,
           funlist = NULL,
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
    
    if (is.null(idata_method)) {
      idata_method <- 'm2'
    }
    
    
    
    # This in plot_conditional_effects_calling if(!eval(full.args$deriv_model)){
    plot_conditional_effects_calling <- FALSE
    syscalls1 <- sys.calls()[[1]]
    syscallsall <- paste(deparse(syscalls1), collapse = "\n")
    for (xc in 1:length(syscallsall)) {
      if(any(grepl('plot_conditional_effects', syscallsall[[xc]]))) {
        plot_conditional_effects_calling <- TRUE
      }
    }
    
    # Checks for newdata and arguments
    # For plot_conditional_effects_calling, newdata is not evaluted
    # For indirectcall i.e.,  model$xcall arguments are passed from the
    # plot_curves() and growthparameters() functions
    
    indirectcall <- FALSE
    if(!plot_conditional_effects_calling) {
      if(!is.null(model$xcall)) {
        arguments <- get_args_(as.list(match.call())[-1], model$xcall)
        full.args <- evaluate_call_args(cargs = arguments, 
                                        fargs = NULL, 
                                        dargs = NULL, 
                                        verbose = verbose)
        full.args$object <- full.args$model
        newdata <- full.args$newdata
        indirectcall <- TRUE
      } else {
        full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                        fargs = formals(), 
                                        dargs = list(...), 
                                        verbose = verbose)
        
        full.args$model <- model
        newdata <- do.call(get.newdata, full.args)
      }
      full.args$newdata <- newdata
    }
    
    
    if(plot_conditional_effects_calling) {
      full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                      fargs = formals(), 
                                      dargs = list(...), 
                                      verbose = verbose)
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
    
    if(!is.null(funlist)) {
      if(!is.list(funlist)) {
        stop("funlist must be a list")
      } else {
        o <- funlist
      }
    }
    
    
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
              "envir", "fullframe")
    
    if(!indirectcall) {
      calling.args <- post_processing_args_sanitize(model = model,
                                                    xcall = match.call(),
                                                    resp = resp,
                                                    envir = envir,
                                                    deriv = deriv, 
                                                    dots = list(...),
                                                    misc = misc,
                                                    verbose = verbose)
    } else if(indirectcall) {
      calling.args <- full.args
    }
    
    calling.args$object <- full.args$model
    if(is.null(calling.args$newdata)) {
      if(!is.null(newdata)) calling.args$newdata <- newdata
    }
    
    
    growthparameters_calling <- FALSE
    syscalls1 <- sys.calls()[[1]]
    syscallsall <- paste(deparse(syscalls1), collapse = "\n")
    for (xc in 1:length(syscallsall)) {
      if(any(grepl('growthparameters', syscallsall[[xc]]))) {
        growthparameters_calling <- TRUE
      }
    }
    if(growthparameters_calling) {
      if(calling.args$re_formula_opt == "V") {
        calling.args$re_formula <- NULL
      } else if(calling.args$re_formula_opt == "v") {
        calling.args$re_formula <- NA
      }
    }
    
    
    
    . <- do.call(predict, calling.args)
    
    
    if(!is.null((eval(full.args$deriv)))) {
      if(eval(full.args$deriv > 0)) { 
        if(!eval(full.args$deriv_model)) {
          full.args$. <- .
          . <- do.call(mapderivqr, full.args)
        } else {
          . <- .
        }
      }
    }
    
    
    
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
    
    
    # fullframe
    full.args$idata_method <- idata_method
    full.args$fullframe <- eval(full.args$fullframe)
    print(full.args$summary)
    if(!is.null(eval(full.args$fullframe))) {
      if(eval(full.args$fullframe)) {
        if(!eval(full.args$summary)) {
          stop("fullframe can not be combined with summary = FALSE")
        }
        if(full.args$idata_method == 'm1') {
          stop("fullframe can not be combined with idata_method = 'm1'")
        }
      }
    }
    if(is.null(eval(full.args$fullframe))) {
      if (!is.na(model$model_info$univariate_by)) {
        if(full.args$idata_method == 'm1') setfullframe <- FALSE
        if(full.args$idata_method == 'm2') setfullframe <- TRUE
      } else {
        setfullframe <- FALSE
      }
    }
    if (!is.na(model$model_info$univariate_by)) {
      if(is.null(full.args$fullframe)) 
        full.args$fullframe <- fullframe <- FALSE
      if(full.args$fullframe & full.args$idata_method == 'm1') 
        setfullframe <- FALSE
      if(full.args$fullframe & full.args$idata_method == 'm2') 
        setfullframe <- TRUE
      if(!full.args$fullframe) 
        setfullframe <- FALSE
      if(setfullframe) {
        uvarby <- model$model_info$univariate_by
        uvarbyresp <- paste0(uvarby, resp)
        uvarbynewdata <- eval(full.args$newdata) %>% 
          dplyr::filter(!!dplyr::sym(uvarbyresp) == 1)
        if(setfullframe) . <- cbind(., uvarbynewdata)
      }
    }
    . 
  }



#' @rdname predict_draws.bgmfit
#' @export
predict_draws <- function(model, ...) {
  UseMethod("predict_draws")
}


