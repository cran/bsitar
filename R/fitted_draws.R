


#' @title Estimate fitted (Expected) values for the Bayesian SITAR model
#' 
#' @description The \strong{fitted_draws()} function is a wrapper around the
#'   [brms::fitted.brmsfit()] function, which allows users to obtain fitted
#'   values (and their summaries) from the posterior draws. For more details,
#'   refer to the documentation for [brms::fitted.brmsfit()]. An alternative
#'   approach is to [get_predictions()] function which is based on the
#'   \pkg{marginaleffects}.
#' 
#' @details The \strong{fitted_draws()} function computes the fitted values from
#'   the posterior draws. While the [brms::fitted.brmsfit()] function from the
#'   \pkg{brms} package can be used to obtain fitted (distance) values when the
#'   outcome (e.g., height) is untransformed, it returns fitted values on the
#'   log or square root scale if the outcome is transformed. In contrast,
#'   \strong{fitted_draws()} returns fitted values on the original scale.
#'   Additionally, \strong{fitted_draws()} computes the first derivative
#'   (velocity) on the original scale, after applying the necessary
#'   back-transformation. Apart from these differences, both
#'   functions—[brms::fitted.brmsfit()] and [fitted_draws()]—operate in the same
#'   manner, allowing users to specify all options available in
#'   [brms::fitted.brmsfit()].
#'   
#' @param deriv An integer indicating whether to estimate the distance curve or
#'   its derivative (velocity curve). The default \code{deriv = 0} is for the
#'   distance curve, while \code{deriv = 1} is for the velocity curve.
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit()]
#'   function. For details on available options, please refer to
#'   \code{brms::fitted.brmsfit()}.
#'   
#' @inheritParams growthparameters.bgmfit
#' @inheritParams brms::fitted.brmsfit
#' 
#' @return An array of predicted mean response values when \code{summarise =
#'   FALSE}, or a \code{data.frame} when \code{summarise = TRUE}. For further
#'   details, refer to [brms::fitted.brmsfit].
#' 
#' @rdname fitted_draws
#' @export
#' 
#' @seealso [brms::fitted.brmsfit()]
#' 
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid time-consuming model estimation, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check if the model fit object 'berkeley_exfit' exists
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Population average distance curve
#' fitted_draws(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' fitted_draws(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' fitted_draws(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' fitted_draws(model, deriv = 1, re_formula = NULL)
#' }
#' 
fitted_draws.bgmfit <-
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
           grid_add = NULL,
           ipts = NULL,
           deriv = 0,
           model_deriv = TRUE,
           summary = TRUE,
           robust = FALSE,
           transform_draws = NULL,
           scale = c("response", "linear"),
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
           xvar = NULL,
           difx = NULL,
           idvar = NULL,
           itransform = NULL,
           newdata_fixed = NULL,
           envir = NULL,
           ...) {
    
    
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- envir
    }
    
    assign_function_to_environment(transform_draws, 'transform_draws', 
                                   envir = NULL)
    
    model$model_info[['transform_draws']] <- transform_draws
    

    if(is.null(dpar)) {
      dpar <- "mu"
    }
    
    model <- getmodel_info(model = model, 
                           dpar = dpar, 
                           resp = resp, 
                           deriv = NULL, 
                           verbose = verbose)
    
   
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
    
    if(is.null(model_deriv)) {
      model_deriv <- TRUE
    }
    
    if (is.null(idata_method)) {
      idata_method <- 'm2'
    }
   
    
    rlang_trace_back <- rlang::trace_back()
    
    check_trace_back.bgmfit <- grepl("plot_conditional_effects", 
                                     rlang_trace_back[[1]])
    if(all(!check_trace_back.bgmfit)) {
      plot_conditional_effects_calling <- FALSE
    } else {
      plot_conditional_effects_calling <- TRUE 
    }
    
    check_trace_back.bgmfit <- grepl("growthparameters", 
                                     rlang_trace_back[[1]])
    if(all(!check_trace_back.bgmfit)) {
      growthparameters_calling <- FALSE
    } else {
      growthparameters_calling <- TRUE 
    }
    
    if(any(grepl("modelbased_growthparameters", rlang_trace_back[[1]]))) {
      growthparameters_calling <- FALSE 
    }
    
    if(any(grepl("plot_curves", rlang_trace_back[[1]]))) {
      growthparameters_calling <- FALSE 
    }
    
    
    ########################################################
    
    # For sigma
    if (deriv > 0) {
      need_velocity_curve <- TRUE
    } else {
      need_velocity_curve <- FALSE
    }
    
    if(deriv == 0) {
      only_distance_curve <- TRUE
    }
    
    if(need_velocity_curve) {
      need_xvar_must <- TRUE
    } else {
      need_xvar_must <- FALSE
    }
    
    ########################################################
    
    indirectcall <- FALSE
    if(!plot_conditional_effects_calling) {
      if(!is.null(model$xcall)) {
        arguments <- get_args_(as.list(match.call())[-1], model$xcall)
        # This when call coming from 'plot_curves()' or 'growthparameters()'
        if(dpar == "sigma") {
          if(!is.null(model$model_info[["which_sigma_model"]])) {
            if(model$model_info[['which_sigma_model']] == "basic") {
              set_sanitize <- FALSE
            } else if(model$model_info[['which_sigma_model']] != "basic") {
              set_sanitize <- FALSE
            }
          } # if(!is.null(model$model_info[["which_sigma_model"]])) {
        } else if(dpar != "sigma") {
          set_sanitize <- FALSE
        } # if(dpar == "sigma") { else...
        
         
        full.args <- evaluate_call_args(cargs = arguments, 
                                        fargs = NULL, 
                                        dargs = NULL, 
                                        sanitize_CustomDoCall_args = set_sanitize,
                                        check_formalArgs  = NULL,
                                        check_formalArgs_exceptions  = NULL,
                                        check_trace_back  = NULL,
                                        envir = parent.frame(),
                                        verbose = verbose)
        full.args$object <- full.args$model
        newdata          <- full.args$newdata
        indirectcall     <- TRUE
      } else {
        full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                        fargs = formals(), 
                                        dargs = list(...), 
                                        sanitize_CustomDoCall_args = FALSE,
                                        check_formalArgs  = NULL,
                                        check_formalArgs_exceptions  = NULL,
                                        check_trace_back  = NULL,
                                        envir = parent.frame(),
                                        verbose = verbose)
        full.args$model <- model
        
      }
    }
    
    xcall_str <- NULL
    if(plot_conditional_effects_calling) {
      full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                      fargs = formals(), 
                                      dargs = list(...), 
                                      sanitize_CustomDoCall_args = TRUE,
                                      check_formalArgs = NULL,
                                      check_formalArgs_exceptions = NULL,
                                      check_trace_back = NULL,
                                      envir = parent.frame(),
                                      verbose = verbose)
      xcall_str           <- full.args$xcall_str
      full.args$xcall_str <- NULL
    }
    
    
    expose_method_set <- model$model_info[['expose_method']]
    
    model$model_info[['expose_method']] <- 'NA' # Over ride method 'R'
    
    # 6.03.2025
    if(is.null(xcall_str)) {
      setxcall_   <- match.call()
    } else {
      setxcall_ <- xcall_str
    }
    
    post_processing_checks_args <- list()
    post_processing_checks_args[['model']]    <- model
    post_processing_checks_args[['xcall']]    <- setxcall_ 
    post_processing_checks_args[['resp']]     <- resp
    post_processing_checks_args[['envir']]    <- envir
    post_processing_checks_args[['deriv']]    <- deriv
    post_processing_checks_args[['all']]      <- FALSE
    post_processing_checks_args[['verbose']]  <- verbose
    post_processing_checks_args[['check_d0']] <- FALSE
    post_processing_checks_args[['check_d1']] <- TRUE
    post_processing_checks_args[['check_d2']] <- FALSE
    
    o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)

    post_processing_checks_args[['all']]      <- TRUE
    oall <- CustomDoCall(post_processing_checks, post_processing_checks_args)
    post_processing_checks_args[['all']]      <- FALSE
    
    if(!is.null(funlist)) {
      if(!is.list(funlist)) {
        stop("funlist must be a list")
      } else {
        o <- funlist
      }
    }
    
    # 6.03.2025
    # see slopes will be mandatory
    check_fun <- FALSE
    if(deriv > 0) {
      available_d1 <- o[['available_d1']]
      if(!available_d1) {
        model_deriv <- FALSE
        call_slopes <- TRUE
        post_processing_checks_args[['deriv']]    <- 0
        o <- CustomDoCall(post_processing_checks, post_processing_checks_args)
      }
      check_fun <- TRUE
    }
    
    
    ############################################
    # if sigma_model == "basic" and all function can be set with deriv 0 -> 1
    # out[['sigma_model_is_ba_set_d0_as_d1_val']] is zero
    if(!is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) {
      if(o[['sigma_model_is_ba_set_d0_as_d1']]) {
        deriv <- o[['sigma_model_is_ba_set_d0_as_d1_val']]
        sigma_model_is_ba_set_d0_as_d1_funs <- 
          o[['sigma_model_is_ba_set_d0_as_d1_funs']]
        for (i in names(sigma_model_is_ba_set_d0_as_d1_funs)) {
          # model$model_info$exefuns[[i]] <- NULL
          # assign(i, sigma_model_is_ba_set_d0_as_d1_funs[[i]], envir = envir)
          model$model_info$exefuns[[i]] <- 
            sigma_model_is_ba_set_d0_as_d1_funs[[i]]
        }
        check_fun <- FALSE
      } # o[['sigma_model_is_ba_set_d0_as_d1']]
    } # if(!is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) {
    

    if(dpar == "sigma") {
      if(deriv > 0) {
        if(!is.null(o[['sigma_model']])) {
          if(o[['sigma_model']] == "ls") {
            # nothing
          } else if(o[['sigma_model']] != "ls") {
            if(!is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) {
              if(!o[['sigma_model_is_ba_set_d0_as_d1']]) {
                check_fun    <- TRUE
                available_d1 <- FALSE
                model_deriv  <- FALSE
                call_slopes  <- TRUE
              }
            } else if(is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) {
              check_fun    <- TRUE
              available_d1 <- FALSE
              model_deriv  <- FALSE
              call_slopes  <- TRUE
            } # if(!is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) else if(is.null
          } # if(o[['sigma_model']] == "ls") { else if(o[['sigma_model']] ...
        } # if(!is.null(o[['sigma_model']])) {
      } # if(deriv > 0) {
    } # if(dpar == "sigma") {
    
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall,
                      usesavedfuns = usesavedfuns,
                      deriv = post_processing_checks_args[['deriv']],
                      envir = envir,
                      model_deriv = model_deriv,
                      ...)

    if(is.null(test)) {
      return(invisible(NULL))
    }
    
    
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
    
    misc <- c("verbose", "usesavedfuns", "clearenvfuns", "envir", "fullframe")
    
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

    
    if(growthparameters_calling) {
      if(!is.null(calling.args$re_formula_opt)) {
        if(calling.args$re_formula_opt == "V") {
          calling.args$re_formula <- NULL
        } else if(calling.args$re_formula_opt == "v") {
          calling.args$re_formula <- NA
        }
      } # if(!is.null(calling.args$re_formula_opt)) {
    } # if(growthparameters_calling) {
    
    
    
    # 6.03.2025
    if(is.null(newdata) & !indirectcall ) {
      calling.args_newdata         <- calling.args
      calling.args_newdata$model   <- calling.args_newdata$object
      calling.args_newdata$newdata <- model$model_info$bgmfit.data 
      calling.args_newdata$dpar    <- dpar
      get.newdata_args <- list()
      for (i in methods::formalArgs(get.newdata)) {
        get.newdata_args[[i]] <- calling.args_newdata[[i]]
      }
      get.newdata_args$ipts <- calling.args$ipts <- 
        set_for_check_ipts(ipts = ipts, nipts = 50, 
                           dpar = dpar, verbose = verbose)
      newdata <- CustomDoCall(get.newdata, get.newdata_args)
      rm('calling.args_newdata')
      rm('get.newdata_args')
      calling.args$newdata <- newdata
    }
    
    
    # 6.03.2025
    if(!indirectcall & !growthparameters_calling) {
      if(is.null(newdata)) {
        calling.args_newdata         <- calling.args
        calling.args_newdata$dpar    <- dpar
        get.newdata_args <- list()
        for (i in methods::formalArgs(get.newdata)) {
          get.newdata_args[[i]] <- calling.args_newdata[[i]]
        }
        newdata <- CustomDoCall(get.newdata, get.newdata_args)
        rm('calling.args_newdata')
        rm('get.newdata_args')
        calling.args$newdata <- newdata
      }
    }
    
    # 6.03.2025
    if(is.null(full.args$newdata)) {
      full.args$newdata <- calling.args$newdata
    }
    
    # set up mesage
    if(check_fun) {
      if(!available_d1) {
        message_for_model_deriv_FALSE <- ""
        message_for_model_deriv_FALSE <- 
          paste0(message_for_model_deriv_FALSE, "\n",
                 "calculating deriv by differentiation of distance curve")
        
        if(is.null(ipts)) {
          message_for_model_deriv_FALSE <- 
            paste0(message_for_model_deriv_FALSE, "\n",
                   "It is strongly recommended not to set'ipts = NULL'", "\n",
                   "Typically, ipts > 100 is needed to get smooth deriv curve")
        }
        if(any(grepl("fitted_draws", rlang_trace_back[[1]])) |
           any(grepl("predict_draws", rlang_trace_back[[1]]))) {
          message_for_model_deriv_FALSE <- 
            paste0(message_for_model_deriv_FALSE, "\n",
                   "A better approach would be use 'get_predictions()' ",
                   "instead")
        } # if(grepl("fitted_draws", rlang_trace_back[[1]]) |
      } # if(!available_d1) {
    } # if(check_fun) {
    
    
    calling.args <- 
      sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                 arguments = calling.args, 
                                 check_formalArgs = NULL,
                                 check_formalArgs_exceptions = c('object'),
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    calling.args$ipts <- ipts <- check_ipts(ipts = calling.args$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)
    
    if(!check_fun) {
      . <- CustomDoCall(fitted, calling.args)
    }
    
    
    if(!check_fun) {
      if(deriv == 0) {
        if(!is.null(calling.args$model$
                    model_info$transform_draws)) {
          . <- calling.args$model$
            model_info$transform_draws(.)
        } else if(!is.null(calling.args$object$
                           model_info$transform_draws)) {
          . <- calling.args$object$
            model_info$transform_draws(.)
        } else {
          if(verbose) message("transform_draws function is NULL, check it")
        }
      } # if(deriv = 0) {
    } # if(!check_fun) {
    
    if(!indirectcall) {
      if(dpar == "sigma") {
        sigma_model <- get_sigmamodel_info(model = model,
                                           newdata = newdata,
                                           dpar = dpar, 
                                           resp = resp, 
                                           what = 'model',
                                           cov = NULL, 
                                           all = FALSE, 
                                           verbose = verbose)
        
        calling.args$model$model_info[['which_sigma_model']] <- 
          model$model_info[['which_sigma_model']] <- sigma_model
        
        if(is.null(transform_draws)) {
          transform_draws <- 
            check_set_transform_draws_sigma(model = model, 
                                            dpar = dpar, 
                                            xvar = xvar, 
                                            resp = resp, 
                                            auto = TRUE,
                                            transform_draws = transform_draws,
                                            itransform = itransform,
                                            verbose = verbose)
          
          # Imp to assign calling.args[['transform_draws']] 
          calling.args[['transform_draws']] <- transform_draws
        }
        
        
        if(sigma_model == "basic") {
          if(!is.null(ipts)) {
            stop("For sigma_model = ",  
                 collapse_comma(sigma_model), ", the ipts should be NULL", 
                 "\n  ", 
                 "Currently, you have set this argument as ipts = ", ipts)
          }
        }
        
        msg_sigma_model_no_xvar <- 
          paste0("Although 'xvar' is strictly not required for estimating 
           distance curve when sigma_model = ",  collapse_comma(sigma_model), 
                 " but still it is better to specify 'xvar' to correctly label
           and plot x-axis. Otherwise x-axis wil be based on the xvar
           from the 'mu' part"
          )
        
        clean_msg_sigma_model_no_xvar <- trimws(gsub("\\s+", " ",
                                                     msg_sigma_model_no_xvar))
        
        
        if(sigma_model != "ls" && !need_xvar_must && !need_velocity_curve) {
          # if(sigma_model == "basic" && !need_velocity_curve) {
          if(is.null(xvar)) {
            if(verbose) {
              message(clean_msg_sigma_model_no_xvar)
            }
          }
        }
        
        if(sigma_model != "ls" && need_velocity_curve) {
          xvar <- check_set_xvar_sigma(model = model, 
                                       dpar = dpar, 
                                       xvar = xvar, 
                                       resp = resp, 
                                       auto = TRUE,
                                       verbose = verbose)
          
          newdata <- set_manual_datagrid(model = model,
                                         newdata = newdata,
                                         resp = resp, 
                                         dpar = NULL, 
                                         idvar = NULL,
                                         xvar = xvar,
                                         difx = difx,
                                         difx_asit = FALSE,
                                         auto = TRUE,
                                         xrange = NULL,
                                         length.out = NULL,
                                         grid_add = grid_add,
                                         grid_type= NULL,
                                         FUN = NULL,
                                         FUN_character = NULL,
                                         FUN_factor = NULL,
                                         FUN_logical = NULL,
                                         FUN_numeric = NULL,
                                         FUN_integer = NULL,
                                         FUN_binary = NULL,
                                         FUN_other = NULL,
                                         verbose = verbose)
          
          calling.args$model$model_info[['xvar_for_sigma_model_basic']] <- xvar
          calling.args[['xvar']] <- xvar
          calling.args[['difx']] <- difx
          calling.args$newdata <- newdata
        } # if(sigma_model == "basic") {
      } # if(dpar == "sigma") {
    } # if(!indirectcall) {
    
    
    if(indirectcall) {
      # plot_curves() and growthparameters()
      if(dpar == "sigma") {
        if(!is.null(model$model_info[["which_sigma_model"]])) {
          if(model$model_info[['which_sigma_model']] == "basic") {
            calling.args[['xvar']] <- 
              xvar <- model$model_info[['xvar_for_sigma_model_basic']]
          } else if(model$model_info[['which_sigma_model']] != "basic") {
            calling.args[['xvar']] <- xvar
          }
        } # if(!is.null(model$model_info[["which_sigma_model"]])) {
      } # if(dpar == "sigma") {
      
      if(dpar != "sigma") {
        calling.args[['xvar']] <- xvar
      } # if(dpar != "sigma") {
      
    } # if(indirectcall) {
  
    if(!is.null(attr(calling.args$newdata, 'difx'))) {
      calling.args[['difx']] <- difx <- attr(calling.args$newdata, 'difx')
    } 
    
    if(!plot_conditional_effects_calling) {
      if(check_fun) {
        if(deriv > 0) {
          if(available_d1) {
            . <- CustomDoCall(fitted, calling.args)
          }
          if(!available_d1) {
            if(verbose) {
              message(message_for_model_deriv_FALSE)
            }
            calling.args_mapderivqr_args <- calling.args
            calling.args_mapderivqr_args[['summary']] <- FALSE
            y0 <- CustomDoCall(fitted, calling.args_mapderivqr_args)
            if(!is.null(calling.args_mapderivqr_args$model$
                        model_info$transform_draws)) {
              y0 <- calling.args_mapderivqr_args$model$
                model_info$transform_draws(y0)
            } else if(!is.null(calling.args_mapderivqr_args$object$
                               model_info$transform_draws)) {
              y0 <- calling.args_mapderivqr_args$object$
                model_info$transform_draws(y0)
            } else {
              stop("transform_draws function is NULL, check it")
            }
            mapderivqr_args <- list()
            mapderivqr_args[['y0']] <- y0
            mapderivqr_args[['model']] <- calling.args[['object']]
            mapderivqr_args[['newdata']] <- calling.args[['newdata']]
            mapderivqr_args[['xvar']] <- calling.args[['xvar']]
            mapderivqr_args[['difx']] <- calling.args[['difx']]
            mapderivqr_args[['idvar']] <- calling.args[['idvar']]
            mapderivqr_args[['levels_id']] <- calling.args[['levels_id']]
            mapderivqr_args[['deriv']] <- calling.args[['deriv']]
            mapderivqr_args[['resp']] <- calling.args[['resp']]
            mapderivqr_args[['probs']] <- calling.args[['probs']]
            mapderivqr_args[['summary']] <- calling.args[['summary']]
            mapderivqr_args[['robust']] <- calling.args[['robust']]
            mapderivqr_args[['dpar']] <- calling.args[['dpar']]
            mapderivqr_args[['verbose']] <- calling.args[['verbose']]
            . <- CustomDoCall(mapderivqr, mapderivqr_args)
          }
        } # if(deriv > 0) {
      } # if(check_fun) {
    } # if(!plot_conditional_effects_calling) {
    
    full.args <-
      sanitize_CustomDoCall_args(what = "CustomDoCall",
                                 arguments = full.args,
                                 # check_formalArgs = plot_conditional_effects.bgmfit,
                                 check_formalArgs_exceptions = c('object'),
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
    full.args$object <- full.args$model
    
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    full.args$ipts <- ipts <- check_ipts(ipts = full.args$ipts, 
                                            nipts = NULL, 
                                            check_fun  = check_fun, 
                                            available_d1 = available_d1, 
                                            xcall = NULL, verbose = verbose)

    if(plot_conditional_effects_calling) {
      if(check_fun) {
        if(deriv > 0) {
          if(available_d1) {
            . <- CustomDoCall(fitted, full.args)
          }
          if(!available_d1) {
            if(verbose) {
              message(message_for_model_deriv_FALSE)
            }
            calling.args_mapderivqr_args <- full.args
            calling.args_mapderivqr_args[['summary']] <- FALSE
            y0 <- CustomDoCall(fitted, calling.args_mapderivqr_args)
            if(!is.null(calling.args_mapderivqr_args$model$
                        model_info$transform_draws)) {
              y0 <- calling.args_mapderivqr_args$model$
                model_info$transform_draws(y0)
            } else if(!is.null(calling.args_mapderivqr_args$object$
                               model_info$transform_draws)) {
              y0 <- calling.args_mapderivqr_args$object$
                model_info$transform_draws(y0)
            } else {
              stop("transform_draws function is NULL, check it")
            }
            mapderivqr_args <- list()
            mapderivqr_args[['y0']] <- y0
            mapderivqr_args[['model']] <- calling.args[['object']]
            mapderivqr_args[['newdata']] <- calling.args[['newdata']]
            mapderivqr_args[['xvar']] <- calling.args[['xvar']]
            mapderivqr_args[['difx']] <- calling.args[['difx']]
            mapderivqr_args[['idvar']] <- calling.args[['idvar']]
            mapderivqr_args[['levels_id']] <- calling.args[['levels_id']]
            mapderivqr_args[['deriv']] <- calling.args[['deriv']]
            mapderivqr_args[['resp']] <- calling.args[['resp']]
            mapderivqr_args[['probs']] <- calling.args[['probs']]
            mapderivqr_args[['summary']] <- calling.args[['summary']]
            mapderivqr_args[['robust']] <- calling.args[['robust']]
            mapderivqr_args[['dpar']] <- calling.args[['dpar']]
            mapderivqr_args[['verbose']] <- calling.args[['verbose']]
            . <- CustomDoCall(mapderivqr, mapderivqr_args)
          }
        } # if(deriv > 0) {
      } # if(check_fun) {
    } # if(plot_conditional_effects_calling) {

   
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
    
    
    full.args$idata_method <- idata_method
    full.args$fullframe <- eval(full.args$fullframe)
    
    if(!is.null(eval(full.args$fullframe))) {
      if(eval(full.args$fullframe)) {
        if (is.null(resp)) {
          resp_rev_ <- resp
        } else if (!is.null(resp)) {
          resp_rev_ <- paste0("_", resp)
        }
        xvar_ <- paste0('xvar', resp_rev_)
        yvar_ <- paste0('yvar', resp_rev_)
        groupvar_ <- paste0('groupvar', resp_rev_)
        if(is.null(xvar)) xvar <- model$model_info[[xvar_]]
        yvar <- model$model_info[[yvar_]]
        hierarchical_ <- paste0('hierarchical', resp_rev_)
        if(is.null(levels_id) & is.null(idvar)) {
          idvar <- model$model_info[[groupvar_]]
          if (!is.null(model$model_info[[hierarchical_]])) {
            idvar <- model$model_info[[hierarchical_]]
          }
        } else if (!is.null(levels_id)) {
          idvar <- levels_id
        } else if (!is.null(idvar)) {
          idvar <- idvar
        }
        xvar  <- xvar
        idvar <- idvar
        if(length(idvar) > 1) idvar <- idvar[1]
        yvar  <- 'yvar'
      } # if(eval(full.args$fullframe)) {
    } # if(!is.null(eval(full.args$fullframe))) {

    
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
      if (!is.na(model$model_info$univariate_by$by)) {
        if(full.args$idata_method == 'm1') setfullframe <- FALSE
        if(full.args$idata_method == 'm2') setfullframe <- TRUE
      } else {
        setfullframe <- FALSE
      }
    }
    
    if (!is.na(model$model_info$univariate_by$by)) {
      if(is.null(full.args$fullframe)) {
        full.args$fullframe <- fullframe <- FALSE
      }
      if(full.args$fullframe & full.args$idata_method == 'm1') {
        setfullframe <- FALSE
      }
      if(full.args$fullframe & full.args$idata_method == 'm2') {
        setfullframe <- TRUE
      }
      if(!full.args$fullframe) {
        setfullframe <- FALSE
      }
      if(setfullframe) {
        uvarby <- model$model_info$univariate_by$by
        uvarbyresp <- paste0(uvarby, resp)
        uvarbynewdata <- eval(full.args$newdata) %>% 
          dplyr::filter(!!dplyr::sym(uvarbyresp) == 1)
        . <- cbind(., uvarbynewdata)
        
        itransform_set <- get_itransform_call(itransform)
        if(any(itransform_set != "")) {
          . <- prepare_transformations(data = ., model = model, 
                                       itransform = itransform_set) 
        } # if(any(itransform_set != "")) {
      } # if(setfullframe) {
    } # if (!is.na(model$model_info$univariate_by$by)) {
    
    if (is.na(model$model_info$univariate_by$by)) {
      if(!is.null(eval(full.args$fullframe))) {
        if(eval(full.args$fullframe)) {
          cbindtonewdata <- eval(full.args$newdata)
          . <- cbind(cbindtonewdata, .)
        }
      }
      
      itransform_set <- get_itransform_call(itransform)
      cbindtonewdata <- eval(full.args$newdata)
      if(any(itransform_set != "")) {
       . <- prepare_transformations(data = ., model = model, 
                                    itransform = itransform_set) 
      } # if(any(itransform_set != "")) {
    } # if (is.na(model$model_info$univariate_by$by)) {
    . 
  } # end fitted_draws



#' @rdname fitted_draws
#' @export
fitted_draws <- function(model, ...) {
  UseMethod("fitted_draws")
}

