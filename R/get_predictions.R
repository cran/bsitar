

#' @title Estimate and plot growth curves for the Bayesian SITAR model
#' 
#' @description The \strong{get_predictions()} function estimates and plots
#'   growth curves (distance and velocity) by using the \pkg{marginaleffects}
#'   package as a back-end. This function computes growth curves
#'   ([marginaleffects::predictions()]), average growth curves
#'   ([marginaleffects::avg_predictions()]), and plots growth them via the
#'   [marginaleffects::plot_predictions()]. See
#'   [here](https://marginaleffects.com/) for details. Note that the
#'   \pkg{marginaleffects} package is highly flexible, and therefore, it is
#'   expected that the user has a strong understanding of its workings.
#'   Furthermore, since \pkg{marginaleffects} is rapidly evolving, the results
#'   obtained from the current implementation should be considered experimental.
#'   Also note that unlike [get_predictions()], the [fitted_draws()] and
#'   [predict_draws()] uses the \pkg{brms} package as a back-end.
#'
#' @details The \strong{get_predictions()} function estimates fitted values (via
#'   [brms::fitted.brmsfit()]) or the posterior draws from the posterior
#'   distribution (via [brms::predict.brmsfit()]) depending on the \code{type}
#'   argument.
#'   
#' @param average A logical indicating whether to internally call the
#'    [marginaleffects::predictions()] or [marginaleffects::avg_predictions()]
#'    function. If \code{FALSE} (default), [marginaleffects::predictions()] is
#'    called; otherwise, [marginaleffects::avg_predictions()] is used when
#'    \code{average = TRUE}.
#'
#' @param plot A logical specifying whether to plot predictions by calling
#'   [marginaleffects::plot_predictions()] (\code{TRUE}) or not (\code{FALSE}).
#'   If \code{FALSE} (default), [marginaleffects::predictions()] or
#'   [marginaleffects::avg_predictions()] are called to compute predictions (see
#'   \code{average} for details). Note that
#'   [marginaleffects::plot_predictions()] allows either \code{condition} or
#'   \code{by} arguments, but not both. Therefore, when the \code{condition}
#'   argument is not \code{NULL}, the \code{by} argument is set to \code{NULL}.
#'   This step is required because \strong{get_predictions()} automatically
#'   assigns the \code{by} argument when the model includes a covariate.
#'   
#' @param mapping_facet A named list that can be used to pass the aesthetic
#'   \code{mapping} and \code{facet} arguments to the \code{ggplot2} object
#'   (default \code{NULL}). The main use of this is to get overlay line plot
#'   instead of separate line plot for each \code{id} variable. This can also be
#'   used to remove a layer. For example, if user wants to remove the confidence
#'   intervals band from the lines i.e., \code{geom_ribbon}, then include
#'   \code{rm_geom_ribbon = TRUE}. Note the prefix \code{'rm_'}. This is a
#'   general approach in which any layer name with  prefix \code{'rm_'} included
#'   in the \code{mapping_facet} will be removed. An example is shown below:\cr
#'   \code{mapping_facet = list(group = 'id', colour = 'id', facet_wrap =
#'   c('study'), rm_geom_ribbon = TRUE)}.
#' 
#' @param deriv An integer to indicate whether to estimate the distance curve or
#'   its derivative (i.e., velocity curve). The \code{deriv = 0} (default) is
#'   for the distance curve, whereas \code{deriv = 1} is for the velocity curve.
#'
#' @param method A character string specifying the computation method: whether
#'   to use the \pkg{marginaleffects} machinery at the post-draw stage, i.e.,
#'   [marginaleffects::comparisons()] (\code{method = 'pkg'}) or to use custom
#'   functions for efficiency and speed (\code{method = 'custom'}, default).
#'   Note that \code{method = 'custom'} is particularly useful when testing
#'   hypotheses. Also, when \code{method = 'custom'},
#'   [marginaleffects::predictions()] is used internally instead of
#'   [marginaleffects::comparisons()].
#' 
#' @param constrats_by A character vector (default \code{NULL}) specifying the
#'   variable(s) by which hypotheses (at the post-draw stage) should be tested.
#'   Note that the variable(s) in \code{constrats_by} should be a subset of the
#'   variables included in the \code{'by'} argument.
#'   
#' @param constrats_at A character vector (default \code{NULL}) specifying the
#'   variable(s) at which hypotheses (at the post-draw stage) should be tested.
#'   \code{constrats_at} is particularly useful when the number of rows in the
#'   estimates is large because \pkg{marginaleffects} does not allow hypotheses
#'   testing when the number of rows exceeds 25.
#'  
#' @inheritParams growthparameters.bgmfit
#' @inheritParams get_growthparameters.bgmfit
#' @inheritParams get_comparisons.bgmfit
#' @inheritParams marginaleffects::predictions
#' @inheritParams marginaleffects::plot_predictions
#' @inheritParams brms::prepare_predictions
#' @inheritParams brms::fitted.brmsfit
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit()] 
#' function. Please see \code{brms::fitted.brmsfit()} for details on 
#' various options available.
#' 
#' @return An array of predicted mean response values. See
#'   [brms::fitted.brmsfit] for details.
#' 
#' @rdname get_predictions
#' @export
#' 
#' @seealso [marginaleffects::predictions()]
#'   [marginaleffects::avg_predictions()]
#'   [marginaleffects::plot_predictions()]
#' 
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
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
#' # Population average distance curve
#' get_predictions(model, deriv = 0, re_formula = NA, draw_ids = 1:2)
#' 
#' # Individual-specific distance curves
#' get_predictions(model, deriv = 0, re_formula = NULL, draw_ids = 1:2)
#' 
#' # Population average velocity curve
#' get_predictions(model, deriv = 1, re_formula = NA, draw_ids = 1:2)
#' 
#' # Individual-specific velocity curves
#' get_predictions(model, deriv = 1, re_formula = NULL, draw_ids = 1:2)
#' }
#' 
get_predictions.bgmfit <-
  function(model,
           resp = NULL,
           dpar = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           newdata = NULL,
           datagrid = NULL,
           re_formula = NA,
           newdata2 = NULL,
           allow_new_levels = FALSE,
           sample_new_levels = "gaussian", 
           parameter = NULL,
           xrange = 1,
           acg_velocity = 0.10,
           digits = 2,
           numeric_cov_at = NULL,
           aux_variables = NULL,
           grid_add = NULL,
           levels_id = NULL,
           avg_reffects = NULL,
           idata_method = NULL,
           ipts = NULL,
           seed = 123,
           future = FALSE,
           future_session = 'multisession',
           future_splits = TRUE,
           future_method = 'future',
           future_re_expose = NULL,
           usedtplyr = FALSE,
           usecollapse = TRUE,
           cores = NULL,
           fullframe = FALSE, 
           average = FALSE, 
           plot = FALSE, 
           mapping_facet = NULL,
           showlegends = NULL, 
           variables = NULL,
           condition = NULL,
           deriv = 0,
           model_deriv = TRUE,
           method = 'custom',
           marginals = NULL,
           pdrawso = FALSE, 
           pdrawsp = FALSE, 
           pdrawsh = FALSE, 
           type = NULL,
           by = NULL,
           conf_level = 0.95,
           transform = NULL,
           transform_draws = NULL,
           byfun = NULL,
           wts = NULL,
           hypothesis = NULL,
           equivalence = NULL,
           eps = NULL,
           constrats_by = NULL,
           constrats_at = NULL,
           reformat = NULL,
           estimate_center = NULL,
           estimate_interval = NULL,
           dummy_to_factor = NULL, 
           verbose = FALSE,
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
    ec_agg <- getOption("marginaleffects_posterior_center")
    ei_agg <- getOption("marginaleffects_posterior_interval")
    if(is.null(ec_agg)) ec_agg <- "mean"
    if(is.null(ei_agg)) ei_agg <- "eti"
    
    lean_ <- getOption("marginaleffects_lean")
    options("marginaleffects_lean" = FALSE)
    on.exit(options("marginaleffects_lean" = lean_), add = TRUE)
    
    insight::check_if_installed('cheapr', prompt = FALSE, stop = FALSE)
    
    try(zz <- insight::check_if_installed(c("marginaleffects"), 
                                          minimum_version = 
                                            get_package_minversion(
                                              'marginaleffects'
                                              ), 
                                          prompt = FALSE,
                                          stop = FALSE))
    
    
    if(!isTRUE(zz)) {
      message("Please install the latest version of the 'marginaleffects' 
              package",
              "\n ",
              "remotes::install_github('vincentarelbundock/marginaleffects')")
      return(invisible(NULL))
    }
    
    if(usedtplyr) {
      try(zz <- insight::check_if_installed(c("dtplyr"), 
                                            minimum_version =
                                              get_package_minversion(
                                                'dtplyr'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
      
      try(zz <- insight::check_if_installed(c("data.table"), 
                                            minimum_version =
                                              get_package_minversion(
                                                'data.table'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
    }
    
    callfuns <- TRUE
    setmarginals <- FALSE
    if(!is.null(marginals)) {
      setmarginals <- TRUE
      if(method == 'custom') callfuns <- FALSE
      if(method == 'pkg')    callfuns <- FALSE
    }
    
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- envir
    }
    
    if(!is.null(transform) & !is.null(transform_draws)) {
      stop("Please specify either transform or transform_draws, not both")
    }
    
    if(is.null(dpar)) {
      dpar <- "mu"
    }
    
    model <- getmodel_info(model = model, 
                           dpar = dpar, 
                           resp = resp, 
                           deriv = NULL, 
                           verbose = verbose)

    
    if(is.null(model$test_mode)) {
      model$test_mode <- FALSE
    }
    
    if(!model$test_mode) {
      unlock_replace_bind(package = "insight", what = "get_data",
                      replacement = custom_get_data.brmsfit, ept_str = T)
      if(verbose) {
        message(" As model[['test_mode']] = FLASE, the full data by the",
                "\n ", 
                "insight::get_data() is extracted via 'custom_get_data.brmsfit'",
                "\n ", 
                "This full data is needed for marginaleffects functions",
                "\n ", 
                "'To over ride this approach, set model[['test_mode']] = TRUE")
      }
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
    
    ndraws_org <- ndraws
    ndraws_exe <- FALSE
    if(!is.null(ndraws)) {
      ndraws_exe <- TRUE
    } else if(is.null(ndraws)) {
      ndraws <- brms::ndraws(model)
      ndraws_exe <- TRUE
    }
    
    if (is.null(resp)) {
      resp_rev_ <- resp
    } else if (!is.null(resp)) {
      resp_rev_ <- paste0("_", resp)
    }
    
    # For sigma
    xvar_      <- paste0('xvar', resp_rev_)
    sigmaxvar_ <- paste0('sigma', xvar_)
    cov_       <- paste0('cov', resp_rev_)
    sigmacov_  <- paste0('sigma', cov_)
    uvarby     <- model$model_info$univariate_by$by
    
    if(dpar == "mu") {
      if(is.null(xvar)) {
        xvar   <- model$model_info[[xvar_]]
      }
      cov    <- model$model_info[[cov_]]
    } else if(dpar == "sigma") {
      if(!is.na(model$model_info[[sigmaxvar_]])) {
        xvar   <- model$model_info[[sigmaxvar_]]
      } else if(is.na(model$model_info[[sigmaxvar_]]) & 
                !is.null(model$model_info[[xvar_]])) {
        xvar   <- model$model_info[[xvar_]]
      }
      cov    <- model$model_info[[sigmacov_]]
    }
    
    check_set_fun <- check_set_fun_transform(model = model, 
                                             which = 'ixfuntransform2',
                                             dpar = dpar, 
                                             resp= resp, 
                                             transform = itransform,
                                             auto = TRUE, 
                                             verbose = verbose)
    
    ifunx_ <- check_set_fun[['setfun']]
    if(check_set_fun[['was_null']]) {
      model$model_info[[check_set_fun[['setfunname']]]] <- ifunx_
    }
    # just added, 
    funx_ <- NULL

    if (is.null(idata_method)) {
      idata_method <- 'm2'
    }
    
    if(idata_method == 'm1') {
      stop("For marginaleffects based functions, the " ,
           " \n",
           " 'idata_method' argument must be either NULL or 'm2'" )
    }
    
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
    drawid <- NULL;
    
    conf <- conf_level
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    # set_names_  <- c('Estimate', 'Est.Error', probtitles)
    set_names_  <- c('Estimate', probtitles)
  
    expose_method_set <- model$model_info[['expose_method']]
    
    model$model_info[['expose_method']] <- 'NA' # Over ride 'R'
    
    setxcall_                                 <- match.call()
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
    
    # Imp to correctly call slopes to marginal effects
    # This will assign d0 and set call_slopes = TRUE
    if(method == "pkg") {
      if(deriv > 0) {
        model_deriv <- FALSE
        post_processing_checks_args[['deriv']] <- 0
      }
    }
    if(method == "custom") {
      if(deriv > 0) {
        if(!model_deriv) {
          post_processing_checks_args[['deriv']] <- 0
        }
      }
    }
    
    method_call <- 'predictions'
    comparison  <- NULL
    model_deriv <- model_deriv
    deriv       <- deriv
    
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
    
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall, 
                      usesavedfuns = usesavedfuns, 
                      deriv = deriv, envir = envir, 
                      model_deriv = model_deriv, 
                      ...)
    
    if(is.null(test)) return(invisible(NULL))
    
    ######################################################################
    ######################################################################
    
    # For sigma
    deriv.org       <- deriv
    model_deriv.org <- model_deriv
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
    
    # if plot, need xvar
    if(!plot) {
      need_xvar_must <- need_xvar_must
    } else {
      need_xvar_must <- TRUE
    }
    
    model$model_info[['difx']] <- difx

    if(dpar == "sigma") {
      sigma_model <- get_sigmamodel_info(model = model,
                                         newdata = newdata,
                                         dpar = dpar, 
                                         resp = resp, 
                                         what = 'model',
                                         cov = NULL, 
                                         all = FALSE, 
                                         verbose = verbose)
      
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
        
        model$model_info[['xvar_for_sigma_model_basic']] <- xvar
      } # if(sigma_model == "basic") {
    } # if(dpar == "sigma") {
    
    
    if(!is.null(transform)) {
      # new check added if(!is.function(transform)) {
      if(!is.function(transform)) {
        if(is.logical(transform)) {
          if(!transform) transform_draws <- 'identity'
        } else if(!is.logical(transform)) {
          if(transform == "exp") transform_draws <- 'exp'
          if(transform == "ln") transform_draws <- 'log'
        }
      } # if(!is.function(transform)) {
    } # if(!is.null(transform)) {
    
    assign_function_to_environment(transform_draws, 'transform_draws',
                                   envir = NULL)
    
    model$model_info[['transform_draws']] <- transform_draws
    
    
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
      if(deriv.org > 0) {
        if(!is.null(o[['sigma_model']])) {
          if(o[['sigma_model']] == "ls") {
            # nothing
          } else if(o[['sigma_model']] != "ls") {
            if(!is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) {
              if(!o[['sigma_model_is_ba_set_d0_as_d1']]) {
                check_fun    <- FALSE # TRUE
                available_d1 <- FALSE
                model_deriv  <- FALSE
                call_slopes  <- TRUE # FALSE # TRUE
              }
            } else if(is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) {
              check_fun    <- FALSE # TRUE
              available_d1 <- FALSE
              model_deriv  <- FALSE
              call_slopes  <- TRUE # FALSE # TRUE
            } # if(!is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) else if(is.null
          } # if(o[['sigma_model']] == "ls") { else if(o[['sigma_model']] ...
        } # if(!is.null(o[['sigma_model']])) {
      } # if(deriv > 0) {
    } # if(dpar == "sigma") {
    
    
    
    ########################################################
    
    
    # If default marginal effects 'dydx', then 
    call_predictions <- TRUE
    call_slopes      <- FALSE
    if(!model_deriv) {
      if(deriv > 0) {
        deriv <- 0
        call_predictions <- FALSE
        call_slopes      <- TRUE
      }
    } # if(!model_deriv) {
    
   
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
    post_processing_checks_args[['deriv']]    <- deriv
    
    
    ######################################################
    # somehow, condition, not by gives correct result result for slope 
    force_condition_and_by_switch_plot <- FALSE
    if(dpar == "sigma") {
      if(deriv.org > 0) {
        if(!is.null(o[['sigma_model']])) {
          if(o[['sigma_model']] == "ls") {
            # nothing
          } else if(o[['sigma_model']] != "ls") {
            force_condition_and_by_switch_plot <- TRUE
            if(is.null(difx)) {
              if(verbose) {
                message("The difx has been set same as variables i.e.," , 
                        "\n ",
                        collapse_comma(variables),
                        "\n ")
              }
              difx <- variables
            }
          }
        } # if(!is.null(o[['sigma_model']])) {
      } # if(deriv.org > 0) {
    } # if(dpar == "sigma") {
    
    
    if(force_condition_and_by_switch_plot) {
      if(is.null(variables) & is.null(difx)) {
        stop("For dpar = 'sigma', please specify 'variables' or 'difx' argument")
      }
    }
    
    if(force_condition_and_by_switch_plot) {
      return_plot      <- plot
      return_plot_est  <- FALSE
      plot             <- TRUE
    } else {
      return_plot      <- plot
    }
    
    ######################################################################
    ######################################################################
    
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
    
    
    if(!is.null(model$xcall)) {
      if(grepl("get_predictions", model$xcall)) {
        xcall <- "get_predictions"
      } else if(grepl("modelbased_growthparameters", model$xcall)) {
        xcall <- "modelbased_growthparameters"
      } else if(grepl("get_growthparameters", model$xcall)) {
        xcall <- "get_growthparameters"
      }
    } else {
      rlang_trace_back <- rlang::trace_back()
      check_trace_back.bgmfit <- grepl(".bgmfit", rlang_trace_back[[1]])
      if(all(!check_trace_back.bgmfit)) {
        # nothing
      } else {
        rlang_trace_back.bgmfit_i <- min(which(check_trace_back.bgmfit == TRUE))
        rlang_trace_back.bgmfit <- rlang_trace_back[[1]][[rlang_trace_back.bgmfit_i]]
        rlang_call_name <- rlang::call_name(rlang_trace_back.bgmfit)
        xcall <- rlang_call_name
      }
    }
    
    
    check_if_package_installed(model, xcall = xcall)
    
    model$xcall <- xcall
    
    
    call_from_modelbased_growthparameters <- FALSE
    call_from_modelbased_growthparameters_nonS3 <- FALSE
    if(xcall == "modelbased_growthparameters.bgmfit" |
       xcall == "modelbased_growthparameters") {
      call_from_modelbased_growthparameters <- TRUE
    }
    
  
    scallstatus <- sys.status()
    
    
    if(xcall == "modelbased_growthparameters_nonS3" |
       xcall == "CustomDoCall") {
      arguments <- get_args_(as.list(match.call())[-1], xcall, xclass = "", 
                             scallstatus = scallstatus)
      call_from_modelbased_growthparameters_nonS3 <- TRUE
    } else {
      arguments <- get_args_(as.list(match.call())[-1], xcall, xclass = NULL, 
                             scallstatus = scallstatus)
    }
    
    arguments$model        <- model
    arguments$usesavedfuns <- usesavedfuns
    
    get.cores_             <- get.cores(arguments$cores)
    
    # 28.09.2024
    if(is.null(get.cores_[['max.cores']])) {
      if(is.null(arguments$cores)) 
        get.cores_[['max.cores']] <- future::availableCores() - 1
    }
    
    arguments$cores <- setincores <-  get.cores_[['max.cores']]
    .cores_ps <- get.cores_[['.cores_ps']]
    
    get_future_args <- get_future_plan_args(future = future, 
                                            future_session = future_session, 
                                            oldfutureplan = future::plan(),
                                            setincores = setincores,
                                            verbose = FALSE)
    if(!is.null(get_future_args)) {
      future_plan_args <- get_future_args[['future_plan_args']]
      setplanis        <- get_future_args[['setplanis']]
      oldfutureplan    <- future::plan()
      do.call(future::plan, future_plan_args)
      on.exit(future::plan(oldfutureplan), add = TRUE)
      # marginaleffects future options
      getmarginaleffects_parallel <- 
        getOption("marginaleffects_parallel")
      getmarginaleffects_parallel_inferences <- 
        getOption("marginaleffects_parallel_inferences")
      options(marginaleffects_parallel = TRUE)
      options(marginaleffects_parallel_inferences = TRUE)
      on.exit(options("marginaleffects_parallel" = getmarginaleffects_parallel), 
              add = TRUE)
      on.exit(options("marginaleffects_parallel_inferences" = 
                        getmarginaleffects_parallel_inferences), 
              add = TRUE)
      # multicore
      if (inherits(future::plan(), "multicore")) {
        multthreadplan <- getOption("future.fork.multithreading.enable")
        options(future.fork.multithreading.enable = TRUE)
        on.exit(options("future.fork.multithreading.enable" = multthreadplan), 
                add = TRUE)
      }
    } else {
      
    } # if(!is.null(get_future_args)) { else {
    
    
   
    
    draw_ids_org <- draw_ids
    draw_ids_exe <- FALSE
    if(!is.null(draw_ids)) {
      draw_ids_exe <- TRUE
      ndraws_exe   <- FALSE
      draw_ids     <- draw_ids
    }
    
    
    
    future_splits_exe <- FALSE
    if(!is.null(future_splits)) {
      future_splits_exe <- TRUE
      if(is.logical(future_splits)) {
        if(future_splits) {
          if(ndraws_exe) {
            chunk_size_den <- setincores
            ndraws_seq <- sample.int(ndraws)
            chunk_size <- ndraws / chunk_size_den
            future_splits_at <- split(ndraws_seq, 
                                      ceiling(seq_along(ndraws_seq)/chunk_size))
            future_splits_at <- unname(future_splits_at)
          } else if(draw_ids_exe) {
            chunk_size_den <- setincores
            ndraws_seq <- draw_ids
            chunk_size <- length(draw_ids) / chunk_size_den
            future_splits_at <- split(ndraws_seq, 
                                      ceiling(seq_along(ndraws_seq)/chunk_size))
            future_splits_at <- unname(future_splits_at)
          }
        } else if(!future_splits) {
          future_splits_exe <- FALSE
          future_splits <- future_splits
        }
          
      } else if(is.list(future_splits)) {
        future_splits_at <- future_splits
        
      } else if(is.vector(future_splits)) {
        if(!is.numeric(future_splits)) {
          stop("future_splits must be a numeric vector of lenghth 2")
        } else if(length(future_splits) == 1) {
          if(draw_ids_exe) ndraws_exe <- FALSE
          if(ndraws_exe) {
            chunk_size_den <- future_splits
            ndraws_seq <- sample.int(ndraws)
            chunk_size <- ndraws / chunk_size_den
            future_splits_at <- split(ndraws_seq, 
                                      ceiling(seq_along(ndraws_seq)/chunk_size))
            future_splits_at <- unname(future_splits_at)
          } else if(draw_ids_exe) {
            chunk_size_den <- future_splits
            ndraws_seq <- draw_ids
            chunk_size <- length(draw_ids) / chunk_size_den
            future_splits_at <- split(ndraws_seq, 
                                      ceiling(seq_along(ndraws_seq)/chunk_size))
            future_splits_at <- unname(future_splits_at)
          } else if(is.null(ndraws_org)) {
            chunk_size_den <- future_splits
            ndraws_seq <- sample.int(ndraws)
            chunk_size <- ndraws / chunk_size_den
            future_splits_at <- split(ndraws_seq, 
                                      ceiling(seq_along(ndraws_seq)/chunk_size))
            future_splits_at <- unname(future_splits_at)
          } else if(!is.null(draw_ids_org)) {
            chunk_size_den <- future_splits
            ndraws_seq <- draw_ids
            chunk_size <- length(draw_ids) / chunk_size_den
            future_splits_at <- split(ndraws_seq, 
                                      ceiling(seq_along(ndraws_seq)/chunk_size))
            future_splits_at <- unname(future_splits_at)
          }
        } else if(length(future_splits) != 2) {
          stop("future_splits must be a numeric vector of lenghth 2")
        } else {
          if(future_splits[2] > future_splits[1]) {
            stop2c("The first element of 'future_splits' should equal to, or 
               greater than the second element")
          }
          future_splits_at <- parallel::splitIndices(future_splits[1], 
                                                     future_splits[2])
        }
      }
    }
   
      
    if(future_splits_exe) {
      if(plot) {
        future_splits_exe <- FALSE
        future_splits     <- NULL
        future_splits_at  <- NULL
        if(verbose) {
          message2c("future_splits can not be used when plot = TRUE. 
                  future_splits set as FALSE")
        }
      } # if(plot) {
      if(method == 'pkg') {
        future_splits_exe <- FALSE
        future_splits     <- NULL
        future_splits_at  <- NULL
        if(verbose) {
          message2c("future_splits can not be used when method = 'pkg'.
               future_splits set as FALSE")
        }
      } # if(method == 'pkg') {
    }
      
      
    
    if(!future_splits_exe) {
      future_splits_exe_future <- FALSE
      future_splits_exe_dofuture <- FALSE
    } else if(future_splits_exe) {
      if(future_method == 'future') {
        future_splits_exe_future <- TRUE
        future_splits_exe_dofuture <- FALSE
      }
      if(future_method == 'dofuture') {
        future_splits_exe_future <- FALSE
        future_splits_exe_dofuture <- TRUE
      }
    }
      
    
    
    re_expose <- FALSE
    if (future) {
      need_future_re_expose_cpp <- FALSE
      if(any(grepl("pstream__",
                   deparse(model$model_info$exefuns[[1]])))) {
        need_future_re_expose_cpp <- TRUE
      }
      
      
      if(is.null(future_re_expose)) {
        if(setplanis == "multisession") {
          if(need_future_re_expose_cpp) {
            re_expose <- TRUE
            if(verbose) {
              message2c("For multisession plan, argument 'future_re_expose' 
                        has been set as TRUE")
            }
          } else if(!need_future_re_expose_cpp) {
            if(verbose) {
              message2c("To speed up the calulations, it is advised to 
                        set 'future_re_expose' = TRUE")
            }
          }
        }
      } else if(!is.null(future_re_expose)) {
        if(future_re_expose) {
          re_expose <- TRUE
        } else if(!future_re_expose) {
          if(!need_future_re_expose_cpp) {
            # if(expose_method_set == "R") {
            if(verbose) {
              message2c("To speed up the calulations, it is advised to 
                        set 'future_re_expose' = TRUE")
            }
          } 
          if(need_future_re_expose_cpp & setplanis == "multisession") {
            # if(expose_method_set != "R") {
            stop2c("For plan multisession, the functions need to be 
                   re_exposed by setting 'future_re_expose' = TRUE")
          }
        }
      }
    } # if (future) {
    
    
   
    if (!future) {
      future_splits_at <- NULL
      future_splits_exe <- FALSE
      future_splits_exe_future <- FALSE
      future_splits_exe_dofuture <- FALSE
    }
      
    
    full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                    # fargs = formals(), 
                                    fargs = arguments, 
                                    dargs = list(...), 
                                    verbose = verbose)
    
    full.args$model       <- model
    full.args$model_deriv <- model_deriv
    full.args$newdata     <- newdata
    
    
    if(!is.null(full.args$hypothesis)) {
      if(method == 'pkg') {
        if(!is.null(full.args$by)) {
          if(is.logical(full.args$by)) {
            # stop("Argument 'by' is required for hypothesis")
          }
        }
      } else if(method == 'custom') {
        if(!is.null(full.args$by)) {
          if(is.logical(full.args$by)) {
            # stop("Argument 'by' is required for hypothesis")
          }
        } else if(is.null(full.args$by)) {
          # stop("Argument 'by' is required for hypothesis")
        }
      }
    }
    
   
    
    valid_hypothesis <- c("pairwise", "reference", "sequential", 
                          "revpairwise", "revreference", "revsequential")
    
    new_valid_hypothesis <- c("number", 
                              "string such as 'a = b'",
                              "formula of the following forms: 
                            ~ rhs, 
                            ~ lhs ~ rhs, 
                            lhs ~ rhs | group")
    
    if(!is.null(full.args$hypothesis)) {
      if(method == 'custom') {
        if(is.language(full.args$hypothesis)) {
          # stop("Argument 'hypothesis' must be one of the following strings: ",
          #      collapse_comma(valid_hypothesis))
        }
      } else if(method == 'pkg') {
        if(is.character(full.args$hypothesis)) {
          stop2c("Argument 'hypothesis' must be one of the following form: ",
                 collapse_comma(new_valid_hypothesis))
        }
      }
    }
    
    
    
    
    
    # 6.03,2025 - > coming from get_dv
    if(call_from_modelbased_growthparameters |
       call_from_modelbased_growthparameters_nonS3) {
      full.args[['peak']] <- NULL
      full.args[['takeoff']] <- NULL
      full.args[['trough']] <- NULL
      full.args[['acgv']] <- NULL
    }
    
    
    full.args <- 
      sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                 arguments = full.args, 
                                 check_formalArgs = get_predictions.bgmfit,
                                 check_formalArgs_exceptions = NULL,
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
    
    full.args$dpar    <- dpar
    
    get.newdata_args <- list()
    for (i in methods::formalArgs(get.newdata)) {
      get.newdata_args[[i]] <- full.args[[i]]
    }
    
    
    get.newdata_args$ipts <- full.args$ipts <- ipts <- 
      set_for_check_ipts(ipts = ipts, nipts = 50, dpar = dpar, verbose = verbose)
    
    full.args$newdata <- newdata <- CustomDoCall(get.newdata, 
                                                 get.newdata_args)
  

    # Interpolation points
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    full.args$ipts <- ipts <- check_ipts(ipts = full.args$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)
    
    if(!is.na(uvarby)) {
      uvarby_ind <- paste0(uvarby, resp)
      varne <- paste0(uvarby, resp)
      if(usedtplyr) {
        newdata <- newdata %>% dtplyr::lazy_dt() %>% 
          dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
      } else if(!usedtplyr) {
        newdata <- newdata %>% 
          dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
      }
    }
    
    full.args$newdata <- newdata
    full.args[["..."]] <- NULL
    
    
    if(!is.null(full.args[['transform_draws']])) {
      full.args[['transform']] <- transform <- full.args[['transform_draws']]
      if(verbose) message("'transform' set based on 'transform_draws'")
    } else if(!is.null(transform_draws)) {
      full.args[['transform']] <- transform <- transform_draws
      if(verbose) message("'transform' set based on 'transform_draws'")
    } else {
      #####
    }
    full.args[['transform']] <- transform <- transform_draws

    
    ########################################################################
    
    full.args$comparison  <- comparison
    full.args$model_deriv <- model_deriv
    # correct_comparison <- TRUE
    # correct_method_call <- TRUE 
    
    allowed_method_call_args <- c("predictions")
    
    # full.args[['comparison']]  <- set_comparisons_method_call_arg[['comparison']]
    # full.args[['method_call']] <- set_comparisons_method_call_arg[['method_call']]
    
    full.args[['comparison']]  <-  comparison
    full.args[['method_call']]  <-  method_call
    
    ########################################################################
    

    predictions_arguments <- full.args
    

    # Drop that not required for marginaleffects::
    exclude_args <- as.character(quote(
      c(
        parameter,
        xrange,
        acg_velocity,
        digits,
        numeric_cov_at,
        acg_asymptote,
        levels_id,
        avg_reffects,
        deriv,
        model_deriv,
        aux_variables, 
        idata_method, 
        # 19.09.2024
        # condition, 
        reformat, 
        dummy_to_factor, 
        expose_function,
        ipts,
        seed,
        future,
        future_session,
        future_splits,
        future_method,
        future_re_expose,
        cores,
        verbose,
        usesavedfuns,
        clearenvfuns,
        envir, 
        fullframe,
        average,
        plot,
        showlegends,
        average,
        estimate_center,
        estimate_interval, 
        reformat,
        method,
        marginals,
        constrats_by,
        constrats_at,
        usedtplyr,
        usecollapse,
        pdrawso,
        pdrawsp,
        pdrawsh,
        funlist,
        itransform,
        newdata_fixed,
        transform_draws,
        xvar,
        difx,
        mapping_facet,
        method_call
      )
    ))[-1]
    
    
    if(method == "pkg") {
      if(!plot) {
      #  exclude_args <- c(exclude_args, 'condition')
        # if(!is.null(predictions_arguments$condition)) {
        #   predictions_arguments$condition <- NULL
        # }
      }
    }
    
    
    if(call_slopes) {
      exclude_args <- c(exclude_args, 'transform', 'byfun')
    }
    
    
    # for get_dv
    if(call_from_modelbased_growthparameters) {
      exclude_args <- c(exclude_args, 'parallel', 'pdraws', 'constrats_subset')
    }
    
    for (exclude_argsi in exclude_args) {
      predictions_arguments[[exclude_argsi]] <- NULL
    }
    
    
    
    set_variables <- setup_variables_var(model, 
                                         variables = variables, 
                                         xvar = xvar, 
                                         eps = eps, 
                                         method = method, 
                                         deriv = deriv, 
                                         model_deriv = model_deriv,
                                         call_predictions = call_predictions,
                                         call_slopes = call_slopes,
                                         difx = difx, 
                                         xcall = xcall,
                                         cov = cov, 
                                         dpar = dpar,
                                         xvar_strict = TRUE, 
                                         switch_plot = FALSE,
                                         verbose = verbose)
    
    
    predictions_arguments$variables <- set_variables
  
    set_group <- setup_by_var(model = model, 
                              by = by, 
                              cov = cov, 
                              xvar = xvar, 
                              dpar = dpar,
                              method = method,
                              plot = plot,
                              condition = condition,
                              deriv = deriv,
                              difx = difx,
                              xcall = xcall,
                              xvar_strict = TRUE,
                              switch_plot = force_condition_and_by_switch_plot,
                              verbose = verbose)
   
    predictions_arguments$by         <- set_group
    
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)

    if(is.null(showlegends)) {
      if(is.null(predictions_arguments$re_formula)) {
        showlegends <- FALSE
      } else {
        showlegends <- TRUE
      }
    }
    
    
    # Set up datagrid
    if(!is.null(datagrid)) {
      if(is.data.frame(datagrid)) {
        set_datagrid <- datagrid
        predictions_arguments$newdata <- set_datagrid
      } else if(is.list(datagrid)) {
        if(is.null(datagrid[['model']]))
          setmodel <- model
        else
          setmodel <- datagrid$model
        if (is.null(datagrid[['newdata']]))
          setnewdata <- newdata
        else
          setnewdata <- datagrid$newdata
        if (is.null(datagrid[['grid_type']]))
          setgrid_type <- "mean_or_mode"
        else
          setgrid_type <- datagrid$grid_type
        if (is.null(datagrid[['by']]))
          setgrid_by <- NULL
        else
          setgrid_by <- datagrid$by
        
        if (is.null(datagrid[['FUN']]))
          setgrid_FUN <- NULL
        else
          setgrid_FUN <- datagrid$FUN
        if (is.null(datagrid[['FUN_character']]))
          setgrid_FUN_character <- NULL
        else
          setgrid_FUN_character <- datagrid$FUN_character
        if (is.null(datagrid[['FUN_factor']]))
          setgrid_FUN_factor <- NULL
        else
          setgrid_FUN_factor <- datagrid$FUN_factor
        if (is.null(datagrid[['FUN_logical']]))
          setgrid_FUN_logical <- NULL
        else
          setgrid_FUN_logical <- datagrid$FUN_logical
        if (is.null(datagrid[['FUN_numeric']]))
          setgrid_FUN_numeric <- NULL
        else
          setgrid_FUN_numeric <- datagrid$FUN_numeric
        if (is.null(datagrid[['FUN_integer']]))
          setgrid_FUN_integer <- NULL
        else
          setgrid_FUN_integer <- datagrid$FUN_integer
        if (is.null(datagrid[['FUN_binary']]))
          setgrid_FUN_binary <- NULL
        else
          setgrid_FUN_binary <- datagrid$FUN_binary
        if (is.null(datagrid[['FUN_other']]))
          setgrid_FUN_other <- NULL
        else
          setgrid_FUN_other <- datagrid$FUN_other
        
        if(is.null(datagrid[[xvar]])) 
          setxvar <- newdata[[xvar]] 
        else 
          setxvar <- datagrid$newdata[[xvar]]
        
        datagrid_arguments <- list(model = setmodel,
                                   newdata = setnewdata,
                                   by = setgrid_by,
                                   grid_type = setgrid_type,
                                   FUN = setgrid_FUN,
                                   FUN_character = setgrid_FUN_character,
                                   FUN_factor = setgrid_FUN_factor,
                                   FUN_logical = setgrid_FUN_logical,
                                   FUN_numeric = setgrid_FUN_numeric,
                                   FUN_integer = setgrid_FUN_integer,
                                   FUN_binary = setgrid_FUN_binary,
                                   FUN_other = setgrid_FUN_other)
        
        datagrid_arguments[[xvar]] <- setxvar
        if(setgrid_type == "mean_or_mode") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- set_group
        } else if(setgrid_type == "balanced") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- NULL
          # correctly set predictions_arguments[['by']]  too 
          predictions_arguments[['by']] <- 
            setdiff(predictions_arguments[['by']], cov)
        }
        set_datagrid <- CustomDoCall(marginaleffects::datagrid, 
                                     datagrid_arguments)
        predictions_arguments$newdata <- set_datagrid
      } else {
        stop("datagrid should be a data frame or named list")
      }
    } else if(is.null(datagrid)) {
      predictions_arguments$newdata <- predictions_arguments$newdata
    }

    predictions_arguments[['datagrid']] <- NULL
    
    # Somehow draw_ids not passed correctly if not specified explicitly as arg
    get_draw_ids <- predictions_arguments[['draw_ids']]
    if(!is.null(get_draw_ids)) {
      if(any(check_is_numeric_like(get_draw_ids))) {
        if(is.character(get_draw_ids)) {
          get_draw_ids <- ept(get_draw_ids)
        }
      }
    }
    if(is.null(eval(get_draw_ids))) {
      set_draw_ids <- NULL
    } else if(is.numeric(eval(get_draw_ids))) {
      set_draw_ids <- get_draw_ids
    } else if(!eval(get_draw_ids)) {
      set_draw_ids <- NULL
    }
   predictions_arguments[['draw_ids']] <- set_draw_ids
   
   
   exclude_args_con_by <- exclude_args
   
  # new
   if(method == "custom") {
     if(plot) {
       if(!is.null(predictions_arguments[['condition']])) {
       #  predictions_arguments[['by']] <- NULL
       }
     } else if(!plot) {
      # predictions_arguments[['condition']] <- NULL
     }
   }
 
   
   ###############################
   if(method == 'pkg') {
     remove_cols_ <- c('tmp_idx', 'predicted_lo', 
                       'predicted_hi', 'predicted', 'rowid')
   } else if(method == 'custom') {
     remove_cols_ <- c('term', 'contrast', 'tmp_idx', 'predicted_lo', 
                       'predicted_hi', 'predicted', 'rowid')
   }
   
   ###############################
   if(dpar == "sigma") {
     if(deriv.org == 0) {
       if(!is.null(variables)) {
         variables <- NULL
       }
       if(!is.null(predictions_arguments[['variables']])) {
         predictions_arguments[['variables']] <- NULL
       }
     }
   }
   ###############################
   
   
   
   out_sf_hy <- NULL
   allowed_methods <- c('pkg', 'custom')
   if(!method %in% allowed_methods) 
     stop("Argument 'method' should be one of the following:",
          "\n ", 
          collapse_comma(allowed_methods)
     )
  
   
   
   # Even for "pkg", we need following
   return_plot_est_pkg <- FALSE
   if(method == "pkg") {
     if(!force_condition_and_by_switch_plot) {
       if(!plot) {
         if(!is.null(predictions_arguments$condition)) {
           if(is.null(predictions_arguments$by)) {
             return_plot_est_pkg <- TRUE
           }
         } else if(is.null(predictions_arguments$condition)) {
           # return_plot_est_pkg <- TRUE
         }
       }
     }
   }
  
   
   return_plot_est_custom <- FALSE
   if(method == "custom") {
     if(!force_condition_and_by_switch_plot) {
       if(!plot) {
         if(!is.null(predictions_arguments$condition)) {
           if(is.null(predictions_arguments$by)) {
             return_plot_est_custom <- TRUE
           } 
         } else if(is.null(predictions_arguments$condition)) {
           # return_plot_est_custom <- TRUE
         }
       }
     }
   }
   
  
  
   
   if(method == 'pkg') {
     # predictions_arguments$by <- NULL
     if(call_predictions) {
       if(!plot) {
         if(return_plot_est_pkg) {
           predictions_arguments[['draw']] <- FALSE
           out_sf <- CustomDoCall(marginaleffects::plot_predictions,
                                  predictions_arguments)
           return(out_sf)
         } # if(return_plot_est_pkg) {
         if(!average) {
           out_sf <- CustomDoCall(marginaleffects::predictions, 
                             predictions_arguments)
         } else if(average) {
           out_sf <- CustomDoCall(marginaleffects::avg_predictions, 
                             predictions_arguments)
         }
       } else if(plot) {
         if(force_condition_and_by_switch_plot) {
           # predictions_arguments <- predictions_arguments.org
           predictions_arguments <- 
             condition_by_switch_fun(arg = predictions_arguments, 
                                     condition = 'condition', by = 'by',
                                     rm_by = TRUE, rm_condition = FALSE,
                                     verbose = verbose)
           if(!return_plot) {
             predictions_arguments[['draw']] <- FALSE
             out_sf <- CustomDoCall(marginaleffects::plot_predictions,
                                    predictions_arguments)
             if(return_plot_est) return(out_sf)
           } # if(!return_plot) {
           if(return_plot) {
             # predictions_arguments[['draw']] <- FALSE
             outp <- CustomDoCall(marginaleffects::plot_predictions,
                                  predictions_arguments)
             
             outp <- edit_mapping_facet(outp = outp, 
                                        by = predictions_arguments$by,
                                        condition = predictions_arguments$condition,
                                        xcall = xcall,
                                        method = method,
                                        mapping_facet = mapping_facet,
                                        showlegends = showlegends,
                                        funx_ = funx_,
                                        ifunx_ = ifunx_,
                                        envir = envir,
                                        which_aes = NULL, 
                                        print = FALSE,
                                        verbose = verbose)
             
             return(outp)
           } # if(!return_plot) {
         } else if(!force_condition_and_by_switch_plot) {
           outp <- CustomDoCall(marginaleffects::plot_predictions,
                                predictions_arguments)
          
           outp <- edit_mapping_facet(outp = outp, 
                                      by = predictions_arguments$by,
                                      condition = predictions_arguments$condition,
                                      xcall = xcall,
                                      method = method,
                                      mapping_facet = mapping_facet,
                                      showlegends = showlegends,
                                      funx_ = funx_,
                                      ifunx_ = ifunx_,
                                      envir = envir,
                                      which_aes = NULL, 
                                      print = FALSE,
                                      verbose = verbose)
          
           return(outp)
         } # else if(!force_condition_and_by_switch_plot) {
       } # else if(plot) {
     } # if(call_predictions) {
     
    
     if(call_slopes) {
       if(!plot) {
         if(return_plot_est_pkg) {
           predictions_arguments[['draw']] <- FALSE
           out_sf <- CustomDoCall(marginaleffects::plot_slopes,
                                  predictions_arguments)
           return(out_sf)
         } # if(return_plot_est_pkg) {
         if(!average) {
           out_sf <- CustomDoCall(marginaleffects::slopes,
                                  predictions_arguments)
         } else if(average) {
           out_sf <- CustomDoCall(marginaleffects::avg_slopes, 
                                  predictions_arguments)
         }
       } else if(plot) {
         if(force_condition_and_by_switch_plot) {
           # predictions_arguments <- predictions_arguments.org
           predictions_arguments <- 
             condition_by_switch_fun(arg = predictions_arguments, 
                                     condition = 'condition', by = 'by',
                                     rm_by = TRUE, rm_condition = FALSE,
                                     verbose = verbose)
           if(!return_plot) {
             predictions_arguments[['draw']] <- FALSE
             out_sf <- CustomDoCall(marginaleffects::plot_slopes,
                                    predictions_arguments)
             if(return_plot_est) return(out_sf)
           } # if(!return_plot) {
           if(return_plot) {
             # predictions_arguments[['draw']] <- FALSE
             outp <- CustomDoCall(marginaleffects::plot_slopes,
                                    predictions_arguments)
             
             outp <- edit_mapping_facet(outp = outp, 
                                        by = predictions_arguments$by,
                                        condition = predictions_arguments$condition,
                                        xcall = xcall,
                                        method = method,
                                        mapping_facet = mapping_facet,
                                        showlegends = showlegends,
                                        funx_ = funx_,
                                        ifunx_ = ifunx_,
                                        envir = envir,
                                        which_aes = NULL, 
                                        print = FALSE,
                                        verbose = verbose)
             
             return(outp)
           } # if(!return_plot) {
         } else if(!force_condition_and_by_switch_plot) {
           outp <- CustomDoCall(marginaleffects::plot_slopes,
                             predictions_arguments)
           
           outp <- edit_mapping_facet(outp = outp, 
                                      by = predictions_arguments$by,
                                      condition = predictions_arguments$condition,
                                      xcall = xcall,
                                      method = method,
                                      mapping_facet = mapping_facet,
                                      showlegends = showlegends,
                                      funx_ = funx_,
                                      ifunx_ = ifunx_,
                                      envir = envir,
                                      which_aes = NULL, 
                                      print = FALSE,
                                      verbose = verbose)
           
           return(outp)
         } # else if(!force_condition_and_by_switch_plot) {
       } # else if(plot) {
     } # if(call_slopes) {
   } # if(method == 'pkg') {

   
   if(method == 'pkg') {
     cols_to_remove_pkg <- c("df")
     out_sf <- out_sf %>% dplyr::select(-dplyr::any_of(cols_to_remove_pkg))
   }
   
   
   
   
   
   ############################################################################
   # method == 'custom' -  
   # below code exact same for get_predictions & get_comparisons
   ############################################################################
   testthat_mode  <- FALSE
   testthat_mode <- testthat_mode
   
   pdrawsp_est <- NULL
   pdrawsh_est <- NULL
   
   if(method == 'custom') {
     if(grepl("get_comparisons", xcall)) {
       predictions_arguments    <- comparisons_arguments
     }
     by                       <- predictions_arguments[['by']] 
     hypothesis_method_custom <- predictions_arguments[['hypothesis']]
     
     custom_method_call       <- full.args[['method_call']]
     
     if(grepl("get_predictions", xcall)) {
       predictions_arguments[['hypothesis']] <- NULL # evaluated later
     }
     
     if(grepl("get_comparisons", xcall)) {
       if(is.null(hypothesis_method_custom)) {
         if(is.null(custom_method_call))  custom_method_call <- 'comparisons'
       } else {
         if(is.null(custom_method_call)) custom_method_call <- 'predictions'
       }
     }
     
      
     if(custom_method_call == 'predictions') {
       predictions_arguments[['cross']]      <- NULL
       predictions_arguments[['method']]     <- NULL
       predictions_arguments[['hypothesis']] <- NULL # evaluated later
     }
     
     if(custom_method_call == 'comparisons') {
       comparisons_arguments[['transform']]     <- NULL
     }
     
      
      if(custom_method_call == 'predictions') {
        if(future_splits_exe) {
          for (i in names(predictions_arguments)) {
            predictions_arguments[[i]] <- eval(predictions_arguments[[i]])
          }
        }
      }

      if(custom_method_call == 'comparisons') {
        if(future_splits_exe) {
          for (i in names(comparisons_arguments)) {
            comparisons_arguments[[i]] <- eval(comparisons_arguments[[i]])
          }
        }
      }
      
     
     
     
     if(plot) {
       if(!force_condition_and_by_switch_plot) {
         if(custom_method_call == 'predictions') {
           if(!is.null(predictions_arguments[['by']])) {
             if(is.logical(predictions_arguments[['by']])) {
               if(!predictions_arguments[['by']]) {
                 predictions_arguments[['by']] <- NULL
               }
             }
           }
         } else if(custom_method_call == 'comparisons') {
           comparisons_arguments <- comparisons_arguments
         }
       } else if(force_condition_and_by_switch_plot) {
         if(custom_method_call == 'predictions') {
           predictions_arguments <- predictions_arguments
         } else if(custom_method_call == 'comparisons') {
           comparisons_arguments <- comparisons_arguments
         }
       } # if(!force_condition_and_by_switch_plot) { else
       
       if(custom_method_call == 'predictions') {
         if(is.null(predictions_arguments[['by']]) & 
            is.null(predictions_arguments[['condition']])) {
           stop2c("One of the `condition` and `by` arguments must 
              be supplied. Currentlu both 'by' and 'both' are NULL")
         } # if(is.null(predictions_arguments[['by']]) & 
       }
       
       if(custom_method_call == 'comparisons') {
         if(is.null(comparisons_arguments[['by']]) & 
            is.null(comparisons_arguments[['condition']])) {
           stop2c("One of the `condition` and `by` arguments must 
              be supplied. Currentlu both 'by' and 'both' are NULL")
         } # if(is.null(comparisons_arguments[['by']]) & 
       }
       
     } # if(plot) {
     
     
     
     # Set up calls
     if(callfuns) {
       if(future_splits_exe_future | 
          future_splits_exe_dofuture |
          !future_splits_exe) {
         # call_predictions
         if(call_predictions) {
           if(!plot) {
             if(!return_plot_est_custom) {
               if(!check_fun) {
                 if(custom_method_call == 'predictions') {
                   if(!average) funcall <- marginaleffects::predictions
                   if( average) funcall <- marginaleffects::avg_predictions
                   funcallargs          <- predictions_arguments
                 } else if(custom_method_call == 'comparisons') {
                   if(!average) funcall <- marginaleffects::comparisons
                   if( average) funcall <- marginaleffects::avg_comparisons
                   funcallargs          <- comparisons_arguments
                 }
               } 
               if(check_fun) {
                 if(deriv > 0) {
                   if(available_d1) {
                     if(custom_method_call == 'predictions') {
                       if(!average) funcall <- marginaleffects::predictions
                       if( average) funcall <- marginaleffects::avg_predictions
                       funcallargs          <- predictions_arguments
                     } else if(custom_method_call == 'comparisons') {
                       if(!average) funcall <- marginaleffects::comparisons
                       if( average) funcall <- marginaleffects::avg_comparisons
                       funcallargs          <- comparisons_arguments
                     }
                   } 
                   if(!available_d1) {
                     if(custom_method_call == 'predictions') {
                       if(!average) funcall <- marginaleffects::slopes
                       if( average) funcall <- marginaleffects::avg_slopes
                       funcallargs          <- predictions_arguments
                     } else if(custom_method_call == 'comparisons') {
                       if(!average) funcall <- marginaleffects::comparisons
                       if( average) funcall <- marginaleffects::avg_comparisons
                       funcallargs          <- comparisons_arguments
                     }
                   }
                 } # if(deriv > 0) {
               } # if(check_fun) {
             } else if(return_plot_est_custom) {
               if(!check_fun) {
                 if(custom_method_call == 'predictions') {
                   if(!average) funcall <- marginaleffects::plot_predictions
                   if( average) funcall <- marginaleffects::plot_predictions
                   funcallargs          <- predictions_arguments
                 } else if(custom_method_call == 'comparisons') {
                   if(!average) funcall <- marginaleffects::plot_comparisons
                   if( average) funcall <- marginaleffects::plot_comparisons
                   funcallargs          <- comparisons_arguments
                 }
               } 
               if(check_fun) {
                 if(deriv > 0) {
                   if(available_d1) {
                     if(custom_method_call == 'predictions') {
                       if(!average) funcall <- marginaleffects::plot_predictions
                       if( average) funcall <- marginaleffects::plot_predictions
                       funcallargs          <- predictions_arguments
                     } else if(custom_method_call == 'comparisons') {
                       if(!average) funcall <- marginaleffects::plot_comparisons
                       if( average) funcall <- marginaleffects::plot_comparisons
                       funcallargs          <- comparisons_arguments
                     }
                   } 
                   if(!available_d1) {
                     if(custom_method_call == 'predictions') {
                       if(!average) funcall <- marginaleffects::plot_slopes
                       if( average) funcall <- marginaleffects::plot_slopes
                       funcallargs          <- predictions_arguments
                     } else if(custom_method_call == 'comparisons') {
                       if(!average) funcall <- marginaleffects::plot_comparisons
                       if( average) funcall <- marginaleffects::plot_comparisons
                       funcallargs          <- comparisons_arguments
                     }
                   }
                 } # if(deriv > 0) {
               } # if(check_fun) {
             } # if(!return_plot_est_custom) { else if(return_plot_est_custom) {
           } else if(plot) {
             if(force_condition_and_by_switch_plot & !return_plot) {
               if(custom_method_call == 'predictions') {
                 if(!average) funcall <- marginaleffects::plot_predictions
                 if( average) funcall <- marginaleffects::plot_predictions
                 funcallargs          <- predictions_arguments
               } else if(custom_method_call == 'comparisons') {
                 if(!average) funcall <- marginaleffects::plot_comparisons
                 if( average) funcall <- marginaleffects::plot_comparisons
                 funcallargs          <- comparisons_arguments
               }
               funcallargs[['draw']] <- FALSE
             } else if(!force_condition_and_by_switch_plot) {
               if(custom_method_call == 'predictions') {
                 if(!average) funcall <- marginaleffects::plot_predictions
                 if( average) funcall <- marginaleffects::plot_predictions
                 funcallargs          <- predictions_arguments
               } else if(custom_method_call == 'comparisons') {
                 if(!average) funcall <- marginaleffects::plot_comparisons
                 if( average) funcall <- marginaleffects::plot_comparisons
                 funcallargs          <- comparisons_arguments
               }
             } # else if(!force_condition_and_by_switch_plot) {
           } # if(!plot) { else if(plot) {
         } # if(call_predictions) {
         
         
         
         # call_slopes
         if(call_slopes) {
           if(!plot) {
             if(!return_plot_est_custom) {
               if(custom_method_call == 'predictions') {
                 if(!average) funcall <- marginaleffects::slopes
                 if( average) funcall <- marginaleffects::avg_slopes
                 funcallargs          <- predictions_arguments
               } else if(custom_method_call == 'comparisons') {
                 if(!average) funcall <- marginaleffects::comparisons
                 if( average) funcall <- marginaleffects::avg_comparisons
                 funcallargs          <- comparisons_arguments
               }
             } else if(return_plot_est_custom) {
               if(custom_method_call == 'predictions') {
                 if(!average) funcall <- marginaleffects::plot_slopes
                 if( average) funcall <- marginaleffects::plot_slopes
                 funcallargs          <- predictions_arguments
               } else if(custom_method_call == 'comparisons') {
                 if(!average) funcall <- marginaleffects::plot_comparisons
                 if( average) funcall <- marginaleffects::plot_comparisons
                 funcallargs          <- comparisons_arguments
               }
             } # if(!return_plot_est_custom) {else if(return_plot_est_custom) {
           } else if(plot) {
             if(force_condition_and_by_switch_plot & !return_plot) {
               if(custom_method_call == 'predictions') {
                 if(!average) funcall <- marginaleffects::plot_slopes
                 if( average) funcall <- marginaleffects::plot_slopes
                 funcallargs          <- predictions_arguments
               } else if(custom_method_call == 'comparisons') {
                 if(!average) funcall <- marginaleffects::plot_comparisons
                 if( average) funcall <- marginaleffects::plot_comparisons
                 funcallargs          <- comparisons_arguments
               }
               funcallargs[['draw']] <- FALSE
             } else if(!force_condition_and_by_switch_plot) {
               if(custom_method_call == 'predictions') {
                 if(!average) funcall <- marginaleffects::plot_slopes
                 if( average) funcall <- marginaleffects::plot_slopes
                 funcallargs          <- predictions_arguments
               } else if(custom_method_call == 'comparisons') {
                 if(!average) funcall <- marginaleffects::plot_comparisons
                 if( average) funcall <- marginaleffects::plot_comparisons
                 funcallargs          <- comparisons_arguments
               }
             } # else if(!force_condition_and_by_switch_plot) {
           } # if(!plot) { else if(plot) {
         } # if(call_slopes) {
         
       } # if(future_splits_exe_future | future_splits_exe_dofuture) {
     } # if(callfuns) {

    
    
     # 6.03.2025
     if(!future_splits_exe & callfuns) {
       if(!plot) {
         if(return_plot_est_custom) {
           funcallargs[['draw']] <- FALSE
           out_sf <- CustomDoCall(funcall, funcallargs)
           out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]
           attr(out_sf$estimate, "label") <- NULL
           return(out_sf)
         } # if(return_plot_est_custom) {
         out <-  CustomDoCall(funcall, funcallargs)
       } else if(plot) {
         if(force_condition_and_by_switch_plot & !return_plot) {
           out <-  CustomDoCall(funcall, funcallargs)
         } else if(!force_condition_and_by_switch_plot) {
           outp <- CustomDoCall(funcall, funcallargs)
           outp <- edit_mapping_facet(outp = outp, 
                                      by = predictions_arguments$by,
                                      condition = predictions_arguments$condition,
                                      xcall = xcall,
                                      method = method,
                                      mapping_facet = mapping_facet,
                                      showlegends = showlegends,
                                      funx_ = funx_,
                                      ifunx_ = ifunx_,
                                      envir = envir,
                                      which_aes = NULL, 
                                      print = FALSE,
                                      verbose = verbose)
           return(outp)
         } # else if(!force_condition_and_by_switch_plot) {
       } # if(!plot) { else if(plot) {
     } # if(!future_splits_exe) {
    
     if(future_splits_exe_future & callfuns) {
       myzfun <- function(x, funcall, funcallargs) {
         funcallargs[['draw_ids']] <- x
         funcallargs[['ndraws']] <- NULL
         `%>%` <- bsitar::`%>%`
         if(re_expose) {
           if(verbose) {
             message2c("need to expose functions for 'multisession'")
           }
           funcallargs[['model']] <- 
             bsitar::expose_model_functions(funcallargs[['model']])
         }
         setenv <- funcallargs[['model']]$model_info$envir
         assign(o[[1]],
                funcallargs[['model']]$model_info[['exefuns']][[o[[2]]]], 
                envir = setenv)
         return(CustomDoCall(funcall, funcallargs))
       }
       if(!plot) {
         if(return_plot_est_custom | !return_plot_est_custom) {
           funcallargs[['draw']] <- FALSE
           out <- future.apply::future_lapply(future_splits_at, 
                                              FUN = myzfun,
                                              funcall = funcall,
                                              funcallargs = funcallargs,
                                              future.envir = parent.frame(),
                                              future.globals = TRUE,
                                              future.seed = TRUE)
         } # if(return_plot_est_custom) {
         if(!return_plot_est_custom) {
           # stop("Future not supported for this call...")
         }
         # stop("Future not supported for this call...")
         # out <-  future.apply::future_lapply(future_splits_at, 
         #                                     FUN = myzfun,
         #                                     funcall = funcall,
         #                                     funcallargs = funcallargs,
         #                                     future.envir = parent.frame(),
         #                                     future.globals = TRUE,
         #                                     future.seed = TRUE)
       } else if(plot) {
         if(force_condition_and_by_switch_plot & !return_plot) {
           out <-  future.apply::future_lapply(future_splits_at, 
                                               FUN = myzfun,
                                               funcall = funcall,
                                               funcallargs = funcallargs,
                                               future.envir = parent.frame(),
                                               future.globals = TRUE,
                                               future.seed = TRUE)
         } else if(!force_condition_and_by_switch_plot) {
           outp <- CustomDoCall(funcall, funcallargs)
           outp <- edit_mapping_facet(outp = outp, 
                                      by = predictions_arguments$by,
                                      condition = predictions_arguments$condition,
                                      xcall = xcall,
                                      method = method,
                                      mapping_facet = mapping_facet,
                                      showlegends = showlegends,
                                      funx_ = funx_,
                                      ifunx_ = ifunx_,
                                      envir = envir,
                                      which_aes = NULL, 
                                      print = FALSE,
                                      verbose = verbose)
           return(outp)
         } # else if(!force_condition_and_by_switch_plot) {
       } # if(!plot) { else if(plot) {
     } # if(future_splits_exe_future) {
     
     
     
     
     if(future_splits_exe_dofuture & callfuns) {
       `%doFuture_function%` <- doFuture::`%dofuture%`
       # somehow .options.future = list(seed = TRUE) not working, so set below
       dofutureplan <- getOption("doFuture.rng.onMisuse")
       options(doFuture.rng.onMisuse = "ignore")
       on.exit(options("doFuture.rng.onMisuse" = dofutureplan), add = TRUE)
       
       if(!plot) {
         if(return_plot_est_custom | !return_plot_est_custom) {
           funcallargs[['draw']] <- FALSE
           out <- foreach::foreach(x = 1:length(future_splits_at),
                                   .options.future = list(seed = TRUE),
                                   .options.future =
                                     list(globals = c('future_splits_at',
                                                      're_expose',
                                                      'o',
                                                      'call_predictions',
                                                      'call_slopes',
                                                      'CustomDoCall',
                                                      'funcallargs',
                                                      'funcall',
                                                      'verbose'))
           ) %doFuture_function% {
             x <- future_splits_at[[x]]
             funcallargs[['draw_ids']] <- x
             funcallargs[['ndraws']] <- NULL
             `%>%` <- bsitar::`%>%`
             if(re_expose) {
               if(verbose) {
                 message2c("need to expose functions for 'multisession'")
               }
               funcallargs[['model']] <- bsitar::expose_model_functions(funcallargs[['model']])
             }
             setenv <- funcallargs[['model']]$model_info$envir
             assign(o[[1]],
                    funcallargs[['model']]$model_info[['exefuns']][[o[[2]]]], 
                    envir = setenv)
             return(CustomDoCall(funcall, funcallargs))
           } # %doFuture_function% {
         } # if(return_plot_est_custom) {
         if(!return_plot_est_custom) {
           # stop("Future not supported for this call...")
         }
       } else if(plot) {
         if(force_condition_and_by_switch_plot & !return_plot) {
           out <- foreach::foreach(x = 1:length(future_splits_at),
                                   .options.future = list(seed = TRUE),
                                   .options.future =
                                     list(globals = c('future_splits_at',
                                                      're_expose',
                                                      'o',
                                                      'call_predictions',
                                                      'call_slopes',
                                                      'CustomDoCall',
                                                      'funcallargs',
                                                      'funcall',
                                                      'verbose'))
           ) %doFuture_function% {
             x <- future_splits_at[[x]]
             funcallargs[['draw_ids']] <- x
             funcallargs[['ndraws']] <- NULL
             `%>%` <- bsitar::`%>%`
             if(re_expose) {
               if(verbose) {
                 message2c("need to expose functions for 'multisession'")
               }
               funcallargs[['model']] <- bsitar::expose_model_functions(funcallargs[['model']])
             }
             setenv <- funcallargs[['model']]$model_info$envir
             assign(o[[1]],
                    funcallargs[['model']]$model_info[['exefuns']][[o[[2]]]], 
                    envir = setenv)
             return(CustomDoCall(funcall, funcallargs))
           } # %doFuture_function% {
         } else if(!force_condition_and_by_switch_plot) {
           # already not supported
           # outp <- CustomDoCall(funcall, funcallargs)
           # outp <- edit_mapping_facet(outp = outp, 
           #                            by = predictions_arguments$by,
           #                            condition = predictions_arguments$condition,
           #                            xcall = xcall,
           #                            method = method,
           #                            mapping_facet = mapping_facet,
           #                            showlegends = showlegends,
           #                            funx_ = funx_,
           #                            ifunx_ = ifunx_,
           #                            envir = envir,
           #                            which_aes = NULL, 
           #                            print = FALSE,
           #                            verbose = verbose)
           # return(outp)
         } # else if(!force_condition_and_by_switch_plot) {
       } # if(!plot) { else if(plot) {
     } # if(future_splits_exe_dofuture & callfuns) {
     
     
     
     
     
     
     
     
     # posterior_draws_function <- function(x, ...) {
     #   out[[x]] %>% 
     #     marginaleffects:: posterior_draws(shape = "long") %>% 
     #     dplyr::mutate(drawid = as.numeric(drawid)) %>% 
     #     dplyr::mutate(drawid = future_splits_at[[x]] [.data[['drawid']]]) %>% 
     #     dplyr::mutate(drawid = as.factor(drawid)) %>% 
     #     dplyr::relocate(drawid, .before = 'draw')
     # }
     # 
     # 
     # consecutive_drawid_function <- function(x, ...) {
     #   x %>% 
     #     dplyr::group_by(drawid) %>% 
     #     dplyr::mutate(drawid = dplyr::cur_group_id()) %>% 
     #     dplyr::mutate(drawid = as.factor(drawid)) %>% 
     #     dplyr::ungroup()
     # }
     # 
     # 
     # posterior_draws_dt <- function(i) {
     #   dt <- as.data.table(marginaleffects::posterior_draws(out[[i]], 
     #                                                        shape = "long"))
     #   dt[, drawid := as.factor(future_splits_at[[i]][as.numeric(drawid)])]
     #   data.table::setcolorder(dt, "drawid")
     #   return(dt)
     # }
     
     posterior_draws_collapse <- function(i) {
       dt <- collapse::qDT(marginaleffects::posterior_draws(out[[i]], 
                                                            shape = "long"))
       draw_idx <- as.numeric(dt$drawid)
       dt[, drawid := as.factor(future_splits_at[[i]][draw_idx])]
       return(dt)
     }
     
     
     # somehow this need consequence number
     if(!future_splits_exe) {
       if(callfuns) {
         if(pdrawso) return(out)
         onex0 <- out %>% marginaleffects::posterior_draws()
       }
     } else if(future_splits_exe) {
       if(callfuns) {
         if(pdrawso) {
           out <- out %>% CustomDoCall(rbind, .)
           return(out)
         }
         # dplyr
         # onex0 <- lapply(1:length(future_splits_at), 
         #                 FUN = posterior_draws_function)
         # onex0 <- onex0 %>% CustomDoCall(rbind, .)
         # onex0 <- consecutive_drawid_functionx(onex0)
         # dplyr data.table
         # onex0 <- data.table::rbindlist(lapply(seq_along(out), 
         #                                       posterior_draws_dt))
         # onex0[, drawid := as.factor(.GRP), by = drawid]
         # dplyr collapse
         onex0 <- collapse::unlist2d(future.apply::future_lapply(
           seq_along(out), posterior_draws_collapse), 
           idcols = FALSE, DT = TRUE)
         onex0[, drawid := as.factor(collapse::groupid(drawid))] 
       }
     }
     
     
     
     marginals_list_consecutive_drawid_function <- function(x, ...) {
       if(x == 1) {
         oux <- out[[x]]
         oux$drawid <- as.numeric(oux$drawid)
       } else {
         maxpre <- max(as.numeric(levels(out[[x-1]]$drawid)))
         oux <- out[[x]]
         oux$drawid <- as.numeric(oux$drawid) + maxpre * x
       }
       return(oux)
     }
     
     
     
     
     if(setmarginals) {
       if(inherits(marginals, 'list')) {
         onex0 <-
           {. <- lapply(1:length(marginals), marginals_list_consecutive_drawid_function)
           list2DF(lapply(setNames(seq_along(.[[1]]), names(.[[1]])), function(i)
             unlist(lapply(., `[[`, i), FALSE, FALSE)))}
         onex0$drawid <- cheapr::factor_(onex0$drawid)
       } else {
         onex0 <- marginals
       }
     }
     
     
     
     
     if(!isFALSE(pdrawsp)) {
       if(!is.character(pdrawsp)) pdrawsp <- "return"
       selectchoicesr <- c("return", 'add') 
       checkmate::assert_choice(pdrawsp, choices = selectchoicesr)
       if(pdrawsp == 'return') {
         return(onex0)
       } else if(pdrawsp == 'add') {
         pdrawsp_est <- onex0
       } else {
         # nothing
       }
     }
     
     
     
     # 16.10.2024
     # when marginals are given, then need to summarise
     if(setmarginals) {
       namesx <- c('estimate', 'conf.low', 'conf.high')
       setdrawidparm_at_ <- c(by, namesx)
       out_sf <- 
         onex0 %>%
         collapse::fgroup_by(by) %>%
         collapse::fsummarise(collapse::mctl(
           get_pe_ci_collapse(.data[['draw']], ec_agg = ec_agg, 
                              ei_agg = ei_agg, na.rm = TRUE, 
                              nthreads = arguments$cores, 
                              conf = conf, probs = probs))
         )  %>% collapse::frename(., setdrawidparm_at_)
     } # if(setmarginals) {
     
     
     
     if(!setmarginals) {
       setdrawidparm <- by
       # print(onex0)
       namesx <- c('estimate', 'conf.low', 'conf.high')
       if(!isFALSE(setdrawidparm)) setdrawidparm_ <- c(setdrawidparm, namesx)
       if( isFALSE(setdrawidparm)) setdrawidparm_ <- c(namesx)
       out_sf <- onex0 %>% 
         collapse::fsubset(., drawid == 1) %>% 
         collapse::fselect(., setdrawidparm_)
     }
     
     
     ######## new 
     if(is.null(hypothesis)) {
       constrats_by <- by
       if(is.null(constrats_by)) {
         constrats_by <- constrats_by
       } else if(!is.null(constrats_by)) {
         if(is.logical(constrats_by)) {
           if(!constrats_by) {
             constrats_by <- NULL
           }
         } else if(!is.logical(constrats_by)) {
           if(is.character(constrats_by)) {
             constrats_by <- constrats_by
           } else {
             stop2c("constrats_by must be character")
           }
         }
       } 
       constrats_by      <- constrats_by
       check_names_onex0 <- c('rowid', 'term', 'contrast')
       set_constrats_mfx <- c()
       for (i in names(onex0)) {
         if(i %in% check_names_onex0) {
           set_constrats_mfx <- c(set_constrats_mfx, i)
         }
       }
       # set_constrats_mfx  <- c('rowid', 'term', 'contrast')
       set_constrats_by  <- c(constrats_by, set_constrats_mfx)
       namesx            <- c('estimate', 'conf.low', 'conf.high')
       setdrawidparm_at_ <- c(set_constrats_by, namesx)
       
       out_sf <- onex0 %>%
         collapse::fgroup_by(set_constrats_by) %>%
         collapse::fsummarise(collapse::mctl(
           get_pe_ci_collapse(.data[['draw']], ec_agg = ec_agg, 
                              ei_agg = ei_agg, na.rm = TRUE, 
                              nthreads = arguments$cores, 
                              conf = conf, probs = probs))
         ) %>%
         collapse::frename(., setdrawidparm_at_) %>% 
         # collapse::roworder('term') %>% 
         # collapse::roworderv(set_constrats_by) %>% 
         setcolorder(., set_constrats_mfx)
       
       if("term" %in% set_constrats_by) {
         out_sf <- collapse::roworder(out_sf, 'term')
       }
     } # if(is.null(hypothesis)) {
     
     
     
     
     if(!is.null(hypothesis)) {
       # For hypothesis
       groupvarshyp1 <- c('drawid')
       groupvarshyp2 <- c('term')
       if(!is.null(constrats_at)) {
         if(!is.character(constrats_at)) 
           stop2c("'constrats_at' must be a character vector")
         for (caxi in constrats_at) {
           if(!caxi %in% names(onex0)) {
             stop2c("Variable '", caxi, ". specified in 'constrats_at' is not in",
                    "\n ", 
                    " the estimates. Note that ", caxi, " should also be included",
                    "\n ", 
                    " in the 'by' argument. The current 'by' argument includes:", 
                    "\n ",
                    collapse_comma(by)
             )
           }
           groupvarshyp1 <- c(caxi, groupvarshyp1)
           groupvarshyp2 <- c(caxi, groupvarshyp2)
         } # for (caxi in names(constrats_at)) {
       } # if(!is.null(constrats_at)) {
       
       # Not needed with new hypothesis 
       # if(is.null(constrats_by)) {
       #   stop2c("Please specify 'constrats_by' argument when testing 
       #          'hypothesis'",
       #          "\n ",
       #          " The available options are: ",
       #          collapse_comma(by))
       # }
       
       if(!is.null(constrats_by)) {
         if(!is.character(constrats_by)) 
           stop2c("'constrats_by' must be a character vector")
         for (caxi in constrats_by) {
           if(!caxi %in% names(onex0)) {
             stop2c("Variable '", caxi, ". specified in 'constrats_by' is not in",
                    "\n ", 
                    " the estimates. Note that ", caxi, " should also be included",
                    "\n ", 
                    " in the 'by' argument. The current 'by' argument includes:", 
                    "\n ",
                    collapse_comma(by)
             )
           }
         } # for (caxi in names(constrats_by)) {
       } # if(!is.null(constrats_by)) {
       
       
       
       
       set_constrats_by <- c(constrats_by, 'draw')
       namesx <- c('estimate', 'conf.low', 'conf.high')
       
       setdrawidparm_at <- c(constrats_at, 'term')
       setdrawidparm_at_ <- c(setdrawidparm_at, namesx)
       
       
       # temhyy <-
       #   onex0 %>% 
       #   collapse::fgroup_by(groupvarshyp1) %>%
       #   collapse::fselect(set_constrats_by) %>% 
       #   collapse::frename('estimate' = 'draw') %>% 
       #   collapse::fsummarise(collapse::qDF(
       #     get_hypothesis_x(.data,
       #                      by = constrats_by,
       #                      hypothesis = hypothesis,
       #                      draws = 'estimate'))) 
       # 
       # out_sf_hy <- 
       #   temhyy %>%
       #   collapse::fgroup_by(groupvarshyp2) %>%
       #   collapse::fsummarise(collapse::mctl(
       #     get_pe_ci_collapse(.data[['estimate']], ec_agg = ec_agg, 
       #                        ei_agg = ei_agg, na.rm = TRUE, 
       #                        nthreads = arguments$cores, 
       #                        conf = conf, probs = probs))
       #   ) %>%
       #   collapse::frename(., setdrawidparm_at_)
       
       
       
       # parameter is just a placeholder
       onex0 <- data.table::setDT(onex0)[, 'parameter' := "apgv"]
       if(!is.null(constrats_by)) {
         hypothesis_by_what <- constrats_by
       } else {
         hypothesis_by_what <- full.args$by
       }
       out_sf_hy <- get_comparison_hypothesis(data = onex0, 
                                              full.args = full.args, 
                                              by = hypothesis_by_what,
                                              evaluate_comparison = FALSE,
                                              evaluate_hypothesis = TRUE,
                                              rope_test = FALSE,
                                              pd_test = FALSE,
                                              get_range_null_form = FALSE,
                                              get_range_null_value = FALSE,
                                              format = FALSE,
                                              verbose = FALSE)
       out_sf_hy <- out_sf_hy[, 'parameter' := NULL]
       out_sf_hy <- data.table::setnames(out_sf_hy, "hypothesis", "term")
       
       # if(!isFALSE(pdrawsh)) {
       #   temhyy <-
       #     onex0 %>%
       #     collapse::fgroup_by(groupvarshyp1) %>%
       #     collapse::fselect(set_constrats_by) %>%
       #     collapse::frename('estimate' = 'draw') %>%
       #     collapse::fsummarise(collapse::qDF(
       #       get_hypothesis_x(.data,
       #                        by = constrats_by,
       #                        hypothesis = hypothesis,
       #                        draws = 'estimate')))
       #   selectchoicesr <- c("return", 'add') 
       #   checkmate::assert_choice(pdrawsh, choices = selectchoicesr)
       #   if(pdrawsh == 'return') {
       #     return(temhyy)
       #   } else if(pdrawsh == 'add') {
       #     pdrawsh_est <- temhyy
       #   } 
       # } # if(!isFALSE(pdrawsh)) {
       
     } # if(!is.null(hypothesis)) {
   } # if(method == 'custom') {
   
   
   
   out_sf <- out_sf %>% data.frame() %>% 
     dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                 ~ round(., digits = digits))) %>%
     data.frame()
   
   
   if(!is.null(pdrawsh_est)) {
     if(usecollapse) {
       pdrawsh_est <- pdrawsh_est %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits)))
     } else {
       pdrawsh_est <- pdrawsh_est %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits)))
     }
   }
   
   
   
   if(!is.null(pdrawsp_est)) {
     if(usecollapse) {
       pdrawsp_est <- pdrawsp_est %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits)))
     } else {
       pdrawsp_est <- pdrawsp_est %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits)))
     }
   }
   
   
   ##############################################################
   ##############################################################
   # prepare_data2
   newdata_before_itransform <- newdata
   itransform_set <- get_itransform_call(itransform = itransform,
                                         model = model, 
                                         newdata = newdata,
                                         dpar = dpar, 
                                         resp = resp,
                                         auto = TRUE,
                                         verbose = verbose)
   
   itransform_set_x_for_sigma_model <- c("varpower", 
                                         "varconstpower",
                                         "varexp", 
                                         "fitted",
                                         "fittedz",
                                         "fittedpower", 
                                         "fittedexp", 
                                         "mean", 
                                         "meanpower", 
                                         "meanexp", 
                                         "residual",
                                         "residualpower",
                                         "residualexp")
   
   if(!is.null(model$model_info[['which_sigma_model']])) {
     sigma_model <- model$model_info[['which_sigma_model']]
     if(sigma_model %in% itransform_set_x_for_sigma_model) {
       if(!is.null(itransform)) {
         itransform_set <- c(itransform_set, 'x')
       }
     }
   }
   
   
   if(any(itransform_set != "")) {
     if(!is.null(out_sf)) {
       out_sf <- prepare_transformations(data = out_sf, model = model,
                                         itransform = itransform_set)
     }
     if(!is.null(out_sf_hy)) {
       out_sf_hy <- prepare_transformations(data = out_sf_hy, model = model,
                                            itransform = itransform_set)
     }
     if(!is.null(pdrawsp_est)) {
       pdrawsp_est <- prepare_transformations(data = pdrawsp_est, model = model,
                                              itransform = itransform_set)
     }
     if(!is.null(pdrawsh_est)) {
       pdrawsh_est <- prepare_transformations(data = pdrawsh_est, model = model,
                                              itransform = itransform_set)
     }
   } # if(any(itransform_set != "")) {
   
   ##############################################################
   ##############################################################
   
   
   
   if(is.null(reformat)) {
     if(is.null(hypothesis) && is.null(equivalence)) {
       reformat <- TRUE
     } else {
       reformat <- FALSE
     }
   }
   
   if (reformat) {
     out_sf <- out_sf %>% 
       dplyr::rename(!!as.symbol(set_names_[1]) := 
                       dplyr::all_of('estimate')) %>% 
       dplyr::rename(!!as.symbol(set_names_[2]) := 
                       dplyr::all_of('conf.low')) %>% 
       dplyr::rename(!!as.symbol(set_names_[3]) := 
                       dplyr::all_of('conf.high')) %>% 
       data.frame()
     
     if(method == 'pkg') {
       remove_cols_ <- c('tmp_idx', 'predicted_lo', 
                         'predicted_hi', 'predicted', 'rowid')
     } else if(method == 'custom') {
       remove_cols_ <- c('term', 'contrast', 'tmp_idx', 'predicted_lo', 
                         'predicted_hi', 'predicted', 'rowid')
     }
    
    
     
     out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]
     row.names(out_sf) <- NULL
   
     
     
     
     if(!is.null(out_sf_hy)) {
       out_sf_hy <- out_sf_hy %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits))) %>% 
         dplyr::rename(!!as.symbol(set_names_[1]) := 
                         dplyr::all_of('estimate')) %>% 
         dplyr::rename(!!as.symbol(set_names_[2]) := 
                         dplyr::all_of('conf.low')) %>% 
         dplyr::rename(!!as.symbol(set_names_[3]) := 
                         dplyr::all_of('conf.high')) %>% 
         dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
         data.frame()
     } # if(!is.null(out_sf_hy)) {
     
     if(!is.null(pdrawsp_est)) {
       pdrawsp_est <- pdrawsp_est %>% 
         dplyr::rename(!!as.symbol(set_names_[1]) := 
                         dplyr::all_of('estimate')) %>% 
         dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
         data.frame()
     } # if(!is.null(pdrawsp_est)) {
     
     
     if(!is.null(pdrawsh_est)) {
       pdrawsh_est <- pdrawsh_est %>% 
         dplyr::rename(!!as.symbol(set_names_[1]) := 
                         dplyr::all_of('estimate')) %>% 
         dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
         data.frame()
     } # if(!is.null(pdrawsh_est)) {
     
   } # if (reformat) {
   
   
   
   ###########################################
   # convert factor variable that do not carry attributes ...
   as_factor_as_character_factor_df <- function(df) {
     as_factor_as_character_factor <- function(x) {
       as.factor(as.character.factor(x))
     }
     df %>% dplyr::mutate_if(is.factor, as_factor_as_character_factor )
   }
   if(!is.null(out_sf)) {
     out_sf <- as_factor_as_character_factor_df(out_sf)
   }
   if(!is.null(out_sf_hy)) {
     out_sf_hy <- as_factor_as_character_factor_df(out_sf_hy)
   }
   if(!is.null(pdrawsp_est)) {
     pdrawsp_est <- as_factor_as_character_factor_df(pdrawsp_est)
   }
   if(!is.null(pdrawsh_est)) {
     pdrawsh_est <- as_factor_as_character_factor_df(pdrawsh_est)
   }
   ###########################################
   
   
   out <- list()
   if(!is.null(out_sf)) {
     out[['estimate']] <- out_sf %>% dplyr::ungroup()
   }
   
   if(!is.null(out_sf_hy)) {
     out[['contrast']] <- out_sf_hy %>% dplyr::ungroup()
   }
   
   if(!is.null(pdrawsp_est)) {
     out[['pdrawsp_est']] <- pdrawsp_est %>% dplyr::ungroup()
   }
   
   if(!is.null(pdrawsh_est)) {
     out[['pdrawsh_est']] <- pdrawsh_est %>% dplyr::ungroup()
   }
   
   if(length(out) == 1) out <- out[[1]]
   
   return(out)
  }



#' @rdname get_predictions
#' @export
get_predictions <- function(model, ...) {
  UseMethod("get_predictions")
}



#' An alias of get_predictions()
#' @rdname get_predictions
#' @export
marginal_draws <- function(model, ...) {
  warning2c(
    "The function `marginal_draws()` has been renamed to 
    `get_predictions()`.\n",
    "Please update your code to use `get_predictions()` instead",
    " of the old function with ``marginal_`` prefix which will be removed in  
    the next release ",
    "The new name better reflect the role of this function and to harmonise the
    naming scheme across the package. In particular, the earlier name with 
    the ``marginal_`` prefix unintentionally suggested that this function
    is used only for marginal inference, whereas they in fact 
    support both marginal and conditional inferences.",
    call. = FALSE
  )
  # .Deprecated("get_predictions")
  UseMethod("get_predictions")
}
