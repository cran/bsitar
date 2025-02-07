

#' @title Estimate growth curves
#' 
#' @description The \strong{marginal_draws()} function estimates and plots
#'   growth curves (distance and velocity) by using the \pkg{marginaleffects}
#'   package as a back-end. This function can compute growth curves (via
#'   [marginaleffects::predictions()]), average growth curves (via
#'   [marginaleffects::avg_predictions()]), or plot growth curves (via
#'   [marginaleffects::plot_predictions()]). Please see
#'   [here](https://marginaleffects.com/) for details. Note that the
#'   \pkg{marginaleffects} package is highly flexible, and therefore, it
#'   is expected that the user has a strong understanding of its workings.
#'   Furthermore, since \pkg{marginaleffects} is rapidly evolving, the
#'   results obtained from the current implementation should be considered
#'   experimental.
#'
#' @details The \strong{marginal_draws()} function estimates fitted values (via
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
#'   This step is required because \strong{marginal_draws()} automatically
#'   assigns the \code{by} argument when the model includes a covariate.
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
#' @inheritParams growthparameters_comparison.bgmfit
#' @inheritParams marginal_comparison.bgmfit
#' @inheritParams marginaleffects::predictions
#' @inheritParams marginaleffects::plot_predictions
#' @inheritParams brms::fitted.brmsfit
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit()] 
#' function. Please see \code{brms::fitted.brmsfit()} for details on 
#' various options available.
#' 
#' @return An array of predicted mean response values. See
#'   [brms::fitted.brmsfit] for details.
#' 
#' @export marginal_draws.bgmfit
#' @export
#' 
#' @seealso [marginaleffects::predictions()]
#'   [marginaleffects::avg_predictions()]
#'   [marginaleffects::plot_predictions()]
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
#' # Population average distance curve
#' marginal_draws(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' marginal_draws(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' marginal_draws(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' marginal_draws(model, deriv = 1, re_formula = NULL)
#' }
#' 
marginal_draws.bgmfit <-
  function(model,
           resp = NULL,
           dpar = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           newdata = NULL,
           datagrid = NULL,
           re_formula = NA,
           allow_new_levels = FALSE,
           sample_new_levels = "gaussian", 
           parameter = NULL,
           xrange = 1,
           acg_velocity = 0.10,
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
           future_splits = NULL,
           future_method = 'future',
           future_re_expose = NULL,
           usedtplyr = FALSE,
           usecollapse = TRUE,
           cores = NULL,
           fullframe = FALSE, 
           average = FALSE, 
           plot = FALSE, 
           showlegends = NULL, 
           variables = NULL,
           condition = NULL,
           deriv = 0,
           deriv_model = TRUE,
           method = 'custom',
           marginals = NULL,
           pdrawso = FALSE, 
           pdrawsp = FALSE, 
           pdrawsh = FALSE, 
           type = NULL,
           by = NULL,
           conf_level = 0.95,
           transform = NULL,
           byfun = NULL,
           wts = NULL,
           hypothesis = NULL,
           equivalence = NULL,
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
      if(method == 'pkg') callfuns <- FALSE
    }
    
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
      ndraws <- ndraws
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
   
    
    # If default marginal effects 'dydx', then 
    call_predictions <- TRUE
    call_slopes      <- FALSE
    if(!deriv_model) {
      if(deriv > 0) {
        deriv <- 0
        call_predictions <- FALSE
        call_slopes      <- TRUE
      }
    } # if(!deriv_model) {
    
    
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
    
    if(!is.null(model$model_info$decomp)) {
      if(model$model_info$decomp == "QR") deriv_model<- FALSE
    }
    
    expose_method_set <- model$model_info[['expose_method']]
    
    model$model_info[['expose_method']] <- 'NA' # Over ride 'R'
    
    
    
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
    
    
    xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
    scall <- sys.calls()
    
    get_xcall <- function(xcall, scall) {
      scall <- scall[[length(scall)]]
      if(any(grepl("marginal_draws", scall, fixed = T)) |
         any(grepl("marginal_draws.bgmfit", scall, fixed = T))) {
        xcall <- "marginal_draws"
      } else {
        xcall <- xcall
      } 
    }
    
    if(!is.null(model$xcall)) {
      if(model$xcall == "marginal_draws") {
        xcall <- "marginal_draws"
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
    
    
    
    
    # 28.09.2024
    if(is.null(get.cores_[['max.cores']])) {
      if(is.null(arguments$cores)) 
        get.cores_[['max.cores']] <- future::availableCores() - 1
    }
    
    arguments$cores <- setincores <-  get.cores_[['max.cores']]
    .cores_ps <- get.cores_[['.cores_ps']]
    
    if (future) {
      getfutureplan <- future::plan()
      if (future_session == 'multisession') {
        set_future_session <- future_session
      } else if (future_session == 'multicore') {
        set_future_session <- future_session
      } else if (future_session == 'sequential') {
        set_future_session <- future_session
      }
      
      setplanis <- set_future_session
      if(set_future_session == 'sequential') {
        future::plan(setplanis)
      } else {
        future::plan(setplanis, workers = setincores)
      }
      on.exit(future::plan(getfutureplan), add = TRUE)
      
      # setplanis <- set_future_session
      # fuplan <- getOption("future.plan")
      # options("future.plan" = setplanis)
      # future::plan(setplanis)
      # # on.exit(options("future.plan" = fuplan), add = TRUE)
      # on.exit(future::plan(fuplan), add = TRUE)
      

      # coresplan <- getOption("mc.cores")
      # options(mc.cores = setincores)
      # on.exit(options("mc.cores" = coresplan), add = TRUE)
      
      if (future_session == 'multicore') {
        multthreadplan <- getOption("future.fork.multithreading.enable")
        options(future.fork.multithreading.enable = TRUE)
        on.exit(options("future.fork.multithreading.enable" = multthreadplan), 
                add = TRUE)
      }
      
      
      # getOption("future.plan")
      #future::availableCores()
      
      # oopts <- options(future.plan = "multisession")
      # on.exit(options(oopts))
    }
    
    
    
    draw_ids_org <- draw_ids
    draw_ids_exe <- FALSE
    if(!is.null(draw_ids)) {
      draw_ids_exe <- TRUE
      draw_ids <- draw_ids
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
          future_splits_at <- parallel::splitIndices(future_splits[1], 
                                                     future_splits[2])
        }
      }
    }
   
      
      if(future_splits_exe) {
        if(plot) {
          stop("future_splits can not be used when plot = TRUE")
        }
        if(method == 'pkg') {
          stop("future_splits can not be used when method = 'pkg'")
        }
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
              message("For multisession plan, argument 'future_re_expose' has been set as TRUE")
            }
          } else if(!need_future_re_expose_cpp) {
            if(verbose) {
              message("To speed up the calulations, it is advised to set future_re_expose = TRUE")
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
              message("To speed up the calulations, it is advised to set future_re_expose = TRUE")
            }
          } 
          if(need_future_re_expose_cpp & setplanis == "multisession") {
            # if(expose_method_set != "R") {
            stop("For plan multisession, the functions need to be re_exposed by setting future_re_expose = TRUE")
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
    
    full.args$model <- model
    full.args$deriv_model <- deriv_model
    
    full.args$newdata <- newdata
    newdata           <- do.call(get.newdata, full.args)
    
    
    
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
        deriv_model,
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
        funlist
      )
    ))[-1]
    
    # don't exclude variables
   # if(call_predictions) exclude_args <- c(exclude_args, 'variables')
    
    if(call_slopes) exclude_args <- c(exclude_args, 'transform', 'byfun')
    
    
    for (exclude_argsi in exclude_args) {
      predictions_arguments[[exclude_argsi]] <- NULL
    }
    
    
    
    
    
    if(call_slopes) {
      if (!is.null(variables)) {
        if (!is.character(variables)) {
          stop("'variables' argument must be a character string such as", 
               "\n ",
               " variables = ", "'", xvar, "'"
          )
        } else {
          set_variables <- variables
          if(!grepl(xvar, variables)) {
            set_variables <- xvar
          } else if(!is.null(set_variables[[xvar]])) {
            
          }
        }
      } else if (is.null(variables)) {
        set_variables <- xvar
      } 
    } # if(call_slopes) {
    
    
    
    # Decide if set by = NULL and then here pick and replace 'by' set_group 
    
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
    
    
    if(call_slopes) predictions_arguments$variables  <- set_variables
    predictions_arguments$by         <- set_group
    
    if(is.null(predictions_arguments$by)) predictions_arguments$by < 'NULL'
    
    
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
                                   FUN_character = setgrid_FUN_character,
                                   FUN_factor = setgrid_FUN_factor,
                                   FUN_logical = setgrid_FUN_logical,
                                   FUN_numeric = setgrid_FUN_numeric,
                                   FUN_integer = setgrid_FUN_integer,
                                   FUN_binary = setgrid_FUN_binary,
                                   FUN_other = setgrid_FUN_other
        )
        datagrid_arguments[[xvar]] <- setxvar
        if(setgrid_type == "mean_or_mode") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- set_group
        } else if(setgrid_type == "balanced") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- NULL
          # correctly set predictions_arguments[['by']]  too 
          predictions_arguments[['by']] <- 
            setdiff(predictions_arguments[['by']], cov)
        }
        set_datagrid <- do.call(marginaleffects::datagrid, datagrid_arguments)
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
    if(is.null(eval(get_draw_ids))) {
      set_draw_ids <- NULL
    } else if(is.numeric(eval(get_draw_ids))) {
      set_draw_ids <- get_draw_ids
    } else if(!eval(get_draw_ids)) {
      set_draw_ids <- NULL
    }
   predictions_arguments[['draw_ids']] <- set_draw_ids
   
   
   
   
   # 19.09.2024
   # For marginal_draws(...,  plot = T), either condition or by allowed
   # Therefore, when plot = T, condition is kept and by dropped, 
   # otherwise by is kept and condition dropped
   
   exclude_args_con_by <- exclude_args
 
   if(plot) {
     if(!is.null(predictions_arguments[['condition']]))
       predictions_arguments[['by']] <- NULL
   } else {
     predictions_arguments[['condition']] <- NULL
   }
   
   
   

   
   
   
   
   
   
   out_sf_hy <- NULL
   allowed_methods <- c('pkg', 'custom')
   if(!method %in% allowed_methods) 
     stop("Argument 'method' should be one of the following:",
          "\n ", 
          collapse_comma(allowed_methods)
     )
   # if(!is.null(predictions_arguments[['by']])) {
   #   checbyx <- predictions_arguments[['by']]
   #   if(all(checbyx == "")) method <- 'pkg'
   #   if(is.logical(checbyx)) {
   #     if(!checbyx) method <- 'pkg'
   #   }
   # }
  
   if(method == 'pkg') {
     if(call_predictions) {
       if(!plot) {
         if(!average) {
           out_sf <- do.call(marginaleffects::predictions, 
                             predictions_arguments)
         } else if(average) {
           out_sf <- do.call(marginaleffects::avg_predictions, 
                             predictions_arguments)
         }
       } else if(plot) {
         . <- do.call(marginaleffects::plot_predictions, predictions_arguments)
         outp <- .
         if(!showlegends) outp <- outp + 
           ggplot2::theme(legend.position = 'none')
         return(outp)
       }
     } # if(call_predictions) {
     
     if(call_slopes) {
       if(!plot) {
         if(!average) {
           out_sf <- do.call(marginaleffects::slopes, predictions_arguments)
         } else if(average) {
           out_sf <- do.call(marginaleffects::avg_slopes, predictions_arguments)
         }
       } else if(plot) {
         . <- do.call(marginaleffects::plot_slopes, predictions_arguments)
         outp <- .
         outp <- outp + ggplot2::theme(legend.position = 'none')
         return(outp)
       }
     } # if(call_slopes) {
   } # if(method == 'pkg') {

    
  
   # This from marginaleffects does not allow na.rm
   # So use stats::quantile instead
   get_etix <- utils::getFromNamespace("get_eti", "marginaleffects")
   get_etix <- stats::quantile
   get_hdix <- utils::getFromNamespace("get_hdi", "marginaleffects")
   get_pe_ci <- function(x, draw = NULL, na.rm = TRUE, ...) {
     if(data.table::is.data.table(x) | is.data.frame(x)) {
       if(is.null(draw)) {
         stop("please specify the 'draw' argument")
       }
       x <- x %>% dplyr::select(dplyr::all_of(draw)) %>% 
         unlist() %>% as.numeric()
     }
     if(ec_agg == "mean") estimate <- mean(x, na.rm = na.rm)
     if(ec_agg == "median") estimate <- median(x, na.rm = na.rm)
     if(ei_agg == "eti") luci = get_etix(x, probs = probs, na.rm = na.rm)
     if(ei_agg == "hdi") luci = get_hdix(x, credMass = conf)
     tibble::tibble(
       estimate = estimate, conf.low = luci[1],conf.high = luci[2]
     )
   }
   
   
   
   get_pe_ci_collapse <- function(x, na.rm = TRUE,...) {
     if(ec_agg == "mean")  estimate <- 
         collapse::fmean(x, 
                         na.rm = na.rm, 
                         nthreads = arguments$cores) 
     
     if(ec_agg == "median") estimate <- 
         collapse::fmedian(x, 
                           na.rm = na.rm, 
                           nthreads = arguments$cores)
     
     if(ei_agg == "eti") luci = collapse::fquantile(x, probs = probs, 
                                                    na.rm = na.rm)
     if(ei_agg == "hdi") luci = get_hdix(x, credMass = conf)
     cbind(estimate, luci[1], luci[2]) 
   }
   
   
   
   pdrawsp_est <- NULL
   pdrawsh_est <- NULL
   # pdraws_est <- NULL

   if(method == 'custom') {
     predictions_arguments[['cross']] <- NULL
     predictions_arguments[['method']] <- NULL
     predictions_arguments[['hypothesis']] <- NULL # hypothesis evaluated later
     by <- predictions_arguments[['by']] 
     
     
     if(future_splits_exe) {
       # Note that since predictions_arguments are passed to multisession, 
       # evaluate each argument
       for (i in names(predictions_arguments)) {
         predictions_arguments[[i]] <- eval(predictions_arguments[[i]])
       }
     }
     
     
     if(!future_splits_exe & callfuns) {
       if(call_predictions) {
         if(!plot) {
           if(!average) {
             out <- do.call(marginaleffects::predictions, predictions_arguments)
           } else if(average) {
             out <- do.call(marginaleffects::avg_predictions, predictions_arguments)
           }
         } else if(plot) {
           out <- do.call(marginaleffects::plot_predictions, 
                          predictions_arguments)
           outp <- out
           if(!showlegends) outp <- outp + 
             ggplot2::theme(legend.position = 'none')
           return(outp)
         }
       } # if(call_predictions) {
       
       if(call_slopes) {
         if(!plot) {
           if(!average) {
             out <- do.call(marginaleffects::slopes, predictions_arguments)
           } else if(average) {
             out <- do.call(marginaleffects::avg_slopes, predictions_arguments)
           }
         } else if(plot) {
           out <- do.call(marginaleffects::plot_slopes, predictions_arguments)
           outp <- out
           outp <- outp + ggplot2::theme(legend.position = 'none')
           return(outp)
         }
       } # if(call_slopes) {
     } # if(!future_splits_exe) {
     
     # print(future_splits_exe)
     # print(future_splits_exe)
     # print(future_splits_exe_future)
     # print(future_splits_exe_dofuture)
     
     
     
     
     if(future_splits_exe_future & callfuns) {
       if(call_predictions) {
         if(!plot) {
           if(!average) {
             myzfun <- function(x) {
               predictions_arguments[['draw_ids']] <- x
               predictions_arguments[['ndraws']] <- NULL
               `%>%` <- bsitar::`%>%`
               if(re_expose) {
                 if(verbose) message("need to expose functions for 'multisession'")
                 predictions_arguments[['model']] <- 
                   bsitar::expose_model_functions(predictions_arguments[['model']])
               }
                 # Re-assign appropriate function
                 setenv <- predictions_arguments[['model']]$model_info$envir
                 assign(o[[1]],
                        predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                        envir = setenv)
               do.call(marginaleffects::predictions, predictions_arguments)
             }
             out <-  future.apply::future_lapply(future_splits_at,
                                                 future.envir = parent.frame(),
                                                 future.globals = TRUE,
                                                 # future.globals = 
                                                 #   c('future_splits_at',
                                                 #     'setplanis',
                                                 #     'verbose',
                                                 #     'predictions_arguments'),
                                                 future.seed = TRUE,
                                                 FUN = myzfun)
           } else if(average) {
             myzfun <- function(x) {
               predictions_arguments[['draw_ids']] <- x
               predictions_arguments[['ndraws']] <- NULL
               `%>%` <- bsitar::`%>%`
               if(re_expose) {
                 if(verbose) message("need to expose functions for 'multisession'")
                 predictions_arguments[['model']] <- 
                   bsitar::expose_model_functions(predictions_arguments[['model']])
               }
                 # Re-assign appropriate function
                 setenv <- predictions_arguments[['model']]$model_info$envir
                 assign(o[[1]],
                        predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                        envir = setenv)
               do.call(marginaleffects::avg_predictions, predictions_arguments)
             }
             out <-  future.apply::future_lapply(future_splits_at,
                                                 future.envir = parent.frame(),
                                                 future.globals = TRUE,
                                                 # future.globals = 
                                                 #   c('future_splits_at',
                                                 #     'setplanis',
                                                 #     'verbose',
                                                 #     'predictions_arguments'),
                                                 future.seed = TRUE,
                                                 FUN = myzfun)
           }
         } else if(plot) {
           out <- do.call(marginaleffects::plot_predictions, 
                          predictions_arguments)
           outp <- out
           if(!showlegends) outp <- outp + 
             ggplot2::theme(legend.position = 'none')
           return(outp)
         }
       } # if(call_predictions) {
       
       if(call_slopes) {
         if(!plot) {
           if(!average) {
             myzfun <- function(x) {
               predictions_arguments[['draw_ids']] <- x
               predictions_arguments[['ndraws']] <- NULL
               `%>%` <- bsitar::`%>%`
               if(re_expose) {
                 if(verbose) message("need to expose functions for 'multisession'")
                 predictions_arguments[['model']] <- 
                   bsitar::expose_model_functions(predictions_arguments[['model']])
               }
                 # Re-assign appropriate function
                 setenv <- predictions_arguments[['model']]$model_info$envir
                 assign(o[[1]],
                        predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                        envir = setenv)
               do.call(marginaleffects::slopes, predictions_arguments)
             }
             out <-  future.apply::future_lapply(future_splits_at,
                                                 future.envir = parent.frame(),
                                                 future.globals = TRUE,
                                                 # future.globals = 
                                                 #   c('future_splits_at',
                                                 #     'setplanis',
                                                 #     'verbose',
                                                 #     'predictions_arguments'),
                                                 future.seed = TRUE,
                                                 FUN = myzfun)
           } else if(average) {
             myzfun <- function(x) {
               predictions_arguments[['draw_ids']] <- x
               predictions_arguments[['ndraws']] <- NULL
               `%>%` <- bsitar::`%>%`
               if(re_expose) {
                 if(verbose) message("need to expose functions for 'multisession'")
                 predictions_arguments[['model']] <- 
                   bsitar::expose_model_functions(predictions_arguments[['model']])
               }
                 # Re-assign appropriate function
                 setenv <- predictions_arguments[['model']]$model_info$envir
                 assign(o[[1]],
                        predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                        envir = setenv)
               do.call(marginaleffects::avg_slopes, predictions_arguments)
             }
             out <-  future.apply::future_lapply(future_splits_at,
                                                 future.envir = parent.frame(),
                                                 future.globals = TRUE,
                                                 # future.globals = 
                                                 # c('future_splits_at',
                                                 #   'setplanis',
                                                 #   'verbose',
                                                 #   'predictions_arguments'),
                                                 # future.seed = TRUE,
                                                 FUN = myzfun)
           }
         } else if(plot) {
           out <- do.call(marginaleffects::plot_slopes, predictions_arguments)
           outp <- out
           outp <- outp + ggplot2::theme(legend.position = 'none')
           return(outp)
         }
       } # if(call_slopes) {
     } # if(future_splits_exe_future) {
     
     
     
     
     
    
     # detach("package:doFuture", unload=TRUE) 
     
     if(future_splits_exe_dofuture & callfuns) {
       `%doFuture_function%` <- doFuture::`%dofuture%`
       # somehow .options.future = list(seed = TRUE) not working, so set below
       dofutureplan <- getOption("doFuture.rng.onMisuse")
       options(doFuture.rng.onMisuse = "ignore")
       on.exit(options("doFuture.rng.onMisuse" = dofutureplan), add = TRUE)
       
       if(call_predictions) {
         if(!plot) {
           if(!average) {
             out <- foreach::foreach(x = 1:length(future_splits_at),
                                        .options.future = list(seed = TRUE),
                                        .options.future =
                                       list(globals = structure(TRUE, add = c('future_splits_at',
                                                        'setplanis',
                                                        'verbose',
                                                        'o', 
                                                        're_expose',
                                                        'predictions_arguments')))
             ) %doFuture_function% {
               x <- future_splits_at[[x]]
               predictions_arguments[['draw_ids']] <- x
               predictions_arguments[['ndraws']] <- NULL
               `%>%` <- bsitar::`%>%`
               if(re_expose) {
                 if(verbose) message("need to expose functions for 'multisession'")
                 predictions_arguments[['model']] <-
                   bsitar::expose_model_functions(predictions_arguments[['model']])
               }
                 # Re-assign appropriate function
                 setenv <- predictions_arguments[['model']]$model_info$envir
                 assign(o[[1]],
                        predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                        envir = setenv)
              do.call(marginaleffects::predictions, predictions_arguments)
             }
           } else if(average) {
             out <- foreach::foreach(x = 1:length(future_splits_at),
                                     .options.future = list(seed = TRUE),
                                     .options.future =
                                       list(globals = c('future_splits_at',
                                                        'setplanis',
                                                        'verbose',
                                                        'o', 
                                                        're_expose',
                                                        'predictions_arguments'))
             ) %doFuture_function% {
               x <- future_splits_at[[x]]
               predictions_arguments[['draw_ids']] <- x
               predictions_arguments[['ndraws']] <- NULL
               `%>%` <- bsitar::`%>%`
               if(re_expose) {
                 if(verbose) message("need to expose functions for 'multisession'")
                 predictions_arguments[['model']] <- 
                   bsitar::expose_model_functions(predictions_arguments[['model']])
               }
                 # Re-assign appropriate function
                 setenv <- predictions_arguments[['model']]$model_info$envir
                 assign(o[[1]],
                        predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                        envir = setenv)
               do.call(marginaleffects::avg_predictions, predictions_arguments)
             }
           }
         } else if(plot) {
           out <- do.call(marginaleffects::plot_predictions, 
                          predictions_arguments)
           outp <- out
           if(!showlegends) outp <- outp + 
             ggplot2::theme(legend.position = 'none')
           return(outp)
         }
       } # if(call_predictions) {
       
       if(call_slopes) {
         if(!plot) {
           if(!average) {
             out <- foreach::foreach(x = 1:length(future_splits_at),
                                     .options.future = list(seed = TRUE),
                                     .options.future =
                                       list(globals = c('future_splits_at',
                                                        'setplanis',
                                                        'verbose',
                                                        'o', 
                                                        're_expose',
                                                        'predictions_arguments'))
             ) %doFuture_function% {
               x <- future_splits_at[[x]]
               predictions_arguments[['draw_ids']] <- x
               predictions_arguments[['ndraws']] <- NULL
               `%>%` <- bsitar::`%>%`
               if(re_expose) {
                 if(verbose) message("need to expose functions for 'multisession'")
                 predictions_arguments[['model']] <- 
                   bsitar::expose_model_functions(predictions_arguments[['model']])
               }
                 # Re-assign appropriate function
                 setenv <- predictions_arguments[['model']]$model_info$envir
                 assign(o[[1]],
                        predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                        envir = setenv)
               do.call(marginaleffects::slopes, predictions_arguments)
             }
           } else if(average) {
             out <- foreach::foreach(x = 1:length(future_splits_at),
                                     .options.future = list(seed = TRUE),
                                     .options.future =
                                       list(globals = c('future_splits_at',
                                                        'setplanis',
                                                        'verbose',
                                                        'o', 
                                                        're_expose',
                                                        'predictions_arguments'))
             ) %doFuture_function% {
               x <- future_splits_at[[x]]
               predictions_arguments[['draw_ids']] <- x
               predictions_arguments[['ndraws']] <- NULL
               `%>%` <- bsitar::`%>%`
               if(re_expose) {
                 if(verbose) message("need to expose functions for 'multisession'")
                 predictions_arguments[['model']] <- 
                   bsitar::expose_model_functions(predictions_arguments[['model']])
               }
                 # Re-assign appropriate function
                 setenv <- predictions_arguments[['model']]$model_info$envir
                 assign(o[[1]],
                        predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                        envir = setenv)
               do.call(marginaleffects::avg_slopes, predictions_arguments)
             }
           }
         } else if(plot) {
           out <- do.call(marginaleffects::plot_slopes, predictions_arguments)
           outp <- out
           outp <- outp + ggplot2::theme(legend.position = 'none')
           return(outp)
         }
       } # if(call_slopes) {
     } # if(future_splits_exe_dofuture) {
     
     
     
     
     posterior_draws_function <- function(x, ...) {
       out[[x]] %>% 
         marginaleffects:: posterior_draws(shape = "long") %>% 
         dplyr::mutate(drawid = as.numeric(drawid)) %>% 
         dplyr::mutate(drawid = future_splits_at[[x]] [.data[['drawid']]]) %>% 
         dplyr::mutate(drawid = as.factor(drawid)) %>% 
         dplyr::relocate(drawid, .before = 'draw')
     }
     
     
     consecutive_drawid_function <- function(x, ...) {
       x %>% 
         dplyr::group_by(drawid) %>% 
         dplyr::mutate(drawid = dplyr::cur_group_id()) %>% 
         dplyr::mutate(drawid = as.factor(drawid)) %>% 
         dplyr::ungroup()
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
           out <- out %>% do.call(rbind, .)
           return(out)
         }
         onex0 <- lapply(1:length(future_splits_at),  FUN = posterior_draws_function)
         onex0 <- onex0 %>% do.call(rbind, .)
         # Note that above onex0 has drawid exact same as splits
         # but somehow, we need consecutive drawid for summarising
         onex0 <- consecutive_drawid_function(onex0)
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
       oux
     }
     
     
     # trs <- lapply(1:length(out),  FUN = marginals_list_consecutive_drawid_function)
     # trs <- trs %>% do.call(rbind, .)
     # trs$drawid <- as.factor(trs$drawid)
     
     # https://stackoverflow.com/questions/19697700/how-to-speed-up-rbind
     # This below is very fast version of combined following two steps 
     # trs <- lapply(1:length(out),  FUN = marginals_list_consecutive_drawid_function)
     # trs <- trs %>% do.call(rbind, .)
     # trs <-
     # {. <- lapply(1:length(out), marginals_list_consecutive_drawid_function)
     #   list2DF(lapply(setNames(seq_along(.[[1]]), names(.[[1]])), function(i)
     #                  unlist(lapply(., `[[`, i), FALSE, FALSE)))}
     
     # trs <- lapply(1:length(out),  FUN = marginals_list_consecutive_drawid_function)
     # trs <- 
     # list2DF(lapply(setNames(seq_along(.[[1]]), names(trs[[1]])), function(i)
     #   unlist(lapply(trs, `[[`, i), FALSE, FALSE)))
     
     # This below is very fast version of below steps 
     # trs$drawid <- as.factor(trs$drawid)
     # trs$drawid <- cheapr::factor_(trs$drawid)

     
     if(setmarginals) {
       if(inherits(marginals, 'list')) {
       # if(is.list(marginals)) {
         # onex0 <- lapply(1:length(marginals),  FUN = marginals_list_consecutive_drawid_function)
         # onex0 <- onex0 %>% do.call(rbind, .)
         # onex0$drawid <- as.factor(onex0$drawid)
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
         
       }
     }
     

     # 16.10.2024
     # w hen marginals are given, then need to summarise
     if(setmarginals) {
       namesx <- c('estimate', 'conf.low', 'conf.high')
       setdrawidparm_at_ <- c(by, namesx)
       out_sf <- 
         onex0 %>%
         collapse::fgroup_by(by) %>%
         collapse::fsummarise(collapse::mctl(
           get_pe_ci_collapse(.data[['draw']]))
         )  %>% collapse::frename(., setdrawidparm_at_)
     } # if(setmarginals) {
     
     
     
     # 4.10.2024
     # The by is setting as set_group instead of by from argument
     # Therefore checking if not FALSE
     
     if(!setmarginals) {
       setdrawidparm <- by
       namesx <- c('estimate', 'conf.low', 'conf.high')
       if(!isFALSE(setdrawidparm)) setdrawidparm_ <- c(setdrawidparm, namesx)
       if( isFALSE(setdrawidparm)) setdrawidparm_ <- c(namesx)
       
       out_sf <- onex0 %>% collapse::fsubset(., drawid == 1) %>%
         collapse::fselect(., setdrawidparm_)
     }
     
     
    

     
     if(!is.null(hypothesis)) {
       # For hypothesis
       groupvarshyp1 <- c('drawid')
       groupvarshyp2 <- c('term')
       if(!is.null(constrats_at)) {
         if(!is.character(constrats_at)) 
           stop("'constrats_at' must be a character vector")
         for (caxi in constrats_at) {
           if(!caxi %in% names(onex0)) {
             stop("Variable '", caxi, ". specified in 'constrats_at' is not in",
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
       
       
       if(is.null(constrats_by)) {
         stop("Please specify 'constrats_by' argument when testing 'hypothesis'",
              "\n ",
              " The available options are: ",
              collapse_comma(by)
         )
       }
       
       if(!is.null(constrats_by)) {
         if(!is.character(constrats_by)) 
           stop("'constrats_by' must be a character vector")
         for (caxi in constrats_by) {
           if(!caxi %in% names(onex0)) {
             stop("Variable '", caxi, ". specified in 'constrats_by' is not in",
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
       
       
       if(nrow(out_sf) > 25) {
         if(is.null(constrats_at)) {
           cat(" Note that the 'marginaleffects' package does not allow" ,
               "\n",
               "'hypothesis' argument when estimates rows are more than 25",
               "\n",
               "To avoid this issue, you can use 'constrats_at' argument",
               "\n"
           )
         }
       }
       
       
       # get_pe_ci_collapse <- function(x, na.rm = TRUE,...) {
       #   if(ec_agg == "mean")  estimate <- 
       #       collapse::fmean(x, 
       #                       na.rm = na.rm, 
       #                       nthreads = arguments$cores) 
       #   
       #   if(ec_agg == "median") estimate <- 
       #       collapse::fmedian(x, 
       #                         na.rm = na.rm, 
       #                         nthreads = arguments$cores)
       #   
       #   if(ei_agg == "eti") luci = collapse::fquantile(x, probs = probs, 
       #                                                  na.rm = na.rm)
       #   if(ei_agg == "hdi") luci = get_hdix(x, credMass = conf)
       #   cbind(estimate, luci[1], luci[2]) 
       # }
       
       set_constrats_by <- c(constrats_by, 'draw')
       namesx <- c('estimate', 'conf.low', 'conf.high')
       
       setdrawidparm_at <- c(constrats_at, 'term')
       setdrawidparm_at_ <- c(setdrawidparm_at, namesx)
       
       
       temhyy <-
         onex0 %>% 
         collapse::fgroup_by(groupvarshyp1) %>%
         collapse::fselect(set_constrats_by) %>% 
         collapse::frename('estimate' = 'draw') %>% 
         collapse::fsummarise(collapse::qDF(
           get_hypothesis_x(.data,
                            by = constrats_by,
                            hypothesis = hypothesis,
                            draws = 'estimate'))) 
       
       out_sf_hy <- 
       temhyy %>%
         collapse::fgroup_by(groupvarshyp2) %>%
         collapse::fsummarise(collapse::mctl(
           get_pe_ci_collapse(.data[['estimate']]))
         ) %>%
         collapse::frename(., setdrawidparm_at_)
       
       
       if(!isFALSE(pdrawsh)) {
         selectchoicesr <- c("return", 'add') 
         checkmate::assert_choice(pdrawsh, choices = selectchoicesr)
         if(pdrawsh == 'return') {
           return(temhyy)
         } else if(pdrawsh == 'add') {
           pdrawsh_est <- temhyy
         } else {
           
         }
       }
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
         # dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
         dplyr::rename(!!as.symbol(set_names_[1]) := 
                         dplyr::all_of('estimate')) %>% 
         # For pdrawsp_est and pdrawsh_est, there are no conf columns, only estimates
         # dplyr::rename(!!as.symbol(set_names_[2]) := 
         #                 dplyr::all_of('conf.low')) %>% 
         # dplyr::rename(!!as.symbol(set_names_[3]) := 
         #                 dplyr::all_of('conf.high')) %>% 
         dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
         data.frame()
     } # if(!is.null(pdrawsp_est)) {
     
     if(!is.null(pdrawsh_est)) {
       pdrawsh_est <- pdrawsh_est %>% 
         # dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
         dplyr::rename(!!as.symbol(set_names_[1]) := 
                         dplyr::all_of('estimate')) %>% 
         # For pdrawsp_est and pdrawsh_est, there are no conf columns, only estimates
         # dplyr::rename(!!as.symbol(set_names_[2]) := 
         #                 dplyr::all_of('conf.low')) %>% 
         # dplyr::rename(!!as.symbol(set_names_[3]) := 
         #                 dplyr::all_of('conf.high')) %>% 
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
   if(!is.null(out_sf)) out_sf <- as_factor_as_character_factor_df(out_sf)
   if(!is.null(out_sf_hy)) out_sf_hy <- as_factor_as_character_factor_df(out_sf_hy)
   if(!is.null(pdrawsp_est)) pdrawsp_est <- as_factor_as_character_factor_df(pdrawsp_est)
   if(!is.null(pdrawsh_est)) pdrawsh_est <- as_factor_as_character_factor_df(pdrawsh_est)
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


#' @rdname marginal_draws.bgmfit
#' @export
marginal_draws <- function(model, ...) {
  UseMethod("marginal_draws")
}


