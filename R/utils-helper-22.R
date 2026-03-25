

#' @title Estimate model based individual growth parameters
#' 
#' @description The \strong{modelbased_growthparameters_call()} function estimates
#'   individual growth parameters by mapping population average estimate of age
#'   of interest (such as age at peak growth velocity or age at take off) on
#'   individual velocity curves defined by individual level random effects.
#'   After the individual level age parameter is found, the corresponding
#'   distance and velocity is also estimated.
#'   
#' @details Since SITAR is a shape-invariant model, each individual curve has a  
#' peak velocity point that can be mapped by knowing the population average age 
#' at peak velocity. This hold true even when a individual lacks measurements at 
#' the expected turning point.
#' 
#' @inheritParams growthparameters.bgmfit
#' @inheritParams get_growthparameters.bgmfit
#' @inheritParams get_comparisons.bgmfit
#' @inheritParams marginaleffects::predictions
#' @inheritParams marginaleffects::plot_predictions
#' @inheritParams brms::fitted.brmsfit
#' 
#' @param ... Additional arguments passed to the function. 
#' 
#' @return A named list of 3 comprising individual level estimate of
#'   \strong{age}, \strong{distance} and \strong{velocity}. Each of the list is
#'   a data frame with one row per individual and six columns.\cr
#'   
#'   \strong{age} \cr \item{id}{subject identifier} \item{Estimate}{subject's
#'   age corresponding to \code{x}.} \item{Est.Error}{SD of Estimate} \item{Q2.5
#'   }{Lower CI} \item{Q97.5}{Upper CI}
#'   \item{missing}{logical flags where TRUE means subject's specified age lies
#'   outside their measurement range}
#'   
#'   \strong{distance} \cr \item{id}{subject identifier}
#'   \item{Estimate}{distance corresponding to subject's age}
#'   \item{Est.Error}{SD of Estimate} \item{Q2.5 }{Lower CI} \item{Q97.5}{Upper
#'   CI}
#'   \item{missing}{logical flags where TRUE means subject's specified age lies
#'   outside their measurement range}
#'   
#'   \strong{velocity} \cr \item{id}{subject identifier}
#'   \item{Estimate}{velocity corresponding to subject's age}
#'   \item{Est.Error}{SD of Estimate} \item{Q2.5 }{Lower CI} \item{Q97.5}{Upper
#'   CI}
#'   \item{missing}{logical flags where TRUE means subject's specified age lies
#'   outside their measurement range}
#' 
#' @keywords internal
#' @noRd
#' 
modelbased_growthparameters_call.bgmfit <-
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
           parallel = FALSE,
           cores = NULL,
           fullframe = FALSE, 
           average = FALSE, 
           plot = FALSE, 
           showlegends = NULL, 
           variables = NULL,
           deriv = NULL,
           model_deriv = NULL,
           method = 'pkg',
           marginals = NULL, 
           pdraws = FALSE, 
           pdrawso = FALSE,
           pdrawsp = FALSE, 
           pdrawsh = FALSE, 
           comparison = "difference",
           type = NULL,
           by = FALSE,
           bys = NULL,
           condition = NULL,
           conf_level = 0.95,
           transform = NULL,
           transform_draws = NULL,
           cross = FALSE,
           wts = NULL,
           hypothesis = NULL,
           equivalence = NULL,
           eps = NULL,
           constrats_by = FALSE,
           constrats_at = FALSE,
           constrats_subset = FALSE,
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
    
    
    # 6.03.2025
    # Also, set 'marginaleffects_lean' to FALSE for 'posterior_draws()' to work
    lean_ <- getOption("marginaleffects_lean")
    options("marginaleffects_lean" = FALSE)
    on.exit(options("marginaleffects_lean" = lean_), add = TRUE)
    
    
    
    try(zz <- insight::check_if_installed(c("marginaleffects"), 
                                          minimum_version = 
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
    
    if(usecollapse) {
      usedtplyrcheck <- usedtplyr
      usedtplyr <- FALSE
      if(verbose & usedtplyrcheck) message("Setting usedtplyr = FALSE because ",
                                           "usecollapse = TRUE")
    } else {
      usedtplyr <- usedtplyr
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
    
    if(usecollapse) {
      try(zz <- insight::check_if_installed(c("collapse"), 
                                            minimum_version =
                                              get_package_minversion(
                                                'collapse'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
      
      try(zz <- insight::check_if_installed(c("doParallel"), 
                                            minimum_version =
                                              get_package_minversion(
                                                'doParallel'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
      
      try(zz <- insight::check_if_installed(c("foreach"), 
                                            minimum_version =
                                              get_package_minversion(
                                                'foreach'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
      
      try(zz <- insight::check_if_installed(c("parallel"), 
                                            minimum_version =
                                              get_package_minversion(
                                                'parallel'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
    }
    
    
    if(method == 'pkg') {
      stop("use method == 'custom' for modelbased_growthparameters_call()")
    }
    
    
    callfuns     <- TRUE
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
    
    # 20.03.2025
    assign_function_to_environment(transform_draws, 'transform_draws', 
                                   envir = NULL)
    model$model_info[['transform_draws']] <- transform_draws
    
    # 20.03.2025
    # Depending on dpar 'mu' or 'sigma', subset model_info
    # This only when set_sigma_manual used to model a b c 
    # Not when a function such as splines::ns etc used in sigma_formula
   
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
    if(is.null(xvar)) xvar   <- model$model_info[[xvar_]]
    cov_   <- paste0('cov', resp_rev_)
    cov    <- model$model_info[[cov_]]
    uvarby <- model$model_info$univariate_by$by
    
    ########################################################
    groupvar_ <- paste0('groupvar', resp_rev_)
    
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
    
    
    
    
    ########################################################
    # funx_ <- paste0('xfuntransform2', resp_rev_)
    
    
    check_set_fun <- check_set_fun_transform(model = model, 
                                             which = 'xfuntransform2',
                                             dpar = dpar, 
                                             resp= resp, 
                                             transform = NULL,
                                             auto = FALSE, 
                                             verbose = verbose)
    
    funx_ <- check_set_fun[['setfun']]
    if(check_set_fun[['was_null']]) {
      model$model_info[[check_set_fun[['setfunname']]]] <- funx_
    }
    
    
    ########################################################
    
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
    
    
    
    ########################################################
    
    
    # Over ride 'ifunx_()' i.e, return xvar on scale used in model fit
    # This is restricted to  when using 'pdrawsp' etc. 
    # use cae ->  get_dv
    # 6.03.205
    # itransform_set <- get_itransform_call(itransform)
    
    itransform_set <- get_itransform_call(itransform = itransform,
                                          model = model, 
                                          newdata = newdata,
                                          dpar = dpar, 
                                          resp = resp,
                                          auto = FALSE,
                                          verbose = verbose)
    
    if(itransform_set == "") {
      if(!isFALSE(pdrawsp)) {
        if(!is.character(pdrawsp)) pdrawsp <- "return"
        selectchoicesr <- c("return", 'add') 
        checkmate::assert_choice(pdrawsp, choices = selectchoicesr)
        if(pdrawsp == 'return' | pdrawsp == 'add') {
          ifunx_ <- function(x)x
        } 
      } # if(!isFALSE(pdrawsp)) {
    } # if(itransform_set == "") {
    ########################################################
    
    
    
    if(is.null(model_deriv)) {
      if(is.null(deriv)) {
        model_deriv <- FALSE
      } else if(deriv == 0) {
        model_deriv <- FALSE
      } else if(deriv == 1) {
        model_deriv <- TRUE
      }
    } else if(!is.null(model_deriv)) {
      if(is.null(deriv) & !model_deriv) {
        deriv <- 0
        model_deriv <- FALSE
      } else if(is.null(deriv) & model_deriv) {
        deriv <- 1
        model_deriv <- TRUE
      }
    }
    
    
    # if(is.null(deriv) & is.null(model_deriv)) {
    #   deriv <- 0
    #   model_deriv <- FALSE
    # } else if(deriv == 0 & is.null(model_deriv)) {
    #   deriv <- 0
    #   model_deriv <- FALSE
    # } else if(deriv == 1 & is.null(model_deriv)) {
    #   deriv <- 1
    #   model_deriv <- TRUE
    # } else if(is.null(deriv) & !model_deriv) {
    #   deriv <- 0
    #   model_deriv <- FALSE
    # } else if(is.null(deriv) & model_deriv) {
    #   deriv <- 1
    #   model_deriv <- TRUE
    # }
    
    
    
    # 15 06 2025
    allowed_methods <- c('pkg', 'custom')
    if(!method %in% allowed_methods) 
      stop("Argument 'method' should be one of the following:",
           "\n ", 
           collapse_comma(allowed_methods)
      )
    
    
    
    if(method == 'custom') {
      deriv <- 1
      model_deriv <- TRUE
      if(verbose) message("For method = 'custom', deriv is set to TRUE")
    }
    
    
    
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
    drawid <- NULL;
    
    draw <- NULL;
    j <- NULL;
    i <- NULL;
    
    
    
    # allowed_parms <- c(
    #   'apgv',
    #   'pgv',
    #   'atgv',
    #   'tgv',
    #   'acgv',
    #   'cgv'
    # )
    # 
    # if(is.null(parameter)) parameter <- c('apgv', 'pgv')
    # 
    # parameter <- base::tolower(parameter)
    # 
    # if (is.null(parameter)) {
    #   parm <- 'apgv' 
    # } else if(length(parameter) == 1 && parameter == 'all') {
    #   parm <- allowed_parms 
    # } else if(length(parameter) == 1) {
    #   parm <- parameter
    # } else if(length(parameter) > 1) {
    #   # parameter <- base::tolower(parameter)
    #   for (parameteri in parameter) {
    #     if(!parameteri %in% allowed_parms) {
    #       allowed_parms_err <- c(allowed_parms, 'all')
    #       stop("Allowed parameter options are ", 
    #            paste(paste0("'", allowed_parms_err, "'"), collapse = ", ")
    #       )
    #     }
    #   }
    #   parm <- parameter
    # }
    # parm <- base::tolower(parm)
    # 
    # 
    # if(length(parm) > 1) {
    #   if(plot) stop("Please specify only one parameter")
    # }
    
    
    allowed_parms <- c(
      'apgv',
      'atgv',
      'acgv'
    )
    
    
    if(is.null(parameter)) {
      parameter <- c('apgv')
    } else if(is.numeric(parameter)) {
      parameter <- parameter
    } else if(is.character(parameter)) {
      parameter <- parameter
    } else if(is.symbol(parameter)) {
      parameter <- deparse(parameter)
    } else if(is.matrix(parameter)) {
      
    } else {
      stop("parameter must be a character string of a numeric value")
    }
    
    
    
    if(length(parameter) > 1) {
      stop("'parameter' must a length one")
    }
    
    
    
    parameter_numeric <- FALSE
    if(is.numeric(parameter)) {
      parameter_numeric <- TRUE
    } else {
      if(check_is_numeric_like(parameter)) {
        parameter_numeric <- TRUE
      }
    }
    
    
    if(is.character(parameter) & !parameter_numeric) {
      parameter <- base::tolower(parameter)
      for (parameteri in parameter) {
        if(!parameteri %in% allowed_parms) {
          allowed_parms_err <- allowed_parms # c(allowed_parms, 'all')
          stop("Allowed parameter options are ", 
               paste(paste0("'", allowed_parms_err, "'"), collapse = ", ")
          )
        }
      }
    } # if(is.character(parameter)) {
    
    
    
    parm <- parameter
    
    
    conf <- conf_level
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', probtitles)
    
    # if(!is.null(model$model_info$decomp)) {
    #   if(model$model_info$decomp == "QR") model_deriv<- FALSE
    # }
    
    expose_method_set <- model$model_info[['expose_method']]
    
    model$model_info[['expose_method']] <- 'NA' # Over ride 'R'
    
    
    setxcall_   <- match.call()
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
    
    # o <- post_processing_checks(model = model,
    #                             xcall = match.call(),
    #                             resp = resp,
    #                             envir = envir,
    #                             deriv = deriv, 
    #                             all = FALSE,
    #                             verbose = verbose, 
    #                             check_d1 = T)
    # 
    # oall <- post_processing_checks(model = model,
    #                                xcall = match.call(),
    #                                resp = resp,
    #                                envir = envir,
    #                                deriv = deriv, 
    #                                all = TRUE,
    #                                verbose = FALSE)
    
    
    
    
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
    
    
    
    ########################################################
    # method = 'pkg' uses 'y <- (hi - lo)...' approach for deriv = 0 or deriv > 0
    # see gparms_fun call_comparison_gparms_fun
    # So it work when assigned d0 or d1
    # However, for method = 'custom', need to switch between
    # marginaleffects::predictions /  marginaleffects::slopes
    ########################################################
    
    ########################################################
    # Decide here -> no 'x', 'log' or 'sqrt'
    # Then need 'model_deriv = FALSE' and 'deriv = 0'
    # This because print(o[[2]]) function d1 is not reliable- see prepare_function
    # This has already been taken care off in prepare_function by excluding d1
    ########################################################
    # This been moved to post_processing_checks
    
    # available_fund1 <- FALSE
    # for (i in names(model$model_info$exefuns)) {
    #   check_funds <- ifelse(grepl("\\d$", i), sub(".*?(\\d+)$", "\\1", i), "")
    #   if(grepl("1", check_funds)) {
    #     available_fund1 <- TRUE
    #   }
    # }
    
    ########################################################
    # If no d1, override and and deriv arguments 
    ########################################################
    # This been moved to post_processing_checks
    # if(!available_fund1) {
    #   deriv       <- 0
    #   model_deriv <- FALSE
    #   if(verbose) {
    #     message("Since no 'd1' found, setting 'model_deriv = FALSE', 'deriv = 0'")
    #   }
    # }
    ########################################################
    
    
    # Below deriv/model_deriv will be over riddent for both 'pkg' and 'custom'
    
    # This borrowed from get_predictions
    call_predictions <- TRUE
    call_slopes      <- FALSE
    
    available_d1 <- o[['available_d1']]
    o_org <- o
    if(!available_d1) {
      deriv       <- 0
      model_deriv <- FALSE
      call_predictions <- FALSE
      call_slopes      <- TRUE
      # re-get o[[2]] as _do
      post_processing_checks_args[['deriv']]    <- 0
      o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
      o_available_d1_F <- o
    }
    
    post_processing_checks_args[['deriv']]    <- deriv
    
    
    # 20.03.2025
    if(!is.null(model$model_info[['sigma_fun_mode']])) {
      sigma_fun_mode <- model$model_info[['sigma_fun_mode']]
      if(dpar == "sigma") {
        if(deriv > 0) {
          if(sigma_fun_mode == "inline") {
            check_fun    <- TRUE
            available_d1 <- FALSE
            model_deriv  <- FALSE
            call_slopes  <- TRUE
          }
        }
      }
    }
    
    
    
    # test <- setupfuns(model = model, resp = resp,
    #                   o = o, oall = oall, 
    #                   usesavedfuns = usesavedfuns, 
    #                   deriv = post_processing_checks_args[['deriv']], 
    #                   envir = envir, 
    #                   model_deriv = model_deriv, 
    #                   ...)
    # 
    # if(is.null(test)) return(invisible(NULL))
    # 
    
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
    
    
    
    # xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
    # scall <- sys.calls()
    # 
    # get_xcall <- function(xcall, scall) {
    #   scall <- scall[[length(scall)]]
    #   if(any(grepl("modelbased_growthparameters_call", scall, fixed = T)) |
    #      any(grepl("modelbased_growthparameters_call.bgmfit", scall, fixed = T))) {
    #     xcall <- "modelbased_growthparameters_call"
    #   } else {
    #     xcall <- xcall
    #   }
    # }
    # 
    # 
    # 
    # if(xcall == "do.call" | xcall == "CustomDoCall") {
    #   zzz <- gsub_space(paste(deparse(sys.calls()[[1]]), collapse = ""))
    #   zzz <- regmatches(zzz, gregexpr("(?<=\\().*?(?=\\))", zzz, perl=T))[[1]]
    #   zzz <- strsplit(zzz, ",")[[1]][1]
    #   xcall <- strsplit(zzz, "\\.")[[1]][1]
    # } else {
    #   if(!is.null(model$xcall)) {
    #     if(model$xcall == "modelbased_growthparameters_call") {
    #       xcall <- "modelbased_growthparameters_call"
    #     }
    #   } else {
    #     scall <- sys.calls()
    #     xcall <- get_xcall(xcall, scall)
    #   }
    # }
    
    
    # if(!is.null(model$xcall)) {
    #   if(model$xcall == "modelbased_growthparameters_call") {
    #     xcall <- "modelbased_growthparameters_call"
    #   }
    # } else {
    #   scall <- sys.calls()
    #   xcall <- get_xcall(xcall, scall)
    # }
    
    
    # xcall <- xcall
    
    # for first TRUE, use min(which(lv == TRUE))
    
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
      
      if (future_session == 'multicore') {
        multthreadplan <- getOption("future.fork.multithreading.enable")
        options(future.fork.multithreading.enable = TRUE)
        on.exit(options("future.fork.multithreading.enable" = multthreadplan), 
                add = TRUE)
      }
    }
    
    
    
    
    draw_ids_org <- draw_ids
    draw_ids_exe <- FALSE
    if(!is.null(draw_ids)) {
      draw_ids_exe <- TRUE
      draw_ids <- draw_ids
      eval_draw_ids <- eval(draw_ids)
    } else if(is.null(draw_ids)) {
      eval_draw_ids <- seq(1, ndraws, 1)
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
              message("To speed up the calulations, it is advised ",
                      "to set 'future_re_expose = TRUE'")
            }
          } 
          if(need_future_re_expose_cpp & setplanis == "multisession") {
            # if(expose_method_set != "R") {
            stop("For plan 'multisession', the functions need to be ",
                 "\n ",
                 "re_exposed by setting 'future_re_expose = TRUE'")
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
                                    fargs = arguments, 
                                    dargs = list(...), 
                                    verbose = verbose)
    
    full.args$model <- model
    full.args$model_deriv <- model_deriv
    
    
    if(is.null(full.args$hypothesis) & is.null(full.args$equivalence)) {
      plot <- plot
    } else {
      plot <- FALSE
      if(verbose & plot) {
        message("Argument plot = TRUE is not allowed when either hypothesis ", 
                "or equivalence is not NULL",
                "\n ",
                "Therefor, setting 'plot = FALSE'") 
      }
    }
    
    
    full.args$newdata <- newdata
    
    full.args <- 
      sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                 arguments = full.args, 
                                 check_formalArgs = NULL,
                                 check_formalArgs_exceptions = NULL,
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
   
    # Interpolation points
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    full.args$ipts <- ipts <- check_ipts(ipts = full.args$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)
    
    full.args$dpar <- dpar
    
    get.newdata_args <- list()
    for (i in methods::formalArgs(get.newdata)) {
      get.newdata_args[[i]] <- full.args[[i]]
    }
    
    newdata        <- CustomDoCall(get.newdata, get.newdata_args)
    
    

    
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
      full.args[['transform']] <- full.args[['transform_draws']]
      if(verbose) message("'transform' set based on 'transform_draws'")
    }
    
    
    comparisons_arguments <- full.args
    
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
        ipts,
        seed,
        future,
        future_session,
        future_splits,
        future_method,
        future_re_expose,
        cores,
        fullframe,
        dummy_to_factor,
        verbose,
        expose_function,
        usesavedfuns,
        clearenvfuns,
        envir,
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
        constrats_subset,
        usedtplyr,
        usecollapse, 
        parallel,
        cores,
        pdraws,
        pdrawso,
        pdrawsp,
        pdrawsh,
        bys,
        funlist,
        itransform,
        newdata_fixed,
        transform_draws,
        incl_autocor,
        internalmethod
      )
    ))[-1]
    
    if(plot) {
      exclude_args <- c(exclude_args, "cross")
    }
    
    for (exclude_argsi in exclude_args) {
      comparisons_arguments[[exclude_argsi]] <- NULL
    }
    
    
    
    if(deriv == 0 & model_deriv) 
      stop("If argument 'model_deriv = TRUE', then 'deriv' should be '1'")
    if(deriv == 1 & !model_deriv) 
      stop("If argument 'model_deriv' = FALSE, then 'deriv' should be '0'")
    
    
    if (!is.null(variables)) {
      if (!is.list(variables)) {
        if(!model_deriv) {
          stop("'variables' argument must be a named list and the first ", 
               "element should be ", xvar, 
               "\n ",
               " specified as follows ",
               "\n ",
               " variables = list(", xvar, "=", "1e-6",")",
               "\n ",
               " where 1e-6 is the default value for the argument 'eps'"
          )
        }
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
    
    
    
    allowed_comparison <- c('difference', 'differenceavg')
    
    if(!comparison %in% allowed_comparison) {
      stop("Allowed comparison options are ", 
           paste(paste0("'", allowed_comparison, "'"), collapse = ", ")
      )
    }
    
    
    if(comparison == 'differenceavg') {
      if(!average) {
        stop("For comparison = 'differenceavg' ", 
             ", the argument 'average' should be TRUE")
      }
      if(is.null(hypothesis)) {
        stop("For comparison = 'differenceavg' ", 
             ", the argument 'hypothesis' is required.",
             " \n",
             "An example of Non-linear hypothesis testing via hypothesis",
             " argument is as follows:",
             " \n",
             "hypothesis = 'b2 - b1 = 0.2'",
             " where b2 and b1 are row indices",
             " \n",
             "(see https://marginaleffects.com/vignettes/comparisons.html)",
             " for more details",
             " \n",
             "Note that user need to set comparison = 'differenceavg'",
             " and average = TRUE when ",
             " \n", 
             "testing non-linear hypothesis as described above"
        )
      }
    }
    
    
    set_group <- setup_by_var(model = model, 
                              by = by, 
                              cov = cov, 
                              xvar = xvar, 
                              dpar = dpar,
                              method = 'pkg',
                              plot = plot,
                              condition = condition,
                              deriv = NULL,
                              difx = NULL,
                              xcall = xcall,
                              xvar_strict = TRUE,
                              switch_plot = FALSE,
                              verbose = verbose)
    
    
    if (acg_velocity >= 1 | acg_velocity <= 0) {
      stop("The acg_velocity should be set between 0.01 and 0.99")
    }
    
    call_comparison_gparms_fun <- function(parm, 
                                           eps, 
                                           by, 
                                           aggregate_by, 
                                           newdata, ...) {
      gparms_fun = function(hi, lo, x, ...) {
        if(deriv == 0) y <- (hi - lo) / eps
        if(deriv > 0)  y <- (hi + lo) / 2
        # Here need to aggregate based on by argument
        if(aggregate_by) {
          try(insight::check_if_installed(c("grDevices", "stats"), stop = FALSE, 
                                          prompt = FALSE))
          xy <- grDevices::xy.coords(x, y)
          xy <- unique(as.data.frame(xy[1:2])[order(xy$x), ])
          if(!isFALSE(by)) {
            if(ec_agg == "mean")   xy <- stats::aggregate(.~x, data=xy, 
                                                          mean, 
                                                          na.action = na.omit,
                                                          drop = TRUE)
            if(ec_agg == "median") xy <- stats::aggregate(.~x, data=xy, 
                                                          median, 
                                                          na.action = na.omit,
                                                          drop = TRUE)
          }
          x <- xy$x
          y <- xy$y
        } # if(aggregate_by) {
        
        
        
        if (parm == 'apgv') {
          out <- sitar::getPeak(x = x, y = y)[1]
          out <- ifunx_(out) # prepare_data2
        } else if (parm == 'pgv') {
          out <- sitar::getPeak(x = x, y = y)[2]
        } else if (parm == 'atgv') {
          out <- sitar::getTakeoff(x = x, y = y)[1]
          out <- ifunx_(out) # prepare_data2
        } else if (parm == 'tgv') {
          out <- sitar::getTakeoff(x = x, y = y)[2]
        } else if (parm == 'acgv') {
          cgv  <- acg_velocity * sitar::getPeak(x = x, y = y)[2]
          vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
          out <-  x[vcgi]
          out <- ifunx_(out) # prepare_data2
        } else if (parm == 'cgv') {
          cgv  <- acg_velocity * sitar::getPeak(x = x, y = y)[2]
          vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
          out <-  y[vcgi]
        } else if (parm == 'xxxx') {
          
        } else {
          stop('parm not valid')
        }
        out <- round(out, digits = digits)
        out
      } # gparms_fun
      
      comparisons_arguments$variables  <- set_variables
      comparisons_arguments$by         <- set_group
      comparisons_arguments$comparison <- gparms_fun
      
      
      # This is for method = 'pkg' - does't matter call_predictions/call_slopes?
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
      
      
      # 19.09.2024
      # For get_predictions(...,  plot = T), either condition or by allowed
      # Therefore, when plot = T, condition is kept and by dropped, 
      # otherwise by is kept and condition dropped
      
      exclude_args_con_by <- exclude_args
      
      # Need to move belwo in outer_call_comparison_gparms_fun
      if(plot) {
        if(!is.null(comparisons_arguments[['condition']]))
          comparisons_arguments[['by']] <- NULL
      } else {
        comparisons_arguments[['condition']] <- NULL
      }
      
      
      if(!plot) {
        if(!average) {
          if(callfuns) out <- do.call(marginaleffects::comparisons, 
                                      comparisons_arguments)
        } else if(average) {
          if(callfuns) out <- do.call(marginaleffects::avg_comparisons, 
                                      comparisons_arguments)
        }
      }
      if(plot) {
        if(isFALSE(set_group)) comparisons_arguments$by <- NULL
        if(is.null(comparisons_arguments[['by']]) &
           is.null(comparisons_arguments[['condition']])) {
          stop("'plot = TRUE' is not yet supported for ",
               "'modelbased_growthparameters_call()'"
          )
        }
        out <- do.call(marginaleffects::plot_comparisons, 
                       comparisons_arguments)
        outp <- out
        if(!showlegends) outp <- outp + ggplot2::theme(legend.position = 'none')
        return(outp)
      }
      
      return(out)
    } # call_comparison_gparms_fun
    
    
    
    #################
    enverr. <- environment()
    outer_call_comparison_gparms_fun <- function(parm, 
                                                 eps,
                                                 by,
                                                 aggregate_by,
                                                 newdata,
                                                 ...) {
      
      assign('err.', FALSE, envir = enverr.)
      
      tryCatch(
        expr = {
          # 13 03 2024
          gout <- call_comparison_gparms_fun(parm = parm, eps = eps,
                                             by = by,
                                             aggregate_by = aggregate_by,
                                             newdata = newdata)
        },
        error = function(e) {
          assign('err.', TRUE, envir = enverr.)
          message(e)
        }
      )
      err. <- get('err.', envir = enverr.)
      
      # This happen when NA in one factor for hypothesis
      if(!exists('gout')) return(invisible(NULL))
      
      if(length(gout) == 0) err. <- TRUE
      
      if (err.) {
        gout <- NULL
      } else if (!err.) {
        gout <- gout
      }
      
      
      if(plot) {
        return(gout)
      }
      
      # If results not available for all factor, theN add NA it
      if(!is.null(gout)) {
        if(!isFALSE(by)) {
          if(is.null(hypothesis) && is.null(equivalence)) {
            goutnames <- names(gout)
            groupvars <-  eval(by)
            newdatajoin <- newdata %>% dplyr::group_by_at(groupvars) %>%
              dplyr::filter(dplyr::row_number() == 1) %>% dplyr::ungroup()
            gout <- newdatajoin %>% dplyr::left_join(., gout, by = groupvars) %>%
              dplyr::select(dplyr::all_of(goutnames)) %>%
              dplyr::select(-dplyr::any_of(c('predicted_lo', 'predicted_hi',
                                             'predicted', 'tmp_idx'))) %>%
              dplyr::mutate(!! as.name('parameter') := parm) %>%
              dplyr::relocate(dplyr::all_of('parameter')) %>%
              data.frame()
          }
        }
      }
      return(gout)
    }
    
    
    
    
    
    ###############
    eval_re_formula <- eval(comparisons_arguments$re_formula)
    if(is.null(eval_re_formula)) {
      aggregate_by <- TRUE
    } else if(is.na(eval_re_formula)) {
      aggregate_by <- FALSE
    }
    
    
    #######################################
    out_sf_hy <- NULL
    allowed_methods <- c('pkg', 'custom')
    if(!method %in% allowed_methods) 
      stop("Argument 'method' should be one of the following:",
           "\n ", 
           collapse_comma(allowed_methods)
      )
    
    # if(!is.null(comparisons_arguments[['by']])) {
    #   checbyx <- comparisons_arguments[['by']]
    #   if(all(checbyx == "")) method <- 'pkg'
    #   if(is.logical(checbyx)) {
    #     if(!checbyx) method <- 'pkg'
    #   }
    # }
    
    
    
    if(method == 'pkg') {
      
      if(plot) {
        out_sf <- outer_call_comparison_gparms_fun(
          parm = parm, eps = eps, 
          by = comparisons_arguments$by,
          aggregate_by = aggregate_by,
          newdata = newdata
        ) 
        return(out_sf)
      }
      
      
      if (length(parm) == 1) {
        out_sf <- outer_call_comparison_gparms_fun(
          parm = parm, eps = eps, 
          by = comparisons_arguments$by,
          aggregate_by = aggregate_by,
          newdata = newdata
        ) 
        if(is.null(out_sf)) return(invisible(NULL))
        if(!"parameter" %in% colnames(out_sf)) {
          out_sf <- out_sf %>% dplyr::mutate(!!as.symbol('parameter') := parm) %>%
            dplyr::relocate(!!as.symbol('parameter')) %>% data.frame() 
        }
      } else if (length(parm) > 1) {
        list_cout <- list()
        for (allowed_parmsi in parm) {
          list_cout[[allowed_parmsi]] <- outer_call_comparison_gparms_fun(
            parm = allowed_parmsi, eps = eps, 
            by = comparisons_arguments$by,
            aggregate_by = aggregate_by,
            newdata = newdata
          )
        }
        
        out_sf <- do.call(rbind, list_cout) %>% data.frame() 
        if(!"parameter" %in% colnames(out_sf)) {
          out_sf <- out_sf %>% tibble::rownames_to_column(., "parameter")
        }
        out_sf <- out_sf %>% 
          dplyr::mutate(!!as.symbol('parameter') := sub("*\\.[0-9]", "", 
                                                        !!as.symbol('parameter')))
      }
    } # if(method == 'pkg') {
    
    
    # Start method = 'custom'
    
    pdrawsp_est <- NULL
    pdrawsh_est <- NULL
    pdraws_est <- NULL
    
    if(method == 'custom') {
      predictions_arguments                 <- comparisons_arguments
      predictions_arguments[['cross']]      <- NULL
      predictions_arguments[['method']]     <- NULL
      predictions_arguments[['hypothesis']] <- NULL # hypothesis evaluated later
      # From get_predictions
      #################################################
      if(call_slopes) {
        predictions_arguments[['transform']]      <- NULL
        predictions_arguments[['byfun']]          <- NULL
        
      }
      
      
      ##########################################################
      # start new layer of parameter_numeric
      ##########################################################
      if(!parameter_numeric) { # start parameter_numeric
        
        
        ########################################################################
        ########################################################################
        ########################################################################
        # This if(call_slopes) { borrowed from get_predictions
        # needed when no d1 and slopes function used to get deriv from distance
        if(call_slopes) {
          predictions_arguments[['comparison']]     <- NULL
          
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
          
        } # end if(call_slopes)  borrowed from get_predictions
        
        ########################################################################
        ########################################################################
        ########################################################################
        ########################################################################
        
        if(!is.null(predictions_arguments[['variables']])) {
          if(!is.list(eval(predictions_arguments[['variables']]))) {
            # 21.09.2024
            # In fact list does not work, variables = c('class') works
            # stop("Argument 'variables' must be a named list")
          }
        }
        
        # Imp, add xvar to the by if missing
        by <- predictions_arguments[['by']]
        
        # if(!any(grepl(xvar, by)))  by <- c(xvar, eval(by))
        
        if(isFALSE(by)) {
          by <- xvar
        } else if(!any(grepl(xvar, by))) {
          by <- c(xvar, eval(by))
        }
        
        by                            <- eval(by)
        predictions_arguments[['by']] <- by 
        
        
        if(future_splits_exe) {
          # Note that since predictions_arguments are passed to multisession, 
          # evaluate each argument
          for (i in names(predictions_arguments)) {
            predictions_arguments[[i]] <- eval(predictions_arguments[[i]])
          }
        }
        
        
        if(!future_splits_exe & callfuns) {
          if(!average) {
            if(call_predictions) {
              out <- do.call(marginaleffects::predictions, predictions_arguments)
            } 
            if(call_slopes) {
              out <- do.call(marginaleffects::slopes, predictions_arguments)
            }
          } else if(average) {
            if(call_predictions) {
              out <- do.call(marginaleffects::avg_predictions, predictions_arguments)
            } 
            if(call_slopes) {
              out <- do.call(marginaleffects::avg_slopes, predictions_arguments)
            }
          }
          
        } # if(!future_splits_exe) {
        
        
        
        if(future_splits_exe_future & callfuns) {
          if(!average) {
            myzfun <- function(x) {
              predictions_arguments[['draw_ids']] <- x
              predictions_arguments[['ndraws']]   <- NULL
              `%>%` <- bsitar::`%>%`
              if(re_expose) {
                if(verbose) message("need to expose functions for 'multisession'")
                predictions_arguments[['model']] <- 
                  bsitar::expose_model_functions(predictions_arguments[['model']])
              }
              # Re-assign appropriate function
              setenv <- predictions_arguments[['model']]$model_info$envir
              
              assign(
                o[[1]],
                predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                envir = setenv
              )
              if(call_predictions) {
                do.call(marginaleffects::predictions, predictions_arguments)
              } 
              if(call_slopes) {
                do.call(marginaleffects::slopes, predictions_arguments)
              }
              # do.call(marginaleffects::predictions, predictions_arguments)
            }
            out <-  future.apply::future_lapply(future_splits_at,
                                                future.envir = parent.frame(),
                                                future.globals = TRUE,
                                                future.seed = TRUE,
                                                FUN = myzfun)
          } else if(average) {
            myzfun <- function(x) {
              predictions_arguments[['draw_ids']] <- x
              predictions_arguments[['ndraws']]   <- NULL
              `%>%` <- bsitar::`%>%`
              if(re_expose) {
                if(verbose) message("need to expose functions for 'multisession'")
                predictions_arguments[['model']] <- 
                  bsitar::expose_model_functions(predictions_arguments[['model']])
              }
              # Re-assign appropriate function
              setenv <- predictions_arguments[['model']]$model_info$envir
              
              assign(
                o[[1]],
                predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                envir = setenv
              )
              if(call_predictions) {
                do.call(marginaleffects::avg_predictions, predictions_arguments)
              } 
              if(call_slopes) {
                do.call(marginaleffects::avg_slopes, predictions_arguments)
              }
              # do.call(marginaleffects::avg_predictions, predictions_arguments)
            }
            out <-  future.apply::future_lapply(future_splits_at,
                                                future.envir = parent.frame(),
                                                future.globals = TRUE,
                                                future.seed = TRUE,
                                                FUN = myzfun)      
          }
        } # if(future_splits_exe_future) {
        
        
        
        
        if(future_splits_exe_dofuture & callfuns) {
          `%doFuture_function%` <- doFuture::`%dofuture%`
          # somehow .options.future = list(seed = TRUE) not working, so set below
          dofutureplan <- getOption("doFuture.rng.onMisuse")
          options(doFuture.rng.onMisuse = "ignore")
          on.exit(options("doFuture.rng.onMisuse" = dofutureplan), add = TRUE)
          if(!average) {
            out <- foreach::foreach(x = 1:length(future_splits_at),
                                    .options.future = list(seed = TRUE),
                                    .options.future =
                                      list(globals = c('future_splits_at',
                                                       'setplanis',
                                                       'verbose',
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
              
              assign(
                o[[1]],
                predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]],
                envir = setenv
              )
              if(call_predictions) {
                do.call(marginaleffects::predictions, predictions_arguments)
              } 
              if(call_slopes) {
                do.call(marginaleffects::slopes, predictions_arguments)
              }
              # do.call(marginaleffects::predictions, predictions_arguments)
            }
          } else if(average) {
            out <- foreach::foreach(x = 1:length(future_splits_at),
                                    .options.future = list(seed = TRUE),
                                    .options.future =
                                      list(globals = c('future_splits_at',
                                                       'setplanis',
                                                       'verbose',
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
              
              assign(
                o[[1]],
                predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
                envir = setenv
              )
              if(call_predictions) {
                do.call(marginaleffects::avg_predictions, predictions_arguments)
              } 
              if(call_slopes) {
                do.call(marginaleffects::avg_slopes, predictions_arguments)
              }
              # do.call(marginaleffects::avg_predictions, predictions_arguments)
            }
          } 
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
            zxdraws <- out %>% marginaleffects::posterior_draws()
          }
        } else if(future_splits_exe) {
          if(callfuns) {
            if(pdrawso) {
              out <- out %>% do.call(rbind, .)
              return(out)
            }
            zxdraws <- lapply(1:length(future_splits_at), 
                              FUN = posterior_draws_function)
            
            zxdraws <- zxdraws %>% do.call(rbind, .)
            # Note that above zxdraws has drawid exact same as splits
            # but somehow, we need consecutive drawid for summarising
            zxdraws <- consecutive_drawid_function(zxdraws)
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
        
        
        
        if(setmarginals) {
          if(inherits(marginals, 'list')) {
            zxdraws <-
              {. <- lapply(1:length(marginals), 
                           marginals_list_consecutive_drawid_function)
              list2DF(lapply(setNames(seq_along(.[[1]]), names(.[[1]])), function(i)
                unlist(lapply(., `[[`, i), FALSE, FALSE)))}
            zxdraws$drawid <- cheapr::factor_(zxdraws$drawid)
          } else {
            zxdraws <- marginals
          }
        }
        
        
        
        
        
        
        
        
        
        
        by_pdraws <- by
        
        # Imp, remove xvar from the by
        by <- base::setdiff(eval(by), eval(xvar)) 
        
        
        
        
        
        getparmsx <- function(x, y, parm = NULL, xvar = NULL, draw = NULL,
                              aggregate_by = FALSE, ...) {
          
          if(data.table::is.data.table(x) | is.data.frame(x)) {
            if(is.null(xvar)) {
              stop("please specify the 'xvar' argument")
            }
            if(is.null(xvar)) {
              stop("please specify the 'draw' argument")
            }
            temx <- x
            x <- temx[[xvar]] %>% unlist() %>% as.numeric()
            y <- temx[[draw]] %>% unlist() %>% as.numeric()
          }
          
          # aggregate_by <- FALSE
          if(aggregate_by) {
            try(insight::check_if_installed(c("grDevices", "stats"), stop = FALSE, 
                                            prompt = FALSE))
            xy <- grDevices::xy.coords(x, y)
            xy <- unique(as.data.frame(xy[1:2])[order(xy$x), ])
            if(!isFALSE(by)) {
              if(ec_agg == "mean")   xy <- stats::aggregate(.~x, data=xy, 
                                                            mean, 
                                                            na.action = na.omit,
                                                            drop = TRUE)
              if(ec_agg == "median") xy <- stats::aggregate(.~x, data=xy, 
                                                            median, 
                                                            na.action = na.omit,
                                                            drop = TRUE)
            }
            x <- xy$x
            y <- xy$y
          } # if(aggregate_by) {
          
          parm_c <- list()
          pgvx <- NULL
          for (parmi in parm) {
            if('apgv' %in% parmi) {
              # prepare_data2
              parm_c[[parmi]] <- sitar::getPeak(x, y)[1] %>% ifunx_() 
            }
            if('pgv' %in% parmi) {
              parm_c[[parmi]] <- pgvx <- sitar::getPeak(x, y)[2]
            }
            if('atgv' %in% parmi) {
              # prepare_data2
              parm_c[[parmi]] <- sitar::getTakeoff(x, y)[1] %>% ifunx_()
            }
            if('tgv' %in% parmi) {
              parm_c[[parmi]] <- sitar::getTakeoff(x, y)[2]
            }
            if('acgv' %in% parmi | 'acgv' %in% parmi) {
              if(is.null(pgvx)) pgvx <- sitar::getPeak(x, y)[2]
              cgv  <- acg_velocity * pgvx
              vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
            }
            # prepare_data2
            if('acgv' %in% parmi) parm_c[[parmi]] <- x[vcgi] %>% ifunx_()
            if('cgv' %in% parmi)  parm_c[[parmi]] <- y[vcgi]
          }
          out <- parm_c %>% do.call(cbind, .) %>% data.frame()
          out
        }
        
        
        if(usedtplyr) {
          getparmsx2                     <- getparmsx
          hypothesisargs                 <- formals(getparmsx2)
          hypothesisargs[['xvar']]       <- xvar
          hypothesisargs[['draw']]       <- 'draw'
          hypothesisargs[['parm']]       <- parm
          formals(getparmsx2)            <- hypothesisargs
          drawidby                       <- c('drawid', by)
          onex0 <- zxdraws %>% dtplyr::lazy_dt() %>%
            dplyr::group_by_at(drawidby) %>% 
            dplyr::group_modify(., getparmsx2, .keep = F) %>% 
            dplyr::ungroup()
          
        } else if(usecollapse) {
          drawidby  <- c('drawid', by)
          drawidby_ <- c(drawidby, 'parameter', 'estimate')
          parmest   <- 'draw'
          
          if(any(c('apgv', 'pgv') %in% parm)) getpest <- TRUE else getpest <- FALSE
          if(any(c('atgv', 'tgv') %in% parm)) gettest <- TRUE else gettest <- FALSE
          if(any(c('acgv', 'cgv') %in% parm)) getcest <- TRUE else getcest <- FALSE
          
          if(getpest) namesp <- cbind('apgv', 'pgv') else namesp <- NULL
          if(gettest) namest <- cbind('atgv', 'tgv') else namest <- NULL
          if(getcest) namesc <- cbind('acgv', 'cgv') else namesc <- NULL
          
          getcgvfunc <- function(x, y, p = NULL, ...) {
            if(is.null(p)) {
              cgv <- acg_velocity * sitar::getPeak(x, y)[2] 
            } else {
              cgv <- acg_velocity * p
            }
            vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
            return(c(x = x[vcgi], y = y[vcgi]))
          }
          
          funx <- function(x,...) {
            if(getpest) {
              dfp <- sitar::getPeak(x[,1], x[,2]) 
              dfp[1] <- ifunx_(dfp[1]) # prepare_data2
            } else {
              dfp <- NULL
            }
            if(gettest) {
              dft <- sitar::getTakeoff(x[,1], x[,2]) 
              dft[1] <- ifunx_(dft[1]) # prepare_data2
            } else {
              dft <- NULL
            }
            if(getcest) {
              dfc <- getcgvfunc(x[,1], x[,2]) 
              dfc[1] <- ifunx_(dfc[1]) # prepare_data2
            } else {
              dfc <- NULL
            }
            cbind(c(namesp, namest, namesc), matrix(c(dfp, dft, dfc)))
          }
          
          onex0 <- zxdraws %>% 
            collapse::fgroup_by(drawidby) %>% 
            collapse::fsummarise(collapse::mctl(funx(cbind(.data[[xvar]], 
                                                           .data[[parmest]])), 
                                                names = F)) %>% 
            collapse::ftransformv(., 'V2', as.numeric) %>% 
            collapse::frename(., drawidby_) %>% 
            collapse::fsubset(., parameter %in% parm)
          
        
          
        } else {
          drawid_c <- list()
          for (drawidi in 1:nlevels(zxdraws$drawid)) {
            drawid_c[[drawidi]] <- zxdraws %>% dplyr::filter(drawid == drawidi) %>%
              dplyr::group_by_at(by) %>%
              dplyr::group_modify(., ~ getparmsx(.x[[xvar]] , .x$draw, parm = parm),
                                  .keep = TRUE) %>%
              dplyr::mutate(drawid = drawidi)
          }
          onex0 <- drawid_c %>% do.call(rbind, .) %>% data.frame()
        }
        
        
        #######################################################
        
        # This from marginaleffects does not allow na.rm
        # So use stats::quantile instead
        get_etix <- utils::getFromNamespace("get_eti", "marginaleffects")
        get_etix <- stats::quantile
        get_hdix <- utils::getFromNamespace("get_hdi", "marginaleffects")
        get_pe_ci <- function(x, draw = NULL, na.rm = TRUE, ...) {
          if(data.table::is.data.table(x) | is.data.frame(x)) {
            if(is.null(draw)) {
              stop("Please specify the 'draw' argument")
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
        
        
      } # end if(!parameter_numeric) {
      ##########################################################
      # end new layer of parameter_numeric
      ##########################################################
      
      
      #################################################################
      #################################################################
      #################################################################
      #################################################################
      #################################################################
      #################################################################
      
      # "marginaleffects"  "get_predictions" "fitted_draws"
      # Using "marginaleffects" because it gives consiatent result future T/F
      
      call_marginaleffects_get_predictions <- "marginaleffects" 
      
      
      # progressr::handlers(global = TRUE)
      # progressr::handlers("progress", "beepr")
      
      xyadj_xyv_warp_fun <- function(ApvX0, newdata, ...) {
        xmin <- NULL;
        xmax <- NULL;
        
        # p <- progressor(ApvX0 = xs)
        
        newdata.in <- newdata
        newdata <- newdata%>% 
          fastplyr::f_group_by(!!as.name(idvar)) %>% 
          fastplyr::f_slice(1) %>% fastplyr::f_ungroup()
        
        # For storing results in my_matrix returned from analyze_data()
        nlevels_idvar <- nlevels(newdata [[idvar]])
        
        # Here nrow  is the number of rows - xyadj, yyadj and vyadj
        my_matrix <- matrix(NA, nrow = 3, ncol = nlevels_idvar)
        
        
        ApvX0[['drawid']] <- eval_draw_ids
        
        ApvX0_for_list_c <- tibble::as_tibble(ApvX0) %>% 
          fastplyr::f_select(-'parameter') %>% 
          fastplyr::f_rename('age' = 'estimate') %>% 
          dplyr::mutate(tempjoinby = 'tempjoinby') %>% 
          dplyr::relocate('age') 
        
        newdata_for_list_c <-tibble::as_tibble(newdata) %>%
          fastplyr::f_select(-'age') %>% 
          dplyr::mutate(tempjoinby = 'tempjoinby') 
        
        newdata_list_c <- list()
        for (i in 1:nrow(ApvX0)) {
          newdata_list_c[[i]] <- fastplyr::f_left_join(newdata_for_list_c, 
                                                       ApvX0_for_list_c[i,], 
                                                       by = "tempjoinby") %>% 
            fastplyr::f_select(-'tempjoinby') 
        }
        
        
        xyadj_curves_args <- list()
        xyadj_curves_args[['model']]  <- model
        xyadj_curves_args[['ndraws']] <- ndraws
        xyadj_curves_args[['tomean']] <- FALSE
        xyadj_curves_args[['get_dv']] <- TRUE
        
        
        fitted_draws_args <- list()
        fitted_draws_args[['model']] <- model
        fitted_draws_args[['ndraws']] <- ndraws
        fitted_draws_args[['draw_ids']] <- draw_ids
        fitted_draws_args[['re_formula']] <- NULL
        fitted_draws_args[['summary']] <- FALSE
        fitted_draws_args[['newdata']] <- newdata 
        fitted_draws_args_d0 <- fitted_draws_args
        fitted_draws_args_d1 <- fitted_draws_args
        fitted_draws_args_d0[['deriv']] <- 0
        fitted_draws_args_d1[['deriv']] <- 1
        
        
        get_predictions_args <- predictions_arguments
        get_predictions_args[['re_formula']] <- NULL
        get_predictions_args[['peak']] <- NULL
        get_predictions_args[['takeoff']] <- NULL
        get_predictions_args[['trough']] <- NULL
        get_predictions_args[['acgv']] <- NULL
        get_predictions_args[['newdata_fixed']] <- 0
        get_predictions_args[['pdrawsp']] <-  "return"
        get_predictions_args[['newdata']] <- newdata 
        
        get_predictions_args[['comparison']] <- NULL
        
        if(call_marginaleffects_get_predictions == "marginaleffects") {
          get_predictions_args[['by']] <- xvar
          get_predictions_args[['newdata_fixed']] <- NULL
          get_predictions_args[['pdrawsp']] <- NULL
          get_predictions_args[['deriv']] <- NULL
        } # if(call_marginaleffects_get_predictions == "marginaleffects") {
        
        
        if(call_marginaleffects_get_predictions == "fitted_draws") {
          get_predictions_args[['summary']] <- FALSE
        }
        
        get_predictions_args_d0 <- get_predictions_args
        get_predictions_args_d1 <- get_predictions_args
        get_predictions_args_d0[['deriv']] <- 0
        get_predictions_args_d1[['deriv']] <- 1
        
        
        if(call_marginaleffects_get_predictions == "marginaleffects") {
          get_predictions_args_d0[['deriv']] <- NULL
          get_predictions_args_d1[['deriv']] <- NULL 
          get_predictions_args_d0[['by']] <- NULL
          get_predictions_args_d0[['variables']] <- NULL
        }
        
        if(call_marginaleffects_get_predictions == "get_predictions") {
          get_predictions_args_d0[['by']] <- NULL
          get_predictions_args_d0[['variables']] <- NULL
          get_predictions_args_d1[['by']] <- NULL
        } # if(call_marginaleffects_get_predictions == "get_predictions") {
        
        
        if(call_marginaleffects_get_predictions == "marginaleffects" |
           call_marginaleffects_get_predictions == "get_predictions") {
          call_predictions2 <- TRUE
          call_slopes2      <- FALSE
          if(available_d1) {
            call_predictions2 <- TRUE
            call_slopes2      <- FALSE
            post_processing_checks_args[['deriv']]    <- 0
            o.2 <- o.2.0 <- do.call(post_processing_checks, 
                                    post_processing_checks_args)
            post_processing_checks_args[['deriv']]    <- 1
            o.2 <- o.2.1 <- do.call(post_processing_checks, 
                                    post_processing_checks_args)
          } else if(!available_d1) {
            call_predictions2 <- FALSE
            call_slopes2      <- TRUE
            post_processing_checks_args[['deriv']]    <- 0
            o.2 <- o.2.0 <- do.call(post_processing_checks, 
                                    post_processing_checks_args)
          }
        } # if(call_marginaleffects_get_predictions == "marginaleffects") {
        
        
        
        
        analyze_data <- function(data, ...) {
          xyadj_curves_args[['newdata']] <- data
          xyadj_curves_args[['draw_ids']] <- data[['drawid']] [1]
          xyadj_dv <- CustomDoCall(xyadj_curves, xyadj_curves_args)
          
          xyadj_dv <- 
            apply(xyadj_dv, 2, 
                  xyadj_curves_args$model$model_info$xfuntransform2) 
          
          # i <-  data[['loopcounter']] [1]
          
          fitted_draws_args_d0[['draw_ids']] <- data[['drawid']] [1]
          fitted_draws_args_d1[['draw_ids']] <- data[['drawid']] [1]
          fitted_draws_args_d0[['newdata']] [['age']]  <- xyadj_dv 
          fitted_draws_args_d1[['newdata']] [['age']]  <- xyadj_dv 
          
          get_predictions_args_d0[['draw_ids']] <- data[['drawid']] [1]
          get_predictions_args_d1[['draw_ids']] <- data[['drawid']] [1]
          get_predictions_args_d0[['newdata']] [['age']]  <- xyadj_dv 
          get_predictions_args_d1[['newdata']] [['age']]  <- xyadj_dv 
          
          
          if(call_marginaleffects_get_predictions == "marginaleffects") {
            if(!average) {
              if(call_predictions2) {
                assign(o[[1]], model$model_info[['exefuns']][[o.2.0[[2]]]], 
                       envir = envir)
                yyadj_dv <- do.call(marginaleffects::predictions, 
                                    get_predictions_args_d0)
                
                assign(o[[1]], model$model_info[['exefuns']][[o.2.1[[2]]]], 
                       envir = envir)
                
                vyadj_dv <- do.call(marginaleffects::predictions, 
                                    get_predictions_args_d1)
              }
              if(call_slopes2) {
                assign(o[[1]], model$model_info[['exefuns']][[o.2.0[[2]]]], 
                       envir = envir)
                
                yyadj_dv <- do.call(marginaleffects::predictions, 
                                    get_predictions_args_d0)
                
                assign(o[[1]], model$model_info[['exefuns']][[o.2.0[[2]]]], 
                       envir = envir)
                
                vyadj_dv <- do.call(marginaleffects::slopes, 
                                    get_predictions_args_d1)
              }
            } else if(average) {
              if(call_predictions2) {
                assign(o[[1]], model$model_info[['exefuns']][[o.2.0[[2]]]], 
                       envir = envir)
                yyadj_dv <- do.call(marginaleffects::avg_predictions, 
                                    get_predictions_args_d0)
                
                assign(o[[1]], model$model_info[['exefuns']][[o.2.1[[2]]]], 
                       envir = envir)
                
                vyadj_dv <- do.call(marginaleffects::avg_predictions, 
                                    get_predictions_args_d1)
              }
              if(call_slopes2) {
                assign(o[[1]], model$model_info[['exefuns']][[o.2.0[[2]]]], 
                       envir = envir)
                
                yyadj_dv <- do.call(marginaleffects::avg_predictions, 
                                    get_predictions_args_d0)
                
                assign(o[[1]], model$model_info[['exefuns']][[o.2.0[[2]]]], 
                       envir = envir)
                
                vyadj_dv <- do.call(marginaleffects::avg_slopes, 
                                    get_predictions_args_d1)
              }
            }
          } # if(call_marginaleffects_get_predictions == "marginaleffects") {
          
          
          
          
          
          if(call_marginaleffects_get_predictions == "get_predictions") {
            if(call_predictions2) {
              assign(o[[1]], model$model_info[['exefuns']][[o.2.0[[2]]]], 
                     envir = envir)
              yyadj_dv <- do.call(get_predictions, 
                                  get_predictions_args_d0)
              
              assign(o[[1]], model$model_info[['exefuns']][[o.2.1[[2]]]], 
                     envir = envir)
              
              vyadj_dv <- do.call(get_predictions, 
                                  get_predictions_args_d1)
            }
            if(call_slopes2) {
              assign(o[[1]], model$model_info[['exefuns']][[o.2.0[[2]]]], 
                     envir = envir)
              
              yyadj_dv <- do.call(get_predictions, 
                                  get_predictions_args_d0)
              
              vyadj_dv <- do.call(get_predictions, 
                                  get_predictions_args_d1)
            }
          } # if(call_marginaleffects_get_predictions == "get_predictions") {
          
          
          
          # if(call_marginaleffects_get_predictions == "get_predictions") {
          #   yyadj_dv <- CustomDoCall(get_predictions, get_predictions_args_d0)
          #   vyadj_dv <- CustomDoCall(get_predictions, get_predictions_args_d1)
          # } # if(call_marginaleffects_get_predictions == "get_predictions") {
          # 
          
          if(call_marginaleffects_get_predictions == "fitted_draws") {
            yyadj_dv <- CustomDoCall(fitted_draws, get_predictions_args_d0)
            vyadj_dv <- CustomDoCall(fitted_draws, get_predictions_args_d1)
          }
          
          
          if(call_marginaleffects_get_predictions == "marginaleffects") {
            yyadj_dv <- yyadj_dv[["estimate"]]
            vyadj_dv <- vyadj_dv[["estimate"]]
          }
          if(call_marginaleffects_get_predictions == "get_predictions") {
            yyadj_dv <- yyadj_dv[["estimate"]]
            vyadj_dv <- vyadj_dv[["estimate"]]
          }
          
          
          
          my_matrix[1, ] <- xyadj_dv
          my_matrix[2, ] <- yyadj_dv
          my_matrix[3, ] <- vyadj_dv
          my_matrix
        } # end analyze_data
        
        if(future) {
          results_list <- 
            future.apply::future_lapply(newdata_list_c, 
                                        future.envir = new.env(),
                                        future.globals = F,
                                        # future.packages = 'bsitar',
                                        future.seed = TRUE,
                                        FUN = analyze_data) 
          
        } else {
          results_list <- lapply(newdata_list_c, analyze_data) 
        }
        
        
        xyvyadj_rows <- do.call(cbind, results_list)
        
        
        xyadj_rows <- matrix(xyvyadj_rows[1,], nrow = length(results_list), 
                             byrow = T)
        yyadj_rows <- matrix(xyvyadj_rows[2,], nrow = length(results_list), 
                             byrow = T)
        vyadj_rows <- matrix(xyvyadj_rows[3,], nrow = length(results_list), 
                             byrow = T)
        
        xyadj_rows <- apply(xyadj_rows, 2, model$model_info$ixfuntransform2)
        
        
        # When only one draw, apply converts xyadj_rows to as vector
        if(length(eval_draw_ids) == 1) {
          xyadj_rows <- matrix(xyadj_rows, nrow = 1)
        }
        
        xyadj_summary <- brms::posterior_summary(xyadj_rows)
        yyadj_summary <- brms::posterior_summary(yyadj_rows)
        vyadj_summary <- brms::posterior_summary(vyadj_rows)
        
        
        
        # data frame and add idavr
        xyadj_summary <- xyadj_summary %>% data.frame() %>% 
          dplyr::mutate(!! as.name(idvar) := newdata[[idvar]]) %>% 
          dplyr::relocate(!! as.name(idvar))
        
        # add missing xvar within the range
        xyadj_summary <- 
          dplyr::inner_join(xyadj_summary,
                            newdata.in %>%
                              dplyr::select(dplyr::all_of(c(idvar, 
                                                            xvar))) %>%
                              dplyr::nest_by(dplyr::pick(1)) %>%
                              dplyr::summarise(xmin = min(dplyr::pick(1)),
                                               xmax = max(dplyr::pick(1)),
                                               .groups = 'drop'),
                            by = idvar) %>%
          dplyr::mutate(missing = !between(dplyr::pick(2) %>% 
                                             dplyr::pull(), xmin, 
                                           xmax)) %>%
          dplyr::select(-c(xmin, xmax))
        
        
        # data frame and add idavr
        yyadj_summary <- yyadj_summary %>% data.frame() %>% 
          dplyr::mutate(!! as.name(idvar) := newdata[[idvar]]) %>% 
          dplyr::relocate(!! as.name(idvar))
        
        yyadj_summary$missing <- xyadj_summary$missing
        
        vyadj_summary <- vyadj_summary %>% data.frame() %>% 
          dplyr::mutate(!! as.name(idvar) := newdata[[idvar]]) %>% 
          dplyr::relocate(!! as.name(idvar))
        
        vyadj_summary$missing <- xyadj_summary$missing
        
        
        out <- list()
        out[['age']]      <- xyadj_summary
        out[['distance']] <- yyadj_summary
        out[['velocity']] <- vyadj_summary
        
        return(out)
      } # end xyadj_xyv_warp_fun
      
      
      
      # parameter_numeric
      if(parameter_numeric) {
        value <- NULL;
        onex0 <- tibble::as_tibble(1:length(eval_draw_ids)) %>% 
          dplyr::mutate(estimate = parameter) %>% 
          dplyr::mutate(parameter = "xxx") %>% # dummy
          dplyr::rename(drawid = value) %>% 
          dplyr::relocate(drawid, parameter, estimate)
      }
      
      
      
      
      xyadj_xyv <- xyadj_xyv_warp_fun(onex0, model$model_info$bgmfit.data)
      
      return(xyadj_xyv)
      
      
      #################################################################
      #################################################################
      #################################################################
      #################################################################
      #################################################################
      #################################################################
      #################################################################
      #################################################################
      
      
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
      
      
      
      
      if(!isFALSE(pdraws)) {
        selectchoicesr <- c("return", "add", "returns", "adds") 
        checkmate::assert_choice(pdraws, choices = selectchoicesr)
        if(pdraws == 'return') {
          return(zxdraws)
        } else if(pdraws == 'add') {
          pdraws_est <- zxdraws
        } else {
          setdrawidparm <- c(by_pdraws)
          namesx <- c('estimate', 'conf.low', 'conf.high')
          setdrawidparm_ <- c(setdrawidparm, namesx)
          # No need to summarise again because estimate already is summary
          # if(usecollapse) {
          #   what_summary <- 'draw' # 'estimate' 'draw'
          #   setdrawidparm <- c(by_pdraws)
          #   namesx <- c('estimate', 'conf.low', 'conf.high')
          #   setdrawidparm_ <- c(setdrawidparm, namesx)
          # 
          #   zxdraws_summary <-
          #     zxdraws %>% collapse::fgroup_by(setdrawidparm) %>%
          #     collapse::fsummarise(collapse::mctl(
          #       get_pe_ci_collapse(.data[[what_summary]]))
          #     ) %>%
          #     collapse::ftransformv(., 'V2', as.numeric) %>%
          #     collapse::frename(., setdrawidparm_)
          # 
          #   row.names(zxdraws_summary) <- NULL
          # }
          
          if(usedtplyr) {
            zxdraws_summary <- zxdraws %>% dplyr::filter(., drawid == 1) %>%
              dplyr::select(., dplyr::all_of(setdrawidparm_))
          } else if(usecollapse) {
            zxdraws_summary <- zxdraws %>% collapse::fsubset(., drawid == 1) %>%
              collapse::fselect(., setdrawidparm_)
          } else {
            zxdraws_summary <- zxdraws %>% dplyr::filter(., drawid == 1) %>%
              dplyr::select(., dplyr::all_of(setdrawidparm_))
          }
          
          
          if(pdraws == 'returns') return(zxdraws_summary)
          if(pdraws == 'adds') pdraws_est <- zxdraws_summary
        }
      } # if(!isFALSE(pdraws)) {
      
      
      if(isFALSE(constrats_by)) {
        constrats_by <- NULL
      } else if(!isFALSE(constrats_by)) {
        if(is.null(constrats_by)) {
          if(is.null(hypothesis)) {
            constrats_by <- NULL
          } else if(!is.null(hypothesis)) {
            if(!isFALSE(by)) constrats_by <- setdiff(by, xvar) 
          }
        } else if(!is.null(constrats_by)) {
          constrats_by <- constrats_by
        }
      }
      
      
      if(!is.null(constrats_by)) {
        if(!is.character(constrats_by)) 
          stop("The 'constrats_by' argument should be a character string")
        for (axi in 1:length(constrats_by)) {
          caxi <- constrats_by[axi]
          if(!is.character(caxi)) {
            stop("The 'constrats_by' argument '", caxi, " should be a character")
          }
          if(!caxi %in% by) {
            stop("The 'constrats_by' argument '", caxi, "' is not available in",
                 " the 'by' argument.",
                 "\n ", 
                 " Note that '", caxi, "' must be included in the 'by' argument.",
                 "\n ",
                 " The current 'by' argument includes:",
                 "\n ",
                 collapse_comma(by)
            )
          }
        }
      }
      
      
      if(isFALSE(constrats_at)) {
        constrats_at <- NULL
      } else if(!isFALSE(constrats_at)) {
        if(is.null(constrats_at)) {
          constrats_at <- list()
          for(byi in by) {
            if(is.numeric(constrats_at[[byi]])) {
              constrats_at[[byi]] <- 'unique'
            }
          }
          if(length(constrats_at) == 0) constrats_at <- NULL
        }
      }
      
      # For hypothesis
      groupvarshyp1 <- c('drawid')
      groupvarshyp2 <- c('term')
      if(!is.null(constrats_at)) {
        for (caxi in names(constrats_at)) {
          if(!caxi %in% names(onex0)) {
            stop("Variable '", caxi, ". specified in 'constrats_at' is not in",
                 "\n ", 
                 " the 'by' arg. Note that ", caxi, " should also be included",
                 "\n ", 
                 " in the 'by' argument. The current 'by' argument includes:", 
                 "\n ",
                 collapse_comma(by)
            )
          }
          allowed_char_constrats_at <- c('max', 'min', 'unique', 'range')
          if(is.character(constrats_at[[caxi]])) {
            if(length(constrats_at[[caxi]]) > 1) {
              stop(caxi, " specified in 'constrats_at' as character should be a 
                 single character:",
                   "\n ", 
                   collapse_comma(allowed_char_constrats_at)
              )
            }
            if(!constrats_at[[caxi]] %in% allowed_char_constrats_at) {
              stop(constrats_at[[caxi]], " specified in 'constrats_at' as 
                 character should be one of",
                   "\n ", 
                   collapse_comma(allowed_char_constrats_at)
              )
            }
            getatval<- paste0(constrats_at[[caxi]], "(", 
                              "onex0", "[['", caxi, "']]", ")" )
            getatval <- eval(parse(text = getatval))                
          } # if(is.character(constrats_at[[caxi]])) {
          
          if(!is.character(constrats_at[[caxi]])) {
            getatval <- constrats_at[[caxi]]
          }
          constrats_at[[caxi]] <- getatval
          groupvarshyp1 <- c(caxi, groupvarshyp1)
          groupvarshyp2 <- c(caxi, groupvarshyp2)
        } # for (caxi in names(constrats_at)) {
        
        for (caxi in names(constrats_at)) {
          onex1 <- base::subset(onex0, onex0[[caxi]] %in% constrats_at[[caxi]])
          if(nrow(onex1) == 0) {
            stop(caxi, " specified in 'constrats_at' has resulted in zero rows",
                 "\n ", 
                 ""
            )
          }
        }
      } # if(!is.null(constrats_at)) {
      
      
      if(is.null(constrats_at)) {
        onex1 <- onex0
      }
      
      #########################33
      if(isFALSE(constrats_subset)) {
        constrats_subset <- NULL
      } 
      
      
      if(!is.null(constrats_subset)) {
        for (caxi in names(constrats_subset)) {
          if(!caxi %in% names(onex0)) {
            stop("Variable '", caxi, ". specified in 'constrats_subset' is not ",
                 " avaialble in the 'by' argument.",
                 "\n ", 
                 " Please include ", caxi, " in the 'by' argument.",
                 "\n ", 
                 " The current 'by' argument includes:", 
                 "\n ",
                 collapse_comma(by)
            )
          }
          getatval <- constrats_subset[[caxi]]
          constrats_subset[[caxi]] <- getatval
          groupvarshyp1 <- c(caxi, groupvarshyp1)
          groupvarshyp2 <- c(caxi, groupvarshyp2)
        } # for (caxi in names(constrats_subset)) {
        
        
        for (caxi in names(constrats_subset)) {
          onex1 <- 
            base::subset(onex1, onex1[[caxi]] %in% constrats_subset[[caxi]])
          if(nrow(onex1) == 0) {
            stop(caxi, " specified in 'constrats_subset' resulted in zero rows",
                 "\n ", 
                 ""
            )
          }
        }
      } # if(!is.null(constrats_subset)) {
      
      if(usedtplyr) {
        get_pe_ci2 <- get_pe_ci
        hypothesisargs <- formals(get_pe_ci2)
        hypothesisargs[['x']]          <- as.name('x')
        hypothesisargs[['draw']]       <- 'draw'
        hypothesisargs[['...']]        <- NULL
        hypothesisargs[['na.rm']]      <- TRUE
        hypothesisargs <- base::append(hypothesisargs, as.name('...') , 
                                       after = 1)
        names(hypothesisargs)[2] <- '...'
        formals(get_pe_ci2) <- hypothesisargs
        parmi_estimate <- 'estimate'
        
        setdrawid     <- c('drawid', by)
        setdrawidparm <- c(by, 'parameter')
        
        out_sf_and_later_hy <-
          onex1 %>% dtplyr::lazy_dt() %>% 
          tidyr::pivot_longer(!dplyr::all_of(setdrawid), names_to = 'parameter', 
                              values_to = "draw")
        
        out_sf <- out_sf_and_later_hy %>% 
          dplyr::select(-dplyr::all_of('drawid')) %>% 
          dplyr::group_by_at(setdrawidparm) %>% 
          dplyr::group_modify(., get_pe_ci2, .keep = F) %>% 
          dtplyr::lazy_dt() %>% 
          dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, 
                                       fixed = TRUE))
      } else if(usecollapse) {
        if(is.null(bys)) {
          setdrawid     <- c('drawid', by)
          setdrawidparm <- c(by, 'parameter')
        } else if(!is.null(bys)) {
          setdrawid     <- c('drawid', bys)
          setdrawidparm <- c(bys, 'parameter')
        }
        
        
        namesx <- c('estimate', 'conf.low', 'conf.high')
        setdrawidparm_ <- c(setdrawidparm, namesx)
        
        out3 <-
          onex1 %>% collapse::fgroup_by(setdrawidparm) %>% 
          collapse::fsummarise(collapse::mctl(
            get_pe_ci_collapse(.data[['estimate']]))
          ) %>% 
          collapse::ftransformv(., 'V2', as.numeric) %>% 
          collapse::frename(., setdrawidparm_) 
        
        row.names(out3) <- NULL
        out_sf <- out3
        
      } else {
        summary_c <- list()
        for (parmi in parm) {
          summary_c[[parmi]] <- onex1 %>% # out2 %>%
            dplyr::reframe(
              dplyr::across(c(dplyr::all_of(parmi)), get_pe_ci, 
                            .unpack = TRUE),
              .by = dplyr::all_of(!! by)
            ) %>% dplyr::rename_with(., ~ gsub(paste0(parmi, "_"), "", .x, 
                                               fixed = TRUE)) %>% 
            dplyr::mutate(parameter = parmi) %>% 
            dplyr::relocate(parameter)
        }
        out3 <- summary_c %>% do.call(rbind, .) %>% data.frame()
        row.names(out3) <- NULL
        out_sf <- out3
      }
      
      
      # Hypothesis
      if(!is.null(hypothesis)) {
        get_hypothesis_x_modify <- function(.x, hypothesis, by, draws, ...) {
          get_hypothesis_x(x = .x, hypothesis = hypothesis, by = by, 
                           draws = draws)
        }
        if(length(tibble::as_tibble(onex1)) > 25) {
          # if(nrow(onex1) > 25) {
          if(is.null(constrats_at)) {
            message("Note that the 'marginaleffects' package does not allow" ,
                    "\n",
                    "hypothesis argument when estimates rows are more than 25",
                    "\n",
                    "To avoid this issue, you can use 'constrats_at' argument",
                    "\n"
            )
          }
        }
        
        
        
        if(usedtplyr) {
          parmi_estimate <- 'estimate'
          get_hypothesis_x_modifyx2 <- get_hypothesis_x
          hypothesisargs <- formals(get_hypothesis_x_modifyx2)
          hypothesisargs[['x']]          <- as.name('x')
          hypothesisargs[['hypothesis']] <- hypothesis
          hypothesisargs[['draws']]      <- as.name(parmi_estimate)
          hypothesisargs[['by']]         <- constrats_by
          hypothesisargs[['newdata']]    <- NULL
          hypothesisargs <- base::append(hypothesisargs, as.name('...') , 
                                         after = 1)
          names(hypothesisargs)[2] <- '...'
          formals(get_hypothesis_x_modifyx2) <- hypothesisargs
          
          get_pe_ci2 <- get_pe_ci
          hypothesisargs <- formals(get_pe_ci2)
          hypothesisargs[['x']]          <- as.name('x')
          hypothesisargs[['draw']]       <- parmi_estimate
          hypothesisargs[['...']]        <- NULL
          hypothesisargs[['na.rm']]      <- TRUE
          hypothesisargs <- base::append(hypothesisargs, as.name('...') , 
                                         after = 1)
          names(hypothesisargs)[2] <- '...'
          formals(get_pe_ci2) <- hypothesisargs
          
          groupvarshyp1 <- c("drawid", 'parameter')
          groupvarshyp2 <- c('term', 'parameter')
          
          out_sf_hy <-
            # onex1 %>% dtplyr::lazy_dt() %>%
            out_sf_and_later_hy %>% 
            dplyr::group_by_at(groupvarshyp1) %>% 
            dplyr::mutate(!! parmi_estimate := eval(parse(text = 'draw'))) %>% 
            dplyr::group_modify(., get_hypothesis_x_modifyx2, .keep = F) %>% 
            dplyr::ungroup() %>% 
            dplyr::select(-dplyr::all_of('drawid')) %>% 
            dtplyr::lazy_dt(.) %>%
            dplyr::group_by_at(groupvarshyp2) %>% 
            dplyr::group_modify(., get_pe_ci2, .keep = F) %>% 
            dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, 
                                         fixed = TRUE)) 
        } else if(usecollapse) {
          get_hypothesis_x_fx <- function(x,...) {
            get_hypothesis_x(x,
                             hypothesis = hypothesis,
                             by = by, 
                             newdata = NULL,
                             draws = 'estimate'
            ) %>% as.matrix()
          }
          
          
          
          setdrawidh     <- c('drawid', 'parameter')
          setdrawidh_  <- c(setdrawidh, 'term', 'estimate')
          temhyy <- 
            onex1 %>% collapse::fgroup_by( setdrawidh ) %>% 
            collapse::fsummarise(collapse::mctl(get_hypothesis_x_fx(.data))) %>% 
            collapse::ftransformv(., 'V2', as.numeric) %>% 
            collapse::frename(., setdrawidh_)
          
          
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
          
          
          
          setdrawidparmh <- c('term', 'parameter')
          namesx <- c('estimate', 'conf.low', 'conf.high')
          setdrawidparm_ <- c(setdrawidparmh, namesx)
          
          out_sf_hy <-
            temhyy %>% collapse::fgroup_by(setdrawidparmh) %>% 
            collapse::fsummarise(collapse::mctl(
              get_pe_ci_collapse(.data[['estimate']]))
            ) %>% 
            collapse::frename(., setdrawidparm_) 
          
          row.names(out_sf_hy) <- NULL
          
        } else {
          parmi_ci_c <- list()
          for (parmi in parm) {
            parmi_estimate <- 'estimate'
            measurevar <- parmi
            parmi_ci_c[[parmi]] <-
              onex1 %>% dplyr::group_by_at(groupvarshyp1) %>% 
              dplyr::mutate(!! parmi_estimate := 
                              eval(parse(text = measurevar))) %>% 
              dplyr::group_modify(., ~get_hypothesis_x_modify(.x,
                                                              hypothesis = 
                                                                hypothesis, 
                                                              by = constrats_by, 
                                                              draws = estimate
              ), .keep = F) %>% 
              dplyr::ungroup() %>% 
              dplyr::reframe(
                dplyr::across(c(dplyr::all_of(parmi_estimate)), get_pe_ci, 
                              .unpack = TRUE),
                .by = dplyr::all_of(!! groupvarshyp2)
              ) %>%
              dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, 
                                           fixed = TRUE)) %>% 
              dplyr::mutate(!! 'parameter' := parmi) %>% 
              dplyr::relocate(dplyr::all_of('parameter'))
          }
          out_sf_hy <- parmi_ci_c %>% do.call(rbind, .) %>% data.frame()
        }
        
        
      } # if(!is.null(hypothesis)) {
    } # if(method == 'custom') {
    
    
    
    
    #######################################
    
    if(!isFALSE(by)) {
      byarrange <- by 
    } else if(isFALSE(by)) {
      byarrange <- NULL
    } 
    
    if(!is.null(byarrange)) {
      if(!is.null(bys)) byarrange <- bys else byarrange <- by
    }
    
    # if(!byarrange) byarrange <- NULL
    
    if(length(byarrange) != 0) {
      out_sf <- out_sf %>% data.frame() %>% 
        dplyr::relocate(dplyr::all_of('parameter')) %>% 
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                    ~ round(., digits = digits))) %>% 
        dplyr::arrange(match(parameter, parm), !!as.name(byarrange)) %>% 
        dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
        data.frame()
    } else if(length(byarrange) == 0) {
      out_sf <- out_sf %>% data.frame() %>% 
        dplyr::relocate(dplyr::all_of('parameter')) %>% 
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                    ~ round(., digits = digits))) %>% 
        dplyr::arrange(match(parameter, parm)) %>% 
        dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
        data.frame()
    }
    
    
    
    if(!is.null(out_sf_hy)) {
      if(usecollapse) {
        out_sf_hy <- out_sf_hy %>% data.frame() %>% 
          dplyr::arrange(match(parameter, parm)) %>% 
          dplyr::relocate(dplyr::all_of('parameter')) %>% 
          dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
          dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                      ~ round(., digits = digits)))
      } else {
        out_sf_hy <- out_sf_hy %>% data.frame() %>% 
          dplyr::relocate(dplyr::all_of('parameter')) %>% 
          dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                      ~ round(., digits = digits)))
      }
    }
    
    
    
    if(!is.null(pdraws_est)) {
      if(usecollapse) {
        pdraws_est <- pdraws_est %>% data.frame() %>% 
          dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                      ~ round(., digits = digits)))
      } else {
        pdraws_est <- pdraws_est %>% data.frame() %>% 
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
    
    
    ##############################################################
    ##############################################################
    # prepare_data2
    itransform_set <- get_itransform_call(itransform)
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
      # 'pdraws_est' is not threre for get_predictions and get_comparisons
      if(!is.null(pdraws_est)) {
        pdraws_est <- prepare_transformations(data = pdraws_est, model = model,
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
      
      remove_cols_ <- c('tmp_idx', 'predicted_lo', 
                        'predicted_hi', 'predicted')
      if(!is.null(full.args$hypothesis) | !is.null(full.args$equivalence)) {
        remove_cols_ <- remove_cols_
      }
      if(is.null(full.args$hypothesis) && is.null(full.args$equivalence)) {
        remove_cols_ <- c('term', 'contrast', remove_cols_)
      }
      out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]
      colnames(out_sf) <- sub("(.)", "\\U\\1", colnames(out_sf), perl = TRUE)
      row.names(out_sf) <- NULL
      attr(out_sf$Parameter, "dimnames") <- NULL
      
      if(!is.null(out_sf_hy)) {
        out_sf_hy <- out_sf_hy %>% 
          dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
          dplyr::rename(!!as.symbol(set_names_[1]) := 
                          dplyr::all_of('estimate')) %>% 
          dplyr::rename(!!as.symbol(set_names_[2]) := 
                          dplyr::all_of('conf.low')) %>% 
          dplyr::rename(!!as.symbol(set_names_[3]) := 
                          dplyr::all_of('conf.high')) %>% 
          dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
          data.frame()
      } # if(!is.null(out_sf_hy)) {
      
      if(!is.null(pdraws_est)) {
        pdraws_est <- pdraws_est %>% 
          # dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
          dplyr::rename(!!as.symbol(set_names_[1]) := 
                          dplyr::all_of('estimate')) %>% 
          # For pdrawsp_est/pdrawsh_est, there are no conf columns, only estimates
          dplyr::rename(!!as.symbol(set_names_[2]) :=
                          dplyr::all_of('conf.low')) %>%
          dplyr::rename(!!as.symbol(set_names_[3]) :=
                          dplyr::all_of('conf.high')) %>%
          dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
          data.frame()
      } # if(!is.null(pdraws_est)) {
      
      if(!is.null(pdrawsp_est)) {
        pdrawsp_est <- pdrawsp_est %>% 
          # dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
          dplyr::rename(!!as.symbol(set_names_[1]) := 
                          dplyr::all_of('estimate')) %>% 
          # For pdrawsp_est/pdrawsh_est, there are no conf columns, only estimates
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
          # For pdrawsp_est/pdrawsh_est, there are no conf columns, only estimates
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
    if(!is.null(out_sf)) {
      out_sf <- as_factor_as_character_factor_df(out_sf)
    }
    if(!is.null(out_sf_hy)) {
      out_sf_hy <- as_factor_as_character_factor_df(out_sf_hy)
    }
    if(!is.null(pdraws_est)) {
      pdraws_est <- as_factor_as_character_factor_df(pdraws_est)
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
    
    if(!is.null(pdraws_est)) {
      out[['pdraws_est']] <- pdraws_est %>% dplyr::ungroup()
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


#' @noRd
#' @exportS3Method modelbased_growthparameters_call bgmfit
modelbased_growthparameters_call <- function(model, ...) {
  UseMethod("modelbased_growthparameters_call")
}


