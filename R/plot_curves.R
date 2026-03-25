

#' @title Plot growth curves for the Bayesian SITAR model
#'
#' @description The \strong{plot_curves()} function visualizes six different
#'   types of growth curves using the \pkg{ggplot2} package. Additionally, it
#'   allows users to create customized plots from the data returned as a
#'   \code{data.frame}. For an alternative approach, the [get_predictions()]
#'   function can be used, which not only estimates adjusted curves but also
#'   enables comparison across groups using the \code{hypotheses} argument.
#'
#' @details The \strong{plot_curves()} function is a generic tool for
#'   visualizing the following six curves:
#'   - Population average distance curve
#'   - Population average velocity curve
#'   - Individual-specific distance curves
#'   - Individual-specific velocity curves
#'   - Unadjusted individual growth curves (i.e., observed growth curves)
#'   - Adjusted individual growth curves (adjusted for the model-estimated 
#'   random effects)
#'   
#'   Internally, \strong{plot_curves()} calls the [growthparameters()] function
#'   to estimate and summarize the distance and velocity curves, as well as to
#'   compute growth parameters such as the age at peak growth velocity (APGV).
#'   The function also calls [fitted_draws()] or [predict_draws()] to make
#'   inferences based on posterior draws. As a result, \strong{plot_curves()}
#'   can plot either fitted or predicted curves. For more details, see
#'   [fitted_draws()] and [predict_draws()] to understand the difference between
#'   fitted and predicted values.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @param opt A character string containing one or more of the following
#'   plotting options:
#'   - 'd': Population average distance curve
#'   - 'v': Population average velocity curve
#'   - 'D': Individual-specific distance curves
#'   - 'V': Individual-specific velocity curves
#'   - 'u': Unadjusted individual-specific distance curves
#'   - 'a': Adjusted individual-specific distance curves (adjusted for 
#'   random effects)
#'   
#'   Note that 'd' and 'D' cannot be specified simultaneously, nor can 'v' and
#'   'V'. Other combinations are allowed, e.g., 'dvau', 'Dvau', 'dVau', etc.
#'
#' @param apv A logical value (default \code{FALSE}) indicating whether to
#'   calculate and plot the age at peak velocity (APGV) when \code{opt} includes
#'   'v' or 'V'.
#'  
#' @param bands A character string containing one or more of the following
#'   options, or \code{NULL} (default), indicating if CI bands should be plotted
#'   around the curves:
#'   - 'd': Band around the distance curve
#'   - 'v': Band around the velocity curve
#'   - 'p': Band around the vertical line denoting the APGV parameter
#'   
#'   The \code{'dvp'} option will include CI bands for distance and velocity
#'   curves, and the APGV.
#'  
#' @param conf A numeric value (default \code{0.95}) specifying the confidence
#'   interval (CI) level for the bands. See [bsitar::growthparameters()] for
#'   more details.
#'
#' @param trim A numeric value (default \code{0}) indicating the number of long
#'   line segments to be excluded from the plot when the option 'u' or 'a' is
#'   selected. See [sitar::plot.sitar] for further details.
#'
#' @param layout A character string defining the plot layout. The default
#'   \code{'single'} layout overlays distance and velocity curves on a single
#'   plot when \code{opt} includes combinations like \code{'dv'}, \code{'Dv'},
#'   \code{'dV'}, or \code{'DV'}. The alternative layout option \code{'facet'}
#'   uses \code{facet_wrap} from \pkg{ggplot2} to map and draw plots when
#'   \code{opt} includes two or more letters.
#'
#' @param linecolor The color of the lines when the layout is \code{'facet'}.
#'   The default is \code{NULL}, which sets the line color to \code{'grey50'}.
#'
#' @param linecolor1 The color of the first line when the layout is
#'   \code{'single'}. For example, in \code{opt = 'dv'}, the distance line is
#'   controlled by \code{linecolor1}. The default \code{NULL} sets
#'   \code{linecolor1} to \code{'orange2'}.
#'
#' @param linecolor2 The color of the second line when the layout is
#'   \code{'single'}. For example, in \code{opt = 'dv'}, the velocity line is
#'   controlled by \code{linecolor2}. The default \code{NULL} sets
#'   \code{linecolor2} to \code{'green4'}.
#'
#' @param label.x An optional character string to label the x-axis. If
#'   \code{NULL} (default), the x-axis label will be taken from the predictor
#'   (e.g., age).
#'
#' @param label.y An optional character string to label the y-axis. If
#'   \code{NULL} (default), the y-axis label will be taken from the plot type
#'   (e.g., distance, velocity). When \code{layout = 'facet'}, the label is
#'   removed, and the same label is used as the title.
#'   
#' @param label.title An optional character string to label the title. Default
#'   \code{NULL}.
#'   
#' @param label.subtitle An optional character string to label the title. Default
#'   \code{NULL}
#'
#' @param legendpos A character string to specify the position of the legend. If
#'   \code{NULL} (default), the legend position is set to 'bottom' for distance
#'   and velocity curves in the \code{'single'} layout. For individual-specific
#'   curves, the legend position is set to \code{'none'} to suppress the legend.
#'
#' @param linetype.apv A character string to specify the type of the vertical
#'   line marking the APGV. Default \code{NULL} sets the linetype to
#'   \code{dotted}.
#'
#' @param linewidth.main A numeric value to specify the line width for distance
#'   and velocity curves. The default \code{NULL} sets the width to 0.35.
#'
#' @param linewidth.apv A numeric value to specify the width of the vertical
#'   line marking the APGV. The default \code{NULL} sets the width to 0.25.
#'
#' @param linetype.groupby A character string specifying the line type for
#'   distance and velocity curves when drawing plots for a model with factor
#'   covariates or individual-specific curves. The default is \code{NA}, which
#'   sets the line type to 'solid' and suppresses legends.
#'
#' @param color.groupby A character string specifying the line color for
#'   distance and velocity curves when drawing plots for a model with factor
#'   covariates or individual-specific curves. The default is \code{NA}, which
#'   suppresses legends.
#'
#' @param band.alpha A numeric value to specify the transparency of the CI bands
#'   around the curves. The default \code{NULL} sets the transparency to 0.4.
#'
#' @param show_age_takeoff A logical value (default \code{TRUE}) to indicate
#'   whether to display the ATGV line(s) on the plot.
#'
#' @param show_age_peak A logical value (default \code{TRUE}) to indicate
#'   whether to display the APGV line(s) on the plot.
#'
#' @param show_age_cessation A logical value (default \code{TRUE}) to indicate
#'   whether to display the ACGV line(s) on the plot.
#'
#' @param show_vel_takeoff A logical value (default \code{FALSE}) to indicate
#'   whether to display the TGV line(s) on the plot.
#'
#' @param show_vel_peak A logical value (default \code{FALSE}) to indicate
#'   whether to display the PGV line(s) on the plot.
#'
#' @param show_vel_cessation A logical value (default \code{FALSE}) to indicate
#'   whether to display the CGV line(s) on the plot.
#'
#' @param returndata A logical value (default \code{FALSE}) to indicate whether
#'   to plot the data or return it as a \code{data.frame}.
#'
#' @param returndata_add_parms A logical value (default \code{FALSE}) to specify
#'   whether to add growth parameters to the returned \code{data.frame}. Ignored
#'   when \code{returndata = FALSE}. Growth parameters are added when the
#'   \code{opt} argument includes 'v' or 'V' and \code{apv = TRUE}.
#'
#' @param aux_variables An optional argument to specify variables passed to the
#'   \code{ipts} argument, useful when fitting location-scale or measurement
#'   error models.
#'
#' @inheritParams growthparameters.bgmfit
#' @inherit brms::fitted.brmsfit params
#' @inherit brms::prepare_predictions.brmsfit params
#'
#' @return A plot object (default) or a \code{data.frame} when \code{returndata
#'   = TRUE}.
#'
#' @rdname plot_curves
#' @export
#' 
#' @importFrom rlang .data
#' @importFrom graphics curve
#'
#' @seealso [growthparameters()] [fitted_draws()] [predict_draws()]
#'
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model is fit to 
#' # the 'berkeley_exdata' and saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether the model fit object 'berkeley_exfit' exists
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Population average distance and velocity curves with default options
#' plot_curves(model, opt = 'dv')
#' 
#' # Individual-specific distance and velocity curves with default options
#' # Note that \code{legendpos = 'none'} will suppress the legend positions. 
#' # This suppression is useful when plotting individual-specific curves
#' 
#' plot_curves(model, opt = 'DV')
#' 
#' # Population average distance and velocity curves with APGV
#' 
#' plot_curves(model, opt = 'dv', apv = TRUE)
#' 
#' # Individual-specific distance and velocity curves with APGV
#' 
#' plot_curves(model, opt = 'DV', apv = TRUE)
#' 
#' # Population average distance curve, velocity curve, and APGV with CI bands
#' # To construct CI bands, growth parameters are first calculated for each  
#' # posterior draw and then summarized across draws. Therefore,summary 
#' # option must be set to FALSE
#' 
#' plot_curves(model, opt = 'dv', apv = TRUE, bands = 'dvp', summary = FALSE)
#' 
#' # Adjusted and unadjusted individual curves
#' # Note ipts = NULL (i.e., no interpolation of predictor (i.e., age) to plot a 
#' # smooth curve). This is because it does not a make sense to interploate data 
#' # when estimating adjusted curves. Also, layout = 'facet' (and not default 
#' # layout = 'single') is used for the ease of visualizing the plotted 
#' # adjusted and unadjusted individual curves. However, these lines can be 
#' # superimposed on each other by setting the set layout = 'single'.
#' # For other plots shown above, layout can be set as 'single' or 'facet'
#' 
#' # Separate plots for adjusted and unadjusted curves (layout = 'facet')
#' plot_curves(model, opt = 'au', ipts = NULL, layout = 'facet')
#' 
#' # Superimposed adjusted and unadjusted curves (layout = 'single')
#' plot_curves(model, opt = 'au', ipts = NULL, layout = 'single')
#' 
#' }
#' 
plot_curves.bgmfit <- function(model,
                               opt = 'dv',
                               apv = FALSE,
                               bands = NULL,
                               conf = 0.95,
                               resp = NULL,
                               dpar = NULL,
                               ndraws = NULL,
                               draw_ids = NULL,
                               newdata = NULL,
                               summary = FALSE,
                               digits = 2,
                               re_formula = NULL,
                               numeric_cov_at = NULL,
                               aux_variables = NULL,
                               grid_add = NULL,
                               levels_id = NULL,
                               avg_reffects = NULL,
                               ipts = NULL,
                               model_deriv = TRUE,
                               xrange = NULL,
                               xrange_search = NULL,
                               takeoff = FALSE,
                               trough = FALSE,
                               acgv = FALSE,
                               acgv_velocity = 0.10,
                               seed = 123,
                               estimation_method = 'fitted',
                               allow_new_levels = FALSE,
                               sample_new_levels = "uncertainty",
                               incl_autocor = TRUE,
                               robust = FALSE,
                               transform_draws = NULL,
                               scale = c("response", "linear"),
                               future = FALSE,
                               future_session = 'multisession',
                               cores = NULL,
                               trim = 0,
                               layout = 'single',
                               linecolor = NULL,
                               linecolor1 = NULL,
                               linecolor2 = NULL,
                               label.x = NULL,
                               label.y = NULL,
                               label.title = NULL,
                               label.subtitle = NULL,
                               legendpos = NULL,
                               linetype.apv = NULL,
                               linewidth.main = NULL,
                               linewidth.apv = NULL,
                               linetype.groupby = NA,
                               color.groupby = NA,
                               band.alpha = NULL,
                               show_age_takeoff = TRUE,
                               show_age_peak = TRUE,
                               show_age_cessation = TRUE,
                               show_vel_takeoff = FALSE,
                               show_vel_peak = FALSE,
                               show_vel_cessation = FALSE,
                               returndata = FALSE,
                               returndata_add_parms = FALSE,
                               parms_eval = FALSE,
                               idata_method = NULL,
                               parms_method = 'getPeak',
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
  
  if(is.null(dpar)) {
    dpar <- "mu"
  }
  
  
  model <- getmodel_info(model = model, 
                         dpar = dpar, 
                         resp = resp, 
                         deriv = NULL, 
                         verbose = verbose)
  
  model$model_info[['model_deriv']] <- model_deriv
  model$model_info[['dpar']]        <- dpar
  
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
  
  
  try(insight::check_if_installed(c("jtools", "ggplot2"), stop = FALSE, 
                                  prompt = FALSE))
  
  
  # Initiate non formalArgs()
  yvar <- NULL;
  # IDvar <- NULL;
  groupbytest <- NULL;
  # IDvar <- NULL;
  subindicatorsi <- NULL;
  dy <- NULL;
  yvar <- NULL;
  subindicatorsi <- NULL;
  # IDvar <- NULL;
  Estimate <- NULL;
  groupby <- NULL;
  groupby_line <- NULL;
  groupby_color <- NULL;
  Parameter <- NULL;
  cov_factor_vars <- NULL;
  Estimate.x <- NULL;
  groupby.x <- NULL;
  groupby_line.x <- NULL;
  groupby_color.x <- NULL;
  Estimate.y <- NULL;
  groupby.y <- NULL;
  groupby_line.y <- NULL;
  groupby_color.y <- NULL;
  groupby_fistr <- NULL;
  cov_vars <- NULL;
  ':=' <- NULL;
  . <- NULL;
  uvarby <- NULL;
  
 
  xcall <- match.call()
  match.call.list.in <- as.list(match.call())[-1]
  
  dots <- list(...)
  if ("peak" %in% names(dots)) {
    if (missing(apv)) {
      match.call.list.in[['apv']] <- dots[["peak"]]
    } else if (!missing(apv)) {
      if(is.symbol(match.call.list.in[['apv']])) {
        match.call.list.in[['apv']] <- eval(match.call.list.in[['apv']])
      }
      warning("both 'apv =' and 'peak =' found, ignoring 'peak='")
    }
  }
  
  
  setxcall_ <- match.call()
  post_processing_checks_args <- list()
  post_processing_checks_args[['model']]    <- model
  post_processing_checks_args[['xcall']]    <- setxcall_
  post_processing_checks_args[['resp']]     <- resp
  post_processing_checks_args[['envir']]    <- envir
  # post_processing_checks_args[['deriv']]    <- deriv
  # post_processing_checks_args[['all']]      <- FALSE
  # post_processing_checks_args[['verbose']]  <- verbose
  post_processing_checks_args[['check_d0']] <- FALSE
  post_processing_checks_args[['check_d1']] <- TRUE
  post_processing_checks_args[['check_d2']] <- FALSE
  
  o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
  
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
  model$xcall            <- xcall
  arguments              <- get_args_(match.call.list.in, xcall)
  arguments$model        <- model
  arguments$dpar         <- dpar
  arguments$usesavedfuns <- usesavedfuns
  
  if(is.null(envir)) {
    arguments$envir <- envir <- parent.frame()
  }
  
  if(is.null(ndraws)) {
    arguments$ndraws <- ndraws <- brms::ndraws(model)
  }
  
  if(is.null(model_deriv)) {
    arguments$model_deriv <- model_deriv <- TRUE
  }
  
  if (is.null(idata_method)) {
    arguments$idata_method <- idata_method <- 'm2'
  }
  
  arguments$deriv <- NULL
  
  arguments$ipts <- ipts <- set_for_check_ipts(ipts = ipts, nipts = 50, 
                                               dpar = dpar, verbose = verbose)
  
  
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', 'Est.Error', probtitles)
  
  cores <- 1
  get.cores_ <- get.cores(arguments$cores)
  arguments$cores <- cores <-  get.cores_[['max.cores']] 
  .cores_ps <- get.cores_[['.cores_ps']]
  
  if (future) {
    if(is.null(cores)) stop2c("Please set the number of cores for 'future' by  
                            using the the 'cores' argument, e.g. cores = 4")
    if (future_session == 'multisession') {
      future::plan('multisession', workers = cores)
    } else if (future_session == 'multicore') {
      future::plan('multicore', workers = cores)
    }
  }
  

  if (is.null(newdata)) {
    newdata <- model$model_info$bgmfit.data
  } else {
    newdata <- newdata
  } 
  
  
  if (opt == 'd' | opt == 'D') {
    only_distance_curve <- TRUE
  } else {
    only_distance_curve <- FALSE
  }
  
  if (grepl("v", opt, ignore.case = F) |
      grepl("V", opt, ignore.case = F)) {
    need_velocity_curve <- TRUE
  } else {
    need_velocity_curve <- FALSE
  }
  
  if(only_distance_curve) {
    need_velocity_curve <- FALSE
  }
  
  if(need_velocity_curve) {
    need_xvar_must <- TRUE
  } else {
    need_xvar_must <- FALSE
  }
  
  if(returndata) {
    need_xvar_must <- need_xvar_must
  } else {
    need_xvar_must <- TRUE
  }
  
  arguments$model$model_info[['difx']] <- difx

  if(dpar == "sigma") {
    sigma_model <- get_sigmamodel_info(model = model,
                                       newdata = newdata,
                                       dpar = dpar, 
                                       resp = resp, 
                                       what = 'model',
                                       cov = NULL, 
                                       all = FALSE, 
                                       verbose = verbose)
    
    arguments$model$model_info[['which_sigma_model']] <- 
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
      
      arguments[['transform_draws']] <- transform_draws
    }
    
    if(sigma_model == "basic") {
      if(!is.null(ipts)) {
        stop2c("For sigma_model = ",  
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
      
      
      
      arguments$model$model_info[['xvar_for_sigma_model_basic']] <- xvar
      arguments$newdata <- newdata
    } # if(sigma_model == "basic") {
  } # if(dpar == "sigma") {
  
  assign_function_to_environment(transform_draws, 'transform_draws',
                                 envir = NULL)

  arguments$model$model_info[['transform_draws']] <-
  model$model_info[['transform_draws']] <- transform_draws
  
  get.newdata_args <- list()
  get.newdata_args[['model']]          <- model
  get.newdata_args[['newdata']]        <- newdata
  get.newdata_args[['xvar']]           <- xvar
  get.newdata_args[['idvar']]          <- idvar
  get.newdata_args[['resp']]           <- resp
  get.newdata_args[['numeric_cov_at']] <- numeric_cov_at
  get.newdata_args[['aux_variables']]  <- aux_variables
  get.newdata_args[['levels_id']]      <- levels_id
  get.newdata_args[['xrange']]         <- xrange
  get.newdata_args[['idata_method']]   <- idata_method
  get.newdata_args[['newdata_fixed']]  <- newdata_fixed
  get.newdata_args[['verbose']]        <- verbose
  get.newdata_args[['ipts']]           <- NULL
  
  get.newdata_args$dpar      <- dpar
  newdata.xyadj <- CustomDoCall(get.newdata, get.newdata_args)
  
  get.newdata_args[['ipts']] <- ipts
  newdata       <- CustomDoCall(get.newdata, get.newdata_args)
  
  arguments$newdata <- newdata
  
  list_c <- attr(newdata, 'list_c')
  for (list_ci in names(list_c)) {
    assign(list_ci, list_c[[list_ci]])
  }
  check__ <- c('xvar', 'yvar', 'idvar', 'cov_vars', 'cov_factor_vars', 
               'cov_numeric_vars', 'groupby_fstr', 'groupby_fistr', 
               'uvarby', 'subindicatorsi')
  
  for (check___ in check__) {
    if(!exists(check___)) assign(check___, NULL)
  }
  
  Xx <- xvar
  Yy <- yvar
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  if (is.null(bands)) {
    bands <- ''
  }
  
  if (grepl("d", opt, ignore.case = F) &
      grepl("D", opt, ignore.case = F)) {
    stop2c(
      "Options 'd' and 'D' can not be specified simultanously",
      "\n ",
      " Please check opt argument which is set as '",
      opt,
      "'",
      "\n ",
      " Either specify option 'd' for the population average distance curve",
      "\n ",
      "Or, else option 'D' for the individual specific distance curves",
      "\n ",
      "\n ",
      " Also note that you can combine option 'd' or 'D' with:",
      "\n ",
      "  - option 'v' (population average velocity curve) or",
      "\n ",
      "   option 'V' (individual specific velocity curves)",
      "\n ",
      " - option 'a' (adjusted curves) and 'u' (unadjusted curves)",
      "\n ",
      "\n ",
      "For example, opt = 'dvau', opt = 'dVau', opt = 'Dvau' or opt = 'DVau"
    )
  }
  
  if (grepl("v", opt, ignore.case = F) &
      grepl("V", opt, ignore.case = F)) {
    stop2c(
      "Options 'v' and 'V' can not be specified simultanously",
      "\n ",
      " Please check opt argument which is set as '",
      opt,
      "'",
      "\n ",
      " Either specify option 'v' for the population average velocity curve",
      "\n ",
      "Or, else option 'V' for the individual specific velocity curves",
      "\n ",
      "\n ",
      " Also note that you can combine option 'v' or 'V' with:",
      "\n ",
      "  - option 'd' (population average distance curve) or",
      "\n ",
      "   option 'D' (individual specific distance curves)",
      "\n ",
      " - option 'a' (adjusted curves) and 'u' (unadjusted curves)",
      "\n ",
      "\n ",
      "For example, opt = 'dvau', opt = 'dVau', opt = 'Dvau' or opt = 'DVau"
    )
  }
  
  if (grepl("p", bands, ignore.case = T) & summary) {
    stop2c(
      "To construct bands (e.g., 95%) around the parameter estimates",
      "\n ",
      " (such as APGV, PGV), they are first calculated for each",
      "\n ",
      " posterior draw and then summarised for the give conf limit.",
      "\n ",
      " Therefore,summary option must be set to FALSE"
    )
  }
  
  
  if (grepl("a", bands, ignore.case = T) & summary) {
    stop2c(
      "To construct bands (e.g., 95%) around the adjusted curve estimates, ",
      "\n ",
      " the summary option must be set to FALSE"
    )
  }
  
  
  if (grepl("a", opt, ignore.case = F) |
      grepl("u", opt, ignore.case = F)) {
    ipts <- NULL
    if(verbose) {
      message("The ipts has been set to NULL i.e., ipts = NULL",
              "\n ",
              "This because it does't a make sense to interploate data when",
           "\n ",
           " estimating adjusted/unadjusted curves")
    }
    
    testdata1 <- model$data %>% dplyr::select(dplyr::all_of(idvar)) %>% 
      droplevels() %>% 
      dplyr::mutate(
        groupbytest = interaction(dplyr::across(dplyr::all_of(idvar)))
        ) %>% 
      dplyr::select(dplyr::all_of(groupbytest)) %>% dplyr::ungroup()
    
    testdata2 <- newdata %>% dplyr::select(dplyr::all_of(idvar)) %>% 
      droplevels() %>% 
      dplyr::mutate(groupbytest = 
                      interaction(dplyr::across(dplyr::all_of(idvar)))) %>% 
      dplyr::select(dplyr::all_of(groupbytest)) %>% dplyr::ungroup()
  }
  
  pv <- FALSE
  if (returndata & nchar(opt) > 1) {
    stop2c(
      "For returndata, please specify only one option at a time",
      "\n ",
      " (out of the total six optiona available, i.e., dvDVau)",
      "\n ",
      " For example, opt = 'd'"
    )
  }
  
  if (!grepl("v", opt, ignore.case = F) &
      !grepl("V", opt, ignore.case = F)) {
    apv <- arguments$apv <- FALSE
    pv <- arguments$pv <- FALSE
  }
  
  if(is.null(arguments$pv)) arguments$pv <- FALSE
  
  if(length(list(...)) != 0) arguments <- c(arguments, list(...))
  
  arguments$draw_ids <- draw_ids
  
  for (i in names(arguments)) {
    if(is.symbol(arguments[[i]])) {
      if(deparse(arguments[[i]]) == "") {
        arguments[[i]] <- NULL
      }
    }
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
  

  d. <- CustomDoCall(growthparameters, arguments)
  
  if(is.null(d.)) return(invisible(NULL))
  
  p.                    <- d.[['parameters']]
  probtitles            <- d.[['probtitles']]
  groupby_str_d         <- d.[['groupby_str_d']]
  groupby_str_v         <- d.[['groupby_str_v']]
  
  p.as.d.out_attr       <- p.
  
  d.[['parameters']]    <- NULL
  d.[['probtitles']]    <- NULL
  d.[['groupby_str_d']] <- NULL
  d.[['groupby_str_v']] <- NULL
  
  d. <- d. %>% CustomDoCall(rbind, .) %>% data.frame()
  row.names(d.) <- NULL
  
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
    d. <- prepare_transformations(data = d., model = model,
                                  itransform = itransform_set)
    
    newdata <- prepare_transformations(data = newdata, model = model,
                                  itransform = itransform_set)
  }
 

  curve.d <- 'distance'
  curve.v <- 'velocity'
  
  name.apv <- "APGV"
  name.pv <- "PGV"
  
  name.atv <- "ATGV"
  name.tv <- "TGV"
  
  name.acv <- "ACGV"
  name.cv <- "CGV"
  
  name.vline <- c()
  if(show_age_takeoff)   name.vline <- c(name.vline, name.atv)
  if(show_age_peak)      name.vline <- c(name.vline, name.apv)
  if(show_age_cessation) name.vline <- c(name.vline, name.acv)
  
  name.hline <- c()
  if(show_vel_takeoff)   name.hline <- c(name.hline, name.tv)
  if(show_vel_peak)      name.hline <- c(name.hline, name.pv)
  if(show_vel_cessation) name.hline <- c(name.hline, name.cv)
  
  
  name.hline <- c()

  x_minimum <- min(newdata[[Xx]])
  x_maximum <- max(newdata[[Xx]])
  x_minimum <- floor(x_minimum)
  x_maximum <- floor(x_maximum)
  
  single_plot_pair_color_dv_au <- c('black', 'red')
  
  if (nchar(opt) > 2) {
    layout <- 'facet'
  }
  
  if (is.null(linecolor)) {
    color_single <- "grey50"
  }
  
  if (is.null(linecolor1)) {
    color.d <- "orange2"
      color.adj <- "orange2"
  }
  if (is.null(linecolor2)) {
    color.v <- "green4"
      color.unadj <- "green4"
  }
  
  if (is.null(label.y)) {
    label.d     <- firstup(curve.d)
    label.v     <- firstup(curve.v)
    label.adj   <- firstup('adjusted')
    label.unadj <- firstup('unadjusted')
  } else {
    label.d     <- label.v     <- label.y
    label.adj   <- label.unadj <- label.y
  }
  
  if (is.null(label.x)) {
    label.x     <- paste0(Xx, "")
    # label.x     <- paste0(firstup(Xx), "")
  }
  
  
  if (is.null(legendpos)) {
    if (grepl("D", opt, ignore.case = F) |
        grepl("V", opt, ignore.case = F)) {
      legendpos <- "none"
    } else {
      legendpos <- "bottom"
      legendpos.adj.unadj <- "topleft"
    }
  } else if (!is.null(legendpos)) {
    legendpos <- legendpos.adj.unadj <- legendpos
  }
  
  if (is.null(linetype.apv)) {
    linetype.apv <- 'dotted'
    linetype.pv <- 'dotted'
  } else {
    linetype.apv <- linetype.pv <- linetype.apv
  }
  
  if (is.null(linewidth.main)) {
    linewidth.main <- 0.5
  }
  
  if (is.null(linewidth.apv)) {
    linewidth.apv <- 0.5
    linewidth.pv <- 0.5
  } else {
    linewidth.apv <- linewidth.pv <- linewidth.apv
  }
  
  if (is.null(band.alpha)) {
    band.alpha <- 0.25
  }
  
  if (grepl("d", opt, ignore.case = T) |
      grepl("v", opt, ignore.case = T)) {
    curves <- unique(d.$curve)
    if (length(curves) == 1) {
      layout <- 'facet'
    } else {
      layout <- layout
    }
    
    if (layout == 'facet')
      color.d <- color.v <- color_single
    
    if (grepl("d", opt, ignore.case = T)) {
      d.o <- d.
      index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
      dist.. <- substr(opt, index_opt, index_opt)
      if (grepl("^[[:upper:]]+$", dist..)) {
        d. <-
          d. %>% 
          dplyr::mutate(
            groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
            )
      } else if (!grepl("^[[:upper:]]+$", dist..)) {
        if (is.null(groupby_str_d))
          d. <- d. %>% dplyr::mutate(groupby = NA)
        if (!is.null(groupby_str_d))
          d. <-
            d. %>% 
            dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
              )
      }
      
      if(is.na(d.[['groupby']][1])) {
        d.$groupby_line <- 'solid'
        d.$groupby_color <- 'black'
      } else {
        d.$groupby_line <- d.$groupby
        d.$groupby_color <- d.$groupby
      }
      
      plot.o.d <- d. %>% dplyr::filter(curve == curve.d) %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate,
            group = groupby,
            linetype = groupby_line,
            color = groupby_color
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::scale_x_continuous(breaks = seq(x_minimum, x_maximum, 1)) +
        ggplot2::labs(x = "", y = "", title = label.d) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      plot.o.d <- set_lines_colors(plot.o.d, length(unique(d.[['groupby']])), 
                                   linetype.groupby = linetype.groupby, 
                                   color.groupby = color.groupby)
    
      if (grepl("d", bands, ignore.case = T)) {
        plot.o.d <- plot.o.d +
          ggplot2::geom_ribbon(
            data = d. %>% dplyr::filter(curve == curve.d),
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '')]],
              ymax = .data[[paste0(probtitles[2], '')]],
              group = groupby,
              linetype = groupby_line,
              color = groupby_color,
              fill = groupby_color
            ),
            alpha = band.alpha
          )
        plot.o.d <- set_lines_colors_ribbon(plot.o.d, guideby = 'color')
      }
      
      d. <- d.o
      if ('curve' %in% names(d.)) {
        d.out <- d. %>% dplyr::select(-dplyr::all_of('curve')) # curve to 'curve'
      } else {
        d.out <- d.
      }
    }
    if (!grepl("d", opt, ignore.case = T)) {
      plot.o.d <- NULL
    }
    
    if (grepl("v", opt, ignore.case = T)) {
      d.o <- d.
      index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
      velc.. <- substr(opt, index_opt, index_opt)
      if (grepl("^[[:upper:]]+$", velc..)) {
        d. <-
          d. %>% 
          dplyr::mutate(
            groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
            )
      } else if (!grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_v))
          d. <- d. %>% dplyr::mutate(groupby = NA)
        if (!is.null(groupby_str_v))
          d. <-
            d. %>% dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
              )
      }
      
      if(is.na(d.[['groupby']][1])) {
        d.$groupby_line <- 'solid'
        d.$groupby_color <- 'black'
      } else {
        d.$groupby_line <- d.$groupby
        d.$groupby_color <- d.$groupby
      }
      
      plot.o.v <- d. %>% dplyr::filter(curve == curve.v) %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate,
            group = groupby,
            linetype = groupby_line,
            color = groupby_color
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::scale_x_continuous(breaks = seq(x_minimum, x_maximum, 1)) +
        # ggplot2::labs(title = label.v) +
        # ggplot2::xlab("") +
        # ggplot2::ylab("") +
        ggplot2::labs(x = "", y = "", title = label.v) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      
      plot.o.v <- set_lines_colors(plot.o.v, length(unique(d.[['groupby']])), 
                                   linetype.groupby = linetype.groupby, 
                                   color.groupby = color.groupby)
      
      if (grepl("v", bands, ignore.case = T)) {
        plot.o.v <- plot.o.v +
          ggplot2::geom_ribbon(
            data = d. %>% dplyr::filter(curve == curve.v),
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '')]],
              ymax = .data[[paste0(probtitles[2], '')]],
              group = groupby,
              linetype = groupby_line,
              color = groupby_color,
              fill = groupby_color
            ),
            alpha = band.alpha
          )
        plot.o.v <- set_lines_colors_ribbon(plot.o.v, guideby = 'color')
      }
      
      if (pv) {
        data_hline <- p. %>% dplyr::filter(Parameter == name.pv)
        plot.o.v <- plot.o.v +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = .data[['Estimate']]),
            linewidth = linewidth.pv,
            linetype = linetype.pv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o.v <- plot.o.v +
            ggplot2::annotate(
              "rect",
              xmin = -Inf,
              xmax = Inf,
              ymin = (data_hline[[paste0(probtitles[1], '')]]),
              ymax = (data_hline[[paste0(probtitles[2], '')]]),
              alpha = band.alpha
            )
        }
      }
      
      if (!is.null(name.vline) & !is.null(p.)  ) {
        data_vline <- 
          dplyr::filter(p., grepl(paste(name.vline, collapse  = "|"),
                                              Parameter))
        plot.o.v <- plot.o.v +
          ggplot2::geom_vline(
            data = data_vline,
            ggplot2::aes(xintercept = .data[['Estimate']]),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        if (grepl("p", bands, ignore.case = T)) {
          plot.o.v <- plot.o.v +
            ggplot2::annotate(
              "rect",
              xmin = (data_vline[[paste0(probtitles[1], '')]]),
              xmax = (data_vline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      } #  if (!is.null(name.vline) & !is.null(p.) ) {
      
      if (!is.null(name.hline) & !is.null(p.) ) {
        data_hline <- 
          dplyr::filter(p., grepl(paste(name.hline, collapse  = "|"), 
                                              Parameter))
        plot.o.v <- plot.o.v +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = t.s.axis$fwd(.data[['Estimate']]) ),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        if (grepl("p", bands, ignore.case = T)) {
          plot.o.v <- plot.o.v +
            ggplot2::annotate(
              "rect",
              xmin = (data_vline[[paste0(probtitles[1], '')]]),
              xmax = (data_vline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      } #  if (!is.null(name.hline)) {
      
      d. <- d.o
      if ('curve' %in% names(d.)) {
        d.out <- d. %>% dplyr::select(-dplyr::all_of('curve'))
      } else {
        d.out <- d.
      }
    }
    
    if (!grepl("v", opt, ignore.case = T)) {
      plot.o.v <- NULL
    }
    
    if (length(curves) > 1 & layout == 'facet') {
      plot.o <- patchwork::wrap_plots(plot.o.d + plot.o.v) %>%
        add_global_label(
          Xlab = label.x,
          Ylab = "",
          size = 5,
          Xgap = 0.08,
          Ygap = 0.04
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = legendpos,
                       legend.direction = 'horizontal')
      
    } else if (length(curves) == 1 & layout == 'facet') {
      if (!is.null(plot.o.d)) {
        plot.o <- plot.o.d +
          ggplot2::labs(x = label.x, y = label.d) +
          ggplot2::theme(plot.title = ggplot2::element_blank()) +
          ggplot2::theme(legend.position = legendpos,
                         legend.direction = 'horizontal')
      } else if (!is.null(plot.o.v)) {
        plot.o <- plot.o.v +
          ggplot2::labs(x = label.x, y = label.v) +
          ggplot2::theme(plot.title = ggplot2::element_blank()) +
          ggplot2::theme(legend.position = legendpos,
                         legend.direction = 'horizontal')
      }
    }
    
    
    if (length(curves) > 1 & layout == 'single') {
      data_d <- subset(d., curve == "distance")
      data_v <- subset(d., curve == "velocity")
      by_join_ <- c(idvar, Xx, groupby_str_d)
      by_join_ <- unique(by_join_)
      data_dv <- dplyr::left_join(data_d, data_v, by = by_join_)
      data_dv.o <- data_dv
      if (grepl("d", opt, ignore.case = T)) {
        index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
        dist.. <- substr(opt, index_opt, index_opt)
        if (grepl("^[[:upper:]]+$", dist..)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
              ) %>%
            dplyr::mutate(
              groupby.x = 
                interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
              )
        } else if (!grepl("^[[:upper:]]+$", dist..)) {
          if (is.null(groupby_str_d)) {
            data_dv <- data_dv %>% dplyr::mutate(groupby = NA) %>%
              dplyr::mutate(groupby.x = NA)
          } else if (!is.null(groupby_str_d)) {
            data_dv <- data_dv %>%
              dplyr::mutate(
                groupby = 
                  interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
                ) %>%
              dplyr::mutate(
                groupby.x = 
                  interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
                )
          }
        }
      }
      
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
        if (grepl("^[[:upper:]]+$", velc..)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
              ) %>%
            dplyr::mutate(groupby.y = groupby)
        } else if (!grepl("^[[:upper:]]+$", velc..)) {
          if (is.null(groupby_str_v)) {
            data_dv <- data_dv %>% dplyr::mutate(groupby = NA) %>%
              dplyr::mutate(groupby.y = NA)
          } else if (!is.null(groupby_str_v)) {
            data_dv <- data_dv %>%
              dplyr::mutate(
                groupby = 
                  interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
                ) %>%
              dplyr::mutate(
                groupby.y = 
                  interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
                )
          }
        }
      }
      
      if (grepl("^[[:upper:]]+$", dist..) &
          !grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_v)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby.x = 
                interaction(dplyr::across(dplyr::all_of(groupby_str_d))),
              groupby.y = NA)
        } else if (!is.null(groupby_str_v)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby.x = 
                interaction(dplyr::across(dplyr::all_of(groupby_str_d))),
              groupby.y =
                interaction(dplyr::across(dplyr::all_of(groupby_str_v))))
        }
      }
      
      if (!grepl("^[[:upper:]]+$", dist..) &
          grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_d)) {
          data_dv <-
            data_dv %>% dplyr::mutate(
              groupby.x = NA,
              groupby.y =
                interaction(dplyr::across(dplyr::all_of(groupby_str_v))))
        } else if (!is.null(groupby_str_d)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby.x = 
                interaction(dplyr::across(dplyr::all_of(groupby_str_d))),
              groupby.y =
                interaction(dplyr::across(dplyr::all_of(groupby_str_v))))
        }
      }
      
      t.s.axis <- with(data_dv, transform.sec.axis(Estimate.x, Estimate.y))
      
      if(is.na(uvarby)) {
        if(is.na(data_dv[['groupby.x']][1])) {
          legendlabs_mult_singel <- c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          data_dv$groupby_line.x <- legendlabs_mult_singel[1]
          data_dv$groupby_color.x <- legendlabs_mult_singel[1]
          data_dv$groupby_line.y <- legendlabs_mult_singel[2]
          data_dv$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          data_dv$groupby_line.x <- data_dv$groupby.x
          data_dv$groupby_color.x <- data_dv$groupby.x
          data_dv$groupby_line.y <- data_dv$groupby.y
          data_dv$groupby_color.y <- data_dv$groupby.y
          legendlabs_mult_mult <- unique(data_dv[['groupby.x']])
        }
      }
      
      
      
      if(!is.na(uvarby)) {
        if(is.null(cov_factor_vars)) {
          legendlabs_mult_singel <- c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          data_dv$groupby_line.x <- legendlabs_mult_singel[1]
          data_dv$groupby_color.x <- legendlabs_mult_singel[1]
          data_dv$groupby_line.y <- legendlabs_mult_singel[2]
          data_dv$groupby_color.y <- legendlabs_mult_singel[2]
          legendlabs_mult_mult <- NULL # unique(data_dv[['groupby.x']])
        } else {
          data_dv$groupby_line.x <- data_dv$groupby.x
          data_dv$groupby_color.x <- data_dv$groupby.x
          data_dv$groupby_line.y <- data_dv$groupby.y
          data_dv$groupby_color.y <- data_dv$groupby.y
          legendlabs_mult_mult <- unique(data_dv[['groupby.x']])
        }
      }
      
      if(length( unique(round(data_dv$Estimate.y, 6))) == 1) {
        stop2c("The velocity estimates are identical over the entire range of x",
             "\n  ", 
             "Can't draw distance and velocity curves together")
      }
      
      plot.o <- data_dv %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate.x,
            group = groupby.x,
            linetype = groupby_line.x,
            colour = groupby_color.x
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = t.s.axis$fwd(Estimate.y),
            group = groupby.y,
            linetype = groupby_line.y,
            colour = groupby_color.y
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::scale_y_continuous(sec.axis =
                                      ggplot2::sec_axis(~ t.s.axis$rev(.),
                                                        name = label.v)) +
        ggplot2::labs(x = label.x, y = label.d, color = "") +
        ggplot2::scale_x_continuous(breaks = seq(x_minimum, x_maximum, 1)) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90))
      
      
      getbuiltingg <- ggplot2::ggplot_build(plot.o)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- 
        rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- 
        rep(get_color_, ngrpanels)
      
      if(!exists('legendlabs_mult_line')) legendlabs_mult_line <- 'solid'
      if(!exists('legendlabs_mult_color')) legendlabs_mult_color <- 'black'
      if(!exists('legendlabs_mult_singel')) legendlabs_mult_singel <- 'solid'
      
      
      
      
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- legendlabs_mult_line
        get_color_ <- legendlabs_mult_color
        legendlabs_ <- legendlabs_mult_singel
      }
      
      plot.o <- plot.o + 
        ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
        ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
      
      
      if (grepl("d", bands, ignore.case = T)) {
        plot.o <- plot.o +
          ggplot2::geom_ribbon(
            data = data_dv,
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '.x')]],
              ymax = .data[[paste0(probtitles[2], '.x')]],
              group = groupby.x,
              linetype = groupby_line.x,
              # color = groupby_color.x,
              fill = groupby_color.x,
            ),
            alpha = band.alpha
          )
      }
      
      if (grepl("v", bands, ignore.case = T)) {
        plot.o <- plot.o +
          ggplot2::geom_ribbon(
            data = data_dv,
            ggplot2::aes(
              ymin = t.s.axis$fwd(.data[[paste0(probtitles[1], '.y')]]),
              ymax = t.s.axis$fwd(.data[[paste0(probtitles[2], '.y')]]),
              group = groupby.y,
              linetype = groupby_line.y,
              # color = groupby_color.y,
              fill = groupby_color.y,
            ),
            alpha = band.alpha
          )
      }
      
      
      if((grepl("d", bands, ignore.case = T) & 
          !grepl("v", bands, ignore.case = T)) |
         !grepl("d", bands, ignore.case = T) & 
         grepl("v", bands, ignore.case = T)
      ) {
        one_band <- TRUE
      } else {
        one_band <- FALSE
      }
      
      if(one_band & ngrpanels == 1) {
        if(grepl("d", bands, ignore.case = T) & 
           !grepl("v", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::scale_fill_manual(values=legendlabs_mult_color[1], 
                                       guide = 'none')
        }
        if(!grepl("d", bands, ignore.case = T) & 
           grepl("v", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::scale_fill_manual(values=legendlabs_mult_color[2], 
                                       guide = 'none')
        }
      }
      
      if(!one_band & ngrpanels == 1) {
        plot.o <- plot.o +
          ggplot2::scale_fill_manual(values=get_color_, guide = 'none')
      }
      
      if (pv) {
        data_hline <- p. %>% dplyr::filter(Parameter == name.pv)
        plot.o <- plot.o +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = t.s.axis$fwd(.data[['Estimate']])),
            linewidth = linewidth.pv,
            linetype = linetype.pv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::annotate(
              "rect",
              xmin = -Inf,
              xmax = Inf,
              ymin = t.s.axis$fwd(data_hline[[paste0(probtitles[1], '')]]),
              ymax = t.s.axis$fwd(data_hline[[paste0(probtitles[2], '')]]),
              alpha = band.alpha
            )
        }
      }
      
      if (!is.null(name.vline) & !is.null(p.)) {
        data_vline <- dplyr::filter(p., grepl(paste(name.vline, 
                                                    collapse  = "|"), 
                                              Parameter))
        plot.o <- plot.o +
          ggplot2::geom_vline(
            data = data_vline,
            ggplot2::aes(xintercept = .data[['Estimate']]),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::annotate(
              "rect",
              xmin = (data_vline[[paste0(probtitles[1], '')]]),
              xmax = (data_vline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      } # if (!is.null(name.vline) & !is.null(p.) ) {
      
      if (!is.null(name.hline) & !is.null(p.) ) {
        data_hline <- dplyr::filter(p., grepl(paste(name.hline, 
                                                    collapse  = "|"), 
                                              Parameter))
        plot.o <- plot.o +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = t.s.axis$fwd(.data[['Estimate']]) ),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::annotate(
              "rect",
              xmin = (data_hline[[paste0(probtitles[1], '')]]),
              xmax = (data_hline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      } # if (!is.null(name.hline)) {
      
      data_dv <- data_dv.o
      if ('curve' %in% names(data_dv)) {
        d.out <- data_dv %>% dplyr::select(-dplyr::all_of('curve'))
      } else {
        d.out <- data_dv
      }
    }
  }
  
  groupby_str_au <- groupby_fistr
  
  if (grepl("a", opt, ignore.case = T) |
      grepl("u", opt, ignore.case = T)) {
    if (!is.null(cov_vars)) {
     # stop("Adjusted curves not yet supported for model with covariate(s)")
    }
    
    if (grepl("a", opt, ignore.case = T)) {
      xyadj_ed <- xyadj_curves(model, 
                               x = NULL,
                               y = NULL,
                               id = NULL,
                               v = NULL,
                               newdata = newdata.xyadj, 
                               ndraws = ndraws,
                               draw_ids = draw_ids,
                               resp = resp, 
                               tomean = TRUE,
                               conf = conf, 
                               robust = robust,
                               summary = summary, 
                               numeric_cov_at = numeric_cov_at,
                               aux_variables = aux_variables,
                               levels_id = levels_id,
                               ipts = ipts,
                               xrange = xrange, 
                               idata_method = idata_method,
                               verbose = verbose,
                               model_deriv = NULL,
                               deriv = NULL, 
                               envir = envir,
                               ...) 
 
      out_a_ <- trimlines_curves(model, 
                                 x = Xx,
                                 y = Yy,
                                 id = idvar,
                                 newdata = xyadj_ed, 
                                 ndraws = ndraws,
                                 draw_ids = draw_ids,
                                 resp = resp, 
                                 level = 0,
                                 trim = trim, 
                                 estimation_method = estimation_method,
                                 verbose = verbose,
                                 model_deriv = NULL,
                                 deriv = NULL, 
                                 envir = envir,
                                 ...)
      
      if(any(itransform_set != "")) {
        out_a_ <- prepare_transformations(data = out_a_, model = model,
                                          itransform = itransform_set)
      }
      
      d.out <- out_a_

      dots <- list(...)
      set_get_dv <- FALSE
      if(!is.null(dots$get_dv)) {
        if(dots$get_dv) {
          if(verbose) message("executing 'get_dv'!")
          set_get_dv <- TRUE
        }
      }
      
      if(set_get_dv) {
        return(out_a_)
      }
      
      if(!is.null(dots$xadj_tmt)) {
        if(dots$xadj_tmt) {
          return(out_a_)
        }
      }
      
      if(!is.null(dots$xadj_tmf)) {
        if(dots$xadj_tmf) {
          return(out_a_)
        }
      }

      out_a_ <-
        out_a_ %>%
        dplyr::mutate(
          groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_au)))
          )
      
      # x_minimum_a_ <- floor(min(out_a_[[Xx]]))
      # x_maximum_a_ <- ceiling(max(out_a_[[Xx]]))
      x_minimum_a_ <- x_minimum
      x_maximum_a_ <- x_maximum
      
      out_a_ <- out_a_[out_a_[[Xx]] >= x_minimum_a_ & 
                         out_a_[[Xx]] <= x_maximum_a_, ]

      out_a_ <- out_a_ %>% dplyr::mutate(groupby.x = groupby, 
                                         groupby.y = groupby.x)
      
      if(is.na(uvarby)) {
        if(is.na(out_a_[['groupby']][1])) {
          legendlabs_mult_singel <- c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_a_$groupby_line.x <- legendlabs_mult_singel[1]
          out_a_$groupby_color.x <- legendlabs_mult_singel[1]
          out_a_$groupby_line.y <- legendlabs_mult_singel[2]
          out_a_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_a_$groupby_line.x <- out_a_$groupby.x
          out_a_$groupby_color.x <- out_a_$groupby.x
          out_a_$groupby_line.y <- out_a_$groupby.y
          out_a_$groupby_color.y <- out_a_$groupby.y
          legendlabs_mult_mult <- unique(out_a_[['groupby']])
        }
      }
      
      
      if(!is.na(uvarby)) {
        if(is.null(cov_factor_vars)) {
          legendlabs_mult_singel <- c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_a_$groupby_line.x <- legendlabs_mult_singel[1]
          out_a_$groupby_color.x <- legendlabs_mult_singel[1]
          out_a_$groupby_line.y <- legendlabs_mult_singel[2]
          out_a_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_a_$groupby_line.x <- out_a_$groupby.x
          out_a_$groupby_color.x <- out_a_$groupby.x
          out_a_$groupby_line.y <- out_a_$groupby.y
          out_a_$groupby_color.y <- out_a_$groupby.y
          legendlabs_mult_mult <- unique(out_a_[['groupby']])
        }
      }
      
      plot.o.a <- out_a_ %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = !!as.name(Yy),
            group = groupby.x # ,
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::labs(x = label.x, y = label.d, color = "") +
        ggplot2::scale_x_continuous(breaks =
                                      seq(x_minimum_a_, x_maximum_a_, 1)) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(y = paste0("Adjusted ", "individual curves")) +
        ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90)) +
        ggplot2::labs(title = label.adj) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      
      getbuiltingg <- ggplot2::ggplot_build(plot.o.a)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- 
        rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- 
        rep(get_color_, ngrpanels)
      
      # These will be carried forward for ribbon also (below)
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- 'solid' # legendlabs_mult_line
        get_color_ <- 'black' # legendlabs_mult_color
        legendlabs_ <- NULL # legendlabs_mult_singel
      }
      
      plot.o.a <- plot.o.a +
        ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
        ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
      
      if (grepl("a", bands, ignore.case = T)) {
        plot.o.a <- plot.o.a +
          ggplot2::geom_ribbon(
            data = out_a_,
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '')]],
              ymax = .data[[paste0(probtitles[2], '')]],
              group = groupby.x,
              linetype = groupby_line.x,
              # color = groupby_color.x,
              fill = groupby_color.x,
            ),
            alpha = band.alpha
          )
      }
      
      
      if (nchar(opt) == 1) {
        plot.o.a <- plot.o.a +
          ggplot2::theme(plot.title = ggplot2::element_blank())
      } else if (nchar(opt) > 1) {
        plot.o.a <- plot.o.a +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
      if (!grepl("u", opt, ignore.case = T)) {
        suppressMessages({
        })
      }
    } else if (!grepl("a", opt, ignore.case = T)) {
      plot.o.a <- NULL
    }
    
    
    if (grepl("u", opt, ignore.case = T)) {
      xyunadj_ed <- xyunadj_curves(model, 
                                   x = NULL,
                                   y = NULL,
                                   id = NULL,
                                   newdata = NULL,
                                   ndraws = ndraws,
                                   draw_ids = draw_ids,
                                   resp = resp, 
                                   verbose = verbose,
                                   model_deriv = NULL,
                                   deriv = NULL, 
                                   envir = envir,
                                   ...)
      
      out_u_ <- trimlines_curves(model, 
                                 x = Xx,
                                 y = Yy,
                                 id = idvar,
                                 newdata = xyunadj_ed, 
                                 ndraws = ndraws,
                                 draw_ids = draw_ids,
                                 resp = resp, 
                                 level = 0,
                                 trim = trim, 
                                 estimation_method = estimation_method,
                                 verbose = verbose,
                                 model_deriv = NULL,
                                 deriv = NULL, 
                                 envir = envir,
                                 ...)
      
      if(any(itransform_set != "")) {
        out_u_ <- prepare_transformations(data = out_u_, model = model,
                                          itransform = itransform_set)
      }
      d.out <- out_u_
      
      out_u_ <-
        out_u_ %>%
        dplyr::mutate(
          groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_au)))
          )
      
      out_u_ <- out_u_ %>% dplyr::mutate(groupby.x = groupby, 
                                         groupby.y = groupby.x)
    
      if(is.na(uvarby)) { 
        if(is.na(out_u_[['groupby']][1])) {
          legendlabs_mult_singel <- c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_u_$groupby_line.x <- legendlabs_mult_singel[1]
          out_u_$groupby_color.x <- legendlabs_mult_singel[1]
          out_u_$groupby_line.y <- legendlabs_mult_singel[2]
          out_u_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_u_$groupby_line.x <- out_u_$groupby.x
          out_u_$groupby_color.x <- out_u_$groupby.x
          out_u_$groupby_line.y <- out_u_$groupby.y
          out_u_$groupby_color.y <- out_u_$groupby.y
          legendlabs_mult_mult <- unique(out_u_[['groupby']])
        }
      }
      
      if(!is.na(uvarby)) { 
        if(is.null(cov_factor_vars)) {
          legendlabs_mult_singel <- c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_u_$groupby_line.x <- legendlabs_mult_singel[1]
          out_u_$groupby_color.x <- legendlabs_mult_singel[1]
          out_u_$groupby_line.y <- legendlabs_mult_singel[2]
          out_u_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_u_$groupby_line.x <- out_u_$groupby.x
          out_u_$groupby_color.x <- out_u_$groupby.x
          out_u_$groupby_line.y <- out_u_$groupby.y
          out_u_$groupby_color.y <- out_u_$groupby.y
          legendlabs_mult_mult <- unique(out_u_[['groupby']])
        }
      }
      
    plot.o.u <- out_u_ %>%
      ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
      ggplot2::geom_line(
        ggplot2::aes(
          y = !!as.name(Yy),
          group = groupby.y
        ),
        linewidth = linewidth.main
      ) +
      ggplot2::labs(x = label.x, y = label.d, color = "") +
      ggplot2::scale_x_continuous(breaks = seq(x_minimum, x_maximum, 1)) +
      jtools::theme_apa(legend.pos = legendpos) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(y = paste0("Unadjusted ", "individual curves")) +
      ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90)) +
      ggplot2::labs(title = label.unadj) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      getbuiltingg <- ggplot2::ggplot_build(plot.o.u)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- 
        rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- 
        rep(get_color_, ngrpanels)
      
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- 'solid' # legendlabs_mult_line
        get_color_ <- 'black' # legendlabs_mult_color
        legendlabs_ <- NULL # legendlabs_mult_singel
      }
      
      plot.o.u <- plot.o.u +
        ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
        ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
      
      
      if (!grepl("a", opt, ignore.case = T)) {
        suppressMessages({
        })
      }
      if (nchar(opt) == 1) {
        plot.o.u <- plot.o.u +
          ggplot2::theme(plot.title = ggplot2::element_blank())
      } else if (nchar(opt) > 1) {
        plot.o.u <- plot.o.u +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
    } else if (!grepl("u", opt, ignore.case = T)) {
      plot.o.u <- NULL
    }
    
    if (grepl("a", opt, ignore.case = T) &
        !grepl("u", opt, ignore.case = T)) {
      plot.o <- plot.o.a
    } else if (!grepl("a", opt, ignore.case = T) &
               grepl("u", opt, ignore.case = T)) {
      plot.o <- plot.o.u
    } else if (grepl("a", opt, ignore.case = T) &
               grepl("u", opt, ignore.case = T)) {
      
      if (layout == 'facet') {
        out_a_u_ <-
          d.out <- out_a_ %>% dplyr::mutate(curve = 'Adjusted') %>%
          dplyr::bind_rows(., out_u_ %>%
                             dplyr::mutate(curve = 'Unadjusted')) %>%
          data.frame()
        
        # x_minimum_a_ <- floor(min(out_a_[[Xx]]))
        # x_maximum_a_ <- ceiling(max(out_a_[[Xx]]))
        x_minimum_a_ <- x_minimum
        x_maximum_a_ <- x_maximum
        
        out_a_u_ <- out_a_u_[out_a_u_[[Xx]] >= x_minimum_a_ & 
                               out_a_u_[[Xx]] <= x_maximum_a_, ]
        
        out_a_u_ <- out_a_u_ %>% dplyr::mutate(groupby.x = groupby, 
                                               groupby.y = groupby.x)
        
        
        if(is.na(uvarby)) {
          if(is.na(out_a_u_[['groupby']][1])) {
            legendlabs_mult_singel <- c('Distance', 'Velocity')
            legendlabs_mult_color <- single_plot_pair_color_dv_au
            legendlabs_mult_line <- c('solid', 'solid')
            out_a_u_$groupby_line.x <- legendlabs_mult_singel[1]
            out_a_u_$groupby_color.x <- legendlabs_mult_singel[1]
            out_a_u_$groupby_line.y <- legendlabs_mult_singel[2]
            out_a_u_$groupby_color.y <- legendlabs_mult_singel[2]
          } else {
            out_a_u_$groupby_line.x <- out_a_u_$groupby.x
            out_a_u_$groupby_color.x <- out_a_u_$groupby.x
            out_a_u_$groupby_line.y <- out_a_u_$groupby.y
            out_a_u_$groupby_color.y <- out_a_u_$groupby.y
            legendlabs_mult_mult <- unique(out_a_u_[['groupby']])
          }
        }
        
        if(!is.na(uvarby)) {
          if(is.null(cov_factor_vars)) {
            legendlabs_mult_singel <- c('Distance', 'Velocity')
            legendlabs_mult_color <- single_plot_pair_color_dv_au
            legendlabs_mult_line <- c('solid', 'solid')
            out_a_u_$groupby_line.x <- legendlabs_mult_singel[1]
            out_a_u_$groupby_color.x <- legendlabs_mult_singel[1]
            out_a_u_$groupby_line.y <- legendlabs_mult_singel[2]
            out_a_u_$groupby_color.y <- legendlabs_mult_singel[2]
          } else {
            out_a_u_$groupby_line.x <- out_a_u_$groupby.x
            out_a_u_$groupby_color.x <- out_a_u_$groupby.x
            out_a_u_$groupby_line.y <- out_a_u_$groupby.y
            out_a_u_$groupby_color.y <- out_a_u_$groupby.y
            legendlabs_mult_mult <- unique(out_a_u_[['groupby']])
          }
        }
        
        plot.o <- out_a_u_ %>%
          ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
          ggplot2::geom_line(
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby.x,
              linetype = groupby_line.x,
              colour = groupby_color.x
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::labs(x = label.x,
                        y = label.d,
                        color = "") +
          ggplot2::scale_x_continuous(breaks =
                                        seq(x_minimum_a_, x_maximum_a_, 1)) +
          jtools::theme_apa(legend.pos = legendpos) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::theme(axis.title.y.right =
                           ggplot2::element_text(angle = 90)) +
          ggplot2::theme(legend.position = "none") +
          ggplot2::labs(y = paste0("Individual curves")) +
          ggplot2::facet_wrap(~ curve, scales = 'free_x')
        
        
        getbuiltingg <- ggplot2::ggplot_build(plot.o)
        get_line_  <- getbuiltingg$data[[1]]["linetype"]
        get_color_ <- getbuiltingg$data[[1]]["colour"]
        get_fill_  <- getbuiltingg$data[[1]]["colour"]
        ngrpanels  <- getbuiltingg$data[[1]]["group"]
        ngrpanels <- length(unique(unlist(ngrpanels)))
        get_line_ <- unique(unlist(get_line_))
        get_color_ <- unique(unlist(get_color_))
        if(length(get_line_) != ngrpanels) get_line_ <- 
          rep(get_line_, ngrpanels)
        if(length(get_color_) != ngrpanels) get_color_ <- 
          rep(get_color_, ngrpanels)
        
        if(ngrpanels > 1) {
          get_line_ <- get_line_
          get_color_ <- get_color_
          legendlabs_ <- legendlabs_mult_mult
        } else if(ngrpanels == 1) {
          get_line_ <- 'solid'  # legendlabs_mult_line
          get_color_ <- 'black' # legendlabs_mult_color
          legendlabs_ <- NULL   # legendlabs_mult_singel
        }
        
        plot.o <- plot.o +
          ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
          ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
        
      }
      
      if (layout == 'single') {
        plot.o <- out_a_ %>%
          ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
          ggplot2::geom_line(
            data = out_u_,
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby,
              # linetype = linetype.groupby, # This does not plot lines TODO
              colour = label.unadj
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::geom_line(
            data = out_a_,
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby,
              # linetype = linetype.groupby, # This does not plot lines TODO
              colour = label.adj
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::labs(x = label.x,
                        y = label.d,
                        color = "") +
          ggplot2::scale_color_manual(values = single_plot_pair_color_dv_au) +
          ggplot2::scale_x_continuous(breaks =
                                        seq(x_minimum_a_, x_maximum_a_, 1)) +
          jtools::theme_apa(legend.pos = legendpos.adj.unadj) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::labs(y = paste0("Individual curves")) +
          ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90))
        
      }
    }
  }
  
  if (nchar(opt) > 2) {
    if (!exists('plot.o.d'))
      plot.o.d <- NULL
    if (!exists('plot.o.v'))
      plot.o.v <- NULL
    if (!exists('plot.o.a'))
      plot.o.a <- NULL
    if (!exists('plot.o.u'))
      plot.o.u <- NULL
    
    suppressMessages({
      if (!is.null(plot.o.d)) {
        plot.o.d <- plot.o.d +
          ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
      if (!is.null(plot.o.v)) {
        plot.o.v <- plot.o.v +
          ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
      if (!is.null(plot.o.a)) {
        plot.o.a <- plot.o.a +
          ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
      if (!is.null(plot.o.u)) {
        plot.o.u <- plot.o.u +
          ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
    })
    
    plot.list <- list(plot.o.d, plot.o.v, plot.o.a, plot.o.u)
    plot.list <- plot.list[lengths(plot.list) != 0]
    
    plot.o <- patchwork::wrap_plots(plot.list,
                                    ncol = 2, nrow = NULL) %>%
      add_global_label(
        Xlab = label.x,
        Ylab = "",
        size = 5,
        Xgap = 0.08,
        Ygap = 0.04)
    plot.o <- plot.o +  patchwork::plot_layout(guides = "collect")
  }
 
  if (!returndata) {
    print(plot.o)
    if (grepl("d", opt, ignore.case = F) |
        grepl("v", opt, ignore.case = F)) {
      if(apv | takeoff | trough | acgv)  print(p.)
    }
    if(!is.null(p.as.d.out_attr)) {
      plot.o[['growthparameters']] <- p.as.d.out_attr
    }
    return(plot.o)
  } else if (returndata) {
    attr(d.out, 'growthparameters') <- p.as.d.out_attr
    if(returndata_add_parms) {
      if(!is.null(p.as.d.out_attr)) {
        d.out <- add_parms_to_curve_data(d.out, 
                                         gpdata = NULL,
                                         Parametername = "Parameter",
                                         parmcols = set_names_,
                                         nonparmcols = groupby_str_v,
                                         byjoincols = groupby_str_v)
      } # if(!is.null(p.as.d.out_attr)) {
    } # if(returndata_add_parms) {
    return(d.out)
  } # else if (returndata) {
  
} # end plot_curves




#' @rdname plot_curves
#' @export
plot_curves <- function(model, ...) {
  UseMethod("plot_curves")
}


