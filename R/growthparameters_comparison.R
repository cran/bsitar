

#' Compare growth parameters
#' 
#'@description The \strong{growthparameters_comparison()} function estimates and
#'  compare growth parameters such as peak growth velocity and the age at peak
#'  growth velocity. This function is a wrapper around the
#'  [marginaleffects::comparisons()] and [marginaleffects::avg_comparisons()].
#'  The [marginaleffects::comparisons()] computes unit-level (conditional)
#'  estimates whereas [marginaleffects::avg_comparisons()] return average
#'  (marginal) estimates. A detailed explanation is available
#'  [here](https://marginaleffects.com). Note that for the current use case,
#'  i.e., to estimate and compare growth parameters, the arguments
#'  \code{variables} and \code{comparion} of [marginaleffects::comparisons()]
#'  and [marginaleffects::avg_comparisons()] are modified (see below).
#'  Furthermore, comparison of growth parameters is performed via the
#'  \code{hypothesis} argument of the [marginaleffects::comparisons()] and
#'  [marginaleffects::avg_comparisons()] functions.
#'
#'
#' @details The \code{growthparameters_comparison} function estimates and
#'   returns the following growth parameters:
#' \itemize{
#'   \item pgv  - peak growth velocity
#'   \item apgv - age at peak growth velocity
#'   \item tgv  - takeoff growth velocity
#'   \item atgv - age at takeoff growth velocity
#'   \item cgv  - cessation growth velocity
#'   \item acgv - age at cessation growth velocity
#' }
#' 
#' The takeoff growth velocity is the lowest velocity just before the peak
#' starts and it indicates the beginning of the pubertal growth spurt. The
#' cessation growth velocity indicates the end of the active pubertal growth
#' spurt and is calculated as some percentage of the peak velocity (\code{pgv}).
#' Typically, a 10 percent of the \code{pgv} is considered as a good indicator
#' of the cessation of the active pubertal growth spurt
#' \insertCite{Anna2022}{bsitar}. The percentage is controlled via the
#' \code{acg_velocity} argument which takes a positive real value bounded
#' between 0 and 1 (default \code{0.1} implying 10 percent). 
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @param datagrid Generate a grid of user-specified values for use in the
#'   \code{newdata} argument in various functions of the \pkg{marginaleffects}
#'   package. This is useful to define where in the predictor space we want to
#'   evaluate the quantities of interest. See [marginaleffects::datagrid()] for
#'   details. The default value for the \code{datagrid} is \code{NULL} implying
#'   that no custom grid is constructed. To set a data grid, the argument should
#'   be a data.frame constructed by using the [marginaleffects::datagrid()]
#'   function, or else a named list which are internally used for setting up the
#'   grid. For the user convenience, we also allow setting an empty list
#'   \code{datagrid = list()} in which case essential arguments such as
#'   \code{model}, \code{newdata} are taken up from the respective arguments
#'   specified elsewhere. Further, the level 1 predictor (such as age) and any
#'   covariate included in the model fit (e.g., gender) are also automatically
#'   inferred from the \code{model} object.
#' 
#' @param parameter A single character string, or a character vector specifying
#'   the growth parameter(s) to be estimated. Options are \code{'tgv'} (takeoff
#'   growth velocity), \code{'atgv'} (age at takeoff growth velocity),
#'   \code{'pgv'} (peak growth velocity), \code{'apgv'} (age at peak growth
#'   velocity), \code{'cgv'} (cessation growth velocity), and \code{'acgv'} (age
#'   at cessation growth velocity), and \code{'all'}. If \code{parameter = NULL}
#'   (default), age at peak growth velocity (\code{'apgv'}) is estimated where
#'   when \code{parameter = 'all'}, all six parameters are estimated. Note that
#'   option \code{'all'} can not be used when argument \code{by} is \code{TRUE}.
#' 
#' @param acg_velocity A real number to set the percentage of peak growth growth
#'   velocity as the cessation velocity when estimating the \code{cgv} and
#'   \code{acgv} growth parameters. The \code{acg_velocity} should be greater
#'   than \code{0} and less than \code{1}. The default \code{acg_velocity =
#'   0.10} indicates that a 10 per cent of the peak growth velocity will be used
#'   to get the cessation velocity and the corresponding age at the cessation
#'   velocity. For example if peak growth velocity estimate is \code{10
#'   mm/year}, then cessation growth velocity is \code{1 mm/year}.
#' 
#' @param digits An integer (default \code{2}) to set the decimal places for the
#'   estimated growth parameters. The \code{digits} is passed on to the
#'   [base::round()] function.
#'   
#' @param average A logical to indicate whether to internally call the
#'    [marginaleffects::comparisons()] or the
#'    [marginaleffects::avg_comparisons()] function. If \code{FALSE} (default),
#'    [marginaleffects::comparisons()] is called otherwise
#'    [marginaleffects::avg_comparisons()] when \code{average = TRUE}.
#'
#' @param plot A logical to specify whether to plot comparisons by calling the
#'   [marginaleffects::plot_comparisons()] function (\code{FALSE}) or not
#'   (\code{FALSE}). If \code{FALSE} (default), then
#'   [marginaleffects::comparisons()] or [marginaleffects::avg_comparisons()]
#'   are called to compute predictions (see \code{average} for details).
#'   
#' @param showlegends An argument to specify whether to show legends
#'   (\code{TRUE}) or not (\code{FALSE}). If \code{NULL} (default), then
#'   \code{showlegends} is internally set to \code{TRUE} if \code{re_formula =
#'   NA}, and \code{FALSE} if \code{re_formula = NULL}. 
#'
#' @param variables For estimating growth parameters in the current use case,
#'   the \code{variables} is the level 1 predictor such as
#'   \code{age}/\code{time}. The \code{variables} is a named list where value is
#'   set via the \code{esp} argument (default 1e-6). If \code{NULL}, the
#'   \code{variables} is set internally by retrieving the relevant information
#'   from the \code{model}. Otherwise, user can define it as follows:
#'   \code{variables = list('x' = 1e-6)} where \code{'x'} is the level 1
#'   predictor. Note that \code{variables = list('age' = 1e-6)} is the default
#'   behavior for the \pkg{marginaleffects} because velocity is typically
#'   calculated by differentiating the distance curve via \code{dydx} approach,
#'   and therefore argument \code{deriv} is automatically set as \code{0} and
#'   \code{deriv_model} as \code{FALSE}. If user want to estimate parameters
#'   based on the model based first derivative, then argument \code{deriv} must
#'   be set as \code{1} and internally argument \code{variables} is defined as
#'   \code{variables = list('age' = 0)} i.e, original level 1 predictor
#'   variable, \code{'x'}. It is important to consider that if default behavior
#'   is used i.e, \code{deriv = 0} and \code{variables = list('x' = 1e-6)}, then
#'   user can not pass additional arguments to the \code{variables} argument. On
#'   the other hand, alternative approach i.e, \code{deriv = 0} and
#'   \code{variables = list('x' = 0)}, additional options can be passed to the
#'   [marginaleffects::comparisons()] and [marginaleffects::avg_comparisons()]
#'   functions.
#'   
#' @param deriv A numeric to specify whether to estimate parameters based on the
#'   differentiation of the distance curve or the model based first derivative.
#'   Please see argument \code{variables} for more details.  
#' 
#' @param comparison For estimating growth parameters in the current use case,
#'   options allowed for the \code{comparison} are \code{'difference'} and
#'   \code{'differenceavg'}. Note that \code{comparison} is a placeholder and is
#'   only used to setup the the internal function that estimates
#'   \code{'parameter'} via [sitar::getPeak()], [sitar::getTakeoff()] and
#'   [sitar::getTrough()] functions to estimate various growth parameters.
#'   Options \code{'difference'} and \code{'differenceavg'} are internally
#'   restructured according to the user specified \code{hypothesis} argument.
#'   
#' @param reformat A logical (default \code{TRUE}) to reformat the  output
#'   returned by the \code{marginaleffects} as a data.frame with column names
#'   re-defined as follows: \code{conf.low} as \code{Q2.5}, and \code{conf.high}
#'   as \code{Q97.5} (assuming that \code{conf_int = 0.95}). Also, following
#'   columns are dropped from the data frame: \code{term}, \code{contrast},
#'   \code{tmp_idx}, \code{predicted_lo}, \code{predicted_hi}, \code{predicted}.
#'   
#' @param estimate_center A character string (default \code{NULL}) to specify
#'   whether to center estimate as \code{'mean'} or as \code{'median'}. Note
#'   that \code{estimate_center} is used to set the global options as follows:
#'   \cr
#'   \code{ options("marginaleffects_posterior_center" = "mean")}, or \cr 
#'   \code{options("marginaleffects_posterior_center" = "median")} \cr
#'   The pre-specified global options are restored on exit via the
#'   [base::on.exit()].
#' 
#' @param estimate_interval A character string (default \code{NULL}) to specify
#'   whether to compute credible intervals as equal-tailed intervals,
#'   \code{'eti'} or highest density intervals, \code{'hdi'}. Note that
#'   \code{estimate_interval} is used to set the global options as follows: \cr
#'   \code{ options("marginaleffects_posterior_interval" = "eti")}, or \cr
#'   \code{ options("marginaleffects_posterior_interval" = "hdi")} \cr
#'   The pre-specified global options are restored on exit via the
#'   [base::on.exit()].
#' 
#' @inheritParams  growthparameters.bgmfit
#' @inheritParams  marginaleffects::comparisons
#' @inheritParams  marginaleffects::avg_comparisons
#' @inheritParams  marginaleffects::plot_comparisons
#' @inheritParams  marginaleffects::datagrid
#' @inheritParams  brms::fitted.brmsfit
#'
#' @return A data frame objects with estimates and CIs for computed parameter(s)
#' 
#' @export growthparameters_comparison.bgmfit
#' @export
#' 
#' @seealso [marginaleffects::comparisons()]
#'   [marginaleffects::avg_comparisons()]
#'   [marginaleffects::plot_comparisons()]
#' 
#' @references
#' \insertAllCited{}
#' 
#' @inherit berkeley author
#' 
#' @examples
#' 
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' growthparameters_comparison(model, parameter = 'apgv', ndraws = 10)
#' }
#' 
growthparameters_comparison.bgmfit <- function(model,
                                   resp = NULL,
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
                                   cores = NULL,
                                   average = FALSE, 
                                   plot = FALSE, 
                                   showlegends = NULL, 
                                   variables = NULL,
                                   deriv = NULL,
                                   deriv_model = NULL,
                                   comparison = "difference",
                                   type = NULL,
                                   by = FALSE,
                                   conf_level = 0.95,
                                   transform = NULL,
                                   cross = FALSE,
                                   wts = NULL,
                                   hypothesis = NULL,
                                   equivalence = NULL,
                                   eps = NULL,
                                   reformat = NULL,
                                   estimate_center = NULL,
                                   estimate_interval = NULL,
                                   dummy_to_factor = NULL, 
                                   verbose = FALSE,
                                   expose_function = FALSE,
                                   usesavedfuns = NULL,
                                   clearenvfuns = NULL,
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
  
  
  
  
  try(zz <- insight::check_if_installed(c("marginaleffects"), 
                                        minversion = 
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
  
  
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- parent.frame()
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
  
  # Note, newdata is not model$data but rather model$model_info$bgmfit.data
  # This is must for univariate_by
  if(is.null(newdata)) {
    #newdata <- model$model_info$bgmfit.data
  }
  
  if(!is.na(uvarby)) {
    uvarby_ind <- paste0(uvarby, resp)
    varne <- paste0(uvarby, resp)
   # newdata <- newdata %>% dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
  }
  
  
  if(is.null(deriv) & is.null(deriv_model)) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(deriv == 0 & is.null(deriv_model)) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(deriv == 1 & is.null(deriv_model)) {
    deriv <- 1
    deriv_model <- TRUE
  } else if(is.null(deriv) & !deriv_model) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(is.null(deriv) & deriv_model) {
    deriv <- 1
    deriv_model <- TRUE
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
  
  
  allowed_parms <- c(
    'apgv',
    'pgv',
    'atgv',
    'tgv',
    'acgv',
    'cgv'
    )
  
  
  if (is.null(parameter)) {
    parm <- 'apgv' 
  } else if(parameter == 'all') {
    parm <- allowed_parms
  } else {
    if(!parameter %in% allowed_parms) {
      allowed_parms_err <- c(allowed_parms, 'all')
      stop("Allowed parameter options are ", 
           paste(paste0("'", allowed_parms_err, "'"), collapse = ", ")
      )
    }
    parm <- parameter
  }
  
  if(length(parm) > 1) {
    if(plot) stop("Please specify only one parameter when plot = TRUE")
  }
  
  
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
                             minversion = get_package_minversion('brms'), 
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
    if(any(grepl("growthparameters_comparison", scall, fixed = T)) |
       any(grepl("growthparameters_comparison.bgmfit", scall, fixed = T))) {
      xcall <- "growthparameters_comparison"
    } else {
      xcall <- xcall
    } 
  }
  
  if(!is.null(model$xcall)) {
    if(model$xcall == "growthparameters_comparison") {
      xcall <- "growthparameters_comparison"
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
  arguments$cores <- setincores <-  get.cores_[['max.cores']]
  .cores_ps <- get.cores_[['.cores_ps']]
  
  if (future) {
    if (future_session == 'multisession') {
      future::plan('multisession', workers = setincores)
    } else if (future_session == 'multicore') {
      future::plan('multicore', workers = setincores)
    }
  }
  
  
  
  full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                  # fargs = formals(), 
                                  fargs = arguments, 
                                  dargs = list(...), 
                                  verbose = verbose)
  
  full.args$model <- model
  full.args$deriv_model <- deriv_model
  
  if(is.null(full.args$hypothesis) && is.null(full.args$equivalence)) {
    plot <- plot
  } else {
    plot <- FALSE
    if(verbose) {
      message("Argument plot = TRUE is not allowed when either hypothesis ", 
              "or equivalence is not NULL",
              "\n ",
              "Therefor, argument setting plot = FALSE") 
    }
  }
  
  full.args$newdata <- newdata
  newdata           <- do.call(get.newdata, full.args)
  
  if(!is.na(uvarby)) {
    uvarby_ind <- paste0(uvarby, resp)
    varne <- paste0(uvarby, resp)
     newdata <- newdata %>% dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
  }
  full.args$newdata <- newdata
  

  # keeping ... cause marginaleffects:: argument is missing, with no default
  full.args[["..."]] <- NULL
  
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
      deriv_model,
      ipts,
      seed,
      future,
      future_session,
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
      reformat
    )
  ))[-1]
  
  if(plot) {
    exclude_args <- c(exclude_args, "cross")
  }
  
  for (exclude_argsi in exclude_args) {
    comparisons_arguments[[exclude_argsi]] <- NULL
  }
  
  
 
  
  if(deriv == 0 & deriv_model) 
    stop("If deriv_model = TRUE, deriv should be 1")
  if(deriv == 1 & !deriv_model) 
    stop("If deriv_model = FALSE, deriv should be 0")
  
  if (!is.null(variables)) {
    if (!is.list(variables)) {
      if(!deriv_model) {
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
          ec_agg <- getOption("marginaleffects_posterior_center")
          if(is.null(ec_agg)) ec_agg <- "mean"
          if(ec_agg == "mean")   xy <- stats::aggregate(.~x, data=xy, 
                                                        mean, drop = TRUE)
          if(ec_agg == "median") xy <- stats::aggregate(.~x, data=xy, 
                                                        median, drop = TRUE)
        }
        x <- xy$x
        y <- xy$y
      } # if(aggregate_by) {
      
      
      
      if (parm == 'apgv') {
        out <- sitar::getPeak(x = x, y = y)[1]
      } else if (parm == 'pgv') {
        out <- sitar::getPeak(x = x, y = y)[2]
      } else if (parm == 'atgv') {
        out <- sitar::getTakeoff(x = x, y = y)[1]
      } else if (parm == 'tgv') {
        out <- sitar::getTakeoff(x = x, y = y)[2]
      } else if (parm == 'acgv') {
        cgv  <- acg_velocity * sitar::getPeak(x = x, y = y)[2]
        vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
        out <-  x[vcgi]
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
        if(is.null(datagrid[['model']])) {
          setmodel <- model 
        } else {
          setmodel <- datagrid$model
        }
        if(is.null(datagrid[['newdata']])) {
          setnewdata <- newdata
        } else {
          setnewdata <- datagrid$newdata
        }
        if(is.null(datagrid[['grid_type']])) {
          setgrid_type <- "mean_or_mode"
        } else {
          setgrid_type <- datagrid$grid_type
        }
        if(is.null(datagrid[[xvar]])) {
          setxvar <- newdata[[xvar]] 
        } else {
          setxvar <- datagrid$newdata[[xvar]]
        }
        datagrid_arguments <- list(model = setmodel,
                                   newdata = setnewdata,
                                   grid_type = setgrid_type)
        datagrid_arguments[[xvar]] <- setxvar
        if(setgrid_type == "mean_or_mode") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- set_group
        } else if(setgrid_type == "balanced") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- NULL
          # correctly set comparisons_arguments[['by']] too 
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
    
    # The datagrid argument is not allowed. It served its purpose by defining 
    # the newdata. So remove it from the arguments
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
    
    suppressWarnings({
      if(!plot) {
        if(!average) {
          out <- do.call(marginaleffects::comparisons, 
                         comparisons_arguments)
        } else if(average) {
          out <- do.call(marginaleffects::avg_comparisons, 
                         comparisons_arguments)
        }
      }
      if(plot) {
        if(isFALSE(set_group)) comparisons_arguments$by <- NULL
        out <- do.call(marginaleffects::plot_comparisons, 
                       comparisons_arguments)
        outp <- out
        if(!showlegends) outp <- outp + ggplot2::theme(legend.position = 'none')
        return(outp)
      }
    })
    
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
      }
    )
    err. <- get('err.', envir = enverr.)
    
    # This happen when NA in one factor for hypothesis
    if(!exists('gout')) return(invisible(NULL))
    
    if(length(gout) == 0) err. <- TRUE
    
    # for (byi in eval(by)) {
    #   if(is.factor(newdata[[byi]])) {
    #     if(length(unique(gout[[byi]])) != 
    #        length(unique(newdata[[byi]])) ) {
    #       err. <- TRUE
    #     }
    #   }
    # }
    
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
  
  
  out_sf <- out_sf %>% 
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                         ~ round(., digits = digits))) %>% 
    dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
    data.frame()
  
  
  if(is.null(reformat)) {
    if(is.null(hypothesis) && is.null(equivalence)) {
      reformat <- TRUE
    } else {
      reformat <- FALSE
    }
  }
  
  if (reformat) {
    out_sf <- out_sf %>% 
      dplyr::rename(!!as.symbol(set_names_[1]) := estimate) %>% 
      dplyr::rename(!!as.symbol(set_names_[2]) := conf.low) %>% 
      dplyr::rename(!!as.symbol(set_names_[3]) := conf.high) 
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
  }
   
  
  return(out_sf)
}




#' @rdname growthparameters_comparison.bgmfit
#' @export
growthparameters_comparison <- function(model, ...) {
  UseMethod("growthparameters_comparison")
}


