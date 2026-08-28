


#' Diagnostic plots for Bayesian SITAR models
#'
#' @description
#' This function generate various model diagnostic plots including
#' residual-based diagnostics, MCMC convergence diagnostics, and posterior
#' predictive checks. When multiple plots are requested, patch-able plots are
#' combined into a single display using \pkg{patchwork}.
#'
#' @param model An object of class \code{bgmfit}
#'
#' @param newdata Optional \code{data.frame} used to compute residual-based
#'   diagnostics. If \code{NULL}, the function attempts to use the data stored
#'   in the fitted model object. Supplying \code{newdata} explicitly is safer
#'   and is recommended for reproducibility.
#' 
#' @param plots Character vector specifying which diagnostics to produce.
#'   Supported values are:
#'   \describe{
#'     \item{`"rvf"`}{residual vs fitted plot showing residuals plotted against 
#'     fitted values derived from posterior draws.}
#'     \item{`"rvp"`}{residual vs predictor (such as age) plot showing 
#'     residuals plotted against predictor values derived from posterior draws.}
#'     \item{`"qqn"`}{Normal Q-Q plot based on summarized posterior residuals.}
#'     \item{`"qqp"`}{Pearson Q-Q plot based on standardized posterior
#'     residuals.}
#'     \item{`"qq"`}{Same as \code{"qqp"}.}
#'     \item{`"pairs"`}{Pair plots for MCMC draws via
#'     [bayesplot::mcmc_pairs()].}
#'     \item{`"acfp"`}{Autocorrelation plots of parameters for MCMC draws via
#'     [brms::mcmc_plot()] with \code{type = "acf"}.}
#'     \item{`"acfr"`}{Autocorrelation plots of residuals via
#'     custom acf_residuals_custom function.}
#'     \item{`"acf"`}{Same as \code{"acfr"}.}
#'     \item{`"trace"`}{Trace plots for selected parameters via
#'     [brms::mcmc_plot()] with \code{type = "trace"}.}
#'     \item{`"dens_overlay"`}{Density overlays for selected parameters via
#'     [brms::mcmc_plot()] with \code{type = "dens_overlay"}.}
#'     \item{`"rhat"`}{Potential scale reduction diagnostics via
#'     [brms::mcmc_plot()] with \code{type = "rhat"}.}
#'     \item{`"rhat_hist"`}{Potential scale reduction diagnostics via
#'     [brms::mcmc_plot()] with \code{type = "rhat_hist"}.}
#'     \item{`"neff"`}{Effective sample size diagnostics via
#'     [brms::mcmc_plot()] with \code{type = "neff"}.}
#'     \item{`"ppc_overlay"`}{Posterior predictive density overlay via
#'     [brms::pp_check()] with \code{type = "dens_overlay"}.}
#'     \item{`"ppc_hist"`}{Posterior predictive histogram via
#'     [brms::pp_check()] with \code{type = "hist"}.}
#'     \item{`"ppc_scatter"`}{Posterior predictive observed-versus-replicated
#'     average scatter plot via [brms::pp_check()] with \code{type =
#'     "scatter_avg"}.}
#'     \item{`"ppc_stat"`}{Posterior predictive check for a summary statistic 
#'     via [brms::pp_check()] with \code{type = "stat"}.}
#'   }
#'
#' @param set_draws Character scalar indicating which fitted-value draw type to
#'   use for the residual-versus-fitted plot. Must be one of:
#'   \describe{
#'     \item{`"epred"`}{Expected posterior predictive values from
#'     [tidybayes::add_epred_draws()]. This is usually the preferred option when
#'     available.}
#'     \item{`"linpred"`}{Linear predictor values from
#'     [tidybayes::add_linpred_draws()].}
#'     \item{`"prediction"`}{Posterior predictive draws from
#'     [tidybayes::add_predicted_draws()].}
#'   }
#'   The \pkg{tidybayes} documentation recommends using expectation-scale
#'   predictions from [tidybayes::add_epred_draws()] when available, because
#'   these correspond to the expectation of the posterior predictive
#'   distribution more directly than the linear predictor in many models.
#'
#' @param resp Optional name of the response variable to plot. This is
#'   primarily useful for multivariate models, where a specific response
#'   must often be selected.
#'
#' @param ndraws Optional integer specifying the number of posterior draws to
#'   use in residual-based diagnostics and posterior predictive checks.
#'
#' @param seed Optional random seed used when subsampling posterior draws.
#'
#' @param re_formula Optional formula passed to prediction
#'   methods to control whether group-level effects are included.
#'
#' @param dpar Optional distributional parameter name passed to the underlying
#'   draw function when relevant.
#'
#' @param variable Optional character vector of parameter names for MCMC
#'   diagnostic plots such as trace, autocorrelation, and density overlays.
#'
#' @param combine Logical; if \code{TRUE} and more than one plot is generated,
#'   compatible plots are combined with [patchwork::wrap_plots()]. If
#'   \code{FALSE}, the function returns the individual plot objects without
#'   combining them.
#'
#' @param ncol Optional integer giving the number of columns to use when
#'   combining multiple plots with \pkg{patchwork}. If \code{NULL}, a simple
#'   default is chosen automatically.
#'
#' @param point_alpha Numeric transparency level used for points in the
#'   residual-versus-fitted plot.
#'   
#' @param bins  Numeric value passed to \code{'stat_bin()'} for
#'   \code{'ppc_hist'} and  \code{'ppc_stat'} plots.
#'
#' @param smooth_se Logical; should the residual smoothing line display a
#'   standard-error ribbon.
#'
#' @param smooth_method Smoother passed to [ggplot2::geom_smooth()] for the
#'   residual-versus-fitted plot. The default is \code{"lm"}.
#'
#' @param ppc_stat Summary statistic to use when \code{plots} includes
#'   \code{"ppc_stat"}. This value is passed to [brms::pp_check()] as
#'    \code{type = "stat", stat = 'ppc_stat'}.
#' 
#' @param rank_overlay Logical; if \code{TRUE}, sets the \pkg{bayesplot} color
#'   scheme to \code{"blue"}. This does not itself create a rank-overlay plot,
#'   but allows a simple style customization for compatible
#'   \pkg{bayesplot}-based outputs.
#'   
#' @param add_plot_df Optional logical (default \code{FALSE}) indicating whether
#'   to add \code{plot_df} to the returned value.
#'   
#' @param wrap_title Optional logical indicating whether to wrap plot title. If
#'   \code{NULL} (default), then \code{wrap_title} is internally set as
#'   \code{TRUE} when three or more plots are requested.
#'   
#' @param title_size Text size of plot title.
#' 
#' @param qq_plot_args A named list of arguments passed on to the
#'   \code{`"qq_plot_pearson"`}. Ignored except when \code{plots = `"qqp"`} or
#'   \code{plots = `"qq"`}. Note that \code{plots = `"qq"`} is internally
#'   evaluated as \code{plots = `"qqp"`}.
#'   
#' @param patch_plot.margin Optional setting the margins for the
#'   \pkg{patchwork}.
#'   
#' @param funlist Currently ignored and serves as a placeholder. It will be
#'   implemented in future updates.
#' @param future Currently ignored and serves as a placeholder.
#' @param future_session Currently ignored and serves as a placeholder.
#' @param future_splits Currently ignored and serves as a placeholder.
#' @param future_method Currently ignored and serves as a placeholder.
#' @param future_re_expose Currently ignored and serves as a placeholder.
#' @param transform Currently ignored and serves as a placeholder.
#' @param transform_draws Currently ignored and serves as a placeholder.
#' @param itransform Currently ignored and serves as a placeholder.
#' @param model_deriv Currently ignored and serves as a placeholder.
#' @param dummy_to_factor Currently ignored and serves as a placeholder.
#' @param xvar Currently ignored and serves as a placeholder.
#' @param difx Currently ignored and serves as a placeholder.
#' @param idvar Currently ignored and serves as a placeholder.
#' @param conf_level Currently ignored and serves as a placeholder.
#' @param plot Currently ignored and serves as a placeholder.
#' @param method Currently ignored and serves as a placeholder.
#' @param deriv Currently ignored and serves as a placeholder.
#' @param idata_method Currently ignored and serves as a placeholder.
#' @param ipts Currently ignored and serves as a placeholder.
#' @param newdata_fixed Currently ignored and serves as a placeholder.
#'   
#' @inheritParams plot_curves
#' @inheritParams add_model_criterion
#' @inheritParams get_predictions
#' @inheritParams brms::mcmc_plot
#' @inheritParams tidybayes::add_epred_draws
#'
#' @returns
#' A list. The returned structure depends on \code{combine} and the number of
#' plots requested:
#' \describe{
#'   \item{If `combine = FALSE` or only one plot is requested}{A named list of
#'   plot objects.}
#'   \item{If `combine = TRUE` and multiple compatible plots are requested}{A
#'   list with components:
#'     \describe{
#'       \item{`combined`}{A patchwork object combining all \code{patchable}
#'       plots.}
#'       \item{`plots`}{A named list containing all individual plot objects.}
#'       \item{`plot_df`}{The joined residual/fitted draws data frame used for
#'       residual-based plots, or \code{NULL} if no residual-based plots were
#'       requested.}
#'     }
#'   }
#' }
#'
#' @details
#' This function provides a single interface to several common Bayesian model
#' diagnostics. Residual-based plots are constructed with \pkg{tidybayes}, using
#' posterior residual draws together with fitted-value draws from
#' [tidybayes::add_epred_draws()], [tidybayes::add_linpred_draws()], or
#' [tidybayes::add_predicted_draws()] MCMC diagnostics are produced through
#' \code{'mcmc_plot'} which is a wrapper around plotting functionality from
#' \pkg{bayesplot}. Posterior predictive checks are produced through
#' [brms::pp_check()]
#' 
#' Residual-based diagnostics are created by joining posterior residual draws
#' from [tidybayes::add_residual_draws()] with fitted-value draws from one of
#' [tidybayes::add_epred_draws()], [tidybayes::add_linpred_draws()], or
#' [tidybayes::add_predicted_draws()]. The resulting joined draws are used for
#' the residual-versus-fitted plot, and the residuals are summarized with
#' [tidybayes::median_qi()] before creating the Q-Q plot. The \pkg{tidybayes}
#' documentation notes that the default output column names depend on the draw
#' function: \code{".epred"}, \code{".linpred"}, \code{".prediction"}, and
#' \code{".residual"}.
#' 
#' MCMC diagnostic plots are delegated to \code{'mcmc_plot'} which provides a
#' convenient wrapper around \pkg{bayesplot} diagnostics such as trace plots,
#' autocorrelation plots, density overlays, R-hat, and effective sample size
#' summaries. 
#'
#' Posterior predictive checks are delegated to [brms::pp_check()], which
#' exposes \pkg{bayesplot} \code{PPC} types such as density overlays,
#' histograms, average scatter checks, and statistic-based checks.
#'
#' @section Interpretation:
#' - For Gaussian models, residual-versus-fitted plots and Q-Q
#'   plots are often useful summaries of residual behavior.
#' - \code{set_draws = "epred"} is usually the best default because expected
#'   predictions are often easier to interpret than linear predictors or raw
#'   posterior predictive draws. 
#'
#' @seealso
#' [brms::pp_check()], [brms::mcmc_plot()],
#' [tidybayes::add_residual_draws()], [tidybayes::add_epred_draws()],
#' [tidybayes::add_linpred_draws()], [tidybayes::add_predicted_draws()],
#' [patchwork::wrap_plots()]
#' 
#' @rdname plot_diagnostics
#' @export
#' 
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid model estimation which can take time, the Bayesian SITAR model fit
#' # to the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' model <- getNsObject(berkeley_exfit)
#' 
#' # Residual vs Fitted plot
#' res1 <- plot_diagnostics(
#'   model = model,
#'   plots = "rvf"
#' )
#' 
#'
#' # MCMC diagnostics for selected parameters
#' res2 <- plot_diagnostics(
#'   model = model,
#'   plots = c("trace", "dens_overlay"),
#'   variable = c("b_a_Intercept", "b_b_Intercept")
#' )
#' res2$combined
#'
#' # Posterior predictive checks
#' res3 <- plot_diagnostics(
#'   model = model,
#'   plots = c("ppc_overlay", "ppc_scatter"),
#'   ndraws = 10,
#'   ppc_stat = "sd"
#' )
#' res3$combined
#'
#' # Use linear predictor instead of expected predictions
#' res4 <- plot_diagnostics(
#'   model = model,
#'   plots = "rvf",
#'   set_draws = "linpred"
#' )
#' }
#'
plot_diagnostics.bgmfit <- function(
    model,
    newdata = NULL,
    plots = c(
      "rvf",
      "rvp",
      "qq",
      "qqn",
      "qqp",
      "pairs",
      "acf",
      "acfp",
      "acfr",
      "trace",
      "dens_overlay",
      "rhat",
      "rhat_hist",
      "neff",
      "ppc_overlay",
      "ppc_hist",
      "ppc_scatter",
      "ppc_stat"
    ),
    set_draws = "epred",
    resp = NULL,
    ndraws = NULL,
    draw_ids = NULL,
    seed = 123,
    re_formula = NULL,
    dpar = NULL,
    variable = NULL,
    combine = TRUE,
    ncol = NULL,
    point_alpha = 0.15,
    bins = 30,
    regex = FALSE,
    smooth_se = FALSE,
    smooth_method = "lm",
    ppc_stat = "mean",
    rank_overlay = FALSE,
    add_plot_df = FALSE,
    expose_function = FALSE,
    usesavedfuns = NULL,
    clearenvfuns = NULL,
    category = ".category",
    wrap_title = NULL,
    title_size = 12,
    each_object = FALSE,
    qq_plot_args = list(draw_ids = NULL, 
                        draw_ids_select = 1, 
                        summary = "mean", 
                        qq_type = 'qq',
                        seed = 123),
    patch_plot.margin = ggplot2::margin(10, 0, 10, 0),
    funlist = NULL,
    future = FALSE,
    future_session = 'multisession',
    future_splits = TRUE,
    future_method = 'future',
    future_re_expose = NULL,
    transform = NULL,
    transform_draws = NULL,
    itransform = NULL,
    model_deriv = NULL,
    dummy_to_factor = NULL, 
    xvar = NULL,
    difx = NULL,
    idvar = NULL,
    levels_id = NULL,
    conf_level = 0.95,
    deriv = 0,
    plot = FALSE,
    method = 'pkg',
    ipts = FALSE,
    idata_method = NULL,
    newdata_fixed = NULL,
    verbose = FALSE,
    envir = NULL,
    ...) {
  
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
  
  xvar_      <- paste0('xvar', resp_rev_)
  sigmaxvar_ <- paste0('sigma', xvar_)
  cov_       <- paste0('cov', resp_rev_)
  sigmacov_  <- paste0('sigma', cov_)
  uvarby     <- model$model_info$univariate_by$by
  yvar_      <- paste0('yvar', resp_rev_)
  yvar       <- model$model_info[[xvar_]]
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
  groupvar_     <- paste0('groupvar', resp_rev_)
  yvar_         <- paste0('yvar', resp_rev_)
  yvar          <- model$model_info[[yvar_]]
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
  funx_ <- NULL
  
  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  
  if(idata_method == 'm1') {
    stop("For marginaleffects based functions, the " ,
         " \n",
         " 'idata_method' argument must be either NULL or 'm2'" )
  }
  
  .draw <- NULL;
  .fitted_value <- NULL;
  .residual <- NULL;
  mean_residual <- NULL;
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
  set_names_  <- c('Estimate', probtitles)
  
  expose_method_set <- model$model_info[['expose_method']]
  model$model_info[['expose_method']] <- 'NA' 
  
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
  post_processing_checks_args[['deriv']] <- 0
  if(is.null(model_deriv)) model_deriv <- TRUE
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
    } 
  } 
  
  if(!is.null(transform)) {
    if(!is.function(transform)) {
      if(is.logical(transform)) {
        if(!transform) transform_draws <- 'identity'
      } else if(!is.logical(transform)) {
        if(transform == "exp") transform_draws <- 'exp'
        if(transform == "ln") transform_draws <- 'log'
      }
    } 
  } 
  
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
    }
  }

  if(dpar == "sigma") {
    if(deriv.org > 0) {
      if(!is.null(o[['sigma_model']])) {
        if(o[['sigma_model']] == "ls") {
          # 
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
          } 
        } 
      } 
    } 
  } 
  
  call_predictions <- TRUE
  call_slopes      <- FALSE
  if(!model_deriv) {
    if(deriv > 0) {
      deriv <- 0
      call_predictions <- FALSE
      call_slopes      <- TRUE
    }
  } 
  
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
 
  force_condition_and_by_switch_plot <- FALSE
  if(dpar == "sigma") {
    if(deriv.org > 0) {
      if(!is.null(o[['sigma_model']])) {
        if(o[['sigma_model']] == "ls") {
          # 
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
      } 
    } 
  } 
  
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
    if(grepl("plot_diagnostics", model$xcall)) {
      xcall <- "plot_diagnostics"
    } else if(grepl("plot_diagnostics", model$xcall)) {
      xcall <- "plot_diagnostics"
    } else if(grepl("plot_diagnostics", model$xcall)) {
      xcall <- "plot_diagnostics"
    }
  } else {
    rlang_trace_back <- rlang::trace_back()
    check_trace_back.bgmfit <- grepl(".bgmfit", rlang_trace_back[[1]])
    if(all(!check_trace_back.bgmfit)) {
      # 
    } else {
      rlang_trace_back.bgmfit_i <- 
        min(which(check_trace_back.bgmfit == TRUE))
      rlang_trace_back.bgmfit <- 
        rlang_trace_back[[1]][[rlang_trace_back.bgmfit_i]]
      rlang_call_name <- 
        rlang::call_name(rlang_trace_back.bgmfit)
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
    if (inherits(future::plan(), "multicore")) {
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
    } 
    if(method == 'pkg') {
      future_splits_exe <- FALSE
      future_splits     <- NULL
      future_splits_at  <- NULL
      if(verbose) {
        message2c("future_splits can not be used when method = 'pkg'.
               future_splits set as FALSE")
      }
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
          if(verbose) {
            message2c("To speed up the calulations, it is advised to 
                        set 'future_re_expose' = TRUE")
          }
        } 
        if(need_future_re_expose_cpp & setplanis == "multisession") {
          stop2c("For plan multisession, the functions need to be 
                   re_exposed by setting 'future_re_expose' = TRUE")
        }
      }
    }
  } 

  if (!future) {
    future_splits_at <- NULL
    future_splits_exe <- FALSE
    future_splits_exe_future <- FALSE
    future_splits_exe_dofuture <- FALSE
  }
  
  full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                  # fargs = arguments,
                                  fargs = formals(),
                                  dargs = list(...), 
                                  sanitize_CustomDoCall_args = FALSE,
                                  check_formalArgs  = NULL,
                                  check_formalArgs_exceptions  = NULL,
                                  check_trace_back  = NULL,
                                  envir = envir,
                                  verbose = verbose)

  full.args$model       <- model
  full.args$model_deriv <- model_deriv
  full.args$newdata     <- newdata

  full.args <- 
    sanitize_CustomDoCall_args(what = "CustomDoCall", 
                               arguments = full.args, 
                               check_formalArgs = plot_diagnostics.bgmfit,
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
  }
  
  full.args$newdata <- newdata
  full.args[["..."]] <- NULL
  
  if(!is.null(full.args[['transform_draws']])) {
    full.args[['transform']] <- transform <- full.args[['transform_draws']]
    if(verbose) message("'transform' set based on 'transform_draws'")
  } else if(!is.null(transform_draws)) {
    full.args[['transform']] <- transform <- transform_draws
    if(verbose) message("'transform' set based on 'transform_draws'")
  } 
  
  full.args[['transform']] <- transform <- transform_draws
  full.args$model_deriv <- model_deriv

  allowed_set_draws <- c("epred", "linpred", "prediction")
  if (length(set_draws) != 1 || !set_draws %in% allowed_set_draws) {
    stop("`set_draws` must be one of: ", 
         paste(allowed_set_draws, collapse = ", "))
  }
  
  missing_pkgs <- c("ggplot2", "dplyr", "patchwork", "bayesplot")[
    !vapply(c("ggplot2", "dplyr", "patchwork", "bayesplot"),
            requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0) {
    stop("Please install required packages: ", 
         paste(missing_pkgs, collapse = ", "))
  }
  
  enverr. <- environment()
  
  needs_tidybayes <- any(c("rvf", "rvp", "qq", "qqn", "qqp") %in% plots)
  if (needs_tidybayes && !requireNamespace("tidybayes", quietly = TRUE)) {
    stop("Please install `tidybayes` for residual-based plots.")
  }

  plots <- unique(plots)
  valid_plots <- c(
    "rvf", "rvp", "qq", "qqn", "qqp", "pairs", "acf", "acfp", "acfr", "trace", 
    "dens_overlay", "rhat", "rhat_hist", "neff", "ppc_overlay", "ppc_hist",
    "ppc_scatter", "ppc_stat"
  )
  
  bad <- setdiff(plots, valid_plots)
  if (length(bad) > 0) {
    stop("Unknown plot type(s): ", paste(bad, collapse = ", "))
  }

  if ("qq" %in% plots) {
    plots <- gsub("^qq$", "qqp", plots, fixed = F)
    plots <- unique(plots)
  }
  
  if ("acf" %in% plots) {
    plots <- gsub("^acf$", "acfr", plots, fixed = F)
    plots <- unique(plots)
  }

  out <- list()
  
  if (any(c("rvf", "rvp", "qqn", "qqp") %in% plots)) {
    draw_fun <- switch(
      set_draws,
      epred = tidybayes::add_epred_draws,
      linpred = tidybayes::add_linpred_draws,
      prediction = tidybayes::add_predicted_draws
    )
    resi_fun <- tidybayes::add_residual_draws
    draw_col <- paste0(".", set_draws)
    full.args[['newdata']] <- full.args[['newdata']] %>% 
      dplyr::mutate(.row = dplyr::row_number())
    full.args[['object']] <- full.args[['model']]
    draw_fun_args <- resi_fun_args <- list()
    resi_fun_args_names <- c('object', 'newdata', 'resp', 'ndraws',
                             'draw_ids', 'seed', 're_formula', 'category')
    draw_fun_args_names <- c('object', 'newdata', 'resp', 'ndraws',
                             'draw_ids', 'seed', 're_formula', 'category', 
                             'dpar')
    for (draw_fun_args_namesi in draw_fun_args_names) {
      draw_fun_args[[draw_fun_args_namesi]] <- full.args[[draw_fun_args_namesi]] 
    }
    for (resi_fun_args_namesi in resi_fun_args_names) {
      resi_fun_args[[resi_fun_args_namesi]] <- full.args[[resi_fun_args_namesi]] 
    }
    resi_fun_args[['value']]   <- ".residual"
    draw_fun_args[['value']] <- draw_col
    fitted_df <- CustomDoCall(draw_fun, draw_fun_args) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(.row, .draw, dplyr::all_of(draw_col)) %>%
      dplyr::rename(.fitted_value = dplyr::all_of(draw_col))
    resid_df <- CustomDoCall(resi_fun, resi_fun_args) 
    plot_df <- dplyr::left_join(
      resid_df,
      fitted_df,
      by = c(".row", ".draw")
    )
    
    if ("rvf" %in% plots) {
      xlab_txt <- switch(
        set_draws,
        epred = "Expected predicted values",
        linpred = "Linear predicted values",
        prediction = "Posterior predicted values"
      )
      out$rvf <- plot_df %>%
        ggplot2::ggplot(ggplot2::aes(x = .fitted_value, y = .residual)) +
        ggplot2::geom_point(alpha = point_alpha, color = "gray40") +
        ggplot2::geom_hline(yintercept = 0, linetype="dashed", color = "red") +
        ggplot2::geom_smooth(
          method = smooth_method,
          formula = y ~ x,
          se = smooth_se,
          color = "blue"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          x = xlab_txt,
          y = "Residuals",
          title = "Residuals vs Fitted"
        )
    }
    
    if ("qqn" %in% plots) {
      qq_df <- plot_df %>%
        tidybayes::median_qi(.residual)
      out$qqn <- qq_df %>%
        ggplot2::ggplot(ggplot2::aes(sample = .residual)) +
        ggplot2::geom_qq() +
        ggplot2::geom_qq_line() +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Q-Q Plot of Residuals") +
        ggplot2::labs(x = "Theoretical Quantiles", y = "Residuals")
    }
    
    if ("qqp" %in% plots) {
      if(!is.list(qq_plot_args)) qq_plot_args <- list()
      qq_plot_args[['model']]        <- model
      qq_plot_args[['data']]         <- full.args[['newdata']]
      qq_plot_args[['resid_draws']]  <- resid_df
      if(is.null(qq_plot_args[['draw_ids']])) 
        qq_plot_args[['draw_ids']] <- 1:pmin(brms::ndraws(model), 10)
      if(is.null(qq_plot_args[['draw_ids_select']])) 
        qq_plot_args[['draw_ids_select']] <- 1
      if(is.null(qq_plot_args[['summary']])) 
        qq_plot_args[['summary']] <- 'mean'
      if(is.null(qq_plot_args[['qq_type']])) 
        qq_plot_args[['qq_type']] <- 'qq'
      if(is.null(qq_plot_args[['seed']])) 
        qq_plot_args[['seed']] <- 123
      out$qqp <- do.call(qq_plot_pearson, qq_plot_args)
    }
  } 
  
  if (any(c( "rvp") %in% plots)) {
    resid_df_mean <- resid_df %>% 
      dplyr::group_by(.row) %>%
      dplyr::summarise(
        mean_residual = mean(.residual) ,
        .groups = "drop"
      ) %>%
      dplyr::bind_cols(full.args[['newdata']] %>% 
                         dplyr::select(dplyr::all_of(xvar))) 
    resid_df_mean[[xvar]] <- ifunx_(resid_df_mean[[xvar]])
    
    out$rvp <- resid_df_mean %>% 
      ggplot2::ggplot(ggplot2::aes(x = .data[['age']], y = mean_residual)) +
      ggplot2::geom_point(alpha = point_alpha, color = "gray40") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_smooth(
        method = smooth_method,
        formula = y ~ x,
        se = smooth_se,
        color = "blue"
        ) +
      ggplot2::labs(
        title = "Residuals vs. Predictor",
        x = xvar,
        y = "Average Residuals"
      ) +
      ggplot2::theme_minimal()
  }

  suppressMessages({
    if ("trace" %in% plots) {
      out$trace <- brms::mcmc_plot(
        model, type = "trace", variable = variable, regex = regex 
      )
      out$trace <- out$trace + 
      ggplot2::labs(title = "Trace plot") 
    }
  })

  if ("dens_overlay" %in% plots) {
    if(brms::nchains(model) > 1) {
      out$dens_overlay <- brms::mcmc_plot(
        model, type = "dens_overlay", variable = variable, regex = regex 
      )
      out$dens_overlay <- out$dens_overlay + 
        ggplot2::labs(title = "Density overlay plot") 
    } else {
      if(verbose) message("dens_overlay failed because of only one chain")
    }
  }

  if ("pairs" %in% plots) {
    insight::check_if_installed('ggpubr')
    brms_pairs_brmsfit <- utils::getFromNamespace('pairs.brmsfit', 'brms')
    assign('err.', FALSE, envir = enverr.)
    tryCatch(
      expr = {
        suppressWarnings({
          pairs_obj <- brms_pairs_brmsfit(
            x = model, variable = variable, regex = regex 
          )
        })
      },
      error = function(e) {
        assign('err.', TRUE, envir = enverr.)
      }
    )
    err. <- get('err.', envir = enverr.)
    if (err.) {
      if(verbose) message("pairs plot failed")
      out$pairs <- NULL
    } else {
      out$pairs <- pairs_obj 
      out$pairs <- ept("ggpubr::as_ggplot(out$pairs)")
      out$pairs <- out$pairs +
        ggplot2::labs(title = "Pairs plot")
    }
  }
  
  if ("acfp" %in% plots) {
    if(brms::ndraws(model) > 20) {
      out$acfp <- brms::mcmc_plot(
        model, type = "acf", variable = variable, regex = regex 
      )
      out$acfp <- out$acfp + 
        ggplot2::labs(title = "Autocorrelation plot") 
    } else {
      if(verbose) message("acf for parameters failed because of < 20 iter")
    }
  }
  
  if ("acfr" %in% plots) {
    out$acfr <- acf_residuals_custom(model,
                                     residual_type = c("raw", "standardized"),
                                     summary = c("mean", "median"),
                                     re_formula = re_formula,
                                     ndraws = ndraws,
                                     draw_ids = draw_ids,
                                     idvar = idvar,
                                     xvar = xvar,
                                     sort_residuals = TRUE,
                                     lag_max = NULL,
                                     plot_type = c("ggplot", "base", "none"),
                                     ci_level = NULL,
                                     acf_cutoffs = NULL,
                                     acf_cutoffs_pm = FALSE,
                                     print = FALSE,
                                     main_title = "ACF of Residuals")
    out$acfr <- out$acfr + 
      ggplot2::labs(title = "Autocorrelation of Residuals") 
  }
  
  if ("rhat" %in% plots) {
    suppressWarnings({ 
      suppressMessages({ 
        out$rhat <- brms::mcmc_plot(model, type = "rhat")
        out$rhat <- out$rhat + 
          ggplot2::labs(title = "Rhat plot") 
      })
    })
  }
  
  if ("rhat_hist" %in% plots) {
    suppressWarnings({ 
      suppressMessages({
      out$rhat_hist <- brms::mcmc_plot(model, type = "rhat_hist", bins  = bins)
      out$rhat_hist <- out$rhat_hist + 
        ggplot2:: scale_x_continuous(labels = scales::label_number(accuracy = 0.001),
                                     limits = 
                                       c(1, max(out$rhat_hist$data$value, 
                                                na.rm = TRUE))) +
        ggplot2:: scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
        ggplot2::labs(title = "Rhat plot") 
    })
      })
  }
  
  if ("neff" %in% plots) {
    suppressWarnings({ 
      suppressMessages({
        out$neff <- brms::mcmc_plot(model, type = "neff")
        out$neff <- out$neff + 
          ggplot2::labs(title = "Effective sample size plot") 
      })
    })   
  }
  
 
  if ("ppc_overlay" %in% plots) {
    if(brms::nchains(model) > 0) {
      out$ppc_overlay <- brms::pp_check(
        model, type = "dens_overlay", ndraws = ndraws, size = 2,  alpha = 4
      ) +
        ggplot2::theme(
          legend.position = "inside",
          legend.position.inside = c(0.15, 0.9)
        ) + 
        ggplot2::labs(title = "Posterior Predictive Check: Density Overlay")
    } else {
      if(verbose) message("ppc_overlay failed because of only one chain")
    }
  }

  if ("ppc_hist" %in% plots) {
    out$ppc_hist <- brms::pp_check(
      model, type = "hist", ndraws = ndraws, bins  = bins
    ) +
      ggplot2::labs(title = "Posterior Predictive Check: Histogram")
  }
  
  if ("ppc_scatter" %in% plots) {
    out$ppc_scatter <- brms::pp_check(
      model, type = "scatter_avg", ndraws = ndraws, size = 1
    ) +
      ggplot2::labs(title = "Posterior Predictive Check: Scatter Average")
  }
  
  suppressMessages({
    if ("ppc_stat" %in% plots) {
      out$ppc_stat <- brms::pp_check(
        model, type = "stat", stat = ppc_stat, ndraws = ndraws, bins  = bins
      ) +
        ggplot2::labs(title = paste("Posterior Predictive Check:", ppc_stat))
    }
  })

  for (namespi in names(out)) {
    out[[namespi]] <-  out[[namespi]] +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size))
  }
  
  if(is.null(wrap_title)) {
    if(length(plots) >=3) wrap_title <- TRUE else wrap_title <- FALSE
  }

  if(!wrap_title) {
    for (namespi in names(out)) {
      out[[namespi]] <-  out[[namespi]] +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
  } else if(wrap_title) {
    insight::check_if_installed('ggtext')
    for (namespi in names(out)) {
      out[[namespi]] <-  out[[namespi]] +
        ept(
        "ggplot2::theme(plot.title = 
                         ggtext::element_textbox_simple(
                           hjust = 0,
                           vjust = 0,
                           margin = ggplot2::margin(b = 10),
                           halign = 0.5,
                           valign = 0.5))"
        ) 
    }
  }

  if (isTRUE(rank_overlay) && requireNamespace("bayesplot", quietly = TRUE)) {
    bayesplot::color_scheme_set("blue")
  }
  
  if (!combine || length(out) <= 1) {
    if(length(out) == 1) return(out[[1]]) else return(out) 
  }

  is_patchable <- vapply(
    out,
    function(x) inherits(x, c("gg", "ggplot")),
    logical(1)
  )
  
  patchable_plots <- out[is_patchable]
  
  if (length(patchable_plots) == 0) {
    return(out)
  }
  
  if (is.null(ncol)) {
    ncol <- if (length(patchable_plots) <= 2) 1 else 2
  }

  patchable_plots <- patchable_plots[plots]
  
  if(is.null(patch_plot.margin)) {
    combined <- patchwork::wrap_plots(patchable_plots, ncol = ncol)
  } else if(!is.null(patch_plot.margin)) {
    combined <- patchwork::wrap_plots(patchable_plots, ncol = ncol) &
      ggplot2::theme(plot.margin = patch_plot.margin) 
  }
  
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
  }
  
  out_all <- list(combined = combined, plots = out)
  if(add_plot_df) {
    out_all[['plot_df']] = if (exists("plot_df")) plot_df else NULL
  }
  
  if(each_object) {
    return(out_all)
  } else {
    if(is.null(combined)) return(out) else return(combined)
  }
  
}



#' @rdname plot_diagnostics
#' @export
plot_diagnostics <- function(model, ...) {
  UseMethod("plot_diagnostics")
}


