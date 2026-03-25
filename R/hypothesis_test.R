

#' @title Comprehensive hypothesis testing framework for the Bayesian SITAR model
#'
#' @description
#' Performs non-linear hypothesis testing via [brms::hypothesis()], Region of
#' Practical Equivalence (\code{ROPE}) via [bayestestR::equivalence_test()],
#' probability of direction (\code{pd}) via [bayestestR::p_direction()], and
#' pairwise and contrast-style comparisons via [marginaleffects::comparisons()],
#' using a unified interface.
#'
#' The \code{ROPE} evaluates whether the parameter values should be accepted or
#' rejected against an explicitly formulated "null hypothesis". It checks the
#' percentage of the posterior credible intervals such as the High Density
#' Intervals (\code{HDI})  defined as the null region (the \code{ROPE}). If this
#' percentage is sufficiently low, the null hypothesis is rejected. If this
#' percentage is sufficiently high, the null hypothesis is accepted.  If the HDI
#' is completely outside the \code{ROPE}, the "null hypothesis" for this
#' parameter is "rejected". If the \code{ROPE} completely covers the HDI, i.e.,
#' all most credible values of a parameter are inside the region of practical
#' equivalence, the null hypothesis is accepted. Else, it is undecided whether
#' to accept or reject the null hypothesis. See [bayestestR::equivalence_test()]
#' for further details.
#' 
#' The \code{PD}, also known as the Maximum Probability of Effect (\code{MPE})
#' is the probability that a parameter (described by its posterior distribution)
#' is strictly positive or negative (whichever is the most probable). Although
#' differently expressed, this index is fairly similar (i.e., is strongly
#' correlated) to the frequentest p-value (see details). See
#' [bayestestR::p_direction()] for further details.
#' 
#' @details
#' This function collects arguments from multiple upstream packages:
#' \itemize{
#'   \item \code{hypothesis_test} arguments (e.g., \code{hypothesis},
#'   \code{class}, \code{scope}, \code{alpha}, \code{robust}, \code{seed}) map
#'   to [brms::hypothesis()] and, in particular,
#'   \code{scope} controls where hypothesis components referenced in
#'   \code{hypothesis} are looked up. Note that to avoid conflicts with the
#'   \code{hypothesis} argument from
#'   [marginaleffects::comparisons()],
#'   the \code{hypothesis} argument for the brms package is renamed as
#'   \code{hypothesis_str}. In other words, the \code{hypothesis_str} (see
#'   below) is used to setup the
#'   \code{hypothesis} argument when calling the
#'   [brms::hypothesis()].
#'   \item
#'   \code{equivalence_test} arguments (e.g., \code{range}, \code{ci},
#'   \code{remove_na}, \code{parameters}) map to
#'   [bayestestR::equivalence_test()]
#'   for ROPE-based practical equivalence decisions.
#'   \item \code{p_direction} arguments (e.g., \code{method}, \code{null},
#'   \code{ci}, \code{remove_na}, \code{parameters}) map to
#'   [bayestestR::p_direction()] for
#'   probability of direction based practical decisions.
#'   \item Additional arguments in \code{...} are forwarded to
#'   [marginaleffects::comparisons()]
#'   via \code{get_growthparameters} to compute contrasts between
#'   predictions under different predictor values.
#' }
#'
#' Users should supply either \code{hypothesis_str} (\pkg{brms} hypothesis),
#' \code{parameters} (\pkg{bayestestR}), or \code{parameter}
#' (growth-parameter selector used by this package), but not multiple
#' overlapping selectors.
#'
#' Supported inputs: \code{bgmfit} objects (via all three packages),
#' posterior data frames (via \pkg{brms} and \pkg{bayestestR} methods), or
#' \code{brmsfit} objects (via brms).
#'
#' @param model A fitted model object of class \code{bgmfit}, or a posterior
#'   data frame. See Details for engine-specific support.
#'   
#' @param by Optional grouping variable(s) used to stratify summaries or
#'   comparisons. The exact meaning depends on the selected engine and methods.
#'   
#' @param hypothesis_str A character vector specifying one or more non-linear
#'   hypothesis concerning parameters of the model. Note that the argument
#'   \code{hypothesis_str} is the renamed version of the \code{hypothesis}
#'   argument for the [brms::hypothesis()]. This renaming is done to avoid
#'   conflicts with the \code{hypothesis} argument for the
#'   [marginaleffects::comparisons()].
#'   
#' @param rope_test Logical; if \code{TRUE}, run ROPE-based practical
#'   equivalence testing via [bayestestR::equivalence_test()]. This is rarely
#'   used because appropriate setting \code{TRUE/FALSE} is assigned
#'   automatically depending on the \code{equivalence_test}. Note that when
#'   \code{equivalence_test} is used to set up the \code{rope_test}, then all
#'   missing arguments to the \code{equivalence_test} list are assigned from the
#'   default or user defined arguments included in the [hypothesis_test()] call.
#'   See \code{equivalence_test} for details.
#'   
#' @param pd_test Logical; if \code{TRUE}, compute probability of direction via
#'   [bayestestR::p_direction()]. This is rarely used because appropriate
#'   setting \code{TRUE/FALSE} is assigned automatically depending on the
#'   \code{p_direction}. Note that when \code{p_direction} is used to set up the
#'   \code{pd_test}, then all missing arguments to the \code{p_direction} list
#'   are assigned from the default or user defined arguments included in the
#'   [hypothesis_test()] call. See \code{p_direction} for details.
#' 
#' @param digits Number of digits to use when printing numeric results. Default
#'   \code{2}.
#' 
#' @param engine Optional engine selector specifying which backend(s) to use:
#'   \itemize{
#'     \item \code{"brms"} for [brms::hypothesis()] non-linear hypothesis
#'     testing
#'     \item \code{"bayestestR"} for [bayestestR::equivalence_test()] (ROPE) and
#'     [bayestestR::p_direction()] (PD) testing
#'     \item \code{"marginaleffects"} for [marginaleffects::comparisons()]
#'     pairwise/contrast comparisons
#'     \item \code{"mbcombo"} for combined \pkg{marginaleffects} and
#'     \pkg{bayestestR} workflows
#'     \item \code{NULL} (default) automatically selects appropriate engine(s)
#'     based on \code{model} class and requested analyses
#'   }
#'   
#' @param probs Probability levels for credible intervals. Defaults to
#'   \code{c(0.025, 0.975)} (95\% credible intervals). Also determines default
#'   values for \code{alpha}, \code{ci}, and \code{conf_level} when those are
#'   \code{NULL} (see below).
#'   
#' @param conf_level Confidence level used in [marginaleffects::comparisons()].
#'   The default is \code{NULL}, in which case \code{conf_level} is inferred
#'   from \code{probs}. With the default \code{probs}, this corresponds to
#'   \code{conf_level = 0.95}, which is also the default in
#'   [marginaleffects::comparisons()]. Note that if \code{probs = NULL}, then
#'   \code{conf_level} is used to determine the \code{probs} as well as the
#'   \code{alpha} and \code{ci} (see below).
#'   
#' @param alpha Alpha level for hypothesis tests used in [brms::hypothesis()].
#'   The default is \code{NULL}, in which case \code{alpha} is inferred from
#'   \code{probs}. With the default \code{probs}, this corresponds to
#'   \code{alpha = 0.05}, which is also the default in [brms::hypothesis()].
#'
#' @param ci Confidence level used in [bayestestR::equivalence_test()] and
#'   [bayestestR::p_direction()]. The default is \code{NULL}, in which case
#'   \code{ci} is inferred from \code{probs}. With the default \code{probs},
#'   this corresponds to \code{ci = 0.95}, which matches the defaults in those
#'   functions.
#'
#' @param evaluate_comparison Logical indicating whether to compute and return
#'   comparisons from [marginaleffects::comparisons()] via the
#'   [get_growthparameters()]. Defaults to \code{NULL}, in which case it is
#'   automatically determined from other arguments. Only
#'   evaluated when \code{model} is class \code{"bsitar"} and \code{engine =
#'   "marginaleffects" or "mbcombo"}.
#' 
#' @param comparison_by Character string specifying the \code{by} variable for
#'   comparisons. Defaults to \code{NULL}, in which case it is inferred from the
#'   \code{by} argument.
#' 
#' @param comparison_range_null Data frame with columns \code{parameter},
#'   \code{by} variables, \code{range} (for [bayestestR::equivalence_test()]),
#'   and \code{null} (for [bayestestR::p_direction()]). Defaults to \code{NULL},
#'   in which case \code{range} and \code{null} values are set automatically.
#'
#' @param comparison_range Like \code{comparison_range_null} but without the
#'   \code{null} column. Used when \code{equivalence_test = TRUE} and
#'   \code{p_direction = FALSE}.
#'
#' @param comparison_null Like \code{comparison_range_null} but without the
#'   \code{range} column. Used when \code{equivalence_test = FALSE} and
#'   \code{p_direction = TRUE}.
#'   
#' @param evaluate_hypothesis Logical indicating whether to compute and return
#'   hypothesis from [marginaleffects::comparisons()] via the
#'   [get_growthparameters()]. Defaults to \code{NULL}, in which case it is
#'   automatically determined from other arguments. Only
#'   evaluated when \code{model} is class \code{"bsitar"} and \code{engine =
#'   "marginaleffects" or "mbcombo"}.
#' 
#' @param hypothesis_by Character string specifying the \code{by} variable for
#'   hypotheses. Defaults to \code{NULL}, in which case it is inferred from the
#'   \code{by} argument.
#' 
#' @param hypothesis_range_null Data frame with columns \code{parameter},
#'   \code{by} variables, \code{range} (for [bayestestR::equivalence_test()]),
#'   and \code{null} (for [bayestestR::p_direction()]). Defaults to \code{NULL},
#'   in which case \code{range} and \code{null} values are set automatically.
#'
#' @param hypothesis_range Like \code{hypothesis_range_null} but without the
#'   \code{null} column. Used when \code{equivalence_test = TRUE} and
#'   \code{p_direction = FALSE}.
#'
#' @param hypothesis_null Like \code{hypothesis_range_null} but without the
#'   \code{range} column. Used when \code{equivalence_test = FALSE} and
#'   \code{p_direction = TRUE}.
#'   
#' 
#' @param rope_as_percent Logical indicating whether to return ROPE results as
#'   percentages (\code{TRUE}, default) or fractions (\code{FALSE}). Only
#'   evaluated when \code{model} is class \code{"bsitar"} and \code{engine} is
#'   \code{"bayestestR"} or \code{"mbcombo"}.
#' 
#' @param pd_as_percent Logical indicating whether to return PD results as
#'   percentages (\code{TRUE}, default) or fractions (\code{FALSE}). Only
#'   evaluated when \code{model} is class \code{"bsitar"} and \code{engine} is
#'   \code{"bayestestR"} or \code{"mbcombo"}.
#' 
#' @param format Logical; if \code{TRUE}, merge interval-bound columns returned
#'   by [bayestestR::equivalence_test()] into compact range columns.
#'   Specifically, \code{ROPE_low} and \code{ROPE_high} are combined into
#'   \code{ROPE_range}, and \code{HDI_low} and \code{HDI_high} are combined into
#'   \code{HDI_range}. This option only modifies results originating from
#'   [bayestestR::equivalence_test()].
#'
#' @param reformat Logical; if \code{TRUE}, post-process the data frame returned
#'   by [marginaleffects::comparisons()] to make it consistent with
#'   posterior-summary conventions. This includes renaming interval columns such
#'   as \code{conf.low} to \code{Q2.5} and \code{conf.high} to \code{Q97.5}, and
#'   optionally dropping columns which are often not needed for reporting (e.g.,
#'   \code{term}, \code{contrast}, etc.).
#'
#' @param get_range_null_form Logical; if \code{TRUE}, return the expected
#'   structure of \code{range} ([bayestestR::equivalence_test()]) and
#'   \code{null} ([bayestestR::p_direction()]) values. User can use this format
#'   to set the \code{range} and \code{null} values that are then passed to the
#'   [bayestestR::equivalence_test()] and [bayestestR::equivalence_test()] via
#'   \code{comparison_range_null} and \code{hypothesis_range_null}.
#'   
#' @param get_range_null_value Logical; if \code{TRUE}, return the default
#'   \code{range} ([bayestestR::equivalence_test()]) and \code{null}
#'   ([bayestestR::p_direction()]) values that are used to set up the
#'   \code{range} and \code{null} values for [bayestestR::equivalence_test()]
#'   and [bayestestR::equivalence_test()].
#'   
#' @param method_call A character string to pass \code{'method'} argument to 
#' the \code{'get_growthparameters'} function. Default \code{'pkg'}. 
#' 
#' @param na.rm Logical; if \code{TRUE} (default), then remove \code{NA} values.
#'  
#' @param plot Logical or a character string \code{"return"}; if \code{TRUE},
#'   plot the hypothesis test results and return the data frame. Set \code{plot
#'   = "return"} to return the plot object itself instead of the data frame.
#'   Note that only [brms::hypothesis()] and [bayestestR::equivalence_test()]
#'   support plotting.
#' 
#' @param ... Additional arguments forwarded to [get_growthparameters()]
#' 
#' 
#' @inheritParams marginaleffects::comparisons
#' @inheritParams brms::hypothesis.brmsfit
#' @inheritParams marginaleffects::hypotheses
#' @inheritParams bayestestR::equivalence_test
#' @inheritParams bayestestR::p_direction
#' @inheritParams get_growthparameters.bgmfit
#' 
#' @inherit brms::hypothesis.brmsfit return description details 
#' seealso sections references note format
#'
#' @return
#' An object (typically a list) containing some or all of the following
#' components, depending on the requested analyses and available methods:
#' \itemize{
#'   \item A [brms::hypothesis()] result object for hypothesis tests.
#'   \item A [bayestestR::equivalence_test()] result for ROPE and practical
#'   equivalence.
#'   \item A [bayestestR::p_direction()] result for probability of direction.
#'   \item A [marginaleffects::comparisons()] result for contrasts between
#'   predictions.
#' }
#'
#' @seealso
#' \itemize{
#'   \item [brms::hypothesis()] for hypothesis testing.
#'   \item [bayestestR::equivalence_test()] for ROPE tests.
#'   \item [bayestestR::p_direction()] for probability of direction tests.
#'   \item [marginaleffects::comparisons()] for prediction contrasts.
#' }
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether the model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Speed up examples by using a subset of posterior draws
#' draw_ids <- 1:5
#' 
#' ## Example 1: Hypothesis testing (brms-style)
#'
#' set_hypothesis <- c("a_Intercept = 0", "b_Intercept > 1")
#'
#' # brms reference
#' hyp_ex1 <- brms::hypothesis(model, hypothesis = set_hypothesis, draw_ids = draw_ids)
#'
#' # hypothesis_test (auto-detects brms engine)
#' hyp_ex11 <- hypothesis_test(model, 
#'                            draw_ids = draw_ids,
#'                            hypothesis_str = set_hypothesis)
#'
#' # If see package is installed (install.packages("see")), then you can plot it 
#' # plot(hyp_ex1)
#' # plot(hyp_ex11)
#' 
#' ## Example 2: ROPE equivalence testing (bayestestR)
#'
#' set_parameters <- c("b_a_Intercept", "b_b_Intercept")
#' set_range <- list(b_a_Intercept = c(100, 150), b_b_Intercept = c(-2, 2))
#'
#' # bayestestR reference - If you installed bayestestR package
#' # hyp_ex1 <- bayestestR::equivalence_test(model, 
#' #                                      parameters = set_parameters, 
#' #                                      range = set_range,
#' #                                      draw_ids = draw_ids)
#'
#' # hypothesis_test (auto-detects bayestestR engine)
#' # If see package is installed (install.packages("see")), then you can plot 
#' # the results by setting plot = TRUE
#' hyp_ex11 <- hypothesis_test(model,
#'                            draw_ids = draw_ids,
#'                            parameters = set_parameters,
#'                            range = set_range,
#'                            plot = FALSE)
#'
#' # print(hyp_ex1)
#' print(hyp_ex11)
#' 
#' # Note: bayestestR often ignores parameters/effects for nonlinear bsitar models.
#' # hypothesis_test correctly respects user-specified parameters.
#' 
#' ## Example 3: p-direction testing (bayestestR)
#'
#' set_parameters <- c("b_a_Intercept", "b_b_Intercept")
#' set_null <- list(b_a_Intercept = 0, b_b_Intercept = 0)
#'
#' # bayestestR reference - If you installed bayestestR package
#' # hyp_ex1 <- bayestestR::p_direction(model, 
#' #                                 parameters = set_parameters,
#' #                                 null = set_null,
#' #                                 draw_ids = draw_ids)
#'
#' # hypothesis_test (auto-detects bayestestR engine)
#' # If see package is installed (install.packages("see")), then you can plot 
#' # the results by setting plot = TRUE
#' hyp_ex11 <- hypothesis_test(model,
#'                            draw_ids = draw_ids,
#'                            parameters = set_parameters,
#'                            null = set_null,
#'                            plot = FALSE,
#'                            equivalence_test = FALSE,
#'                            p_direction = TRUE)
#'
#' # print(hyp_ex1)
#' print(hyp_ex11)
#' 
#' ## Example 4: Growth parameter hypothesis testing (marginaleffects)
#'
#' set_parameter <- c("apgv", "pgv")
#' set_range <- list(apgv = 1, pgv = c(0, 0.5))
#' set_null <- list(apgv = 1, pgv = 0)
#'
#' # Directly using get_growthparameters()
#' hyp_ex1 <- get_growthparameters(model, 
#'                                     parameter = set_parameter,
#'                                     draw_ids = draw_ids)
#' # Uncomment for ROPE/p-direction tests:
#' # equivalence_test = list(range = set_range),
#' # p_direction = list(null = set_null)
#' 
#' ## Example 5: hypothesis_test() with growth parameters
#'
#' hyp_ex1 <- hypothesis_test(model,
#'                           parameter = set_parameter,
#'                           equivalence_test = list(range = set_range),
#'                           p_direction = list(null = set_null),
#'                           draw_ids = draw_ids)
#' 
#' ## Example 6: Pass posterior draws from get_growthparameters()
#'
#' # pdrawsp = TRUE returns draws for downstream hypothesis testing
#' hyp_draws <- get_growthparameters(model,
#'                                       parameter = set_parameter,
#'                                       pdrawsp = TRUE,
#'                                       equivalence_test = list(range = set_range),
#'                                       p_direction = list(null = set_null),
#'                                       draw_ids = draw_ids)
#'
#' # Use draws directly in hypothesis_test()
#' hyp_ex1 <- hypothesis_test(hyp_draws,
#'                           parameter = set_parameter,
#'                           equivalence_test = list(range = set_range),
#'                           p_direction = list(null = set_null),
#'                           draw_ids = draw_ids)
#' 
#' }
#'
#' @rdname hypothesis_test
#' @export
#'
#' @inherit berkeley author
#' 
hypothesis_test.bgmfit <- function(model,
                                   by = NULL,
                                   hypothesis = NULL,
                                   hypothesis_str = NULL,
                                   parameters = NULL,
                                   parameter = NULL, 
                                   class = "b",
                                   group = "",
                                   scope = c("standard", "ranef", "coef"),
                                   robust = FALSE,
                                   seed = NULL,
                                   equivalence_test = NULL,
                                   p_direction = NULL,
                                   rope_test = NULL,
                                   pd_test = NULL,
                                   range = "default",
                                   method = "direct",
                                   method_call = "pkg",
                                   null = NULL,
                                   as_p = FALSE,
                                   remove_na = TRUE,
                                   rvar_col = NULL,
                                   effects = "fixed",
                                   component = "conditional",
                                   evaluate_comparison = NULL,
                                   comparison_by = NULL,
                                   comparison_range_null = NULL,
                                   comparison_range = NULL,
                                   comparison_null = NULL,
                                   evaluate_hypothesis = NULL,
                                   hypothesis_by = NULL,
                                   hypothesis_range_null = NULL,
                                   hypothesis_range = NULL,
                                   hypothesis_null = NULL,
                                   rope_as_percent = TRUE,
                                   pd_as_percent = TRUE,
                                   format = NULL,
                                   reformat = TRUE,
                                   estimate_center = NULL,
                                   estimate_interval = NULL,
                                   na.rm = TRUE,
                                   alpha = 0.05,
                                   ci = 0.95,
                                   conf_level = NULL,
                                   probs = c(0.025, 0.975),
                                   get_range_null_form = NULL,
                                   get_range_null_value = NULL,
                                   plot = FALSE,
                                   digits = 2,
                                   engine = NULL,
                                   usesavedfuns = FALSE,
                                   verbose = FALSE,
                                   envir = NULL,
                                   ...) {
  
  Parameter <- NULL;
  Effects <- NULL;
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- envir
  }
  
  if(!is.null(by)) {
    if(!is.character(by)) {
      stop("Argument 'by' must be a character string or a character vector")
    }
  }
  
  if(is.null(probs) & is.null(conf_level)) {
    stop2c("Please specify either 'probs' or 'conf_level'")
  } else if(!is.null(conf_level)) {
    conf <- conf_level
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  } else if(!is.null(probs)) {
    conf <- ci <- conf_level <- probs[2] - probs[1] 
    probs <- probs
  }
  ci <- conf
  alpha <- 1-conf
  
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', probtitles)
  
  
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
  
  ##########################################
  
  allowed_engine <- c('brms', 'marginaleffects', 'bayestestR', 'mbcombo')
  if(!is.null(engine)) {
    if(grepl("^marg", engine)) {
      engine <- 'marginaleffects'
    } else if(grepl("^brm", engine)) {
      engine <- 'brms'
    } else if(grepl("^bayestestR", engine)) {
      engine <- 'bayestestR'
    } else if(grepl("^bayes", engine)) {
      engine <- 'bayestestR'
    } else if(grepl("^mb", engine)) {
      engine <- 'mbcombo'
    } else if(!engine %in% allowed_engine) {
      stop("Argument engine must be one of the following: ",
           collapse_comma(allowed_engine))
    } # else if(!engine %in% allowed_engine) {
  } 
  
  obj_model <- obj_df <- obj_mfx <- obj_mfx_matrix <- obj_model_pseudo <- FALSE
  obj_marginaleffects <- FALSE
  if(inherits(model, 'bgmfit')) {
    if(!inherits(model, 'data.frame') & !inherits(model, 'data.table')) {
      obj_model <- TRUE
    } else if(inherits(model, 'data.frame') | inherits(model, 'data.table')) {
      if(inherits(model, 'mfx_draws') | inherits(model, 'mfx')) {
        obj_mfx <- TRUE
      } else {
        obj_df <- TRUE
      }
      model_attr_orig <- attr(model, 'class') 
      xcall <- 'hypothesis_test.bgmfit'
      obj_model_pseudo <- TRUE
    }
  } else if(!inherits(model, 'bgmfit')) {
    if(inherits(model, 'data.frame') | inherits(model, 'data.table')) {
      # New layes added 
      # key_names_mfx <- c("drawid",  "parameter", "draw")
      # if(all(key_names_mfx %in% attr(model, 'names'))) {
      #   obj_mfx <- TRUE
      # } else 
      key_names_mfx <- c("drawid",  "parameter", "draw")
      if(all(key_names_mfx %in% attr(model, 'names'))) {
        obj_mfx <- TRUE
      } else if(inherits(model, 'mfx_draws') | inherits(model, 'mfx')) {
        obj_mfx <- TRUE
      } else {
        obj_df <- TRUE
      }
      model_attr_orig <- attr(model, 'class') 
      # attr(model, 'class') <- c(model_attr_orig, 'bgmfit')
      xcall <- 'hypothesis_test.bgmfit'
      obj_model_pseudo <- TRUE
    } else if(inherits(model, 'matrix')) {
      # get_draws(shape = "DxP")
      obj_mfx_matrix <- TRUE
    }
  } else if(!inherits(model, 'bgmfit')) {
    if(inherits(model, 'comparisons') |
       inherits(model, 'predictions') |
       inherits(model, 'slopes') |
       inherits(model, 'marginaleffects')) {
      obj_marginaleffects <- TRUE
    }
  }
 
  ##############################################################################
  if(!is.null(parameter)) {
    allowed_parms      <- c('apgv', 'pgv', 'atgv', 'tgv', 'acgv', 'cgv')
    allowed_parms_size <- c('spgv', 'stgv', 'scgv')
    check_set_parm_out <- check_set_parm(parameter = parameter,
                                         allowed_parms = allowed_parms,
                                         allowed_parms_size = allowed_parms_size,
                                         default_parms = c('apgv'),
                                         setpreparms = FALSE,
                                         plot = plot,
                                         verbose = FALSE)
    
    eout_check_set_parm_out <- list2env(check_set_parm_out)
    for (eoutii in names(eout_check_set_parm_out)) {
      if(!is.null(eout_check_set_parm_out[[eoutii]])) {
        assign(eoutii, eout_check_set_parm_out[[eoutii]])
      }
    }
    if(!exists('parm'))               parm <- NULL
    if(!exists('sat_ptc'))            sat_ptc <- NULL
    if(!exists('parameter'))          parameter <- NULL
    if(!exists('parameter_arg'))      parameter_arg <- NULL
    if(!exists('parameter_sat'))      parameter_sat <- NULL
    if(!exists('string_sat'))         string_sat <- NULL
    if(!exists('numeric_sat'))        numeric_sat <- NULL
    if(!exists('string_numeric_sat')) string_numeric_sat <- NULL
    # For get_comparison_hypothesis
    parms_sat_elements <- list()
    parms_sat_elements[['sat_ptc']]            <- sat_ptc
    parms_sat_elements[['parameter_sat']]      <- parameter_sat
    parms_sat_elements[['string_sat']]         <- string_sat
    parms_sat_elements[['numeric_sat']]        <- numeric_sat
    parms_sat_elements[['string_numeric_sat']] <- string_numeric_sat
  } else if(is.null(parameter)) {
    parms_sat_elements <- NULL
  }
  
  
  ##############################################################################
  
  # Override if user sets range null not as equivalence_test and p_direction
  # lists but as open arguments
  if(!is.null(parameter)) {
    if(is.null(equivalence_test)) {
      if(!is.null(range))  {
        if(is.character(range)) {
          # if(range != "default") equivalence_test <- list(range = range)
        } else {
          equivalence_test <- list(range = range)
        }
      }
    }
    if(is.null(p_direction)) {
      if(!is.null(null))  p_direction <- list(null = null)
    }
  }
  
  # Note:
  # The missing arguments for equivalence_test will be later updated from the
  # hypothesis_test() arguments via bayestestR_equivalence_test_df_args - 1301
  # The missing arguments for p_direction will be later updated from the
  # hypothesis_test() arguments via bayestestR_p_direction_df_args - 1301
  if(is.null(rope_test)) {
    if(!is.null(equivalence_test)) {
      if(is.logical(equivalence_test)) {
        rope_test <- equivalence_test
      } else if(is.list(equivalence_test)) {
        rope_test <- TRUE
      }
    }
  }
  
  if(is.null(pd_test)) {
    if(!is.null(p_direction)) {
      if(is.logical(p_direction)) {
        pd_test <- p_direction
      } else if(is.list(p_direction)) {
        pd_test <- TRUE
      }
    }
  }
  
  if(is.null(rope_test) & is.null(pd_test)) {
    if(obj_model) {
      if(!is.null(parameters)) {
        rope_test  <- TRUE
        pd_test <- FALSE
        if(verbose) {
          stop2c("For the model object, argument 'parameters' is specified but  
                 both 'rope_test' and 'pd_test' are NULL. The  
                 expected engine appears to be 'bayestestR', so the default for 
                 'rope_test' changes from NULL to TRUE")
        }
      }
    }
  }
  
  
  if(is.null(rope_test)) {
    rope_test <- FALSE 
  }
  if(is.null(pd_test)) {
    pd_test <- FALSE 
  }
  

  if(pd_test) {
    if(is.null(p_direction)) {
      if(is.null(null)) {
        if(verbose) {
          message2c("For 'pd_test = TRUE', the default 'null' for 'pd test' 
                  changed to '0'")
        }
        p_direction <- list(null = 0)
      } else if(!is.null(null)) {
        p_direction <- list(null = null)
      }
    }
  }
  
  if(rope_test) {
    if(is.null(equivalence_test)) {
      if(is.null(range)) {
        if(verbose) {
          message2c("For 'rope_test = TRUE', the default 'range' for 'rope test' 
                  changed to 'default'")
        }
        equivalence_test <- list(range = 'default')
      } else if(!is.null(range)) {
        equivalence_test <- list(range = range)
      }
    }
  }
  
  # Extract draws from the obj_marginaleffects
  if(obj_marginaleffects) {
    if(verbose) {
      message2c("draws extrcated from the obj_marginaleffects")
    }
    allowed_engine_obj_marginaleffects <-c('brms', 'marginaleffects', 'mbcombo')
    if(!is.null(engine)) {
      if(!engine %in% allowed_engine_obj_marginaleffects) {
        stop("Argument engine must be one of the following: ",
             collapse_comma(allowed_engine))
      } # if(!engine %in% allowed_engine) {
      if(engine == 'brms') {
        model <- marginaleffects::get_draws(model, shape = "DxP")
        obj_mfx_matrix <- TRUE
        obj_model <- obj_df <- obj_mfx <- obj_model_pseudo <- FALSE
        obj_marginaleffects <- FALSE
      } else if(engine == 'marginaleffects' | engine == 'mbcombo') {
        model <- marginaleffects::get_draws(model, shape = "long")
        obj_mfx <- TRUE
        obj_model <- obj_df <- obj_mfx_matrix <- obj_model_pseudo <- FALSE
        obj_marginaleffects <- FALSE
      } # if(engine == 'brms') { else if(engine == 'marginaleffects' |
    } else if(is.null(engine)) {
      model <- marginaleffects::get_draws(model, shape = "long")
      obj_mfx <- TRUE
      obj_model <- obj_df <- obj_mfx_matrix <- obj_model_pseudo <- FALSE
      obj_marginaleffects <- FALSE
      if(!rope_test & !pd_test) {
        engine <- 'marginaleffects' 
      } else if(rope_test | pd_test) {
        engine <- 'mbcombo' 
      }
    }
  } # if(obj_marginaleffects) {
  
  check_messages_1 <- 
    paste0(
      "Please specify either 'hypothesis_str', 'parameters', or 'parameter' 
      - not all three or any combination of them. ",
      
      "The argument 'hypothesis_str' (used by brms) refers to the hypothesis
      components (See ?brms::hypothesis_str) ",
      
      " whereas 'parameters' (used by bayestestR) is synonymous with  
      hypothesis_str. (see ?bayestestR::equivalence_test) ",
      
      "The 'parameter', on the other hand, refers to the specific growth
       parameter such as apgv and pgv  (see ?get_growthparameters for
      details)")
  
  # Set default
  # This when everything NULL, and default setting up the hypothesis
  if(is.null(parameters) & is.null(parameter) & is.null(hypothesis_str)) {
    if(obj_model) {
      parameter <- c('apgv', 'pgv')
      engine <- 'marginaleffects'
      if(is.null(hypothesis)) {
        # if(!is.null(by)) hypothesis <- ~ pairwise
        if(!is.null(by)) hypothesis <- as.formula(~ pairwise)
      }
      if(is.null(rope_test)) rope_test <- FALSE
      if(is.null(pd_test)) pd_test <- FALSE
    }
  }
  
  if(!is.null(parameters) & !is.null(parameter) & !is.null(hypothesis_str)) {
    if(is.null(engine)) {
      stop2c(check_messages_1)
    } else if(!is.null(engine)) {
      if(engine == 'brms') {
        parameter <- NULL
        parameters <- NULL
        hypothesis <- NULL
      } else if(engine == 'bayestestR') {
        hypothesis_str <- NULL
        parameter <- NULL
      } else if(engine == 'marginaleffects') {
        hypothesis_str <- NULL
        parameters <- NULL
      } else if(engine == 'mbcombo') {
        hypothesis_str <- NULL
        parameters <- NULL
      }
    } # if(is.null(engine)) { else if(!is.null(engine)) {
  } else if(!is.null(parameters) & !is.null(parameter) &
            is.null(hypothesis_str)) {
    if(is.null(engine)) {
      stop2c(check_messages_1)
    } else if(!is.null(engine)) {
      if(engine == 'brms') {
        if(is.null(hypothesis_str)) stop2c("For engine = ",
                                           engine, 
                                           " 'hypothesis_str' is required")
        parameter <- NULL
        parameters <- NULL
      } 
    } # if(is.null(engine)) { else if(!is.null(engine)) {
  } else if(!is.null(parameters) & is.null(parameter) & 
            !is.null(hypothesis_str)) {
    if(is.null(engine)) {
      stop2c(check_messages_1)
    } else if(!is.null(engine)) {
      if(engine == 'marginaleffects' | engine == 'mbcombo') {
        if(is.null(parameter)) stop2c("For engine = ", 
                                      engine, " 'parameter' is required")
        hypothesis_str <- NULL
        parameters <- NULL
      } 
    } # if(is.null(engine)) { else if(!is.null(engine)) {
  } else if(is.null(parameters) & !is.null(parameter) & 
            !is.null(hypothesis_str)) {
    if(is.null(engine)) {
      stop2c(check_messages_1)
    } else if(!is.null(engine)) {
      if(engine == 'bayestestR') {
        if(obj_model) {
          if(is.null(parameters)) {
            stop2c("For engine = ", engine, " and model object, the argument
                                          'parameters' is required")
          }
        }
        hypothesis_str <- NULL
        parameter <- NULL
      } 
    } # if(is.null(engine)) { else if(!is.null(engine)) {
  }
  
  if(!is.null(engine)) {
    if(engine == 'brms' | 
       engine == 'bayestestR' | 
       engine == 'marginaleffects' |
       engine == 'mbcombo') {
      if(!inherits(model, 'bgmfit')) {
        stop2c("You have set engine = ", collapse_comma(engine)," but the object 
           is not of class 'bsitar'. The class of object is ",
               collapse_comma(attr(model, 'class')))
      }
    } 
  }
  
  get_engine_obj_args_names  <- expression(engine, 
                                           hypothesis_str,
                                           parameters, 
                                           parameter,
                                           hypothesis,
                                           obj_model,
                                           obj_df,
                                           obj_mfx,
                                           obj_mfx_matrix,
                                           obj_model_pseudo,
                                           obj_marginaleffects
  )
  get_engine_obj_args_names <- as.character(get_engine_obj_args_names)
  
  get_engine_obj_args_all_c      <- c()
  get_engine_obj_args_internal_c <- c()
  get_engine_obj_args_print_c    <- c()
  for (i in get_engine_obj_args_names) {
    zzz <- deparse(get(i))
    zzz <- gsub("\"", "", zzz)
    zzz <- paste0(i, " = ", zzz)
    get_engine_obj_args_all_c <- c(get_engine_obj_args_all_c, zzz)
    if(grepl("^obj_", i)) {
      get_engine_obj_args_internal_c <- c(get_engine_obj_args_internal_c, zzz)
    } else if(!grepl("^obj_", i)) {
      get_engine_obj_args_print_c <- c(get_engine_obj_args_print_c, zzz)
    }
  } # for (i in get_engine_obj_args) {
  
  get_engine_obj_args_all <- paste(get_engine_obj_args_all_c,
                                   collapse = "; ")
  get_engine_obj_args_internal <- paste(get_engine_obj_args_internal_c,
                                        collapse = "; ")
  get_engine_obj_args_print <- paste(get_engine_obj_args_print_c,
                                     collapse = "; ")
  
  
  if(!is.logical(verbose)) {
    if(verbose == 1) {
      print('Arguments to build objects')
      print(get_engine_obj_args_print)
    }
  }
  
  
  # obj_model -> Check and set up engine and hypothesis_str/parameters/parameter
  obj_model_brms <- obj_model_bayestestR <- FALSE
  obj_model_marginaleffects <- obj_model_mbcombo <- FALSE
  if(obj_model) {
    allowed_engine_obj_model <- c('brms', 'marginaleffects', 
                                  'bayestestR', 'mbcombo')
    if(is.null(engine)) {
      if(!is.null(hypothesis_str) & !is.null(parameters) & 
         !is.null(parameter)) {
        stop2c("For model object, specify either 'hypothesis_str' or 
               'parameters' or 'parameter, not more than one of them")
      } else if(is.null(hypothesis_str) & is.null(parameters) & 
                is.null(parameter)) {
        stop2c("For model object, one of the 'hypothesis_str' or 
               'parameters' or 'parameter must be specified")
      } else if(!is.null(hypothesis_str)) {
        engine <- 'brms'
        obj_model_brms <- TRUE
      } else if(!is.null(parameters)) {
        engine <- 'bayestestR'
        obj_model_bayestestR <- TRUE
      } else if(!is.null(parameter)) {
        if(!rope_test & !pd_test) {
          engine <- 'marginaleffects'
          obj_model_marginaleffects <- TRUE
        } else if(rope_test | pd_test) {
          engine <- 'mbcombo'
          obj_model_mbcombo <- TRUE
        }
      } 
    } else if(!is.null(engine)) {
      if(!engine %in% allowed_engine_obj_model) {
        stop("Argument engine must be one of the following: ",
             collapse_comma(allowed_engine))
      } # if(!engine %in% allowed_engine) {
      if(engine == 'brms') {
        if(is.null(hypothesis_str)) {
          stop2c("For engine = 'brms', 'hypothesis_str' must be specified")
        }
        parameter <- NULL
        parameters <- NULL
        obj_model_brms <- TRUE
      } else if(engine == 'bayestestR') {
        if(is.null(parameters)) {
          stop2c("For engine = ", engine, " and model object, the argument
                                          'parameters' is required")
        }
        parameter <- NULL
        hypothesis_str <- NULL
        obj_model_bayestestR <- TRUE
      } else if(engine == 'marginaleffects' | engine == 'mbcombo') {
        if(is.null(parameter)) {
          stop2c("For engine = 'marginaleffects' and engine = 'mbcombo' with 
               model object, 'parameter' must be specified")
        }
        parameters <- NULL
        hypothesis_str <- NULL
        if(!rope_test & !pd_test) {
          obj_model_marginaleffects <- TRUE
        } else if(rope_test | pd_test) {
          obj_model_mbcombo <- TRUE
        }
      } # else if(engine == 'marginaleffects' | engine == 'mbcombo') {
    } # if(is.null(engine)) { else if(!is.null(engine)) {
  } # if(obj_df) {
  
 
  
  # obj_df & obj_mfx_matrix
  # -> Check and set up engine and hypothesis_str / parameters / parameter
  obj_df_brms <- obj_df_bayestestR <- FALSE
  if(obj_df | obj_mfx_matrix) {
    allowed_engine_obj_df <- c('brms', 'bayestestR')
    if(is.null(engine)) {
      parameter <- NULL
      if(!is.null(hypothesis_str) & !is.null(parameters)) {
        stop2c("For data frame object, specify either 'hypothesis_str' or 
               'parameters', not both")
      } else if(is.null(hypothesis_str) & is.null(parameters)) {
        stop2c("For data frame object, either 'hypothesis_str' or 
               'parameters' must be specified")
      } else if(!is.null(hypothesis_str)) {
        engine <- 'brms'
        obj_df_brms <- TRUE
      } else if(!is.null(parameters)) {
        engine <- 'bayestestR'
        obj_df_bayestestR <- TRUE
      } 
    } else if(!is.null(engine)) {
      if(!engine %in% allowed_engine_obj_df) {
        stop("Argument engine must be one of the following: ",
             collapse_comma(allowed_engine))
      } # if(!engine %in% allowed_engine) {
      if(engine == 'brms') {
        if(is.null(hypothesis_str)) {
          stop2c("For engine = 'brms', 'hypothesis_str' must be specified")
        }
        parameter <- NULL
        parameters <- NULL
        obj_df_brms <- TRUE
      } else if(engine == 'bayestestR') {
        if(is.null(hypothesis_str)) {
          stop2c("For engine = 'bayestestR' and data frame object, 'parameters'  
               must be specified")
        }
        parameter <- NULL
        hypothesis_str <- NULL
        obj_df_bayestestR <- TRUE
      }
    } # if(is.null(engine)) { else if(!is.null(engine)) {
  } # if(obj_df) {
  
  
  
  # obj_mfx -> Check and set up engine and hypothesis_str/parameters/parameter
  obj_mfx_bayestestR <- FALSE
  obj_mfx_marginaleffects <- obj_mfx_mbcombo <- FALSE
  if(obj_mfx) {
    allowed_engine_obj_mfx <- c('marginaleffects', 'mbcombo')
    if(is.null(engine)) {
      if(!is.null(hypothesis_str) & !is.null(parameters)) {
        stop2c("For mfx object, 'hypothesis_str' or  
               'parameters' are not allowed")
      } else if(is.null(parameter)) {
        stop2c("For mfx object, 'parameter must be specified")
      } else if(!is.null(parameter)) {
        if(!rope_test & !pd_test) {
          engine <- 'marginaleffects'
          obj_mfx_marginaleffects <- TRUE
        } else if(rope_test | pd_test) {
          engine <- 'mbcombo'
          obj_mfx_mbcombo <- TRUE
        }
      } 
    } else if(!is.null(engine)) {
      if(!engine %in% allowed_engine_obj_mfx) {
        stop("Argument engine must be one of the following: ",
             collapse_comma(allowed_engine))
      } # if(!engine %in% allowed_engine) {
      if(engine == 'marginaleffects' | engine == 'mbcombo') {
        if(is.null(parameter)) {
          stop2c("For engine = 'marginaleffects' and engine = 'mbcombo' with 
               mfx object, 'parameter' must be specified")
        }
        parameters <- NULL
        hypothesis_str <- NULL
        if(!rope_test & !pd_test) {
          obj_mfx_marginaleffects <- TRUE
        } else if(rope_test | pd_test) {
          obj_mfx_mbcombo <- TRUE
        }
      } # else if(engine == 'marginaleffects' | engine == 'mbcombo') {
    } # if(is.null(engine)) { else if(!is.null(engine)) {
  } # if(obj_mfx) {
  
  # Available choices and executions
  set_engine_obj_args_names  <- expression(obj_model_brms,
                                           obj_model_bayestestR,
                                           obj_model_marginaleffects,
                                           obj_model_mbcombo,
                                           obj_df_brms,
                                           obj_df_bayestestR,
                                           obj_mfx_bayestestR,
                                           obj_mfx_marginaleffects,
                                           obj_mfx_mbcombo)
  
  set_engine_obj_args_names <- as.character(set_engine_obj_args_names)
  
  set_engine_obj_args_all_c      <- c()
  set_engine_obj_args_internal_c <- c()
  set_engine_obj_args_print_c    <- c()
  for (i in set_engine_obj_args_names) {
    zzz <- deparse(get(i))
    zzz <- gsub("\"", "", zzz)
    zzz <- paste0(i, " = ", zzz)
    set_engine_obj_args_all_c <- c(set_engine_obj_args_all_c, zzz)
    set_engine_obj_args_print_c <- c(set_engine_obj_args_print_c, zzz)
  }
  
  set_engine_obj_args_all <- paste(set_engine_obj_args_all_c,
                                   collapse = "; ")
  set_engine_obj_args_internal <- paste(set_engine_obj_args_internal_c,
                                        collapse = "; ")
  set_engine_obj_args_print <- paste(set_engine_obj_args_print_c,
                                     collapse = "; ")
  
  
  if(!is.logical(verbose)) {
    if(verbose == 1) {
      print('Objects created')
      print(set_engine_obj_args_print)
    }
  }
  
  
  ##############################################################################
  # Get and collect arguments
  ##############################################################################
  
  if(obj_model | obj_model_pseudo ) {
    if(!is.null(model$xcall)) {
      if(grepl("get_growthparameters", model$xcall)) {
        xcall <- "get_growthparameters"
      } else if(grepl("modelbased_growthparameters", model$xcall)) {
        xcall <- "modelbased_growthparameters"
      } 
    } else {
      rlang_trace_back <- rlang::trace_back()
      check_trace_back.bgmfit <- grepl(".bgmfit", rlang_trace_back[[1]])
      if(all(!check_trace_back.bgmfit)) {
        # nothing
      } else {
        rlang_trace_back.bgmfit_i <- min(which(check_trace_back.bgmfit == TRUE))
        rlang_trace_back.bgmfit <- 
          rlang_trace_back[[1]][[rlang_trace_back.bgmfit_i]]
        rlang_call_name <- rlang::call_name(rlang_trace_back.bgmfit)
        xcall <- rlang_call_name
      }
    }
    
    # Imp - when model_marginaleffects, xcall must be 'get_growthparameters
    if(obj_model_marginaleffects | obj_model_mbcombo ) {
      xcall <- 'get_growthparameters'
    }
    
    
    check_if_package_installed(model, xcall = xcall)
    # model$xcall <- xcall
    arguments <- get_args_(as.list(match.call())[-1], xcall)
    arguments$model <- model
    arguments$usesavedfuns <- usesavedfuns
    # CustomDoCall 
    arguments <- sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                            arguments = arguments, 
                                            check_formalArgs = NULL,
                                            check_trace_back = NULL,
                                            envir = parent.frame())
    
    get.cores_ <- get.cores(arguments$cores)
    # 28.09.2024
    if(is.null(get.cores_[['max.cores']])) {
      if(is.null(arguments$cores)) 
        get.cores_[['max.cores']] <- future::availableCores() - 1
    }
    arguments$cores <- setincores <-  get.cores_[['max.cores']]
    .cores_ps <- get.cores_[['.cores_ps']]
    
    arguments[['ci']] <- ci
    arguments[['conf_level']] <- conf_level
    arguments[['alpha']] <- alpha
    
    arguments[['estimate_center']] <- estimate_center
    arguments[['estimate_interval']] <- estimate_interval
    arguments[['ec_agg']] <- ec_agg
    arguments[['ei_agg']] <- ei_agg
    arguments[['nthreads']] <- arguments$cores 
    
    if(is.null(arguments$pd_test)) {
      arguments$pd_test <- pd_test
    }
    if(is.null(arguments$rope_test)) {
      arguments$rope_test <- rope_test
    }
    
    
    if(is.null(arguments$p_direction)) {
      arguments$p_direction <- p_direction
    }
    if(is.null(arguments$equivalence_test)) {
      arguments$equivalence_test <- equivalence_test
    }
    
    if(NullFALSE(arguments$p_direction)) {
      if(is.list(arguments$p_direction)) {
        arguments$null <- arguments$p_direction$null
      }
    }
    
    if(NullFALSE(arguments$equivalence_test)) {
      if(is.list(arguments$equivalence_test)) {
        arguments$range <- arguments$equivalence_test$range
      }
    }
    
    arguments$parms_sat_elements <- parms_sat_elements
    
    if(obj_model_pseudo) {
      attr(model, 'class') <- model_attr_orig
    }
    
  } # if(obj_model | obj_model_pseudo) {
  
  
  
  ##############################################################################
  # Build arguments
  ##############################################################################
  
  full.arguments <-  arguments
  
  brms_hypothesis_args_df_names <- 
    c('x', 'hypothesis_str', 'class', 'group', 'scope', 'alpha',
      'robust', 'seed')
  
  brms_hypothesis_args_model_names <- 
    c('x', 'hypothesis_str', 'alpha', 'robust')
  
  bayestestR_equivalence_test_df_names <- 
    c('x', 'range', 'ci', 'rvar_col', 'verbose')
  
  bayestestR_equivalence_test_model_names <- 
    c('x', 'range', 'ci', 'effects', 'component', 'parameters', 'verbose')
  
  bayestestR_p_direction_df_names <- 
    c('x', 'method', 'null', 'as_p', 'remove_na', 'rvar_col')
  
  bayestestR_p_direction_model_names <- 
    c('x', 'method', 'null', 'effects', 'component', 'parameters', 
      'as_p', 'remove_na')
  
  bayestestR_p_direction_model_get_predicted_names <- 
    c('x', 'method', 'null', 'as_p', 'remove_na', 'use_iterations', 'verbose')
  
  
  eqpd_df<-c('format', 'get_form', 'get_value', 'digits', 
             'as_percent', 'inline')
  
  bayestestR_equivalence_test_df_names <- 
    c(bayestestR_equivalence_test_df_names, eqpd_df)
  bayestestR_p_direction_df_names <- 
    c(bayestestR_p_direction_df_names, eqpd_df)
  
  bayestestR_equivalence_test_df_names_ALL <- 
    c(brms_hypothesis_args_df_names,
      brms_hypothesis_args_model_names,
      bayestestR_equivalence_test_df_names,
      bayestestR_equivalence_test_model_names,
      bayestestR_p_direction_df_names,
      bayestestR_p_direction_model_names,
      bayestestR_p_direction_model_get_predicted_names)
  
  brms_hypothesis_args_df_args                    <- list()
  brms_hypothesis_args_model_args                 <- list()
  bayestestR_equivalence_test_df_args             <- list()
  bayestestR_equivalence_test_model_args          <- list()
  bayestestR_p_direction_df_args                  <- list()
  bayestestR_p_direction_model_args               <- list()
  bayestestR_p_direction_model_get_predicted_args <- list()
  
  for (i in brms_hypothesis_args_df_names) {
    if(i != 'x') {
      brms_hypothesis_args_df_args[[i]] <- full.arguments[[i]]
    }
  }
  for (i in brms_hypothesis_args_model_names) {
    if(i != 'x') {
      brms_hypothesis_args_model_args[[i]] <- full.arguments[[i]]
    }
  }
  for (i in bayestestR_equivalence_test_df_names) {
    if(i != 'x') {
      bayestestR_equivalence_test_df_args[[i]] <- full.arguments[[i]]
    }
  }
  for (i in bayestestR_equivalence_test_model_names) {
    if(i != 'x') {
      bayestestR_equivalence_test_model_args[[i]] <- full.arguments[[i]]
    }
  }
  for (i in bayestestR_p_direction_df_names) {
    if(i != 'x') {
      bayestestR_p_direction_df_args[[i]] <- full.arguments[[i]]
    }
  }
  for (i in bayestestR_p_direction_model_names) {
    if(i != 'x') {
      bayestestR_p_direction_model_args[[i]] <- full.arguments[[i]]
    }
  }
  for (i in bayestestR_p_direction_model_get_predicted_names) {
    if(i != 'x') {
      bayestestR_p_direction_model_get_predicted_args[[i]]<- full.arguments[[i]]
    }
  }
  
  
  
  
  
  if(!rope_test) {
    bayestestR_equivalence_test_df_args <- NULL
  }
  
  if(!pd_test) {
    bayestestR_p_direction_df_args <- NULL
  }
  
  
  ##############################################################################
  ##############################################################################
  
  # This when everything NULL, and default setting up the hypothesis
  obj_model_marginaleffects_direct_only <- FALSE
  if(obj_model_marginaleffects &
     !obj_model_mbcombo & 
     !obj_mfx_marginaleffects & 
     !obj_mfx_bayestestR & 
     !obj_mfx_mbcombo) {
    obj_model_marginaleffects_direct_only <- TRUE
  }
  
  ##############################################################################
  # Call marginal_* for obj_model_marginaleffects & obj_model_mbcombo
  # Evaluate all marginal_* mbcombo_*
  ##############################################################################
  
  # If there is no need for *_mbcombo but user requests evaluate_comparison
  # Then pass on to the get_comparison_hypothesis()
  if(!rope_test & !pd_test) {
    if(!is.null(hypothesis)) {
      if(NullFALSE(evaluate_comparison)) {
        obj_model_marginaleffects_direct_only <- FALSE
        obj_model_mbcombo <- TRUE
      }
    }
  }
  
  if(obj_model_marginaleffects | obj_model_mbcombo) {
    dots <- list(...)
    arguments_dots <- c(arguments, list(...))
    get_growthparameters_bgmfit_names <- 
      methods::formalArgs(get_growthparameters.bgmfit)
    get_growthparameters_bgmfit_args <- list()
    for (i in get_growthparameters_bgmfit_names) {
      if(i != "...") {
        get_growthparameters_bgmfit_args[[i]] <- arguments_dots[[i]]
      }
    }
    
    get_growthparameters_bgmfit_args_in <- 
      get_growthparameters_bgmfit_args
    
    get_growthparameters_bgmfit_args[['method']] <- method_call
    get_growthparameters_bgmfit_args[['equivalence_test']] <- NULL
    get_growthparameters_bgmfit_args[['p_direction']] <- NULL
    
    if(obj_model_marginaleffects_direct_only & !obj_model_mbcombo) {
      get_growthparameters_bgmfit_args[['pdrawsp']] <- FALSE
      get_growthparameters_bgmfit_args[['hypothesis']] <- hypothesis
    } else if(!obj_model_marginaleffects_direct_only & obj_model_mbcombo) {
      get_growthparameters_bgmfit_args[['pdrawsp']] <- TRUE
      # get_growthparameters_bgmfit_args_in <- 
      #   get_growthparameters_bgmfit_args[['hypothesis']]
      get_growthparameters_bgmfit_args[['hypothesis']] <- NULL
    } 
    
    
    for (i in get_growthparameters_bgmfit_names) {
      assign(i, get_growthparameters_bgmfit_args[[i]])
    }
    
    if(exists('resp')) {
      if(is.null(resp)) {
        resp <- NULL
      }
    } else {
      resp <- NULL
    }
    
    if(is.null(dpar)) {
      dpar <- "mu"
    }
    
    model <- getmodel_info(model = model, 
                           dpar = dpar, 
                           resp = resp, 
                           deriv = NULL, 
                           verbose = FALSE)
    
    get_growthparameters_bgmfit_args$model <- model
    get_growthparameters_bgmfit_args$usesavedfuns <- TRUE
    get_growthparameters_bgmfit_args$verbose <- FALSE
    
    data_draws <- CustomDoCall(get_growthparameters, 
                               get_growthparameters_bgmfit_args)
    
    if(obj_model_marginaleffects_direct_only) {
    } 
    
    if(obj_model_marginaleffects_direct_only & !obj_model_mbcombo) {
      return(data_draws)
    } else if(!obj_model_marginaleffects_direct_only & obj_model_mbcombo) {
      get_growthparameters_bgmfit_args[['hypothesis']] <-
        get_growthparameters_bgmfit_args_in[['hypothesis']]
    } 
    
    if(!data.table::is.data.table(data_draws)) {
      data_draws <- data.table::as.data.table(data_draws) 
    }
    if (!("draw" %in% names(data_draws)) && ("estimate" %in% names(data_draws))) {
      data_draws <- data.table::setnames(data_draws, "estimate", "draw")
    }
    # Convert character columns to factors (skip if already factor)
    char_cols <- names(data_draws)[sapply(data_draws, is.character)]
    data_draws[, (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]
    # return(data_draws)
  } # if(obj_model_marginaleffects | obj_model_mbcombo) {
  
  
  ##############################################################################
  # Evaluate all marginal_* mbcombo_*
  ##############################################################################
  
  if(obj_mfx_marginaleffects | obj_mfx_bayestestR | obj_mfx_mbcombo |
     obj_model_marginaleffects | obj_model_mbcombo) {
    
    if(obj_model_marginaleffects | obj_model_mbcombo) {
      data_draws <- data_draws
    } else {
      data_draws <- model
    }
    
    full.args                       <- arguments
    full.args[['equivalence_test']] <- bayestestR_equivalence_test_df_args
    full.args[['p_direction']]      <- bayestestR_p_direction_df_args
    
    # This call need only for marginal_* when deciding inline or not
    eqpdargs <- set_up_equivalence_test_p_direction_args(
      inbound_arguments = full.args, checking_inline = FALSE, 
      xcall = xcall, verbose = FALSE)
    
    
    
    full.args <- eqpdargs[['inbound_arguments']]
    check_equivalence_test_full.args <- 
      eqpdargs[['check_equivalence_test_full.args']]
    check_p_direction_full.args <- eqpdargs[['check_p_direction_full.args']]
    
    if(exists('get_growthparameters_bgmfit_args_in')) {
      p_direction <- get_growthparameters_bgmfit_args_in[['p_direction']]
      equivalence_test <- 
        get_growthparameters_bgmfit_args_in[['equivalence_test']]
    }
    
    if(is.null(equivalence_test)) {
      rope_test <- eqpdargs[['rope_test']]
    } else {
      rope_test <- rope_test
    }
    
    
    if(is.null(p_direction)) {
      pd_test <- eqpdargs[['pd_test']]
    } else {
      pd_test <- pd_test
    }
    
    if(is.null(get_range_null_form)) {
      get_range_null_form <- eqpdargs[['get_range_null_form']]
    }
    if(is.null(get_range_null_value)) {
      get_range_null_value <- eqpdargs[['get_range_null_value']]
    }
    if(is.null(format)) {
      format <- eqpdargs[['format']]
    }
    
    
    if(is.null(evaluate_hypothesis)) {
      if(exists('get_growthparameters_bgmfit_args')) {
        hypothesis <- get_growthparameters_bgmfit_args[['hypothesis']]
      }
      if(is.null(hypothesis)) {
        evaluate_hypothesis <- FALSE
      } else if(!is.null(hypothesis)) {
        evaluate_hypothesis <- TRUE
      }
    }
    
    
    # Note get_comparison_hypothesis() computes both 'comparisons' and 
    # 'hypothesis' but here we return only one, either 'comparisons' or
    # 'hypothesis'. Can make it genera;lize later by adding and argument later
    # that controls this behavior. Now this is done via arg evaluate_comparison
    
    if(is.null(evaluate_comparison)) {
      if(evaluate_hypothesis) {
        evaluate_comparison <- FALSE
      } else {
        evaluate_comparison <- TRUE
      }
    }
    
    full.args[['evaluate_hypothesis']] <- evaluate_hypothesis
    full.args[['evaluate_comparison']] <- evaluate_comparison
    
    
    get_comparison_hypothesis_args <- list()
    get_comparison_hypothesis_args[['data']] <- data_draws
    get_comparison_hypothesis_args[['full.args']] <- full.args
    get_comparison_hypothesis_args[['by']] <- full.args$by
    get_comparison_hypothesis_args[['evaluate_comparison']] <- 
      evaluate_comparison
    get_comparison_hypothesis_args[['evaluate_hypothesis']] <- 
      evaluate_hypothesis
    
    get_comparison_hypothesis_args[['rope_test']] <- rope_test
    get_comparison_hypothesis_args[['pd_test']] <- pd_test
    
    get_comparison_hypothesis_args[['get_range_null_form']] <- 
      get_range_null_form
    get_comparison_hypothesis_args[['get_range_null_value']] <- 
      get_range_null_value
    
    get_comparison_hypothesis_args[['comparison_by']] <- 
      comparison_by
    get_comparison_hypothesis_args[['comparison_range_null']] <- 
      comparison_range_null
    get_comparison_hypothesis_args[['comparison_range']] <- 
      comparison_range
    get_comparison_hypothesis_args[['comparison_null']] <- 
      comparison_null
    get_comparison_hypothesis_args[['hypothesis_by']] <- 
      hypothesis_by
    get_comparison_hypothesis_args[['hypothesis_range_null']] <- 
      hypothesis_range_null
    get_comparison_hypothesis_args[['hypothesis_range']] <- 
      hypothesis_range
    get_comparison_hypothesis_args[['hypothesis_null']] <- 
      hypothesis_null
    get_comparison_hypothesis_args[['rope_as_percent']] <- 
      rope_as_percent
    get_comparison_hypothesis_args[['pd_as_percent']] <- 
      pd_as_percent
    
    get_comparison_hypothesis_args[['format']] <- format
    get_comparison_hypothesis_args[['verbose']] <- verbose
    
    get_comparison_hypothesis_args[['parms_sat_elements']] <- parms_sat_elements
    
    
    out <- CustomDoCall(get_comparison_hypothesis, 
                        get_comparison_hypothesis_args)
    
    out <- DT_to_data_frames(out)
    
    if(is.null(reformat)) {
      out <- marginalstyle_reformat(out = out, 
                                    set_names_= set_names_)
    } else if(!is.null(reformat)) {
      if(reformat) out <- marginalstyle_reformat(out = out, 
                                                 set_names_= set_names_)
    }
    
    
  } # if(obj_mfx_marginaleffects | obj_mfx_bayestestR | obj_mfx_mbcombo) {
  
  
  ##############################################################################
  # Evaluate obj_model_brms
  ##############################################################################
  
  if(obj_model_brms) {
    brms_hypothesis_args_model_args[['x']] <- model
    brms_hypothesis_args_model_args[['hypothesis']] <- 
      brms_hypothesis_args_model_args[['hypothesis_str']]
    brms_hypothesis_args_model_args[['hypothesis_str']] <- NULL
    out <- do.call(brms::hypothesis, brms_hypothesis_args_model_args)
  }
  
  ##############################################################################
  # Evaluate obj_df_brms
  ##############################################################################
  
  if(obj_df_brms) {
    brms_hypothesis_args_df_args[['x']] <- model
    brms_hypothesis_args_df_args[['hypothesis']] <- 
      brms_hypothesis_args_df_args[['hypothesis_str']]
    brms_hypothesis_args_df_args[['hypothesis_str']] <- NULL
    out <- do.call(brms::hypothesis, brms_hypothesis_args_df_args)
  }
  
  
  ##############################################################################
  # Evaluate obj_model_bayestestR
  ##############################################################################
  
  param_names_bayestestR <- NULL
  effects_names_bayestestR <- NULL
  if(obj_model_bayestestR) {
    bayestestR_equivalence_test_model_args[['x']] <- model
    bayestestR_p_direction_model_args[['x']]      <- model
    
    if(rope_test & pd_test) {
      stop2("For engine = 'bayestestR', specify either 
            'rope_test' or 'pd_test', not both")
    } else if(!rope_test & !pd_test) {
      stop2("For engine = 'bayestestR', either 'rope_test'
            or 'pd_test' must be TRUE")
    } else if(rope_test) {
      param_names_bayestestR <-bayestestR_equivalence_test_model_args$parameters
      effects_names_bayestestR <- bayestestR_equivalence_test_model_args$effects
      out <- do.call(bayestestR::equivalence_test, 
                     bayestestR_equivalence_test_model_args)
    } else if(pd_test) {
      param_names_bayestestR <- bayestestR_p_direction_model_args$parameters
      effects_names_bayestestR <- bayestestR_p_direction_model_args$effects
      out <- do.call(bayestestR::p_direction, 
                     bayestestR_p_direction_model_args)
    } 
    
    # For nonlinear brms models, somehow parameters and effects are not filtered
    if(bayestestR_p_direction_model_args$effects == "fixed") {
      out <- out %>% 
        dplyr::filter(Parameter %in% param_names_bayestestR) %>% 
        dplyr::filter(Effects %in% effects_names_bayestestR)
    }
    
    param_names_bayestestR_cd_msg <- 
      paste0("The bayestestR output for the specified parameters %s is NULL. 
           Parameters must use full names with 'b_' prefix for fixed effects 
           (e.g., 'b_a_Intercept'). Run `brms::variables(model)` on your 
           fitted bsitar model to list all valid names.")
    
    if(nrow(out) == 0) {
      param_names_bayestestR_cd <- deparse(param_names_bayestestR)
      warning2c(sprintf(param_names_bayestestR_cd_msg,
                        param_names_bayestestR_cd))
    } else if(nrow(out) != length(param_names_bayestestR)) {
      param_names_bayestestR_cd <- setdiff(param_names_bayestestR, 
                                           out$Parameter)
      param_names_bayestestR_cd <- deparse(param_names_bayestestR_cd)
      warning2c(sprintf(param_names_bayestestR_cd_msg, 
                        param_names_bayestestR_cd))
    }
    # This "object_name" is must for plot
    attr(out, "object_name") <- 'model'
  } # if(obj_model_bayestestR) {
  
  ##############################################################################
  # Evaluate obj_df_bayestestR
  ##############################################################################
  
  if(obj_df_bayestestR) {
    if(plot) {
      stop2c("For engine = 'bayestestR', 
             plot is allowed only for model object, not for the data frame")
    }
    bayestestR_equivalence_test_df_args[['x']] <- model
    bayestestR_p_direction_df_args[['x']]      <- model
    
    if(rope_test & pd_test) {
      stop2("For engine = 'bayestestR', specify either 
            'rope_test' or 'pd_test', not both")
    } else if(!rope_test & !pd_test) {
      stop2("For engine = 'bayestestR', either 'rope_test'
            or 'pd_test' must be TRUE")
    } else if(rope_test) {
      out <- do.call(bayestestR::equivalence_test, 
                     bayestestR_equivalence_test_df_args)
    } else if(pd_test) {
      out <- do.call(bayestestR::p_direction, 
                     bayestestR_p_direction_df_args)
    } 
  } # if(obj_df_bayestestR) {
  
  
  ##############################################################################
  # plot - supported for model based brms and bayestestR
  ##############################################################################
  
  # patchwork::wrap_plots(plots, ncol = 2)  
  # cowplot::plot_grid(plotlist = plots, ncol = 2)
  # do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
  
  if(isTRUE(plot) | plot == "return") {
    if(!is.null(param_names_bayestestR)) {
      nplots <- length(param_names_bayestestR)
      if(nplots > 1) {
        plot_out_list <- list()
        for (i in param_names_bayestestR) {
          plot_out_list[[i]] <- out %>% dplyr::filter(Parameter == i) %>% 
            plot() 
        }
        # plot_out <- do.call(gridExtra::grid.arrange, 
        #                     c(plot_out_list, nrow = nplots))
        plot_out <- patchwork::wrap_plots(plot_out_list, nrow = nplots) +
          patchwork::plot_layout(guides = "keep")
      }
    } else if(nplots == 1) {
      plot_out <- plot(out) + jtools::theme_apa(legend.pos = 'bottom')
    }
  } else if(!is.logical(plot) & plot != "return") {
    stop2c("The plot should be either logical TRUE/FALSE, or a string 'return'")
  }
  
  
  if(isTRUE(plot)) {
    invisible(suppressMessages(print(plot_out)))
    # grid::grid.draw(plot_out)
  } else {
    if(plot == 'return') return(plot_out)
  }
  
  
  
  return(out)
}




#' @rdname hypothesis_test
#' @export
hypothesis_test <- function(model, ...) {
  UseMethod("hypothesis_test")
}


#' @noRd
#' @exportS3Method hypothesis_test data.frame
hypothesis_test.data.frame <- hypothesis_test.bgmfit

#' @noRd
#' @exportS3Method hypothesis_test data.table
hypothesis_test.data.table <- hypothesis_test.bgmfit



