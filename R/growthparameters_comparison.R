

#' @title Estimate and Compare Growth Parameters
#' 
#' @description The \strong{growthparameters_comparison()} function estimates
#' and compares growth parameters, such as peak growth velocity and the age at
#' peak growth velocity. This function serves as a wrapper around
#' [marginaleffects::comparisons()] and [marginaleffects::avg_comparisons()].
#' The [marginaleffects::comparisons()] function computes unit-level
#' (conditional) estimates, whereas [marginaleffects::avg_comparisons()]
#' returns average (marginal) estimates. A detailed explanation is available
#' [here](https://marginaleffects.com). Note that for the current use case—
#' estimating and comparing growth parameters—the arguments \code{variables} and
#' \code{comparison} in [marginaleffects::comparisons()] and
#' [marginaleffects::avg_comparisons()] are modified (see below).
#' Furthermore, comparisons of growth parameters are performed via the
#' \code{hypothesis} argument of both the [marginaleffects::comparisons()] and
#' [marginaleffects::avg_comparisons()] functions. Please note that the
#' \pkg{marginaleffects} package is highly flexible, and users are expected to
#' have a strong understanding of its workings. Additionally, since the
#' \pkg{marginaleffects} package is rapidly evolving, results obtained from the
#' current implementation should be considered experimental.
#'
#' @details The \code{growthparameters_comparison} function estimates and
#' returns the following growth parameters:
#' \itemize{
#'   \item \code{pgv} - peak growth velocity
#'   \item \code{apgv} - age at peak growth velocity
#'   \item \code{tgv} - takeoff growth velocity
#'   \item \code{atgv} - age at takeoff growth velocity
#'   \item \code{cgv} - cessation growth velocity
#'   \item \code{acgv} - age at cessation growth velocity
#' }
#' 
#' The takeoff growth velocity is the lowest velocity just before the peak
#' begins, indicating the start of the pubertal growth spurt. The cessation
#' growth velocity marks the end of the active pubertal growth spurt and is
#' calculated as a certain percentage of the peak velocity (\code{pgv}).
#' Typically, 10 percent of \code{pgv} is considered a good indicator for the
#' cessation of the active pubertal growth spurt \insertCite{Anna2022}{bsitar}.
#' This percentage is controlled via the \code{acg_velocity} argument, which
#' accepts a positive real value bounded between 0 and 1 (with the default value
#' being \code{0.1}, indicating 10 percent).
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @param datagrid A grid of user-specified values to be used in the
#' \code{newdata} argument of various functions in the \pkg{marginaleffects}
#' package. This allows you to define the regions of the predictor space
#' where you want to evaluate the quantities of interest. See
#' [marginaleffects::datagrid()] for more details. By default, the
#' \code{datagrid} is set to \code{NULL}, meaning no custom grid is constructed.
#' To set a custom grid, the argument should either be a data frame created
#' using [marginaleffects::datagrid()], or a named list, which is internally
#' used for constructing the grid. For convenience, you can also pass an empty
#' list \code{datagrid = list()}, in which case essential arguments like
#' \code{model} and \code{newdata} are inferred from the respective arguments
#' specified elsewhere. Additionally, the first-level predictor (such as age)
#' and any covariates included in the model (e.g., gender) are automatically
#' inferred from the \code{model} object.
#' 
#' @param parameter A single character string or a character vector specifying
#'  the growth parameter(s) to be estimated. Options include \code{'tgv'}
#'  (takeoff growth velocity), \code{'atgv'} (age at takeoff growth velocity),
#'  \code{'pgv'} (peak growth velocity), \code{'apgv'} (age at peak growth
#'  velocity), \code{'cgv'} (cessation growth velocity), \code{'acgv'} (age at
#'  cessation growth velocity), and \code{'all'}. If \code{parameter = NULL}
#'  (default), age at peak growth velocity (\code{'apgv'}) is estimated. When
#'  \code{parameter = 'all'}, all six parameters are estimated. Note that the
#'  \code{'all'} option cannot be used when the \code{by} argument is set to
#'  \code{TRUE}.
#'
#' @param acg_velocity A real number specifying the percentage of peak growth
#'   velocity to be used as the cessation velocity when estimating the
#'   \code{cgv} and \code{acgv} growth parameters. The \code{acg_velocity}
#'   should be greater than \code{0} and less than \code{1}. The default value
#'   of \code{acg_velocity = 0.10} indicates that 10 percent of the peak growth
#'   velocity will be used to calculate the cessation growth velocity and the
#'   corresponding age at cessation velocity. For example, if the peak growth
#'   velocity estimate is \code{10 mm/year}, then the cessation growth velocity
#'   will be \code{1 mm/year}.
#' 
#' @param digits An integer (default \code{2}) specifying the number of decimal
#'   places to round the estimated growth parameters. The \code{digits} value is
#'   passed to the [base::round()] function.
#'
#' @param average A logical value indicating whether to internally call the
#'   [marginaleffects::comparisons()] or the
#'   [marginaleffects::avg_comparisons()] function. If \code{FALSE} (default),
#'   [marginaleffects::comparisons()] is called, otherwise
#'   [marginaleffects::avg_comparisons()] is used when \code{average = TRUE}.
#'
#' @param plot A logical value specifying whether to plot comparisons by calling
#'   the [marginaleffects::plot_comparisons()] function (\code{TRUE}) or not
#'   (\code{FALSE}). If \code{FALSE} (default), then
#'   [marginaleffects::comparisons()] or [marginaleffects::avg_comparisons()]
#'   are called to compute predictions (see the \code{average} argument for
#'   details).
#'
#' @param showlegends A logical value to specify whether to show legends
#'   (\code{TRUE}) or not (\code{FALSE}). If \code{NULL} (default), the value of
#'   \code{showlegends} is internally set to \code{TRUE} if \code{re_formula =
#'   NA}, and \code{FALSE} if \code{re_formula = NULL}.
#'
#' @param variables A named list specifying the level 1 predictor, such as
#'   \code{age} or \code{time}, used for estimating growth parameters in the
#'   current use case. The \code{variables} list is set via the \code{esp}
#'   argument (default value is \code{1e-6}). If \code{variables} is
#'   \code{NULL}, the relevant information is retrieved internally from the
#'   \code{model}. Alternatively, users can define \code{variables} as a named
#'   list, e.g., \code{variables = list('x' = 1e-6)} where \code{'x'} is the
#'   level 1 predictor. By default, \code{variables = list('age' = 1e-6)} in the
#'   \pkg{marginaleffects} package, as velocity is usually computed by
#'   differentiating the distance curve using the \code{dydx} approach. When
#'   using this default, the argument \code{deriv} is automatically set to
#'   \code{0} and \code{deriv_model} to \code{FALSE}. If parameters are to be
#'   estimated based on the model's first derivative, \code{deriv} must be set
#'   to \code{1} and \code{variables} will be defined as \code{variables =
#'   list('age' = 0)}. Note that if the default behavior is used (\code{deriv =
#'   0} and \code{variables = list('x' = 1e-6)}), additional arguments cannot be
#'   passed to \code{variables}. In contrast, when using an alternative approach
#'   (\code{deriv = 0} and \code{variables = list('x' = 0)}), additional options
#'   can be passed to the [marginaleffects::comparisons()] and
#'   [marginaleffects::avg_comparisons()] functions.
#'
#' @param method A character string indicating whether to compute estimates
#'   using the \code{'marginaleffects'} package (\code{method = 'pkg'}) or
#'   custom functions for efficiency and speed (\code{method = 'custom'},
#'   default). The \code{method = 'pkg'} option is only suitable for simple
#'   cases and should be used with caution. \code{method = 'custom'} is the
#'   preferred option because it allows for simultaneous estimation of multiple
#'   parameters (e.g., \code{'apgv'} and \code{'pgv'}). This method works during
#'   the post-draw stage, supports multiple parameter comparisons via the
#'   \code{hypothesis} argument, and allows users to add or return draws (see
#'   \code{pdraws} for details). If \code{method = 'pkg'}, the \code{by}
#'   argument must not contain the predictor (e.g., \code{age}), and
#'   \code{variables} must either be \code{NULL} (which defaults to
#'   \code{list(age = 1e-6)}) or a list with factor variables like
#'   \code{variables = list(class = 'pairwise')} or \code{variables = list(age =
#'   1e-6, class = 'pairwise')}. With \code{method = 'custom'}, the \code{by}
#'   argument can include predictors, which will be ignored, and
#'   \code{variables} should not contain predictors, but can accept factor
#'   variables as a vector (e.g., \code{variables = c('class')}). Using
#'   \code{method = 'custom'} is strongly recommended for better performance and
#'   flexibility.
#'
#' @param marginals A \code{list}, \code{data.frame}, or \code{tibble} returned
#'   by the \pkg{marginaleffects} functions (default \code{NULL}). This is only
#'   evaluated when \code{method = 'custom'}. The \code{marginals} can be the
#'   output from \pkg{marginaleffects} functions or posterior draws from
#'   \code{marginaleffects::posterior_draws()}. The \code{marginals} argument is
#'   primarily used for internal purposes.
#'
#' @param constrats_by A character vector (default \code{FALSE}) specifying the
#'   variable(s) by which estimates and contrasts (during the post-draw stage)
#'   should be computed using the \code{hypothesis} argument. The variable(s) in
#'   \code{constrats_by} should be a subset of those specified in the \code{by}
#'   argument. If \code{constrats_by = NULL}, it will copy all variables from
#'   \code{by}, except for the level-1 predictor (e.g., \code{age}). To disable
#'   this automatic behavior, use \code{constrats_by = FALSE}. This argument is
#'   evaluated only when \code{method = 'custom'} and \code{hypothesis} is not
#'   \code{NULL}.
#'
#' @param constrats_at A named list (default \code{FALSE}) to specify the values
#'   at which estimates and contrasts should be computed during the post-draw
#'   stage using the \code{hypothesis} argument. The values can be specified as
#'   \code{'max'}, \code{'min'}, \code{'unique'}, or \code{'range'} (e.g.,
#'   \code{constrats_at = list(age = 'min')}) or as numeric values or vectors
#'   (e.g., \code{constrats_at = list(age = c(6, 7))}). If \code{constrats_at =
#'   NULL}, level-1 predictors (e.g., \code{age}) are automatically set to their
#'   unique values (i.e., \code{constrats_at = list(age = 'unique')}). To turn
#'   off this behavior, use \code{constrats_at = FALSE}. Note that
#'   \code{constrats_at} only affects data subsets prepared via
#'   [marginaleffects::datagrid()] or the \code{newdata} argument. The argument
#'   is evaluated only when \code{method = 'custom'} and \code{hypothesis} is
#'   not \code{NULL}.
#'
#' @param constrats_subset A named list (default \code{FALSE}) to filter the
#'   estimates at which contrasts are computed using the \code{hypothesis}
#'   argument. This is similar to \code{constrats_at}, except that
#'   \code{constrats_subset} filters based on a character vector of variable
#'   names (e.g., \code{constrats_subset = list(id = c('id1', 'id2'))}) rather
#'   than numeric values. The argument is evaluated only when \code{method =
#'   'custom'} and \code{hypothesis} is not \code{NULL}.
#'
#' @param deriv A numeric value specifying whether to estimate parameters based
#'   on the differentiation of the distance curve or the model's first
#'   derivative. Please refer to the \code{variables} argument for more details.
#' 
#' @param comparison A character string specifying the comparison type for
#'   growth parameter estimation. Options are \code{'difference'} and
#'   \code{'differenceavg'}. This argument sets up the internal function for
#'   estimating parameters using [sitar::getPeak()], [sitar::getTakeoff()], and
#'   [sitar::getTrough()] functions. These options are restructured according to
#'   the user-specified \code{hypothesis} argument.
#'
#' @param reformat A logical (default \code{TRUE}) indicating whether to
#'   reformat the output returned by \code{marginaleffects} as a data frame.
#'   Column names are redefined as \code{conf.low} to \code{Q2.5} and
#'   \code{conf.high} to \code{Q97.5} (assuming \code{conf_int = 0.95}).
#'   Additionally, some columns (\code{term}, \code{contrast}, etc.) are dropped
#'   from the data frame.
#'
#' @param estimate_center A character string (default \code{NULL}) specifying
#'   how to center estimates: either \code{'mean'} or \code{'median'}. This
#'   option sets the global options as follows:
#'   \code{options("marginaleffects_posterior_center" = "mean")} or
#'   \code{options("marginaleffects_posterior_center" = "median")}. These global
#'   options are restored upon function exit using [base::on.exit()].
#'
#' @param estimate_interval A character string (default \code{NULL}) to specify
#'   the type of credible intervals: \code{'eti'} for equal-tailed intervals or
#'   \code{'hdi'} for highest density intervals. This option sets the global
#'   options as follows: \code{options("marginaleffects_posterior_interval" =
#'   "eti")} or \code{options("marginaleffects_posterior_interval" = "hdi")},
#'   and is restored on exit using [base::on.exit()].
#'
#' @param usedtplyr A logical (default \code{FALSE}) indicating whether to use
#'   the \pkg{dtplyr} package for summarizing the draws. This package uses
#'   \pkg{data.table} as a back-end. It is useful when the data has a large
#'   number of observations. For typical use cases, it does not make a
#'   significant performance difference. The \code{usedtplyr} argument is
#'   evaluated only when \code{method = 'custom'}.
#'
#' @param usecollapse A logical (default \code{FALSE}) to indicate whether to
#'   use the \pkg{collapse} package for summarizing the draws.
#' 
#' @param parallel A logical (default \code{FALSE}) indicating whether to use
#'   parallel computation (via \pkg{doParallel} and \pkg{foreach}) when
#'   \pkg{usecollapse = TRUE}. When \code{parallel = TRUE},
#'   [parallel::makeCluster()] sets the cluster type as \code{"PSOCK"}, which
#'   works on all operating systems, including \code{Windows}. If you want to
#'   use a faster option for Unix-based systems, you can set \code{parallel =
#'   "FORK"}, but note that it is not compatible with \code{Windows} systems.
#' 
#' @param cores A positive integer (default \code{1}) specifying the number of
#'   CPU cores to use when \code{parallel = TRUE}. To automatically detect the
#'   number of available cores, set \code{cores = NULL}. This is useful for
#'   optimizing performance when working with large datasets.
#'
#' @param pdraws A character string (default \code{FALSE}) that indicates
#'   whether to return the raw posterior draws. Options include:
#'   \itemize{
#'     \item \code{'return'}: returns the raw draws,
#'     \item \code{'add'}: adds the raw draws to the final return object,
#'     \item \code{'returns'}: returns the summary of the raw draws,
#'     \item \code{'adds'}: adds the summary of raw draws to the final return
#'     object.
#'   }
#'   The \code{pdraws} are the velocity estimates for each posterior sample. For
#'   more details, see [marginaleffects::posterior_draws()].
#'
#' @param pdrawso A character string (default \code{FALSE}) to indicate whether
#'   to return the original posterior draws for parameters. Options include:
#'   \itemize{
#'     \item \code{'return'}: returns the original posterior draws,
#'     \item \code{'add'}: adds the original posterior draws to the outcome.
#'   }
#'   When \code{pdrawso = TRUE}, the default behavior is \code{pdrawso =
#'   'return'}. Note that the posterior draws are returned before calling
#'   [marginaleffects::posterior_draws()].
#' 
#' @param pdrawsp A character string (default \code{FALSE}) to indicate whether
#'   to return the posterior draws for parameters. Options include:
#'   \itemize{
#'     \item \code{'return'}: returns the posterior draws for parameters,
#'     \item \code{'add'}: adds the posterior draws to the outcome.
#'   }
#'   When \code{pdrawsp = TRUE}, the default behavior is \code{pdrawsp =
#'   'return'}. The \code{pdrawsp} represent the parameter estimates for each of
#'   the posterior samples, and the summary of these are the estimates returned.
#' 
#' @param pdrawsh A character string (default \code{FALSE}) to indicate whether
#'   to return the posterior draws for parameter contrasts. Options include:
#'   \itemize{
#'     \item \code{'return'}: returns the posterior draws for contrasts.
#'   }
#'   The summary of posterior draws for parameters is the default returned
#'   object. The \code{pdrawsh} represent the contrast estimates for each of the
#'   posterior samples, and the summary of these are the contrast returned.
#' 
#' @param bys A character string (default \code{NULL}) specifying the variables
#'   over which the parameters need to be summarized. If \code{bys} is not
#'   \code{NULL}, the summary statistics will be calculated for each unique
#'   combination of the specified variables.
#'
#' @param future_splits A list (default \code{NULL}) that can be an unnamed
#'   numeric list, a logical value, or a numeric vector of length 1 or 2. It is
#'   used to split the processing of posterior draws into smaller subsets for
#'   parallel computation.
#'   - If passed as a list (e.g., \code{future_splits = list(1:6, 7:10)}), 
#'   each sequence of
#'   numbers is passed to the \code{draw_ids} argument.
#'   - If passed as a numeric vector (e.g., \code{future_splits = c(10, 2)}), 
#'   the first element
#'   specifies the number of draws (see \code{draw_ids}) and the second element
#'   indicates the number of splits. The splits are created using
#'   [parallel::splitIndices()].
#'   - If passed as a numeric vector of length 1, the first element is 
#'   internally set as the
#'   number of draws (\code{ndraws} or \code{draw_ids}) depending on which one
#'   is not \code{NULL}.
#'   - If \code{TRUE}, a numeric vector for \code{future_splits} is created 
#'   based on the number
#'   of draws (\code{ndraws}) and the number of cores (\code{cores}).
#'   - If \code{FALSE}, \code{future_splits} is ignored.
#'   The use case for \code{future_splits} is to save memory and improve
#'   performance, especially on \code{Linux} systems when \code{future::plan()}
#'   is set to \code{multicore}. Note: on Windows systems, R processes may not
#'   be freed automatically when using \code{'multisession'}. In such cases, the
#'   R processes can be interrupted using [installr::kill_all_Rscript_s()].
#'
#' @param future_method A character string (default \code{'future'}) to specify
#'   the method for parallel computation. Options include:
#'   \itemize{
#'     \item \code{'future'}: Uses [future::future()] along with
#'     [future.apply::future_lapply()] for parallel execution.
#'     \item \code{'foreach'}: Uses [foreach::foreach()] with the
#'     \code{'dofuture'} function from the \code{doFuture} package for parallel
#'     execution.
#'   }
#'
#' @param future_re_expose A logical (default \code{NULL}) to indicate whether
#'   to re-expose \code{Stan} functions when \code{future = TRUE}. This is
#'   especially relevant when [future::plan()] is set to \code{'multisession'},
#'   as already exposed C++ \code{Stan} functions cannot be passed across
#'   multiple sessions.
#'
#'   - When \code{future_re_expose = NULL} (the default), \code{future_re_expose} 
#'   is automatically set to \code{TRUE} for the \code{'multisession'} plan.
#'   - It is advised to explicitly set \code{future_re_expose = TRUE} for speed 
#'   gains when using parallel processing with \code{future = TRUE}.
#'
#' 
#' @inheritParams  growthparameters.bgmfit
#' @inheritParams  marginaleffects::comparisons
#' @inheritParams  marginaleffects::avg_comparisons
#' @inheritParams  marginaleffects::plot_comparisons
#' @inheritParams  marginaleffects::datagrid
#' @inheritParams  brms::fitted.brmsfit
#'
#' @return A data frame object with estimates and credible intervals (CIs) for
#'   the computed parameter(s). The returned data frame includes the parameter
#'   estimates, along with lower and upper bounds of the credible intervals,
#'   typically labeled as \code{Q2.5} and \code{Q97.5}, assuming a 95%
#'   confidence level. The specific columns may vary depending on the
#'   computation method and the parameters being estimated.
#' 
#' @import data.table
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
#' # Note that since no covariate is part of the model fit, the below example 
#' # doesn't make sense and is included here only for the purpose of completeness.
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
                                   average = FALSE, 
                                   plot = FALSE, 
                                   showlegends = NULL, 
                                   variables = NULL,
                                   deriv = NULL,
                                   deriv_model = NULL,
                                   method = 'custom',
                                   marginals = NULL, 
                                   pdraws = FALSE, 
                                   pdrawso = FALSE,
                                   pdrawsp = FALSE, 
                                   pdrawsh = FALSE, 
                                   comparison = "difference",
                                   type = NULL,
                                   by = FALSE,
                                   bys = NULL,
                                   conf_level = 0.95,
                                   transform = NULL,
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
                                   envir = NULL, ...
                                   ) {
  
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
   
  
  
  # 15 06 2025
  allowed_methods <- c('pkg', 'custom')
  if(!method %in% allowed_methods) 
    stop("Argument 'method' should be one of the following:",
         "\n ", 
         collapse_comma(allowed_methods)
    )
  
  
  
  if(method == 'custom') {
    deriv <- 1
    deriv_model <- TRUE
    if(verbose) message("for method = 'custom', deriv set to TRUE")
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
  
  allowed_parms <- c(
    'apgv',
    'pgv',
    'atgv',
    'tgv',
    'acgv',
    'cgv'
    )
  
  if(is.null(parameter)) parameter <- 'apgv'
  
  parameter <- base::tolower(parameter)
  
  if (is.null(parameter)) {
    parm <- 'apgv' 
  } else if(length(parameter) == 1 && parameter == 'all') {
    parm <- allowed_parms 
  } else if(length(parameter) == 1) {
    parm <- parameter
  } else if(length(parameter) > 1) {
    # parameter <- base::tolower(parameter)
    for (parameteri in parameter) {
      if(!parameteri %in% allowed_parms) {
        allowed_parms_err <- c(allowed_parms, 'all')
        stop("Allowed parameter options are ", 
             paste(paste0("'", allowed_parms_err, "'"), collapse = ", ")
        )
      }
    }
    parm <- parameter
  }
  parm <- base::tolower(parm)
  
  
  
  if(length(parm) > 1) {
    if(plot) stop("Please specify only one parameter when plot = TRUE")
  }
  
  
  conf <- conf_level
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
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
                                  fargs = arguments, 
                                  dargs = list(...), 
                                  verbose = verbose)
  
  full.args$model <- model
  full.args$deriv_model <- deriv_model
  
 
  if(is.null(full.args$hypothesis) & is.null(full.args$equivalence)) {
    plot <- plot
  } else {
    plot <- FALSE
    if(verbose & plot) {
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
      future_splits,
      future_method,
      future_re_expose,
      cores,
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
      funlist
    )
  ))[-1]
  
  if(plot) {
    exclude_args <- c(exclude_args, "cross")
  }
  
  for (exclude_argsi in exclude_args) {
    comparisons_arguments[[exclude_argsi]] <- NULL
  }

  
  if(deriv == 0 & deriv_model) 
    stop("If argument 'deriv_model' = TRUE, the argument 'deriv' should be 1")
  if(deriv == 1 & !deriv_model) 
    stop("If argument 'deriv_model' = FALSE, the argument 'deriv' should be 0")
  
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
  
  # set_variables <- variables
  # print(variables)

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
    # For marginal_draws(...,  plot = T), either condition or by allowed
    # Therefore, when plot = T, condition is kept and by dropped, 
    # otherwise by is kept and condition dropped
    
    exclude_args_con_by <- exclude_args
    
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
  
  
  pdrawsp_est <- NULL
  pdrawsh_est <- NULL
  pdraws_est <- NULL
  
  if(method == 'custom') {
    predictions_arguments <- comparisons_arguments
    predictions_arguments[['cross']] <- NULL
    predictions_arguments[['method']] <- NULL
    predictions_arguments[['hypothesis']] <- NULL # hypothesis evaluated later
    
    
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
    
    by <- eval(by)
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
        out <- do.call(marginaleffects::predictions, predictions_arguments)
      } else if(average) {
        out <- do.call(marginaleffects::avg_predictions, predictions_arguments)
      }
      # out <- out %>% marginaleffects::posterior_draws()
      # if(pdrawso) return(out)
      # zxdraws <- out # %>% marginaleffects::posterior_draws()
    } # if(!future_splits_exe) {
    
    
    
    if(future_splits_exe_future & callfuns) {
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
        zxdraws <- lapply(1:length(future_splits_at),  FUN = posterior_draws_function)
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
          {. <- lapply(1:length(marginals), marginals_list_consecutive_drawid_function)
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
        if('apgv' %in% parmi) parm_c[[parmi]] <- sitar::getPeak(x, y)[1]
        if('pgv' %in% parmi) parm_c[[parmi]] <- pgvx <- sitar::getPeak(x, y)[2]
        if('atgv' %in% parmi) parm_c[[parmi]] <- sitar::getTakeoff(x, y)[1]
        if('tgv' %in% parmi) parm_c[[parmi]] <- sitar::getTakeoff(x, y)[2]
        if('acgv' %in% parmi | 'acgv' %in% parmi) {
          if(is.null(pgvx)) pgvx <- sitar::getPeak(x, y)[2]
          cgv  <- acg_velocity * pgvx
          vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
        }
        if('acgv' %in% parmi) parm_c[[parmi]] <- x[vcgi]
        if('cgv' %in% parmi) parm_c[[parmi]] <- y[vcgi]
      }
      out <- parm_c %>% do.call(cbind, .) %>% data.frame()
      out
    }
   
    
    if(usedtplyr) {
      getparmsx2 <- getparmsx
      hypothesisargs <- formals(getparmsx2)
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
        if(getpest) dfp <- sitar::getPeak(x[,1], x[,2]) else dfp <- NULL
        if(gettest) dft <- sitar::getTakeoff(x[,1], x[,2]) else dft <- NULL
        if(getcest) dfc <- getcgvfunc(x[,1], x[,2]) else dfc <- NULL
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
        drawid_c[[drawidi]] <-  zxdraws %>% dplyr::filter(drawid == drawidi) %>%
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
        # Actually no need to summarise again because estimate already is summary
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
        # For pdrawsp_est and pdrawsh_est, there are no conf columns, only estimates
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
  if(!is.null(pdraws_est)) pdraws_est <- as_factor_as_character_factor_df(pdraws_est)
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




#' @rdname growthparameters_comparison.bgmfit
#' @export
growthparameters_comparison <- function(model, ...) {
  UseMethod("growthparameters_comparison")
}


