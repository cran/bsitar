

#' @title Estimate and compare growth parameters for the Bayesian SITAR model
#' 
#' @description The \strong{get_growthparameters()} function estimates and
#'   compares growth parameters, such as peak growth velocity and the age at
#'   peak growth velocity. This function serves as a wrapper around
#'   [marginaleffects::comparisons()] and [marginaleffects::avg_comparisons()].
#'   The [marginaleffects::comparisons()] function computes unit-level
#'   (conditional) estimates, whereas [marginaleffects::avg_comparisons()]
#'   returns average (marginal) estimates. A detailed explanation is available
#'   [here](https://marginaleffects.com). Note that for the current use case—
#'   estimating and comparing growth parameters—the arguments \code{variables}
#'   and \code{comparison} in [marginaleffects::comparisons()] and
#'   [marginaleffects::avg_comparisons()] are modified (see below). Furthermore,
#'   comparisons of growth parameters are performed via the \code{hypothesis}
#'   argument of both the [marginaleffects::comparisons()] and
#'   [marginaleffects::avg_comparisons()] functions. Please note that the
#'   \pkg{marginaleffects} package is highly flexible, and users are expected to
#'   have a strong understanding of its workings. Additionally, since the
#'   \pkg{marginaleffects} package is rapidly evolving, results obtained from
#'   the current implementation should be considered experimental.
#'
#' @details The \code{get_growthparameters} function estimates and
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
#' @param parameter A single character string or character vector specifying the
#'   growth parameter(s) to estimate. Must be one of the following options,
#'   ordered by growth curve progression:
#'   
#'   \itemize{
#'     \item \code{'apgv'}: Age at peak growth velocity - age corresponding to
#'     the peak growth rate (i.e. maximum growth velocity) during adolescent
#'     growth spurt (default when \code{NULL})
#'     \item \code{'pgv'}: Peak growth velocity - maximum growth rate during 
#'       adolescent growth spurt
#'     \item \code{'atgv'}: Age at takeoff growth velocity - age where initial 
#'       rapid growth begins accelerating toward peak
#'     \item \code{'tgv'}: Takeoff growth velocity - initial high growth rate 
#'       at curve start
#'     \item \code{'acgv'}: Age at cessation growth velocity - age where growth 
#'       rate approaches near zero and growth effectively stops
#'     \item \code{'cgv'}: Cessation growth velocity - final low growth rate 
#'       approaching maturity
#'     
#'     \item \code{'spgv'}: Size (length/height) at peak growth velocity age,
#'     \code{'apgv'}
#'     \item \code{'stgv'}: Size at takeoff growth velocity age, \code{'atgv'}
#'     \item \code{'scgv'}: Size at cessation growth velocity age, \code{'acgv'}
#'     
#'     \item \code{'satXX'}: Size at specific age XX years of predictor
#'     \code{xvar} (e.g., \code{'sat15'} = size at 15 years, \code{'sat12.5'} =
#'     size at 12.5 years)
#'     
#'     \item \code{'all'}: All six primary velocity/age parameters
#'     (\code{'apgv','pgv','atgv','tgv','acgv','cgv'}). Cannot be used when
#'     \code{by = TRUE}.
#'   }
#'
#'   Multiple parameters can be requested simultaneously as a vector, e.g.,
#'   \code{c('apgv', 'pgv', 'spgv')}. Custom \code{'satXX'} format requires no
#'   spaces between \code{'sat'} prefix and numeric age.
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
#'   \code{0} and \code{model_deriv} to \code{FALSE}. If parameters are to be
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
#' @param marginals A \code{list}, \code{data.frame}, or \code{tibble} of
#'   velocity returned by the \pkg{marginaleffects} functions (default
#'   \code{NULL}). This is only evaluated when \code{method = 'custom'}. The
#'   \code{marginals} can be the output from \pkg{marginaleffects} functions or
#'   posterior draws from \code{marginaleffects::posterior_draws()}. The
#'   \code{marginals} argument is primarily used for internal purposes only.
#'   
#' @param preparms A \code{list}, \code{data.frame}, or \code{tibble} of pre
#'   computed parameters (default \code{NULL}). This is only evaluated when
#'   \code{method = 'custom'}. The \code{preparms} argument is primarily used
#'   for internal purposes only.
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
#'   \code{constrats_at = list(age = 'min')}) or as numeric vectors
#'   (e.g., \code{constrats_at = list(age = c(6, 7))}). If \code{constrats_at =
#'   NULL}, level-1 predictors (e.g., \code{age}) are automatically set to their
#'   unique values (i.e., \code{constrats_at = list(age = 'unique')}). To turn
#'   off this behavior, use \code{constrats_at = FALSE}. Note that
#'   \code{constrats_at} only affects data subsets prepared via
#'   [marginaleffects::datagrid()] or the \code{newdata} argument. The argument
#'   is evaluated only when \code{method = 'custom'} and \code{hypothesis} is
#'   not \code{NULL}.
#'
#' @param constrats_subset A named list (default \code{FALSE}) to subset draws
#'   for which contrasts should be computed via the \code{hypothesis} argument.
#'   This is similar to \code{constrats_at}, except that \code{constrats_subset}
#'   filters draws based on the character vector of variable names (e.g.,
#'   \code{constrats_subset = list(id = c('id1', 'id2'))}) rather than numeric
#'   values. The argument is evaluated only when \code{method = 'custom'} and
#'   \code{hypothesis} is not \code{NULL}.
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
#' @param equivalence_test A named arguments list (default \code{NULL}) passed
#'   to [bayestestR::equivalence_test()] for ROPE-based hypothesis testing. Core
#'   arguments:
#'   \itemize{
#'     \item \code{range = "default"}: ROPE bounds (numeric vector or
#'     \code{"default"}). Defaults values that are automatically set up for
#'     \code{range} are: \code{c(-1, 1)} for age parameters (\code{apgv},
#'     \code{atgv}, \code{acgv}); \code{c(-0.5, 0.5)} for velocity parameters
#'     (\code{pgv}, \code{tgv}, \code{cgv}). \item \code{ci}: Credible interval
#'     level. Default \code{0.95}. \item \code{rvar_col}: Random variable column
#'     name for long-format data. Default \code{NULL}. \item \code{verbose}:
#'     Print progress messages if \code{TRUE}. Default \code{FALSE}.
#'   }
#'   
#'   Additional package-specific controls:
#'   \itemize{
#'     \item \code{format}: If \code{TRUE}, merge
#'     \code{ROPE_low}/\code{ROPE_high} into \code{ROPE_range} and
#'     \code{HDI_low}/\code{HDI_high} into \code{HDI_range}. \item
#'     \code{reformat}: If \code{TRUE}, standardize
#'     [marginaleffects::comparisons()] output (\code{conf.low} → \code{Q2.5},
#'     \code{conf.high} → \code{Q97.5}) and drop reporting columns (\code{term},
#'     \code{contrast}, etc.).
#'     \item \code{get_range_null_form}: If \code{TRUE}, return expected
#'     \code{range} and \code{null} structure for manual specification via
#'     \code{comparison_range_null} and \code{hypothesis_range_null}.
#'     \item \code{digits}: Number of digits to use when printing numeric
#'     results. Default \code{2}.
#'     \item \code{as_percent}: Logical indicating whether to return ROPE 
#'     results as percentages (\code{TRUE}, default) or fractions (\code{FALSE}). 
#'     Only evaluated when \code{model} is class \code{"bsitar"} and
#'     \code{engine} is \code{"bayestestR"} or \code{"mbcombo"}.
#'     \item \code{na.rm}: If \code{TRUE} (default), remove \code{NA} values.
#'     \item \code{inline}: Internal use only; executes custom equivalence
#'     function (not for users).
#'   }
#' 
#' @param p_direction A named arguments list (default \code{NULL}) passed to
#'   [bayestestR::p_direction()] for Probability of Direction (\code{pd})
#'   computation. Commonly used elements:
#'   \itemize{
#'     \item \code{method}: Computation method. Default \code{"direct"}. \item
#'     \code{null}: Null value (default \code{0} for all growth parameters).
#'     \item \code{as_p}: Return p-value-like scale. Default \code{FALSE}. \item
#'     \code{remove_na}: Remove missing values. Default \code{TRUE}.
#'   }
#'   
#'   Additional package-specific controls:
#'   \itemize{
#'     \item \code{get_range_null_form}: If \code{TRUE}, return expected
#'     \code{range} and \code{null} structure for manual specification via
#'     \code{comparison_range_null} and \code{hypothesis_range_null}.
#'     \item \code{digits}: Number of digits to use when printing numeric
#'     results. Default \code{2}.
#'     \item \code{as_percent}: Logical indicating whether to return PD 
#'     results as percentages (\code{TRUE}, default) or fractions (\code{FALSE}).
#'     Only evaluated when \code{model} is class \code{"bsitar"} and
#'     \code{engine} is \code{"bayestestR"} or \code{"mbcombo"}.
#'     \item \code{na.rm}: If \code{TRUE} (default), remove \code{NA} values.
#'     \item \code{inline}: Internal use only; executes custom equivalence
#'     function (not for users).
#'   }
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
#'   over which the parameters are summarized. The summary statistics are
#'   computed for each unique combination of variables specified.
#'
#' @param future_splits A list, numeric vector, or logical value (default
#'   \code{TRUE}). Controls how posterior draws are partitioned into smaller
#'   subsets for parallel computation. This helps manage memory and improve
#'   performance, particularly on Linux when using
#'   \code{future::plan("multicore")}.
#'
#'   Input options:
#'   \itemize{
#'     \item \strong{List}: (e.g., \code{list(1:6, 7:10)}). Each element is a
#'     sequence of numbers passed to \code{draw_ids} to process in a separate
#'     chunk.
#'     \item \strong{Numeric vector of length 2}: (e.g., \code{c(100, 4)}). The
#'     first element is the total number of draws, and the second is the number
#'     of splits. Indices are generated via [parallel::splitIndices()].
#'     \item \strong{Numeric vector of length 1}: Specifies the total number of
#'     draws. Splits are calculated automatically.
#'     \item \strong{\code{TRUE}}: Automatically creates splits based on
#'     \code{ndraws} and \code{cores}. If both \code{ndraws} and \code{draw_ids}
#'     are specified, \code{draw_ids} takes precedence. This is consistent with
#'     post-processing functions in the \code{bsitar} and \code{brms} packages.
#'     \item \strong{\code{FALSE}}: Splitting is disabled.
#'   }
#'   \strong{Note}: On Windows with \code{future::plan("multisession")},
#'   background R processes may not free resources automatically. If needed, use
#'   [installr::kill_all_Rscript_s()] to terminate them.
#'
#' @param future_method A character string (default \code{'future'}) specifying
#'   the parallel computation backend. Options include:
#'   \itemize{
#'     \item \code{'future'}: Uses [future.apply::future_lapply()] for
#'     execution.
#'     \item \code{'dofuture'}: Uses [foreach::foreach()] via the
#'     \code{doFuture} package.
#'   }
#'
#' @param future_re_expose A logical value (default \code{NULL}) indicating whether 
#'   to re-expose internal \code{Stan} functions when \code{future = TRUE}. 
#'   This is critical when [future::plan()] is set to \code{"multisession"}, 
#'   as compiled C++ functions cannot be exported across distinct R sessions.
#'   \itemize{
#'     \item If \code{NULL} (default), it is automatically set to \code{TRUE} 
#'       when the plan is \code{"multisession"}.
#'     \item Explicitly setting this to \code{TRUE} is recommended for improved 
#'       performance during parallel execution.
#'   }
#' 
#' @inheritParams  get_predictions.bgmfit
#' @inheritParams  growthparameters.bgmfit
#' @inheritParams  marginaleffects::comparisons
#' @inheritParams  marginaleffects::avg_comparisons
#' @inheritParams  marginaleffects::plot_comparisons
#' @inheritParams  marginaleffects::datagrid
#' @inheritParams  brms::prepare_predictions
#' @inheritParams  brms::fitted.brmsfit
#'
#' @return A \code{data.frame} containing posterior estimates and credible
#'   intervals (CIs). The output includes the parameter estimates and their
#'   corresponding bounds (e.g., \code{Q2.5} and \code{Q97.5} for a 95%
#'   interval). The specific column names and structure may vary depending on
#'   the computation method.
#' 
#' @import data.table
#' 
#' @rdname get_growthparameters
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
#' get_growthparameters(model, parameter = 'apgv', ndraws = 10)
#' }
#' 
get_growthparameters.bgmfit <- function(model,
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
                                        parallel = FALSE,
                                        cores = NULL,
                                        fullframe = FALSE, 
                                        average = FALSE, 
                                        plot = FALSE, 
                                        mapping_facet = NULL,
                                        showlegends = NULL, 
                                        variables = NULL,
                                        condition = NULL,
                                        deriv = NULL,
                                        model_deriv = NULL,
                                        method = 'custom',
                                        marginals = NULL, 
                                        preparms = NULL,
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
                                        transform_draws = NULL,
                                        cross = FALSE,
                                        wts = NULL,
                                        hypothesis = NULL,
                                        equivalence = NULL,
                                        equivalence_test = NULL,
                                        p_direction = NULL,
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
    message2c("Please install the latest version of the 
              'marginaleffects' package",
              "\n ",
              "remotes::install_github('vincentarelbundock/marginaleffects')")
    return(invisible(NULL))
  }
  
  if(usecollapse) {
    usedtplyrcheck <- usedtplyr
    usedtplyr <- FALSE
    if(verbose & usedtplyrcheck) message2c("Setting usedtplyr = FALSE because ",
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
  
  callfuns     <- TRUE
  setmarginals <- FALSE
  setpreparms  <- FALSE
  
  if(!is.null(marginals) & !is.null(preparms)) {
    stop2c("Please specify either marginals or preparms, not both")
  }
  
  if(!is.null(marginals)) {
    setmarginals <- TRUE
    if(method == 'custom') callfuns <- FALSE
    if(method == 'pkg')    callfuns <- FALSE
  }
  
  if(!is.null(preparms)) {
    setmarginals <- TRUE
    setpreparms  <- TRUE
    if(method == 'custom') callfuns <- FALSE
    if(method == 'pkg')    callfuns <- FALSE
  }
  
  
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- envir
  }
  
  if(!is.null(transform) & !is.null(transform_draws)) {
    stop2c("Please specify either transform or transform_draws, not both")
  }
  
  
  
  if(is.null(dpar)) {
    dpar <- "mu"
  }
  
  model <- getmodel_info(model = model, 
                         dpar = dpar, 
                         resp = resp, 
                         deriv = NULL, 
                         verbose = verbose)
  
  # For sigma
  # Run this to get full data via modified get_data() for insight
  # See 'custom_get_data.brmsfit' in utils-helper-1
  if(!model$test_mode) {
    unlock_replace_bind(package = "insight", what = "get_data",
                        replacement = custom_get_data.brmsfit, 
                        ept_str = T)
    if(verbose) {
      message2c(" As model[['test_mode']] = FLASE, the full data by the",
                "\n ", 
                "insight::get_data() is extracted via 'custom_get_data.brmsfit'",
                "\n ", 
                "This full data is needed for marginaleffects functions",
                "\n ", 
                "'To over ride this approach, set model[['test_mode']] = TRUE")
    }
    
    unlock_replace_bind(package = "marginaleffects", what = "get_ci_draws",
                        replacement = custom_get_ci_draws, 
                        ept_str = T)
    
    # unlock_replace_bind(package = "marginaleffects", what = "equivalence",
    #                     replacement = custom_marginaleffects_equivalence, 
    #                     ept_str = T)
  } # if(!model$test_mode) {
  
  
  
  # Run this to get full data via modified get_data() for insight
  # See 'custom_get_data.brmsfit' in utils-helper-1
  if(!model$test_mode) {
    unlock_replace_bind(package = "insight", what = "get_data",
                        replacement = custom_get_data.brmsfit, ept_str = T)
    if(verbose) {
      message2c(" As model[['test_mode']] = FLASE, the full data by the",
                "\n ", 
                "insight::get_data() is extracted via 'custom_get_data.brmsfit'",
                "\n ", 
                "This full data is needed for marginaleffects functions",
                "\n ", 
                "'To over ride this approach, set model[['test_mode']] = TRUE")
    }
  } # if(!model$test_mode) {
  
  
  
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
  
  # This below for utilit heper 14 which is alos used in modelbased_
  # For sigma
  xvar_      <- paste0('xvar', resp_rev_)
  sigmaxvar_ <- paste0('sigma', xvar_)
  cov_       <- paste0('cov', resp_rev_)
  sigmacov_  <- paste0('sigma', cov_)
  uvarby     <- model$model_info$univariate_by$by
  if(is.null(uvarby)) uvarby <- NA 
  
  
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
  } # if(dpar == "mu") { else if(dpar == "sigma") {
  
  groupvar_     <- paste0('groupvar', resp_rev_)
  yvar_         <- paste0('yvar', resp_rev_)
  yvar          <- model$model_info[[yvar_]]
  hierarchical_ <- paste0('hierarchical', resp_rev_)
  
  
  if(is.null(levels_id) & is.null(idvar)) {
    idvar <- model$model_info[[groupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      idvar <- model$model_info[[hierarchical_]]
    }
    # 29.08.2025 - re assign idvar to groupvar_ if hierarchical_
    model$model_info[[groupvar_]] <- idvar # idvar[1]
  } else if (!is.null(levels_id)) {
    idvar <- levels_id
  } else if (!is.null(idvar)) {
    idvar <- idvar
  }
  
  cov_       <- paste0('cov', resp_rev_)
  sigmacov_  <- paste0('sigma', cov_)
  
  # When no random effects and hierarchical, IDvar <- NULL problem 02 03 2024
  if(is.null(idvar)) {
    if(is.null(idvar)) {
      if(!is.null(model$model_info[['idvars']])) {
        idvar <- model$model_info[['idvars']]
      }
    }
  }
  
  
  
  # # For sigma
  # xvar_      <- paste0('xvar', resp_rev_)
  # sigmaxvar_ <- paste0('sigma', xvar_)
  # cov_       <- paste0('cov', resp_rev_)
  # sigmacov_  <- paste0('sigma', cov_)
  # uvarby     <- model$model_info$univariate_by$by
  # 
  # if(dpar == "mu") {
  #   if(is.null(xvar)) {
  #     xvar   <- model$model_info[[xvar_]]
  #   }
  #   cov    <- model$model_info[[cov_]]
  # } else if(dpar == "sigma") {
  #   
  #   if(!is.na(model$model_info[[sigmaxvar_]])) {
  #     xvar   <- model$model_info[[sigmaxvar_]]
  #   } else if(is.na(model$model_info[[sigmaxvar_]]) & 
  #             !is.null(model$model_info[[xvar_]])) {
  #     xvar   <- model$model_info[[xvar_]]
  #   }
  #   
  #   cov    <- model$model_info[[sigmacov_]]
  # } # if(dpar == "mu") { else if(dpar == "sigma") {
  # 
  
  ########################################################
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
  funx_ <- NULL
  
  
  ########################################################
  
  itransform_set <- get_itransform_call(itransform = itransform,
                                        model = model, 
                                        newdata = newdata,
                                        dpar = dpar, 
                                        resp = resp,
                                        auto = TRUE,
                                        verbose = verbose)
  
  if(all(itransform_set == "")) {
    if(!isFALSE(pdrawsp)) {
      if(!is.character(pdrawsp)) pdrawsp <- "return"
      selectchoicesr <- c("return", 'add') 
      checkmate::assert_choice(pdrawsp, choices = selectchoicesr)
      if(pdrawsp == 'return' | pdrawsp == 'add') {
        ifunx_ <- function(x)x
      } 
    } # if(!isFALSE(pdrawsp)) {
  } # if(itransform_set == "") {
  
  
  
  # For sigma
  if(is.null(deriv)) {
    deriv <- 1
  } else {
    if(deriv == 0) {
      if(dpar == "sigma") {
        stop2c("Please set 'deriv = 1', or 'NULL'")
      }
    }
  }
  
  
  deriv.org       <- deriv
  model_deriv.org <- model_deriv
  # For growthparameetrs, always needed - i.e., need_velocity_curve = TRUE
  need_velocity_curve <- TRUE
  need_xvar_must      <- TRUE
  
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
  
  
  
  
  # 15 06 2025
  allowed_methods <- c('pkg', 'custom')
  if(!method %in% allowed_methods) 
    stop2c("Argument 'method' should be one of the following:",
           "\n ", 
           collapse_comma(allowed_methods)
    )
  
  
  
  if(method == 'custom') {
    deriv <- 1
    model_deriv <- TRUE
    if(verbose) {
      if(!setpreparms) {
        message2c(" For method = 'custom', deriv is set to TRUE.\n")
      }
    }
  }
  
  
  
  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  
  if(idata_method == 'm1') {
    stop2c("For marginaleffects based functions, the " ,
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
  
  
  
  # sat_ptc will be size at apgv/atgv/acgv whereas numeric_sat at 12, 13 etc 
  allowed_parms      <- c('apgv', 'pgv', 'atgv', 'tgv', 'acgv', 'cgv')
  allowed_parms_size <- c('spgv', 'stgv', 'scgv')
  check_set_parm_out <- check_set_parm(parameter = parameter,
                                       allowed_parms = allowed_parms,
                                       allowed_parms_size = allowed_parms_size,
                                       default_parms = c('apgv', 'pgv'),
                                       setpreparms = setpreparms,
                                       plot = plot,
                                       verbose = FALSE)
  eout_check_set_parm_out <- list2env(check_set_parm_out)
  for (eoutii in names(eout_check_set_parm_out)) {
    # if(!is.null(eout_check_set_parm_out[[eoutii]])) {
    #   assign(eoutii, eout_check_set_parm_out[[eoutii]])
    # }
    assign(eoutii, eout_check_set_parm_out[[eoutii]])
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
 
  
  peak_parm_TF <- takeoff_parm_TF <- cessation_parm_TF <- sat_parm_TF <- FALSE
  if('apgv' %in% parm | 'pgv' %in% parm | 'spgv' %in% parm) {
    peak_parm_TF <- TRUE
  }
  if('atgv' %in% parm | 'tgv' %in% parm | 'stgv' %in% parm) {
    takeoff_parm_TF <- TRUE
  }
  if('acgv' %in% parm | 'cgv' %in% parm | 'scgv' %in% parm) {
    cessation_parm_TF <- TRUE
  }
  if(!is.null(parameter_sat)) {
    sat_parm_TF <- TRUE
  }
  
  
  conf <- conf_level
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', probtitles)
  
  
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
  
  
  
  if(!is.null(funlist)) {
    if(!is.list(funlist)) {
      stop2c("'funlist' must be a list")
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
  
  
  
  call_predictions <- TRUE
  call_slopes      <- FALSE
  available_d1 <- o[['available_d1']]
  if(!available_d1) {
    deriv            <- 0
    model_deriv      <- FALSE
    call_predictions <- FALSE
    call_slopes      <- TRUE
    # re-get o[[2]] as _do
    post_processing_checks_args[['deriv']]    <- 0
    o    <- CustomDoCall(post_processing_checks, 
                         post_processing_checks_args)
  }
  
  post_processing_checks_args[['deriv']]    <- deriv
  
  
  ######################################################################
  ######################################################################
  
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
          } # if(!is.null(o[['sigma_model_is_ba_set_d0_as_d1']])) else if(
        } # if(o[['sigma_model']] == "ls") { else if(o[['sigma_model']] ...
      } # if(!is.null(o[['sigma_model']])) {
    } # if(deriv > 0) {
  } # if(dpar == "sigma") {
  
  
  ######################################################
  
  model$model_info[['difx']] <- difx
  # arguments$model$model_info[['difx']] <- difx
  
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
      # if(sigma_model == "basic" && !need_velocity_curve) {
      if(is.null(xvar)) {
        if(verbose) {
          message2c(clean_msg_sigma_model_no_xvar)
        }
      }
    }
    
    if(sigma_model != "ls" && need_velocity_curve) {
      # if(sigma_model == "basic" && need_velocity_curve) {
      # for deriv > 0, imp each id to have enough data points
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
              message2c("The difx has been set same as variables i.e.," , 
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
      stop2c("For dpar = 'sigma', please specify 'variables' 
             or 'difx' argument")
    }
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
  
  
  check_if_package_installed(model, xcall = xcall)
  
  model$xcall <- xcall
  
  
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
  
  
  
  get_future_args <- get_future_plan_args(future = future, 
                                          future_session = future_session, 
                                          oldfutureplan = future::plan(),
                                          setincores = setincores,
                                          verbose = FALSE)
  
  if(is.null(get_future_args)) get_future_args <- NULL
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
    } # if (inherits(future::plan(), "multicore")) {
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
        stop2c("future_splits must be a numeric vector of lenghth 2")
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
        stop2c("future_splits must be a numeric vector of lenghth 2")
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
                      set future_re_expose = TRUE")
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
            message2c("To speed up the calulations, it is advised ",
                      "to set 'future_re_expose = TRUE'")
          }
        } 
        if(need_future_re_expose_cpp & setplanis == "multisession") {
          stop2c("For plan 'multisession', the functions need to be ",
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
  
  full.args$model       <- model
  full.args$model_deriv <- model_deriv
  full.args$newdata     <- newdata
  
  if(!is.null(full.args$hypothesis)) {
    if(method == 'pkg') {
      if(!is.null(full.args$by)) {
        if(is.logical(full.args$by)) {
          stop2c("Argument 'by' is required for hypothesis")
        }
      }
    } else if(method == 'custom') {
      if(!is.null(full.args$by)) {
        if(is.logical(full.args$by)) {
          stop2c("Argument 'by' is required for hypothesis")
        }
      } else if(is.null(full.args$by)) {
        stop2c("Argument 'by' is required for hypothesis")
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
        # stop2c("Argument 'hypothesis' must be one of the following strings: ",
        #      collapse_comma(valid_hypothesis))
      }
    } else if(method == 'pkg') {
      if(is.character(full.args$hypothesis)) {
        stop2c("Argument 'hypothesis' must be one of the following form: ",
               collapse_comma(new_valid_hypothesis))
      }
    }
  }
  
  
  
  if(is.null(full.args$hypothesis) & is.null(full.args$equivalence)) {
    plot <- plot
  } else {
    plot <- FALSE
    if(verbose & plot) {
      message2c("Argument plot = TRUE is not allowed when either hypothesis ", 
                "or equivalence is not NULL",
                "\n ",
                "Therefor, setting 'plot = FALSE'") 
    }
  }
  
  
  
  
  full.args <- 
    sanitize_CustomDoCall_args(what = "CustomDoCall", 
                               arguments = full.args, 
                               check_formalArgs = 
                                 get_growthparameters.bgmfit,
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
  
  
  
  if(!exists('check_fun'))    check_fun    <- FALSE
  if(!exists('available_d1')) available_d1 <- FALSE
  
  
  
  if(!setpreparms) {
    full.args$ipts <- ipts <- check_ipts(ipts = full.args$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)
  }
  
  
  
  
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
    if(verbose) message2c("'transform' set based on 'transform_draws'")
  } else if(!is.null(transform_draws)) {
    full.args[['transform']] <- transform <- transform_draws
    if(verbose) message2c("'transform' set based on 'transform_draws'")
  } else {
    #####
  }
  full.args[['transform']] <- transform <- transform_draws
  
  
  comparisons_arguments <- full.args
  
  eqpdargs <- set_up_equivalence_test_p_direction_args(
    inbound_arguments = comparisons_arguments, 
    checking_inline = TRUE,
    xcall = xcall,
    verbose = FALSE)
  
  comparisons_arguments <- 
    eqpdargs[['inbound_arguments']]
  check_equivalence_test_full.args <- 
    eqpdargs[['check_equivalence_test_full.args']]
  check_p_direction_full.args <- 
    eqpdargs[['check_p_direction_full.args']]
  rope_test <- eqpdargs[['rope_test']]
  pd_test <- eqpdargs[['pd_test']]
  get_range_null_form <- eqpdargs[['get_range_null_form']]
  get_range_null_value <- eqpdargs[['get_range_null_value']]
  
  format <- eqpdargs[['format']]
  
  
  
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
      preparms,
      constrats_by,
      constrats_at,
      constrats_subset,
      usedtplyr,
      usecollapse, 
      parallel,
      cores,
      fullframe,
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
      internalmethod,
      xvar,
      difx,
      mapping_facet
    )
  ))[-1]
  
  if(plot) {
    exclude_args <- c(exclude_args, "cross")
  }
  
  for (exclude_argsi in exclude_args) {
    comparisons_arguments[[exclude_argsi]] <- NULL
  }
  
  
  set_variables <- setup_variables_var(model, 
                                       variables = variables, 
                                       xvar = xvar, 
                                       eps = eps, 
                                       method = method, 
                                       deriv = deriv, 
                                       model_deriv = model_deriv,
                                       difx = difx, 
                                       xcall = xcall,
                                       cov = cov, 
                                       dpar = dpar,
                                       xvar_strict = TRUE, 
                                       switch_plot = FALSE,
                                       verbose = verbose)
  
  
  allowed_comparison <- c('difference', 'differenceavg')
  
  if(!comparison %in% allowed_comparison) {
    stop2c("Allowed comparison options are ", 
           paste(paste0("'", allowed_comparison, "'"), collapse = ", ")
    )
  }
  
  
  if(comparison == 'differenceavg') {
    if(!average) {
      stop2c("For comparison = 'differenceavg' ", 
             ", the argument 'average' should be TRUE")
    }
    if(is.null(hypothesis)) {
      stop2c("For comparison = 'differenceavg' ", 
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
  
  
  # Skipping the below stop()
  # by could be 'id' only, not necessarily include xvar
  if(model$xcall == "modelbased_growthparameters" |
     model$xcall == "modelbased_growthparameters.bgmfit" |
     model$xcall == "get_growthparameters" |
     model$xcall == "get_growthparameters.bgmfit") {
    xvar_strict <- FALSE
  } else {
    xvar_strict <- TRUE
  }
  
  
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
                            xvar_strict = xvar_strict,
                            switch_plot = force_condition_and_by_switch_plot,
                            verbose = verbose)
  
  
  
  
  
  if (acg_velocity >= 1 | acg_velocity <= 0) {
    stop2c("The 'acg_velocity' should be set between 0.01 and 0.99")
  }
  
  
  
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
                                 FUN_other = setgrid_FUN_other
      )
      datagrid_arguments[[xvar]] <- setxvar
      if(setgrid_type == "mean_or_mode") {
        if(!isFALSE(set_group)) datagrid_arguments[['by']] <- set_group
      } else if(setgrid_type == "balanced") {
        if(!isFALSE(set_group)) datagrid_arguments[['by']] <- NULL
        comparisons_arguments[['by']] <- NULL
      }
      set_datagrid <- CustomDoCall(marginaleffects::datagrid, 
                                   datagrid_arguments)
      comparisons_arguments$newdata <- set_datagrid
    } else {
      stop2c("'datagrid' should be a data frame or named list")
    }
  } else if(is.null(datagrid)) {
    comparisons_arguments$newdata <- comparisons_arguments$newdata
  }
  
  comparisons_arguments[['datagrid']] <- NULL
  
  # Somehow draw_ids not passed correctly if not specified explicitly as arg
  get_draw_ids <- comparisons_arguments[['draw_ids']]
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
  comparisons_arguments[['draw_ids']] <- set_draw_ids

  gparms_fun = function(hi, lo, x, deriv, parm, by, ec_agg, ei_agg) {
    if(deriv == 0) y <- (hi - lo) / eps
    if(deriv > 0)  y <- (hi + lo) / 2
    
    # Here need to aggregate based on by argument
    if(aggregate_by) {
      try(insight::check_if_installed(c("grDevices", "stats"), stop = FALSE, 
                                      prompt = FALSE))
      xy_list <- grDevices::xy.coords(x, y)[1:2]
      tempxy <- collapse::qDF(xy_list)
      tempxy <- collapse::roworder(collapse::funique(tempxy), x)
      if (!isFALSE(by)) {
        if (ec_agg %in% c("mean", "median")) {
          agg_fun <- if(ec_agg == "mean") collapse::fmean else collapse::fmedian
          tempxy <- collapse::collap(collapse::na_omit(tempxy), ~ x, 
                                     FUN = agg_fun)
        }
      }
      x <- tempxy$x
      y <- tempxy$y
    } # if(aggregate_by) {
    parm_out <- list()
    tempout_p_exists <- FALSE
    for (i in parm) {
      if(i == 'apgv' | i == 'pgv') {
        tempout_p <- sitar::getPeak(x = x, y = y)
        if(i == 'apgv') parm_out[[i]] <- ifunx_(tempout_p[1])
        if(i ==  'pgv') parm_out[[i]] <- tempout_p[2]
        tempout_p_exists <- TRUE
      }
      if(i == 'atgv' | i == 'tgv') {
        tempout_t <- sitar::getTakeoff(x = x, y = y)
        if(i == 'atgv') parm_out[[i]] <- ifunx_(tempout_t[1])
        if(i ==  'tgv') parm_out[[i]] <- tempout_t[2]
      }
      if(i == 'acgv' | i == 'cgv') {
        if( tempout_p_exists) tempout_c <- tempout_p[2]
        if(!tempout_p_exists) tempout_c <- sitar::getPeak(x = x, y = y)[2]
        cgv  <- acg_velocity * tempout_c
        vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
        if(i == 'acgv') parm_out[[i]] <- ifunx_(x[vcgi])
        if(i ==  'cgv') parm_out[[i]] <- y[vcgi]
      }
    } # for (i in parm) {
    # for spgv, 
    if(length(parm) == 1) parm_out <- parm_out[[1]] 
    return(parm_out)
  } # gparms_fun
  
  
  
  
  
  
  call_comparison_gparms_fun <- function(parm, 
                                         eps, 
                                         comparisons_arguments, 
                                         aggregate_by, 
                                         newdata,
                                         ec_agg,
                                         ei_agg,
                                         enverr.,
                                         verbose) {
    formals(comparisons_arguments$comparison)[['by']] <- 
      comparisons_arguments[['by']]
    formals(comparisons_arguments$comparison)[['deriv']]  <- deriv
    formals(comparisons_arguments$comparison)[['parm']]   <- parm
    formals(comparisons_arguments$comparison)[['ec_agg']] <- ec_agg
    formals(comparisons_arguments$comparison)[['ei_agg']] <- ei_agg
    # For multisession, need to set all options()
    options("marginaleffects_posterior_center" = ec_agg)
    options("marginaleffects_posterior_interval" = ei_agg)
    
    
    suppresswar_args <-
      suppresswar_equivalence_test_p_direction_args(
        inbound_arguments = comparisons_arguments, parm = parm, verbose = FALSE)
    
    comparisons_arguments <- suppresswar_args[['inbound_arguments']]
    suppresswar           <- suppresswar_args[['suppresswar']]
    
    
    call_funcall <- function(funcall, args, suppresswar) {
      if(suppresswar) {
        suppressWarnings({return(CustomDoCall(funcall, args))})
      } else {
        return(CustomDoCall(funcall, args))
      }
    } # call_funcall
    
    if(!plot) {
      if(callfuns) {
        if(!average) setfuncall <- marginaleffects::comparisons
        if( average) setfuncall <- marginaleffects::avg_comparisons
      }
      out <- call_funcall(setfuncall, comparisons_arguments, suppresswar)
    } else if(plot) {
      if(isFALSE(set_group)) comparisons_arguments$by <- NULL
      if(is.null(comparisons_arguments[['by']]) &
         is.null(comparisons_arguments[['condition']])) {
        stop2c("For 'plot = TRUE', the 'by' argument should be specified")
      }
      setfuncall <- marginaleffects::plot_comparisons
      outp <- call_funcall(setfuncall, comparisons_arguments, suppresswar)
      # outp <- edit_mapping_facet() -> not relevenat here
      return(outp)
    } # if(!plot) { else if(plot) {
    return(out)
  } # call_comparison_gparms_fun
  
  
  
  #################
  enverr. <- environment()
  outer_call_comparison_gparms_fun <- function(parm, 
                                               eps,
                                               comparisons_arguments,
                                               aggregate_by,
                                               newdata,
                                               ec_agg,
                                               ei_agg,
                                               keep_mfx_draws,
                                               attr_mfx_draws,
                                               enverr.,
                                               verbose) {
    assign('err.', FALSE, envir = enverr.)
    # tryCatch needed when NA in one factor for hypothesis
    gout <- call_comparison_gparms_fun(parm = parm, eps = eps,
                                       comparisons_arguments = 
                                         comparisons_arguments,
                                       aggregate_by = aggregate_by,
                                       newdata = newdata,
                                       ec_agg = ec_agg,
                                       ei_agg = ei_agg,
                                       enverr. = enverr.,
                                       verbose = verbose)
    
    tryCatch(
      expr = {
        gout <- call_comparison_gparms_fun(parm = parm, eps = eps,
                                           comparisons_arguments = 
                                             comparisons_arguments,
                                           aggregate_by = aggregate_by,
                                           newdata = newdata,
                                           ec_agg = ec_agg,
                                           ei_agg = ei_agg,
                                           enverr. = enverr.,
                                           verbose = verbose)
      },
      error = function(e) {
        assign('err.', TRUE, envir = enverr.)
        if(verbose) message(e)
      }
    )
    err. <- get('err.', envir = enverr.)
    if(!exists('gout')) return(invisible(NULL))
    if(length(gout) == 0) err. <- TRUE
    if (err.) gout <- NULL
    if(plot) return(gout)
    # Create attr item here because next step deforms the mfx@draws 
    attr_gout_mfx <- attr_gout_mfx_draws <- NULL
    if(attr_mfx_draws) {
      if(keep_mfx_draws == 'mfx') {
        attr_gout_mfx <- gout
      } else if(keep_mfx_draws == 'mfx_draws') {
        attr_gout_mfx_draws <- marginaleffects::get_draws(gout)
      }
    }
    
    # If results not available for all factor, then add NA it
    # But this will deform the mfx@draws 
    by <- comparisons_arguments$by
    if(is.null(by)) by <- FALSE
    if(!is.null(gout)) {
      if(!isFALSE(by)) {
        if(is.null(hypothesis) && is.null(equivalence)) {
          goutnames <- names(gout)
          groupvars <-  eval(by)
          newdatajoin <- collapse::funique(newdata, cols = groupvars)
          gout_joined <- collapse::join(newdatajoin, gout, on = groupvars, 
                                        how = "left", verbose = FALSE)
          cols_to_keep <- setdiff(goutnames, c('predicted_lo', 'predicted_hi', 
                                               'predicted', 'tmp_idx'))
          gout <- gout_joined %>%  collapse::get_vars(cols_to_keep) %>% 
            collapse::ftransform(parameter = parm) # %>% collapse::qDF()  
        }
      }
    }
    
    if(!is.null(attr_gout_mfx)) {
      attr(gout, 'mfx') <- attr_gout_mfx
    } else if(!is.null(attr_gout_mfx_draws)) {
      attr(gout, 'mfx_draws') <- attr_gout_mfx_draws
    }
    return(gout)
  } # outer_call_comparison_gparms_fun
  
  
  
  
  ###############
  eval_re_formula <- eval(comparisons_arguments$re_formula)
  if(is.null(eval_re_formula)) {
    aggregate_by <- TRUE
  } else if(is.na(eval_re_formula)) {
    aggregate_by <- FALSE
  }
  if(is.null(showlegends)) {
    if(is.null(comparisons_arguments$re_formula)) {
      showlegends <- FALSE
    } else {
      showlegends <- TRUE
    }
  }
  out_sf_hy <- NULL
  allowed_methods <- c('pkg', 'custom')
  if(!method %in% allowed_methods) {
    stop2c("Argument 'method' should be one of the following:",
           "\n ", 
           collapse_comma(allowed_methods))
  }
  
  
  #############################################################################
  # method == 'pkg'
  #############################################################################
  
  if(method == 'pkg') {
    comparisons_arguments$variables  <- set_variables
    comparisons_arguments$by         <- set_group
    comparisons_arguments$comparison <- gparms_fun
    # For method = 'pkg' - doesn't matter call_predictions/call_slopes?
    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    comparisons_arguments$by <- setdiff(comparisons_arguments$by, xvar)
    if(is_emptyx(comparisons_arguments$by)) {
      comparisons_arguments$by <- NULL
    }
    for (i in names(comparisons_arguments$variables)) {
      if(i != xvar) comparisons_arguments$variables[[i]] <- NULL
    }
    
    if(!is.null(comparisons_arguments[['equivalence_test']])) {
      comparisons_arguments[['equivalence_test']][['by']] <- 
        comparisons_arguments[['by']]
      comparisons_arguments[['equivalence_test']][['hypothesis']] <- 
        comparisons_arguments[['hypothesis']]
      comparisons_arguments[['equivalence_test']][['equivalence']] <- 
        comparisons_arguments[['equivalence']]
    }
    
    if(!is.null(comparisons_arguments[['p_direction']])) {
      comparisons_arguments[['p_direction']][['by']] <- 
        comparisons_arguments[['by']]
    }
    
    # if gout is altered by adding NA, then attr_mfx_draws will allows attr
    list_mfx_draws <- TRUE # keep mfx_draws for post processing
    keep_mfx_draws <- 'mfx_draws' # If mfx 'gout', if mfx_draws -> get_draws()
    attr_mfx_draws <- list_mfx_draws
    
    if(plot) {
      out_sf <- outer_call_comparison_gparms_fun(
        parm = parm, eps = eps, 
        comparisons_arguments = comparisons_arguments,
        aggregate_by = aggregate_by,
        newdata = newdata,
        ec_agg = ec_agg,
        ei_agg = ei_agg,
        keep_mfx_draws = keep_mfx_draws,
        attr_mfx_draws = attr_mfx_draws,
        enverr. = enverr.,
        verbose = verbose) 
      return(out_sf)
    } else if(!plot) {
      # 1. Define a worker function
      run_comparison_worker <- function(p) {
        # Run the external function
        res <- outer_call_comparison_gparms_fun(
          parm = p, 
          eps = eps, 
          comparisons_arguments = comparisons_arguments,
          aggregate_by = aggregate_by,
          newdata = newdata,
          ec_agg = ec_agg,
          ei_agg = ei_agg,
          keep_mfx_draws = keep_mfx_draws,
          attr_mfx_draws = attr_mfx_draws,
          enverr. = enverr.,
          verbose = verbose)
        # Handle NULL results immediately
        if (is.null(res)) return(NULL)
        if (!"parameter" %in% names(res)) {
          res <- collapse::ftransform(res, 'parameter' = p)
          res <- collapse::colorder(res, parameter)
        }
        return(res)
      }
      
      
      
      # make method = 'pkg' also future ready
      if(future_splits_exe_future & callfuns) {
        
      }
      
      # out <-  future.apply::future_lapply(future_splits_at,
      #                                     future.envir = parent.frame(),
      #                                     future.globals = TRUE,
      #                                     future.seed = TRUE,
      #                                     FUN = myzfun)

      get_mfx_draws_fun <- function(x, list_cout, parm, 
                                    keep_mfx_draws, verbose) {
        parmi <- parm[x]
        if(keep_mfx_draws == "mfx_draws") {
          out <- attr(list_cout[[x]], "mfx_draws")
          if(is.null(out)) {
            if(verbose) warning2c("No draws for parameter ", parmi)
            return(NULL)
          }
          attr(out, 'class') <- c(attr(out, 'class'), "mfx_draws")
        } else if(keep_mfx_draws == "mfx") {
          out <- marginaleffects::get_draws(attr(list_cout[[x]], "mfx"))
          if(is.null(out)) {
            if(verbose) warning2c("No draws for parameter ", parmi)
            return(NULL)
          }
          attr(out, 'class') <- c(attr(out, 'class'), "mfx")
        }
        attr(out, 'class') <- c(attr(out, 'class'), "bgmfit")
        return(out)
      }
      
      # 2. Execution Logic
      if (length(parm) == 1) {
        out_sf <- run_comparison_worker(parm)
        if(list_mfx_draws) {
          draws_list <- lapply(1:length(parm), get_mfx_draws_fun, list(out_sf), 
                               parm, keep_mfx_draws, verbose)
          names(draws_list) <- parm
        } # if(list_mfx_draws) {
        draws_list_dt <- data.table::rbindlist(
          draws_list, 
          idcol = "parameter")[, 
                               parameter := factor(parameter)]
        
        attr(draws_list_dt, 'class') <- c(attr(draws_list_dt, 'class'), 
                                          c("bgmfit", 'mfx_draws'))
        
        if (is_emptyx(out_sf)) {
          if(verbose) {
            message2c("All draws are NA for the specified parameter: ",
                      collapse_comma(parm))
          }
          return(NA)
        }
        out_sf <- collapse::qDF(out_sf) 
      } else { # Multiple parameter case
        if (isTRUE(future)) {
          list_cout <- future.apply::future_lapply(parm, 
                                                   run_comparison_worker,
                                                   future.envir = enverr.,
                                                   future.stdout = TRUE,
                                                   future.globals = TRUE,
                                                   future.packages = NULL,
                                                   future.scheduling = 1.0,
                                                   future.chunk.size = Inf,
                                                   future.seed = FALSE)
        } else {
          list_cout <- lapply(parm, run_comparison_worker)
        }
        
        if(list_mfx_draws) {
          draws_list <- lapply(1:length(parm), get_mfx_draws_fun, list_cout, 
                               parm, keep_mfx_draws, verbose)
          names(draws_list) <- parm
          col_names_vars <- c('draw', 'estimate', 'conf.low', 'conf.high')
          draws_list <- populate_na_elements(draws_list, 
                                             col_names = col_names_vars)
          draws_list_dt <- 
            data.table::rbindlist(draws_list, 
                                  idcol = "parameter")[, 
                                                       parameter := 
                                                         factor(parameter)]
          
          attr(draws_list_dt, 'class') <- c(attr(draws_list_dt, 'class'), 
                                            c("bgmfit", 'mfx_draws'))
        }
        out_sf <- data.table::rbindlist(list_cout, fill = TRUE)
        if (is_emptyx(out_sf)) {
          if(verbose) {
            message2c("All draws are NA for the specified parameter: ",
                      collapse_comma(parm))
          }
          return(NA)
          # out_sf <- NA
        } else {
          out_sf <- collapse::qDF(out_sf)
        }
      } # if (length(parm) == 1) { else { # Multiple parameter case
      remove_it_pkg <- c('term', 'contrast')
      cols_to_keep <- setdiff(names(out_sf), remove_it_pkg)
      out_sf <- collapse::get_vars(out_sf, cols_to_keep)
      row.names(out_sf) <- NULL
    } # if(plot) { else if(!plot) {
  } # if(method == 'pkg') {
  

  if(method == 'pkg') {
    draws_list_dt <- 
      get_size_from_age_draws (age_draws_dt = draws_list_dt,
                               model = model,
                               xvar = xvar,
                               by = comparisons_arguments[['by']],
                               draw_ids = comparisons_arguments[['draw_ids']],
                               sat = numeric_sat,
                               sat_name = string_sat,
                               parameter = sat_ptc,
                               vat = NULL,
                               re_formula = re_formula,
                               parameter_name = 'parameter',
                               draw_name = 'draw',
                               drawid_name = 'drawid',
                               drawindex_name = 'drawindex',
                               before = NULL,
                               after = NULL,
                               first = FALSE,
                               last = TRUE,
                               skip_absent=FALSE)
    
    # If size is estimated, then trigger check_equivalence_test_full.args
    if(!is.null(sat_ptc) | !is.null(numeric_sat)) {
      parm <- unique(draws_list_dt[['parameter']])
      check_equivalence_test_full.args <- TRUE
    }
  } # if(method == 'pkg') {
  
  
  
  # pdrawsp can be used in get_comparison_hypothesis()
  if(method == 'pkg') {
    if(!isFALSE(pdrawsp)) {
      if(!is.character(pdrawsp)) pdrawsp <- "return"
      selectchoicesr <- c("return", 'add') 
      checkmate::assert_choice(pdrawsp, choices = selectchoicesr)
      if(pdrawsp == 'return') {
        full.args_by <- full.args$by
        if(!is.null(full.args_by)) {
          if(is.logical(full.args_by)) {
            if(!full.args_by) full.args_by <- NULL
          }
        }
        selectcols_custom <- c('parameter', 'drawid', 'draw', full.args_by)
        draws_list_dt <- draws_list_dt[, mget(selectcols_custom)]
        return(draws_list_dt)
      } else if(pdrawsp == 'add') {
        pdrawsp_est <- onex0
      } else {
        
      }
    }
  } # if(method == 'pkg') {
  
  
  if(method == 'pkg') {
    if(check_equivalence_test_full.args | check_p_direction_full.args) {
      if(is.null(full.args$by)) {
        selectcols_custom <- c('parameter', 'drawid', 'draw',  full.args$by)
      } else if(is.logical(full.args$by)) {
        if(!full.args$by) {
          selectcols_custom <- c('parameter', 'drawid', 'draw')
        } else if(full.args$by) {
          selectcols_custom <- c('parameter', 'drawid', 'draw')
        }
      } else if(!is.logical(full.args$by)) {
        checkmate::assert_character(full.args$by)
        selectcols_custom <- c('parameter', 'drawid', 'draw',  full.args$by)
      }
      draws_list_dt     <- draws_list_dt[, mget(selectcols_custom)]
      out_eqpt <- get_comparison_hypothesis(data = draws_list_dt, 
                                            full.args = full.args, 
                                            by = full.args$by,
                                            evaluate_comparison = TRUE,
                                            evaluate_hypothesis = TRUE,
                                            rope_test = rope_test,
                                            pd_test = pd_test,
                                            get_range_null_form = 
                                              get_range_null_form, 
                                            get_range_null_value = 
                                              get_range_null_value,
                                            parms_sat_elements = 
                                              parms_sat_elements,
                                            format = format,
                                            verbose = FALSE)
      
      if(!is.null(string_sat)) {
        out_eqpt <- out_eqpt[parameter == string_sat, 
                             parameter := string_numeric_sat]
      }
      
      out_eqpt <- DT_to_data_frames(out_eqpt)
      
      if(is.null(reformat)) {
        out_eqpt <- marginalstyle_reformat(out = out_eqpt, 
                                           set_names_ = set_names_)
      } else if(!is.null(reformat)) {
        if(reformat) out_eqpt <- marginalstyle_reformat(out = out_eqpt, 
                                                        set_names_ = set_names_)
      }
      return(out_eqpt)
    }
  } # if(method == 'pkg') {
  
  
  #############################################################################
  # End method == 'pkg'
  #############################################################################
  
  
  
  
  ##################################################################
  ##################################################################
  ##################################################################
  
  getparmsx <- function(x, y, parm = NULL, xvar = NULL, draw = NULL,
                        aggregate_by = FALSE, ...) {
    
    if(data.table::is.data.table(x) | is.data.frame(x)) {
      if(is.null(xvar)) {
        stop2c("Please specify the 'xvar' argument")
      }
      if(is.null(xvar)) {
        stop2c("Please specify the 'draw' argument")
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
    parm_out <- list()
    tempout_p_exists <- FALSE
    for (i in parm) {
      if(i == 'apgv' | i == 'pgv') {
        tempout_p <- sitar::getPeak(x = x, y = y)
        if(i == 'apgv') parm_out[[i]] <- ifunx_(tempout_p[1])
        if(i ==  'pgv') parm_out[[i]] <- tempout_p[2]
        tempout_p_exists <- TRUE
      }
      if(i == 'atgv' | i == 'tgv') {
        tempout_t <- sitar::getTakeoff(x = x, y = y)
        if(i == 'atgv') parm_out[[i]] <- ifunx_(tempout_t[1])
        if(i ==  'tgv') parm_out[[i]] <- tempout_t[2]
      }
      if(i == 'acgv' | i == 'cgv') {
        if( tempout_p_exists) tempout_c <- tempout_p[2]
        if(!tempout_p_exists) tempout_c <- sitar::getPeak(x = x, y = y)[2]
        cgv  <- acg_velocity * tempout_c
        vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
        if(i == 'acgv') parm_out[[i]] <- ifunx_(x[vcgi])
        if(i ==  'cgv') parm_out[[i]] <- y[vcgi]
      }
    } # for (i in parm) {
    # if(length(parm) == 1) parm_out <- parm_out[[1]] 
    # return(parm_out)
    out <- parm_out %>% CustomDoCall(cbind, .) %>% data.frame()
    return(out)
  } # getparmsx <- function
  
  ##################################################################
  
  
  
  # new 
  if(force_condition_and_by_switch_plot) {
    if(method == 'custom') {
      force_condition_and_by_switch_plot_arg <- comparisons_arguments
      force_condition_and_by_switch_plot_arg$difx <- full.args$difx
      force_condition_and_by_switch_plot_arg$variables <- full.args$variables
      force_condition_and_by_switch_plot_arg$by <- full.args$by
      
      force_condition_and_by_switch_plot_arg$transform <- FALSE
      force_condition_and_by_switch_plot_arg$pdrawso   <- TRUE
      marginals <- CustomDoCall(get_predictions, 
                                force_condition_and_by_switch_plot_arg)
      marginals <- marginals %>% marginaleffects::posterior_draws()
      setmarginals <- TRUE
      callfuns     <- FALSE
    } # if(method == 'custom') {
  } # if(force_condition_and_by_switch_plot) {
  
  
  if(force_condition_and_by_switch_plot) {
    summarise_over_x <- difx
  } else {
    summarise_over_x <- xvar
  }
  
  
  ###############################################
  
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
    
    ########################################################################
    if(call_slopes) {
      predictions_arguments[['comparison']]     <- NULL
      if(call_slopes) {
        if (!is.null(variables)) {
          if (!is.character(variables)) {
            stop2c("'variables' argument must be a character string such as", 
                   "\n ",
                   " variables = ", "'", summarise_over_x, "'"
            )
          } else {
            set_variables <- variables
            if(!summarise_over_x %in% variables) {
              if(is.null(difx)) {
                set_variables <- summarise_over_x
              } else if(!is.null(difx)) {
                set_variables <- difx
              }
              # set_variables <- summarise_over_x
            } else  {
              # if(!is.null(set_variables[[summarise_over_x]]))
            }
          }
        } else if (is.null(variables)) {
          if(is.null(difx)) {
            set_variables <- summarise_over_x
          } else if(!is.null(difx)) {
            set_variables <- difx
          }
        } 
      } # if(call_slopes) {
      
      # Decide if set by = NULL and then here pick and replace 'by' set_group 
      if(is.null(by)) {
        if(is.null(cov)) {
          set_group <- FALSE
        } else if(!is.null(cov)) {
          set_group <- cov
          if (!set_group %in% cov) {
            stop2c('by must be one of the ', cov)
          } 
        }
      } else if(!is.null(by)) {
        if (!isFALSE(by)) {
          set_group <- by
        } else if (isFALSE(by)) {
          set_group <- FALSE
        }
      }
      
      
      predictions_arguments$variables  <- set_variables
      predictions_arguments$by         <- set_group
      if(is.null(predictions_arguments$by)) {
        predictions_arguments$by < 'NULL'
      }
      assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    } # end if(call_slopes)
    
    ########################################################################
    
    
    
    # Imp, add xvar to 'by' if missing
    by <- predictions_arguments[['by']]
    
    if(isFALSE(by)) {
      by <- summarise_over_x
    } else if(!any(grepl(summarise_over_x, by))) {
      by <- c(summarise_over_x, eval(by))
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
        if(call_predictions) {
          out <- CustomDoCall(marginaleffects::predictions, 
                              predictions_arguments)
        } 
        if(call_slopes) {
          out <- CustomDoCall(marginaleffects::slopes, predictions_arguments)
        }
      } else if(average) {
        if(call_predictions) {
          out <- CustomDoCall(marginaleffects::avg_predictions, 
                              predictions_arguments)
        } 
        if(call_slopes) {
          out <- CustomDoCall(marginaleffects::avg_slopes, 
                              predictions_arguments)
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
            if(verbose) message2c("need to expose functions for 'multisession'")
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
            out <- CustomDoCall(marginaleffects::predictions, 
                                predictions_arguments)
          } 
          if(call_slopes) {
            out <- CustomDoCall(marginaleffects::slopes, 
                                predictions_arguments)
          }
          return(out)
        } # end myzfun
        out <-  future.apply::future_lapply(future_splits_at,
                                            future.envir = parent.frame(),
                                            future.globals = TRUE,
                                            future.seed = TRUE,
                                            future.packages = ('bsitar'),
                                            FUN = myzfun)
      } else if(average) {
        myzfun <- function(x) {
          predictions_arguments[['draw_ids']] <- x
          predictions_arguments[['ndraws']]   <- NULL
          # model <- predictions_arguments[['model']] 
          `%>%` <- bsitar::`%>%`
          if(re_expose) {
            if(verbose) message2c("need to expose functions for 'multisession'")
            predictions_arguments[['model']] <- 
              bsitar::expose_model_functions(predictions_arguments[['model']])
          }
          # Re-assign appropriate function
          setenv <- predictions_arguments[['model']]$model_info$envir
          environment(predictions_arguments) <- setenv
          assign(
            o[[1]],
            predictions_arguments[['model']]$model_info[['exefuns']][[o[[2]]]], 
            envir = setenv
          )
          if(call_predictions) {
            out <- CustomDoCall(marginaleffects::avg_predictions,
                                predictions_arguments)
          } 
          if(call_slopes) {
            out <- CustomDoCall(marginaleffects::avg_slopes, 
                                predictions_arguments)
          }
          return(out)
        } # end myzfun
        out <-  future.apply::future_lapply(future_splits_at,
                                            future.envir = parent.frame(),
                                            future.globals = TRUE,
                                            future.seed = TRUE,
                                            FUN = myzfun)
      } # if(!average) { else if(average) {
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
                                                   're_expose',
                                                   'o',
                                                   'call_predictions',
                                                   'call_slopes',
                                                   'CustomDoCall',
                                                   'verbose',
                                                   'predictions_arguments'))
        ) %doFuture_function% {
          x <- future_splits_at[[x]]
          predictions_arguments[['draw_ids']] <- x
          predictions_arguments[['ndraws']] <- NULL
          `%>%` <- bsitar::`%>%`
          if(re_expose) {
            if(verbose) message2c("need to expose functions for 'multisession'")
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
            out <- CustomDoCall(marginaleffects::predictions, 
                                predictions_arguments)
          } 
          if(call_slopes) {
            out <- CustomDoCall(marginaleffects::slopes, 
                                predictions_arguments)
          }
          return(out)
        }
      } else if(average) {
        out <- foreach::foreach(x = 1:length(future_splits_at),
                                .options.future = list(seed = TRUE),
                                .options.future =
                                  list(globals = c('future_splits_at',
                                                   're_expose',
                                                   'o',
                                                   'call_predictions',
                                                   'call_slopes',
                                                   'CustomDoCall',
                                                   'verbose',
                                                   'predictions_arguments'))
        ) %doFuture_function% {
          x <- future_splits_at[[x]]
          predictions_arguments[['draw_ids']] <- x
          predictions_arguments[['ndraws']] <- NULL
          `%>%` <- bsitar::`%>%`
          if(re_expose) {
            if(verbose) message2c("need to expose functions for 'multisession'")
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
            out <- CustomDoCall(marginaleffects::avg_predictions, 
                                predictions_arguments)
          } 
          if(call_slopes) {
            out <- CustomDoCall(marginaleffects::avg_slopes, 
                                predictions_arguments)
          }
          return(out)
        }
      } 
    } # if(future_splits_exe_dofuture) {
    
    
    
    
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
    
    
    # somehow this need consecutive number
    if(!future_splits_exe) {
      if(callfuns) {
        if(pdrawso) return(out)
        zxdraws <- out %>% marginaleffects::posterior_draws()
      }
    } else if(future_splits_exe) {
      if(callfuns) {
        if(pdrawso) {
          out <- out %>% CustomDoCall(rbind, .)
          return(out)
        }
        
        zxdraws <- collapse::unlist2d(future.apply::future_lapply(
          seq_along(out), posterior_draws_collapse), 
          idcols = FALSE, DT = TRUE)
        zxdraws[, drawid := as.factor(collapse::groupid(drawid))] 
        
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
    
    
    if(setpreparms) {
      if(inherits(preparms, 'list')) {
        zxdraws <-
          {. <- lapply(1:length(preparms), 
                       marginals_list_consecutive_drawid_function)
          list2DF(lapply(setNames(seq_along(.[[1]]), names(.[[1]])), function(i)
            unlist(lapply(., `[[`, i), FALSE, FALSE)))}
        zxdraws$drawid <- cheapr::factor_(zxdraws$drawid)
      } else {
        zxdraws <- preparms
      }
    }
    
    
    
    by_pdraws <- by
    
    # Imp, remove xvar from the by
    by <- base::setdiff(eval(by), eval(summarise_over_x)) 
    
    
    # For pre computed parameters, the below is not required
    if(!setpreparms) {
      if(usedtplyr) {
        getparmsx2                     <- getparmsx
        hypothesisargs                 <- formals(getparmsx2)
        hypothesisargs[['xvar']]       <- summarise_over_x
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
        # drawidby_ <- c(drawidby, 'parameter', 'estimate')
        drawidby_ <- c(drawidby, 'parameter', 'draw')
        parmest   <- 'draw'
        
        if(any(c('apgv','pgv') %in% parm)) getpest <- TRUE else getpest <- FALSE
        if(any(c('atgv','tgv') %in% parm)) gettest <- TRUE else gettest <- FALSE
        if(any(c('acgv','cgv') %in% parm)) getcest <- TRUE else getcest <- FALSE
        
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
          # if(is.null(dfp) & is.null(dft) & is.null(dfc)) return(x)
          cbind(c(namesp, namest, namesc), matrix(c(dfp, dft, dfc)))
        }
        
        onex0 <- zxdraws %>% 
          collapse::fgroup_by(drawidby) %>% 
          collapse::fsummarise(collapse::mctl(
            funx(cbind(.data[[summarise_over_x]], 
                       .data[[parmest]])), 
            names = F)) %>% 
          collapse::ftransformv(., 'V2', as.numeric) %>% 
          collapse::frename(., drawidby_) %>% 
          collapse::fsubset(., parameter %in% parm)
        
      } else {
        drawid_c <- list()
        for (drawidi in 1:nlevels(zxdraws$drawid)) {
          drawid_c[[drawidi]] <-  zxdraws %>% 
            dplyr::filter(drawid == drawidi) %>%
            dplyr::group_by_at(by) %>%
            dplyr::group_modify(., ~ getparmsx(.x[[summarise_over_x]] , 
                                               .x$draw, parm = parm),
                                .keep = TRUE) %>%
            dplyr::mutate(drawid = drawidi)
        }
        onex0 <- drawid_c %>% CustomDoCall(rbind, .) %>% data.frame()
      }
    } else if(setpreparms) {
      onex0 <- zxdraws
    } # if(!setpreparms) {
    
    
    
    
    #######################################################
    
    get_etix <- utils::getFromNamespace("get_eti", "marginaleffects")
    get_etix <- stats::quantile
    get_hdix <- utils::getFromNamespace("get_hdi", "marginaleffects")
    get_pe_ci <- function(x, draw = NULL, na.rm = TRUE, ...) {
      if(data.table::is.data.table(x) | is.data.frame(x)) {
        if(is.null(draw)) {
          stop2c("Please specify the 'draw' argument")
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
    
    

    onex0 <- 
      get_size_from_age_draws (age_draws_dt = onex0,
                               model = model,
                               xvar = xvar,
                               by = comparisons_arguments[['by']],
                               draw_ids = comparisons_arguments[['draw_ids']],
                               sat = numeric_sat,
                               sat_name = string_sat,
                               parameter = sat_ptc,
                               vat = NULL,
                               re_formula = re_formula,
                               parameter_name = 'parameter',
                               draw_name = 'draw',
                               drawid_name = 'drawid',
                               drawindex_name = 'drawindex',
                               before = NULL,
                               after = NULL,
                               first = FALSE,
                               last = TRUE,
                               skip_absent=FALSE)
    
    # If size is estimated, then trigger check_equivalence_test_full.args
    if(!is.null(sat_ptc) | !is.null(numeric_sat)) {
      parm <- unique(onex0[['parameter']])
      check_equivalence_test_full.args <- TRUE
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
          if(!isFALSE(by)) constrats_by <- setdiff(by, summarise_over_x) 
        }
      } else if(!is.null(constrats_by)) {
        constrats_by <- constrats_by
      }
    }
    
    
    if(!is.null(constrats_by)) {
      if(!is.character(constrats_by)) 
        stop2c("The 'constrats_by' argument should be a character string")
      for (axi in 1:length(constrats_by)) {
        caxi <- constrats_by[axi]
        if(!is.character(caxi)) {
          stop2c("The 'constrats_by' argument '", caxi," should be a character")
        }
        if(!caxi %in% by) {
          stop2c("The 'constrats_by' argument '", caxi, "' is not available in",
                 " the 'by' argument.",
                 "\n ", 
                 " Note that '", caxi, "' must be included in 'by' argument.",
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
      } else if(!is.null(constrats_at)) {
        if(!is.list(constrats_at)) stop2c("'constrats_at' must be a named list")
      }
    }
    
    
    # For hypothesis
    groupvarshyp1 <- c('drawid')
    groupvarshyp2 <- c('term')
    if(!is.null(constrats_at)) {
      for (caxi in names(constrats_at)) {
        if(!caxi %in% names(onex0)) {
          stop2c("Variable '", caxi, ". specified in 'constrats_at' is not in",
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
            stop2c(caxi, " specified in 'constrats_at' as character should be a 
                 single character:",
                   "\n ", 
                   collapse_comma(allowed_char_constrats_at)
            )
          }
          if(!constrats_at[[caxi]] %in% allowed_char_constrats_at) {
            stop2c(constrats_at[[caxi]], " specified in 'constrats_at' as 
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
          stop2c(caxi, " specified in 'constrats_at' has resulted in zero rows",
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
          stop2c("Variable '", caxi, ". specified in 'constrats_subset' is ",
                 " not avaialble in the 'by' argument.",
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
          stop2c(caxi, " specified in 'constrats_subset' resulted in zero rows",
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
      # parmi_estimate <- 'estimate'
      parmi_estimate <- 'draw'
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
          get_pe_ci_collapse(.data[['draw']],
                             # get_pe_ci_collapse(.data[['estimate']],
                             ec_agg = ec_agg, 
                             ei_agg = ei_agg, na.rm = TRUE, 
                             nthreads = arguments$cores, 
                             conf = conf, probs = probs))
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
      out3 <- summary_c %>% CustomDoCall(rbind, .) %>% data.frame()
      row.names(out3) <- NULL
      out_sf <- out3
    }
    
    
    # check_equivalence_test_full.args | check_p_direction_full.args
    
    # Hypothesis
    if(!is.null(hypothesis) | check_equivalence_test_full.args) {
      # if(!is.null(hypothesis)) {
      get_hypothesis_x_modify <- function(.x, hypothesis, by, draws, ...) {
        get_hypothesis_x(x = .x, hypothesis = hypothesis, by = by, 
                         draws = draws)
      }
      if(length(tibble::as_tibble(onex1)) > 25) {
        # if(nrow(onex1) > 25) {
        if(is.null(constrats_at)) {
          message2c("Note that the 'marginaleffects' package does not allow" ,
                    "\n",
                    "hypothesis argument when estimates rows are more than 25",
                    "\n",
                    "To avoid this issue, you can use 'constrats_at' argument",
                    "\n"
          )
        }
      }
      
      
      
      if(usedtplyr) {
        # parmi_estimate <- 'estimate'
        parmi_estimate <- 'draw'
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
        
        # if(!isFALSE(pdrawsh)) {
        #   setdrawidh          <- NULL
        #   setdrawidh_         <- NULL
        #   get_hypothesis_x_fx <- NULL
        #   temhyy <- 
        #     onex1 %>% collapse::fgroup_by( setdrawidh ) %>% 
        #   collapse::fsummarise(collapse::mctl(get_hypothesis_x_fx(.data))) %>%
        #     collapse::ftransformv(., 'V1', as.numeric) %>% 
        #     collapse::frename(., setdrawidh_)
        #   
        #   # if(!isFALSE(pdrawsh)) {
        #   selectchoicesr <- c("return", 'add') 
        #   checkmate::assert_choice(pdrawsh, choices = selectchoicesr)
        #   if(pdrawsh == 'return') {
        #     return(temhyy)
        #   } else if(pdrawsh == 'add') {
        #     pdrawsh_est <- temhyy
        #   } else {
        #     
        #   }
        # } # if(!isFALSE(pdrawsh)) {
        
        
        if(!data.table::is.data.table(onex1)) {
          onex1 <- data.table::as.data.table(onex1) 
        }
        
        onex1 <- clean_draws(onex1, 
                             variable = 'draw', 
                             group = 'parameter', 
                             verbose = FALSE)
        
        if(!is.null(constrats_by)) {
          hypothesis_by_what <- constrats_by
        } else {
          hypothesis_by_what <- full.args$by
        }
        
        
        
        # Here we can set evaluate_comparison = TRUE, that will be estimate
        out_sf_hy <- get_comparison_hypothesis(data = onex1, 
                                               full.args = full.args, 
                                               by = hypothesis_by_what,
                                               evaluate_comparison = TRUE,
                                               evaluate_hypothesis = TRUE,
                                               rope_test = rope_test,
                                               pd_test = pd_test,
                                               get_range_null_form = 
                                                 get_range_null_form, 
                                               get_range_null_value = 
                                                 get_range_null_value,
                                               parms_sat_elements = 
                                                 parms_sat_elements,
                                               format = format,
                                               verbose = FALSE)
        
        out_sf_hy <- out_sf_hy[['hypothesis']]
       
        # out_sf_hy <- data.table::setnames(out_sf_hy, "hypothesis", "term")
        out_sf_hy <- DT_to_data_frames(out_sf_hy)
        row.names(out_sf_hy) <- NULL
      } else {
        parmi_ci_c <- list()
        for (parmi in parm) {
          # parmi_estimate <- 'estimate'
          parmi_estimate <- 'draw'
          measurevar <- parmi
          parmi_ci_c[[parmi]] <-
            onex1 %>% dplyr::group_by_at(groupvarshyp1) %>% 
            dplyr::mutate(!! parmi_estimate := 
                            eval(parse(text = measurevar))) %>% 
            dplyr::group_modify(., ~get_hypothesis_x_modify(.x,
                                                            hypothesis = 
                                                              hypothesis, 
                                                            by = constrats_by, 
                                                            draws = draw
                                                            # draws = estimate
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
        out_sf_hy <- parmi_ci_c %>% CustomDoCall(rbind, .) %>% data.frame()
      }
    } # if(!is.null(hypothesis)) {
  } # if(method == 'custom') {
  
  #######################################
  
  if(!isFALSE(by)) {
    byarrange <- by 
  } else if(isFALSE(by)) {
    byarrange <- NULL
  } 
  
  if(method == 'pkg') {
    byarrange <- NULL
  }
  
  if(!is.null(byarrange)) {
    if(!is.null(bys)) byarrange <- bys else byarrange <- by
  }
  
  
  ##############################################################
  ##############################################################
  
  if(!is.null(string_sat)) {
    if(!is.null(out_sf)) {
      out_sf <- rename_vector_in_column_dt(out_sf, 
                                           column = 'parameter', 
                                           it = string_sat,
                                           by = string_numeric_sat)
    }
    if(!is.null(out_sf_hy)) {
      out_sf_hy <- rename_vector_in_column_dt(out_sf_hy, 
                                              column = 'parameter', 
                                              it = string_sat,
                                              by = string_numeric_sat)
    }
    if(!is.null(pdraws_est)) {
      pdraws_est <- rename_vector_in_column_dt(pdraws_est,
                                               column = 'parameter', 
                                               it = string_sat,
                                               by = string_numeric_sat)
    }
    if(!is.null(pdrawsp_est)) {
      pdrawsp_est <- rename_vector_in_column_dt(pdrawsp_est, 
                                                column = 'parameter', 
                                                it = string_sat,
                                                by = string_numeric_sat)
    }
    if(!is.null(pdrawsh_est)) {
      pdrawsh_est <- rename_vector_in_column_dt(pdrawsh_est, 
                                                column = 'parameter', 
                                                it = string_sat,
                                                by = string_numeric_sat)
    }
  } # if(!is.null(string_sat)) {
  
  
  ##############################################################
  ##############################################################
  
  
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
        dplyr::rename(!!as.symbol(set_names_[1]) := 
                        dplyr::all_of('estimate')) %>% 
        dplyr::rename(!!as.symbol(set_names_[2]) :=
                        dplyr::all_of('conf.low')) %>%
        dplyr::rename(!!as.symbol(set_names_[3]) :=
                        dplyr::all_of('conf.high')) %>%
        dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
        data.frame()
    } # if(!is.null(pdraws_est)) {
    
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
    out[['estimate']] <- NULL
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



#' @rdname get_growthparameters
#' @export
get_growthparameters <- function(model, ...) {
  UseMethod("get_growthparameters")
}

#' An alias of 'get_growthparameters()'
#' @rdname get_growthparameters
#' @export
#' 
get_gp <- function(model, ...) {
  UseMethod("get_growthparameters")
}


#' An alias of 'get_growthparameters()'
#' @rdname get_growthparameters
#' @export
growthparameters_comparison <- function(model, ...) {
  warning2c(
    "The function `growthparameters_comparison()` has been renamed to 
    `get_growthparameters()`.\n",
    "Please update your code to use `get_growthparameters()` instead",
    " of the old function with which will be removed in  
    the next release ",
    call. = FALSE
  )
  # .Deprecated("growthparameters_comparison")
  UseMethod("get_growthparameters")
}



#' An alias of get_growthparameters()
#' @rdname get_growthparameters
#' @export
marginal_growthparameters <- function(model, ...) {
  warning2c(
    "The function `marginal_growthparameters()` has been renamed to 
    `get_growthparameters()`.\n",
    "Please update your code to use `get_growthparameters()` instead",
    " of the old function with ``marginal_`` prefix which will be removed in  
    the next release ",
    "The new name better reflect the role of this function and to harmonise the
    naming scheme across the package. In particular, the earlier name with 
    the ``marginal_`` prefix unintentionally suggested that this function
    is used only for marginal inference, whereas they in fact 
    support both marginal and conditional inferences.",
    call. = FALSE
  )
  # .Deprecated("get_growthparameters")
  UseMethod("get_growthparameters")
}

