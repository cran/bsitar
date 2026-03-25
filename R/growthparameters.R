

#' @title Estimate growth parameters for the Bayesian SITAR model
#'
#' @description The \strong{growthparameters()} function estimates both
#'   population-average and individual-specific growth parameters (e.g., age at
#'   peak growth velocity). It also provides measures of uncertainty, including
#'   standard errors (SE) and credible intervals (CIs). For a more advanced
#'   analysis, consider using the [get_growthparameters()] function, which not
#'   only estimates adjusted parameters but also enables comparisons of these
#'   parameters across different groups.The [get_growthparameters()] function is
#'   based on the \pkg{marginaleffects} package.
#'
#' @details The \strong{growthparameters()} function internally calls either the
#'   [fitted_draws()] or the [predict_draws()] function to estimate
#'   first-derivative growth parameters for each posterior draw. The estimated
#'   growth parameters include:
#'   - Age at Peak Growth Velocity (APGV)
#'   - Peak Growth Velocity (PGV)
#'   - Age at Takeoff Growth Velocity (ATGV)
#'   - Takeoff Growth Velocity (TGV)
#'   - Age at Cessation of Growth Velocity (ACGV)
#'   - Cessation Growth Velocity (CGV)
#'
#'   APGV and PGV are estimated using the [sitar::getPeak()] function, while
#'   ATGV and TGV are estimated using the [sitar::getTakeoff()] function. The
#'   [sitar::getTrough()] function is employed to estimate ACGV and CGV. The
#'   parameters from each posterior draw are then summarized to provide
#'   estimates along with uncertainty measures (SEs and CIs).
#'
#'   Please note that estimating cessation and takeoff growth parameters may not
#'   be possible if there are no distinct pre-peak or post-peak troughs in the
#'   data.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @param newdata An optional data frame for estimation. If \code{NULL}
#'   (default), \code{newdata} is retrieved from the \code{model}.
#'
#' @param resp A character string (default \code{NULL}) to specify the response
#'   variable when processing posterior draws for \code{univariate_by} and
#'   \code{multivariate} models. See [bsitar::bsitar()] for details on
#'   \code{univariate_by} and \code{multivariate} models.
#'
#' @param ndraws A positive integer indicating the number of posterior draws to
#'   use in estimation. If \code{NULL} (default), all draws are used.
#'
#' @param draw_ids An integer specifying the specific posterior draw(s) to use
#'   in estimation (default \code{NULL}).
#'
#' @param summary A logical value indicating whether only the estimate should be
#'   computed (\code{TRUE}), or whether the estimate along with SE and CI should
#'   be returned (\code{FALSE}, default). Setting \code{summary} to \code{FALSE}
#'   will increase computation time. Note that \code{summary = FALSE} is
#'   required to obtain correct estimates when \code{re_formula = NULL}.
#'
#' @param robust A logical value to specify the summary options. If \code{FALSE}
#'   (default), the mean is used as the measure of central tendency and the
#'   standard deviation as the measure of variability. If \code{TRUE}, the
#'   median and median absolute deviation (MAD) are applied instead. Ignored if
#'   \code{summary} is \code{FALSE}.
#'
#' @param transform_draws A function applied to individual draws from the
#'   posterior distribution before computing summaries (default \code{NULL}).
#'   The argument \code{transform_draws} is derived from the
#'   [marginaleffects::predictions()] function and should not be confused with
#'   the \code{transform} argument from the deprecated
#'   [brms::posterior_predict()] function. It's important to note that for both
#'   [marginaleffects::predictions()] and [marginaleffects::avg_predictions()],
#'   the \code{transform_draws} argument takes precedence over the
#'   \code{transform} argument. Note that when \code{transform_draws = NULL}, an
#'   attempt is made to automatically set \code{transform_draws = 'exp'} for
#'   \code{dpar = 'sigma'}. User can set  \code{transform_draws = FALSE} to turn
#'   off this automatic assignment of \code{'exp'} to the
#'   \code{transform_draws}. It is also important to set \code{transform_draws =
#'   FALSE} when computing the first derivative (velocity) for \code{dpar =
#'   'sigma'}.
#'
#' @param re_formula Option to indicate whether or not to include
#'   individual/group-level effects in the estimation. When \code{NA} (default),
#'   individual-level effects are excluded, and population average growth
#'   parameters are computed. When \code{NULL}, individual-level effects are
#'   included in the computation, and the resulting growth parameters are
#'   individual-specific. In both cases (\code{NA} or \code{NULL}), continuous
#'   and factor covariates are appropriately included in the estimation.
#'   Continuous covariates are set to their means by default (see
#'   \code{numeric_cov_at} for details), while factor covariates remain
#'   unaltered, allowing for the estimation of covariate-specific population
#'   average and individual-specific growth parameters.
#'
#' @param peak A logical value (default \code{TRUE}) indicating whether to
#'   calculate the age at peak velocity (APGV) and the peak velocity (PGV)
#'   parameters.
#'
#' @param takeoff A logical value (default \code{FALSE}) indicating whether to
#'   calculate the age at takeoff velocity (ATGV) and the takeoff growth
#'   velocity (TGV) parameters.
#'
#' @param trough A logical value (default \code{FALSE}) indicating whether to
#'   calculate the age at cessation of growth velocity (ACGV) and the cessation
#'   of growth velocity (CGV) parameters.
#'
#' @param acgv A logical value (default \code{FALSE}) indicating whether to
#'   calculate the age at cessation of growth velocity from the velocity curve.
#'   If \code{TRUE}, the age at cessation of growth velocity (ACGV) and the
#'   cessation growth velocity (CGV) are calculated based on the percentage of
#'   the peak growth velocity, as defined by the \code{acgv_velocity} argument
#'   (see below). The \code{acgv_velocity} is typically set at 10 percent of the
#'   peak growth velocity. ACGV and CGV are calculated along with the
#'   uncertainty (SE and CI) around the ACGV and CGV parameters.
#'
#' @param acgv_velocity The percentage of the peak growth velocity to use when
#'   estimating \code{acgv}. The default value is \code{0.10}, i.e., 10 percent
#'   of the peak growth velocity.
#'
#' @param estimation_method A character string specifying the estimation method
#'   when calculating the velocity from the posterior draws. The \code{'fitted'}
#'   method internally calls [bsitar::fitted_draws()], while the
#'   \code{'predict'} method calls [bsitar::predict_draws()]. See
#'   [brms::fitted.brmsfit()] and [brms::predict.brmsfit()] for details.
#'
#' @param numeric_cov_at An optional (named list) argument to specify the value
#'   of continuous covariate(s). The default \code{NULL} option sets the
#'   continuous covariate(s) to their mean. Alternatively, a named list can be
#'   supplied to manually set these values. For example, \code{numeric_cov_at =
#'   list(xx = 2)} will set the continuous covariate variable 'xx' to 2. The
#'   argument \code{numeric_cov_at} is ignored when no continuous covariates are
#'   included in the model.
#'
#' @param levels_id An optional argument to specify the \code{ids} for the
#'   hierarchical model (default \code{NULL}). It is used only when the model is
#'   applied to data with three or more levels of hierarchy. For a two-level
#'   model, \code{levels_id} is automatically inferred from the model fit. For
#'   models with three or more levels, \code{levels_id} is inferred from the
#'   model fit under the assumption that hierarchy is specified from the lowest
#'   to the uppermost level, i.e., \code{id} followed by \code{study}, where
#'   \code{id} is nested within \code{study}. However, it is not guaranteed that
#'   \code{levels_id} is sorted correctly, so it is better to set it manually
#'   when fitting a model with three or more levels of hierarchy.
#'
#' @param avg_reffects An optional argument (default \code{NULL}) to calculate
#'   (marginal/average) curves and growth parameters, such as APGV and PGV. If
#'   specified, it must be a named list indicating the \code{over} (typically a
#'   level 1 predictor, such as age), \code{feby} (fixed effects, typically a
#'   factor variable), and \code{reby} (typically \code{NULL}, indicating that
#'   parameters are integrated over the random effects). For example,
#'   \code{avg_reffects = list(feby = 'study', reby = NULL, over = 'age')}.
#'
#' @param aux_variables An optional argument to specify the variable(s) that can
#'   be passed to the \code{ipts} argument (see below). This is useful when
#'   fitting location-scale models and measurement error models. If
#'   post-processing functions throw an error such as \code{variable 'x' not
#'   found in either 'data' or 'data2'}, consider using \code{aux_variables}.
#'   
#' @param grid_add An optional argument to specify the variable(s) that can be
#'   passed to the [marginaleffects::datagrid()]. This is useful when fitting
#'   location-scale models and measurement error models. If
#'   post-processing functions throw an error such as \code{variable 'x' not
#'   found in either 'data' or 'data2'}, consider using \code{grid_add}.
#'   Note that unlike \code{aux_variables} which are passed to the internal data
#'   functions such as \code{'get.newdata'}, the \code{grid_add} are passed to
#'   the [marginaleffects::datagrid()].
#'
#' @param ipts An integer specifying the number of points for interpolating the
#'   predictor variable (e.g., age) to generate smooth curves for predictions
#'   and plots. This value is used as the \code{length.out} argument for
#'   [seq()], controlling the smoothness of distance and velocity curves without
#'   altering the predictor range.
#'   \describe{
#'     \item{\code{NULL} (the default)}{Engages automatic behavior based on the
#'       \code{dpar} argument: it is internally set to \code{50} for the mean
#'       response (\code{dpar = 'mu'}), ensuring smooth curves, while for the
#'       distributional parameters (e.g., \code{dpar = 'sigma'}), no
#'       interpolation is performed.}
#'     \item{An integer (e.g., \code{100})}{Explicitly sets the number of
#'       interpolation points.}
#'     \item{\code{FALSE}}{Disables all interpolation, forcing predictions to
#'       be made only at the original data points of the predictor variable.}
#'   }
#'   
#'   This argument affects the following post-processing functions: \cr
#'   [fitted_draws()], [predict_draws()], [growthparameters()], [plot_curves()], 
#'   [get_predictions()], [get_comparisons()], and 
#'   [get_growthparameters()].
#'
#' @param model_deriv A logical value specifying whether to estimate the
#'   velocity curve from the derivative function or by differentiating the
#'   distance curve. Set \code{model_deriv = TRUE} for functions that require
#'   the velocity curve, such as \code{growthparameters()} and
#'   \code{plot_curves()}. Set it to \code{NULL} for functions that use the
#'   distance curve (i.e., fitted values), such as \code{loo_validation()} and
#'   \code{plot_ppc()}.
#'
#' @param conf A numeric value (default \code{0.95}) to compute the confidence
#'   interval (CI). Internally, \code{conf} is translated into paired
#'   probability values as \code{c((1 - conf)/2, 1 - (1 - conf)/2)}. For
#'   \code{conf = 0.95}, this computes a 95% CI where the lower and upper limits
#'   are named \code{Q.2.5} and \code{Q.97.5}, respectively.
#'
#' @param xrange An integer to set the predictor range (e.g., age) when
#'   executing the interpolation via \code{ipts}. By default, \code{NULL} sets
#'   the individual-specific predictor range. Setting \code{xrange = 1} applies
#'   the same range for individuals within the same higher grouping variable
#'   (e.g., study). Setting \code{xrange = 2} applies an identical range across
#'   the entire sample. Alternatively, a numeric vector (e.g., \code{xrange =
#'   c(6, 20)}) can be provided to set the range within the specified values.
#'
#' @param xrange_search A vector of length two or a character string
#'   \code{'range'} to set the range of the predictor variable (\code{x}) within
#'   which growth parameters are searched. This is useful when there is more
#'   than one peak and the user wants to summarize the peak within a specified
#'   range of the \code{x} variable. The default value is \code{xrange_search =
#'   NULL}.
#'
#' @param digits An integer (default \code{2}) to set the decimal places for
#'   rounding the results using the [base::round()] function.
#'
#' @param seed An integer (default \code{123}) that is passed to the estimation
#'   method to ensure reproducibility.
#'
#' @param future A logical value (default \code{FALSE}) indicating whether to
#'   perform parallel computations. If \code{TRUE}, posterior summaries are
#'   computed in parallel using [future.apply::future_sapply()].
#' 
#' @param future_session A character string or a named list specifying the
#'   parallel plan when \code{future = TRUE}.
#'   \itemize{
#'     \item \strong{Character string}: Defaults to \code{"multisession"}.
#'       Use \code{"multicore"} for forking (not supported on Windows).
#'     \item \strong{Named list}: Useful for advanced plans like
#'       [future.mirai::mirai_cluster()]. The list must contain:
#'       \itemize{
#'         \item \code{future_session}: The planner function (e.g.,
#'           \code{mirai_cluster}).
#'         \item Additional named arguments passed to the planner (e.g.,
#'           \code{daemons} for [mirai::daemons()]).
#'       }
#'       Example: \code{list(future_session = mirai_cluster, daemons =
#'       list(...))}.
#'   }
#' 
#' @param cores An integer specifying the number of cores for parallel
#'   execution.
#'   \itemize{
#'     \item If \code{NULL} (default), the number of cores is set to 
#'       \code{future::availableCores() - 1}.
#'     \item On non-Windows systems, this can be controlled globally via the 
#'       \code{mc.cores} option.
#'   }
#'
#' @param parms_eval A logical value to specify whether or not to compute growth
#'   parameters on the fly. This is for internal use only and is mainly needed
#'   for compatibility across internal functions.
#'
#' @param idata_method A character string specifying the interpolation method
#'   (default \code{NULL}). The number of interpolation points is controlled by
#'   the \code{ipts} argument.
#'
#'   Available options:
#'   \itemize{
#'     \item \code{'m1'}: Adapted from the \pkg{iapvbs} package (documented
#'     [here](https://rdrr.io/github/Zhiqiangcao/iapvbs/src/R/exdata.R)). This
#'     method internally constructs the data frame based on the model
#'     configuration. \emph{Note:} This method may fail if the model includes
#'     covariates (especially in \code{univariate_by} models).
#'     \item \code{'m2'}: Based on the \pkg{JMbayes} package (documented
#'     [here](https://github.com/drizopoulos/JMbayes/blob/master/R/dynPred_lme.R)).
#'     This method uses the exact data frame from the model fit
#'     (\code{fit$data}) and is generally more robust.
#'   }
#'   
#'   If \code{idata_method = NULL}, \code{'m2'} is automatically selected. 
#'   It is recommended to use \code{'m2'} if \code{'m1'} encounters errors 
#'   with covariate-dependent models.
#'
#' @param parms_method A character string specifying the method used when
#'   evaluating \code{parms_eval}. The default method is \code{getPeak}, which
#'   uses the [sitar::getPeak()] function from the \code{sitar} package.
#'   Alternatively, \code{findpeaks} uses the \code{findpeaks} function from the
#'   \code{pracma} package. This parameter is for internal use and ensures
#'   compatibility across internal functions.
#'
#' @param verbose A logical argument (default \code{FALSE}) to specify whether
#'   to print information collected during the setup of the object(s).
#'
#' @param fullframe A logical value indicating whether to return a
#'   \code{fullframe} object in which \code{newdata} is bound to the summary
#'   estimates. Note that \code{fullframe} cannot be used with \code{summary =
#'   FALSE}, and it is only applicable when \code{idata_method = 'm2'}. A
#'   typical use case is when fitting a \code{univariate_by} model. This option
#'   is mainly for internal use.
#'
#' @param dummy_to_factor A named list (default \code{NULL}) to convert dummy
#'   variables into a factor variable. The list must include the following
#'   elements:
#'   - \code{factor.dummy}: A character vector of dummy variables to be 
#'   converted to factors.
#'   - \code{factor.name}: The name for the newly created factor variable 
#'   (default is \code{'factor.var'} if \code{NULL}).
#'   - \code{factor.level}: A vector specifying the factor levels.
#'   If \code{NULL}, levels are taken from \code{factor.dummy}.
#'   If \code{factor.level} is provided, its length must match
#'   \code{factor.dummy}.
#'
#' @param expose_function A logical value or \code{NULL} (default \code{FALSE}).
#'   Controls whether Stan functions are exposed for post-processing.
#'   \itemize{
#'     \item \code{TRUE}: Explicitly exposes Stan functions saved from the model
#'     fit. This is required when calculating \code{fit criteria} or
#'     \code{bayes_R2} during model optimization.
#'     \item \code{FALSE} (default): Stan functions are not exposed.
#'     \item \code{NULL}: The setting is inherited from the original
#'     \code{bsitar()} model fit configuration.
#'   }
#'   \strong{Note}: In the [optimize_model()] function, the default is
#'   \code{NULL} (inheriting behavior), whereas other post-processing functions
#'   default to \code{FALSE}.
#'
#' @param usesavedfuns A logical value (default \code{NULL}) indicating whether
#'   to use already exposed and saved Stan functions. This is typically set
#'   automatically based on the \code{expose_functions} argument from the
#'   [bsitar()] call. Manual specification of \code{usesavedfuns} is rarely
#'   needed and is intended for internal testing, as improper use can lead to
#'   unreliable estimates.
#'
#' @param clearenvfuns A logical value indicating whether to clear the exposed
#'   Stan functions from the environment (\code{TRUE}) or not (\code{FALSE}). If
#'   \code{NULL}, \code{clearenvfuns} is set based on the value of
#'   \code{usesavedfuns}: \code{TRUE} if \code{usesavedfuns = TRUE}, or
#'   \code{FALSE} if \code{usesavedfuns = FALSE}.
#'
#' @param funlist A list (default \code{NULL}) specifying function names. This
#'   is rarely needed, as required functions are typically retrieved
#'   automatically. A use case for \code{funlist} is when \code{sigma_formula},
#'   \code{sigma_formula_gr}, or \code{sigma_formula_gr_str} use an external
#'   function (e.g., \code{poly(age)}). The \code{funlist} should include
#'   function names defined in the \code{globalenv()}. For functions needing
#'   both distance and velocity curves (e.g., \code{plot_curves(..., opt =
#'   'dv')}), \code{funlist} must include two functions: one for the distance
#'   curve and one for the velocity curve.
#'   
#' @param xvar A character string (default \code{NULL}) specifying the
#'   \code{'x'} variable. Rarely used because \code{xvar} is inferred
#'   internally. A use case is when conflicting variables exist (e.g.,
#'   \code{sigma_formula}) and user wants to set a specific variable as
#'   \code{'x'}.
#'   
#' @param difx A character string (default \code{NULL}) specifying the
#'   \code{'x'} variable that should be used for manual differentiation of the
#'   distance curve. Internally, the \code{xvar} is set as \code{difx} if
#'   specified. The argument \code{difx} is evaluated only when \code{dpar =
#'   'sigma'}, ignored otherwise. Note that argument \code{xvar} itself is rarely
#'   used because \code{xvar} is inferred internally. A use case is when
#'   conflicting variables exist (e.g., \code{sigma_formula}) and user wants to
#'   set a specific variable as \code{'x'}.
#'   
#' @param idvar A character string (default \code{NULL}) specifying the
#'   \code{'id'} variable. Rarely used because \code{idvar} is inferred
#'   internally.
#'
#' @param itransform A character string (default \code{NULL}) indicating the
#'   variables names that are reverse transformed. Options are  \code{c("x",
#'   "y", "sigma")}. The \code{itransform} is primarily used to get the
#'   \code{xvar} variable at original scale i.e., \code{itransform = 'x'}. To
#'   turn of all transformations, use \code{itransform = ""}. when
#'   \code{itransform = NULL}, the appropriate transformation for \code{xvar} is
#'   selected automatically. Note that when no match for \code{xvar} is found in
#'   the \code{data,frame}, the \code{itransform} will be ignored within the 
#'   calling function, \code{'prepare_transformations()'}.
#'   
#' @param newdata_fixed An indicator to specify whether to check data 
#'   format and structure for the user provided \code{newdata}, and apply needed
#'   \code{prepare_data2} and \code{prepare_transformations}
#'   (\code{newdata_fixed = NULL}, default), return user provided \code{newdata}
#'   (\code{newdata = TRUE}) as it is without checking for the data format or
#'   applying \code{prepare_data2} and \code{prepare_transformations}
#'   (\code{newdata_fixed = 0}), check for the data format and if needed,
#'   prepare data format using \code{prepare_data2} (\code{newdata_fixed = 1}),
#'   or apply \code{prepare_transformations} only assuming that data format is
#'   correct (\code{newdata_fixed = 2}). It is strongly recommended that user
#'   either leave the \code{newdata = NULL} and \code{newdata_fixed = NULL} in
#'   which case data used in the model fitting is automatically retrieved and
#'   checked for the required data format and transformations, and if needed,
#'   \code{prepare_data2} and \code{prepare_transformations} are applied
#'   internally. The other flags provided for  \code{newdata_fixed = 0, 1, 2}
#'   are mainly for the internal use during post-processing.
#' 
#' @param envir The environment used for function evaluation. The default is
#'   \code{NULL}, which sets the environment to \code{parent.frame()}. Since
#'   most post-processing functions rely on \pkg{brms}, it is recommended to set
#'   \code{envir = globalenv()} or \code{envir = .GlobalEnv}, especially for
#'   derivatives like velocity curves.
#'
#' @param ... Additional arguments passed to the \code{brms::fitted.brmsfit()}
#'   and \code{brms::predict()} functions.
#'
#' @return A data frame with either five columns (when \code{summary = TRUE}) or
#'   two columns (when \code{summary = FALSE}, assuming \code{re_formual =
#'   NULL}). The first two columns, common to both scenarios, are
#'   \code{'Parameter'} and \code{'Estimate'}, representing the growth parameter
#'   (e.g., APGV, PGV) and its estimate. When \code{summary = TRUE}, three
#'   additional columns are included: \code{'Est.Error'} and two columns
#'   representing the lower and upper bounds of the confidence intervals, named
#'   \code{Q.2.5} and \code{Q.97.5} (for the 95% CI). If \code{re_formual =
#'   NULL}, an additional column with individual identifiers (e.g., \code{id})
#'   is included.
#' 
#' @inherit brms::fitted.brmsfit params
#' 
#' @rdname growthparameters
#' @export
#'
#' @importFrom rlang .data
#' 
#' @inherit brms::prepare_predictions.brmsfit params 
#'
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR Model 
#' 
#' # To avoid mode estimation, which takes time, the Bayesian SITAR model fit 
#' # to the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check if the model fit object 'berkeley_exfit' exists and load it
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Population average age and velocity during the peak growth spurt
#' growthparameters(model, re_formula = NA)
#' 
#' # Population average age and velocity during the take-off and peak 
#' # growth spurt (APGV, PGV, ATGV, TGV)
#' growthparameters(model, re_formula = NA, peak = TRUE, takeoff = TRUE)
#' 
#' # Individual-specific age and velocity during the take-off and peak
#' # growth spurt (APGV, PGV, ATGV, TGV)
#' growthparameters(model, re_formula = NULL, peak = TRUE, takeoff = TRUE)
#' }
#' 
growthparameters.bgmfit <- function(model,
                               newdata = NULL,
                               resp = NULL,
                               dpar = NULL,
                               ndraws = NULL,
                               draw_ids = NULL,
                               summary = FALSE,
                               robust = FALSE,
                               transform_draws = NULL,
                               scale = c("response", "linear"),
                               re_formula = NA,
                               peak = TRUE,
                               takeoff = FALSE,
                               trough = FALSE,
                               acgv = FALSE,
                               acgv_velocity = 0.10,
                               estimation_method = 'fitted',
                               allow_new_levels = FALSE,
                               sample_new_levels = "uncertainty",
                               incl_autocor = TRUE,
                               numeric_cov_at = NULL,
                               levels_id = NULL,
                               avg_reffects = NULL,
                               aux_variables = NULL,
                               grid_add = NULL,
                               ipts = NULL,
                               model_deriv = TRUE,
                               conf = 0.95,
                               xrange = NULL,
                               xrange_search = NULL,
                               digits = 2,
                               seed = 123,
                               future = FALSE,
                               future_session = 'multisession',
                               cores = NULL,
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
  
  set_ipts_null <- FALSE
  
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
  
  if(is.null(ndraws)) {
    ndraws <- brms::ndraws(model)
  }
  
  if(is.null(model_deriv)) {
    model_deriv <- TRUE
  }
  
  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  
  
  # Initiate non formalArgs()
  # xvar <- NULL;
  acgv_asymptote <- NULL;
  apv <- NULL;
  Parameter <- NULL;
  # IDvar <- NULL;
  groupby_fistr <- NULL;
  groupby_fstr <- NULL;
  subindicatorsi <- NULL;
  Estimate <- NULL;
  ':=' <- NULL;
  . <- NULL;
  XXi <- NULL;
  uvarby <- NULL;
  
  msg_nrow_data_1 <- 
    paste0(" Subsetting dataframe for unique length of predictor variable ", 
           'xvar', 
           "\n ", 
           "resulted in a single row of dataframe.",
           "\n ",
           "If this is intentional, then ignore this message.",
           "\n ",
           "The data were subsetted for a unique combination of following:\n",
           'groupby_fstr_xvars',
           "\n "
    )
  
  msg_xvar_all_same <- 
    paste0("The number of unique values for predictor variable ", 
           'xvar', 
           "\n ", 
           "is one i.e., all values are exact same.",
           "\n ",
           "If this is intentional, then ignore this message")
  
  
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
  
  
  dots <- list(...)
  set_get_dv <- FALSE
  if(!is.null(dots$get_dv)) {
    if(dots$get_dv) {
      if(verbose) message("executing 'get_dv'!")
      set_get_dv <- TRUE
    }
  }
  
  if(set_get_dv) {
    ifunx_ <- function(x)x
  }
  
  setxcall_ <- match.call()
  post_processing_checks_args <- list()
  post_processing_checks_args[['model']]    <- model
  post_processing_checks_args[['xcall']]    <- setxcall_
  post_processing_checks_args[['resp']]     <- resp
  post_processing_checks_args[['envir']]    <- envir
  post_processing_checks_args[['verbose']]  <- verbose
  post_processing_checks_args[['check_d0']] <- FALSE
  post_processing_checks_args[['check_d1']] <- TRUE
  post_processing_checks_args[['check_d2']] <- FALSE
  
  oo    <- CustomDoCall(post_processing_checks, post_processing_checks_args)

  if(!is.null(model$xcall)) {
    if(grepl("plot_curves", model$xcall)) {
      xcall <- "plot_curves"
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
  
  model$xcall <- xcall
  
  check_if_package_installed(model, xcall = xcall)
  
  arguments              <- get_args_(as.list(match.call())[-1], xcall)
  arguments$model        <- model
  arguments$usesavedfuns <- usesavedfuns
  
  if(xcall != 'plot_curves') {
    need_velocity_curve <- TRUE
    need_xvar_must      <- TRUE
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
  } # if(xcall != 'plot_curves') {
  
  
  
  
  if(xcall == 'plot_curves') {
    arguments$plot <- TRUE
  } else {
    arguments$plot <- FALSE
  }
 
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', 'Est.Error', probtitles)
  
  get.cores_ <- get.cores(arguments$cores)
  arguments$cores <- setincores <-  get.cores_[['max.cores']] 
  .cores_ps <- get.cores_[['.cores_ps']]
  
  if (future) {
    if(is.null(cores)) stop2c("Please set the number of cores for 'future' by  
                            using the the 'cores' argument, e.g. cores = 4")
    if (arguments$future_session == 'multisession') {
      future::plan('multisession', workers = setincores)
    } else if (arguments$future_session == 'multicore') {
      future::plan('multicore', workers = setincores)
    }
  }
  
  call_posterior_summary <- function(dat) {
    if (!robust) {
      . <- posterior::summarise_draws(dat,
                                      ~ mean(.x, na.rm = T),
                                      ~ sd(.x, na.rm = T),
                                      ~ quantile(.x,
                                                 probs = probs,
                                                 na.rm = T),
                                      .cores = .cores_ps)[-1] %>%
        data.frame() %>% stats::setNames(set_names_)
    } else if (robust) {
      . <- posterior::summarise_draws(
        dat,
        ~ median(.x, na.rm = T),
        ~ mad(.x, na.rm = T),
        ~ quantile(.x, probs = probs, na.rm = T),
        .cores = .cores_ps
      )[-1] %>%
        data.frame() %>% stats::setNames(set_names_)
    }
    as.data.frame(.)
  }
  

  summarise_gp <-
    function(.x,
             probs,
             future,
             cores,
             peak,
             takeoff,
             trough,
             acgv,
             xrange_search,
             summary,
             robust,
             ...) {
      
      Xnames <- names(.x)[grepl("^P._D.", names(.x))]
      .x <- .x %>% data.frame()
      
      if(!is.null(xrange_search)) {
        if(length(xrange_search) == 1) {
          if(is.symbol(xrange_search)) xrange_search <- deparse(xrange_search)
          if(xrange_search == 'range') {
            ullimits <- range(.x[[xvar]])
            .x <- .x %>% dplyr::mutate(XXi := eval(parse(text = xvar))) %>% 
              dplyr::filter(XXi %in% (ullimits[1]:ullimits[2]))
          }
        } else if(length(xrange_search) == 2) {
            ullimits <- xrange_search
            .x <- .x %>% dplyr::mutate(XXi := eval(parse(text = xvar))) %>% 
              dplyr::filter(XXi %in% (ullimits[1]:ullimits[2]))
            } else {
              stop2c("argument xrange_search should be either 
                 'range' or vector of length 2")
          }
      }
     
      
      if (peak) {
        if (future)
          out_1 <-
            future.apply::future_sapply(Xnames, function(x)
              sitar::getPeak(.x[[xvar]], as.numeric(.x[[x]])))
        if (!future)
          out_1 <-
            sapply(Xnames, function(x)
              sitar::getPeak(.x[[xvar]], as.numeric(.x[[x]])))
        out_1 <- t(out_1)
        out_1[,1] <- ifunx_(out_1[,1])
        colnames(out_1) <- c("APGV", "PGV")
      } else if (!peak) {
        out_1 <- NULL
      }
      
      if (takeoff) {
        if (future)
          out_2 <-
            future.apply::future_sapply(Xnames, function(x)
              sitar::getTakeoff(.x[[xvar]], as.numeric(.x[[x]])))
        if (!future)
          out_2 <-
            sapply(Xnames, function(x)
              sitar::getTakeoff(.x[[xvar]], as.numeric(.x[[x]])))
        out_2 <- t(out_2)
        out_2[,1] <- ifunx_(out_2[,1])
        colnames(out_2) <- c("ATGV", "TGV")
      } else if (!takeoff) {
        out_2 <- NULL
      }
      
      if (trough) {
        if (future)
          out_3 <-
            future.apply::future_sapply(Xnames, function(x)
              sitar::getTrough(.x[[xvar]], as.numeric(.x[[x]])))
        if (!future)
          out_3 <-
            sapply(Xnames, function(x)
              sitar::getTrough(.x[[xvar]], as.numeric(.x[[x]])))
        out_3 <- t(out_3)
        out_3[,1] <- ifunx_(out_3[,1])
        colnames(out_3) <- c("ACGV", "CGV")
      } else if (!trough) {
        out_3 <- NULL
      }
      
      if (acgv) {
        if(!is.null(acgv_asymptote)) stop2c("Currently acgv_asymptote not 
                                           supported. Please use acgv_velocity")
        if(!is.null(acgv_asymptote) & !is.null(acgv_velocity) ) {
          stop2c("Specify either acgv_asymptote or acgv_velocity but not both")
        }
        if(is.null(acgv_asymptote) & is.null(acgv_velocity) ) {
          stop2c("Specify either acgv_asymptote or acgv_velocity")
        }
        
        set_get___fun <- function(x) {
          if (!peak) {
            pkkk <- sitar::getPeak(.x[[xvar]], as.numeric(.x[[x]]))
            pkkk[,1] <- ifunx_(pkkk[,1])
          }
          if ( peak | !is.null(apv)) { # if ( peak | apv) {
            pkkk <- out_1
          }
          get__ <- as.numeric(.x[[x]])
          set_x_for_afo <- .x[[xvar]]
          tempbind <- cbind(get__, set_x_for_afo) %>% 
            data.frame() %>% 
            filter(set_x_for_afo > pkkk [1])
          get__ <- tempbind[,1]
          set_x_for_afo <-  tempbind[,2]
          if(length(get__) != 0) {
            get___pct <- max(get__) * acgv_velocity
            target.index <- 
              which(abs(get__ - get___pct) == min(abs(get__ - get___pct)))
            get_asymptote_pct_x <- set_x_for_afo[target.index]
            if(length(get_asymptote_pct_x) > 1) {
              get_asymptote_pct_x <- mean(get_asymptote_pct_x)
            }
            
            if(length(get_asymptote_pct_x) == 0) {
              get_asymptote_pct_x <- NA
            }

            out_3_temp <- c(get_asymptote_pct_x, get___pct)
          } else if(length(get__) == 0) {
            out_3_temp <- c(NA, NA)
          }
          out_3_temp
        }
        if (future)
          out_3 <-
            future.apply::future_sapply(Xnames, function(x)
              set_get___fun(x))
        if (!future)
          out_3 <-
            sapply(Xnames, function(x)
              set_get___fun(x))
         out_3 <- t(out_3)
         out_3[,1] <- ifunx_(out_3[,1])
         colnames(out_3) <- c("ACGV", "CGV")
      } else if (!acgv) {
        out_3 <- NULL
      }
      
      if(set_get_dv) {
        if(peak)    return(out_1)
        if(takeoff) return(out_2)
        if(trough)  return(out_3)
        if(acgv)    return(out_3)
      }
      

      xframe <- out_1
      if (exists('out_2'))
        xframe <- cbind(xframe, out_2)
      if (exists('out_3'))
        xframe <- cbind(xframe, out_3)
      xframe <- xframe %>% as.matrix()
      pnames <- colnames(xframe)[!grepl("^P._D.", colnames(xframe))]
      if (!robust) {
        o_ <- posterior::summarise_draws(
          xframe,
          ~ mean(.x, na.rm = T),
          ~ sd(.x, na.rm = T),
          ~ quantile(.x, probs = probs, na.rm = T),
          .cores = .cores_ps
        )[-1] %>%
          data.frame() %>% stats::setNames(set_names_) %>%
          dplyr::mutate(Parameter = pnames) %>%
          dplyr::relocate(Parameter)
      } else if (robust) {
        o_ <- posterior::summarise_draws(
          xframe,
          ~ median(.x, na.rm = T),
          ~ mad(.x, na.rm = T),
          ~ quantile(.x, probs = probs, na.rm = T),
          .cores = .cores_ps
        )[-1] %>%
          data.frame() %>% stats::setNames(set_names_) %>%
          dplyr::mutate(Parameter = pnames) %>%
          dplyr::relocate(Parameter)
      }
      
      if (summary) {
        . <- o_[, 1:2]
      } else if (!summary) {
        . <-  o_
      }
      .
    }
  
  
  multiNewVar <- function(df, df2, varname){
    df %>% dplyr::mutate(., !!varname := df2[[varname]])
  }
  
  get_growthparameters <-
    function(out_v_,
             newdata,
             groupby_str,
             summary,
             digits,
             ...) {
      dots <- list(...)
      set_get_dv <- FALSE
      if(!is.null(dots$get_dv)) {
        if(dots$get_dv) {
          set_get_dv <- TRUE
        }
      }
      
      if(set_get_dv) {
        if(sum(c(peak, takeoff, trough, acgv), na.rm = TRUE) > 1) {
          stop2c("For 'get_dv', set only one as TRUE: ",
               "'peak', 'takeoff', 'trough', or 'acgv'")
        }
      }
      
      if (!summary) {
        out__ <- out_v_ %>%
          data.frame() %>%
          stats::setNames(paste0('P._D.', names(.))) %>%
          dplyr::mutate(!!xvar := newdata[[xvar]])
        for(i in idvar) {
          out__ <- out__ %>% multiNewVar(df=., df2 = newdata, varname=i)
        } 
      } else if (summary) {
        out__ <- out_v %>% data.frame() %>%
          dplyr::select(1) %>%  data.frame() %>%
          stats::setNames(paste0('P._D.', names(.))) %>%
          dplyr::mutate(!!xvar := newdata[[xvar]])
        for(i in idvar) {
          out__ <- out__ %>% multiNewVar(df=., df2 = newdata, varname=i)
        } 
      }
      
      if (!is.null(groupby_str)) {
        out__ <-
          cbind(out__, newdata %>% 
                  dplyr::select(dplyr::all_of(groupby_str))) %>%
          data.frame()
        out__ <-
          out__ %>% dplyr::group_by(dplyr::across(dplyr::all_of(groupby_str)))
      } else if (is.null(groupby_str)) {
        out__ <- cbind(out__, newdata)
      }

      if (is.null(re_formula)) {
        if (!summary) {
          parameters <- out__ %>%
            dplyr::group_modify(
              ~ summarise_gp(
                .x,
                probs = probs,
                future = future,
                cores = setincores,
                peak = peak,
                takeoff = takeoff,
                trough = trough,
                acgv = acgv,
                xrange_search = xrange_search,
                summary = summary,
                robust = robust, set_get_dv = set_get_dv
              )
            ) 
        } else if (summary) {
          parameters <- out__ %>%
            dplyr::group_modify(
              ~ summarise_gp(
                .x,
                probs = probs,
                future = future,
                cores = setincores,
                peak = peak,
                takeoff = takeoff,
                trough = trough,
                acgv = acgv,
                xrange_search = xrange_search,
                summary = summary,
                robust = robust, set_get_dv = set_get_dv
              )
            ) 
        }
      } else if (!is.null(re_formula)) {
        if (!summary) {
          parameters <- out__ %>%
            dplyr::group_modify(
              ~ summarise_gp(
                .x,
                probs = probs,
                future = future,
                cores = setincores,
                peak = peak,
                takeoff = takeoff,
                trough = trough,
                acgv = acgv,
                xrange_search = xrange_search,
                summary = summary,
                robust = robust, set_get_dv = set_get_dv
              )
            ) 
        } else if (summary) {
          parameters <- out__ %>%
            dplyr::group_modify(
              ~ summarise_gp(
                .x,
                probs = probs,
                future = future,
                cores = setincores,
                peak = peak,
                takeoff = takeoff,
                trough = trough,
                acgv = acgv,
                xrange_search = xrange_search,
                summary = summary,
                robust = robust, set_get_dv = set_get_dv
              )
            ) 
        }
      }
      if(set_get_dv) {
        return(unname(parameters[,1]))
      }
      parameters <- parameters %>% dplyr::ungroup()
      parameters <- parameters %>% 
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                  ~ round(., digits = digits)))
      return(parameters)
    } # end get_growthparameters
  
  
  
  get_avg_over <- function(raw_re, newdata, by, probs, robust) {
    raw_re_c <- c()
    getitEstimate <- getitarray <- NULL
    for (i in 1:(dim(raw_re)[1])) {
      getitEstimate <- raw_re[i,]
      raw_re_c[i] <- cbind(newdata, getitEstimate) %>% data.frame() %>% 
        dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
        dplyr::summarise(getitEstimate = mean(getitEstimate), 
                         .groups = 'drop') %>%
        dplyr::ungroup() %>%
        dplyr::select(dplyr::all_of(getitEstimate)) 
    }
    
    getitarray <- array(unlist(raw_re_c), 
                        dim=c(length(raw_re_c[[1]]), length(raw_re_c)  ))
    getitarray <- t(getitarray)
    getitarray
  }
  
  
  
  if (arguments$plot) {
    arguments <- sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                            arguments = arguments, 
                                            check_formalArgs = NULL,
                                            check_trace_back = NULL,
                                            envir = parent.frame())
    
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('check_fun')) available_d1 <- FALSE
    
    arguments$ipts <- ipts <- check_ipts(ipts = arguments$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)
    
    out_summary <- list()
    if(!is.null(arguments$...)) {
      arguments <- c(arguments, list(arguments$...))
    }
    
    for (argumentsi in names(arguments)) {
      if (length(arguments[[argumentsi]]) != 0) {
        if (argumentsi != "...") {
          arguments[[argumentsi]] <- eval(arguments[[argumentsi]])
        } 
      }
    }
    
    arguments[which(names(arguments) %in% "")] <- NULL
    
    list2env(arguments, envir = parent.env(new.env()))
    
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', 'Est.Error', probtitles)
    
    get.cores_ <- get.cores(arguments$cores)
    arguments$cores <- setincores <-  get.cores_[['max.cores']] 
    .cores_ps <- get.cores_[['.cores_ps']]
    
    if (future) {
      if(is.null(cores)) stop2c("Please set the number of cores for 'future' by  
                            using the the 'cores' argument, e.g. cores = 4")
      if (arguments$future_session == 'multisession') {
        future::plan('multisession', workers = setincores)
      } else if (arguments$future_session == 'multicore') {
        future::plan('multicore', workers = setincores)
      }
    }
    
    list_c <- attr(newdata, 'list_c')
    if(is.null(list_c)) {
      newdata <- get.newdata(model, 
                             newdata = newdata, 
                             xvar = xvar,
                             idvar = idvar,
                             resp = resp, 
                             dpar = dpar,
                             numeric_cov_at = numeric_cov_at,
                             aux_variables = aux_variables,
                             levels_id = levels_id,
                             ipts = ipts,
                             xrange = xrange,
                             idata_method = idata_method,
                             dummy_to_factor = dummy_to_factor,
                             newdata_fixed = newdata_fixed,
                             verbose = verbose)
      list_c <- attr(newdata, 'list_c')
    } else {
      list_c <- list_c
    }
    
  
    
    for (list_ci in names(list_c)) {
      assign(list_ci, list_c[[list_ci]])
    }
    check__ <- c('xvar', 'yvar', 'idvar', 'cov_vars', 'cov_factor_vars', 
                 'cov_numeric_vars', 'groupby_fstr', 'groupby_fistr', 
                 'uvarby', 'subindicatorsi')
    
    for (check___ in check__) {
      if(!exists(check___)) assign(check___, NULL)
    }
    
    newdata___ <- newdata
    
    if(!is.null(avg_reffects)) {
      if (grepl("d", opt, ignore.case = T)) {
        index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
        dist.. <- substr(opt, index_opt, index_opt)
        if (dist.. != "" & grepl("^[[:upper:]]+$", dist..) )
          stop2c("use option 'd' (and not 'D') with avg_reffects" )
      }
      if (grepl("v", opt, ignore.case = T) ) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
        if (velc.. != "" & grepl("^[[:upper:]]+$", velc..) )
          stop2c("use option 'v' (and not 'V') with avg_reffects" )
      }
    }
    
    if(is.null(avg_reffects)) {
      if (grepl("d", opt, ignore.case = T)) {
        index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
        dist.. <- substr(opt, index_opt, index_opt)
      } else if (!grepl("d", opt, ignore.case = T)) {
        dist.. <- ""
      }
      
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
      } else if (!grepl("v", opt, ignore.case = T)) {
        velc.. <- ""
      }
      
      if (apv & velc.. == "") {
        stop2c("You have set apv = TRUE but your opt argument",
             "\n ",
             "contains no 'v' or 'V' option")
      }
      
      if (dist.. != "") {
        if (grepl("^[[:upper:]]+$", dist..)) {
          groupby_str_d <- groupby_fistr
        } else  if (!grepl("^[[:upper:]]+$", dist..)) {
          groupby_str_d <- groupby_fstr
        }
        if (identical(groupby_str_d, character(0)))
          groupby_str_d <- NULL
        groupby_str_d <- groupby_str_d
      } else {
        groupby_str_d <- NULL
      }
      
      
      
      if (velc.. != "") {
        if (grepl("^[[:upper:]]+$", velc..)) {
          groupby_str_v <- groupby_fistr
        } else  if (!grepl("^[[:upper:]]+$", velc..)) {
          groupby_str_v <- groupby_fstr
        }
        if (identical(groupby_str_v, character(0)))
          groupby_str_v <- NULL
      } else {
        groupby_str_v <- NULL
      }
      
      
      
      if (dist.. != "") {
        newdata <- newdata___
        if (grepl("^[[:upper:]]+$", dist..)) {
          arguments$re_formula <- NULL
        } else if (!grepl("^[[:upper:]]+$", dist..)) {
          arguments$re_formula <- NA
          if (!is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(groupby_fstr, xvar)
          } else if (is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(xvar)
          }
     
          if(!is.na(uvarby)) {
            groupby_fstr_xvars_forx <- c(uvarby, groupby_fstr_xvars)
          } else {
            groupby_fstr_xvars_forx <- groupby_fstr_xvars
          }
          
          if(dpar != "sigma") {
            tempidname <- idvar[1] 
            newdata <- newdata %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars_forx))) %>%
              dplyr::slice(1) %>% # droplevels() %>% 
              dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr))) %>%
              dplyr::mutate(!! as.name(tempidname) := dplyr::first( ept(tempidname) )) %>% 
              dplyr::ungroup()
          } # if(dpar != "sigma") {
          
          if(nrow(newdata) == 1) {
            if(verbose) {
              old     <- c('xvar', 'groupby_fstr_xvars_forx')
              new     <- sapply(old, function(x) collapse_comma(ept(x)) )
              message(multi_gsub_exact(msg_nrow_data_1, old, new) )
              old <- new <- NULL
            }
          } else if(nrow(newdata) > 1) {
            if(length(unique(newdata[[xvar]])) == 1)
              if(verbose) {
                old     <- c('xvar')
                new     <- sapply(old, function(x) collapse_comma(ept(x)) )
                message(multi_gsub_exact(msg_xvar_all_same, old, new) )
                old <- new <- NULL
              }
          }
        }
        
        
        arguments$newdata    <- newdata
        arguments$deriv      <- 0
        arguments$probs      <- probs
        arguments$fullframe  <- NULL 
        arguments$itransform <- ""
        if(set_ipts_null) {
          arguments$ipts <- NULL 
        }
        
        if (estimation_method == 'fitted') {
          out_d_ <- CustomDoCall(fitted_draws, arguments)
        } else if (estimation_method == 'predict') {
          out_d_ <- CustomDoCall(predict_draws, arguments)
        }
        
        if(is.null(out_d_)) return(invisible(NULL))
        
        if (!summary) {
          out_d <- call_posterior_summary((out_d_))
        } else if (summary) {
          out_d <- out_d_
        }
        

        if(!is.na(model$model_info$univariate_by$by)) {
          newdata <- newdata %>%
            dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
            droplevels()
        }
        
        out_summary[['distance']] <-  
          cbind(newdata,
                out_d %>% data.frame() %>%
                  dplyr::mutate(curve = 'distance')) %>%
          data.frame()
        
      } else if (dist.. == "") {
        out_summary[['distance']] <- NULL
      }
      
      
      if (velc.. != "") {
        newdata <- newdata___
        if (grepl("^[[:upper:]]+$", velc..)) {
          arguments$re_formula <- NULL
        } else if (!grepl("^[[:upper:]]+$", velc..)) {
          arguments$re_formula <- NA
          if (!is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(groupby_fstr, xvar)
          } else if (is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(xvar)
          }
          
          if(!is.na(uvarby)) {
            groupby_fstr_xvars_forx <- c(uvarby, groupby_fstr_xvars)
          } else {
            groupby_fstr_xvars_forx <- groupby_fstr_xvars
          }
          
          if(dpar != "sigma") {
            tempidname <- idvar[1] 
            newdata <- newdata %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars_forx))) %>%
              dplyr::slice(1) %>% # droplevels() %>% 
              dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr))) %>%
              dplyr::mutate(!! as.name(tempidname) := dplyr::first( ept(tempidname) )) %>% 
              dplyr::ungroup()
          } # if(dpar != "sigma") {
          
          if(nrow(newdata) == 1) {
            if(verbose) {
              old     <- c('xvar', 'groupby_fstr_xvars_forx')
              new     <- sapply(old, function(x) collapse_comma(ept(x)) )
              message(multi_gsub_exact(msg_nrow_data_1, old, new) )
              old <- new <- NULL
            }
          } else if(nrow(newdata) > 1) {
            if(length(unique(newdata[[xvar]])) == 1)
              if(verbose) {
                old     <- c('xvar')
                new     <- sapply(old, function(x) collapse_comma(ept(x)) )
                message(multi_gsub_exact(msg_xvar_all_same, old, new) )
                old <- new <- NULL
              }
          }
        }
        
        arguments$newdata    <- newdata
        arguments$deriv      <- 1
        arguments$probs      <- probs
        arguments$fullframe  <- NULL 
        arguments$itransform <- ""
        if(set_ipts_null) {
          arguments$ipts <- NULL 
        }
        
        if (estimation_method == 'fitted') {
          out_v_ <- CustomDoCall(fitted_draws, arguments)
        } else if (estimation_method == 'predict') {
          out_v_ <- CustomDoCall(predict_draws, arguments)
        }
        
        if(is.null(out_v_)) return(invisible(NULL))
        
        out_v__apv_ <- out_v_
        if (!summary) {
          out_v <- call_posterior_summary((out_v_))
        } else if (summary) {
          out_v <- out_v_
        }
        
        if(!is.na(model$model_info$univariate_by$by)) {
          newdata <- newdata %>%
            dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
            droplevels()
        }
        
        out_summary[['velocity']] <-
          cbind(newdata,
                out_v %>% data.frame() %>%
                  dplyr::mutate(curve = 'velocity')) %>%
          data.frame()
      } else if (velc.. == "") {
        out_summary[['velocity']] <- NULL
      }
      
      
      if (apv | takeoff | trough | acgv) {
        if(!is.na(model$model_info$univariate_by$by)) {
          newdata <- newdata %>%
            dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
            droplevels()
        }
        
        out_v__apv_ <- t(out_v__apv_)
        out_summary[['parameters']] <-
          get_growthparameters(out_v__apv_, newdata, groupby_str_v, summary, 
                               digits, ...)
      }
      out_summary[['groupby_str_d']] <- groupby_str_d
      out_summary[['groupby_str_v']] <- groupby_str_v
      out_summary[['probtitles']] <- probtitles
    } # if(is.null(avg_reffects)) {
    
    
    
    
    if(!is.null(avg_reffects)) {
      if (grepl("d", opt, ignore.case = T)) {
        index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
        dist.. <- substr(opt, index_opt, index_opt)
      } else if (!grepl("d", opt, ignore.case = T)) {
        dist.. <- ""
      }
      
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
      } else if (!grepl("v", opt, ignore.case = T)) {
        velc.. <- ""
      }
      
      if (apv & velc.. == "") {
        stop2c("You have set apv = TRUE but your opt argument",
             "\n ",
             "contains no 'v' or 'V' option")
      }
      
      
      if (dist.. != "") {
        groupby_str_d <- groupby_fstr
        if (identical(groupby_str_d, character(0)))
          groupby_str_d <- NULL
        groupby_str_d <- groupby_str_d
      } else {
        groupby_str_d <- NULL
      }
      
      
      if (velc.. != "") {
        groupby_str_v <- groupby_fstr
        if (identical(groupby_str_v, character(0)))
          groupby_str_v <- NULL
      } else {
        groupby_str_v <- NULL
      }
      
      summary_org <- arguments$summary
      arguments$summary <- FALSE
      arguments$re_formula <- NULL
      
      groupby_str_d <- c(groupby_str_d, avg_reffects[['feby']])
      groupby_str_v <- c(groupby_str_v, avg_reffects[['feby']])
      
      if (dist.. != "") {
        newdata <- newdata___
        if (grepl("^[[:upper:]]+$", dist..)) {
        } else if (!grepl("^[[:upper:]]+$", dist..)) {
          if (!is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(groupby_fstr, xvar)
          } else if (is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(xvar)
          }
        }
        
        arguments$newdata <- newdata
        arguments$deriv <- 0
        arguments$probs <- probs
        arguments$fullframe <- NULL 
        arguments$itransform <- ""
        if(set_ipts_null) {
          arguments$ipts <- NULL 
        }
        
        if (estimation_method == 'fitted') {
          out_d_ <- CustomDoCall(fitted_draws, arguments)
        } else if (estimation_method == 'predict') {
          out_d_ <- CustomDoCall(predict_draws, arguments)
        }
        
        if(is.null(out_d_)) return(invisible(NULL))
        
        arguments$summary <- summary_org
        
        if(!is.na(model$model_info$univariate_by$by)) {
          newdata <- newdata %>%
            dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
            droplevels()
        }
        
        selectby <- avg_reffects[['reby']]
        selectover <- avg_reffects[['over']]
        selectby_over <- c(selectby, selectover)
        out_d_ <- get_avg_over(out_d_, newdata = newdata, by = selectby_over,
                               probs = probs, robust = robust)
        
        out_d <- brms::posterior_summary(out_d_, probs = probs, robust = robust)
  
        
        newdata <- newdata %>%
          dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), 
                          .keep_all = TRUE) %>% 
          dplyr::arrange(!! as.name(selectby_over)) %>% 
          droplevels()
        
  
        
        out_summary[['distance']] <-  
          cbind(newdata,
                out_d %>% data.frame() %>%
                  dplyr::mutate(curve = 'distance')) %>%
          data.frame()
        
      } else if (dist.. == "") {
        out_summary[['distance']] <- NULL
      }
      
      
      summary_org <- arguments$summary
      arguments$summary <- FALSE
      arguments$re_formula <- NULL
      
      if (velc.. != "") {
        newdata <- newdata___
        if (grepl("^[[:upper:]]+$", velc..)) {
        } else if (!grepl("^[[:upper:]]+$", velc..)) {
          if (!is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(groupby_fstr, xvar)
          } else if (is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(xvar)
            
          }
        }
        
        arguments$newdata <- newdata
        arguments$deriv <- 1
        arguments$probs <- probs
        arguments$fullframe <- NULL 
        arguments$itransform <- ""
        if(set_ipts_null) {
          arguments$ipts <- NULL 
        }
        
        if (estimation_method == 'fitted') {
          out_v_ <- CustomDoCall(fitted_draws, arguments)
        } else if (estimation_method == 'predict') {
          out_v_ <- CustomDoCall(predict_draws, arguments)
        }
        
        if(is.null(out_v_)) return(invisible(NULL))
        
        arguments$summary <- summary_org
        
        if(!is.na(model$model_info$univariate_by$by)) {
          newdata <- newdata %>%
            dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
            droplevels()
        }
        
        selectby <- avg_reffects[['reby']]
        selectover <- avg_reffects[['over']]
        selectby_over <- c(selectby, selectover)
        out_v_ <- get_avg_over(out_v_, newdata = newdata, by = selectby_over,
                               probs = probs, robust = robust)
        
        out_v <- brms::posterior_summary(out_v_, probs = probs, robust = robust)
       
        newdata <- newdata %>%
          dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), 
                          .keep_all = TRUE) %>% 
          dplyr::arrange(!! as.name(selectby_over)) %>% 
          droplevels()
        
        out_summary[['velocity']] <-
          cbind(newdata,
                out_v %>% data.frame() %>%
                  dplyr::mutate(curve = 'velocity')) %>%
          data.frame()
      } else if (velc.. == "") {
        out_summary[['velocity']] <- NULL
      }
      
      if (apv) {
        newdata <- newdata %>%
          dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), 
                          .keep_all = TRUE) %>% 
          dplyr::arrange(!! as.name(selectby_over)) %>% 
          droplevels()
        
        out_v_ <- t(out_v_)
        out_summary[['parameters']] <-
          get_growthparameters(out_v_, newdata, groupby_str_v, summary, 
                               digits, ...) # out_v_
      }
      out_summary[['groupby_str_d']] <- groupby_str_d
      out_summary[['groupby_str_v']] <- groupby_str_v
      out_summary[['probtitles']] <- probtitles
    } # if(!is.null(avg_reffects)) {
    return(out_summary)
  } # if (arguments$plot) {
  
 
  if (!arguments$plot) {
    if(!is.null(arguments$deriv)) {
      stop2c("argument 'deriv' is not allowed for 'growthparameters()'")
    }
    
    
    arguments$ipts <- ipts <- set_for_check_ipts(ipts = ipts, nipts = 50, 
                                                 dpar = dpar, verbose = verbose)
    
    
    arguments <- sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                            arguments = arguments, 
                                            check_formalArgs = NULL,
                                            check_trace_back = NULL,
                                            envir = parent.frame())
    
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    arguments$ipts <- ipts <- check_ipts(ipts = arguments$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)
    
    if(set_get_dv) {
      ipts <- NULL
      arguments$summary <- summary <- FALSE
      if(verbose) {
        message("setting 'ipts = NULL' for 'get_dv'")
        message("setting 'summary = FALSE' for 'get_dv'")
      }
      if(!is.null(arguments$avg_reffects)) {
        stop2c("'avg_reffects' must be NULL for get_dv = TRUE")
      }
    } # if (!arguments$plot) {
    
    newdata <- get.newdata(model, 
                           newdata = newdata, 
                           xvar = xvar,
                           idvar = idvar,
                           resp = resp, 
                           dpar = dpar,
                           numeric_cov_at = numeric_cov_at,
                           aux_variables = aux_variables,
                           levels_id = levels_id,
                           ipts = ipts,
                           xrange = xrange,
                           idata_method = idata_method,
                           dummy_to_factor = dummy_to_factor,
                           newdata_fixed = newdata_fixed,
                           verbose = verbose)
    

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
    
    if(is.null(arguments$re_formula)) {
      opt <- 'V'
    }
    if(!is.null(arguments$re_formula)) {
      opt <- 'v'
    }
    
    if(!is.null(avg_reffects)) {
      if (grepl("v", opt, ignore.case = T) ) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
        if (velc.. != "" & grepl("^[[:upper:]]+$", velc..) )
          stop2c("use option 'v' (and not 'V') with avg_reffects" )
      }
    }
    
    if(is.null(avg_reffects)) {
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
      } else if (!grepl("v", opt, ignore.case = T)) {
        velc.. <- ""
      }
      
      if (velc.. != "") {
        if (grepl("^[[:upper:]]+$", velc..)) {
          groupby_str_v <- groupby_fistr
        } else  if (!grepl("^[[:upper:]]+$", velc..)) {
          groupby_str_v <- groupby_fstr
        }
        if (identical(groupby_str_v, character(0)))
          groupby_str_v <- NULL
      } else {
        groupby_str_v <- NULL
      }
      
      if (grepl("^[[:upper:]]+$", velc..)) {
        arguments$re_formula <- NULL
      } else if (!grepl("^[[:upper:]]+$", velc..)) {
        arguments$re_formula <- NA
        if (!is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(groupby_fstr, xvar)
        } else if (is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(xvar)
        }
        
        
        if(!is.na(uvarby)) {
          groupby_fstr_xvars_forx <- c(uvarby, groupby_fstr_xvars)
        } else {
          groupby_fstr_xvars_forx <- groupby_fstr_xvars
        }
        
        if(dpar != "sigma") {
          tempidname <- idvar[1] 
          newdata <- newdata %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars_forx))) %>%
            dplyr::slice(1) %>% # droplevels() %>% 
            dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr))) %>%
            dplyr::mutate(!! as.name(tempidname) := dplyr::first( ept(tempidname) )) %>% 
            dplyr::ungroup()
        } # if(dpar != "sigma") {
        
        if(nrow(newdata) == 1) {
          if(verbose) {
            old     <- c('xvar', 'groupby_fstr_xvars_forx')
            new     <- sapply(old, function(x) collapse_comma(ept(x)) )
            message(multi_gsub_exact(msg_nrow_data_1, old, new) )
            old <- new <- NULL
          }
        } else if(nrow(newdata) > 1) {
          if(length(unique(newdata[[xvar]])) == 1)
            if(verbose) {
              old     <- c('xvar')
              new     <- sapply(old, function(x) collapse_comma(ept(x)) )
              message(multi_gsub_exact(msg_xvar_all_same, old, new) )
              old <- new <- NULL
            }
        }
      } # else if (!grepl("^[[:upper:]]+$", velc..)) {
      
    
      arguments$newdata <- newdata
      arguments$deriv <- 1
      arguments$probs <- probs
      arguments$re_formula_opt <- opt
      arguments$fullframe      <- NULL 
      arguments$itransform     <- ""
      if(set_ipts_null) {
        arguments$ipts <- NULL 
      } 

      if (estimation_method == 'fitted') {
        out_v_ <- CustomDoCall(fitted_draws, arguments)
      } else if (estimation_method == 'predict') {
        out_v_ <- CustomDoCall(predict_draws, arguments)
      }
      
     if(is.null(out_v_)) return(invisible(NULL))

      out_v__apv_ <- out_v_
      if (!summary) {
        out_v <- call_posterior_summary((out_v_))
      } else if (summary) {
        out_v <- out_v_
      }
      
      if(!is.na(model$model_info$univariate_by$by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
          droplevels()
      }
      out_v__apv_ <- t(out_v__apv_)
      parameters <-
        get_growthparameters(out_v__apv_, newdata, groupby_str_v, summary, 
                             digits, ...)
    } # if(is.null(avg_reffects)) {
    
    
    if(!is.null(avg_reffects)) {
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
      } else if (!grepl("v", opt, ignore.case = T)) {
        velc.. <- ""
      }
      
      summary_org <- arguments$summary
      arguments$summary <- FALSE
      arguments$re_formula <- NULL
      
      groupby_str_d <- avg_reffects[['feby']]
      groupby_str_v <- avg_reffects[['feby']]
      
      if (grepl("^[[:upper:]]+$", velc..)) {
      } else if (!grepl("^[[:upper:]]+$", velc..)) {
        if (!is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(groupby_fstr, xvar)
        } else if (is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(xvar)
          
        }
      }
      
      arguments$newdata <- newdata
      arguments$deriv <- 1
      arguments$probs <- probs
      arguments$re_formula_opt <- opt
      arguments$fullframe      <- NULL 
      arguments$itransform     <- "" 
      if(set_ipts_null) {
        arguments$ipts <- NULL 
      }
      
      if (estimation_method == 'fitted') {
        out_v_ <- CustomDoCall(fitted_draws, arguments)
      } else if (estimation_method == 'predict') {
        out_v_ <- CustomDoCall(predict_draws, arguments)
      }
      
      if(is.null(out_v_)) {
        return(invisible(NULL))
      }
      
      arguments$summary <- summary_org

      selectby <- avg_reffects[['reby']]
      selectover <- avg_reffects[['over']]
      selectby_over <- c(selectby, selectover)
      out_v_ <- get_avg_over(out_v_, newdata = newdata, by = selectby_over,
                             probs = probs, robust = robust)
      
      
      out_v <- out_v_
      
      if(!is.na(model$model_info$univariate_by$by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
          droplevels()
      }
      
      newdata <- newdata %>%
        dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), 
                        .keep_all = TRUE) %>% 
        dplyr::arrange(!! as.name(selectby_over)) %>% 
        droplevels()
      
      if(!is.na(model$model_info$univariate_by$by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
          droplevels()
      }
      
      
      newdata <- newdata %>%
        dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), 
                        .keep_all = TRUE) %>% 
        dplyr::arrange(!! as.name(selectby_over)) %>% 
        droplevels()
     
      out_v_ <- t(out_v_)
      parameters <-
        get_growthparameters(out_v_, newdata, groupby_str_v, summary, 
                             digits, ...)
    } # if(!is.null(avg_reffects)) {
    
    newdata_before_itransform <- newdata
    
    itransform_set <- get_itransform_call(itransform = itransform,
                                          model = model, 
                                          newdata = newdata,
                                          dpar = dpar, 
                                          resp = resp,
                                          auto = FALSE,
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
      parameters <- prepare_transformations(data = parameters, model = model,
                                        itransform = itransform_set)
    }
    return(parameters)
  } # if (!arguments$plot) {
  
} # end growthparameters



#' @rdname growthparameters
#' @export
growthparameters <- function(model, ...) {
  UseMethod("growthparameters")
}


