


#' @title Fit Bayesian SITAR Model
#'
#' @description The \strong{bsitar()} function fits the Bayesian version of the
#'   Super Imposition by Translation and Rotation (\emph{SITAR}) model. The
#'   \emph{SITAR} model is a nonlinear mixed-effects model that has been widely
#'   used to summarize growth processes (such as height and weight) from early
#'   childhood through adulthood.
#'
#'   The frequentest version of the \emph{SITAR} model can be fit using the
#'   already available R package, \pkg{sitar} \insertCite{R-sitar}{bsitar}.
#'   However, the \pkg{bsitar} package offers an enhanced Bayesian
#'   implementation that improves modeling capabilities. In addition to the
#'   conventional univariate analysis (i.e., modeling a single outcome),
#'   \pkg{bsitar} also supports:
#'   
#'   \itemize{
#'     \item Univariate-by-subgroup analysis: This allows for simultaneous
#'     modeling of a single outcome across different subgroups defined by a
#'     factor variable (e.g., gender). The advantage is that posterior draws for
#'     each subgroup are part of a single model object, enabling comparisons of
#'     coefficients across groups and testing of various hypotheses.
#'     \item Multivariate analysis: This approach involves simultaneous joint
#'     modeling of two or more outcomes, allowing for more comprehensive growth
#'     modeling.
#'   }
#'   
#'   The Bayesian implementation in \pkg{bsitar} provides a more flexible and
#'   robust framework for growth curve modeling compared to the frequentist
#'   approach.
#'  
#' @details The \emph{SITAR} is a shape-invariant nonlinear mixed-effects growth
#' curve model that fits a population average (i.e., mean) curve to the data and
#' aligns each individual's growth trajectory to the underlying population curve
#' via a set of (typically) three random effects: \code{size}, \code{timing},
#' and \code{intensity}. Additionally, a slope parameter can be included as a
#' random effect to estimate the variability in adult growth rate (see
#' [sitar::sitar()] for details).
#' 
#' The concept of a shape-invariant model (SIM) was first introduced by
#' \insertCite{Lindstrom1995;textual}{bsitar}, and later used by
#' \insertCite{Beath2007;textual}{bsitar} to model infant growth data (birth to
#' 2 years). The current version of the \emph{SITAR} model was developed by
#' \insertCite{Cole2010;textual}{bsitar} and has been extensively used for
#' modeling growth data 
#' \insertCite{@see @nembidzaneUsingSITARMethod2020 and @Sandhu2020}{bsitar}.
#' 
#' The frequentist version of the \emph{SITAR} model can be fit using the
#' already available R package, \pkg{sitar} \insertCite{R-sitar}{bsitar}. The
#' framework of the Bayesian implementation of the \emph{SITAR} model in the
#' \pkg{bsitar} package is similar to the \pkg{sitar} package, with the main
#' difference being that \pkg{sitar} uses [splines::ns()] to construct the
#' B-splines based natural cubic spline design matrix, whereas \pkg{bsitar}
#' implements a different strategy to create natural cubic splines. The
#' \pkg{bsitar} offers three different types of splines: \pkg{nsp}, \pkg{nsk},
#' and \pkg{rcs}. Both \pkg{nsp} and \pkg{nsk} use the B-splines basis to
#' generate the natural cubic spline design matrix as implemented in
#' [splines2::nsp()] and [splines2::nsk()], whereas \pkg{rcs} is based on the
#' truncated power basis approach (see
#' \insertCite{harrell2001regression;textual}{bsitar} and
#' \insertCite{R-Hmisc;textual}{bsitar} for details) to construct the spline
#' design matrix. While all approaches produce the same growth curves, the
#' model-estimated spline coefficients differ from each other. 
#' 
#' A key challenge in spline based regression methods is to select the number
#' and location of knots. The standard approach is to place knots by a regular
#' sequence of quantiles between the outer boundaries. A regression curve can
#' easily be fitted to the sample using a relatively high number of knots. The
#' problem is then over fitting, where a regression model has a good fit to the
#' given sample but does not generalize well to other samples. A low knot count
#' is thus preferred. However, the standard knot selection process can lead to
#' under performance in the sparser regions of the predictor variable,
#' especially when using a low number of knots. It can also lead to over fitting
#' in the denser regions. The \pkg{sitar} package offers an option (see
#' \code{knots_selection} argument) to use a search algorithm that implements a
#' backward method for knot selection that shows reduced prediction error and
#' a lower information criterion (IC) or cross-validation scores compared to the
#' standard knot selection process in simulation experiments. See
#' \insertCite{Arnes2023}{bsitar} for details. The approach suggested by
#' \insertCite{Arnes2023}{bsitar} is implemented via option \code{method = 'bs'}
#' when setting up the \code{knots_selection} argument. We strongly recommend
#' that users carefully consider the placement of knots—whether using the
#' approach suggested here or another method of their choice—to ensure optimal
#' model performance.
#' 
#' Like \pkg{sitar}, the \insertCite{Cole2010}{bsitar}, the \pkg{bsitar} package
#' fits the \emph{SITAR} model with (usually) three random effects: size
#' (parameter \code{a}), timing (parameter \code{b}), and intensity (parameter
#' \code{c}). Additionally, there is a slope parameter (parameter \code{d}) that
#' models the variability in the adult slope of the growth curve (see
#' [sitar::sitar()] for details).
#' 
#' Note that author of the \pkg{sitar} package \insertCite{Cole2010}{bsitar}
#' enforces the inclusion of \code{d} parameter as a random effect only,
#' excluding it from the fixed structure of the model. However, the \pkg{bsitar}
#' package allows inclusion of the \code{d} parameter in both the fixed and/or
#' random effects structures of the \emph{SITAR} model.
#' 
#' For the three-parameter version of the \emph{SITAR} model (default), the
#' fixed effects structure (i.e., population average trajectory) is specified as
#' \code{fixed = 'a+b+c'}, and the random effects structure, capturing the
#' deviation of individual trajectories from the population average curve, is
#' specified as \code{random = 'a+b+c'}.
#' 
#' The \pkg{bsitar} package offers flexibility in model specification. For
#' example:
#' \itemize{
#'   \item A fixed-effect version of the \emph{SITAR} model can be fit by
#'   setting \code{random = ''}. \item The fixed-effect structure can include a
#'   subset of parameters,
#'   such as size and timing (\code{fixed = 'a+b'}) or size and intensity
#'   (\code{fixed = 'a+c'}).
#'   \item For a four-parameter version of the \emph{SITAR} model, parameter
#'   \code{d} is included in the \code{fixed} and/or \code{random} arguments.
#' }
#' 
#' The \pkg{bsitar} package internally depends on the \pkg{brms} 
#' package \insertCite{@see @R-brms; @brms2021}{bsitar}, which fits a wide 
#' range of hierarchical linear and nonlinear regression models, including 
#' multivariate models. The \pkg{brms} package itself depends on **Stan** for 
#' full Bayesian 
#' inference \insertCite{@see @teamStanReferenceManual; @gelman2015}{bsitar}. 
#' Like \pkg{brms}, the \pkg{bsitar} package allows flexible prior 
#' specifications based on user's knowledge of growth processes (e.g., timing 
#' and intensity of growth spurts).
#' 
#' The \pkg{brms} package uses a combination of \code{normal} and
#' \code{student_t} distributions for regression coefficients, group-level
#' random effects, and the distributional parameter (\code{sigma}), while
#' \pkg{rstanarm} uses \code{normal} distributions for regression coefficients
#' and group-level random effects, but sets \code{exponential} for the
#' distributional parameter (\code{sigma}). By default, \pkg{bsitar} uses
#' \code{normal} distributions for all parameters, including regression
#' coefficients, standard deviations of group-level random effects, and the
#' distributional parameter. Additionally, \pkg{bsitar} provides flexibility in
#' choosing scale parameters for location-scale distributions (such as
#' \code{normal} and \code{student_t}).
#' 
#' The \pkg{bsitar} package also allows three types of model specifications:
#' \code{'univariate'}, \code{'univariate_by'}, and \code{'multivariate'}:
#' \itemize{
#'   \item \code{'univariate'} fits a single model to an outcome variable.
#'   \item \code{'univariate_by'} fits two or more sub-models to an outcome
#'   defined by a factor variable (e.g., sex).
#'   \item \code{'multivariate'} fits a joint model to multiple outcomes with
#'   shared random effects.
#' }
#' 
#' The \pkg{bsitar} package offers full flexibility in specifying predictors,
#' degrees of freedom for design matrices, priors, and initial values. The
#' package also allows users to specify options in a user-friendly manner (e.g.,
#' \code{univariate_by = sex} is equivalent to \code{univariate_by = 'sex'}).
#' 
#'@param x Predictor variable (typically age in years). For a \code{univariate}
#'  model, \code{x} is a single variable. For \code{univariate_by} and
#'  \code{multivariate} models, \code{x} can either be the same for all
#'  sub-models, or different for each sub-model. For example, when fitting a
#'  bivariate model, \code{x = list(x1, x2)} specifies that \code{x1} is the
#'  predictor variable for the first sub-model, and \code{x2} is the predictor
#'  for the second sub-model. To use \code{x1} as a common predictor variable
#'  for both sub-models, you can specify \code{x = list(x1)} or simply \code{x =
#'  x1}.
#'
#'@param y Response variable (e.g., repeated height measurements). For
#'  \code{univariate} and \code{univariate_by} models, \code{y} is specified as
#'  a single variable. In the case of a \code{univariate_by} model, the response
#'  vector for each sub-model is created and named internally based on the
#'  factor levels of the variable used to define the sub-model. For example,
#'  specifying \code{univariate_by = sex} creates response vectors \code{Female}
#'  and \code{Male} when \code{Female} is the first level and \code{Male} is the
#'  second level of the \code{sex} variable. In a \code{multivariate} model, the
#'  response variables are provided as a list, such as \code{y = list(y1, y2)},
#'  where \code{y1} is the response variable for the first sub-model, and
#'  \code{y2} is the response for the second sub-model. Note that for the
#'  \code{multivariate} model, the data are not stacked; instead, response
#'  vectors are separate variables in the \code{data} and must be of equal
#'  length.
#'
#' @param id A factor variable uniquely identifying the groups (e.g.,
#'   individuals) in the data frame. For \code{univariate_by} and
#'   \code{multivariate} models, the \code{id} can be the same (typically) for
#'   all sub-models, or different for each sub-model (see the \code{x} argument
#'   for details on setting different arguments for sub-models).
#'
#' @param data A data frame containing variables such as \code{x}, \code{y},
#'   \code{id}, etc.
#'
#' @param df Degrees of freedom for the natural cubic spline design matrix
#'   (default \code{4}). The \code{df} is used internally to construct knots
#'   (quantiles of the \code{x} distribution) for the spline design matrix. For
#'   \code{univariate_by} and \code{multivariate} models, the \code{df} can be
#'   the same across sub-models (e.g., \code{df = 4}) or different for each
#'   sub-model, such as \code{df = list(4, 5)}, where \code{df = 4} applies to
#'   the first sub-model and \code{df = 5} applies to the second sub-model.
#'
#' @param knots A numeric vector specifying the knots for the natural cubic
#'   spline design matrix (default \code{NULL}). Note that you cannot specify
#'   both \code{df} and \code{knots} at the same time, nor can both be
#'   \code{NULL}. In other words, either \code{df} or \code{knots} must be
#'   specified. Like \code{df}, the \code{knots} can be the same for all
#'   sub-models, or different for each sub-model when fitting
#'   \code{univariate_by} and \code{multivariate} models (see \code{df} for
#'   details).
#'   
#' @param knots_selection A named list to control knot optimization during model
#'   fitting. This feature explores a larger candidate set (typically \code{df +
#'   4}) and selects an optimal subset that minimizes the chosen information
#'   criterion. The following elements are available:
#'   \describe{
#'   \item{\code{select}}{Selection strategy for knots and/or degrees of
#'   freedom. Options are \code{"knots"}, \code{"df"}, or \code{"both"}.
#'   \code{"knots"} (fully integrated into the \pkg{bsitar} workflow)
#'   automatically optimizes knot locations for a given \code{df} and passes the
#'   result to internal functions. The \code{"both"} option first determines the
#'   optimal degree of freedom, then selects knot locations using this
#'   \code{df}, which may result in a different (optimized) \code{df} from the
#'   user's initial choice. Users should note this potential change.
#'   Alternatively, \code{select = "df"} can be combined with \code{return =
#'   TRUE} to recover the optimal degree of freedom, which may then be reused as
#'   an argument in [bsitar::bsitar()].}
#'   \item{\code{all_scores}}{Logical (default \code{FALSE}). When \code{TRUE}
#'   and \code{select = "df"}, returns the information criterion value for each
#'   evaluated \code{df}. Ignored for other selection types.}
#'   \item{\code{plot_all_scores}}{Logical (default \code{FALSE}). When
#'   \code{TRUE} and \code{select = "df"} with \code{all_scores = TRUE}, plots
#'   the criterion for each \code{df}. Ignored otherwise.}
#'   \item{\code{nsearch}}{Number of candidate knots in the search space. If
#'   \code{NULL} (default), set automatically to \code{df + 4} to provide
#'   sufficient coverage for optimization.}
#'   \item{\code{criteria}}{Model selection criterion. Options include
#'   \code{"AIC"}, \code{"BIC"}, and cross-validation (\code{"CV"}). Lower
#'   values are preferred. Defaults to \code{"AIC"}. Note: Cross-validation
#'   often requires considerably more computation than \code{"AIC"} or
#'   \code{"BIC"}.}
#'   \item{\code{cvk}}{Number of folds for cross-validation (default 10). Only
#'   used if \code{criteria = "CV"}.}
#'   \item{\code{cviter}}{Number of cross-validation iterations to average
#'   (default 100). Only used if \code{criteria = "CV"}.}
#'   \item{\code{bkrange}}{Logical (default \code{TRUE}). If \code{TRUE}, uses
#'   the full predictor range for boundary knots. Otherwise, boundary knots are
#'   set via quantiles (see [Hmisc::rcspline.eval()]). This option is relevant
#'   only for \code{stype = "rcs"} and overrides any later \code{bkrange} value
#'   globally for \code{knots_selection}.}
#'   \item{\code{fix_bknots}}{Logical (default \code{TRUE}). If \code{TRUE},
#'   fixes boundary knots during internal optimization. If \code{FALSE},
#'   boundary knots are also candidates for optimization. Applies only when
#'   \code{stype = "rcs"}; see notes on \code{bkrange}.}
#'   \item{\code{kspace}}{Knot spacing approach:
#'   \itemize{
#'     \item \code{"un"}: Uniform quantile-based spacing (default)
#'     \item \code{"nu"}: Non-uniform spacing (algorithm-selected). Note: the
#'     actual number of knots may differ from the specified \code{df}.
#'   }}
#'   \item{\code{method}}{Algorithm for knot optimization:
#'   \itemize{
#'     \item \code{"bs"}: built-in systematic workflow for knot optimization
#'     (default)
#'     \item \code{"rs"}: Recursive search based on user-defined function. Not
#'     implemented yet.
#'   }}
#'   \item{\code{when}}{Stage of knot optimization relative to predictor 
#'   centering:
#'   \itemize{
#'     \item \code{"bc"}: Before centering (raw \code{x} values; default)
#'     \item \code{"ac"}: After centering (i.e., \code{x - xoffset})
#'   }}
#'   \item{\code{what}}{Diagnostics or plotting output type. 
#'   Default \code{"knots"}. Available options:
#'   \itemize{
#'     \item \code{"knots"}: Optimized internal and boundary knots
#'     \item \code{"df"}: Selected degree of freedom
#'     \item \code{"plot1"}: Data and curve (original knots)
#'     \item \code{"plot2"}: Data and curve (optimized knots)
#'     \item \code{"plot3"}: Data and curves for both knot sets
#'     \item \code{"plot4"}: Data, curve, and confidence bands (original knots)
#'     \item \code{"plot5"}: Data, curve, and confidence bands (optimized knots)
#'     \item \code{"plot6"}: Data, curve, and confidence bands (both knots)
#'     \item \code{"plot7"}: Residual plot (original knots)
#'     \item \code{"plot8"}: Residual plot (optimized knots)
#'     \item \code{"plot9"}: Residual plots (comparison)
#'   }}
#'   \item{\code{return}}{Logical. If \code{TRUE}, returns the diagnostic plot
#'   object (as specified by \code{what}). Default is \code{FALSE}.}
#'   \item{\code{print}}{Logical. If \code{TRUE}, displays the diagnostic plot
#'   (as specified by \code{what}). Default is \code{FALSE}.}
#'   }
#'   
#'   Note that when both \code{return} and \code{print} are \code{FALSE}, the
#'   optimized knots are passed to the [bsitar::bsitar()] function for model
#'   fitting. If \code{return = TRUE}, the knots are returned immediately, and
#'   model fitting is not performed. When \code{return = FALSE} and \code{print
#'   = TRUE}, the selected knots are displayed via a plot, and the model fitting
#'   proceeds using [bsitar::bsitar()].
#'
#'   Example usage:
#'   \code{
#'   bsitar(
#'     x = age, y = height, id = id, data = berkeley_exdata,
#'     knots_selection = list(
#'       select = "knots", nsearch = NULL, criteria = "AIC",
#'       bkrange = TRUE, fix_bknots = TRUE, method = "bs",
#'       when = "bc", what = "plot3", return = FALSE, print = TRUE
#'     ),
#'     seed = 123)
#'   }
#'   
#' @param fixed A character string specifying the fixed effects structure
#'   (default \code{'a+b+c'}). For \code{univariate_by} and \code{multivariate}
#'   models, you can specify different fixed effect structures for each
#'   sub-model. For example, \code{fixed = list('a+b+c', 'a+b')} implies that
#'   the fixed effects structure for the first sub-model is \code{'a+b+c'}, and
#'   for the second sub-model it is \code{'a+b'}.
#' 
#' @param random A character string specifying the random effects structure
#'   (default \code{'a+b+c'}). The approach to setting the \code{random} is the
#'   same as for the fixed effects structure (see \code{fixed}).
#' 
#' @param xoffset An optional character string or numeric value to set the
#'   origin of the predictor variable, \code{x} (i.e., centering of \code{x}).
#'   Available options include:
#'  - \code{'mean'}: The mean of \code{x} (i.e., \code{mean(x)}),
#'  - \code{'max'}: The maximum value of \code{x} (i.e., \code{max(x)}),
#'  - \code{'min'}: The minimum value of \code{x} (i.e., \code{min(x)}),
#'  - \code{'apv'}: Age at peak velocity estimated from the velocity curve 
#'  derived from a simple linear model fit to the data,
#'  - Any real number (e.g., \code{xoffset = 12}). Note that the value should be
#'  on the original scale of the \code{x} variable because this number will be
#'  automatically transformed based on the \code{xfunxoffset} which is set same
#'  as \code{xfun} unless \code{xfunxoffset} is explicitly set as \code{FALSE}.
#'  See \code{xfunxoffset} for details.
#' The default is \code{xoffset = 'mean'}.
#' 
#' For \code{univariate_by} and \code{multivariate} models, \code{xoffset} can
#' be the same or different for each sub-model (see \code{x} for details on
#' setting different arguments for sub-models). If \code{xoffset} is a numeric
#' value, it will be transformed internally (e.g., \code{log} or \code{sqrt})
#' depending on the \code{xfun} argument. Similarly, when \code{xoffset} is
#' \code{'mean'}, \code{'min'}, or \code{'max'}, these values are calculated
#' after applying the \code{log} or \code{sqrt} transformation to \code{x}.
#' 
#' @param bstart An optional character string or numeric value to set the origin
#'   of the fixed effect parameter \code{b}. The \code{bstart} argument is used
#'   to establish the location parameter for location-scale based priors (such
#'   as \code{normal()}) via the \code{b_prior_beta} argument, and/or the
#'   initial value via the \code{b_init_beta} argument. The available options
#'   for \code{bstart} are:
#'  - \code{'mean'}: The mean of \code{x} (i.e., \code{mean(x)}),
#'  - \code{'min'}: The minimum value of \code{x} (i.e., \code{min(x)}),
#'  - \code{'max'}: The maximum value of \code{x} (i.e., \code{max(x)}),
#'  - \code{'apv'}: Age at peak velocity estimated from the velocity curve 
#'  derived from a simple linear model fit to the data,
#'  - Any real number (e.g., \code{bstart = 12}).
#'  
#' The default is \code{bstart = 'xoffset'} (i.e., the same value as
#' \code{xoffset}). For \code{univariate_by} and \code{multivariate} models,
#' \code{bstart} can be the same for all sub-models (typically), or different
#' for each sub-model (refer to \code{x} for details on setting different
#' arguments for sub-models).
#' 
#' @param cstart An optional character string or numeric value to set the origin
#'   of the fixed effect parameter \code{c}. The \code{cstart} argument is used
#'   to establish the location parameter for location-scale based priors (such
#'   as \code{normal()}) via the \code{c_prior_beta} argument, and/or the
#'   initial value via the \code{c_init_beta} argument. The available options
#'   for \code{cstart} are:
#'  - \code{'pv'}: Peak velocity estimated from the velocity curve derived 
#'  from the simple linear model fit to the data,
#'  - Any real number (e.g., \code{cstart = 1}).
#'  
#' Note that since parameter \code{c} is estimated on the exponential scale, the
#' \code{cstart} should be adjusted accordingly. The default \code{cstart} is
#' \code{'0'} (i.e., \code{cstart = '0'}). For \code{univariate_by} and
#' \code{multivariate} models, \code{cstart} can be the same for all sub-models
#' (typically), or different for each sub-model (refer to \code{x} for details
#' on setting different arguments for sub-models).
#' 
#' @param xfun An optional character string specifying the transformation of the
#'   predictor variable \code{x}. The default value is \code{NULL}, indicating
#'   that no transformation is applied and the model is fit to the data with the
#'   original scale of \code{x}. The available transformation options are:
#'  - \code{'log'}: Logarithmic transformation,
#'  - \code{'sqrt'}: Square root transformation.
#'  
#' For \code{univariate_by} and \code{multivariate} models, the \code{xfun} can
#' be the same for all sub-models (typically), or different for each sub-model
#' (refer to \code{x} for details on setting different arguments for
#' sub-models).
#' 
#' @param yfun An optional character string specifying the transformation of the
#'   response variable \code{y}. The default value is \code{NULL}, indicating
#'   that no transformation is applied and the model is fit to the data with the
#'   original scale of \code{y}. The available transformation options are:
#'  - \code{'log'}: Logarithmic transformation,
#'  - \code{'sqrt'}: Square root transformation.
#'  
#' For \code{univariate_by} and \code{multivariate} models, the \code{yfun} can
#' be the same for all sub-models (typically), or different for each sub-model
#' (refer to \code{x} for details on setting different arguments for
#' sub-models).
#' 
#' @param xfunxoffset Transformation applied to \code{xoffset} for \code{x}
#'   variable (default \code{TRUE}). See \code{xfun} for details. The default
#'   \code{xfunxoffset = TRUE} sets its value to same as \code{xfun}. Users
#'   rarely need to specify \code{xfunxoffset} themselves. One potential
#'   application is when a user wants to turn it off by setting
#'   \code{xfunxoffset = FALSE}. Note that \code{xfunxoffset} is called only
#'   when \code{xoffset} is a numeric values and not data based such as
#'   \code{xoffset = mean}. This is because even when \code{xfunxoffset} is not
#'   \code{TRUE}, \code{xoffset} is still automatically adjusted by \code{xfun},
#'   as it is inferred from the transformed \code{x} variable.
#' 
#' @param bound An optional real number specifying the value by which the span
#'   of the predictor variable \code{x} should be extended (default is
#'   \code{0.04}). This extension can help in modeling edge cases. For more
#'   details, refer to the \pkg{sitar} package documentation. For
#'   \code{univariate_by} and \code{multivariate} models, the \code{bound} can
#'   be the same for all sub-models (typically), or different for each sub-model
#'   (refer to \code{x} for details on setting different arguments for
#'   sub-models).
#' 
#' @param stype A character string or a named list specifying the spline type to
#'   be used. The available options are:
#'  - \code{'rcs'} (default): Constructs the spline design matrix using the truncated
#'   power basis (Harrell's method), implemented in [Hmisc::rcspline.eval()].
#'  - \code{'nsk'}: Implements a B-spline based natural cubic spline method,
#'   similar to [splines2::nsk()].
#'  - \code{'nsp'}: Implements a B-spline based natural cubic spline method,
#'   similar to [splines2::nsp()].
#' 
#' The \code{'rcs'} method uses a truncated power basis, whereas \code{'nsk'}
#' and \code{'nsp'} are B-spline-based methods. Unlike [splines2::nsp()] and
#' [splines2::nsk()], which normalize the spline basis by default, \code{'nsk'}
#' and \code{'nsp'} return the non-normalized version of the spline. If
#' normalization is desired, the user can specify \code{normalize = TRUE} in a
#' list. For example, to use a normalized \code{'nsp'}, one can specify
#' \code{stype = list(type = 'nsp', normalize = TRUE)}.
#' 
#' For more details, see [Hmisc::rcspline.eval()], [splines2::nsk()], and
#' [splines2::nsp()].
#' 
#'@param terms_rhs An optional character string (default \code{NULL}) specifying
#'  terms on the right-hand side of the response variable, but before the
#'  formula tilde sign \code{~}. The \code{terms_rhs} is used when fitting a
#'  measurement error model.
#'
#'  For example, when fitting a model with measurement error in the response
#'  variable, the formula in [brms::brmsformula()] could be specified as
#'  \code{brmsformula(y | mi(sdy) ~ ...)}. In this case, \code{mi(sdy)} is
#'  passed to the formula via \code{terms_rhs = 'mi(sdy)'}.
#'
#'  For a \code{multivariate} model, each outcome can have its own measurement
#'  error variable. For instance, the \code{terms_rhs} can be specified as a
#'  list: \code{terms_rhs = list(mi(sdy1), mi(sdy2))}.
#'
#'  Note that [brms::brmsformula()] does not allow combining \code{mi()} with
#'  the \code{subset()} formulation which used in fitting \code{univariate_by}
#'  models.
#'  
#'@param a_formula Formula for the fixed effect parameter, \code{a} (default
#'  \code{~ 1}). Users can specify different formulas when fitting
#'  \code{univariate_by} and \code{multivariate} models.
#'
#'  For example, \code{a_formula = list(~1, ~1 + cov)} specifies that the
#'  \code{a_formula} for the first sub-model includes only an intercept, while
#'  the second sub-model includes both an intercept and a covariate \code{cov}.
#'  The covariate(s) can be either continuous or factor variables. For factor
#'  covariates, dummy variables are created internally using
#'  [stats::model.matrix()].
#'
#'  The formula can include a combination of continuous and factor variables, as
#'  well as their interactions.
#'
#'@param b_formula Formula for the fixed effect parameter, \code{b} (default
#'  \code{~ 1}). See \code{a_formula} for details on how to specify the formula.
#'  The behavior and structure of \code{b_formula} are similar to
#'  \code{a_formula}.
#'
#'@param c_formula Formula for the fixed effect parameter, \code{c} (default
#'  \code{~ 1}). See \code{a_formula} for details on how to specify the formula.
#'  The behavior and structure of \code{c_formula} are similar to
#'  \code{a_formula}.
#'
#'@param d_formula Formula for the fixed effect parameter, \code{d} (default
#'  \code{~ 1}). See \code{a_formula} for details on how to specify the formula.
#'  The behavior and structure of \code{d_formula} are similar to
#'  \code{a_formula}.
#'  
#'@param s_formula Formula for the fixed effect parameter, \code{s} (default
#'  \code{~ 1}). The \code{s_formula} sets up the spline design matrix.
#'  Typically, covariates are not included in the \code{s_formula} to limit the
#'  population curve to a single curve for the entire data. In fact, the
#'  \pkg{sitar} package does not provide an option to include covariates in the
#'  \code{s_formula}. However, the \pkg{bsitar} package allows the inclusion of
#'  covariates. In such cases, the user must justify the modeling of separate
#'  curves for each category when the covariate is a factor variable.
#'
#' @param a_formula_gr Formula for the random effect parameter, \code{a}
#'   (default \code{~ 1}). Similar to \code{a_formula}, users can specify
#'   different formulas when fitting \code{univariate_by} and
#'   \code{multivariate} models. The formula can include continuous and/or
#'   factor variables, including their interactions as covariates (see
#'   \code{a_formula} for details). In addition to setting up the design matrix
#'   for the random effect parameter \code{a}, users can define the group
#'   identifier and the correlation structure for random effects using the
#'   vertical bar \code{||} notation. For example, to include only an intercept
#'   for the random effects \code{a}, \code{b}, and \code{c}, you can specify:
#' 
#' \code{a_formula_gr = ~1}, \code{b_formula_gr = ~1}, \code{c_formula_gr = ~1}.
#' 
#' To specify the group identifier (e.g., \code{id}) and an unstructured
#' correlation structure, use the vertical bar notation:
#' 
#' \code{a_formula_gr = ~ (1|i|id)} \cr
#' \code{b_formula_gr = ~ (1|i|id)} \cr
#' \code{c_formula_gr = ~ (1|i|id)} \cr
#' 
#' Here, \code{i} within the vertical bars is a placeholder, and a common
#' identifier (e.g., \code{i}) shared across the random effect formulas will
#' model them as unstructured correlated random effects. For more details on
#' this vertical bar approach, please see \code{[brms::brm()]}.
#' 
#' An alternative approach to specify the group identifier and correlation
#' structure is through the \code{group_by} argument. To achieve the same setup
#' as described above with the vertical bar approach, users can define the
#' formula part as:
#' 
#' \code{a_formula_gr = ~1}, \code{b_formula_gr = ~1}, \code{c_formula_gr = ~1}, 
#' 
#' and use \code{group_by} as \code{group_by = list(groupvar = id, cor = un)},
#' \cr where \code{id} specifies the group identifier and \code{un} sets the
#' unstructured correlation structure. See the \code{group_by} argument for more
#' details.
#'
#' @param b_formula_gr Formula for the random effect parameter, \code{b}
#'   (default \code{~ 1}). Similar to \code{a_formula_gr}, user can specify
#'   different formulas when fitting \code{univariate_by} and
#'   \code{multivariate} models. The formula can include continuous and/or
#'   factor variable(s), including their interactions as covariates (see
#'   \code{a_formula_gr} for details). In addition to setting up the design
#'   matrix for the random effect parameter \code{b}, the user can set up the
#'   group identifier and the correlation structure for random effects via the
#'   vertical bar \code{||} approach. For example, consider only an intercept
#'   for the random effects \code{a}, \code{b}, and \code{c} specified as \cr
#'   \code{a_formula_gr = ~1}, \cr
#'   \code{b_formula_gr = ~1}, and \cr 
#'   \code{c_formula_gr = ~1}. \cr
#'   To specify the group identifier (e.g., \code{id}) and an unstructured
#'   correlation structure, the formula argument can be specified as: \cr
#'   \code{a_formula_gr = ~ (1|i|id)} \cr \code{b_formula_gr = ~ (1|i|id)} \cr
#'   \code{c_formula_gr = ~ (1|i|id)} \cr where \code{i} within the vertical
#'   bars \code{||} is just a placeholder. A common identifier (i.e., \code{i})
#'   shared across random effect formulas are modeled as unstructured
#'   correlated. For more details on the vertical bar approach, please see
#'   [brms::brm()].
#' 
#' @param c_formula_gr Formula for the random effect parameter, \code{c}
#'   (default \code{~ 1}). See \code{b_formula_gr} for details.
#'
#' @param d_formula_gr Formula for the random effect parameter, \code{d}
#'   (default \code{~ 1}). See \code{b_formula_gr} for details.
#'  
#' @param a_formula_gr_str Formula for the random effect parameter, \code{a}
#'   (default \code{NULL}), used when fitting a hierarchical model with three or
#'   more levels of hierarchy. For example, a model applied to data that
#'   includes repeated measurements (level 1) on individuals (level 2), which
#'   are further nested within growth studies (level 3). 
#'   
#'   For \code{a_formula_gr_str} argument, only the vertical bar approach (see
#'   \code{a_formula_gr}) can be used to define the group identifiers and
#'   correlation structure. An example of setting up a formula for a three-level
#'   model with random effect parameters \code{a}, \code{b}, and \code{c} is as
#'   follows: \cr 
#'   \code{a_formula_gr_str = ~ (1|i|id:study) + (1|i2|study)} \cr 
#'   \code{b_formula_gr_str = ~ (1|i|id:study) + (1|i2|study)} \cr 
#'   \code{c_formula_gr_str = ~ (1|i|id:study) + (1|i2|study)} \cr 
#'   
#'   In this example, \code{|i|} and \code{|i2|} set up unstructured correlation
#'   structures for the random effects at the individual and study levels,
#'   respectively. Note that \code{|i|} and \code{|i2|} must be distinct, as
#'   random effect parameters cannot be correlated across different levels of
#'   hierarchy.
#'   
#'   Additionally, users can specify models with any number of hierarchical
#'   levels and include covariates in the random effect formula.
#'   
#' @param b_formula_gr_str Formula for the random effect parameter, \code{b} 
#'   (default \code{NULL}), used when fitting a hierarchical model with three 
#'   or more levels of hierarchy. For details, see \code{a_formula_gr_str}.
#'
#' @param c_formula_gr_str Formula for the random effect parameter, \code{c} 
#'   (default \code{NULL}), used when fitting a hierarchical model with three 
#'   or more levels of hierarchy. For details, see \code{a_formula_gr_str}.
#'
#' @param d_formula_gr_str Formula for the random effect parameter, \code{d} 
#'   (default \code{NULL}), used when fitting a hierarchical model with three 
#'   or more levels of hierarchy. For details, see \code{a_formula_gr_str}.
#'   
#' @param d_adjusted A logical indicator to adjust the scale of the predictor
#'   variable \code{x} when fitting the model with the random effect parameter
#'   \code{d}. The coefficient of parameter \code{d} is estimated as a linear
#'   function of \code{x}, i.e., \code{d * x}. If \code{FALSE} (default), the
#'   original \code{x} is used. When \code{d_adjusted = TRUE}, \code{x} is
#'   adjusted for the timing (\code{b}) and intensity (\code{c}) parameters as
#'  \code{x} - \code{b}) * \code{exp(c)} i.e., \code{d * ((x-b)*exp(c))}. The
#'  adjusted scale of \code{x} reflects individual developmental age rather than
#'  chronological age. This makes d more sensitive to the timing of puberty in
#'  individuals. See [sitar::sitar()] function for details.
#'
#' @param sigma_formula Formula for the fixed effect distributional parameter,
#'   \code{sigma}. The \code{sigma_formula} sets up the fixed effect design
#'   matrix, which may include continuous and/or factor variables (and their
#'   interactions) as covariates for the distributional parameter. In other
#'   words, setting up the covariates for \code{sigma_formula} follows the same
#'   approach as for other fixed parameters, such as \code{a} (see
#'   \code{a_formula} for details). Note that \code{sigma_formula} estimates the
#'   \code{sigma} parameter on the \code{log} scale. By default,
#'   \code{sigma_formula} is \code{NULL}, as the [brms::brm()] function itself
#'   models \code{sigma} as a residual standard deviation (\code{RSD}) parameter
#'   on the link scale. The \code{sigma_formula}, along with the arguments
#'   \code{sigma_formula_gr} and \code{sigma_formula_gr_str}, allows
#'   \code{sigma} to be estimated as a random effect. The setup for fixed and
#'   random effects for \code{sigma} is similar to the approach used for other
#'   parameters such as \code{a}, \code{b}, and \code{c}.
#'   
#'   It is important to note that an alternative way to set up the fixed effect
#'   design matrix for the distributional parameter \code{sigma} is to use the
#'   \code{dpar_formula} argument. The advantage of \code{dpar_formula} over
#'   \code{sigma_formula} is that it allows users to specify both linear and
#'   nonlinear formulations using the [brms::lf()] and [brms::nlf()] syntax.
#'   These functions offer more flexibility, such as centering the predictors
#'   and enabling or disabling cell mean centering by excluding the intercept
#'   via \code{0 + } formulation. However, a disadvantage of the
#'   \code{dpar_formula} approach is that random effects cannot be included for
#'   \code{sigma}.
#'   
#'   \code{sigma_formula} and \code{dpar_formula} cannot be specified together.
#'   When either \code{sigma_formula} or \code{dpar_formula} is used, the
#'   default estimation of \code{RSD} by [brms::brm()] is automatically turned
#'   off.
#'   
#'   Users can specify an external function, such as \code{poly},
#'   \code{splines::ns} etc. There are two ways to use external function. \cr 
#'   
#'  \itemize{
#'   \item A conventional approach of using routine function such as \cr
#'   \code{sigma_formula = ~ 1 + splines::ns(age, df = 3)} \cr
#'   \item An external function with a single argument (the predictor), e.g.,
#'   \code{sigma_formula = poly(age)}. Additional arguments are specified
#'   externally. For example, to set the degree of the polynomial to 3, a copy
#'   of the \code{poly} function can be created and  modified as follows: \cr
#'   \code{mypoly = poly; formals(mypoly)[['degree']] <- 3; mypoly(age)}. \cr An
#'   advantage of this approach is that spline coefficient names are small. \cr
#'  }
#'  
#' @param sigma_formula_gr Formula for the random effect parameter, \code{sigma}
#'   (default \code{NULL}). See \code{a_formula_gr} for details. Similar to
#'   \code{sigma_formula}, external functions such as \code{poly},
#'   \code{splines::ns} etc. can be used. For further details, please refer to
#'   the description of \code{sigma_formula}.
#' 
#' @param sigma_formula_gr_str Formula for the random effect parameter,
#'   \code{sigma}, when fitting a hierarchical model with three or more levels
#'   of hierarchy. See \code{a_formula_gr_str} for details. As with
#'   \code{sigma_formula}, external functions can be used. For example, \cr
#'   \code{sigma_formula_gr_str = (1 + splines::ns(age, df = 3) | 11 | gr(id, by
#'   = NULL))} For further details, please refer to the description of
#'   \code{sigma_formula}.
#'  
#' @param sigma_formula_manual A custom formula for advanced modeling of the
#'   distributional parameter \code{sigma} using the [brms::nlf()] and
#'   [brms::lf()] functions. This is an advanced and experimental feature for
#'   specifying complex variance structures. The \code{sigma_formula_manual} can
#'   be can be a character string or a list of formulas.
#' 
#'   \strong{Use Cases:} The primary use cases for \code{sigma_formula_manual}
#'   are:
#'   \enumerate{
#'     \item \strong{Basic variance modelling:} A simple formula for
#'     modelling heteroscedasticity (\code{sigma}).
#'     \item \strong{Advanced variance modelling:} Calculating variance weights
#'       where the variance of residuals depends on a specified co variate,
#'       replicating functionality from the \pkg{nlme} package.
#'     \item \strong{Location-Scale Models:} Applying a full `SITAR`-like model
#'     to the scale (\code{sigma}) parameter in addition to the location
#'     (\code{mu}).
#'   }
#' 
#'  \strong{Basic variance modelling:} 
#'  A simple example of modeling \code{sigma} directly:
#'   \preformatted{
#'   sigma_formula_manual = list(nlf(sigma ~ z) + lf(z ~ 1 + (1 | gr(id))))
#'   }
#'   Note that use can specify any external function in the linear predictor 
#'   part of the formula \code{'lf()'} in order to model nonlinear trajectories.
#'   For example, user can include [splines2::nsk()] in the model as follows:
#'   \preformatted{
#'   sigma_formula_manual = list(nlf(sigma ~ z) + 
#'   lf(z ~ 1 + splines2::nsk(age) + (1 | gr(id))))
#'   }
#'   Automatic Prior Assignment:\cr Priors for these linear predictors are
#'   assigned automatically for both mean and group-level random effects. The
#'   function uses priors specified via arguments like
#'   \code{'sigma_prior_beta'}, \code{'sigma_cov_prior_beta'},
#'   \code{'sigma_prior_sd'}, etc. These are the same arguments otherwise used
#'   for setting priors on parameters defined by \code{'sigma_formula'},
#'   \code{'sigma_formula_gr'}, and \code{'sigma_formula_gr_str'}.
#'   
#' \strong{Advanced variance modelling:} This approach models heteroscedasticity
#' using an explicit variance function. The \pkg{'bsitar'} package provides six
#' different methods for variance modeling, five of which are implemented in the
#' \pkg{'nlme'} package. The variance modeling approach uses a helper function
#' \code{'sigmavarfun'} (short hand form, \code{'vf'}) that sets up the
#' appropriate formulation for variance modeling within the \pkg{'bsitar'}
#' package. 
#' 
#' For multivariate models, the \code{'sigmavarfun'} is appropriately renamed to
#' indicate response specific \code{'sigmavarfun'} by adding \code{'_response'}
#' as prefix. The documentation below provides a detailed overview of the
#' available methods and their usage.
#'
#' \strong{Available Methods} (short hands in parentheses):
#' \itemize{
#'   \item \code{'varpower'} (\code{'vp'}): 
#'   Implements [nlme::varPower()].
#'   \item \code{'varconstpower'} (\code{'cp'}): 
#'   Implements [nlme::varConstPower()].
#'   \item \code{'varexp'} (\code{'ve'}): 
#'   Implements [nlme::varExp()] with \code{'form ~ x'} as follows:\cr
#'   \code{'nlme::varExp(form ~ x)'}
#'   \item \code{'fitted'} (\code{'fi'}): 
#'   Implements [nlme::varExp()] with \code{'form ~ fitted(.)'} as follows:\cr
#'   \code{'nlme::varExp(form ~ fitted(.))'}
#'   \item \code{'residual'} (\code{'re'}): 
#'   Implements [nlme::varExp()] with \code{'form ~ resid(.)'} as follows:\cr
#'   \code{'nlme::varExp(form ~ resid(.))'}
#'   \item \code{'mean'} (\code{'me'}): A sixth method based on an example from
#'   the \pkg{brms} reference manual, which models the square root of the fitted
#'   values as follows:\cr
#'    \code{'nlme::varExp(form ~ sqrt(fitted(.)))'}
#' }
#'
#' Below are examples showing how to use \code{'sigmavarfun'} to specify each of
#' the six variance models. We encourage the use of short hand form, \code{'vf'}
#' to avoid any errors in correctly spelling the full form \code{'sigmavarfun'}
#'
#' \strong{1. varpower:}
#' \preformatted{
#'   nlf(sigma ~ vf(param1, param2, predictor), method = 'vp') +
#'   lf(param1 + param2 ~ 1)
#' }
#'
#' \strong{2. varConstPower:}
#' \preformatted{
#'   nlf(sigma ~ vf(param1, param2, param3, predictor), method = cp') + 
#'   lf(param1 + param2 + param3 ~ 1)
#' }
#'
#' \strong{3. varExp:}
#' \preformatted{
#'   nlf(sigma ~ vf(param1, param2, predictor), method = 've') +
#'   lf(param1 + param2 ~ 1)
#' }
#'
#' \strong{4. fitted:}
#' \preformatted{
#'   nlf(sigma ~ vf(param1, param2, identity()), method = 'fi') +
#'   lf(param1 + param2 ~ 1)
#' }
#' 
#' \strong{5. residual:}
#' \preformatted{
#'   nlf(sigma ~ vf(param1, param2, identity(), resp), method = 're') + 
#'   lf(param1 + param2 ~ 1)
#' }
#' 
#' \strong{6. mean:}
#' \preformatted{
#'   nlf(sigma ~ vf(param1, param2, identity()), method = 'me') +
#'   lf(param1 + param2 ~ 1)
#' }
#' 
#' \strong{Internal Predictor Transformations:}
#' The function applies internal transformations based on the chosen method:
#' \itemize{
#'   \item For \code{'varpower'} and \code{'varConstPower'}, the predictor is
#'     transformed to \code{log(abs(predictor))}.
#'   \item For \code{'varexp'}, the predictor is not transformed.
#'   \item For \code{'fitted'} and \code{'residual'}, \code{identity()} is
#'     internally set to \code{fitted(.)}.
#'   \item For \code{'mean'}, \code{identity()} is internally set to
#'     \code{sqrt(fitted(.))}.
#' }
#' 
#' Note that the default \code{'fitted'} (see above, \code{4. fitted:})
#' internally gets coded as \code{method = 'fittedexp'} (or short hand
#' \code{method = 'fe'}) which sets up the variance form similar to the
#' \code{varExp}. To set up the variance form similar to the \code{varpower},
#' please use \code{method = 'fittedpower'} (or short hand \code{method =
#' 'fp'}). Another form, which models the non zero values for the
#' \code{'fitted'}, can be specified as \code{method = 'fittedz'} (or short hand
#' \code{method = 'fz'}).
#' 
#' Similar to the \code{fitted} (please above), the default \code{method =
#' 'residual'} will set up the variance form similar to the \code{varExp} by
#' internally coding the method argument as \code{method = 'residualexp'} (i.e.,
#' \code{method = 're'}). The \code{power} form can be set as \code{method =
#' 'residualpower'} (or, \code{method = 'rp'}).
#' 
#' Like \code{fitted} and \code{residual} (please above), the default
#' \code{method = 'mean'} will set up the variance form similar to the
#' \code{varExp} by internally coding the method argument as \code{method =
#' 'meanexp'} (i.e., \code{method = 'me'}). The \code{power} form can be set as
#' \code{method = 'meanpower'} (or, \code{method = 'mp'}).
#'
#' \strong{Covariates and Random Effects:}
#' The linear predictor, \code{'lf()'}, can be extended to include covariates
#' and group-level random effects. For example, \code{'lf(param1 + param2 ~ 1)'}
#' can be expanded as follows:
#' \preformatted{
#'   lf(param1 + param2 ~ 1 + covariate + (1 || gr(id, by = groupid)))
#' }
#'
#' \strong{Automatic Prior Assignment:} 
#' Priors for these linear predictors are assigned automatically for both mean
#' and group-level random effects. The function uses priors specified via
#' arguments like \code{'sigma_prior_beta'}, \code{'sigma_cov_prior_beta'},
#' \code{'sigma_prior_sd'}, etc. These are the same arguments otherwise used for
#' setting priors on parameters defined by \code{'sigma_formula'},
#' \code{'sigma_formula_gr'}, and \code{'sigma_formula_gr_str'}.
#'
#' To disable this automatic prior assignment, you can add the argument
#' \code{prior = 'self'} to the `nlf()` function. For example:
#' \preformatted{
#'   nlf(..., method = 'vp', prior = 'self')
#' }
#'
#'  \strong{Note:} If user disables the automatic prior assignment, it is
#'  advised to first get the required prior structure by calling the [bsitar()]
#'  with the \code{'get_prior = TRUE'} argument. The relevant portions of this
#'  structure can then be edited and added back to the model via the
#'  \code{'add_self_priors'} argument in [bsitar()].
#'  
#'  For \code{multivariate} model, the \code{sigma_formula_manual} can be used
#'  set up different form and predictor variables for each outcome separately by
#'  using the list approach as shown below: \cr
#'   \code{sigma_formula_manual = list(
#'   nlf(sigma ~ vf(param1, param2, identity()), method = 'fi') +
#'   lf(param1 + param2 ~ 1) ,
#'   nlf(sigma ~ vf(param1, param2, identity()), method = 'fi') +
#'   lf(param1 + param2 ~ 1)
#'   )}. \cr
#'   Note that formulation name \code{'nlf()'} should be same across all
#'   outcomes. The appropriate response assignment is done internally. Also,
#'   parameter should not contain dots or underscores.
#'
#' \strong{Location-Scale Models:} Another important use case for
#' \code{sigma_formula_manual} is modeling \code{sigma} in a location scale
#' \code{SITAR} model, where the \code{SITAR} formula can be applied to the
#' scale parameter (\code{sigma}) similar to the one used for modelling the
#' location parameter \code{mu}. An example is:
#'   
#'  \code{nlf(sigma ~ ls(x, sigmaa, sigmab, sigmac, sigmas1,
#'  sigmas2, sigmas3, sigmas4), loop = FALSE) +
#'  lf(sigmaa ~ 1+(1 |110| gr(id, by = NULL))+(1 |330| gr(study, by = NULL))) +
#'  lf(sigmab ~ 1+(1 |110| gr(id, by = NULL))+(1 |330| gr(study, by = NULL))) +
#'  lf(sigmac ~ 1+(1 |110| gr(id, by = NULL))+(1 |330| gr(study, by = NULL))) +
#'  lf(sigmas1 + sigmas2 + sigmas3 + sigmas4 ~ 1)}.
#'   
#'  Here, \code{ls}, which is an abbreviation for location-scale, is a
#'  placeholder and gets replaced by the the name of actual function that is
#'  used for modelling the scale \code{sigma} part of the model. For example, if
#'  the name of the \code{mu} function is \code{sigmaSITARFun}, then the name of
#'  the  \code{sigma} function would be \code{sigmaSITARFun} and hence \code{ls}
#'  gets replaced as \code{sigmaSITARFun}. All functions and corresponding names
#'  are created internally.
#'  
#'  The first argument \code{x} of the function is again a placeholder for the
#'  actual predictor variable defined for the \code{sigma} modelling via the
#'  \code{sigmax} argument.In other words, the \code{x} gets replaced by the
#'  \code{sigmax} argument evaluated. For example, when \code{sigmax = age},
#'  then the \code{x} gets replaced by \code{age}. Like \code{ls}, internal 
#'  renaming of \code{x} is handled automatically.
#'  
#'  The growth parameters \code{sigmaa}, \code{sigmab}, \code{sigmac} should
#'  match the input used for the \code{sigmafixed} and \code{sigmarandom}
#'  argument. In other words, either \code{sigmafixed} or \code{sigmarandom}
#'  must include the form for the corresponding growth parameter defined in the
#'  \code{ls} function. For example, when either or both \code{sigmafixed} and
#'  \code{sigmarandom} are defined as \code{a+b+c}, then all three growth
#'  parameters \code{sigmaa}, \code{sigmab}, \code{sigmac} should be part of the
#'  \code{ls} function. However, when only sunset of parameters are specified via
#'  the \code{sigmafixed} and \code{sigmarandom} form, say for example
#'  \code{a+b}, the \code{ls} function should exclude \code{sigmac} parameter.
#'  
#'  Similarly, the number of spline parameters are based on the \code{sigmadf}
#'  argument. For example, when \code{sigmadf = 4}, then four spline parameters
#'  \code{sigmas1, ..., sigmas4} are included, as shown above in the example.
#'  
#'  The other relevant information that is passed to the \code{ls} function can
#'  be set via the \code{sigmaknots}, \code{sigmaxoffset}, \code{sigmaxfun}, and
#'  \code{sigmabound} arguments.
#' 
#'  Note that for the \code{location-scale} model, priors must be set up
#'  manually using the \code{add_self_priors} argument. To see which priors are
#'  required, the user can run the code with \code{get_priors = TRUE}. Also note
#'  that the default initial values for \code{location-scale} model are random.
#'   
#' @param sigmax Predictor for the distributional parameter \code{sigma}. See
#'   \code{x} for details. Ignored if \code{sigma_formula_manual = NULL}.
#'   
#' @param sigmaid A factor variable uniquely identifying the groups (e.g.,
#'   individuals) in the data frame. If \code{NULL} (default), then
#'   \code{sigmaid} is set same as \code{id}. For \code{univariate_by} and
#'   \code{multivariate} models, the \code{sigmaid} can be the same (typically)
#'   for all sub-models, or different for each sub-model (see the \code{id}
#'   argument for details on setting different arguments for sub-models).
#'   Ignored if \code{sigma_formula_manual = NULL}.
#' 
#' @param sigmadf Degree of freedom for the spline function used for
#'   \code{sigma}. See \code{df} for details. Ignored if
#'   \code{sigma_formula_manual = NULL}.
#' 
#' @param sigmaknots Knots for the spline function used for \code{sigma}. See
#'   \code{knots} for details. Ignored if \code{sigma_formula_manual = NULL}.
#'   
#' @param sigmafixed Fixed effect formula for the \code{sigma} structure. See
#'   \code{fixed} for details. Ignored if \code{sigma_formula_manual = NULL}.
#'
#' @param sigmarandom Random effect formula for the \code{sigma} structure. See
#'   \code{random} for details. Ignored if \code{sigma_formula_manual = NULL}.
#'   Currently not used even when \code{sigma_formula_manual} is specified.
#'
#' @param sigmaxoffset Offset for the \code{x} in the \code{sigma} structure.
#'   See \code{xoffset} for details. Ignored if \code{sigma_formula_manual =
#'   NULL}.
#'
#' @param sigmabstart Starting value for the \code{b} parameter in the
#'   \code{sigma} structure. See \code{bstart} for details. Ignored if
#'   \code{sigma_formula_manual = NULL}. Currently not used even when
#'   \code{sigma_formula_manual} is specified.
#'
#' @param sigmacstart Starting value for the \code{c} parameter in the
#'   \code{sigma} structure. See \code{cstart} for details. Ignored if
#'   \code{sigma_formula_manual = NULL}. Currently not used even when
#'   \code{sigma_formula_manual} is specified.
#'
#' @param sigmaxfun Transformation of \code{sigmax} variable for \code{sigma}
#'   structure. See \code{xfun} for details. Ignored if
#'   \code{sigma_formula_manual = NULL}.
#' 
#' @param sigmabound Bounds for the \code{x} in the \code{sigma} structure. See
#'   \code{bound} for details. Ignored if \code{sigma_formula_manual = NULL}.
#'
#' @param sigmaxfunxoffset Transformation applied to \code{sigmaxoffset} for
#'   \code{sigmax} variable (default \code{TRUE}). See \code{sigmaxfun} for
#'   details. The default \code{sigmaxfunxoffset = TRUE} sets its value to same
#'   as \code{sigmaxfun}. Users rarely need to specify \code{sigmaxfunxoffset}
#'   themselves. One potential application is when a user wants to turn it off
#'   by setting \code{sigmaxfunxoffset = FALSE}. Note that
#'   \code{sigmaxfunxoffset} is called only when \code{sigmaxoffset} is a
#'   numeric values and not data based such as \code{sigmaxoffset = mean}. This
#'   is because even when \code{sigmaxfunxoffset} is not \code{TRUE},
#'   \code{sigmaxoffset} is still automatically adjusted by \code{sigmaxfun}, as
#'   it is inferred from the transformed \code{sigmax} variable.
#' 
#' @param dpar_formula Formula for the distributional fixed effect parameter,
#'   \code{sigma} (default \code{NULL}). See \code{sigma_formula} for details.
#'   
#' @param autocor_formula Formula to set up the autocorrelation structure of 
#'   residuals (default \code{NULL}). Allowed autocorrelation structures include:
#'   \itemize{
#'   \item autoregressive moving average (\code{arma}) of order \code{p} and 
#'     \code{q}, specified as \code{autocor_formula = ~arma(p = 1, q = 1)}.
#'   \item autoregressive (\code{ar}) of order \code{p}, specified as 
#'     \code{autocor_formula = ~ar(p = 1)}.
#'   \item moving average (\code{ma}) of order \code{q}, specified as 
#'     \code{autocor_formula = ~ma(q = 1)}.
#'   \item unstructured (\code{unstr}) over time (and individuals), specified as
#'     \code{autocor_formula = ~unstr(time, id)}.
#'   }
#'   See [brms::brm()] for further details on modeling the autocorrelation 
#'   structure of residuals.
#'   
#' @param family Family distribution (default \code{gaussian}) and the link
#'   function (default \code{identity}). See [brms::brm()] for details on
#'   available distributions and link functions, and how to specify them. For
#'   \code{univariate_by} and \code{multivariate} models, the \code{family} can
#'   be the same for all sub-models (e.g., \code{family = gaussian()}) or
#'   different for each sub-model, such as \code{family = list(gaussian(),
#'   student())}, which sets \code{gaussian} distribution for the first
#'   sub-model and \code{student_t} distribution for the second. Note that the
#'   \code{family} argument is ignored if \code{custom_family} is specified
#'   (i.e., if \code{custom_family} is not \code{NULL}).
#' 
#' @param custom_family Specifies custom families (i.e., response distribution).
#'   Default is \code{NULL}. For details, see [brms::custom_family()]. Note that
#'   user-defined Stan functions must be exposed by setting
#'   \code{expose_functions = TRUE}.
#'   
#' @param custom_formula Specifies custom formula. Default is \code{NULL}. For
#'   details, see [brms::brmsformula()]. Currently ignored.
#'
#' @param custom_prior Specifies custom prior Default is \code{NULL}. For
#'   details, see [brms::prior()]. Currently ignored. It is primarily designed
#'   to support setting custom prior for \code{custom_formula}. Note that
#'   \code{custom_prior} is different from the \code{set_self_priors} which
#'   evaluated throughout the call.
#' 
#' @param custom_stanvars Allows the preparation and passing of user-defined
#'   variables to be added to Stan's program blocks (default \code{NULL}). This
#'   is primarily useful when defining a \code{custom_family}. For more details
#'   on specifying \code{stanvars}, see [brms::custom_family()]. Note that
#'   \code{custom_stanvars} are passed directly without conducting any sanity
#'   checks.
#' 
#' @param group_arg Specify arguments for group-level random effects. The
#'   \code{group_arg} should be a named list that may include \code{groupvar},
#'   \code{dist}, \code{cor}, and \code{by} as described below:
#'   \itemize{
#'   \item \code{groupvar} specifies the subject identifier. If \code{groupvar =
#'   NULL} (default), \code{groupvar} is automatically assigned based on the
#'   \code{id} argument. \item \code{dist} specifies the distribution from which
#'   the random effects are drawn (default \code{gaussian}). Currently,
#'   \code{gaussian} is the only available distribution (as per the
#'   [brms::brm()] documentation).
#'   \item \code{by} can be used to estimate a separate variance-covariance
#'   structure (i.e., standard deviation and correlation parameters) for random
#'   effect parameters (default \code{NULL}). If specified, the variable used
#'   for \code{by} must be a factor variable. For example, \code{by = 'sex'}
#'   estimates separate variance-covariance structures for males and females.
#'   \item \code{cor} specifies the covariance (i.e., correlation) structure for
#'   random effect parameters. The default covariance is unstructured (\code{cor
#'   = un}) for all model types (i.e., \code{univariate}, \code{univariate_by},
#'   and \code{multivariate}). The alternative correlation structure available
#'   for \code{univariate} and \code{univariate_by} models is \code{diagonal},
#'   which estimates only the variance parameters (standard deviations), while
#'   setting the covariance (correlation) parameters to zero. For
#'   \emph{multivariate} models, options include \code{un}, \code{diagonal}, and
#'   \code{un_s}. The \code{un} structure models a full unstructured
#'   correlation, meaning that the group-level random effects across response
#'   variables are drawn from a joint multivariate normal distribution with
#'   shared correlation parameters. The \code{cor = diagonal} option estimates
#'   only variance parameters for each sub-model, while setting the correlation
#'   parameters to zero. The \code{cor = un_s} option allows for separate
#'   estimation of unstructured variance-covariance parameters for each response
#'   variable. 
#'   }
#'   
#'   Note that it is not necessary to define all or any of these options
#'   (\code{groupvar}, \code{dist}, \code{cor}, or \code{by}), as they will
#'   automatically be set to their default values if unspecified. Additionally,
#'   only \code{groupvar} from the \code{group_arg} argument is passed to the
#'   \emph{univariate_by} and \emph{multivariate} models, as these models have
#'   their own additional options specified via the \code{univariate_by} and
#'   \code{multivariate} arguments. Lastly, the \code{group_arg} is ignored when
#'   random effects are specified using the vertical bar \code{||} approach (see
#'   \code{a_formula_gr} for details) or when fitting a hierarchical model with
#'   three or more levels of hierarchy (see \code{a_formula_gr_str} for
#'   details).
#'   
#' @param sigma_group_arg Specify arguments for modeling distributional-level
#'   random effects for \code{sigma}. The setup for \code{sigma_group_arg}
#'   follows the same approach as described for group-level random effects (see
#'   \code{group_arg} for details).
#'
#' @param univariate_by Set up the univariate-by-subgroup model fitting (default
#'   \code{NULL}) via a named list with the following elements:
#'   \itemize{
#'   \item \code{by} (optional, character string): Specifies the factor variable
#'   used to define the sub-models (default \code{NA}).
#'   \item \code{cor} (optional, character string): Defines the correlation
#'   structure. Options include \code{un} (default) for a full unstructured
#'   variance-covariance structure and \code{diagonal} for a structure with only
#'   variance parameters (i.e., standard deviations) and no covariance (i.e.,
#'   correlations set to zero).
#'   \item \code{terms} (optional, character string): Specifies the method for
#'   setting up the sub-models. Options are \code{'subset'} (default) and
#'   \code{'weights'}. See \code{brms::`addition-terms`} for more details.
#'   }
#' 
#' @param multivariate Set up the multivariate model fitting (default
#'  \code{NULL}) using a named list with the following elements:
#'  \itemize{
#'  \item \code{mvar} (logical, default \code{FALSE}): Indicates whether to fit
#'  a multivariate model.
#'  \item \code{cor} (optional, character string): Specifies the correlation
#'  structure for group-level random effects. Available options are:
#'  \itemize{
#'  \item \code{"un"} (default): Models a full unstructured correlation, where
#'  group-level random effects across response variables are drawn from a joint
#'  multivariate normal distribution with shared correlation parameters.
#'  \item \code{"diagonal"}: Estimates only the variance parameters for each
#'  sub-model, with the correlation parameters set to zero.
#'  \item \code{"un_s"}: Estimates unstructured variance-covariance parameters
#'  separately for each response variable.
#'  }
#'  \item \code{rescor} (logical, default \code{TRUE}): Indicates whether to
#'  estimate the residual correlation between response variables.
#'  \item \code{rcorr_by} (character string, default \code{NULL}): The variable
#'  by which separate residual correlations between response variables are
#'  estimated. This must be a factor variable present in the \code{data}.
#'  Ignored when \code{rescor = FALSE}.
#'  \item \code{rcorr_gr} (character string, default \code{NULL}): The grouping
#'  variable (typically a level 2 variable like \code{id}) for which residual
#'  correlations between response variables are estimated. This must be a factor
#'  variable present in the \code{data}. Ignored when \code{rescor = FALSE}.
#'  \item \code{rcorr_method} (character string, default \code{NULL}): Specifies
#'  the method for estimating residual correlations: \code{"lkj"} (uses an LKJ
#'  distribution) or \code{"cde"} (uses the Cholesky decomposition of the
#'  covariance matrix with a uniform prior of \code{-1, 1} for each off-diagonal
#'  element). If \code{NULL}, \code{"lkj"} is set as the default. Ignored when
#'  \code{rescor = FALSE}.
#'  \item \code{rcorr_prior} (numeric vector, default \code{NULL}): Used to
#'  specify the LKJ prior for each level of \code{rcorr_by} when
#'  \code{rcorr_method = "lkj"}. The length of \code{rcorr_prior} should be
#'  either one (to use the same prior for all levels) or equal to the number of
#'  levels in \code{rcorr_by}. For \code{rcorr_method = 'cde'}, uniform priors
#'  \code{-1,1} are assigned to all correlation parameters. The argument
#'  \code{rcorr_prior} is ignored when \code{rescor = FALSE}.
#'  }
#'   
#' @param a_prior_beta Specify priors for the fixed effect parameter, \code{a}.
#'   (default \code{normal(ymean, ysd, autoscale = FALSE)}). The following key
#'   points are applicable for all prior specifications. For full details, see
#'   [brms::prior()]:
#'   \itemize{
#'   \item Allowed distributions: \code{normal}, \code{student_t},
#'   \code{cauchy}, \code{lognormal}, \code{uniform}, \code{exponential},
#'   \code{gamma}, and \code{inv_gamma} (inverse gamma).
#'   
#'   \item For each distribution, upper and lower bounds can be set via
#'   \code{lb} and \code{ub} (default \code{NA}).
#'   
#'   \item Location-scale based distributions (such as \code{normal},
#'   \code{student_t}, \code{cauchy}, and \code{lognormal}) have an
#'   \code{autoscale} option (default \code{FALSE}). This option multiplies the
#'   scale parameter by a numeric value. While \pkg{brms} typically uses a
#'   scaling factor of 1.0 or 2.5, the \pkg{bsitar} package allows any real
#'   number to be used (e.g., \code{autoscale = 5.0}).
#'   
#'   \item For location-scale distributions, \code{fxl} (\code{function
#'   location}) and \code{fxs} (\code{function scale}) are available to apply
#'   transformations to the location and scale parameters. For example, setting
#'   \code{normal(2, 5, fxl = 'log', fxs = 'sqrt')} translates to
#'   \code{normal(log(2), sqrt(5))}.
#'
#'   \item \code{fxls} (\code{function location scale}) transforms both location
#'   and scale parameters. The transformation applies when both parameters are
#'   involved, as in the log-transformation for normal priors:
#'   \code{log_location = log(location / sqrt(scale^2 / location^2 + 1))},
#'   \code{log_scale = sqrt(log(scale^2 / location^2 + 1))}. This can be
#'   specified as a character string or a list of functions.
#'   
#'   \item For strictly positive distributions like \code{exponential},
#'   \code{gamma}, and \code{inv_gamma}, the lower bound (\code{lb}) is
#'   automatically set to zero.
#'   
#'   \item For uniform distributions, the option \code{addrange} widens the
#'   prior range symmetrically. For example, \code{uniform(a, b, addrange = 5)}
#'   adjusts the range to \code{uniform(a-5, b+5)}.
#'   
#'   \item For exponential distributions, the rate parameter is evaluated as the
#'   inverse of the specified value. For instance, \code{exponential(10.0)} is
#'   internally treated as \code{exponential(1.0 / 10.0)} =
#'   \code{exponential(0.1)}.
#'   
#'   \item Users do not need to specify each option explicitly, as missing
#'   options will automatically default to their respective values. For example,
#'   \code{a_prior_beta = normal(location = 5, scale = 1)} is equivalent to
#'   \code{a_prior_beta = normal(5, 1)}.
#'   
#'   \item For \code{univariate_by} and \code{multivariate} models, priors can
#'   either be the same for all submodels (e.g., \code{a_prior_beta = normal(5,
#'   1)}) or different for each submodel (e.g., \code{a_prior_beta =
#'   list(normal(5, 1), normal(10, 5))}).
#'   
#'   \item For location-scale distributions, the location parameter can be
#'   specified as the mean (\code{ymean}) or median (\code{ymedian}) of the
#'   response variable, and the scale parameter can be specified as the standard
#'   deviation (\code{ysd}) or median absolute deviation (\code{ymad}).
#'   Alternatively, coefficients from a simple linear model can be used (e.g.,
#'   \code{lm(y ~ age)}).
#'   
#'   Example prior specifications include: 
#'   \code{a_prior_beta = normal(ymean, ysd)}, 
#'   \code{a_prior_beta = normal(ymedian, ymad)}, 
#'   \code{a_prior_beta = normal(lm, ysd)}.
#'   
#'   Note that options such as \code{ymean}, \code{ymedian}, \code{ysd},
#'   \code{ymad}, and \code{lm} are available only for the fixed effect
#'   parameter \code{a}, not for other parameters like \code{b}, \code{c}, or
#'   \code{d}.
#'   }
#'   
#' @param b_prior_beta Specify priors for the fixed effect parameter, \code{b}.
#'   The default prior is \code{normal(0, 2.0, autoscale = FALSE)}. For full
#'   details on prior specification, please refer to \code{a_prior_beta}.
#'   
#'   \itemize{
#'   \item Allowed distributions include \code{normal}, \code{student_t},
#'   \code{cauchy}, \code{lognormal}, \code{uniform}, \code{exponential},
#'   \code{gamma}, and \code{inv_gamma}. \item You can set upper and lower
#'   bounds (\code{lb}, \code{ub}) as needed (default is \code{NA}). \item The
#'   \code{autoscale} option controls scaling of the prior’s scale parameter. By
#'   default, this is set to \code{FALSE}. \item Further customization and
#'   transformations can be applied, similar to the \code{a_prior_beta}
#'   specification. 
#'   }
#'   
#' @param c_prior_beta Specify priors for the fixed effect parameter, \code{c}.
#'   The default prior is \code{normal(0, 1.0, autoscale = FALSE)}. For full
#'   details on prior specification, please refer to \code{a_prior_beta}.
#'   
#'   \itemize{
#'   \item Allowed distributions include \code{normal}, \code{student_t},
#'   \code{cauchy}, \code{lognormal}, \code{uniform}, \code{exponential},
#'   \code{gamma}, and \code{inv_gamma}. \item Upper and lower bounds
#'   (\code{lb}, \code{ub}) can be set as necessary (default is \code{NA}).
#'   \item The \code{autoscale} option is also available for scaling the prior's
#'   scale parameter (default \code{FALSE}). \item Similar to
#'   \code{a_prior_beta}, further transformations or customization can be
#'   applied.
#'   }
#'
#' @param d_prior_beta Specify priors for the fixed effect parameter, \code{d}.
#'   The default prior is \code{normal(0, 1.0, autoscale = FALSE)}. For full
#'   details on prior specification, please refer to \code{a_prior_beta}. Note
#'   that to set the scale of location-scale based priors, the user can set
#'   scale as standard deviation \code{'xsd'} or the median absolute deviation
#'   \code{'xmad'} of the predictor variable \code{'x'}.
#'   
#'   \itemize{
#'   \item Allowed distributions include \code{normal}, \code{student_t},
#'   \code{cauchy}, \code{lognormal}, \code{uniform}, \code{exponential},
#'   \code{gamma}, and \code{inv_gamma}. \item The option to set upper and lower
#'   bounds (\code{lb}, \code{ub}) is available (default is \code{NA}). \item
#'   \code{autoscale} allows scaling of the prior’s scale parameter and is
#'   \code{FALSE} by default. \item For more advanced transformations or
#'   customization, similar to \code{a_prior_beta}, these options are available.
#'   }
#'   
#' @param s_prior_beta Specify priors for the fixed effect parameter, \code{s}
#'   (i.e., spline coefficients). The default prior is \code{normal('lm', 'lm',
#'   autoscale = FALSE)}. The general approach is similar to the one described
#'   for other fixed effect parameters (see \code{a_prior_beta} for details).
#'   Key points to note:
#'   \itemize{
#'   \item When using location-scale based priors with 'lm' (e.g.,
#'   \code{s_prior_beta = normal(lm, 'lm')}), the location parameter is set from
#'   the spline coefficients obtained from the simple linear model fit, and the
#'   scale parameter is based on the standard deviation of the spline design
#'   matrix. The location and scale parameters are typically set to \code{lm}
#'   (default), and \code{autoscale} is set to \code{FALSE}. \item For
#'   location-scale based priors, the option \code{sethp} (logical, default
#'   \code{FALSE}) is available to define hierarchical priors. Setting
#'   \code{sethp = TRUE} alters the prior setup to use hierarchical priors:
#'   \code{s ~ normal(0, 'lm')} becomes \code{s ~ normal(0, 'hp')}, where
#'   \code{'hp'} is defined as \code{hp ~ normal(0, 'lm')}. The scale for the
#'   hierarchical prior is automatically taken from the \code{s} parameter, and
#'   it can also be defined
#'   using the same \code{sethp} option. For example, \code{s_prior_beta =
#'   normal(0, 'lm', sethp = cauchy)} will result in \code{s ~ normal(0, 'lm')},
#'   \code{hp ~ cauchy(0, 'lm')
#'   }.
#'   \item For \code{uniform} priors, you can use the option \code{addrange} to
#'   symmetrically expand the prior range.
#'   }
#'   It has been observed that location-scale based prior distributions (such as
#'   \code{normal}, \code{student_t}, and \code{cauchy}) typically work well for
#'   spline coefficients.
#'
#'  Note that the scale parameter for above \code{s_prior_beta = normal(lm, lm)}
#'  (which is same as \code{s_prior_beta = normal(lm,lm1)}) is derived from the
#'  standard deviation of the outcome and the spline design matrix as \cr
#'  \code{sd(y)/sd(X)} where \code{y} is the outcome and \code{X} design matrix. The other variants of scale parameters are: \cr \code{s_prior_beta =
#'  normal(lm,lm2)}) for which the scale parameter is defined as:
#'  \code{lm_se/sd(X)} where \code{lm_se} is the vector of standard error
#'  obtained from the linear model fit and \code{X} is the design matrix. \cr
#'  \code{s_prior_beta = normal(lm,lm3)}) for which the scale parameter is
#'  defined as \code{lm_se} where \code{lm_se} is the vector of standard error
#'  obtained from the linear model fit.
#'
#' @param a_cov_prior_beta Specify priors for the covariate(s) included in the
#'   fixed effect parameter, \code{a} (default \code{normal(0, 5.0, autoscale =
#'   FALSE)}). The approach for specifying priors is similar to
#'   \code{a_prior_beta}, with a few differences:
#'   \itemize{
#'   \item The options \code{'ymean'}, \code{'ymedian'}, \code{'ysd'}, and
#'   \code{'ymad'} are not allowed for \code{a_cov_prior_beta}. \item The
#'   \code{'lm'} option for the location parameter allows the covariate
#'   coefficient(s) to be obtained from a simple linear model fit to the data.
#'   Note that the \code{'lm'} option is only allowed for
#'   \code{a_cov_prior_beta} and not for covariates in other fixed or random
#'   effect parameters. \item Separate priors can be specified for submodels
#'   when fitting \code{univariate_by} and \code{a_prior_beta} models (see
#'   \code{a_prior_beta} for details).
#'   }
#'   
#' @param b_cov_prior_beta Specify priors for the covariate(s) included in the
#'   fixed effect parameter, \code{b} (default \code{normal(0, 1.0, autoscale =
#'   FALSE)}). See \code{a_cov_prior_beta} for details.
#'   
#' @param c_cov_prior_beta Specify priors for the covariate(s) included in the
#'   fixed effect parameter, \code{c} (default \code{normal(0, 0.1, autoscale =
#'   FALSE)}). See \code{a_cov_prior_beta} for details.
#'   
#' @param d_cov_prior_beta Specify priors for the covariate(s) included in the
#'   fixed effect parameter, \code{d} (default \code{normal(0, 1.0, autoscale =
#'   FALSE)}). See \code{a_cov_prior_beta} for details.
#'   
#' @param s_cov_prior_beta Specify priors for the covariate(s) included in the
#'   fixed effect parameter, \code{s} (default \code{normal(0, 10.0, autoscale =
#'   FALSE)}). As described in \code{s_formula}, the \emph{SITAR} model does not
#'   allow covariates in the spline design matrix. If covariates are specified
#'   (see \code{s_formula}), the approach to setting priors for the covariates
#'   in parameter \code{s} is the same as for \code{a} (see
#'   \code{a_cov_prior_beta}). For location-scale based priors, the option
#'   \code{'lm'} sets the location parameter based on spline coefficients
#'   obtained from fitting a simple linear model to the data.
#'   
#' @param a_prior_sd Specify priors for the random effect parameter, \code{a}.
#'   (default \code{normal(0, 'ysd', autoscale = FALSE)}). The prior is applied
#'   to the standard deviation (the square root of the variance), not the
#'   variance itself. The approach for setting the prior is similar to
#'   \code{a_prior_beta}, with the location parameter always set to zero. The
#'   lower bound is automatically set to \code{0} by \code{brms::brm()}. For
#'   \code{univariate_by} and \code{multivariate} models, priors can be the same
#'   or different for each submodel (see \code{a_prior_beta}).
#'
#' @param b_prior_sd Specify priors for the random effect parameter, \code{b}.
#'   (default \code{normal(0, 2.0, autoscale = FALSE)}). See \code{a_prior_sd}
#'   for details.
#' 
#' @param c_prior_sd Specify priors for the random effect parameter, \code{c}.
#'   (default \code{normal(0, 1.0, autoscale = FALSE)}). See \code{a_prior_sd}
#'   for details.
#' 
#' @param d_prior_sd Specify priors for the random effect parameter, \code{d}.
#'   (default \code{normal(0, 1.0, autoscale = FALSE)}). See \code{a_prior_sd}
#'   for details.
#' 
#' @param a_cov_prior_sd Specify priors for the covariate(s) included in the
#'   random effect parameter, \code{a}. (default \code{normal(0, 5.0, autoscale
#'   = FALSE)}). The approach is the same as described for
#'   \code{a_cov_prior_beta}, except that no pre-defined options (e.g.,
#'   \code{'lm'}) are allowed.
#'
#' @param b_cov_prior_sd Specify priors for the covariate(s) included in the
#'   random effect parameter, \code{b}. (default \code{normal(0, 1.0, autoscale
#'   = FALSE)}). See \code{a_cov_prior_sd} for details.
#'
#' @param c_cov_prior_sd Specify priors for the covariate(s) included in the
#'   random effect parameter, \code{c}. (default \code{normal(0, 0.1, autoscale
#'   = FALSE)}). See \code{a_cov_prior_sd} for details.
#'
#' @param d_cov_prior_sd Specify priors for the covariate(s) included in the
#'   random effect parameter, \code{d}. (default \code{normal(0, 1.0, autoscale
#'   = FALSE)}). See \code{a_cov_prior_sd} for details.
#'
#' @param a_prior_sd_str Specify priors for the random effect parameter,
#'   \code{a}, when fitting a hierarchical model with three or more levels of
#'   hierarchy. (default \code{NULL}). The approach is the same as described for
#'   \code{a_prior_sd}.
#'
#' @param b_prior_sd_str Specify priors for the random effect parameter,
#'   \code{b}, when fitting a hierarchical model with three or more levels of
#'   hierarchy. (default \code{NULL}). The approach is the same as described for
#'   \code{a_prior_sd_str}.
#'
#' @param c_prior_sd_str Specify priors for the random effect parameter,
#'   \code{c}, when fitting a hierarchical model with three or more levels of
#'   hierarchy. (default \code{NULL}). The approach is the same as described for
#'   \code{a_prior_sd_str}.
#'
#' @param d_prior_sd_str Specify priors for the random effect parameter,
#'   \code{d}, when fitting a hierarchical model with three or more levels of
#'   hierarchy. (default \code{NULL}). The approach is the same as described for
#'   \code{a_prior_sd_str}.
#'
#' @param a_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'   random effect parameter, \code{a}, when fitting a hierarchical model with
#'   three or more levels of hierarchy. (default \code{NULL}). The approach is
#'   the same as described for \code{a_cov_prior_sd}.
#'
#' @param b_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'   random effect parameter, \code{b}, when fitting a hierarchical model with
#'   three or more levels of hierarchy. (default \code{NULL}). The approach is
#'   the same as described for \code{a_cov_prior_sd_str}.
#'
#' @param c_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'   random effect parameter, \code{c}, when fitting a hierarchical model with
#'   three or more levels of hierarchy. (default \code{NULL}). The approach is
#'   the same as described for \code{a_cov_prior_sd_str}.
#'
#' @param d_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'   random effect parameter, \code{d}, when fitting a hierarchical model with
#'   three or more levels of hierarchy. (default \code{NULL}). The approach is
#'   the same as described for \code{a_cov_prior_sd_str}.
#' 
#' @param sigma_prior_beta Specify priors for the fixed effect distributional
#'   parameter, \code{sigma}. (default \code{normal(0, 1.0, autoscale =
#'   FALSE)}). The approach is similar to that for \code{a_prior_beta}.
#'
#' @param sigma_cov_prior_beta Specify priors for the covariate(s) included in
#'   the fixed effect distributional parameter, \code{sigma}. (default
#'   \code{normal(0, 0.5, autoscale = FALSE)}). Follows the same approach as
#'   \code{a_cov_prior_beta}.
#'
#' @param sigma_prior_sd Specify priors for the random effect distributional
#'   parameter, \code{sigma}. (default \code{normal(0, 0.25, autoscale =
#'   FALSE)}). Same approach as \code{a_prior_sd}.
#'
#' @param sigma_cov_prior_sd Specify priors for the covariate(s) included in the
#'   random effect distributional parameter, \code{sigma}. (default
#'   \code{normal(0, 0.15, autoscale = FALSE)}). Follows the same approach as
#'   \code{a_cov_prior_sd}.
#'
#' @param sigma_prior_sd_str Specify priors for the random effect distributional
#'   parameter, \code{sigma}, when fitting a hierarchical model with three or
#'   more levels of hierarchy. (default \code{NULL}). Same approach as
#'   \code{a_prior_sd_str}.
#'
#' @param sigma_cov_prior_sd_str Specify priors for the covariate(s) included in
#'   the random effect distributional parameter, \code{sigma}, when fitting a
#'   hierarchical model with three or more levels of hierarchy. (default
#'   \code{NULL}). Follows the same approach as \code{a_cov_prior_sd_str}.
#'
#' @param rsd_prior_sigma Specify priors for the residual standard deviation
#'   parameter \code{sigma} (default \code{normal(0, 'ysd', autoscale =
#'   FALSE)}). Evaluated when both \code{dpar_formula} and \code{sigma_formula}
#'   are \code{NULL}. For location-scale based distributions, user can specify
#'   standard deviation (\code{ysd}) or the median absolute deviation
#'   (\code{ymad}) of outcome as the scale parameter. Also, residual standard
#'   deviation from the linear mixed model (\code{nlme::lme()}) or the linear
#'   model (\code{base::lm()}) fitted to the data. These are specified as
#'   \code{'lme_rsd'} and \code{'lm_rsd'}, respectively. Note that if
#'   \code{nlme::lme()} fails to converge, the option \code{'lm_rsd'} is set
#'   automatically. The argument \code{rsd_prior_sigma} is evaluated when both
#'   \code{dpar_formula} and \code{sigma_formula} are set to \code{NULL}.
#'
#' @param dpar_prior_sigma Specify priors for the fixed effect distributional
#'   parameter \code{sigma} (default \code{normal(0, 'ysd', autoscale =
#'   FALSE)}). Evaluated when \code{sigma_formula} is \code{NULL}. See
#'   \code{rsd_prior_sigma} for details.
#'
#' @param dpar_cov_prior_sigma Specify priors for the covariate(s) included in
#'   the fixed effect distributional parameter \code{sigma}. (default
#'   \code{normal(0, 1.0, autoscale = FALSE)}). Evaluated when
#'   \code{sigma_formula} is \code{NULL}.
#' 
#' @param autocor_prior_acor Specify priors for the autocorrelation parameters 
#'   when fitting a model with \code{'arma'}, \code{'ar'}, or \code{'ma'} 
#'   autocorrelation structures (see \code{autocor_formula}). The only allowed 
#'   distribution is \code{uniform}, bounded between -1 and +1 (default 
#'   \code{uniform(-1, 1, autoscale = FALSE)}). For the unstructured residual 
#'   correlation structure, use \code{autocor_prior_unstr_acor}.
#' 
#' @param autocor_prior_unstr_acor Specify priors for the autocorrelation 
#'   parameters when fitting a model with the unstructured (\code{'un'}) 
#'   autocorrelation structure (see \code{autocor_formula}). The only allowed 
#'   distribution is \code{lkj} (default \code{lkj(1)}). See \code{gr_prior_cor} 
#'   for details on setting up the \code{lkj} prior.
#'
#' @param gr_prior_cor Specify priors for the correlation parameter(s) of 
#'   group-level random effects (default \code{lkj(1)}). The only allowed 
#'   distribution is \code{lkj}, specified via a single parameter \code{eta} 
#'   (see \code{brms::prior()} for details).
#'
#' @param gr_prior_cor_str Specify priors for the correlation parameter(s) of
#'   group-level random effects when fitting a hierarchical model with three or
#'   more levels of hierarchy (default \code{lkj(1)}). Same as
#'   \code{gr_prior_cor}.
#'
#' @param sigma_prior_cor Specify priors for the correlation parameter(s) of
#'   distributional random effects \code{sigma} (default \code{lkj(1)}). The
#'   only allowed distribution is \code{lkj} (see \code{gr_prior_cor} for
#'   details). Note that \code{brms::brm()} does not currently allow different
#'   \code{lkj} priors for the group level and distributional random effects
#'   sharing the same group identifier (\code{id}).
#'
#' @param sigma_prior_cor_str Specify priors for the correlation parameter(s) of
#'   distributional random effects \code{sigma} when fitting a hierarchical
#'   model with three or more levels of hierarchy (default \code{lkj(1)}). Same
#'   as \code{sigma_prior_cor}.
#'
#' @param mvr_prior_rescor Specify priors for the residual correlation parameter
#'   when fitting a multivariate model (default \code{lkj(1)}). The only allowed
#'   distribution is \code{lkj} (see \code{gr_prior_cor} for details).
#' 
#' @param init Initial values for the sampler. Options include:
#'  \itemize{
#'    \item \code{'random'} (default): \strong{Stan} randomly generates initial
#'    values for each parameter within a range defined by \code{init_r} (see
#'    below), or between -2 and 2 in unconstrained space if \code{init_r = NULL}.
#'    \item \code{'0'}: All parameters are initialized to zero. 
#'    \item \code{'prior'}: Initializes parameters based on the specified prior.
#'    \item \code{NULL}: Initial values are provided by the corresponding init
#'    arguments defined below.
#'  }
#'  
#'  Note that \code{init = NULL} assigns initials for fixed effects, and
#'  variance co variance parameters (\code{vcov_init_0 = FALSE}) based on the
#'  individual setting for each parameter. If you want to initiate all
#'  parameters as \code{'random'}, then you must set \code{init = 'random'}
#'  which will be translated to \code{init = NULL} argument for \code{init =
#'  'rstan'} and \code{init = 'cmdstanr'}.
#'
#' @param init_r A positive real value specifying the range for random initial
#'   values (default \code{0.5}. This argument is used only when \code{init =
#'   'random'}. Note that the default setting for \code{Stan} is \code{2.0} to
#'   assign random initials between a range \code{-2.0, 2.0} on the
#'   unconstrained parameter space.
#'
#' @param a_init_beta Initial values for the fixed effect parameter, \code{a}
#'  (default \code{'random'}). Available options include:
#'  \itemize{
#'    \item \code{'0'}: Initializes the parameter to zero. \item
#'    \code{'random'}: Initializes with random values within a specified range.
#'    \item \code{'prior'}: Uses values drawn from the prior distribution. \item
#'    \code{'ymean'}: Initializes with the mean of the response variable. \item
#'    \code{'ymedian'}: Initializes with the median of the response variable.
#'    \item \code{'lm'}: Initializes with the coefficients from a simple linear
#'    model fitted to the data.
#'  }
#'  
#'  Note that options \code{'ymean'}, \code{'ymedian'}, and \code{'lm'} are only
#'  available for the fixed effect parameter \code{a}. For \code{univariate_by}
#'  and \code{multivariate} models, initial values can be the same across
#'  submodels (e.g., \code{a_init_beta = '0'}) or different for each submodel
#'  (e.g., \code{list(a_init_beta = '0', a_init_beta = 'lm')}).
#'
#' @param b_init_beta Initial values for the fixed effect parameter, \code{b}
#'   (default \code{'random'}). See \code{a_init_beta} for details on available
#'   options.
#' 
#' @param c_init_beta Initial values for the fixed effect parameter, \code{c}
#'   (default \code{'random'}). See \code{a_init_beta} for details on available
#'   options.
#'
#' @param d_init_beta Initial values for the fixed effect parameter, \code{d}
#'   (default \code{'random'}). See \code{a_init_beta} for details on available
#'   options.
#' 
#' @param s_init_beta Initial values for the fixed effect parameter, \code{s} 
#'  (default \code{'random'}). Available options include:
#'  \itemize{
#'    \item \code{'0'}: Initializes the parameter to zero. \item
#'    \code{'random'}: Initializes with random values within a specified range.
#'    \item \code{'prior'}: Uses values drawn from the prior distribution. \item
#'    \code{'lm'}: Initializes with the coefficients from a simple linear model
#'    fitted to the data.
#'  }
#'
#' @param a_cov_init_beta Initial values for the covariate(s) included in the
#'   fixed effect parameter, \code{a} (default \code{'random'}). Available
#'   options include:
#'  \itemize{
#'    \item \code{'0'}: Initializes the covariates to zero. \item
#'    \code{'random'}: Initializes with random values within a specified range.
#'    \item \code{'prior'}: Uses values drawn from the prior distribution. \item
#'    \code{'lm'}: Initializes with the coefficients from a simple linear model
#'    fitted to the data.
#'  }
#'  
#'  Note that the \code{'lm'} option is only available for
#'  \code{a_cov_init_beta} and not for covariates in other parameters such as
#'  \code{b}, \code{c}, or \code{d}.
#'
#' @param b_cov_init_beta Initial values for the covariate(s) included in the
#'   fixed effect parameter, \code{b} (default \code{'random'}). See
#'   \code{a_cov_init_beta} for details.
#' 
#' @param c_cov_init_beta Initial values for the covariate(s) included in the
#'   fixed effect parameter, \code{c} (default \code{'random'}). See
#'   \code{a_cov_init_beta} for details.
#'
#' @param d_cov_init_beta Initial values for the covariate(s) included in the
#'   fixed effect parameter, \code{d} (default \code{'random'}). See
#'   \code{a_cov_init_beta} for details.
#'
#' @param s_cov_init_beta Initial values for the covariate(s) included in the
#'   fixed effect parameter, \code{s} (default \code{'lm'}). See
#'   \code{a_cov_init_beta} for details. The option \code{'lm'} sets the spline
#'   coefficients obtained from a simple linear model fitted to the data.
#'   However, note that \code{s_cov_init_beta} serves as a placeholder and is
#'   not evaluated, as covariates are not allowed for the \code{s} parameter.
#'   For more details on covariates for \code{s}, refer to \code{s_formula}.
#' 
#' @param a_init_sd Initial value for the standard deviation of the group-level
#'   random effect parameter, \code{a} (default \code{'random'}). Available
#'   options are:
#' 
#'  \itemize{
#'    \item \code{'random'}: Initializes with random values within a specified
#'    range. 
#'    \item \code{'prior'}: Uses values drawn from the prior distribution. 
#'    \item \code{'ysd'}: Sets the standard deviation (\code{sd}) of the
#'    response variable as the initial value.
#'    \item \code{'ymad'}: Sets the median absolute deviation (\code{mad}) of
#'    the response variable as the initial value.
#'    \item \code{'lme_sd_a'}: Sets the initial value based on the standard
#'    deviation of the random intercept obtained from a linear mixed model
#'    (\code{nlme::lme()}) fitted to the data. If \code{nlme::lme()} fails to
#'    converge, the option \code{'lm_sd_a'} will be used automatically.
#'    \item \code{'lm_sd_a'}: Sets the square root of the residual variance
#'    obtained from a simple linear model applied to the data as the initial
#'    value.
#'  }
#'  
#'  Note that the options \code{'ysd'}, \code{'ymad'}, \code{'lme_sd_a'}, and
#'  \code{'lm_sd_a'} are available only for the random effect parameter \code{a}
#'  and not for other group-level random effects.
#'  
#'  Additionally, when fitting \code{univariate_by} and \code{multivariate}
#'  models, the user can set the same initial values for all sub-models, or
#'  different initial values for each sub-model.
#'  
#' @param b_init_sd Initial value for the standard deviation of the group-level
#'   random effect parameter, \code{b} (default \code{'random'}). Refer to
#'   \code{a_init_sd} for available options and details.
#'
#' @param c_init_sd Initial value for the standard deviation of the group-level
#'   random effect parameter, \code{c} (default \code{'random'}). Refer to
#'   \code{a_init_sd} for available options and details.
#'  
#' @param d_init_sd Initial value for the standard deviation of the group-level
#'   random effect parameter, \code{d} (default \code{'random'}). Refer to
#'   \code{a_init_sd} for available options and details.
#'  
#' @param a_cov_init_sd Initial values for the covariate(s) included in the
#'   random effect parameter \code{a} (default \code{'random'}). Available
#'   options include:
#'  \itemize{
#'    \item \code{'random'}: Random initialization.
#'    \item \code{'prior'}: Uses prior distribution values.
#'  }
#'  
#' @param b_cov_init_sd Initial values for the covariate(s) included in the
#'   random effect parameter \code{b} (default \code{'random'}). Refer to
#'   \code{a_cov_init_sd} for available options and details.
#'  
#' @param c_cov_init_sd Initial values for the covariate(s) included in the
#'   random effect parameter \code{c} (default \code{'random'}). Refer to
#'   \code{a_cov_init_sd} for available options and details.
#'  
#' @param d_cov_init_sd Initial values for the covariate(s) included in the
#'   random effect parameter \code{d} (default \code{'random'}). Refer to
#'   \code{a_cov_init_sd} for available options and details.
#'  
#' @param sigma_init_beta Initial values for the fixed effect distributional
#'   parameter \code{sigma} (default \code{'random'}). Available options
#'   include:
#'  \itemize{
#'    \item \code{'random'}: Random initialization.
#'    \item \code{'prior'}: Uses prior distribution values.
#'  }
#'  
#' @param sigma_cov_init_beta Initial values for the covariate(s) included in
#'   the fixed effect distributional parameter \code{sigma} (default
#'   \code{'random'}). Refer to \code{sigma_init_beta} for available options and
#'   details.
#'  
#' @param sigma_init_sd Initial value for the standard deviation of the
#'   distributional random effect parameter \code{sigma} (default
#'   \code{'random'}). The approach is the same as described earlier for the
#'   group-level random effect parameters such as \code{a} (See \code{a_init_sd}
#'   for details).
#'
#' @param sigma_cov_init_sd Initial values for the covariate(s) included in the
#'   distributional random effect parameter \code{sigma} (default
#'   \code{'random'}). The approach is the same as described for
#'   \code{a_cov_init_sd} (See \code{a_cov_init_sd} for details).
#' 
#' @param gr_init_cor Initial values for the correlation parameters of
#'   group-level random effects parameters (default \code{'random'}). Allowed
#'   options are:
#'  \itemize{
#'    \item \code{'random'}: Random initialization.
#'    \item \code{'prior'}: Uses prior distribution values.
#'    \item \code{'prior'}: A vector of length equal to the number of lower 
#'    triangle elements. For example, the initials for a model with three random 
#'    effects parameters can be specified as specified as
#'    \code{gr_init_cor = list(c(0.5. 0.5, 0.5))}
#'  }
#'  
#'  Note that when \code{vcov_init_0 = TRUE}, the \code{gr_init_cor} will be 
#'  set as '0'.
#'  
#' @param sigma_init_cor Initial values for the correlation parameters of
#'   distributional random effects parameter \code{sigma} (default
#'   \code{'random'}). Allowed options are:
#'  \itemize{
#'    \item \code{'random'}: Random initialization.
#'    \item \code{'prior'}: Uses prior distribution values.
#'  }
#'  
#' @param rsd_init_sigma Initial values for the residual standard deviation
#'   parameter, \code{sigma} (default \code{'random'}). Options available are:
#'  \itemize{
#'    \item \code{'0'}: Initializes the residual standard deviation to zero.
#'    \item \code{'random'}: Random initialization of the residual standard
#'    deviation.
#'    \item \code{'prior'}: Initializes the residual standard deviation based on
#'    prior distribution values.
#'    \item \code{'lme_rsd'}: Sets the initial value based on the standard
#'    deviation of residuals obtained from the linear mixed model
#'    (\code{nlme::lme()}) fitted to the data.
#'    \item \code{'lm_rsd'}: Sets the initial value as the square root of the
#'    residual variance from the simple linear model fitted to the data.
#'  }
#'  
#'  Note that if \code{nlme::lme()} fails to converge, the option
#'  \code{'lm_rsd'} is set automatically. The argument \code{rsd_init_sigma} is
#'  evaluated when both \code{dpar_formula} and \code{sigma_formula} are set to
#'  \code{NULL}.
#'  
#' @param dpar_init_sigma Initial values for the distributional parameter
#'   \code{sigma} (default \code{'random'}). The approach and available options
#'   are the same as described for \code{rsd_init_sigma}. This argument is
#'   evaluated only when \code{dpar_formula} is not \code{NULL}.
#'
#' @param dpar_cov_init_sigma Initial values for the covariate(s) included in
#'   the distributional parameter \code{sigma} (default \code{'random'}).
#'   Allowed options are \code{'0'}, \code{'random'}, and \code{'prior'}.
#'
#' @param autocor_init_acor Initial values for the autocorrelation parameter
#'   (see \code{autocor_formula} for details). Allowed options are \code{'0'},
#'   \code{'random'}, and \code{'prior'} (default \code{'random'}).
#'
#' @param autocor_init_unstr_acor Initial values for unstructured residual
#'   autocorrelation parameters (default \code{'random'}). Allowed options are
#'   \code{'0'}, \code{'random'}, and \code{'prior'}. The approach for setting
#'   initials for \code{autocor_init_unstr_acor} is the same as for
#'   \code{gr_init_cor}.
#'
#' @param mvr_init_rescor Initial values for the residual correlation parameter
#'   when fitting a \code{multivariate} model (default \code{'random'}). Allowed
#'   options are \code{'0'}, \code{'random'}, and \code{'prior'}.
#'
#' @param r_init_z Initial values for the standardized group-level random effect
#'   parameters (default \code{'random'}). These parameters are part of the
#'   Non-Centered Parameterization (NCP) approach used in the [brms::brm()].
#'
#' @param vcov_init_0 A logical to set initial values for variance (standard
#'   deviation) and covariance (correlation) parameters to zero (when
#'   \code{vcov_init_0 = TRUE}). This allows for setting custom initial values
#'   for the fixed effects parameters while keeping the variance-covariance
#'   parameters at zero. When \code{vcov_init_0 = FALSE} (default), then
#'   variance-covariance parameters are assigned random initial values unless
#'   each individual parameter has it own initial values setting (e.g.,
#'   \code{a_init_sd = 0}). Note that \code{vcov_init_0} is ignored when global
#'   initial values are assigned for all parameters via \code{init} argument.
#'   
#' @param jitter_init_beta A named list or numeric value to add a small amount
#'   of noise to an initial value or vector of initial values for the population
#'   level parameters.
#'   
#'   When \code{jitter_init_beta} is specified as a numeric value, it is treated
#'   as the percentage of perturbation applied to the initials. This value must
#'   be between \code{0} and \code{100}. Internally, the percentage is converted
#'   to a proportion (\code{percentage / 100}) and passed as the \code{amount}
#'   argument to the [base::jitter()] function. The \code{factor} argument is
#'   kept at its default value, \code{1}.
#'
#'   The default, \code{jitter_init_beta = NULL}, means no perturbation is
#'   applied, so the same initial values are used for all chains. For mild
#'   perturbation, you might use a value such as \code{jitter_init_beta = 10}.
#'
#'   Note that jitter is applied proportionally to the specified initial value,
#'   not as an absolute amount. For example, if the initial value is \code{100},
#'   setting \code{jitter_init_beta = 0.1} causes the perturbed value to fall
#'   within the range \code{90} to \code{110}. Conversely, if the initial value
#'   is \code{10}, the perturbed value will be within \code{9} to \code{11}.
#'
#'   If \code{jitter_init_beta} is provided as a named list, these elements are
#'   passed to the [base::jitter()] function. In addition to the \code{factor}
#'   and \code{amount} arguments, you may specify a \code{percent} argument,
#'   which is handled in the same way as a single numeric value. If both
#'   \code{percent} and \code{factor} are provided in the list, the effective
#'   perturbation is their product.
#'
#'   To use the default behavior of [base::jitter()], supply an empty
#'   \code{list()}, which will then be populated with the defaults:
#'   \code{list(..., factor = 1, amount = NULL)}. Please refer to the
#'   [base::jitter()] documentation for further details on how the \code{factor}
#'   and \code{amount} arguments affect the perturbation.
#'   
#'
#' @param jitter_init_sd A named list or numeric value to add a small amount of
#'   noise to an initial value or vector of initial values for the standard
#'   deviation of random effect parameters.For \code{jitter_init_sd} a
#'   reasonable option of setting one percent as the perturbed value
#'   \code{jitter_init_sd = 1} has been found to work well during early testing.
#'   See \code{jitter_init_beta} for details on various options available to
#'   perturb the initials.
#'
#' @param jitter_init_cor A named list or numeric value to add a small amount of
#'   noise to an initial value or vector of initial values for the correlations
#'   of random effect parameters.For \code{jitter_init_cor} a reasonable option
#'   of setting one percent as the perturbed value \code{jitter_init_sd = 0.1}
#'   has been found to work well during early testing. See
#'   \code{jitter_init_beta} for details on various options available to perturb
#'   the initials.
#' 
#' 
#' @param prior_data An optional argument (a named list, default \code{NULL})
#'   that can be used to pass information to the prior arguments for each
#'   parameter (e.g., \code{a_prior_beta}). The \code{prior_data} is
#'   particularly helpful when passing a long vector or matrix as priors. These
#'   vectors and matrices can be created in the R framework and then passed
#'   using the \code{prior_data}. For example, to pass a vector of location and
#'   scale parameters when setting priors for covariate coefficients (with 10
#'   dummy variables) included in the fixed effects parameter \code{a}, the
#'   following steps can be used:
#'  \itemize{
#'  \item Create the named objects \code{prior_a_cov_location} and 
#'  \code{prior_a_cov_scale} in the R environment:
#'  \code{prior_a_cov_location <- rnorm(n = 10, mean = 0, sd = 1)} 
#'  \code{prior_a_cov_scale <- rep(5, 10)}.
#'  \item Specify these objects in the \code{prior_data} list:
#'  \code{prior_data = list(prior_a_cov_location = prior_a_cov_location, 
#'  prior_a_cov_scale = prior_a_cov_scale)}.
#'  \item Use the \code{prior_data} objects to set up the priors:
#'  \code{a_cov_prior_beta = normal(prior_a_cov_location, prior_a_cov_scale)}.
#'  }
#'  
#' @param init_data An optional argument (a named list, default \code{NULL})
#'   that can be used to pass information to the initial arguments. The approach
#'   is identical to how \code{prior_data} is handled (as described above).
#'
#' @param init_custom Specify a custom initialization object (a named list). The
#'   named list is directly passed to the \code{init} argument without verifying
#'   the dimensions or name matching. If initial values are set for some
#'   parameters via parameter-specific arguments (e.g., \code{a_init_beta = 0}),
#'   \code{init_custom} will only be passed to those parameters that do not have
#'   initialized values. To override this behavior and use all of
#'   \code{init_custom} values regardless of parameter-specific initials, set
#'   \code{init = 'custom'}.
#'
#' @param expose_function An optional argument (logical, default \code{FALSE})
#'   to indicate whether to expose the Stan function used in model fitting.
#'
#' @param get_stancode An optional argument (logical, default \code{FALSE}) to
#'   retrieve the Stan code (see \code{[brms::stancode()]} for details).
#'
#' @param get_standata An optional argument (logical, default \code{FALSE}) to
#'   retrieve the Stan data (see \code{[brms::standata()]} for details).
#'
#' @param get_formula An optional argument (logical, default \code{FALSE}) to
#'   retrieve the model formula (see \code{[brms::brmsformula()]} for details).
#'
#' @param get_stanvars An optional argument (logical, default \code{FALSE}) to
#'   retrieve the Stan variables (see \code{[brms::stanvar()]} for details).
#'
#' @param get_priors An optional argument (logical, default \code{FALSE}) to
#'   retrieve the priors (see \code{[brms::get_prior()]} for details). Note 
#'   that \code{get_priors = TRUE} will return priors based on the final code.
#'   In case user want to return the basic default priors, then it can be 
#'   achieved by setting \code{get_priors = "default"}. 
#'
#' @param get_priors_eval An optional argument (logical, default \code{FALSE})
#'   to retrieve the priors specified by the user.
#'
#' @param get_init_eval An optional argument (logical, default \code{FALSE}) to
#'   retrieve the initial values specified by the user.
#'
#' @param validate_priors An optional argument (logical, default \code{FALSE})
#'   to validate the specified priors (see \code{[brms::validate_prior()]} for
#'   details).
#'
#' @param set_self_priors An optional argument (default \code{NULL}) to manually
#'   specify the priors. \code{set_self_priors} is passed directly to
#'   \code{[brms::brm()]} without performing any checks.
#'
#' @param add_self_priors An optional argument (default \code{NULL}) to append
#'   part of the prior object. This is for internal use only.
#'
#' @param set_replace_priors An optional argument (default \code{NULL}) to
#'   replace part of the prior object. This is for internal use only.
#'
#' @param set_same_priors_hierarchy An optional argument (default \code{NULL})
#'   to replace part of the prior object. This is for internal use only.
#'
#' @param outliers An optional argument (default \code{NULL}) to remove
#'   outliers. This should be a named list passed directly to
#'   \code{[sitar::velout()]} and \code{[sitar::zapvelout()]} functions. This is
#'   for internal use only.
#'
#' @param unused An optional formula defining variables that are unused in the
#'   model but should still be stored in the model's data frame. Useful when
#'   variables are needed during post-processing.
#'
#' @param chains The number of Markov chains (default 4).
#'
#' @param iter The total number of iterations per chain, including warmup
#'   (default 2000).
#' 
#' @param warmup A positive integer specifying the number of warmup (aka
#'   burn-in) iterations. This also specifies the number of iterations used for
#'   stepsize adaptation, so warmup draws should not be used for inference. The
#'   number of warmup iterations should not exceed \code{iter}, and the default
#'   is \code{iter/2}.
#'
#' @param thin A positive integer specifying the thinning interval. Set
#'   \code{thin > 1} to save memory and computation time if \code{iter} is
#'   large. Thinning is often used in cases with high autocorrelation of MCMC
#'   draws. An indication of high autocorrelation is poor mixing of chains
#'   (i.e., high \code{rhat} values) despite the model recovering parameters
#'   well. A useful diagnostic to check for autocorrelation of MCMC draws is the
#'   \code{mcmc_acf} function from the \pkg{bayesplot} package.
#'
#' @param cores Number of cores to be used when executing the chains in
#'   parallel. See [brms::brm()] for details. Unlike [brms::brm()], which
#'   defaults the \code{cores} argument to \code{cores=getOption("mc.cores",
#'   1)}, the default \code{cores} in the \pkg{bsitar} package is
#'   \code{cores=getOption("mc.cores", 'optimize')}, which optimizes the
#'   utilization of system resources. The maximum number of cores that can be
#'   deployed is calculated as the maximum number of available cores minus 1.
#'   When the number of available cores exceeds the number of chains (see
#'   \code{chains}), then the number of cores is set equal to the number of
#'   chains.
#'  
#'  Another option is to set \code{cores} as \code{getOption("mc.cores",
#'  'maximise')}, which sets the number of cores to the maximum number of cores
#'  available on the system regardless of the number of chains specified.
#'  Alternatively, the user can specify \code{cores} in the same way as
#'  [brms::brm()] with \code{getOption("mc.cores", 1)}.
#'  
#'  These options can be set globally using \code{options(mc.cores = x)}, where
#'  \code{x} can be \code{'optimize'}, \code{'maximise'}, or \code{1}. The
#'  \code{cores} argument can also be directly specified as an integer (e.g.,
#'  \code{cores = 4}).
#'  
#' @param backend A character string specifying the package to be used when
#'   executing the Stan model. The available options are \code{"rstan"} (the
#'   default) or \code{"cmdstanr"}. The backend can also be set globally for the
#'   current \R session using the \code{"brms.backend"} option. See
#'   [brms::brm()] for more details.
#'
#'@param threads Number of threads to be used in within-chain parallelization.
#'  Note that unlike the [brms::brm()] which sets the \code{threads} argument as
#'  \code{getOption("brms.threads", NULL)} implying that no within-chain
#'  parallelization is used by default, the \pkg{bsitar} package, by default,
#'  sets \code{threads} as \code{getOption("brms.threads", 'optimize')} to
#'  utilize the available resources from the modern computing systems. The
#'  number of threads per chain is set as the maximum number of cores available
#'  minus 1. Another option is to set \code{threads} as
#'  \code{getOption("brms.threads", 'maximise')} which set the number threads
#'  per chains same as the  maximum number of cores available. User can also set
#'  the \code{threads} similar to the \code{brms} i.e.,
#'  \code{getOption("brms.threads", NULL)}. All these three options can be set
#'  globally as \code{options(brms.threads = x}) where x can be
#'  \code{'optimize'}, \code{'maximise'} or \code{NULL}.
#'  Alternatively, the number of threads can be set directly as \code{threads
#'  = threading(x)} where \code{X} is an integer. Other arguments that can be
#'  passed to the \code{threads} are \code{grainsize} and the \code{static}. See
#'  [brms::brm()] for further details on within-chain parallelization.
#'  
#' @param opencl The platform and device IDs of the OpenCL device to use for GPU
#'   support during model fitting. If you are unsure about the IDs of your
#'   OpenCL device, \code{c(0,0)} is typically the default that should work. For
#'   more details on how to find the correct platform and device IDs, refer to
#'   [brms::opencl()]. This parameter can also be set globally for the current
#'   \R session using the \code{"brms.opencl"} option.
#' 
#' @param normalize Logical flag indicating whether normalization constants
#'   should be included in the Stan code (default is \code{TRUE}). If set to
#'   \code{FALSE}, normalization constants are omitted, which may increase
#'   sampling efficiency. However, this requires Stan version >= 2.25. Note that
#'   setting \code{normalize = FALSE} will disable some post-processing
#'   functions, such as [brms::bridge_sampler()]. This option can be controlled
#'   globally via the \code{brms.normalize} option.
#'
#' @param algorithm A character string specifying the estimation method to use. 
#'  Available options are:
#'  \itemize{
#'    \item \code{"sampling"} (default): Markov Chain Monte Carlo (MCMC) method.
#'    \item \code{"meanfield"}: Variational inference with independent normal
#'    distributions.
#'    \item \code{"fullrank"}: Variational inference with a multivariate normal
#'    distribution.
#'    \item \code{"fixed_param"}: Sampling from fixed parameter values.
#'  }
#'  This parameter can be set globally via the \code{"brms.algorithm"} option
#'  (see \code{\link{options}} for more details).
#'
#' @param control A named \code{list} to control the sampler's behavior. The
#'   default settings are the same as those in [brms::brm()], with one
#'   exception: the \code{max_treedepth} has been increased from 10 to 12 to
#'   better explore the typically challenging posterior geometry in nonlinear
#'   models. However, the \code{adapt_delta}, which is often increased for
#'   nonlinear models, retains its default value of 0.8 to avoid unnecessarily
#'   increasing sampling time. For full details on control parameters and their
#'   default values, refer to [brms::brm()].
#' 
#' @param pathfinder_args A named \code{list} of arguments passed to the
#'   \code{'pathfinder'} algorithm. This is used to set
#'   \code{'pathfinder'}-based initial values for the \code{'MCMC'} sampling.
#'   Note that \code{'pathfinder_args'} currently only works when \code{backend
#'   = "cmdstanr"}. If \code{pathfinder_args} is not \code{NULL} and the user
#'   specifies \code{backend = "rstan"}, the backend will automatically be
#'   changed to \code{cmdstanr}.
#'
#' @param pathfinder_init A logical value (default \code{FALSE}) indicating
#'   whether to use initial values from the \code{'pathfinder'} algorithm when
#'   fitting the final model (i.e., \code{'MCMC'} sampling). Note that
#'   \code{'pathfinder_args'} currently works only when \code{backend =
#'   "cmdstanr"}. If \code{pathfinder_args} is not \code{NULL} and the user
#'   specifies \code{backend = "rstan"}, the backend will automatically switch
#'   to \code{cmdstanr}. The arguments passed to the \code{'pathfinder'}
#'   algorithm are specified via \code{'pathfinder_args'}; if
#'   \code{'pathfinder_args'} is \code{NULL}, the default arguments from
#'   \code{'cmdstanr'} will be used.
#'   
#' @param data_custom A \code{data.frame} object (default \code{NULL}). This is
#'   mainly for internal testing and not be used for routine model fitting.
#'   
#' @param genquant_xyadj A logical (default \code{NULL}) indicating whether to
#'   generate \code{xyadj} in the \code{generated quantities} block of the
#'   \code{Stan}.
#' 
#' @param sample_prior A character string indicating whether to draw samples
#'   from the priors in addition to the posterior draws. Options are \code{"no"}
#'   (the default), \code{"yes"}, and \code{"only"}. These prior draws can be
#'   used for various purposes, such as calculating Bayes factors for point
#'   hypotheses via [brms::hypothesis()]. Note that improper priors (including
#'   the default improper priors used by \code{brm}) are not sampled. For proper
#'   priors, see [brms::set_prior()]. Also, prior draws for the overall
#'   intercept are not obtained by default for technical reasons. See
#'   [brms::brmsformula()] for instructions on obtaining prior draws for the
#'   intercept. If \code{sample_prior} is set to \code{"only"}, draws will be
#'   taken solely from the priors, ignoring the likelihood, which allows you to
#'   generate draws from the prior predictive distribution. In this case, all
#'   parameters must have proper priors.
#'
#' @param save_pars An object generated by [brms::save_pars()] that controls
#'   which parameters should be saved in the model. This argument does not
#'   affect the model fitting process itself but provides control over which
#'   parameters are retained in the final output.
#'
#' @param drop_unused_levels A logical value indicating whether unused factor
#'   levels in the data should be dropped. The default is \code{TRUE}.
#'
#' @param stan_model_args A \code{list} of additional arguments passed to
#'   \code{\link[rstan:stan_model]{rstan::stan_model}} when using the
#'   \code{backend = "rstan"} or \code{backend = "cmdstanr"}. This allows
#'   customization of how models are compiled.
#' 
#' @param refresh An integer specifying the frequency of printing every nth
#'   iteration. By default, \code{NULL} indicates that the refresh rate will be
#'   automatically set by [brms::brm()]. Setting \code{refresh} is especially
#'   useful when \code{thin} is greater than \code{1}, in which case the refresh
#'   rate is recalculated as (\code{refresh} * \code{thin}) / \code{thin}.
#'
#' @param silent A verbosity level between \code{0} and \code{2}. When set to
#'   \code{1} (the default), most informational messages from the compiler and
#'   sampler are suppressed. Setting it to \code{2} suppresses even more
#'   messages. The sampling progress is still printed. To turn off all printing,
#'   set \code{refresh = 0}. Additionally, when using \code{backend = "rstan"},
#'   you can prevent the opening of additional progress bars by setting
#'   \code{open_progress = FALSE}.
#' 
#' @param seed An integer or \code{NA} (default) specifying the seed for random
#'   number generation, ensuring reproducibility of results. If set to
#'   \code{NA}, \pkg{Stan} will randomly select the seed.
#'
#' @param save_model A character string or \code{NULL} (default). If provided,
#'   the Stan code for the model will be saved in a text file with the name
#'   corresponding to the string specified in \code{save_model}.
#'
#' @param fit An instance of class \code{brmsfit} from a previous fit (default
#'   is \code{NA}). If a \code{brmsfit} object is provided, the compiled model
#'   associated with the fitted result is reused, and any arguments that modify
#'   the model code or data are ignored. It is generally recommended to use the
#'   \code{\link[brms:update.brmsfit]{update}} method for this purpose, rather
#'   than directly passing the \code{fit} argument.
#'
#' @param file Either \code{NULL} or a character string. If a character string
#'   is provided, the fitted model object is saved using \code{\link{saveRDS}}
#'   in a file named after the string supplied in \code{file}. The \code{.rds}
#'   extension is automatically added. If the specified file already exists, the
#'   existing model object is loaded and returned instead of refitting the
#'   model. To overwrite an existing file, you must manually remove the file or
#'   specify the \code{file_refit} argument. The file name is stored within the
#'   \code{brmsfit} object for later use.
#'
#' @param file_refit Modifies when the fit stored via the \code{file} argument
#'   is re-used. This can be set globally for the current \R session via the
#'   \code{"brms.file_refit"} option (see \code{\link{options}}). The possible
#'   options are:
#'  \itemize{
#'    \item \code{"never"} (default): The fit is always loaded if it exists, and
#'    fitting is skipped.
#'    \item \code{"always"}: The model is always refitted, regardless of
#'    existing fits.
#'    \item \code{"on_change"}: The model is refitted only if the model, data,
#'    algorithm, priors, \code{sample_prior}, \code{stanvars}, covariance
#'    structure, or similar parameters have changed.
#'  }
#'  
#'  If you believe a false positive occurred, you can use
#'  \code{[brms::brmsfit_needs_refit()]} to investigate why a refit is deemed
#'  necessary. A refit will not be triggered for changes in additional
#'  parameters of the fit (e.g., initial values, number of iterations, control
#'  arguments). A known limitation is that a refit will be triggered if
#'  within-chain parallelization is switched on/off.
#'  
#'
#' @param future Logical; If \code{TRUE}, the \pkg{\link[future:future]{future}}
#'   package is used for parallel execution of the chains. In this case, the
#'   \code{cores} argument will be ignored. The execution type is controlled via
#'   \code{\link[future:plan]{plan}} (see the examples section below). This
#'   argument can be set globally for the current \R session via the
#'   \code{"future"} option.
#'   
#' @param sum_zero Currently ignored. Placeholder for future development.
#' 
#' @param global_args Currently ignored. Placeholder for future development.
#'
#' @param parameterization A character string specifying the type of
#'   parameterization to use for drawing group-level random effects. Options
#'   are: \code{'ncp'} for Non-Centered Parameterization (NCP), and \code{'cp'}
#'   for Centered Parameterization (CP).
#'  
#'  The NCP is generally recommended when the likelihood is weak (e.g., few
#'  observations per individual) and is the default approach (and only option)
#'  in [brms::brm()].
#'  
#'  The CP parameterization is typically more efficient when a relatively large
#'  number of observations are available across individuals. We consider a
#'  'relatively large number' as at least 10 repeated measurements per
#'  individual. If there are fewer than 10, NCP is used automatically. This
#'  behavior applies only when \code{parameterization = NULL}. To explicitly set
#'  CP parameterization, use \code{parameterization = 'cp'}.
#'  
#'  Note that since [brms::brm()] does not support CP, the \code{stancode}
#'  generated by [brms::brm()] is edited internally before fitting the model
#'  using [rstan::rstan()] or \code{"cmdstanr"}, depending on the chosen
#'  \code{backend}. Therefore, CP parameterization is considered experimental
#'  and may fail if the structure of the generated \code{stancode} changes in
#'  future versions of [brms::brm()].
#'  
#' @param verbose An optional argument (logical, default \code{FALSE}) to
#'   indicate whether to print information collected during setting up the model
#'   formula priors, and initials. As an example, the user might be interested
#'   in knowing the response variables created for the sub model when fitting a
#'   univariate-by-subgroup model. This information can then be used in setting
#'   the desired order of options passed to each such model such as \code{df},
#'   \code{prior}, \code{initials} etc.
#'   
#' @param ... Further arguments passed to [brms::brm()]. This may include
#'   additional arguments that are either passed directly to the underlying
#'   model fitting function, or used for internal purposes. Specifically, the
#'   \code{...} can also be used to pass arguments used for testing and
#'   debugging, such as: \code{match_sitar_a_form},
#'   \code{sigmamatch_sitar_a_form}, \code{displayit}, \code{setcolh},
#'   \code{setcolb}, \code{decomp} and \code{smat}. Note that for all arguments
#'   (such as intercept) to correctly pass to the respective functions, use
#'   \code{smat} and not \code{stype} when setting up the spline arguments.
#'  
#'  These internal arguments are typically not used in regular model fitting but
#'  can be relevant for certain testing scenarios or advanced customization.
#'  Users are generally not expected to interact with these unless working on
#'  debugging or testing specific features of the model fitting process.
#'
#' @return An object of class \code{brmsfit, bsitar}, which contains the
#'   posterior draws, model coefficients, and other useful information related
#'   to the model fitting. This object includes details such as the fitted
#'   model, the data used, prior distributions, and any other relevant outputs
#'   from the Stan model fitting process. The resulting object can be used for
#'   further analysis, diagnostics, and post-processing, including model summary
#'   statistics, predictions, and visualizations.
#'
#'@export
#'
#'@inheritParams brms::brm
#'
#'@importFrom stats as.formula coef df dist filter fitted gaussian lm mad median
#'  model.matrix predict quantile rbeta sd setNames smooth.spline rnorm runif
#'  rcauchy rexp rlnorm rgamma rlnorm loess na.omit residuals complete.cases
#'  deriv formula update
#' 
#'@importFrom rlang .data
#'
#'@importFrom utils combn head installed.packages packageVersion tail data
#'
#'@importFrom Rdpack reprompt
#'
#'@import brms
#'
#' @note The package is under continuous development, and new models,
#'   post-processing features, and improvements are being actively worked on.
#'   Keep an eye on future releases for additional functionality and updates to
#'   enhance model fitting, diagnostics, and analysis capabilities.
#'
#' @seealso [brms::brm()] [brms::brmsformula()] [brms::prior()]
#'
#' @references
#'  \insertAllCited{}
#'  
#' @inherit berkeley author
#' 
#' @examples
#' \donttest{
#' # Below, we fit a SITAR model to a subset of the Berkley height data, 
#' # specifically the data for 70 girls between the ages of 8 and 18.  
#' # This subset is used as an example in the vignette for the 'sitar' package.
#' #
#' # The original Berkley height data contains repeated growth measurements for
#' # 66 boys and 70 girls (ages 0-21). For this example, we use a subset of the 
#' # data for 70 girls aged 8 to 18 years.
#' #
#' # For details on the full Berkley height dataset, refer to 'sitar' package
#' # documentation (help file: ?sitar::berkeley). Further details on the subset
#' # of the data used here can be found in the vignette ('Fitting_models_with_SITAR', 
#' # package = 'sitar').
#' 
#' # Load the 'berkeley_exdata' that has been pre-saved
#' berkeley_exdata <- getNsObject(berkeley_exdata)
#' 
#' # Fit frequentist SITAR model with df = 3 using the sitar package 
#' 
#' model_ml <- sitar::sitar(x = age, y = height, id = id, 
#'                           df = 3, 
#'                           data = berkeley_exdata, 
#'                           xoffset = 'mean',
#'                           fixed = 'a+b+c', 
#'                           random = 'a+b+c',
#'                           a.formula = ~1, 
#'                           b.formula = ~1, 
#'                           c.formula = ~1
#'                           )
#' 
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid time-consuming model estimation, the Bayesian SITAR model fit has 
#' # been saved as an example fit ('berkeley_exfit'). This model was fit using 
#' # 2 chains (2000 iterations per chain) with thinning set to 6 for memory  
#' # efficiency. Users are encouraged to refit the model using default settings 
#' # (4 chains, 2000 iterations per chain, thin = 1) as suggested by the Stan
#' # team. Note that with thinning set to 6 (thin = 6), only one sixth of total  
#' # draws will be saved and hence the effective sample size is expected to be 
#' # small.
#' 
#' # Check if the pre-saved model 'berkeley_exfit' exists
#' # berkeley_exfit <- bsitar:::berkeley_exfit
#' 
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#'  
#' if(exists('berkeley_exfit')) {
#'   model <- berkeley_exfit
#' } else {
#'   # Fit model with default priors
#'   # Refer to the documentation for prior on each parameter
#'   model <- bsitar(x = age, y = height, id = id, 
#'                   df = 3, 
#'                   data = berkeley_exdata,
#'                   xoffset = 'mean', 
#'                   fixed = 'a+b+c', 
#'                   random = 'a+b+c',
#'                   a_formula = ~1, 
#'                   b_formula = ~1, 
#'                   c_formula = ~1, 
#'                   threads = brms::threading(NULL),
#'                   chains = 2, cores = 2, iter = 1000, thin = 6)
#'                   
#' }
#' 
#' # Generate model summary
#' summary(model)
#' 
#' # Compare model summary with the frequentist SITAR model
#' print(model_ml)
#' 
#' # Check model fit via posterior predictive checks using plot_ppc.
#' # This function is based on pp_check from the 'brms' package.
#' plot_ppc(model, ndraws = NULL)
#' 
#' # Plot distance and velocity curves using plot_conditional_effects.
#' # This function works like conditional_effects from the 'brms' package,
#' # with the added option to plot velocity curves.
#' 
#' # Distance curve
#' plot_conditional_effects(model, deriv = 0)
#' 
#' # Velocity curve
#' plot_conditional_effects(model, deriv = 1)
#' 
#' # Plot distance and velocity curves along with parameter estimates using 
#' # plot_curves (similar to plot.sitar from the sitar package).
#' plot_curves(model, apv = TRUE)
#' 
#' # Compare plots with the frequentist SITAR model
#' plot(model_ml)
#' }
#'
bsitar <- function(x,
                   y,
                   id,
                   data,
                   df = 4,
                   knots = NA,
                   knots_selection = NULL,
                   fixed = 'a + b + c',
                   random = 'a + b + c',
                   xoffset = 'mean',
                   bstart = xoffset,
                   cstart = 0,
                   xfun = NULL,
                   yfun = NULL,
                   xfunxoffset = TRUE, 
                   bound = 0.04,
                   stype = 'nsk',
                   terms_rhs = NULL,
                   a_formula = ~ 1,
                   b_formula = ~ 1,
                   c_formula = ~ 1,
                   d_formula = ~ 1,
                   s_formula = ~ 1,
                   a_formula_gr = ~ 1,
                   b_formula_gr = ~ 1,
                   c_formula_gr = ~ 1,
                   d_formula_gr = ~ 1,
                   a_formula_gr_str = NULL,
                   b_formula_gr_str = NULL,
                   c_formula_gr_str = NULL,
                   d_formula_gr_str = NULL,
                   d_adjusted = FALSE,
                   sigma_formula = NULL,
                   sigma_formula_gr = NULL,
                   sigma_formula_gr_str = NULL,
                   sigma_formula_manual = NULL,
                   sigmax = NULL,
                   sigmaid  = NULL,
                   sigmadf = 4,
                   sigmaknots = NA,
                   sigmafixed = NULL,
                   sigmarandom = "",
                   sigmaxoffset = 'mean',
                   sigmabstart = 'sigmaxoffset',
                   sigmacstart = 0,
                   sigmaxfun = NULL,
                   sigmabound = 0.04,
                   sigmaxfunxoffset = TRUE,
                   dpar_formula = NULL,
                   autocor_formula = NULL,
                   family = gaussian(),
                   custom_family = NULL,
                   custom_formula = NULL,
                   custom_prior = NULL,
                   custom_stanvars  = NULL,
                   group_arg = list(
                     groupvar = NULL,
                     by = NULL,
                     cor = 'un',
                     cov = NULL,
                     dist = 'gaussian'
                   ),
                   sigma_group_arg = list(
                     groupvar = NULL,
                     by = NULL,
                     cor = un,
                     cov = NULL,
                     dist = 'gaussian'
                   ),
                   univariate_by = list(by = NA, 
                                        cor = 'un', 
                                        terms = 'subset'),
                   multivariate = list(mvar = FALSE,
                                       cor = 'un',
                                       rescor = TRUE,
                                       rcorr_by = NULL,
                                       rcorr_gr = NULL,
                                       rcorr_method = NULL,
                                       rcorr_prior = NULL),
                   a_prior_beta = normal(ymean, ysd, autoscale = FALSE),
                   b_prior_beta = normal(0, 2, autoscale = FALSE),
                   c_prior_beta = normal(0, 1, autoscale = FALSE),
                   d_prior_beta = normal(0, 1.0, autoscale = FALSE),
                   s_prior_beta = normal(lm, lm, autoscale = FALSE),
                   a_cov_prior_beta = normal(0, 5.0, autoscale = FALSE),
                   b_cov_prior_beta = normal(0, 1.0, autoscale = FALSE),
                   c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                   d_cov_prior_beta = normal(0, 1.0, autoscale = FALSE),
                   s_cov_prior_beta = normal(lm, lm, autoscale = FALSE),
                   a_prior_sd = normal(0, ysd, autoscale = FALSE),
                   b_prior_sd = normal(0, 2, autoscale = FALSE),
                   c_prior_sd = normal(0, 1, autoscale = FALSE),
                   d_prior_sd = normal(0, 1.0, autoscale = FALSE),
                   a_cov_prior_sd = normal(0, 5.0, autoscale = FALSE),
                   b_cov_prior_sd = normal(0, 1.0, autoscale = FALSE),
                   c_cov_prior_sd = normal(0, 0.1, autoscale = FALSE),
                   d_cov_prior_sd = normal(0, 1.0, autoscale = FALSE),
                   a_prior_sd_str = NULL,
                   b_prior_sd_str = NULL,
                   c_prior_sd_str = NULL,
                   d_prior_sd_str = NULL,
                   a_cov_prior_sd_str = NULL,
                   b_cov_prior_sd_str = NULL,
                   c_cov_prior_sd_str = NULL,
                   d_cov_prior_sd_str = NULL,
                   sigma_prior_beta = normal(0, 1, autoscale = FALSE),
                   sigma_cov_prior_beta = normal(0, 0.5, autoscale = FALSE),
                   sigma_prior_sd = normal(0, 0.25, autoscale = FALSE),
                   sigma_cov_prior_sd = normal(0, 0.15, autoscale = FALSE),
                   sigma_prior_sd_str = NULL,
                   sigma_cov_prior_sd_str = NULL,
                   rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                   dpar_prior_sigma = normal(0, ysd, autoscale = FALSE),
                   dpar_cov_prior_sigma = normal(0, 1, autoscale = FALSE),
                   autocor_prior_acor = uniform(-1, 1, autoscale = FALSE),
                   autocor_prior_unstr_acor = lkj(1),
                   gr_prior_cor = lkj(1),
                   gr_prior_cor_str = lkj(1),
                   sigma_prior_cor = lkj(1),
                   sigma_prior_cor_str = lkj(1),
                   mvr_prior_rescor = lkj(1),
                   init = NULL,
                   init_r = 0.5,
                   a_init_beta = 'lm',
                   b_init_beta = 0,
                   c_init_beta = 0,
                   d_init_beta = 0,
                   s_init_beta = 'lm',
                   a_cov_init_beta = 0,
                   b_cov_init_beta = 0,
                   c_cov_init_beta = 0,
                   d_cov_init_beta = 0,
                   s_cov_init_beta = 0,
                   a_init_sd = 'random',
                   b_init_sd = 'random',
                   c_init_sd = 'random',
                   d_init_sd = 'random',
                   a_cov_init_sd = 'random',
                   b_cov_init_sd = 'random',
                   c_cov_init_sd = 'random',
                   d_cov_init_sd = 'random',
                   sigma_init_beta = 'random',
                   sigma_cov_init_beta = 'random',
                   sigma_init_sd = 'random',
                   sigma_cov_init_sd = 'random',
                   gr_init_cor = 'random',
                   sigma_init_cor = 'random',
                   rsd_init_sigma = 'random',
                   dpar_init_sigma = 'random',
                   dpar_cov_init_sigma = 'random',
                   autocor_init_acor = 'random',
                   autocor_init_unstr_acor = 'random',
                   mvr_init_rescor = 'random',
                   r_init_z = 'random',
                   vcov_init_0 = FALSE,
                   jitter_init_beta = NULL,
                   jitter_init_sd = NULL,
                   jitter_init_cor = NULL,
                   prior_data = NULL,
                   init_data = NULL,
                   init_custom = NULL,
                   expose_function = FALSE,
                   get_stancode = FALSE,
                   get_standata = FALSE,
                   get_formula = FALSE,
                   get_stanvars = FALSE,
                   get_priors = FALSE,
                   get_priors_eval = FALSE,
                   get_init_eval = FALSE,
                   validate_priors = FALSE,
                   set_self_priors = NULL,
                   add_self_priors = NULL,
                   set_replace_priors = NULL,
                   set_same_priors_hierarchy = FALSE,
                   outliers = NULL, 
                   unused = NULL,
                   chains = 4,
                   iter = 2000,
                   warmup = floor(iter / 2),
                   thin = 1,
                   cores = getOption("mc.cores", "optimize"),
                   backend = getOption("brms.backend", "rstan"),
                   threads = getOption("brms.threads", "optimize"),
                   opencl = getOption("brms.opencl", NULL),
                   normalize = getOption("brms.normalize", TRUE),
                   algorithm = getOption("brms.algorithm", "sampling"),
                   control = list(adapt_delta = 0.9, max_treedepth = 15),
                   empty = FALSE,
                   rename = TRUE,
                   pathfinder_args = NULL,
                   pathfinder_init = FALSE,
                   data2 = NULL,
                   data_custom = NULL,
                   genquant_xyadj  = FALSE,
                   sample_prior = 'no',
                   save_pars = NULL,
                   drop_unused_levels = TRUE,
                   stan_model_args = list(),
                   refresh = NULL,
                   silent = 1,
                   seed = 123,
                   save_model = NULL,
                   fit = NA,
                   file = NULL,
                   file_compress = TRUE,
                   file_refit = getOption("brms.file_refit", "never"),
                   future = getOption("future", FALSE),
                   sum_zero = FALSE,
                   global_args = NULL,
                   parameterization = 'ncp',
                   verbose = FALSE,
                   ...) {
  
  
  mcall <- match.call()
  no_default_args <- c("x", "y", "id", "data", "...")
  
  if(is.null(global_args)) {
    global_args <- FALSE
  } else if(!is.null(global_args)) {
    if(is.logical(global_args)) {
      global_args <- global_args
    } else {
      stop2c("'global_args' must be either NULL or a logical")
    }
  }
  
  if(global_args) {
    mcall <- mcall_dictionary(mcall, envir = NULL, xenvir = NULL, 
                                        exceptions = no_default_args)
  }
  
  if(!global_args) {
    call_eval_globals_in_mcall <- TRUE
  } else {
    call_eval_globals_in_mcall <- FALSE
  }

  # 25.01.2026
  set_eval_globals_in_mcall <- FALSE
  if(call_eval_globals_in_mcall) {
    if(called_via_do_call()) {
      # nothing
    } else if(called_via_CustomDoCall()) {
      stop("use 'do.call', not the 'CustomDoCall' for bsitar()")
    } else {
      set_eval_globals_in_mcall <- TRUE
      set_eval_globals_in_mcall_names <- names(mcall)
      set_exceptions <- c("data", "...")
      set_exceptions <- c(set_exceptions, 'family')
      mcall <- eval_globals_in_mcall(mcall, exceptions = set_exceptions) 
    }
  } # if(call_eval_globals_in_mcall) {
  
  
  if(!'init' %in% names(mcall)) {
    mcall <-rlang::call_modify(mcall, init = init)
  }
  if(!'init_r' %in% names(mcall)) {
    mcall <-rlang::call_modify(mcall, init_r = init_r)
  }
  
  
  mcall_ <- mcall
  
   data_check_for_modifications <- FALSE
   if(is.null(mcall_$data)) {
     stop2c("Data argument must be specified")
   } else if(is.language(mcall_$data)) {
     data_check_for_modifications <- TRUE
     data_name_str_check <- deparse(mcall_$data)
   } else if(is.symbol(mcall_$data)) {
     data_check_for_modifications <- TRUE
     data_name_str_check <- deparse(mcall_$data)
   } else if(is.data.frame(mcall_$data)) {
     data_check_for_modifications <- FALSE
   } else if(tibble::is.tibble(mcall_$data)) {
     data_check_for_modifications <- FALSE
   }

   if(data_check_for_modifications) {
     data_name_str_check <- deparse(mcall_$data)
     data_name_str_check <- gsub("\"", "", data_name_str_check, fixed = T)
     data_name_str       <- check_forpipe(data_name_str_check, return = 'name')
     data_name_pipe      <- check_forpipe(data_name_str_check, return = 'logical')
     data_name_str_attr  <- check_forpipe(data_name_str_check, return = 'attr')
   } else {
     data_name_pipe <- FALSE
   }

   if(data_name_pipe) {
     assign(data_name_str_attr, eval(mcall_$data) )
     mcall_$data <- as.symbol(data_name_str_attr)
   }
   
   # This when data argument is a data frame. e.g., during update_model()
   if(is.data.frame(mcall_$data) | tibble::is_tibble(mcall_$data)) {
     data_name_str_attr <- 'data'
     assign(data_name_str_attr, mcall_$data)
     mcall_$data <- as.symbol(data_name_str_attr)
   }
 
   threads_char <- function(threads_, chains = 1) {
     if(is_emptyx(chains)) chains <- 1
     if(!is.character(threads_)) {
       stop2c("threads_ must be a string of length 1")
     }
     if(length(threads_) > 1) {
       stop2c("threads_ must be a string of length 1")
     }
     allowed_temp_threads_char <- c('optimize', 'maximise')
     if(!threads_ %in% allowed_temp_threads_char) {
       stop2c("Allowed string options for 'threads' are: ",
              collapse_comma(allowed_temp_threads_char))
     }
     if( is.character(threads_) & threads_ == "maximise") {
       max.threads <- 
         as.numeric(future::availableCores(methods = "system", omit = 0))
       if(max.threads < 1) max.threads <- 1
     } else if( is.character(threads_) & threads_ == "optimize") {
       max.threads <- 
         as.numeric(future::availableCores(methods = "system", omit = 1))
       if(max.threads < 1) max.threads <- 1
       max.threads <- floor(max.threads /  chains)
     }
     return(max.threads)
   }
   
   
   
  # Check and allow setting threads as NULL or integer
   mcall_threads_ <- eval(mcall$threads) # For CustomDoCall
  deparse_mcall_threads_check <- paste(deparse(mcall_threads_), collapse = "")
  deparse_sub_mcall_threads_check <- deparse(substitute(mcall_threads_))
  deparse_sub_mcall_threads_check <- paste(deparse_sub_mcall_threads_check, 
                                           collapse = "")
  deparse_sub_mcall_threads_check <- gsub_space(deparse_sub_mcall_threads_check)
 
  if(grepl("getOption",  deparse_sub_mcall_threads_check)) {
    mcall_threads_ <- mcall_threads_ 
  } else if(!grepl("threading",  deparse_mcall_threads_check)) {
    if(!is.list(mcall_threads_)) {
      temp_threads_ <- mcall_threads_
      if(is.null(temp_threads_)) {
        temp_threads_ <- temp_threads_
      } else if(is.numeric(temp_threads_)) {
        temp_threads_ <- as.integer(temp_threads_)
      } else if(is.integer(temp_threads_)) {
        temp_threads_ <- temp_threads_
      } else if(is.character(temp_threads_)) {
        temp_threads_ <- threads_char(temp_threads_, 
                                      chains = eval(mcall$chains))
        
      } else {
        
        stop2c("Argument 'threads' must be 'NULL', a string, or an 'integer'")
      } # if(is.null(temp_threads_)) { else if else if....
      mcall_threads_  <- brms::threading()
      if(!is.null(temp_threads_)) {
        mcall_threads_$threads <- temp_threads_
      }
    } # if(!is.list(mcall_threads_)) {
  } # if(grepl("getOption",  deparse_sub_mcall_threads_check)) { else 
  
  mcall$threads <- mcall_threads_
  
  newcall_checks <- c('save_pars')
  if(!is.null(mcall$threads)) {
    if(is.list(mcall$threads)) {
      newcall_checks <- newcall_checks
    } else {
      newcall_checks <- c(newcall_checks, 'threads')
    }
  }
  
  newcall <- check_brms_args(mcall, newcall_checks)
  mcall <- mcall_ <- newcall
  
  # Check and set Alias argument for a b c ... formula
  dots_allias <- list(...)
  collect_dot_names <- c()
  for (ia in letters[1:26]) {
    set_name_dot <- paste0(ia, ".", 'formula')
    set_name_uns <- paste0(ia, "_", 'formula')
    collect_dot_names <- c(collect_dot_names, set_name_dot)
    if (set_name_dot %in% names(dots_allias)) {
      if (eval(bquote(missing(.(set_name_uns)))) ) { 
        mcall[[set_name_uns]] <- dots_allias[[set_name_dot]]
      } else if (!eval(bquote(missing(.(set_name_uns)))) ) { 
        err_msg <- paste0("both '", set_name_uns, "' and '" , 
                          set_name_dot, "' found, ignoring '",
                          set_name_dot, "'")
        if(verbose) warning(err_msg)
      }
    }
  }
  
  # Remove  'collect_dot_names' to avoid conflict with 'brms' dot arguments
  for (collect_dot_namesi in collect_dot_names) {
    if(!is.null(mcall[[collect_dot_namesi]])) 
      mcall[[collect_dot_namesi]] <- NULL
  }
  
  # Check and set Alias argument for d_adjusted (SITAR)
  collect_dot_names <- c()
  for (ia in letters[4]) {
    set_name_dot <- paste0(ia, ".", 'adjusted')
    set_name_uns <- paste0(ia, "_", 'adjusted')
    collect_dot_names <- c(collect_dot_names, set_name_dot)
    if (set_name_dot %in% names(dots_allias)) {
      if (eval(bquote(missing(.(set_name_uns)))) ) { 
        mcall[[set_name_uns]] <- dots_allias[[set_name_dot]]
      } else if (!eval(bquote(missing(.(set_name_uns)))) ) { 
        err_msg <- paste0("both '", set_name_uns, "' and '" , 
                          set_name_dot, "' found, ignoring '",
                          set_name_dot, "'")
        if(verbose) warning(err_msg)
      }
    }
  }
  
  # Remove  'collect_dot_names' to avoid conflict with 'brms' dot arguments
  for (collect_dot_namesi in collect_dot_names) {
    if(!is.null(mcall[[collect_dot_namesi]])) 
      mcall[[collect_dot_namesi]] <- NULL
  }
  
  # Clear alias argument for formula and adjusted
  rm(dots_allias)
  mcall <- mcall_ <- mcall
  
  # Problem with rethinking occurs during 'expose_model_function()'
  if("rethinking" %in% (.packages())){
    message2c("Package 'rethinking' detached and unloaded ato avoid conflict",
            " \nwith the rstan version ", utils::packageVersion('rstan'))
    detach("package:rethinking", unload=TRUE) 
  }
  
  if(utils::packageVersion('rstan') < "2.26") {
    if(expose_function) {
      stop2c("Argument 'expose_function' not supported ",
             "for the rstan version ",
             utils::packageVersion('rstan'))
    }
  }
  
 quote_random_as_init_arg <- function(temp_init_call_in, mcall,...) {
   if(is.null(temp_init_call_in)) {
     temp_init_call_c <- temp_init_call_in
   } 
   if(is.symbol(temp_init_call_in) | is.numeric(temp_init_call_in)) {
     if(!is.character(temp_init_call_in)) {
       temp_init_call_c <- deparse(temp_init_call_in)
     } else {
       temp_init_call_c <- temp_init_call_in
     }
   } else if(is.character(temp_init_call_in)) {
     temp_init_call_c <- temp_init_call_in
   }
   if(is.language(temp_init_call_in)) {
     temp_init_call <- deparse(temp_init_call_in)
     temp_init_call <- gsub("[[:space:]]", "", temp_init_call)
     temp_init_call <- regmatches(temp_init_call, 
                                  gregexpr("(?<=\\().*?(?=\\))", 
                                           temp_init_call, perl=T))[[1]]
     if(length(temp_init_call) != 0) {
       temp_init_call <- strsplit(temp_init_call, ",")[[1]]
       temp_init_call_c <- c()
       for (temp_init_calli in temp_init_call) {
         if(!grepl("\"", temp_init_calli)) {
           temp_init_call2 <- deparse(temp_init_calli)
         } else {
           temp_init_call2 <- temp_init_calli
         }
         temp_init_call_c <- c(temp_init_call_c, temp_init_call2)
       }
       temp_init_call_c <- gsub("[[:space:]]", "", temp_init_call_c)
       
       
       temp_init_call_c <- paste0("list(", 
                                  paste(temp_init_call_c, 
                                        collapse = ",") , ")")
       temp_init_call_c <- str2lang(temp_init_call_c)
     } else if(is.symbol(temp_init_call_in)) { 
       temp_init_call_c <- deparse(substitute(temp_init_call_c))
       temp_init_call_c <- gsub("\"" , "", temp_init_call_c)
     } else {
       temp_init_call_c <- mcall$init
     }
   } 
   temp_init_call_c
 }
  
 
 if(!is.null(mcall$init)) {
   mcall$init <- quote_random_as_init_arg(mcall$init, mcall)
 } else if(is.null(mcall$init)) {
   mcall$init <- "NULL"
 }
 
  for (inxc in letters[1:26]) {
    what_inxc <- paste0(inxc, "_", "init", "_", "beta", "")
    if(!is.null(mcall[[what_inxc]])) mcall[[what_inxc]] <- 
        quote_random_as_init_arg(mcall[[what_inxc]], mcall)
    what_inxc <- paste0(inxc, "_", "cov", "_", "init", "_", "beta", "")
    if(!is.null(mcall[[what_inxc]])) mcall[[what_inxc]] <- 
      quote_random_as_init_arg(mcall[[what_inxc]], mcall)
    what_inxc <- paste0(inxc, "_", "init", "_", "sd", "")
    if(!is.null(mcall[[what_inxc]])) mcall[[what_inxc]] <- 
      quote_random_as_init_arg(mcall[[what_inxc]], mcall)
    what_inxc <- paste0(inxc, "_", "cov", "_", "init", "_", "sd", "")
    if(!is.null(mcall[[what_inxc]])) mcall[[what_inxc]] <- 
      quote_random_as_init_arg(mcall[[what_inxc]], mcall)
  }
  
  what_inxc <- paste0('sigma', "_", "init", "_", "beta", "")
  if(!is.null(mcall[[what_inxc]])) mcall[[what_inxc]] <- 
    quote_random_as_init_arg(mcall[[what_inxc]], mcall)
  
  what_inxc <- paste0('sigma', "_", "cov", "_", "init", "_", "beta", "")
  if(!is.null(mcall[[what_inxc]])) mcall[[what_inxc]] <- 
    quote_random_as_init_arg(mcall[[what_inxc]], mcall)
 
  # Initiate non methods::formalArgs()
  a <- b <- c <- d <- e <- f <- g <- h <- i <- NULL;
  sitar <- NULL;
  mean <- NULL;
  xoffset <- NULL;
  gaussian <- NULL;
  un <- NULL;
  normal <- NULL;
  student_t <- NULL;
  cauchy <- NULL;
  lognormal <- NULL;
  uniform <- NULL;
  exponential <- NULL;
  gamma <- NULL;
  inv_gamma <- NULL;
  lkj <- NULL;
  idsi <- NULL;
  dfsi <- NULL;
  knotssi <- NULL;
  d_formulasi <- NULL;
  ysi <- NULL;
  a_formula_gr_strsi <- NULL;
  b_formula_gr_strsi <- NULL;
  c_formula_gr_strsi <- NULL;
  d_formula_gr_strsi <- NULL;
  e_formula_gr_strsi <- NULL;
  f_formula_gr_strsi <- NULL;
  g_formula_gr_strsi <- NULL;
  h_formula_gr_strsi <- NULL;
  i_formula_gr_strsi <- NULL;
  s_formula_gr_strsi <- NULL;
  d_adjustedsi <- NULL;
  xsi <- NULL;
  xfunsi <- NULL;
  yfunsi <- NULL;
  boundsi <- NULL;
  xoffsetsi <- NULL;
  bstartsi <- NULL;
  cstartsi <- NULL;
  apvsi <- NULL;
  pvsi <- NULL;
  group_arg_groupvar <- NULL;
  multivariate_rescor <- NULL; # why ?
  univariate_by_by <- NULL;
  sigma_arg_groupvar <- NULL;
  a_init_betasi <- NULL;
  b_init_betasi <- NULL;
  c_init_betasi <- NULL;
  d_init_betasi <- NULL;
  e_init_betasi <- NULL;
  f_init_betasi <- NULL;
  g_init_betasi <- NULL;
  h_init_betasi <- NULL;
  i_init_betasi <- NULL;
  s_init_betasi <- NULL;
  a_cov_init_betasi <- NULL;
  b_cov_init_betasi <- NULL;
  c_cov_init_betasi <- NULL;
  d_cov_init_betasi <- NULL;
  e_cov_init_betasi <- NULL;
  f_cov_init_betasi <- NULL;
  g_cov_init_betasi <- NULL;
  h_cov_init_betasi <- NULL;
  i_cov_init_betasi <- NULL;
  s_cov_init_betasi <- NULL;
  a_init_sdsi <- NULL;
  b_init_sdsi <- NULL;
  c_init_sdsi <- NULL;
  d_init_sdsi <- NULL;
  e_init_sdsi <- NULL;
  f_init_sdsi <- NULL;
  g_init_sdsi <- NULL;
  h_init_sdsi <- NULL;
  i_init_sdsi <- NULL;
  s_init_sdsi <- NULL;
  a_cov_init_sdsi <- NULL;
  b_cov_init_sdsi <- NULL;
  c_cov_init_sdsi <- NULL;
  d_cov_init_sdsi <- NULL;
  e_cov_init_sdsi <- NULL;
  f_cov_init_sdsi <- NULL;
  g_cov_init_sdsi <- NULL;
  h_cov_init_sdsi <- NULL;
  i_cov_init_sdsi <- NULL;
  s_cov_init_sdsi <- NULL;
  sigma_init_betasi <- NULL;
  sigma_cov_init_betasi <- NULL;
  sigma_init_sdsi <- NULL;
  sigma_cov_init_sdsi <- NULL;
  rsd_init_sigmasi <- NULL;
  dpar_init_sigmasi <- NULL;
  dpar_cov_init_sigmasi <- NULL;
  autocor_init_acorsi <- NULL;
  autocor_init_unstr_acorsi <- NULL;
  gr_init_corsi <- NULL;
  sigma_init_corsi <- NULL;
  mvr_init_rescorsi <- NULL;
  r_init_zsi <- NULL;
  a_prior_betasi <- NULL;
  b_prior_betasi <- NULL;
  c_prior_betasi <- NULL;
  d_prior_betasi <- NULL;
  e_prior_betasi <- NULL;
  f_prior_betasi <- NULL;
  g_prior_betasi <- NULL;
  h_prior_betasi <- NULL;
  i_prior_betasi <- NULL;
  s_prior_betasi <- NULL;
  a_cov_prior_betasi <- NULL;
  b_cov_prior_betasi <- NULL;
  c_cov_prior_betasi <- NULL;
  d_cov_prior_betasi <- NULL;
  e_cov_prior_betasi <- NULL;
  f_cov_prior_betasi <- NULL;
  g_cov_prior_betasi <- NULL;
  h_cov_prior_betasi <- NULL;
  i_cov_prior_betasi <- NULL;
  s_cov_prior_betasi <- NULL;
  a_prior_sdsi <- NULL;
  b_prior_sdsi <- NULL;
  c_prior_sdsi <- NULL;
  d_prior_sdsi <- NULL;
  e_prior_sdsi <- NULL;
  f_prior_sdsi <- NULL;
  g_prior_sdsi <- NULL;
  h_prior_sdsi <- NULL;
  i_prior_sdsi <- NULL;
  s_prior_sdsi <- NULL;
  a_cov_prior_sdsi <- NULL;
  b_cov_prior_sdsi <- NULL;
  c_cov_prior_sdsi <- NULL;
  d_cov_prior_sdsi <- NULL;
  e_cov_prior_sdsi <- NULL;
  f_cov_prior_sdsi <- NULL;
  g_cov_prior_sdsi <- NULL;
  h_cov_prior_sdsi <- NULL;
  i_cov_prior_sdsi <- NULL;
  s_cov_prior_sdsi <- NULL;
  gr_prior_corsi <- NULL;
  sigma_prior_corsi <- NULL;
  sigma_prior_betasi <- NULL;
  sigma_cov_prior_betasi <- NULL;
  sigma_prior_sdsi <- NULL;
  sigma_cov_prior_sdsi <- NULL;
  rsd_prior_sigmasi <- NULL;
  dpar_prior_sigmasi <- NULL;
  dpar_cov_prior_sigmasi <- NULL;
  autocor_prior_acorsi <- NULL;
  autocor_prior_unstr_acorsi <- NULL;
  mvr_prior_rescorsi <- NULL;
  initsi <- NULL;
  hierarchical_gr_names <- NULL;
  sigma_hierarchical_gr_names <- NULL;
  lb <- NULL;
  ub <- NULL;
  init_rsi <- NULL;
  `:=` <- NULL;
  . <- NULL;
  s_formulasi <- NULL;
  ymean <- NULL;
  ysd <- NULL;
  lm <- NULL;
  checks. <- NULL;
  NoccPI <- NULL;
  NoccAI <- NULL;
  xs <- NULL;
  ids <- NULL;
  dfs <- NULL;
  XXi <- NULL;
  onepic <- NULL;
  temp1 <- NULL;
  temp2 <- NULL;
  sigma_formula_manualsi <- NULL;
  sigmaxsi <- NULL;
  sigmadfsi <- NULL;
  sigmaknotssi <- NULL;
  sigmafixedsi <- NULL;
  sigmarandomsi <- NULL;
  sigmaxoffsetsi <- NULL;
  sigmaxfunsi <- NULL;
  sigmayfunsi <- NULL;
  sigmaboundsi <- NULL;
  sigmabstartsi <- NULL;
  sigmacstartsi <- NULL;
  sigmaids <- NULL;
  sigmaxoffset <- NULL;
  sigmadfs <- NULL;
  ixfuntransformsi <- NULL;
  iyfuntransformsi <- NULL;
  sigmaixfuntransformsi <- NULL;
  nsp <- NULL;
  nsk <- NULL;
  rcs <- NULL;
  sigma_formulasi <- NULL;
  sigma_formula_grsi <- NULL;
  sigma_formula_gr_strsi <- NULL;
  dpar <- NULL;
  nlpar <- NULL;
  sigma_formula_manual_prior_via_sigma_formula <- NULL;
  
  ##############################################################
  # cp -> multi_normal_cholesky_lpdf / multi_normal_lpdf
  ##############################################################
  
  if(is.null(parameterization)) {
    parameterization <- 'cp'
  }
  
  if(parameterization == 'cp') {
    cp_via <- "multi_normal_cholesky_lpdf"
  } else if(parameterization == 'multi_normal_cholesky_lpdf' |
            parameterization == 'multi_normal_cholesky' |
            parameterization == 'normal_cholesky' |
            parameterization == 'normal_chol') {
    cp_via <- "multi_normal_cholesky_lpdf"
    parameterization <- 'cp'
  } else if(parameterization == 'multi_normal_lpdf' |
            parameterization == 'multi_normal' |
            parameterization == 'normal') {
    cp_via <- "multi_normal_lpdf"
    parameterization <- 'cp'
  } else if(parameterization == 'ncp') {
    cp_via <- NULL
  } else {
    cp_via <- NULL
  }
  
  ##############################################################
  ##############################################################
  
  enclose_c_list_elemnts_with_quotes <- function(xcall, elements) {
    xcall_x_c <- list()
    for (x in elements) {
      xcall_x <- deparse(xcall[[x]])
      xcall_x <- paste0(xcall_x, collapse = " ")
      xcall_x <- gsub_space(xcall_x)
      xcall_x <- gsub("\"", "", xcall_x)
      str_c       <- FALSE
      str_list    <- FALSE
      if(grepl("c\\(", xcall_x, fixed = F)) {
        pater_add <- "c("
        str_c     <- TRUE
      } else if(grepl("list\\(",xcall_x, fixed = F)) {
        pater_add <- "list("
        str_list  <- TRUE
      } else {
        # stop2c(x ," must be a list or a vector")
      }
      if(str_c | str_list) {
        xcall_x_str_get <- xcall_x
        if(grepl("\\(", xcall_x_str_get)) {
          xcall_x_str_get <- regmatches(xcall_x, gregexpr("(?<=\\().*(?=\\))",
                                                          xcall_x, perl=T))[[1]]
        } else {
          xcall_x_str_get <- xcall_x
        }
      } else if(!str_c & !str_list) {
        xcall_x_str_get <- xcall_x
      }
      xcall_x_str_get <- strsplit(xcall_x_str_get, ",")[[1]]
      xcall_x_i_str_c <- list()
      xcall_x_i_str_i_get <- list()
      for (i in 1:length(xcall_x_str_get)) {
        xcall_x_i_str_i <- xcall_x_str_get[i]
        if(grepl("^function\\(", xcall_x_i_str_i)) {
          xcall_x_i_str_i_get <- str2lang(xcall_x_i_str_i)
          xcall_x_i_str_c[[i]] <- xcall_x_i_str_i_get
        } else {
          if(grepl("\\(", xcall_x_i_str_i)) {
            xcall_x_i_str_i <- regmatches(xcall_x_i_str_i, 
                                          gregexpr("(?<=\\().*(?=\\))",
                                                   xcall_x_i_str_i, perl=T))[[1]]
          } else {
            xcall_x_i_str_i <- xcall_x_i_str_i
          }
          xcall_x_i_str_i <- paste0(sprintf("'%s'", xcall_x_i_str_i), 
                                    collapse = ",")
          
          xcall_x_i_str_i_get <-  str2lang(xcall_x_i_str_i)
          
          
          xcall_x_i_str_c[[i]] <- xcall_x_i_str_i_get
        }
      } # for (i in length(xcall_x_i_str)) {
      if(!is.null(xcall_x_i_str_c)) {
        xcall[[x]] <- xcall_x_i_str_c
      }
    } # for (x in elements) {
    return(xcall) 
  } # enclose_c_list_elemnts_with_quotes
  
  ##############################################################
  ##############################################################

 # Enclose primitive functions with quote "" if specified as c() or list()
  enclose_c_list_elemnts_with_quotes_these <- c("xfun",
                                                "yfun",
                                                "sigmaxfun",
                                                "xfunxoffset",
                                                "sigmaxfunxoffset")
  
 
  mcall <- enclose_c_list_elemnts_with_quotes(mcall,
                                       enclose_c_list_elemnts_with_quotes_these)

  
  ##############################################################
  ##############################################################
  
  # for terms_rhs
  quote_elements <- function(call, element,
                             return_call = TRUE, 
                             strict_list = FALSE) {
    expr <- call[[element]]
  
    if (is.symbol(expr)) {
      expr <- get(as.character(expr), envir = parent.frame())
  
      if (!is.language(expr)) {
        stop2c("Symbol must refer to a language object (e.g., quote(c(...)))")
      }
    }
  
    if(strict_list) {
      allowed_lang <- 'list'
    } else {
      allowed_lang <- c("c", "list")
    }
    
    quoted <- deparse(substitute(expr)) 
    
    if(grepl("^list\\(", quoted, fixed = F) | 
       grepl("^c\\(", quoted, fixed = F)) {
      if(strict_list) {
        if(grepl("^c\\(", quoted, fixed = F)) {
          stop2c("Argument ",  collapse_comma(element), " must be a list(...)")
        } # if(grepl("^c\\(", quoted, fixed = F)) {
      } # if(strict_list) {
    } else {
      quoted <- deparse(substitute(expr))
      quoted <- gsub("\"", "", quoted)
      return(quoted)
    }
  
    constructor <- expr[[1]]
    args <- as.list(expr)[-1]
  
    quoted <- lapply(args, function(arg) {
      if (is.character(arg)) {
        # Already a string — quote it once
        sprintf("'%s'", arg)
      } else {
        # Expression — deparse it into code, then quote it
        sprintf("'%s'", paste(deparse(arg), collapse = ""))
      }
    })
    quoted <- gsub("\'", "", quoted)
    if(return_call) {
      # Reconstruct the call (as a call to c(...) or list(...))
      out <- as.call(c(constructor, quoted))
    } else {
      out <- quoted
    }
    return(out)
  }
  
  
  mcall$terms_rhs <- quote_elements(call = mcall, 
                                    element = 'terms_rhs', 
                                    return_call = TRUE, 
                                    strict_list = TRUE)
 
  ##############################################################
  ##############################################################
  
  get_pos_val_unnamed_list_element <- function(my_list) {
    unnamed_elements_info <- list()
    for (i in seq_along(my_list)) {
      if (is.null(names(my_list)[i]) || names(my_list)[i] == "") {
        unnamed_elements_info[[length(unnamed_elements_info) + 1]] <- list(
          position = i,
          value = my_list[[i]]
        )
      }
    }
    collect_pos_val <- c()
    for (i in 1:length(unnamed_elements_info)) {
      pos <- unnamed_elements_info[[i]][[1]]
      val <- unnamed_elements_info[[i]][[2]]
      pos_val <- paste0("position = ", pos, "; value = ", val)
      collect_pos_val <- c(collect_pos_val, pos_val)
    }
    collect_pos_val <- paste(collect_pos_val, collapse = "\n")
    return(collect_pos_val)
  }
  
  check_each_element_of_list_named <- function(call, 
                                               element, 
                                               assert_names = NULL) {
    x <- call[[element]] %>% eval()
    
    if(is.null(x)) {
      return(NULL)
    }
    
    if(!is.list(x)) {
      stop2c("'x' must be a list")
    }
    
    x_name    <- element
    all_names <- remove_empty_string_from_vector(x)
    out       <- length(all_names) == length(x)
    
    if(!is.null(assert_names)) {
      names_mismatch <- setdiff(assert_names, all_names)
    }
    
    n_missing_names <- length(x) - length(all_names) 
    
    if(!out) {
      stop2c("Argument ", collapse_comma(x_name), " must be a named list with ",
           "\n  ",
           "each element having a unique name.", 
           "\n  ",
           "Allowed names are: ", 
           # "\n  ",
             collapse_comma(assert_names), 
           "\n  ",
           "The position and value of elements without name(s) are: ", 
           "\n  ",
           get_pos_val_unnamed_list_element(x))
    }
    
    if(is.null(x[['select']]))   x[['select']]   <- 'knots'
    if(is.null(x[['when']]))     x[['when']]     <- 'bc'
    if(is.null(x[['criteria']])) x[['criteria']] <- 'AIC'
    if(is.null(x[['nsearch']]))  x[['nsearch']]  <- NULL
    
    if(is.null(x[['what']])) {
      x[['what']]   <- 'none'
    } else {
      if(isFALSE(x[['what']])) {
        x[['what']]   <- 'none'
      }
    }
    
    if(is.null(x[['return']])) {
      x[['return']]   <- FALSE
    } 
    
    if(is.null(x[['print']])) {
      x[['print']]   <- FALSE
    } 
    
    if(is.null(x[['bkrange']])) {
      x[['bkrange']]   <- FALSE
    } 
    
    if(is.null(x[['fix_bknots']])) {
      x[['fix_bknots']]   <- TRUE
    } 
    
    if(is.null(x[['method']])) {
      x[['method']]   <- 'bs'
    } else if(x[['method']] == 'rs') {
      stop2c("Only 'bs' is allowed as a method for knots_selection")
    }
    
    if(is.null(x[['all_scores']])) {
      x[['all_scores']]   <- FALSE
    } 
    
    if(is.null(x[['plot_all_scores']])) {
      x[['plot_all_scores']]   <- FALSE
    } 
    
    if(is.null(x[['kspace']])) {
      x[['kspace']]   <- 'un'
    } 
    
    if(is.null(x[['cvk']])) {
      x[['cvk']]   <- 10
    } 
    
    if(is.null(x[['cviter']])) {
      x[['cviter']]   <- 100
    } 
    
    if(out) {
      if(!is.null(assert_names)) {
        if(!is.character(assert_names)) {
          stop2c("'assert_names' must be a character or character vector")
        } else if(length(all_names) != length(assert_names)) {
          # 
        } else {
          # names_mismatch <- setdiff(assert_names, all_names)
          if(!is_emptyx(names_mismatch)) {
            stop2("mismatch in names. Following name is missing: ", 
                  collapse_comma(names_mismatch))
          }
        }
      }
    }
    
    if(!grepl('stats::', x[['criteria']])) {
      if(x[['criteria']] == "AIC" | x[['criteria']] == "BIC") {
        x[['criteria']] <- paste0('stats::', x[['criteria']])
      } else if(x[['criteria']] == "CV") {
        x[['criteria']] <- deparse(x[['criteria']])
      }
    }  
    
    return(x)
  }
  
  
  knots_selection_assert_names <- c('select', 
                                   'criteria', 
                                   'nsearch',
                                   'bkrange',
                                   'fix_bknots',
                                   'all_scores',
                                   'plot_all_scores',
                                   'cvk',
                                   'cviter',
                                   'kspace',
                                   'method',
                                   'when', 
                                   'what',
                                   'print',
                                   'return')
 
  
  if(!is.null(mcall[['knots_selection']])) {
    mcall[['knots_selection']] <- knots_selection <- 
      check_each_element_of_list_named(call = mcall, 
                                       element = 'knots_selection', 
                                       assert_names = 
                                         knots_selection_assert_names)
  }
  

  ##############################################################
  ##############################################################
  
  enverr. <- environment()
  for (i in names(mcall)[-1]) {
    no_default_args_plus_family <- c(no_default_args, "family")
    if (!i %in% no_default_args_plus_family) {
      assign('err.', FALSE, envir = enverr.)
      tryCatch(
        expr = {
          if (is.function(suppressWarnings(eval(mcall[[i]])))) {
            checks. <- deparse_0(mcall[[i]])
          } else {
            suppressWarnings(checks. <- eval(mcall[[i]]))
          }
        },
        error = function(e) {
          assign('err.', TRUE, envir = enverr.)
        }
      )
      err. <- get('err.', envir = enverr.)
      if(length(checks.) == 0) err. <- TRUE
      if (err.) {
        mcall[[i]] <- mcall[[i]]
      } else if (!err.) {
        if (is.list(checks.)) {
          if (is.list(checks.[[1]])) {
            mcall[[i]] <- checks.[[1]]
          } else if (!is.list(checks.[[1]])) {
            if (is.list(checks.)) {
              if (is.symbol(mcall[[i]]))
                mcall[[i]] <- deparse_0(mcall[[i]]) # for set_self_priors
              # suppressWarnings 14 01 2024
                suppressWarnings(mcall[[i]] <- eval(mcall[[i]]))
                temp       <- str2lang(deparse_0((mcall[[i]])))
                mcall[[i]] <- temp
            } else if (!is.list(checks.)) {
              mcall[[i]] <- checks.
            }
          }
        } else {
          mcall[[i]] <-  checks.
        }
      }
    } else if (i %in% no_default_args_plus_family) {
      mcall[[i]] <-  mcall[[i]]
    }
  }
  
  
  arguments <- as.list(mcall)[-1]

  match.call.defaults <- function(...) {
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))
    for(i in setdiff(names(formals), names(call)))
      call[i] <- list( formals[[i]] )
    match.call(sys.function(sys.parent()), call)
  }
  
  call.full <- match.call.defaults()
  # 25.01.2026
  if(set_eval_globals_in_mcall) {
    set_eval_globals_in_mcall_names <- 
      remove_empty_string_from_vector(set_eval_globals_in_mcall_names)
    set_exceptions <- setdiff(names(call.full), 
                              set_eval_globals_in_mcall_names)
    set_exceptions <- c(set_exceptions, "data", "...")
    call.full <- eval_globals_in_mcall(call.full,
                                       exceptions = set_exceptions) 
  }
  
  
  # call.full <- match.call.defaults()
  call.full <- call.full[-length(call.full)]
 
  for (call.fulli in names(call.full)) {
    if(call.fulli != "") {
      if(call.fulli == 'family' & 
         is.language(call.full[[call.fulli]])) {
        call.full[[call.fulli]] <- deparse(call.full[[call.fulli]])
      } else if(call.fulli == 'stan_model_args')  { 
        if(length(eval(call.full[[call.fulli]])) == 0) {
          call.full[[call.fulli]] <- NULL
        } else {
          call.full[[call.fulli]] <- call.full[[call.fulli]] 
        }
      } else {
        
      }
    } else {
      #
    }
  }
  

  f_funx_arg <- formals(bsitar)
  nf_funx_arg_names <-
    intersect(names(arguments), names(f_funx_arg))
  arguments <-
    c(arguments, f_funx_arg[names(f_funx_arg) %!in% nf_funx_arg_names])
  
  # New -> don't get default when user specified init/init_r
  arguments$init <- mcall$init
  arguments$init_r <- mcall$init_r
  
  
  
  # add_sigma_by_mu - add_sigma_by_ls
  # 'sigma_formula_manual' is used to set sigma by mu ('add_sigma_by_mu') and 
  # location scale mode ('add_sigma_by_ls'). For 'add_sigma_by_mu' no need to add
  # the following substitute and deparse but are must for 'add_sigma_by_ls'
  # Even if we do when 'add_sigma_by_mu', it reamins unaffected.
  # Anothe important difference is that for 'add_sigma_by_mu', the setdepar0sgub
  # should not include ; but for 'add_sigma_by_ls', it must  
  # have 'sigma_formula_manual'. But if we add list(), then can be excluded.
  
  # Also, if(is.language(sigma_formula_manual_fun)) check is mandatory only for'
  # 'add_sigma_by_ls'

  sigma_formula_manual_fun     <- substitute(sigma_formula_manual)
  sigma_formula_manual_fun_str <- deparse(sigma_formula_manual_fun)
  sigma_formula_manual_fun_str <-paste(gsub_space(sigma_formula_manual_fun_str),
                                        collapse = "")
  
  count_number_nlf <- gregexpr("brms::nlf\\(|nlf\\(", 
                               sigma_formula_manual_fun_str)[[1]][1]

  if(count_number_nlf > 1) {
    if(is.language(sigma_formula_manual_fun)) {
      if(!grepl("^list\\(", sigma_formula_manual_fun_str)) {
        stop2c("Argument 'sigma_formula_manual' must be a list or a string")
      }
    }
  }
  
  sigma_formula_manual <- sigma_formula_manual_fun_str
  
  
  setdepar0sgub <- c("sigma_formula", "sigma_formula_gr")
  if(count_number_nlf > 1) {
    setdepar0sgub <- setdepar0sgub
  } else {
    setdepar0sgub <- c(setdepar0sgub, "sigma_formula_manual")
  }
  
 
  
  for (argumentsnamesi in names(arguments)) {
    if(argumentsnamesi %in% setdepar0sgub) {
      if(!is.null(arguments[[argumentsnamesi]])) {
        if(grepl("^list", arguments[[argumentsnamesi]] )[1]) {
          argumentsnamesi_li <- as.list(arguments[[argumentsnamesi]])
          argumentsnamesi_li <- argumentsnamesi_li[-1]
          argumentsnamesi_c_list <- list()
          for (ixz in 1:length(argumentsnamesi_li)) {
            argumentsnamesi_c_list[[ixz]] <- argumentsnamesi_li[[ixz]] %>% 
              deparse_0() %>% gsub_space()
          }
          arguments[[argumentsnamesi]] <- argumentsnamesi_c_list %>% unlist()
        } else {
          arguments[[argumentsnamesi]] <- arguments[[argumentsnamesi]] %>% 
            deparse_0() %>% gsub_space()
        }
      }
    }
  }
  

  familyzzzx <- arguments$family
  if(grepl("^c\\(", deparse_0(familyzzzx), fixed = FALSE)) {
    stop2c("Argument family should be a list() and not a vector 'c()'")
  } else if(grepl("^\\(", deparse_0(familyzzzx), fixed = FALSE)) {
    familyzzzx <- paste0("list", "", deparse_0(familyzzzx) , "")
    familyzzzx <- str2lang(familyzzzx)
  } else if(!is.list(familyzzzx) & !grepl("list", familyzzzx)[[1]]) {
    familyzzzx <- paste0("list", "(", deparse_0(familyzzzx) , ")")
    familyzzzx <- str2lang(familyzzzx)
  } else {
    familyzzzx <- familyzzzx 
  }
  familyzzzx2 <- familyzzzx
  arguments$family <- familyzzzx
  

  checks_start_names <- c('bstart', 'cstart', 'apv', 'pv')
  for (checks_start_namesi in checks_start_names) {
    if(checks_start_namesi %in% names(mcall_)) {
      if(is.null(mcall_[[checks_start_namesi]])) {
        arguments[[checks_start_namesi]] <- 'NULL'
      }
    }
  }
  
  
  checks_start_names <- c('sigmabstart', 'sigmacstart')
  for (checks_start_namesi in checks_start_names) {
    if(checks_start_namesi %in% names(mcall_)) {
      if(is.null(mcall_[[checks_start_namesi]])) {
        arguments[[checks_start_namesi]] <- 'NULL'
      }
    }
  }
  
  
  # Override when restricting to abcd
  override_select_model <- TRUE # FALSE
  if(override_select_model) arguments$select_model <- select_model <- 'sitar'
  
  

  # 01.07/2025
  getdotslist <- list(...)
  
# getdotslist decomp
# QR_Xmat -> not important. Internally set as Qc
# QR_center -> center spl matrix before applying QR decom
# QR_complete -> whether to get complete QR matrix. Don't change, must be FALSE
# QR_flip -> flip negative to positive of Q R matrices. This matches R with Stan
# QR_scale -> scaling, default sqrt(N-1) when QR_scale = NULL
# qr_gq -> controls wheter to add v_sx vectors to gen quant
  
  QR_decomp_args <- getdotslist[['decomp']]
  if(is.null(QR_decomp_args)) {
    decomp <- NULL
    QR_Xmat      <- NULL
    QR_center    <- NULL
    QR_complete  <- NULL
    QR_flip      <- NULL
    QR_scale     <- NULL
    QR_gq        <- FALSE
  } else if(!is.null(QR_decomp_args)) {
    if(is.list(QR_decomp_args)) {
      if(is.null(QR_decomp_args[['decomp']])) {
        decomp      <- 'QR'
      } else {
        decomp      <- QR_decomp_args[['decomp']]
      }
      if(is.null(QR_decomp_args[['Xmat']])) {
        QR_Xmat      <- 'Qc'
      } else {
        QR_Xmat      <- QR_decomp_args[['Xmat']]
      }
      if(is.null(QR_decomp_args[['center']])) {
        QR_center      <- FALSE
      } else {
        QR_center      <- QR_decomp_args[['center']]
      }
      if(is.null(QR_decomp_args[['complete']])) {
        QR_complete      <- FALSE
      } else {
        QR_complete      <- QR_decomp_args[['complete']]
      }
      if(is.null(QR_decomp_args[['flip']])) {
        QR_flip      <- TRUE
      } else {
        QR_flip      <- QR_decomp_args[['flip']]
      }
      if(is.null(QR_decomp_args[['scale']])) {
        QR_scale      <- "sqrt(N-1)"
      } else {
        QR_scale      <- QR_decomp_args[['scale']]
      }
      if(is.null(QR_decomp_args[['gq']])) {
        QR_gq      <- FALSE
      } else {
        QR_gq      <- QR_decomp_args[['gq']]
      }
    } else if(!is.list(QR_decomp_args)) {
      if(!is.character(QR_decomp_args)) {
        stop2c("argument 'decomp' must be a list or character")
      } else {
        if(QR_decomp_args != "QR") {
          stop2c("only 'QR' decomp allowed")
        }
        decomp       <- QR_decomp_args
        QR_Xmat      <- 'Qc'
        QR_center    <- FALSE
        QR_complete  <- FALSE
        QR_flip      <- TRUE
        QR_scale     <- NULL
        QR_gq        <- FALSE
      }
    } # else if(!is.list(QR_decomp_args)) {
  } # else if(!is.null(QR_decomp_args)) {
  
  
  add_rcsfunmatqrinv_genquant <- QR_gq
  
 
  if(is.character(arguments$select_model)) {
    select_model <- arguments$select_model
  } else if(is.symbol(arguments$select_model)) {
    select_model <- deparse(arguments$select_model)
  } else if(!is.character(arguments$select_model) |
            !is.symbol(arguments$select_model)
            ) {
    stop2c("The argument 'select_model' must be a symbol or 
         single character string")
  }
  
  # For editing scode (if required for later use TODO)
  select_model_edit <- select_model
  
  if(select_model == 'logistic1e') select_model <- 'logistic1'
  if(select_model == 'logistic2e') select_model <- 'logistic2'
  if(select_model == 'logistic3e') select_model <- 'logistic3'
  
  # For default prior setting (if required for later use TODO)
  select_model_arg <- select_model
  
  
  if(select_model == 'pb')       select_model <- 'pb1'
  if(select_model == 'logistic') select_model <- 'logistic1'
  if(select_model == 'sitar3')   select_model <- 'sitar'
  
  if(select_model == 'rcs' | 
     select_model == 'rcsf' | 
     select_model == 'rcsfr') {
    if(select_model == 'rcsf') {
      rcs_add_re_spline <- FALSE 
    } else {
      rcs_add_re_spline <- TRUE
    }
    select_model <- 'rcs'
  }
  
  # For ns() based SITAR, a intercept is matched if rcs based s1 is adjusted as
  # A=a-(s1*min(knots))
  # We keeping same form for mu (match_sitar_a_form = TRUE) but not for sigma
  # Note below that these can be controlled via ... dots
  
  # 01.07/2025 -> getdotslist moved up to allow passing decomp 
  # 24.08.2024
  # getdotslist <- list(...)
  
  # Note that for all argumnets such as intercept to correctly pass to the 
  # use smat and not stype
  
  # spline types supported are 'rcs', 'nsp' and 'nsk', 'bsp', 'isp', 'msp'
  # The argument . is exposed that allows setting spline type as string 
  # However for developmental purposes, an additional option is allowed that 
  # pass a named list ('smat') via ... that allows a more elaborate control on 
  # various aspects of splines. These are currently tested and will be exposed
  # later. These are 
  # 1. type - a character tring to set spline type i.e, rcs, nsp and nsk
  # 2. normalize - a logical (T/F) to specify to normalize H matrix
  # 3. centerval -a real number to center the intercept a given value
  # 4. intercept -a logical (T/F) to specify whether or not to return matrix 
  # with intercept
  # 5. preH - a logical (T/F) to use precomputed H matrix to compute within 
  # the function 
  # 6. include - a logical (T/F) to indicate if .stan splines be included 
  # via '#include' or 
  # read it and include as it is in the function block 
  # 7. sfirst - a logical (T/F) to indicate whether to include only the first
  # elements from s1, s2, s2,.. vectors in functions block e.g., s1[1], s2[1]
  # This sfirst should be used only when splines have no covariate or random effect
  # This typically is the case of sitar model
  # 8. sparse - a logical (T/F) to indicate whether to use sparse in the function
  # block where Spl * s vector is used. For this sparse, sfirst need to be T.
  # Again, this is helpful in case of sitar model.
  # Note that according to stan documentation, sparse = T speeds up the computation
  # only when sparsity > 90 %. The sparcity in sitar model is '0'% except when 
  # type = 'rcs' and decomp = NULL.
  # 9. check_sparsity - a logical (T/F) to check the sparsity (%) in the 
  # function block where Spl * s vector is used. For this sparse, both 'sfirst'
  # 'sparse' need to be T, also chains = 1, and iter = 2. (since this is printed)
  # Again, this is helpful in case of sitar model
  # 10. degree - an integers - for rcs, then 'nk' is df + 1
  # 11. bkrange - 
  #     should 'rcs' boundary knots set using xrange (TRUE) or quntile (FALSE). 
  # 12. fix_bknots 
  #     should 'rcs' boundary knots be fixed (TRUE) from qunatile or not (FALSE)  
  # 13. return - 
  #     should the plot showing the knots returned (TRUE) or not (FALSE, default). 
  # 14. print 
  #     should the plot showing the knots returned (TRUE) or not (FALSE, default). 
  # 15. what 
  #     what to print i.e, which plot. see knots_selection argument  
  # 16. when 
  #     when to print / return, before centering 'bc', or after centering 'ac'
  
  # Although stype has when option, somehow it does not work
  # Therefore, let the default 'bc' keep working ....
  # if(is.null(knots_selection)) {
  #   if(smat_when == "bc") {
  
  
  allowed_spline_type <- c('rcs', 'nsp', 'nsk', 'bsp', 'msp', 'isp', 'moi')
  allowed_spline_type_exception_msg <- 
    paste("The options available are:", 
          paste(paste(paste0("'", allowed_spline_type, "'"), collapse =", "), 
                collapse =", ")
    )
  
  # 5.06.2025 -> this needed for CustomDoCall in update_model
  quote_allowed_spline_type <- function(aaax, allowed_spline_type) {
    for (ix in allowed_spline_type) {
      gsub_it <- ix
      gsub_by <- paste0("'", ix, "'")
      aaax <- gsub(gsub_it, gsub_by, aaax, fixed = TRUE)
      if(grepl("\"", aaax, fixed = TRUE)) {
        aaax <- gsub("\"", "", aaax, fixed = TRUE)
      }
    }
    return(aaax)
  }
  
  
  # For CustomDoCall
  if(!is.null(stype)) {
     if(is.symbol(stype)) stype <- deparse(stype)
  }

  stype_temp_str <- deparse(substitute(stype))
  stype_temp_str <- paste0(gsub_space(stype_temp_str), collapse = " ")
  
  stype_temp_str <- gsub_quote1(stype_temp_str)
  
  # Handle list[[stype]]
  if(grepl("\\[\\[", stype_temp_str)) {
    stype_temp_str <- stype
  } else if(grepl("\\$", stype_temp_str)) {
    stype_temp_str <- stype
  } else {
    stype_temp_str <- stype_temp_str
  }
  
  if(grepl("^list\\(", stype_temp_str)) {
    stype_temp_str <-  quote_allowed_spline_type(stype_temp_str, 
                                                 allowed_spline_type)
    stype <- ept(stype_temp_str)
    stype[['type']] <- stype[[1]]
  } else {
    stype <- stype_temp_str
    stype <- gsub("\"", "", stype)
  }
  
  allowed_smat_options <- c('type', 
                            'centerval',
                            'intercept',
                            'degree',
                            'normalize', 
                            'derivs', 
                            'preH', 
                            'include',
                            "sfirst",
                            "sparse", 
                            "check_sparsity", 
                            "bkrange",
                            "fix_bknots",
                            "what",
                            'when', 
                            "return",
                            "print")
  
  # Note that to pass all arguments correctly, use 'smat' and not 'stype'
  # later, after testing, the 'smat' will be inferred from 'stype'
  spline_type_via_stype <- FALSE
  if(!is.null(getdotslist[['smat']])) {
    spline_type <- getdotslist[['smat']]
    if(!is.list(spline_type)) {
      stop2c("'smat' set via '...' must be a named list")
    }
    checknamessmat <- names(spline_type)
    checknamessmat <- checknamessmat[nzchar(checknamessmat)] 
    if(length(checknamessmat) != length(spline_type)) {
      stop2c("Each element of 'smat' set via '...' must be named.",
           "\n  ", 
           "Allowed options are:",
           "\n  ", 
           collapse_comma(allowed_smat_options))
    } 
    for (checknamessmati in checknamessmat) {
      if(!checknamessmati %in% allowed_smat_options)
        stop2c("Option  '", checknamessmati, "' is invalid for 'smat'",
             "\n  ", 
             "Allowed options are:",
             "\n  ", 
             collapse_comma(allowed_smat_options))
    }
  } else {
    spline_type <- stype
    spline_type_via_stype <- TRUE
  } 
  
  if(any(spline_type == "NULL")) {
    spline_type_via_stype <- FALSE
  }
  
  if(is.null(getdotslist[['smat']]) & !spline_type_via_stype) {
    spline_type <- 'rcs'
    # if(verbose) message2c("'rcs' set as default spline type")
  }
    
  # While testing, expose only type and normalize from the 'stype' 
  allowed_spline_type_list_names_c <- c('type', 
                                        'normalize')
  
  allowed_spline_type_list_names_msg <- 
    paste("argument 'spline_type' must be a named list, allowed names are:\n", 
          paste(paste(paste0("'", allowed_spline_type_list_names_c, "'"), 
                      collapse =", "), 
                collapse =", ")
    )
    
  spline_type_list <- list()
  if(is.null(spline_type)) {
    spline_type_list[['type']]        <- NULL
    spline_type_list[['centerval']]   <- 0
    spline_type_list[['degree']]      <- 3L
    spline_type_list[['intercept']]   <- FALSE
    spline_type_list[['normalize']]   <- TRUE
    spline_type_list[['derivs']]      <- FALSE
    spline_type_list[['preH']]        <- TRUE
    spline_type_list[['include']]     <- TRUE
    spline_type_list[['sfirst']]      <- FALSE
    spline_type_list[['bkrange']]     <- TRUE
    spline_type_list[['fix_bknots']]  <- TRUE
    spline_type_list[['what']]        <- NULL
    spline_type_list[['when']]        <- NULL
    spline_type_list[['return']]      <- FALSE
    spline_type_list[['print']]       <- FALSE
    spline_type_list[['sparse']]      <- FALSE
    spline_type_list[['check_sparsity']]      <- FALSE
  } else if(!is.null(spline_type)) {
    if(is.list(spline_type)) {
      if(is.null(spline_type[['normalize']])) {
        spline_type[['normalize']]   <- TRUE
      } 
      if(is.null(spline_type[['preH']])) {
        spline_type[['preH']]   <- TRUE
      } 
      if(is.null(spline_type[['sfirst']])) {
        spline_type[['sfirst']]   <- FALSE
      } 
      if(is.null(spline_type[['sparse']])) {
        spline_type[['sparse']]   <- FALSE
      } 
      if(is.null(spline_type[['check_sparsity']])) {
        spline_type[['check_sparsity']]   <- FALSE
      } 
      if(is.null(spline_type[['bkrange']])) {
        spline_type[['bkrange']]   <- TRUE
      } 
      if(is.null(spline_type[['fix_bknots']])) {
        spline_type[['fix_bknots']]   <- TRUE
      } 
      if(is.null(spline_type[['what']])) {
        spline_type[['what']]   <- NULL
      } 
      if(is.null(spline_type[['when']])) {
        spline_type[['when']]   <- NULL
      } 
      if(is.null(spline_type[['return']])) {
        spline_type[['return']]   <- FALSE
      } 
      if(is.null(spline_type[['print']])) {
        spline_type[['print']]   <- FALSE
      } 
      if(length(spline_type) > 0) {
        # if only type specified and unnamed, name it
        if(is.null(names(spline_type))) { 
          if(length(spline_type) == 1) {
            names(spline_type) <- 'type'
          } else if(length(spline_type) == 2) {
            # if only type and normalize specified and unnamed, name them
            # This is only when spline type set via stype argument and not ...
            if(spline_type_via_stype) {
              names(spline_type) <- c('type', 'normalize')
              if(verbose) {
                message2c("stype arguments named as 'type', 'normalize'")
              }
            } else {
              stop2c(allowed_spline_type_list_names_msg)
            }
          } else if(length(spline_type) == 3) {
            if(spline_type_via_stype) {
              names(spline_type) <- c('type', 'normalize', "preH")
              if(verbose) {
                message2c("stype arguments named as 'type', 
                          'normalize', 'preH'")
              }
            } else {
              stop2c(allowed_spline_type_list_names_msg)
            }
          } else {
            stop2c(allowed_spline_type_list_names_msg)
          }
        }
        if(!is.null(spline_type[['type']])) {
          if(!is.character(spline_type[['type']])) {
            stop2c(paste0(spline_type[['type']], 
                          " must be a character string"))
          } else {
            spline_type_list[['type']] <- spline_type[['type']]
          }
        } else if(is.null(spline_type[['type']])) {
          # 
        }
        
        # change check message same as 'centerval' for other
        if(!is.null(spline_type[['degree']])) {
          if(!is.numeric(spline_type[['degree']])) {
            stop2c("Argument 'degree' must be a numeric value",
                 " but instead specified as ", 
                 "'", paste0(spline_type[['degree']], "'"))
          } else {
            spline_type_list[['degree']] <- spline_type[['degree']]
          }
        } else if(is.null(spline_type[['degree']])) {
          spline_type_list[['degree']] <- 3
        }
        
        
        
        
        if(!is.null(spline_type[['intercept']])) {
          if(!is.logical(as.logical(spline_type[['intercept']]))) {
            stop2c(paste0(spline_type[['intercept']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['intercept']] <- spline_type[['intercept']]
          }
        } else if(is.null(spline_type[['intercept']])) {
          spline_type_list[['intercept']] <- FALSE
        }
        
        if(!is.null(spline_type[['normalize']])) {
          if(!is.logical(as.logical(spline_type[['normalize']]))) {
            stop2c(paste0(spline_type[['normalize']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['normalize']] <- spline_type[['normalize']]
          }
        } else if(is.null(spline_type[['normalize']])) {
          spline_type_list[['normalize']] <- FALSE
          if(verbose) message2c(paste0("'", FALSE ,
                                     "' set as default spline normalize"))
        }
        
        if(!is.null(spline_type[['derivs']])) {
          if(!is.integer(as.logical(spline_type[['derivs']]))) {
            stop2c("Argument 'derivs' must be an integer",
                 " but instead specified as ", 
                 "'", paste0(spline_type[['degree']], "'"))
          } else {
            spline_type_list[['derivs']] <-  spline_type[['derivs']] 
          }
        } else if(is.null(spline_type[['derivs']])) {
          spline_type_list[['derivs']]    <- 0
        }
        
        if(!is.null(spline_type[['preH']])) {
          if(!is.logical(as.logical(spline_type[['preH']]))) {
            stop2c(paste0(spline_type[['preH']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['preH']]    <-  spline_type[['preH']] 
          }
        } else if(is.null(spline_type[['preH']])) {
          spline_type_list[['preH']]    <- FALSE
        }
        
        if(!is.null(spline_type[['include']])) {
          if(!is.logical(as.logical(spline_type[['include']]))) {
            stop2c(paste0(spline_type[['include']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['include']]    <-  spline_type[['include']] 
          }
        } else if(is.null(spline_type[['include']])) {
          spline_type_list[['include']]    <- TRUE
        }
        
        if(!is.null(spline_type[['path']])) {
          if(!is.character(spline_type[['path']])) {
            stop2c(paste0(spline_type[['path']]," must be a character string"))
          } else {
            spline_type_list[['path']]    <-  spline_type[['path']] 
          }
        } else if(is.null(spline_type[['path']])) {
          spline_type_list[['path']]    <- NULL
        }
        
        # change check message same as 'centerval' for other
        if(!is.null(spline_type[['centerval']])) {
          if(!is.numeric(spline_type[['centerval']])) {
            stop2c("Argument 'centerval' must be a numeric value",
                 " but instead specified as ", 
                 "'", paste0(spline_type[['centerval']], "'"))
          } else {
            spline_type_list[['centerval']] <- spline_type[['centerval']]
          }
        } else if(is.null(spline_type[['centerval']])) {
          spline_type_list[['centerval']] <- 0
        }
        
        
        if(!is.null(spline_type[['sfirst']])) {
          if(!is.logical(as.logical(spline_type[['sfirst']]))) {
            stop2c(paste0(spline_type[['sfirst']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['sfirst']] <- spline_type[['sfirst']]
          }
        } else if(is.null(spline_type[['sfirst']])) {
          spline_type_list[['sfirst']] <- FALSE
        }
        
        if(!is.null(spline_type[['bkrange']])) {
          if(!is.logical(as.logical(spline_type[['bkrange']]))) {
            stop2c(paste0(spline_type[['bkrange']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['bkrange']] <- spline_type[['bkrange']]
          }
        } else if(is.null(spline_type[['bkrange']])) {
          spline_type_list[['bkrange']] <- TRUE
        }
        
        if(!is.null(spline_type[['fix_bknots']])) {
          if(!is.logical(as.logical(spline_type[['fix_bknots']]))) {
            stop2c(paste0(spline_type[['fix_bknots']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['fix_bknots']] <- spline_type[['fix_bknots']]
          }
        } else if(is.null(spline_type[['fix_bknots']])) {
          spline_type_list[['fix_bknots']] <- TRUE
        }
        
        if(!is.null(spline_type[['what']])) {
          if(!is.character(spline_type[['what']])) {
            stop2c(paste0(spline_type[['what']], 
                        " must be a NULL or a character string"))
          } else {
            spline_type_list[['what']] <- spline_type[['what']]
          }
        } else if(is.null(spline_type[['what']])) {
          spline_type_list[['what']] <- NULL
        }
        
        if(!is.null(spline_type[['when']])) {
          if(!is.character(spline_type[['when']])) {
            stop2c(paste0(spline_type[['when']], 
                        " must be a NULL or a character string"))
          } else {
            spline_type_list[['when']] <- spline_type[['when']]
          }
        } else if(is.null(spline_type[['when']])) {
          spline_type_list[['when']] <- NULL
        }
        
        if(!is.null(spline_type[['return']])) {
          if(!is.logical(as.logical(spline_type[['return']]))) {
            stop2c(paste0(spline_type[['return']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['return']] <- spline_type[['return']]
          }
        } else if(is.null(spline_type[['return']])) {
          spline_type_list[['return']] <- TRUE
        }
        
        if(!is.null(spline_type[['print']])) {
          if(!is.logical(as.logical(spline_type[['print']]))) {
            stop2c(paste0(spline_type[['print']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['print']] <- spline_type[['print']]
          }
        } else if(is.null(spline_type[['print']])) {
          spline_type_list[['print']] <- TRUE
        }
        
        if(!is.null(spline_type[['sparse']])) {
          if(!is.logical(as.logical(spline_type[['sparse']]))) {
            stop2c(paste0(spline_type[['sparse']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['sparse']] <- spline_type[['sparse']]
          }
        } else if(is.null(spline_type[['sparse']])) {
          spline_type_list[['sparse']] <- FALSE
        }
        
        if(!is.null(spline_type[['check_sparsity']])) {
          if(!is.logical(as.logical(spline_type[['check_sparsity']]))) {
            stop2c(paste0(spline_type[['check_sparsity']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['check_sparsity']] <- 
              spline_type[['check_sparsity']]
          }
        } else if(is.null(spline_type[['check_sparsity']])) {
          spline_type_list[['check_sparsity']] <- FALSE
        }
        
      } else if(length(spline_type) == 0) {
        spline_type_list[['type']]        <- NULL
        spline_type_list[['intercept']]   <- FALSE
        spline_type_list[['centerval']]   <- 0
        spline_type_list[['degree']]      <- 3L
        spline_type_list[['normalize']]   <- FALSE
        spline_type_list[['derivs']]      <- FALSE
        spline_type_list[['preH']]        <- FALSE
        spline_type_list[['include']]     <- TRUE
        spline_type_list[['path']]        <- NULL
        spline_type_list[['sfirst']]      <- FALSE
        spline_type_list[['bkrange']]     <- TRUE
        spline_type_list[['fix_bknots']]  <- TRUE
        spline_type_list[['what']]        <- NULL
        spline_type_list[['when']]        <- NULL
        spline_type_list[['return']]      <- FALSE
        spline_type_list[['print']]       <- FALSE
        spline_type_list[['sparse']]      <- FALSE
        spline_type_list[['check_sparsity']]  <- FALSE
      } # if(length(spline_type) > 0) {
    } else if(!is.list(spline_type)) { 
      if(is.character(spline_type)) {
        spline_type_list[['type']]        <- spline_type
        spline_type_list[['intercept']]   <- FALSE
        spline_type_list[['centerval']]   <- 0
        spline_type_list[['degree']]      <- 3L
        spline_type_list[['normalize']]   <- TRUE
        spline_type_list[['derivs']]      <- FALSE
        spline_type_list[['preH']]        <- TRUE
        spline_type_list[['include']]     <- TRUE
        spline_type_list[['path']]        <- NULL
        spline_type_list[['sfirst']]      <- FALSE
        spline_type_list[['bkrange']]     <- TRUE
        spline_type_list[['fix_bknots']]  <- TRUE
        spline_type_list[['what']]        <- NULL
        spline_type_list[['when']]        <- NULL
        spline_type_list[['return']]      <- FALSE
        spline_type_list[['print']]       <- FALSE
        spline_type_list[['sparse']]      <- FALSE
        spline_type_list[['check_sparsity']]  <- FALSE
      } else if(!is.character(spline_type)) {
        stop2c('augument spline_type must be a 
               character string or a named list')
      } # if(is.character(spline_type)) {
    } # else if(!is.null(spline_type)) {
  } # if(is.null(spline_type)) {
  
  
  smat <- spline_type_list[['type']] 

  # This to check spline type set using the ... smat
  if(!smat %in% allowed_spline_type)
   stop2c(paste0("The spline type must be a character string.", 
               "\n  ", allowed_spline_type_exception_msg))

  
  # Except for rcs, match_sitar_a_form should be FALSE
  if(smat == 'rcs') {
    # QR
  } else if(smat == 'ns') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  } else if(smat == 'nsp') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  } else if(smat == 'nsk') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  } else if(smat == 'bsp') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  } else if(smat == 'msp') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  } else if(smat == 'isp') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  } 
  
  
  
  SplinefunxPre     <- 'GS'
  Splinefunxsuf     <- '_call'
  SplinefunxR       <- paste0(SplinefunxPre, "_", smat, Splinefunxsuf)
  SplinefunxStan    <- paste0(SplinefunxR, "_", 'stan')
  
 
  if((smat == 'nsp' | smat == 'nsk' |
      smat == 'bsp' | smat == 'msp' | 
      smat == 'isp') & 
     !spline_type_via_stype) {
    smat_intercept    <- as.integer(spline_type_list[['intercept']])
    smat_centerval    <- as.numeric(spline_type_list[['centerval']])
    smat_degree       <- as.integer(spline_type_list[['degree']])
    smat_normalize    <- as.integer(spline_type_list[['normalize']])
    smat_derivs       <- as.integer(spline_type_list[['derivs']])
    smat_preH         <- as.integer(spline_type_list[['preH']])
    smat_include_stan <- as.integer(spline_type_list[['include']])
    smat_include_path <- spline_type_list[['path']]
    smat_sfirst       <- as.integer(spline_type_list[['sfirst']])
    smat_bkrange      <- as.integer(spline_type_list[['bkrange']])
    smat_fix_bknots   <- as.integer(spline_type_list[['fix_bknots']])
    smat_what         <- (spline_type_list[['what']])
    smat_when         <- (spline_type_list[['when']])
    smat_return       <- as.integer(spline_type_list[['return']])
    smat_print        <- as.integer(spline_type_list[['print']])
    smat_sparse       <- as.integer(spline_type_list[['sparse']])
    smat_check_sparsity <- as.integer(spline_type_list[['check_sparsity']])
  } else if((smat == 'rcs') & !spline_type_via_stype) {
    smat_intercept    <- as.integer(spline_type_list[['intercept']])
    smat_centerval    <- as.numeric(spline_type_list[['centerval']])
    smat_degree       <- as.integer(spline_type_list[['degree']])
    smat_normalize    <- as.integer(spline_type_list[['normalize']])
    smat_derivs       <- as.integer(spline_type_list[['derivs']])
    smat_preH         <- as.integer(spline_type_list[['preH']])
    smat_include_stan <- as.integer(spline_type_list[['include']])
    smat_include_path <- spline_type_list[['path']]
    smat_sfirst       <- as.integer(spline_type_list[['sfirst']])
    smat_bkrange      <- as.integer(spline_type_list[['bkrange']])
    smat_fix_bknots   <- as.integer(spline_type_list[['fix_bknots']])
    smat_what         <-  (spline_type_list[['what']])
    smat_when         <-  (spline_type_list[['when']])
    smat_return       <- as.integer(spline_type_list[['return']])
    smat_print        <- as.integer(spline_type_list[['print']])
    smat_sparse       <- as.integer(spline_type_list[['sparse']]) 
    smat_check_sparsity <- as.integer(spline_type_list[['check_sparsity']])
  } else if((smat == 'nsp' | smat == 'nsk' |
             smat == 'bsp' | smat == 'msp' | 
             smat == 'isp') & 
            spline_type_via_stype) {
    smat_intercept    <- 0
    smat_centerval    <- 0
    smat_degree       <- 3L
    smat_normalize    <- as.integer(spline_type_list[['normalize']])
    smat_derivs       <- 0
    smat_preH         <- as.integer(spline_type_list[['preH']])
    smat_include_stan <- 0
    smat_include_path <- NULL
    smat_sfirst       <- as.integer(spline_type_list[['sfirst']])
    smat_bkrange      <- as.integer(spline_type_list[['bkrange']])
    smat_fix_bknots   <- as.integer(spline_type_list[['fix_bknots']])
    smat_what         <- (spline_type_list[['what']])
    smat_when         <- (spline_type_list[['when']])
    smat_return       <- as.integer(spline_type_list[['return']])
    smat_print        <- as.integer(spline_type_list[['print']])
    smat_sparse       <- as.integer(spline_type_list[['sparse']]) 
    smat_check_sparsity <- as.integer(spline_type_list[['check_sparsity']])
  } else if((smat == 'rcs') & spline_type_via_stype) {
    smat_intercept    <- 0
    smat_centerval    <- 0
    smat_degree       <- 3L
    smat_normalize    <- as.integer(spline_type_list[['normalize']])
    smat_derivs       <- 0
    smat_preH         <- as.integer(spline_type_list[['preH']])
    smat_include_stan <- 0
    smat_include_path <- NULL
    smat_sfirst       <- as.integer(spline_type_list[['sfirst']])
    smat_bkrange      <- as.integer(spline_type_list[['bkrange']])
    smat_fix_bknots   <- as.integer(spline_type_list[['fix_bknots']])
    smat_what         <- (spline_type_list[['what']])
    smat_when         <- (spline_type_list[['when']])
    smat_return       <- as.integer(spline_type_list[['return']])
    smat_print        <- as.integer(spline_type_list[['print']])
    smat_sparse       <- as.integer(spline_type_list[['sparse']]) 
    smat_check_sparsity <- as.integer(spline_type_list[['check_sparsity']])
  } else {
    # allow further checks - for later use
  }
  
  
  if(smat_sparse) {
    if(!smat_sfirst) {
      stop2c("If 'smat_sparse = TRUE', then 
             'smat_sfirst' must also be set as 'TRUE'")
    }
  }
  
  

  if(smat_check_sparsity) {
    if(!smat_sfirst | !smat_sparse) {
      stop2c("If 'check_sparsity = TRUE', then both
      'smat_sfirst' and 'smat_sparse'
           must also be set as 'TRUE'")
    }
    if(arguments$chains > 1) stop2c("'chains' must be 
                                    set as '1' when check_sparsity = TRUE'")
    if(arguments$iter > 2) stop2c("'iter' must be set as
                                  '2' when check_sparsity = TRUE'")
  }
  
  
 # if(smat == 'rcs') {
 #   smat_bkrange <- smat_bkrange
 # } else {
 #   # if(as.logical(smat_bkrange)) {
 #   #   stop2c("Argument 'bkrange' can only be used when smat is 'rcs', and not ", 
 #   #        collapse_comma(smat),
 #   #        "\n  ", 
 #   #        "Please check the 'stype' argument and correct it")
 #   # }
 #   # smat_bkrange <- NULL
 # }
  
  
  
  # over ride match_sitar_a_form to match smat_bkrange
  if(smat == 'rcs') {
    if(!as.logical(smat_bkrange) | !as.logical(smat_fix_bknots)) {
      temp_match_sitar_a_form <- getdotslist[['match_sitar_a_form']]
      if(is.null(temp_match_sitar_a_form)) {
        getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
        if(verbose) {
          message2c("The 'match_sitar_a_form' has been set as 'FALSE' because 
                  'bkrange' was set FALSE")
        }
      }
    }
  } # if(smat == 'rcs') {
  
  
  
  if(smat == 'nsp' | smat == 'nsk' | smat == 'rcs') {
    if(smat_degree != 3) {
      stop2c("'nsp', 'nsk' and 'rcs', the degree must be 3")
    }
  }
     

  if(smat_include_stan == 1) {
    # stop2c("Please set smat_include_stan = 0")
    smat_include_stan <- 0
   if(verbose) message2c("'smat_include_stan' is set to '0'")
  }
 
 
   if(verbose) {
     message2c(paste0("setting spline type as '", smat, "'"))
     if(smat != "rcs") {
       message2c(paste0("setting degree for spline type '",
                      spline_type_list[['degree']], "' as: ", smat_degree))
       message2c(paste0("setting intercept for spline type '",
                      spline_type_list[['type']], "' as: ", smat_intercept))
       message2c(paste0("setting normalize for spline type '",
                      spline_type_list[['type']], "' as: ", smat_normalize))
       message2c(paste0("setting centerval for spline type '",
                      spline_type_list[['type']], "' as: ", smat_centerval))
     }
   }
   

  
  
  # Handle fast_nsk
  if(is.null(getdotslist[['fast_nsk']])) {
    if(arguments$backend == "cmdstanr") {
      fast_nsk <- 2L
    } else if(arguments$backend == "rstan" | 
              eval(getOption("brms.backend", "rstan")) == "rstan") {
      if(utils::packageVersion('rstan') > '2.35.0') {
        fast_nsk <- 2L
      } else {
        fast_nsk <- 0L
      }
    } else {
      fast_nsk <- 0L
    }
  } else if(!is.null(getdotslist[['fast_nsk']])) {
    fast_nsk <- as.integer(getdotslist[['fast_nsk']])
    getdotslist[['fast_nsk']] <- NULL
  }
  
  
  # if smat != 'nsk' | 'nsp' -> set fast_nsk <- 0L
  if(smat == 'nsp') {
    fast_nsk <- fast_nsk
  } else if(smat == 'nsk') {
    fast_nsk <- fast_nsk
  } else if(smat == 'bsp') {
    fast_nsk <- fast_nsk
  } else if(smat == 'isp') {
    fast_nsk <- fast_nsk
  } else if(smat == 'msp') {
    fast_nsk <- fast_nsk
  } else if(smat == 'rcs') {
    fast_nsk <- 0L
  } else {
    fast_nsk <- 0L
  }
  

  if(smat == 'isp') {
    smat_moi <- TRUE
  } else {
    smat_moi <- FALSE
  }
  
  
 

  # 24.08.2024
  if(is.null(getdotslist[['match_sitar_a_form']])) {
    if(is.null(decomp)) {
      match_sitar_a_form <- TRUE
    } else if(!is.null(decomp)) {
      if(decomp == 'QR') {
        getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
      } # if(decomp == 'QR') {
    } # if(!is.null(decomp)) {
    # match_sitar_a_form <- TRUE
  } else {
    match_sitar_a_form <- getdotslist[['match_sitar_a_form']]
  }
  
  if(verbose) {
    message2c("The 'a' form for sitar model has been set as ", 
            "'", match_sitar_a_form, "'",
            " (see '...$match_sitar_a_form')")
  }
  
  
  
  # 24.08.2024
  if(is.null(getdotslist[['sigmamatch_sitar_a_form']])) {
    sigmamatch_sitar_a_form <- FALSE
  } else {
    sigmamatch_sitar_a_form <- getdotslist[['sigmamatch_sitar_a_form']]
  }
  
  
  
  # 24.08.2024
  if(is.null(getdotslist[['sigmad_adjusted']])) {
    sigmad_adjusted <- FALSE
  } else {
    sigmad_adjusted <- getdotslist[['sigmad_adjusted']]
  }
  
  
  # New
  # Now fixed/random form will decide d form, so setting universally FALSE
  getdotslist[['match_sitar_d_form']] <- NULL
  if(!is.null(getdotslist[['match_sitar_d_form']])) {
    match_sitar_d_form <- getdotslist[['match_sitar_d_form']] 
  }
  
  if(!is.null(getdotslist[['match_sitar_d_form']])) {
    if(getdotslist[['match_sitar_d_form']]) {
      # arguments$select_model <- select_model <- 'sitar4r'
    }
  }

  # 24.08.2024
  if(is.null(getdotslist[['match_sitar_d_form']])) { 
    match_sitar_d_form <- FALSE
    if(grepl('sitar4', select_model)) {
      if(select_model == 'sitar4fr') match_sitar_d_form <- FALSE
      if(select_model == 'sitar4f')  match_sitar_d_form <- FALSE
      if(select_model == 'sitar4r')  match_sitar_d_form <- TRUE
      if(select_model == 'sitar4')   match_sitar_d_form <- FALSE
      sitar_nparms <- 4
      select_model <- 'sitar'
    } else if(select_model == 'sitar') {
      sitar_nparms <- 3
      select_model <- 'sitar'
      match_sitar_d_form <- FALSE
    } else {
      match_sitar_d_form <- FALSE
    }
  } else {
    if(select_model == 'sitar4r')  {
      match_sitar_d_form <- getdotslist[['match_sitar_d_form']]
    } else {
      if(getdotslist[['match_sitar_d_form']]) {
        # stop2c("match_sitar_d_form = TRUE only allowed for 
        #        sitar model 'sitar4r'")
      }
    }
  }
  
  
  
 
  sitar_models    <- c('sitar', 'sitar3', 'sitar4', 
                       'sitar4f', 'sitar4fr', 'sitar4r')
  pb_models       <- c('pb1', 'pb2', 'pb3')
  logistic_models <- c('logistic1', 'logistic2', 'logistic3')
  rcs_models      <- c('rcs', 'rcsf', 'rcsfr')
  
  allowed_model_names <- c(sitar_models, 
                           pb_models, 
                           logistic_models, 
                           rcs_models)
  
  sitar_models    <- paste0(sitar_models, collapse=", ")
  pb_models       <- paste0(pb_models, collapse=", ")
  logistic_models <- paste0(logistic_models, collapse=", ")
  rcs_models      <- paste0(rcs_models, collapse=", ")
  
  all_models <- paste0('NLME models: ', "\n  ", 
                       '  SITAR: ',  sitar_models, "\n  ", 
                       '  PB: ', pb_models, "\n  ", 
                       '  LOGISTIC: ', logistic_models, "\n  ",
                       'LME models: ', "\n  ", 
                       '  RCS: ', rcs_models)
  
 
  if(!select_model_arg %in% allowed_model_names) {
    stop2c("Currently supported models (via 'select_model' argument) are:",
         "\n ",
         " ", all_models
         )
  }
  
  
  for (ip in names(arguments)) {
    if (grepl("_init_", ip)) {
      d_mcall_ <- deparse_0(mcall_[[ip]])
      if(is.symbol(mcall_[[ip]])) {
        arguments[[ip]] <- d_mcall_
      } else if(grepl("c(", d_mcall_, fixed = T)) {
        arguments[[ip]] <- gsub("c(", "list(", d_mcall_, fixed = T)
      } else if(grepl("list(", d_mcall_, fixed = T)) {
        arguments[[ip]] <- mcall_[[ip]]
      }
    }
  }
  
  remove_spaces <- c('a_formula_gr_str', 'b_formula_gr_str', 
                     'c_formula_gr_str', 'd_formula_gr_str',
                     'e_formula_gr_str', 'f_formula_gr_str', 
                     'g_formula_gr_str', 'h_formula_gr_str', 
                     'i_formula_gr_str', 's_formula_gr_str', 
                     'sigma_formula_gr_str')
  
  for (ip in remove_spaces) {
    arguments[[ip]] <-  gsub_space(arguments[[ip]] )
  }
  
  brms_arguments_list <-
    c(
      'chains',
      'iter',
      'warmup',
      'thin',
      'cores',
      'backend',
      'threads',
      'opencl',
      'normalize',
      'algorithm',
      'control',
      'empty',
      'rename',
      'sample_prior',
      'save_pars',
      'drop_unused_levels',
      'stan_model_args',
      'refresh',
      'silent',
      'seed',
      'save_model',
      'fit',
      'file',
      'file_refit',
      'file_compress',
      'future',
      "data2"
    )
  
  

  if(is.numeric(arguments$cores)) {
   oldopts <- options(mc.cores = arguments$cores)
   on.exit(options(oldopts))
  }
  
  iter <-  arguments$iter
  warmup <-  arguments$warmup <- eval(arguments$warmup)
 
 
  brms_arguments <- list()
  for (brms_arguments_listi in brms_arguments_list) {
    brms_arguments[[brms_arguments_listi]] <-
      arguments[[brms_arguments_listi]]
    arguments[[brms_arguments_listi]] <- NULL
  }
  
  
  brms_arguments <- mget(brms_arguments_list)
  
  # Set path for s files
  if(smat == 'nsp' | smat == 'nsk') {
    if(smat_include_stan) {
      if(is.null(brms_arguments$stan_model_args)) {
        brms_arguments$stan_model_args <- list()
        brms_arguments$stan_model_args[['include_paths']] <- "."
        if(verbose) 
          message2c("path for .stan file(s) set to '.' via 'stan_model_args'")
      } else if(!is.null(brms_arguments$stan_model_args)) {
        if(is.list(brms_arguments$stan_model_args)) {
          if(is.null(brms_arguments$stan_model_args[['include_paths']])) {
            brms_arguments$stan_model_args[['include_paths']] <- "."
            if(verbose) 
              message2c("path for .stan file(s) set to '.' via 'stan_model_args'")
          }
        }
      } # if(is.null(brms_arguments$stan_model_args)) {
    } # if(smat_include_stan) {
  } # if(smat == 'nsp' | smat == 'nsk') {

  
  if(smat == 'nsp' | smat == 'nsk') {
    if(smat_include_stan) {
      if(is.null( brms_arguments$stan_model_args[['include_paths']])) {
        stop2c("Please specify path for .stan file(s) via 'stan_model_args'")
      }
    }
  }
  
  
  if (eval(brms_arguments$backend) != "rstan" &
      eval(brms_arguments$backend) != "mock" &
      eval(brms_arguments$backend) != "cmdstanr") {
    stop2c("The backend argument must be either 'rstan', 'mock', or 'cmdstanr'",
         "\n ",
         "\ Please check it which you have specified as: ", 
         eval(brms_arguments$backend))
  }
  
  
  if(is.null(getdotslist[['displayit']])) { 
    displayit <- 'col'
  } else {
    displayit <- getdotslist[['displayit']]
  }
  
  
  if(is.null(getdotslist[['setcolh']])) { 
    setcolh <- 47
  } else {
    setcolh <- getdotslist[['setcolh']]
  }
  
  
  if(is.null(getdotslist[['setcolb']])) { 
    setcolb <- 3
  } else {
    setcolb <- getdotslist[['setcolb']]
  }
  
 
  # Quote unquoted character (e.g., sex to 'sex')
  list_to_quoted_if_not <- function(x) {
    splitmvar <- x
    splitmvar <- gsub("\\s", "", splitmvar)
    splitmvar <- paste(splitmvar, collapse = "")
    splitmvar_w <-
      gsub("[\\(\\)]", "", regmatches(splitmvar, gregexpr("\\(.*?\\)",
                                                          splitmvar))[[1]])
    splitmvar_w2 <- strsplit(splitmvar_w, ",")[[1]]
    splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
    gsubs_c_counter <- 0
    for (i in splitmvar_w3) {
      gsubs_c_counter <- gsubs_c_counter + 1
      if (gsubs_c_counter < max(length(splitmvar_w3))) {
        pattern <- paste0(i, "=", "\\s*(.*?)\\s*", ",")
      } else {
        pattern <- paste0(i, "=", "\\s*(.*?)\\s*", ")")
      }
      majors <- regmatches(splitmvar, regexec(pattern, splitmvar))
      majors2 <- majors[[1]][2]
      majors2 <- majors[[1]][2]
      if (grepl("^T$", majors2)) {
        majors2 <- gsub("^T$", "TRUE", majors2)
      }
      if (grepl("^F$", majors2)) {
        majors2 <- gsub("^F$", "FALSE", majors2)
      }
      majors2 <- gsub("\"", "", majors2)
      majors3 <- paste0("\"", majors2, "\"")
      if (gsubs_c_counter == 1) {
        splitmvar2 <- gsub(noquote(majors2), majors3, 
                           splitmvar, fixed = FALSE)
      } else {
        # splitmvar2 <- gsub(noquote(majors2), majors3, splitmvar2, fixed = F)
        splitmvar2 <- gsub(paste0('\\<', noquote(majors2), '\\>'), majors3, 
                           splitmvar2, fixed = FALSE)
      }
    }
    
    for (i in 1:length(splitmvar_w3))
      splitmvar2 <- gsub("\"\"", "\"", splitmvar2)
    splitmvar3 <- eval(parse(text = splitmvar2))
    zzz <- splitmvar3
    
    for (z in names(splitmvar3)) {
      assign('err.', FALSE, envir = enverr.)
      tryCatch(
        expr = {
          eval(parse(text = zzz[[z]]), envir = parent.frame())
        },
        error = function(e) {
          assign('err.', TRUE, envir = enverr.)
        }
      )
      err. <- get('err.', envir = enverr.)
      if (!err.) {
        # if brms::brmsfamily(family), eval eliminates family 16 1. 2024
        if(z != "family") c_c_ <- eval(parse(text = zzz[[z]]))
        if(z == "family") c_c_ <- zzz[[z]] 
        checkclass <- class(c_c_)
        if (checkclass == "NULL")
          checkclass_ <- NULL
        else
          checkclass_ <- NA
        if (is.logical(c_c_) | is.null(checkclass_))
          zzz[[z]] <- c_c_
      } else {
        zzz[[z]] <- zzz[[z]]
      }
    }
    return(zzz)
  }
  
  
  list_to_quoted_if_not_si <- function(xx) {
    xx.o <- xx
    prefix_ <- strsplit(xx, "\\(")[[1]][1]
    prefix_by <- "list"
    xx <- gsub(paste0("^", prefix_, ""), prefix_by, xx)
    if (sub("\\).*", "", sub(".*\\(", "", xx)) != "") {
      xxx <- list_to_quoted_if_not(xx)
      xxx <- gsub("\"" , "'", deparse_0(xxx))
      xxx <- gsub(paste0("^", prefix_by, ""), prefix_, xxx)
      xxx <- gsub("\\s", "", xxx)
    } else {
      xxx <- xx.o
    }
    xxx
  }
  
  
  list_to_quoted_if_not_si_lf <- function(xx) {
    xx.o <- xx
    prefix_ <- strsplit(xx, "\\(")[[1]][1]
    prefix_by <- "list"
    xx <- gsub(paste0("^", prefix_, ""), prefix_by, xx)
    if (sub("\\).*", "", sub(".*\\(", "", xx)) != "") {
      xxt <- sub("\\).*", "", sub(".*\\(", "", xx))
      xxtf <-
        strsplit(xxt, ",")[[1]][grepl("~", strsplit(xxt, ",")[[1]])]
      xxtnf <-
        strsplit(xxt, ",")[[1]][!grepl("~", strsplit(xxt, ",")[[1]])]
      xxtf <- gsub("\\s", "", xxtf)
      xxtnf <- gsub("\\s", "", xxtnf)
      xx <-
        paste0(prefix_by, "(", paste(xxtnf, collapse = ","), ")")
      xxx <- list_to_quoted_if_not(xx)
      xxx <- gsub("\"" , "'", deparse_0(xxx))
      xxx <- gsub(paste0("^", prefix_by, ""), prefix_, xxx)
      xxx <-
        gsub(paste0(prefix_, "\\("),
             paste0(prefix_, "(", xxtf, ","),
             xxx)
      xxx <- gsub("\\s", "", xxx)
    } else {
      xxx <- xx.o
    }
    xxx
  }
  
  
 
  # set multivariate arguments
  if (gsub("\\s", "",
           paste(deparse(substitute(multivariate)), collapse = "")) == "NULL" |
      gsub("\\s", "",
           paste(deparse(substitute(multivariate)), collapse = "")) == "NA" |
      gsub("\\s", "",
           paste(deparse(substitute(multivariate)), collapse = "")) == "FALSE" |
      gsub("\\s", "",
           paste(deparse(substitute(multivariate)), collapse = "")) == "F") {
    multivariate <- list()
    multivariate$mvar <- FALSE
  } else if (gsub("\\s", "",
                  paste(deparse(substitute(multivariate)),
                        collapse = "")) == "TRUE" |
             gsub("\\s", "",
                  paste(deparse(substitute(multivariate)),
                        collapse = "")) == "T") {
    multivariate <- list()
    multivariate$mvar <- TRUE
    
  } else if (!grepl("^list", gsub("\\s", "", paste(deparse(
    substitute(multivariate)
  ), collapse = ""))) &
  !is.null(gsub("\\s", "", paste(deparse(
    substitute(multivariate)
  ), collapse = "")))) {
    if (is.symbol(substitute(multivariate))) {
      if(!grepl("^list", deparse_0(multivariate))) { # For CustomDoCall
        multivariate <- gsub("\\s", "", paste(deparse(substitute(multivariate)), 
                                              collapse = ""))
        if (multivariate == "T") multivariate <- eval(parse(text = multivariate))
        multivariate <- as.list(multivariate)
        names(multivariate) <- 'mvar'
      } # For CustomDoCall
      
    } else if (is.character(substitute(multivariate))) {
      multivariate <- multivariate
      multivariate <- as.list(multivariate)
      names(multivariate) <- 'mvar'
    }
  }
  if (grepl("^list", gsub("\\s", "", paste(deparse(
    substitute(multivariate)
  ), collapse = ""))) &
  length(strsplit(gsub("\\s", "", paste(
    deparse(substitute(multivariate)), collapse = ""
  )), ",")[[1]]) == 1) {
    if (!is.null(gsub("\\s", "", paste(deparse(
      substitute(multivariate)
    ), collapse = "")))) {
      if (is.language(substitute(multivariate))) {
        ttt <-
          gsub("\\s", "", paste(deparse(substitute(
            multivariate
          )), collapse = ""))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop2c("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("mvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        multivariate <- list_to_quoted_if_not(ttt)
        
      } else if (grepl("^list", multivariate)) {
        ttt <- deparse_0(as.name(substitute(multivariate)))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop2c("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("mvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        multivariate <- list_to_quoted_if_not(ttt)
        for (multivariatei in 1:length(multivariate)) {
          if (!is.null(multivariate[[multivariatei]])) {
            multivariate[[multivariatei]] <-
              gsub("'", "", multivariate[[multivariatei]])
          }
        }
      } else {
        if (!is.null(multivariate)) {
          if (!grepl("^list", multivariate)) {
            multivariate <- multivariate
            multivariate <- as.list(multivariate)
            names(multivariate) <- 'mvar'
          }
        } else if (is.null(multivariate)) {
          multivariate <- as.list(multivariate)
        }
      }
    }
  }
  if (length(strsplit(gsub("\\s", "", paste(
    deparse(substitute(multivariate)), collapse = ""
  )), ",")[[1]]) > 1) {
    ttt <-
      gsub("\\s", "", paste(deparse(substitute(multivariate)), collapse = ""))
    multivariate <- list_to_quoted_if_not(ttt)
  }
  
  
  

  # Set univariate_by arguments
  if (gsub("\\s", "",
           paste(deparse(substitute(univariate_by)), 
                 collapse = "")) == "NULL" |
      gsub("\\s", "",
           paste(deparse(substitute(univariate_by)), 
                 collapse = "")) == "NA" |
      gsub("\\s", "",
           paste(deparse(substitute(univariate_by)), 
                 collapse = "")) == "FALSE" |
      gsub("\\s", "",
           paste(deparse(substitute(univariate_by)), 
                 collapse = "")) == "F") {
    univariate_by <- list()
    univariate_by$by <- NA
  } else if (!grepl("^list", gsub("\\s", "", paste(deparse(
    substitute(univariate_by)
  ), collapse = ""))) &
  !is.null(gsub("\\s", "", paste(deparse(
    substitute(univariate_by)
  ), collapse = "")))) {
    if (is.symbol(substitute(univariate_by))) {
      if(!grepl("^list", deparse_0(multivariate))) { # For CustomDoCall
        univariate_by <- gsub("\\s", "", paste(deparse(substitute(univariate_by)),
                                               collapse = ""))
        univariate_by <- univariate_by
        univariate_by <- as.list(univariate_by)
        names(univariate_by) <- 'by'
      }
      
    } else if (is.character(substitute(univariate_by))) {
      univariate_by <- univariate_by
      univariate_by <- as.list(univariate_by)
      names(univariate_by) <- 'by'
    }
  }
  if (grepl("^list", gsub("\\s", "", paste(deparse(
    substitute(univariate_by)
  ), collapse = ""))) &
  length(strsplit(gsub("\\s", "", paste(
    deparse(substitute(univariate_by)), collapse = ""
  )), ",")[[1]]) == 1) {
    if (!is.null(gsub("\\s", "", paste(deparse(
      substitute(univariate_by)
    ), collapse = "")))) {
      if (is.language(substitute(univariate_by))) {
        ttt <-
          gsub("\\s", "", paste(deparse(substitute(
            univariate_by
          )), collapse = ""))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop2c("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("by=", strsplit(temp, "=")[[1]]), ttt)
        }
        univariate_by <- list_to_quoted_if_not(ttt)
        
      } else if (grepl("^list", univariate_by)) {
        ttt <- deparse_0(as.name(substitute(univariate_by)))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop2c("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("by=", strsplit(temp, "=")[[1]]), ttt)
        }
        univariate_by <- list_to_quoted_if_not(ttt)
        for (univariate_byi in 1:length(univariate_by)) {
          if (!is.null(univariate_by[[univariate_byi]])) {
            univariate_by[[univariate_byi]] <-
              gsub("'", "", univariate_by[[univariate_byi]])
          }
        }
      } else {
        if (!is.null(univariate_by)) {
          if (!grepl("^list", univariate_by)) {
            univariate_by <- univariate_by
            univariate_by <- as.list(univariate_by)
            names(univariate_by) <- 'by'
          }
        } else if (is.null(univariate_by)) {
          univariate_by <- as.list(univariate_by)
        }
      }
    }
  }
  if (length(strsplit(gsub("\\s", "", paste(
    deparse(substitute(univariate_by)), collapse = ""
  )), ",")[[1]]) > 1) {
    ttt <-
      gsub("\\s", "", paste(deparse(substitute(univariate_by)), 
                            collapse = ""))
    univariate_by <- list_to_quoted_if_not(ttt)
  }
  
  
  
  
  # Set group_arg arguments 
  if (!paste(deparse(substitute(group_arg)), collapse = "") == "NULL"  &
      !any(grepl("^list", gsub("\\s", "", paste(
        deparse(substitute(group_arg)), collapse = ""
      )))) &
      any(gsub("\\s", "", paste(deparse(
        substitute(group_arg)
      ), collapse = "")) == "NULL")) {
    group_arg <- list()
    group_arg$groupvar <- NULL
  } else if (!any(grepl("^list", gsub("\\s", "", paste(
    deparse(substitute(group_arg)), collapse = ""
  )))) &
  any(gsub("\\s", "", paste(deparse(
    substitute(group_arg)
  ), collapse = "")) != "NULL")) {
    if (paste(deparse(substitute(group_arg)), collapse = "") == "T" |
        paste(deparse(substitute(group_arg)), collapse = "") == "TRUE" |
        paste(deparse(substitute(group_arg)), collapse = "") == "F" |
        paste(deparse(substitute(group_arg)), collapse = "") == "FALSE" |
        paste(deparse(substitute(group_arg)), collapse = "") == "NA") {
      stop2c("group_arg should be either NULL or a character",
           " denoting the group idetifier")
    }
    if (is.symbol(substitute(group_arg))) {
      group_arg <-
        gsub("\\s", "", paste(deparse(substitute(group_arg)), collapse = ""))
      group_arg <- group_arg
      group_arg <- as.list(group_arg)
      names(group_arg) <- 'groupvar'
    } else if (is.character(substitute(group_arg))) {
      group_arg <- group_arg
      group_arg <- as.list(group_arg)
      names(group_arg) <- 'groupvar'
    }
  }
  if (any(grepl("^list", gsub("\\s", "",
                              paste(
                                deparse(substitute(group_arg)),
                                collapse = ""
                              )))) &
      length(strsplit(gsub("\\s", "",
                           paste(
                             deparse(substitute(group_arg)),
                             collapse = ""
                           )), ",")[[1]]) == 1) {
    if (!is.null(gsub("\\s", "", paste(deparse(
      substitute(group_arg)
    ),
    collapse = "")))) {
      if (is.language(substitute(group_arg))) {
        ttt <- gsub("\\s", "", paste(deparse(substitute(group_arg)),
                                     collapse = ""))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "T" |
            temp == "TRUE" |
            temp == "F" |
            temp == "FALSE") {
          stop2c(
            "group_arg should be either NULL or a character",
            " denoting the group idetifier"
          )
        }
        if (temp == "") {
          stop2c("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("groupvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        group_arg <- list_to_quoted_if_not(ttt)
      } else if (grepl("^list", group_arg)) {
        ttt <- deparse_0(as.name(substitute(group_arg)))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop2c("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("groupvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        group_arg <- list_to_quoted_if_not(ttt)
        for (group_argi in 1:length(group_arg)) {
          if (!is.null(group_arg[[group_argi]])) {
            group_arg[[group_argi]] <- gsub("'", "", group_arg[[group_argi]])
          }
        }
      } else {
        if (!is.null(group_arg)) {
          if (!grepl("^list", gsub("\\s", "",
                                   paste(
                                     deparse(substitute(group_arg)),
                                     collapse = ""
                                   )))) {
            group_arg <- group_arg
            group_arg <- as.list(group_arg)
            names(group_arg) <- 'groupvar'
          }
        } else if (is.null(group_arg)) {
          group_arg <- as.list(group_arg)
        }
      }
    } else if (is.null(group_arg)) {
      group_arg <- list()
      group_arg$groupvar <- NULL
    }
  }
  if (any(grepl("^list", gsub("\\s", "",
                              paste(
                                deparse(substitute(group_arg)),
                                collapse = ""
                              )))) &
      length(strsplit(gsub("\\s", "",
                           paste(
                             deparse(substitute(group_arg)),
                             collapse = ""
                           )), ",")[[1]]) > 1) {
    ttt <-
      gsub("\\s", "", paste(deparse(substitute(group_arg)), collapse = ""))
    group_arg <- list_to_quoted_if_not(ttt)
  }
  if (length(group_arg) == 0) {
    group_arg <- list()
    group_arg$groupvar <- NULL
  }
  if (!is.null(group_arg$groupvar) &
      !is.character(group_arg$groupvar)) {
    stop2c("group_arg should be either NULL or a character",
         " denoting the group idetifier")
  }
  
 
  # Set up sigma_group_arg arguments 
  if (!paste(deparse(substitute(sigma_group_arg)), collapse = "") == "NULL"  &
      !any(grepl("^list", gsub("\\s", "", paste(
        deparse(substitute(sigma_group_arg)), collapse = ""
      )))) &
      any(gsub("\\s", "", paste(deparse(
        substitute(sigma_group_arg)
      ), collapse = "")) == "NULL")) {
    sigma_group_arg <- list()
    sigma_group_arg$groupvar <- NULL
  } else if (!any(grepl("^list", gsub("\\s", "", paste(
    deparse(substitute(sigma_group_arg)), collapse = ""
  )))) &
  any(gsub("\\s", "", paste(deparse(
    substitute(sigma_group_arg)
  ), collapse = "")) != "NULL")) {
    if (paste(deparse(substitute(sigma_group_arg)), collapse = "") == "T" |
        paste(deparse(substitute(sigma_group_arg)), collapse = "") == "TRUE" |
        paste(deparse(substitute(sigma_group_arg)), collapse = "") == "F" |
        paste(deparse(substitute(sigma_group_arg)), collapse = "") == "FALSE" |
        paste(deparse(substitute(sigma_group_arg)), collapse = "") == "NA") {
      stop2c("sigma_group_arg should be either NULL or a character",
           " denoting the group idetifier")
    }
    if (is.symbol(substitute(sigma_group_arg))) {
      sigma_group_arg <-
        gsub("\\s", "", 
             paste(deparse(substitute(sigma_group_arg)), collapse = ""))
      sigma_group_arg <- sigma_group_arg
      sigma_group_arg <- as.list(sigma_group_arg)
      names(sigma_group_arg) <- 'groupvar'
    } else if (is.character(substitute(sigma_group_arg))) {
      sigma_group_arg <- sigma_group_arg
      sigma_group_arg <- as.list(sigma_group_arg)
      names(sigma_group_arg) <- 'groupvar'
    }
  }
  if (any(grepl("^list", gsub("\\s", "",
                              paste(
                                deparse(substitute(sigma_group_arg)),
                                collapse = ""
                              )))) &
      length(strsplit(gsub("\\s", "",
                           paste(
                             deparse(substitute(sigma_group_arg)),
                             collapse = ""
                           )), ",")[[1]]) == 1) {
    if (!is.null(gsub("\\s", "", paste(deparse(
      substitute(sigma_group_arg)
    ),
    collapse = "")))) {
      if (is.language(substitute(sigma_group_arg))) {
        ttt <- gsub("\\s", "", paste(deparse(substitute(sigma_group_arg)),
                                     collapse = ""))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "T" |
            temp == "TRUE" |
            temp == "F" |
            temp == "FALSE") {
          stop2c(
            "sigma_group_arg should be either NULL or a character",
            " denoting the group idetifier"
          )
        }
        if (temp == "") {
          stop2c("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("groupvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        sigma_group_arg <- list_to_quoted_if_not(ttt)
      } else if (grepl("^list", sigma_group_arg)) {
        ttt <- deparse_0(as.name(substitute(sigma_group_arg)))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop2c("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("groupvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        sigma_group_arg <- list_to_quoted_if_not(ttt)
        for (sigma_group_argi in 1:length(sigma_group_arg)) {
          if (!is.null(sigma_group_arg[[sigma_group_argi]])) {
            sigma_group_arg[[sigma_group_argi]] <- 
              gsub("'", "", sigma_group_arg[[sigma_group_argi]])
          }
        }
      } else {
        if (!is.null(sigma_group_arg)) {
          if (!grepl("^list", gsub("\\s", "",
                                   paste(
                                     deparse(substitute(sigma_group_arg)),
                                     collapse = ""
                                   )))) {
            sigma_group_arg <- sigma_group_arg
            sigma_group_arg <- as.list(sigma_group_arg)
            names(sigma_group_arg) <- 'groupvar'
          }
        } else if (is.null(sigma_group_arg)) {
          sigma_group_arg <- as.list(sigma_group_arg)
        }
      }
    } else if (is.null(sigma_group_arg)) {
      sigma_group_arg <- list()
      sigma_group_arg$groupvar <- NULL
    }
  }
  if (any(grepl("^list", gsub("\\s", "",
                              paste(
                                deparse(substitute(sigma_group_arg)),
                                collapse = ""
                              )))) &
      length(strsplit(gsub("\\s", "",
                           paste(
                             deparse(substitute(sigma_group_arg)),
                             collapse = ""
                           )), ",")[[1]]) > 1) {
    ttt <- gsub("\\s", "", paste(deparse(substitute(sigma_group_arg)), 
                                 collapse = ""))
    sigma_group_arg <- list_to_quoted_if_not(ttt)
  }
  if (length(sigma_group_arg) == 0) {
    sigma_group_arg <- list()
    sigma_group_arg$groupvar <- NULL
  }
  if (!is.null(sigma_group_arg$groupvar) &
      !is.character(sigma_group_arg$groupvar)) {
    stop2c("sigma_group_arg should be either NULL or a character",
         " denoting the group idetifier")
  }
  
  
  
  # Add defaults to univariate_by, multivariate, and group_arg arguments
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    univariate_by$by <- gsub("\\s", "", univariate_by$by)
  }
  if (identical(univariate_by$by, character(0))) {
    univariate_by$by <- NA
  }
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (univariate_by$by == "" |
        univariate_by$by == FALSE | is.null(univariate_by$by)) {
      univariate_by$by <- NA
    }
    if (univariate_by$by == TRUE) {
      stop2c(
        "For univeriate-by-subgroup model fitting (via univariate_by argument)",
        "\n ",
        "argument 'by' should be a variable name, '', NULL, or FALSE"
      )
    }
  }
  
  
  
  if (multivariate$mvar &
      !(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    stop2c(
      "You have set multivariate as TRUE and also specified ",
      "\n ",
      " univeriate-by-subgroup model (see univariate_by argument)",
      "\n ",
      " Please specify either multivariate or univariate_by argument"
    )
  }
  
  
  if (is.symbol(arguments[["y"]]) |
      is.character(arguments[["y"]])) {
    nys <- length(arguments[["y"]])
  } else {
    nys <- length(arguments[["y"]]) - 1
  }
  
  
  if (multivariate$mvar & nys == 1) {
    stop2c(
      "You have set multivariate as TRUE but provided only one outcome ",
      "\n ",
      " Please set y as list or vector of multiple outcomes such as ",
      "\n ",
      " list(outcome1, outcome2) or y = c(outcome1, outcome2)"
    )
  }
  
  if (!multivariate$mvar & nys > 1) {
    stop2c(
      "You have set multivariate as FALSE but provided more than one outcome",
      "\n ",
      " Please set y as a symbol / list / vector of single outcome such as",
      "\n ",
      " y = outcome, y = list(outcome1) or y = c(outcome1)"
    )
  }
  
  if (!(is.na(univariate_by$by) |
        univariate_by$by == "NA") & nys > 1) {
    stop2c(
      "You have specified univariate_by model for ",
      univariate_by$by,
      "for which ",
      "\n ",
      " only one outcome varibale should be specified but have provided ",
      nys,
      " outcomes",
      "\n ",
      " Please set y as a symbol / list / vector of single outcome ",
      "\n ",
      " such as y = outcome, y = list(outcome1) or y = c(outcome1)"
    )
  }
  
  if (multivariate$mvar) {
    if (is.null(multivariate$cor))
      multivariate$cor <- "un"
    if (is.null(multivariate$rescor))
      multivariate$rescor <- TRUE
    if (is.null(multivariate$rcorr_by))
      multivariate$rcorr_by <- NULL
    if (is.null(multivariate$rcorr_gr))
      multivariate$rcorr_gr <- NULL
    if (is.null(multivariate$rcorr_prior))
      multivariate$rcorr_prior <- NULL
    if (is.null(multivariate$rcorr_method))
      multivariate$rcorr_method <- NULL
  }
  if (!multivariate$mvar) {
    if (is.null(multivariate$cor))
      multivariate$cor <- "un"
    if (is.null(multivariate$rescor))
      multivariate$rescor <- FALSE # TRUE -> FALSE
    if (is.null(multivariate$rcorr_by))
      multivariate$rcorr_by <- NULL
    if (is.null(multivariate$rcorr_gr))
      multivariate$rcorr_gr <- NULL
    if (is.null(multivariate$rcorr_prior))
      multivariate$rcorr_prior <- NULL
    if (is.null(multivariate$rcorr_method))
      multivariate$rcorr_method <- NULL
  }
  
  
  
  if (is.na(univariate_by$by)) {
    if (is.null(univariate_by$cor))
      univariate_by$cor <- "un"
  }
  if (!is.na(univariate_by$by)) {
    if (is.null(univariate_by$cor))
      univariate_by$cor <- "un"
  }
  
  if (is.na(univariate_by$by)) {
    if (is.null(univariate_by$terms))
      univariate_by$terms <- "subset"
  }
  if (!is.na(univariate_by$by)) {
    if (is.null(univariate_by$terms))
      univariate_by$terms <- "subset"
  }
  
  if (is.null(group_arg$groupvar))
    group_arg$groupvar <- NULL
  if (is.null(group_arg$by))
    group_arg$by <- NULL
  if (is.null(group_arg$cor))
    group_arg$cor <- "un"
  if (is.null(group_arg$dist))
    group_arg$dist <- "gaussian"
  
  if (is.null(sigma_group_arg$groupvar))
    sigma_group_arg$groupvar <- NULL
  if (is.null(sigma_group_arg$by))
    sigma_group_arg$by <- NULL
  if (is.null(sigma_group_arg$cor))
    sigma_group_arg$cor <- "un"
  if (is.null(sigma_group_arg$dist))
    sigma_group_arg$dist <- "gaussian"
  
  
  
  multivariate$verbose <-
    univariate_by$verbose <- group_arg$verbose <- verbose
  
  sigma_group_arg$verbose <- verbose
  
  
  # Temporary placeholder for the number of response for univariate_by
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    temp_ <- univariate_by$by
    if (!temp_ %in% colnames(data)) {
      stop2c(
        paste(
          "\nvariable '",
          temp_,
          "' used for setting univariate_by-univariate submodels is missing"
        )
      )
    }
    if (!is.factor(data[[temp_]])) {
      stop2c("The 'univariate_by' variable '", 
             temp_, "' should be a factor variable")
    }
    nlevtemp_ <- nlevels(data[[temp_]])
    nys <- nlevtemp_
  }
  
  
  # Perform checks and set-up the 'to convert arguments' 
  to_list_if_not <- function(.x, nys, arguments, ...) {
    if (nys == 1) {
      if (!is.symbol(arguments[[.x]]) & !is.character(arguments[[.x]])) {
        arguments[[.x]] <- deparse_0(arguments[[.x]])
      } else {
        arguments[[.x]] <- arguments[[.x]]
      }
      if (is.symbol(arguments[[.x]]) &
          !is.character(arguments[[.x]])) {
        arguments[[.x]] <- deparse_0(arguments[[.x]])
      } else {
        arguments[[.x]] <- arguments[[.x]]
      }
      if (!is.character(.x)) {
        .xx <- eval(parse(text = .x))
      } else {
        .xx <- .x
      }
    }
    if (nys > 1) {
      .xx <- .x
    }
    
    if (!is.character(.xx) & !is.list(.xx)) {
      .xx <- deparse_0(.xx)
    } else {
      .xx <- .xx
    }
    . <- lapply(.xx, function(x)
      if (is.list(x))
        x <- x
      else
        x <- list(x)[[1]])
    assign(.x, ., envir = parent.frame())
  }
  
  eval_c_list_args <- function(.x, nys, arguments, return_args = FALSE, ...) {
    if (is.language(arguments[[.x]]) &
        (strsplit(deparse_0(arguments[[.x]]), "\\(")[[1]][1] !=
         "c" &
         strsplit(deparse_0(arguments[[.x]]), "\\(")[[1]][1] != "list")) {
      arguments[[.x]] <- deparse(arguments[[.x]])
    } else {
      arguments[[.x]] <- (arguments[[.x]])
    }
    if (is.logical(arguments[[.x]])) {
      arguments[[.x]] <- deparse_0(arguments[[.x]])
    } else {
      arguments[[.x]] <- arguments[[.x]]
    }
    if (is.numeric(arguments[[.x]])) {
      arguments[[.x]] <- deparse_0(arguments[[.x]])
    } else {
      arguments[[.x]] <- (arguments[[.x]])
    }
    .xo <- .x
    .x <- arguments[[.x]]
    fun_ <- function(.x) {
      if (!is.character(.x))
        .x <- deparse_0(.x)
      else
        .x <- .x
      .x <- gsub("\\s", "", .x)
    }
    if (is.symbol(arguments[[.xo]]))
      .x <- deparse_0(.x)
    else
      .x <- .x
    if (is.symbol(arguments[[.xo]]))
      args_s <- mapply(fun_, .x)
    if (!is.symbol(arguments[[.xo]]))
      args_s <- mapply(fun_, .x)[-1]
    if (is.character(arguments[[.xo]]))
      args_s <- mapply(fun_, .x)
    if (length(args_s) > 1)
      args_s <- mapply(fun_, .x)[-1]
    attr(args_s, "names") <- NULL
    if (length(args_s) < nys)
      args_s <- rep(args_s, nys)
    if (length(args_s) > nys)
      args_s <- args_s[1:nys]
    if(return_args) {
      return(args_s)
    }
    assign(paste0(.xo, "s"), args_s, envir = parent.frame())
  }
  
  

  getArgNames <-
    function(value)
      methods::formalArgs(deparse_0(substitute(value)[[1]]))
  
  convert_to_list <- getArgNames(bsitar())
  
  # enverr. <- parent.frame()
  for (ip in convert_to_list) {
    if (grepl("_init_", ip)) {
      assign('err.', FALSE, envir = enverr.)
      tryCatch(
        expr = {
          out <- suppressWarnings(ept(ip))
        },
        error = function(e) {
          assign('err.', TRUE, envir = enverr.)
        }
      )
      err. <- get('err.', envir = enverr.)
      if (!err.) {
        if (length(out) > 1 & !is.list(out)) {
          stop2c(
            "Initials specified as vector [e.g, c(1, 2)] but must be a list, ",
            "\n ",
            " Note, initials can also be specified by using a single character",
            "\n ",
            " such as 0, random, or an object defined in the init_data",
            "\n ",
            " please check the following init arg: ",
            ip
          )
        }
      }
    }
  }
  
  
  # Convert arguments to the required format for setting sub-options 
  single_args <- c(
    "data",
    "group_arg",
    "sigma_group_arg",
    "univariate_by",
    "multivariate",
    "prior_data",
    "init_data",
    "init_custom",
    "jitter_init_beta",
    "jitter_init_sd",
    "jitter_init_cor",
    "expose_function",
    "verbose",
    "normalize",
    "seed",
    "brms_arguments",
    "get_stancode",
    "get_standata",
    "get_formula",
    "get_stanvars",
    "get_priors",
    "get_priors_eval",
    "validate_priors",
    "get_init_eval",
    "set_self_priors",
    "add_self_priors",
    "set_replace_priors",
    "set_same_priors_hierarchy",
    "outliers",
    "select_model",
    "decomp",
    "parameterization",
    "custom_family",
    "custom_formula",
    "custom_prior",
    "custom_stanvars",
    'pathfinder_args',
    'pathfinder_init',
    'data_custom',
    'genquant_xyadj',
    "fast_nsk",
    "global_args",
    "sum_zero",
    "...",
    "data2")
  
  for (i in convert_to_list) {
    if (!i %in% single_args) {
      to_list_if_not(i, nys, arguments)
    }
  }
  
  # This will assign paste0(..., 's') -> return_args  = FALSE
  for (i in convert_to_list) {
    if (!i %in% single_args) {
      eval_c_list_args(i, nys, arguments, return_args = FALSE)
      # eval_c_list_args(i, nys, arguments)
    }
  }
  
  # Create list for later use in data_custom -> return_args  = TRUE
  # 24.02.2025
  eval_c_list_args_data_custom <- list()
  for (i in convert_to_list) {
    if (!i %in% single_args) {
      eval_c_list_args_data_custom[[i]] <- 
        eval_c_list_args(i, nys, arguments, return_args = TRUE)
    }
  }
  
  less_args <- extra_args <- c()
  outcomes_l <- paste0(" (", paste(ys, collapse = ", "), ")")
  for (i in convert_to_list) {
    .x <- i
    if (is.call(arguments[[.x]]))
      nl <- length(arguments[[.x]]) - 1
    if (!is.call(arguments[[.x]]))
      nl <- length(arguments[[.x]])
    if (nl > nys)
      extra_args <- c(extra_args, .x)
    if (nl < nys)
      less_args <- c(less_args, .x)
  }
  
  
  if (verbose) {
    setmsgtxt <- paste0("\n Preparing data")
    if (displayit == 'msg') {
      message2c(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  if(is.list(xfuns) & length(xfuns) == 0) {
    xfuns <- rep('NULL', length(ys))
  }
  if(is.list(yfuns) & length(yfuns) == 0) {
    yfuns <- rep('NULL', length(ys))
  }
  
  if(is.list(sigmaxfuns) & length(sigmaxfuns) == 0) {
    sigmaxfuns <- rep('NULL', length(ys))
  }
  
  # For model_info
  xfuns_user      <- xfuns
  yfuns_user      <- yfuns
  sigmaxfuns_user <- sigmaxfuns
  
  if(!is.null(outliers)) {
    if(is.null(outliers$remove))    outliers$remove <- TRUE
    if(is.null(outliers$icode))     outliers$icode <- c(4,5,6)
    if(is.null(outliers$limit))     outliers$limit <- 5
    if(is.null(outliers$velpower))  outliers$velpower <- 0.5
    if(is.null(outliers$lag))       outliers$lag <- 1
    if(is.null(outliers$linearise)) outliers$linearise <- FALSE
    if(is.null(outliers$verbose))   outliers$verbose <- FALSE
  }
  
  check_sigmax_not_c <- arguments[['sigmax']] %>% deparse()
  check_sigmax_not_c <- paste(gsub_space(check_sigmax_not_c), collapse = "")
  if(grepl("c(", check_sigmax_not_c, fixed = T)) {
    stop2c("Argument 'sigmax' should be a list")
  }
  
  if(is_emptyx(sigmaxs)) {
    sigmaxs <- NA
  }
  
  # Now if sigmax -> sigmaxs = FALSE, then no xs will be set as sigmax
  sigmaxs <- check_and_replace_sort_to_full(str = sigmaxs,
                                   x = c("T", "F", "FALSE", "NULL", "NA"),
                                   # what = c("TRUE", "NA", "NA", "TRUE", "NA"), 
                                   what = c("TRUE", "NA", "NA", "NA", "NA"), 
                                   allowed_left = "(^|[^[:alnum:]])",
                                   allowed_right = "($|[^[:alnum:]])"
                                   )
  
  # prepare_data2, when 'univariate_by', first run is to get names
  prepare_data_args <- list()
  prepare_data_args[['data']]          <- data
  prepare_data_args[['xvar']]          <- xs
  prepare_data_args[['yvar']]          <- ys
  prepare_data_args[['idvar']]         <- ids
  prepare_data_args[['univariate_by']] <- univariate_by
  prepare_data_args[['multivariate']]  <- multivariate
  prepare_data_args[['outliers']]      <- outliers
  prepare_data_args[['envir']]         <- enverr.
  prepare_data_args[['sigmaxvar']]     <- sigmaxs
  prepare_data_args[['nys']]           <- nys
  prepare_data_args[['subset']]        <- FALSE
  prepare_data_args[['returnys']]      <- FALSE
  prepare_data_args[['verbose']]       <- FALSE
  prepare_data_args[['displayit']]     <- displayit
  prepare_data_args[['setcolb']]       <- setcolb

  # Imp 
  # When sigmax is not specified in the call, then prepare_data2 will
  # automatically generate sigmax based on xs by adding sigma as prefix
  
  # .org.in will be used for data_custom
  data.org.in    <- data
  xs.org.in      <- xs
  ys.org.in      <- ys
  ids.org.in     <- ids
  sigmaxs.org.in <- sigmaxs
  
  data          <- CustomDoCall(prepare_data2, prepare_data_args)
  xs            <- attr(data, "xs")
  ys            <- attr(data, "ys")
  ids           <- attr(data, "ids")
  sigmaxs       <- attr(data, "sigmaxs")
  subindicators <- attr(data, "subindicators")
  
  check_variable_numeric_exists(data, c(xs, ys))
  
  ###########################################################
  ###########################################################
  dataout <- priorlist <- NULL

  bflist <- list()
  initialslist <- initialslist_s <- 
    prior_stanvarlist <- auxillary_stanvarlist <-
    data_stanvarlist <- bflist
  
  set_model_sigma_by_mu_fun_str_c <- sigmatau_strsi_c <- list()
  sigmaspfncname_c <- list()
  add_identityfun_c <- list()
  

  funlist <- c()
  xoffsetvaluelist <- xoffsetnamelist <- knotsvaluelist <- funlist
  knotsnamelist <- spfun_collect <- funlist
  fixedvaluelist <- fixednamelist <- funlist
  randomvaluelist <- randomnamelist <- groupvarvaluelist <- funlist
  yvarvaluelist <- ynamelist <- covvaluelist <- covnamelist <- funlist
  groupvarnamelist <- xvarvaluelist <- xnamelist <- funlist
  hierarchicalvarnamelist <- hierarchicalvarvaluelist <- funlist
  
  idvarvaluelist <- idnamelist <- funlist
  sigmaidvarvaluelist <- sigmaidnamelist <- funlist

  sigmaxoffsetvaluelist <- sigmaxoffsetnamelist <- funlist
  sigmaxvarvaluelist    <- sigmaxnamelist <- funlist
  sigmacovnamelist      <- sigmacovvaluelist <- funlist
  setsigmaxvarvaluelist <- setsigmaxvarnamelist <- funlist

  sigma_groupvarnamelist <- sigma_groupvarvaluelist <- funlist
  sigma_hierarchicalvarnamelist <- sigma_hierarchicalvarvaluelist <- funlist
  
  funlist_r <- funlist_rnamelist <- funlist_rvaluelist <- list()
  include_fun_nameslist <- funlist_r
  include_fun_nameslist_rnamelist <- funlist_r
  include_fun_nameslist_rvaluelist <- funlist_r

  gq_funs <- list()
  spfncname_c <- c()
  d_adjustedvaluelist <- d_adjustednamelist <- funlist
  SplineCallvaluelist <- SplineCallnamelist <- funlist

  ####
  sigmaspfun_collect <- funlist
  sigmafunlist_r <- sigmafunlist_rnamelist <- sigmafunlist_rvaluelist <- list()
  sigmafunlist <- funlist
  sigmafunlist_r <- funlist_r
  sigmagq_funs <- gq_funs
  sigmafixedvaluelist <- sigmafixednamelist <- funlist
  sigmarandomvaluelist <-sigmarandomnamelist <-sigmagroupvarvaluelist <-funlist
  sigmad_adjustedvaluelist <- sigmad_adjustednamelist <- funlist
  
  sigmavarspfncname_c <- list()
  sigmavarfunlist_rnamelist <- sigmavarfunlist_rvaluelist <- list()
  sigmavarfunlist <- funlist
  sigmavarfunlist_r <- funlist_r
  sigmavargq_funs <- gq_funs
  
  sigmabasicfunlist_rnamelist <- sigmabasicfunlist_rvaluelist <- list()
  sigmabasicgq_funs <- gq_funs
  
  sigmabasicfunlist <- funlist
  sigmabasicfunlist_r <- funlist_r
  
  sigmabasicfunnamevaluelist    <- sigmabasicfunnamenamelist <- funlist
  sigmabasicfunattrvaluelist    <- sigmabasicfunattrnamelist <- funlist
  
  sigmamodelnamevaluelist    <- sigmamodelnamenamelist <- funlist
  
  xfunvaluelist      <- xfunnamelist      <- funlist
  yfunvaluelist      <- yfunnamelist      <- funlist
  sigmaxfunvaluelist <- sigmaxfunnamelist <- funlist

  xfuntransformvaluelist  <- xfuntransformnamelist  <- funlist
  ixfuntransformvaluelist <- ixfuntransformnamelist <- funlist
  
  xfuntransform2valuelist  <- xfuntransform2namelist  <- funlist
  ixfuntransform2valuelist <- ixfuntransform2namelist <- funlist

  yfuntransformvaluelist <- yfuntransformnamelist <- funlist
  iyfuntransformvaluelist <- iyfuntransformnamelist <- funlist

  sigmaxfuntransformvaluelist <- sigmaxfuntransformnamelist <- funlist
  sigmaixfuntransformvaluelist <- sigmaixfuntransformnamelist <- funlist
  
  
  sigmaxfuntransform2valuelist <- sigmaxfuntransform2namelist <- funlist
  sigmaixfuntransform2valuelist <- sigmaixfuntransform2namelist <- funlist

  ###########################################################
  ###########################################################
  
  # Start loop over response
  for (ii in 1:length(ys)) {
    if (nys > 1)
      resp <- ys[ii]
    else
      resp <- ""
    
    subindicatorsi <- subindicators[ii]
    
    # For multivariate and univariate_by models, over ride smat_preH to FALSE
    # This ensures that appropriate function is constructed for different df
    if(nys > 1) {
      smat_preH <- 0
      if(verbose) {
        message2c("The 'smat_preH' is set as 'FLASE' for multivariate model")
      }
    }
    
    # Define function names, moved here up now
    # Let spfncname name be common without response, which later added
    spfncname        <- paste0(toupper(select_model), "", 'Fun')
    spfncname_common <- spfncname
    
    for (i in convert_to_list) {
      if (!i %in% single_args) {
        assign(paste0(i, "s", "i"), eval(parse(text = paste0(i, "s")))[ii])
      }
    }
    
    if (is.null(group_arg$groupvar)) {
      group_arg$groupvar <- idsi
    }
      
    if (is.null(sigma_group_arg$groupvar))
      sigma_group_arg$groupvar <- idsi
    if (!is.numeric(ept(dfsi)) & !is.numeric(ept(knotssi))) {
      stop2c("Either 'df' or 'knots' must be specified")
    }
    if (is.numeric(ept(dfsi)) & is.numeric(ept(knotssi))) {
      # stop2c("Both 'df' and 'knots' are specified. Specify one of them\n")
      dfsi <- 'NULL'
      if(verbose) {
        message2c("The user specified knots are used, hence",
                " the df argument ignored")
      }
    }
    
    if (!is.numeric(ept(sigmadfsi)) & !is.numeric(ept(sigmaknotssi))) {
      stop2c("Either df or knots must be specified for sigma")
    }
    if (is.numeric(ept(sigmadfsi)) & is.numeric(ept(sigmaknotssi))) {
      # stop2c("Both df and knots specified. Specify one of them\n")
      dfsi <- 'NULL'
      if(verbose) {
        message2c("The user specified knots are used for sigma, hence",
                " the df argument ignored")
      }
    }
    
    for (agsxi in letters[1:26]) {
      if(is.null(arguments[[paste0(agsxi, "_", "formula" , "")]])) {
        assign(paste0(agsxi, "_", "formula",        "si") , NULL)
        assign(paste0(agsxi, "_", "formula_gr",     "si") , NULL)
        assign(paste0(agsxi, "_", "formula_gr_str", "si") , NULL)
        assign(paste0(agsxi, "_", "prior_beta",     "si") , NULL)
        assign(paste0(agsxi, "_", "cov_prior_beta", "si") , NULL)
        assign(paste0(agsxi, "_", "prior_sd",       "si") , NULL)
        assign(paste0(agsxi, "_", "cov_prior_sd",   "si") , NULL)
        assign(paste0(agsxi, "_", "init_beta",      "si") , NULL)
        assign(paste0(agsxi, "_", "cov_init_beta",  "si") , NULL)
        assign(paste0(agsxi, "_", "init_sd",        "si") , NULL)
        assign(paste0(agsxi, "_", "cov_init_sd",    "si") , NULL)
      }
    }
    
    validate_fixed_random_parms <- function(fixedsi, 
                                            randomsi, 
                                            allowed_parm_letters, 
                                            select_model,
                                            match_sitar_d_form) {
      
      parm_letters_fixed <- strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
      parm_letters_fixed <- sort(parm_letters_fixed)
      parm_letters_fixed <- parm_letters_fixed[1:length(allowed_parm_letters)]
      parm_letters_fixed <- parm_letters_fixed[!is.na(parm_letters_fixed)]
      
      parm_letters_random <- strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
      parm_letters_random <- sort(parm_letters_random)
      parm_letters_random <- parm_letters_random[1:length(allowed_parm_letters)]

      if(select_model == 'pb1' | 
         select_model == 'pb2' | 
         select_model == 'pb3' |
         select_model == 'logistic1' |
         select_model == 'logistic2' |
         select_model == 'logistic3' 
         ) {
        if(length(parm_letters_fixed) != length(allowed_parm_letters))
          stop2c("For model '", select_model, "'", ", 
               the number of parameters must be ",
               length(allowed_parm_letters),
               " \n ", 
               "(parameters ", 
               paste(paste0("'", allowed_parm_letters, "'"), collapse = " "),
               ")"
          )
      }
      
      if(select_model == 'sitar') {
        if(length(parm_letters_fixed) > length(allowed_parm_letters))
          stop2c("For model '", select_model, "'", ", 
               the maximum number of parameters is ",
               length(allowed_parm_letters),
               " \n ", 
               "(parameters ", 
               paste(paste0("'", allowed_parm_letters, "'"), collapse = " "),
               ")"
          )
      }
      
      get_parm_letters <- parm_letters_fixed
      sub_parm_letters_fixed <- intersect(allowed_parm_letters, 
                                          get_parm_letters)
      sub_parm_letters <- sub_parm_letters_fixed
      inv_parm_letters <- get_parm_letters[!get_parm_letters %in% 
                                             sub_parm_letters]
      inv_parm_letters <- sort(inv_parm_letters)
      if(length(inv_parm_letters) > 0) {
        see_what_formual <- paste0(" Please see and correctly ", "'", 
                                   'fixed', "'", " argument")
        not_allowed_parsm <- paste(paste0("'", inv_parm_letters, "'"), 
                                   collapse = " ")
        msg_1 <- paste0("Parameter ", not_allowed_parsm, 
                        " not allowed for ", "'", select_model, "'", " model")
        msg_2 <- paste0(" Allowed parameters are ", 
                        paste(paste0("'", 
                                     allowed_parm_letters, "'"), 
                              collapse = " "))
        stop2c(msg_1, "\n ", msg_2, " \n ", see_what_formual)
      }
      
      get_parm_letters <- parm_letters_random
      sub_parm_letters_random <- intersect(allowed_parm_letters, 
                                           get_parm_letters)
      sub_parm_letters <- sub_parm_letters_random
      inv_parm_letters <- get_parm_letters[!get_parm_letters %in% 
                                             sub_parm_letters]
      inv_parm_letters <- sort(inv_parm_letters)
      if(length(inv_parm_letters) > 0) {
        see_what_formual <- paste0(" Please see and correctly ", "'", 
                                   'random', "'", " argument")
        not_allowed_parsm <- paste(paste0("'", inv_parm_letters, "'"), 
                                   collapse = " ")
        msg_1 <- paste0("Parameter ", 
                        not_allowed_parsm, " not allowed for ", "'", 
                        select_model, "'", " model" )
        msg_2 <- paste0(" Allowed parameters are ", 
                        paste(paste0("'", 
                                     allowed_parm_letters, "'"), 
                              collapse = " "))
        stop2c(msg_1, "\n ", msg_2, " \n ", see_what_formual)
      }
      
      sub_parm_letters_fixed_random <- intersect(parm_letters_fixed, 
                                                 parm_letters_random)
      inv_parm_letters_fixed_random <- 
        parm_letters_random[!parm_letters_random %in% parm_letters_fixed]
      inv_parm_letters_fixed_random <- sort(inv_parm_letters_fixed_random)
      
      if(length(inv_parm_letters_fixed_random) > 0) {
        not_allowed_parsm <- paste(paste0("'", 
                                          inv_parm_letters_fixed_random, "'"), 
                                   collapse = " ")
        msg_mismatch_fixed_random_str <- 
          paste0(
            "Parameter ", not_allowed_parsm, " included in the random part of ",
            "\n ",
            " the model but missing from the fixed effects.",
            "\n ",
            " Please check and correctly specify the",
            "\n ",
            "'fixed'/'random' arguments for the ",
            toupper(select_model), " model."
          )
        
        if(select_model == 'sitar') {
          if(!match_sitar_d_form) stop2c(msg_mismatch_fixed_random_str)
        } else {
          stop2c(msg_mismatch_fixed_random_str)
        }
      } # if(length(inv_parm_letters_fixed_random) > 0) {
      
      sub_parm_letters_fixed <- sort(sub_parm_letters_fixed)
      sub_parm_letters_random <- sort(sub_parm_letters_random)
      
      out_fixed <- paste(sub_parm_letters_fixed, collapse = "+")
      out_random <- paste(sub_parm_letters_random, collapse = "+")
      list(fixed = out_fixed, random = out_random)
    } # validate_fixed_random_parms
    
    # Over ride when restricting to abcd
    if(override_select_model) {
      if(grepl("d", fixedsi) & grepl("d", randomsi)) {
        sitar_nparms <- 4
        # match_sitar_d_form <- FALSE
      } else if(grepl("d", fixedsi) & !grepl("d", randomsi)) {
        sitar_nparms <- 4
        # match_sitar_d_form <- FALSE
      } else if(!grepl("d", fixedsi) & grepl("d", randomsi)) {
        sitar_nparms <- 4
        match_sitar_d_form <- TRUE
      } else if(!grepl("d", fixedsi) & !grepl("d", randomsi)) {
        sitar_nparms <- 3
        # match_sitar_d_form <- FALSE
      }
    }
    
    # Model specific number of fixed and random parameters
    allowed_parm_letters <- NULL
    # covers all sitar models
    if(grepl("^sitar", select_model)) {
      allowed_parm_letters <- letters[1:sitar_nparms]
    }
    # if(select_model == 'sitar') allowed_parm_letters <- letters[1:sitar_nparms]
    if(select_model == 'pb1')   allowed_parm_letters <- letters[1:5]
    if(select_model == 'pb2')   allowed_parm_letters <- letters[1:6]
    if(select_model == 'pb3')   allowed_parm_letters <- letters[1:6]
    if(select_model == 'logistic1')   allowed_parm_letters <- letters[1:3]
    if(select_model == 'logistic2')   allowed_parm_letters <- letters[1:6]
    if(select_model == 'logistic3')   allowed_parm_letters <- letters[1:9]
    
    if(select_model == 'rcs')   allowed_parm_letters <- letters[1]
    
    fixedsi_randomsi <- 
      validate_fixed_random_parms(
        fixedsi = fixedsi, 
        randomsi = randomsi,
        allowed_parm_letters = allowed_parm_letters, 
        select_model = select_model,
        match_sitar_d_form = match_sitar_d_form)
    
    fixedsi <- fixedsi_randomsi[['fixed']]
    randomsi <- fixedsi_randomsi[['random']]
    
    abc_fixedsi  <- strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
    abc_randomsi <- strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
    
    # Covariate not allowed when matching to sitar 'd' form
    if(select_model == 'sitar') {
      if (!match_sitar_d_form) {
        if (!grepl("d", fixedsi, fixed = T) &
            grepl("d", randomsi, fixed = T)) {
          stop2c(
            "Parameter 'd' is missing in the fixed effects part of the model ",
            "\n ",
            " but specified in the random effects part of the model ",
            "\n ",
            " Either include 'd' in the fixed effects too or else ",
            "\n ",
            " remove it from the random effect part of the model"
          )
        }
      }
      
      # New d as random cov
      d_as_random_only_cov <- FALSE
      if ((grepl("d", fixedsi, fixed = T) |
           grepl("d", randomsi, fixed = T)) &
          (!grepl("^~1$", d_formulasi) |
           !grepl("^~1$", d_formula_grsi))) {
        d_as_random_only_cov <- TRUE
      }

      if (match_sitar_d_form) {
        if ((grepl("d", fixedsi, fixed = T) |
             grepl("d", randomsi, fixed = T)) &
            (!grepl("^~1$", d_formulasi) |
             !grepl("^~1$", d_formula_grsi) |
             !grepl("^~1$", d_formula_gr_strsi))) {
          # stop2c(
          #   "Parameter 'd' is missing in the fixed effects part of the model ",
          #   "\n ",
          #   " but specified in the random effects part of the model ",
          #   "\n ",
          #   " (This is to match with the 'sitar' package's formulation)",
          #   "\n ",
          #   " For this current formulation ",
          #   "\n ",
          #   " covariate(s) are not allowed")
        }
      } # if (match_sitar_d_form) {
    } # if(select_model == 'sitar') {
    
    if(select_model == 'sitar') {
      if(!any(grepl('s', fixedsi))) fixedsi <- paste0(fixedsi, "+", "s")
    }
    
    # New
    d_as_random_only <- FALSE
    if(select_model == 'sitar') {
      if(!grepl("d", fixedsi, fixed = T) & 
         grepl("d", randomsi, fixed = T)) {
        d_as_random_only <- TRUE
      }
      
      if(d_as_random_only) {
        if(d_as_random_only) d_formulasi <- "~0"
        match_sitar_d_form <- TRUE
      }
      
      # Set all _gr to NULL if parameter is not random
      for (i in letters[1:26]) {
        if(!grepl(i, randomsi, fixed = T)) {
          assign(paste0(i, '_formula_grsi'), NULL)
          assign(paste0(i, 'formula_gr_strsi'), NULL)
        }
      }
      
    } # if(select_model == 'sitar') {
    
    
    
    # for (i in letters[1:26]) {
    #   if(!grepl(i, randomsi, fixed = T)) {
    #     assign(paste0(i, '_formula_grsi'), NULL)
    #     assign(paste0(i, 'formula_gr_strsi'), NULL)
    #   }
    # }
    
   
    
    
    
    
    if(select_model == 'rcs') {
      if(!any(grepl('s', fixedsi))) fixedsi <- paste0(fixedsi, "+", "s")
      if(!any(grepl('s', randomsi)) & rcs_add_re_spline) {
          randomsi <- paste0(randomsi, "+", "s")
      }
      if(any(grepl('s', randomsi)) & !rcs_add_re_spline) {
        stop2c("you have specified select_model = 'rcsf' (i.e., no random",
             "\n ",
             "spline effects) but your random argument have parameter 's'.",
             "\n ",
             "Please check your 'select_model' and 'random' arguments")
      }
    }
    
    ###########################################################################
    # sigma_formula_manual
    ###########################################################################
    # This below for sigma_formula_manual
    # The objective is set sigma_formulasi and sigma_formula_gr_strsi 
    # by extracting relevant portion. These 'sigma_formulasi' and 
    # 'sigma_formula_gr_strsi' formuale are used for prior setting and covars
    # Also, 'add_default_args_to_nlf_lf' can be used for 'dpar_formuala'
    # Add missing parameters to the sigma_formula_manual
    # This check might be needed for dpar_formual
    
    # set_model_sigma_by_fz -> nlme::varExp( form ~ fitted(.)) - negzero
    # set_model_sigma_by_fp -> nlme::varPower( form ~ fitted(.))
    # set_model_sigma_by_fe -> nlme::varExp( form ~ fitted(.))
    # set_model_sigma_by_ve -> nlme::varExp()
    # set_model_sigma_by_vp -> nlme::varPower()
    # set_model_sigma_by_mp -> mean described in brms manual - sqrt(varPower)
    # set_model_sigma_by_me -> mean described in brms manual - sqrt(varExp)
    # set_model_sigma_by_ls -> location scale via sigmafunction()
    # set_model_sigma_by_rp -> nlme::varPower( form ~ fitted(.)) residual
    # set_model_sigma_by_re -> nlme::varExp( form ~ fitted(.)) residual
    
    sigma_formula_manualsi <- paste(gsub_space(sigma_formula_manualsi), 
                                    collapse = "")
    sigma_formula_manualsi <- gsub("\"" , "'", 
                                   sigma_formula_manualsi, fixed = T)
    
   
    add_arg_to_sigma_formula_manual <- function(x, arg, what) {
      gsub_it <- replace_string_part(x = x, 
                                     start = "nlf(", 
                                     end = ")",  replace = "",
                                     extract = T, cat_str = FALSE, 
                                     exclude_start = FALSE, 
                                     exclude_end = T)
      gsub_by <- paste0(gsub_it, ",", paste0(arg, "=", what), "")
      out <- gsub(gsub_it, gsub_by, x, fixed = T)
      return(out)
    }
    
    sigma_formula_manualsi_set <- FALSE
    if(sigma_formula_manualsi != "NULL" |
       grepl("nlf(", sigma_formula_manualsi, fixed = TRUE)) {
      if(grepl("method=", sigma_formula_manualsi)) {
        sigma_formula_manualsi_set <- TRUE
      } else if(grepl("method=none", sigma_formula_manualsi)) {
        sigma_formula_manualsi_set <- FALSE
      } else if(grepl("method=no", sigma_formula_manualsi)) {
        sigma_formula_manualsi_set <- FALSE
      } else {
        sigma_formula_manualsi <-
          add_arg_to_sigma_formula_manual(x = sigma_formula_manualsi,
                                          arg = "method",
                                          what = "basic")
        sigma_formula_manualsi_set <- TRUE
        # sigma_formula_manualsi_set <- FALSE
      }
    }
    

    sigma_formula_manualsi <- paste0(gsub_space(sigma_formula_manualsi), 
                                     collapse = "")
    
    # Note that brms complains of duplicate names  Hmisc::rcspline.ev... 
    # when T/TRUE is used in the below form. 
    # The F/FALSE does not result in the same error as T/TRUE
    # allowed_left could be "="
    sigma_formula_manualsi <- 
    check_and_replace_sort_to_full(str = sigma_formula_manualsi,
                           x = c("T", "F"),
                           what = c("TRUE", "FALSE"), 
                           allowed_left = "(^|[^[:alnum:]])", # = 
                           allowed_right = "($|[^[:alnum:]])")
    
    
    if(!sigma_formula_manualsi_set) {
      set_model_sigma_by_ba <- FALSE
      set_model_sigma_by_fz <- FALSE
      set_model_sigma_by_fp <- FALSE
      set_model_sigma_by_fe <- FALSE
      set_model_sigma_by_ve <- FALSE
      set_model_sigma_by_vp <- FALSE
      set_model_sigma_by_cp <- FALSE
      set_model_sigma_by_mp <- FALSE
      set_model_sigma_by_me <- FALSE
      set_model_sigma_by_ls <- FALSE
      set_model_sigma_by_rp <- FALSE
      set_model_sigma_by_re <- FALSE
      expose_sigma_ls_fun   <- FALSE
      add_identityfun       <- FALSE
      sigmamodelsi          <- NULL
      sigmatau_strsi        <- NULL
      sigmabasicfunnamesi   <- NULL
      sigmabasicfunattrsi   <- NULL
      sigmavarspfncname_temp  <- NULL
      sigmavarspfncname <- NULL
      sigma_formula_manual_prior_via_sigma_formula <- FALSE
    }
   
    
    if(sigma_formula_manualsi_set) {
      set_model_sigma_by_ba <- FALSE
      set_model_sigma_by_fz <- FALSE
      set_model_sigma_by_fp <- FALSE
      set_model_sigma_by_fe <- FALSE
      set_model_sigma_by_ve <- FALSE
      set_model_sigma_by_vp <- FALSE
      set_model_sigma_by_cp <- FALSE
      set_model_sigma_by_mp <- FALSE
      set_model_sigma_by_me <- FALSE
      set_model_sigma_by_ls <- FALSE
      set_model_sigma_by_rp <- FALSE
      set_model_sigma_by_re <- FALSE
      expose_sigma_ls_fun   <- FALSE
      add_identityfun       <- FALSE
      sigmamodelsi          <- NULL
      sigmatau_strsi        <- NULL
      sigmabasicfunnamesi   <- NULL
      sigmabasicfunattrsi   <- NULL
      sigmavarspfncname_temp  <- NULL
      sigmavarspfncname <- NULL
      sigma_formula_manual_prior_via_sigma_formula <- FALSE
      ##########################################################################
      # Get sigma method and prior arg
      # for nys > 1, get below ars once only otherwise it will lead to ""
      # another option is to create list.. [[ii]] but we are replacing priors
      # for sigma_formula_manualsi at once for all outcomes. Therefore it is 
      # assumed that same behaviors is expected
      method_nlf_custom_arg_full <- c("varpower", 
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
                                      "residualexp",
                                      "ls")

      allowed_nlf_custom_arg_short <- c("vp", 
                                        "cp",
                                        "ve", 
                                        "fz",
                                        "fe",
                                        "fp",
                                        "mp",
                                        "me",
                                        "rp", 
                                        "re",
                                        "ls")

      method_nlf_custom_arg_full   <- c("basic", method_nlf_custom_arg_full)
      allowed_nlf_custom_arg_short <- c("ba", allowed_nlf_custom_arg_short)
      
      method_nlf_custom_arg <- c(method_nlf_custom_arg_full, 
                                 allowed_nlf_custom_arg_short)
      
      
      prior_nlf_custom_arg  <- c("self", "auto")
      
      
      get_nlf_newstr <- 
        get_nlf_custom_arg(str = sigma_formula_manualsi,
                           search = "method",
                           allowed_nlf_custom_arg = method_nlf_custom_arg,
                           clean = TRUE)
      sigma_formula_manualsi <- get_nlf_newstr[1]
      nlf_sigma_method_arg   <- get_nlf_newstr[2]
      
      
      get_nlf_newstr <- 
        get_nlf_custom_arg(str = sigma_formula_manualsi,
                           search = "prior",
                           allowed_nlf_custom_arg = prior_nlf_custom_arg,
                           clean = TRUE)
      sigma_formula_manualsi <- get_nlf_newstr[1]
      nlf_sigma_prior_arg    <- get_nlf_newstr[2]
      
      nlf_sigma_set_check    <- TRUE
      if(is.na(nlf_sigma_prior_arg)) {
        nlf_sigma_prior_arg <- ""
      } else if(!is.na(nlf_sigma_prior_arg)) {
        if(nlf_sigma_prior_arg == "self") nlf_sigma_set_check <- FALSE
        if(nlf_sigma_prior_arg != "self") nlf_sigma_set_check <- TRUE
      }

      ##########################################################################
      if(grepl("identity(", sigma_formula_manualsi, fixed = TRUE)) {
        add_identityfun <- TRUE
      } else {
        add_identityfun <- FALSE
      }
      add_identityfun_c[[ii]] <- add_identityfun
      ##########################################################################
  
      if(nlf_sigma_method_arg == "basic") {
        set_model_sigma_by_ba <- TRUE
      } else if(nlf_sigma_method_arg == "varpower") {
        set_model_sigma_by_vp <- TRUE
      } else if(nlf_sigma_method_arg == "varconstpower") {
        set_model_sigma_by_cp <- TRUE
      } else if(nlf_sigma_method_arg == "varexp") {
        set_model_sigma_by_ve <- TRUE
      } else if(nlf_sigma_method_arg == "fittedz") {
        set_model_sigma_by_fz <- TRUE
      } else if(nlf_sigma_method_arg == "fitted" |
                nlf_sigma_method_arg == "fittedpower") {
        set_model_sigma_by_fp <- TRUE
      } else if(nlf_sigma_method_arg == "fittedexp") {
        set_model_sigma_by_fe <- TRUE
      } else if(nlf_sigma_method_arg == "mean" |
                nlf_sigma_method_arg == "meanpower") {
        set_model_sigma_by_mp <- TRUE
      } else if(nlf_sigma_method_arg == "meanexp") {
        set_model_sigma_by_me <- TRUE
      } else if(nlf_sigma_method_arg == "residual" |
                nlf_sigma_method_arg == "residualpower") {
        set_model_sigma_by_rp <- TRUE
      } else if(nlf_sigma_method_arg == "residualexp") {
        set_model_sigma_by_re <- TRUE
      } else if(nlf_sigma_method_arg == "ls") {
        set_model_sigma_by_ls <- TRUE
      } else if(nlf_sigma_method_arg == "none") {
        set_model_sigma_by_no <- TRUE
      } else {
        stop2c("method not found in nlf() for sigma_formula_manual")
      }
      
      sigmamodelsi <- nlf_sigma_method_arg
      
      ##########################################################################
      check_if_varname_exact(str = sigma_formula_manualsi,
                             x = xsi,
                             allowed_left = "._",
                             allowed_right = "._")
      
      if(nlf_sigma_method_arg == "ls") {
        if(is.na(sigmaxsi)) {
          stop2c("For location scale model 'ls', the argument 
                 'sigmax' is needed")
        }
      }
      
      ##########################################################################
      getsetform0tilde <- replace_string_part(x = sigma_formula_manualsi,
                            start = "nlf(sigma",
                            end = "(",
                            extract = TRUE, 
                            cat_str = FALSE,
                            exclude_start = FALSE, 
                            exclude_end = FALSE) 
      
      if(grepl("~1+", getsetform0tilde)) {
        stop2c("sigma_formula_manual must not contain intercept as sigma~1+...")
      }
      
      ##########################################################################
      if(set_model_sigma_by_fz |
         set_model_sigma_by_fp |
         set_model_sigma_by_fe |
         set_model_sigma_by_ve | 
         set_model_sigma_by_vp |
         set_model_sigma_by_cp |
         set_model_sigma_by_mp |
         set_model_sigma_by_me |
         set_model_sigma_by_rp |
         set_model_sigma_by_re ) {
        
        # replace placeholder vf with sigmavarfun
        sigmavarspfncname_temp      <- "sigmavarfun"
        sigmavarspfncname_org <- 
          replace_string_part(x = sigma_formula_manualsi,
                              start = "~",
                              end = "(",
                              replace = "",
                              extract = TRUE, 
                              cat_str = FALSE,
                              exclude_start = TRUE, 
                              exclude_end = TRUE) 
        
        sigma_formula_manualsi <- gsub(sigmavarspfncname_org, 
                                       sigmavarspfncname_temp,
                                       sigma_formula_manualsi, fixed = T)
        
        sigmavarspfncname_common    <- sigmavarspfncname_temp
        sigmavarspfncname           <- sigmavarspfncname_temp
        if(nys > 1) {
          sigmavarspfncname      <- paste0(ysi, "_", sigmavarspfncname_temp)
        } 
        # new 
        sigma_formula_manualsi <- gsub(sigmavarspfncname_temp, 
                                       sigmavarspfncname,
                                       sigma_formula_manualsi, fixed = T)
        
      }
      
      
      ##########################################################################
      
      if(set_model_sigma_by_ve |
         set_model_sigma_by_vp |
         set_model_sigma_by_cp) {
        if(grepl('identity()', sigma_formula_manualsi, fixed = T) |
           grepl('()', sigma_formula_manualsi, fixed = T)) {
          stop2c("Expecting predictor as the last argument 
                 for variance function:",
               "\n  ",
               collapse_comma(sigmavarspfncname_org))
        }
      }
      
      ##########################################################################
      if(set_model_sigma_by_ba |
         set_model_sigma_by_fz |
         set_model_sigma_by_fp |
         set_model_sigma_by_fe |
         set_model_sigma_by_ve | 
         set_model_sigma_by_vp |
         set_model_sigma_by_cp |
         set_model_sigma_by_mp |
         set_model_sigma_by_me |
         set_model_sigma_by_rp |
         set_model_sigma_by_re ) {
        sigma_formula_manual_prior_via_sigma_formula <- TRUE
      }
      ##########################################################################
      # add_sigma_by_ls - Extract sigmaspfncname
      if(set_model_sigma_by_ls) {
        expose_sigma_ls_fun <- TRUE
        sigma_formula_manualsi_str_full <- 
          replace_string_part(x = sigma_formula_manualsi,
                              start = "~",
                              end = "(",
                              replace = "",
                              extract = T)
        sigmaspfncname_temp <- gsub("~", "", sigma_formula_manualsi_str_full, 
                                    fixed = T)
        sigmaspfncname_temp <- gsub("(", "", sigmaspfncname_temp, 
                                    fixed = T)
        
        if(!grepl("^sigma", sigmaspfncname_temp)) {
          sigmaspfncname_temp.org <- sigmaspfncname_temp
          sigmaspfncname_temp     <- paste0("sigma", spfncname)
          sigma_formula_manualsi  <- gsub(sigmaspfncname_temp.org, 
                                          sigmaspfncname_temp,
                                          sigma_formula_manualsi, fixed = T)
        }
        
        # replace x predictor using sigmaxsi
        placeholderx <- 
        replace_string_part(x = sigma_formula_manualsi,
                            start = paste0(sigmaspfncname_temp, "("),
                            end = ",",
                            extract = TRUE, 
                            cat_str = FALSE,
                            exclude_start = TRUE, 
                            exclude_end = TRUE) 
        
        sigma_formula_manualsi  <- sub(placeholderx, 
                                        sigmaxsi,
                                        sigma_formula_manualsi, fixed = T)

        
        sigmaspfncname_common  <- sigmaspfncname_temp
        sigmaspfncname_c[[ii]] <- sigmaspfncname_temp
        if(nys > 1) {
          sigmaspfncname         <- paste0(ysi, "_", sigmaspfncname_temp)
          sigma_formula_manualsi <- gsub(sigmaspfncname_temp, sigmaspfncname,
                                         sigma_formula_manualsi, fixed = T)
        } else {
          sigmaspfncname <- sigmaspfncname_temp
        }
      }
      
      
      
      ##########################################################################
      # Now check for distinct predictor for mu and sigma
      check_if_varname_exact(str = sigma_formula_manualsi,
                             x = xsi,
                             allowed_left = "._",
                             allowed_right = "._")
      
      ##########################################################################
      # Add missing parameters to the sigma_formula_manual
      sigma_formula_manualsi <- 
        add_default_args_to_nlf_lf(str = sigma_formula_manualsi, 
                                   nys = nys, 
                                   ysi = ysi, 
                                   check = FALSE,
                                   extract_covar = FALSE,
                                   extract_nlpar = FALSE, 
                                   data_varnames = colnames(data),
                                   verbose = FALSE)
      ##########################################################################
      set_model_sigma_by_mu_fun_str_full <- 
        replace_string_part(x = sigma_formula_manualsi,
                            start = "~",
                            end = ")",
                            replace = "",
                            extract = T)
      
      set_model_sigma_by_mu_fun_str_full <- 
        gsub("~", "", set_model_sigma_by_mu_fun_str_full, fixed = T)
      
      set_model_sigma_by_mu_fun_str <- set_model_sigma_by_mu_fun_str_full
      ##########################################################################
      
      
      ##########################################################################
      # check sigmavarspfncname_common
      if(set_model_sigma_by_fz |
         set_model_sigma_by_fp |
         set_model_sigma_by_fe |
         set_model_sigma_by_ve | 
         set_model_sigma_by_vp |
         set_model_sigma_by_cp |
         set_model_sigma_by_mp |
         set_model_sigma_by_me |
         set_model_sigma_by_rp |
         set_model_sigma_by_re ) {
        # write message 
        msg_for_setting_sigma_var_function <-
          paste0(
            "The sigma formulation for variance modeling must be ",
            "specified via '", sigmavarspfncname_common, "()'.",
            "\n\n",
            "The bsitar package provides six different methods for variance ",
            "modeling, five of which are implemented in the nlme package. ",
            "These methods are:",
            "\n  ",
            "'nlme::varPower()'",
            "\n  ",
            "'nlme::varConstPower()'",
            "\n  ",
            "'nlme::varExp(form ~ x)'",
            "\n  ",
            "'nlme::varExp(form ~ fitted(.))'",
            "\n  ",
            "'nlme::varExp(form ~ resid(.))'",
            "\n\n",
            "These are specified using the 'method' argument within 'nlf()', ",
            "for example: nlf(..., method = 'xx') where 'xx' is the method (see below).",
            "\n\n",
            "The 'method' argument for each of the five nlme approaches ",
            "is as follows (short hands in parentheses):",
            "\n  ",
            "'varpower' ('vp')",
            "\n  ",
            "'varconstpower' ('cp')",
            "\n  ",
            "'varexp' ('ve')",
            "\n  ",
            "'fitted' ('fi')",
            "\n  ",
            "'residual' ('re')",
            "\n\n",
            "In addition, bsitar provides a sixth method based on an ",
            "example from the brms reference manual, which models the ",
            "square root of the fitted values.",
            "\n",
            "This can be specified with method = 'mean' or its short hand 'me'.",
            "\n\n",
            "Below are examples showing how to use '", sigmavarspfncname_common, "' ",
            "to specify each of the six variance models:",
            "\n\n",
            "1. varpower:",
            "\n  ",
            "nlf(sigma ~ vf(param1, param2, predictor), method = 'vp') +",
            "\n  ",
            "lf(param1 + param2 ~ 1)",
            "\n\n",
            "2. varConstPower:",
            "\n  ",
            "nlf(sigma ~ vf(param1, param2, param3, predictor), method = 'cp') +",
            "\n  ",
            "lf(param1 + param2 + param3 ~ 1)",
            "\n\n",
            "3. varExp:",
            "\n  ",
            "nlf(sigma ~ vf(param1, param2, predictor), method = 've') +",
            "\n  ",
            "lf(param1 + param2 ~ 1)",
            "\n\n",
            "4. fitted:",
            "\n  ",
            "nlf(sigma ~ vf(param1, param2, identity()), method = 'fi') +",
            "\n  ",
            "lf(param1 + param2 ~ 1)",
            "\n\n",
            "5. residual:",
            "\n  ",
            "nlf(sigma ~ vf(param1, param2, identity(), response), method = 're') +",
            "\n  ",
            "lf(param1 + param2 ~ 1)",
            "\n\n",
            "6. mean:",
            "\n  ",
            "nlf(sigma ~ vf(param1, param2, identity()), method = 'me') +",
            "\n  ",
            "lf(param1 + param2 ~ 1)",
            "\n\n",
            "Internal Predictor Transformations:",
            "\n  ",
            "- For 'varpower' and 'varConstPower', the predictor is ",
            "transformed to log(abs(predictor)).",
            "\n  ",
            "- For 'varexp', the predictor is not transformed.",
            "\n  ",
            "- For 'fitted' and 'residual', identity() is internally ",
            "set to fitted(.).",
            "\n  ",
            "- For 'mean', identity() is internally set to sqrt(fitted(.)).",
            "\n\n",
            "The linear predictor, lf(), can be extended to include ",
            "covariates and group-level random effects.",
            "\n",
            "For example, 'lf(param1 + param2 ~ 1)' can become:",
            "\n  ",
            "'lf(param1 + param2 ~ 1 + covariate + (1 || gr(id, by = groupid)))'",
            "\n\n",
            "Automatic Prior Assignment:",
            "\n",
            "Priors for these linear predictors are assigned automatically. ",
            "This applies to both mean and group-level random effects. ",
            "The function uses priors specified via arguments like ",
            "'sigma_prior_beta', 'sigma_cov_prior_beta', 'sigma_prior_sd', etc. ",
            "These are the same arguments otherwise used for setting priors on ",
            "parameters defined by 'sigma_formula', 'sigma_formula_gr', and ",
            "'sigma_formula_gr_str'.",
            "\n\n",
            "To disable this automatic prior assignment, you can add ",
            "the argument prior = 'self' to the nlf() function. ",
            "For example:",
            "\n  ",
            "'nlf(..., method = 'vp', prior = 'self')'",
            "\n\n",
            "Note: If you disable automatic prior assignment, it is advised ",
            "to first get the required prior structure by running 'bsitar' ",
            "with the 'get_prior = TRUE' argument. The relevant portions ",
            "of this structure can then be edited and added back to the ",
            "model via the 'add_self_priors' argument in 'bsitar'. Please ",
            "see the documentation for 'add_self_priors' for more details."
          )
        # check and display message as error
        if(!grepl(paste0(sigmavarspfncname_common, "("), 
                  set_model_sigma_by_mu_fun_str, fixed = TRUE)) {
          stop2c(msg_for_setting_sigma_var_function)
        } 
      }
      # End check sigmavarspfncname_common
      
      
      set_model_sigma_by_mu_fun_str_c[[ii]] <- set_model_sigma_by_mu_fun_str
      
      
      sigmatau_strsi <- 
        add_default_args_to_nlf_lf(str = sigma_formula_manualsi, 
                                   nys = nys, 
                                   ysi = ysi, 
                                   check = FALSE,
                                   extract_covar = FALSE,
                                   extract_nlpar = TRUE, 
                                   data_varnames = colnames(data),
                                   verbose = FALSE)
      
      sigmatau_strsi_c[[ii]] <- sigmatau_strsi
      
      
      if(set_model_sigma_by_ba) {
        # replace functions with :: / ::: with _ in 'sigma_formula_manualsi'
        getouttemp<- get_function_names_code_from_string(sigma_formula_manualsi)
        sigma_formula_manualsi <- getouttemp[['str']]
        sigmabasicfunnamesi    <- getouttemp[['name']]
        sigmabasicfunattrsi    <- getouttemp[['attr']]
        # Also, assign those functions to the environment 
        package_env <- as.environment("package:bsitar")
        if(length(getouttemp[['code']] != 0)) {
          intx <- 0
          for (funi in 1:length(getouttemp[['code']])) {
            intx <- intx + 1
            if(sigmabasicfunattrsi[intx] %in% allowed_namespace_for_sigma_d1()) {
              assign(gsub("<-.*$", "", getouttemp[['code']][funi]),
                     ept(getouttemp[['code']][funi]) , envir = package_env )
            }
          }
        }
      }
      
      ##########################################################################
      if(sigma_formula_manual_prior_via_sigma_formula) {
        # Split into sigma_formula and sigma_formula_gr
        sigma_formula_manualsi_for_parms <- sigma_formula_manualsi
        sigma_formula_manualsi_for_parms <- 
          gsub("+lf", ",lf", sigma_formula_manualsi_for_parms, fixed = T)
        sigma_formula_manualsi_for_parms <- 
          paste0("c(", sigma_formula_manualsi_for_parms, ")")
        get_lf_part_sigma_formula_manualsi <- 
          ept(sigma_formula_manualsi_for_parms)
        nthtau <- length(get_lf_part_sigma_formula_manualsi)
        get_lf_part_sigma_formula_manualsi_form <- 
          get_lf_part_sigma_formula_manualsi[[nthtau]]
        sigma_formulasi <- get_lf_part_sigma_formula_manualsi_form %>% deparse()
        sigma_formulasi <- paste0(gsub_space(sigma_formulasi), collapse = "")
        sigma_formulasi_check <- strsplit(sigma_formulasi, "+(", fixed = T)[[1]]
        sigma_formulasi <- sigma_formulasi_check[1]
        if(length(strsplit(sigma_formulasi, "~", fixed = T)[[1]]) > 1) {
          sigma_formulasi <- strsplit(sigma_formulasi, "~", fixed = T)[[1]][-1]
        } else {
          sigma_formulasi <- sigma_formulasi
        }
        sigma_formulasi <- paste0("~", sigma_formulasi)
        if(length(sigma_formulasi_check) == 1) {
          sigma_formula_gr_strsi <- NULL
        } else {
          sigma_formula_gr_strsi <- sub(".*?\\+\\(", "", sigma_formulasi_check)
          sigma_formula_gr_strsi <- 
            sigma_formula_gr_strsi[2:length(sigma_formula_gr_strsi)]
          sigma_formula_gr_strsi <- paste(sigma_formula_gr_strsi, collapse= "+")
          sigma_formula_gr_strsi <- 
            paste0(gsub_space(sigma_formula_gr_strsi), collapse = "")
          sigma_formula_gr_strsi <- paste0("(", sigma_formula_gr_strsi)
        }
        if(is.null(sigma_formula_gr_strsi)) {
          sigma_formula_gr_strsi <- "NULL"
        }
      } # if(sigma_formula_manual_prior_via_sigma_formula) {
      ##########################################################################
    } # if(sigma_formula_manualsi_set) {
    
    
    ###########################################################################
    # end of sigma_formula_manual - # add_sigma_by_mu
    ###########################################################################
    # see if above function can be used
    # Add missing parameters to the dpar_formula
    if (!is.null(dpar_formulasi)) {
      if (grepl("^1$", dpar_formulasi)) {
        dpar_formulasi <- paste0("lf(", "sigma", "~", dpar_formulasi, ")")
      } else if (grepl("^~1", dpar_formulasi)) { 
        dpar_formulasi <- paste0("lf(", "sigma", dpar_formulasi, ")")
      } else if (grepl("^sigma~1", dpar_formulasi)) { 
        dpar_formulasi <- paste0("lf(", "", dpar_formulasi, ")")
      } else {
        dpar_formulasi <- dpar_formulasi
      }
      if (grepl("lf\\(", dpar_formulasi) |
          grepl("nlf\\(", dpar_formulasi)) {
        if (grepl("^lf\\(", dpar_formulasi) &
            !grepl("nlf\\(", dpar_formulasi)) {
          lf_list <- c('flist',
                       'dpar',
                       'resp',
                       'center',
                       'cmc',
                       'sparse',
                       'decomp') #
        } else if (!grepl("^lf\\(", dpar_formulasi) &
                   grepl("^nlf\\(", dpar_formulasi)) {
          lf_list <- c('flist', 'dpar', 'resp', 'loop ') 
        }
        lf_list_c <- c()
        for (lf_listi in lf_list) {
          if (!grepl(lf_listi, dpar_formulasi)) {
            if (lf_listi == 'center') {
              o. <- paste0(lf_listi, "=", 'TRUE')
            } else if (lf_listi == 'cmc') {
              o. <- paste0(lf_listi, "=", 'TRUE')
            } else if (lf_listi == 'resp') {
              if (nys > 1) {
                o. <- paste0(lf_listi, "=", paste0("'", ysi, "'"))
                # o. <- paste0(lf_listi, "=", 'NULL')
              } else {
                o. <- paste0(lf_listi, "=", 'NULL')
              }
            } else {
              o. <- paste0(lf_listi, "=", 'NULL')
            }
            lf_list_c <- c(lf_list_c, o.)
          }
        }
        lf_list_c <- paste(lf_list_c, collapse = ",")
        if (lf_list_c != "")
          lf_list_c <- paste0(",", lf_list_c)
        dpar_formulasi <- gsub(")$", lf_list_c, dpar_formulasi)
        dpar_formulasi <- paste0(dpar_formulasi, ")")
      }
    }
    
    
    # Check for higher level model and update level 2 random formula
    f_checks_gr_gr_str <- function(a, b) {
      if(!is.null(a)) {
        gr_st_id <- sub(".*\\|", "", a) 
        a_ <- paste0("'", deparse_0(substitute(a)), "'")
        b_ <- paste0("'", deparse_0(substitute(b)), "'")
        b_out <- NULL
        if(is.null(b[[1]])) {
          if(grepl(":", gr_st_id, fixed = T) | 
             grepl("/", gr_st_id, fixed = T)) {
            stop2c("Models beyound two levels of hierarchy are not supported yet",
                 "\n ",
                 "An alternative to argument ", a_, " is to use ",
                 "\n ",
                 "argument ", b_, " to directly pass on the ", 
                 "\n ",
                 "random formula to the brms and then either accept",
                 "\n ",
                 "default priors placed by brms for those varinace covarinace", 
                 "\n ",
                 "or else use get_prios to place priors manually the pass ",
                 "\n ",
                 "by using argument 'set_self_priors'"
            )
          }
        } else if(!is.null(b[[1]])) {
          b_out <- b
        }
        out <- b_out
      } # if(!is.null(a)) {
      if(is.null(a)) {
        out <- NULL
      }
      out
    } # f_checks_gr_gr_str
    
    test_gr_sr_str_function <- function(x_grsi, x_gr_strsi) {
      if(!is.null(x_grsi)) {
        if(x_gr_strsi != 'NULL') {
          if(x_grsi == 'NULL') {
            x_grsi <- "~1"
          } else if(x_grsi != "~1") {
            if(verbose) {
              message2c("Argument '", 
                      substitute(x_grsi), "' changed from '", 
                      x_grsi , "' to  '~1'.")
              message2c("Instead of '", 
                      substitute(x_grsi), " = ", x_grsi, "', 
                    the covariates are now specified as '", 
                      substitute(x_gr_strsi), " = ", x_grsi, "'")
            }
            x_grsi <- "~1"
            
          } else {
            x_grsi <- x_grsi
          }
        } 
        out <- x_grsi
      } # if(!is.null(x_grsi)) {
      if(is.null(x_grsi)) {
        out <- NULL
      }
      return(out)
    } 
    
    # Over ride when restricting to abcd
    if(!exists('s_formula_grsi')) s_formula_grsi <- NULL
    
    a_formula_grsi <- 
      test_gr_sr_str_function(a_formula_grsi, a_formula_gr_strsi)
    b_formula_grsi <- 
      test_gr_sr_str_function(b_formula_grsi, b_formula_gr_strsi)
    c_formula_grsi <- 
      test_gr_sr_str_function(c_formula_grsi, c_formula_gr_strsi)
    d_formula_grsi <- 
      test_gr_sr_str_function(d_formula_grsi, d_formula_gr_strsi)
    e_formula_grsi <- 
      test_gr_sr_str_function(e_formula_grsi, e_formula_gr_strsi)
    f_formula_grsi <- 
      test_gr_sr_str_function(f_formula_grsi, f_formula_gr_strsi)
    g_formula_grsi <- 
      test_gr_sr_str_function(g_formula_grsi, g_formula_gr_strsi)
    h_formula_grsi <- 
      test_gr_sr_str_function(h_formula_grsi, h_formula_gr_strsi)
    i_formula_grsi <- 
      test_gr_sr_str_function(i_formula_grsi, i_formula_gr_strsi)
    s_formula_grsi <- 
      test_gr_sr_str_function(s_formula_grsi, s_formula_gr_strsi)
    
    a_fcgs_out <- f_checks_gr_gr_str(a_formula_grsi, a_formula_gr_strsi)
    b_fcgs_out <- f_checks_gr_gr_str(b_formula_grsi, b_formula_gr_strsi)
    c_fcgs_out <- f_checks_gr_gr_str(c_formula_grsi, c_formula_gr_strsi)
    d_fcgs_out <- f_checks_gr_gr_str(d_formula_grsi, d_formula_gr_strsi)
    e_fcgs_out <- f_checks_gr_gr_str(e_formula_grsi, e_formula_gr_strsi)
    f_fcgs_out <- f_checks_gr_gr_str(f_formula_grsi, f_formula_gr_strsi)
    g_fcgs_out <- f_checks_gr_gr_str(g_formula_grsi, g_formula_gr_strsi)
    h_fcgs_out <- f_checks_gr_gr_str(h_formula_grsi, h_formula_gr_strsi)
    i_fcgs_out <- f_checks_gr_gr_str(i_formula_grsi, i_formula_gr_strsi)
    s_fcgs_out <- f_checks_gr_gr_str(s_formula_grsi, s_formula_gr_strsi)
    
    sigma_formula_grsi_NULL <- sigma_formula_gr_strsi_NULL <- FALSE
    if (is.null(sigma_formula_grsi[[1]][1]) |
        sigma_formula_grsi == "NULL") {
      sigma_formula_grsi_NULL <- TRUE
    }
    
    if (is.null(sigma_formula_gr_strsi[[1]][1]) |
        sigma_formula_gr_strsi == "NULL") {
      sigma_formula_gr_strsi_NULL <- TRUE
    }
    
    if(randomsi == "") {
      if (!sigma_formula_grsi_NULL |
          !sigma_formula_gr_strsi_NULL) {
        stop2c("Random effect for parameter 'sigma' are not allowed",
             " \n ",
             " if no group level random effect is specified i.e., random = ''",
             " \n ",
             " Therefore, 
             please set argument sigma_formula_gr/sigma_formula_gr_str to NULL")
      }
    }
    
    sigma_formula_grsi <- test_gr_sr_str_function(sigma_formula_grsi, 
                                                  sigma_formula_gr_strsi)
   
    
    if(sigma_formula_gr_strsi != 'NULL') {
      if(!grepl("^~", sigma_formula_gr_strsi)) {
        sigma_formula_gr_strsi <- paste0("~", sigma_formula_gr_strsi)
      }
    }
    if(is.null(sigma_formula_gr_strsi[[1]])) {
      sigma_formula_gr_strsi <- 'NULL'
    }
    
    if(sigma_formula_grsi != 'NULL') {
      if(!grepl("^~", sigma_formula_grsi)) {
        sigma_formula_grsi <- paste0("~", sigma_formula_grsi)
      }
    }
    if(is.null(sigma_formula_grsi[[1]])) {
      sigma_formula_grsi <- 'NULL'
    }
    
    sigma_fcgs_out <- f_checks_gr_gr_str(sigma_formula_grsi, 
                                         sigma_formula_gr_strsi)
    
    if(!is.null(a_fcgs_out)) {
      if(a_formula_grsi == "~1" & !is.null(a_formula_gr_strsi[[1]])) {
        a_formula_grsi <- strsplit(a_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(b_fcgs_out)) {
      if(b_formula_grsi == "~1" & !is.null(b_formula_gr_strsi[[1]])) {
        b_formula_grsi <- strsplit(b_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(c_fcgs_out)) {
      if(c_formula_grsi == "~1" & !is.null(c_formula_gr_strsi[[1]])) {
        c_formula_grsi <- strsplit(c_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(d_fcgs_out)) {
      if(!is.null(d_formula_grsi[[1]]) & !is.null(d_formula_gr_strsi[[1]])) {
        if(d_formula_grsi == "~1" & !is.null(d_formula_gr_strsi[[1]])) {
          d_formula_grsi <- strsplit(d_formula_gr_strsi, 
                                     "+(", fixed = T)[[1]][1]
        }
      }
    }
    
    if(!is.null(e_fcgs_out)) {
      if(e_formula_grsi == "~1" & !is.null(e_formula_gr_strsi[[1]])) {
        e_formula_grsi <- strsplit(e_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(f_fcgs_out)) {
      if(f_formula_grsi == "~1" & !is.null(f_formula_gr_strsi[[1]])) {
        f_formula_grsi <- strsplit(f_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(g_fcgs_out)) {
      if(g_formula_grsi == "~1" & !is.null(g_formula_gr_strsi[[1]])) {
        g_formula_grsi <- strsplit(g_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(h_fcgs_out)) {
      if(h_formula_grsi == "~1" & !is.null(h_formula_gr_strsi[[1]])) {
        h_formula_grsi <- strsplit(h_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(i_fcgs_out)) {
      if(i_formula_grsi == "~1" & !is.null(i_formula_gr_strsi[[1]])) {
        i_formula_grsi <- strsplit(i_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(s_fcgs_out)) {
      if(s_formula_grsi == "~1" & !is.null(s_formula_gr_strsi[[1]])) {
        s_formula_grsi <- strsplit(s_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    if(!is.null(sigma_fcgs_out) & sigma_fcgs_out != 'NULL') {
      if(sigma_formula_grsi == "~1" & !is.null(sigma_formula_gr_strsi[[1]])) {
        sigma_formula_grsi <- strsplit(sigma_formula_gr_strsi, 
                                       "+(", fixed = T)[[1]][1]
      }
    }
    
    # 24.08.2024
    gsub_paranth_formula_grsi <- function(x) {
      if(!is.null(x)) {
        x_formula_grsi <- x
        strpartstrx <- strsplit(x_formula_grsi, "|", fixed = T)[[1]]
        strpartstrx_form <- strpartstrx[1]
        strpartstrx_form <-  gsub("~(", "~",  strpartstrx_form, fixed = T)
        if(length(strpartstrx) > 1 ) {
          strpartstrx_grpa <- strpartstrx[2:length(strpartstrx)]
          strpartstrx_grpa <- gsub("[()]", "", strpartstrx_grpa)
          strpartstrx_grpa2 <- paste0("", strpartstrx_grpa, collapse = "|")
          x_formula_grsi <- paste0(strpartstrx_form, "|", strpartstrx_grpa2)
        } else {
          x_formula_grsi <- strpartstrx_form
        }
      } else if(is.null(x)) {
        x_formula_grsi <- NULL
      }
      if(!is.null(x_formula_grsi)) {
        if(grepl("^\\(", x_formula_grsi)) {
          x_formula_grsi <- gsub("^\\(", "", x_formula_grsi)
        }
      }
      return(x_formula_grsi)
    }
    
    a_formula_grsi <- gsub_paranth_formula_grsi(a_formula_grsi)
    b_formula_grsi <- gsub_paranth_formula_grsi(b_formula_grsi)
    c_formula_grsi <- gsub_paranth_formula_grsi(c_formula_grsi)
    if(!is.null(d_formula_grsi))  d_formula_grsi <- 
      gsub_paranth_formula_grsi(d_formula_grsi)
    e_formula_grsi <- gsub_paranth_formula_grsi(e_formula_grsi)
    f_formula_grsi <- gsub_paranth_formula_grsi(f_formula_grsi)
    g_formula_grsi <- gsub_paranth_formula_grsi(g_formula_grsi)
    h_formula_grsi <- gsub_paranth_formula_grsi(h_formula_grsi)
    i_formula_grsi <- gsub_paranth_formula_grsi(i_formula_grsi)
    
    sigma_formula_grsi <- gsub_paranth_formula_grsi(sigma_formula_grsi)
    
    set_higher_levels <- TRUE
    if(is.null(a_fcgs_out) & 
       is.null(b_fcgs_out) & 
       is.null(c_fcgs_out) & 
       is.null(d_fcgs_out) &
       is.null(e_fcgs_out) &
       is.null(f_fcgs_out) &
       is.null(g_fcgs_out) &
       is.null(h_fcgs_out) &
       is.null(i_fcgs_out) &
       is.null(s_fcgs_out)
       ) {
      set_higher_levels <- FALSE
    }
    
    sigma_set_higher_levels <- TRUE
    if(is.null(sigma_fcgs_out) | sigma_fcgs_out == 'NULL') {
      sigma_set_higher_levels <- FALSE
    }
    
    check_formuals <- c(
        "a_formulasi",
        "b_formulasi",
        "c_formulasi",
        "d_formulasi",
        "e_formulasi",
        "f_formulasi",
        "g_formulasi",
        "h_formulasi",
        "i_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        "e_formula_grsi",
        "f_formula_grsi",
        "g_formula_grsi",
        "h_formula_grsi",
        "i_formula_grsi",
        "s_formula_grsi",
        "sigma_formulasi",
        "sigma_formula_grsi")
    
    check_formuals_t_f <- c(check_formuals, 
                            "dpar_formulasi",
                            'sigma_formula_gr_strsi')
    
    # new
    for (check_formualsi in check_formuals_t_f) {
      if(!is.null(ept(check_formualsi)) & length(ept(check_formualsi)) !=0 ) {
        assign(check_formualsi, replace_t_f_to_full(ept(check_formualsi)))
      }
    }
    
    for (check_formualsi in check_formuals) {
      if(!is.null(ept(check_formualsi)) & length(ept(check_formualsi)) !=0 ) {
        if (!grepl("~1", ept(check_formualsi)) &
            !grepl("~0", ept(check_formualsi))) {
          check_formualsi_with1 <-
            gsub("^~", "~1+", ept(check_formualsi), fixed = F)
          if(!grepl("^~", ept(check_formualsi))) {
            if(!grepl("^sigma", check_formualsi))
              check_formualsi_with1 <- paste0("~", check_formualsi_with1)
          }
          assign(check_formualsi, check_formualsi_with1)
        }
      } # if(!is.null(ept(check_formualsi))) {
      if(is.null(ept(check_formualsi)) | length(ept(check_formualsi)) ==0 ) {
        assign(check_formualsi, NULL)
      }
    } 
    
    if (is.null(sigma_formula_gr_strsi[[1]][1]) |
        sigma_formula_gr_strsi == "NULL") {
      sigma_formula_gr_strsi <- NULL
    }
    
    if (is.null(sigma_formula_grsi[[1]][1]) |
        sigma_formula_grsi == "NULL") {
      sigma_formula_grsi <- NULL
    }
    
    if (is.null(dpar_formulasi[[1]][1]) |
        dpar_formulasi == "NULL") {
      dpar_formulasi <- NULL
    }
    
    if (is.null(autocor_formulasi[[1]][1]) |
        autocor_formulasi == "NULL") {
      autocor_formi <- NULL
    } else {
      autocor_formulasi <- gsub("\"", "", autocor_formulasi)
      if(!grepl("^~", autocor_formulasi)) {
        stop2c('autocor_formula argument should be a formula. E.g.,',
             "\n ",
             " autocor_formula = ~arms(p=1,q=1)",
             "\n ", 
             " It seems you forgot to add '~' before the autocor structure")
      }
      autocor_formi <- autocor_formulasi
    } 
    
    if(!is.null(autocor_formi)) {
      tempunstx <- autocor_formi # '~unstr(time=visit, patient)'
      tempunstx <- gsub("[[:space:]]", "", tempunstx)
      if(grepl("unstr(", tempunstx, fixed = T)) {
        tempunstx_1 <- regmatches(tempunstx, gregexpr("(?<=\\().*?(?=\\))", 
                                                      tempunstx, perl=T))[[1]]
        tempunstx_2 <- strsplit(tempunstx_1, ",")[[1]][1]
        if(grepl("time=", tempunstx_2, fixed = T)) {
          tempunstx_3 <- sub(".*time=", "", tempunstx_2) 
        } else if(!grepl("time=", tempunstx_2, fixed = T)) {
          tempunstx_3 <- tempunstx_2
        }
        cortimeNlags_var <- tempunstx_3
      } # if(grepl("unstr(", tempunstx, fixed = T)) {
      
      if(!grepl("unstr(", tempunstx, fixed = T)) {
        cortimeNlags_var <- NULL
      }
    } 
      
    if(is.null(autocor_formi)) {
      cortimeNlags_var <- NULL
    }
   
    
    if (is.null(familysi[[1]][1]) |
        familysi == "NULL") {
      familysi <- NULL
    }
    
   # For backward compatibility if model fit using family = gaussian()
   if (!is.null(familysi)) {
     if(familysi == "gaussian()") {
       familysi <- "brms::brmsfamily(family = gaussian)"
     } else {
       if(!grepl('family(', familysi, fixed = T)) {
         # stop2c("The 'family' argument must be specified by explicitly using the",
         #      "\n  ",
         #      "family(....) or brms::brmsfamily(....) form",
         #      "\n  ",
         #      "For example, to specify the 'gaussian()' family, please use:",
         #      "\n  ",
         #      "brmsfamily(family='gaussian', link='identity', link_sigma='log')"
         #      )
       }
     }
   }
   
    if (!is.null(familysi)) {
      familysi_check <- familysi
      if(grepl('brmsfamily', familysi_check) &
          grepl('brms::', familysi_check)) {
        familysi_check <- familysi_check
      } else if( grepl('brmsfamily', familysi_check) &
                 !grepl('brms::', familysi_check)) {
        familysi_check <- paste0('brms::', familysi_check)
      } else if( grepl('family', familysi_check) &
                 !grepl('brms::', familysi_check)) {
        # familysi_check <- paste0('brms::brmsfamily', familysi_check)
        familysi_check <- gsub('family', 'brms::brmsfamily', 
                               familysi_check, fixed = T)
      } else if( grepl('\\(', familysi_check) &
                 !grepl('family', familysi_check)) {
        familysi_check_w <- strsplit(familysi_check, "[^a-zA-Z]+")[[1]]
        familysi_check_w <- collapse_comma(familysi_check_w)
        familysi_check_w <- paste(familysi_check_w, collapse = ",")
        familysi_check_w <- paste0('brms::brmsfamily', "(", familysi_check_w, ")")
        familysi_check <- familysi_check_w
      } else if(!grepl('brmsfamily', familysi_check) & 
                !grepl('family', familysi_check)) {
         # stop2c("Argument family should be specified as brmsfamily(family,...)",
         #      "\n ", 
         #      "where family is the name of family such as 'gaussian' and",
         #      "\n ", 
         #      "... are the family specific argument such as link and link_sigma",
         #      "\n ", 
         #      "For example, 'gaussian' family can be set explicitly as follows:",
         #      "\n ", 
         #      "brmsfamily('gaussian', link = 'identity', link_sigma = 'log')"
         #      )
        familysi_check <- paste0("brms::brmsfamily(", "'", 
                                 familysi_check, 
                                 "'", ")")
      }
      familysi <- familysi_check
    } # if (!is.null(familysi)) {
  
   
    
    if (!is.null(familysi)) {
      familysi_temp <- list_to_quoted_if_not_si(familysi)
      if(familysi_temp == "NA" | is.na(familysi_temp)) {
        familysi <- familysi
      } else {
        familysi <- familysi_temp
      }
    }
    
    familysi <- gsub_space(familysi)
    
    if (!is.null(dpar_formulasi)) {
      if (grepl("^lf\\(", dpar_formulasi) |
          grepl("^nlf\\(", dpar_formulasi)) {
      } else {
        dpar_formulasi <- dpar_formulasi
      }
    }
    

    N_J_all <- length(unique(data[[idsi]]))
  
    ##########################
    # add_sigma_by_ls
    setsigmaxvarsi <- FALSE
    if(set_model_sigma_by_ls) {
      setsigmaxvarsi <- TRUE
    } 
  
    
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
      datai <- data %>%
        dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>%
        droplevels()
      if (!subindicatorsi %in% colnames(datai)) {
        stop2c("variable ", subindicatorsi, " not in the dataframe")
      }
      if (!xsi %in% colnames(datai))
        stop2c("variable ", xsi, " not in the dataframe")
      if (!idsi %in% colnames(datai))
        stop2c("variable ", idsi, " not in the dataframe")
    } # if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    
    if ((is.na(univariate_by$by) | univariate_by$by == "NA")) {
      datai <- data %>% droplevels()
      if (!ysi %in% colnames(datai))
        stop2c("variable ", ysi, " not in the dataframe")
      if (!xsi %in% colnames(datai))
        stop2c("variable ", xsi, " not in the dataframe")
      if (!idsi %in% colnames(datai))
        stop2c("variable ", idsi, " not in the dataframe")
    }
   
    # 28 01 2024
    drop_na_vars <- c(xsi, ysi, idsi)
    datai <- datai %>% tidyr::drop_na(., dplyr::any_of(drop_na_vars))
    
    check_variable_numeric_exists(datai, c(xsi, ysi))
   
    if(!is.null(cortimeNlags_var)) {
      if(!is.factor(datai[[cortimeNlags_var]])) {
        datai[[cortimeNlags_var]] <- as.factor(datai[[cortimeNlags_var]])
        datai[[cortimeNlags_var]] <- droplevels(datai[[cortimeNlags_var]])
      } else {
        datai[[cortimeNlags_var]] <-  datai[[cortimeNlags_var]]
      }
      cortimeNlags <- nlevels(unique(datai[[tempunstx_3]]))
    } else if(is.null(cortimeNlags_var)) {
      cortimeNlags <- NULL
    }
    

    if(is.null(parameterization)) {
      checkoccs <- datai %>% 
        dplyr::filter(!is.na(ysi)) %>% 
        droplevels() %>% 
        dplyr::mutate(nid = dplyr::n_distinct(idsi)) %>%
        dplyr::group_by_at(idsi) %>% 
        dplyr::mutate(NoccPI = dplyr::row_number()) %>% 
        dplyr::mutate(NoccAI = max(NoccPI)) %>% 
        dplyr::ungroup()
      
      if(min(checkoccs$NoccAI) >= 10) {
        parameterization = 'cp'
      } else {
        parameterization = 'ncp'
      }
    }
    
    fit_edited_scode <- FALSE
    if(select_model == 'logistic1e' |
       select_model == 'logistic2e' |
       select_model == 'logistic3e' |
       parameterization == 'cp') {
      fit_edited_scode <- TRUE
    }
    if(sum_zero) {
      fit_edited_scode <- TRUE
    }
    
    # add_sigma_by_mu
    if(set_model_sigma_by_fz |
       set_model_sigma_by_fp |
       set_model_sigma_by_fe | 
       set_model_sigma_by_mp | 
       set_model_sigma_by_me | 
       set_model_sigma_by_rp |
       set_model_sigma_by_re) {
      fit_edited_scode <- TRUE
    }
    
    # add_rescor_by
    set_rescor_by <- FALSE
    if (nys > 1) {
      if(multivariate$mvar) {
        if(multivariate$rescor) {
          if(!is.null(multivariate$rcorr_by)) fit_edited_scode <- TRUE
          if(!is.null(multivariate$rcorr_by)) set_rescor_by    <- TRUE
        } # if(!is.null(multivariate$rcorr_by)) {
      } # if(multivariate$mvar) {
    } # if (nys > 1) {
 
    ######################################################################
    ######################################################################
    
    # Refactor to use function() for transformations
    # This will allow using optimize_x = list(function(x) log(x + 3/4))
    # Note that instead of calling log(data[[xsi]]), 'xfuntransformsi' be used 
    # Check if xfunsi, yfunsi and sigmaxfunsi
    # Should sigmayfunsi be the sigma_link family
    
    familysi_ept      <- ept(familysi)
    family_link_sigma <- familysi_ept[['link_sigma']]
    
    remove_sigma_parameter <- FALSE
    if(is.null(family_link_sigma)) {
      remove_sigma_parameter <- TRUE
    }
    
    # If there is no sigma parameter such as family gamma
    if(!remove_sigma_parameter) {
      # Leave it as identity
      if(family_link_sigma == "log") {
        sigmayfunsi <- NULL # 'log'
      } else {
        sigmayfunsi <- NULL
      }
    }
    
    set_xfunsi      <- check_if_arg_set(xfunsi)
    set_yfunsi      <- check_if_arg_set(yfunsi)
    set_sigmaxfunsi <- check_if_arg_set(sigmaxfunsi)
    set_sigmayfunsi <- check_if_arg_set(sigmayfunsi)
   
    if (!set_xfunsi) {
      xfuntransformsi <- function(x)x
      assign('xfuntransformsi', xfuntransformsi, envir = enverr.)
    } else if (set_xfunsi) {
      if(xfunsi == "log") {
        xfuntransformsi <- function(x)log(x)
      } else if(xfunsi == "sqrt") {
        xfuntransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(xfunsi))) {
        xfuntransformsi <- ept(xfunsi)
      } else {
        stop2c(paste0(
          "The xfun argument must be either a string ('log' or 'sqrt'),", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('xfuntransformsi', xfuntransformsi, envir = enverr.)
    }
    
    if (!set_yfunsi) {
      yfuntransformsi <- function(x)x
      assign('yfuntransformsi', yfuntransformsi, envir = enverr.)
    } else if (set_yfunsi) {
      if(yfunsi == "log") {
        yfuntransformsi <- function(x)log(x)
      } else if(yfunsi == "sqrt") {
        yfuntransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(yfunsi))) {
        yfuntransformsi <- ept(yfunsi)
      } else {
        stop2c(paste0("The yfun argument must be a string ('log' or 'sqrt'),", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('yfuntransformsi', yfuntransformsi, envir = enverr.)
    }
    
    
    if (!set_sigmaxfunsi) {
      sigmaxfuntransformsi <- function(x)x
      assign('sigmaxfuntransformsi', sigmaxfuntransformsi, envir = enverr.)
    } else if (set_sigmaxfunsi) {
      if(sigmaxfunsi == "log") {
        sigmaxfuntransformsi <- function(x)log(x)
      } else if(sigmaxfunsi == "sqrt") {
        sigmaxfuntransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(sigmaxfunsi))) {
        sigmaxfuntransformsi <- ept(sigmaxfunsi)
      } else {
        stop2c(paste0("The xfun argument must be either 'log' or 'sqrt',", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('sigmaxfuntransformsi', sigmaxfuntransformsi, envir = enverr.)
    }
    
    if (!set_sigmayfunsi) {
      sigmayfuntransformsi <- function(x)x
      assign('sigmayfuntransformsi', sigmayfuntransformsi, envir = enverr.)
    } else if (set_sigmayfunsi) {
      if(sigmayfunsi == "log") {
        sigmayfuntransformsi <- function(x)log(x)
      } else if(sigmayfunsi == "sqrt") {
        sigmayfuntransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(sigmayfunsi))) {
        sigmayfuntransformsi <- ept(sigmayfunsi)
      } else {
        stop2c(paste0("The xfun argument must be either 'log' or 'sqrt',", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('sigmayfuntransformsi', sigmayfuntransformsi, envir = enverr.)
    }
    
    
    if (!is.null(sigmaxoffsetsi[[1]][1]) & sigmaxoffsetsi != "NULL") {
      sigmaxoffsetsi <- sigmaxoffsetsi
    } else {
      sigmaxoffsetsi <- xoffsetsi
      if(verbose) message2c("xoffset for sigma is set same as for mu")
    }
    
    # add_sigma_by_ls - sigmaxoffsetsi only relevant for location scale model
    if(!set_model_sigma_by_ls) {
      sigmaxoffsetsi <- 0
    }
    
    # Check if xfunsi, yfunsi and sigmaxfunsi
    # These are offset specific funs - if NULL, then default xfun will be 
    # applied to the offset - same for sigmaoffset
    set_xfunxoffsetsi      <- check_if_arg_set(xfunxoffsetsi)
    set_sigmaxfunxoffsetsi <- check_if_arg_set(sigmaxfunxoffsetsi)
    
    
    if (!set_xfunxoffsetsi) {
      xfunxoffsettransformsi <- function(x)x
      assign('xfunxoffsettransformsi', xfunxoffsettransformsi, envir = enverr.)
    } else if(set_xfunxoffsetsi) {
      if(xfunxoffsetsi == "T" | xfunxoffsetsi == "TRUE") {
        if(verbose) {
          message2c("Argument 'xfunxoffset' set same as 'xfun'")
        }
        xfunxoffsettransformsi <- xfuntransformsi
      } else if(xfunxoffsetsi == "log") {
        xfunxoffsettransformsi <- function(x)log(x)
      } else if(xfunxoffsetsi == "sqrt") {
        xfunxoffsettransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(xfunxoffsetsi))) {
        xfunxoffsettransformsi <- ept(xfunxoffsetsi)
      } else {
        stop2c(paste0(
          "The xfunxoffset argument must be either a string ('log' or 'sqrt'),", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('xfunxoffsettransformsi', xfunxoffsettransformsi, envir = enverr.)
    }
    
    if (!set_sigmaxfunxoffsetsi) {
      sigmaxfunxoffsettransformsi <- function(x)x
      assign('sigmaxfunxoffsettransformsi', sigmaxfunxoffsettransformsi, 
             envir = enverr.)
    } else if (set_sigmaxfunxoffsetsi) {
      if(sigmaxfunxoffsetsi == "T" | sigmaxfunxoffsetsi == "TRUE") {
        if(verbose) {
          message2c("Argument 'sigmaxfunxoffset' set same as 'sigmaxfun'")
        }
        sigmaxfunxoffsettransformsi <- sigmaxfuntransformsi
      } else if(sigmaxfunxoffsetsi == "log") {
        sigmaxfunxoffsettransformsi <- function(x)log(x)
      } else if(sigmaxfunxoffsetsi == "sqrt") {
        sigmaxfunxoffsettransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(sigmaxfunxoffsetsi))) {
        sigmaxfunxoffsettransformsi <- ept(sigmaxfunxoffsetsi)
      } else {
        stop2c(paste0("The 'sigmaxfunxoffset' must be a
                      string ('log' / 'sqrt'),", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('sigmaxfunxoffsettransformsi', sigmaxfunxoffsettransformsi, 
             envir = enverr.)
    }
    
    #############################################################
    ############################################################
    # check functions 
    funs_c <- c(xfuntransformsi, 
                yfuntransformsi,
                sigmaxfuntransformsi,
                sigmayfuntransformsi,
                xfunxoffsettransformsi,
                sigmaxfunxoffsettransformsi)
    
    names(funs_c) <- c("xfuntransformsi",
                       "yfuntransformsi",
                       "sigmaxfuntransformsi",
                       "sigmayfuntransformsi",
                       "xfunxoffsettransformsi",
                       "sigmaxfunxoffsettransformsi")
    
    for (i in 1:length(funs_c)) {
      set_fun <- funs_c[[i]]
      set_fun_name <- names(funs_c)[i]
      assign(set_fun_name, check_and_rename_funs_args_to_x(set_fun, 
                                                           checkname = 'x') )
    }
   
    ############################################################
    ############################################################
    
    # Assign reverse functions also
    assign("ixfuntransformsi",  
           inverse_transform(base::body(xfuntransformsi)), 
           envir = enverr.)
    
    assign("iyfuntransformsi",  
           inverse_transform(base::body(yfuntransformsi)), 
           envir = enverr.)
    
    assign("sigmaixfuntransformsi",  
           inverse_transform(base::body(sigmaxfuntransformsi)), 
           envir = enverr.)
    
    assign("sigmaiyfuntransformsi",  
           inverse_transform(base::body(sigmayfuntransformsi)), 
           envir = enverr.)
    
    
    if(!set_model_sigma_by_ls) {
      if(!is.na(sigmaxsi)) {
        datai[[sigmaxsi]] <- NULL
        sigmaxsi          <- NA
        sigmaxs[ii]       <- NA
      }
    }
    
    
    #################################################################
    #################################################################
    check_for_validy_of_prepare_transformations    <- TRUE
    prepare_transformations_args                   <- list()
    prepare_transformations_args[['data']]         <- datai
    prepare_transformations_args[['xvar']]         <- xsi
    prepare_transformations_args[['yvar']]         <- ysi
    prepare_transformations_args[['sigmaxvar']]    <- sigmaxsi
    prepare_transformations_args[['xfun']]         <- xfuntransformsi
    prepare_transformations_args[['yfun']]         <- yfuntransformsi
    prepare_transformations_args[['sigmaxfun']]    <- sigmaxfuntransformsi
    prepare_transformations_args[['ixfun']]        <- FALSE
    prepare_transformations_args[['iyfun']]        <- FALSE
    prepare_transformations_args[['sigmaixfun']]   <- FALSE
    prepare_transformations_args[['xoffset']]      <- NULL
    prepare_transformations_args[['sigmaxoffset']] <- NULL
    prepare_transformations_args[['model']]        <- NULL
    prepare_transformations_args[['envir']]        <- NULL
    prepare_transformations_args[['verbose']]      <- FALSE
    prepare_transformations_args[['transform']]    <- ""
    prepare_transformations_args[['itransform']]   <- ""
    
    # Here need to transform x, y, and sigmax 
    # x y must because accordingly xoffset and lm matrix etc created
    if(check_for_validy_of_prepare_transformations) {
      check_for_validy_of_prepare_transformations_0 <- datai
    }
    
    datai <- CustomDoCall(prepare_transformations, prepare_transformations_args)
    
    if(check_for_validy_of_prepare_transformations) {
      check_for_validy_of_prepare_transformations_1 <- datai
    }
    
   
    if (is.numeric(ept(knotssi))) {
      knots <- ept(knotssi)
    }
    
    if (is.numeric(ept(dfsi))) {
      knots <- unname(gkn(datai[[xsi]], ept(dfsi), ept(boundsi)))
      if(verbose) {
        message2c("For '", smat, "' knots are created internally based on 'df'",
                "\n ",
                " Note that knots are constructed and then adjusted by xoffset",
                "\n ",
                " such as knots - xoffset") 
      }
    }
    
    knots_from_gkn <- knots
    
    #################################################################
    #################################################################
    
    # Over ride smat_bkrange with knots_selection[['bkrange']]
    if(!is.null(mcall[['knots_selection']])) {
      knots_selection_bkrange    <- mcall[['knots_selection']][['bkrange']]
      knots_selection_fix_bknots <- mcall[['knots_selection']][['fix_bknots']]
      knots_selection_what       <- mcall[['knots_selection']][['what']]
      knots_selection_when       <- mcall[['knots_selection']][['when']]
      knots_selection_method     <- mcall[['knots_selection']][['method']]
      knots_selection_return     <- mcall[['knots_selection']][['return']]
      knots_selection_print      <- mcall[['knots_selection']][['print']]
      knots_selection_bkrange    <- knots_selection_bkrange %>% as.logical()
      knots_selection_fix_bknots <- knots_selection_fix_bknots %>% as.logical()
      knots_selection_return     <- knots_selection_return %>% as.logical()
      knots_selection_print      <- knots_selection_print %>% as.logical()
      knots_selection_what       <- knots_selection_what
      knots_selection_when       <- knots_selection_when
      knots_selection_method     <- knots_selection_method
      
      smat_bkrange               <- knots_selection_bkrange
      smat_fix_bknots            <- knots_selection_fix_bknots
      smat_what                  <- knots_selection_what
      smat_when                  <- knots_selection_when
      smat_return                <- knots_selection_return
      smat_print                 <- knots_selection_print
      if(knots_selection_bkrange) {
        set_bknots <- NULL
      } else {
        set_bknots <- checkgetiknotsbknots(knots, 'bknots')
      }
    } else if(is.null(mcall[['knots_selection']])) {
      set_bknots <- checkgetiknotsbknots(knots, 'bknots')
      if(as.logical(smat_bkrange)) set_bknots <- NULL
      knots_selection_fix_bknots <- TRUE
      if(is_emptyx(smat_what)) smat_what <- 'plot1'
      if(is_emptyx(smat_when)) smat_when <- 'bc'
    }

    # get_knost_from_df_arg used later for plot, so keep it out here
    get_knost_from_df_arg <- list()
    get_knost_from_df_arg[['x']]         <-  datai[[xsi]]
    get_knost_from_df_arg[['knots']]     <-  NULL
    get_knost_from_df_arg[['bknots']]    <-  set_bknots
    get_knost_from_df_arg[['df']]        <-  ept(dfsi)
    get_knost_from_df_arg[['degree']]    <-  smat_degree
    get_knost_from_df_arg[['intercept']] <-  smat_intercept
    get_knost_from_df_arg[['derivs']]    <-  smat_derivs
    get_knost_from_df_arg[['centerval']] <-  smat_centerval
    get_knost_from_df_arg[['normalize']] <-  smat_normalize
    get_knost_from_df_arg[['preH']]      <-  smat_preH
    get_knost_from_df_arg[['sfirst']]    <-  smat_sfirst
    get_knost_from_df_arg[['sparse']]    <-  smat_sparse
    get_knost_from_df_arg[['bound']]     <-  ept(boundsi)
    get_knost_from_df_arg[['xoffset']]   <-  NULL # don't set xoffset here
    get_knost_from_df_arg[['bkrange']]   <-  smat_bkrange
    get_knost_from_df_arg[['fix_bknots']]<-  smat_fix_bknots
    get_knost_from_df_arg[['smat']]      <-  smat
    
    
    if(knotssi == "NA" | is.na(knotssi)) {
      knots <- do.call(get_knost_from_df, get_knost_from_df_arg)
      if(verbose) {
        message2c("For '",smat,"' knots are created internally based on the 'df'",
                "\n ",
                " The boundary knots are adjusted for bounds. The full knots ",
                "\n ",
                " i.e, internal as well as the boundary knots are then ",
                "\n ",
                " adjusted for the xoffset i.e., knots - xoffset")
      }
    } # if(knotssi == "NA" | is.na(knotssi)) {
    
   
    
    
    #################################################################
    #################################################################
    
    knots_from_gkn      <- knots_from_gkn
    knots_from_new_funs <- knots
    
    knots_maxdp          <- max(get_decimal_places(knots_from_gkn))
    knots_from_gkn       <- round(knots_from_gkn, knots_maxdp)
    knots_from_new_funs  <- round(knots_from_new_funs, knots_maxdp)

    knots_from_new_funs_msg <- 
      paste0("knots changed from earlier version: ", 
             "\n  ",
             "Earlier: ",
             paste(deparse(knots_from_gkn), collapse = ", "), 
             "\n  ",
             "Now:     ",
             paste(deparse(knots_from_new_funs), collapse = ", "))
    
    if(is.null(mcall[['knots_selection']])) {
      if(!identical(knots_from_new_funs, knots_from_gkn)) {
        if(smat == 'rcs') {
          if(verbose) message2c(knots_from_new_funs_msg)
        } else if(smat == 'nsp' | smat == 'nsk') {
          stop2c(knots_from_new_funs_msg)
        } else {
          # stop2c(knots_from_new_funs_msg)
        }
      } # if(!identical(knots_from_new_funs, knots_from_gkn)) {
    } # if(is.null(mcall[['knots_selection']])) {
    
    
    if(!is.null(mcall[['knots_selection']])) {
      if(!identical(knots_from_new_funs, knots_from_gkn)) {
        if(verbose) message2c(knots_from_new_funs_msg)
      } # if(!identical(knots_from_new_funs, knots_from_gkn)) {
    } # if(!is.null(mcall[['knots_selection']])) {
    
    #################################################################
    #################################################################
    
    knots_selection      <- mcall[['knots_selection']] %>% eval()
    
    if(!is.null(knots_selection)) {
      if(is.null(knots_selection[['nsearch']])) {
        knots_selection[['nsearch']] <- length(knots) - 2 + 4
      } else {
        knots_selection[['nsearch']] <- knots_selection[['nsearch']]
      }
    } # if(!is.null(knots_selection)) {
    
   
    if(!is.null(knots_selection)) {
      knots_selection_arg <- get_knost_from_df_arg
      knots_selection_arg[['dataset']] <- datai
      knots_selection_arg[['dependent']] <- ysi
      knots_selection_arg[['independents']] <- xsi
      knots_selection_arg[['target_nknots']] <- length(knots) - 2
      
      knots_selection_arg[['max_nknots']] <- length(knots) - 2
      
      knots_selection_arg[['print']] <- knots_selection[['print']]
      knots_selection_arg[['return']] <- knots_selection[['return']]
      knots_selection_arg[['select']] <- knots_selection[['select']]
      knots_selection_arg[['initial_nknots']] <- knots_selection[['nsearch']]
      knots_selection_arg[['cost_fn']]<- knots_selection[['criteria']]
      knots_selection_arg[['icr_fn']] <- knots_selection[['criteria']]
      knots_selection_arg[['all_scores']] <- knots_selection[['all_scores']]
      knots_selection_arg[['all_knots']]  <- FALSE
      knots_selection_arg[['method']]     <- knots_selection_method
      knots_selection_arg[['cvk']]     <- knots_selection[['cvk']]
      knots_selection_arg[['cviter']]  <- knots_selection[['cviter']]
      knots_selection_arg[['kspace']]  <- knots_selection[['kspace']]
      knots_selection_arg[['plot_all_scores']] <- 
        knots_selection[['plot_all_scores']]
      
      knots_selection_arg[['smat']]      <- smat
      knots_selection_arg[['knots']]     <-  NULL
      knots_selection_arg[['bknots']]    <- checkgetiknotsbknots(knots,'bknots')
      knots_selection_arg[['bkrange']]   <- knots_selection_bkrange
      knots_selection_arg[['fix_bknots']]<- knots_selection_fix_bknots
      knots_selection_arg[['df']]        <-  df
      knots_selection_arg[['degree']]    <-  smat_degree
      knots_selection_arg[['intercept']] <-  smat_intercept
      knots_selection_arg[['derivs']]    <-  smat_derivs
      knots_selection_arg[['centerval']] <-  smat_centerval
      knots_selection_arg[['normalize']] <-  smat_normalize
      knots_selection_arg[['preH']]      <-  smat_preH
      knots_selection_arg[['sfirst']]    <-  smat_sfirst
      knots_selection_arg[['sparse']]    <-  smat_sparse
      knots_selection_arg[['verbose']]   <-  verbose
      # Some check for knots_selection_arg
      if(knots_selection_arg[['select']] == 'df') {
        if(!knots_selection_arg[['return']]) {
          stop2c("please use return = TRUE when select = 'df'")
        }
        knots_selection_arg[['return_df']] <- TRUE
      } else {
        knots_selection_arg[['return_df']] <- FALSE
      }
      
      if(knots_selection_arg[['select']] == 'df') {
        knots_selection_arg[['search_df']] <- TRUE
      } else if(knots_selection_arg[['select']] == 'both') {
        knots_selection_arg[['search_df']] <- TRUE
      } else if(knots_selection_arg[['select']] == 'knots') {
        knots_selection_arg[['search_df']] <- FALSE
      } else {
        stop2c("For argument 'knots_selection', the option 'select' must be 
             either 'df', 'knots', or 'both'. The current choice ", 
             collapse_comma(knots_selection_arg[['select']]),
             " is invalid")
      }
      # End of Some check for knots_selection_arg
      
    } # iif(!is.null(knots_selection)) {
   
    

    if(!is.null(knots_selection)) {
      if(knots_selection[['when']] == 'bc') {
        knots_old <- knots
        model_old <- NULL
        model_new <- do.call(get_choose_model, knots_selection_arg)
        if(knots_selection_arg[['return_df']]) {
          return(model_new)
        }
        if(knots_selection_arg[['kspace']] == "un") {
          knots_new <- model_new$ns$knot_placements$fullknots
        } 
        if(knots_selection_arg[['kspace']] == "nu") {
          knots_new <- model_new$ns_nu$knot_placements$fullknots
        }
        getmaxdp  <- max(get_decimal_places(knots_old))
        knots_new <- round(knots_new, getmaxdp)
        knots     <- knots_new # These knots will be passed to bsitar
        #######
        plot_object <- get_print_return_obj(knots = knots_old, 
                                            model = model_old,
                                            model_new = model_new, 
                                            knots_new = knots_new, 
                                            nys = nys,
                                            xsi = xsi, 
                                            ysi = ysi,
                                            select = 
                                              knots_selection_arg[['select']],
                                            kspace = 
                                              knots_selection_arg[['kspace']],
                                            what = knots_selection_what,
                                            print = knots_selection_print, 
                                            list_arg = knots_selection_arg,
                                            return = knots_selection_return,
                                            verbose = verbose)
        if(knots_selection_print) {
          print(plot_object)
        }
        if(knots_selection_return) {
          if(knots_selection_what == 'knots') {
            retout <- list(knots_old = knots_old,
                           knots_new = knots_new)
            return(retout)
          } 
          return(plot_object)
        }
        #######
      } # if(knots_selection[['when']] == 'bc') {
    } # if(!is.null(knots_selection)) {
    

    

    # Although stype has when option, somehow it does not work
    # Therefore, let the default 'bc' keep working 
    if(is.null(knots_selection)) {
      if(smat_when == "bc") {
        smat_knots_plot_arg <- get_knost_from_df_arg
        smat_knots_plot_arg[['dataset']] <- datai
        smat_knots_plot_arg[['dependent']] <- ysi
        smat_knots_plot_arg[['independents']] <- xsi
        plot_object <- get_print_return_obj(knots = knots, 
                                            model = NULL,
                                            model_new = NULL, 
                                            knots_new = NULL, 
                                            nys = nys,
                                            xsi = xsi, 
                                            ysi = ysi,
                                            select = 
                                              knots_selection_arg[['select']],
                                            kspace = 
                                              knots_selection_arg[['kspace']],
                                            what = smat_what,
                                            print = smat_print, 
                                            list_arg = smat_knots_plot_arg,
                                            return = smat_return,
                                            verbose = verbose)
        
        if(smat_print) {
          print(plot_object)
        }
        if(smat_return) {
          return(plot_object)
        }
      } # if(smat_when == "bc") {
    } # if(is.null(knots_selection)) {
    

    
    #################################################################
    #################################################################
    
    knots_from_gkn       <- knots_from_gkn
    knots_from_new_funs  <- knots_from_new_funs
    knots_from_selection <- knots # knots are from latest gkn -funs - selection
    
    knots_maxdp          <- max(get_decimal_places(knots_from_gkn))
    knots_from_gkn       <- round(knots_from_gkn, knots_maxdp)
    knots_from_new_funs  <- round(knots_from_new_funs, knots_maxdp)
    knots_from_selection <- round(knots_from_selection, knots_maxdp)
    
    knots_from_new_funs_selection_msg <- 
      paste0("knots summary from different versions: ", 
             "\n  ",
             "From gkn:       ",
             paste(deparse(knots_from_gkn), collapse = ", "), 
             "\n  ",
             "From df:        ",
             paste(deparse(knots_from_new_funs), collapse = ", "), 
             "\n  ",
             "From selection: ",
             paste(deparse(knots_from_selection), collapse = ", "))
    
    if(!is.null(mcall[['knots_selection']])) {
        if(verbose) message2c(knots_from_new_funs_selection_msg)
    } # if(!is.null(mcall[['knots_selection']])) {
    
   
    
    #################################################################
    #################################################################
    
    if(smat == 'rcs') {
      if(length(knots) <= 2) {
        stop2c("For '",smat,"' the minimum number of knots should be '3'")
      }
    } else {
      if(length(knots) <= 1) {
        # stop2c("For '",smat,"' the minimum number of knots should be '2'")
      }
    }
    
    #################################################################
    #################################################################
 
    if(select_model == "sitar") {
      if (match_sitar_d_form) {
        if (length(knots) > 2) {
          itemp <- strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
          itemp <- itemp[!grepl("d", itemp)]
          fixedsi <- paste(itemp, collapse = "+")
        }
      }
    }
    
    nabci   <- length(strsplit(gsub("\\+", " ", fixedsi), " ")[[1]])
    nabcrei <- length(strsplit(gsub("\\+", " ", randomsi), " ")[[1]])
   
    ##################################################################
    ##################################################################
    
    if (xoffsetsi == 'NA' | xoffsetsi == '') {
      if(grepl('sitar', select_model) | grepl('rcs', select_model)) {
        xoffsetsi <- 'mean'
      } else {
        xoffsetsi <- 0
      }
    }
    
    if(bstartsi == 'xoffset') {
      bstartsi <- xoffsetsi
    }
    
    xoffset <- eval_xoffset_bstart_args(x = xsi, 
                                        y = ysi, 
                                        knots = knots, 
                                        data = datai, 
                                        eval_arg = xoffsetsi, 
                                        offsetfunsi = xfunxoffsettransformsi, 
                                        smat = smat,
                                        degree = smat_degree,
                                        intercept = smat_intercept, 
                                        derivs = smat_derivs,
                                        centerval = smat_centerval,
                                        normalize = smat_normalize,
                                        preH = smat_preH,
                                        sfirst = smat_sfirst, 
                                        sparse = smat_sparse,
                                        arg = 'xoffset',
                                        dpar = "mu",
                                        verbose = verbose)
    
    
    bstart <- eval_xoffset_bstart_args(x = xsi, 
                                       y = ysi, 
                                       knots = knots, 
                                       data = datai, 
                                       eval_arg = bstartsi, 
                                       offsetfunsi = xfunxoffsettransformsi, 
                                       smat = smat,
                                       degree = smat_degree,
                                       intercept = smat_intercept, 
                                       derivs = smat_derivs,
                                       centerval = smat_centerval,
                                       normalize = smat_normalize,
                                       preH = smat_preH,
                                       sfirst = smat_sfirst, 
                                       sparse = smat_sparse,
                                       arg = 'bstart',
                                       dpar = "mu",
                                       verbose = verbose)
    
    bstart <- bstart - xoffset
    
    cstart <- eval_xoffset_cstart_args(x = xsi, 
                                       y = ysi, 
                                       knots = knots, 
                                       data = datai, 
                                       eval_arg = cstartsi, 
                                       smat = smat,
                                       degree = smat_degree,
                                       intercept = smat_intercept, 
                                       derivs = smat_derivs,
                                       centerval = smat_centerval,
                                       normalize = smat_normalize,
                                       preH = smat_preH,
                                       sfirst = smat_sfirst, 
                                       sparse = smat_sparse,
                                       arg = 'cstart',
                                       dpar = "mu",
                                       verbose = verbose)
    
    
    xoffset      <- round(xoffset, 8)
    
    if(is_emptyx(xoffset)) {
      xoffset      <- 0
    }
    
    if(smat == 'bsp' |  smat == 'msp' |  smat == 'isp') {
      if(!is.null(knots)) {
        knots        <- knots - xoffset
        knots        <- round(knots, 8)
        nknots       <- length(knots)
        df           <- length(knots) - 1
      } else if(is.null(knots)) {
        knots        <- NULL
        knots        <- NULL
        nknots       <- NULL
        df           <- NULL
      }
    } else {
      knots        <- knots - xoffset
      knots        <- round(knots, 8)
      nknots       <- length(knots)
      df           <- length(knots) - 1
    }
   
    
    ##########################################################################
    ##########################################################################
   
    # when sigmaxsi i.e, ls model not evaluated, the sigmaxoffset remains a
    # a string 'sigmaxoffset' as list. Therefore, set it to NULL or 0
    # Bettwe would be assign numeric values such as 0/10 etc to 
    # sigmaxoffsetsi and not sigmaxoffset
    # same is true for xoffset
    # The xoffset and sigmaxoffset are used in model_info values
    # TODO: CORRECT THE ABOVE
    
    if(is.na(sigmaxsi) |  sigmaxsi != "NA") {
      sigmaxoffset <- sigmaxoffsetsi
    }

    # should be ger _ls
    if(!is.na(sigmaxsi) &  sigmaxsi != "NA") {
      if (is.numeric(ept(sigmaknotssi))) {
        sigmaknots <- ept(sigmaknotssi)
        # sigmaknotssi should take precedence over sigmadf
        # Since sigmadf is automatically set as df, need to shut it off
        sigmadfsi  <- "NA"
      }
      
      if (is.numeric(ept(sigmadfsi))) {
        sigmaknots <- (unname(gkn(datai[[sigmaxsi]], 
                                  ept(sigmadfsi), ept(sigmaboundsi))))
      }
      
      if(sigmabstartsi == 'sigmaxoffset') {
        sigmabstartsi <- sigmaxoffsetsi
      }
      
      if(sigmaxoffsetsi == "apv") {
        stop2c("xoffset can not be 'apv' for sigma")
      }
      
      sigmaxoffset <- eval_xoffset_bstart_args(x = sigmaxsi, 
                                               y = ysi, 
                                               knots = sigmaknots, 
                                               data = datai, 
                                               eval_arg = sigmaxoffsetsi, 
                                               offsetfunsi = 
                                                 sigmaxfunxoffsettransformsi, 
                                               smat = smat,
                                               degree = smat_degree,
                                               intercept = smat_intercept, 
                                               derivs = smat_derivs,
                                               centerval = smat_centerval,
                                               normalize = smat_normalize,
                                               preH = smat_preH,
                                               sfirst = smat_sfirst, 
                                               sparse = smat_sparse,
                                               arg = 'xoffset',
                                               dpar = "sigma",
                                               verbose = verbose)
      
    
      
      sigmabstart <- eval_xoffset_bstart_args(x = sigmaxsi, 
                                              y = ysi, 
                                              knots = sigmaknots, 
                                              data = datai, 
                                              eval_arg = sigmabstartsi, 
                                              offsetfunsi = 
                                                sigmaxfunxoffsettransformsi, 
                                              smat = smat,
                                              degree = smat_degree,
                                              intercept = smat_intercept, 
                                              derivs = smat_derivs,
                                              centerval = smat_centerval,
                                              normalize = smat_normalize,
                                              preH = smat_preH,
                                              sfirst = smat_sfirst, 
                                              sparse = smat_sparse,
                                              arg = 'bstart',
                                              dpar = "sigma",
                                              verbose = verbose)
      
      sigmabstart <- sigmabstart - sigmaxoffset
      
      sigmacstart <- eval_xoffset_cstart_args(x = sigmaxsi, 
                                              y = ysi, 
                                              knots = sigmaknots, 
                                              data = datai, 
                                              eval_arg = sigmacstartsi, 
                                              smat = smat,
                                              degree = smat_degree,
                                              intercept = smat_intercept, 
                                              derivs = smat_derivs,
                                              centerval = smat_centerval,
                                              normalize = smat_normalize,
                                              preH = smat_preH,
                                              sfirst = smat_sfirst, 
                                              sparse = smat_sparse,
                                              arg = 'cstart',
                                              dpar = "sigma",
                                              verbose = verbose)
      
      
      
      
      sigmaxoffset <- round(sigmaxoffset, 8)
      
      if(is_emptyx(sigmaxoffset)) {
        sigmaxoffset      <- 0
      }
      
      
      if(smat == 'bsp' |  smat == 'msp' |  smat == 'isp') {
        if(!is.null(knots)) {
          sigmaknots        <- sigmaknots - sigmaxoffset
          sigmaknots        <- round(sigmaknots, 8)
          sigmanknots       <- length(sigmaknots)
          sigmadf           <- length(sigmaknots) - 1
        } else if(is.null(knots)) {
          sigmaknots        <- NULL
          sigmaknots        <- NULL
          sigmanknots       <- NULL
          sigmadf           <- NULL
        }
      } else {
        sigmaknots        <- sigmaknots - sigmaxoffset
        sigmaknots        <- round(sigmaknots, 8)
        sigmanknots       <- length(sigmaknots)
        sigmadf           <- length(sigmaknots) - 1
      }
    } # if(!is.na(sigmaxsi) &  sigmaxsi != "NA") {
    
    
    ##########################################################################
    ##########################################################################
    
    
    
    
    
    
    ##########################################################################
    ##########################################################################
    
    
    
    ######################################################################
    ######################################################################
    
    # Once offsets defined, make a copy of funs and ifun with offsets included
    xfuntransform2si        <- xfuntransformsi
    bodyoffun               <- deparse(body(xfuntransform2si))
    addtobodyoffun          <- paste0("-", xoffset)
    bodyoffun2              <- paste0(bodyoffun, addtobodyoffun)
    body(xfuntransform2si)  <- str2lang(bodyoffun2)
    ixfuntransform2si       <- inverse_transform(base::body(xfuntransform2si))
    
    
    if(!is.na(sigmaxsi) &  sigmaxsi != "NA") {
      sigmaxfuntransform2si   <- sigmaxfuntransformsi
      bodyoffun               <- deparse(body(sigmaxfuntransform2si))
      addtobodyoffun          <- paste0("-", sigmaxoffset)
      bodyoffun2              <- paste0(bodyoffun, addtobodyoffun)
      body(sigmaxfuntransform2si)   <- str2lang(bodyoffun2)
      sigmaixfuntransform2si  <- 
        inverse_transform(base::body(sigmaxfuntransform2si))
    } else {
      sigmaxfuntransform2si <- NULL
      sigmaixfuntransform2si <- NULL
    }
    
    
    
    ######################################################################
    ######################################################################
    # Once xoffset, sigmaxoffset and lm matrix etc created, 
    # inverse transform x and sigma but leave y as such 
    # Note that xoffset and sigmaxoffset still NULL
    prepare_transformations_args[['data']]         <- datai
    prepare_transformations_args[['xvar']]         <- xsi
    prepare_transformations_args[['yvar']]         <- NULL
    prepare_transformations_args[['sigmaxvar']]    <- sigmaxsi
    prepare_transformations_args[['xfun']]         <- xfuntransformsi
    prepare_transformations_args[['yfun']]         <- NULL
    prepare_transformations_args[['sigmaxfun']]    <- sigmaxfuntransformsi
    prepare_transformations_args[['ixfun']]        <- TRUE
    prepare_transformations_args[['iyfun']]        <- FALSE
    prepare_transformations_args[['sigmaixfun']]   <- TRUE
    prepare_transformations_args[['xoffset']]      <- NULL
    prepare_transformations_args[['sigmaxoffset']] <- NULL
    prepare_transformations_args[['transform']]    <- ""
    prepare_transformations_args[['itransform']]   <- ""

    
    datai <- CustomDoCall(prepare_transformations, prepare_transformations_args)
    
    if(check_for_validy_of_prepare_transformations) {
      check_for_validy_of_prepare_transformations_2 <- datai
    }
    
    ######################################################################
    ######################################################################
    # Now re transform 'x' and 'sigmax' with xoffset and 'sigmaxoffset'
    # but leave y as such 
    # Note below that 'xoffset' and 'sigmaxoffset' are not NULL
    
     prepare_transformations_args[['data']]         <- datai
     prepare_transformations_args[['xvar']]         <- xsi
     prepare_transformations_args[['yvar']]         <- NULL
     prepare_transformations_args[['sigmaxvar']]    <- sigmaxsi
     prepare_transformations_args[['xfun']]         <- xfuntransformsi
     prepare_transformations_args[['yfun']]         <- NULL
     prepare_transformations_args[['sigmaxfun']]    <- sigmaxfuntransformsi
     prepare_transformations_args[['ixfun']]        <- FALSE
     prepare_transformations_args[['iyfun']]        <- FALSE
     prepare_transformations_args[['sigmaixfun']]   <- FALSE
     prepare_transformations_args[['xoffset']]      <- xoffset
     prepare_transformations_args[['sigmaxoffset']] <- sigmaxoffset
     prepare_transformations_args[['transform']]    <- ""
     prepare_transformations_args[['itransform']]   <- ""
    
    datai <- CustomDoCall(prepare_transformations, prepare_transformations_args)
    
    if(check_for_validy_of_prepare_transformations) {
      check_for_validy_of_prepare_transformations_3 <- datai
    }
    
    #################################################################
    #################################################################
    # This below just for final checks to confirm transformations work fine 
    
    if(check_for_validy_of_prepare_transformations) {
      prepare_transformations_args[['data']]         <- 
        check_for_validy_of_prepare_transformations_3
      prepare_transformations_args[['ixfun']]        <- TRUE
      prepare_transformations_args[['sigmaixfun']]   <- TRUE
      prepare_transformations_args[['transform']]    <- ""
      prepare_transformations_args[['itransform']]   <- ""
      
      check_for_validy_of_prepare_transformations_4 <- 
        CustomDoCall(prepare_transformations, prepare_transformations_args)
      
      prepare_transformations_args[['data']]         <- 
        check_for_validy_of_prepare_transformations_4
      prepare_transformations_args[['ixfun']]        <- FALSE
      prepare_transformations_args[['sigmaixfun']]   <- FALSE
      prepare_transformations_args[['transform']]    <- ""
      prepare_transformations_args[['itransform']]   <- ""
      
      check_for_validy_of_prepare_transformations_5 <- 
        CustomDoCall(prepare_transformations, prepare_transformations_args)
      
      # Ignore outcome ysi because it remains transformed
      check_for_validy_of_prepare_transformations_0 <- 
        check_for_validy_of_prepare_transformations_0 %>% 
        # dplyr::select(-ysi)
        dplyr::select(-dplyr::all_of(ysi))
      check_for_validy_of_prepare_transformations_4 <- 
        check_for_validy_of_prepare_transformations_4 %>% 
        dplyr::select(-dplyr::all_of(ysi))
      
      check_for_validy_of_prepare_transformations_3 <- 
        check_for_validy_of_prepare_transformations_3 %>% 
        dplyr::select(-dplyr::all_of(ysi))
      check_for_validy_of_prepare_transformations_5 <- 
        check_for_validy_of_prepare_transformations_5 %>% 
        dplyr::select(-dplyr::all_of(ysi))
      
      # 
      if(!isTRUE(all.equal(check_for_validy_of_prepare_transformations_0,
                    check_for_validy_of_prepare_transformations_4))) {
        stop2c("Something wrong with 'prepare_transformations'")
      }
      if(!isTRUE(all.equal(check_for_validy_of_prepare_transformations_3,
                    check_for_validy_of_prepare_transformations_5))) {
        stop2c("Something wrong with 'prepare_transformations'")
      }
    } # end if(check_for_validy_of_prepare_transformations) {
    
   
    ##################################################################
    ##################################################################
    
    iknots <- checkgetiknotsbknots(knots, 'iknots')
    bknots <- checkgetiknotsbknots(knots, 'bknots')
    
    SplineCall <- substitute(TEMPNAME(x = datai[[xsi]],
                                      knots = iknots,
                                      bknots = bknots,
                                      degree = smat_degree,
                                      intercept = smat_intercept,
                                      derivs = smat_derivs,
                                      centerval = smat_centerval,
                                      normalize = smat_normalize,
                                      preH = smat_preH,
                                      sfirst = smat_sfirst,
                                      sparse = smat_sparse))
    
    
    if(smat == 'rcs') {
      SplineCall[[1]] <- quote(GS_rcs_call)
    } else if(smat == 'nsp') {
      SplineCall[[1]] <- quote(GS_nsp_call)
    } else if(smat == 'nsk') {
      SplineCall[[1]] <- quote(GS_nsk_call)
    } else if(smat == 'bsp') {
      SplineCall[[1]] <- quote(GS_bsp_call)
    } else if(smat == 'msp') {
      SplineCall[[1]] <- quote(GS_msp_call)
    } else if(smat == 'isp') {
      SplineCall[[1]] <- quote(GS_isp_call)
    }
    
    mat_s           <- eval(SplineCall)
    SplineCall[[2]] <- quote(x)
    

    #################################################################
    #################################################################
    
    if(!is.null(knots_selection)) {
      if(knots_selection[['when']] == 'ac') {
        dataset_temp_knots_selection <- datai
        dataset_temp_knots_selection[[xsi]] <- datai[[xsi]]
        knots_selection_arg[['dataset']] <- dataset_temp_knots_selection
        knots_selection_arg[['dependent']] <- ysi
        knots_selection_arg[['independents']] <- xsi
        knots_selection_arg[['bknots']] <- checkgetiknotsbknots(knots, 'bknots')
        knots_old <- knots
        model_old <- NULL
        model_new <- do.call(get_choose_model, knots_selection_arg)
        if(knots_selection_arg[['return_df']]) {
          return(model_new)
        }
        if(knots_selection_arg[['kspace']] == "un") {
          knots_new <- model_new$ns$knot_placements$fullknots
        } 
        if(knots_selection_arg[['kspace']] == "nu") {
          knots_new <- model_new$ns_nu$knot_placements$fullknots
        }
        getmaxdp  <- max(get_decimal_places(knots_old))
        knots_new <- round(knots_new, getmaxdp)
        knots     <- knots_new # These knots will be passed to bsitar
        #######
        plot_object <- get_print_return_obj(knots = knots_old, 
                                            model = model_old,
                                            model_new = model_new, 
                                            knots_new = knots_new, 
                                            nys = nys,
                                            xsi = xsi, 
                                            ysi = ysi,
                                            select = 
                                              knots_selection_arg[['select']],
                                            kspace = 
                                              knots_selection_arg[['kspace']],
                                            what = knots_selection_what,
                                            print = knots_selection_print, 
                                            list_arg = knots_selection_arg,
                                            return = knots_selection_return,
                                            verbose = verbose)
        if(knots_selection_print) {
          print(plot_object)
        }
        if(knots_selection_return) {
          if(knots_selection_what == 'knots') {
            retout <- list(knots_old = knots_old,
                           knots_new = knots_new)
            return(retout)
          } 
          return(plot_object)
        }
        #######
      } # if(knots_selection[['when']] == 'ac') {
    } # if(!is.null(knots_selection)) {
    
    
    #################################################################
    #################################################################
    
    # Although stype has when option, somehow it does not work
    # Therefore, let the default 'bc' keep working 
    if(is.null(knots_selection)) {
      if(smat_when == "ac") {
        smat_knots_plot_arg <- get_knost_from_df_arg
        smat_knots_plot_arg[['dataset']] <- datai
        smat_knots_plot_arg[['dependent']] <- ysi
        smat_knots_plot_arg[['independents']] <- xsi
        plot_object <- get_print_return_obj(knots = knots, 
                                            model = NULL,
                                            model_new = NULL, 
                                            knots_new = NULL, 
                                            nys = nys,
                                            xsi = xsi, 
                                            ysi = ysi,
                                            what = smat_what,
                                            print = smat_print, 
                                            list_arg = smat_knots_plot_arg,
                                            return = smat_return,
                                            verbose = verbose)
        
        if(smat_print) {
          print(plot_object)
        }
        if(smat_return) {
          return(plot_object)
        }
      } # if(smat_when == "ac") {
    } # if(is.null(knots_selection)) {
   
   
   
    #################################################################
    #################################################################
    
    # After setting xfuntransformsi etc.., 
    # re-assign xfunsi etc... based on xfuntransformsi, yfuntransformsi etc.., 
    
    xfunsi      <- strsplit(gsub_space(deparse(body(xfuntransformsi))), 
                                "\\(")[[1]][1]
    yfunsi      <- strsplit(gsub_space(deparse(body(yfuntransformsi))),
                                "\\(")[[1]][1]
    sigmaxfunsi <- strsplit(gsub_space(deparse(body(sigmaxfuntransformsi))), 
                                "\\(")[[1]][1]
    
    sigmayfunsi <- strsplit(gsub_space(deparse(body(sigmayfuntransformsi))), 
                            "\\(")[[1]][1]
    
    if(xfunsi == "")      xfunsi      <- NULL
    if(yfunsi == "")      yfunsi      <- NULL
    if(sigmaxfunsi == "") sigmaxfunsi <- NULL
    if(sigmayfunsi == "") sigmayfunsi <- NULL
    
    # sqrt is same as ^0.5
    if(grepl("\\^0.5$", xfunsi)) {
      xfunsi      <- "sqrt"
    }
    if(grepl("\\^0.5$", yfunsi)) {
      yfunsi      <- "sqrt"
    }
    if(grepl("\\^0.5$", sigmaxfunsi)) {
      sigmaxfunsi <- "sqrt"
    }
    if(grepl("\\^0.5$", sigmayfunsi)) {
      sigmayfunsi <- "sqrt"
    }
    
    xfunxoffsetsi <- 
      strsplit(gsub_space(deparse(body(xfuntransformsi))), 
                            "\\(")[[1]][1]
    
    sigmaxfunxoffsetsi <- 
      strsplit(gsub_space(deparse(body(sigmaxfuntransformsi))), 
                            "\\(")[[1]][1]
    
    
    if(xfunxoffsetsi == "")      xfunxoffsetsi      <- NULL
    if(sigmaxfunxoffsetsi == "") sigmaxfunxoffsetsi <- NULL
    
    xfunvalue             <- xfunsi
    yfunvalue             <- yfunsi
    sigmaxfunvalue        <- sigmaxfunsi
    sigmayfunvalue        <- sigmayfunsi
    xfunxoffsetvalue      <- xfunxoffsetsi
    sigmaxfunxoffsetvalue <- sigmaxfunxoffsetsi
    
    #################################################################
    #################################################################
    
    includefunnameslistname  <- 'include_fun_names'
    funlist_r_name           <- 'funlist_r'  
    sigmafunlist_r_name      <- 'sigmafunlist_r'
    SplineCall_name          <- "SplineCall"
    sigmavarfunlist_r_name   <- 'sigmavarfunlist_r'
    sigmabasicfunlist_r_name <- 'sigmabasicfunlist_r'
   
    getxname     <- "getX"
    getknotsname <- "getKnots"
    getpreHname  <- "getpreH"
    
    SplineCallname  <- "SplineCall"
    
    if (nys > 1) {
      spfncname       <- paste0(ysi, "_", spfncname)
      getxname        <- paste0(ysi, "_", getxname)
      getknotsname    <- paste0(ysi, "_", getknotsname)
      getpreHname     <- paste0(ysi, "_", getpreHname)
      SplineCallname  <- paste0(ysi, "_", SplineCallname)
    }  
    
    SplineCallnamelist[[ii]]  <- SplineCall_name
    SplineCallvaluelist[[ii]] <- SplineCall
   
    spfncname_c <- c(spfncname_c, spfncname)
    
    spfun_collect <- c(spfun_collect, 
                       c(spfncname, paste0(spfncname, "_", c("d1", "d2")))
                       )
    
    # QR, also 'decomp_editcode' = TRUE conflicts with init != random
    decomp_editcode <- FALSE
    if(select_model == 'rcs') {
      decomp_editcode <- FALSE
    }
    
    # This to check s covs - re
    checkscovsi <-  getcovlist(s_formulasi)
    
    # This can be set to TRUE for RCS
    if(!is.null(checkscovsi)) {
      add_b_Qr_genquan_s_coef <- FALSE
    } else {
      add_b_Qr_genquan_s_coef <- FALSE
    }
    
    #################################################################
    #################################################################
    
    SbasisN <- ncol(mat_s)
    
    internal_function_args_names <-
      c(
        "fixedsi",
        "randomsi",
        "spfncname",
        "getxname",
        "getknotsname",
        "getpreHname",
        "match_sitar_a_form",
        'match_sitar_d_form',
        'd_as_random_only',
        "d_adjustedsi",
        'xfunsi',
        'yfunsi',
        'xoffset',
        'sigmaxoffset',
        'brms_arguments',
        'select_model',
        'decomp', 
        'decomp_editcode',
        'nys',
        'checkscovsi',
        'add_b_Qr_genquan_s_coef',
        'add_rcsfunmatqrinv_genquant',
        "verbose",
        "smat",
        "smat_degree",
        "smat_intercept",
        "smat_derivs",
        "smat_centerval",
        "smat_normalize",
        "smat_preH",
        "smat_sfirst",
        "smat_sparse",
        "smat_moi",
        "smat_check_sparsity",
        "smat_include_stan",
        "smat_include_path",
        "SplinefunxPre",
        "Splinefunxsuf",
        "SplinefunxR",
        "SplinefunxStan",
        "ii",
        "set_xfunsi",
        "set_yfunsi",
        "set_sigmaxfunsi",
        "set_xfunxoffsetsi",
        "set_sigmaxfunxoffsetsi",
        "xfuntransformsi",
        "yfuntransformsi",
        "sigmaxfuntransformsi",
        "xfunxoffsettransformsi",
        "sigmaxfunxoffsettransformsi",
        "fast_nsk",
        "sum_zero",
        "QR_Xmat",
        "QR_center",
        "QR_complete",
        "QR_flip",
        "QR_scale",
        "SbasisN",
        "sigmaxsi",
        'familysi',
        'sigmaxfunsi',
        'sigmayfunsi',
        'sigmayfuntransformsi',
        'family_link_sigma'
        )
    
    internal_function_args <- list()
    internal_function_args <- mget(internal_function_args_names)
    
    if (verbose) {
      if (ii == 1) {
        setmsgtxt <- paste0("\n Preparing function")
        if (displayit == 'msg') {
          message2c(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
    }
    
    get_s_r_funs <- prepare_function_nsp_rcs(
      x = xsi,
      y = ysi,
      id = idsi,
      knots = knots,
      nknots = nknots,
      data = datai,
      internal_function_args = internal_function_args)
    
    funlist[ii]     <- get_s_r_funs[['rcsfun']]
    funlist_r[[ii]] <- get_s_r_funs[['r_funs']]
    gq_funs[[ii]]   <- get_s_r_funs[['gq_funs']]
    
    include_fun_nameslist[[ii]] <- get_s_r_funs[['include_fun_names']]
    
    
    ##########################################################################
    # prepare sigma function
    ##########################################################################    
    # why ?
    if(!set_model_sigma_by_ls) {
      sigmad_adjustedsi <- 'NULL'
    }
    
    # Define sigma function
    if(set_model_sigma_by_ls) {
      sigmagetxname       <- paste0("sigma", getxname)
      sigmagetknotsname   <- paste0("sigma", getknotsname)
      sigmagetpreHname    <- paste0("sigma", getpreHname)
      if (nys > 1) {
        sigmagetxname     <- paste0(ysi, "_", sigmagetxname)
        sigmagetknotsname <- paste0(ysi, "_", sigmagetknotsname)
        sigmagetpreHname  <- paste0(ysi, "_", sigmagetpreHname)
      } 

      sigmaspfun_collect <-
        c(sigmaspfun_collect, c(sigmaspfncname, 
                                paste0(sigmaspfncname, "_",  c("d1", "d2"))
                                ))
      
      sigmadecomp_editcode <- FALSE
      if(select_model == 'rcs') {
        sigmadecomp_editcode <- FALSE
      }
      
      # This to check s covs - re
      sigmacheckscovsi <-  getcovlist(s_formulasi)
      
      # This can be set to TRUE for RCS
      if(!is.null(sigmacheckscovsi)) {
        sigmaadd_b_Qr_genquan_s_coef <- FALSE
        # if(select_model == 'rcs') add_b_Qr_genquan_s_coef <- TRUE
      } else {
        sigmaadd_b_Qr_genquan_s_coef <- FALSE
      }
      
      # This controls whether to add scode for genquant block for QR model
      # Relevant in prepare_function
      sigmaadd_rcsfunmatqrinv_genquant <- FALSE # TRUE
      
      # copy internal_function_args but later replace them by sigma args
      sigmainternal_function_args <- internal_function_args
      
   
      sigmamatch_sitar_d_form          <- match_sitar_d_form
      sigmad_adjustedsi                <- d_adjustedsi
      sigmayfunsi                      <- sigmayfunsi
      sigmabrms_arguments              <- brms_arguments
      sigmaselect_model                <- select_model
      sigmadecomp                      <- decomp
      sigmadecomp_editcode             <- decomp_editcode
      sigmanys                         <- nys
      sigmacheckscovsi                 <- checkscovsi
      sigmaadd_b_Qr_genquan_s_coef     <- add_b_Qr_genquan_s_coef
      sigmaadd_rcsfunmatqrinv_genquant <- add_rcsfunmatqrinv_genquant
      sigmainternal_function_args[['fixedsi']] <- sigmafixedsi
      sigmainternal_function_args[['randomsi']] <- sigmarandomsi
      sigmainternal_function_args[['spfncname']] <- sigmaspfncname
      sigmainternal_function_args[['getxname']] <- sigmagetxname
      sigmainternal_function_args[['getknotsname']] <- sigmagetknotsname
      sigmainternal_function_args[['match_sitar_a_form']] <- 
        sigmamatch_sitar_a_form
      sigmainternal_function_args[['match_sitar_d_form']] <- 
        sigmamatch_sitar_d_form
      sigmainternal_function_args[['d_adjustedsi']] <- 
        sigmad_adjustedsi
      sigmainternal_function_args[['sigmaxfunsi']] <- 
        sigmaxfunsi
      sigmainternal_function_args[['sigmayfunsi']] <- 
        sigmayfunsi
      sigmainternal_function_args[['sigmaxoffset']] <- 
        sigmaxoffset
      sigmainternal_function_args[['brms_arguments']] <- 
        sigmabrms_arguments
      sigmainternal_function_args[['select_model']] <- 
        sigmaselect_model
      sigmainternal_function_args[['decomp']] <- 
        sigmadecomp
      sigmainternal_function_args[['decomp_editcode']] <- 
        sigmadecomp_editcode
      sigmainternal_function_args[['nys']] <- 
        sigmanys
      sigmainternal_function_args[['checkscovsi']] <- 
        sigmacheckscovsi
      sigmainternal_function_args[['add_b_Qr_genquan_s_coef']] <- 
        sigmaadd_b_Qr_genquan_s_coef
      sigmainternal_function_args[['add_rcsfunmatqrinv_genquant']] <- 
        sigmaadd_rcsfunmatqrinv_genquant
      
      sigmainternal_function_args[['dpar_function']] <- "sigma"
      
      if (verbose) {
        if (ii == 1) {
          setmsgtxt <- paste0("\n Preparing function for sigma")
          if (displayit == 'msg') {
            message2c(setmsgtxt)
          } else if (displayit == 'col') {
            col <- setcolh
            cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
          }
        }
      }
      
      if (is.null(sigmaidsi[[1]][1]) |
          sigmaidsi == "NULL") {
        sigmaidsi <- idsi
      }

      # These are copied from the mu part
      # Note that idsi is not used and ysi is placeholder
      sigmaysi   <- ysi
      sigmadatai <- datai
      sigmaxsi   <- sigmaxsi
      
      sigmaget_s_r_funs <-
        prepare_function_nsp_rcs(
          x = sigmaxsi,
          y = sigmaysi,
          id = sigmaidsi,
          knots = sigmaknots,
          nknots = sigmanknots,
          data = sigmadatai,
          internal_function_args = sigmainternal_function_args
        )
      
      sigmafunlist[ii]     <- sigmaget_s_r_funs[['rcsfun']]
      sigmafunlist_r[[ii]] <- sigmaget_s_r_funs[['r_funs']]
      sigmagq_funs[[ii]]   <- sigmaget_s_r_funs[['gq_funs']]
      
      include_fun_nameslist[[ii]] <- c(include_fun_nameslist[[ii]], 
                                       sigmaget_s_r_funs[['include_fun_names']])
      
    } # if(set_model_sigma_by_ls) {
    
    ##########################################################################
    # prepare sigma function - end
    ##########################################################################  
    
    ############################################################################
    # add define vf()
    ############################################################################
    
    if(set_model_sigma_by_fz |
       set_model_sigma_by_fp |
       set_model_sigma_by_fe |
       set_model_sigma_by_ve | 
       set_model_sigma_by_vp |
       set_model_sigma_by_cp |
       set_model_sigma_by_mp |
       set_model_sigma_by_me |
       set_model_sigma_by_rp |
       set_model_sigma_by_re ) {
      sigmavarget_s_r_funs <- list()
      
      if(any(unlist(add_identityfun_c))) {
        add_identityfun_stan_fun <- TRUE
      } else {
        add_identityfun_stan_fun <- FALSE
      }
      
      add_identityfun_stan_fun_scode <- 
        "vector identity(vector x) {
          return (x);
        }"
      
      add_absifel_stan_fun_scode <- 
        "vector negzero(vector x) {
        int N = size(x);
        vector[N] out;
        for (i in 1:N) {
          if(x[i] < 0.0) {
            out[i] = 0.0;
          } else {
            out[i] = (x[i]);
          }
        }
        return (out);
        }"
      

      if(set_model_sigma_by_fz) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor) {
          return (parameter0 + (negzero(parameter1 .* predictor)));
        }"
      } else if(set_model_sigma_by_fp) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor) {
          return (parameter0 + (parameter1 .* log(abs(predictor))));
        }"
      } else if(set_model_sigma_by_fe) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor) {
          return (parameter0 + (parameter1 .* predictor));
        }"
      } else if(set_model_sigma_by_ve) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor) {
          return (parameter0 + (parameter1 .* predictor));
        }"
      } else if(set_model_sigma_by_vp) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor) {
          return (parameter0 + (parameter1 .* log(abs(predictor))));
        }"
      } else if(set_model_sigma_by_cp) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector parameter2, 
        vector predictor) {
          return (parameter0 + (parameter1 .* (parameter2 + log(abs(predictor)))));
        }"
      } else if(set_model_sigma_by_mp) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor) {
          return (parameter0 + (parameter1 .* log(abs(sqrt(predictor)))));
        }"
      } else if(set_model_sigma_by_me) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor) {
          return (parameter0 + (parameter1 .* sqrt(predictor)));
        }"
      } else if(set_model_sigma_by_rp) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor, vector y) {
          return (parameter0 + (parameter1 .* log(abs(sqrt(y-predictor)))));
        }"
      } else if(set_model_sigma_by_re) {
        add_sigmavarfun_stan_fun_scode <- "
        vector xxx_set_name_xxx(vector parameter0, vector parameter1, vector predictor, vector y) {
          return (parameter0 + (parameter1 .* sqrt( abs( (y-predictor)  ) )  ));
        }"
      }
      
      
      add_sigmavarfun_stan_fun_scode <- 
        gsub("xxx_set_name_xxx", sigmavarspfncname,
             add_sigmavarfun_stan_fun_scode, fixed = TRUE)
      
      # set identity and negzero, which are common without ys, null if ii > 1
      if(ii > 1) {
        add_identityfun_stan_fun_scode <- NULL
        add_absifel_stan_fun_scode     <- NULL
      }
      
      collect_sigmavar_stan_fun_scode_all_c <- c()
      
      add_identityfun_stan_fun_scode_all <- ""
      if(add_identityfun_stan_fun) {
        collect_sigmavar_stan_fun_scode_all_c <- 
          c(collect_sigmavar_stan_fun_scode_all_c,
            add_identityfun_stan_fun_scode)
      }
      
      if(set_model_sigma_by_fz) {
        collect_sigmavar_stan_fun_scode_all_c <- 
          c(collect_sigmavar_stan_fun_scode_all_c,
            add_absifel_stan_fun_scode)
      }
      
      
      collect_sigmavar_stan_fun_scode_all_c <- 
        c(collect_sigmavar_stan_fun_scode_all_c,
          add_sigmavarfun_stan_fun_scode)
      
      add_identityfun_stan_fun_scode_all <- 
        paste0(collect_sigmavar_stan_fun_scode_all_c, collapse = "\n") 
      
      #########################################################################
      #########################################################################
      # Extract r funs from code
      
      add_identityfun_stan_fun_scode_all_c_r_funs_c <- c()
      add_to_include_fun_nameslist_c <- c()
      for (i in 1:length(collect_sigmavar_stan_fun_scode_all_c)) {
        set_str <- collect_sigmavar_stan_fun_scode_all_c[i]
        set_str <- trimws(set_str, which = c("both", "left", "right"), 
                          whitespace = "[ \t\r\n]")
        
        if(grepl("^\n", set_str)) {
          set_start <- "\n"
        } else if(!grepl("^\n", set_str)) {
          set_start <- ""
        }
        set_str_name <- replace_string_part(x = set_str, start = set_start, 
                                            end = "(",  replace = "",
                                            extract = T, cat_str = FALSE, 
                                            exclude_start = T, exclude_end = T)
        set_SplinefunxStan <- strsplit(set_str_name, " ")[[1]][-1]
        set_SplinefunxR    <- gsub(paste0("", "_", "stan"), "", 
                                   set_SplinefunxStan, fixed = T)
        
        getsighamrfun <- 
          extract_r_fun_from_scode_sigma(set_str,
                                         SplinefunxStan = set_SplinefunxStan,
                                         SplinefunxR = set_SplinefunxR) 
        
        add_to_include_fun_nameslist_c <- 
          c(add_to_include_fun_nameslist_c, set_SplinefunxR)
        
        add_identityfun_stan_fun_scode_all_c_r_funs_c <- 
          c(add_identityfun_stan_fun_scode_all_c_r_funs_c, getsighamrfun)
      }
      
      #########################################################################
      #########################################################################
      
      add_identityfun_stan_fun_scode_all <- 
        paste0("\n", add_identityfun_stan_fun_scode_all)
      add_identityfun_stan_fun_scode_all <- 
        remove_spaces_and_tabs(add_identityfun_stan_fun_scode_all)
      
      add_identityfun_stan_fun_scode_all_c_r_funs_c <- 
        remove_spaces_and_tabs(add_identityfun_stan_fun_scode_all_c_r_funs_c)
      
      
      sigmavarget_s_r_funs[['rcsfun']] <- add_identityfun_stan_fun_scode_all
      sigmavarget_s_r_funs[['r_funs']] <- 
        add_identityfun_stan_fun_scode_all_c_r_funs_c
      sigmavarget_s_r_funs[['gq_funs']] <- NULL
      sigmavarget_s_r_funs[['include_fun_names']] <- 
        add_to_include_fun_nameslist_c
      
      sigmavarfunlist[ii]     <- sigmavarget_s_r_funs[['rcsfun']]
      sigmavarfunlist_r[[ii]] <- sigmavarget_s_r_funs[['r_funs']]
      sigmavargq_funs[[ii]]   <- sigmavarget_s_r_funs[['gq_funs']]
      
      # now adding this, because it stores outcome specific vf()
      include_fun_nameslist[[ii]] <-
        c(include_fun_nameslist[[ii]],
          sigmavarget_s_r_funs[['include_fun_names']])
      
    } # if(set_model_sigma_by_fz | ......
    
    ############################################################################
    # end add define vf()
    ############################################################################
    
    
    ############################################################################
    # end add define sigma basic functions
    ############################################################################
   
    if(set_model_sigma_by_ba) {
      sigmabasicget_s_r_funs <- list()
      # This getouttemp has been collected above, search below line 
      # replace functions with :: / ::: with _ in sigma_formula_manualsi
      # getouttemp <- get_function_names_code_from_string(sigma_formula_manualsi)
      sigmabasicget_s_r_funs[['rcsfun']] <- NULL
      sigmabasicget_s_r_funs[['r_funs']] <-  getouttemp[['code']]
      sigmabasicget_s_r_funs[['gq_funs']] <- NULL
      sigmabasicget_s_r_funs[['include_fun_names']] <- getouttemp[['name']]
      sigmabasicfunlist[ii]     <- sigmabasicget_s_r_funs[['rcsfun']]
      sigmabasicfunlist_r[[ii]] <- sigmabasicget_s_r_funs[['r_funs']]
      sigmabasicgq_funs[[ii]]   <- sigmabasicget_s_r_funs[['gq_funs']]
      include_fun_nameslist[[ii]] <- 
        c(include_fun_nameslist[[ii]], 
          sigmabasicget_s_r_funs[['include_fun_names']])
    }
    
    
    ############################################################################
    # end add define sigma basic functions
    ############################################################################
   
    # if decomp = QR
    # make mat_s as Q, so that correct lm based initials
    if(!is.null(decomp)) {
      if(decomp == 'QR') {
        QR_decomp_R_out <- QR_decomp_R(X = mat_s, 
                                       center = QR_center, 
                                       complete = QR_complete, 
                                       flip = QR_flip, 
                                       scale = QR_scale)
        mat_s <- QR_decomp_R_out[['Q']]
      }
    }
    

    #################################################
    internal_formula_args_names <-
      c(
        "a_formulasi",
        "b_formulasi",
        "c_formulasi",
        "d_formulasi",
        "e_formulasi",
        "f_formulasi",
        "g_formulasi",
        "h_formulasi",
        "i_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        "e_formula_grsi",
        "f_formula_grsi",
        "g_formula_grsi",
        "h_formula_grsi",
        "i_formula_grsi",
        "s_formula_grsi",
        "terms_rhssi",
        "sigma_formulasi",
        "sigma_formula_grsi",
        "sigma_formula_manualsi",
        "dpar_formulasi",
        "autocor_formi",
        "subindicatorsi",
        "fixedsi",
        "randomsi",
        "univariate_by",
        "multivariate",
        "group_arg",
        "sigma_group_arg",
        "df",
        "mat_s",
        "spfncname",
        'nys',
        'ysi',
        'familysi',
        'custom_family',
        'xfunsi',
        'xoffset',
        'match_sitar_a_form',
        'match_sitar_d_form',
        'd_as_random_only',
        "a_formula_gr_strsi",
        "b_formula_gr_strsi",
        "c_formula_gr_strsi",
        "d_formula_gr_strsi",
        "e_formula_gr_strsi",
        "f_formula_gr_strsi",
        "g_formula_gr_strsi",
        "h_formula_gr_strsi",
        "i_formula_gr_strsi",
        "s_formula_gr_strsi",
        "sigma_formula_gr_strsi",
        "set_higher_levels",
        "sigma_set_higher_levels",
        "select_model",
        "verbose",
        "unusedsi",
        "smat",
        "smat_degree",
        "smat_intercept",
        "smat_derivs",
        "smat_centerval",
        "smat_normalize",
        "smat_preH",
        "smat_sfirst",
        "smat_sparse",
        "smat_moi",
        "smat_check_sparsity",
        "smat_include_stan",
        "smat_include_path",
        "SplinefunxPre",
        "Splinefunxsuf",
        "SplinefunxR",
        "SplinefunxStan",
        "set_model_sigma_by_ba",
        "set_model_sigma_by_ls",
        "sigma_formula_manualsi_set",
        "sigma_formula_manual_prior_via_sigma_formula",
        "SbasisN",
        'sigmaxfunsi',
        'sigmayfunsi',
        'family_link_sigma'
      )

    internal_formula_args <- list()
    internal_formula_args <- mget(internal_formula_args_names)
    
    if (verbose) {
      if (ii == 1) {
        setmsgtxt <- paste0("\n Preparing formula")
        if (displayit == 'msg') {
          message2c(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
    }
    
    formula_bf <-
      prepare_formula(
        x = xsi,
        y = ysi,
        id = idsi,
        knots = knots,
        nknots = nknots,
        data = datai,
        internal_formula_args = internal_formula_args)
    
    
    formula_bf_to_check_loop <- formula_bf
    
    list_out <- attr(formula_bf, "list_out")
   
    attributes(formula_bf) <- NULL
    
    eout <- list2env(list_out)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
   
    group_arg$groupvar  <- group_arg_groupvar
    multivariate$rescor <- multivariate_rescor  # why ?
    univariate_by$by    <- univariate_by_by
    covariates_         <- covariates_
    sigmacovariates_    <- sigmacovariates_
    set_higher_levels   <- set_higher_levels
    
    # Check covariate and gr(..., by) are distinct
    abc_grby    <- list_out[['abc_grby']]
    sigma_grby  <- list_out[['sigma_grby']]
    abc_check_grby_covariates_ <- intersect(abc_grby, covariates_)
    sigma_check_grby_covariates_ <- intersect(sigma_grby, sigmacovariates_)
    abc_check_grby_covariates_msg <- 
    paste0("The names of covariate(s) and the 'by' variable in ",
    "gr(..., by=) for '%s' must not be same. The variables can ",
    "be exactly same but have different names. For example, if ",
    "covariate(s) in the fixed effects '%s' are 'covname', then ",
    " make a copy of this variables such as 'covnameid' in the data ",
    "frame and then use it as a by variable gr(..., by=covnameid)")
    
    if(!is_emptyx(abc_check_grby_covariates_)) {
      stop2(sprintf(
        paste0(abc_check_grby_covariates_msg), 'a, b, c', 'a, b, c'))
    }
    if(!is_emptyx(sigma_check_grby_covariates_)) {
      stop2(sprintf(
          paste0(abc_check_grby_covariates_msg),'sigma', 'sigma'))
    }
    # End Check covariate and gr(..., by) are distinct
   
    
    sigma_set_higher_levels  <- sigma_set_higher_levels
    sigma_group_arg$groupvar <- sigma_arg_groupvar

    lm_val_list <-
      names(eout)[grep(pattern = "^lm_|^lme_", names(eout))]
    lm_val_list <- sort(lm_val_list)
    
    lm_val_list_not <-
      names(eout)[!names(eout) %in%
                    names(eout)[grep(pattern = "^lm_|^lme_", 
                                     names(eout))]]
    lm_val_list_not <- sort(lm_val_list_not)
   
    
    cov_list_names <- ls()[grepl(pattern = "_cov", ls())]
    cov_list_names <-
      cov_list_names[!grepl(pattern = "_init_", cov_list_names)]
    cov_list_names <-
      cov_list_names[!grepl(pattern = "_prior_", cov_list_names)]
    cov_list_names <-
      cov_list_names[!grepl(pattern = "^lm_", cov_list_names)]
    cov_list_names <- sort(cov_list_names)
    
    bflist[[ii]] <- formula_bf
    
    
    ######################################################
    
    loess_fit <- paste0("loess(", ysi, "~", xsi, ",", 'datai', ")")
    loess_fitx <- eval(parse(text = loess_fit))
    
    ymean   <- mean(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
    ymedian <- median(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
    if(select_model == 'sitar' | select_model == 'rcs') {
      ymax  <- max(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
      ymin  <- min(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
      ymaxs <- NULL
    } else if(select_model != 'sitar' & select_model != 'rcs') {
      ymax <- round(max(predict(loess_fitx)), 2)
      ymin  <- min(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
      ymaxs <- round(ymax * 0.95, 2)
    }
    
    ysd     <- sd(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
    ymad    <- mad(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
    xsd     <- sd(datai[[xsi]], na.rm = TRUE) %>% round(., 2)
    
    
    # This for logistic3 model
    ymeanxmin_ysdxmin <- 
      datai %>% dplyr::mutate(XXi := eval(parse(text = xsi))) %>% 
      dplyr::filter(XXi %in% 
                      (min(XXi):min(XXi)+0)) %>% 
      dplyr::mutate(onepic = 1) %>% dplyr::group_by(onepic) %>% 
      dplyr::mutate(temp1 = mean(eval(parse(text = ysi)))) %>% 
      dplyr::mutate(temp2 = sd(eval(parse(text = ysi)))) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(temp1, temp2) %>% 
      dplyr::filter(dplyr::row_number() == 1) %>% 
      unlist() %>% as.numeric()
    
    ymeanxmin <- round(ymeanxmin_ysdxmin[1], 2)
    ysdxmin   <- round(ymeanxmin_ysdxmin[2], 2)
    
    ymeanxmax_ysdxmax <- 
      datai %>% dplyr::mutate(XXi := eval(parse(text = xsi))) %>% 
      dplyr::filter(XXi %in% 
                      (max(XXi):max(XXi)+0)) %>% 
      dplyr::mutate(onepic = 1) %>% dplyr::group_by(onepic) %>% 
      dplyr::mutate(temp1 = mean(eval(parse(text = ysi)))) %>% 
      dplyr::mutate(temp2 = sd(eval(parse(text = ysi)))) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(temp1, temp2) %>% 
      dplyr::filter(dplyr::row_number() == 1) %>% 
      unlist() %>% as.numeric()
    
    ymeanxmax <- round(ymeanxmax_ysdxmax[1], 2)
    ysdxmax   <- round(ymeanxmax_ysdxmax[2], 2)
    
    ymeanxmid_ysdxmid <- 
      datai %>% dplyr::mutate(XXi := eval(parse(text = xsi))) %>% 
      dplyr::filter(XXi %in% 
                      (((max(XXi)-min(XXi))/2):
                         ((max(XXi)-min(XXi))/1.4))) %>% 
      dplyr::mutate(onepic = 1) %>% dplyr::group_by(onepic) %>% 
      dplyr::mutate(temp1 = mean(eval(parse(text = ysi)))) %>% 
      dplyr::mutate(temp2 = sd(eval(parse(text = ysi)))) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(temp1, temp2) %>% 
      dplyr::filter(dplyr::row_number() == 1) %>% 
      unlist() %>% as.numeric()
    
    ymeanxmid <- round(ymeanxmid_ysdxmid[1], 2)
    ysdxmid   <- round(ymeanxmid_ysdxmid[2], 2)
    
    ymeanxmidxmaxdiff <- ymeanxmax - ymeanxmid
    ysdxmidxmaxdiff   <- (ysdxmax + ysdxmin) / 2
    
    ymeanxmidxmaxdiff <- round(ymeanxmidxmaxdiff, 2)
    ysdxmidxmaxdiff   <- round(ysdxmidxmaxdiff, 2)
    
    if(is.na(ymeanxmidxmaxdiff)) ymeanxmidxmaxdiff <- (ymeanxmax + ymeanxmin)/2
    
    # Add missing arguments when restricting to abcd
    if(is.null(pvsi))   pvsi  <- list(NULL)
    if(is.null(apvsi))  apvsi <- list(NULL)
    
    if (!is.null(pvsi[[1]][1]) & pvsi != "NULL") {
      setpv <- (eval(parse(text = pvsi)))
      if(grepl("sitar", select_model)) cstart <- setpv
      if(grepl("pb", select_model))    cstart <- setpv / 5.0
      if(grepl("pb", select_model))    dstart <- setpv
    } else if (!is.null(cstartsi[[1]][1]) & cstartsi != "NULL") {
      setpv <- (cstart)
      if(grepl("sitar", select_model)) cstart <- setpv
      if(grepl("pb", select_model))    cstart <- setpv / 5.0
      if(grepl("pb", select_model))    dstart <- setpv
    } else {
      if(grepl("sitar", select_model)) cstart <- 0.01
      if(grepl("pb", select_model))    cstart <- 0.01
      if(grepl("pb", select_model))    dstart <- 0.01 * 5.0
    }
    
    if (!is.null(apvsi[[1]][1]) & apvsi != "NULL") {
      setapv <- eval(parse(text = apvsi))
      if(grepl("sitar", select_model)) bstart <- setapv
      if(grepl("pb", select_model))    estart <- bstart
    } else if (!is.null(bstartsi[[1]][1]) & bstartsi != "NULL") {
      if(grepl("sitar", select_model)) bstart <- bstart
      if(grepl("pb", select_model))    estart <- bstart
    } else {
      if(grepl("sitar", select_model)) bstart <- 13
      if(grepl("pb", select_model))    estart <- 13
    }
    
    # Set missing start values to 0
    for (gwatxi in letters[1:26]) {
      gwatx__ <- paste0(gwatxi, 'start')
      if(!exists(gwatx__)) assign(gwatx__, 0)
    }
    
    # TODO
    # acov_sd etc setting numeric but later can be worked out to infer from
    # some meaningful way - until then, these will not be used in anywhere
    # and are included here just as placeholders
    lm_a_cov_sd <- 10
    lm_b_cov_sd <- 1
    lm_c_cov_sd <- 0.5
    
    vcov_init_0e <- eval(parse(text =  "vcov_init_0si" ))
    vcov_init_0e <- eval(parse(text =  vcov_init_0e ))
    
    # If vcov_init_0e = TRUE, override random options for below elements
    # Note that these are placeholders, actual setting to zero done later
    
    if(vcov_init_0e) {
      for (inxc in letters[1:26]) {
        assign(paste0(inxc, "_", "init", "_", "sd", "si"), '0')
        assign(paste0(inxc, "_", "cov", "_", "init", "_", "sd", "si"), '0')
      }
      assign('gr_init_corsi',         '0')
      assign('sigma_init_corsi',      '0')
      assign('rsd_init_sigmasi',      '0')
      assign('sigma_init_sdsi',       '0')
      assign('sigma_cov_init_sdsi',   '0')
      assign('dpar_init_sigmasi',     '0')
      assign('dpar_cov_init_sigmasi', '0')
      assign('mvr_init_rescorsi',     '0')
      assign('r_init_zsi',            '0')
    }
    
    
    # Check random initials when fitting 'multivariate' or 'univariate_by' model
    if(nys > 1) {
      param_cov_intit_mismatch_msg <- 
        "If any nonlinear parameter (i.e., a, b, c, d, s) is assigned 'random'
    initials, the same 'random' initials should be used for both the intercept 
    and the corresponding covariate coefficients. Currently, this condition is 
    not met for the nonlinear parameter "
      
      param_cov_intit_mismatch_msg2 <- 
        "This strict consistency in assigning 'random' initials is necessary 
    because each nonlinear parameter is vectorized and therefore cannot 
    contain empty elements."
      
      for (i in letters) {
        for (j in c('beta', 'sd')) {
          param_intit_name     <- paste0(i, "_", "init", "_", j)
          param_cov_intit_name <- paste0(i, "_", "cov", "_", "init", "_", j)
          param_intit_str      <- paste0(param_intit_name, "si")
          param_cov_intit_str  <- paste0(param_cov_intit_name, "si")
          param_intit          <- ept(param_intit_str)
          param_cov_intit      <- ept(param_cov_intit_str)
          # if formula is ~ 0+..., then 'random' conflict corrected internally
          # see if below check can be suppressed when ~0+
          # param_formula_name   <- paste0(i, "_", "formula", "", "")
          if(!is.null(param_intit)) {
            if(param_intit == "random" | param_cov_intit == "random") {
              if(!identical(param_intit, param_cov_intit)) {
                param_cov_intit_mismatch_display <-
                  paste0(param_cov_intit_mismatch_msg, " ", collapse_comma(i),
                         " where the intercept (",
                         collapse_comma(param_intit_name),
                         ") is set to ",  collapse_comma(param_intit),
                         " while the covariate coefficient (",
                         collapse_comma(param_cov_intit_name),
                         ") is set to ",  collapse_comma(param_cov_intit),
                         ".")
                param_cov_intit_mismatch_display <- 
                  paste0(param_cov_intit_mismatch_display, 
                         param_cov_intit_mismatch_msg2)
                stop2c(param_cov_intit_mismatch_display)
              }
            } # if(param_intit == "random" | param_cov_intit == "random") {
          } # if(!is.null(param_intit)) {
        } # for (j in c('beta', 'sd')) {
      } # for (i in letters[1:20]) {
    } # if(nys > 1) {
    
    
    prior_data_internal_names <-
      c(
        lm_val_list,
        "ymean",
        "ymedian",
        "ymax",
        "ymaxs",
        "ymin",
        "ysd",
        "ymad",
        "xsd",
        'ymeanxmin', 
        'ysdxmin', 
        'ymeanxmax', 
        'ysdxmax', 
        'ymeanxmid', 
        'ysdxmid',
        'ymeanxmidxmaxdiff',
        'ysdxmidxmaxdiff',
        "lm_a_cov_sd",
        "lm_b_cov_sd",
        "lm_c_cov_sd",
        "bstart",
        "cstart",
        "dstart",
        "estart"
        )
    
  
    prior_args_internal_names <-
      c(
        lm_val_list_not,
        cov_list_names,
        "a_formulasi",
        "b_formulasi",
        "c_formulasi",
        "d_formulasi",
        "e_formulasi",
        "f_formulasi",
        "g_formulasi",
        "h_formulasi",
        "i_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "fixedsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        "e_formula_grsi",
        "f_formula_grsi",
        "g_formula_grsi",
        "h_formula_grsi",
        "i_formula_grsi",
        "s_formula_grsi",
        "sigma_formulasi",
        "sigma_formula_grsi",
        "sigma_formula_gr_strsi",
        "sigma_formula_manualsi",
        "autocor_formi",
        "randomsi",
        "nabci",
        "nabcrei",
        "univariate_by",
        "multivariate",
        "group_arg",
        "sigma_group_arg",
        "initsi",
        "df",
        "idsi",
        "sigmaidsi",
        "ys",
        "resp",
        "ii",
        "nys",
        "N_J_all",
        "dpar_formulasi",
        "normalize",
        "seed",
        'match_sitar_a_form',
        'match_sitar_d_form',
        'cortimeNlags_var',
        'cortimeNlags',
        "select_model",
        "verbose",
        "SbasisN"
        )
    
    
    prior_data_internal <- list()
    prior_data_internal <- mget(prior_data_internal_names)
    
    
    prior_args_internal <- list()
    prior_args_internal <- mget(prior_args_internal_names)
    
    
    if (!is.null(prior_data[[1]])) {
      get_common_names_lists <-
        intersect(names(prior_data_internal), names(prior_data))
      ttt_n1 <- paste(names(prior_data_internal), collapse = ", ")
      ttt_nn2 <- paste(get_common_names_lists, collapse = ", ")
      if (!identical(get_common_names_lists, character(0))) {
        stop2c(
          "Names in prior_data list should not be following reserved names:",
          "\n ",
          ttt_n1,
          "\n ",
          " Please change the following name(s) ",
          "\n ",
          ttt_nn2
        )
      }
    }
    
    
    
    init_data_internal <- prior_data_internal
    init_args_internal <- prior_args_internal
    
    # Add if(!is.null(a_init_betasi)).. when restricting to abcd
    # check and set default initials (class = b)
    if(!is.null(a_init_betasi)) a_init_betasi <- 
      set_default_inits(select_model_arg, a_init_betasi)
    if(!is.null(b_init_betasi)) b_init_betasi <- 
      set_default_inits(select_model_arg, b_init_betasi)
    if(!is.null(c_init_betasi)) c_init_betasi <- 
      set_default_inits(select_model_arg, c_init_betasi)
    if(!is.null(d_init_betasi)) d_init_betasi <- 
      set_default_inits(select_model_arg, d_init_betasi)
    if(!is.null(e_init_betasi)) e_init_betasi <- 
      set_default_inits(select_model_arg, e_init_betasi)
    if(!is.null(f_init_betasi)) f_init_betasi <- 
      set_default_inits(select_model_arg, f_init_betasi)
    if(!is.null(g_init_betasi)) g_init_betasi <- 
      set_default_inits(select_model_arg, g_init_betasi)
    if(!is.null(h_init_betasi)) h_init_betasi <- 
      set_default_inits(select_model_arg, h_init_betasi)
    if(!is.null(i_init_betasi)) i_init_betasi <- 
      set_default_inits(select_model_arg, i_init_betasi)
    
    
    init_arguments <-
      list(
        a_init_beta = a_init_betasi,
        b_init_beta = b_init_betasi,
        c_init_beta = c_init_betasi,
        d_init_beta = d_init_betasi,
        e_init_beta = e_init_betasi,
        f_init_beta = f_init_betasi,
        g_init_beta = g_init_betasi,
        h_init_beta = h_init_betasi,
        i_init_beta = i_init_betasi,
        s_init_beta = s_init_betasi,
        a_cov_init_beta = a_cov_init_betasi,
        b_cov_init_beta = b_cov_init_betasi,
        c_cov_init_beta = c_cov_init_betasi,
        d_cov_init_beta = d_cov_init_betasi,
        e_cov_init_beta = e_cov_init_betasi,
        f_cov_init_beta = f_cov_init_betasi,
        g_cov_init_beta = g_cov_init_betasi,
        h_cov_init_beta = h_cov_init_betasi,
        i_cov_init_beta = i_cov_init_betasi,
        s_cov_init_beta = s_cov_init_betasi,
        a_init_sd = a_init_sdsi,
        b_init_sd = b_init_sdsi,
        c_init_sd = c_init_sdsi,
        d_init_sd = d_init_sdsi,
        e_init_sd = e_init_sdsi,
        f_init_sd = f_init_sdsi,
        g_init_sd = g_init_sdsi,
        h_init_sd = h_init_sdsi,
        i_init_sd = i_init_sdsi,
        s_init_sd = s_init_sdsi,
        a_cov_init_sd = a_cov_init_sdsi,
        b_cov_init_sd = b_cov_init_sdsi,
        c_cov_init_sd = c_cov_init_sdsi,
        d_cov_init_sd = d_cov_init_sdsi,
        e_cov_init_sd = e_cov_init_sdsi,
        f_cov_init_sd = f_cov_init_sdsi,
        g_cov_init_sd = g_cov_init_sdsi,
        h_cov_init_sd = h_cov_init_sdsi,
        i_cov_init_sd = i_cov_init_sdsi,
        s_cov_init_sd = s_cov_init_sdsi,
        sigma_init_beta = sigma_init_betasi,
        sigma_cov_init_beta = sigma_cov_init_betasi,
        sigma_init_sd = sigma_init_sdsi,
        sigma_cov_init_sd = sigma_cov_init_sdsi,
        rsd_init_sigma = rsd_init_sigmasi,
        dpar_init_sigma = dpar_init_sigmasi,
        dpar_cov_init_sigma = dpar_cov_init_sigmasi,
        autocor_init_acor = autocor_init_acorsi,
        autocor_init_unstr_acor = autocor_init_unstr_acorsi,
        gr_init_cor = gr_init_corsi,
        sigma_init_cor = sigma_init_corsi,
        mvr_init_rescor = mvr_init_rescorsi,
        r_init_z = r_init_zsi
      )
    
    
    if (verbose) {
      if (ii == 1) {
        setmsgtxt <- paste0("\n Preparing priors and initials")
        if (displayit == 'msg') {
          message2c(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
    }
    
   
   
    # Add if(!is.null(a_prior_betasi)).. when restricting to abcd
    # check and set default priors (class = b)
    if(!is.null(a_prior_betasi)) a_prior_betasi <- 
      set_default_priors(select_model_arg, a_prior_betasi)
    if(!is.null(b_prior_betasi)) b_prior_betasi <- 
      set_default_priors(select_model_arg, b_prior_betasi)
    if(!is.null(c_prior_betasi)) c_prior_betasi <- 
      set_default_priors(select_model_arg, c_prior_betasi)
    if(!is.null(d_prior_betasi)) d_prior_betasi <- 
      set_default_priors(select_model_arg, d_prior_betasi)
    if(!is.null(e_prior_betasi)) e_prior_betasi <- 
      set_default_priors(select_model_arg, e_prior_betasi)
    if(!is.null(f_prior_betasi)) f_prior_betasi <- 
      set_default_priors(select_model_arg, f_prior_betasi)
    if(!is.null(g_prior_betasi)) g_prior_betasi <- 
      set_default_priors(select_model_arg, g_prior_betasi)
    if(!is.null(h_prior_betasi)) h_prior_betasi <- 
      set_default_priors(select_model_arg, h_prior_betasi)
    if(!is.null(i_prior_betasi)) i_prior_betasi <- 
      set_default_priors(select_model_arg, i_prior_betasi)
    
    # Add if(!is.null(a_prior_sdsi)).. when restricting to 'abcd'
    # check and set default priors (class = sd)
    if(!is.null(a_prior_sdsi)) a_prior_sdsi <- 
      set_default_priors(select_model_arg, a_prior_sdsi)
    if(!is.null(b_prior_sdsi)) b_prior_sdsi <- 
      set_default_priors(select_model_arg, b_prior_sdsi)
    if(!is.null(c_prior_sdsi)) c_prior_sdsi <- 
      set_default_priors(select_model_arg, c_prior_sdsi)
    if(!is.null(d_prior_sdsi)) d_prior_sdsi <- 
      set_default_priors(select_model_arg, d_prior_sdsi)
    if(!is.null(e_prior_sdsi)) e_prior_sdsi <- 
      set_default_priors(select_model_arg, e_prior_sdsi)
    if(!is.null(f_prior_sdsi)) f_prior_sdsi <- 
      set_default_priors(select_model_arg, f_prior_sdsi)
    if(!is.null(g_prior_sdsi)) g_prior_sdsi <- 
      set_default_priors(select_model_arg, g_prior_sdsi)
    if(!is.null(h_prior_sdsi)) h_prior_sdsi <- 
      set_default_priors(select_model_arg, h_prior_sdsi)
    if(!is.null(i_prior_sdsi)) i_prior_sdsi <- 
      set_default_priors(select_model_arg, i_prior_sdsi)
   
    set_priors_initials_agrs <- list()

    set_priors_initials_agrs $ a_prior_beta <- a_prior_betasi
    set_priors_initials_agrs $ b_prior_beta <- b_prior_betasi
    set_priors_initials_agrs $ c_prior_beta <- c_prior_betasi
    set_priors_initials_agrs $ d_prior_beta <- d_prior_betasi
    set_priors_initials_agrs $ e_prior_beta <- e_prior_betasi
    set_priors_initials_agrs $ f_prior_beta <- f_prior_betasi
    set_priors_initials_agrs $ g_prior_beta <- g_prior_betasi
    set_priors_initials_agrs $ h_prior_beta <- h_prior_betasi
    set_priors_initials_agrs $ i_prior_beta <- i_prior_betasi
    set_priors_initials_agrs $ s_prior_beta <- s_prior_betasi
    
    set_priors_initials_agrs $ a_cov_prior_beta <- a_cov_prior_betasi
    set_priors_initials_agrs $ b_cov_prior_beta <- b_cov_prior_betasi
    set_priors_initials_agrs $ c_cov_prior_beta <- c_cov_prior_betasi
    set_priors_initials_agrs $ d_cov_prior_beta <- d_cov_prior_betasi
    set_priors_initials_agrs $ e_cov_prior_beta <- e_cov_prior_betasi
    set_priors_initials_agrs $ f_cov_prior_beta <- f_cov_prior_betasi
    set_priors_initials_agrs $ g_cov_prior_beta <- g_cov_prior_betasi
    set_priors_initials_agrs $ h_cov_prior_beta <- h_cov_prior_betasi
    set_priors_initials_agrs $ i_cov_prior_beta <- i_cov_prior_betasi
    set_priors_initials_agrs $ s_cov_prior_beta <- s_cov_prior_betasi
    
    set_priors_initials_agrs $ a_prior_sd <- a_prior_sdsi
    set_priors_initials_agrs $ b_prior_sd <- b_prior_sdsi
    set_priors_initials_agrs $ c_prior_sd <- c_prior_sdsi
    set_priors_initials_agrs $ d_prior_sd <- d_prior_sdsi
    set_priors_initials_agrs $ e_prior_sd <- e_prior_sdsi
    set_priors_initials_agrs $ f_prior_sd <- f_prior_sdsi
    set_priors_initials_agrs $ g_prior_sd <- g_prior_sdsi
    set_priors_initials_agrs $ h_prior_sd <- h_prior_sdsi
    set_priors_initials_agrs $ i_prior_sd <- i_prior_sdsi
    set_priors_initials_agrs $ s_prior_sd <- s_prior_sdsi
    
    set_priors_initials_agrs $ a_cov_prior_sd <- a_cov_prior_sdsi
    set_priors_initials_agrs $ b_cov_prior_sd <- b_cov_prior_sdsi
    set_priors_initials_agrs $ c_cov_prior_sd <- c_cov_prior_sdsi
    set_priors_initials_agrs $ d_cov_prior_sd <- d_cov_prior_sdsi
    set_priors_initials_agrs $ e_cov_prior_sd <- e_cov_prior_sdsi
    set_priors_initials_agrs $ f_cov_prior_sd <- f_cov_prior_sdsi
    set_priors_initials_agrs $ g_cov_prior_sd <- g_cov_prior_sdsi
    set_priors_initials_agrs $ h_cov_prior_sd <- h_cov_prior_sdsi
    set_priors_initials_agrs $ i_cov_prior_sd <- i_cov_prior_sdsi
    set_priors_initials_agrs $ s_cov_prior_sd <- s_cov_prior_sdsi
    
    set_priors_initials_agrs $ gr_prior_cor         <- gr_prior_corsi
    set_priors_initials_agrs $ sigma_prior_cor      <- sigma_prior_corsi
    set_priors_initials_agrs $ sigma_prior_beta     <- sigma_prior_betasi
    set_priors_initials_agrs $ sigma_cov_prior_beta <- sigma_cov_prior_betasi
    
    set_priors_initials_agrs $ sigma_prior_sd      <- sigma_prior_sdsi
    set_priors_initials_agrs $ sigma_cov_prior_sd  <- sigma_cov_prior_sdsi
    set_priors_initials_agrs $ rsd_prior_sigma     <- rsd_prior_sigmasi
    set_priors_initials_agrs $ dpar_prior_sigma    <- dpar_prior_sigmasi
    
  
    set_priors_initials_agrs $ dpar_cov_prior_sigma     <- 
      dpar_cov_prior_sigmasi
    set_priors_initials_agrs $ autocor_prior_acor       <- autocor_prior_acorsi
    set_priors_initials_agrs $ autocor_prior_unstr_acor <- 
      autocor_prior_unstr_acorsi
    set_priors_initials_agrs $ mvr_prior_rescor         <- mvr_prior_rescorsi
    set_priors_initials_agrs $ prior_data               <- prior_data
    set_priors_initials_agrs $ prior_data_internal      <- prior_data_internal
    set_priors_initials_agrs $ prior_args_internal      <- prior_args_internal
    set_priors_initials_agrs $ init_arguments           <- init_arguments
    set_priors_initials_agrs $ init_data                <- init_data
    set_priors_initials_agrs $ init_data_internal       <- init_data_internal
    set_priors_initials_agrs $ init_args_internal       <- init_args_internal
    set_priors_initials_agrs $ custom_order_prior_str   <- ""
    
    set_priors_initials_agrs $ d_as_random_only         <- d_as_random_only
  
    bpriors <- CustomDoCall(set_priors_initials, set_priors_initials_agrs)
    
    stanvar_priors <- attr(bpriors, "stanvars")
    
    initials <- attr(bpriors, "initials")
    
    # check and add hierarchical prior (for 3 level and more)
    # First, sd
    set_class_what <- 'sd'
    set_org_priors_initials_agrs_what <- set_priors_initials_agrs
    set_randomsi_higher_levsl <- strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
    
    check_sigma_str <- 
      eval(parse(text = paste0('sigma', "covcoefnames_gr_str")))
    
    if(!is.null(check_sigma_str)) {
      set_randomsi_higher_levsl <- c(set_randomsi_higher_levsl, 'sigma')
    }
    
    evaluate_higher_level_sd_priors <- function(set_nlpar_, 
                                                set_class,
                                                set_prior,
                                                set_cov_prior,
                                                org_priors_initials_agrs,
                                                set_env,
                                                ...) {
      
      custom_order_prior_str <- c(paste0(set_nlpar_, "_prior_sd"),
                                  paste0(set_nlpar_, "_cov_prior_sd"))
     
      
      eval_what <- eval(parse(text = paste0(set_nlpar_, 
                                            "covcoefnames_gr_str")), 
                        envir = set_env_what)
      if(!is.null(eval_what)) {
        gr_str_id   <- eval(parse(text = paste0(set_nlpar_, 
                                                "covcoefnames_gr_str_id")), 
                            envir = set_env_what)
        gr_str_coef <- eval(parse(text = paste0(set_nlpar_, 
                                                "covcoefnames_gr_str")), 
                            envir = set_env_what)
        gr_str_form <- eval(parse(text = paste0(set_nlpar_, 
                                                "covcoefnames_gr_str_form")), 
                            envir = set_env_what)
        gr_str_ncov <- eval(parse(text = paste0(set_nlpar_, 
                                                "ncov_gr_str")), 
                            envir = set_env_what)
        temp_gr_str_priors <- list()
        temp_gr_str_stanvars <- c()
        temp_gr_str_inits <- c()
        set_priors_initials_agrs_str <- org_priors_initials_agrs 
        # this for adding _prior_cor 
        
       counter_start_from_one_for_prior <- 0
       
       # 24.08.2024
       # now after 24.08.2024 update, 1:length(eval_what) needed, why?
       for (istrx in 1:length(eval_what)) {
        # for (istrx in 2:length(eval_what)) {
          counter_start_from_one_for_prior <- 
            counter_start_from_one_for_prior + 1
          if(set_nlpar_ == 'sigma') {
            assign('sigma_arg_groupvar', gr_str_id[[istrx]], 
                   envir = set_env_what)
          } else {
            assign('group_arg_groupvar', gr_str_id[[istrx]], 
                   envir = set_env_what)
          }
          assign( paste0(set_nlpar_, "_formula_grsi"), 
                  gr_str_form[[istrx]], envir = set_env_what)
          assign( paste0(set_nlpar_, "covcoefnames_gr"), 
                  gr_str_coef[[istrx]], envir = set_env_what)
          assign( paste0(set_nlpar_, "ncov_gr"), 
                  gr_str_coef[[istrx]], envir = set_env_what)
          
          prior_args_internal_str <- list()
          prior_args_internal_str <- mget(prior_args_internal_names, 
                                          envir = set_env_what)
          set_priors_initials_agrs_str $ prior_args_internal <- 
            prior_args_internal_str

          set_priors_initials_agrs_str $ custom_order_prior_str <- 
            custom_order_prior_str
       
          set_priors_initials_agrs_str [[paste0(set_nlpar_, 
                                                "_prior_sd")]]  <- 
            set_prior[counter_start_from_one_for_prior]
          set_priors_initials_agrs_str [[paste0(set_nlpar_, 
                                                "_cov_prior_sd")]] <- 
            set_cov_prior[counter_start_from_one_for_prior]

          bpriors_str <- CustomDoCall(set_priors_initials, 
                                 set_priors_initials_agrs_str, 
                                 envir = set_env_what)

          stanvars_str <- attr(bpriors_str, "stanvars")
          initials_str <- attr(bpriors_str, "initials")
          temp_gr_str_stanvars <- c(temp_gr_str_stanvars, stanvars_str)
          temp_gr_str_priors[[istrx]] <- bpriors_str
        }
        temp_gr_str_priors <- temp_gr_str_priors %>% CustomDoCall(rbind, .)
        out <- list(temp_gr_str_priors = temp_gr_str_priors,
                    temp_gr_str_stanvars = temp_gr_str_stanvars,
                    temp_gr_str_inits = temp_gr_str_inits)
      }
      return(out)
    } 
    
    
    temp_gr_str_priors_sd <- list()
    temp_gr_str_stanvars_sd <-  temp_gr_str_inits_sd <- c()
    for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
      set_nlpar_what <- set_randomsi_higher_levsli
      set_env_what   <- environment()
      
      if(set_nlpar_what == "sigma") {
        n_higher_str <- length(eval(parse(text = 
                                            paste0(set_nlpar_what, "_",
                                                   "hierarchical_gr_names")),
                                    envir = set_env_what))
      } else {
        n_higher_str <- length(eval(parse(text = 
                                            paste0("",
                                                   "hierarchical_gr_names")),
                                    envir = set_env_what))
      }
      n_higher_str   <- n_higher_str - 1
     
      
      if(n_higher_str > 0) {
        set_assign_prior_what <- '_prior'
        check_prior_ifp <- 
          extract_prior_str_lv(ept(paste0(set_nlpar_what, 
                                          paste0(set_assign_prior_what, 
                                                 "_", set_class_what, 
                                                 "_strsi"))))
        check_prior_ifp_true_false <- FALSE
        if(length(check_prior_ifp) == 1) {
          if(!is.null(ept(paste0(set_nlpar_what, 
                                 paste0(set_assign_prior_what, "_", 
                                        set_class_what, "_strsi") ))[[1]][1])) {
            check_prior_ifp_true_false <- TRUE 
          } else {
            check_prior_ifp_true_false <- FALSE
          }
          if(ept(paste0(set_nlpar_what, 
                        paste0(set_assign_prior_what, "_", 
                               set_class_what, "_strsi") )) != "NULL") {
            check_prior_ifp_true_false <- TRUE
          } else {
            check_prior_ifp_true_false <- FALSE
          }
        } else if(length(check_prior_ifp) > 1) {
          check_prior_ifp_true_false <- TRUE
        }
        if(check_prior_ifp_true_false) {
          assign(paste0(set_nlpar_what, paste0(set_assign_prior_what, "_",
                                               set_class_what, "si")),  
                 check_prior_ifp)
        }
        set_prior_what <- ept(paste0(set_nlpar_what, 
                                     paste0(set_assign_prior_what, "_", 
                                            set_class_what, "si") ))
        
        
        paste_message <- paste("Length of prior elements for random effect ",
             "'", set_nlpar_what, "'",
             " \n",
             "  specified by using the argument ", 
             "'", paste0(set_nlpar_what, 
                         paste0(set_assign_prior_what, "_", 
                                set_class_what, "_str")), "'",  " ",
             " \n",
             "  should be one or same as the levels of hierarchy minus one.",
             " \n",
             "  (minus one because the prior for the second level of hierarchy",
             " \n",
             "  is taken from the ", 
             "'", paste0(set_nlpar_what, paste0(set_assign_prior_what, 
                                                "_", set_class_what, "")), "'"
        )
        
        if(length(set_prior_what) > 1 & 
           length(set_prior_what) != n_higher_str) {
          stop2c(paste_message)
        } else if(length(set_prior_what) == 1) {
          set_prior_what <- rep(set_prior_what, n_higher_str)
        }
        paste_message <- NULL

        
        set_assign_prior_what <- '_cov_prior'
        check_prior_ifp <- 
          extract_prior_str_lv(ept(paste0(set_nlpar_what, 
                                          paste0(set_assign_prior_what, 
                                                 "_", set_class_what, 
                                                 "_strsi"))))
        check_prior_ifp_true_false <- FALSE
        if(length(check_prior_ifp) == 1) {
          if(!is.null(ept(paste0(set_nlpar_what, 
                                 paste0(set_assign_prior_what, "_", 
                                        set_class_what, "_strsi") ))[[1]][1])) {
            check_prior_ifp_true_false <- TRUE 
          } else {
            check_prior_ifp_true_false <- FALSE
          }
          if(ept(paste0(set_nlpar_what, 
                        paste0(set_assign_prior_what, "_", 
                               set_class_what, "_strsi") )) != "NULL") {
            check_prior_ifp_true_false <- TRUE
          } else {
            check_prior_ifp_true_false <- FALSE
          }
        } else if(length(check_prior_ifp) > 1) {
          check_prior_ifp_true_false <- TRUE
        }
        if(check_prior_ifp_true_false) {
          assign(paste0(set_nlpar_what, paste0(set_assign_prior_what, "_",
                                               set_class_what, "si")),  
                 check_prior_ifp)
        }
        set_cov_prior_what <- ept(paste0(set_nlpar_what, 
                                         paste0(set_assign_prior_what, "_", 
                                                set_class_what, "si") ))
        
        if(length(set_cov_prior_what) > 1 & 
           length(set_cov_prior_what) != n_higher_str) {
          stop2c("Length of prior elements for random effect parameter ",
               "'", set_nlpar_what, "'",
               " \n",
               "  specified by using the argument ", 
               "'", paste0(set_nlpar_what, 
                           paste0(set_assign_prior_what, "_", 
                                  set_class_what, "_str")), "'",  " ",
               " \n",
               "  should be one or same as the levels of hierarchy minus one.",
               " \n",
               "  (minus one because prior for the first level of hierarchy",
               " \n",
               "  is taken from the ", 
               "'", paste0(set_nlpar_what, 
                           paste0(set_assign_prior_what, "_", 
                                  set_class_what, "")), "'"
          )
        } else if(length(set_prior_what) == 1) {
          set_cov_prior_what <- rep(set_cov_prior_what, n_higher_str)
        }
       
       
        out2 <- evaluate_higher_level_sd_priors(set_nlpar_ = set_nlpar_what, 
                                        set_class  = set_class_what,
                                        set_prior = set_prior_what,
                                        set_cov_prior = set_cov_prior_what,
                                        set_env = set_env_what,
                                        org_priors_initials_agrs = 
                                          set_org_priors_initials_agrs_what)
      
      temp_gr_str_priors_sd[[set_randomsi_higher_levsli]] <- 
        out2 $ temp_gr_str_priors
      temp_gr_str_stanvars_sd <- 
        c(temp_gr_str_stanvars_sd, out2 $ temp_gr_str_stanvars)
      temp_gr_str_inits_sd <- 
        c(temp_gr_str_inits_sd,    out2 $ temp_gr_str_inits)
      } 
    } 
    
    higher_level_priors <- temp_gr_str_priors_sd %>% CustomDoCall(rbind, .)
    bpriors             <- rbind(bpriors, higher_level_priors)
    
    if(length(temp_gr_str_stanvars_sd) > 0) {
      stanvar_priors_c <- temp_gr_str_stanvars_sd_c <- c()
      for (i in 1:length(stanvar_priors)) {
        stanvar_priors_c <- c(stanvar_priors_c, stanvar_priors[i])
      }
      for (i in 1:length(temp_gr_str_stanvars_sd)) {
        temp_gr_str_stanvars_sd_c <- c(temp_gr_str_stanvars_sd_c, 
                                       temp_gr_str_stanvars_sd[i])
      }
      stanvar_priors <- c(stanvar_priors_c, temp_gr_str_stanvars_sd_c)
    } 
    
    if(length(temp_gr_str_inits_sd) > 0) {
      initials_c <- temp_gr_str_inits_sd_c <- c()
      for (i in 1:length(initials)) {
        initials_c <- c(initials_c, initials[i])
      }
      for (i in 1:length(temp_gr_str_inits_sd)) {
        temp_gr_str_inits_sd_c <- c(temp_gr_str_inits_sd_c, 
                                    temp_gr_str_inits_sd[i])
      }
      initials <- c(initials_c, temp_gr_str_inits_sd) 
    } 
    
    # Now, cor priors    
    # Adding cor priors is tricky because of complex |x| formulations
    set_class_what <- 'cor'
    set_org_priors_initials_agrs_what <- set_priors_initials_agrs
    set_randomsi_higher_levsl <- 'gr'
    check_sigma_str <- eval(parse(text = paste0('sigma', 
                                                "covcoefnames_gr_str")))
    if(!is.null(check_sigma_str)) {
      set_randomsi_higher_levsl <- c(set_randomsi_higher_levsl, 'sigma')
    }
    
    
    evaluate_higher_level_corr_priors <- function(set_nlpar_, 
                                                set_class,
                                                set_prior,
                                                id_higher_str = id_higher_str,
                                                corr_higher_str_tf = 
                                                  corr_higher_str_tf,
                                                org_priors_initials_agrs,
                                                set_env,
                                                ...) {
      
      custom_order_prior_str <- c(paste0(set_nlpar_, "_prior_cor"))
      
      temp_gr_str_priors <- list()
      temp_gr_str_stanvars <- c()
      temp_gr_str_inits <- c()

      set_priors_initials_agrs_str <- org_priors_initials_agrs 
      
      gr_str_id <- id_higher_str
      
      counter_start_from_one_for_prior <- 0
      
      # 24.08.2024
      # Somehow now after 24.08.2024, 2:length(eval_what) needed, why?
      
      for (istrx in 1:length(gr_str_id)) {
      # for (istrx in 2:length(gr_str_id)) {
        counter_start_from_one_for_prior <- 
          counter_start_from_one_for_prior + 1
        get_corr_higher_str_tf <- corr_higher_str_tf[istrx]
                  if(get_corr_higher_str_tf) {
          if(set_nlpar_ == 'sigma') {
            assign('sigma_arg_groupvar', gr_str_id[istrx], 
                   envir = set_env_what)
            set_priors_initials_agrs_str $ sigma_prior_cor <- 
              set_prior[counter_start_from_one_for_prior]
          } else {
            assign('group_arg_groupvar', gr_str_id[istrx], 
                   envir = set_env_what)
            set_priors_initials_agrs_str $ gr_prior_cor <- 
              set_prior[counter_start_from_one_for_prior]
          }
          prior_args_internal_str <- list()
          prior_args_internal_str <- mget(prior_args_internal_names, 
                                          envir = set_env_what)
          set_priors_initials_agrs_str $ prior_args_internal <- 
            prior_args_internal_str
          
          set_priors_initials_agrs_str $ custom_order_prior_str <-
            custom_order_prior_str
          
          bpriors_str <- CustomDoCall(set_priors_initials, 
                                 set_priors_initials_agrs_str, 
                                 envir = set_env_what)
          stanvars_str <- attr(bpriors_str, "stanvars")
          initials_str <- attr(bpriors_str, "initials")
          temp_gr_str_stanvars <- c(temp_gr_str_stanvars, stanvars_str)
          temp_gr_str_inits    <- c(temp_gr_str_inits, initials_str)
          temp_gr_str_priors[[istrx]] <- bpriors_str
          
          # 26 12 2023
          bpriors_str_checks <- bpriors_str
          attr(bpriors_str_checks, "stanvars") <- NULL
          attr(bpriors_str_checks, "initials") <- NULL
          bpriors_str_checks <- bpriors_str_checks
          attributes(bpriors_str_checks) <- NULL
          if(!is.list(bpriors_str_checks)) {
            if(bpriors_str_checks == "") {
              temp_gr_str_priors[[istrx]] <- temp_gr_str_stanvars <- NULL
              temp_gr_str_inits <- NULL
            }
          }
        } # if(get_corr_higher_str_tf) {
        if(!get_corr_higher_str_tf) {
          temp_gr_str_priors[[istrx]] <- temp_gr_str_stanvars <- NULL
            temp_gr_str_inits <- NULL
        }
        
      } 
      
      temp_gr_str_priors <- temp_gr_str_priors %>% CustomDoCall(rbind, .)
      out <- list(temp_gr_str_priors = temp_gr_str_priors,
                  temp_gr_str_stanvars = temp_gr_str_stanvars,
                  temp_gr_str_inits = temp_gr_str_inits)
      
    return(out)
    } # end evaluate_higher_level_corr_priors
    

    temp_gr_str_priors_corr <- list()
    temp_gr_str_stanvars_corr <-  temp_gr_str_inits_corr <- c()

    for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
      set_nlpar_what <- set_randomsi_higher_levsli
      set_env_what   <- environment()
      id_higher_str  <- eval(parse(text = paste0(set_nlpar_what, 
                                                 "_str_unique_id")), 
                             envir = set_env_what)
      
      n_higher_str   <- length(id_higher_str)
      if(set_nlpar_what == "sigma") {
        n_higher_str <- 
          length(eval(parse(text = paste0(set_nlpar_what, "_",
                                          "hierarchical_gr_names")),
                                    envir = set_env_what))
      } else {
        n_higher_str <- 
          length(eval(parse(text = paste0("",
                                          "hierarchical_gr_names")),
                                    envir = set_env_what))
      }
      n_higher_str   <- n_higher_str - 1
      corr_higher_str_tf <- eval(parse(text = paste0(set_nlpar_what, 
                                                     "_str_corr_tf")),
                                 envir = set_env_what)
      if(n_higher_str > 0) {
        set_assign_prior_what <- '_prior'
        check_prior_ifp <- 
          extract_prior_str_lv(ept(paste0(set_nlpar_what, 
                                          paste0(set_assign_prior_what, 
                                                 "_", set_class_what, 
                                                 "_strsi"))))
        check_prior_ifp_true_false <- FALSE
        if(length(check_prior_ifp) == 1) {
          if(!is.null(ept(paste0(set_nlpar_what, 
                                 paste0(set_assign_prior_what, "_", 
                                        set_class_what, "_strsi") ))[[1]][1])) {
            check_prior_ifp_true_false <- TRUE 
          } else {
            check_prior_ifp_true_false <- FALSE
          }
          if(ept(paste0(set_nlpar_what, 
                        paste0(set_assign_prior_what, "_", 
                               set_class_what, "_strsi") )) != "NULL") {
            check_prior_ifp_true_false <- TRUE
          } else {
            check_prior_ifp_true_false <- FALSE
          }
        } else if(length(check_prior_ifp) > 1) {
          check_prior_ifp_true_false <- TRUE
        }
        if(check_prior_ifp_true_false) {
          assign(paste0(set_nlpar_what, paste0(set_assign_prior_what, "_",
                                               set_class_what, "si")),  
                 check_prior_ifp)
        }
        set_prior_cor_what <- ept(paste0(set_nlpar_what, 
                                     paste0(set_assign_prior_what, "_", 
                                            set_class_what, "si") ))
        
        if(length(set_prior_cor_what) > 1 & 
           length(set_prior_cor_what) != n_higher_str) {
          stop2c("Length of prior elements for random effect parameter ",
               "'", set_prior_cor_what, "'",
               " \n",
               "  specified by using the argument ", 
               "'", paste0(set_prior_cor_what, paste0(set_assign_prior_what, 
                                                      "_", 
                                                      set_class_what, "_str")), 
               "'",  " ",
               " \n",
               "  should be one or same as the levels of hierarchy minus one.",
               " \n",
               "  (minus one because prior for the second level of hierarchy",
               " \n",
               "  is taken from the ", 
               "'", paste0(set_prior_cor_what, 
                           paste0(set_assign_prior_what, "_", 
                                  set_class_what, "")), "'"
          )
        } else if(length(set_prior_cor_what) == 1) {
          set_prior_cor_what <- rep(set_prior_cor_what, n_higher_str)
        }
       
        out2 <- evaluate_higher_level_corr_priors(
          set_nlpar_ = set_nlpar_what, 
          set_class  = set_class_what,
          set_prior = set_prior_cor_what,
          id_higher_str = id_higher_str,
          corr_higher_str_tf = corr_higher_str_tf,
          set_env = set_env_what,
          org_priors_initials_agrs = set_org_priors_initials_agrs_what)
       
        temp_gr_str_priors_corr[[set_randomsi_higher_levsli]] <- 
          out2 $ temp_gr_str_priors
        temp_gr_str_stanvars_corr <- 
          c(temp_gr_str_stanvars_corr, out2 $ temp_gr_str_stanvars)
        temp_gr_str_inits_corr <- 
          c(temp_gr_str_inits_corr,    out2 $ temp_gr_str_inits)
      } 
    } # for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
    
    higher_level_priors_corr <- temp_gr_str_priors_corr %>% 
      CustomDoCall(rbind, .)
    bpriors <- rbind(bpriors, higher_level_priors_corr)
   
    if(length(temp_gr_str_stanvars_corr) > 0) {
      stanvar_priors_c <- temp_gr_str_stanvars_corr_c <- c()
      for (i in 1:length(stanvar_priors)) {
        stanvar_priors_c <- c(stanvar_priors_c, stanvar_priors[i])
      }
      for (i in 1:length(temp_gr_str_stanvars_corr)) {
        temp_gr_str_stanvars_corr_c <- c(temp_gr_str_stanvars_corr_c, 
                                         temp_gr_str_stanvars_corr[i])
      }
      stanvar_priors <- c(stanvar_priors_c, temp_gr_str_stanvars_corr_c)
    } 
    
    if(length(temp_gr_str_inits_corr) > 0) {
      initials_c <- temp_gr_str_inits_corr_c <- c()
      for (i in 1:length(initials)) {
        initials_c <- c(initials_c, initials[i])
      }
      for (i in 1:length(temp_gr_str_inits_corr)) {
        temp_gr_str_inits_corr_c <- c(temp_gr_str_inits_corr_c, 
                                      temp_gr_str_inits_corr[i])
      }
      initials <- c(initials_c, temp_gr_str_inits_corr) 
    } 
   
    
    priorlist <- rbind(priorlist, bpriors)
    stanvar_priors_names <- names(stanvar_priors)
    
    if(!"stanvars" %in% attr(stanvar_priors, 'class')) {
      attr(stanvar_priors, 'class') <- c("stanvars", 
                                         attr(stanvar_priors, 'class'))
    }
    
    prior_stanvarlist[[ii]] <- stanvar_priors 
    
    scode_auxillary <- attr(bpriors, "scode_auxillary")
    auxillary_stanvarlist[[ii]] <- scode_auxillary
    
    initialslist[[ii]]   <- initials
    initialslist_s[[ii]] <- initsi
    

    
    #################################################################
    #################################################################
    
    xvar_name                    <- "xvar"
    yvar_name                    <- "yvar"
    idvar_name                   <- "idvar"
    sigmaidvar_name              <- "sigmaidvar"
    cov_name                     <- "cov"
    sigmacov_name                <- "sigmacov"
    xfun_name                    <- "xfun"
    yfun_name                    <- "yfun"
    sigmaxfun_name               <- "sigmaxfun"
    xfuntransform_name           <- "xfuntransform"
    ixfuntransform_name          <- "ixfuntransform"
    xfuntransform2_name          <- "xfuntransform2"
    ixfuntransform2_name         <- "ixfuntransform2"
    yfuntransform_name           <- "yfuntransform"
    iyfuntransform_name          <- "iyfuntransform"
    sigmaxfuntransform_name      <- "sigmaxfuntransform"
    sigmaixfuntransform_name     <- "sigmaixfuntransform"
    sigmaxfuntransform2_name     <- "sigmaxfuntransform2"
    sigmaixfuntransform2_name    <- "sigmaixfuntransform2"
    xoffset_name                 <- "xoffset"
    knots_name                   <- "knots"
    nknots_name                  <- "nknots"
    fixed_name                   <- "fixed"
    random_name                  <- "random"
    groupvar_name                <- "groupvar"
    sigma_groupvar_name          <- "sigma_groupvar"
    hierarchical_name            <- "hierarchical"
    sigma_hierarchical_name      <- "sigma_hierarchical"
    d_adjusted_name              <- "d_adjusted"
    sigmafixed_name              <- "sigmafixed"
    sigmarandom_name             <- "sigmarandom"
    sigmad_adjusted_name         <- "sigmad_adjusted"
    sigmaxoffset_name            <- "sigmaxoffset"
    sigmaxvar_name               <- "sigmaxvar"
    setsigmaxvar_name            <- "setsigmaxvar"
    
    sigmabasicfunname_name      <- "sigmabasicfunname"
    sigmabasicfunattr_name      <- "sigmabasicfunattr"
    sigmamodelname_name         <- "sigmamodel"
    
    # these paste0(..., 's') will be combined across ys
    # these need to be defined only once
    # Corresponding values 'xfun_names_val' ... are defined outside the loop
    if(ii == 1) {
      xvar_names                <- paste0(xvar_name,                "s")
      yvar_names                <- paste0(yvar_name,                "s")
      idvar_names               <- paste0(idvar_name,               "s")
      sigmaidvar_names          <- paste0(sigmaidvar_name,          "s")
      sigmaxvar_names           <- paste0(sigmaxvar_name,          "s")
      cov_names                 <- paste0(cov_name,                 "s")
      sigmacov_names            <- paste0(sigmacov_name,            "s")
      xfun_names                <- paste0(xfun_name,                "s")
      yfun_names                <- paste0(yfun_name,                "s")
      sigmaxfun_names           <- paste0(sigmaxfun_name,           "s")
      xfuntransform_names       <- paste0(xfuntransform_name,       "s")
      ixfuntransform_names      <- paste0(ixfuntransform_name,      "s")
      xfuntransform2_names      <- paste0(xfuntransform2_name,      "s")
      ixfuntransform2_names     <- paste0(ixfuntransform2_name,     "s")
      yfuntransform_names       <- paste0(yfuntransform_name,       "s")
      iyfuntransform_names      <- paste0(iyfuntransform_name,      "s")
      sigmaxfuntransform_names  <- paste0(sigmaxfuntransform_name,  "s")
      sigmaixfuntransform_names <- paste0(sigmaixfuntransform_name, "s")
      sigmaxfuntransform2_names <- paste0(sigmaxfuntransform2_name, "s")
      sigmaixfuntransform2_names<- paste0(sigmaixfuntransform2_name,"s")
      xoffset_names             <- paste0(xoffset_name,             "s")
      sigmaxoffset_names        <- paste0(sigmaxoffset_name,        "s")
      setsigmaxvar_names        <- paste0(setsigmaxvar_name,        "s")
      sigmabasicfunname_names   <- paste0(sigmabasicfunname_name,  "s")
      sigmabasicfunattr_names   <- paste0(sigmabasicfunattr_name,  "s")
      sigmamodelname_names      <- paste0(sigmamodelname_name,  "s")
    }
 
    if (nys > 1) {
      xvar_name                 <- paste0(xvar_name,                 "_", ysi)
      yvar_name                 <- paste0(yvar_name,                 "_", ysi)
      idvar_name                <- paste0(idvar_name,                "_", ysi)
      sigmaidvar_name           <- paste0(sigmaidvar_name,           "_", ysi)
      cov_name                  <- paste0(cov_name,                  "_", ysi)
      sigmacov_name             <- paste0(sigmacov_name,             "_", ysi)
      xfun_name                 <- paste0(xfun_name,                 "_", ysi)
      yfun_name                 <- paste0(yfun_name,                 "_", ysi)
      sigmaxfun_name            <- paste0(sigmaxfun_name,            "_", ysi)
      xfuntransform_name        <- paste0(xfuntransform_name,        "_", ysi)
      ixfuntransform_name       <- paste0(ixfuntransform_name,       "_", ysi)
      xfuntransform2_name       <- paste0(xfuntransform2_name,       "_", ysi)
      ixfuntransform2_name      <- paste0(ixfuntransform2_name,      "_", ysi)
      yfuntransform_name        <- paste0(yfuntransform_name,        "_", ysi)
      iyfuntransform_name       <- paste0(iyfuntransform_name,       "_", ysi)
      sigmaxfuntransform_name   <- paste0(sigmaxfuntransform_name,   "_", ysi)
      sigmaixfuntransform_name  <- paste0(sigmaixfuntransform_name,  "_", ysi)
      sigmaxfuntransform2_name  <- paste0(sigmaxfuntransform2_name,  "_", ysi)
      sigmaixfuntransform2_name <- paste0(sigmaixfuntransform2_name, "_", ysi)
      xoffset_name              <- paste0(xoffset_name,              "_", ysi)
      knots_name                <- paste0(knots_name,                "_", ysi)
      nknots_name               <- paste0(nknots_name,               "_", ysi)
      fixed_name                <- paste0(fixed_name,                "_", ysi)
      random_name               <- paste0(random_name,               "_", ysi)
      groupvar_name             <- paste0(groupvar_name,             "_", ysi)
      sigma_groupvar_name       <- paste0(sigma_groupvar_name,       "_", ysi)
      hierarchical_name         <- paste0(hierarchical_name,         "_", ysi)
      sigma_hierarchical_name   <- paste0(sigma_hierarchical_name,   "_", ysi)
      d_adjusted_name           <- paste0(d_adjusted_name,           "_", ysi)
      sigmafixed_name           <- paste0(sigmafixed_name,           "_", ysi)
      sigmarandom_name          <- paste0(sigmarandom_name,          "_", ysi)
      sigmad_adjusted_name      <- paste0(sigmad_adjusted_name,      "_", ysi)
      sigmaxoffset_name         <- paste0(sigmaxoffset_name,         "_", ysi)
      sigmaxvar_name            <- paste0(sigmaxvar_name,            "_", ysi)
      setsigmaxvar_name         <- paste0(setsigmaxvar_name,         "_", ysi)
      sigmabasicfunname_name    <- paste0(sigmabasicfunname_name,    "_", ysi)
      sigmabasicfunattr_name    <- paste0(sigmabasicfunattr_name,    "_", ysi)
      sigmamodelname_name       <- paste0(sigmamodelname_name,       "_", ysi)
    }
    
    xnamelist[[ii]]                        <- xvar_name
    xvarvaluelist[[ii]]                    <- xsi
    ynamelist[[ii]]                        <- yvar_name
    yvarvaluelist[[ii]]                    <- ysi
    idnamelist[[ii]]                       <- idvar_name
    idvarvaluelist[[ii]]                   <- idsi
    sigmaidnamelist[[ii]]                  <- sigmaidvar_name
    sigmaidvarvaluelist[[ii]]              <- sigmaidsi
    covnamelist[[ii]]                      <- cov_name
    covvaluelist[[ii]]                     <- covariates_
    knotsnamelist[[ii]]                    <- knots_name
    knotsvaluelist[[ii]]                   <- knots
    xfunvaluelist[[ii]]                    <- xfunvalue
    xfunnamelist[[ii]]                     <- xfun_name
    yfunvaluelist[[ii]]                    <- yfunvalue
    yfunnamelist[[ii]]                     <- yfun_name
    sigmaxfunvaluelist[[ii]]               <- sigmaxfunvalue
    sigmaxfunnamelist[[ii]]                <- sigmaxfun_name
    xfuntransformvaluelist[[ii]]           <- xfuntransformsi
    xfuntransformnamelist[[ii]]            <- xfuntransform_name
    ixfuntransformvaluelist[[ii]]          <- ixfuntransformsi
    ixfuntransformnamelist[[ii]]           <- ixfuntransform_name
    xfuntransform2valuelist[[ii]]          <- xfuntransform2si
    xfuntransform2namelist[[ii]]           <- xfuntransform2_name
    ixfuntransform2valuelist[[ii]]         <- ixfuntransform2si
    ixfuntransform2namelist[[ii]]          <- ixfuntransform2_name
    yfuntransformvaluelist[[ii]]           <- yfuntransformsi
    yfuntransformnamelist[[ii]]            <- yfuntransform_name
    iyfuntransformvaluelist[[ii]]          <- iyfuntransformsi
    iyfuntransformnamelist[[ii]]           <- iyfuntransform_name
    sigmaxfuntransformvaluelist[[ii]]      <- sigmaxfuntransformsi
    sigmaxfuntransformnamelist[[ii]]       <- sigmaxfuntransform_name
    sigmaxfuntransform2valuelist[[ii]]     <- sigmaxfuntransform2si
    sigmaxfuntransform2namelist[[ii]]      <- sigmaxfuntransform2_name
    sigmaixfuntransformvaluelist[[ii]]     <- sigmaixfuntransformsi
    sigmaixfuntransformnamelist[[ii]]      <- sigmaixfuntransform_name
    sigmaixfuntransform2valuelist[[ii]]    <- sigmaixfuntransform2si
    sigmaixfuntransform2namelist[[ii]]     <- sigmaixfuntransform2_name
    include_fun_nameslist_rnamelist[[ii]]  <- includefunnameslistname
    include_fun_nameslist_rvaluelist[[ii]] <- unlist(include_fun_nameslist)
    funlist_rnamelist[[ii]]                <- funlist_r_name
    funlist_rvaluelist[[ii]]               <- unlist(funlist_r)
    sigmafunlist_rnamelist[[ii]]           <- sigmafunlist_r_name
    sigmafunlist_rvaluelist[[ii]]          <- unlist(sigmafunlist_r)
    
    sigmavarfunlist_rnamelist[[ii]]        <- sigmavarfunlist_r_name
    sigmavarfunlist_rvaluelist[[ii]]       <- unlist(sigmavarfunlist_r)
    
    sigmabasicfunlist_rnamelist[[ii]]      <- sigmabasicfunlist_r_name
    sigmabasicfunlist_rvaluelist[[ii]]     <- unlist(sigmabasicfunlist_r)
    
    xoffsetnamelist[[ii]]                  <- xoffset_name
    xoffsetvaluelist[[ii]]                 <- xoffset # xoffsetsi
    sigmaxoffsetnamelist[[ii]]             <- sigmaxoffset_name
    sigmaxoffsetvaluelist[[ii]]            <- sigmaxoffset # sigmaxoffsetsi
    groupvarnamelist[[ii]]                 <- groupvar_name
    groupvarvaluelist[[ii]]                <- group_arg_groupvar
    sigma_groupvarnamelist[[ii]]           <- sigma_groupvar_name
    sigma_groupvarvaluelist[[ii]]          <- sigma_arg_groupvar
    hierarchicalvarnamelist[[ii]]          <- hierarchical_name
    hierarchicalvarvaluelist[[ii]]         <- hierarchical_gr_names
    sigma_hierarchicalvarnamelist[[ii]]    <- sigma_hierarchical_name
    sigma_hierarchicalvarvaluelist[[ii]]   <- sigma_hierarchical_gr_names
    sigmacovnamelist[[ii]]                 <- sigmacov_name
    sigmacovvaluelist[[ii]]                <- sigmacovariates_
    d_adjustednamelist[[ii]]               <- d_adjusted_name
    d_adjustedvaluelist[[ii]]              <- ept(d_adjustedsi)
    sigmad_adjustednamelist[[ii]]          <- sigmad_adjusted_name
    sigmad_adjustedvaluelist[[ii]]         <- ept(sigmad_adjustedsi)
    sigmaxnamelist[[ii]]                   <- sigmaxvar_name
    sigmaxvarvaluelist[[ii]]               <- sigmaxsi
    
    setsigmaxvarnamelist[[ii]]             <- setsigmaxvar_name
    setsigmaxvarvaluelist[[ii]]            <- setsigmaxvarsi
    
    
    sigmabasicfunnamenamelist[[ii]]        <- sigmabasicfunname_name
    sigmabasicfunnamevaluelist[[ii]]       <- sigmabasicfunnamesi
    
    sigmabasicfunattrnamelist[[ii]]        <- sigmabasicfunattr_name
    sigmabasicfunattrvaluelist[[ii]]       <- sigmabasicfunattrsi
    
    sigmamodelnamenamelist[[ii]]           <- sigmamodelname_name
    sigmamodelnamevaluelist[[ii]]          <- sigmamodelsi
    
    fixednamelist[[ii]]        <- fixed_name
    fixedvaluelist[[ii]]       <- abc_fixedsi
    
    sigmafixednamelist[[ii]]   <- sigmafixed_name
    sigmafixedvaluelist[[ii]]  <- strsplit(gsub("\\+", " ", sigmafixedsi), 
                                           " ")[[1]]
    
    randomnamelist[[ii]]       <- random_name
    randomvaluelist[[ii]]      <- abc_randomsi
    
    sigmarandomnamelist[[ii]]  <- sigmarandom_name
    sigmarandomvaluelist[[ii]] <- strsplit(gsub("\\+", " ",  sigmarandomsi), 
                                           " ")[[1]]
    
    
    #####################################################################
    #####################################################################
    # add data scode
    x_Naux_str           <- 'Naux'
    if (nys == 1) {
      N_name    <- 'N'
      Naux_name <- x_Naux_str
      xsi_name  <- xsi
      ysi_name  <- ysi
    } else if (nys > 1) {
      Naux_name <- paste0(x_Naux_str, "_", ysi)
      N_name    <- paste0('N', "_", ysi)
      xsi_name  <- paste0(xsi, "_", ysi)
      ysi_name  <- paste0(ysi, "_", ysi)
    }
    
    x_xoffset          <- xoffset
    name_xoffset_name  <- xoffset_name
    scode_xoffset_name <- paste0("real ", name_xoffset_name, ";")
    
    x_nknots          <- nknots
    name_nknots_name  <- nknots_name
    scode_nknots_name <- paste0("int ", name_nknots_name, ";")
    
    x_knots          <- knots
    name_knots_name  <- knots_name
    scode_knots_name <- paste0("vector[", name_nknots_name, "] ", 
                               name_knots_name, ";")
    
    
    data_stanvarlist[[ii]] <- 
      brms::stanvar(x = x_xoffset,
                    name = name_xoffset_name,
                    scode = scode_xoffset_name,
                    block = "data",
                    position = "start",
                    pll_args = NULL) +
      brms::stanvar(x = x_nknots,
                    name = name_nknots_name,
                    scode = scode_nknots_name,
                    block = "data",
                    position = "start",
                    pll_args = NULL) +
      brms::stanvar(x = x_knots,
                    name = name_knots_name,
                    scode = scode_knots_name,
                    block = "data",
                    position = "start",
                    pll_args = NULL) 
    
    #####################################################################
    #####################################################################
    
    # restoring original data
    # Just before leaving the loop, restore all inverse transformations
    prepare_transformations_args[['data']]         <- datai
    prepare_transformations_args[['xvar']]         <- xsi
    prepare_transformations_args[['yvar']]         <- ysi
    prepare_transformations_args[['sigmaxvar']]    <- sigmaxsi
    prepare_transformations_args[['xfun']]         <- xfuntransformsi
    prepare_transformations_args[['yfun']]         <- yfuntransformsi
    prepare_transformations_args[['sigmaxfun']]    <- sigmaxfuntransformsi
    prepare_transformations_args[['ixfun']]        <- TRUE
    prepare_transformations_args[['iyfun']]        <- TRUE
    prepare_transformations_args[['sigmaixfun']]   <- TRUE
    prepare_transformations_args[['xoffset']]      <- xoffset
    prepare_transformations_args[['sigmaxoffset']] <- sigmaxoffset
    prepare_transformations_args[['transform']]    <- ""
    prepare_transformations_args[['itransform']]   <- ""
    
    datai <- CustomDoCall(prepare_transformations, prepare_transformations_args)
   
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA"))
      dataout <- rbind(dataout, datai)
    else
      dataout <- datai

    if (!(is.na(univariate_by$by) | univariate_by$by == "NA"))
      uvarbyTF <- TRUE
    else
      uvarbyTF <- FALSE
  }  # End of the loop
  # End of Start loop over response i.e. ii ...
  
  #######################################################################
  #######################################################################
  
  
  if (verbose) {
    if (multivariate$mvar) {
      setmsgtxt <-
        paste0(
          "\n Combining formula, function, priors and initials",
          "\n ",
          "for multivariate model fitting"
        )
    }
    if (!is.na(univariate_by$by)) {
      setmsgtxt <-
        paste0(
          "\n Combining formula, function, priors and initials",
          "\n ",
          "for univariate-by-subgroup model fitting"
        )
    }
    if (displayit == 'msg') {
      message2c(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  #######################################################################
  #######################################################################
  
  xvar_names_val                <- xvarvaluelist %>% unlist()
  yvar_names_val                <- yvarvaluelist %>% unlist()
  idvar_names_val               <- idvarvaluelist %>% unlist()
  sigmaidvar_names_val          <- sigmaidvarvaluelist %>% unlist()
  cov_names_val                 <- covvaluelist %>% unlist()
  xfun_names_val                <- xfunvaluelist %>% unlist()
  yfun_names_val                <- yfunvaluelist %>% unlist()
  xfuntransform_names_val       <- xfuntransformvaluelist %>% unlist()
  xfuntransform2_names_val      <- xfuntransform2valuelist %>% unlist()
  yfuntransform_names_val       <- yfuntransformvaluelist %>% unlist()
  ixfuntransform_names_val      <- ixfuntransformvaluelist %>% unlist()
  ixfuntransform2_names_val     <- ixfuntransform2valuelist %>% unlist()
  iyfuntransform_names_val      <- iyfuntransformvaluelist %>% unlist()
  xoffset_names_val             <- xoffsetvaluelist %>% unlist()
  
  # unlke bewlo sigmas infom this will be ratined always
  setsigmaxvar_names_val        <- setsigmaxvarvaluelist %>% unlist()
  
  sigmaxvar_names_val           <- sigmaxvarvaluelist %>% unlist()
  sigmacov_names_val            <- sigmacovvaluelist %>% unlist()
  sigmaxfun_names_val           <- sigmaxfunvaluelist %>% unlist()
  sigmaxfuntransform_names_val  <- sigmaxfuntransformvaluelist %>% unlist()
  sigmaxfuntransform2_names_val <- sigmaxfuntransform2valuelist %>% unlist()
  sigmaixfuntransform_names_val <- sigmaixfuntransformvaluelist %>% unlist()
  sigmaixfuntransform2_names_val<- sigmaixfuntransform2valuelist %>% unlist()
  sigmaxoffset_names_val        <- sigmaxoffsetvaluelist %>% unlist()
  
 
  
  #######################################################################
  #######################################################################
  
  dataout_restoted.org.in <- dataout
  
  # Now we are out of loop, re transform 
  prepare_transformations_args[['data']]         <- dataout
  prepare_transformations_args[['xvar']]         <- xvarvaluelist %>% unlist()
  prepare_transformations_args[['yvar']]         <- yvarvaluelist %>% unlist()
  prepare_transformations_args[['sigmaxvar']]    <- 
    sigmaxvarvaluelist %>% unlist()
  prepare_transformations_args[['xfun']]         <- 
    xfuntransformvaluelist %>% unlist()
  prepare_transformations_args[['yfun']]         <- 
    yfuntransformvaluelist %>% unlist()
  prepare_transformations_args[['sigmaxfun']]    <- 
    sigmaxfuntransformvaluelist %>% unlist()
  prepare_transformations_args[['ixfun']]        <- FALSE
  prepare_transformations_args[['iyfun']]        <- FALSE
  prepare_transformations_args[['sigmaixfun']]   <- FALSE
  prepare_transformations_args[['xoffset']]      <- xoffsetvaluelist
  prepare_transformations_args[['sigmaxoffset']] <- sigmaxoffsetvaluelist
  prepare_transformations_args[['transform']]    <- ""
  prepare_transformations_args[['itransform']]   <- ""
  
  dataout <- CustomDoCall(prepare_transformations, prepare_transformations_args)
  
  
  #######################################################################
  #######################################################################
  
  # assemble 'bformula'
  bflist_c_list <- list()
  bflist_c <- c()
  for (il in 1:length(bflist)) {
    bflist_c_list[[il]] <- ept(bflist[[il]])
    bflist_c <- c(bflist_c, paste0("bflist_c_list[[", il, "]]"))
  }
  bformula <- ept(paste(bflist_c, collapse = "+"))
  
  if (nys > 1) {
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
      bformula <- bformula + brms::set_rescor(FALSE)
    }
    if (multivariate$mvar && multivariate$rescor) {
      bformula <- bformula + brms::set_rescor(TRUE)
    }
    if (multivariate$mvar && !multivariate$rescor) {
      bformula <- bformula + brms::set_rescor(FALSE)
    }
  }
  
  ###################################################################
  ###################################################################
  # assembe brmsdata and brmspriors
  brmsdata   <- dataout
  brmspriors <- priorlist
  # IMP - brms does not allow different lb for sd params (e.e, all to be NA)
  # Error: Conflicting boundary information for coefficients of class 'sd'.
  # Because prior function automatically sets lb 0 for positive priors 
  # such as exponential the following is need (again done at line 4753 )
  
  brmspriors <- brmspriors %>% 
    dplyr::mutate(lb = dplyr::if_else(class == 'sd', NA, lb))
  brmspriors <- brmspriors %>% 
    dplyr::mutate(ub = dplyr::if_else(class == 'sd', NA, ub))
  
  
  #######################################################################
  #######################################################################
  # some checks for sigma var modelling - except _ba / _ls / _no
  #######################################################################
  #######################################################################
  
  if(set_model_sigma_by_fz |
     set_model_sigma_by_fp |
     set_model_sigma_by_fe |
     set_model_sigma_by_ve | 
     set_model_sigma_by_vp |
     set_model_sigma_by_cp |
     set_model_sigma_by_mp |
     set_model_sigma_by_me |
     set_model_sigma_by_rp |
     set_model_sigma_by_re ) {
  
    # check - 1
    ithx <- 0
    for (outrespbames in ys) {
      ithx <- ithx + 1
      if(nys == 1) {
        base_forms   <- bformula
        base_pforms  <- base_forms[['pforms']]
        sigma_forms  <- base_pforms[['sigma']]
      } else if(nys > 1) {
        base_forms   <- bformula[['forms']][[outrespbames]]
        base_pforms  <- base_forms[['pforms']]
        sigma_forms  <- base_pforms[['sigma']]
      }
      sigma_forms_  <- sigma_forms[[3]]
      sigma_forms_vars <- all.vars(sigma_forms_)
      sigmatau_strsi_i <- sigmatau_strsi_c[[ithx]]
      all_form_pars <- base_pforms %>% names()
      all_forms_vars_ <- setdiff(sigmatau_strsi_i, all_form_pars)
      sigma_forms_vars_not_allowed <- setdiff(all_forms_vars_,
                                              colnames(brmsdata))
      if(length(sigma_forms_vars_not_allowed) > 0) {
        stop2c("You have included the following variable(s) in the formula ", 
             "\n  ", 
             "which are are not present in the data: ",
             "\n ", 
             collapse_comma(sigma_forms_vars_not_allowed),
             "\n  ",
             ". Please check the relevant part of the code shown below: ",
             "\n ", 
             paste0(gsub_space(deparse(sigma_forms)), collapse = ""))
      }
    }
  
    
    # check - 2
    ithx <- 0
    for (outrespbames in ys) {
      ithx <- ithx + 1
      sigmatau_strsi_i <- sigmatau_strsi_c[[ithx]]
      if(set_model_sigma_by_fz |
         set_model_sigma_by_fp | 
         set_model_sigma_by_fe | 
         set_model_sigma_by_ve | 
         set_model_sigma_by_vp |
         set_model_sigma_by_mp |
         set_model_sigma_by_me |
         set_model_sigma_by_rp |
         set_model_sigma_by_re ) {
        length_of_sigma_var_nlpars <- 2
      } else if(set_model_sigma_by_cp) {
        length_of_sigma_var_nlpars <- 3
      }
      if(length(sigmatau_strsi_i) != length_of_sigma_var_nlpars) {
        stop2c("For sigma variance method ", 
               collapse_comma(nlf_sigma_method_arg), 
             ", the number of nlpar should be ", 
             length_of_sigma_var_nlpars,
             "\n ",
             " but the number of specified nlpar is ", length(sigmatau_strsi_i),
             "\n ",
             collapse_comma(sigmatau_strsi_i)
        )
      }
    }
    
    # check - 3
    if(nys > 1) {
      # check for _ls
      if(length(unique(unlist(sigmaspfncname_c))) > 1) {
        stop2c("The name of 'sigma' functions defined for 
               modelling location-scale",
             "\n model should be same across all responses.",
             "\n Currently specified names are: ", 
             collapse_comma(unique(unlist(sigmaspfncname_c))),
             "\n Also note that function name should have 'sigma' as prefix",
             "\n without any underscore such as 'sigmafun'. In case prefix",
             "\n 'sigma' is missing, this will be added internally i.e., if ",
             "\n function name is 'fun', then it will be renamed as 'sigmafun'"
        )
      }
      # check for _mu
      if(! all_inner_lengths_equal_in_list(sigmatau_strsi_c) ) {
        stop2c("The number of 'nlpar' parameters defined for 
               modelling sigma as a",
             "\n  function of mean should be same across all responses.",
             "\n  (the response specific renaming is done internally)",
             "\n  Currently specified names are: ", 
             collapse_comma(sigmatau_strsi_c)
        )
      }
      
      if(!all_elements_identical_in_list(sigmatau_strsi_c)) {
        stop2c("The names of 'nlpar' parameters defined for modelling sigma as a",
             "\n  function of mean should be same across all responses.",
             "\n  (the response specific renaming is done internally)",
             "\n  Currently specified names are: ", 
             collapse_comma(unique(unlist(sigmatau_strsi_c)))
        )
      }
    } # if(nys > 1) {
    
    
    
    # check - 4 - checks and update brmsprior 
    if(sigma_formula_manual_prior_via_sigma_formula & 
       !is.null(sigmatau_strsi)) {
      warn_sigma_self_prior_msg <- 
        paste0(" There are custom nlpar parameters for the distributional",
               "\n  ",
               "parameters sigma (", 
               collapse_comma(sigmatau_strsi), ")",
               "\n  ",
               "For each of these nlpar parameter, no prior is set",
               "\n  ",
               "as you have set prior=self in the nlf() function call.",
               "\n  ",
               "Thus, defaut priors set by 'brms' will be used automatically",
               "\n  ",
               "If you want to use priors genearted by the 'bsitar()' for the",
               "\n  ",
               "distributional parameters sigma, then please remove prior=self",
               "\n  ",
               "from the nlf(). In that case, same priors will be set for each",
               "\n  ",
               "nlpar parameter across all outcomes, In case you stick",
               "\n  ",
               "to use the  nlf(..., prior=self), then you can supply custom",
               "\n  ",
               "priors that will be added to the prior object",
               "\n  ",
               "These priors can be set using 'add_self_priors = x'",
               "\n  ",
               "call where x is the prior object with your custom priors")
      
      user_prompt_msg <- 
        paste0(" There are more than one nlpar parameters for the  ",
               "\n ",
               "distributional parameters sigma (", 
               collapse_comma(sigmatau_strsi), ")",
               "\n ",
               "For each of these nlpar parameters, same priors will be set.",
               "\n ",
               "If you want to use self defined priors, then exit and use ",
               "\n ",
               "prior=self in the nlf() function such as nlf(..., prior=self)",
               "\n ",
               "In that case, use 'add_self_priors = xx' in bsitar() call",
               "\n ",
               "where xx is the prior object with your custom priors") 
    }
    
    
    # check - 5 update brmspriors
    priorobject <- brmspriors
    
    check_prompt <- FALSE
    check_verbose <- verbose
    # add_sigma_by_mu
    if(sigma_formula_manual_prior_via_sigma_formula & 
       !is.null(sigmatau_strsi)) {
      set_user_prompt <- check_prompt
      if(length(sigmatau_strsi) > 1) {
        if(set_user_prompt) {
          message2c(user_prompt_msg)
          if (user_prompt()) {
            #
          } else {
            return(invisible(NULL))
          }
        } # if(set_user_prompt) {
      } # if(length(sigmatau_strsi) > 1) {
      priorobject_no_sigam   <- priorobject %>% dplyr::filter(dpar != 'sigma')
      if(nlf_sigma_prior_arg != 'self') {
        priorobject_only_sigam <- priorobject %>% dplyr::filter(dpar == 'sigma')
        priorobject_only_sigam_sigmatau_strsi_c <- list()
        counter_sigmatau_strsi <- 0
        for (i in sigmatau_strsi) {
          counter_sigmatau_strsi <- counter_sigmatau_strsi + 1
          if(nlf_sigma_method_arg == "fitted") {
            if(counter_sigmatau_strsi == length(sigmatau_strsi)) {
              # wanted but 
              #Prior argument 'coef' may not be specified when using boundaries.
              # set_lb = "0"
              set_lb = ""
            } else {
              set_lb = ""
            }
          } else {
            set_lb = ""
          }
          priorobject_sigmatau_strsi <- 
            priorobject_only_sigam %>% 
            dplyr::mutate(nlpar = dplyr::if_else(dpar == "sigma", 
                                                 i, nlpar)) %>% 
            dplyr::mutate(class = dplyr::if_else(dpar == "sigma" & 
                                                   class != 'sd', "b", 
                                                 class)) %>% 
            dplyr::mutate(dpar = dplyr::if_else(dpar == "sigma", "", dpar)) %>% 
            dplyr::mutate(lb = dplyr::if_else(nlpar == i, set_lb, lb))
          
          priorobject_only_sigam_sigmatau_strsi_c[[i]] <- 
            priorobject_sigmatau_strsi
        }
        priorobject <- rbind(priorobject_no_sigam, 
                             do.call(rbind, 
                                     priorobject_only_sigam_sigmatau_strsi_c))
      } else if(nlf_sigma_prior_arg == 'self') {
        priorobject <- priorobject_no_sigam
        if(check_verbose) {
          message2c(warn_sigma_self_prior_msg)
        }
      } # else if(nlf_sigma_prior_arg == 'self') {
    } # if(sigma_formula_manual_prior_via_sigma_formula
    
    
    brmspriors <- priorobject
  } # if(set_model_sigma_by_fz |.....
  
  #######################################################################
  # end some checks for sigma var modelling - except _ba / _ls / _no
  #######################################################################
  
  
  
  #######################################################################
  #######################################################################
  # some checks for sigma var modelling - for _ba / _ls / _no
  #######################################################################
  #######################################################################

  if(set_model_sigma_by_ba) {
    # check - 5 update brmspriors
    priorobject <- brmspriors
    check_prompt <- FALSE
    check_verbose <- verbose
    # add_sigma_by_mu
    if(sigma_formula_manual_prior_via_sigma_formula & !is.null(sigmatau_strsi)) {
      set_user_prompt <- check_prompt
      if(length(sigmatau_strsi) > 1) {
        if(set_user_prompt) {
          message2c(user_prompt_msg)
          if (user_prompt()) {
            #
          } else {
            return(invisible(NULL))
          }
        } # if(set_user_prompt) {
      } # if(length(sigmatau_strsi) > 1) {
      priorobject_no_sigam   <- priorobject %>% dplyr::filter(dpar != 'sigma')
      if(nlf_sigma_prior_arg != 'self') {
        priorobject_only_sigam <- priorobject %>% dplyr::filter(dpar == 'sigma')
        priorobject_only_sigam_sigmatau_strsi_c <- list()
        counter_sigmatau_strsi <- 0
        for (i in sigmatau_strsi) {
          counter_sigmatau_strsi <- counter_sigmatau_strsi + 1
          if(nlf_sigma_method_arg == "fitted") {
            if(counter_sigmatau_strsi == length(sigmatau_strsi)) {
              # wanted but 
              # Prior argument 'coef' may not be specified when using boundaries.
              # set_lb = "0"
              set_lb = ""
            } else {
              set_lb = ""
            }
          } else {
            set_lb = ""
          }
          priorobject_sigmatau_strsi <- 
            priorobject_only_sigam %>% 
            dplyr::mutate(nlpar = dplyr::if_else(dpar == "sigma", 
                                                 i, nlpar)) %>% 
            dplyr::mutate(class = dplyr::if_else(dpar == "sigma" & 
                                                   class != 'sd', "b", 
                                                 class)) %>% 
            dplyr::mutate(dpar = dplyr::if_else(dpar == "sigma", "", dpar)) %>% 
            dplyr::mutate(lb = dplyr::if_else(nlpar == i, set_lb, lb))
          
          priorobject_only_sigam_sigmatau_strsi_c[[i]] <- 
            priorobject_sigmatau_strsi
        }
        priorobject <- rbind(priorobject_no_sigam, 
                             do.call(rbind, 
                                     priorobject_only_sigam_sigmatau_strsi_c))
      } else if(nlf_sigma_prior_arg == 'self') {
        priorobject <- priorobject_no_sigam
        if(check_verbose) {
          message2c(warn_sigma_self_prior_msg)
        }
      } # else if(nlf_sigma_prior_arg == 'self') {
    } # if(sigma_formula_manual_prior_via_sigma_formula
    
    brmspriors <- priorobject
  } # if(set_model_sigma_by_ba) {
  
  #######################################################################
  #######################################################################
  # end some checks for sigma var modelling - for _ba / _ls / _no
  #######################################################################
  #######################################################################
  
  
  
  ###################################################################
  ###################################################################
  # assemble stan functions code
  
  fun_scode <- paste(funlist, collapse = "\n")
  
  for (j in 1:length(setsigmaxvar_names_val)) {
    if(setsigmaxvar_names_val[j]) {
      fun_scode <- paste(fun_scode, sigmafunlist[j], collapse = "\n")
    }
  }
   # sigmavar
  if(length(sigmavarfunlist) > 0) {
    for (j in 1:length(sigmavarfunlist)) {
        fun_scode <- paste(fun_scode, sigmavarfunlist[[j]], collapse = "\n")
    }
  }
  
  # brms::stanvar code without functions {} block enclosure
  bstanvars <- brms::stanvar(scode = fun_scode, block = "function")
  
  # Now add functions { to fun_scode will be used in expose_model_function
  fun_scode <- paste0("functions {", "\n", fun_scode, "\n", "}")
  
  ###################################################################
  ###################################################################
  
  if (length(data_stanvarlist) != 0) {
    data_stanvarlistlist <- c()
    for (i in 1:nys) {
      data_stanvarlistlist[i] <- paste0("data_stanvarlist[[", i, "]]")
    }
    bstanvars <-
      bstanvars + eval(parse(text = paste(data_stanvarlistlist, 
                                          collapse = "+")))
  }
  
  prior_stanvarlistlist <- c()
  for (i in 1:nys) {
    prior_stanvarlistlist[i] <- paste0("prior_stanvarlist[[", i, "]]")
  }
  
  bstanvars <-
    bstanvars + eval(parse(text = paste(prior_stanvarlistlist, 
                                        collapse = "+")))
  
  
  if (length(auxillary_stanvarlist) != 0) {
    auxillary_stanvarlistlist <- c()
    for (i in 1:nys) {
      auxillary_stanvarlistlist[i] <-
        paste0("auxillary_stanvarlist[[", i, "]]")
    }
    bstanvars <-
      bstanvars + eval(parse(text = paste(
        auxillary_stanvarlistlist, collapse = "+"
      )))
  }
  
  ############################################################################
  
  # add_rescor_by
  Rescor_by_levels <- NULL
  if (set_rescor_by) {
      if(!is.null(multivariate$rcorr_by)) {
        
        if(!is.null(multivariate$rcorr_gr)) {
          Rescor_gr_id  <- multivariate$rcorr_gr 
        } else {
          if(verbose) message2c("Rescor_gr_id set as group_arg$groupvar")
          Rescor_gr_id <- group_arg$groupvar # 'id'
        }
        Rescor_gr_id_integer     <- as.integer(brmsdata[[Rescor_gr_id]])
        
        Rescor_by_id   <- multivariate$rcorr_by
        
        if(!is.factor(brmsdata[[Rescor_by_id]])) {
          stop2c("The variable ", "'", multivariate$rcorr_by, "'",
               " set as 'rcorr_by' must be a factor variable")
        } else if(is.factor(brmsdata[[Rescor_by_id]])) {
          if(nlevels(brmsdata[[Rescor_by_id]]) == 1) {
            stop2c("The variable ", "'", multivariate$rcorr_by, "'",
                 " set as 'rcorr_by' must be a factor variable",
                 "\n  ",
                 "with at least two levels")
          }
        }
        
        Rescor_by_levels <- levels(brmsdata[[Rescor_by_id]])
        Rescor_by_levels <- paste0(Rescor_by_id, Rescor_by_levels)
        
        Rescor_by_id_integer <- brmsdata %>% 
          dplyr:: group_by(!!as.name(Rescor_gr_id)) %>%
          dplyr:: filter(dplyr::row_number() == 1) %>% 
          dplyr::pull(Rescor_by_id) %>% 
          as.vector() 
        
        Rescor_by_id_integer     <- as.integer(as.factor(Rescor_by_id_integer))
        Rescor_by_id_integer_max <- max(Rescor_by_id_integer)
        
        if(!is.null(multivariate$rcorr_method)) {
          rcorr_method_choices <- c('lkj', 'cde')
          if(!multivariate$rcorr_method %in% rcorr_method_choices) {
            stop2c("The residual correletion method specified as ", "'",
                 multivariate$rcorr_method,"'",
                 " is invalid",
                 "\n  ", 
                 "The available residual correletion methods are:",
                 "\n  ", 
                 "'lkj': models residual correletions using LKJ prior",
                 "\n  ", 
                 "'cde': models residual correletions using the Cholesky",
                 "\n   ", 
                 "decomposition of covariance matrix. In this method, uniform",
                 "\n   ", 
                 "priors (-1,1) are assigned to each correlation parameter"
                 )
          }
          Rescor_method  <- multivariate$rcorr_method 
        } else {
          Rescor_method <- 'lkj'
        }
        
        if(!is.null(multivariate$rcorr_prior)) {
          Rescor_prior  <- multivariate$rcorr_prior 
          if(length(Rescor_prior) == 1) {
            Rescor_prior <- rep(Rescor_prior, Rescor_by_id_integer_max) %>% 
              as.vector()
          } else {
            if(length(Rescor_prior) != Rescor_by_id_integer_max) {
              stop2c("lenght of 'rcorr_prior' must be either 1 or same",
                   "as levels of 'rescor_by'")
            }
          }
        } else {
          Rescor_prior <- rep(1, Rescor_by_id_integer_max) %>% as.vector()
        }     
        # create stanvars 
        Rescor_by_stanvars <- 
          brms::stanvar(x = Rescor_by_id_integer_max, 
                        name = 'Rescor_Nby',
                        scode = "int Rescor_Nby;", 
                        block = 'data') + 
          brms::stanvar(x = Rescor_gr_id_integer, 
                        name = 'Rescor_gr_id',
                        scode = "array[N] int<lower=1, upper=N_1> Rescor_gr_id;", 
                        block = 'data') + 
          brms::stanvar(x = Rescor_by_id_integer, 
                        name = 'Rescor_by_id',
                        scode = "array[N_1] int<lower=1, upper=Rescor_Nby> Rescor_by_id;",  
                        block = 'data') 
        if(Rescor_method == 'lkj') {
          Rescor_by_stanvars <- Rescor_by_stanvars + 
            brms::stanvar(x = Rescor_prior, 
                          name = 'Rescor_prior',
                          scode = "vector[Rescor_Nby] Rescor_prior;",  
                          block = 'data') 
        }
        if(Rescor_method == 'cde') {
          Rescor_by_stanvars <- Rescor_by_stanvars + 
            brms::stanvar(x = (nys * (nys - 1) / 2), 
                          name = 'N_rhos',
                          scode = "int N_rhos;",  
                          block = 'data') 
        }
        bstanvars <- bstanvars + Rescor_by_stanvars
        # end create stanvars 
      } # if(!is.null(multivariate$rescor_by)) {
  } # if (set_rescor_by) {
  

  
  # rescor_by rescor_gr rescor_method rescor_lkj   multivariate$mvar
  if (is.list(initialslist) & length(initialslist) == 0) {
    brmsinits <- NULL
  } else if (is.list(initialslist) & length(initialslist) > 0) {
    clistlist <- c()
    for (i in 1:length(initialslist)) {
      clistlist <- c(clistlist, ept(paste0("initialslist[[", i, "]]")))
    }
    brmsinits <- clistlist
  }
  
  
  if (!is.null(brmsinits)) {
    if (multivariate$mvar & multivariate$cor == "un") {
     
      c_it <- "sd_"
      brmsinits_names <- names(brmsinits)
      brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', 
                                                brmsinits_names)]
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      temppp <- unlist(unname(temppp))
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      brmsinits[[keys[1]]] <- temppp %>%
        unname()
      
      c_it <- "L_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      t_names <- l_comb <- d_comb <- c()
      for (lnamei in names(temppp)) {
        t <- temppp[[lnamei]]
        l <- t[lower.tri(t)]
        names(l) <-
          apply(combn(colnames(t), 2), 2, paste, collapse = "_")
        d_comb <- c(d_comb, ncol(t))
        l_comb <- c(l_comb, l)
        t_names <- c(t_names, colnames(t))
      }
      
      create_cor_mat <- function(n, cor = NULL) {
        n_elements <- n
        m <- diag(n_elements)
        m_upper <- m_lower <- matrix(0, n_elements, n_elements)
        nc <- n_elements * (n_elements - 1) / 2
        if (is.null(cor)) {
          x <- rep(0, nc)
        } else {
          x <- cor
          if (length(x) != nc) {
            stop2c("length of correlation vector must be ",
                 nc,
                 "\n, ",
                 ", but found ",
                 length(x))
          }
        }
        m_lower[lower.tri(m_lower, diag = FALSE)] <- x
        m_upper <- t(m_lower)
        M <- m_lower + m + m_upper
        return(M)
      }
      
     
      if(!is.null(t_names)) { # 17.02.2025
        tt_names <- apply(combn(t_names, 2), 2, paste, collapse = "_")
        tt_dims <- sum(d_comb)
        tt_nc <- (tt_dims * (tt_dims - 1) / 2)
        tt_12 <- create_cor_mat(tt_dims, rep(0, tt_nc))
        colnames(tt_12) <- rownames(tt_12) <- t_names
        tt_ll <- tt_12[lower.tri(tt_12)]
        names(tt_ll) <- apply(combn(colnames(tt_12), 2), 2, paste, 
                              collapse = "_")
        tt_ll[names(l_comb)] <- l_comb
        tt_ll[!names(tt_ll) %in% names(l_comb)] <- 0
        brmsinits[[keys[1]]] <- create_cor_mat(tt_dims, tt_ll)
      } # 17.02.2025
      
      c_it <- "z_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      brmsinits[[keys[1]]] <- CustomDoCall(rbind, temppp)
    } else if (multivariate$mvar &
               (multivariate$cor == "un" | multivariate$cor ==
                "un_s") &
               !any(grepl("^L_", names(brmsinits)))) {
     
      c_it <- "sd_"
      brmsinits_names <- names(brmsinits)
      brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', 
                                                brmsinits_names)]
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      temppp <- unlist(unname(temppp))
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      for (sdi in 1:length(temppp)) {
        brmsinits[[paste0(c_it, sdi)]] <- temppp[sdi] %>% unname()
      }
    }  
    
    
    # keep only one Lrescor
    if (multivariate$mvar & multivariate$rescor) {
      c_it <- "Lrescor_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      if(!is_emptyx(temppp)) { # 17.02.2025
        brmsinits[["Lrescor"]] <- temppp[[1]]
      } # 17.02.2025
    }
    

    # add_rescor_by
    if (set_rescor_by) {
      matrix_to_repated_as_array <- brmsinits[["Lrescor"]]
      # if 'matrix_to_repated_as_array == NULL', it means random init values
      if(!is.null(matrix_to_repated_as_array)) {
        array_of_matrices_array <- array(
          data = rep(matrix_to_repated_as_array, Rescor_by_id_integer_max),
          dim = c(nrow(matrix_to_repated_as_array),
                  ncol(matrix_to_repated_as_array),
                  Rescor_by_id_integer_max))
        # This is how rstan restructures Lrescor initials
        reordered_array <- aperm(array_of_matrices_array, perm = c(3, 2, 1))
        brmsinits[["Lrescor"]] <- reordered_array 
      } # if(!is.null(matrix_to_repated_as_array)) {
    } # if (set_rescor_by) {
    
    

    if ((multivariate$mvar & multivariate$cor == "diagonal") |
        (!is.na(univariate_by$by) & 
         univariate_by$cor == "diagonal") |
        group_arg$cor == "diagonal" |
        sigma_group_arg$cor == "diagonal") {
     
      c_it <- "sd_"
      brmsinits_names <- names(brmsinits)
      brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', 
                                                brmsinits_names)]
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      temppp <- unlist(unname(temppp))
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      xxx <- temppp
      ilc <- list()
      ilc_c <- 0
      for (nysi_ in 1:nys) {
        for (il in letters[1:26]) {
          ilc_c <- ilc_c + 1
          na <- paste0("^", il, nysi_)
          nb <- paste0("^", il, "cov", nysi_)
          nanb <- paste0(na, "|", nb)
          if (length(xxx[grepl(nanb, names(xxx))]) > 0) {
            ilc[[ilc_c]] <- xxx[grepl(nanb, names(xxx))]
          }
        }
      }
      
      if(!is_emptyx(ilc)) {
        ilc <- ilc[lengths(ilc) != 0]
        names(ilc) <- paste0("sd_", 1:length(ilc))
      }
      
      
      for (sdi in names(ilc)) {
        brmsinits[[sdi]] <- ilc[[sdi]]
      }
      
      c_it <- "z_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      
      if(!is_emptyx(ilc)) {
        for (zi in 1:length(ilc)) {
          brmsinits[[paste0(c_it, zi)]] <- temppp[[zi]]
        }
      }
      #
    }
  }  
  
  
  # For multivariate, it makes sense to keep initials for betas only otherwise
  # dimensional mismatch
  if (!is.null(brmsinits) & length(initialslist) != nys) {
    if (multivariate$mvar & multivariate$cor == "un") {
      c_it_names <- c("sd_", "L_", "z_", "Lrescor")
      for (c_it in c_it_names) {
        brmsinits_names <- names(brmsinits)
        brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', 
                                                  brmsinits_names)]
        keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
        brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      }
    }
  }
  
  
  if (all(sapply("random", grepl, initialslist_s))) {
    brmsinits <- "random"
    brmsinits_r <- ept(init_rsi)
    brmsinits_ <- NULL
  } else if (!all(sapply(NULL, grepl, ept(initialslist_s)  ))) {
    # } else if (all(sapply("0", grepl, initialslist_s))) {
    if (all(sapply(0, grepl, ept(initialslist_s)  ))) {
      brmsinits_r <- ept(init_rsi)
      brmsinits <- "0"
      brmsinits_ <- NULL
    }
  } else {
    brmsinits <- brmsinits
    brmsinits_r <- ept(init_rsi) # NULL
    brmsinits_ <- ""
  }
  


  check_set_init_r <- FALSE # new
  if(initialslist_s[[1]][1] == "NULL") { # new
    brmsinits <- brmsinits
  } else if (is.null(ept(initialslist_s)[[1]][1]) | # new else if
             ept(initialslist_s) == "NULL") {
    brmsinits <- brmsinits
  } else {
    if(all(brmsinits != "random")) {
      brmsinits <- brmsinits_r <- ept(initialslist_s)[[1]][1]
    } else if(all(brmsinits == "random")) {
      brmsinits <- brmsinits_r <- NULL
      check_set_init_r <- TRUE
    }
    brmsinits_ <- NULL
  }
  
 
  
  
  # New to set init_r for 'cmdstanr' and 'rstan' when init = random
  if(check_set_init_r) {
    if(is.character(init_rsi)) {
      brmsinits_r <- ept(init_rsi)
      if(!is.null(brmsinits_r)) {
        brmsinits <- brmsinits_r
      }
    } else if(!is.character(init_rsi)) {
      if(!is.null(brmsinits_r)) {
        brmsinits <- brmsinits_r
      }
    }
  } # if(check_set_init_r) {

  
  if(!is.null(brmsinits)) {
    if(!is.list(brmsinits)) {
      if(brmsinits == 0) brmsinits <- '0'
    }
  }
  
  
  # print(brmsinits)
  # print(brmsinits_r)
  # stop()
 
 
  for (inm in names(brmsinits)) {
    if (is.matrix(brmsinits[[inm]])) {
      colnames(brmsinits[[inm]]) <- rownames(brmsinits[[inm]]) <- NULL
      t__ <- brmsinits[[inm]]
      if (!is.null(attr(t__, "dimnames")))
        attr(t__, "dimnames") <- NULL
      brmsinits[[inm]] <- t__
    }
    if (is.vector(brmsinits[[inm]])) {
      t__ <- brmsinits[[inm]]
      if (!is.null(attr(t__, "names")))
        attr(t__, "names") <- NULL
      brmsinits[[inm]] <- t__
    }
  }
  
  
  
  if (!is.null(brmsinits_)) {
    eval_inits_fun <-
      function(inits,
               jitter_init_beta,
               jitter_init_sd,
               jitter_init_cor,
               digits) {
        
        #####################################################################
        # Define get_jitter_list, jitter_x and jitter_mat Functions
        #####################################################################
        get_jitter_list <- function(what_to_jitter) {
          if(is.null(what_to_jitter)) {
            what_to_jitter      <- what_to_jitter
            what_to_jitter_list <- NULL
          } else if(!is.null(what_to_jitter)) {
            if(is.list(what_to_jitter)) {
              what_to_jitter_list <- list()
              if(is.null(what_to_jitter[['factor']])) {
                what_to_jitter_list[['factor']] <- 1
              } else {
                what_to_jitter_list[['factor']] <- what_to_jitter[['factor']]
              }
              if(is.null(what_to_jitter[['amount']])) {
                what_to_jitter_list[['amount']] <- NULL
              } else {
                what_to_jitter_list[['amount']] <- what_to_jitter[['amount']]
              }
              if(is.null(what_to_jitter[['percent']])) {
                what_to_jitter_list[['percent']] <- NULL
              } else {
                what_to_jitter_list[['percent']] <- what_to_jitter[['percent']]
              }
            } else if(!is.list(what_to_jitter)) {
              if(is.numeric(what_to_jitter)) {
                # This will set up jitter_amount <- abs(x[i]) * a
                what_to_jitter_list <- list()
                what_to_jitter_list[['factor']] <- 1
                what_to_jitter_list[['amount']] <- NULL 
                what_to_jitter_list[['percent']] <- what_to_jitter 
              } else if(!is.numeric(what_to_jitter)) {
                stop2c("The 'what_to_jitter' argument must be either NULL,  
                         a names list, or a numeric value defining percentrage")
              }
            } # if(is.list(what_to_jitter)){else if(!is.list(what_to_jitter)){
          } # if(is.null(what_to_jitter)) {else if(!is.null(what_to_jitter)) {
          return(what_to_jitter_list)
        }
        
        jitter_x <- function(x, what_to_jitter_list, digits) {
          set_jitter_factor  <- what_to_jitter_list[['factor']]
          set_jitter_amount  <- what_to_jitter_list[['amount']]
          set_jitter_percent <- what_to_jitter_list[['percent']]
          if(!is.null(set_jitter_amount) & !is.null(set_jitter_percent)) {
            stop2c("Please specify either amount or percent for jitter, 
                   not both")
          }
          if(is.null(set_jitter_percent)) {
            set_jitter_prop <- NULL
          } else if(!is.null(set_jitter_percent)) {
            set_jitter_prop <- set_jitter_percent / 100
          }
          
          if(is.null(set_jitter_factor)) {
            jitter_factor <- 1 
          } else if(!is.null(set_jitter_factor)) {
            jitter_factor <- set_jitter_factor
          }
          x <- unname(x)
          col <- c()
          for (i in 1:length(x)) {
            if(!is.null(set_jitter_amount)) {
              jitter_amount <- set_jitter_amount 
            } else if(!is.null(set_jitter_prop)) {
              jitter_amount <- abs(x[i]) * set_jitter_prop
            } else {
              jitter_amount <- NULL # default, as in base::jitter
            }
            col <- c(col, jitter(x[i], factor = jitter_factor, 
                                 amount = jitter_amount))
          }
          col <- round(col, digits)
          return(col)
        }
        
        jitter_mat <- function(x, what_to_jitter_list, digits) {
          set_jitter_factor  <- what_to_jitter_list[['factor']]
          set_jitter_amount  <- what_to_jitter_list[['amount']]
          set_jitter_percent <- what_to_jitter_list[['percent']]
          if(!is.null(set_jitter_amount) & !is.null(set_jitter_percent)) {
            stop2c("Please specify either amount or percent for jitter,
                   not both")
          }
          if(is.null(set_jitter_percent)) {
            set_jitter_prop <- NULL
          } else if(!is.null(set_jitter_percent)) {
            set_jitter_prop <- set_jitter_percent / 100
          }
          
          if(is.null(set_jitter_factor)) {
            jitter_factor <- 1 
          } else if(!is.null(set_jitter_factor)) {
            jitter_factor <- set_jitter_factor
          }
          mat_out <- x
          x <- x[lower.tri(x)]
          col <- c()
          for (i in 1:length(x)) {
            if(!is.null(set_jitter_amount)) {
              jitter_amount <- set_jitter_amount 
            } else if(!is.null(set_jitter_prop)) {
              jitter_amount <- abs(x[i]) * set_jitter_prop
            } else {
              jitter_amount <- NULL # default, as in base::jitter
            }
            col <- c(col, jitter(x[i], factor = jitter_factor, 
                                 amount = jitter_amount))
          }
          col <- round(col, digits)
          col <- ifelse(col > 1, 1, col)
          col <- ifelse(col < -1, 1, col)
          mat_out[lower.tri(mat_out)] <-
            mat_out[upper.tri(mat_out)] <- col
          return(mat_out)
        }
        
        
        #####################################################################
        # End Define get_jitter_list, jitter_x and jitter_mat Functions
        #####################################################################
        
        if (is.character(jitter_init_beta)) {
          jitter_init_beta <- ept(jitter_init_beta)
        }
        if (is.character(jitter_init_sd)) {
          jitter_init_sd <- ept(jitter_init_sd)
        }
        if (is.character(jitter_init_cor)) {
          jitter_init_cor <- ept(jitter_init_cor)
        }
          
        jitter_init_beta_list <- get_jitter_list(jitter_init_beta)
        jitter_init_sd_list   <- get_jitter_list(jitter_init_sd)
        jitter_init_cor_list  <- get_jitter_list(jitter_init_cor)
        eval_inits <- c()
        for (i_init in names(inits)) {
          if (grepl("^b_", i_init)) {
            if (!is.null(jitter_init_beta)) {
              values_i <-
                jitter_x(x = inits[[i_init]], 
                         what_to_jitter_list = jitter_init_beta_list, 
                         digits = digits)
            } else if (is.null(jitter_init_beta)) {
              values_i <- inits[[i_init]]
              values_i <- round(values_i, digits)
            }
            eval_inits[[i_init]] <- values_i
          } else if (grepl("^sd_", i_init)) {
            if (!is.null(jitter_init_sd)) {
              values_i <- jitter_x(x = inits[[i_init]], 
                                   what_to_jitter_list = jitter_init_sd_list, 
                                   digits)
              values_i <- abs(values_i)
              values_i <-
                ifelse(values_i <= 0, values_i + 0.01, values_i)
            } else if (is.null(jitter_init_sd)) {
              values_i <- inits[[i_init]]
              values_i <- abs(round(values_i, digits))
              values_i <-
                ifelse(values_i <= 0, values_i + 0.01, values_i)
            }
            eval_inits[[i_init]] <- values_i
          } else if (grepl("^L_", i_init)) {
            if (!is.null(jitter_init_cor)) {
              values_i <- jitter_mat(inits[[i_init]], jitter_init_cor, digits)
            } else if (is.null(jitter_init_cor)) {
              values_i <- inits[[i_init]]
              values_i <- round(values_i, digits)
            }
            eval_inits[[i_init]] <- values_i
          } else {
            eval_inits[[i_init]] <- inits[[i_init]]
          }
        }  # for(i_init in names(inits)) {
        eval_inits
        return(eval_inits)
      }
    
    
    if(is.null(set_self_priors)) {
      temp_prior <- brmspriors
    } else if(!is.null(set_self_priors)) {
      temp_prior <- set_self_priors
    }
    
    # 20.03.2025 - moved to 'final_scode'
    # but added support for returning tempriorstr if get_priors == "default"
    if(!is.logical(get_priors)) {
      if(get_priors == "default") {
        tempriorstr <- brms::get_prior(formula = bformula,
                                       stanvars = bstanvars,
                                       prior = temp_prior,
                                       data = brmsdata)
        return(tempriorstr)
      }
    } # if(!is.logical(get_priors)) {
    
    ################################################################
    
    for (j in 1:length(setsigmaxvar_names_val)) {
      if(any(setsigmaxvar_names_val[j])) {
        temp_prior <- temp_prior %>% dplyr::filter(class != "sigma")
      }
    }
    
    # add_sigma_by_mu
    if(sigma_formula_manual_prior_via_sigma_formula & 
       !is.null(sigmatau_strsi)) {
      set_user_prompt <- FALSE # already prompted at level of brmsprior
      if(length(sigmatau_strsi) > 1) {
        if(set_user_prompt) {
          message2c(user_prompt_msg)
          if (user_prompt()) {
            #
          } else {
            return(invisible(NULL))
          }
        } # if(set_user_prompt) {
      } # if(length(sigmatau_strsi) > 1) {
      priorobject <- temp_prior
      priorobject_no_sigam   <- priorobject %>% dplyr::filter(dpar != 'sigma')
      if(nlf_sigma_prior_arg != 'self') {
        priorobject_only_sigam <- priorobject %>% dplyr::filter(dpar == 'sigma')
        priorobject_only_sigam_sigmatau_strsi_c <- list()
        counter_sigmatau_strsi <- 0
        for (i in sigmatau_strsi) {
          counter_sigmatau_strsi <- counter_sigmatau_strsi + 1
          if(nlf_sigma_method_arg == "fitted") {
            if(counter_sigmatau_strsi == length(sigmatau_strsi)) {
              # wanted but 
              # Prior argument 'coef' may not be specified when using boundaries.
              # set_lb = "0"
            } else {
              set_lb = ""
            }
          } else {
            set_lb = ""
          }
          priorobject_sigmatau_strsi <- 
            priorobject_only_sigam %>% 
            dplyr::mutate(nlpar = dplyr::if_else(dpar == "sigma", 
                                                 i, nlpar)) %>% 
            dplyr::mutate(class = dplyr::if_else(dpar == "sigma" & 
                                                   class != 'sd', "b", 
                                                 class)) %>% 
            dplyr::mutate(dpar = dplyr::if_else(dpar == "sigma", "", dpar)) %>% 
            dplyr::mutate(lb = dplyr::if_else(nlpar == i, set_lb, lb))
          
          priorobject_only_sigam_sigmatau_strsi_c[[i]] <- 
            priorobject_sigmatau_strsi
        }
        priorobject <- rbind(priorobject_no_sigam, 
                             do.call(rbind, 
                                     priorobject_only_sigam_sigmatau_strsi_c))
      } else if(nlf_sigma_prior_arg == 'self') {
        priorobject <- priorobject_no_sigam
        if(check_verbose) {
          message2c(warn_sigma_self_prior_msg)
        }
      } # else if(nlf_sigma_prior_arg == 'self') {
      temp_prior <- priorobject
    } # if(sigma_formula_manual_prior_via_sigma_formula
    
    
    if(remove_sigma_parameter) {
      temp_prior <- temp_prior %>% dplyr::filter(class != 'sigma')
    }
      
   
    temp_stancode2 <- brms::make_stancode(formula = bformula,
                                    stanvars = bstanvars,
                                    prior = temp_prior,
                                    data = brmsdata)
    temp_standata2 <- brms::make_standata(formula = bformula,
                                    stanvars = bstanvars,
                                    prior = temp_prior,
                                    data = brmsdata)
    
   
    move_from_model_to_qq_for_bqinv <- 
      function(temp_stancode2x, 
               section = 'model',
               replacemuby = NULL,
               spfncname_c = NULL,
               spfncname_c_vector = NULL,
               decomp_editcode = decomp_editcode) {
      regex_for_section <- 
        paste(".*(",section,"\\s*\\{.*?\\priors including constants).*", 
              sep = '')
      filtered_stan_code <- gsub(temp_stancode2x, pattern = regex_for_section, 
                                 replacement = "\\1")
      filtered_stan_code <- gsub('if (!prior_only) {', '', 
                                 filtered_stan_code, fixed = T)
      filtered_stan_code <- gsub('vector[N] mu;', '', 
                                 filtered_stan_code, fixed = T)
      filtered_stan_code <- gsub("target.*", '', filtered_stan_code, fixed = F)
      
      
      editedcode2 <- filtered_stan_code
      clines_p <- strsplit(filtered_stan_code, "\n")[[1]]
      for (il in clines_p) {
        editedcode2 <- gsub(pattern = "//", replacement = "//", 
                            x = editedcode2, fixed = T)
        editedcode2 <- gsub(pattern = "//[^\\\n]*", replacement = "", 
                            x = editedcode2)
        editedcode2 <- gsub(paste0(il, ""), "", editedcode2, fixed = T)
      }
      
      zz <- strsplit(editedcode2, "\n")[[1]]
      zz_c <- c()
      for (iz in 1:length(zz)) {
        if(!is_emptyx(gsub_space(zz[iz]))) {
          zz_in <- zz[iz]
          zz_in <- paste0("  ", zz_in)
          zz_c <- c(zz_c, zz_in)
        }
      }
      
      zz_c <- paste(zz_c, collapse = "\n")
      
      if(!decomp_editcode) {
        htx <- zz_c
        htx <- strsplit(htx, "\n")[[1]]
        nonmulines <- mulines <- c()
        for (htxi in htx) {
          tzx <- gsub("[[:space:]]", "", htxi)
          if(!grepl('mu', tzx)) {
            nonmulines <- c(nonmulines, htxi)
          } else if(grepl('mu', tzx)) {
            mulines <- c(mulines, htxi)
          }
        }
        nonmulines <- paste0(nonmulines, collapse = "\n")

        lines_mu_subs <- c()
        for (htxi in mulines) {
            htxi_ <- gsub("[[:space:]]", "", htxi)
            if(grepl("=", htxi)) {
              htxi_ <- gsub("=(", "=", htxi_, fixed = T)
              htxi_ <- gsub("));", ");", htxi_, fixed = T)
              htxi_name <- strsplit(htxi_, "(", fixed = T)[[1]][1]
              htxi_2 <- regmatches(htxi_, gregexpr("(?<=\\().*?(?=\\))", 
                                                   htxi_, perl=T))[[1]]
              htxi_2 <- strsplit(htxi_2, ",", fixed = T)[[1]]
              htxi_3 <- c()
              nlp_s_number <- c()
              for (htxi_2i in htxi_2) {
                if(grepl('nlp_', htxi_2i) & grepl('_s', htxi_2i)) {
                  nlpsparms <- TRUE
                } else {
                  nlpsparms <- FALSE
                }
                if(!nlpsparms) htxi_3 <- c(htxi_3, htxi_2i)
                if(nlpsparms) nlp_s_number <- c(nlp_s_number, htxi_2i)
              }
              htxi_4 <- paste(htxi_3, collapse = ",")
              htxi_5 <- paste0(htxi_name, "(", htxi_4, ");")
              htxi_5 <- gsub("mu", " XR_inv", htxi_5)
              htxi_5 <- gsub("=", " = ", htxi_5)
              htxi_5 <- gsub(spfncname, paste0(spfncname, 
                                                    'QRsmatinv'), htxi_5)
              npsn   <- length(nlp_s_number)
              htxi_6 <- paste0("matrix[", npsn, ", ", npsn, "] ", htxi_5)
              lines_mu_subs <- c(lines_mu_subs, htxi_6)
            } # if(grepl("=", htxi)) {
        } # for (htxi in mulines) {
        lines_mu_subs <- paste0(lines_mu_subs, collapse = "\n")
        zz_c <- paste0(nonmulines, "\n", "  ", lines_mu_subs)
        zz_c2 <- zz_c
      }
        
      if(decomp_editcode) {
        zz_c_ <- strsplit(zz_c, "\n", fixed = T)[[1]]
        zz_c2 <- c()
        for (zz_ci in 1:length(zz_c_)) {
          if(!grepl('mu=', gsub_space(zz_c_[zz_ci])) ) {
            g <- zz_c_[zz_ci]
          }
          zz_c2 <- c(zz_c2, g)
        }
        zz_c2 <- paste(zz_c2, collapse = "\n")
      }
      return(zz_c2)
    }
    
  
    if(!is.null(decomp)) {
      if(add_rcsfunmatqrinv_genquant ) {
        temp_stancode_gqinv <- 
          brms::make_stancode(formula = bformula,
                              stanvars = bstanvars,
                              prior = brmspriors,
                              threads = brms::threading(NULL),
                              data = brmsdata)
        
        gq_funs_2 <- list()
        for (gq_funslen in 1:length(gq_funs)) {
          gq_funs_2[[gq_funslen]] <- gsub(gsub("\\;.*", "", 
                                               gq_funs[[gq_funslen]]),
                                          "", gq_funs[[gq_funslen]], 
                                          fixed = T)
        }
        
        qgcode <- 
          move_from_model_to_qq_for_bqinv(
            temp_stancode_gqinv, 
            replacemuby = NA, 
            spfncname_c = spfncname_c,
            spfncname_c_vector = spfncname_c_vector,
            decomp_editcode = decomp_editcode)
        
        
        qgcode <- gsub("\n  }\n  }", "\n  }", qgcode) 
        
        gq_funs_2 <- paste(unlist(gq_funs_2), collapse = "\n")
        gq_funs_2 <- paste0(qgcode, '\n', gq_funs_2)
        
        bstanvars <- bstanvars + 
          brms::stanvar(scode = gq_funs_2, block = "genquant", 
                        position = "end")
      } # if(add_rcsfunmatqrinv_genquant ) {
      
      
      
      if(decomp_editcode & add_rcsfunmatqrinv_genquant) {
        spfncname_c_mat <- c()
        spfncname_c_vector <- c()
        for (spfncname_ci in spfncname_c) {
          spfncname_ci2 <- gsub(spfncname_ci, 
                                paste0(spfncname_ci, 'QRsmat'), spfncname_ci)
          spfncname_invmat <- gsub(spfncname_ci, 
                                   paste0(spfncname_ci, 'QRsmatinv'), 
                                   spfncname_ci)
          waht_C <- 'C_1'
          waht_Cby <- paste0(waht_C, 'X')
          spfncname_c_vector <- c(spfncname_c_vector, waht_C)
          tdcode <- paste0('matrix[', 'N', ',', 4, ']', " ", waht_Cby, " = ",
                           spfncname_ci2, "(", waht_C, ",",
                           'rep_vector(0.0, num_elements(', waht_C, '))',
                           ");")
          waht_Cby_inv <- 'XR_inv'
          tdcode_inv <- paste0('matrix[', 4, ',', 4, ']', " ", 
                               waht_Cby_inv, " = ",
                               spfncname_invmat, "(", waht_C, ",",
                               'rep_vector(0.0, num_elements(', waht_C, '))',
                               ");")
          
          tdcode <- paste0(tdcode, "\n", tdcode_inv)
          spfncname_c_mat <- c(spfncname_c_mat, tdcode)
        }
        spfncname_c_mat2 <- paste(spfncname_c_mat, collapse = "\n")
        bstanvars <- bstanvars + brms::stanvar(scode = spfncname_c_mat2, 
                                               block = "tdata", 
                                               position = "end")
      } # if(decomp_editcode & add_rcsfunmatqrinv_genquant) {
    } # if(!is.null(decomp)) {
    
    
    
    if(vcov_init_0e) {
      initialsx2 <- brmsinits
      for (initialsi in names(initialsx2)) {
        if(grepl("sd_", initialsi)) {
          if(!grepl("sd_nu", initialsi, fixed = T)) {
            initialsx2[[initialsi]] <- NULL
            newinits <- set_init_gr_effects(temp_stancode2, 
                                            temp_standata2,
                                            parameterization = parameterization,
                                            what = 'sd')
            initialsx2 <- c(initialsx2, newinits)
          }
        }
        if(grepl("L_", initialsi)) {
          initialsx2[[initialsi]] <- NULL
          newinits <- set_init_gr_effects(temp_stancode2, 
                                          temp_standata2, 
                                          parameterization = parameterization,
                                          what = 'L')
          initialsx2 <- c(initialsx2, newinits)
        }
        if(grepl("z_", initialsi)) {
          initialsx2[[initialsi]] <- NULL
          newinits <- set_init_gr_effects(temp_stancode2, 
                                          temp_standata2, 
                                          parameterization = parameterization,
                                          what = 'z')
          initialsx2 <- c(initialsx2, newinits)
        }
      }
      uni_name <- unique(names(initialsx2))
      initialsx2 <- initialsx2[uni_name] 
      brmsinits <- initialsx2
    } 
    
    
    
    
    if(vcov_init_0e) {
      if(parameterization == 'cp') {
        initialsx2 <- brmsinits
        temp_stancode2cp <- edit_scode_ncp_to_cp_new(temp_stancode2, 
                                                 genq_only = FALSE, 
                                                 normalize = normalize, 
                                                 cp_via = cp_via)
        
        newinits <- set_init_gr_effects(temp_stancode2cp, 
                                        temp_standata2, 
                                        parameterization = parameterization,
                                        what = 'r')
        initialsx2 <- c(initialsx2, newinits)
        uni_name <- unique(names(initialsx2))
        initialsx2 <- initialsx2[uni_name] 
        brmsinits <- initialsx2
      }
    }
    
    brmsinits <- lapply(1:brms_arguments$chains, function(id) {
      eval_inits_fun(
        inits = brmsinits,
        jitter_init_beta = jitter_init_beta,
        jitter_init_sd = jitter_init_sd,
        jitter_init_cor = jitter_init_cor,
        digits = 4)
    })
  }
  
  
  


  # Add stanvars for logistic3e
  if(select_model_edit == 'logistic3e') {
    temp_stancode_logistic3e <- brms::make_stancode(formula = bformula,
                                          stanvars = bstanvars,
                                          prior = brmspriors,
                                          data = brmsdata)
    temp_standata_logistic3e <- brms::make_standata(formula = bformula,
                                                    stanvars = bstanvars,
                                                    prior = brmspriors,
                                                    data = brmsdata)
    
    
    check_p_dimes <- c()
    for (clines_tpi in names(temp_standata_logistic3e)) {
      for (igr in 1:9) {
        if(grepl(paste0("^K", "_"), clines_tpi) & 
           grepl(paste0("_", letters[igr]), clines_tpi)) {
           temd <- temp_standata_logistic3e[[clines_tpi]]
          check_p_dimes <- c(check_p_dimes, temd)
        }
      }
    }
    
    if(!all(check_p_dimes==check_p_dimes[1])) {
      stop2c('All parameters must have the same number of parameters')
    }
    
    
    check_p_attr1 <- c()
    for (clines_tpi in names(temp_standata_logistic3e)) {
      for (igr in 1:9) {
        if(grepl(paste0("^X", "_"), clines_tpi) & 
           grepl(paste0("_", letters[igr]), clines_tpi)) {
          temd <- temp_standata_logistic3e[[clines_tpi]]
          temd2 <- attr(temd, "assign")[1]
          check_p_attr1 <- c(check_p_attr1, temd2)
        }
      }
    }
    
    
    if(check_p_dimes[1] > 1) {
      if(any(check_p_attr1 == 0)) {
        stop2c('All parameters must have the covariate form as ~0+')
      }
    }
    
    
    outlogistic3e <- edit_scode_for_logistic3(temp_stancode_logistic3e, 
                                              normalize = normalize)
    
    bstanvars <- bstanvars + brms::stanvar(scode = outlogistic3e$pcode, 
                                           block = "parameters", 
                                           position = "start")
    
    bstanvars <- bstanvars + brms::stanvar(scode = outlogistic3e$fcode, 
                                           block = "functions")
    
    edit_ncov  <- as.integer(check_p_dimes[1])
    edit_npar  <- as.integer(3)
    edit_min_d <- array(rep(0, edit_ncov), dim = edit_ncov)
    edit_max_d <- array(rep(500,  edit_ncov), dim = edit_ncov)
    edit_min_v <- array(rep(0, edit_ncov), dim = edit_ncov)
    edit_max_v <- array(rep(2.5,    edit_ncov), dim = edit_ncov)
    edit_min_t <- array(rep(0,    edit_ncov), dim = edit_ncov)
    edit_max_t <- array(rep(18,   edit_ncov), dim = edit_ncov)
    bstanvars <- bstanvars + brms::stanvar(x = edit_ncov, name = 'Kedit',
                                           block = "data")
    bstanvars <- bstanvars + brms::stanvar(x = edit_npar, name = 'Cedit',
                                           block = "data")
    bstanvars <- bstanvars + brms::stanvar(x = edit_min_d, name = 'min_d',
                                           block = "data")
    bstanvars <- bstanvars + brms::stanvar(x = edit_max_d, name = 'max_d',
                                           block = "data")
    bstanvars <- bstanvars + brms::stanvar(x = edit_min_v, name = 'min_v',
                                           block = "data")
    bstanvars <- bstanvars + brms::stanvar(x = edit_max_v, name = 'max_v',
                                           block = "data")
    bstanvars <- bstanvars + brms::stanvar(x = edit_min_t, name = 'min_t',
                                           block = "data")
    bstanvars <- bstanvars + brms::stanvar(x = edit_max_t, name = 'max_t',
                                           block = "data")
    
  } # if(select_model_edit == 'logistic3e') {
  
  
  # Set brm arguments
  setup_brms_args <-
    function(formula,
             prior,
             stanvars,
             data,
             init_set,
             init_str,
             init_r,
             seed,
             verbose,
             setarguments,
             brmsdots) {
      
      exc_args <- c("formula", "prior", "stanvars", "init", "data")
      if (eval(setarguments$backend) == "rstan")
        exc_args <- c(exc_args, "stan_model_args")
      for (exc_argsi in exc_args) {
        if (exc_argsi %in% names(setarguments))
          setarguments[[exc_argsi]] <- NULL
      }
      setarguments$formula <- formula
      setarguments$prior <- prior
      setarguments$stanvars <- stanvars
      setarguments$data <- data
      
      if (eval(setarguments$backend) == "cmdstanr") {
        if (all(sapply("0", grepl, init_str))) {
          setarguments$init <- 0
          custom_init <- FALSE
        } else if (all(sapply("random", grepl, init_str))) {
          setarguments$init <- NULL
          custom_init <- FALSE
        } else {
          setarguments$init <- init_set
          custom_init <- TRUE
        }
        if (!custom_init & !is.null(init_r)) {
          setarguments$init <- init_r
        }
      }
      
      if (eval(setarguments$backend) == "rstan") {
        if (all(sapply("0", grepl, init_str))) {
          setarguments$init <- "0"
          custom_init <- FALSE
        } else if (all(sapply("random", grepl, init_str))) {
          setarguments$init <- "random"
          custom_init <- FALSE
        } else {
          setarguments$init <- init_set
          custom_init <- TRUE
        }
        if (!custom_init & !is.null(init_r)) {
          setarguments$init_r <- init_r
        }
      }
      
      if (eval(setarguments$backend) == "rstan" | 
          eval(setarguments$backend) == "cmdstanr") {
        if (is.null(eval(setarguments$control))) {
          setarguments$control <- list(adapt_delta = 0.8, max_treedepth = 15)
        }
        if (is.na(eval(setarguments$seed))) {
          setarguments$seed <- seed
        }
        
        cores_   <- eval(setarguments$cores)
        threads_ <- eval(setarguments$threads)
        
        
        if(!is.null(getOption('mc.cores'))) {
          # cores_ <- NULL
        }
        
      if(is.null(cores_) & is.null(getOption('mc.cores'))) {
        max.cores <- 1 # getOption("mc.cores", 1) -> from ?brms::brm
      }
        
       if(!is.null(cores_)) {
         if(cores_ == "maximise") {
           max.cores <- 
             as.numeric(future::availableCores(methods = "system", omit = 0))
           if(max.cores < 1) max.cores <- 1
         } else if(cores_ == "optimize") {
           max.cores <- 
             as.numeric(future::availableCores(methods = "system", omit = 1))
           if(max.cores < 1) max.cores <- 1
           if(max.cores > eval(setarguments$chains)) {
             max.cores <- eval(setarguments$chains)
           }
         } else if(check_is_numeric_like(cores_)) {
           max.cores <- eval(cores_)
         } else {
           stop2c("Argument cores must be either an integer or a character
                  'maximise' or 'optimize'")
         }
       } else if(!is.null(getOption('mc.cores'))) {
         if(is.null(cores_)) {
           max.cores <- getOption('mc.cores')
         } else if(!is.null(cores_)) {
           if(cores_ != "maximise" & cores_ != "optimize") {
             max.cores <- getOption('mc.cores')
           } else if(!check_is_numeric_like(cores_)) {
             max.cores <- getOption('mc.cores')
           }
         } # if(is.null(cores_)) { else if(!is.null(cores_)) {
         # max.cores <- getOption('mc.cores')
       } # if(!is.null(cores_)) { else if(!is.null(getOption('mc.cores'))) {
        
       setarguments$cores <-  max.cores
        
        
        if(!is.list(threads_)) {
          if(is.null(threads_)) {
            threads_ <- deparse(threads_)
          }
          if( is.character(threads_) & threads_ == "maximise") {
            max.threads <- threads_char(threads_, 
                                        chains = eval(setarguments$chains))
          } else if( is.character(threads_) & threads_ == "optimize") {
            max.threads <- threads_char(threads_, 
                                        chains = eval(setarguments$chains))
          } else if(!is.null(getOption('brms.threads')) &
                    (is.character(threads_) & threads_ != "maximise") &
                    (is.character(threads_) & threads_ != "optimize")) {
            max.threads <- getOption('brms.threads')
          } else if(is.null(getOption('brms.threads')) &
                    (is.character(threads_) & threads_ != "maximise") &
                    (is.character(threads_) & threads_ != "optimize")) {
            max.threads <- getOption('brms.threads')
          } else {
            max.threads <- eval(setarguments$cores)
          }
          setarguments$threads <-  brms::threading(max.threads)
        }
      } 
      
      if (eval(setarguments$backend) == "cmdstanr") {
        if (is.list(eval(setarguments$stan_model_args)) &
            eval(length(setarguments$stan_model_args)) == 0) {
          setarguments$stan_model_args <- list(
            # pedantic = FALSE,
            
            # Setting this leads to error or multiple --O1 stanflag
            
            # Setting this results in some compilation error - Eigen
            # stanc_options = list("O1")
            
            # , cpp_options = list(#'STAN_CPP_OPTIMS=true',
            #                      # 'CXXFLAGS = -O2',
            #                      # 'STANCFLAGS+= --warn-pedantic --O0',
            #                      'STAN_NO_RANGE_CHECKS=true'
            #                      )
            )
        }
      }
      
      
      if(verbose) {
        message2c(setarguments$stan_model_args)
      }
      
      
      
      if (eval(setarguments$backend) == "rstan" & 
          packageVersion("rstan") < "2.26.1") {
        setarguments$threads <- setarguments$threads 
      }
      
      if (eval(setarguments$backend) == "mock") {
        max.threads <- getOption('brms.threads')
        setarguments$threads <- brms::threading(max.threads)
        max.cores <- getOption('mc.cores')
        setarguments$cores <-  max.cores
      }
      
      if (length(brmsdots) > 0) {
        setarguments <- c(setarguments, brmsdots)
      }
      return(setarguments)
    }
  
  
  
  if (verbose) {
    setmsgtxt <- paste0("\n Setting-up brms arguments")
    if (displayit == 'msg') {
      message2c(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  
  
  brmsdots_ <- list(...)
  # 24.08.2024
  getdotslistnames <- c("match_sitar_a_form", "match_sitar_d_form",
                         "sigmamatch_sitar_a_form", "displayit", 
                        "setcolh", "setcolb", "smat", "decomp")
  
  for (getdotslisti in getdotslistnames) {
    brmsdots_[[getdotslisti]] <- NULL
  }
  
  for (collect_dot_namesi in collect_dot_names) {
    if(!is.null(brmsdots_[[collect_dot_namesi]])) 
      brmsdots_[[collect_dot_namesi]] <- NULL
  }
 
  if(!is.null(custom_stanvars)) {
    bstanvars <- bstanvars + custom_stanvars
  }
  
  
  brm_args <-
    setup_brms_args(
      formula = bformula,
      prior = brmspriors,
      stanvars = bstanvars,
      data = brmsdata,
      init = brmsinits,
      init_str = initialslist_s,
      init_r = brmsinits_r,
      seed = seed,
      verbose = verbose,
      setarguments = brms_arguments,
      brmsdots = brmsdots_)
  
  # 27.02.2025
  # when fitting univariate_by, the subset is found in brm_args, why?
  # Just drop it 
  brm_args$subset <- NULL
  
  # 14.05.2025
  brm_args$fast_nsk <- NULL
  
  if(!is.null(custom_family)) {
    brm_args$family <- custom_family
  }
  
  if (verbose) {
    setmsgtxt <- paste0("\n Fitting model")
    if (displayit == 'msg') {
      message2c(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  cat("\n")
  insert_new_priors <- function(setdf_1, setdf_2) {
    cc <- zz <- list()
    for (i in 1:nrow(setdf_1)) {
      getx <- setdf_1[i,]
      zz[[i]] <- setdf_2 %>% 
        dplyr::mutate(prior = dplyr::if_else(.data$class == getx[['class']] &
                                               .data$coef == getx[['coef']] &
                                               .data$group == getx[['group']] &
                                               .data$resp == getx[['resp']] &
                                               .data$dpar == getx[['dpar']] &
                                               .data$nlpar == getx[['nlpar']],
                                             .data$getx$prior,
                                             .data$setdf_2$prior)) %>% 
        dplyr::filter(.data$class == getx[['class']] &
                        .data$coef == getx[['coef']] &
                        .data$group == getx[['group']] &
                        .data$resp == getx[['resp']] &
                        .data$dpar == getx[['dpar']] &
                        .data$nlpar == getx[['nlpar']])
      
      cc[[i]] <- setdf_2 %>% 
        dplyr::mutate(prior = dplyr::if_else(.data$class != getx[['class']] &
                                               .data$coef != getx[['coef']] &
                                               .data$group != getx[['group']] &
                                               .data$resp != getx[['resp']] &
                                               .data$dpar != getx[['dpar']] &
                                               .data$nlpar != getx[['nlpar']] ,
                                             .data$setdf_2$prior,
                                             .data$getx$prior)) %>% 
        
        dplyr::filter(.data$class != getx[['class']] & 
                        .data$coef != getx[['coef']] &
                        .data$group != getx[['group']] &
                        .data$resp != getx[['resp']] &
                        .data$dpar != getx[['dpar']] &
                        .data$nlpar != getx[['nlpar']])
      
    }
    p1 <- cc %>% CustomDoCall(rbind, .)
    p2 <- zz %>% CustomDoCall(rbind, .)
    p1p2 <- rbind(p1, p2)
    return(p1p2)
  }
  
  if(set_higher_levels) {
    brmspriors_sdcor <- brmspriors %>% 
      dplyr::filter(.data$class == 'sd' | .data$class == 'cor')
    brmspriors_sdcor_gr <- brmspriors_sdcor$group
    
    brmsfit_sdcor <- CustomDoCall(brms::get_prior, brm_args) %>% 
      dplyr::filter(.data$class == 'sd' | .data$class == 'cor')
    
    brmsfit_sdcor_prior_gr <- brmsfit_sdcor %>% 
      dplyr::filter(!.data$group %in%  brmspriors_sdcor_gr)
    
    brmspriors_brmsfit_sdcor <- brmspriors %>% 
      dplyr::bind_rows(., brmsfit_sdcor_prior_gr)
    
    brmspriors <- brmspriors_brmsfit_sdcor
  }
  
  brm_args$prior <- brmspriors
  
  if(!is.null(set_self_priors) & 
     !is.null(add_self_priors) & 
     !is.null(set_replace_priors)) {
    stop2c("Amongst 'set_self_priors', 
         'add_self_priors' and 'set_replace_priors' arguments,",
         "\n ",
         " only one can be specified at a time")
  }
  
  if(get_priors & get_priors_eval & validate_priors & 
     get_stancode & get_standata & get_formula & get_stanvars) {
    stop2c("Amongst 'get_priors' 'get_priors_eval', 'validate_priors' ",
         "\n ",
         "'get_stancode', 'get_standata', 'get_formula', 'get_stanvars' ",
         "\n ",
         " arguments, only one can be set to TRUE at a time")
  }
  
  lbbb_ <- ubbb_ <- NULL
  tempprior_hold <- brmspriors # brm_args$prior 
  setpriornamesorder <- colnames(tempprior_hold)
  tempprior_hold$lbbb_ <- tempprior_hold$lb
  tempprior_hold$ubbb_ <- tempprior_hold$ub
  tempprior_hold$lb <- tempprior_hold$ub <- NULL
  tempprior_hold <- tempprior_hold %>% 
    dplyr::mutate(lbbb_ = dplyr::if_else(class == 'sd', NA, lbbb_))
  tempprior_hold <- tempprior_hold %>% 
    dplyr::mutate(ubbb_ = dplyr::if_else(class == 'sd', NA, ubbb_))
  tempprior_hold$lb <- tempprior_hold$lbbb_
  tempprior_hold$ub <- tempprior_hold$ubbb_
  tempprior_hold$lbbb_ <- tempprior_hold$ubbb_ <- NULL
  tempprior_hold <- tempprior_hold %>% 
    dplyr::relocate(dplyr::all_of(setpriornamesorder))
  brmspriors <-   tempprior_hold
  
  if(!is.null(set_self_priors) & 
     is.null(add_self_priors) &
     is.null(set_replace_priors)) {
    brmspriors <- set_self_priors
  }
  

  if(is.null(set_self_priors) & is.null(set_replace_priors)) {
    brmspriors <- brmspriors
  }
  
  # 24.08.2024
  if(!is.null(add_self_priors)) {
    add_self_priors <- add_self_priors %>%  
      dplyr::filter(source == 'user' & coef != "")
    
    if (is.null(sigma_formula_manualsi[[1]][1])) {
      brmspriors_toadd <- brmspriors
    } else if(sigma_formula_manualsi == "NULL") {
      brmspriors_toadd <- brmspriors
    } else {
      brmspriors_toadd <-  brmspriors %>% 
        dplyr::filter(., !grepl('sigma', .data[['nlpar']]))
    }
    
    brmspriors <- brmspriors_toadd %>% dplyr::bind_rows(., add_self_priors)
  } # if(!is.null(add_self_priors)) {
  
  
  ################################################################
  
  for (j in 1:length(setsigmaxvar_names_val)) {
    if(any(setsigmaxvar_names_val[j])) {
      brmspriors <- brmspriors %>% dplyr::filter(class != "sigma")
    }
  }
  
  ################################################################
 
  brm_args$prior <- brmspriors
  
  decomp_escode2<- function(temp_stancode2x) {
    htx <- strsplit(temp_stancode2x, "\n")[[1]]
    lines_all <- c()
    for (htxi in 1:length(htx)) {
      htxi_ <- htx[htxi]
      if(grepl('^mu', gsub("[[:space:]]", "", htxi_))) {
        htxi_ <- gsub("[[:space:]]", "", htxi_)
        htxi_ <- gsub("=(", "=", htxi_, fixed = T)
        htxi_ <- gsub("));", ");", htxi_, fixed = T)
        htxi_name <- strsplit(htxi_, "(", fixed = T)[[1]][1]
        if(grepl('[', htxi_, fixed = T)) {
          htxi_c1 <- strsplit(htxi_, "[", fixed = T)[[1]][1]
        } else {
          htxi_c1 <- strsplit(htxi_, ",", fixed = T)[[1]][1]
        }
        htxi_others <- strsplit(htxi_, ",", fixed = T)[[1]][-1]
        htxi_others <- paste(htxi_others, collapse = ", ")
        htxi_c1 <- paste0(htxi_c1, 'X')
        htxi_c1 <- strsplit(htxi_c1, "(", fixed = T)[[1]][2]
        htxi_name <- paste0(htxi_name, 'X')
        if(grepl('[', htxi_, fixed = T)) {
          htxi_c1 <- paste0(htxi_c1, '[start:end,]')
        } else {
          htxi_c1 <- htxi_c1
        }
        htxi_final <- paste0(htxi_name, "(", htxi_c1, ", ", htxi_others)
      } else {
        htxi_final <- htxi_
      }
      lines_all <- c(lines_all, htxi_final)
      lines_all <- paste(lines_all, collapse = "\n")
    }
    dvciit <- paste0('data vector', " ", 'C_1')
    dvciby <- paste0('data matrix', " ", 'C_1', 'X')
    lines_all <- gsub(dvciit, dvciby, lines_all, fixed = T)
    dvciit <- ', C_1'
    dvciby <-  paste0(dvciit, 'X')
    lines_all <- gsub(dvciit, dvciby, lines_all, fixed = T)
    return(lines_all)
  }
  
  ####################################################################

  # stringr::str_extract
  stringr_str_extract_base_regexpr <- function(text, 
                                               pattern, 
                                               replacement = NULL, 
                                               all = FALSE,
                                               returnmatch = FALSE) {
    if(!all) matches   <- regexpr(pattern = pattern, text = text)
    if( all) matches   <- gregexpr(pattern = pattern, text = text)
    setitscode_i_temp1 <- regmatches(text, m = matches)
    if(returnmatch) return(setitscode_i_temp1)
    if(is_emptyx(setitscode_i_temp1)) {
      setitscode_i <- text
    } else {
      setitscode_i_temp2 <- paste0(setitscode_i_temp1, replacement)
      setitscode_i_temp3 <- gsub(setitscode_i_temp1, setitscode_i_temp2, 
                                 text) # setitscode_i_temp
      setitscode_i <- setitscode_i_temp3
    }
    return(setitscode_i)
  } 
  
  
  string_patterns_replacements <- function(string, 
                                           patterns, 
                                           replacements) {
    for (i in seq_along(patterns))
      string <- gsub(patterns[i], replacements[i], string, perl=T)
    string
  } 
  
  
  
  build_missing_via_fun <- function(what, whatabc, bere,
                                    setitscode_i_rep_vector_exa,
                                    outcomes, param = 'fixed',
                                    checkfire = NULL) {
    
    check_for_ipaterns <- paste0(paste0("_", whatabc, "$"), 
                                 collapse = "|")
    
    grepl_defined_fixed_abci_xxx <- c()
    for (grepl_defined_fixed_abci in what) {
      if(grepl(check_for_ipaterns, grepl_defined_fixed_abci)) {
        grepl_defined_fixed_abci_xxx <- c(grepl_defined_fixed_abci_xxx, 
                                          grepl_defined_fixed_abci)
      }
    }
    
    if(length(outcomes) == 1) {
      find_uniqie_outcomes <- "nlp"
    } else {
      find_uniqie_outcomes <- gsub("yvar_", "nlp_", names(outcomes))
    }
    
    find_uniqie_outcomes2 <- c()
    for (find_uniqie_outcomesi in find_uniqie_outcomes) {
      for (i in whatabc) {
        dhhhc <- paste0(find_uniqie_outcomesi, "_", i)
        find_uniqie_outcomes2 <- c(find_uniqie_outcomes2, dhhhc)
      }
    }
    whatabc_all_defined <- unique(find_uniqie_outcomes2)
    
    whatabc_not_defined <- setdiff(whatabc_all_defined,
                                   grepl_defined_fixed_abci_xxx)
    
    # create fixed and be rep_vector if parm not defined in random
    if(!is.null(checkfire)) {
      if(param != "random") {
        checkfire2 <- checkfire
        names(checkfire2) <- gsub("random", "nlp", names(checkfire2))
        whatabc_not_defined2 <- c()
        for (f in find_uniqie_outcomes) {
          for (w in whatabc_not_defined) {
            if(grepl(f, w)) {
              checkfire3 <- checkfire2[[f]]
              j <- strsplit(w , "_(?!.*_)", perl=TRUE)[[1]][2]
              for (i in checkfire3) {
                if(grepl(i, j)) {
                  whatabc_not_defined2 <- c(whatabc_not_defined2, w)
                } # if(grepl(checkfire3x, jjj)) {
              } # for (i in checkfire3) {
            } # if(grepl(f, w)) {
          } # for (w in whatabc_not_defined) {
        } # for (f in find_uniqie_outcomes) {
        whatabc_not_defined3 <- setdiff(whatabc_not_defined,
                                        whatabc_not_defined2)
        whatabc_not_defined <- whatabc_not_defined3
      } # if(param != "random") {
    } # if(!is.null(checkfire)) {
    
    
    # this because rep_vector() not defined earlier for random parameters
    if(param == "random") {
      whatabc_not_defined <- whatabc_all_defined
    }
  
    grepl_defined_fixed_abci_xxx2 <- c()
    for (whatabc_not_definedi in whatabc_not_defined) {
      for (setitscode_i_rep_vector_exai in setitscode_i_rep_vector_exa) { 
        defined_random <- 
          stringr_str_extract_base_regexpr(setitscode_i_rep_vector_exai, 
                                           "nlp_\\w+", 
                                           paste_x_Naux_str, all = T,
                                           returnmatch = T)[[1]]
        j <- strsplit(defined_random , "_(?!.*_)", perl=TRUE)[[1]][1]
        jj <- strsplit(whatabc_not_definedi , "_(?!.*_)", perl=TRUE)[[1]][1]
        if(grepl(jj, j)) {
          build_missing <- gsub(defined_random, whatabc_not_definedi, 
                                setitscode_i_rep_vector_exai)
        } # if(grepl(jj, j)) {
      } # for (setitscode_i_rep_vector_exai
      grepl_defined_fixed_abci_xxx2 <- c(grepl_defined_fixed_abci_xxx2, 
                                         build_missing)
    } # for (whatabc_not_definedi in whatabc_not_defined) {
    
    grepl_defined_fixed_abci_xxx3 <- c()
    for (grepl_defined_fixed_abci_xxx2i in grepl_defined_fixed_abci_xxx2) {
      grepl_defined_fixed_abci_xxx3_be <- 
        x_gsubit_gsubby(grepl_defined_fixed_abci_xxx2i,
                        gsubit = "nlp_", gsubby = bere,
                        pasteit = TRUE, fixed = FALSE)
      grepl_defined_fixed_abci_xxx3 <- c(grepl_defined_fixed_abci_xxx3,
                                         grepl_defined_fixed_abci_xxx3_be)
    }
    return(grepl_defined_fixed_abci_xxx3)
  } # end build_missing_via_fun
  
  
  
  check_d_defined_fun <- function(what, outcomes, param = NULL, lookfor = "d") {
    if(!is.null(param)) {
      add_parm <- paste0("_", param)
    } else {
      add_parm <- NULL
    }
    
    if(length(outcomes) == 1) {
      get_unique_length_ys <- "nlp" # paste0("nlp_", outcomes[[1]])
    } else {
      get_unique_length_ys <- gsub("yvar_", "nlp_", names(outcomes))
    }
    
    lookfor_str <- paste0("_", lookfor, "$")
    x_c <- c()
    for (j in what) {
      if(grepl(lookfor_str, j)) {
        x <- j
      } else {
        x <- NULL
      }
      x_c <- c(x_c, x)
    }
    # get_unique_length_ys <- sort(get_unique_length_ys)
    # x_c <- sort(x_c)
    if(!is.null(x_c)) {
      xxz_c <- c()
      for (i in get_unique_length_ys) {
        if(any(grepl(i, x_c))) {
          xxz <- paste0("int ", i, add_parm, "_", lookfor, 
                        "_", 'defined', " = 1;")
        } else if(!all(grepl(i, x_c))) {
          xxz <- paste0("int ", i, add_parm, "_", lookfor, 
                        "_", 'defined', " = 0;")
        }
        xxz_c <- c(xxz_c, xxz)
      }
    } else if(is.null(x_c)) {
      xxz_c <- c()
      for (i in get_unique_length_ys) {
        xxz <- paste0("int ", i, add_parm, "_", 
                      lookfor, "_", 'defined', " = 0;")
        xxz_c <- c(xxz_c, xxz)
      }
    }
    return(xxz_c)
  } # end check_d_defined_fun
  
  
  ####################################################################
  ####################################################################
  
  strsplit_type <- function(x,
                            split,
                            type = "remove",
                            perl = FALSE,
                            ...) {
    if (type == "remove") {
      # use base::strsplit
      out <- base::strsplit(x = x, split = split, perl = perl, ...)
    } else if (type == "before") {
      # split before the delimiter and keep it
      out <- base::strsplit(x = x,
                            split = paste0("(?<=.)(?=", split, ")"),
                            perl = TRUE,
                            ...)
    } else if (type == "after") {
      # split after the delimiter and keep it
      out <- base::strsplit(x = x,
                            split = paste0("(?<=", split, ")"),
                            perl = TRUE,
                            ...)
    } else {
      # wrong type input
      stop2c("type must be remove, after or before!")
    }
    return(out)
  }
  
  
  
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  
  if(genquant_xyadj) {
    make_stancode_custom_args <- brm_args
    make_standata_custom_args <- brm_args
    # Make sure threadus NULL, otherwise it will create grainsize etic
    make_stancode_custom_args[['threads']] <- brms::threading(NULL)
    make_standata_custom_args[['threads']] <- brms::threading(NULL)
    
    # check_for_validy_of_prepare_transformations makes sense only if same data
    if(is.null(data_custom)) {
      check_for_validy_of_prepare_transformations <- 
        check_for_validy_of_prepare_transformations
    } else {
      check_for_validy_of_prepare_transformations <- FALSE
    }
    
    # mandatory to make Naux variables - this is must for univariate_by
    use_data_custom <- FALSE
    if(is.null(data_custom)) {
      data_custom <- data.org.in
      use_data_custom <- FALSE
    } else {
      data_custom <- data_custom
      use_data_custom <- TRUE
    }
    
    
    ####################################################################
    ####################################################################
    
    data_custom.org.in                   <- data_custom
    prepare_data_args <- list()
    prepare_data_args[['data']]          <- data_custom
    prepare_data_args[['xvar']]          <- xs.org.in
    prepare_data_args[['yvar']]          <- ys.org.in
    prepare_data_args[['idvar']]         <- ids.org.in
    prepare_data_args[['univariate_by']] <- univariate_by
    prepare_data_args[['multivariate']]  <- multivariate
    prepare_data_args[['outliers']]      <- outliers
    prepare_data_args[['envir']]         <- enverr.
    prepare_data_args[['sigmaxvar']]     <- sigmaxs.org.in
    prepare_data_args[['nys']]           <- nys
    prepare_data_args[['subset']]        <- FALSE
    prepare_data_args[['returnys']]      <- FALSE
    prepare_data_args[['verbose']]       <- FALSE
    prepare_data_args[['displayit']]     <- displayit
    prepare_data_args[['setcolb']]       <- setcolb

    data_custom_data  <- CustomDoCall(prepare_data2, prepare_data_args)
    
    if(check_for_validy_of_prepare_transformations) {
      check_for_validy_of_prepare_transformations_0_custom <- data_custom_data
    }
    
    # xvarvaluelist etc must not be list but a name vector, expecially when ifun
    prepare_transformations_args[['data']]         <- data_custom_data
    prepare_transformations_args[['xvar']]         <- xvarvaluelist %>% unlist()
    prepare_transformations_args[['yvar']]         <- yvarvaluelist %>% unlist()
    prepare_transformations_args[['sigmaxvar']]    <- 
      sigmaxvarvaluelist %>% unlist()
    prepare_transformations_args[['xfun']]         <- xfuntransformvaluelist 
    prepare_transformations_args[['yfun']]         <- yfuntransformvaluelist
    prepare_transformations_args[['sigmaxfun']]    <- 
      sigmaxfuntransformvaluelist
    prepare_transformations_args[['ixfun']]        <- FALSE
    prepare_transformations_args[['iyfun']]        <- FALSE
    prepare_transformations_args[['sigmaixfun']]   <- FALSE
    prepare_transformations_args[['xoffset']]      <- xoffsetvaluelist
    prepare_transformations_args[['sigmaxoffset']] <- sigmaxoffsetvaluelist
    prepare_transformations_args[['transform']]    <- ""
    prepare_transformations_args[['itransform']]   <- ""
    
    data_custom_data <- CustomDoCall(prepare_transformations, 
                                     prepare_transformations_args)
    
    if(check_for_validy_of_prepare_transformations) {
      check_for_validy_of_prepare_transformations_3_custom <- data_custom_data
    }
    
    
    if(check_for_validy_of_prepare_transformations) {
      prepare_transformations_args[['data']]         <- data_custom_data
      prepare_transformations_args[['yfun']]         <- NULL
      prepare_transformations_args[['ixfun']]        <- TRUE
      prepare_transformations_args[['sigmaixfun']]   <- TRUE
      prepare_transformations_args[['transform']]    <- ""
      prepare_transformations_args[['itransform']]   <- ""
      
      check_for_validy_of_prepare_transformations_4_custom <- 
        CustomDoCall(prepare_transformations, prepare_transformations_args)
      
      prepare_transformations_args[['data']]         <- 
        check_for_validy_of_prepare_transformations_4_custom
      prepare_transformations_args[['ixfun']]        <- FALSE
      prepare_transformations_args[['sigmaixfun']]   <- FALSE
      prepare_transformations_args[['transform']]    <- ""
      prepare_transformations_args[['itransform']]   <- ""
      
      check_for_validy_of_prepare_transformations_5_custom <- 
        CustomDoCall(prepare_transformations, prepare_transformations_args)
      
      
      # Ignore outcome ysi becuase it remains transformed
      check_for_validy_of_prepare_transformations_0_custom <- 
        check_for_validy_of_prepare_transformations_0_custom %>% 
        # dplyr::select(-ys)
        dplyr::select(-dplyr::all_of(ys))
      check_for_validy_of_prepare_transformations_4_custom <- 
        check_for_validy_of_prepare_transformations_4_custom %>% 
        dplyr::select(-dplyr::all_of(ys))
      check_for_validy_of_prepare_transformations_3_custom <- 
        check_for_validy_of_prepare_transformations_3_custom %>% 
        dplyr::select(-dplyr::all_of(ys))
      check_for_validy_of_prepare_transformations_5_custom <- 
        check_for_validy_of_prepare_transformations_5_custom %>% 
        dplyr::select(-dplyr::all_of(ys))
      
      if(!isTRUE(all.equal(check_for_validy_of_prepare_transformations_0_custom,
                           check_for_validy_of_prepare_transformations_4_custom
                           ))) {
        stop2c("Something wrong with 'prepare_transformations' for data_custom")
      }
      
      if(!isTRUE(all.equal(check_for_validy_of_prepare_transformations_3_custom,
                           check_for_validy_of_prepare_transformations_5_custom
                           ))) {
        stop2c("Something wrong with 'prepare_transformations' for data_custom")
      }
      
    } # if(check_for_validy_of_prepare_transformations) {
    
   
    
    ################################################################
    
    for (j in 1:length(setsigmaxvar_names_val)) {
      if(!setsigmaxvar_names_val[j]) {
        data_custom_data <- data_custom_data %>% 
          dplyr::select(-sigmaxvarvaluelist[[j]])
      }
    }
    
    ################################################################
    
    if(check_for_validy_of_prepare_transformations) {
      # compare with original
      if(!isTRUE(all.equal(dataout,
                           data_custom_data))) {
        stop2c("Something wrong with 'prepare_transformations' for data_custom")
      }
    }
    
    data_custom_data_ys            <- attr(data_custom_data, "ys")
    data_custom_data_subindicators <- attr(data_custom_data, "subindicators")
    
    ####################################################################
    ####################################################################
    
    if(!is.null(data_custom_data)) {
      make_stancode_custom_args[['data']] <- data_custom_data
      make_standata_custom_args[['data']] <- data_custom_data
    }
    
    
    make_standata_custom_args_exclude <- c('stanvars', 'prior_only')
    for (i in make_standata_custom_args_exclude) {
      make_standata_custom_args[[i]] <- NULL
    }
    
    data_patterns_search_for_replace <- c("Z_\\w+", 
                                          "X_\\w+", 
                                          "N_\\w+", 
                                          "K_\\w+",
                                          "M_\\w+",
                                          "J_\\w+",
                                          "C_\\w+", 
                                          "Y_\\w+", # multivariate
                                          "\\bnresp\\b", # multivariate
                                          "\\bnrescor\\b", # multivariate
                                          "\\bN\\b",
                                          "\\bY\\b")
    
    data_patterns_search_for_replace_bar <- 
      paste(data_patterns_search_for_replace, 
            collapse = "|")
    
    
    if(is.null(data_custom)) {
      paste_x_Naux_str <- ""
    } else if(!is.null(data_custom)) {
      paste_x_Naux_str <- paste0("_", x_Naux_str)
      make_stancode_custom_args[['return_data']] <- TRUE
      gqdata_stanvarlist_data_si <- CustomDoCall(make_stancode_custom, 
                                            make_stancode_custom_args)
      
      
      data_custom_standata <- CustomDoCall(make_standata_custom, 
                                      make_standata_custom_args)
      
      data_custom_standata[["prior_only"]] <- NULL
      
      data_custom_standata_code <- 
        strsplit(gqdata_stanvarlist_data_si, "\n")[[1]]
      
      data_custom_standata_code2 <- data_custom_standata_code
      
      data_custom_standata_code2 <- c()
      for (i in 1:length(data_custom_standata_code)) {
        tweasj <- normalize_stancode_custom(data_custom_standata_code[i])
        if(!is_emptyx(tweasj)) {
          data_custom_standata_code2 <- c(data_custom_standata_code2, tweasj)
        }
      }
      
      add_data2_c <- list()
      collect_datax <- list()
      collect_scodex <- collect_stanvarname <- c()
      for (j in data_custom_standata_code2) {
        get_i <- tail(strsplit(j, split=" ")[[1]],1)
        get_i2 <- gsub(";", "", get_i)
        # get_i3 <- paste0("^", get_i2, "$")
        get_i3 <- paste0("\\b", get_i2, "\\b")
        get_jdat <- data_custom_standata[grepl(get_i3, 
                                               names(data_custom_standata))]
        
        get_jstanvarname <- 
          stringr_str_extract_base_regexpr(get_i2, 
                                           data_patterns_search_for_replace_bar, 
                                           paste_x_Naux_str, all = T,
                                           returnmatch = F)[[1]]
        
        addxt <- 
          stringr_str_extract_base_regexpr(j, 
                                           data_patterns_search_for_replace_bar, 
                                           paste_x_Naux_str, all = T,
                                           returnmatch = T)[[1]]
        
        gsubby <- paste0(addxt, paste_x_Naux_str)
        
        get_jcode <- string_patterns_replacements(j, addxt, gsubby)
        
        add_data2_c[[get_jstanvarname]] <-
          brms::stanvar(x = get_jdat[[1]],
                        name = get_jstanvarname,
                        scode = get_jcode,
                        block = "data",
                        position = "start",
                        pll_args = NULL)
      }
      
      add_data2_stanvarlistlist <- c()
      for (i in 1:length(add_data2_c)) {
        add_data2_stanvarlistlist[i] <- paste0("add_data2_c[[", i, "]]")
      }
      add_data2stanvars <- eval(parse(text = paste(add_data2_stanvarlistlist,
                                                   collapse = "+")))
      
      brm_args$stanvars <- brm_args$stanvars + add_data2stanvars
    } # if(is.null(data_custom)) { else if(!is.null(data_custom)) {
    
    
    ####################################################################
    ####################################################################
    
    make_stancode_custom_args[['return_gq']] <- TRUE
    gqdata_stanvarlist_si <- CustomDoCall(make_stancode_custom, 
                                     make_stancode_custom_args)
    
    gq_custom_standata_code <- strsplit(gqdata_stanvarlist_si, "\n")[[1]]
    
    add_data2_c <- c()
    add_data1_c <- c()
    
    grepl_random_abc <- c()
    grepl_random_abc_expected <- letters[1:4]
    add_missing_random_abc <- c()
    grepl_fixed_abc <- c()
    grepl_fixed_abc_expected <- letters[1:4]
    add_missing_fixed_abc <- c()
    setitscode_i_rep_vector_exa <- c()
    
    grepl_defined_fixed_abc <- c()
    grepl_defined_random_abc <- c()
    counterx <- 0
    for (i in 1:length(gq_custom_standata_code)) {
      counterx <- counterx + 1
      setitx_i <- deparse(i)
      setitscode_i_temp <- gq_custom_standata_code[counterx]
      
      # Z_ -> random effect design matrix
      setitscode_i <- stringr_str_extract_base_regexpr(setitscode_i_temp, 
                                                       " Z_\\w+", 
                                                       paste_x_Naux_str)
      # X_ -> fixed effect design matrix
      setitscode_i <- stringr_str_extract_base_regexpr(setitscode_i, 
                                                       " X_\\w+", 
                                                       paste_x_Naux_str)
      # N_ -> outcome specific number of observations
      setitscode_i <- stringr_str_extract_base_regexpr(setitscode_i, 
                                                       "N_\\w+", 
                                                       paste_x_Naux_str)
      # N_ -> outcome specific predictor, age
      setitscode_i <- stringr_str_extract_base_regexpr(setitscode_i, 
                                                       "C_\\w+", 
                                                       paste_x_Naux_str)
      
      # N -> single N denoting total number of observations
      setitscode_i <- stringr_str_extract_base_regexpr(setitscode_i,
                                                       "\\bN\\b", 
                                                       paste_x_Naux_str)
      
      
      if(grepl("rep_vector", setitscode_i)) {
        if(grepl("vector", setitscode_i) & 
           grepl("nlp", setitscode_i) & 
           grepl("_a", setitscode_i)) {
          setitscode_i_rep_vector_exa <- c(setitscode_i_rep_vector_exa, 
                                           setitscode_i)
        }
      }
      
      if(grepl("rep_vector", setitscode_i)) {
        setitscode_i2 <- x_gsubit_gsubby(setitscode_i,
                                         gsubit = "nlp_", gsubby = "be_",
                                         pasteit = TRUE, fixed = FALSE)
        add_data2_c <- c(add_data2_c, setitscode_i2)
      }
      
      
      if(grepl("+=", setitscode_i, fixed = T)) {
        if(grepl("X_", setitscode_i, fixed = T)) {
          setitscode_i2 <- x_gsubit_gsubby(setitscode_i,
                                           gsubit = "nlp_", gsubby = "be_",
                                           pasteit = TRUE, fixed = FALSE)
          add_data2_c <- c(add_data2_c, setitscode_i2)
          defined_fixed <- 
            stringr_str_extract_base_regexpr(setitscode_i, 
                                             "nlp_\\w+", 
                                             paste_x_Naux_str, all = T,
                                             returnmatch = T)[[1]]
          grepl_defined_fixed_abc <- c(grepl_defined_fixed_abc, defined_fixed)
        }
      }
      
      if(grepl("{", setitscode_i, fixed = T)) {
        add_data2_c <- c(add_data2_c, setitscode_i)
      }
      if(grepl("+=", setitscode_i, fixed = T)) {
        if(grepl("Z_", setitscode_i, fixed = T)) {
          setitscode_i2 <- x_gsubit_gsubby(setitscode_i,
                                           gsubit = "nlp_", gsubby = "re_",
                                           pasteit = TRUE, fixed = FALSE)
          add_data2_c <- c(add_data2_c, setitscode_i2)
          defined_random <- 
            stringr_str_extract_base_regexpr(setitscode_i, 
                                             "nlp_\\w+", 
                                             paste_x_Naux_str, all = T,
                                             returnmatch = T)[[1]]
          grepl_defined_random_abc <- c(grepl_defined_random_abc, 
                                        defined_random)
        }
      }
      if(grepl("}", setitscode_i, fixed = T)) {
        add_data2_c <- c(add_data2_c, setitscode_i)
      }
      
      add_data1_c <- c(add_data1_c, setitscode_i)
    } # end for (i in 1:length(gq_custom_standata_code)) {
    
    
    randomvaluelist_forcheck        <- randomvaluelist
    names(randomvaluelist_forcheck) <- unlist(randomnamelist)
  
    
    yvarvaluelist_forcheck        <- yvarvaluelist
    names(yvarvaluelist_forcheck) <- unlist(ynamelist)
    
    # add missing
    build_missing_via_fun_org <- 
      build_missing_via_fun(grepl_defined_fixed_abc, letters[1:4], "",
                            setitscode_i_rep_vector_exa,
                            outcomes = yvarvaluelist_forcheck,
                            param = 'fixed', 
                            checkfire = randomvaluelist_forcheck)
    
    build_missing_via_fun_fixed <- 
      build_missing_via_fun(grepl_defined_fixed_abc, letters[1:4], "be_",
                            setitscode_i_rep_vector_exa,
                            outcomes = yvarvaluelist_forcheck,
                            param = 'fixed', 
                            checkfire = randomvaluelist_forcheck)
    
    build_missing_via_fun_random <- 
      build_missing_via_fun(grepl_defined_random_abc, letters[1:4], "re_",
                            setitscode_i_rep_vector_exa,
                            outcomes = yvarvaluelist_forcheck,
                            param = 'random', 
                            checkfire = NULL)
    
    d_not_defined_tf <- FALSE
    d_not_defined <- paste0("_", "d", "$")
    for (i in build_missing_via_fun_fixed) {
      if(grepl(d_not_defined, i)) {
        d_not_defined_tf <- TRUE
      }
    }
    for (i in build_missing_via_fun_random) {
      if(grepl(d_not_defined, i)) {
        d_not_defined_tf <- TRUE
      }
    }
    
    if(d_not_defined_tf) {
      d_not_defined_scode <- paste0("int d_not_defined = 1;")
    } else {
      d_not_defined_scode <- paste0("int d_not_defined = 0;")
    }
    
    
    build_missing_via_fun_org_si <-
      paste(build_missing_via_fun_org, collapse = "\n")
    
    build_missing_via_fun_fixed_si <-
      paste(build_missing_via_fun_fixed, collapse = "\n")
    
    build_missing_via_fun_random_si <-
      paste(build_missing_via_fun_random, collapse = "\n")
    
    add_missing_gqdata_stanvarlist_si <-
      paste0(build_missing_via_fun_org_si, "\n",
             build_missing_via_fun_fixed_si, "\n",
             build_missing_via_fun_random_si)
    
    
    # check d define 
    check_if_d_fixed <- check_d_defined_fun(grepl_defined_fixed_abc,
                                            outcomes = yvarvaluelist_forcheck,
                                            param = NULL, 
                                            lookfor = "d")
    
    check_if_d_fixed <- x_gsubit_gsubby(check_if_d_fixed,
                                        gsubit = "nlp_", gsubby = "be_",
                                        pasteit = TRUE, fixed = FALSE)
    
    
    check_if_d_random <- check_d_defined_fun(grepl_defined_random_abc,
                                             outcomes = yvarvaluelist_forcheck,
                                             param = NULL, 
                                             lookfor = "d")
    
    check_if_d_random <- x_gsubit_gsubby(check_if_d_random,
                                         gsubit = "nlp_", gsubby = "re_",
                                         pasteit = TRUE, fixed = FALSE)
    
    
    
    check_if_d_fixed_random <- c(check_if_d_fixed, check_if_d_random)
    check_if_d_fixed_random <- paste(check_if_d_fixed_random, collapse = "\n")
    
    add_missing_gqdata_stanvarlist_si <-
      paste0(add_missing_gqdata_stanvarlist_si, "\n",
             check_if_d_fixed_random)
    
    
    add_data1_c_pasted <- paste(add_data1_c, collapse = "\n")
    add_data2_c_pasted <- paste(add_data2_c, collapse = "\n")
    
    add_data12_c_pasted <- paste0(add_data1_c_pasted, "\n", 
                                  add_data2_c_pasted)
    
    # add all mising rep_vector at the beginin
    add_data12_c_pasted <- paste0(add_missing_gqdata_stanvarlist_si, "\n",
                                  add_data12_c_pasted)
    
    
    add_data2stanvars <- 
      brms::stanvar(x = NULL,
                    name = add_data12_c_pasted,
                    scode = add_data12_c_pasted,
                    block = "genquant",
                    position = "end",
                    pll_args = NULL)
    
    #######################################################################
    #######################################################################
    
    gqdata_stanvarlist_si <- add_data2stanvars
    
    # we have come out of ii loop over ysi, so recreate elements temporarily
    add_gqdata_stanvarlist_si_c <- c()
    for (xnss in ys) {
      if(nys == 1) {
        name_xyadj_tomeanx_t <- "tomeanx_true"
        name_xyadj_tomeanx_f <- "tomeanx_false"
        name_xyadj_tomeany_t <- "tomeany_true"
        name_xyadj_tomeany_f <- "tomeany_false"
        x_N_name <- paste0("N", "", "")
        Naux_name_define <- paste0(x_N_name, paste_x_Naux_str)
      } else if(nys > 1) {
        name_xyadj_tomeanx_t <- paste0("tomeanx_true", "_", xnss)
        name_xyadj_tomeanx_f <- paste0("tomeanx_false", "_", xnss)
        name_xyadj_tomeany_t <- paste0("tomeany_true", "_", xnss)
        name_xyadj_tomeany_f <- paste0("tomeany_false", "_", xnss)
        x_N_name <- paste0("N", "_", xnss)
        Naux_name_define <- paste0(x_N_name, paste_x_Naux_str)
      }
      
      
      add_gqdata_stanvarlist_si <- 
        paste0(
          paste0("vector[", Naux_name_define, "] ", name_xyadj_tomeanx_t, ";"),
          "\n", 
          paste0("vector[", Naux_name_define, "] ", name_xyadj_tomeanx_f, ";"),
          "\n", 
          paste0("vector[", Naux_name_define, "] ", name_xyadj_tomeany_t, ";"),
          "\n", 
          paste0("vector[", Naux_name_define, "] ", name_xyadj_tomeany_f, ";")
        )
      add_gqdata_stanvarlist_si_c <- c(add_gqdata_stanvarlist_si_c,
                                       add_gqdata_stanvarlist_si)
    }
    
    add_gqdata_stanvarlist_si <- paste(add_gqdata_stanvarlist_si_c, 
                                       collapse = "\n")
    
    
    # we have come out of ii loop over ysi, so recreate here temporarily
    define_add_gqdata_stanvarlist_si_c <- c()
    for (xnss in ys) {
      if(nys == 1) {
        name_xyadj_tomeanx_t <- "tomeanx_true"
        name_xyadj_tomeanx_f <- "tomeanx_false"
        name_xyadj_tomeany_t <- "tomeany_true"
        name_xyadj_tomeany_f <- "tomeany_false"
        xsi_as_C_name        <-  "C_1"
        ysi_as_Y_name        <-  "Y"
      } else if(nys > 1) {
        name_xyadj_tomeanx_t <- paste0("tomeanx_true", "_", xnss)
        name_xyadj_tomeanx_f <- paste0("tomeanx_false", "_", xnss)
        name_xyadj_tomeany_t <- paste0("tomeany_true", "_", xnss)
        name_xyadj_tomeany_f <- paste0("tomeany_false", "_", xnss)
        xsi_as_C_name        <- paste0("C", "_", xnss, "_", "1")
        ysi_as_Y_name        <- paste0("Y", "_", xnss)
      }
      
      xsi_as_C_name <- paste0(xsi_as_C_name, paste_x_Naux_str)
      ysi_as_Y_name <- paste0(ysi_as_Y_name, paste_x_Naux_str)
      
      # (x - .data$b) * exp(.data$c) + xoffset
      define_tomeanx_t <- 
        paste0(
          "(", xsi_as_C_name, " - nlp_re_b)", 
          " .* exp(nlp_re_b)", " + ", 
          "(", "nlp_be_b", " + ", xoffset_name, ")"
        )
      
      # x / exp(.data$c) + .data$b + xoffset
      define_tomeanx_f <- 
        paste0(
          "", xsi_as_C_name, "", 
          " ./ exp(nlp_re_b)", " + ",
          "(", "nlp_be_b", " + ", xoffset_name, ")"
        )
      
      
      # y - .data$a -.data$d*if_else(.data$d.adjusted,.data$x.adj - xoffset,x) 
      # (X[ia] - knots[ja] > 0 ? X[ia] - knots[ja] : 0);
      define_tomeany_t <- 
        paste0(
          "(", ysi_as_Y_name, " - nlp_re_b)", 
          " - nlp_re_d", " .* ", 
          "(",
          # "mean(nlp_re_d) != 0 ? ", 
          "nlp_be_d_defined == 1 || nlp_re_d_defined == 1 ? ", 
          "(",
          name_xyadj_tomeanx_t, " - ", paste0("(nlp_be_b", " + ", 
                                              xoffset_name, ")"),
          ")",
          " : ",
          xsi_as_C_name , 
          ")"
        ) 
      
      
      # y + .data$a + .data$d*if_else(.data$d.adjusted,.data$x.adj - xoffset,x) 
      define_tomeany_f <- 
        paste0(
          "(", ysi_as_Y_name, " + nlp_re_b)", 
          " + nlp_re_d", " .* ", 
          "(",
          # "mean(nlp_re_d) != 0 ? ", 
          "nlp_be_d_defined == 1 || nlp_re_d_defined == 1 ? ", 
          "(",
          name_xyadj_tomeanx_f, " - ", paste0("(nlp_be_b", " + ", 
                                              xoffset_name, ")"),
          ")",
          " : ",
          xsi_as_C_name , 
          ")"
        ) 
      
      
      define_tomeany_tf_all <- 
        paste0(paste0(name_xyadj_tomeanx_t, " = ", define_tomeanx_t, ";"),
               "\n",
               paste0(name_xyadj_tomeanx_f, " = ", define_tomeanx_f, ";"), 
               "\n", 
               paste0(name_xyadj_tomeany_t, " = ", define_tomeany_t, ";"),
               "\n",
               paste0(name_xyadj_tomeany_f, " = ", define_tomeany_f, ";"))
      
      if(nys > 1) {
        gsubit <- "nlp_re"
        gsubby <- paste0(gsubit, "_", xnss)
        define_tomeany_tf_all <- gsub(gsubit, gsubby, define_tomeany_tf_all, 
                                      fixed = T)
        
        gsubit <- "nlp_be"
        gsubby <- paste0(gsubit, "_", xnss)
        define_tomeany_tf_all <- gsub(gsubit, gsubby, define_tomeany_tf_all, 
                                      fixed = T)
      }
      
      define_add_gqdata_stanvarlist_si_c <- 
        c(define_add_gqdata_stanvarlist_si_c, define_tomeany_tf_all)
    }
    
    define_add_gqdata_stanvarlist_si <- 
      paste(define_add_gqdata_stanvarlist_si_c, collapse = "\n")
    
    
    
    start_add_gqdata_stanvarlist_si <- paste0(add_gqdata_stanvarlist_si,
                                              "\n",
                                              "{")
    
    end_add_gqdata_stanvarlist_si <- paste0("\n",
                                            define_add_gqdata_stanvarlist_si,
                                            "\n",
                                            "}",
                                            "\n")
    
    
    add_start_add_gqdata_stanvarlist_si <-
      brms::stanvar(x = NULL,
                    name = 'start_add_gqdata_stanvarlist_si',
                    scode = start_add_gqdata_stanvarlist_si,
                    block = "genquant",
                    position = "end",
                    pll_args = NULL)
    
    
    add_end_add_gqdata_stanvarlist_si <-
      brms::stanvar(x = NULL,
                    name = 'end_add_gqdata_stanvarlist_si',
                    scode = end_add_gqdata_stanvarlist_si,
                    block = "genquant",
                    position = "end",
                    pll_args = NULL)
    
    
    brm_args$stanvars <- 
      brm_args$stanvars + 
      add_start_add_gqdata_stanvarlist_si + 
      add_data2stanvars + 
      add_end_add_gqdata_stanvarlist_si
    
  } # end if(genquant_xyadj) {
  
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################

  
  if(get_priors) {
    tempriorstr <- CustomDoCall(brms::get_prior, brm_args)
    return(tempriorstr)
  }

  
  if(remove_sigma_parameter) {
    brm_args[['prior']] <- brm_args[['prior']] %>% 
      dplyr::filter(class != 'sigma')
  }

  scode_final  <- CustomDoCall(brms::make_stancode, brm_args)
  sdata        <- CustomDoCall(brms::make_standata, brm_args)
  
  
  if(parameterization == 'cp') {
    scode_final <- edit_scode_ncp_to_cp_new(scode_final, 
                                        genq_only = FALSE, 
                                        normalize = normalize, 
                                        cp_via = cp_via)
  } else if(parameterization == 'ncp') {
    scode_final <- scode_final
  }
  
  
  if(select_model_edit == 'logistic3e') {
    outedit_ <- edit_scode_for_logistic3(scode_final, 
                                         normalize = normalize)
    
    scode_final <- outedit_$editedcode
  }
  
  
  
  if(!is.null(decomp)) {
    if(decomp_editcode) scode_final <- decomp_escode2(scode_final)
  }
  

  
  get_priors_eval_numeric <- TRUE
  if(get_priors_eval & get_priors_eval_numeric) {
    get_priors_eval_out <- priors_to_textdata(spriors = brm_args$prior,
                                                  sdata = sdata,
                                              raw = TRUE)
  }
  
  if(get_priors_eval & !get_priors_eval_numeric) {
    get_priors_eval_out <- brm_args$prior
  }
  
  
  
  full_custom <- FALSE
  if(is.null(initsi) | initsi == 'NULL') {
    init_custom <- init_custom
  } else if(initsi[[1]] == 'custom') {
    if(!is.null(init_custom)) full_custom <- TRUE
  }
  
  
  
  
  exe_model_fit <- TRUE
  if(get_stancode |
     get_standata |
     get_formula |
     get_stanvars |
     get_priors |
     get_priors_eval |
     validate_priors |
     get_init_eval) {
    exe_model_fit <- FALSE
  }
  
  
  
  
  if(!fit_edited_scode) {
   if(!exe_model_fit) {
      if(get_priors) {
        return(CustomDoCall(brms::get_prior, brm_args))
      } else if(get_standata) {
        return(CustomDoCall(brms::make_standata, brm_args))
      } else if(get_stancode) {
        return(scode_final)
      } else if(get_priors_eval) {
        return(get_priors_eval_out)
      } else if(validate_priors) {
        return(CustomDoCall(brms::validate_prior, brm_args))
      } else if(get_init_eval) {
        return(brm_args$init)
      } else if(get_formula) {
        return(brm_args$formula)
      } else if(get_stanvars) {
        return(brm_args$stanvars)
      }
     } # if(!fit_edited_scode) {
  } # if(!exe_model_fit) {
  
  
  # over ride exe_model_fit TRUE
  # This to not by pass below return statements
  
  fit_edited_scode_exe_model_fit <- exe_model_fit
  exe_model_fit                  <- TRUE

  

  if(exe_model_fit) {
    if(brm_args$backend == "rstan") {
      if(length(brm_args$init) == 1) {
        if(brm_args$init == "0") {
          init_custom <- NULL
        } else if(brm_args$init == "random") {
          init_custom <- NULL
        } else {
          init_custom <- init_custom
        }
      } else {
        init_custom <- init_custom
      }
    }
   

    if(brm_args$backend == "cmdstanr") {
      if(is.null(brm_args$init)) {
        init_custom <- NULL
      } else if(!is.list(brm_args$init) & length(brm_args$init) == 1) {
        if(brm_args$init == "0") {
          init_custom <- NULL
        } else if(brm_args$init == "random") {
          init_custom <- NULL
        } else if(brm_args$init == 0) {
          init_custom <- NULL
        } else {
          init_custom <- init_custom
        }
      } else if(!is.null(brm_args$init)) {
        init_custom <- init_custom
      }
    } 
    
    
    
    if(!is.null(init_custom)) {
      init_fun <- function(chain_id = 1) init_custom
      if(!is.list(init_custom[[1]])) {
        init_custom <- 
          lapply(1:brm_args$chains, function(id) init_fun(chain_id = id))
      } else if(is.list(init_custom[[1]]) & length(init_custom) == 1) {
        init_custom <- rep(init_custom, brm_args$chains)
      } else {
        if(length(init_custom) != length(brm_args$init)) {
          stop2c("Custom initials specified via 'init_custom' argument must",
               "\n ", 
               " be a single named list (e.g., custom_init = list(x= 2,xx=5)) ",
               "\n ", 
               " or else a list of list matching the number of chains")
        }
      }
      
      
      if(!full_custom) {
        new_init_append <- list()
        init_old <- brm_args$init
        init_append <- init_custom
        
        additional_init_names <- 
          setdiff(names(init_append[[1]]), names(init_old[[1]]))
        for (ilen in 1:length(init_old)) {
          additional_init <- init_append[[ilen]][additional_init_names]
          new_init_append[[ilen]] <- c(init_old[[ilen]], additional_init)
        }
        brm_args$init <- new_init_append
      } else if(full_custom) {
        brm_args$init <- init_custom
      } # if(!full_custom) {
      
    } # if(!is.null(init_custom)) {
    


    # This when all lists of list NULL (e.g., when all init args random)
    if(length(brm_args$init[lengths(brm_args$init) != 0]) == 0) {
      if(brm_args$backend == 'cmdstanr') brm_args$init <- NULL
      if(brm_args$backend == 'rstan')    brm_args$init <- 'random'
    }
 
    
    # Set refresh based on thin argument when thin > 1
    if(!is.null(brm_args$refresh) & brm_args$thin > 1) {
      brm_args$refresh <- 
        ceiling((brm_args$refresh * brm_args$thin) / brm_args$thin)
    } else {
      brm_args$refresh <- brm_args$refresh
    }
    
    # Get and evaluate file argument
    # This to save object 'file' at the end with model info
    get_file          <- brm_args$file
    get_file_refit    <- brm_args$file_refit
    get_file_compress <- brm_args$file_compress
    get_write_brmsfit <- utils::getFromNamespace("write_brmsfit", "brms")
    get_read_brmsfit  <- utils::getFromNamespace("read_brmsfit", "brms")
    
    get_file_refit_options <- 
      utils::getFromNamespace("file_refit_options", "brms")
    
    get_file_refit <- match.arg(get_file_refit, get_file_refit_options())
    if (!is.null(get_file) && get_file_refit == "never") {
      x_from_file <- get_read_brmsfit(get_file)
      if (!is.null(x_from_file)) {
        return(x_from_file)
      }
    }
    
    if (!is.null(get_file) && get_file_refit == "on_change") {
      x_from_file <- get_read_brmsfit(get_file)
      if (!is.null(x_from_file)) {
        needs_refit <- brms::brmsfit_needs_refit(
          x_from_file, scode = brms::stancode(x_from_file), 
          sdata = brms::standata(x_from_file),
          data = x_from_file$data, 
          algorithm = brm_args$algorithm, silent = brm_args$silent
        )
        if (!needs_refit) {
          return(x_from_file)
        }
      }
    }
    
    # Set it to NULL to avoid re saving later
    brm_args$file <- NULL
    
    brm_args <- sanitize_algorithm_args(args = brm_args,
                                        algorithm = brm_args$algorithm,
                                        verbose = FALSE)
    
    
    
    if(is.logical(pathfinder_init)) {
      if(!pathfinder_init) pathfinder_init <- FALSE
      if( pathfinder_init) pathfinder_init <- TRUE
    } else if(!is.logical(pathfinder_init)) {
      if(pathfinder_init == "pathfinder") {
        pathfinder_init <- TRUE
      }
    }
    
    if(brm_args$backend == "cmdstanr" | !is.null(pathfinder_args) | 
       pathfinder_init) {
      clinenumber <- getOption("cmdstanr.print_line_numbers")
      options("cmdstanr.print_line_numbers" = TRUE)
      on.exit(options("cmdstanr.print_line_numbers" = clinenumber), add = TRUE)
      
      cmaxrows <- getOption("cmdstanr_max_rows")
      options("cmdstanr_max_rows" = 20)
      on.exit(options("cmdstanr_max_rows" = cmaxrows), add = TRUE)
      
      cwarninits <- getOption("cmdstanr_warn_inits")
      options("cmdstanr_warn_inits" = FALSE)
      on.exit(options("cmdstanr_warn_inits" = cwarninits), add = TRUE)
    }
   
  
    # 23.05.2025
    # For stanc_options "01", compilation hangs for cmdstanr when thread != NULL
    # This does't happen for brms e,g,, 
    # fit6 <- brm(bf(y ~ x, sigma ~ 0 + x), data = data_het,
    #             stan_model_args = list(stanc_options = list("O1")),
    #             backend = 'cmdstanr')
    
    # Investigate it but for now set 'stanc_options' as NULL
    
    if(is.null(brm_args$threads$threads)) {
      if(brm_args$backend == "cmdstanr") {
        # brm_args$stan_model_args <- list()
        brm_args$stan_model_args$stanc_options <- NULL
      }
    }
    
    
    if(parameterization == "ncp") {
      if(sum_zero) {
        scode_final_sum_zero <- scode_final_sum_zero_int <- scode_final
        for (i in 1:100) {
          char_to_search <- paste0("z_", i)
          found_char <- grepl(char_to_search, scode_final_sum_zero, fixed = T)
          if(found_char) {
            gsub_it <- "transformed data {"
            gsub_by_add <-   paste0("  real<lower=0> sigmauz_", i, 
                                    " = sqrt(N_", i, " * inv(N_", i, 
                                    " - 1.0))", ";")
            
            gsub_by <- paste0(gsub_it, "\n", gsub_by_add)
            scode_final_sum_zero <- gsub(gsub_it, gsub_by, 
                                         scode_final_sum_zero, fixed = T)
            
            gsub_it <- paste0("matrix[M_", i, ", N_", i, "] z_", i, ";")
            gsub_by <- paste0("array[M_", i, "] sum_to_zero_vector[N_", i,
                              "] uz_", i, ";")
            scode_final_sum_zero <- gsub(gsub_it, gsub_by, 
                                         scode_final_sum_zero, fixed = T)
            gsub_it <- paste0("matrix[N_", i, ", M_", i, "] r_", i, ";")
            gsub_by_add <- paste0("  matrix[M_", i, ", N_", i, "] z_", i, ";", 
                                  "\n", 
                                  paste0("  for(i in 1:M_", i, ") z_", i, 
                                         "[i, ] = uz_", i, "[i][]'", 
                                         ";")) #note[]
            gsub_by <- paste0(gsub_it, "\n", gsub_by_add)
            scode_final_sum_zero <- gsub(gsub_it, gsub_by, 
                                         scode_final_sum_zero, fixed = T)
            if(brm_args$normalize) {
              gsub_it <- paste0("target += std_normal_lpdf(to_vector(z_", i,
                                "))", ";")
            } else {
              gsub_it <- paste0("target += std_normal_lupdf(to_vector(z_", i, 
                                "))", ";")
            }
            gsub_by <- paste0("for(i in 1:M_", i, ") target += normal_lpdf(uz_", 
                              i, 
                              "[i][] | 0.0, sigmauz_", i, ")", 
                              ";") # note []
            scode_final_sum_zero <- gsub(gsub_it, gsub_by, 
                                         scode_final_sum_zero, fixed = T)
          } # if(found_char) {
        } # for (i in 1:10) {
        
   
        
        if(is.null(brm_args$init)) {
          brm_args$init <- brm_args$init
        } else if(!is.list(brm_args$init)) {
          if(brm_args$init == 0) brm_args$init <- brm_args$init
        } else {
          for (variablexx in 1:length(brm_args$init)) {
            for (i in 1:100) {
              char_to_search <- paste0("z_", i)
              found_char <- grepl(char_to_search, 
                                  scode_final_sum_zero_int, fixed = T)
              if(found_char) {
                brm_args$init[[variablexx]] [[paste0("uz_", i)]] <- 
                  brm_args$init[[variablexx]] [[paste0("z_", i)]]
                
                brm_args$init[[variablexx]] [[paste0("z_", i)]] <- NULL
              } # if(found_char) {
            } # for (i in 1:10) {
          } # for (variablexx in 1:length(brm_args$init)) {
        }
        
        brm_args$sum_zero <-  NULL
        brms_arguments$sum_zero <-  NULL
        scode_final <- scode_final_sum_zero
      } # if(sum_zero) {
    } # if(parameterization == "ncp") {
    
    
    
    
    if(parameterization == "cp") {
      if(sum_zero) {
        scode_final_sum_zero <- scode_final_sum_zero_int <- scode_final
        for (i in 1:100) {
          char_to_search <- paste0("r_", i)
          found_char <- grepl(char_to_search, scode_final_sum_zero, fixed = T)
          if(found_char) {
            gsub_it <- "transformed data {"
            gsub_by_add <-   paste0("  real<lower=0> sigmauz_", i, 
                                    " = sqrt(N_", i, " * inv(N_", i, 
                                    " - 1.0))", ";")
            
            gsub_by <- paste0(gsub_it, "\n", gsub_by_add)
            scode_final_sum_zero <- gsub(gsub_it, gsub_by, 
                                         scode_final_sum_zero, fixed = T)
            
            gsub_it <- paste0("matrix[N_", i, ", M_", i, "] r_", i, ";")
            gsub_by <- paste0("array[M_", i, "] sum_to_zero_vector[N_", i,
                              "] uz_", i, ";")
            scode_final_sum_zero <- gsub(gsub_it, gsub_by, 
                                         scode_final_sum_zero, fixed = T)
            gsub_it <- paste0("real lprior = 0;")
            gsub_by_add <- paste0("  matrix[N_", i, ", M_", i, "] r_", i, ";", 
                                  "\n", 
                                  paste0("  for(i in 1:M_", i, ") r_", i, 
                                         "[, i] = uz_", i, "[i]", ";"),
                                  "\n",
                                  paste0("  for(i in 1:M_", i, ") ", 
                                         "lprior += log_sum_exp(uz_", i, 
                                         "[i]", ")", ";")) # note []
            gsub_by <- paste0(gsub_it, "\n", gsub_by_add)
            scode_final_sum_zero <- gsub(gsub_it, gsub_by, 
                                         scode_final_sum_zero, fixed = T)
            if(cp_via == "multi_normal_cholesky_lpdf") {
              pasetitsd <- paste0(" .* ", "rep_vector(sigmauz_", i, ", 
                                  ", "M_", i, ")")
              gsub_it <- paste0("diag_pre_multiply(sd_", i)
              gsub_by <- paste0(gsub_it, pasetitsd)
            } else if(cp_via == "multi_normal_lpdf") {
              pasetitsd <- paste0(" .* ", "rep_vector(sigmauz_", i, 
                                  ", ", "M_", i, ")")
              gsub_it <- 
                paste0("quad_form_diag(multiply_lower_tri_self_transpose(L_", 
                                i, "), ", "sd_", i)
              gsub_by <- paste0(gsub_it, pasetitsd)
            }
            scode_final_sum_zero <- gsub(gsub_it, gsub_by, scode_final_sum_zero, 
                                         fixed = T)
          } # if(found_char) {
        } # for (i in 1:10) {
        
  
        
        if(is.null(brm_args$init)) {
          brm_args$init <- brm_args$init
        } else if(!is.list(brm_args$init)) {
          if(brm_args$init == 0) brm_args$init <- brm_args$init
        } else {
          for (variablexx in 1:length(brm_args$init)) {
            for (i in 1:100) {
              char_to_search <- paste0("r_", i)
              found_char <- grepl(char_to_search, 
                                  scode_final_sum_zero_int, fixed = T)
              if(found_char) {
                if(!is.null(brm_args$init[[variablexx]] [[char_to_search]])) {
                  brm_args$init[[variablexx]] [[paste0("uz_", i)]] <-
                    brm_args$init[[variablexx]] [[char_to_search]] %>% t()
                }
                brm_args$init[[variablexx]] [[char_to_search]] <- NULL
                brm_args$init[[variablexx]] [[paste0("z_", i)]] <- NULL
              } # if(found_char) {
            } # for (i in 1:10) {
          } # for (variablexx in 1:length(brm_args$init)) {
        }
        
        brm_args$sum_zero <-  NULL
        brms_arguments$sum_zero <-  NULL
        scode_final <- scode_final_sum_zero
      } # if(sum_zero) {
    } # if(parameterization == "cp") {
    
    
    
    # add_rescor_by
    if (set_rescor_by) {
      scode_final <- 
        edit_stancode_for_multivariate_rescor_by(stan_code = scode_final, 
                                    threads = brm_args$threads$threads,
                                    normalize = brm_args$normalize,
                                    corr_method = Rescor_method) 
    }
    
    
    # add_sigma_by_mu or fitted(.)
    if(set_model_sigma_by_fz |
       set_model_sigma_by_fp |
       set_model_sigma_by_fe |
       set_model_sigma_by_ve | 
       set_model_sigma_by_vp |
       set_model_sigma_by_cp |
       set_model_sigma_by_mp |
       set_model_sigma_by_me |
       set_model_sigma_by_rp |
       set_model_sigma_by_re ) {
      ithx <- 0
      for (outrespbames in ys) {
        ithx <- ithx + 1
        if(nys == 1) {
          base_forms   <- bformula
          base_pforms  <- base_forms[['pforms']]
          sigma_forms  <- base_pforms[['sigma']]
        } else if(nys > 1) {
          base_forms   <- bformula[['forms']][[outrespbames]]
          base_pforms  <- base_forms[['pforms']]
          sigma_forms  <- base_pforms[['sigma']]
        }
        sigma_has_loop <- attr(sigma_forms, "loop")
        
        if(nys == 1) {
          if(sigma_has_loop) {
            gsub_it_start     <- paste0("sigma", "", "", "[n]")
            onlum <- paste0("mu", "", "", "[n]")
          } else {
            # this " =" restricts to the relevant portion of the string
            gsub_it_start     <- paste0("sigma", "", "", " =") 
            onlum <- paste0("mu", "", "", "")
          }
        } else if(nys > 1) {
          if(sigma_has_loop) {
            gsub_it_start     <- paste0("sigma", "_", outrespbames, "[n]")
            onlum <- paste0("mu", "_", outrespbames, "[n]")
          } else {
            # this " =" restrict to relevant portion
            gsub_it_start     <- paste0("sigma", "_", outrespbames, " =") 
            onlum <- paste0("mu", "_", outrespbames, "")
          }
        }
        
        
        set_model_sigma_by_mu_fun_str <- set_model_sigma_by_mu_fun_str_c[[ithx]]
        sigmatau_strsi <- sigmatau_strsi_c[[ithx]]
        
        # best, just keep ()
        if(set_model_sigma_by_fz |
           set_model_sigma_by_fp |
           set_model_sigma_by_fe |
           set_model_sigma_by_mp |
           set_model_sigma_by_me) {
          gsub_it_end <- "()"
          if(!grepl(gsub_it_end, set_model_sigma_by_mu_fun_str, fixed = T)) {
            stop2c("For 'fitted' and 'mean' var function, the predictor must ",
                 "be a set of empty paranthesis '()' or 'identity()'",
                 "\n  ", 
                 "Please check the following relevant part of the code",
                 "\n  ", 
                 set_model_sigma_by_mu_fun_str)
          }
        } else {
          gsub_it_end <- ")"
        }
        
        extract_sigma_by_mean_o <- replace_string_part(x = scode_final,
                                                       start = gsub_it_start,
                                                       end =  gsub_it_end,
                                                       replace = "",
                                                       extract = T)
        
        if(set_model_sigma_by_fz |
           set_model_sigma_by_fp |
           set_model_sigma_by_fe |
           set_model_sigma_by_mp |
           set_model_sigma_by_me |
           set_model_sigma_by_rp |
           set_model_sigma_by_re) {
          extract_sigma_by_mean <- gsub(paste0("", 
                                               "(", ""  ,")"),
                                        paste0("", 
                                               "(", onlum  ,")"),
                                        extract_sigma_by_mean_o, fixed = T)
        } else {
          extract_sigma_by_mean <- extract_sigma_by_mean_o
        }
        
        if(!sigma_has_loop) {
          extract_sigma_by_mean <- gsub(paste0("", " ", "*"), 
                                        paste0("", " ", ".*"), 
                                        extract_sigma_by_mean, fixed = T)
        }
        
        scode_final <- gsub(extract_sigma_by_mean_o, extract_sigma_by_mean,
                            scode_final, fixed = T)
      } # for (outrespbames in ys) {
    } # if(set_model_sigma_by_fz | ...) {
    
    

    if(!is.null(custom_formula)) {
      if(!brms::is.brmsformula(custom_formula)) {
        stop2c("The 'custom_formula' must be of class 'brmsformula'")
      }
      brm_args[['formula']] <- custom_formula
    }
    
    if(!is.null(custom_prior)) {
      if(!brms::is.brmsprior(custom_prior)) {
        stop2c("The 'custom_prior' must be of class 'brmsprior'")
      }
      brm_args[['prior']] <- custom_prior
    }
    
   
    # stop()
    # do.call(brms::make_stancode, brm_args)
    # fit_edited_scode <- TRUE
    # replace_it <- "ptarget += normal_lpdf(Y[start:end] | mu, sigma)"
    # replace_by <-
    #   "ptarget += normal_lpdf(Y[start:end] | mu, sigma) +
    #   normal_lccdf(10 | mu, sigma) +
    # normal_lcdf(300 | mu, sigma);"
    # 
    #  scode_final <- bsitar:::replace_string_part(scode_final, replace_it, ";", replace_by)
    # 
    # print(scode_final)
    # 
    # do.call(brms::make_stancode, brm_argsxx) 
    
    
    # This for model_infor
    brm_args_prior <- brm_args$prior
    
    if(!fit_edited_scode_exe_model_fit & fit_edited_scode) {
      if(get_priors) {
        return(brm_args$prior)
      } else if(get_standata) {
        return(sdata)
      } else if(get_stancode) {
        return(scode_final)
      } else if(get_priors_eval) {
        return(brm_args$prior)
      } else if(validate_priors) {
        return(CustomDoCall(brms::validate_prior, brm_args))
      } else if(get_init_eval) {
        return(brm_args$init)
      } else if(get_formula) {
        return(brm_args$formula)
      } else if(get_stanvars) {
        return(brm_args$stanvars)
      }
    } 
   
    
    if(fit_edited_scode) {
      if(verbose) message2c("Fitting model via edited stancode...")
      if(brm_args$backend == "cmdstanr") {
         brmsfit <- brms_via_cmdstanr(scode_final, sdata,  
                                      brm_args, brms_arguments,
                                      pathfinder_args = pathfinder_args,
                                      pathfinder_init = pathfinder_init,
                                      Rescor_by_levels = Rescor_by_levels,
                                      verbose = verbose)
      }
      if(brm_args$backend == "rstan") {
        brmsfit  <- brms_via_rstan(scode_final, sdata,
                                   brm_args, brms_arguments,
                                   Rescor_by_levels = Rescor_by_levels,
                                   verbose = verbose)
      }
    } else if(!fit_edited_scode) {
      if(!is.null(pathfinder_args) | pathfinder_init) {
        brmsfit <- brms_via_cmdstanr(scode_final, sdata, brm_args,
                                     brms_arguments,
                                     pathfinder_args = pathfinder_args,
                                     pathfinder_init = pathfinder_init,
                                     Rescor_by_levels = Rescor_by_levels,
                                     verbose = verbose)
      } else {
        brmsfit <- CustomDoCall(brms::brm, brm_args)
      }
    } # if(fit_edited_scode) {
    
    if(brm_args$backend == "mock") {
      brmsfit <- CustomDoCall(brms::brm, brm_args)
    }
    
    # Add class attributes and the model info for post-processing
    attr(brmsfit, 'class') <- c(attr(brmsfit, 'class'), 'bgmfit')
    
    ##############################################################
    ##############################################################
    # restore sigma sqrt form
    if(set_model_sigma_by_fz |
       set_model_sigma_by_fp |
       set_model_sigma_by_fe |
       set_model_sigma_by_mp |
       set_model_sigma_by_me |
       set_model_sigma_by_rp |
       set_model_sigma_by_re) {
      function_restore_mu_sigam_form <- function(model, 
                                                 ys, 
                                                 set_model_sigma_by_mu_fun_str_c
                                                 ) {
        ithx <- 0
        for (outrespbames in ys) {
          ithx <- ithx + 1
          if(nys == 1) {
            base_forms   <- bformula
            base_pforms  <- base_forms[['pforms']]
            sigma_forms  <- base_pforms[['sigma']]
          } else if(nys > 1) {
            base_forms   <- bformula[['forms']][[outrespbames]]
            base_pforms  <- base_forms[['pforms']]
            sigma_forms  <- base_pforms[['sigma']]
          }
          base_mu_fun  <- base_forms[['formula']][[3]]
          base_mu_fun  <- base_mu_fun %>% deparse()
          edit_attr    <- attributes(sigma_forms)
          str_edit     <- sigma_forms
          str_edit     <- str_edit %>% deparse()
          set_model_sigma_by_mu_fun_str <- 
            set_model_sigma_by_mu_fun_str_c[[ithx]]
          str_edit <- gsub(paste0("", 
                                    "(", ""  ,")"), 
                             paste0("", 
                                    "(", base_mu_fun  ,")"),
                           str_edit, fixed = T)
          
          sigma_forms <- str2lang(str_edit)
          attributes(sigma_forms) <- edit_attr
          if(nys == 1) {
            model[['formula']][['pforms']][['sigma']] <- sigma_forms
          } else if(nys > 1) {
        model[['formula']][['forms']][[outrespbames]][['pforms']][['sigma']] <- 
              sigma_forms
          }
        } 
        return(model[['formula']])
      } 
      
      function_restore_mu_sigam_form_new <- 
        function_restore_mu_sigam_form(model = brmsfit,
                                       ys = ys, 
                                       set_model_sigma_by_mu_fun_str_c = 
                                         set_model_sigma_by_mu_fun_str_c)
      brmsfit$formula <- function_restore_mu_sigam_form_new
    } # if(set_model_sigma_by_fi | ...) {
    

    
    ##############################################################
    ##############################################################
    
    model_info <- list()
    
    model_info[['fit_edited_scode']]  <- fit_edited_scode
    
    if(fit_edited_scode) {
      model_info[['emodel']]           <- scode_final
    }
    
    model_info[['parameterization']] <- parameterization
    model_info[['d_adjusted']]       <- d_adjusted
    
    for (i in 1:length(funlist_rnamelist)) {
      model_info[[funlist_rnamelist[[i]]]] <- funlist_rvaluelist[[i]]
    }
    
    for (i in 1:length(include_fun_nameslist_rnamelist)) {
      model_info[[include_fun_nameslist_rnamelist[[i]]]] <- 
        include_fun_nameslist_rvaluelist[[i]]
    }
    
    for (i in 1:length(xoffsetnamelist)) {
      model_info[[xoffsetnamelist[[i]]]] <- xoffsetvaluelist[[i]]
    }
    
    for (i in 1:length(knotsnamelist)) {
      model_info[[knotsnamelist[[i]]]] <- knotsvaluelist[[i]]
    }
    
    for (i in 1:length(fixednamelist)) {
      model_info[[fixednamelist[[i]]]] <- fixedvaluelist[[i]]
    }
    
    for (i in 1:length(randomnamelist)) {
      model_info[[randomnamelist[[i]]]] <- randomvaluelist[[i]]
    }
    
    for (i in 1:length(xfunnamelist)) {
      model_info[[xfunnamelist[[i]]]] <- xfunvaluelist[[i]]
    }
    
    for (i in 1:length(yfunnamelist)) {
      model_info[[yfunnamelist[[i]]]] <- yfunvaluelist[[i]]
    }
   
    for (i in 1:length(groupvarnamelist)) {
      model_info[[groupvarnamelist[[i]]]] <- groupvarvaluelist[[i]]
    }
    
    for (i in 1:length(hierarchicalvarnamelist)) {
      model_info[[hierarchicalvarnamelist[[i]]]] <- 
        hierarchicalvarvaluelist[[i]]
    }
    
    for (i in 1:length(xnamelist)) {
      model_info[[xnamelist[[i]]]] <- xvarvaluelist[[i]]
    }
    
    for (i in 1:length(sigmaxnamelist)) {
      model_info[[sigmaxnamelist[[i]]]] <- sigmaxvarvaluelist[[i]]
    }
    
    for (i in 1:length(ynamelist)) {
      model_info[[ynamelist[[i]]]] <- yvarvaluelist[[i]]
    }
    
    for (i in 1:length(covnamelist)) {
      model_info[[covnamelist[[i]]]] <- covvaluelist[[i]]
    }
    
    if(!is.na(univariate_by$by)) {
      model_info[['subindicators']] <- subindicators
    } 
    
    for (i in 1:length(d_adjustednamelist)) {
      model_info[[d_adjustednamelist[[i]]]] <- d_adjustedvaluelist[[i]]
    }
    
    for (i in 1:length(xfuntransformnamelist)) {
      model_info[[xfuntransformnamelist[[i]]]] <- xfuntransformvaluelist[[i]]
    }
    
    for (i in 1:length(xfuntransform2namelist)) {
      model_info[[xfuntransform2namelist[[i]]]] <- xfuntransform2valuelist[[i]]
    }
    
    for (i in 1:length(yfuntransformnamelist)) {
      model_info[[yfuntransformnamelist[[i]]]] <- yfuntransformvaluelist[[i]]
    }
    
    # Inverse funs are created internally
    for (i in 1:length(ixfuntransformnamelist)) {
      model_info[[ixfuntransformnamelist[[i]]]] <- ixfuntransformvaluelist[[i]]
    }
    
    for (i in 1:length(ixfuntransform2namelist)) {
      model_info[[ixfuntransform2namelist[[i]]]] <- ixfuntransform2valuelist[[i]]
    }

    for (i in 1:length(iyfuntransformnamelist)) {
      model_info[[iyfuntransformnamelist[[i]]]] <- iyfuntransformvaluelist[[i]]
    }
    
    for (i in 1:length(SplineCallnamelist)) {
      model_info[[SplineCallnamelist[[i]]]] <- SplineCallvaluelist[[i]]
    }
    

    model_info[['StanFun_name']]  <- spfncname_common
    model_info[['multivariate']]  <- multivariate
    model_info[['univariate_by']] <- univariate_by
    model_info[['nys']] <- nys
    model_info[['dfs']] <- dfs
    model_info[['xfuns_user']] <- xfuns_user
    model_info[['yfuns_user']] <- yfuns_user
    model_info[['outliers']] <- outliers
    model_info[['bgmfit.data']] <- data.org.in
    model_info[['select_model']] <- select_model
    model_info[['decomp']] <- decomp
    model_info[['fun_scode']] <- fun_scode
    model_info[['envir']] <- enverr.
    
    # The brms_arguments_list required in update_model()
    model_info[['brms_arguments_list']] <- brms_arguments_list
    
    # Full call constructed
    model_info[['call.full.bgmfit']] <- call.full
    # The call by user
    model_info[['call.bgmfit']] <- mcall_
    
    if(!exists('sigmaspfncname_common')) sigmaspfncname_common <- NULL
    if(!exists('sigmaselect_model')) sigmaselect_model <- NULL
    if(!exists('sigmadecomp')) sigmadecomp <- NULL
    if(length(sigmafunlist_rvaluelist) == 0) {
      sigmafunlist_rvaluelist <- as.list(rep("", nys))
    }
    
    if(length(sigmavarfunlist_rvaluelist) == 0) {
      sigmavarfunlist_rvaluelist <- as.list(rep("", nys))
    }
    
    if(!exists('sigmavarspfncname_common')) {
      sigmavarspfncname_common <- NULL
    }
    
    model_info[['sigmavarStanFun_name']] <- sigmavarspfncname_common 
    
    model_info[['sigmaStanFun_name']]    <- sigmaspfncname_common 
    model_info[['sigmaxs']] <- sigmaxs
    model_info[['sigmaids']] <- sigmaids
    model_info[['sigmadfs']] <- sigmadfs
    model_info[['sigmaxfuns_user']] <- sigmaxfuns_user
    model_info[['sigmaselect_model']] <- sigmaselect_model
    model_info[['sigmadecomp']] <- sigmadecomp
    model_info[['sigmad_adjusted']] <- sigmad_adjusted
    
    for (i in 1:length(sigmad_adjustednamelist)) {
      model_info[[sigmad_adjustednamelist[[i]]]] <- 
        sigmad_adjustedvaluelist[[i]]
    }
    for (i in 1:length(sigmafunlist_rnamelist)) {
      model_info[[sigmafunlist_rnamelist[[i]]]] <- 
        sigmafunlist_rvaluelist[[i]]
    }
    
    for (i in 1:length(sigmavarfunlist_rnamelist)) {
      model_info[[sigmavarfunlist_rnamelist[[i]]]] <- 
        sigmavarfunlist_rvaluelist[[i]]
    }
    
    for (i in 1:length(sigmaxfunnamelist)) {
      model_info[[sigmaxfunnamelist[[i]]]] <- sigmaxfunvaluelist[[i]]
    }
    
    for (i in 1:length(sigmacovnamelist)) {
      model_info[[sigmacovnamelist[[i]]]] <- sigmacovvaluelist[[i]]
      # drop NA
      model_info[[sigmacovnamelist[[i]]]] <- 
        model_info[[sigmacovnamelist[[i]]]][!is.na(
          model_info[[sigmacovnamelist[[i]]]])]
    }
    
    for (i in 1:length(sigmaxfuntransformnamelist)) {
      model_info[[sigmaxfuntransformnamelist[[i]]]] <- 
        sigmaxfuntransformvaluelist[[i]]
    }
    
    for (i in 1:length(sigmaxfuntransform2namelist)) {
      model_info[[sigmaxfuntransform2namelist[[i]]]] <- 
        sigmaxfuntransform2valuelist[[i]]
    }
    
    for (i in 1:length(sigmaixfuntransformnamelist)) {
      model_info[[sigmaixfuntransformnamelist[[i]]]] <- 
        sigmaixfuntransformvaluelist[[i]]
    }
    
    for (i in 1:length(sigmaixfuntransform2namelist)) {
      model_info[[sigmaixfuntransform2namelist[[i]]]] <- 
        sigmaixfuntransform2valuelist[[i]]
    }
    
    for (i in 1:length(sigmaxoffsetnamelist)) {
      model_info[[sigmaxoffsetnamelist[[i]]]] <- 
        sigmaxoffsetvaluelist[[i]]
    }
    
    for (i in 1:length(setsigmaxvarnamelist)) {
      model_info[[setsigmaxvarnamelist[[i]]]] <-
        setsigmaxvarvaluelist[[i]]
    }
    
    for (i in 1:length(sigmamodelnamenamelist)) {
      model_info[[sigmamodelnamenamelist[[i]]]] <-
        sigmamodelnamevaluelist[[i]]
    }
    model_info[['sigmamodel_all']] <- unlist(sigmamodelnamevaluelist)
    
    
    if(set_model_sigma_by_ba) {
      for (i in 1:length(sigmabasicfunlist_rnamelist)) {
        model_info[[sigmabasicfunlist_rnamelist[[i]]]] <- 
          sigmabasicfunlist_rvaluelist[[i]]
      }
      for (i in 1:length(sigmabasicfunnamenamelist)) {
        model_info[[sigmabasicfunnamenamelist[[i]]]] <-
          sigmabasicfunnamevaluelist[[i]]
      }
      for (i in 1:length(sigmabasicfunattrnamelist)) {
        model_info[[sigmabasicfunattrnamelist[[i]]]] <-
          sigmabasicfunattrvaluelist[[i]]
      }
    } # if(set_model_sigma_by_ba) {
    
    ##############################################################
    ##############################################################
    # these paste0(..., 's') will be combined across ys
    model_info[[xvar_names]]                <- xvar_names_val
    model_info[[yvar_names]]                <- yvar_names_val
    model_info[[idvar_names]]               <- idvar_names_val
    model_info[[sigmaidvar_names]]          <- sigmaidvar_names_val
    model_info[[sigmaxvar_names]]           <- sigmaxvar_names_val
    model_info[[cov_names]]                 <- cov_names_val
    model_info[[sigmacov_names]]            <- sigmacov_names_val
    model_info[[xfun_names]]                <- xfun_names_val
    model_info[[yfun_names]]                <- yfun_names_val
    model_info[[sigmaxfun_names]]           <- sigmaxfun_names_val
    
    model_info[[xfuntransform_names]]       <- xfuntransform_names_val
    model_info[[xfuntransform2_names]]      <- xfuntransform2_names_val
    model_info[[ixfuntransform_names]]      <- ixfuntransform_names_val
    model_info[[ixfuntransform2_names]]     <- ixfuntransform2_names_val
    
    model_info[[yfuntransform_names]]       <- yfuntransform_names_val
    model_info[[iyfuntransform_names]]      <- iyfuntransform_names
    
    model_info[[sigmaxfuntransform_names]]  <- sigmaxfuntransform_names_val
    model_info[[sigmaxfuntransform2_names]] <- sigmaxfuntransform2_names_val
    model_info[[sigmaixfuntransform_names]] <- sigmaixfuntransform_names_val
    model_info[[sigmaixfuntransform2_names]]<- sigmaixfuntransform2_names_val
    model_info[[xoffset_names]]             <- xoffset_names_val
    model_info[[sigmaxoffset_names]]        <- sigmaxoffset_names_val
    
    model_info[[setsigmaxvar_names]]        <- setsigmaxvar_names_val
    model_info[['genquant_xyadj']]          <- genquant_xyadj
    
    model_info[['prior']]                   <- brm_args_prior
   
    
    ##############################################################
    ##############################################################
    
    brmsfit$model_info           <- model_info
    environment(brmsfit$formula) <- enverr.
    

    # Now message moved to the expose_model_functions()
    if (expose_function & !brm_args$empty) {
      # if (verbose) {
      #   setmsgtxt <-
      #     paste0("\n Exposing Stan functions for post-processing\n")
      #   if (displayit == 'msg') {
      #     message2c(setmsgtxt)
      #   } else if (displayit == 'col') {
      #     col <- setcolh
      #     cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
      #   }
      # }
      # if (!verbose) {
      #   setmsgtxt <-
      #     paste0("\n Exposing Stan functions for post-processing..\n")
      #   message2c(setmsgtxt)
      # }
      
      brmsfit <- expose_model_functions(model = brmsfit, 
                                      scode = fun_scode,
                                      expose = TRUE, 
                                      select_model = NULL,
                                      returnobj = TRUE,
                                      vectorize = FALSE,
                                      verbose = TRUE,
                                      # sigmafun = expose_sigma_ls_fun,
                                      backend = FALSE,
                                      path = NULL,
                                      envir = NULL)
      brmsfit$model_info[['expose_method']] <- 'S'
    } 
    
    # if (!expose_function) {
    if (!expose_function & !brm_args$empty) {
      brmsfit <- expose_model_functions(model = brmsfit,
                                      scode = fun_scode,
                                      expose = FALSE,
                                      select_model = select_model,
                                      returnobj = TRUE,
                                      vectorize = FALSE,
                                      verbose = TRUE,
                                      # sigmafun = expose_sigma_ls_fun,
                                      backend = FALSE,
                                      path = NULL,
                                      envir = NULL)
      brmsfit$model_info[['expose_method']] <- 'R'
    } 
    
    if (verbose) {
      setmsgtxt <- paste0("\nModel Fitting complete")
      if (displayit == 'msg') {
        message2c(setmsgtxt)
      } else if (displayit == 'col') {
        col <- setcolh
        cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
      }
    }
    
    if (!is.null(get_file)) {
      brmsfit <- get_write_brmsfit(brmsfit, get_file, 
                               compress = get_file_compress)
    }
   
    # 20.03.2025
    # This needed for insight::get_data
    attr(brmsfit$data, "data_name") <- data_name_str_attr
    return(brmsfit)
  } # exe_model_fit
  
  
} # End bsitar()


