


#' @title Fit Bayesian SITAR Growth Curve Model
#'
#' @description The \strong{bsitar()} function fits the Bayesian version of the
#'   Super Imposition by Translation and Rotation (\emph{SITAR}) model. The
#'   \emph{SITAR} model is a nonlinear mixed-effects model that has been widely
#'   used to summarize growth processes (such as height and weight) from early
#'   childhood through adulthood.
#'
#'   The frequentist version of the \emph{SITAR} model can be fit using the
#'   already available R package, \pkg{sitar} \insertCite{R-sitar}{bsitar}.
#'   However, the \pkg{bsitar} package offers an enhanced Bayesian
#'   implementation that improves modeling capabilities. In addition to
#'   univariate analysis (i.e., modeling a single outcome), \pkg{bsitar} also
#'   supports:
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
#' Like the \pkg{sitar} package \insertCite{Cole2010}{bsitar}, the \pkg{bsitar}
#' package fits the \emph{SITAR} model with (usually) three random effects: size
#' (parameter \code{a}), timing (parameter \code{b}), and intensity (parameter
#' \code{c}). Additionally, there is a slope parameter (parameter \code{d}) that
#' models the variability in the adult slope of the growth curve (see
#' [sitar::sitar()] for details).
#' 
#' Note that the author of the \pkg{sitar} package \insertCite{Cole2010}{bsitar}
#' enforces the inclusion of the \code{d} parameter as a random effect only,
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
#' The \pkg{sitar} package internally depends on the \pkg{brms} 
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
#'  - Any real number (e.g., \code{xoffset = 12}).
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
#'  - \code{'rcs'}: Constructs the spline design matrix using the truncated 
#'  power basis (Harrell's method), implemented in [Hmisc::rcspline.eval()].
#'  - \code{'nks'}: Implements a B-spline based natural cubic spline method, 
#'  similar to [splines2::nsk()].
#'  - \code{'nsp'}: Implements a B-spline based natural cubic spline method, 
#'  similar to [splines2::nsp()].
#' The default is \code{'nsp'}.
#' 
#' The \code{'rcs'} method uses a truncated power basis, whereas \code{'nks'}
#' and \code{'nsp'} are B-spline-based methods. Unlike [splines2::nsp()] and
#' [splines2::nsk()], which normalize the spline basis by default, \code{'nks'}
#' and \code{'nsp'} return the non-normalized version of the spline. If
#' normalization is desired, the user can specify \code{normalize = TRUE} in a
#' list. For example, to use a normalized \code{'nsp'}, one can specify
#' \code{stype = list(type = 'nsp', normalize = TRUE)}.
#' 
#' For more details, see [Hmisc::rcspline.eval()], [splines2::nsk()], and
#' [splines2::nsp()].
#' 
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
#'  the \code{subset()} formulation used in fitting \code{univariate_by} models.
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
#' where \code{id} specifies the group identifier and \code{un} sets the
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
#'   for the random effects \code{a}, \code{b}, and \code{c} specified as
#'   \code{a_formula_gr = ~1}, \code{b_formula_gr = ~1} and \code{c_formula_gr =
#'   ~1}. To specify the group identifier (e.g., \code{id}) and an unstructured
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
#'   Users can specify an external function, such as \code{poly}, but only with
#'   a single argument (the predictor), i.e., \code{poly(age)}. Additional
#'   arguments must be specified externally. For example, to set the degree of
#'   the polynomial to 3, a copy of the \code{poly} function can be created and
#'   modified as follows: \cr \code{mypoly = poly; formals(mypoly)[['degree']]
#'   <- 3; mypoly(age)}.
#'
#' @param sigma_formula_gr Formula for the random effect parameter, \code{sigma}
#'   (default \code{NULL}). See \code{a_formula_gr} for details. Similar to
#'   \code{sigma_formula}, external functions such as \code{poly} can be used.
#'   For further details, please refer to the description of
#'   \code{sigma_formula}.
#' 
#' @param sigma_formula_gr_str Formula for the random effect parameter,
#'   \code{sigma}, when fitting a hierarchical model with three or more levels
#'   of hierarchy. See \code{a_formula_gr_str} for details. As with
#'   \code{sigma_formula}, external functions such as \code{poly} can be used.
#'   For further details, please refer to the description of
#'   \code{sigma_formula}.
#'  
#' @param sigma_formula_manual Formula for the random effect parameter,
#'  \code{sigma}, provided as a character string that explicitly uses the
#'  [brms::nlf()] and [brms::lf()] functions (default \code{NULL}). An example
#'  is: \cr
#'  \code{nlf(sigma ~ z) + lf(z ~ 1 + age + (1 + age |55| gr(id, by = NULL)))}.
#'   
#'  Another use case for \code{sigma_formula_manual} is modeling a
#'  location-scale model in the \code{SITAR} framework, where the same
#'  \code{SITAR} formula can be used to model the scale (\code{sigma}). An
#'  example is: \cr
#'   
#'  \code{nlf(sigma ~ sigmaSITARFun(logage, sigmaa, sigmab, sigmac, sigmas1,
#'  sigmas2, sigmas3, sigmas4), loop = FALSE) +
#'  lf(sigmaa ~ 1 + (1 |110| gr(id, by = NULL))+(1 |330| gr(study, by = NULL))) +
#'  lf(sigmab ~ 1 + (1 |110| gr(id, by = NULL))+(1 |330| gr(study, by = NULL))) +
#'  lf(sigmac ~ 1 + (1 |110| gr(id, by = NULL))+(1 |330| gr(study, by = NULL))) +
#'  lf(sigmas1 + sigmas2 + sigmas3 + sigmas4 ~ 1)}.
#'   
#'  Here, \code{sigmaSITARFun} (and all other required sub-functions) are
#'  defined through the \code{sigmax}, \code{sigmadf}, \code{sigmaknots},
#'  \code{sigmafixed}, \code{sigmarandom}, \code{sigmaxoffset},
#'  \code{sigmaxfun}, and \code{sigmabound} arguments. Ensure the
#'  \code{sigma_formula_manual} code matches the \code{sigmaSITARFun} function
#'  created by these arguments.
#' 
#'   Note that for \code{sigma_formula_manual}, priors must be set up manually
#'   using the \code{add_self_priors} argument. To see which priors are
#'   required, the user can run the code with \code{get_priors = TRUE}.
#'   Additionally, no initial values are defined, so initial values for these
#'   parameters should be set to either \code{0} or \code{random}.
#'   
#' @param sigmax Predictor for the distributional parameter \code{sigma}. See
#'   \code{x} for details. Ignored if \code{sigma_formula_manual = NULL}.
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
#' @param sigmaxfun Transformation function for \code{x} in the \code{sigma}
#'   structure. See \code{xfun} for details. Ignored if
#'   \code{sigma_formula_manual = NULL}.
#'
#' @param sigmabound Bounds for the \code{x} in the \code{sigma} structure. See
#'   \code{bound} for details. Ignored if \code{sigma_formula_manual = NULL}.
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
#' @param multivariate Set up the multivariate model fitting (default \code{NULL}) 
#'   using a named list with the following elements:
#'   \itemize{
#'   \item \code{mvar} (logical, default \code{FALSE}): Indicates whether to fit
#'   a multivariate model.
#'   \item \code{cor} (optional, character string): Specifies the correlation
#'   structure. Available options are:
#'   \itemize{
#'   \item \code{un} (default): Models a full unstructured correlation, where
#'   group-level random effects across response variables are drawn from a joint
#'   multivariate normal distribution with shared correlation parameters. \item
#'   \code{diagonal}: Estimates only the variance parameters for each sub-model,
#'   with the correlation parameters set to zero. \item \code{un_s}: Estimates
#'   unstructured variance-covariance parameters separately for each response
#'   variable.
#'   }
#'   \item \code{rescor} (logical, default \code{TRUE}): Indicates whether to
#'   estimate the residual correlation between response variables.
#'   }
#'   
#' @param a_prior_beta Specify priors for the fixed effect parameter, \code{a}.
#'   (default \code{normal(lm, ysd, autoscale = TRUE)}). The following key
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
#'   The default prior is \code{normal(0, 1.5, autoscale = FALSE)}. For full
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
#'   The default prior is \code{normal(0, 0.5, autoscale = FALSE)}. For full
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
#'   details on prior specification, please refer to \code{a_prior_beta}.
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
#'   (i.e., spline coefficients). The default prior is \code{normal(0, 'lm',
#'   autoscale = TRUE)}. The general approach is similar to the one described
#'   for other fixed effect parameters (see \code{a_prior_beta} for details).
#'   Key points to note:
#'   \itemize{
#'   \item When using location-scale based priors with 'lm' (e.g.,
#'   \code{s_prior_beta = normal(lm, 'lm')}), the location parameter is set from
#'   the spline coefficients obtained from the simple linear model fit, and the
#'   scale parameter is based on the standard deviation of the spline design
#'   matrix. The location parameter is typically set to 0 (default), and
#'   \code{autoscale} is set to \code{TRUE}. \item For location-scale based
#'   priors, the option \code{sethp} (logical, default \code{FALSE}) is
#'   available to define hierarchical priors. Setting \code{sethp = TRUE} alters
#'   the prior setup to use hierarchical priors: \code{s ~ normal(0, 'lm')}
#'   becomes \code{s ~ normal(0, 'hp')}, where \code{'hp'} is defined as
#'   \code{hp ~ normal(0, 'lm')}. The scale for the hierarchical prior is
#'   automatically taken from the \code{s} parameter, and it can also be defined
#'   using the same \code{sethp} option. For example, \code{s_prior_beta =
#'   normal(0, 'lm', sethp = cauchy)} will result in \code{s ~ normal(0, 'lm')},
#'   \code{hp ~ cauchy(0, 'lm')
#'   }.
#'   \item For \code{uniform} priors, you can use the option \code{addrange} to
#'   symmetrically expand the prior range.
#'   }
#'   It is observed that location-scale based prior distributions (such as
#'   \code{normal}, \code{student_t}, and \code{cauchy}) typically work well for
#'   spline coefficients.
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
#'   (default \code{normal(0, 1.0, autoscale = FALSE)}). See \code{a_prior_sd}
#'   for details.
#' 
#' @param c_prior_sd Specify priors for the random effect parameter, \code{c}.
#'   (default \code{normal(0, 0.25, autoscale = FALSE)}). See \code{a_prior_sd}
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
#'    \item \code{'0'}: All parameters are initialized to zero. \item
#'    \code{'random'}: \strong{Stan} randomly generates initial values for each
#'    parameter within a range defined by \code{init_r} (see below), or between
#'    -2 and 2 in unconstrained space if \code{init_r = NULL}. \item
#'    \code{'prior'}: Initializes parameters based on the specified prior. \item
#'    \code{NULL} (default): Initial values are provided by the corresponding
#'    init arguments defined below.
#'  }
#'
#' @param init_r A positive real value that defines the range for the random
#'   generation of initial values (default \code{NULL}). This argument is used
#'   only when \code{init = 'random'}.
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
#'  }
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
#' @param vcov_init_0 A logical (default \code{FALSE}) to set initial values for
#'   variance (standard deviation) and covariance (correlation) parameters to
#'   zero. This allows for setting custom initial values for the fixed effects
#'   parameters while keeping the variance-covariance parameters at zero.
#' 
#' @param jitter_init_beta A proportion (between 0 and 1) to perturb the initial
#'   values for fixed effect parameters. The default is \code{NULL}, which means
#'   that the same initial values are used across all chains. A sensible option
#'   might be \code{jitter_init_beta = 0.1}, which mildly perturbs the initial
#'   values. Note that the jitter is applied as a proportion of the specified
#'   initial value, not an absolute amount. For example, if the initial value is
#'   \code{100}, setting \code{jitter_init_beta = 0.1} means the perturbed
#'   initial value will be within the range \code{90} to \code{110}. Conversely,
#'   if the initial value is \code{10}, the perturbed value will fall within the
#'   range \code{9} to \code{11}.
#'
#' @param jitter_init_sd A proportion (between 0 and 1) to perturb the initial
#'   values for the standard deviation of random effect parameters. The default
#'   is \code{NULL}, which means the same initial values are used across all
#'   chains. A reasonable option might be \code{jitter_init_sd = 0.01}, which
#'   was found to work well during early testing.
#'
#' @param jitter_init_cor A proportion (between 0 and 1) to perturb the initial
#'   values for the correlation parameters of random effects. The default is
#'   \code{NULL}, which means the same initial values are used across all
#'   chains. An option of setting \code{jitter_init_cor = 0.001} was found to be
#'   effective during early testing.
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
#'@param verbose An optional argument (logical, default \code{FALSE}) to
#'  indicate whether to print information collected during setting up the model
#'  formula priors, and initials. As an example, the user might be interested in
#'  knowing the response variables created for the sub model when fitting a
#'  univariate-by-subgroup model. This information can then be used in setting
#'  the desired order of options passed to each such model such as \code{df},
#'  \code{prior}, \code{initials} etc.
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
#'   retrieve the priors (see \code{[brms::get_prior()]} for details).
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
#' @param ... Further arguments passed to [brms::brm()]. This can include
#'   additional arguments that are either passed directly to the underlying
#'   model fitting function or used for internal purposes. Specifically, the
#'   \code{...} can also be used to pass arguments used for testing and
#'   debugging, such as: \code{match_sitar_a_form}, \code{match_sitar_d_form},
#'   \code{sigmamatch_sitar_a_form}, \code{displayit}, \code{setcolh},
#'   \code{setcolb}.
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
#'@importFrom methods formalArgs
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
#' 
#' \donttest{
#' 
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
#' # 2 chains (2000 iterations per chain) with thinning set to 5 for memory  
#' # efficiency. Users are encouraged to refit the model using default settings 
#' # (4 chains, 2000 iterations per chain, thin = 1) as suggested by the Stan team.
#' # Note that with thinning set to 5 (thin = 5), only one fifth of total draws 
#' # will be saved and hence the effective sample size is expected to be small.
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
#'                   chains = 2, cores = 2, iter = 2000, thin = 5)
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
#' plot_ppc(model, ndraws = 100)
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
                   fixed = a + b + c,
                   random = a + b + c,
                   xoffset = mean,
                   bstart = xoffset,
                   cstart = 0,
                   xfun = NULL,
                   yfun = NULL,
                   bound = 0.04,
                   stype = nsp,
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
                   sigmadf = 4,
                   sigmaknots = NA,
                   sigmafixed = a + b + c,
                   sigmarandom = "",
                   sigmaxoffset = mean,
                   sigmabstart = sigmaxoffset,
                   sigmacstart = 0,
                   sigmaxfun = NULL,
                   sigmabound = 0.04,
                   dpar_formula = NULL,
                   autocor_formula = NULL,
                   family = gaussian(),
                   custom_family = NULL,
                   custom_stanvars  = NULL,
                   group_arg = list(
                     groupvar = NULL,
                     by = NULL,
                     cor = un,
                     cov = NULL,
                     dist = gaussian
                   ),
                   sigma_group_arg = list(
                     groupvar = NULL,
                     by = NULL,
                     cor = un,
                     cov = NULL,
                     dist = gaussian
                   ),
                   univariate_by = list(by = NA, cor = un, terms = subset),
                   multivariate = list(mvar = FALSE,
                                       cor = un,
                                       rescor = TRUE),
                   a_prior_beta = normal(lm, ysd, autoscale = TRUE),
                   b_prior_beta = normal(0, 1.5, autoscale = FALSE),
                   c_prior_beta = normal(0, 0.5, autoscale = FALSE),
                   d_prior_beta = normal(0, 1.0, autoscale = TRUE),
                   s_prior_beta = normal(lm, lm, autoscale = TRUE),
                   a_cov_prior_beta = normal(0, 5.0, autoscale = FALSE),
                   b_cov_prior_beta = normal(0, 1.0, autoscale = FALSE),
                   c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                   d_cov_prior_beta = normal(0, 1.0, autoscale = FALSE),
                   s_cov_prior_beta = normal(lm, lm, autoscale = TRUE),
                   a_prior_sd = normal(0, ysd, autoscale = FALSE),
                   b_prior_sd = normal(0, 1.0, autoscale = FALSE),
                   c_prior_sd = normal(0, 0.25, autoscale = FALSE),
                   d_prior_sd = normal(0, 1.0, autoscale = TRUE),
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
                   init_r = NULL,
                   a_init_beta = lm,
                   b_init_beta = 0,
                   c_init_beta = 0,
                   d_init_beta = random,
                   s_init_beta = lm,
                   a_cov_init_beta = random,
                   b_cov_init_beta = random,
                   c_cov_init_beta = random,
                   d_cov_init_beta = random,
                   s_cov_init_beta = random,
                   a_init_sd = random,
                   b_init_sd = random,
                   c_init_sd = random,
                   d_init_sd = random,
                   a_cov_init_sd = random,
                   b_cov_init_sd = random,
                   c_cov_init_sd = random,
                   d_cov_init_sd = random,
                   sigma_init_beta = random,
                   sigma_cov_init_beta = random,
                   sigma_init_sd = random,
                   sigma_cov_init_sd = random,
                   gr_init_cor = random,
                   sigma_init_cor = random,
                   rsd_init_sigma = random,
                   dpar_init_sigma = random,
                   dpar_cov_init_sigma = random,
                   autocor_init_acor = random,
                   autocor_init_unstr_acor = random,
                   mvr_init_rescor = random,
                   r_init_z = random,
                   vcov_init_0 = FALSE,
                   jitter_init_beta = NULL,
                   jitter_init_sd = NULL,
                   jitter_init_cor = NULL,
                   prior_data = NULL,
                   init_data = NULL,
                   init_custom = NULL,
                   verbose = FALSE,
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
                   control = list(adapt_delta = 0.8, max_treedepth = 15),
                   empty = FALSE,
                   rename = TRUE,
                   pathfinder_args = NULL,
                   pathfinder_init = FALSE,
                   sample_prior = "no",
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
                   parameterization = 'ncp',
                   ...) {
  
  # Note
  # Need to work on data argument when using 'sigma_formula_manual' via nlf/lf
  # This is to deal with how 'data' be prepared, used, and stored
  # This is needed to set up separate x for mu and sigma
  # some work done on 20.09.2024
  # The specific areas to look further for are:
  # 'prepare_data' 'data.org.in' 'sigmaxsi' 'setsigma_formula_manual'
  
  
  
  mcall <- match.call()
  
  mcall <- mcall_ <- mcall_dictionary(mcall, envir = NULL, xenvir = NULL)
  
  newcall_checks <- c('threads', 'save_pars')
  
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
  
  # These must be removed to avoid conflict with 'brms' dot arguments
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
  
  # These must be remove to avoid conflict with brms dot arguments
  for (collect_dot_namesi in collect_dot_names) {
    if(!is.null(mcall[[collect_dot_namesi]])) 
      mcall[[collect_dot_namesi]] <- NULL
  }
  
  # Clear alias argument for formula and adjusted
  rm(dots_allias)
  
  mcall <- mcall_ <- mcall
  
  no_default_args <- c("x", "y", "id", "data", "...")
  
  # Problem with rethinking occurs during the expose_model_function
  if("rethinking" %in% (.packages())){
    message("Package 'rethinking' detached and unloaded ato avoid conflict",
            " \nwith the rstan version ", utils::packageVersion('rstan'))
    detach("package:rethinking", unload=TRUE) 
  }
  
  if(utils::packageVersion('rstan') < "2.26") {
    if(expose_function) stop("Argument 'expose_function' not allowed ",
                             "for this rstan version ",
                             utils::packageVersion('rstan'))
  }
  
  
 quote_random_as_init_arg <- function(temp_init_call_in, mcall,...) {
   if(is.null(temp_init_call_in)) temp_init_call_c <- temp_init_call_in
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
 } # quote_random_as_init_arg
  
 
 
  mcall$init <- quote_random_as_init_arg(mcall$init, mcall)
  
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
 
  
  # Initiate non formalArgs()
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
  multivariate_rescor <- NULL;
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
  sigmaboundsi <- NULL;
  sigmabstartsi <- NULL;
  sigmacstartsi <- NULL;
  sigmaids <- NULL;
  sigmaxoffset <- NULL;
  sigmadfs <- NULL;
  ixfuntransformsi <- NULL;
  iyfuntransformsi <- NULL;
  isigmaxfuntransformsi <- NULL;
  nsp <- NULL;
  nsk <- NULL;
  rcs <- NULL;

  enverr. <- environment()
  for (i in names(mcall)[-1]) {
    no_default_args_plus_family <- c(no_default_args, "family")
    if (!i %in% no_default_args_plus_family) {
      assign('err.', FALSE, envir = enverr.)
      tryCatch(
        expr = {
          # suppressWarnings 14 01 2024
          if (is.function(suppressWarnings(eval(mcall[[i]])))) {
            checks. <- deparse_0(mcall[[i]])
          } else {
            # suppressWarnings 14 01 2024
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
  
  
  setdepar0sgub <- c("sigma_formula", 
                     "sigma_formula_gr", 
                     "sigma_formula_manual")
  
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
  if(grepl("^c\\(", deparse_0(familyzzzx), fixed = F)) {
    stop("Argument family should be a list() and not a vector 'c()'")
  } else if(grepl("^\\(", deparse_0(familyzzzx), fixed = F)) {
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
  
  # Override when restricting to rcs
  if(select_model != 'rcs') decomp <- NULL
  
  if(!is.null(decomp)) {
    if(select_model != 'rcs') 
      stop("Decomposition (decomp = 'QR') is allowed only for the RCS model")
  }
  
  
  if(is.character(arguments$select_model)) {
    select_model <- arguments$select_model
  } else if(is.symbol(arguments$select_model)) {
    select_model <- deparse(arguments$select_model)
  } else if(!is.character(arguments$select_model) |
            !is.symbol(arguments$select_model)
            ) {
    stop("The argument 'select_model' must be a symbol or 
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
  
  
  # 24.08.2024
  getdotslist <- list(...)
  
  # spline types supported are 'rcs', 'nsp' and 'nsk'
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
  
  
  allowed_spline_type <- c('rcs', 'nsp', 'nsk')
  allowed_spline_type_exception_msg <- 
    paste("The options available are:", 
          paste(paste(paste0("'", allowed_spline_type, "'"), collapse =", "), 
                collapse =", ")
    )
  

  
  stype_temp_str <- deparse(substitute(stype))
  if(grepl("^list\\(", stype_temp_str)) {
    # stype[[1]] <- deparse(substitute(stype[[1]]))
    # stype[[1]] <- gsub("\"", "", stype[[1]])
    stype <- stype
  } else {
    stype <- stype_temp_str
    stype <- gsub("\"", "", stype)
  }
  

  spline_type_via_stype <- FALSE
  if(!is.null(getdotslist[['smat']])) {
    spline_type <- getdotslist[['smat']]
  } else {
    spline_type <- stype
    spline_type_via_stype <- TRUE
  } 
  
  if(any(spline_type == "NULL")) {
    spline_type_via_stype <- FALSE
  }
  
  if(is.null(getdotslist[['smat']]) & !spline_type_via_stype) {
    spline_type <- 'rcs'
    # if(verbose) message("'rcs' set as default spline type")
  }
    

  # Only expose type and normalize for stype 
  allowed_spline_type_list_names_c <- c('type', 
                                        # 'intercept', 
                                        # 'centerval', 
                                        # 'preH',
                                        # 'include',
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
    spline_type_list[['intercept']]   <- FALSE
    spline_type_list[['normalize']]   <- FALSE
    spline_type_list[['derivs']]      <- FALSE
    spline_type_list[['preH']]        <- FALSE
    spline_type_list[['include']]     <- TRUE
  } else if(!is.null(spline_type)) {
    if(is.list(spline_type)) {
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
              if(verbose) message("stype arguments named as 'type', 'normalize'")
            } else {
              stop(allowed_spline_type_list_names_msg)
            }
          } else {
            stop(allowed_spline_type_list_names_msg)
          }
        }
        if(!is.null(spline_type[['type']])) {
          if(!is.character(spline_type[['type']])) {
            stop(paste0(spline_type[['type']], " must be a character string"))
          } else {
            spline_type_list[['type']] <- spline_type[['type']]
          }
        } else if(is.null(spline_type[['type']])) {
          # 
        }
        
        if(!is.null(spline_type[['intercept']])) {
          if(!is.logical(as.logical(spline_type[['intercept']]))) {
            stop(paste0(spline_type[['intercept']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['intercept']] <- spline_type[['intercept']]
          }
        } else if(is.null(spline_type[['intercept']])) {
          spline_type_list[['intercept']] <- FALSE
        }
        
        if(!is.null(spline_type[['normalize']])) {
          if(!is.logical(as.logical(spline_type[['normalize']]))) {
            stop(paste0(spline_type[['normalize']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['normalize']] <- spline_type[['normalize']]
          }
        } else if(is.null(spline_type[['normalize']])) {
          spline_type_list[['normalize']] <- FALSE
          if(verbose) message(paste0("'", FALSE ,
                                     "' set as default spline normalize"))
        }
        
        if(!is.null(spline_type[['derivs']])) {
          if(!is.logical(as.logical(spline_type[['derivs']]))) {
            stop(paste0(spline_type[['derivs']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['derivs']] <-  spline_type[['derivs']] 
          }
        } else if(is.null(spline_type[['derivs']])) {
          spline_type_list[['derivs']]    <- FALSE
        }
        
        if(!is.null(spline_type[['preH']])) {
          if(!is.logical(as.logical(spline_type[['preH']]))) {
            stop(paste0(spline_type[['preH']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['preH']]    <-  spline_type[['preH']] 
          }
        } else if(is.null(spline_type[['preH']])) {
          spline_type_list[['preH']]    <- FALSE
        }
        
        if(!is.null(spline_type[['include']])) {
          if(!is.logical(as.logical(spline_type[['include']]))) {
            stop(paste0(spline_type[['include']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['include']]    <-  spline_type[['include']] 
          }
        } else if(is.null(spline_type[['include']])) {
          spline_type_list[['include']]    <- TRUE
        }
        
        if(!is.null(spline_type[['path']])) {
          if(!is.character(spline_type[['path']])) {
            stop(paste0(spline_type[['path']], " must be a character string"))
          } else {
            spline_type_list[['path']]    <-  spline_type[['path']] 
          }
        } else if(is.null(spline_type[['path']])) {
          spline_type_list[['path']]    <- NULL
        }
        
        
        if(!is.null(spline_type[['centerval']])) {
          if(!is.numeric(spline_type[['centerval']])) {
            stop(paste0(spline_type[['centerval']], 
                        " must be logical i.e., TRUE/FALSE"))
          } else {
            spline_type_list[['centerval']] <- spline_type[['centerval']]
          }
        } else if(is.null(spline_type[['centerval']])) {
          spline_type_list[['centerval']] <- 0
        }
        
      } else if(length(spline_type) == 0) {
        spline_type_list[['type']]        <- NULL
        spline_type_list[['intercept']]   <- FALSE
        spline_type_list[['centerval']]   <- 0
        spline_type_list[['normalize']]   <- FALSE
        spline_type_list[['derivs']]      <- FALSE
        spline_type_list[['preH']]        <- FALSE
        spline_type_list[['include']]     <- TRUE
        spline_type_list[['path']]        <- NULL
      } # if(length(spline_type) > 0) {
    } else if(!is.list(spline_type)) {
      if(is.character(spline_type)) {
        spline_type_list[['type']]        <- spline_type
        spline_type_list[['intercept']]   <- FALSE
        spline_type_list[['centerval']]   <- 0
        spline_type_list[['normalize']]   <- FALSE
        spline_type_list[['derivs']]      <- FALSE
        spline_type_list[['preH']]        <- FALSE
        spline_type_list[['include']]     <- TRUE
        spline_type_list[['path']]        <- NULL
      } else if(!is.character(spline_type)) {
        stop('augument spline_type must be a character string or a named list')
      } # if(is.character(spline_type)) {
    } # else if(!is.null(spline_type)) {
  } # if(is.null(spline_type)) {
  
  
  
  smat <- spline_type_list[['type']] 
  
   # This to check spline type set using the ... smat
   if(!smat %in% allowed_spline_type)
     stop(paste0("The spline type must be a character string.", 
                 "\n  ",
                 allowed_spline_type_exception_msg)
     )

  if(smat == 'rcs') {
    
  } else if(smat == 'ns') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  } else if(smat == 'nsp') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  } else if(smat == 'nsk') {
    getdotslist[['match_sitar_a_form']] <- match_sitar_a_form <- FALSE
  }
  
  
  
  if((smat == 'nsp' | smat == 'nsk') & spline_type_via_stype) {
    smat_intercept    <- as.integer(spline_type_list[['intercept']])
    smat_centerval    <- as.numeric(spline_type_list[['centerval']])
    smat_normalize    <- as.integer(spline_type_list[['normalize']])
    smat_derivs       <- as.integer(spline_type_list[['derivs']])
    smat_preH         <- as.integer(spline_type_list[['preH']])
    smat_include_stan <- as.integer(spline_type_list[['include']])
    smat_include_path <- spline_type_list[['path']]
    SplinefunxPre  <- 'GS'
    Splinefunxsuf  <- '_call'
    SplinefunxR    <- paste0(SplinefunxPre, "_", smat, Splinefunxsuf)
    SplinefunxStan <- paste0(SplinefunxR, "_", 'stan')
    # when spline type set via 'stype', set normalize = F & include T
  } else if((smat == 'nsp' | smat == 'nsk') & !spline_type_via_stype) {
    smat_intercept <- 0
    smat_centerval <- 0
    smat_normalize <- as.integer(spline_type_list[['normalize']])
    smat_derivs    <- 0
    smat_preH      <- 0
    smat_include_stan <- 0
    smat_include_path <- NULL
    SplinefunxPre  <- 'GS'
    Splinefunxsuf  <- '_call'
    SplinefunxR    <- paste0(SplinefunxPre, "_", smat, Splinefunxsuf)
    SplinefunxStan <- paste0(SplinefunxR, "_", 'stan')
  } else if(smat == 'rcs') { # placeholder to assign these values to envir 
    smat_intercept <- 0
    smat_centerval <- 0
    smat_normalize <- 0
    smat_derivs    <- 0
    smat_preH      <- 0
    smat_include_stan <- 0
    smat_include_path <- NULL
    SplinefunxPre  <- NULL
    Splinefunxsuf  <- NULL
    SplinefunxR    <- NULL
    SplinefunxStan <- NULL
  } else {
    # allow further checks - for later use
  }
   
  # TODO work on smat_preH smat_include_stan to male them compatible
  
  # 'smat_preH' is not allowed because adding two #include does not
  # work in package
  # Hence preH is added to main .stan files 
  # i.e., over riding smat = list(preH = 1)
  if(smat_preH == 1) {
    # stop("Please set preH = 0")
    smat_preH <- 0
   # if(verbose) message("'preH' is set to '0'")
  }
  
  # 'smat_include_stan' is also not working i.e., even single #include also not
  # working in package
  # Hence 'smat_include_stan' is set to 0 and will be pasted to .stan files
  # i.e., over riding smat = list(smat_include_stan = 1)
  if(smat_include_stan == 1) {
    # stop("Please set smat_include_stan = 0")
    smat_include_stan <- 0
   # if(verbose) message("'smat_include_stan' is set to '0'")
  }
 
   
   if(verbose) {
     message(paste0("setting spline type as '", smat, "'"))
     if(smat != "rcs") {
       message(paste0("setting intercept for spline type '",
                      spline_type_list[['type']], "' as: ", smat_intercept))
       message(paste0("setting normalize for spline type '",
                      spline_type_list[['type']], "' as: ", smat_normalize))
       message(paste0("setting centerval for spline type '",
                      spline_type_list[['type']], "' as: ", smat_centerval))
     }
   }
   

  # 24.08.2024
  if(is.null(getdotslist[['match_sitar_a_form']])) {
    match_sitar_a_form <- TRUE
  } else {
    match_sitar_a_form <- getdotslist[['match_sitar_a_form']]
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
      if(getdotslist[['match_sitar_d_form']])
        stop("match_sitar_d_form = TRUE only allowed for sitar model 'sitar4r'")
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
    stop("Currently supported models (via 'select_model' argument) are:",
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
      'future'
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
          message("path for .stan file(s) set to '.' via 'stan_model_args'")
      } else if(!is.null(brms_arguments$stan_model_args)) {
        if(is.list(brms_arguments$stan_model_args)) {
          if(is.null(brms_arguments$stan_model_args[['include_paths']])) {
            brms_arguments$stan_model_args[['include_paths']] <- "."
            if(verbose) 
              message("path for .stan file(s) set to '.' via 'stan_model_args'")
          }
        }
      } # if(is.null(brms_arguments$stan_model_args)) {
    } # if(smat_include_stan) {
  } # if(smat == 'nsp' | smat == 'nsk') {

  
  if(smat == 'nsp' | smat == 'nsk') {
    if(smat_include_stan) {
      if(is.null( brms_arguments$stan_model_args[['include_paths']])) {
        stop("Please specify path for .stan file(s) via 'stan_model_args'")
      }
    }
  }
  
  
  if (eval(brms_arguments$backend) != "rstan" &
      eval(brms_arguments$backend) != "mock" &
      eval(brms_arguments$backend) != "cmdstanr") {
    stop("The backend argument must be either 'rstan', 'mock', or 'cmdstanr'",
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
                           splitmvar, fixed = F)
      } else {
        # splitmvar2 <- gsub(noquote(majors2), majors3, splitmvar2, fixed = F)
        splitmvar2 <- gsub(paste0('\\<', noquote(majors2), '\\>'), majors3, 
                           splitmvar2, fixed = F)
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
      multivariate <-
        gsub("\\s", "", paste(deparse(substitute(multivariate)), 
                              collapse = ""))
      if (multivariate == "T")
        multivariate <- eval(parse(text = multivariate))
      multivariate <- multivariate
      multivariate <- as.list(multivariate)
      names(multivariate) <- 'mvar'
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
          stop("empty list")
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
          stop("empty list")
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
      univariate_by <-
        gsub("\\s", "", paste(deparse(substitute(
          univariate_by
        )), collapse = ""))
      univariate_by <- univariate_by
      univariate_by <- as.list(univariate_by)
      names(univariate_by) <- 'by'
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
          stop("empty list")
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
          stop("empty list")
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
      stop("group_arg should be either NULL or a character",
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
          stop(
            "group_arg should be either NULL or a character",
            " denoting the group idetifier"
          )
        }
        if (temp == "") {
          stop("empty list")
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
          stop("empty list")
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
    stop("group_arg should be either NULL or a character",
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
      stop("sigma_group_arg should be either NULL or a character",
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
          stop(
            "sigma_group_arg should be either NULL or a character",
            " denoting the group idetifier"
          )
        }
        if (temp == "") {
          stop("empty list")
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
          stop("empty list")
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
    stop("sigma_group_arg should be either NULL or a character",
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
      stop(
        "For univeriate-by-subgroup model fitting (via univariate_by argument)",
        "\n ",
        "argument 'by' should be a variable name, '', NULL, or FALSE"
      )
    }
  }
  
  
  
  if (multivariate$mvar &
      !(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    stop(
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
    stop(
      "You have set multivariate as TRUE but provided only one outcome ",
      "\n ",
      " Please set y as list or vector of multiple outcomes such as ",
      "\n ",
      " list(outcome1, outcome2) or y = c(outcome1, outcome2)"
    )
  }
  
  if (!multivariate$mvar & nys > 1) {
    stop(
      "You have set multivariate as FALSE but provided more than one outcome",
      "\n ",
      " Please set y as a symbol / list / vector of single outcome such as",
      "\n ",
      " y = outcome, y = list(outcome1) or y = c(outcome1)"
    )
  }
  
  if (!(is.na(univariate_by$by) |
        univariate_by$by == "NA") & nys > 1) {
    stop(
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
  }
  if (!multivariate$mvar) {
    if (is.null(multivariate$cor))
      multivariate$cor <- "un"
    if (is.null(multivariate$rescor))
      multivariate$rescor <- TRUE
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
      stop(
        paste(
          "\nvariable",
          temp_,
          "used for setting univariate_by-univariate submodels is missing"
        )
      )
    }
    if (!is.factor(data[[temp_]])) {
      stop(temp_, "should be a factor variable")
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
  
  eval_c_list_args <- function(.x, nys, arguments, ...) {
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
    assign(paste0(.xo, "s"), args_s, envir = parent.frame())
  }
  
  
  getArgNames <-
    function(value)
      formalArgs(deparse_0(substitute(value)[[1]]))
  
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
          stop(
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
    "custom_stanvars",
    'pathfinder_args',
    'pathfinder_init',
    "..."
  )
  
  
  
  for (i in convert_to_list) {
    if (!i %in% single_args) {
      to_list_if_not(i, nys, arguments)
    }
  }
  
  for (i in convert_to_list) {
    if (!i %in% single_args) {
      eval_c_list_args(i, nys, arguments)
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
      message(setmsgtxt)
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

  if(!is.null(outliers)) {
    if(is.null(outliers$remove))    outliers$remove <- TRUE
    if(is.null(outliers$icode))     outliers$icode <- c(4,5,6)
    if(is.null(outliers$limit))     outliers$limit <- 5
    if(is.null(outliers$velpower))  outliers$velpower <- 0.5
    if(is.null(outliers$lag))       outliers$lag <- 1
    if(is.null(outliers$linearise)) outliers$linearise <- FALSE
    if(is.null(outliers$verbose))   outliers$verbose <- FALSE
  }
  
  
  
  data.org.in <- data
  uvarby <- NULL
  data <- prepare_data(data = data,
                       x = xs, 
                       y = ys, 
                       id = ids,
                       uvarby = univariate_by$by, 
                       mvar = multivariate$mvar,
                       xfuns = xfuns, 
                       yfuns = yfuns,
                       outliers = outliers,
                       subset = FALSE,
                       envir = enverr.)
  
  ys <- attr(data, "ys")
  subindicators <- attr(data, "subindicators")
  
  if (!is.na(univariate_by$by) & univariate_by$verbose) {
    resvcts_ <- levels(data[[univariate_by$by]])
    resvcts <- paste0(resvcts_, collapse = " ")
    setmsgtxt <- paste0(
      "\n For univariate-by-subgroup model fitting for variable '",
      uvarby,
      "'",
      " (specified via 'univariate_by' argument)",
      "\n ",
      resvcts,
      " response vectors created based on the factor levels",
      "\n\n ",
      "Please check corresponding arguments list.",
      " E.g, df = list(4, 5) denotes that\n df = 4 is for ",
      resvcts_[1],
      ", and  df = 5 is for ",
      resvcts_[2],
      " (and similalry knots, priors, initials etc)",
      "\n\n ",
      "If it does't correspond correctly, then either reverse the list ",
      "arguments\n such as df = list(5, 4),",
      " or else reverse sort the order of factor levels"
    )
    if (displayit == 'msg') {
      message(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolb
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  

  dataout <- priorlist <- NULL
  
  bflist <- list()
  bflist <- initialslist <- initialslist_s <- initsilist <- bflist
  blanketinitslist <- prior_stanvarlist <- auxillary_stanvarlist <- bflist
  
  funlist <- c()
  xoffsetvaluelist <- xoffsetnamelist <- knotsvaluelist <- funlist
  knotsnamelist <- spfun_collect <- xfunvaluelist <- xfunnamelist <- funlist
  yfunvaluelist <- yfunnamelist <- yyfunvaluelist <- yyfunnamelist <- funlist
  xxfunvaluelist <- xxfunnamelist <- fixedvaluelist <- fixednamelist <- funlist
  randomvaluelist <- randomnamelist <- groupvarvaluelist <- funlist
  yvarvaluelist <- ynamelist <- covvaluelist <- covnamelist <- funlist
  groupvarnamelist <- xvarvaluelist <- xnamelist <- funlist
  hierarchicalvarnamelist <- hierarchicalvarvaluelist <- funlist
  
  sigmacovnamelist <- sigmacovvaluelist <- funlist
  
  sigma_groupvarnamelist <- sigma_groupvarvaluelist <- funlist
  sigma_hierarchicalvarnamelist <- sigma_hierarchicalvarvaluelist <- funlist
  
  funlist_r <- funlist_rnamelist <- funlist_rvaluelist <- list()
  
  gq_funs <- list() 
  spfncname_c <- c()
  
  d_adjustedvaluelist <- d_adjustednamelist <- funlist
  
  
  ####
  sigmaspfncname_c <- c()
  sigmaspfun_collect <- funlist
  
  sigmafunlist_r <- sigmafunlist_rnamelist <- sigmafunlist_rvaluelist <- list()
  
  sigmafunlist <- funlist
  sigmafunlist_r <- funlist_r
  sigmagq_funs <- gq_funs
  
  sigmaxfunvaluelist <- sigmaxfunnamelist <- funlist
  sigmaxxfunvaluelist <- sigmaxxfunnamelist <- funlist
  sigmafixedvaluelist <- sigmafixednamelist <- funlist
  sigmarandomvaluelist <- sigmarandomnamelist <- sigmagroupvarvaluelist <- funlist
  
  sigmad_adjustedvaluelist <- sigmad_adjustednamelist <- funlist
  
  xfuntransformvaluelist <- xfuntransformnamelist <- funlist
  ixfuntransformvaluelist <- ixfuntransformnamelist <- funlist
  
  yfuntransformvaluelist <- yfuntransformnamelist <- funlist
  iyfuntransformvaluelist <- iyfuntransformnamelist <- funlist
  
  sigmaxfuntransformvaluelist <- sigmaxfuntransformnamelist <- funlist
  isigmaxfuntransformvaluelist <- isigmaxfuntransformnamelist <- funlist
 
  
  
  # Start loop over response
  for (ii in 1:length(ys)) {
    if (nys > 1)
      resp <- ys[ii]
    else
      resp <- ""
    subindicatorsi <- subindicators[ii]
    
 
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
      stop("either df or knots must be specified")
    }
    if (is.numeric(ept(dfsi)) & is.numeric(ept(knotssi))) {
      # stop("both df and knots specified. Specify one of them\n")
      dfsi <- 'NULL'
      if(verbose) {
        message("The user specified knots are used, hence",
                " the df argument ignored")
      }
    }
    
    
    
    if (!is.numeric(ept(sigmadfsi)) & !is.numeric(ept(sigmaknotssi))) {
      stop("either df or knots must be specified for sigma")
    }
    if (is.numeric(ept(sigmadfsi)) & is.numeric(ept(sigmaknotssi))) {
      # stop("both df and knots specified. Specify one of them\n")
      dfsi <- 'NULL'
      if(verbose) {
        message("The user specified knots are used for sigma, hence",
                " the df argument ignored")
      }
    }
    
    
    

    for (agsxi in letters[1:26]) {
      if(is.null(arguments[[paste0(agsxi, "_", "formula" , "")]])) {
        assign(paste0(agsxi, "_", "formula" , "si") , NULL)
        assign(paste0(agsxi, "_", "formula_gr" , "si") , NULL)
        assign(paste0(agsxi, "_", "formula_gr_str" , "si") , NULL)
        assign(paste0(agsxi, "_", "prior_beta" , "si") , NULL)
        assign(paste0(agsxi, "_", "cov_prior_beta" , "si") , NULL)
        assign(paste0(agsxi, "_", "prior_sd" , "si") , NULL)
        assign(paste0(agsxi, "_", "cov_prior_sd" , "si") , NULL)
        assign(paste0(agsxi, "_", "init_beta", "si" ) , NULL)
        assign(paste0(agsxi, "_", "cov_init_beta" , "si") , NULL)
        assign(paste0(agsxi, "_", "init_sd", "si" ) , NULL)
        assign(paste0(agsxi, "_", "cov_init_sd" , "si") , NULL)
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
          stop("For model '", select_model, "'", ", 
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
          stop("For model '", select_model, "'", ", 
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
        stop(msg_1, "\n ", msg_2, " \n ", see_what_formual)
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
        stop(msg_1, "\n ", msg_2, " \n ", see_what_formual)
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
          if(!match_sitar_d_form) stop(msg_mismatch_fixed_random_str)
        } else {
          stop(msg_mismatch_fixed_random_str)
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
        match_sitar_d_form <- FALSE
      } else if(grepl("d", fixedsi) & !grepl("d", randomsi)) {
        sitar_nparms <- 4
        match_sitar_d_form <- FALSE
      } else if(!grepl("d", fixedsi) & grepl("d", randomsi)) {
        sitar_nparms <- 4
        match_sitar_d_form <- TRUE
      } else if(!grepl("d", fixedsi) & !grepl("d", randomsi)) {
        sitar_nparms <- 3
        match_sitar_d_form <- FALSE
      }
    }
    
    
    
    # Model specific number of fixed and random parameters
    allowed_parm_letters <- NULL
    if(select_model == 'sitar') allowed_parm_letters <- letters[1:sitar_nparms]
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
    
    
    # Covariate not allowed when matching to sitar 'd' form
    if(select_model == 'sitar') {
      if (!match_sitar_d_form) {
        if (!grepl("d", fixedsi, fixed = T) &
            grepl("d", randomsi, fixed = T)) {
          stop(
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
      
      
      if (match_sitar_d_form) {
        if ((grepl("d", fixedsi, fixed = T) |
             grepl("d", randomsi, fixed = T)) &
            (!grepl("^~1$", d_formulasi) |
             !grepl("^~1$", d_formula_grsi))) {
          stop(
            "Parameter 'd' is missing in the fixed effects part of the model ",
            "\n ",
            " but specified in the random effects part of the model ",
            "\n ",
            " (This is to match with the 'sitar' package's formulation)",
            "\n ",
            " For this formulation (i.e., 'd' is missing in the fixed effects)",
            "\n ",
            " covariate(s) are not allowed"
          )
        }
      }
    } # if(select_model == 'sitar') {
    
    

    if(select_model == 'sitar') {
      if(!any(grepl('s', fixedsi))) fixedsi <- paste0(fixedsi, "+", "s")
    }
    
    
    if(select_model == 'rcs') {
      if(!any(grepl('s', fixedsi))) fixedsi <- paste0(fixedsi, "+", "s")
      if(!any(grepl('s', randomsi)) & rcs_add_re_spline) {
          randomsi <- paste0(randomsi, "+", "s")
      }
      if(any(grepl('s', randomsi)) & !rcs_add_re_spline) {
        stop("you have specified select_model = 'rcsf' (i.e., no random",
             "\n ",
             "spline effects) but your random argument have parameter 's'.",
             "\n ",
             "Please check your 'select_model' and 'random' arguments")
      }
    }
    
    
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
            stop("Models beyound two levels of hierarchy are not supported yet",
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
              message("Argument '", 
                      substitute(x_grsi), "' changed from '", 
                      x_grsi , "' to  '~1'.")
              message("Instead of '", 
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
      out
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
        stop("Random effect for parameter 'sigma' are not allowed",
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
        # strpartstrx_form2 <- strsplit(strpartstrx_form, "(", fixed = T)[[1]] [1]
        # strpartstrx_form <- paste0(, collapse = "")
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
      
      x_formula_grsi
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
    
   
    
    check_formuals <-
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
        "sigma_formulasi",
        "sigma_formula_grsi"
      )
    
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
        stop('autocor_formula argument should be a formula. E.g.,',
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
        familysi_check <- paste0('brms::brmsfamily', familysi_check)
      } else if(!grepl('brmsfamily', familysi_check) & 
                !grepl('family', familysi_check)) {
        stop("Argument family should be specified as brmsfamily(family,...)")
      }
      familysi <- familysi_check
    }
    
    
    if (!is.null(familysi)) {
      familysi <- list_to_quoted_if_not_si(familysi)
    }
    
    
    
    
    
    if (!is.null(dpar_formulasi)) {
      if (grepl("^lf\\(", dpar_formulasi) |
          grepl("^nlf\\(", dpar_formulasi)) {
      } else {
        dpar_formulasi <- dpar_formulasi
      }
    }
    

    N_J_all <- length(unique(data[[idsi]]))
    
    
    ##########################
    setsigma_formula_manual <- FALSE
    if (is.null(sigma_formula_manualsi[[1]][1])) {
      setsigma_formula_manual <- FALSE
    } else if(sigma_formula_manualsi == "NULL") {
      setsigma_formula_manual <- FALSE
    } else {
      setsigma_formula_manual <- TRUE
    }
    
    
    # 20.09.2024
    if (!is.null(sigmaxsi[[1]][1]) & sigmaxsi != "NULL") {
      if(identical(sigmaxsi, xsi)) {
        if(verbose) {
          message("Since both ", sigmaxsi, " and ", xsi, " are identical, ",
                  "\n ", 
                  "the ", sigmaxsi, " has been renamed as ", 
                  paste0("sigma", xsi)
          )
        }
        sigmaxsi <- paste0("sigma", xsi)
      } else {
        sigmaxsi <- sigmaxsi
      }
      data[[sigmaxsi]] <- data[[xsi]]
    } 
    
    # Even if not modelling location scale model, create data[[sigmaxsi]]
    # Need to thing how to shut it off compleletely
    if (is.null(sigmaxsi[[1]][1]) | sigmaxsi == "NULL") {
      sigmaxsi <- paste0("sigma", xsi) 
      # if(verbose) {
      #   message("The predictor for distrubutional parameter (i.e., sigma) is set same as mu i.e, ", xsi, ".",
      #           "\n ", 
      #           "However, it has been renamed as ", 
      #           paste0("sigma", xsi)
      #   )
      # }
      data[[sigmaxsi]] <- data[[xsi]]
    }

    
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
      sortbylayer <- NA
      data <- data %>%
        dplyr::mutate(sortbylayer =
                        forcats::fct_relevel(!!as.name(univariate_by$by),
                                             (levels(
                                               !!as.name(univariate_by$by)
                                             )))) %>%
        dplyr::arrange(sortbylayer) %>%
        dplyr::mutate(!!as.name(idsi) := factor(!!as.name(idsi),
                                                levels = 
                                                  unique(!!as.name(idsi)))) %>% 
        dplyr::select(-sortbylayer)
      
      
      
      # if (!is.null(sigmaxsi[[1]][1]) & sigmaxsi != "NULL") {
      #   if(identical(sigmaxsi, xsi)) {
      #     if(verbose) {
      #       message("Since both ", sigmaxsi, " and ", xsi, " are identical, ",
      #               "\n ", 
      #               "the ", sigmaxsi, " has been renamed as ", 
      #               paste0("sigma", xsi)
      #               )
      #     }
      #     sigmaxsi <- paste0("sigma", sigmaxsi)
      #   }
      #   sigmaxsi <- sigmaxsi
      #   data[[sigmaxsi]] <- data[[sigmaxsi]]
      # } else {
      #   sigmaxsi <- paste0("sigma", xsi) 
      #   data[[sigmaxsi]] <- data[[xsi]]
      #   if(verbose) message("predictor for sigma is set same as for mu")
      # }
      # 
     
      
      datai <- data %>%
        dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>%
        droplevels()
      if (!subindicatorsi %in% colnames(datai)) {
        stop("variable ", subindicatorsi, " not in the dataframe")
      }
      if (!xsi %in% colnames(datai))
        stop("variable ", xsi, " not in the dataframe")
      if (!idsi %in% colnames(datai))
        stop("variable ", idsi, " not in the dataframe")
    }
    
    if ((is.na(univariate_by$by) | univariate_by$by == "NA")) {
      datai <- data %>%
        droplevels()
      if (!ysi %in% colnames(datai))
        stop("variable ", ysi, " not in the dataframe")
      if (!xsi %in% colnames(datai))
        stop("variable ", xsi, " not in the dataframe")
      if (!idsi %in% colnames(datai))
        stop("variable ", idsi, " not in the dataframe")
    }
    
    
    # 28 01 2024
    datai <- datai %>% tidyr::drop_na()
    
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
    
    
    # Refactor to use function() for transformations
    # This will allow using optimize_x = list(function(x) log(x + 3/4))
    # Note that instead of calling log(data[[xsi]]), 'xfuntransformsi' will be used 
   
    # Check if xfunsi set
    set_xfunsi <- check_if_arg_set(xfunsi)
    
    # Note that y transformation is done within the prepare_data function,
    # see data <- prepare_data(
    
    # Check if yfunsi set
    set_yfunsi <- check_if_arg_set(yfunsi)
    
    # Check if sigmaxfunsi called
    set_sigmaxfunsi <- check_if_arg_set(sigmaxfunsi)
   

    if (!set_xfunsi) {
      xfuntransformsi <- function(x)x
      assign('xfuntransformsi', xfuntransformsi, envir = enverr.)
    }
    
    if (set_xfunsi) {
      if(xfunsi == "log") {
        xfuntransformsi <- function(x)log(x)
      } else if(xfunsi == "sqrt") {
        xfuntransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(xfunsi))) {
        xfuntransformsi <- ept(xfunsi)
      } else {
        stop(paste0("The xfun argument must be either a string ('log' or 'sqrt'),", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('xfuntransformsi', xfuntransformsi, envir = enverr.)
    }
    
    
    if (!set_yfunsi) {
      yfuntransformsi <- function(x)x
      assign('yfuntransformsi', yfuntransformsi, envir = enverr.)
    }
    
    if (set_yfunsi) {
      if(yfunsi == "log") {
        yfuntransformsi <- function(x)log(x)
      } else if(yfunsi == "sqrt") {
        yfuntransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(yfunsi))) {
        yfuntransformsi <- ept(yfunsi)
      } else {
        stop(paste0("The xfun argument must be either a string ('log' or 'sqrt'),", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('yfuntransformsi', yfuntransformsi, envir = enverr.)
    }
    
    
    if (!is.null(sigmaxoffset[[1]][1]) & sigmaxoffset != "NULL") {
      sigmaxoffset <- sigmaxoffset
    } else {
      sigmaxoffset <- xoffset
      if(verbose) message("xoffset for sigma is set same as for mu")
    }
    
    
    if (!is.null(sigmaxoffsetsi[[1]][1]) & sigmaxoffsetsi != "NULL") {
      sigmaxoffsetsi <- sigmaxoffsetsi
    } else {
      sigmaxoffsetsi <- xoffsetsi
      if(verbose) message("xoffset for sigma is set same as for mu")
    }
    
    
    
    check_for_nan_inf <- function(x) {
      suppressWarnings({
        nan_inf <- FALSE
        if(is.infinite(x)) {
          nan_inf <- TRUE
        } else if(is.na(x)) {
          nan_inf <- TRUE
        } else if(is.nan(x)) {
          nan_inf <- TRUE
        }
      })
      nan_inf
    }
    
    
    datai[[xsi]] <- xfuntransformsi(datai[[xsi]])
    
    if(check_is_numeric_like(xoffsetsi)) {
      zm <- as.numeric(xoffsetsi)
      if(check_for_nan_inf(xfuntransformsi(zm))) {
        if(verbose) message("'xoffset' value '", zm, 
                            "' can not be transformed to", 
                           " match the xfun based transformation of x " ,"")
      } else {
        xoffsetsi <- round(xfuntransformsi(zm), 2)
        if(verbose) message("'xoffset' value '", zm, "' transformed to '", 
                            xoffsetsi ,"'")
      }
      rm('zm')
    }
   
    
    if (set_xfunsi) {
      xfunvalue <- xfunsi
    } else {
      xfunvalue <- NULL
    }
    
    if (set_yfunsi) {
      yfunvalue <- yfunsi
    } else {
      yfunvalue <- NULL
    }
    
  
    if (!set_sigmaxfunsi) {
      sigmaxfuntransformsi <- function(x)x
      assign('sigmaxfuntransformsi', sigmaxfuntransformsi, envir = enverr.)
    }
    
    
    if (set_sigmaxfunsi) {
      if(sigmaxfunsi == "log") {
        sigmaxfuntransformsi <- function(x)log(x)
      } else if(sigmaxfunsi == "sqrt") {
        sigmaxfuntransformsi <- function(x)sqrt(x)
      } else  if(is.function(ept(sigmaxfunsi))) {
        sigmaxfuntransformsi <- ept(sigmaxfunsi)
      } else {
        stop(paste0("The xfun argument must be either 'log' or 'sqrt',", 
                    "\n  ",
                    "or a function such as function(x)log(x)"))
      }
      assign('sigmaxfuntransformsi', sigmaxfuntransformsi, envir = enverr.)
    }
    
    
    
    if (!is.null(sigmaxfunsi[[1]][1]) & sigmaxfunsi != "NULL") {
      datai[[sigmaxsi]] <- sigmaxfuntransformsi(datai[[sigmaxsi]])
      if(check_is_numeric_like(sigmaxoffsetsi)) {
        zm <- as.numeric(sigmaxoffsetsi)
        if(check_for_nan_inf(sigmaxfuntransformsi(zm))) {
          if(verbose) message("'sigmaxoffset' value '", zm, 
                              "' can not be transformed to", 
                              " match the xfun based transformation of x" ,"")
        } else {
          sigmaxoffsetsi <- round(sigmaxfuntransformsi(zm), 2)
          if(verbose) message("'xoffset' value '", zm, "' transformed to '", 
                              sigmaxoffsetsi ,"'")
        }
        rm('zm')
      }
    }
    
    
    # Assign reverse functions also
    assign("ixfuntransformsi",  sitar::ifun(base::body(xfuntransformsi)), 
           envir = enverr.)
    
    
    assign("iyfuntransformsi",  sitar::ifun(base::body(yfuntransformsi)), 
           envir = enverr.)
    
    
    assign("isigmaxfuntransformsi",  sitar::ifun(base::body(sigmaxfuntransformsi)), 
           envir = enverr.)
    
    
    
    xfuntransformvalue  <- xfuntransformsi
    ixfuntransformvalue <- ixfuntransformsi
    
    yfuntransformvalue  <- yfuntransformsi
    iyfuntransformvalue <- iyfuntransformsi
    
    
    sigmaxfuntransformvalue  <- sigmaxfuntransformsi
    isigmaxfuntransformvalue <- isigmaxfuntransformsi
    
    
    
    if (set_sigmaxfunsi) {
      sigmaxfunvalue <- sigmaxfunsi
    } else {
      sigmaxfunvalue <- NULL
    }
    
    
    
    if (nys == 1) {
      xfun_name <- "xfun"
      yfun_name <- "yfun"
      xxfun_name <- "xvar_xfun"
      yyfun_name <- "yvar_yfun"
      sigmaxfun_name <- "sigmaxfun"
      sigmayfun_name <- "sigmayfun"
      xfuntransform_name <- "xfuntransform"
      ixfuntransform_name <- "ixfuntransform"
      yfuntransform_name <- "yfuntransform"
      iyfuntransform_name <- "iyfuntransform"
      sigmaxfuntransform_name <- "sigmaxfuntransform"
      isigmaxfuntransform_name <- "isigmaxfuntransform"
      
    } else if (nys > 1) {
      xfun_name <- paste0("xfun", "_", ysi)
      yfun_name <- paste0("yfun", "_", ysi)
      xxfun_name <- paste0("xvar_xfun", "_", ysi)
      yyfun_name <- paste0("yvar_yfun", "_", ysi)
      sigmaxfun_name <- paste0("sigmaxfun", "_", ysi)
      sigmayfun_name <- paste0("sigmayfun", "_", ysi)
      xfuntransform_name <- paste0("xfuntransform", "_", ysi)
      ixfuntransform_name <- paste0("ixfuntransform", "_", ysi)
      yfuntransform_name <- paste0("yfuntransform", "_", ysi)
      iyfuntransform_name <- paste0("iyfuntransform", "_", ysi)
      sigmaxfuntransform_name <- paste0("sigmaxfuntransform", "_", ysi)
      isigmaxfuntransform_name <- paste0("isigmaxfuntransform", "_", ysi)
    }
    
    
    xfunvaluelist[[ii]] <- xfunvalue
    xfunnamelist[[ii]] <- xfun_name
    
    yfunvaluelist[[ii]] <- yfunvalue
    yfunnamelist[[ii]] <- yfun_name
    
    sigmaxfunvaluelist[[ii]] <- sigmaxfunvalue
    sigmaxfunnamelist[[ii]] <- sigmaxfun_name
    
    xfuntransformvaluelist[[ii]] <- xfuntransformvalue
    xfuntransformnamelist[[ii]]  <- xfuntransform_name
    
    ixfuntransformvaluelist[[ii]] <- ixfuntransformvalue
    ixfuntransformnamelist[[ii]]  <- ixfuntransform_name
    
    yfuntransformvaluelist[[ii]] <- yfuntransformvalue
    yfuntransformnamelist[[ii]]  <- yfuntransform_name
    
    iyfuntransformvaluelist[[ii]] <- iyfuntransformvalue
    iyfuntransformnamelist[[ii]]  <- iyfuntransform_name
    
    sigmaxfuntransformvaluelist[[ii]] <- sigmaxfuntransformvalue
    sigmaxfuntransformnamelist[[ii]]  <- sigmaxfuntransform_name
    
    isigmaxfuntransformvaluelist[[ii]] <- isigmaxfuntransformvalue
    isigmaxfuntransformnamelist[[ii]]  <- isigmaxfuntransform_name
    
    
    if (set_xfunsi) {
      xxfunvaluelist[[ii]] <- paste0(xfunsi, "(", xsi, ")")
    } else {
      xxfunvaluelist[[ii]] <- NULL
    }
    
    if (set_sigmaxfunsi) {
      sigmaxxfunvaluelist[[ii]] <- paste0(sigmaxfunsi, "(", sigmaxsi, ")")
    } else {
      sigmaxxfunvaluelist[[ii]] <- NULL
    }
    
    
    if (set_yfunsi) {
      yyfunvaluelist[[ii]] <- paste0(yfunsi, "(", ysi, ")")
    } else {
      yyfunvaluelist[[ii]] <- NULL
    }
    
    xxfunnamelist[[ii]] <- xxfun_name
    yyfunnamelist[[ii]] <- yyfun_name # xxfun_name
    
    sigmaxxfunnamelist[[ii]] <- sigmaxfun_name
    
    gkn <- function(x, df, bounds) {
      c(min(x) - bounds * (max(x) - min(x)),
        quantile(x, (1:(df - 1)) / df, na.rm = TRUE), # 28 01 2024
        max(x) +
          bounds * (max(x) - min(x)))
    }
    
    if (is.numeric(ept(knotssi))) {
      knots <- ept(knotssi)
    }
    if (is.numeric(ept(dfsi))) {
      knots <- (unname(gkn(datai[[xsi]], ept(dfsi), ept(boundsi))))
    }
    
   
    
    if (is.numeric(ept(sigmaknotssi))) {
      sigmaknots <- ept(sigmaknotssi)
    }
    if (is.numeric(ept(sigmadfsi))) {
      sigmaknots <- (unname(gkn(datai[[sigmaxsi]], 
                                ept(sigmadfsi), ept(sigmaboundsi))))
    }
    
    
    if(select_model == "sitar") {
      if (match_sitar_d_form) {
        if (length(knots) > 2) {
          itemp <- strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
          itemp <- itemp[!grepl("d", itemp)]
          fixedsi <- paste(itemp, collapse = "+")
        }
      }
    }
    
  
    
    nabci <- length(strsplit(gsub("\\+", " ", fixedsi), " ")[[1]])
    nabcrei <-
      length(strsplit(gsub("\\+", " ", randomsi), " ")[[1]])
   
    
    eval_xoffset_bstart_args <-
      function(x, y, knots, data, eval_arg, xfunsi, arg = 'xoffset') {
        if (eval_arg == "mean") {
          eval_arg.o <- mean(data[[x]])
        } else if (eval_arg == "min") {
          eval_arg.o <- min(data[[x]])
        } else if (eval_arg == "max") {
          eval_arg.o <- max(data[[x]])
        } else if (eval_arg == "apv") {
          # mat_s <- make_spline_matrix(data[[x]], knots)
          if(smat == 'rcs') {
            mat_s <- make_spline_matrix(data[[x]], knots)
          } else if(smat == 'nsp') {
            iknots <- knots[2:(length(knots)-1)]
            bknots <- c(knots[1], knots[length(knots)])
            mat_s <- GS_nsp_call(x = data[[x]], knots = iknots, bknots = bknots, 
                                intercept = smat_intercept, derivs = smat_derivs, 
                                centerval = smat_centerval, 
                                normalize = smat_normalize,
                                preH = smat_preH)
          } else if(smat == 'nsk') {
            iknots <- knots[2:(length(knots)-1)]
            bknots <- c(knots[1], knots[length(knots)])
            mat_s <- GS_nsk_call(x = data[[x]], knots = iknots, bknots = bknots, 
                                intercept = smat_intercept, derivs = smat_derivs, 
                                centerval = smat_centerval, 
                                normalize = smat_normalize,
                                preH = smat_preH)
          }
          lmform <- as.formula(paste0(y, "~1+", "mat_s"))
          lmfit <- lm(lmform, data = data)
          eval_arg.o <- sitar::getPeak(data[[x]],
                                       predict(smooth.spline(data[[x]],
                                                             fitted(lmfit)),
                                               data[[x]], deriv = 1)$y)[1]
          if(is.na(eval_arg.o)) {
            stop(arg, " specified as '", eval_arg, "' returned NA.",
                 "\n ",
                 " Please change ", arg,
                 " argument to 'mean' or a numeric value.")
          }
        } else {
          eval_arg.o <- ept(eval_arg)
        }
        out <- as.numeric(eval_arg.o)
        out <- round(out, 3)
        return(out)
      }
    
    
    eval_xoffset_cstart_args <-
      function(x, y, knots, data, eval_arg, xfunsi) {
        if (eval_arg == "pv") {
          # mat_s <- make_spline_matrix(data[[x]], knots)
          if(smat == 'rcs') {
            mat_s <- make_spline_matrix(data[[x]], knots)
          } else if(smat == 'nsp') {
            iknots <- knots[2:(length(knots)-1)]
            bknots <- c(knots[1], knots[length(knots)])
            mat_s <- GS_nsp_call(x = data[[x]], 
                                 knots = iknots, bknots = bknots, 
                                intercept = smat_intercept, 
                                derivs = smat_derivs, 
                                centerval = smat_centerval, 
                                normalize = smat_normalize,
                                preH = smat_preH)
          } else if(smat == 'nsk') {
            iknots <- knots[2:(length(knots)-1)]
            bknots <- c(knots[1], knots[length(knots)])
            mat_s <- GS_nsk_call(x = data[[x]], 
                                 knots = iknots, bknots = bknots, 
                                 intercept = smat_intercept, 
                                 derivs = smat_derivs, 
                                 centerval = smat_centerval, 
                                 normalize = smat_normalize,
                                 preH = smat_preH)
          }
          lmform <- as.formula(paste0(y, "~1+", "mat_s"))
          lmfit <- lm(lmform, data = data)
          eval_arg.o <- sitar::getPeak(data[[x]],
                                       predict(smooth.spline(data[[x]],
                                                             fitted(lmfit)),
                                               data[[x]], deriv = 1)$y)[2]
          if(is.na(eval_arg.o)) {
            stop("cstart specified as '", eval_arg, "' returned NA.",
                 "\n ",
                 " Please change cstart argument to 'mean' or a numeric value.")
          }
        } else {
          eval_arg.o <- ept(eval_arg)
        }
        out <- as.numeric(eval_arg.o)
        out <- round(out, 3)
        return(out)
      }
    
    
    
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
    
    
    xoffset <-
      eval_xoffset_bstart_args(xsi, ysi, knots, datai, 
                               xoffsetsi, xfunsi, arg = 'xoffset')
    
    bstart <-
      eval_xoffset_bstart_args(xsi, ysi, knots, datai, 
                               bstartsi, xfunsi, arg = 'bstart')
    
    bstart <- bstart - xoffset
    
    cstart <-
      eval_xoffset_cstart_args(xsi, ysi, knots, datai, cstartsi, xfunsi)
    
    
    xoffset <- round(xoffset, 8)
    datai[[xsi]] <- datai[[xsi]] - xoffset
    knots <- knots - xoffset
    knots <- round(knots, 8)
    nknots <- length(knots)
    df <- length(knots) - 1
   
    # mat_s <- make_spline_matrix(datai[[xsi]], knots)
    
    if(smat == 'rcs') {
      mat_s <- make_spline_matrix(datai[[xsi]], knots)
    } else if(smat == 'nsp') {
      iknots <- knots[2:(length(knots)-1)]
      bknots <- c(knots[1], knots[length(knots)])
      mat_s <- GS_nsp_call(x = datai[[xsi]], 
                           knots = iknots, bknots = bknots, 
                          intercept = smat_intercept, 
                          derivs = smat_derivs, 
                          centerval = smat_centerval, 
                          normalize = smat_normalize,
                          preH = smat_preH)
    } else if(smat == 'nsk') {
      iknots <- knots[2:(length(knots)-1)]
      bknots <- c(knots[1], knots[length(knots)])
      mat_s <- GS_nsk_call(x = datai[[xsi]], 
                           knots = iknots, bknots = bknots, 
                           intercept = smat_intercept, 
                           derivs = smat_derivs, 
                           centerval = smat_centerval, 
                           normalize = smat_normalize,
                           preH = smat_preH)
    }
    
    

    
    if(sigmabstartsi == 'sigmaxoffset') {
      sigmabstartsi <- sigmaxoffsetsi
    }
    
    if(sigmaxoffsetsi == "apv") stop("xoffset can not be apv for sigma")
    
    sigmaxoffset <-
      eval_xoffset_bstart_args(sigmaxsi, ysi, sigmaknots, datai, 
                               sigmaxoffsetsi, sigmaxfunsi, arg = 'offset')
    
    sigmabstart <-
      eval_xoffset_bstart_args(sigmaxsi, ysi, sigmaknots, datai, 
                               sigmabstartsi, sigmaxfunsi, arg = 'bstart')
    
    sigmabstart <- sigmabstart - sigmaxoffset
    
    sigmacstart <-eval_xoffset_cstart_args(sigmaxsi, ysi, sigmaknots, datai, 
                               sigmacstartsi, sigmaxfunsi)
    
    sigmaxoffset <- round(sigmaxoffset, 8)
    datai[[sigmaxsi]] <- datai[[sigmaxsi]] - sigmaxoffset
    sigmaknots <- sigmaknots - sigmaxoffset
    sigmaknots <- round(sigmaknots, 8)
    sigmanknots <- length(sigmaknots)
    sigmadf <- length(sigmaknots) - 1
    
    
    # Define names for Stan functions 
    SplineFun_name  <- paste0(toupper(select_model), "", 'Fun') # "DefFun" 
    getX_name       <- "getX"
    getKnots_name   <- "getKnots"
    
    getpreH_name   <- "getpreH"
    
    if (nys > 1) {
      spfncname <- paste0(ysi, "_", SplineFun_name)
      getxname <- paste0(ysi, "_", getX_name)
      getknotsname <- paste0(ysi, "_", getKnots_name)
      getpreHname <- paste0(ysi, "_", getKnots_name)
    } else if (nys == 1) {
      spfncname <- SplineFun_name
      getxname <- getX_name
      getknotsname <- getKnots_name
      getpreHname <- getpreH_name
    }
    
    spfncname_c <- c(spfncname_c, spfncname)
    
    spfun_collect <-
      c(spfun_collect, c(spfncname, paste0(spfncname, "_", c("d1",
                                                             "d2"))))
    
    decomp_editcode <- FALSE
    if(select_model == 'rcs') {
      decomp_editcode <- FALSE
    }
    
    # For QR decomp to pass to prepare_function 
    # This to check s covs - re
    checkscovsi <-  getcovlist(s_formulasi)
    
    # This can be set to TRUE for RCS
    if(!is.null(checkscovsi)) {
      add_b_Qr_genquan_s_coef <- FALSE
    } else {
      add_b_Qr_genquan_s_coef <- FALSE
    }
    
    # This control whether to add scode for genquant block for QR model
    # Relevant in both and prepare_function
    add_rcsfunmatqrinv_genquant <- FALSE # TRUE
    

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
        "d_adjustedsi",
        'xfunsi',
        'yfunsi',
        'xoffset',
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
        "smat_intercept",
        "smat_derivs",
        "smat_centerval",
        "smat_normalize",
        "smat_preH",
        "smat_include_stan",
        "smat_include_path",
        "SplinefunxPre",
        "Splinefunxsuf",
        "SplinefunxR",
        "SplinefunxStan"
      )
    
    internal_function_args <- list()
    internal_function_args <- mget(internal_function_args_names)
    
    if (verbose) {
      if (ii == 1) {
        setmsgtxt <- paste0("\n Preparing function")
        if (displayit == 'msg') {
          message(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
    }
    
  
    
    if(smat == 'rcs') {
      get_s_r_funs <- 
        prepare_function(
          x = xsi,
          y = ysi,
          id = idsi,
          knots = knots,
          nknots = nknots,
          data = datai,
          internal_function_args = internal_function_args
        )
    } else if(smat == 'nsp') {
      get_s_r_funs <- 
        prepare_function_nsp(
          x = xsi,
          y = ysi,
          id = idsi,
          knots = knots,
          nknots = nknots,
          data = datai,
          internal_function_args = internal_function_args
        )
    } else if(smat == 'nsk') {
      get_s_r_funs <- 
        prepare_function_nsp(
          x = xsi,
          y = ysi,
          id = idsi,
          knots = knots,
          nknots = nknots,
          data = datai,
          internal_function_args = internal_function_args
        )
    }
    
    
    
    
    funlist[ii] <- get_s_r_funs[['rcsfun']]
    funlist_r[[ii]] <- get_s_r_funs[['r_funs']]
    gq_funs[[ii]] <- get_s_r_funs[['gq_funs']]
    
    include_fun_names <- get_s_r_funs[['include_fun_names']]
    
    
    #################################################
    #################################################
    #################################################
    #################################################
    # moved up at 42..
    # setsigma_formula_manual <- FALSE
    # if (is.null(sigma_formula_manualsi[[1]][1])) {
    #   setsigma_formula_manual <- FALSE
    # } else if(sigma_formula_manualsi == "NULL") {
    #   setsigma_formula_manual <- FALSE
    # } else {
    #   setsigma_formula_manual <- TRUE
    # }
    
    
    # Define sigma function
    if(setsigma_formula_manual) {
      # Define names for Stan functions 
      sigmaSplineFun_name  <- paste0("sigma", SplineFun_name)
      sigmagetX_name       <- paste0("sigma", getX_name)
      sigmagetKnots_name   <- paste0("sigma", getKnots_name)
      sigmagetpreH_name    <- paste0("sigma", getpreH_name)
      
      # sigmaSplineFun_name  <- paste0("sigma", toupper(select_model), "", 'Fun')
      # sigmagetX_name       <- paste0("sigma", "getX")
      # sigmagetKnots_name   <- paste0("sigma", "getKnots")
      
      
      if (nys > 1) {
        sigmaspfncname <- paste0(ysi, "_", sigmaSplineFun_name)
        sigmagetxname <- paste0(ysi, "_", sigmagetX_name)
        sigmagetknotsname <- paste0(ysi, "_", sigmagetKnots_name)
        sigmagetpreHname <- paste0(ysi, "_", sigmagetpreH_name)
        
      } else if (nys == 1) {
        sigmaspfncname <- sigmaSplineFun_name
        sigmagetxname <- sigmagetX_name
        sigmagetknotsname <- sigmagetKnots_name
        sigmagetpreHname <- sigmagetpreH_name
      }
      
      sigmaspfncname_c <- c(sigmaspfncname_c, sigmaspfncname)
      
      sigmaspfun_collect <-
        c(sigmaspfun_collect, c(sigmaspfncname, 
                                paste0(sigmaspfncname, "_", 
                                       c("d1",
                                         "d2"))))
      
      sigmadecomp_editcode <- FALSE
      if(select_model == 'rcs') {
        sigmadecomp_editcode <- FALSE
      }
      
      
      # For QR decomp to pass to prepare_function 
      # This to check s covs - re
      sigmacheckscovsi <-  getcovlist(s_formulasi)
      
      # This can be set to TRUE for RCS
      if(!is.null(sigmacheckscovsi)) {
        sigmaadd_b_Qr_genquan_s_coef <- FALSE
        # if(select_model == 'rcs') add_b_Qr_genquan_s_coef <- TRUE
      } else {
        sigmaadd_b_Qr_genquan_s_coef <- FALSE
      }
      
      # This control whether to add scode for genquant block for QR model
      # Relevant in both and prepare_function
      sigmaadd_rcsfunmatqrinv_genquant <- FALSE # TRUE
      
      
      # copy internal_function_args but later replace them by sigma args
      sigmainternal_function_args <- internal_function_args
      
      # These are defined for sigma via bsitar()
      # These are copied from the mu part
      sigmamatch_sitar_d_form          <- match_sitar_d_form
      sigmad_adjustedsi                <- d_adjustedsi
      sigmayfunsi                      <- yfunsi
      sigmabrms_arguments              <- brms_arguments
      sigmaselect_model                <- select_model
      sigmadecomp                      <- decomp
      sigmadecomp_editcode             <- decomp_editcode
      sigmanys                         <- nys
      sigmacheckscovsi                 <- checkscovsi
      sigmaadd_b_Qr_genquan_s_coef     <- add_b_Qr_genquan_s_coef
      sigmaadd_rcsfunmatqrinv_genquant <- add_rcsfunmatqrinv_genquant
      
      
      sigmainternal_function_args[['fixedsi']] <- 
        sigmafixedsi
      sigmainternal_function_args[['randomsi']] <- 
        sigmarandomsi
      sigmainternal_function_args[['spfncname']] <- 
        sigmaspfncname
      sigmainternal_function_args[['getxname']] <- 
        sigmagetxname
      sigmainternal_function_args[['getknotsname']] <- 
        sigmagetknotsname
      sigmainternal_function_args[['match_sitar_a_form']] <- 
        sigmamatch_sitar_a_form
      sigmainternal_function_args[['match_sitar_d_form']] <- 
        sigmamatch_sitar_d_form
      sigmainternal_function_args[['d_adjustedsi']] <- 
        sigmad_adjustedsi
      sigmainternal_function_args[['xfunsi']] <- 
        sigmaxfunsi
      sigmainternal_function_args[['yfunsi']] <- 
        sigmayfunsi
      sigmainternal_function_args[['xoffset']] <- 
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

      

      
      
      if (verbose) {
        if (ii == 1) {
          setmsgtxt <- paste0("\n Preparing function for sigma")
          if (displayit == 'msg') {
            message(setmsgtxt)
          } else if (displayit == 'col') {
            col <- setcolh
            cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
          }
        }
      }
      
      
      # These are defined for sigma via bsitar()
      # sigmaxsi <- xsi
      # sigmaknots <- knots
      # sigmanknots <- nknots
      
      # These are copied from the mu part
      # Note that idsi is not used and ysi is placeholder
      sigmaysi <- ysi
      sigmaidsi <- idsi
      sigmadatai <- datai
      
      sigmaxs <- sigmaxs
      sigmaids <- sigmaids
      
      
      sigmaget_s_r_funs <-
        prepare_function_sigma(
          x = sigmaxsi,
          y = sigmaysi,
          id = sigmaidsi,
          knots = sigmaknots,
          nknots = sigmanknots,
          data = sigmadatai,
          internal_function_args = sigmainternal_function_args
        )
      
      sigmafunlist[ii] <- sigmaget_s_r_funs[['rcsfun']]
      sigmafunlist_r[[ii]] <- sigmaget_s_r_funs[['r_funs']]
      sigmagq_funs[[ii]] <- sigmaget_s_r_funs[['gq_funs']]
      
    } # if(setsigma_formula_manual) {
    
    
    if(!setsigma_formula_manual) {
      sigmad_adjustedsi <- 'NULL'
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
        "smat_intercept",
        "smat_derivs",
        "smat_centerval",
        "smat_normalize",
        "smat_preH",
        "smat_include_stan",
        "smat_include_path",
        "SplinefunxPre",
        "Splinefunxsuf",
        "SplinefunxR",
        "SplinefunxStan"
      )
    
    
    internal_formula_args <- list()
    internal_formula_args <- mget(internal_formula_args_names)
    
    if (verbose) {
      if (ii == 1) {
        setmsgtxt <- paste0("\n Preparing formula")
        if (displayit == 'msg') {
          message(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
    }
    
   
    
    if(smat == 'rcs') {
      formula_bf <-
        prepare_formula(
          x = xsi,
          y = ysi,
          id = idsi,
          knots = knots,
          nknots = nknots,
          data = datai,
          internal_formula_args = internal_formula_args
        )
    } else if(smat == 'nsp') {
      formula_bf <-
        prepare_formula(
          x = xsi,
          y = ysi,
          id = idsi,
          knots = knots,
          nknots = nknots,
          data = datai,
          internal_formula_args = internal_formula_args
        )
    } else if(smat == 'nsk') {
      formula_bf <-
        prepare_formula(
          x = xsi,
          y = ysi,
          id = idsi,
          knots = knots,
          nknots = nknots,
          data = datai,
          internal_formula_args = internal_formula_args
        )
    }
    
    
    
    
    list_out <- attr(formula_bf, "list_out")
   
    attributes(formula_bf) <- NULL
    
    
    eout <- list2env(list_out)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
   
    group_arg$groupvar <- group_arg_groupvar
    multivariate$rescor <- multivariate_rescor
    univariate_by$by <- univariate_by_by
    covariates_ <- covariates_
    covariates_sigma_ <- covariates_sigma_
    set_higher_levels <- set_higher_levels
    
    sigma_set_higher_levels <- sigma_set_higher_levels
    
    
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
    # formula_bf_sigma <-
    #   prepare_formula_sigma(
    #     x = xsi,
    #     y = ysi,
    #     id = idsi,
    #     knots = knots,
    #     nknots = nknots,
    #     data = datai,
    #     internal_formula_args = internal_formula_args
    #   )

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
        "verbose"
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
        stop(
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
          message(setmsgtxt)
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
    
    
    # Add if(!is.null(a_prior_sdsi)).. when restricting to abcd
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
    
    
    bpriors <- do.call(set_priors_initials, set_priors_initials_agrs)
    
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
       # Somehow now after 24.08.2024, 1:length(eval_what) needed, why?
       
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

          bpriors_str <- do.call(set_priors_initials, 
                                 set_priors_initials_agrs_str, 
                                 envir = set_env_what)

          stanvars_str <- attr(bpriors_str, "stanvars")
          initials_str <- attr(bpriors_str, "initials")
          temp_gr_str_stanvars <- c(temp_gr_str_stanvars, stanvars_str)
          temp_gr_str_priors[[istrx]] <- bpriors_str
        }
        temp_gr_str_priors <- temp_gr_str_priors %>% do.call(rbind, .)
        out <- list(temp_gr_str_priors = temp_gr_str_priors,
                    temp_gr_str_stanvars = temp_gr_str_stanvars,
                    temp_gr_str_inits = temp_gr_str_inits)
      }
      out
    } 
    
    
    temp_gr_str_priors_sd <- list()
    temp_gr_str_stanvars_sd <-  temp_gr_str_inits_sd <- c()
    for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
      set_nlpar_what <- set_randomsi_higher_levsli
      set_env_what   <- environment()
      # n_higher_str   <- length(eval(parse(text = paste0(set_nlpar_what,
      #                                                   "covcoefnames_gr_str")),
      #                               envir = set_env_what))
      
      # 24.08.2024
      # Somehow now after 24.08.2024, n_higher_str <- n_higher_str needed, why?
      
      # n_higher_str   <- n_higher_str - 1
      
      if(set_nlpar_what == "sigma") {
        n_higher_str <- length(eval(parse(text = paste0(set_nlpar_what, "_",
                                                        "hierarchical_gr_names")),
                                    envir = set_env_what))
      } else {
        n_higher_str <- length(eval(parse(text = paste0("",
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
          stop(paste_message)
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
          stop("Length of prior elements for random effect parameter ",
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
    
    higher_level_priors <- temp_gr_str_priors_sd %>% do.call(rbind, .)
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
            
            bpriors_str <- do.call(set_priors_initials, 
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
        
        temp_gr_str_priors <- temp_gr_str_priors %>% do.call(rbind, .)
       
        out <- list(temp_gr_str_priors = temp_gr_str_priors,
                    temp_gr_str_stanvars = temp_gr_str_stanvars,
                    temp_gr_str_inits = temp_gr_str_inits)
        
      out
    }
    

    temp_gr_str_priors_corr <- list()
    temp_gr_str_stanvars_corr <-  temp_gr_str_inits_corr <- c()

    for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
      set_nlpar_what <- set_randomsi_higher_levsli
      set_env_what   <- environment()

      id_higher_str  <- eval(parse(text = paste0(set_nlpar_what, 
                                                 "_str_unique_id")), 
                             envir = set_env_what)
      
      n_higher_str   <- length(id_higher_str)
      
      # 24.08.2024
      # Somehow now after 24.08.2024, 2:length(eval_what) needed, why?
      
      # n_higher_str   <- n_higher_str - 1
      
      if(set_nlpar_what == "sigma") {
        n_higher_str <- length(eval(parse(text = paste0(set_nlpar_what, "_",
                                                        "hierarchical_gr_names")),
                                    envir = set_env_what))
      } else {
        n_higher_str <- length(eval(parse(text = paste0("",
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
          stop("Length of prior elements for random effect parameter ",
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
          org_priors_initials_agrs = set_org_priors_initials_agrs_what
          )
       
        temp_gr_str_priors_corr[[set_randomsi_higher_levsli]] <- 
          out2 $ temp_gr_str_priors
        temp_gr_str_stanvars_corr <- 
          c(temp_gr_str_stanvars_corr, out2 $ temp_gr_str_stanvars)
        temp_gr_str_inits_corr <- 
          c(temp_gr_str_inits_corr,    out2 $ temp_gr_str_inits)
      } 
      
    } 
    
    
    
    higher_level_priors_corr <- temp_gr_str_priors_corr %>% do.call(rbind, .)
    bpriors                  <- rbind(bpriors, higher_level_priors_corr)
   
    
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
    
    initialslist[[ii]] <- initials
    initialslist_s[[ii]] <- initsi
    
    if (initsi == "random") {
      blanketinits <- "random"
    } else if (initsi == "0") {
      blanketinits <- "0"
    } else {
      blanketinits <- "no"
    }
    
    blanketinitslist[[ii]] <- blanketinits
    
    
    if (nys == 1) {
      xoffset_name <- "xoffset"
      knots_name <- "knots"
      fixed_name <- "fixed"
      random_name <- "random"
      groupvar_name <- "groupvar"
      sigma_groupvar_name <- "sigma_groupvar"
      hierarchical_name <- "hierarchical"
      sigma_hierarchical_name <- "sigma_hierarchical"
      xvar_name <- "xvar"
      yvar_name <- "yvar"
      cov_name <- "cov"
      cov_name_sigma <- "cov_sigma"
      d_adjusted_name <- "d_adjusted"
      sigmafixed_name <- "sigmafixed"
      sigmarandom_name <- "sigmarandom"
      sigmad_adjusted_name <- "sigmad_adjusted"
    } else if (nys > 1) {
      xoffset_name <- paste0("xoffset", "_", ysi)
      knots_name <- paste0("knots", "_", ysi)
      fixed_name <- paste0("fixed", "_", ysi)
      random_name <- paste0("random", "_", ysi)
      groupvar_name <- paste0("groupvar", "_", ysi)
      sigma_groupvar_name <- paste0("sigma_groupvar", "_", ysi)
      hierarchical_name <- paste0("hierarchical", "_", ysi)
      sigma_hierarchical_name <- paste0("sigma_hierarchical", "_", ysi)
      xvar_name <- paste0("xvar", "_", ysi)
      yvar_name <- paste0("yvar", "_", ysi)
      cov_name <- paste0("cov", "_", ysi)
      cov_name_sigma <- paste0("cov_sigma", "_", ysi)
      d_adjusted_name <- paste0("d_adjusted", "_", ysi)
      sigmafixed_name <- paste0("sigmafixed", "_", ysi)
      sigmarandom_name <- paste0("sigmarandom", "_", ysi)
      sigmad_adjusted_name <- paste0("sigmad_adjusted", "_", ysi)
    }
    
    
    funlist_r_name <- 'funlist_r'
    funlist_rnamelist[[ii]] <- funlist_r_name
    funlist_rvaluelist[[ii]] <- funlist_r %>% unlist()
  
    
    sigmafunlist_r_name <- 'sigmafunlist_r'
    sigmafunlist_rnamelist[[ii]] <- sigmafunlist_r_name
    sigmafunlist_rvaluelist[[ii]] <- sigmafunlist_r %>% unlist()
    
    
    xoffsetnamelist[[ii]] <- xoffset_name
    xoffsetvaluelist[[ii]] <- xoffset
    
    knotsnamelist[[ii]] <- knots_name
    knotsvaluelist[[ii]] <- knots
    
    fixednamelist[[ii]] <- fixed_name
    fixedvaluelist[[ii]] <-
      strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
    
    sigmafixednamelist[[ii]] <- sigmafixed_name
    sigmafixedvaluelist[[ii]] <-
      strsplit(gsub("\\+", " ", sigmafixedsi), " ")[[1]]
    
    randomnamelist[[ii]] <- random_name
    randomvaluelist[[ii]] <-
      strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
    
    sigmarandomnamelist[[ii]] <- sigmarandom_name
    sigmarandomvaluelist[[ii]] <-
      strsplit(gsub("\\+", " ", sigmarandomsi), " ")[[1]]
    
    
    groupvarnamelist[[ii]] <- groupvar_name
    groupvarvaluelist[[ii]] <- group_arg_groupvar
    
    
    sigma_groupvarnamelist[[ii]] <- sigma_groupvar_name
    sigma_groupvarvaluelist[[ii]] <- sigma_arg_groupvar
    
    hierarchicalvarnamelist[[ii]] <- hierarchical_name
    hierarchicalvarvaluelist[[ii]] <- hierarchical_gr_names
    
    sigma_hierarchicalvarnamelist[[ii]] <- sigma_hierarchical_name
    sigma_hierarchicalvarvaluelist[[ii]] <- sigma_hierarchical_gr_names
    
    xnamelist[[ii]] <- xvar_name
    xvarvaluelist[[ii]] <- xsi
    
    ynamelist[[ii]] <- yvar_name
    yvarvaluelist[[ii]] <- ysi
    
    covnamelist[[ii]] <- cov_name
    covvaluelist[[ii]] <- covariates_
    
    sigmacovnamelist[[ii]] <- cov_name_sigma
    sigmacovvaluelist[[ii]] <- covariates_sigma_
    
    d_adjustednamelist[[ii]] <- d_adjusted_name
    d_adjustedvaluelist[[ii]] <- ept(d_adjustedsi)
    
    
    sigmad_adjustednamelist[[ii]] <- sigmad_adjusted_name
    sigmad_adjustedvaluelist[[ii]] <- ept(sigmad_adjustedsi)
    
    
    # restoring original data
    
    datai[[xsi]] <- ixfuntransformsi(datai[[xsi]] + xoffset) 
    
    # if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
    #   if (xfunsi == "log") {
    #     datai[[xsi]] <- exp(datai[[xsi]] + xoffset)
    #   } else if (xfunsi == "sqrt") {
    #     datai[[xsi]] <- (datai[[xsi]] + xoffset) ^ 2
    #   }
    # } else if (is.null(xfunsi[[1]][1]) | xfunsi == "NULL") {
    #   datai[[xsi]] <- (datai[[xsi]] + xoffset)
    # }
    
    
    
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA"))
      dataout <- rbind(dataout, datai)
    else
      dataout <- datai

    # 20.09.2024
    # remove sigmaxsi if not using
    if(!setsigma_formula_manual) datai[[sigmaxsi]] <- NULL
    
      
      
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA"))
      uvarbyTF <- TRUE
    else
      uvarbyTF <- FALSE
  }  
  
  
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
      message(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  
  
  brmsdata <- dataout
  brmspriors <- priorlist
  
  # brmspriorsx <<- brmspriors

  # IMP - brms does not allow different lb for sd params (e.e, all to be NA)
  # Error: Conflicting boundary information for coefficients of class 'sd'.
  # Because prior function automatically sets lb 0 for positive priors 
  # such as exponential the following is need (again done at line 4753 )
  
  brmspriors <- brmspriors %>% 
    dplyr::mutate(lb = dplyr::if_else(class == 'sd', NA, lb))
  brmspriors <- brmspriors %>% 
    dplyr::mutate(ub = dplyr::if_else(class == 'sd', NA, ub))
  
  
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
  
  
  if(!setsigma_formula_manual) {
    fun_scode <- paste(funlist, collapse = "\n")
  }
  
  if(setsigma_formula_manual) {
    fun_scode <- paste(funlist, sigmafunlist, collapse = "\n")
  }
  
  fun_scode <- paste0("functions {", "\n", fun_scode, "\n", "}")

  if(!setsigma_formula_manual) {
    bstanvars <-
      brms::stanvar(scode = paste(funlist, collapse = "\n"), block = "function")
  }
  
  if(setsigma_formula_manual) {
    bstanvars <-
      brms::stanvar(scode = paste(funlist, sigmafunlist, collapse = "\n"), 
                    block = "function")
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
            stop("length of correlation vector must be ",
                 nc,
                 "\n, ",
                 ", but found ",
                 length(x))
          }
        }
        m_lower[lower.tri(m_lower, diag = FALSE)] <- x
        m_upper <- t(m_lower)
        M <- m_lower + m + m_upper
        M
      }
      
      tt_names <- apply(combn(t_names, 2), 2, paste, collapse = "_")
      tt_dims <- sum(d_comb)
      tt_nc <- (tt_dims * (tt_dims - 1) / 2)
      tt_12 <- create_cor_mat(tt_dims, rep(0, tt_nc))
      colnames(tt_12) <- rownames(tt_12) <- t_names
      tt_ll <- tt_12[lower.tri(tt_12)]
      names(tt_ll) <-
        apply(combn(colnames(tt_12), 2), 2, paste, collapse = "_")
      tt_ll[names(l_comb)] <- l_comb
      tt_ll[!names(tt_ll) %in% names(l_comb)] <- 0
      brmsinits[[keys[1]]] <- create_cor_mat(tt_dims, tt_ll)
  
      c_it <- "z_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      brmsinits[[keys[1]]] <- do.call(rbind, temppp)
    } else if (multivariate$mvar &
               (multivariate$cor == "un" | multivariate$cor ==
                "un_s") &
               !any(grepl("^L_", names(brmsinits))))
    {
     
      c_it <- "sd_"
      brmsinits_names <- names(brmsinits)
      brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', 
                                                brmsinits_names)]
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      temppp <- unlist(unname(temppp))
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      for (sdi in 1:length(temppp)) {
        brmsinits[[paste0(c_it, sdi)]] <- temppp[sdi] %>%
          unname()
      }
    }  
    
    # keep only one Lrescor
    if (multivariate$mvar & multivariate$rescor) {
      c_it <- "Lrescor_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      brmsinits[["Lrescor"]] <- temppp[[1]]
    }
    
    if ((multivariate$mvar &
         multivariate$cor == "diagonal") |
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
  } else if (all(sapply("0", grepl, initialslist_s))) {
    brmsinits <- "0"
    brmsinits_r <- ept(init_rsi)
    brmsinits_ <- NULL
  } else {
    brmsinits <- brmsinits
    brmsinits_r <- NULL
    brmsinits_ <- ""
  }
  
  
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
        if (is.character(jitter_init_beta))
          jitter_init_beta <- ept(jitter_init_beta)
        if (is.character(jitter_init_sd))
          jitter_init_sd <- ept(jitter_init_sd)
        if (is.character(jitter_init_cor))
          jitter_init_cor <- ept(jitter_init_cor)
        
        if (!is.null(jitter_init_beta) &
            !is.numeric(jitter_init_beta)) {
          stop("Argument jitter_init_beta should be NULL or a numeric value")
        }
        if (!is.null(jitter_init_sd) &
            !is.numeric(jitter_init_sd)) {
          stop("Argument jitter_init_sd should be NULL or a numeric value")
        }
        if (!is.null(jitter_init_cor) &
            !is.numeric(jitter_init_cor)) {
          stop("Argument jitter_init_cor should be NULL or a numeric value")
        }
        
        jitter_x <- function(x, a, digits) {
          x <- unname(x)
          col <- c()
          for (i in 1:length(x)) {
            amount <- abs(x[i]) * a
            col <- c(col, jitter(x[i], factor = 1, amount = amount))
          }
          col <- round(col, digits)
          col
        }
        
        jitter_mat <- function(x, a, digits) {
          mat_out <- x
          x <- x[lower.tri(x)]
          col <- c()
          for (i in 1:length(x)) {
            amount <- abs(x[i]) * a
            col <- c(col, jitter(x[i], factor = 1, amount = amount))
          }
          col <- round(col, digits)
          col <- ifelse(col > 1, 1, col)
          col <- ifelse(col < -1, 1, col)
          mat_out[lower.tri(mat_out)] <-
            mat_out[upper.tri(mat_out)] <- col
          return(mat_out)
        }
        
        eval_inits <- c()
        for (i_init in names(inits)) {
          if (grepl("^b_", i_init)) {
            if (!is.null(jitter_init_beta)) {
              values_i <-
                jitter_x(inits[[i_init]], jitter_init_beta, digits = digits)
            } else if (is.null(jitter_init_beta)) {
              values_i <- inits[[i_init]]
              values_i <- round(values_i, digits)
            }
            eval_inits[[i_init]] <- values_i
          } else if (grepl("^sd_", i_init)) {
            if (!is.null(jitter_init_sd)) {
              values_i <- jitter_x(inits[[i_init]], jitter_init_sd, digits)
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
      
    
    
    # 24.08.2024
    if(get_priors) {
      tempriorstr <- brms::get_prior(formula = bformula,
                               stanvars = bstanvars,
                               prior = temp_prior,
                               data = brmsdata)
      return(tempriorstr)
    }

    

    # 24.08.2024
    if(setsigma_formula_manual) {
      temp_prior <- temp_prior %>% dplyr::filter(class != "sigma")
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
        ##########
        
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
              htxi_5 <- gsub(SplineFun_name, paste0(SplineFun_name, 
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
        temp_stancode2cp <- edit_scode_ncp_to_cp(temp_stancode2, 
                                                 genq_only = FALSE, 
                                                 normalize = normalize)
        
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
        digits = 4
      )
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
      stop('All parameters must have the same number of parameters')
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
        stop('All parameters must have the covariate form as ~0+')
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
        
        cores_ <- eval(setarguments$cores)
        threads_ <- eval(setarguments$threads)
        
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
        } else if(!is.null(getOption('mc.cores')) &
                  cores_ != "maximise" &
                  cores_ != "optimize") {
          max.cores <- getOption('mc.cores')
        } else {
          max.cores <- eval(setarguments$cores)
        }
        setarguments$cores <-  max.cores
        
        
        
        if(!is.list(threads_)) {
          if( is.character(threads_) & threads_ == "maximise") {
            max.threads <- 
              as.numeric(future::availableCores(methods = "system", omit = 0))
            if(max.threads < 1) max.threads <- 1
          } else if( is.character(threads_) & threads_ == "optimize") {
            max.threads <- 
              as.numeric(future::availableCores(methods = "system", omit = 1))
            if(max.threads < 1) max.threads <- 1
            max.threads <- floor(max.threads /  eval(setarguments$chains))
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
            stanc_options = list("O1")
            # , cpp_options = list(#'STAN_CPP_OPTIMS=true',
            #                      # 'CXXFLAGS = -O2',
            #                      # 'STANCFLAGS+= --warn-pedantic --O0',
            #                      'STAN_NO_RANGE_CHECKS=true'
            #                      )
            )
        }
      }
      
      
      if(verbose) {
        message(setarguments$stan_model_args)
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
      message(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  
  
  brmsdots_ <- list(...)
  
  # 24.08.2024
  getdotslistnames <- c("match_sitar_a_form",
                         "match_sitar_d_form",
                         "sigmamatch_sitar_a_form",
                         "displayit", "setcolh", "setcolb",
                        "smat"
                        )
  
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
      seed,
      verbose,
      setarguments = brms_arguments,
      brmsdots = brmsdots_
    )
 
  if(!is.null(custom_family)) {
    brm_args$family <- custom_family
  }
  
  if (verbose) {
    setmsgtxt <- paste0("\n Fitting model")
    if (displayit == 'msg') {
      message(setmsgtxt)
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
    p1 <- cc %>% do.call(rbind, .)
    p2 <- zz %>% do.call(rbind, .)
    p1p2 <- rbind(p1, p2)
    p1p2
  }
  
  
  
  
  if(set_higher_levels) {
    brmspriors_sdcor <- brmspriors %>% 
      dplyr::filter(.data$class == 'sd' | .data$class == 'cor')
    brmspriors_sdcor_gr <- brmspriors_sdcor$group
    
    brmsfit_sdcor <- do.call(brms::get_prior, brm_args) %>% 
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
    stop("Amongst 'set_self_priors', 
         'add_self_priors' and 'set_replace_priors' arguments,",
         "\n ",
         " only one can be specified at a time")
  }
  
  if(get_priors & get_priors_eval & validate_priors & 
     get_stancode & get_standata & get_formula & get_stanvars) {
    stop("Amongst 'get_priors' 'get_priors_eval', 'validate_priors' ",
         "\n ",
         "'get_stancode', 'get_standata', 'get_formula', 'get_stanvars' ",
         "\n ",
         " arguments, only one can be set to TRUE at a time")
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
  
  
  

  
  
  
  

  # 24.08.2024
  if(setsigma_formula_manual) {
    brmspriors <- brmspriors %>% dplyr::filter(class != "sigma")
  }
  
 
  
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
  
  
  
  scode_final  <- do.call(brms::make_stancode, brm_args)
  sdata  <- do.call(brms::make_standata, brm_args)
  
  
  if(parameterization == 'cp') {
    scode_final <- edit_scode_ncp_to_cp(scode_final, 
                                        genq_only = FALSE, 
                                        normalize = normalize)
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
  
  
  if(!exe_model_fit) {
    if(get_priors) {
      return(do.call(brms::get_prior, brm_args))
    } else if(get_standata) {
      return(do.call(brms::make_standata, brm_args))
    } else if(get_stancode) {
      return(scode_final)
    } else if(get_priors_eval) {
      return(get_priors_eval_out)
    } else if(validate_priors) {
      return(do.call(brms::validate_prior, brm_args))
    } else if(get_init_eval) {
      return(brm_args$init)
    } else if(get_formula) {
      return(brm_args$formula)
    } else if(get_stanvars) {
      return(brm_args$stanvars)
    }
  } 
  
  
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
          stop("Custom initials specified via 'init_custom' argument must",
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
 
    
    # Set refresh based on thin argument
    if(!is.null(brm_args$refresh) & brm_args$thin > 1) {
      brm_args$refresh <- 
        ceiling((brm_args$refresh * brm_args$thin) / brm_args$thin)
    }
    
    brm_args$refresh <- NULL
    
    
    # Get and evaluate file argument
    # This to save object 'file' at the end with model info
    get_file <- brm_args$file
    get_file_refit <- brm_args$file_refit
    
    get_file_compress <- brm_args$file_compress
    
    get_write_brmsfit          <-
      utils::getFromNamespace("write_brmsfit", "brms")
    
    get_read_brmsfit          <-
      utils::getFromNamespace("read_brmsfit", "brms")
    
    get_file_refit_options          <-
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

    if(brm_args$backend == "cmdstanr" |
       !is.null(pathfinder_args) | 
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
    
    
    

    if(fit_edited_scode) {
      if(brm_args$backend == "cmdstanr") {
         brmsfit <- brms_via_cmdstanr(scode_final, sdata, brm_args, 
                                      brms_arguments,
                                      pathfinder_args = pathfinder_args,
                                      pathfinder_init = pathfinder_init)
      }
      if(brm_args$backend == "rstan") {
        brmsfit  <- brms_via_rstan(scode_final, sdata, brm_args, 
                                   brms_arguments)
      }
    } else if(!fit_edited_scode) {
      if(!is.null(pathfinder_args) | pathfinder_init) {
        brmsfit <- brms_via_cmdstanr(scode_final, sdata, brm_args,
                                     brms_arguments,
                                     pathfinder_args = pathfinder_args,
                                     pathfinder_init = pathfinder_init)
      } else {
        brmsfit <- do.call(brms::brm, brm_args)
      }
    } # if(fit_edited_scode) {
    
    
    
    
    if(brm_args$backend == "mock") {
      brmsfit <- do.call(brms::brm, brm_args)
    }
    

    # Add class attributes and the model info for post-processing
    attr(brmsfit, 'class') <- c(attr(brmsfit, 'class'), 'bgmfit')
    
    model_info <- list()
    model_info[['emodel']] <- scode_final
    model_info[['parameterization']] <- parameterization
    model_info[['d_adjusted']] <- d_adjusted
    
    for (i in 1:length(funlist_rnamelist)) {
      model_info[[funlist_rnamelist[[i]]]] <- funlist_rvaluelist[[i]]
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
    
    for (i in 1:length(xxfunnamelist)) {
      model_info[[xxfunnamelist[[i]]]] <- xxfunvaluelist[[i]]
    }
    
    for (i in 1:length(yyfunnamelist)) {
      model_info[[yyfunnamelist[[i]]]] <- yyfunvaluelist[[i]]
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
    for (i in 1:length(ixfuntransformnamelist)) {
      model_info[[ixfuntransformnamelist[[i]]]] <- ixfuntransformvaluelist[[i]]
    }
    
    for (i in 1:length(yfuntransformnamelist)) {
      model_info[[yfuntransformnamelist[[i]]]] <- yfuntransformvaluelist[[i]]
    }
    for (i in 1:length(iyfuntransformnamelist)) {
      model_info[[iyfuntransformnamelist[[i]]]] <- iyfuntransformvaluelist[[i]]
    }
    
    
    for (i in 1:length(sigmaxfuntransformnamelist)) {
      model_info[[sigmaxfuntransformnamelist[[i]]]] <- 
        sigmaxfuntransformvaluelist[[i]]
    }
    for (i in 1:length(isigmaxfuntransformnamelist)) {
      model_info[[isigmaxfuntransformnamelist[[i]]]] <- 
        isigmaxfuntransformvaluelist[[i]]
    }
    
    
    model_info[['StanFun_name']] <- SplineFun_name
    model_info[['multivariate']] <- multivariate$mvar
    model_info[['univariate_by']] <- univariate_by$by
    model_info[['nys']] <- nys
    model_info[['ys']] <- ys
    model_info[['xs']] <- xs
    model_info[['ids']] <- ids
    model_info[['dfs']] <- dfs
    model_info[['xfuns']] <- xfuns
    model_info[['yfuns']] <- yfuns
    model_info[['outliers']] <- outliers
    model_info[['bgmfit.data']] <- data.org.in
    model_info[['call.full.bgmfit']] <- call.full
    model_info[['call.bgmfit']] <- mcall_
    model_info[['brms_arguments_list']] <- brms_arguments_list
    model_info[['select_model']] <- select_model
    model_info[['decomp']] <- decomp
    model_info[['fun_scode']] <- fun_scode
    model_info[['envir']] <- enverr.
    model_info[['include_fun_names']] <- include_fun_names
    
    if(setsigma_formula_manual) {
      model_info[['sigmaStanFun_name']] <- sigmaSplineFun_name
      model_info[['sigmaxs']] <- sigmaxs
      model_info[['sigmaids']] <- sigmaids
      model_info[['sigmadfs']] <- sigmadfs
      model_info[['sigmaxfuns']] <- sigmaxfuns
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
      for (i in 1:length(sigmaxfunnamelist)) {
        model_info[[sigmaxfunnamelist[[i]]]] <- sigmaxfunvaluelist[[i]]
      }
      for (i in 1:length(sigmaxxfunnamelist)) {
        model_info[[sigmaxxfunnamelist[[i]]]] <- sigmaxxfunvaluelist[[i]]
      }
      for (i in 1:length(sigmacovnamelist)) {
        model_info[[sigmacovnamelist[[i]]]] <- sigmacovvaluelist[[i]]
      }
    } # if(setsigma_formula_manual) {
    

    brmsfit$model_info <- model_info
    
    environment(brmsfit$formula) <- enverr.
    
    # Now message moved to the expose_model_functions()
    if (expose_function & !brm_args$empty) {
      # if (verbose) {
      #   setmsgtxt <-
      #     paste0("\n Exposing Stan functions for post-processing\n")
      #   if (displayit == 'msg') {
      #     message(setmsgtxt)
      #   } else if (displayit == 'col') {
      #     col <- setcolh
      #     cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
      #   }
      # }
      # if (!verbose) {
      #   setmsgtxt <-
      #     paste0("\n Exposing Stan functions for post-processing..\n")
      #   message(setmsgtxt)
      # }
      
      brmsfit <- expose_model_functions(model = brmsfit, 
                                      scode = fun_scode,
                                      expose = TRUE, 
                                      select_model = NULL,
                                      returnobj = TRUE,
                                      vectorize = FALSE,
                                      verbose = TRUE,
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
                                      envir = NULL)
      brmsfit$model_info[['expose_method']] <- 'R'
    } 
    
    if (verbose) {
      setmsgtxt <- paste0("\nModel Fitting complete")
      if (displayit == 'msg') {
        message(setmsgtxt)
      } else if (displayit == 'col') {
        col <- setcolh
        cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
      }
    }
    
    if (!is.null(get_file)) {
      brmsfit <- get_write_brmsfit(brmsfit, get_file, 
                               compress = get_file_compress)
    }
    
    return(brmsfit)
  } # exe_model_fit
  
}

