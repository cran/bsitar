


#' An internal function to set priors and initials values
#'
#'@description The \code{set_priors_initials}) sets priors and initials values
#'  which are passed from the [bsitar::bsitar()] function to
#'  \code{set_priors_initials}. For univariate-by- subgroup model (specified by
#'  using the \code{univariate_by}) and multivariate model (specified by using
#'  the \code{multivariate}), each argument is automatically matched with the
#'  sub-model(s).
#'
#'@param a_prior_beta Specify priors for the fixed effect parameter, \code{a}.
#'  See [bsitar::bsitar()] for details.
#'
#'@param b_prior_beta Specify priors for the fixed effect parameter, \code{b}.
#'  See [bsitar::bsitar()] for details.
#'
#'@param c_prior_beta Specify priors for the fixed effect parameter, \code{c}.
#'  See [bsitar::bsitar()] for details.
#'
#'@param d_prior_beta Specify priors for the fixed effect parameter, \code{d}.
#'  See [bsitar::bsitar()] for details.
#'
#'@param e_prior_beta Specify priors for the fixed effect parameter, \code{e}.
#'  See [bsitar::bsitar()] for details.
#'
#'@param f_prior_beta Specify priors for the fixed effect parameter, \code{f}.
#'  See [bsitar::bsitar()] for details.
#'
#'@param s_prior_beta Specify priors for the fixed effect parameter, \code{s}.
#'  See [bsitar::bsitar()] for details.
#'
#'@param a_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{a}. See [bsitar::bsitar()] for details.
#'
#'@param b_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{b}. See [bsitar::bsitar()] for details.
#'
#'@param c_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{c}. See [bsitar::bsitar()] for details.
#'
#'@param d_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{d}. See [bsitar::bsitar()] for details.
#'
#'@param e_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{e}. See [bsitar::bsitar()] for details.
#'
#'@param f_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{f}. See [bsitar::bsitar()] for details.
#'
#'@param s_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{s}. See [bsitar::bsitar()] for details.
#'
#'@param a_prior_sd Specify prior on the standard deviation (sd) of random
#'  effect parameter, \code{a}. See [bsitar::bsitar()] for details.
#'
#'@param b_prior_sd Specify prior on the standard deviation (sd) of random
#'  effect parameter, \code{b}. See [bsitar::bsitar()] for details.
#'
#'@param c_prior_sd Specify prior on the standard deviation (sd) of random
#'  effect parameter, \code{c}. See [bsitar::bsitar()] for details.
#'
#'@param d_prior_sd Specify prior on the standard deviation (sd) of random
#'  effect parameter, \code{d}. See [bsitar::bsitar()] for details.
#'
#'@param e_prior_sd Specify prior on the standard deviation (sd) of random
#'  effect parameter, \code{e}. See [bsitar::bsitar()] for details.
#'
#'@param f_prior_sd Specify prior on the standard deviation (sd) of random
#'  effect parameter, \code{f}. See [bsitar::bsitar()] for details.
#'
#'@param a_cov_prior_sd Specify prior on the standard deviation (sd) for the
#'  covariate(s) included in the random effect parameter, \code{a}. See
#'  [bsitar::bsitar()] for details.
#'
#'@param b_cov_prior_sd Specify prior on the standard deviation (sd) for the
#'  covariate(s) included in the random effect parameter, \code{b}. See
#'  [bsitar::bsitar()] for details.
#'
#'@param c_cov_prior_sd Specify prior on the standard deviation (sd) for the
#'  covariate(s) included in the random effect parameter, \code{c}. See
#'  [bsitar::bsitar()] for details.
#'
#'@param d_cov_prior_sd Specify prior on the standard deviation (sd) for the
#'  covariate(s) included in the random effect parameter, \code{d}. See
#'  [bsitar::bsitar()] for details.
#'
#'@param e_cov_prior_sd Specify prior on the standard deviation (sd) for the
#'  covariate(s) included in the random effect parameter, \code{w}. See
#'  [bsitar::bsitar()] for details.
#'
#'@param f_cov_prior_sd Specify prior on the standard deviation (sd) for the
#'  covariate(s) included in the random effect parameter, \code{f}. See
#'  [bsitar::bsitar()] for details.
#'
#'@param sigma_prior_beta Specify prior on the the distributional fixed effect
#'  parameter, \code{sigma}. See [bsitar::bsitar()] for details.
#'
#'@param sigma_cov_prior_beta Specify prior on the covariate(s) included in the
#'  distributional fixed effect parameter, \code{sigma}.  See [bsitar::bsitar()]
#'  for details.
#'
#'@param sigma_prior_sd Specify prior on the standard deviation (sd) of
#'  distributional random effect parameter, \code{sigma}. See [bsitar::bsitar()]
#'  for details.
#'
#'@param sigma_cov_prior_sd Specify prior on the standard deviation (sd) of
#'  covariate(s) included in the distributional random effect parameter,
#'  \code{sigma}. See [bsitar::bsitar()] for details.
#'
#'@param gr_prior_cor Specify prior on the correlation for group level random
#'  effect parameters. See [bsitar::bsitar()] for details.
#'
#'@param sigma_prior_cor Specify prior on the correlation for distribution level
#'  random effect parameters. See [bsitar::bsitar()] for details.
#'
#'@param rsd_prior_sigma Specify prior on the residual standared deviation
#'  parameter, \code{sigma}, See [bsitar::bsitar()] for details,
#'
#'@param dpar_prior_sigma Specify prior on the distributional parameter,
#'  \code{sigma} (which is same as residual standared deviation for Gaussian
#'  distribution). See [bsitar::bsitar()] for details,
#'
#'@param dpar_cov_prior_sigma Specify prior for the covariate(s) included in the
#'  distributional parameter, \code{sigma} (which is same as residual standard
#'  deviation for Gaussian distribution).See [bsitar::bsitar()] for details,
#'
#'@param autocor_prior_acor Specify priors on the the autocorrelation parameters
#'  \code{ar}, \code{ma} and \code{arma}. See [bsitar::bsitar()] for details,
#'
#'@param autocor_prior_unstr_acor Specify priors on the the unstructured
#'  autocorrelation parameter. See [bsitar::bsitar()] for details,
#'
#'@param mvr_prior_rescor Specify priors on the the residual correlation
#'  parameter for multivariate model. See [bsitar::bsitar()] for details,
#'
#'@param prior_data An optional argument (as named list) that can pass value for
#'  prior. See [bsitar::bsitar()] for details,
#'
#'@param prior_data_internal An internal data frame (named list) used to pass on
#'  the relevant information on priors from the [bsitar::bsitar()] function to
#'  the \code{set_priors_initials}.
#'
#'@param prior_args_internal An internal argument list that is passed from the
#'  [bsitar::bsitar()] function to the \code{set_priors_initials} and is used
#'  for setting the priors.
#'
#'@param init_arguments A list containing all the init arguments specified in
#'  the [bsitar::bsitar()] function and now passed on to the
#'  \code{set_priors_initials}.
#'
#'@param init_data An optional data argument (named list) used to pass initial
#'  values. See [bsitar::bsitar()] function, \code{prior_data} for details.
#'
#'@param init_data_internal An internal data frame (named list) to pass on the
#'  relevant information on initials from the [bsitar::bsitar()] function to the
#'  \code{set_priors_initials}.
#'
#'@param init_args_internal An internal argument list that is passed from the
#'  [bsitar::bsitar()] function to the \code{set_priors_initials} and is used
#'  for setting the initials.
#'
#'@param custom_order_prior_str An internal argument that is passed from the
#'  [bsitar::bsitar()] function to the \code{set_priors_initials} when setting
#'  the priors for the model with hierarchy level 3 and beyond. See
#'  [bsitar::bsitar()] for details,
#'
#'@return An object of class \code{brmsprior} (See \code{brmsprior}). In
#'  addition to the priors, the returned object contains a list of initial
#'  values.
#'  
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'  
#' @keywords internal
#' @noRd
#'
set_priors_initials <- function(a_prior_beta,
                                b_prior_beta,
                                c_prior_beta,
                                d_prior_beta,
                                e_prior_beta,
                                f_prior_beta,
                                g_prior_beta,
                                h_prior_beta,
                                i_prior_beta,
                                s_prior_beta,
                                a_cov_prior_beta,
                                b_cov_prior_beta,
                                c_cov_prior_beta,
                                d_cov_prior_beta,
                                e_cov_prior_beta,
                                f_cov_prior_beta,
                                g_cov_prior_beta,
                                h_cov_prior_beta,
                                i_cov_prior_beta,
                                s_cov_prior_beta,
                                a_prior_sd,
                                b_prior_sd,
                                c_prior_sd,
                                d_prior_sd,
                                e_prior_sd,
                                f_prior_sd,
                                g_prior_sd,
                                h_prior_sd,
                                i_prior_sd,
                                s_prior_sd,
                                a_cov_prior_sd,
                                b_cov_prior_sd,
                                c_cov_prior_sd,
                                d_cov_prior_sd,
                                e_cov_prior_sd,
                                f_cov_prior_sd,
                                g_cov_prior_sd,
                                h_cov_prior_sd,
                                i_cov_prior_sd,
                                s_cov_prior_sd,
                                gr_prior_cor,
                                sigma_prior_cor,
                                sigma_prior_beta,
                                sigma_cov_prior_beta,
                                sigma_prior_sd,
                                sigma_cov_prior_sd,
                                rsd_prior_sigma,
                                dpar_prior_sigma,
                                dpar_cov_prior_sigma,
                                autocor_prior_acor,
                                autocor_prior_unstr_acor,
                                mvr_prior_rescor,
                                prior_data             = NULL,
                                prior_data_internal    = NULL,
                                prior_args_internal    = NULL,
                                init_arguments         = NULL,
                                init_data              = NULL,
                                init_data_internal     = NULL,
                                init_args_internal     = NULL,
                                custom_order_prior_str = NULL) {
  
  # Initiate non formalArgs()
  resp <- NULL;
  autocor_formi <- NULL;
  randomsi <- NULL;
  sigma_formula_grsi <- NULL;
  a_formulasi <- NULL;
  b_formulasi <- NULL;
  c_formulasi <- NULL;
  d_formulasi <- NULL;
  e_formulasi <- NULL;
  f_formulasi <- NULL;
  g_formulasi <- NULL;
  h_formulasi <- NULL;
  i_formulasi <- NULL;
  s_formulasi <- NULL;
  select_model <- NULL;
  fixedsi <- NULL;
  a_formula_grsi <- NULL;
  b_formula_grsi <- NULL;
  c_formula_grsi <- NULL;
  d_formula_grsi <- NULL;
  e_formula_grsi <- NULL;
  f_formula_grsi <- NULL;
  g_formula_grsi <- NULL;
  h_formula_grsi <- NULL;
  i_formula_grsi <- NULL;
  s_formula_grsi <- NULL;
  acovcoefnames <- NULL;
  bcovcoefnames <- NULL;
  ccovcoefnames <- NULL;
  dcovcoefnames <- NULL;
  ecovcoefnames <- NULL;
  fcovcoefnames <- NULL;
  gcovcoefnames <- NULL;
  hcovcoefnames <- NULL;
  icovcoefnames <- NULL;
  scovcoefnames <- NULL;
  acovcoefnames_gr <- NULL;
  bcovcoefnames_gr <- NULL;
  ccovcoefnames_gr <- NULL;
  dcovcoefnames_gr <- NULL;
  ecovcoefnames_gr <- NULL;
  fcovcoefnames_gr <- NULL;
  gcovcoefnames_gr <- NULL;
  hcovcoefnames_gr <- NULL;
  icovcoefnames_gr <- NULL;
  scovcoefnames_gr <- NULL;
  dpar_formulasi <- NULL;
  sigma_formulasi <- NULL;
  sigmacovcoefnames <- NULL;
  sigmacovcoefnames_gr <- NULL;
  nabcrei <- NULL;
  dparcovcoefnames <- NULL;
  sigma_arg_groupvar <- NULL;
  group_arg_groupvar <- NULL;
  autocor_formi <- NULL;
  acovcoefnames <- NULL;
  bcovcoefnames <- NULL;
  ccovcoefnames <- NULL;
  dcovcoefnames <- NULL;
  ecovcoefnames <- NULL;
  fcovcoefnames <- NULL;
  gcovcoefnames <- NULL;
  hcovcoefnames <- NULL;
  icovcoefnames <- NULL;
  scovcoefnames <- NULL;
  acovcoefnames_gr <- NULL;
  bcovcoefnames_gr <- NULL;
  ccovcoefnames_gr <- NULL;
  dcovcoefnames_gr <- NULL;
  ecovcoefnames_gr <- NULL;
  fcovcoefnames_gr <- NULL;
  gcovcoefnames_gr <- NULL;
  hcovcoefnames_gr <- NULL;
  icovcoefnames_gr <- NULL;
  scovcoefnames_gr <- NULL;
  sigmacovcoefnames <- NULL;
  sigmacovcoefnames_gr <- NULL;
  nabcrei <- NULL;
  resp <- NULL;
  nys <- NULL;
  fixedsi <- NULL;
  randomsi <- NULL;
  nabci <- NULL;
  nabcrei <- NULL;
  ii <- NULL;
  N_J_all <- NULL;
  dpar_formulasi <- NULL;
  initsi <- NULL;
  seed <- NULL;
  cortimeNlags_var <- NULL;
  cortimeNlags <- NULL;
  verbose <- NULL;
  sigma_formulasi <- NULL;
  s_formulasi <- NULL;
  sigma_formula_grsi <- NULL;
  nys <- NULL;
  cortimeNlags <- NULL;
  . <- NULL;
  d_adjustedsi <- NULL;
  
  
  
  eout <- list2env(prior_data_internal)
  for (eoutii in names(eout)) {
    assign(eoutii, eout[[eoutii]])
  }
  
  
  eout <- list2env(prior_args_internal)
  for (eoutii in names(eout)) {
    assign(eoutii, eout[[eoutii]])
  }
  
  
  
  
  # Depending on select_model, assign null values to all not part of the model
  for (set_randomsi_higher_levsli in c(letters[1:20])) {
    set_nlpar_what <- set_randomsi_higher_levsli
    if(!exists(paste0(set_randomsi_higher_levsli, 'form'))) {
      assign(paste0(set_nlpar_what, '_prior_beta'), NULL)
      assign(paste0(set_nlpar_what, '_cov_prior_beta'), NULL)
      assign(paste0(set_nlpar_what, '_prior_sd'), NULL)
      assign(paste0(set_nlpar_what, '_cov_prior_sd'), NULL)
    } else if(exists(paste0(set_randomsi_higher_levsli, 'form'))) {
      if(is.null(ept(paste0(set_randomsi_higher_levsli, 'form')))) {
        assign(paste0(set_nlpar_what, '_prior_beta'), NULL)
        assign(paste0(set_nlpar_what, '_cov_prior_beta'), NULL)
        assign(paste0(set_nlpar_what, '_prior_sd'), NULL)
        assign(paste0(set_nlpar_what, '_cov_prior_sd'), NULL)
      } # if(is.null(ept(paste0('f', 'form')))) {
    }
  } # for (set_randomsi_higher_levsli in c(letters[1:20])) {
  
  
  
  normalize <- ept(normalize)
  
  if (resp == "")
    resp_ <- ""
  if (resp != "")
    resp_ <- paste0("_", resp)
  
  
  if (is.null(autocor_formi)) {
    autocor_prior_acor <- NULL
    autocor_prior_unstr_acor <- NULL
  }
  
  if (!is.null(autocor_formi)) {
    if(grepl("unstr(", autocor_formi, fixed = T)) 
      autocor_prior_acor <- NULL
    if(!grepl("unstr(", autocor_formi, fixed = T)) 
      autocor_prior_unstr_acor <- NULL
  }
  
  
  getArgNames <-
    function(value)
      formalArgs(deparse(substitute(value)[[1]]))
  
  mvar <- multivariate$mvar
  
  if (!is.null(group_arg$cor)) {
    if (group_arg$cor == "un")
      abccorr <- TRUE
    if (group_arg$cor == "diagonal")
      abccorr <- FALSE
  } else {
    group_arg$cor <- "un"
    abccorr <- TRUE
  }
  
  if (group_arg$cor != "diagonal" & !group_arg$cor != "diagonal") {
    #  abccorr <- FALSE
  }
  
  if(group_arg$cor == "no") {
    abccorr <- FALSE
  }
  
  if (!is.null(group_arg$by)) {
    group_arg$by <- group_arg$by
  } else {
    group_arg$by <- NULL
  }
  
  if (!is.null(group_arg$cov)) {
    group_arg$cov <- group_arg$cov
  } else {
    group_arg$cov <- NULL
  }
  
  if (!is.null(group_arg$dist)) {
    group_arg$dist <- group_arg$dist
  } else {
    group_arg$dist <- 'gaussian'
  }
  
  if (!is.null(group_arg$verbose)) {
    group_arg$verbose <- group_arg$verbose
  } else {
    group_arg$verbose <- FALSE
  }
  
  # this on 9 5 23 to accomodate random = ''
  if(randomsi == "") gr_prior_cor <- NULL
  if(!abccorr)       gr_prior_cor <- NULL
  
  
  
  
  
  
  
  if (!is.null(sigma_group_arg$cor)) {
    if (sigma_group_arg$cor == "un")
      sigmacorr <- TRUE
    if (sigma_group_arg$cor == "diagonal")
      sigmacorr <- FALSE
  } else {
    sigma_group_arg$cor <- "un"
    sigmacorr <- TRUE
  }
  
  if (!is.null(sigma_group_arg$by)) {
    sigma_group_arg$by <- sigma_group_arg$by
  } else {
    sigma_group_arg$by <- NULL
  }
  
  if (!is.null(sigma_group_arg$cov)) {
    sigma_group_arg$cov <- sigma_group_arg$cov
  } else {
    sigma_group_arg$cov <- NULL
  }
  
  if (!is.null(sigma_group_arg$dist)) {
    sigma_group_arg$dist <- sigma_group_arg$dist
  } else {
    sigma_group_arg$dist <- 'gaussian'
  }
  
  if (!is.null(sigma_group_arg$verbose)) {
    sigma_group_arg$verbose <- sigma_group_arg$verbose
  } else {
    sigma_group_arg$verbose <- FALSE
  }
  
  
  # this on 9 5 23 to accomodate no random sigam 
  if(is.null(sigma_formula_grsi)) sigma_prior_cor <- NULL
  if(!sigmacorr)                  sigma_prior_cor <- NULL
  
  
  
  # If group and sigma ids are same, then sigma_prior_cor NULL otherwise 
  # duplicate priors error
  
  if(identical(sigma_group_arg$groupvar,
               group_arg$groupvar)) {
    sigma_prior_cor <- NULL
  }
  
  
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (!is.null(univariate_by$cor)) {
      if (univariate_by$cor == "un")
        uvarabccorr <- TRUE
      if (univariate_by$cor == "diagonal")
        uvarabccorr <- FALSE
    } else {
      univariate_by$cor <- "un"
      uvarabccorr <- TRUE
    }
    if (is.null(univariate_by$verbose))
      univariate_by$verbose <- FALSE
  }
  
  if ((is.na(univariate_by$by) | univariate_by$by == "NA")) {
    univariate_by$cor <- "un"
    uvarabccorr <- FALSE
    univariate_by$verbose <- FALSE
  }
  
  
  
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (!is.null(univariate_by$cor)) {
      if (univariate_by$cor == "un")
        uvarsigmacorr <- TRUE
      if (univariate_by$cor == "diagonal")
        uvarsigmacorr <- FALSE
    } else {
      univariate_by$cor <- "un"
      uvarsigmacorr <- TRUE
    }
    if (is.null(univariate_by$verbose))
      univariate_by$verbose <- FALSE
  }
  
  if ((is.na(univariate_by$by) | univariate_by$by == "NA")) {
    univariate_by$cor <- "un"
    uvarsigmacorr <- FALSE
    univariate_by$verbose <- FALSE
  }
  
  
  
  
  if (multivariate$mvar) {
    if (!is.null(multivariate$cor)) {
      if (multivariate$cor == "un")
        mvarccorr <- "WB"
      if (multivariate$cor == "un_s")
        mvarccorr <- "W"
      if (multivariate$cor == "diagonal")
        mvarccorr <- "none"
    } else {
      multivariate$cor <- "un"
      mvarccorr <- "WB"
    }
    if (!is.null(multivariate$rescor)) {
      multivariate$rescor <- multivariate$rescor
    } else {
      multivariate$rescor <- TRUE
    }
    if (is.null(multivariate$verbose))
      multivariate$verbose <- FALSE
  }
  
  if (!multivariate$mvar) {
    multivariate$cor <- "none"
    multivariate$rescor <- FALSE
    multivariate$verbose <- FALSE
  }
  
  
  if (!multivariate$rescor) {
    mvr_prior_rescor <- NULL
  }
  
  
  
  getcovlist <- function(x) {
    if (is.character(x))
      x <- x
    else
      x <- deparse(x)
    x <- gsub("~", "", gsub("\\s", "", x))
    x <- strsplit(x, "+", fixed = T)[[1]]
    if (length(x) == 1)
      x <- NULL
    else
      x <- x[-1]
    return(x)
  }
  
  if (is.null(getcovlist(a_formulasi)))
    a_cov_prior_beta <- NULL
  if (is.null(getcovlist(b_formulasi)))
    b_cov_prior_beta <- NULL
  if (is.null(getcovlist(c_formulasi)))
    c_cov_prior_beta <- NULL
  if (is.null(getcovlist(d_formulasi)))
    d_cov_prior_beta <- NULL
  if (is.null(getcovlist(e_formulasi)))
    e_cov_prior_beta <- NULL
  if (is.null(getcovlist(f_formulasi)))
    f_cov_prior_beta <- NULL
  if (is.null(getcovlist(g_formulasi)))
    g_cov_prior_beta <- NULL
  if (is.null(getcovlist(h_formulasi)))
    h_cov_prior_beta <- NULL
  if (is.null(getcovlist(i_formulasi)))
    i_cov_prior_beta <- NULL
  if (is.null(getcovlist(s_formulasi)))
    s_cov_prior_beta <- NULL
  
  
  if(select_model != 'sitar' & select_model != 'rcs') {
    s_prior_beta <- s_cov_prior_beta <- NULL
  }
  
  
  if (!grepl("a", fixedsi, fixed = T))
    a_prior_beta <- NULL
  if (!grepl("b", fixedsi, fixed = T))
    b_prior_beta <- NULL
  if (!grepl("c", fixedsi, fixed = T))
    c_prior_beta <- NULL
  if (!grepl("d", fixedsi, fixed = T))
    d_prior_beta <- NULL
  if (!grepl("e", fixedsi, fixed = T))
    e_prior_beta <- NULL
  if (!grepl("f", fixedsi, fixed = T))
    f_prior_beta <- NULL
  if (!grepl("g", fixedsi, fixed = T))
    g_prior_beta <- NULL
  if (!grepl("h", fixedsi, fixed = T))
    h_prior_beta <- NULL
  if (!grepl("i", fixedsi, fixed = T))
    i_prior_beta <- NULL
  if (!grepl("s", fixedsi, fixed = T))
    s_prior_beta <- NULL
  
  if (!grepl("a", randomsi, fixed = T))
    a_prior_sd <- NULL
  if (!grepl("b", randomsi, fixed = T))
    b_prior_sd <- NULL
  if (!grepl("c", randomsi, fixed = T))
    c_prior_sd <- NULL
  if (!grepl("d", randomsi, fixed = T))
    d_prior_sd <- NULL
  if (!grepl("e", randomsi, fixed = T))
    e_prior_sd <- NULL
  if (!grepl("f", randomsi, fixed = T))
    f_prior_sd <- NULL
  if (!grepl("g", randomsi, fixed = T))
    g_prior_sd <- NULL
  if (!grepl("h", randomsi, fixed = T))
    h_prior_sd <- NULL
  if (!grepl("i", randomsi, fixed = T))
    i_prior_sd <- NULL
  if (!grepl("s", randomsi, fixed = T))
    s_prior_sd <- NULL
  
  
  if (is.null(getcovlist(a_formula_grsi)))
    a_cov_prior_sd <- NULL
  if (is.null(getcovlist(b_formula_grsi)))
    b_cov_prior_sd <- NULL
  if (is.null(getcovlist(c_formula_grsi)))
    c_cov_prior_sd <- NULL
  if (is.null(getcovlist(d_formula_grsi)))
    d_cov_prior_sd <- NULL
  if (is.null(getcovlist(e_formula_grsi)))
    e_cov_prior_sd <- NULL
  if (is.null(getcovlist(f_formula_grsi)))
    f_cov_prior_sd <- NULL
  if (is.null(getcovlist(g_formula_grsi)))
    g_cov_prior_sd <- NULL
  if (is.null(getcovlist(h_formula_grsi)))
    h_cov_prior_sd <- NULL
  if (is.null(getcovlist(i_formula_grsi)))
    i_cov_prior_sd <- NULL
  if (is.null(getcovlist(s_formula_grsi)))
    s_cov_prior_sd <- NULL
  
  if (!is.null(a_cov_prior_beta))
    ancov <- length(acovcoefnames)
  else
    ancov <- NULL
  if (!is.null(b_cov_prior_beta))
    bncov <- length(bcovcoefnames)
  else
    bncov <- NULL
  if (!is.null(c_cov_prior_beta))
    cncov <- length(ccovcoefnames)
  else
    cncov <- NULL
  if (!is.null(d_cov_prior_beta))
    dncov <- length(dcovcoefnames)
  else
    dncov <- NULL
  if (!is.null(e_cov_prior_beta))
    encov <- length(ecovcoefnames)
  else
    encov <- NULL
  if (!is.null(f_cov_prior_beta))
    fncov <- length(fcovcoefnames)
  else
    fncov <- NULL
  if (!is.null(g_cov_prior_beta))
    gncov <- length(gcovcoefnames)
  else
    gncov <- NULL
  if (!is.null(h_cov_prior_beta))
    hncov <- length(hcovcoefnames)
  else
    hncov <- NULL
  if (!is.null(i_cov_prior_beta))
    incov <- length(icovcoefnames)
  else
    incov <- NULL
  if (!is.null(s_cov_prior_beta))
    sncov <- length(scovcoefnames)
  else
    sncov <- NULL
  
  
  # if (!is.null(s_cov_prior_beta))
  #   sncov <- length(scovcoefnames)
  # else
  #   sncov <- NULL
  
  if (!is.null(a_cov_prior_sd))
    ancov_gr <- length(acovcoefnames_gr)
  else
    ancov_gr <- NULL
  if (!is.null(b_cov_prior_sd))
    bncov_gr <- length(bcovcoefnames_gr)
  else
    bncov_gr <- NULL
  if (!is.null(c_cov_prior_sd))
    cncov_gr <- length(ccovcoefnames_gr)
  else
    cncov_gr <- NULL
  if (!is.null(d_cov_prior_sd))
    dncov_gr <- length(dcovcoefnames_gr)
  else
    dncov_gr <- NULL
  if (!is.null(e_cov_prior_sd))
    encov_gr <- length(ecovcoefnames_gr)
  else
    encov_gr <- NULL
  if (!is.null(f_cov_prior_sd))
    fncov_gr <- length(fcovcoefnames_gr)
  else
    fncov_gr <- NULL
  if (!is.null(g_cov_prior_sd))
    gncov_gr <- length(gcovcoefnames_gr)
  else
    gncov_gr <- NULL
  if (!is.null(h_cov_prior_sd))
    hncov_gr <- length(hcovcoefnames_gr)
  else
    hncov_gr <- NULL
  if (!is.null(i_cov_prior_sd))
    incov_gr <- length(icovcoefnames_gr)
  else
    incov_gr <- NULL
  if (!is.null(s_cov_prior_sd))
    sncov_gr <- length(scovcoefnames_gr)
  else
    sncov_gr <- NULL
  
  
  
  if (!is.null(dpar_formulasi)) {
    if (!grepl("lf\\(", dpar_formulasi) &
        !grepl("nlf\\(", dpar_formulasi)) {
      dpar_covi_mat_form <- dpar_formulasi
    } else {
      dpar_covi_mat_form <-
        gsub("\\(|)", "", strsplit(dpar_formulasi, "~")[[1]][2])
      dpar_covi_mat_form <- paste0("~", dpar_covi_mat_form)
    }
    dpar_covi_mat_form <- strsplit(dpar_covi_mat_form, ",")[[1]][1]
  }
  
  if (is.null(dpar_formulasi)) {
    dpar_covi_mat_form <- NULL
  }
  
  # rsd_prior_sigma and dpar_prior are mutually exclusive
  
  if (!is.null(dpar_formulasi)) {
    rsd_prior_sigma <- NULL
  } else {
    dpar_prior_sigma <- dpar_cov_prior_sigma <- NULL
  }
  
  
  
  if(!is.null(sigma_formulasi[1]) & sigma_formulasi != 'NULL') {
    rsd_prior_sigma <- dpar_prior_sigma <- dpar_cov_prior_sigma <- NULL
  }
  
  if(is.null(sigma_formulasi[1]) | sigma_formulasi == 'NULL') {
    sigma_prior_beta <- NULL
  }
  if (is.null(getcovlist(sigma_formulasi))) {
    sigma_cov_prior_beta <- NULL
  }
  
  
  
  if(is.null(sigma_formula_grsi)) {
    sigma_prior_sd <- NULL
  }
  if (is.null(getcovlist(sigma_formula_grsi))) {
    sigma_cov_prior_sd <- NULL
  }
  
  
  
  if (!is.null(sigma_cov_prior_beta))
    sigmancov <- length(sigmacovcoefnames)
  else
    sigmancov <- NULL
  
  if (!is.null(sigma_cov_prior_sd))
    sigmancov_gr <- length(sigmacovcoefnames_gr)
  else
    sigmancov_gr <- NULL
  
  
  
  
  if(!is.null(a_formulasi)) {
    if (grepl("~0", a_formulasi, fixed = T)) {
      a_form_0 <- TRUE
      a_cov_prior_beta <- NULL
    } else {
      a_form_0 <- FALSE
    }
  } else {
    a_form_0 <- FALSE
  }
  
  if(!is.null(b_formulasi)) {
    if (grepl("~0", b_formulasi, fixed = T)) {
      b_form_0 <- TRUE
      b_cov_prior_beta <- NULL
    } else {
      b_form_0 <- FALSE
    }
  } else {
    b_form_0 <- FALSE
  }
  
  
  if(!is.null(c_formulasi)) {
    if (grepl("~0", c_formulasi, fixed = T)) {
      c_form_0 <- TRUE
      c_cov_prior_beta <- NULL
    } else {
      c_form_0 <- FALSE
    }
  } else {
    c_form_0 <- FALSE
  }
  
  
  if(!is.null(d_formulasi)) {
    if (grepl("~0", d_formulasi, fixed = T)) {
      d_form_0 <- TRUE
      d_cov_prior_beta <- NULL
    } else {
      d_form_0 <- FALSE
    }
  } else {
    d_form_0 <- FALSE
  }
  
  
  if(!is.null(e_formulasi)) {
    if (grepl("~0", e_formulasi, fixed = T)) {
      e_form_0 <- TRUE
      e_cov_prior_beta <- NULL
    } else {
      e_form_0 <- FALSE
    }
  } else {
    e_prior_beta <- e_cov_prior_beta <- NULL
    e_prior_sd   <- e_cov_prior_sd <- NULL
    e_form_0 <- FALSE
  }
  
  
  if(!is.null(f_formulasi)) {
    if (grepl("~0", f_formulasi, fixed = T)) {
      f_form_0 <- TRUE
      f_cov_prior_beta <- NULL
    } else {
      f_form_0 <- FALSE
    }
  } else {
    f_form_0 <- FALSE
  }
  
  if(!is.null(g_formulasi)) {
    if (grepl("~0", g_formulasi, fixed = T)) {
      g_form_0 <- TRUE
      g_cov_prior_beta <- NULL
    } else {
      g_form_0 <- FALSE
    }
  } else {
    g_form_0 <- FALSE
  }
  
  
  if(!is.null(h_formulasi)) {
    if (grepl("~0", h_formulasi, fixed = T)) {
      h_form_0 <- TRUE
      h_cov_prior_beta <- NULL
    } else {
      h_form_0 <- FALSE
    }
  } else {
    h_form_0 <- FALSE
  }
  
  
  if(!is.null(i_formulasi)) {
    if (grepl("~0", i_formulasi, fixed = T)) {
      i_form_0 <- TRUE
      i_cov_prior_beta <- NULL
    } else {
      i_form_0 <- FALSE
    }
  } else {
    i_form_0 <- FALSE
  }
  
  
  if(!is.null(s_formulasi)) {
    if (grepl("~0", s_formulasi, fixed = T)) {
      s_form_0 <- TRUE
      s_cov_prior_beta <- NULL
    } else {
      s_form_0 <- FALSE
    }
  } else {
    s_form_0 <- FALSE
  }
  
  
  
  
  if(!is.null(a_formula_grsi)) {
    if (grepl("~0", a_formula_grsi, fixed = T)) {
      a_form_0_gr <- TRUE
      a_cov_prior_sd <- NULL
    } else {
      a_form_0_gr <- FALSE
    }
  } else {
    a_form_0_gr <- FALSE
  }
  
  
  if(!is.null(b_formula_grsi)) {
    if (grepl("~0", b_formula_grsi, fixed = T)) {
      b_form_0_gr <- TRUE
      b_cov_prior_sd <- NULL
    } else {
      b_form_0_gr <- FALSE
    }
  } else {
    b_form_0_gr <- FALSE
  }
  
  
  if(!is.null(c_formula_grsi)) {
    if (grepl("~0", c_formula_grsi, fixed = T)) {
      c_form_0_gr <- TRUE
      c_cov_prior_sd <- NULL
    } else {
      c_form_0_gr <- FALSE
    }
  } else {
    c_form_0_gr <- FALSE
  }
  
  
  
  if(!is.null(d_formula_grsi)) {
    if (grepl("~0", d_formula_grsi, fixed = T)) {
      d_form_0_gr <- TRUE
      d_cov_prior_sd <- NULL
    } else {
      d_form_0_gr <- FALSE
    }
  } else {
    d_form_0_gr <- FALSE
  }
  
  
  if(!is.null(e_formula_grsi)) {
    if (grepl("~0", e_formula_grsi, fixed = T)) {
      e_form_0_gr <- TRUE
      e_cov_prior_sd <- NULL
    } else {
      e_form_0_gr <- FALSE
    }
  } else {
    e_form_0_gr <- FALSE
  }
  
  
  if(!is.null(f_formula_grsi)) {
    if (grepl("~0", f_formula_grsi, fixed = T)) {
      f_form_0_gr <- TRUE
      f_cov_prior_sd <- NULL
    } else {
      f_form_0_gr <- FALSE
    }
  } else {
    f_form_0_gr <- FALSE
  }
  
  
  
  if(!is.null(g_formula_grsi)) {
    if (grepl("~0", g_formula_grsi, fixed = T)) {
      g_form_0_gr <- TRUE
      g_cov_prior_sd <- NULL
    } else {
      g_form_0_gr <- FALSE
    }
  } else {
    g_form_0_gr <- FALSE
  }
  
  
  if(!is.null(h_formula_grsi)) {
    if (grepl("~0", h_formula_grsi, fixed = T)) {
      h_form_0_gr <- TRUE
      h_cov_prior_sd <- NULL
    } else {
      h_form_0_gr <- FALSE
    }
  } else {
    h_form_0_gr <- FALSE
  }
  
  
  if(!is.null(i_formula_grsi)) {
    if (grepl("~0", i_formula_grsi, fixed = T)) {
      i_form_0_gr <- TRUE
      i_cov_prior_sd <- NULL
    } else {
      i_form_0_gr <- FALSE
    }
  } else {
    i_form_0_gr <- FALSE
  }
  
  
  if(!is.null(s_formula_grsi)) {
    if (grepl("~0", s_formula_grsi, fixed = T)) {
      s_form_0_gr <- TRUE
      s_cov_prior_sd <- NULL
    } else {
      s_form_0_gr <- FALSE
    }
  } else {
    s_form_0_gr <- FALSE
  }
  
  
  
  
  
  
  
  
  
  
  
  if (!is.null(dpar_formulasi)) {
    if (grepl("^~1$", dpar_covi_mat_form, fixed = F)) {
      dpar_intercept_only <- TRUE
    } else {
      dpar_intercept_only <- FALSE
    }
    if (!is.null(dpar_covi_mat_form) &
        grepl("~0", dpar_covi_mat_form, fixed = T)) {
      dpar_form_0 <- TRUE
      dpar_cov_prior_sigma <- NULL
    } else {
      dpar_form_0 <- FALSE
      if (grepl("^~1$", dpar_covi_mat_form, fixed = F)) {
        dpar_cov_prior_sigma <- NULL
      }
    }
  }
  
  if (is.null(dpar_formulasi)) {
    dpar_form_0 <- FALSE
    dpar_cov_prior_sigma <- NULL
  }
  
  
  
  
  
  
  if (grepl("~0", sigma_formulasi, fixed = T)) {
    sigma_form_0 <- TRUE
    sigma_cov_prior_beta <- NULL
  } else {
    sigma_form_0 <- FALSE
  }
  
  if(is.null(sigma_formulasi)) sigma_form_0 <- FALSE
  
  
  if(!is.null(sigma_formula_grsi)) {
    if (grepl("~0", sigma_formula_grsi, fixed = T)) {
      sigma_form_0_gr <- TRUE
      sigma_cov_prior_sd <- NULL
    } else {
      sigma_form_0_gr <- FALSE
    }
  }
  
  if(is.null(sigma_formula_grsi)) sigma_form_0_gr <- FALSE
  
  
  
  
  
  
  
  
  
  
  
  if ((is.na(univariate_by$by) |
       univariate_by$by == "NA") & !mvar & !abccorr) {
    gr_prior_cor <- NULL
  }
  
  if (is.null(randomsi[[1]]))
    gr_prior_cor <- NULL
  if (nabcrei == 1)
    gr_prior_cor <- NULL
  
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (uvarabccorr) {
      gr_prior_cor <- gr_prior_cor
    } else {
      gr_prior_cor <- NULL
    }
  }
  
  if (mvar && mvarccorr == "none") {
    gr_prior_cor <- NULL
  }
  
  
  
  
  
  
  
  
  if ((is.na(univariate_by$by) |
       univariate_by$by == "NA") & !mvar & !sigmacorr) {
    sigma_prior_cor <- NULL
  }
  
  if(is.null(sigmancov_gr)) {
    sigma_prior_cor <- NULL
  }
  if(!is.null(sigmancov_gr)) {
    #if (length(sigmancov_gr) == 1) sigma_prior_cor <- NULL
  }
  
  
  
  # Note that currently brms does not allow setting separate cor prior for sigma
  # So, either set sigma_prior_cor <- NULL for all or else set group = '
  
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (uvarsigmacorr) {
      sigma_prior_cor <- sigma_prior_cor
    } else {
      sigma_prior_cor <- NULL
    }
  }
  
  if (mvar && mvarccorr == "none") {
    sigma_prior_cor <- NULL
  }
  
  
  # evaluate prior arguments
  eval_prior_args <- function(x, ...) {
    x_org <- x
    if (grepl("_beta", x))
      class <- 'b'
    if (grepl("_sd", x))
      class <- 'sd'
    if (grepl("rsd_", x) & grepl("_sigma", x))
      class <- 'sigma'
    if (grepl("_cor$", x))
      class <- 'cor'
    if (grepl("dpar", x) &
        grepl("_sigma", x))
      class <- '' # no class sigma if sigma ~ .
    
    if (grepl("_rescor$", x))
      class <- 'rescor'
    
    nlpar <- ""
    if (grepl("a_", x))
      nlpar <- 'a'
    if (grepl("b_", x))
      nlpar <- 'b'
    if (grepl("c_", x))
      nlpar <- 'c'
    if (grepl("d_", x))
      nlpar <- 'd'
    if (grepl("e_", x))
      nlpar <- 'e'
    if (grepl("f_", x))
      nlpar <- 'f'
    if (grepl("g_", x))
      nlpar <- 'g'
    if (grepl("h_", x))
      nlpar <- 'h'
    if (grepl("i_", x))
      nlpar <- 'i'
    if (grepl("s_", x))
      nlpar <- 's'
    dpar <- ""
    if (grepl("dpar", x) &
        grepl("_sigma", x) & !grepl("rsd_", x))
      dpar <- "dpar"
    
    sigma_dpar <- ""
    if (grepl("sigma", x) &
        grepl("sigma_", x) & !grepl("rsd_", x) & !grepl("_sigma", x) & 
        !grepl("dpar", x)) {
      sigma_dpar <- "sigma"
      nlpar <- ""
      cov_nlpar <- ""
      # class <- 'b'
    }
    
    
    
    cov_dpar <- ""
    if (grepl("dpar_cov", x))
      dpar <- paste0(dpar, "_cov")
    
    dpar_cov <- dpar
    
    
    cov_sigma_dpar <- ""
    if (grepl("sigma_cov", x)) {
      sigma_dpar <- paste0('sigma', "")
      cov_sigma_dpar <- paste0('sigma', "_cov")
      nlpar <- ""
      cov_nlpar <- ""
      # class <- 'b'
    }
    
    
    
    cov_nlpar <- ""
    if (grepl("a_cov", x))
      cov_nlpar <- 'a'
    if (grepl("b_cov", x))
      cov_nlpar <- 'b'
    if (grepl("c_cov", x))
      cov_nlpar <- 'c'
    if (grepl("d_cov", x))
      cov_nlpar <- 'd'
    if (grepl("e_cov", x))
      cov_nlpar <- 'e'
    if (grepl("f_cov", x))
      cov_nlpar <- 'f'
    if (grepl("g_cov", x))
      cov_nlpar <- 'g'
    if (grepl("h_cov", x))
      cov_nlpar <- 'h'
    if (grepl("i_cov", x))
      cov_nlpar <- 'i'
    if (grepl("s_cov", x))
      cov_nlpar <- 's'
    
    if (!is.null(dpar_prior_sigma) |
        !is.null(dpar_cov_prior_sigma)) {
      dparncov <- length(dparcovcoefnames) - 1
    } else {
      dparncov <- NULL
    }
    
    
    
    if(sigma_dpar == 'sigma' | cov_sigma_dpar == 'sigma_cov') {
      group <- sigma_arg_groupvar
    } else {
      group <- group_arg_groupvar
    }
    
    
    get_acorclass <- function(autocor_formi2) {
      if(grepl("arma\\(", autocor_formi2, fixed = F)) {
        acorclass <- 'arma'
      } else if(grepl("ar\\(", autocor_formi2, fixed = F)) {
        acorclass <- 'ar'
      } else if(grepl("ma\\(", autocor_formi2, fixed = F)) {
        acorclass <- 'ma'
      } else if(grepl("cosy\\(", autocor_formi2, fixed = F)) {
        acorclass <- 'cosy'
      } else if(grepl("car\\(", autocor_formi2, fixed = F)) {
        acorclass <- 'car'
      } else if(grepl("lagsar\\(", autocor_formi2, fixed = F)) {
        acorclass <- 'lagsar'
      } else if(grepl("errorsar\\(", autocor_formi2, fixed = F)) {
        acorclass <- 'errorsar'
      } else if(grepl("unstr\\(", autocor_formi2, fixed = F)) {
        acorclass <- 'Lcortime'
      }
      acorclass
    }
    
    
    setautocorr <- FALSE
    if (grepl("autocor_", x) & grepl("_acor", x)) {
      setautocorr <- TRUE
      # this does not work anymore, so using the function next
      # acorclass <- gsub("~", "", autocor_formi, fixed = T)
      # acorclass <- strsplit(acorclass, "\\(")[[1]][1]
      acorclass <- get_acorclass(autocor_formi)
      # stop()
      if (acorclass == "arma")
        class <- 'arma'
      if (acorclass == "ar")
        class <- 'ar'
      if (acorclass == "ma")
        class <- 'ma'
      if (acorclass == "cosy")
        class <- 'cosy'
      if (acorclass == "car")
        class <- 'car'
      if (acorclass == "lagsar")
        class <- 'lagsar'
      if (acorclass == "errorsar")
        class <- 'errorsar'
      if (acorclass == "Lcortime")
        class <- 'Lcortime'
    }
    
    
    if(setautocorr) {
      if(class == 'Lcortime') {
        acorclassname <- 'unstr' 
      } else {
        acorclassname <- acorclass
      }
      
      allowed_acor_classes <- c('ar', 'ma', 'arma', 'unstr')
      
      if(!acorclassname %in% allowed_acor_classes) {
        stop("Allowed autocorrelation classes are ", 
             paste(allowed_acor_classes, collapse = ", ")) 
      }
    }
    
    
    
    if (class == "b" & nlpar == 'a') {
      if (a_form_0) {
        nrep_of_parms <- length(acovcoefnames)
      } else {
        if (!grepl("a_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("a_cov", x)) {
          nrep_of_parms <- length(acovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'b') {
      if (b_form_0) {
        nrep_of_parms <- length(bcovcoefnames)
      } else {
        if (!grepl("b_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("b_cov", x)) {
          nrep_of_parms <- length(bcovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'c') {
      if (c_form_0) {
        nrep_of_parms <- length(ccovcoefnames)
      } else {
        if (!grepl("c_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("c_cov", x)) {
          nrep_of_parms <- length(ccovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'd') {
      if (d_form_0) {
        nrep_of_parms <- length(dcovcoefnames)
      } else {
        if (!grepl("d_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("d_cov", x)) {
          nrep_of_parms <- length(dcovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'e') {
      if (e_form_0) {
        nrep_of_parms <- length(ecovcoefnames)
      } else {
        if (!grepl("e_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("e_cov", x)) {
          nrep_of_parms <- length(ecovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'f') {
      if (f_form_0) {
        nrep_of_parms <- length(fcovcoefnames)
      } else {
        if (!grepl("f_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("f_cov", x)) {
          nrep_of_parms <- length(fcovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'g') {
      if (g_form_0) {
        nrep_of_parms <- length(gcovcoefnames)
      } else {
        if (!grepl("g_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("g_cov", x)) {
          nrep_of_parms <- length(gcovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'h') {
      if (h_form_0) {
        nrep_of_parms <- length(hcovcoefnames)
      } else {
        if (!grepl("h_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("h_cov", x)) {
          nrep_of_parms <- length(hcovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'i') {
      if (i_form_0) {
        nrep_of_parms <- length(icovcoefnames)
      } else {
        if (!grepl("i_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("i_cov", x)) {
          nrep_of_parms <- length(icovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 's') {
      if (s_form_0) {
        nrep_of_parms <- df * length(scovcoefnames)
      } else {
        if (!grepl("s_cov", x)) {
          nrep_of_parms <- df
        } else if (grepl("s_cov", x)) {
          nrep_of_parms <- df * (length(scovcoefnames) - 1)
        }
      }
    } else if (class == "sd" & nlpar == 'a') {
      if (a_form_0_gr) {
        nrep_of_parms <- length(acovcoefnames_gr)
      } else {
        if (!grepl("a_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("a_cov", x)) {
          nrep_of_parms <- length(acovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'b') {
      if (b_form_0_gr) {
        nrep_of_parms <- length(bcovcoefnames_gr)
      } else {
        if (!grepl("b_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("b_cov", x)) {
          nrep_of_parms <- length(bcovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'c') {
      if (c_form_0_gr) {
        nrep_of_parms <- length(ccovcoefnames_gr)
      } else {
        if (!grepl("c_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("c_cov", x)) {
          nrep_of_parms <- length(ccovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'd') {
      if (d_form_0_gr) {
        nrep_of_parms <- length(dcovcoefnames_gr)
      } else {
        if (!grepl("d_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("d_cov", x)) {
          nrep_of_parms <- length(dcovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'e') {
      if (e_form_0_gr) {
        nrep_of_parms <- length(ecovcoefnames_gr)
      } else {
        if (!grepl("e_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("e_cov", x)) {
          nrep_of_parms <- length(ecovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'f') {
      if (f_form_0_gr) {
        nrep_of_parms <- length(fcovcoefnames_gr)
      } else {
        if (!grepl("f_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("f_cov", x)) {
          nrep_of_parms <- length(fcovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'g') {
      if (g_form_0_gr) {
        nrep_of_parms <- length(gcovcoefnames_gr)
      } else {
        if (!grepl("g_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("g_cov", x)) {
          nrep_of_parms <- length(gcovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'h') {
      if (h_form_0_gr) {
        nrep_of_parms <- length(hcovcoefnames_gr)
      } else {
        if (!grepl("h_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("h_cov", x)) {
          nrep_of_parms <- length(hcovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'i') {
      if (i_form_0_gr) {
        nrep_of_parms <- length(icovcoefnames_gr)
      } else {
        if (!grepl("i_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("i_cov", x)) {
          nrep_of_parms <- length(icovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 's') {
      if (s_form_0_gr) {
        nrep_of_parms <- df * length(scovcoefnames_gr)
      } else {
        if (!grepl("s_cov", x)) {
          nrep_of_parms <- df
        } else if (grepl("s_cov", x)) {
          nrep_of_parms <- df * (length(scovcoefnames_gr) - 1)
        }
      }
    } else if (class == "sigma" &
               (class != 'b' |
                class != 'cor') & is.null(dparncov)) {
      nrep_of_parms <- 1
    } else if (class == "" &
               (class != 'b' |
                class != 'cor') & !is.null(dparncov)) {
      if (dpar_form_0) {
        nrep_of_parms <- length(dparcovcoefnames)
      } else {
        if (!grepl("dpar_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("dpar_cov", x)) {
          nrep_of_parms <- length(dparcovcoefnames) - 1
        }
      }
    } else if (setautocorr &
               class == "" &
               (class != 'b' |
                class != 'cor') & is.null(dparncov)) {
      nrep_of_parms <- 1
    } else if (class == "cor" &
               (class != 'b' |
                class != 'sigma') & is.null(dparncov)) {
      nrep_of_parms <- 1
    } else {
      nrep_of_parms <- 1
    }
    
    
    if(setautocorr) {
      tempzxxx <- autocor_formi 
      tempzxxx <- gsub("[[:space:]]", "", tempzxxx)
      if(grepl("p=", tempzxxx, fixed = F)) {
        acor_dim_p <- sub(",.*", "", sub(".*p=", "", tempzxxx) )  
        acor_dim_p <- gsub(")" , "", acor_dim_p, fixed = T)
      } else if(!grepl("p=", tempzxxx, fixed = F)) {
        acor_dim_p <- '1'
      }
      if(grepl("q=", tempzxxx, fixed = F)) {
        acor_dim_q <- sub(",.*", "", sub(".*q=", "", tempzxxx) )  
        acor_dim_q <- gsub(")" , "", acor_dim_q, fixed = T) 
      } else if(!grepl("q=", tempzxxx, fixed = F)) {
        acor_dim_q <- '1'
      }
      
      if(acorclass == 'ar') {
        if(!grepl("p=", tempzxxx, fixed = F)) 
          stop("Please specify arguments by names e.g., p=1")
      } else if(acorclass == 'ma') {
        if(!grepl("q=", tempzxxx, fixed = F))
          stop("Please specify arguments by names e.g., q=1")
      } else if(acorclass == 'arma') {
        if(!grepl("p=", tempzxxx, fixed = F)) 
          stop("Please specify arguments by names e.g., p=1")
        if(!grepl("q=", tempzxxx, fixed = F))
          stop("Please specify arguments by names e.g., q=1")
      }
      
      nrep_of_parms_p <- eval(parse(text = acor_dim_p))
      nrep_of_parms_q <- eval(parse(text = acor_dim_q))
    }
    
    if(!setautocorr) {
      nrep_of_parms_p <- nrep_of_parms_q <- 1
    }
    
    
    
    if (class == "b" & sigma_dpar == 'sigma' & 
        !grepl("rsd_", x) & !grepl("_sigma", x) & !grepl("dpar", x)) {
      if (sigma_form_0) {
        nrep_of_parms <- length(sigmacovcoefnames)
      } else {
        if (!grepl("sigma_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("sigma_cov", x)) {
          nrep_of_parms <- length(sigmacovcoefnames) - 1
        }
      }
    } 
    
    
    
    if (class == "b" & sigma_dpar == 'sigma') {
      if (sigma_form_0) {
        nrep_of_parms <- length(sigmacovcoefnames)
      } else {
        if (!grepl("sigma_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("sigma_cov", x)) {
          nrep_of_parms <- length(sigmacovcoefnames) - 1
        }
      }  
    }
    
    
    if (class == "sd" & sigma_dpar == 'sigma' & 
        !grepl("rsd_", x) & !grepl("_sigma", x) & !grepl("dpar", x)) {
      if (sigma_form_0_gr) {
        nrep_of_parms <- length(sigmacovcoefnames_gr)
      } else {
        if (!grepl("sigma_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("sigma_cov", x)) {
          nrep_of_parms <- length(sigmacovcoefnames_gr) - 1
        }
      }
    }
    
    
  
    
    get_priors_parms <- function(x,
                                 prior_data,
                                 prior_data_internal,
                                 resp,
                                 nlpar,
                                 dpar,
                                 sigma_dpar,
                                 class,
                                 acorclass,
                                 get_priors_parms_args) {
      if (!is.null(get_priors_parms_args)) {
        eout <- list2env(get_priors_parms_args)
        for (eoutii in names(eout)) {
          assign(eoutii, eout[[eoutii]])
        }
      }
      
      x <- gsub("\"", "", gsub("\\s", "", x))
      
      if (resp == "") {
        resp_ <- ""
      } else if (resp != "") {
        resp_ <- paste0("_", resp)
      }
      
      prior_argument <- x
      
      zz <- prior_str_arg <- eval(parse(text = x))
      zz <- strsplit(zz, "\\(")[[1]]
      dist <- zz[1]
      
      
      #################
      
      list_names <-
        c(
          'prior_str_arg',
          'dist',
          'resp',
          'resp_',
          'nlpar',
          'dpar',
          'class',
          'cov_nlpar',
          'cov_dpar',
          'ancov',
          'bncov',
          'cncov',
          'dncov',
          'encov',
          'fncov',
          'gncov',
          'hncov',
          'incov',
          'sncov',
          'ancov_gr',
          'bncov_gr',
          'cncov_gr',
          'dncov_gr',
          'encov_gr',
          'fncov_gr',
          'gncov_gr',
          'hncov_gr',
          'incov_gr',
          'sncov_gr',
          'dparncov',
          'nabci',
          'nabcrei',
          'fixedsi',
          'randomsi',
          'df',
          'setautocorr',
          'nrep_of_parms',
          'nrep_of_parms_p',
          'nrep_of_parms_q',
          'N_J_all',
          'ii',
          'nys',
          'a_form_0',
          'b_form_0',
          'c_form_0',
          'd_form_0',
          'e_form_0',
          'f_form_0',
          'g_form_0',
          'h_form_0',
          'i_form_0',
          's_form_0',
          'a_form_0_gr',
          'b_form_0_gr',
          'c_form_0_gr',
          'd_form_0_gr',
          'e_form_0_gr',
          'f_form_0_gr',
          'g_form_0_gr',
          'h_form_0_gr',
          'i_form_0_gr',
          's_form_0_gr',
          'sigma_form_0',
          'sigma_form_0_gr',
          "sigma_dpar",
          "cov_sigma_dpar",
          'sigmancov',
          'sigmancov_gr',
          'sigma_dpar',
          'sigma_group_arg',
          
          'group_arg_groupvar',
          'group',
          
          'dpar_form_0',
          'dpar_covi_mat_form',
          'dpar_formulasi',
          'univariate_by',
          'multivariate',
          'group_arg',
          'setautocorr',
          'acorclass',
          'initsi',
          'normalize',
          'seed',
          'cortimeNlags_var',
          'cortimeNlags',
          'verbose'
        )
      
      
      prior_internal_args <- mget(list_names)
      
      # set to NULL as appropriate to match priors
      for (ip in names(init_arguments)) {
        init_ip <- ip
        if (grepl("_init_", init_ip) &
            !grepl("r_init_z", init_ip)) {
          prior_ip <- gsub("_init_", "_prior_", init_ip)
          # if(!is.null(eval(parse(text = prior_ip)))) {
          xx   <- deparse(prior_ip)
          if (is.null(ept(prior_ip)))
            init_arguments[[ip]] <- 'NULL'
          # }
        }
      }
      
      if (ept(nabcrei) == 0)
        init_arguments[['r_init_z']] <- 'NULL'
      
      out_p_str <- prepare_priors(
        prior_argument,
        prior_data,
        prior_data_internal,
        prior_internal_args,
        init_arguments,
        init_data,
        init_data_internal,
        init_args_internal
      )
      
      
      stanvars_data_in <- out_p_str$stanvars_data
      prior_str_arg    <- out_p_str$prior_str_arg
      lowerbound       <- out_p_str$lowerbound
      upperbound       <- out_p_str$upperbound
      initial_in       <- out_p_str$initial_out
      
      return(
        list(
          dist = dist,
          lowerbound = lowerbound,
          upperbound = upperbound,
          define_ = prior_str_arg,
          stanvars_data_in = stanvars_data_in,
          initial_in = initial_in
        )
      )
    }
    
    
    if (resp != "") {
      for (i in names(prior_data_internal)) {
        names(prior_data_internal)[names(prior_data_internal) == i] <-
          paste0(i, "_", resp)
      }
    }
    
    x_name <- deparse(substitute(x_org))
    
    
    get_priors_parms_args <- list(
      df = df,
      nys = nys,
      univariate_by = univariate_by,
      multivariate = multivariate,
      group_arg = group_arg,
      fixedsi = fixedsi,
      randomsi = randomsi,
      nabci = nabci,
      nabcrei = nabcrei,
      ii = ii,
      N_J_all = N_J_all,
      cov_nlpar = cov_nlpar,
      cov_dpar = cov_dpar,
      nrep_of_parms = nrep_of_parms,
      nrep_of_parms_p = nrep_of_parms_p,
      nrep_of_parms_q = nrep_of_parms_q,
      ancov = ancov,
      bncov = bncov,
      cncov = cncov,
      dncov = dncov,
      encov = encov,
      fncov = fncov,
      gncov = gncov,
      hncov = hncov,
      incov = incov,
      sncov = sncov,
      ancov_gr = ancov_gr,
      bncov_gr = bncov_gr,
      cncov_gr = cncov_gr,
      dncov_gr = dncov_gr,
      encov_gr = encov_gr,
      fncov_gr = fncov_gr,
      gncov_gr = gncov_gr,
      hncov_gr = hncov_gr,
      incov_gr = incov_gr,
      sncov_gr = sncov_gr,
      dparncov = dparncov,
      a_form_0 = a_form_0,
      b_form_0 = b_form_0,
      c_form_0 = c_form_0,
      d_form_0 = d_form_0,
      e_form_0 = e_form_0,
      f_form_0 = f_form_0,
      g_form_0 = g_form_0,
      h_form_0 = h_form_0,
      i_form_0 = i_form_0,
      s_form_0 = s_form_0,
      a_form_0_gr = a_form_0_gr,
      b_form_0_gr = b_form_0_gr,
      c_form_0_gr = c_form_0_gr,
      d_form_0_gr = d_form_0_gr,
      e_form_0_gr = e_form_0_gr,
      f_form_0_gr = f_form_0_gr,
      g_form_0_gr = g_form_0_gr,
      h_form_0_gr = h_form_0_gr,
      i_form_0_gr = i_form_0_gr,
      s_form_0_gr = s_form_0_gr,
      sigma_form_0 = sigma_form_0,
      sigma_form_0_gr = sigma_form_0_gr,
      cov_sigma_dpar = cov_sigma_dpar,
      sigmancov = sigmancov,
      sigmancov_gr = sigmancov_gr,
      sigma_group_arg = sigma_group_arg,
      group_arg_groupvar = group_arg_groupvar,
      group = group,
      
      dpar_form_0 = dpar_form_0,
      dpar_covi_mat_form = dpar_covi_mat_form,
      dpar_formulasi = dpar_formulasi,
      setautocorr = setautocorr,
      initsi = initsi,
      init_arguments = init_arguments,
      init_data = init_data,
      init_data_internal = init_data_internal,
      init_args_internal = init_args_internal,
      normalize = normalize,
      seed = seed,
      cortimeNlags_var = cortimeNlags_var,
      cortimeNlags = cortimeNlags,
      verbose = verbose
    )
    
    
    if (setautocorr) {
      if (acorclass == 'arma')
        acorclassclasses <- c("ar", "ma")
      if (acorclass == 'ar')
        acorclassclasses <- c("ar")
      if (acorclass == 'ma')
        acorclassclasses <- c("ma")
      
      
      if (acorclass == 'Lcortime')
        acorclassclasses <- c("Lcortime")
      
      priors_arma_c_define <- list()
      for (acorclassi in acorclassclasses) {
        priors_parms <- get_priors_parms(
          x_name,
          prior_data = prior_data,
          prior_data_internal = prior_data_internal,
          resp = resp,
          nlpar = nlpar,
          dpar = dpar,
          sigma_dpar = sigma_dpar,
          class = acorclassi,
          acorclass = acorclass,
          get_priors_parms_args = get_priors_parms_args
        )
        
        priors_arma_c_define[[acorclassi]] <- priors_parms
      }
    } else {
      priors_parms <- get_priors_parms(
        x_name,
        prior_data = prior_data,
        prior_data_internal = prior_data_internal,
        resp = resp,
        nlpar = nlpar,
        dpar = dpar,
        sigma_dpar = sigma_dpar,
        class = class,
        get_priors_parms_args = get_priors_parms_args
      )
    }
    
    
    define_ <- priors_parms$define_
    
    dist <- priors_parms$dist
    
    lowerbound <- priors_parms$lowerbound
    upperbound <- priors_parms$upperbound
    
    stanvars_data_in <- priors_parms$stanvars_data_in
    initial_in       <- priors_parms$initial_in
    
    
    if (class == 'b') {
      # need to remove lb and ub if specifying coef, otherwise
      # Error: Argument 'coef' may not be specified when using boundaries.
      # mnf <- paste0(nlpar, "_form_0")
      if(nlpar != '') {
        mnf <- paste0(nlpar, "_form_0")
        mnc <- paste0("cov_nlpar")
      }
      if(sigma_dpar == 'sigma' | cov_sigma_dpar == 'sigma_cov') {
        mnf <- paste0('sigma', "_form_0")
        mnc <- paste0("cov_sigma_dpar")
      }
      
      
      
      if (all(is.na(lowerbound)) |
          all(is.na(upperbound))) {
        if (nlpar == 'a')
          coef <- acovcoefnames
        if (nlpar == 'b')
          coef <- bcovcoefnames
        if (nlpar == 'c')
          coef <- ccovcoefnames
        if (nlpar == 'd')
          coef <- dcovcoefnames
        if (nlpar == 'e')
          coef <- ecovcoefnames
        if (nlpar == 'f')
          coef <- fcovcoefnames
        if (nlpar == 'g')
          coef <- gcovcoefnames
        if (nlpar == 'h')
          coef <- hcovcoefnames
        if (nlpar == 'i')
          coef <- icovcoefnames
        if (nlpar == 's') {
          coef <- scovcoefnames
          if (!s_form_0 & !is.null(sncov))
            coef <- coef[1]
        }
        
        
        # brms does not allow Intercept as coef name for dpar sigma with ~1
        # But this only when covaritae missing 
        if (sigma_dpar == 'sigma') {
          dpar <- sigma_dpar
          if(ept(mnf)) {
            coef <- sigmacovcoefnames
          }
          if(!ept(mnf)) {
            if (nlpar == '' & sigma_dpar != '' & 
                length(sigmacovcoefnames) == 1 &
                sigmacovcoefnames[1] == "Intercept" ) {
              coef <- ""
              class <- sigmacovcoefnames
            }
            if (nlpar == '' & sigma_dpar != '' & grepl("+", 
                                                       sigma_formulasi, 
                                                       fixed = T)
            ) {
              coef <- ""
              class <- 'Intercept'
            }
          }

          
          
        } else if (!all(is.na(lowerbound)) & !all(is.na(upperbound))) {
          if (nlpar == 'a')
            coef <- rep("", length(acovcoefnames))
          if (nlpar == 'b')
            coef <- rep("", length(bcovcoefnames))
          if (nlpar == 'c')
            coef <- rep("", length(ccovcoefnames))
          if (nlpar == 'd')
            coef <- rep("", length(dcovcoefnames))
          if (nlpar == 'e')
            coef <- rep("", length(ecovcoefnames))
          if (nlpar == 'f')
            coef <- rep("", length(fcovcoefnames))
          if (nlpar == 'g')
            coef <- rep("", length(gcovcoefnames))
          if (nlpar == 'h')
            coef <- rep("", length(hcovcoefnames))
          if (nlpar == 'i')
            coef <- rep("", length(icovcoefnames))
          if (nlpar == 's') {
            coef <- rep("", length(scovcoefnames))
            if (!s_form_0 & !is.null(sncov))
              coef <- coef[1]
          }
        }
      }
      
      
      
      
      # nlpar a b c d e f - betas also sigma
      if (ept(mnf) & cov_nlpar == "" & cov_sigma_dpar == "") {
        if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
          define_ <- unique(define_)
          lowerbound <- unique(lowerbound)
          upperbound <- unique(upperbound)
          setcoef <- ""
        } else {
          setcoef <- coef
        }
        
        
        
        priors_ <-
          brms::prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = setcoef,
            resp = resp,
            dpar = dpar,
            lb = lowerbound,
            ub = upperbound
          )
      }
      
      if (!ept(mnf) & cov_nlpar == "" & cov_sigma_dpar == "") {
        if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
          define_ <- unique(define_)
          lowerbound <- unique(lowerbound)
          upperbound <- unique(upperbound)
          setcoef <- ""
        } else {
          if (ept(mnc) == "") {
            setcoef <- coef[1]
          } else {
            setcoef <- coef
          }
        }
        priors_ <-
          brms::prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = setcoef,
            resp = resp,
            dpar = dpar,
            lb = lowerbound,
            ub = upperbound
          )
      }
      
      
      
      
      # nlpar s - betas
      if (nlpar == 's' & cov_nlpar == "") {
        nlpar <- paste0(nlpar, 1:df)
        if (grepl("~1", s_formulasi, fixed = T)) {
          if (all(coef == "")) {
            priors_ <- brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
          } else {
            priors_ <-   brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          }
        }
        if (grepl("~0", s_formulasi, fixed = T)) {
          nlpar <- rep(nlpar ,
                       times = 1,
                       each = length(scovcoefnames))
          if (all(coef == "")) {
            priors_ <- brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
          } else {
            priors_ <-   brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          }
        }
      }
      
      
      
      # nlpar cov a - betas
      if (!a_form_0) {
        if (class == 'b' & grepl("a_cov", x) & !is.null(a_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- acovcoefnames
          } else {
            coef <- acovcoefnames[2:length(acovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          
        }
      }
      
      # nlpar cov b - betas
      if (!b_form_0) {
        if (class == 'b' & grepl("b_cov", x) & !is.null(b_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- bcovcoefnames
          } else {
            coef <- bcovcoefnames[2:length(bcovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      # nlpar cov c - betas
      if (!c_form_0) {
        if (class == 'b' & grepl("c_cov", x) & !is.null(c_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- ccovcoefnames
          } else {
            coef <- ccovcoefnames[2:length(ccovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov d - betas
      if (!d_form_0) {
        if (class == 'b' & grepl("d_cov", x) & !is.null(d_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- dcovcoefnames
          } else {
            coef <- dcovcoefnames[2:length(dcovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov e - betas
      if (!e_form_0) {
        if (class == 'b' & grepl("e_cov", x) & !is.null(e_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- ecovcoefnames
          } else {
            coef <- ecovcoefnames[2:length(ecovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      
      # nlpar cov f - betas
      if (!f_form_0) {
        if (class == 'b' & grepl("f_cov", x) & !is.null(f_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- fcovcoefnames
          } else {
            coef <- fcovcoefnames[2:length(fcovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov g - betas
      if (!g_form_0) {
        if (class == 'b' & grepl("g_cov", x) & !is.null(g_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- gcovcoefnames
          } else {
            coef <- gcovcoefnames[2:length(gcovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      
      # nlpar cov h - betas
      if (!h_form_0) {
        if (class == 'b' & grepl("h_cov", x) & !is.null(h_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- hcovcoefnames
          } else {
            coef <- hcovcoefnames[2:length(hcovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov i - betas
      if (!i_form_0) {
        if (class == 'b' & grepl("i_cov", x) & !is.null(i_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- icovcoefnames
          } else {
            coef <- icovcoefnames[2:length(icovcoefnames)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      # sigma cov - betas
      if (!grepl("~0", sigma_formulasi, fixed = T) ) {
        class_org <- class
        if (grepl("sigma_cov", x) & !is.null(sigma_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- sigmacovcoefnames
          } else {
            coef <- sigmacovcoefnames[2:length(sigmacovcoefnames)]
            class <- c( rep('b', length(sigmacovcoefnames[-1])))
          }
          
          
          dpar <- sigma_dpar
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          
        }
        class <- class_org
      }
      
      
      
      # nlpar cov s - betas
      if (!s_form_0) {
        if (class == 'b' & grepl("s_cov", x) & !is.null(s_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- scovcoefnames
          } else {
            coef <- scovcoefnames[2:length(scovcoefnames)]
          }
          nlpar <- paste0(nlpar, 1:df)
          nlpar <- rep(nlpar, times = length(coef), each = 1)
          coef <- rep(coef , times = 1, each = df)
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          
        }
      }
      
      
    } # if(class == 'b')
    
    
    
    
    if (class == 'sd') {
      if(nlpar != '') {
        mnf <- paste0(nlpar, "_form_0_gr")
        mnc <- paste0("cov_nlpar")
      }
      if(sigma_dpar == 'sigma' | cov_sigma_dpar == 'sigma_cov') {
        mnf <- paste0('sigma', "_form_0_gr")
        mnc <- paste0("cov_sigma_dpar")
      }
      
      
      
      if (nlpar == 'a')
        coef <- acovcoefnames_gr
      if (nlpar == 'b')
        coef <- bcovcoefnames_gr
      if (nlpar == 'c')
        coef <- ccovcoefnames_gr
      if (nlpar == 'd')
        coef <- dcovcoefnames_gr
      if (nlpar == 'e')
        coef <- ecovcoefnames_gr
      if (nlpar == 'f')
        coef <- fcovcoefnames_gr
      if (nlpar == 'g')
        coef <- gcovcoefnames_gr
      if (nlpar == 'h')
        coef <- hcovcoefnames_gr
      if (nlpar == 'i')
        coef <- icovcoefnames_gr
      
      if (nlpar == 's') {
        coef <- rep("", length(scovcoefnames_gr))
        if (!s_form_0_gr & !is.null(sncov_gr))
          coef <- coef[1]
      }
      
      if (sigma_dpar == 'sigma') {
        dpar <- sigma_dpar
        coef <- sigmacovcoefnames_gr
      }
      
      
      # nlpar a b c d e f - sd
      
      if (ept(mnf) & cov_nlpar == "") {
        if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
          define_ <- unique(define_)
          lowerbound <- unique(lowerbound)
          upperbound <- unique(upperbound)
          setcoef <- ""
        } else {
          setcoef <- coef
        }
        priors_ <-
          brms::prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = setcoef,
            group = group,
            resp = resp,
            dpar = dpar,
            lb = lowerbound,
            ub = upperbound
          )
      }
      
      if (!ept(mnf) & cov_nlpar == "") {
        if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
          define_ <- unique(define_)
          lowerbound <- unique(lowerbound)
          upperbound <- unique(upperbound)
          setcoef <- ""
        } else {
          if (ept(mnc) == "") {
            setcoef <- coef[1]
          } else {
            setcoef <- coef[-1]
          }
        }
        priors_ <-
          brms::prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = setcoef,
            group = group,
            resp = resp,
            dpar = dpar,
            lb = lowerbound,
            ub = upperbound
          )
      }
      
      
      
      if (nlpar == 's' & cov_nlpar == "") {
        nlpar <- paste0(nlpar, 1:df)
        if (grepl("~1", s_formula_grsi, fixed = T)) {
          if (all(coef == "")) {
            priors_ <- brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
          } else {
            priors_ <-   brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              group = group,
              resp = resp,
              dpar = dpar
            )
          }
        }
        if (grepl("~0", s_formula_grsi, fixed = T)) {
          nlpar <- rep(nlpar ,
                       times = 1,
                       each = length(scovcoefnames_gr))
          if (all(coef == "")) {
            priors_ <- brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
          } else {
            priors_ <-   brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              group = group,
              resp = resp,
              dpar = dpar
            )
          }
        }
      }
      
      
      # nlpar cov a - sd
      if (!a_form_0_gr) {
        if (class == 'sd' & grepl("a_cov", x) & !is.null(a_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- acovcoefnames_gr
          } else {
            coef <- acovcoefnames_gr[2:length(acovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      # nlpar cov b - sd
      if (!b_form_0_gr) {
        if (class == 'sd' & grepl("b_cov", x) & !is.null(b_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- bcovcoefnames_gr
          } else {
            coef <- bcovcoefnames_gr[2:length(bcovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      # nlpar cov c - sd
      if (!c_form_0_gr) {
        if (class == 'sd' & grepl("c_cov", x) & !is.null(c_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- ccovcoefnames_gr
          } else {
            coef <- ccovcoefnames_gr[2:length(ccovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      
      # nlpar cov d - sd
      if (!d_form_0_gr) {
        if (class == 'sd' & grepl("d_cov", x) & !is.null(d_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- dcovcoefnames_gr
          } else {
            coef <- dcovcoefnames_gr[2:length(dcovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov e - sd
      if (!e_form_0_gr) {
        if (class == 'sd' & grepl("e_cov", x) & !is.null(e_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- ecovcoefnames_gr
          } else {
            coef <- ecovcoefnames_gr[2:length(ecovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      
      # nlpar cov f - sd
      if (!f_form_0_gr) {
        if (class == 'sd' & grepl("f_cov", x) & !is.null(f_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- fcovcoefnames_gr
          } else {
            coef <- fcovcoefnames_gr[2:length(fcovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      
      # nlpar cov g - sd
      if (!g_form_0_gr) {
        if (class == 'sd' & grepl("g_cov", x) & !is.null(g_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- gcovcoefnames_gr
          } else {
            coef <- gcovcoefnames_gr[2:length(gcovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov h - sd
      if (!h_form_0_gr) {
        if (class == 'sd' & grepl("h_cov", x) & !is.null(h_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- hcovcoefnames_gr
          } else {
            coef <- hcovcoefnames_gr[2:length(hcovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov i - sd
      if (!h_form_0_gr) {
        if (class == 'sd' & grepl("h_cov", x) & !is.null(h_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- hcovcoefnames_gr
          } else {
            coef <- hcovcoefnames_gr[2:length(hcovcoefnames_gr)]
          }
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      
      # nlpar cov s - betas
      if (!s_form_0_gr) {
        if (class == 'b' & grepl("s_cov", x) & !is.null(s_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- scovcoefnames_gr
          } else {
            coef <- scovcoefnames_gr[2:length(scovcoefnames_gr)]
          }
          nlpar <- paste0(nlpar, 1:df)
          nlpar <- rep(nlpar, times = length(coef), each = 1)
          coef <- rep(coef , times = 1, each = df)
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          
        }
      }
      
      
      
      
      # sigma cov - sd
      if(!is.null(sigma_formula_grsi)) {
        if (!grepl("~0", sigma_formula_grsi, fixed = T)) {
          if (class == 'sd' & grepl("sigma_cov", x) & 
              !is.null(sigma_cov_prior_sd)) {
            if (ept(mnf)) {
              coef <- sigmacovcoefnames_gr
            } else {
              coef <- sigmacovcoefnames_gr[2:length(sigmacovcoefnames_gr)]
            }
            dpar <- sigma_dpar
            priors_ <-
              brms::prior_string(
                define_,
                class = class,
                nlpar = nlpar,
                group = group,
                coef = coef,
                resp = resp,
                dpar = dpar
              )
          }
        }
      }
      
      
      
    } # end if (class == 'sd') {
    
    
    
    # correlation priors (lkj)
    
    # currently brms does not allow setting separate ljk prior for
    # subset and multivariate
    # removing resp leads to duplicate priors, so need to set only once
    
    # Note that currently brms not allowing setting separate cor prior for sigma
    # So, set sigma_prior_cor <- NULL (line 810) for all, or else set group = '
    # But then this won't let assign different prior - dup stanvars
    # So keeping group 
    
    if (class == 'cor') {
      # if(sigma_dpar == "sigma") group <- ""
      if (ii == 1) {
        priors_ <-  brms::prior_string(define_,
                                       class = class,
                                       group = group,
                                       dpar = dpar) # resp = resp,
      } else {
        priors_ <- ""
      }
    }
    
    
    if (class == 'rescor') {
      if (ii == 1) {
        priors_ <-  brms::prior_string(define_,
                                       class = class,
                                       group = "",
                                       dpar = "")
      } else {
        priors_ <- ""
      }
    }
    
    
    
    # residual standard deviation (sigma) prior
    
    if (class == 'sigma' & dpar == "") {
      priors_ <-  brms::prior_string(
        define_,
        class = class,
        lb = lowerbound,
        ub = upperbound,
        resp = resp,
        dpar = dpar
      )
      
    }
    
    # residual standard deviation (sigma) prior - dpar_formula formulation
    
    if (class == "" & !is.null(dpar_formulasi)) {
      # need to remove lb and ub if specifying coef, otherwise
      # Error: Argument 'coef' may not be specified when using boundaries.
      dpar <- 'sigma'
      if (!is.null(dpar_covi_mat_form) &
          !grepl("~1$", dpar_covi_mat_form, fixed = F)) {
        class <- 'b'
        mnf <- paste0('dpar', "_form_0")
        mnc <- paste0("dpar_cov")
        
        if (all(is.na(lowerbound)) |
            all(is.na(upperbound))) {
          if (grepl("^lf\\(", dpar_formulasi)) {
            if (grepl("cmc=F", dpar_formulasi) |
                grepl("cmc=FALSE", dpar_formulasi)) {
              coef <- c("", dparcovcoefnames[2:length(dparcovcoefnames)])
              define_ <- c("", define_[2:length(define_)])
            } else {
              coef <- dparcovcoefnames
            }
          }
          if (!grepl("^lf\\(", dpar_formulasi) |
              !grepl("^nlf\\(", dpar_formulasi)) {
            coef <- dparcovcoefnames
          }
        } else if (!all(is.na(lowerbound)) & !all(is.na(upperbound))) {
          coef <- rep("", length(dparcovcoefnames))
        }
        
        if (ept(mnf)) {
          if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
            define_ <- unique(define_)
            lowerbound <- unique(lowerbound)
            upperbound <- unique(upperbound)
            setcoef <- ""
          } else {
            setcoef <- coef
          }
          
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = setcoef,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
        }
        
        if (!ept(mnf)) {
          if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
            define_ <- unique(define_)
            lowerbound <- unique(lowerbound)
            upperbound <- unique(upperbound)
            setcoef <- ""
          } else {
            if (ept(mnc) != "") {
              setcoef <- coef[1]
            } else {
              setcoef <- coef
            }
          }
          
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = setcoef,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
        }
        
        
        
        
        # residual standard deviation (sigma) covariate prior -
        # dpar_formula formulation
        
        if (!is.null(dpar_covi_mat_form) &
            grepl("~1", dpar_covi_mat_form, fixed = T) &
            !grepl("~1$", dpar_covi_mat_form, fixed = T) &
            !is.null(dpar_cov_prior_sigma)) {
          if (grepl("dpar", x) & !grepl("dpar_cov", x)) {
            if (grepl("^lf\\(", dpar_formulasi)) {
              if (grepl("center=T", dpar_formulasi) |
                  grepl("center=TRUE", dpar_formulasi)) {
                class <- dparcovcoefnames[1]
                coef  <- ""
              } else {
                class <- 'b'
                coef  <- dparcovcoefnames[1]
              }
            } else if (!grepl("^lf\\(", dpar_formulasi) |
                       !grepl("^nlf\\(", dpar_formulasi)) {
              class <- dparcovcoefnames[1]
              coef  <- ""
            }
          }
          
          if (!grepl("dpar", x) & grepl("dpar_cov", x)) {
            class <- 'b'
            coef <- dparcovcoefnames
          }
          
          if (class == 'b') {
            if (grepl("center=T", dpar_formulasi) |
                grepl("center=TRUE", dpar_formulasi)) {
              coef <- coef[-1]
            } else {
              coef <- coef
            }
          }
          if (!grepl("center=T", dpar_formulasi) &
              !grepl("center=TRUE", dpar_formulasi)) {
            if (grepl("dpar_cov", x))
              coef <- coef[-1]
          }
          
          priors_ <-
            brms::prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      
      if (!is.null(dpar_covi_mat_form) &
          grepl("~1$", dpar_covi_mat_form, fixed = F)) {
        if (grepl("center=T", dpar_formulasi) |
            grepl("center=TRUE", dpar_formulasi)) {
          class <- dparcovcoefnames[1]
          coef  <- ""
        } else {
          class <- 'b'
          coef  <- dparcovcoefnames[1]
        }
        
        priors_ <-
          brms::prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = coef,
            resp = resp,
            dpar = dpar
          )
      }
    }
    
    
    # autocorrelation priors
    
    if (setautocorr) {
      coef <- ""
      if (acorclass == 'arma') {
        acorclassclasses <- c("ar", "ma")
        priors_arma_c <- list()
        stanvars_data_in_c <- list()
        for (acorclassi in acorclassclasses) {
          define_ <- priors_arma_c_define[[acorclassi]]$define_
          lowerbound <-
            priors_arma_c_define[[acorclassi]]$lowerbound
          upperbound <-
            priors_arma_c_define[[acorclassi]]$upperbound
          priors_temp <-  brms::prior_string(
            define_,
            class = acorclassi,
            lb = lowerbound,
            ub = upperbound,
            coef = coef,
            resp = resp,
            dpar = dpar
          )
          priors_arma_c[[acorclassi]] <- priors_temp
          stanvars_data_in_c[[acorclassi]] <-
            priors_arma_c_define[[acorclassi]]$stanvars_data_in
        }
        priors_ <- priors_arma_c %>% do.call(rbind, .)
        stanvars_data_in <- stanvars_data_in_c %>% do.call(rbind, .)
      } else {
        priors_ <-  brms::prior_string(
          define_,
          class = class,
          lb = lowerbound,
          ub = upperbound,
          coef = coef,
          resp = resp,
          dpar = dpar
        )
      }
    }
    out_pr <-
      list(
        priors_ = priors_,
        stanvars_data_in = stanvars_data_in,
        initial_in = initial_in
      )
    return(out_pr)
  } 
  
  
  
  # use following custom order
  # This order ensures that corresponding initial arguments are matched
  # with the sequence of prior argument evaluation
  
  custom_order_prior <- c(
    'a_prior_beta',
    'a_cov_prior_beta',
    'b_prior_beta',
    'b_cov_prior_beta',
    'c_prior_beta',
    'c_cov_prior_beta',
    'd_prior_beta',
    'd_cov_prior_beta',
    'e_prior_beta',
    'e_cov_prior_beta',
    'f_prior_beta',
    'f_cov_prior_beta',
    'g_prior_beta',
    'g_cov_prior_beta',
    'h_prior_beta',
    'h_cov_prior_beta',
    'i_prior_beta',
    'i_cov_prior_beta',
    's_prior_beta',
    's_cov_prior_beta',
    'a_prior_sd',
    'a_cov_prior_sd',
    'b_prior_sd',
    'b_cov_prior_sd',
    'c_prior_sd',
    'c_cov_prior_sd',
    'd_prior_sd',
    'd_cov_prior_sd',
    'e_prior_sd',
    'e_cov_prior_sd',
    'f_prior_sd',
    'f_cov_prior_sd',
    'g_prior_sd',
    'g_cov_prior_sd',
    'h_prior_sd',
    'h_cov_prior_sd',
    'i_prior_sd',
    'i_cov_prior_sd',
    's_prior_sd',
    's_cov_prior_sd',
    'sigma_prior_beta',
    'sigma_cov_prior_beta',
    'sigma_prior_sd',
    'sigma_cov_prior_sd',
    
    'gr_prior_cor',
    'sigma_prior_cor',
    
    'rsd_prior_sigma',
    'dpar_prior_sigma',
    'dpar_cov_prior_sigma',
    'autocor_prior_acor',
    'autocor_prior_unstr_acor',
    'mvr_prior_rescor'
  )
  
  
  custom_order_prior_str_evaluating <- FALSE
  if(any(custom_order_prior_str != "")) {
    custom_order_prior <- custom_order_prior_str
    custom_order_prior_str_evaluating <- TRUE
  }
  
  
  
  stanvars_data_5 <- list()
  initial_in_data <- list()
  c_priors <- list()
  for (ip in custom_order_prior) {
    if (grepl("_prior_", ip)) {
      if (!is.null(eval(parse(text = ip)))) {
        xx   <- deparse(ip)
        pget <- eval_prior_args(eval(parse(text = xx)))
        zz   <- pget$priors_
        stanvars_data_5[[ip]] <- pget$stanvars_data_in
        initial_in_data[[ip]] <- pget$initial_in
        c_priors[[ip]] <- zz
      }
    }
  }
  
  evaluated_priors <- c_priors %>% do.call(rbind, .)
  

  newlist <- c()
  for (i in 1:length(stanvars_data_5)) {
    ttt <- stanvars_data_5[[i]][1:length(stanvars_data_5[[i]])]
    newlist <- c(newlist, unname(ttt))
  }
  
  svardatalistlist <- c()
  for (istanvardata in 1:length(newlist)) {
    svardatalistlist[istanvardata] <-
      paste0("newlist[[", istanvardata, "]]")
  }
  
  stanvars <-
    eval(parse(text = paste(svardatalistlist, collapse = "+")))
  
  attr(evaluated_priors, 'stanvars') <- stanvars
  
  initial_in_datazz <- initial_in_data
  
  if (is.list(initial_in_datazz) & length(initial_in_datazz) == 0) {
    initial_in_datazz <- NULL
  }
  

  if(!is.null(initsi[[1]])) {
    if(initsi[[1]] == 'random') {
      initial_in_datazz <- NULL
      combined_inits <- NULL
    }
  }
  
  
  
  if (!is.null(initial_in_datazz)) {
    if (!is.null(gr_prior_cor) | !is.null(sigma_prior_cor) ) {
      list_ck <- list_ck_ <- list()
      list_ck_rescor <- list()
      ik_j <- ik_j_ <- 0
      what_not_to_flatten <- "^L_|^z_|Intercept_sigma|Lrescor"
      what_not_to_flatten2 <- "^L_|^z_"
      if(length(initial_in_datazz$a_prior_beta) == 0) {
        # added |b_s when a_init is random 
        # 06 11 2023 detected -> setting a_init to random also set s_init random
        # This should be paste0(what_not_to_flatten, "|b_a") 
        # and paste0(what_not_to_flatten, "|b_s")
        what_not_to_flatten <- paste0(what_not_to_flatten, "|b_a")
        what_not_to_flatten2 <- paste0(what_not_to_flatten2, "|b_s")
      }
      for (ik in 1:length(initial_in_datazz)) {
        ik_j <- ik_j + 1
        ik_names <- names(initial_in_datazz[[ik]])
        if (!any(grepl(what_not_to_flatten, ik_names))) {
          list_ck[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck[[ik]]) <- ik_names
        } else if (any(grepl(what_not_to_flatten2, ik_names))) {
          mn <- 0
          for (ikl in 1:length(grepl(what_not_to_flatten2, ik_names))) {
            mn <- mn + 1
            ik_j_ <- ik_j_ + 1
            list_ck_[[ik_j_]] <- initial_in_datazz[[ik_j]][[mn]]
          }
          names(list_ck_) <- ik_names
        } else if (grepl("^Intercept_sigma", ik_names)) {
          list_ck[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck[[ik]]) <- ik_names
        } else if (multivariate$mvar & multivariate$rescor &
                   grepl("^Lrescor", ik_names)) {
          list_ck_rescor[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck_rescor[[ik]]) <- ik_names
        }
      }
      list_ck <- list_ck[lengths(list_ck) != 0]
      keys    <- unique(unlist(lapply(list_ck, names)))
      list_ck <-
        setNames(do.call(mapply, c(FUN = c, lapply(
          list_ck, `[`, keys
        ))), keys)
      combined_inits <- c(list_ck, list_ck_)
    }
    

    
    if (is.null(gr_prior_cor) & is.null(sigma_prior_cor) ) {
      list_ck <- list_ck_z <- list_ck_sd <- list()
      list_ck_rescor <- list()
      ik_j <- ik_j_ <- 0
      for (ik in 1:length(initial_in_datazz)) {
        ik_j <- ik_j + 1
        ik_names <- names(initial_in_datazz[[ik]])
        if (!any(grepl("^L_|^z_|Intercept_sigma|Lrescor", ik_names))) {
          list_ck[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck[[ik]]) <- ik_names
        } else if (any(grepl("^L_|^z_", ik_names))) {
          mn <- 0
          for (ikl in 1:length(grepl("^L_|^z_", ik_names))) {
            mn <- mn + 1
            ik_j_ <- ik_j_ + 1
            if (is.matrix(initial_in_datazz[[ik_j]][[mn]]))
              list_ck_z[[ik_j_]] <- initial_in_datazz[[ik_j]][[mn]]
            if (!is.matrix(initial_in_datazz[[ik_j]][[mn]]))
              list_ck_sd[[ik_j_]] <- initial_in_datazz[[ik_j]][[mn]]
          }
          
          names(list_ck_z) <- ik_names[2]
          names(list_ck_sd) <- ik_names[1]
        } else if (grepl("^Intercept_sigma", ik_names)) {
          list_ck[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck[[ik]]) <- ik_names
        } else if (multivariate$mvar & multivariate$rescor &
                   grepl("^Lrescor", ik_names)) {
          list_ck_rescor[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck_rescor[[ik]]) <- ik_names
        }
      }
      list_ck <- list_ck[lengths(list_ck) != 0]
      keys    <- unique(unlist(lapply(list_ck, names)))
      list_ck <-
        setNames(do.call(mapply, c(FUN = c, lapply(
          list_ck, `[`, keys
        ))), keys)
      list_ck_sd <- list_ck_sd[lengths(list_ck_sd) != 0]
      # this on 9 5 23 to accomodate random = ''
      if(length(list_ck_sd) != 0) {
        for (list_ck_sd_i in 1:length(list_ck_sd)) {
          if (length(list_ck_sd[[list_ck_sd_i]]) > 1) {
            nami_ <-
              paste0(names(list_ck_sd[[list_ck_sd_i]][1]),
                     "cov",
                     2:length(list_ck_sd[[list_ck_sd_i]]) - 1)
            names(list_ck_sd[[list_ck_sd_i]]) <-
              c(names(list_ck_sd[[list_ck_sd_i]][1]), nami_)
          }
        }
      }
      names(list_ck_sd) <-
        rep(names(list_ck_sd[1]), length(list_ck_sd))
      list_ck_sd2 <- list_ck_sd
      list_ck_z <- list_ck_z[lengths(list_ck_z) != 0]
      list_ck_z2 <- list()
      # this on 9 5 23 to accomodate random = ''
      if(length(list_ck_z) != 0) {
        for (list_ck_i in 1:length(list_ck_z)) {
          addelemnt <-
            strsplit(gsub("\\+", " ", randomsi), " ")[[1]][list_ck_i]
          list_ck_z2[[paste0("z", "_", addelemnt, resp_, list_ck_i)]] <-
            list_ck_z[[list_ck_i]]
          attr(list_ck_z2[[paste0("z", "_", addelemnt, list_ck_i)]], "names") <-
            NULL
        }
      }
      combined_inits <- c(list_ck, list_ck_sd2, list_ck_z2)
    }
    
    # Don't let it evaluate when custom_order_prior_str != "" 
    # i.e,  when evaluating when hierarchy priors 
    
    if (multivariate$mvar & 
        multivariate$rescor & 
        !custom_order_prior_str_evaluating) {
      list_ck_rescor <- list_ck_rescor[lengths(list_ck_rescor) != 0]
      list_ck_rescor <- list_ck_rescor[[1]]
      combined_inits <- c(combined_inits, list_ck_rescor)
    }
    
    
    # convert vector of 's' initials to named individual (s1, s2)
    if(select_model == "sitar" | select_model == 'rcs') {
      # Don't let when evaluating _str higher custom order
      if("s_prior_beta" %in% custom_order_prior) {
        first_loop <- TRUE  
      } else {
        first_loop <- FALSE  
      }
      if(first_loop) {
        nlpar_s_init <- paste0('_s', 1:df)
        if (grepl("~0", s_formulasi, fixed = T)) {
          nlpar_s_init <-
            rep(nlpar_s_init ,
                times = 1,
                each = length(scovcoefnames))
        } else if (!grepl("~0", s_formulasi, fixed = T)) {
          nlpar_s_init <- rep(nlpar_s_init , times = length(scovcoefnames))
        }
        
        
        subset_sparms <-
          combined_inits[grepl(".*_s$", names(combined_inits))]
        
        subset_sparms_name <- names(subset_sparms)
        if(length(subset_sparms) != 0) subset_sparms_numeric <- 
          subset_sparms[[1]]
        if(length(subset_sparms) == 0) subset_sparms_numeric <- NULL
        if(!is.null(subset_sparms_numeric)) {
          subset_sparms2 <- list()
          subset_sparms2names <- c()
          
          for (subset_sparmsi in 1:length(subset_sparms_numeric)) {
            subset_sparms_namei <-
              gsub("_s", nlpar_s_init[subset_sparmsi], subset_sparms_name)
            subset_sparms2[[subset_sparms_namei]] <-
              subset_sparms_numeric[subset_sparmsi]
            subset_sparms2names <-
              c(subset_sparms2names, subset_sparms_namei)
          }
          names(subset_sparms_numeric) <- subset_sparms2names
          subset_sparms3 <- list()
          for (isi in 1:df) {
            subset_sparms3[[paste0("b", resp_, "_s", isi)]] <-
              subset_sparms_numeric[grep(paste0("b", resp_, "_s", isi),
                                         names(subset_sparms_numeric))]
          }
          subset_sparms <- subset_sparms3
          subset_sparms <-
            subset_sparms[!names(subset_sparms) %in% subset_sparms_name]
          combined_inits <-
            append(combined_inits, subset_sparms, after = grep(
              paste0("^", subset_sparms_name, "$"),
              names(combined_inits)
            ))
          combined_inits <-
            combined_inits[!names(combined_inits) %in% 
                             paste0("",
                                    subset_sparms_name, "")]
          initials <- combined_inits
        } # if(!is.null(subset_sparms_numeric)) {
        if(is.null(subset_sparms_numeric)) {
          initials <- combined_inits
        }
      } # if(first_loop) {
      
      if(!first_loop) initials <- combined_inits
      
    } # if(select_model == "sitar") {
    
    if(select_model != "sitar" & select_model != 'rcs') initials <- combined_inits
    
  } # if (!is.null(initial_in_datazz)) {
  
  
  
  # Mean all initals random
  if(length(combined_inits) == 0) initials <- NULL
  
  if (is.null(initial_in_datazz)) {
    initials <- NULL
  }
  
  ###################3
  
  stanvar_priors_names <- names(stanvars)
  getaux <- "tau"
  stanvar_priors_names_c <- c()
  for (stanvar_priors_namesi in stanvar_priors_names) {
    t <-
      stanvar_priors_namesi[grep(paste0(getaux, '_scale', resp_), 
                                 stanvar_priors_namesi)]
    t <- gsub(paste0('_scale', resp_), "", t, fixed = T)
    stanvar_priors_names_c <- c(stanvar_priors_names_c, t)
  }
  
  add_tau <- list()
  for (stanvar_priors_names_ci in stanvar_priors_names_c) {
    fstandat <-
      unlist(stanvars)[grep(paste0(
        stanvar_priors_names_ci,
        paste0('_scale', resp_, ".sdata")
      ), names(unlist(stanvars)))] %>% as.numeric()
    add_tau[[paste0(stanvar_priors_names_ci, resp_)]] <-
      rep(1, length(fstandat))
  }
  if (length(add_tau) == 0)
    add_tau <- NULL
  
  getaux <- "nu"
  stanvar_priors_names_c <- c()
  for (stanvar_priors_namesi in stanvar_priors_names) {
    t <-
      stanvar_priors_namesi[grep(paste0(getaux, '_scale', resp_), 
                                 stanvar_priors_namesi)]
    t <- gsub(paste0('_scale', resp_), "", t, fixed = T)
    stanvar_priors_names_c <- c(stanvar_priors_names_c, t)
  }
  add_nu <- list()
  for (stanvar_priors_names_ci in stanvar_priors_names_c) {
    add_nu[[paste0(stanvar_priors_names_ci, resp_)]] <-  5
  }
  if (length(add_nu) == 0)
    add_nu <- NULL
  
  initials <- c(initials, add_tau, add_nu)
  
  ################
  revSubstr <- function(x_) {
    x__ <- substr(x_, start = 1, stop = 3)
    x___ <- paste0(rev(strsplit(x__, "_")[[1]]), collapse = "_")
    x___ <- gsub(x__, x___, x_)
    x___
  }
  
  tau_nu_init_list <- c(add_tau, add_nu)
  
  if (length(tau_nu_init_list) != 0) {
    names_tau_nu_parms <- names(tau_nu_init_list)
    names_tau_nu_parmsi_c <- c()
    for (names_tau_nu_parmsi in names_tau_nu_parms) {
      plength <- length(tau_nu_init_list[[names_tau_nu_parmsi]])
      revstr <- revSubstr(names_tau_nu_parmsi)
      if (!grepl("^b_b", names_tau_nu_parmsi, fixed = F)) {
        o <-
          paste0("vector[",
                 plength,
                 "]",
                 " ",
                 revstr,
                 " = ",
                 names_tau_nu_parmsi,
                 ";")
        names_tau_nu_parmsi_c <- c(names_tau_nu_parmsi_c, o)
      }
    }
    names_tau_nu_parmsi_c <- names_tau_nu_parmsi_c
    names_tau_nu_parmsi_cc <-
      paste(names_tau_nu_parmsi_c, collapse = "\n")
    
    scode_auxillary <-
      brms::stanvar(scode = names_tau_nu_parmsi_cc,
                    block = "genquant",
                    position = 'end')
  } else if (length(tau_nu_init_list) == 0) {
    scode_auxillary <- NULL
  }
  
  ##################
  out_listx <- initials
  for (ili in names(initials)) {
    if(length(out_listx[[ili]]) == 1) {
      out_listx[[ili]] <- out_listx[[ili]]
      # here also need array for a b c d e
      # but not for sigma when sigma ~ not used but default rsd formulation 
      if(nys == 1) sigma_par_name_rsd <- "sigma"
      if(nys > 1) sigma_par_name_rsd <- paste0('sigma', resp_)
      if(ili != sigma_par_name_rsd) {
        out_listx[[ili]] <- array(out_listx[[ili]], 
                                  dim = length(out_listx[[ili]]))
      }
    } else if(length(out_listx[[ili]]) > 1 & is.vector(out_listx[[ili]])) {
      out_listx[[ili]] <- array(out_listx[[ili]], 
                                dim = length(out_listx[[ili]]))
    }
    
    if(is.na(ili)) ili <- "xxxxxxxxxxxxxx"
    # for ar and ma, it is always vector , so array
    if(ili == 'ar' | ili == 'ma') {
      out_listx[[ili]] <- array(out_listx[[ili]], 
                                dim = length(out_listx[[ili]]))
    }
  }
  
  initials <- out_listx
  
  # when sigma  formula is ~1+.., then first element is Intercept_sigma and the 
  # remaining are b_sigma
  
  initialsx <- out_listx
  if(!sigma_form_0) {
    if(nys == 1) {
      sigma_par_name <- 'b_sigma'
      Intercept_sigma <- 'Intercept_sigma'
    } else if(nys > 1) {
      sigma_par_name <- paste0('b_sigma', resp_)
      Intercept_sigma <- paste0('Intercept_sigma', resp_)
    }
    if(!is.null(initialsx[[sigma_par_name]])) {
      g_sigma_i <- initialsx[[sigma_par_name]]
      initialsx[[Intercept_sigma]] <- g_sigma_i[1]
      if(length(g_sigma_i) > 1) {
        initialsx[[sigma_par_name]] <-  
          array(g_sigma_i[2:length(g_sigma_i)], 
                dim = length(g_sigma_i[2:length(g_sigma_i)]))
      } # if(length(g_sigma_i) > 1) {
    }
  }
  
  
  initials <- initialsx
  
  
  
  # Re create symmetric square Lcortime which is flattened to a vector
  # This could have been done above like L_|z_|Lrescor etc but this is same
  if(!is.null(initials[['Lcortime']])) {
    NC_dims         <- ept(cortimeNlags) %>% as.numeric()
    initials[['Lcortime']] <- matrix(initials[['Lcortime']], 
                                     nrow = NC_dims, 
                                     ncol = NC_dims)
  } # if(!is.null(initials[['Lcortime']])) {
  
  
  
  attr(evaluated_priors, 'initials') <- initials
  attr(evaluated_priors, 'scode_auxillary') <- scode_auxillary
  
  return(evaluated_priors)
}












#' An internal function to prepare priors for Bayesian SITAR growth curve model
#' 
#' For \code{univariate_by} and \code{multivariate} models (see [bsitar::bsitar()])
#' each argument is automatically matched with the sub model. 
#'
#' @param prior_argument A list containing the prior arguments specified in  
#' the [bsitar::bsitar()] function and then passed from the 
#' [bsitar::set_priors_initials()] function to the \code{prepare_priors}. 
#' 
#' @param prior_data An optional argument (a named list) specified in  
#' the [bsitar::bsitar()] function and then passed from the 
#' [bsitar::set_priors_initials()] function to the \code{prepare_priors}. 
#' The \code{prior_data} used to pass value for priors.
#'  See [bsitar::bsitar()] function, \code{prior_data} for details.
#'  
#' @param prior_data_internal An internal argument (as named list) specified in  
#' the [bsitar::bsitar()] function and then passed from the 
#' [bsitar::set_priors_initials()] function to the \code{prepare_priors}. 
#' 
#' @param prior_internal_args An internal argument list that is passed from the 
#' [bsitar::set_priors_initials()] function to the \code{set_priors_initials} 
#' and is used in setting the priors.
#' 
#' @param init_arguments A list containing the initial arguments specified in  
#' the [bsitar::bsitar()] function and then passed from the 
#' [bsitar::set_priors_initials()] function to the \code{prepare_priors}. 
#' 
#' @param init_data An optional argument (as named list) specified in  
#' the [bsitar::bsitar()] function and then passed from the 
#' [bsitar::set_priors_initials()] function to the \code{prepare_priors}. 
#' The \code{init_data} is used for setting the initials.
#' 
#' @param init_data_internal An internal argument (as named list) specified in  
#' the [bsitar::bsitar()] function and then passed from the 
#' [bsitar::set_priors_initials()] function to the \code{prepare_priors}. 
#' 
#' @param init_args_internal An internal argument list that is passed from the 
#' [bsitar::set_priors_initials()] function to the \code{set_priors_initials} 
#' and is used in setting the initials.
#'
#' @return An object of class \code{brmsprior}. See \code{brmsprior} function 
#' for more details. 
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
prepare_priors <- function(prior_argument,
                           prior_data,
                           prior_data_internal,
                           prior_internal_args,
                           init_arguments,
                           init_data,
                           init_data_internal,
                           init_args_internal) {
  
  
  # Initiate non formalArgs()
  nlpar <- NULL;
  verbose <- NULL;
  cov_nlpar <- NULL;
  dpar <- NULL;
  cov_dpar <- NULL;
  fixedsi <- NULL;
  a_form_0 <- NULL;
  resp_ <- NULL;
  nrep_of_parms <- NULL;
  b_form_0 <- NULL;
  c_form_0 <- NULL;
  d_form_0 <- NULL;
  e_form_0 <- NULL;
  f_form_0 <- NULL;
  g_form_0 <- NULL;
  h_form_0 <- NULL;
  i_form_0 <- NULL;
  s_form_0 <- NULL;
  s_form_0_gr <- NULL;
  sncov <- NULL;
  ancov <- NULL;
  bncov <- NULL;
  cncov <- NULL;
  dncov <- NULL;
  encov <- NULL;
  fncov <- NULL;
  gncov <- NULL;
  hncov <- NULL;
  incov <- NULL;
  sncov <- NULL;
  sigma_dpar <- NULL;
  sigma_form_0 <- NULL;
  cov_sigma_dpar <- NULL;
  sigmancov <- NULL;
  randomsi <- NULL;
  ancov_gr <- NULL;
  bncov_gr <- NULL;
  cncov_gr <- NULL;
  dncov_gr <- NULL;
  encov_gr <- NULL;
  fncov_gr <- NULL;
  gncov_gr <- NULL;
  hncov_gr <- NULL;
  incov_gr <- NULL;
  sncov_gr <- NULL;
  sigmancov_gr <- NULL;
  dparncov <- NULL;
  setautocorr <- NULL;
  group <- NULL;
  normalize <- NULL;
  initsi <- NULL;
  
  eout <- list2env(prior_internal_args)
  for (eoutii in names(eout)) {
    assign(eoutii, eout[[eoutii]])
  }
  
  
  eout <- list2env(init_arguments)
  for (eoutii in names(eout)) {
    assign(eoutii, eout[[eoutii]])
  }
  
  
  empty_sufx <- NULL
  
  # set to 'character(0)' to avoid overhead of reduce sum
  # set to '' to mimic default behavious whihc adds pll_args for data block
  # if change_default_data_pll_args <- FALSE, then nothing changed, i.e., 
  # default brms behavious
  
  change_default_data_pll_args <- TRUE
  set_data_pll_args <- 'character(0)'
  
  stanvars_data <- list()
  
  eout <- list2env(prior_data_internal)
  for (eoutii in names(eout)) {
    assign(eoutii, eout[[eoutii]])
  }
  
  # evaluate user defined prior_data
  if (!is.null(prior_data[[1]])) {
    eout <- list2env(prior_data)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  
  # get elements
  get_within_fist_last_paranthesese <- function(x__) {
    x__ <- sub('\\(', '[', x__)
    x__ <- sub("\\)([^)]*)$", "]\\1", x__)
    x__ <-
      gsub("[\\[\\]]", "", regmatches(x__, gregexpr("\\[.*?\\]", x__))[[1]])
    x__ <- gsub("\\[|\\]", "", x__)
    x__
  }
  
  gsub_comma_within_paranthesese <-
    function(x__, replace_comma_by) {
      tt <-
        gsub("[\\(\\)]", "", regmatches(x__, gregexpr("\\(.*?\\)", x__))[[1]])
      tt2 <- gsub(",", replace_comma_by, tt, fixed = T)
      j <- 0
      for (i in tt) {
        j <- j + 1
        x__ <- gsub(tt[j], tt2[j], x__, fixed = T)
      }
      x__
    }
  
  sep_indicator <- "_"
  p_str_in <- gsub("\\s", "", prior_str_arg)
  splitmvar <- p_str_in
  splitmvar <- gsub("\\s", "", splitmvar)
  splitmvar <- paste(splitmvar, collapse = "")
  splitmvar_w <- get_within_fist_last_paranthesese(splitmvar)
  
  # This for flat prior when no () i.e, flat instead of flat()
  if (identical(splitmvar_w, character(0))) {
    splitmvar_w <- ""
  }
  
  
  splitmvar_w <-
    gsub_comma_within_paranthesese(splitmvar_w, "_comma_")
  splitmvar_w2 <- strsplit(splitmvar_w, ",")[[1]]
  splitmvar_w2 <- gsub("_comma_" , ",", splitmvar_w2)
  splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
  
  ept <- function(x)
    eval(parse(text = x), envir = parent.frame())
  
  
  # extract sethp distribution
  get_sethp_arg <- splitmvar_w2[grepl("^sethp=", splitmvar_w2)]
  sethp_tt <- sub('.*=', '', get_sethp_arg)
  sethp_tt <- paste0("'", sethp_tt, "'")
  get_sethp_arg <-
    paste0(sub("=[^=]+$", "=", get_sethp_arg),  sethp_tt)
  sethp_dist <- ept(get_sethp_arg)
  sethp_dist <- gsub("\"" , "", sethp_dist)
  
  
  if (!is.null(sethp_dist) &
      (sethp_dist == "T" | sethp_dist == "TRUE")) {
    sethp_dist <- "normal"
    prior_str_arg <-
      gsub(paste0("=", sethp_dist),
           paste0("=", "TRUE"),
           prior_str_arg)
    splitmvar_w2 <- gsub(sethp_dist, "TRUE", splitmvar_w2)
  }
  
  
  if (!is.null(sethp_dist) &
      (
        sethp_dist != "''" &
        sethp_dist != "NA" &
        sethp_dist != "F" &
        sethp_dist != "FALSE" &
        sethp_dist != ""
      )) {
    if (sethp_dist == "normal" |
        sethp_dist == "cauchy" |
        sethp_dist == "student_t" |
        sethp_dist == "student_nu" |
        sethp_dist == "exponential") {
      prior_str_arg <-
        gsub(paste0("=", sethp_dist),
             paste0("=", "TRUE"),
             prior_str_arg)
      splitmvar_w2 <-
        gsub(sethp_dist, "TRUE", splitmvar_w2)
    } else {
      stop(
        "Hierarchial distribution (i.e,, sethp = distribution) can only",
        "\n",
        " be normal, cauchy, student_nu, student_t, exponential",
        "\n",
        " if sethp = NA, then priors are set directly and not as hierarchial",
        "\n",
        " if sethp = TRUE, then priors are set as hierarchial",
        "\n",
        " with default normal distribution"
      )
    }
  }
  
  
  dist <- strsplit(gsub("\\s", "", prior_str_arg), "\\(")[[1]][1]
  
  if(dist != 'flat') {
    
    add_missing_mandate_names <- function(x, testi, testi2) {
      j = 0
      out_c <- c()
      for (i in 1:length(x)) {
        j <- j + 1
        if (!x[j] %in% testi[j]) {
          out <- paste0(x[j], "=", testi[j])
        } else {
          out <- testi2[j]
        }
        out_c <- c(out_c, out)
      }
      lengthii <- length(testi)
      testi3 <- testi2[length(out_c) + 1:length(testi2)]
      out_c <- c(out_c, testi3)
      out_c <- head(out_c, lengthii)
      out_c
    }
    
    
    error_handle1 <- function(i, dist, testi) {
      if (!i %in% testi) {
        stop(
          "for ",
          dist,
          " distribution,",
          " parameter '",
          i,
          "' is missing",
          "\n  Please specify ",
          dist,
          " distribution as: ",
          paste0(dist, "(", paste(set_str_names, collapse = ", "), ")")
        )
      }
    }
    
    
    if (dist == "normal" |
        dist == "cauchy" |
        dist == "lognormal") {
      set_str_names <- c("location", "scale")
      if (length(splitmvar_w2) < length(set_str_names))
        stop(
          "please sepecify minimum required ",
          length(set_str_names),
          " parameters, i.e., ",
          paste(set_str_names, collapse = ", ")
        )
      splitmvar_w2 <-
        add_missing_mandate_names(set_str_names, splitmvar_w3, splitmvar_w2)
      splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
      for (i in set_str_names)
        error_handle1(i, dist, splitmvar_w3)
    }
    
    if (dist == "gamma") {
      set_str_names <- c("shape", "scale")
      if (length(splitmvar_w2) < length(set_str_names))
        stop(
          "please sepecify minimum required ",
          length(set_str_names),
          " parameters, i.e., ",
          paste(set_str_names, collapse = ", ")
        )
      splitmvar_w2 <-
        add_missing_mandate_names(set_str_names, splitmvar_w3, splitmvar_w2)
      splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
      for (i in set_str_names)
        error_handle1(i, dist, splitmvar_w3)
    }
    
    
    if (dist == "inv_gamma") {
      set_str_names <- c("shape", "scale")
      if (length(splitmvar_w2) < length(set_str_names))
        stop(
          "please sepecify minimum required ",
          length(set_str_names),
          " parameters, i.e., ",
          paste(set_str_names, collapse = ", ")
        )
      splitmvar_w2 <-
        add_missing_mandate_names(set_str_names, splitmvar_w3, splitmvar_w2)
      splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
      for (i in set_str_names)
        error_handle1(i, dist, splitmvar_w3)
    }
    
    if (dist == "uniform") {
      set_str_names <- c("lower", "upper")
      if (length(splitmvar_w2) < length(set_str_names))
        stop(
          "please sepecify minimum required ",
          length(set_str_names),
          " parameters, i.e., ",
          paste(set_str_names, collapse = ", ")
        )
      splitmvar_w2 <-
        add_missing_mandate_names(set_str_names, splitmvar_w3, splitmvar_w2)
      splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
      for (i in set_str_names)
        error_handle1(i, dist, splitmvar_w3)
    }
    
    if (dist == "exponential") {
      set_str_names <- c("rate")
      if (length(splitmvar_w2) < length(set_str_names))
        stop(
          "please sepecify minimum required ",
          length(set_str_names),
          " parameters, i.e., ",
          paste(set_str_names, collapse = ", ")
        )
      splitmvar_w2 <-
        add_missing_mandate_names(set_str_names, splitmvar_w3, splitmvar_w2)
      splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
      for (i in set_str_names)
        error_handle1(i, dist, splitmvar_w3)
    }
    
    if (dist == "student_t") {
      set_str_names <- c("df", "location", "scale")
      if (length(splitmvar_w2) < length(set_str_names))
        stop(
          "please sepecify minimum required ",
          length(set_str_names),
          " parameters, i.e., ",
          paste(set_str_names, collapse = ", ")
        )
      splitmvar_w2 <-
        add_missing_mandate_names(set_str_names, splitmvar_w3, splitmvar_w2)
      splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
      for (i in set_str_names)
        error_handle1(i, dist, splitmvar_w3)
    }
    
    if (dist == "student_nu") {
      set_str_names <- c("nu_shape", "nu_scale", "location", "scale")
      if (length(splitmvar_w2) < length(set_str_names))
        stop(
          "please sepecify minimum required ",
          length(set_str_names),
          " parameters, i.e., ",
          paste(set_str_names, collapse = ", ")
        )
      splitmvar_w2 <-
        add_missing_mandate_names(set_str_names, splitmvar_w3, splitmvar_w2)
      splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
      for (i in set_str_names)
        error_handle1(i, dist, splitmvar_w3)
    }
    if (dist == "lkj") {
      set_str_names <- c("eta")
      if (length(splitmvar_w2) < length(set_str_names))
        stop(
          "please sepecify minimum required ",
          length(set_str_names),
          " parameters, i.e., ",
          paste(set_str_names, collapse = ", ")
        )
      splitmvar_w2 <-
        add_missing_mandate_names(set_str_names, splitmvar_w3, splitmvar_w2)
      splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
      for (i in set_str_names)
        error_handle1(i, dist, splitmvar_w3)
    }
    
    
    vacoublary_prior_parnames <- c(
      "location",
      "scale",
      "df",
      "nu_shape",
      "nu_scale",
      "rate",
      "shape",
      "lower",
      "upper",
      "eta",
      "lb",
      "ub",
      "autoscale",
      "addrange",
      "sethp"
    )
    
    
    optional_prior_names <-
      c("lb", "ub", "autoscale", "addrange", "sethp")
    
    # add missing optional_prior_names
    missing_optional_prior_names <-
      optional_prior_names[!optional_prior_names %in% splitmvar_w3]
    
    if (!identical(missing_optional_prior_names, character(0))) {
      splitmvar_w2 <-
        c(splitmvar_w2,
          paste0(missing_optional_prior_names, "=", "NA"))
      splitmvar_w3 <-
        c(splitmvar_w3,
          paste0(missing_optional_prior_names, "", ""))
    }
    
    
    
    min_par_names <-
      names(Filter(function(check__)
        check__ > 0, colSums(
          sapply(splitmvar_w3, grepl, vacoublary_prior_parnames)
        )))
    incorrect_names <- splitmvar_w3[!splitmvar_w3 %in% min_par_names]
    
    if (!identical(incorrect_names, character(0))) {
      ttt_n1 <- paste(incorrect_names, collapse = ", ")
      ttt_nn2 <- paste(vacoublary_prior_parnames, collapse = ", ")
      stop(
        "\nFollowing prior parameter name(s) are incorrect/misspelled:\n ",
        ttt_n1,
        "\n",
        "Available prior parameter names are:\n",
        ttt_nn2,
        "\n",
        "For ",
        dist,
        " prior distribution, the mandatory parameter are: ",
        paste(set_str_names, collapse = ", "),
        "\n",
        paste0(dist, "(", paste(set_str_names, collapse = ", "), ")")
      )
    }
    
    
    # now make all enclosed in ''
    tt <- sub('.*=', '', splitmvar_w2)
    tt <- paste0("'", tt, "'")
    
    x <- parms_ <- paste0(splitmvar_w3, "=", tt)
    
    # not needed, directly using name_parameter towards the end
    collect_name_parameter <- c()
    
    for (i in 1:length(x)) {
      pname_ <- substr(x[i], 1, regexpr("\\=", x[i]) - 1)
      x_i <- gsub("\"" , "", x[i])
      x_i <- eval(parse(text = x_i))
      
      sethp <-
        isTRUE(ept(gsub("\"" , "", splitmvar_w2[grepl("sethp", splitmvar_w2)])))
      
      if (ept(sethp) & (dist != "normal" &
                        dist != "cauchy" &
                        dist != "student_nu" &
                        dist != "student_t")) {
        stop(
          "Hierarchical priors are supported only for normal, cauchy, student_nu",
          "\n",
          " and student_t distributins.",
          "\n",
          " Please set sethp = NA or sethp = FALSE for the '",
          dist,
          "\n",
          "'distributins specified for nlpar '",
          nlpar,
          "', class '",
          class,
          "'"
        )
      }
      
      
      # Get scale_factor to multiply with scale parameters
      
      err. <- FALSE
      tryCatch(
        expr = {
          check_for_autoscale <-
            ept(gsub("\"" , "", splitmvar_w2[grepl("autoscale", splitmvar_w2)]))
          ! is.na(check_for_autoscale)
          ! is.logical(check_for_autoscale)
          ! is.numeric(check_for_autoscale)
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (err.) {
        stop("scale factor set by autoscale can only be",
             "\n",
             " NA / TRUE / FLASE or a numeric value")
      } else {
        check_for_autoscale <-
          ept(gsub("\"" , "", splitmvar_w2[grepl("autoscale", splitmvar_w2)]))
        if (!is.na(check_for_autoscale) &
            !is.logical(check_for_autoscale) &
            !is.numeric(check_for_autoscale)) {
          stop("scale factor set by autoscale can only be",
               "\n",
               "NA / TRUE / FLASE or a numeric value")
        }
      }
      
      if (is.na(check_for_autoscale) | !(check_for_autoscale)) {
        scale_factor <- 1
      } else if (check_for_autoscale &
                 !is.numeric(check_for_autoscale)) {
        scale_factor <- 2.5
        if (verbose)
          message("scale factor for autoscale option set to 2.5")
      } else if (is.numeric(check_for_autoscale)) {
        scale_factor <- check_for_autoscale
      }
      
      # Get addrange to add to uniform prior range
      
      if (is.na(ept(gsub("\"" , "", 
                         splitmvar_w2[grepl("addrange", splitmvar_w2)])))) {
        addrange <- 0
      } else {
        addrange <-
          ept(gsub("\"" , "", splitmvar_w2[grepl("addrange", splitmvar_w2)]))
      }
      
      # assigning required parameters
      
      if (pname_ %in% set_str_names) {
        check_evalation_of_numeric_pdata_obj <-
          function(prior_argument,
                   p_str_in,
                   eit,
                   x,
                   pname_,
                   dist,
                   nlpar,
                   class,
                   allowed_parm_options,
                   splitmvar_w2) {
            whatin <-
              sub("=[^=]+$", "", splitmvar_w2[grepl(eit, splitmvar_w2)])
            const_msg <- 
              paste0(" - a numeric value (e.g., 2) or a charater string such as",
                     "\n",
                     "xxx with xxx defined in the use-specified 'prior_data'",
                     "\n",
                     "argument e.g., prior_data = list(xxx = 2)"
              )
            
            if (!is.null(allowed_parm_options)) {
              allowed_parm_options <-
                paste0(allowed_parm_options, collapse = ", ")
              allowed_parm_options <-
                paste0(" - ", allowed_parm_options)
              const_msg <- paste0(allowed_parm_options, "\n", const_msg)
            }
            err. <- FALSE
            tryCatch(
              expr = {
                out <- ept(eit)
              },
              error = function(e) {
                err. <<- TRUE
              }
            )
            if (err.) {
              if (class == 'b' | class == 'sd') {
                stop(
                  "\nFor nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  ", you have specified '",
                  eit,
                  "' as ",
                  pname_,
                  " for the ",
                  dist,
                  " distribution",
                  "\n" ,
                  " But '",
                  eit,
                  "' is not found in the 'prior_data_internal'",
                  "\n" ,
                  " or use-specified 'prior_data' argument",
                  "\n ",
                  " [see specified prior argument: ",
                  prior_argument,
                  " = ",
                  p_str_in,
                  "]",
                  "\n" ,
                  "Avilable ",
                  pname_,
                  " parameter options are:" ,
                  "\n" ,
                  const_msg
                )
              } else if (class == 'sigma') {
                stop(
                  "\nFor residual standard deviation parameter i.e., ",
                  "class ",
                  class,
                  ", you have specified '",
                  eit,
                  "' as ",
                  pname_,
                  " for the ",
                  dist,
                  " distribution",
                  "\n" ,
                  " But '",
                  eit,
                  "' is not found in the 'prior_data_internal'",
                  "\n" ,
                  " or use-specified 'prior_data' argument",
                  "\n ",
                  " [see specified prior argument: ",
                  prior_argument,
                  " = ",
                  p_str_in,
                  "]",
                  "\n" ,
                  "Avilable ",
                  pname_,
                  " parameter options are:" ,
                  "\n" ,
                  const_msg
                )
              } else if (class == '' &
                         grepl("dpar_", prior_argument) &
                         !grepl("dpar_cov", prior_argument)) {
                stop(
                  "\nFor for distributional Intercept parameter i.e., ",
                  "Intercept_sigma ",
                  ", you have specified '",
                  eit,
                  "' as ",
                  pname_,
                  " for the ",
                  dist,
                  " distribution",
                  "\n" ,
                  " But '",
                  eit,
                  "' is not found in the 'prior_data_internal'",
                  "\n" ,
                  " or use-specified 'prior_data' argument",
                  "\n ",
                  " [see specified prior argument: ",
                  prior_argument,
                  " = ",
                  p_str_in,
                  "]",
                  "\n" ,
                  "Avilable ",
                  pname_,
                  " parameter options are:" ,
                  "\n" ,
                  const_msg
                )
              }
              
            }
          }
        
        
        
        
        allowed_parm_options <- NULL
        
        if (grepl("^location$", pname_)) {
          if (class == "b") {
            if (nlpar == "a" & cov_nlpar == "") {
              allowed_parm_options <- c("lm", "ymean", 
                                        "ymedian", "ymax",
                                        "ymeanxmax", "ymeanxmin")
            } else if (nlpar == "b") {
              allowed_parm_options <- c("lm", "ymaxs", "bstart")
            } else if (nlpar == "c") {
              allowed_parm_options <- c("lm", "cstart")
            } else if (nlpar == "d") {
              allowed_parm_options <- c("lm", "dstart", "ymeanxmid")
            } else if (nlpar == "e") {
              allowed_parm_options <- c("lm", "estart")
            } else if (nlpar == "f") {
              allowed_parm_options <- NULL
            } else if (nlpar == "g") {
              allowed_parm_options <- c("ymeanxmax", "ymeanxmidxmaxdiff")
            } else if (nlpar == "h") {
              allowed_parm_options <- NULL
            } else if (nlpar == "i") {
              allowed_parm_options <- c("lm", "estart", "bstart")
            } else if (nlpar == "s") {
              allowed_parm_options <- c("lm")
            } else if (cov_nlpar != "") {
              allowed_parm_options <- c("lm")
            } else {
              allowed_parm_options <- NULL
            }
          }
          
          
          if (class == "sd") {
            allowed_parm_options <- NULL
          }
          if (class == "cor") {
            allowed_parm_options <- NULL
          }
          if (class == "sigma") {
            allowed_parm_options <- NULL
          }
          if (class == '' &
              grepl("dpar_", prior_argument) &
              !grepl("dpar_cov", prior_argument)) {
            allowed_parm_options <- NULL
          }
          if (class == "") {
            allowed_parm_options <- NULL
          }
          if (dpar != "") {
            allowed_parm_options <- NULL
          }
          if (cov_dpar != "") {
            allowed_parm_options <- NULL
          }
          
          if (!is.null(allowed_parm_options)) {
            allowed_init_options_beta <- allowed_parm_options
          } else {
            allowed_init_options_beta <- NULL
          }
        }
        allowed_parm_options <- c("ysd", "ymad")
        ##########
        if (grepl("^scale$", pname_) &
            dist != "gamma" & dist != "inv_gamma") {
          if (class == "b" | class == "sd") {
            if (nlpar == "a" & cov_nlpar == "") {
              allowed_parm_options <- c("ysd", "ymad", "lme_sd_a",
                                        "ysdxmax", "ysdxmin")
            } else if (nlpar == "b") {
              allowed_parm_options <- c("lm", "ymaxs", "bstart")
            } else if (nlpar == "c") {
              allowed_parm_options <- c("lm", "cstart")
            } else if (nlpar == "d") {
              allowed_parm_options <- c("ysd", "ymad", "ysdxmid")
            } else if (nlpar == "e") {
              allowed_parm_options <- c("lm", "estart")
            } else if (nlpar == "f") {
              allowed_parm_options <- NULL
            } else if (nlpar == "g") {
              allowed_parm_options <- c("ysdxmax", "ysdxmidxmaxdiff")
            } else if (nlpar == "h") {
              allowed_parm_options <- NULL
            } else if (nlpar == "i") {
              allowed_parm_options <- NULL
            } else if (nlpar == "s") {
              allowed_parm_options <- c("sdx")
            } else if (cov_nlpar != "") {
              allowed_parm_options <- c("lm")
            } else {
              allowed_parm_options <- NULL
            }
          }
          
          
          
          
          if (class == "cor") {
            allowed_parm_options <- NULL
          }
          
          if (class == "sigma") {
            allowed_parm_options <- c("ysd", "ymad", "lme_rsd", "lm_rsd")
          }
          
          if (class == "") {
            allowed_parm_options <- NULL
          }
          
          if (dpar != "") {
            allowed_parm_options <- NULL
          }
          
          if (cov_dpar != "") {
            allowed_parm_options <- NULL
          }
          
          if (class == '' &
              grepl("dpar_", prior_argument) &
              !grepl("dpar_cov", prior_argument)) {
            allowed_parm_options <- c("ysd", "ymad", "lme_rsd", "lm_rsd")
          }
          
          if (!is.null(allowed_parm_options)) {
            allowed_init_options_sd <- allowed_parm_options
          } else {
            allowed_init_options_sd <- NULL
          }
        }
        
        
        if (grepl("^scale$", pname_) &
            dist == "gamma" | dist == "inv_gamma") {
          allowed_parm_options <- NULL
        }
        
        if (grepl("^shape$", pname_)) {
          allowed_parm_options <- NULL
        }
        
        if (grepl("^rate$", pname_)) {
          if (class == "sigma") {
            allowed_parm_options <-
              c("ysd", "ymad", "lme_sd_a", "lme_rsd", "lm_rsd")
          } else {
            allowed_parm_options <- NULL
          }
          
          if (!is.null(allowed_parm_options)) {
            allowed_init_options_rate <- allowed_parm_options
          } else {
            allowed_init_options_rate <- NULL
          }
          
        } # if (grepl("^rate$", pname_)) {
        
        
        if (grepl("^shape$", pname_)) {
          allowed_init_options_shape <- NULL # 22 4 2023
          allowed_init_options_scale <- NULL # 22 4 2023
        }
        
        if (grepl("^scale$", pname_)) {
          allowed_init_options_shape <- NULL # 22 4 2023
          allowed_init_options_scale <- NULL # 22 4 2023
        }
        
        
        if (grepl("^df$", pname_)) {
          allowed_parm_options <- NULL
        }
        
        if (grepl("^eta$", pname_)) {
          allowed_parm_options <- NULL
        }
        
        
        if (!exists('allowed_init_options_beta'))
          allowed_init_options_beta <- NULL
        if (!exists('allowed_init_options_sd'))
          allowed_init_options_sd <- NULL
        if (!exists('allowed_init_options_rate'))
          allowed_init_options_rate <- NULL
        
        if (!exists('allowed_init_options_shape'))
          allowed_init_options_shape <- NULL
        if (!exists('allowed_init_options_scale'))
          allowed_init_options_scale <- NULL
        
        
        # set location parameter -> for normal, log normal, cauchy, studdent_t
        
        if (grepl("^location$", pname_)) {
          # location nlpar a (class b)
          if (nlpar == "a" & class == "b" & grepl("a", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (a_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymean", empty_sufx)) {
              eit <-  gsub("ymean", paste0("ymean", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymax", empty_sufx)) {
              eit <-  gsub("ymax", paste0("ymax", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymaxs", empty_sufx)) {
              eit <-  gsub("ymaxs", paste0("ymaxs", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              eit <-  gsub("ymedian", paste0("ymedian", resp_), x_i)
              evaluated_parameter <- ept(eit) 
            } else if (x_i == paste0("ymeanxmin", empty_sufx)) {
              eit <-  gsub("ymeanxmin", paste0("ymeanxmin", resp_), x_i)
              evaluated_parameter <- ept(eit) 
            } else if (x_i == paste0("ymeanxmax", empty_sufx)) {
              eit <-  gsub("ymeanxmax", paste0("ymeanxmax", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, ymax, ymeanxmin, ymeanxmax,
                a numeric value (e.g., 2)",
                  "\n",
                  " or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "\n", 
                  " e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar b (class b)
          if (nlpar == "b" & class == "b" & grepl("b", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (b_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
              if (verbose)
                message("location parameter specified as lm for nlpar",
                        nlpar,
                        " is set as 0")
            } else if (x_i == paste0("ymax", empty_sufx)) {
              eit <-  gsub("ymax", paste0("ymax", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymaxs", empty_sufx)) {
              eit <-  gsub("ymaxs", paste0("ymaxs", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymean", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              stop("option ymedian as location parameter not alloweed for nlpar ",
                   nlpar)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the",
                  "prior_data e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar c (class b)
          if (nlpar == "c" & class == "b" & grepl("c", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (c_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
              if (verbose)
                message("location parameter specified as lm for nlpar",
                        nlpar,
                        " is set as 0")
            } else if (x_i == paste0("ymean", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("cstart", empty_sufx)) {
              eit <-  gsub("cstart", paste0("cstart", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar d (class b)
          if (nlpar == "d" & class == "b" & grepl("d", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (d_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
              if (verbose)
                message("location parameter specified as lm for nlpar",
                        nlpar,
                        " is set as 0")
            } else if (x_i == paste0("ymeanxmid", empty_sufx)) {
              eit <-  gsub("ymeanxmid", paste0("ymeanxmid", resp_), x_i)
              evaluated_parameter <- ept(eit) 
            } else if (x_i == paste0("ymeanxmid", empty_sufx)) {
              eit <-  gsub("ymeanxmid", paste0("ymeanxmid", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymean", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("dstart", empty_sufx)) {
              eit <-  gsub("dstart", paste0("dstart", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar e (class b)
          if (nlpar == "e" & class == "b" & grepl("e", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (e_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
              if (verbose)
                message("location parameter specified as lm for nlpar",
                        nlpar,
                        " is set as 0")
            } else if (x_i == paste0("ymean", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar) 
            } else if (x_i == paste0("estart", empty_sufx)) {
              eit <-  gsub("estart", paste0("estart", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar f (class b)
          if (nlpar == "f" & class == "b" & grepl("f", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (f_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
              if (verbose)
                message("location parameter specified as lm for nlpar",
                        nlpar,
                        " is set as 0")
            } else if (x_i == paste0("ymean", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar) 
            } else if (x_i == paste0("estart", empty_sufx)) {
              eit <-  gsub("estart", paste0("estart", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar g (class b)
          if (nlpar == "g" & class == "b" & grepl("g", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (g_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
              if (verbose)
                message("location parameter specified as lm for nlpar",
                        nlpar,
                        " is set as 0")
            } else if (x_i == paste0("ymeanxmax", empty_sufx)) {
              eit <-  gsub("ymeanxmax", paste0("ymeanxmax", resp_), x_i)
              evaluated_parameter <- ept(eit) 
            } else if (x_i == paste0("ymeanxmidxmaxdiff", empty_sufx)) {
              eit <-  gsub("ymeanxmidxmaxdiff", 
                           paste0("ymeanxmidxmaxdiff", resp_), x_i)
              evaluated_parameter <- ept(eit) 
            } else if (x_i == paste0("ymax", empty_sufx)) {
              eit <-  gsub("ymax", paste0("ymax", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymean", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar) 
            } else if (x_i == paste0("estart", empty_sufx)) {
              eit <-  gsub("estart", paste0("estart", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar h (class b)
          if (nlpar == "h" & class == "b" & grepl("h", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (h_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
              if (verbose)
                message("location parameter specified as lm for nlpar",
                        nlpar,
                        " is set as 0")
            } else if (x_i == paste0("ymean", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar) 
            } else if (x_i == paste0("estart", empty_sufx)) {
              eit <-  gsub("estart", paste0("estart", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location nlpar i (class b)
          if (nlpar == "i" & class == "b" & grepl("i", fixedsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (i_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
              if (verbose)
                message("location parameter specified as lm for nlpar",
                        nlpar,
                        " is set as 0")
            } else if (x_i == paste0("ymean", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              stop("option ymean as location parameter not alloweed for nlpar ",
                   nlpar) 
            } else if (x_i == paste0("estart", empty_sufx)) {
              eit <-  gsub("estart", paste0("estart", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar s (class b) - sitar
          if (nlpar == "s" & class == "b") {
            if (x_i == paste0("lm", empty_sufx)) {
              if (s_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2) or",
                  "a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            # checks
            if (nlpar == "s" & !is.null(sncov)) {
              if (length(evaluated_parameter) == 1) {
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              } else if (length(evaluated_parameter) == df) {
                repeach <- nrep_of_parms / df
                evaluated_parameter <-
                  rep(
                    evaluated_parameter,
                    times = 1,
                    each = repeach,
                    length.out = nrep_of_parms
                  )
              }
            } else {
              if (length(evaluated_parameter) < nrep_of_parms)
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              if (length(evaluated_parameter) > nrep_of_parms)
                stop(
                  "prior elements for nlpar ",
                  nlpar, ", class ",  class,
                  " are greater than the parameter dimensions"
                )
            }
          }
          
          
          
          # location nlpar a cov (class b)
          if (cov_nlpar == "a" & class == "b" & !is.null(ancov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!a_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, ymax a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar b cov (class b)
          if (cov_nlpar == "b" & class == "b" & !is.null(bncov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!b_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, ymaxs a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar c cov (class b)
          if (cov_nlpar == "c" & class == "b" & !is.null(cncov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!c_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar d cov (class b)
          if (cov_nlpar == "d" & class == "b" & !is.null(dncov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!d_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar e cov (class b)
          if (cov_nlpar == "e" & class == "b" & !is.null(encov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!e_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location nlpar f cov (class b)
          if (cov_nlpar == "f" & class == "b" & !is.null(fncov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!f_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar g cov (class b)
          if (cov_nlpar == "g" & class == "b" & !is.null(gncov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!g_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar h cov (class b)
          if (cov_nlpar == "h" & class == "b" & !is.null(hncov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!h_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location nlpar i cov (class b)
          if (cov_nlpar == "i" & class == "b" & !is.null(incov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!i_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar s cov (class b)
          if (cov_nlpar == "s" & !is.null(sncov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!s_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            # checks
            if (nlpar == "s" & !is.null(sncov)) {
              if (length(evaluated_parameter) == 1) {
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              } else if (length(evaluated_parameter) == df) {
                repeach <- nrep_of_parms / df
                evaluated_parameter <-
                  rep(evaluated_parameter,
                      times = 1,
                      each = repeach)
              } else {
                #
              }
            } else {
              if (length(evaluated_parameter) < nrep_of_parms)
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              if (length(evaluated_parameter) > nrep_of_parms)
                stop("prior elements for nlpar ",
                     nlpar, ", class ",  class,
                     " are greater than the parameter dimensions"
                )
            }
          }
          
          
          
          
          # location sigma (class b)
          
          if (nlpar == "" & class == "b" & sigma_dpar == 'sigma') {
            if (x_i == paste0("lm", empty_sufx)) {
              if (sigma_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymean", empty_sufx)) {
              eit <-  gsub("ymean", paste0("ymean", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else if (x_i == paste0("ymedian", empty_sufx)) {
              eit <-  gsub("ymedian", paste0("ymedian", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for distributional ",
                  sigma_dpar,
                  ", class ",
                  class,
                  " are:a numeric value (e.g., 2)",
                  "\n",
                  " or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "\n", 
                  " e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for distributional ",
                   sigma_dpar,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location sigma cov (class b)
          if (cov_sigma_dpar != "" & class == "b" & cov_sigma_dpar == 'sigma_cov' & 
              !is.null(sigmancov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!sigma_form_0) {
                # lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                # lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              # eit <-  gsub("lm", lm_gsubby, x_i)
              # evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for distributional ",
                  sigma_dpar,
                  ", class ",
                  class,
                  " are:\n a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for distributional ",
                   sigma_dpar,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location a b c d e random effects
          
          # location nlpar a (class sd, typically 0)
          
          if (nlpar == "a" & class == "sd" & grepl("a", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar b (class sd, typically 0)
          if (nlpar == "b" & class == "sd" & grepl("b", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          # location nlpar c (class sd, typically 0)
          if (nlpar == "c" & class == "sd" & grepl("c", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          # location nlpar d (class sd, typically 0)
          if (nlpar == "d" & class == "sd" & grepl("d", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar e (class sd, typically 0)
          if (nlpar == "e" & class == "sd" & grepl("e", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar f (class sd, typically 0)
          if (nlpar == "f" & class == "sd" & grepl("f", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar g (class sd, typically 0)
          if (nlpar == "g" & class == "sd" & grepl("g", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar h (class sd, typically 0)
          if (nlpar == "h" & class == "sd" & grepl("h", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location nlpar i (class sd, typically 0)
          if (nlpar == "i" & class == "sd" & grepl("i", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          
          # location nlpar s (class sd, typically 0)
          if (nlpar == "s" & class == "sd") {
            if (x_i == paste0("lm", empty_sufx)) {
              if (s_form_0_gr) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2) or",
                  "a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            # checks
            if (nlpar == "s" & !is.null(sncov_gr)) {
              if (length(evaluated_parameter) == 1) {
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              } else if (length(evaluated_parameter) == df) {
                repeach <- nrep_of_parms / df
                evaluated_parameter <-
                  rep(
                    evaluated_parameter,
                    times = 1,
                    each = repeach,
                    length.out = nrep_of_parms
                  )
              }
            } else {
              if (length(evaluated_parameter) < nrep_of_parms)
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              if (length(evaluated_parameter) > nrep_of_parms)
                stop(
                  "prior elements for nlpar ",
                  nlpar, ", class ",  class,
                  " are greater than the parameter dimensions"
                )
            }
          }
          
          
          
          
          # location nlpar a cov (class sd, typically 0)
          if (cov_nlpar == "a" & class == "sd" & !is.null(ancov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar b cov (class sd, typically 0)
          if (cov_nlpar == "b" & class == "sd" & !is.null(bncov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar c cov (class sd, typically 0)
          if (cov_nlpar == "c" & class == "sd" & !is.null(cncov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # location nlpar d cov (class sd, typically 0)
          if (cov_nlpar == "d" & class == "sd" & !is.null(dncov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location nlpar e cov (class sd, typically 0)
          if (cov_nlpar == "e" & class == "sd" & !is.null(encov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location nlpar f cov (class sd, typically 0)
          if (cov_nlpar == "f" & class == "sd" & !is.null(fncov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location nlpar g cov (class sd, typically 0)
          if (cov_nlpar == "g" & class == "sd" & !is.null(gncov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar h cov (class sd, typically 0)
          if (cov_nlpar == "h" & class == "sd" & !is.null(hncov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location nlpar i cov (class sd, typically 0)
          if (cov_nlpar == "i" & class == "sd" & !is.null(incov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          
          # location nlpar i cov (class sd, typically 0)
          if (cov_nlpar == "s" & !is.null(sncov_gr)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!s_form_0_gr) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            # checks
            if (nlpar == "s" & !is.null(sncov_gr)) {
              if (length(evaluated_parameter) == 1) {
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              } else if (length(evaluated_parameter) == df) {
                repeach <- nrep_of_parms / df
                evaluated_parameter <-
                  rep(evaluated_parameter,
                      times = 1,
                      each = repeach)
              } else {
                #
              }
            } else {
              if (length(evaluated_parameter) < nrep_of_parms)
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              if (length(evaluated_parameter) > nrep_of_parms)
                stop("prior elements for nlpar ",
                     nlpar, ", class ",  class,
                     " are greater than the parameter dimensions"
                )
            }
          }
          
          
          
          
          # location sigma (class sd, typically 0)
          if (nlpar == "" & class == "sd" & sigma_dpar == 'sigma') {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for distributional ",
                sigma_dpar,
                ", class ",
                class,
                " are:\n a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for distributional ",
                   sigma_dpar,
                   " are greater than the parameter dimensions")
          }
          
          
          # location sigma cov (class sd, typically 0)
          if (cov_sigma_dpar != "" & class == "sd" & cov_sigma_dpar == 'sigma_cov' & 
              !is.null(sigmancov_gr)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for distributional ",
                sigma_dpar,
                ", class ",
                class,
                " are:\n a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for distributional ",
                   sigma_dpar,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # location rsd param
          if (class == "sigma") {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # location dpar param ~
          if (!is.null(dparncov) & class == "") {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "location parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                "or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
        } # if(grepl("^location$", i))
        
        
        
        
        
        # set scale parameter -> for normal, log normal, cauchy, studdent_t
        
        if (grepl("^scale$", pname_)) {
          # scale a b c d e fixed effects
          
          # scale nlpar a (class b)
          if (nlpar == "a" & class == "b" & grepl("a", fixedsi)) {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("lme_sd_a", empty_sufx)) {
              eit <-  gsub("lme_sd_a", paste0("lme_sd_a", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ysdxmin", empty_sufx)) {
              eit <-  gsub("ysdxmin", paste0("ysdxmin", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit) 
            } else if (x_i == paste0("ysdxmin", empty_sufx)) {
              eit <-  gsub("ysdxmin", paste0("ysdxmin", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n ysd, ysd, lme_sd_a, or a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale nlpar b (class b)
          if (nlpar == "b" & class == "b" & grepl("b", fixedsi)) {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("lme_sd_a", empty_sufx)) {
              eit <-  gsub("lme_sd_a", paste0("lme_sd_a", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n ysd, ysd, lme_sd_a, or a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale nlpar c (class b)
          if (nlpar == "c" & class == "b" & grepl("c", fixedsi)) {
            if (x_i == paste0("vsd", empty_sufx)) {
              eit <-  gsub("vsd", paste0("vsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("vmad", empty_sufx)) {
              eit <-  gsub("vmad", paste0("vmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data", 
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale nlpar d (class b)
          if (nlpar == "d" & class == "b" & grepl("d", fixedsi)) {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("dsd", empty_sufx)) {
              eit <-  gsub("dsd", paste0("vsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("dmad", empty_sufx)) {
              eit <-  gsub("dmad", paste0("dmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("ysdxmid", empty_sufx)) {
              eit <-  gsub("ysdxmid", paste0("ysdxmid", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit) 
            } else if (x_i == paste0("ysdxmid", empty_sufx)) {
              eit <-  gsub("ysdxmid", paste0("ysdxmid", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale nlpar e (class b)
          if (nlpar == "e" & class == "b" & grepl("e", fixedsi)) {
            if (x_i == paste0("xsd", empty_sufx)) {
              eit <-  gsub("xsd", paste0("xsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("xmad", empty_sufx)) {
              eit <-  gsub("xmad", paste0("xmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # scale nlpar f (class b)
          if (nlpar == "f" & class == "b" & grepl("f", fixedsi)) {
            if (x_i == paste0("xsd", empty_sufx)) {
              eit <-  gsub("xsd", paste0("xsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("xmad", empty_sufx)) {
              eit <-  gsub("xmad", paste0("xmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar g (class b)
          if (nlpar == "g" & class == "b" & grepl("g", fixedsi)) {
            if (x_i == paste0("xsd", empty_sufx)) {
              eit <-  gsub("xsd", paste0("xsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("xmad", empty_sufx)) {
              eit <-  gsub("xmad", paste0("xmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("ysdxmax", empty_sufx)) {
              eit <-  gsub("ysdxmax", paste0("ysdxmax", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit) 
            } else if (x_i == paste0("ysdxmidxmaxdiff", empty_sufx)) {
              eit <-  gsub("ysdxmidxmaxdiff", 
                           paste0("ysdxmidxmaxdiff", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar h (class b)
          if (nlpar == "h" & class == "b" & grepl("h", fixedsi)) {
            if (x_i == paste0("xsd", empty_sufx)) {
              eit <-  gsub("xsd", paste0("xsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("xmad", empty_sufx)) {
              eit <-  gsub("xmad", paste0("xmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar i (class b)
          if (nlpar == "i" & class == "b" & grepl("i", fixedsi)) {
            if (x_i == paste0("xsd", empty_sufx)) {
              eit <-  gsub("xsd", paste0("xsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("xmad", empty_sufx)) {
              eit <-  gsub("xmad", paste0("xmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar s (class b) - sitar
          if (nlpar == "s" & class == "b") {
            if (x_i == paste0("lm", empty_sufx)) {
              if (s_form_0) {
                lm_gsubby <- paste0("lm", "_", 'sdx', "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", 'sdx', "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("stau", empty_sufx)) {
              evaluated_parameter <- rep(NA, nrep_of_parms)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data", 
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            # checks
            if (nlpar == "s" & !is.null(sncov)) {
              if (length(evaluated_parameter) == 1) {
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              } else if (length(evaluated_parameter) == df) {
                repeach <- nrep_of_parms / df
                evaluated_parameter <-
                  rep(evaluated_parameter,
                      times = 1,
                      each = repeach)
              }
            } else {
              if (length(evaluated_parameter) < nrep_of_parms)
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              if (length(evaluated_parameter) > nrep_of_parms)
                stop("prior elements for nlpar ",
                     nlpar, ", class ",  class,
                     " are greater than the parameter dimensions"
                )
            }
          }
          
          
          
          # scale nlpar a cov (class b)
          if (cov_nlpar == "a" & class == "b" & !is.null(ancov)) {
            if (x_i == paste0("sdacov", empty_sufx)) {
              eit <-  gsub("sdacov", paste0("acov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale nlpar b cov (class b)
          if (cov_nlpar == "b" & class == "b" & !is.null(bncov)) {
            if (x_i == paste0("sdbcov", empty_sufx)) {
              eit <-  gsub("sdbcov", paste0("bcov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          } # if (cov_nlpar == "b" & class == "b" & !is.null(bncov)) {
          
          
          
          # scale nlpar c cov (class b)
          if (cov_nlpar == "c" & class == "b" & !is.null(cncov)) {
            if (x_i == paste0("sdccov", empty_sufx)) {
              eit <-  gsub("sdccov", paste0("ccov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar d cov (class b)
          if (cov_nlpar == "d" & class == "b" & !is.null(dncov)) {
            if (x_i == paste0("sddcov", empty_sufx)) {
              eit <-  gsub("sddcov", paste0("dcov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar e cov (class b)
          if (cov_nlpar == "e" & class == "b" & !is.null(encov)) {
            if (x_i == paste0("sdecov", empty_sufx)) {
              eit <-  gsub("sdecov", paste0("ecov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # scale nlpar f cov (class b)
          if (cov_nlpar == "f" & class == "b" & !is.null(fncov)) {
            if (x_i == paste0("sdfcov", empty_sufx)) {
              eit <-  gsub("sdfcov", paste0("fcov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # scale nlpar g cov (class b)
          if (cov_nlpar == "g" & class == "b" & !is.null(gncov)) {
            if (x_i == paste0("sdfcov", empty_sufx)) {
              eit <-  gsub("sdfcov", paste0("gcov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar h cov (class b)
          if (cov_nlpar == "h" & class == "b" & !is.null(hncov)) {
            if (x_i == paste0("sdfcov", empty_sufx)) {
              eit <-  gsub("sdfcov", paste0("hcov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar i cov (class b)
          if (cov_nlpar == "f" & class == "b" & !is.null(incov)) {
            if (x_i == paste0("sdfcov", empty_sufx)) {
              eit <-  gsub("sdfcov", paste0("icov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale nlpar s cov (class b) - sitar
          if (cov_nlpar == "s" & class == "b" & !is.null(sncov)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!s_form_0) {
                lm_gsubby <- paste0("lm", "_", 'sdx', "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", 'sdx', "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              }
            }
            # checks
            if (nlpar == "s" & !is.null(sncov)) {
              if (length(evaluated_parameter) == 1) {
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              } else if (length(evaluated_parameter) == df) {
                repeach <- nrep_of_parms / df
                evaluated_parameter <-
                  rep(evaluated_parameter,
                      times = 1,
                      each = repeach)
              } else {
                #
              }
            } else {
              if (length(evaluated_parameter) < nrep_of_parms)
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              if (length(evaluated_parameter) > nrep_of_parms)
                stop("prior elements for nlpar ",
                     nlpar, ", class ",  class,
                     " are greater than the parameter dimensions"
                )
            }
          }
          
          
          
          
          
          
          # scale sigma (class b)
          if (nlpar == "" & class == "b" & sigma_dpar == "sigma") {
            if (x_i == paste0("vsd", empty_sufx)) {
              eit <-  gsub("vsd", paste0("vsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("vmad", empty_sufx)) {
              eit <-  gsub("vmad", paste0("vmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for distributional ",
                  sigma_dpar,
                  ", class ",
                  class,
                  " are:\n a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data", 
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for distributional ",
                   sigma_dpar,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale sigma cov (class b)
          if (cov_sigma_dpar != "" & class == "b" & cov_sigma_dpar == 'sigma_cov' & 
              !is.null(sigmancov)) {
            if (x_i == paste0("sdacov", empty_sufx)) {
              eit <-  gsub("sdacov", paste0("acov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for distributional ",
                  sigma_dpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements fordistributional ",
                   sigma_dpar,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          
          
          
          
          # scale a b c d e f random effects
          
          # scale a (class sd)
          if (nlpar == "a" & class == "sd" & grepl("a", randomsi)) {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("lme_sd_a", empty_sufx)) {
              eit <-  gsub("lme_sd_a", paste0("lme_sd_a", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ysdxmax", empty_sufx)) {
              eit <-  gsub("ysdxmax", paste0("ysdxmax", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit) 
            } else if (x_i == paste0("ysdxmin", empty_sufx)) {
              eit <-  gsub("ysdxmin", paste0("ysdxmin", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n ysd, ysd, lme_sd_a, or a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale b (class sd)
          if (nlpar == "b" & class == "sd" & grepl("b", randomsi)) {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("lme_sd_a", empty_sufx)) {
              eit <-  gsub("lme_sd_a", paste0("lme_sd_a", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n ysd, ysd, lme_sd_a, or a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale c (class sd)
          if (nlpar == "c" & class == "sd" & grepl("c", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "scale parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale d (class sd)
          if (nlpar == "d" & class == "sd" & grepl("d", randomsi)) {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ysdxmid", empty_sufx)) {
              eit <-  gsub("ysdxmid", paste0("ysdxmid", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit) 
            } else if (x_i == paste0("ysdxmidxmaxdiff", empty_sufx)) {
              eit <-  gsub("ysdxmidxmaxdiff", 
                           paste0("ysdxmidxmaxdiff", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit) 
            } else if (x_i == paste0("ysdxmax", empty_sufx)) {
              eit <-  gsub("ysdxmax", paste0("ysdxmax", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n ysd, ysd, lme_sd_a, or a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale e (class sd)
          if (nlpar == "e" & class == "sd" & grepl("e", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "scale parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # scale f (class sd)
          if (nlpar == "f" & class == "sd" & grepl("f", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "scale parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale g (class sd)
          if (nlpar == "g" & class == "sd" & grepl("g", randomsi)) {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ysdxmid", empty_sufx)) {
              eit <-  gsub("ysdxmid", paste0("ysdxmid", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit) 
            } else if (x_i == paste0("ysdxmidxmaxdiff", empty_sufx)) {
              eit <-  gsub("ysdxmidxmaxdiff", 
                           paste0("ysdxmidxmaxdiff", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit) 
            } else if (x_i == paste0("ysdxmax", empty_sufx)) {
              eit <-  gsub("ysdxmax", paste0("ysdxmax", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n ysd, ysd, lme_sd_a, or a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale h (class sd)
          if (nlpar == "h" & class == "sd" & grepl("h", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "scale parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale i (class sd)
          if (nlpar == "i" & class == "sd" & grepl("i", randomsi)) {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "scale parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale nlpar s (class sd) - sitar
          if (nlpar == "s" & class == "sd" & grepl("s", randomsi)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (s_form_0_gr) {
                lm_gsubby <- paste0("lm", "_", 'sdx', "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", 'sdx', "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("stau", empty_sufx)) {
              evaluated_parameter <- rep(NA, nrep_of_parms)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "location parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data", 
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            # checks
            if (nlpar == "s" & !is.null(sncov_gr)) {
              if (length(evaluated_parameter) == 1) {
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              } else if (length(evaluated_parameter) == df) {
                repeach <- nrep_of_parms / df
                evaluated_parameter <-
                  rep(evaluated_parameter,
                      times = 1,
                      each = repeach)
              }
            } else {
              if (length(evaluated_parameter) < nrep_of_parms)
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              if (length(evaluated_parameter) > nrep_of_parms)
                stop("prior elements for nlpar ",
                     nlpar, ", class ",  class,
                     " are greater than the parameter dimensions"
                )
            }
          }
          
          
          
          
          # scale a cov (class sd)
          if (cov_nlpar == "a" & class == "sd" & !is.null(ancov_gr)) {
            if (x_i == paste0("sdacov", empty_sufx)) {
              eit <-  gsub("sdacov", paste0("acov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale b cov (class sd)
          if (cov_nlpar == "b" & class == "sd" & !is.null(bncov_gr)) {
            if (x_i == paste0("sdbcov", empty_sufx)) {
              eit <-  gsub("sdbcov", paste0("bcov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale c cov (class sd)
          if (cov_nlpar == "c" & class == "sd" & !is.null(cncov_gr)) {
            if (x_i == paste0("sdccov", empty_sufx)) {
              eit <-  gsub("sdccov", paste0("ccov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale d cov (class sd)
          if (cov_nlpar == "c" & class == "sd" & !is.null(dncov_gr)) {
            if (x_i == paste0("sdccov", empty_sufx)) {
              eit <-  gsub("sdccov", paste0("ccov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # scale e cov (class sd)
          if (cov_nlpar == "e" & class == "sd" & !is.null(encov_gr)) {
            if (x_i == paste0("sdecov", empty_sufx)) {
              eit <-  gsub("sdecov", paste0("ecov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # scale f cov (class sd)
          if (cov_nlpar == "f" & class == "sd" & !is.null(fncov_gr)) {
            if (x_i == paste0("sdfcov", empty_sufx)) {
              eit <-  gsub("sdfcov", paste0("fcov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale g cov (class sd)
          if (cov_nlpar == "g" & class == "sd" & !is.null(fncov_gr)) {
            if (x_i == paste0("sdfcov", empty_sufx)) {
              eit <-  gsub("sdfcov", paste0("gcov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # scale h cov (class sd)
          if (cov_nlpar == "h" & class == "sd" & !is.null(hncov_gr)) {
            if (x_i == paste0("sdfcov", empty_sufx)) {
              eit <-  gsub("sdfcov", paste0("hcov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          # scale i cov (class sd)
          if (cov_nlpar == "i" & class == "sd" & !is.null(incov_gr)) {
            if (x_i == paste0("sdfcov", empty_sufx)) {
              eit <-  gsub("sdfcov", paste0("icov_sd_gr", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          
          # scale nlpar s cov (class sd) - sitar
          if (cov_nlpar == "s" & class == "sd" & !is.null(sncov_gr)) {
            if (x_i == paste0("lm", empty_sufx)) {
              if (!s_form_0_gr) {
                lm_gsubby <- paste0("lm", "_", 'sdx', "_", "cov", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", 'sdx', "_", "cov", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              }
            }
            # checks
            if (nlpar == "s" & !is.null(sncov_gr)) {
              if (length(evaluated_parameter) == 1) {
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              } else if (length(evaluated_parameter) == df) {
                repeach <- nrep_of_parms / df
                evaluated_parameter <-
                  rep(evaluated_parameter,
                      times = 1,
                      each = repeach)
              } else {
                #
              }
            } else {
              if (length(evaluated_parameter) < nrep_of_parms)
                evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
              if (length(evaluated_parameter) > nrep_of_parms)
                stop("prior elements for nlpar ",
                     nlpar, ", class ",  class,
                     " are greater than the parameter dimensions"
                )
            }
          }
          
          
          
          
          
          # scale sigma (class sd)
          if (nlpar == "" & class == "sd" & sigma_dpar == "sigma") {
            if (x_i == paste0("vsd", empty_sufx)) {
              eit <-  gsub("vsd", paste0("vsd", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else if (x_i == paste0("vmad", empty_sufx)) {
              eit <-  gsub("vmad", paste0("vmad", resp_), x_i)
              evaluated_parameter <- 1 * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- 1 * ept(eit)
              } else {
                stop(
                  "scale parameter options for distributional ",
                  sigma_dpar,
                  ", class ",
                  class,
                  " are:\n a numeric value (e.g., 2) or a charater like zzz",
                  "\n with zzz defined in the prior_data", 
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for distributional ",
                   sigma_dpar,
                   " are greater than the parameter dimensions")
          }
          
          
          
          # scale sigma cov (class sd)
          if (cov_sigma_dpar != "" & class == "sd" & sigma_dpar == "sigma" & 
              !is.null(sigmancov)) {
            if (x_i == paste0("sdacov", empty_sufx)) {
              eit <-  gsub("sdacov", paste0("acov_sd", resp_), x_i)
              evaluated_parameter <- ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- ept(eit)
              } else {
                stop(
                  "scale parameter options for distributional ",
                  sigma_dpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements fordistributional ",
                   sigma_dpar,
                   " are greater than the parameter dimensions")
          }
          
          
          
          
          
          # scale sigma (class sd)
          if (class == "sigma") {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("lme_rsd", empty_sufx)) {
              eit <-  gsub("lme_rsd", paste0("lme_rsd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("lm_rsd", empty_sufx)) {
              eit <-  gsub("lm_rsd", paste0("lm_rsd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          
          
          # scale dpar (class sd) sigma ~ 
          if (!is.null(dparncov) & class == "") {
            if (x_i == paste0("ysd", empty_sufx)) {
              eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("ymad", empty_sufx)) {
              eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("lme_rsd", empty_sufx)) {
              eit <-  gsub("lme_rsd", paste0("lme_rsd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else if (x_i == paste0("lm_rsd", empty_sufx)) {
              eit <-  gsub("lm_rsd", paste0("lm_rsd", resp_), x_i)
              evaluated_parameter <- scale_factor * ept(eit)
            } else {
              check_evalation_of_numeric_pdata_obj(
                prior_argument,
                p_str_in,
                x_i,
                x,
                pname_,
                dist,
                nlpar,
                class,
                allowed_parm_options,
                splitmvar_w2
              )
              if (is.numeric(eval(parse(text = x_i))) |
                  !is.null(eval(parse(text = x_i)))) {
                eit <- x_i
                evaluated_parameter <- scale_factor * ept(eit)
              } else {
                stop(
                  "scale parameter options for nlpar ",
                  nlpar,
                  ", class ",
                  class,
                  " are:\n lm, ymean, ymedian, a numeric value (e.g., 2)",
                  "or a charater such as zzz",
                  "\n with zzz defined in the prior_data",
                  "e.g., prior_data = list(zzz = 2)"
                )
              }
            }
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
        }
        
        
        
        
        # set degree of freedom df parameters -> for student_t
        
        if (grepl("^df$", pname_)) {
          check_evalation_of_numeric_pdata_obj(
            prior_argument,
            p_str_in,
            x_i,
            x,
            pname_,
            dist,
            nlpar,
            class,
            allowed_parm_options,
            splitmvar_w2
          )
          if (is.numeric(eval(parse(text = x_i))) |
              !is.null(eval(parse(text = x_i)))) {
            eit <- x_i
            evaluated_parameter <- ept(eit)
          } else {
            stop(
              "df parameter options for nlpar ",
              nlpar,
              ", class ",
              class,
              " are:\n zzzz, a numeric value (e.g., 2) or a charater such as zzz",
              "\n with zzz defined in the prior_data",
              "e.g., prior_data = list(zzz = 2)"
            )
          }
          # checks
          if (nlpar == "s" & !is.null(sncov)) {
            if (length(evaluated_parameter) == 1) {
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            } else if (length(evaluated_parameter) == df) {
              repeach <- nrep_of_parms / df
              evaluated_parameter <-
                rep(evaluated_parameter,
                    times = 1,
                    each = repeach)
            }
          } else {
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
        }
        
        
        # set nu_shape parameters -> for student_nu
        
        if (grepl("^nu_shape$", pname_)) {
          check_evalation_of_numeric_pdata_obj(
            prior_argument,
            p_str_in,
            x_i,
            x,
            pname_,
            dist,
            nlpar,
            class,
            allowed_parm_options,
            splitmvar_w2
          )
          if (is.numeric(eval(parse(text = x_i))) |
              !is.null(eval(parse(text = x_i)))) {
            eit <- x_i
            evaluated_parameter <- ept(eit)
          } else {
            stop(
              "df parameter options for nlpar ",
              nlpar,
              ", class ",
              class,
              " are:\n zzzz, a numeric value (e.g., 2) or a charater such as zzz",
              "\n with zzz defined in the prior_data",
              "e.g., prior_data = list(zzz = 2)"
            )
          }
          # checks
          if (nlpar == "s" & !is.null(sncov)) {
            if (length(evaluated_parameter) == 1) {
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            } else if (length(evaluated_parameter) == df) {
              repeach <- nrep_of_parms / df
              evaluated_parameter <-
                rep(evaluated_parameter,
                    times = 1,
                    each = repeach)
            }
          } else {
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
        }
        
        
        # set nu_scale parameters -> for student_nu
        
        if (grepl("^nu_scale$", pname_)) {
          check_evalation_of_numeric_pdata_obj(
            prior_argument,
            p_str_in,
            x_i,
            x,
            pname_,
            dist,
            nlpar,
            class,
            allowed_parm_options,
            splitmvar_w2
          )
          if (is.numeric(eval(parse(text = x_i))) |
              !is.null(eval(parse(text = x_i)))) {
            eit <- x_i
            evaluated_parameter <- ept(eit)
          } else {
            stop(
              "df parameter options for nlpar ",
              nlpar,
              ", class ",
              class,
              " are:\n zzzz, a numeric value (e.g., 2) or a charater such as zzz",
              "\n with zzz defined in the prior_data",
              "e.g., prior_data = list(zzz = 2)"
            )
          }
          # checks
          if (nlpar == "s" & !is.null(sncov)) {
            if (length(evaluated_parameter) == 1) {
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            } else if (length(evaluated_parameter) == df) {
              repeach <- nrep_of_parms / df
              evaluated_parameter <-
                rep(evaluated_parameter,
                    times = 1,
                    each = repeach)
            }
          } else {
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
        }
        
        
        
        
        # set rate parameters -> for exponential
        
        if (grepl("^rate$", pname_)) {
          if (x_i == paste0("ysd", empty_sufx)) {
            eit <-  gsub("ysd", paste0("ysd", resp_), x_i)
            evaluated_parameter <- 1 / (1 * ept(eit))
          } else if (x_i == paste0("ymad", empty_sufx)) {
            eit <-  gsub("ymad", paste0("ymad", resp_), x_i)
            evaluated_parameter <- 1 / (1 * ept(eit))
          } else if (x_i == paste0("lme_rsd", empty_sufx)) {
            eit <-  gsub("lme_rsd", paste0("lme_rsd", resp_), x_i)
            evaluated_parameter <- scale_factor * ept(eit)
          } else if (x_i == paste0("lm_rsd", empty_sufx)) {
            eit <-  gsub("lm_rsd", paste0("lm_rsd", resp_), x_i)
            evaluated_parameter <- scale_factor * ept(eit)
          } else if (x_i == paste0("lme_sd_a", empty_sufx)) {
            eit <-  gsub("lme_sd_a", paste0("lme_sd_a", resp_), x_i)
            evaluated_parameter <- scale_factor * ept(eit)
          } else {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- 1 / ept(eit)
            }
          }
          # checks
          if (nlpar == "s" & !is.null(sncov)) {
            if (length(evaluated_parameter) == 1) {
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            } else if (length(evaluated_parameter) == df) {
              repeach <- nrep_of_parms / df
              evaluated_parameter <-
                rep(evaluated_parameter,
                    times = 1,
                    each = repeach)
            }
          } else {
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
        }
        
        
        
        # set shape parameter -> for gamma inv_gamma (scale already covered)
        
        if (grepl("^shape$", pname_)) {
          check_evalation_of_numeric_pdata_obj(
            prior_argument,
            p_str_in,
            x_i,
            x,
            pname_,
            dist,
            nlpar,
            class,
            allowed_parm_options,
            splitmvar_w2
          )
          if (is.numeric(eval(parse(text = x_i))) |
              !is.null(eval(parse(text = x_i)))) {
            eit <- x_i
            evaluated_parameter <- ept(eit)
          } else {
            stop(
              "df parameter options for nlpar ",
              nlpar,
              ", class ",
              class,
              " are:\n zzzz, a numeric value (e.g., 2) or a charater such as zzz",
              "\n with zzz defined in the prior_data", 
              "e.g., prior_data = list(zzz = 2)"
            )
          }
          # checks
          if (nlpar == "s" & !is.null(sncov)) {
            if (length(evaluated_parameter) == 1) {
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            } else if (length(evaluated_parameter) == df) {
              repeach <- nrep_of_parms / df
              evaluated_parameter <-
                rep(evaluated_parameter,
                    times = 1,
                    each = repeach)
            }
          } else {
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
        }
        
        
        
        # set lower upper parameters -> for uniform
        
        if (grepl("^lower$", pname_)) {
          if (x_i == paste0("lm", empty_sufx)) {
            if (nlpar == "a" & class == "b" & grepl("a", fixedsi)) {
              if (a_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)[1]
            }
            if (nlpar == "s") {
              if (s_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            }
          } else {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "lower parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the", 
                "prior_data e.g., prior_data = list(zzz = 2)"
              )
            }
          }
          # checks
          if (nlpar == "s" & !is.null(sncov)) {
            if (length(evaluated_parameter) == 1) {
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            } else if (length(evaluated_parameter) == df) {
              repeach <- nrep_of_parms / df
              evaluated_parameter <-
                rep(evaluated_parameter,
                    times = 1,
                    each = repeach)
            }
          } else {
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          evaluated_parameter <- evaluated_parameter - addrange
          evaluated_parameter_lower <- evaluated_parameter
        }
        
        
        
        if (grepl("^upper$", pname_)) {
          if (x_i == paste0("lm", empty_sufx)) {
            if (nlpar == "a" & class == "b" & grepl("a", fixedsi)) {
              if (a_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)[1]
            }
            if (nlpar == "s") {
              if (s_form_0) {
                lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
              } else {
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              }
              eit <-  gsub("lm", lm_gsubby, x_i)
              evaluated_parameter <- ept(eit)
            }
          } else {
            check_evalation_of_numeric_pdata_obj(
              prior_argument,
              p_str_in,
              x_i,
              x,
              pname_,
              dist,
              nlpar,
              class,
              allowed_parm_options,
              splitmvar_w2
            )
            if (is.numeric(eval(parse(text = x_i))) |
                !is.null(eval(parse(text = x_i)))) {
              eit <- x_i
              evaluated_parameter <- ept(eit)
            } else {
              stop(
                "upper parameter options for nlpar ",
                nlpar,
                ", class ",
                class,
                " are:\n lm, a numeric value (e.g., 2) or a charater such as zzz",
                "\n with zzz defined in the prior_data",
                "e.g., prior_data = list(zzz = 2)"
              )
            }
          }
          # checks
          if (nlpar == "s" & !is.null(sncov)) {
            if (length(evaluated_parameter) == 1) {
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            } else if (length(evaluated_parameter) == df) {
              repeach <- nrep_of_parms / df
              evaluated_parameter <-
                rep(evaluated_parameter,
                    times = 1,
                    each = repeach)
            }
          } else {
            if (length(evaluated_parameter) < nrep_of_parms)
              evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
            if (length(evaluated_parameter) > nrep_of_parms)
              stop("prior elements for nlpar ",
                   nlpar, ", class ",  class,
                   " are greater than the parameter dimensions")
          }
          evaluated_parameter <- evaluated_parameter + addrange
          evaluated_parameter_upper <- evaluated_parameter
        }
        
        
        
        # set eta parameter -> for lkj - also for mvr rescor
        
        if (grepl("^eta$", pname_)) {
          check_evalation_of_numeric_pdata_obj(
            prior_argument,
            p_str_in,
            x_i,
            x,
            pname_,
            dist,
            nlpar,
            class,
            allowed_parm_options,
            splitmvar_w2
          )
          if (is.numeric(eval(parse(text = x_i))) |
              !is.null(eval(parse(text = x_i)))) {
            eit <- x_i
            evaluated_parameter <- ept(eit)
          } else {
            stop(
              "df parameter options for nlpar ",
              nlpar,
              ", class ",
              class,
              " are:\n zzzz, a numeric value (e.g., 2) or a charater such as zzz",
              "\n with zzz defined in the prior_data",
              "e.g., prior_data = list(zzz = 2)"
            )
          }
          if (length(evaluated_parameter) < nrep_of_parms)
            evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
          if (length(evaluated_parameter) > nrep_of_parms)
            stop("prior elements for nlpar ",
                 nlpar, ", class ",  class,
                 " are greater than the parameter dimensions")
        }
        
        
        # set autocr parameter -> ar ma arma
        
        if (setautocorr & class != "b") {
          check_evalation_of_numeric_pdata_obj(
            prior_argument,
            p_str_in,
            x_i,
            x,
            pname_,
            dist,
            nlpar,
            class,
            allowed_parm_options,
            splitmvar_w2
          )
          if (is.numeric(eval(parse(text = x_i))) |
              !is.null(eval(parse(text = x_i)))) {
            eit <- x_i
            evaluated_parameter <- ept(eit)
          } else {
            stop(
              "df parameter options for nlpar ",
              nlpar,
              ", class ",
              class,
              " are:\n zzzz, a numeric value (e.g., 2) or a charater such as zzz",
              "\n with zzz defined in the prior_data",
              "e.g., prior_data = list(zzz = 2)"
            )
          }
          if (length(evaluated_parameter) < nrep_of_parms)
            evaluated_parameter <- rep(evaluated_parameter, nrep_of_parms)
          if (length(evaluated_parameter) > nrep_of_parms)
            stop("prior elements for nlpar ",
                 nlpar, ", class ",  class,
                 " are greater than the parameter dimensions")
        }
        
        
        
        # make unique names
        
        if (nlpar != "")
          prefix <- nlpar
        if (nlpar != "" &
            cov_nlpar != "")
          prefix <- paste0(cov_nlpar, "_", "cov")
        if (class == 'cor')
          prefix <- 'lkj'
        if (class == 'rescor')
          prefix <- 'lkj'
        if (class == 'sigma')
          prefix <- 'sigma'
        if (setautocorr)
          prefix <- class
        if (dpar != "")
          prefix <- dpar
        
        
        if (sigma_dpar == 'sigma' ) {
          prefix <- 'sigma'
        }
        
        if (cov_sigma_dpar == 'sigma_cov' ) {
          prefix <- 'sigma_cov'
        }
        
        
        
        
        
        if (setautocorr) {
          add_cla_to_name <- NULL
        } else if (class == "" |
                   class == "sigma" | dpar != "" & !is.null(dparncov)) {
          add_cla_to_name <- NULL
        } else {
          add_cla_to_name <- paste0(sep_indicator, class)
        }
        
        
        
        # This is required to set unique stanvar names for higher level sd 
        
        if (class == 'sd' | class == 'cor') {
          group_arg_groupvar_paste <-  group # group_arg_groupvar
          group_arg_groupvar_paste <- gsub(":", "_", group_arg_groupvar_paste)
          suppressWarnings({
            set_str_names[i] <- paste0(set_str_names[i], "_", group_arg_groupvar_paste)
          })
        }
        
        
        # This is required to set unique stanvar names forsigma_prior_cor and gr_prior_cor
        
        if (class == 'cor') {
          group_arg_groupvar_paste <-  group
          if (sigma_dpar == 'sigma') lkj_arg_cor_paste <-  'sigma'
          if (sigma_dpar != 'sigma') lkj_arg_cor_paste <-  'gr'
          set_str_names[i] <- paste0(set_str_names[i], "_", 
                                     lkj_arg_cor_paste)
        }
        
        
        
        
        
        
        name_parameter <-
          paste0(prefix,
                 add_cla_to_name,
                 sep_indicator,
                 set_str_names[i], 
                 resp_)
        
        # name_parameter <- paste0(name_parameter, add_gr_id)
        
        assign(name_parameter, evaluated_parameter)
        
        if (change_default_data_pll_args) {
          stanvars_data[[name_parameter]] <-
            brms::stanvar(
              eval(parse(text = name_parameter)),
              name = name_parameter,
              block = "data",
              pll_args = ept(set_data_pll_args)
            )
        } else {
          stanvars_data[[name_parameter]] <-
            brms::stanvar(eval(parse(text = name_parameter)),
                          name = name_parameter, block = "data")
        }
        
        collect_name_parameter <-
          c(collect_name_parameter, name_parameter)
      }
      
      
      
      
      # assigning bounds
      
      name_lb <-
        paste0(prefix, 
               add_cla_to_name, 
               sep_indicator, 
               "lb",
               resp_)
      
      name_ub <-
        paste0(prefix, 
               add_cla_to_name, 
               sep_indicator, 
               "ub", 
               resp_)
      
      
      
      if (grepl("^lb$", pname_)  & !grepl("b_", x_i, fixed = T) ) {
        # if (grepl("^lb$", pname_)) {
        if (!(is.na(eval(parse(text = x_i))) |
              eval(parse(text = x_i)) == "NA")) {
          lowerbound <- eval(parse(text = x_i))
          if (length(lowerbound < length(evaluated_parameter))) {
            lowerbound <- rep(lowerbound, length(evaluated_parameter))
          }
          assign(name_lb, lowerbound)
          if (change_default_data_pll_args) {
            stanvars_data[[name_lb]] <- brms::stanvar(
              eval(parse(text = name_lb)),
              name = name_lb,
              block = "data",
              pll_args = ept(set_data_pll_args)
            )
          } else {
            stanvars_data[[name_lb]] <- 
              brms::stanvar(eval(parse(text = name_lb)),
                            name = name_lb, block = "data")
          }
        } else {
          lowerbound <- NA
          if (length(lowerbound < length(evaluated_parameter))) {
            lowerbound <- rep(lowerbound, length(evaluated_parameter))
          }
          assign(name_lb, rep(lowerbound, length(evaluated_parameter)))
        }
      } # if (grepl("^lb$", pname_)) {
      
      # added on 2/8/2023 to allow parameter as lb ub bound 
      
      if (grepl("^lb$", pname_)  & grepl("b_", x_i, fixed = T) ) {
        set_x_i_p_bound <- x_i
        if(resp_ != "") {
          set_x_i_p_bound <- gsub("b_", paste0("b", resp_, "_"), 
                                  set_x_i_p_bound, fixed = T)
        } 
        lowerbound <- set_x_i_p_bound
        assign(name_lb, lowerbound)
      }
      
      
      if (grepl("^ub$", pname_) & !grepl("b_", x_i, fixed = T) ) {
        # if (grepl("^ub$", pname_)) {
        if (!(is.na(eval(parse(text = x_i))) |
              eval(parse(text = x_i)) == "NA")) {
          upperbound <- eval(parse(text = x_i))
          if (length(upperbound < length(evaluated_parameter))) {
            upperbound <- rep(upperbound, length(evaluated_parameter))
          }
          assign(name_ub, upperbound)
          if (change_default_data_pll_args) {
            stanvars_data[[name_ub]] <- brms::stanvar(
              eval(parse(text = name_ub)),
              name = name_ub,
              block = "data",
              pll_args = ept(set_data_pll_args)
            )
          } else {
            stanvars_data[[name_ub]] <- 
              brms::stanvar(eval(parse(text = name_ub)), 
                            name = name_ub, block = "data")
          }
        } else {
          upperbound <- NA
          if (length(upperbound < length(evaluated_parameter))) {
            upperbound <- rep(upperbound, length(evaluated_parameter))
          }
          assign(name_ub, rep(upperbound, length(evaluated_parameter)))
        }
      } # if (grepl("^ub$", pname_)) {
      
      # added on 2/8/2023 to allow parameter as lb ub bound 
      
      if (grepl("^ub$", pname_)  & grepl("b_", x_i, fixed = T) ) {
        set_x_i_p_bound <- x_i
        if(resp_ != "") {
          set_x_i_p_bound <- gsub("b_", paste0("b", resp_, "_"), 
                                  set_x_i_p_bound, fixed = T)
        } 
        upperbound <- set_x_i_p_bound
        assign(name_ub, upperbound)
      }
      
      
      if (grepl("^lb$", pname_)) {
        if (dist == "lognormal" |
            dist == "gamma" | dist == "inv_gamma" | dist == "exponential") {
          if (all(is.na(lowerbound) |
                  lowerbound == "NA"))
            lowerbound <- name_lb
          assign(name_lb, rep(0, length(evaluated_parameter)))
          if (change_default_data_pll_args) {
            stanvars_data[[name_lb]] <- brms::stanvar(
              eval(parse(text = name_lb)),
              name = name_lb,
              block = "data",
              pll_args = ept(set_data_pll_args)
            )
          } else {
            stanvars_data[[name_lb]] <- 
              brms::stanvar(eval(parse(text = name_lb)),
                            name = name_lb, block = "data")
          }
        }
      }
    } # end of loop for (i in 1:length(x)) {
    
    
    
    
    if (dist == "uniform" ) {
      if (all(is.na(lowerbound) |
              lowerbound == "NA"))
        lowerbound <- name_lb
      assign(name_lb, evaluated_parameter_lower)
      
      if (change_default_data_pll_args) {
        stanvars_data[[name_lb]] <- brms::stanvar(
          eval(parse(text = name_lb)),
          name = name_lb,
          block = "data",
          pll_args = ept(set_data_pll_args)
        )
      } else {
        stanvars_data[[name_lb]] <- 
          brms::stanvar(eval(parse(text = name_lb)),
                        name = name_lb, block = "data")
      }
      
      
      if (all(is.na(upperbound) |
              upperbound == "NA"))
        upperbound <- name_ub
      assign(name_ub, evaluated_parameter_upper)
      
      if (change_default_data_pll_args) {
        stanvars_data[[name_ub]] <- brms::stanvar(
          eval(parse(text = name_ub)),
          name = name_ub,
          block = "data",
          pll_args = ept(set_data_pll_args)
        )
      } else {
        stanvars_data[[name_ub]] <- 
          brms::stanvar(eval(parse(text = name_ub)),
                        name = name_ub, block = "data")
      }
      
      
      if ((identical(evaluated_parameter_lower, evaluated_parameter_upper))) {
        stop(
          "the lower and upper parameters for uniform distribution are identical",
          "\n This could be because of the same values used for lower and upper",
          "\n", 
          " parameters with addrange set as '0",
          "\n Either change the lower and upper parameter values or set",
          "\n", 
          " addrange other than zero"
        )
      }
      
      for (i in 1:length(evaluated_parameter_lower)) {
        if (evaluated_parameter_lower[i] >= evaluated_parameter_upper[i]) {
          stop(
            "lower parameter value at position '",
            i,
            "' ",
            evaluated_parameter_lower[i],
            "should be less than the specified upper parameter value ",
            evaluated_parameter_lower[i]
          )
        }
      }
    }
    
    
    # name_parameter
    if (dist != "student_nu") {
      if (nrep_of_parms != 1) {
        prior_str_arg_out_c <- c()
        for (i in 1:nrep_of_parms) {
          if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
            tt <- paste0(collect_name_parameter, collapse = ", ")
          } else {
            tt <- paste0(collect_name_parameter, "[", i, "]", collapse = ", ")
          }
          prior_str_arg_out_c <- c(prior_str_arg_out_c, tt)
        }
        prior_str_arg_out <- prior_str_arg_out_c
        prior_str_arg_out <- paste0(dist, "(", prior_str_arg_out, ")")
      } else {
        prior_str_arg_out <-
          paste0(dist,
                 "(",
                 paste(collect_name_parameter, collapse = ", "),
                 ")")
      }
    }
    
    
    if (dist == "student_nu") {
      collect_name_parameter_copy <- collect_name_parameter
      collect_name_parameter <-
        c(gsub("_shape", "", collect_name_parameter[1]),
          collect_name_parameter[3:4])
      if (nrep_of_parms != 1) {
        prior_str_arg_out_c <- c()
        for (i in 1:nrep_of_parms) {
          if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
            tt <- paste0(collect_name_parameter, collapse = ", ")
          } else {
            tt <- paste0(collect_name_parameter, "[", i, "]", collapse = ", ")
          }
          prior_str_arg_out_c <- c(prior_str_arg_out_c, tt)
        }
        prior_str_arg_out <- prior_str_arg_out_c
        prior_str_arg_out <-
          paste0("student_t", "(", prior_str_arg_out, ")")
      } else {
        prior_str_arg_out <-
          paste0("student_t",
                 "(",
                 paste(collect_name_parameter, collapse = ", "),
                 ")")
      }
    }
    
    
    if (dist == "student_nu") {
      student_nu_left <-
        gsub("_shape", "", collect_name_parameter_copy[1])
      student_nu_right1 <- collect_name_parameter_copy[1]
      student_nu_right2 <- collect_name_parameter_copy[2]
      stanvars_data[[paste0(student_nu_left, "_", "", "_", "pblock")]] <-
        brms::stanvar(
          scode = paste0(
            "vector<lower=1>[",
            nrep_of_parms ,
            "] ",
            student_nu_left,
            ";"
          ),
          block = "parameter",
          position = "end"
        ) # , pll_args = paste0("vector"," ", hptau_nu)
      
      t_or_lp <- 'target'
      dist_student_nu <- "gamma"
      if (nrep_of_parms != 1) {
        for (i in 1:nrep_of_parms) {
          if (normalize) {
            dist_student_1 <-
              "lpdf"
            dist_student_2 <- "lccdf"
            svarblock <- 'model' # 'tparameters'
            define_studentdisttype_1 <-
              paste0(dist_student_nu, "_", dist_student_1)
            define_studentdisttype_2 <-
              paste0(dist_student_nu, "_", dist_student_2)
            define_studentscode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_studentdisttype_1,
                "(",
                paste0(student_nu_left, "[", i, "]"),
                " | ",
                paste0(student_nu_right1, "[", i, "]") ,
                ", ",
                paste0(student_nu_right2, "[", i, "]") ,
                ")"
              )
            define_studentscode_2 <-
              paste0(
                nrep_of_parms,
                " * ",
                define_studentdisttype_2,
                "(",
                1,
                " | ",
                paste0(student_nu_right1, "[", i, "]") ,
                ", ",
                paste0(student_nu_right2, "[", i, "]") ,
                ")"
              )
            define_studentscode <-
              paste0(define_studentscode_1,
                     "\n    - ",
                     define_studentscode_2,
                     ";")
          } else {
            dist_student_1 <- "lupdf"
            svarblock <- 'model'
            define_studentdisttype_1 <-
              paste0(dist_student_nu, "_", dist_student_1)
            define_studentscode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_studentdisttype_1,
                "(",
                paste0(student_nu_left, "[", i, "]"),
                " | ",
                paste0(student_nu_right1, "[", i, "]")  ,
                ", ",
                paste0(student_nu_right2, "[", i, "]") ,
                ")"
              )
            define_studentscode <- paste0(define_studentscode_1, ";")
          }
          stanvars_data[[paste0(student_nu_left, "_", i, "_", "mblock")]] <-
            brms::stanvar(scode = define_studentscode,
                          block = svarblock,
                          position = "end")
          
        }
      }
      
      
      if (nrep_of_parms == 1) {
        if (normalize) {
          dist_student_1 <-
            "lpdf"
          dist_student_2 <- "lccdf"
          svarblock <- 'model' # 'tparameters'
          define_studentdisttype_1 <-
            paste0(dist_student_nu, "_", dist_student_1)
          define_studentdisttype_2 <-
            paste0(dist_student_nu, "_", dist_student_2)
          define_studentscode_1 <-
            paste0(
              t_or_lp,
              " += ",
              define_studentdisttype_1,
              "(",
              student_nu_left,
              " | ",
              student_nu_right1,
              ", ",
              student_nu_right2,
              ")"
            )
          define_studentscode_2 <-
            paste0(
              nrep_of_parms,
              " * ",
              define_studentdisttype_2,
              "(",
              1,
              " | ",
              student_nu_right1,
              ", ",
              student_nu_right2,
              ")"
            )
          define_studentscode <-
            paste0(define_studentscode_1,
                   "\n    - ",
                   define_studentscode_2,
                   ";")
        } else {
          dist_student_1 <- "lupdf"
          svarblock <- 'model'
          define_studentdisttype_1 <-
            paste0(dist_student_nu, "_", dist_student_1)
          define_studentscode_1 <-
            paste0(
              t_or_lp,
              " += ",
              define_studentdisttype_1,
              "(",
              student_nu_left,
              " | ",
              student_nu_right1,
              ", ",
              student_nu_right2,
              ")"
            )
          define_studentscode <- paste0(define_studentscode_1, ";")
        }
        stanvars_data[[paste0(student_nu_left, "_", "", "_", "mblock")]] <-
          brms::stanvar(scode = define_studentscode,
                        block = svarblock,
                        position = "end")
        
      }
    }
    
    
    
    
    
    
    if (!is.null(stanvars_data[[name_lb]])) {
      if (nrep_of_parms == 1) {
        lowerbound <- name_lb
      } else {
        if (any(is.na(lowerbound))) {
          lowerbound <- paste0(name_lb, "[", 1:nrep_of_parms, "]")
        } else {
          lowerbound <- name_lb
        }
      }
    }
    
    if (!is.null(stanvars_data[[name_ub]])) {
      if (nrep_of_parms == 1) {
        upperbound <- name_ub
      } else {
        if (any(is.na(upperbound))) {
          upperbound <- paste0(name_ub, "[", 1:nrep_of_parms, "]")
        } else {
          upperbound <- name_ub
        }
      }
    }
    
    
    
    
    if (ept(sethp)) {
      original_scale <-
        paste0(prefix, add_cla_to_name, sep_indicator, "scale", resp_)
      tauid <- 'tau'
      hptau <-
        paste0(prefix, add_cla_to_name, sep_indicator, tauid, resp_)
      prior_str_arg_out <-
        gsub(original_scale, hptau, prior_str_arg_out, fixed = T)
      
      if (dist == "student_t") {
        original_df <-
          paste0(prefix, add_cla_to_name, sep_indicator, "df", resp_)
        original_df_val <-
          ept(ept(
            paste0(prefix, add_cla_to_name, sep_indicator, "df", resp_)
          ))
        hptau_df <-
          paste0(
            prefix,
            add_cla_to_name,
            sep_indicator,
            paste0(tauid, sep_indicator, 'df'),
            resp_
          )
        if (length(original_df_val) < nrep_of_parms)
          original_df_val <- rep(original_df_val, nrep_of_parms)
        assign(hptau_df, original_df_val)
      } else {
        hptau_df <-
          paste0(
            prefix,
            add_cla_to_name,
            sep_indicator,
            paste0(tauid, sep_indicator, 'df'),
            resp_
          )
        assign(hptau_df, rep(3, nrep_of_parms))
      }
      
      if (sethp_dist == "student_nu") {
        hptau_nu <-
          paste0(
            prefix,
            add_cla_to_name,
            sep_indicator,
            paste0(tauid, sep_indicator, 'nu'),
            resp_
          )
        hptau_nu_shape <-
          paste0(
            prefix,
            add_cla_to_name,
            sep_indicator,
            paste0(tauid, sep_indicator, 'nu_shape'),
            resp_
          )
        hptau_nu_scale <-
          paste0(
            prefix,
            add_cla_to_name,
            sep_indicator,
            paste0(tauid, sep_indicator, 'nu_scale'),
            resp_
          )
        assign(hptau_nu_shape, rep(2, nrep_of_parms))
        assign(hptau_nu_scale, rep(0.1, nrep_of_parms))
      }
      
      hptau_scale    <-
        paste0(
          prefix,
          add_cla_to_name,
          sep_indicator,
          paste0(tauid, sep_indicator, 'scale'),
          resp_
        )
      hptau_location <-
        paste0(
          prefix,
          add_cla_to_name,
          sep_indicator,
          paste0(tauid, sep_indicator, 'location'),
          resp_
        )
      
      if (sethp_dist == "exponential") {
        hptau_rate    <-
          paste0(
            prefix,
            add_cla_to_name,
            sep_indicator,
            paste0(tauid, sep_indicator, 'rate'),
            resp_
          )
      }
      
      
      t_or_lp <- 'target'
      if (sethp_dist == "normal" |
          sethp_dist == "cauchy" |
          sethp_dist == "student_t" |
          sethp_dist == "student_nu" |
          sethp_dist == "exponential") {
        if (sethp_dist == "normal") {
          if (normalize) {
            dist_1 <-
              "lpdf"
            dist_2 <- "lccdf"
            svarblock <- 'model' # 'tparameters'
            define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
            define_disttype_2 <- paste0(sethp_dist, "_", dist_2)
            define_scode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_disttype_1,
                "(",
                hptau,
                " | ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode_2 <-
              paste0(
                nrep_of_parms,
                " * ",
                define_disttype_2,
                "(",
                0,
                " | ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode <-
              paste0(define_scode_1, "\n    - ", define_scode_2, ";")
          } else {
            dist_1 <- "lupdf"
            svarblock <- 'model'
            define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
            define_scode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_disttype_1,
                "(",
                hptau,
                " | ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode <- paste0(define_scode_1, ";")
          }
          stanvars_data[[paste0(tauid, "_", "mblock")]] <-
            brms::stanvar(scode = define_scode,
                          block = svarblock,
                          position = "end")
        }
        if (sethp_dist == "cauchy") {
          if (normalize) {
            dist_1 <-
              "lpdf"
            dist_2 <- "lccdf"
            svarblock <- 'model' # 'tparameters'
            define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
            define_disttype_2 <- paste0(sethp_dist, "_", dist_2)
            define_scode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_disttype_1,
                "(",
                hptau,
                " | ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode_2 <-
              paste0(
                nrep_of_parms,
                " * ",
                define_disttype_2,
                "(",
                0,
                " | ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode <-
              paste0(define_scode_1, "\n    - ", define_scode_2, ";")
          } else {
            dist_1 <- "lupdf"
            svarblock <- 'model'
            define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
            define_scode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_disttype_1,
                "(",
                hptau,
                " | ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode <- paste0(define_scode_1, ";")
          }
          stanvars_data[[paste0(tauid, "_", "mblock")]] <-
            brms::stanvar(scode = define_scode,
                          block = svarblock,
                          position = "end")
        }
        if (sethp_dist == "student_t") {
          if (normalize) {
            dist_1 <-
              "lpdf"
            dist_2 <- "lccdf"
            svarblock <- 'model' # 'tparameters'
            define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
            define_disttype_2 <- paste0(sethp_dist, "_", dist_2)
            define_scode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_disttype_1,
                "(",
                hptau,
                " | ",
                hptau_df,
                ", ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode_2 <-
              paste0(
                nrep_of_parms,
                " * ",
                define_disttype_2,
                "(",
                0,
                " | ",
                hptau_df,
                ", ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode <-
              paste0(define_scode_1, "\n    - ", define_scode_2, ";")
          } else {
            dist_1 <- "lupdf"
            svarblock <- 'model'
            define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
            define_scode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_disttype_1,
                "(",
                hptau,
                " | ",
                hptau_df,
                ", ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode <- paste0(define_scode_1, ";")
          }
          stanvars_data[[paste0(tauid, "_", "mblock")]] <-
            brms::stanvar(scode = define_scode,
                          block = svarblock,
                          position = "end")
        }
        if (sethp_dist == "student_nu") {
          sethp_dist_student_t <- "student_t"
          if (normalize) {
            dist_1 <-
              "lpdf"
            dist_2 <- "lccdf"
            svarblock <- 'model' # 'tparameters'
            define_disttype_1 <-
              paste0(sethp_dist_student_t, "_", dist_1)
            define_disttype_2 <-
              paste0(sethp_dist_student_t, "_", dist_2)
            define_scode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_disttype_1,
                "(",
                hptau,
                " | ",
                hptau_nu,
                ", ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode_2 <-
              paste0(
                nrep_of_parms,
                " * ",
                define_disttype_2,
                "(",
                1,
                " | ",
                hptau_nu,
                ", ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode <-
              paste0(define_scode_1, "\n    - ", define_scode_2, ";")
          } else {
            dist_1 <- "lupdf"
            svarblock <- 'model'
            define_disttype_1 <-
              paste0(sethp_dist_student_t, "_", dist_1)
            define_scode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_disttype_1,
                "(",
                hptau,
                " | ",
                hptau_nu,
                ", ",
                hptau_location,
                ", ",
                hptau_scale,
                ")"
              )
            define_scode <- paste0(define_scode_1, ";")
          }
          stanvars_data[[paste0(tauid, "_", "mblock")]] <-
            brms::stanvar(scode = define_scode,
                          block = svarblock,
                          position = "end")
          
          # nu gamma
          sethp_dist_student_nu <- "gamma"
          if (normalize) {
            dist_student_1 <-
              "lpdf"
            dist_student_2 <- "lccdf"
            svarblock <- 'model' # 'tparameters'
            define_studentdisttype_1 <-
              paste0(sethp_dist_student_nu, "_", dist_student_1)
            define_studentdisttype_2 <-
              paste0(sethp_dist_student_nu, "_", dist_student_2)
            define_studentscode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_studentdisttype_1,
                "(",
                hptau_nu,
                " | ",
                hptau_nu_shape,
                ", ",
                hptau_nu_scale,
                ")"
              )
            define_studentscode_2 <-
              paste0(
                nrep_of_parms,
                " * ",
                define_studentdisttype_2,
                "(",
                1,
                " | ",
                hptau_nu_shape,
                ", ",
                hptau_nu_scale,
                ")"
              )
            define_studentscode <-
              paste0(define_studentscode_1,
                     "\n    - ",
                     define_studentscode_2,
                     ";")
          } else {
            dist_student_1 <- "lupdf"
            svarblock <- 'model'
            define_studentdisttype_1 <-
              paste0(sethp_dist_student_nu, "_", dist_student_1)
            define_studentscode_1 <-
              paste0(
                t_or_lp,
                " += ",
                define_studentdisttype_1,
                "(",
                hptau_nu,
                " | ",
                hptau_nu_shape,
                ", ",
                hptau_nu_scale,
                ")"
              )
            define_studentscode <- paste0(define_studentscode_1, ";")
          }
          stanvars_data[[paste0(tauid, "_", "nu", "_", "mblock")]] <-
            brms::stanvar(scode = define_studentscode,
                          block = svarblock,
                          position = "end")
        }
        
        if (sethp_dist == "exponential") {
          if (normalize) {
            dist_1 <-
              "lpdf"
            dist_2 <- "lccdf"
            svarblock <- 'model' # 'tparameters'
            define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
            define_disttype_2 <- paste0(sethp_dist, "_", dist_2)
            define_scode_1 <-
              paste0(t_or_lp,
                     " += ",
                     define_disttype_1,
                     "(",
                     hptau,
                     " | ",
                     hptau_rate,
                     ")")
            define_scode_2 <-
              paste0(nrep_of_parms,
                     " * ",
                     define_disttype_2,
                     "(",
                     0,
                     " | ",
                     hptau_rate,
                     ")")
            define_scode <-
              paste0(define_scode_1, "\n    - ", define_scode_2, ";")
          } else {
            dist_1 <- "lupdf"
            svarblock <- 'model'
            define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
            define_scode_1 <-
              paste0(t_or_lp,
                     " += ",
                     define_disttype_1,
                     "(",
                     hptau,
                     " | ",
                     hptau_rate,
                     ")")
            define_scode <- paste0(define_scode_1, ";")
          }
          stanvars_data[[paste0(tauid, "_", "mblock")]] <-
            brms::stanvar(scode = define_scode,
                          block = svarblock,
                          position = "end")
        }
      } else {
        sethp_dist <- 'normal'
        if (normalize) {
          dist_1 <-
            "lpdf"
          dist_2 <- "lccdf"
          svarblock <- 'model' # 'tparameters'
          define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
          define_disttype_2 <- paste0(sethp_dist, "_", dist_2)
          define_scode_1 <-
            paste0(
              t_or_lp,
              " += ",
              define_disttype_1,
              "(",
              hptau,
              " | ",
              hptau_location,
              ", ",
              hptau_scale,
              ")"
            )
          define_scode_2 <-
            paste0(
              nrep_of_parms,
              " * ",
              define_disttype_2,
              "(",
              0,
              " | ",
              hptau_location,
              ", ",
              hptau_scale,
              ")"
            )
          define_scode <-
            paste0(define_scode_1, "\n    - ", define_scode_2, ";")
        } else {
          dist_1 <- "lupdf"
          svarblock <- 'model'
          define_disttype_1 <- paste0(sethp_dist, "_", dist_1)
          define_scode_1 <-
            paste0(
              t_or_lp,
              " += ",
              define_disttype_1,
              "(",
              hptau,
              " | ",
              hptau_location,
              ", ",
              hptau_scale,
              ")"
            )
          define_scode <- paste0(define_scode_1, ";")
        }
        stanvars_data[[paste0(tauid, "_", "mblock")]] <-
          brms::stanvar(scode = define_scode,
                        block = svarblock,
                        position = "end")
      }
      
      stanvars_data[[paste0(tauid, "_", "pblock")]] <-
        brms::stanvar(
          scode = paste0("vector<lower=0>[", nrep_of_parms , "] ", hptau, ";"),
          block = "parameter",
          position = "end"
        ) 
      
      
      if (sethp_dist == "student_nu") {
        stanvars_data[[paste0(tauid, "_", "nu", "_", "pblock")]] <-
          brms::stanvar(
            scode = paste0("vector<lower=1>[", 
                           nrep_of_parms , "] ", hptau_nu, ";"),
            block = "parameter",
            position = "end"
          ) 
      }
      
      # add data stanvars
      if (sethp_dist == "normal" |
          sethp_dist == "cauchy" |
          sethp_dist == "student_nu" |
          sethp_dist == "student_t") {
        if (sethp_dist == "student_t") {
          if (change_default_data_pll_args) {
            stanvars_data[[hptau_df]] <- brms::stanvar(
              ept(hptau_df),
              name = hptau_df,
              block = "data",
              pll_args = ept(set_data_pll_args)
            )
          } else {
            stanvars_data[[hptau_df]] <- brms::stanvar(ept(hptau_df),
                                                       name = hptau_df, block = "data")
          }
        }
        
        
        if (sethp_dist == "student_nu") {
          if (change_default_data_pll_args) {
            stanvars_data[[hptau_nu_shape]] <- brms::stanvar(
              ept(hptau_nu_shape),
              name = hptau_nu_shape,
              block = "data",
              pll_args = ept(set_data_pll_args)
            )
            
            stanvars_data[[hptau_nu_scale]] <-
              brms::stanvar(
                ept(hptau_nu_scale),
                name = hptau_nu_scale,
                block = "data",
                pll_args = ept(set_data_pll_args)
              )
          } else {
            stanvars_data[[hptau_nu_shape]] <- brms::stanvar(ept(hptau_nu_shape),
                                                             name = hptau_nu_shape,
                                                             block = "data")
            
            stanvars_data[[hptau_nu_scale]] <-
              brms::stanvar(ept(hptau_nu_scale),
                            name = hptau_nu_scale,
                            block = "data")
          }
        }
        
        
        if (change_default_data_pll_args) {
          stanvars_data[[hptau_location]] <- brms::stanvar(
            rep(0, nrep_of_parms),
            name = hptau_location,
            block = "data",
            pll_args = ept(set_data_pll_args)
          )
        } else {
          stanvars_data[[hptau_location]] <- brms::stanvar(rep(0, nrep_of_parms),
                                                           name = hptau_location,
                                                           block = "data")
        }
        
        if (change_default_data_pll_args) {
          stanvars_data[[hptau_scale]]    <-
            brms::stanvar(
              eval(parse(text = original_scale)),
              name = hptau_scale,
              block = "data",
              pll_args = ept(set_data_pll_args)
            )
        } else {
          stanvars_data[[hptau_scale]]    <-
            brms::stanvar(eval(parse(text = original_scale)),
                          name = hptau_scale, block = "data")
        }
      }
      
      
      if (sethp_dist == "exponential") {
        if (change_default_data_pll_args) {
          stanvars_data[[hptau_rate]] <- brms::stanvar(
            1 / ept(original_scale),
            name = hptau_rate,
            block = "data",
            pll_args = ept(set_data_pll_args)
          )
        } else {
          stanvars_data[[hptau_rate]] <- brms::stanvar(1 / ept(original_scale),
                                                       name = hptau_rate,
                                                       block = "data")
        }
      }
      stanvars_data[[original_scale]] <- NULL
    }
    
    
  } # end if(dist != 'flat') 
  
  
  if(dist == 'flat') {
    prior_str_arg_out <- ""
    lowerbound <- NA
    upperbound <- NA
    stanvars_data <- NULL
    allowed_init_options_beta <- NULL
    allowed_init_options_sd <- NULL
    allowed_init_options_rate <- NULL
    allowed_init_options_shape <- NULL
    allowed_init_options_scale <- NULL
  }
  
  
  # initials
  if (initsi != "random") {
    # parm <- nlpar
    if(sigma_dpar == 'sigma')  {
      parm <- sigma_dpar
    } else {
      parm <- nlpar
    }
    stanvars_datazz <- stanvars_data
    pstrarg <- prior_str_arg_out
    
    init_internal_args_names <- c(
      'parm',
      'class',
      'dpar',
      'sigma_dpar',
      'resp_',
      'dist',
      'lowerbound',
      'upperbound',
      'allowed_init_options_beta',
      'allowed_init_options_sd',
      'allowed_init_options_rate',
      'allowed_init_options_shape',
      'allowed_init_options_scale',
      'stanvars_datazz',
      'pstrarg',
      'initsi',
      'init_arguments',
      'init_data',
      'init_data_internal',
      'init_args_internal',
      'prior_data',
      'prior_data_internal',
      'prior_internal_args',
      'splitmvar_w2',
      'seed'
    )
    
    
    init_internal_args <- mget(init_internal_args_names)
    init_argument <- gsub("_prior_", "_init_", prior_argument)
    initial_out <-
      prepare_initials(init_argument = init_argument, init_internal_args)
  } else {
    initial_out <- NULL
  }
  
  
  
  
  return(
    list(
      prior_str_arg = prior_str_arg_out,
      lowerbound = lowerbound,
      upperbound = upperbound,
      stanvars_data = stanvars_data,
      initial_out = initial_out
    )
  )
}








#' An internal function to prepare initials for Bayesian SITAR growth curve 
#' model
#' 
#' For \code{univariate_by} and \code{multivariate} models (see [bsitar::bsitar()])
#' each argument is automatically matched with the sub model.
#'
#' @param init_argument A list containing the prior arguments specified in  
#' the [bsitar::bsitar()] function and then passed from the 
#' [bsitar::set_priors_initials()] function to the \code{prepare_priors}. 
#' 
#' @param init_internal_args An internal argument (as named list) specified in  
#' the [bsitar::bsitar()] function and then passed from the 
#' [bsitar::set_priors_initials()] function to the \code{prepare_priors}. 
#'
#' @return A list of initial values. 
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
prepare_initials <- function(init_argument,
                             init_internal_args) {
  
  
  
  ##############################################
  # Initiate non formalArgs()
  ##############################################
  seed <- NULL;
  initsi <- NULL;
  parm <- NULL;
  randomsi <- NULL;
  init_data <- NULL;
  prior_argument <- NULL;
  pstrarg <- NULL;
  resp_ <- NULL;
  allowed_init_options_beta <- NULL;
  nlpar <- NULL;
  stanvars_datazz <- NULL;
  splitmvar_w2 <- NULL;
  ii <- NULL;
  allowed_init_options_sd <- NULL;
  allowed_init_options_rate <- NULL;
  allowed_init_options_shape <- NULL;
  nabcrei <- NULL;
  N_J_all <- NULL;
  nys <- NULL;
  dpar <- NULL;
  ndparcov <- NULL;
  acorclass <- NULL;
  nrep_of_parms_p <- NULL;
  nrep_of_parms_q <- NULL;
  cortimeNlags <- NULL;
  gr_init_cor <- NULL;
  
  
  if (!is.null(init_internal_args)) {
    eout <- list2env(init_internal_args)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  if (!is.null(init_internal_args$init_argument)) {
    eout <- list2env(init_internal_args$init_argument)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  if (!is.null(init_internal_args$prior_data)) {
    eout <- list2env(init_internal_args$prior_data)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  if (!is.null(init_internal_args$prior_data_internal)) {
    eout <- list2env(init_internal_args$prior_data_internal)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  if (!is.null(init_internal_args$prior_internal_args)) {
    eout <- list2env(init_internal_args$prior_internal_args)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  
  if (!is.null(init_internal_args$init_data)) {
    eout <- list2env(init_internal_args$init_data)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  if (!is.null(init_internal_args$init_data_internal)) {
    eout <- list2env(init_internal_args$init_data_internal)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  
  
  set.seed(seed)
  
  # For standardized group level effects as implemented in the centerized 
  # parametrisation approach used in the brms package. 
  init_argument_z <- "r_init_z"
  
  if (initsi == "0") {
    assign(init_argument, "0")
  }
  if (initsi == "prior") {
    if(system.file(package='extraDistr') == "") {
      stop("For prior based initials (i.e., init = 'prior), 
           package 'extraDistr' is required. Please install 'extraDistr'"
      )
    }
    assign(init_argument, "prior")
    assign(init_argument_z, "prior")
  }
  
  check_form_0       <- paste0(parm, "_", 'form_0')
  nparcov            <- paste0(parm, "ncov")
  check_form_0_gr    <- paste0(parm, "_", 'form_0_gr')
  nparcov_gr         <- paste0(parm, "ncov_gr")
  check_sigma_form_0 <- paste0('sigma', "_", 'form_0')
  
  if(nparcov == 'nsigmacov') nparcov <- 'nsigmacov'
  
  abcrandomelements <-
    strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
  
  
  if (!is.null(init_data[[1]])) {
    eout <- list2env(init_data)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  
  
  
  abcrandomelements_c <- c()
  count_ <- 0
  for (abcrandomelements_i in abcrandomelements) {
    count_ <- count_ + 1
    if (exists(paste0(abcrandomelements_i, "n",  "cov_gr")) &
        length(ept(paste0(abcrandomelements_i, "n",  "cov_gr"))) > 0) {
      i__ <- ept(paste0(abcrandomelements_i, "n",  "cov_gr"))
      nb <- rep(paste0(abcrandomelements_i, "cov"), i__ - 1)
      nb <- c(nb, paste0(nb, 1:length(nb)))
      abcrandomelements_c <- c(abcrandomelements_c, nb)
    } else {
      abcrandomelements_c <- c(abcrandomelements_c, abcrandomelements_i)
    }
  }
  
  abcrandomelements <- abcrandomelements_c
  abcrandomelements <- unique(abcrandomelements)
  
  c_t_nabcri <- length(abcrandomelements)
  # nabcri <- length(c_t_nabcri)
  
  
  suffix <-
    strsplit(init_argument, "_")[[1]][length(strsplit(init_argument, "_")[[1]])]
  
  
  # define function to check validity of the initials options specified
  check_evalation_of_numeric_init_obj <-
    function(eit,
             check,
             x,
             pname_,
             dist,
             nlpar,
             class,
             allowed_init_options,
             splitmvar_w2) {
      whatin <-
        sub("=[^=]+$", "", splitmvar_w2[grepl(eit, splitmvar_w2)])
      const_msg <-
        paste0(
          " - a numeric value (e.g., 2) or a charater string such as",
          "\n",
          "xxx with xxx defined in the use-specified 'prior_data'",
          "\n",
          "argument e.g., prior_data = list(xxx = 2)"
        )
      if (!is.null(allowed_init_options)) {
        allowed_init_options <-
          paste0("random,",
                 " 0, ",
                 paste0(allowed_init_options, collapse = ", "))
        allowed_init_options <- paste0(" - ", allowed_init_options)
        const_msg <- paste0(allowed_init_options, "\n", const_msg)
      } else {
        const_msg <- paste0("random,", " 0, ", "Or ", const_msg)
      }
      err. <- FALSE
      tryCatch(
        expr = {
          out <- ept(eit)
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (eit == 'NULL' | eit == 'random')
        err. <- FALSE
      if (err.) {
        if (check == 'args') {
          if (class == 'b' | class == 'sd') {
            stop(
              "\nFor nlpar ",
              nlpar,
              ", class ",
              class,
              ", you have specified '",
              eit,
              "' as an initial argument",
              "\n" ,
              " But '",
              eit,
              "' is not found in the 'init_data_internal'",
              "\ n",
              " or use-specified 'init_data' argument",
              "\n ",
              " [see specified init argument: ",
              init_argument,
              "",
              "",
              "]",
              "\n" ,
              "Avilable  options are:" ,
              "\n" ,
              const_msg
            )
          } else if (class == 'sigma') {
            stop(
              "\nFor residual standard deviation parameter i.e., ",
              "class ",
              class,
              ", you have specified '",
              eit,
              "' as an initial argument",
              "\n" ,
              " But '",
              eit,
              "' is not found in the 'init_data_internal'",
              "\ n",
              " or use-specified 'init_data' argument",
              "\n ",
              " [see specified init argument: ",
              init_argument,
              "",
              "",
              "]",
              "\n" ,
              "Avilable  options are:" ,
              "\n" ,
              const_msg
            )
          } else if (class == '' &
                     grepl("dpar_", prior_argument) &
                     !grepl("dpar_cov", prior_argument)) {
            stop(
              "\nFor for distributional Intercept parameter i.e., ",
              "Intercept_sigma ",
              ", you have specified '",
              eit,
              "' as an initial argument",
              "\n" ,
              " But '",
              eit,
              "' is not found in the 'init_data_internal'",
              "\ n",
              " or use-specified 'init_data' argument",
              "\n ",
              " [see specified init argument: ",
              init_argument,
              "",
              "",
              "]",
              "\n" ,
              "Avilable  options are:" ,
              "\n" ,
              const_msg
            )
          }
          
        }
        if (check == 'dist') {
          stop(
            "For nlpar ",
            nlpar,
            ", class ",
            class,
            ", you have specified '",
            eit,
            "' as ",
            "\n" ,
            pname_,
            " for the ",
            dist,
            " distribution. But '",
            eit,
            "' is not found (check init_data!)",
            "\n" ,
            "Avilable options for the ",
            pname_,
            " parameter of " ,
            dist,
            " distribution ",
            "\n" ,
            "for nlpar ",
            nlpar,
            ", class ",
            class,
            " are:",
            "\n" ,
            const_msg
          )
        }
        rm(err.)
      }
    }
  
  
  # define function to create correlation matrix
  create_cor_mat <- function(n, cor = NULL) {
    n_elements <- n
    m <- diag(n_elements)
    m_upper <- m_lower <- matrix(0, n_elements, n_elements)
    nc <- n_elements * (n_elements - 1) / 2
    if (is.null(cor)) {
      x <- rep(0, nc)
    } else {
      x <- cor
      if (length(x) != nc)
        stop("length of correlation vector must be ",
             nc,
             ", but found ",
             length(x))
    }
    m_lower[lower.tri(m_lower, diag = FALSE)] <- x
    m_upper <- t(m_lower)
    M <- m_lower + m + m_upper
    M
  }
  
  
  # define function evaluation initials based on priors
  eval_prior_based_init <-
    function(dist,
             class,
             lowerbound,
             upperbound,
             length_args,
             ...) {
      if (dist == "normal") {
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_itb    <- paste0('', ")")
        z_replace_byn    <-
          paste0("extraDistr::rtnorm", "\\(", length_args, ", ")
        if (any(is.na(lowerbound)))
          lowerbound <- '-Inf'
        if (any(is.na(upperbound)))
          upperbound <- 'Inf'
        z_replace_byb   <-
          paste0(",",  paste(lowerbound, upperbound, sep = ",")   , ")")
        init_str_arg_out_init <-
          gsub(z_replace_itn, z_replace_byn, pstrarg)
        init_str_arg_out_init <-
          gsub(z_replace_itb, z_replace_byb, init_str_arg_out_init)
        init_str_arg_out_init <-
          gsub("\\s", "", init_str_arg_out_init)
      } else if (dist == "cauchy") {
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_itb    <- paste0('', ")")
        if (is.na(ept(lowerbound)) & is.na(ept(upperbound))) {
          z_replace_byn    <- paste0("rcauchy", "\\(", length_args, ", ")
          init_str_arg_out_init <-
            gsub(z_replace_itn, z_replace_byn, pstrarg)
        } else if (ept(lowerbound) == 0 & is.na(upperbound)) {
          z_replace_byn    <-
            paste0("extraDistr::rhcauchy", "\\(", length_args, ", ")
          init_str_arg_out_init <-
            gsub(z_replace_itn, z_replace_byn, pstrarg)
          init_str_arg_out_init    <-
            paste(strsplit(init_str_arg_out_init, ",")[[1]][-2], collapse = ",")
          init_str_arg_out_init <-
            gsub("\\s", "", init_str_arg_out_init)
        } else {
          stop(
            "For ",
            dist,
            " distribution prior based initials,",
            "\n ",
            " allowed options are unbounded distribution or half ",
            dist,
            " (i.e, lb = 0)",
            "\n ",
            " please check following initial argument: ",
            name_initialsi
          )
        }
      } else if (dist == "student_t") {
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_itb    <- paste0('', ")")
        if (is.na(ept(lowerbound)) & is.na(ept(upperbound))) {
          z_replace_byn    <-
            paste0("extraDistr::rlst", "\\(", length_args, ", ")
          init_str_arg_out_init <-
            gsub(z_replace_itn, z_replace_byn, pstrarg)
        } else if (ept(lowerbound) == 0 & is.na(upperbound)) {
          z_replace_byn    <-
            paste0("extraDistr::rht", "\\(", length_args, ", ")
          init_str_arg_out_init <-
            gsub(z_replace_itn, z_replace_byn, pstrarg)
          init_str_arg_out_init    <-
            paste(strsplit(init_str_arg_out_init, ",")[[1]][-3], collapse = ",")
          init_str_arg_out_init <-
            gsub("\\s", "", init_str_arg_out_init)
        } else {
          stop(
            "For ",
            dist,
            " distribution prior based initials,",
            "\n ",
            " allowed options are unbounded distribution or half ",
            dist,
            " (i.e, lb = 0)",
            "\n ",
            " please check following initial argument: ",
            name_initialsi
          )
        }
      } else if (dist == "gamma") {
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_byn    <-
          paste0("rgamma", "\\(", length_args, ", ")
        init_str_arg_out_init <-
          gsub(z_replace_itn, z_replace_byn, pstrarg)
      } else if (dist == "lognormal ") {
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_byn    <-
          paste0("rlnorm", "\\(", length_args, ", ")
        init_str_arg_out_init <-
          gsub(z_replace_itn, z_replace_byn, pstrarg)
      } else if (dist == "exponential") {
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_byn    <- paste0("rexp", "\\(", length_args, ", ")
        init_str_arg_out_init <-
          gsub(z_replace_itn, z_replace_byn, pstrarg)
      } else if (dist == "inv_gamma") {
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_byn    <-
          paste0("extraDistr::rinvgamma", "\\(", length_args, ", ")
        init_str_arg_out_init <-
          gsub(z_replace_itn, z_replace_byn, pstrarg)
      } else if (dist == "uniform") {
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_byn    <-
          paste0("runif", "\\(", length_args, ", ")
        init_str_arg_out_init <-
          gsub(z_replace_itn, z_replace_byn, pstrarg)
      } else if (dist == "xxxxxxxxx") {
        # placeholder for fututre expansion
      }
      
      if (class == 'cor' | class == 'rescor') {
        setlkjcorr <- function (K, eta = 1) {
          # set number of realizations
          n <- 100
          stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
          stopifnot(eta > 0)
          #if (K == 1) return(matrix(1, 1, 1))
          f <- function() {
            alpha <- eta + (K - 2) / 2
            r12 <- 2 * rbeta(1, alpha, alpha) - 1
            R <- matrix(0, K, K)
            R[1, 1] <- 1
            R[1, 2] <- r12
            R[2, 2] <- sqrt(1 - r12 ^ 2)
            if (K > 2)
              for (m in 2:(K - 1)) {
                alpha <- alpha - 0.5
                y <- rbeta(1, m / 2, alpha)
                z <- rnorm(m, 0, 1)
                z <- z / sqrt(crossprod(z)[1])
                R[1:m, m + 1] <- sqrt(y) * z
                R[m + 1, m + 1] <- sqrt(1 - y)
              }
            return(crossprod(R))
          }
          R <- replicate(n , f())
          if (dim(R)[3] == 1) {
            R <- R[, , 1]
          } else {
            R <- aperm(R, c(3, 1, 2))
          }
          return(R)
        }
        
        z_replace_itn    <- paste0(dist, "\\(")
        z_replace_byn    <-
          paste0("setlkjcorr", "\\(", length_args, ", ")
        init_str_arg_out_init <-
          gsub(z_replace_itn, z_replace_byn, pstrarg)
        Rcor <- ept(init_str_arg_out_init)
        c_r_list <- list()
        for (i in 1:dim(Rcor)[1])
          c_r_list[[i]] <-  Rcor[i:i, , ]
        Rcoravg <- Reduce("+", c_r_list) / length(c_r_list)
        outcor <- Rcoravg
      }
      
      if (class != 'cor' & class != 'rescor') {
        out  <- ept(init_str_arg_out_init)
        if (dist == "exponential")
          out <- 1 / out
      } else if (class == 'cor' | class == 'rescor') {
        out  <- outcor
      }
      
      out
    }
  
  
  
  ################################################
  # class b - beta
  if ((class == 'b' & suffix == 'beta') |
      class == 'b' & suffix == 'beta' & ept("sigma_dpar") == "sigma") {
    
    if(ept("sigma_dpar") == "sigma") {
      name_parm <- paste0(class, "_", parm, resp_)
    } else {
      name_parm <- paste0(class, resp_, "_", parm)
    }
    
    # name_parm <- paste0(class, resp_, "_", parm)
    suffix <- 'beta'
    
    allowed_init_options <- allowed_init_options_beta
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    lowerbound <- lowerbound
    upperbound <- upperbound
    
    out_list <- list_collect <- list()
    start_cnt <- 0
    for (name_initialsi in init_argument) {
      if (grepl("_cov", name_initialsi))
        eit_cov <- TRUE
      else
        eit_cov <- FALSE
      
      start_cnt <- start_cnt + 1
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            evaluated_init <- rep(0, nrep_of_parms)
          } else if (ept(name_initialsi) == 'lm') {
            if (ept(check_form_0)) {
              lm_gsubby <- paste0("lm", "_", nlpar, "_", "all", resp_)
            }
            if (!ept(check_form_0)) {
              if (!eit_cov)
                lm_gsubby <- paste0("lm", "_", nlpar, "", "", resp_)
              if (eit_cov)
                lm_gsubby <-
                  paste0("lm", "_", nlpar, "_", "cov", resp_)
            }
            
            evaluated_init <- ept(lm_gsubby) %>% unname()
          } else if (ept(name_initialsi) == 'ymean') {
            eit <- gsub("ymean",
                        paste0("ymean", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ymax') {
            eit <- gsub("ymax",
                        paste0("ymax", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ymaxs') {
            eit <- gsub("ymaxs",
                        paste0("ymaxs", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'median') {
            eit <- gsub("median",
                        paste0("median", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'bstart') {
            eit <- gsub("bstart",
                        paste0("bstart", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'cstart') {
            eit <- gsub("cstart",
                        paste0("cstart", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'dstart') {
            eit <- gsub("dstart",
                        paste0("dstart", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'estart') {
            eit <- gsub("estart",
                        paste0("estart", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
            
          } else if (ept(name_initialsi) == 'ymeanxmin') {
            eit <- gsub("ymeanxmin",
                        paste0("ymeanxmin", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ymeanxmax') {
            eit <- gsub("ymeanxmax",
                        paste0("ymeanxmax", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ymeanxmid') {
            eit <- gsub("ymeanxmid",
                        paste0("ymeanxmid", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'prior') {
            stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
            stanvname_cnt <- 0
            for (stanvname_i in stanvname_) {
              stanvname_cnt <- stanvname_cnt + 1
              stanvname_cnt_name <-
                stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
              assign(stanvname_cnt_name,
                     stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              length_args <-
                length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
            }
            evaluated_init <-
              eval_prior_based_init(
                dist = dist,
                class = class,
                lowerbound = lowerbound,
                upperbound = upperbound,
                length_args = length_args
              )
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = dist,
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              evaluated_init <- ept(name_initialsi)
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        }
      } else if (ept(name_initialsi) == 'NULL') {
        evaluated_init <- NULL
      }
      list_collect[[name_initialsi]] <- evaluated_init
    }
    
    if (!is.null(list_collect[[name_initialsi]])) {
      tempv <- list_collect %>% unlist() %>% unname()
      tempv <- as.numeric(tempv)
      out_list[[name_parm]] <- tempv
    } else {
      out_list[[name_parm]] <- NULL
    }
  }
  
  
  
  if (class == 'b' & suffix == 'beta' & ept("sigma_dpar") != "sigma") {
    if(ept(check_form_0)) {
      if(length(out_list[[name_parm]]) == 1) {
        out_list[[name_parm]] <- rep(out_list[[name_parm]], ept(nparcov))
      }
    }
  }
  
  
  # new addition on 27 5 2023
  # like a b c d e beta when ~ 0 +..., init for sigma are rep of sigma_init_beta 
  if (class == 'b' & suffix == 'beta' & ept("sigma_dpar") == "sigma") {
    if(ept(ept("init_argument")) != 'random' ) {
      if(ept("sigma_form_0")) {
        out_list[[name_parm]] <- rep(out_list[[name_parm]], ept(nparcov))
      } else if(!ept("sigma_form_0")) {
        if(eit_cov) {
          addcovsigma <- unlist(ept(ept(ept("init_argument"))))
          addcovsigma_n <- ept(nparcov) - 1
          if(addcovsigma_n == length(addcovsigma)) {
            addcovsigma <- addcovsigma
          } else if(length(addcovsigma) == 1) {
            addcovsigma <- rep(addcovsigma, addcovsigma_n)
          }
          if(!is.null(ept(nparcov))) out_list[[name_parm]] <-addcovsigma
        }
      } 
    } 
  } 
  
  
  ################################################
  # class sd - sd
  if (class == 'sd' & suffix == 'sd' & ept("sigma_dpar") != "sigma" |
      class == 'sd' & suffix == 'sd' & ept("sigma_dpar") == "sigma"
  ) {
    name_parm <- paste0('sd', "_", ii)
    suffix <- 'sd'
    
    if (dist == 'normal' |
        dist == 'cauchy' |
        dist == 'student_t' | dist == 'student_nu') {
      allowed_init_options <- allowed_init_options_sd
    } else if (dist == 'exponential') {
      allowed_init_options <- allowed_init_options_rate
    } else if (dist == 'gamma') {
      allowed_init_options <- allowed_init_options_shape
    } else {
      allowed_init_options <- NULL
    }
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    lowerbound <- 0
    upperbound <- upperbound
    
    out_list <- list_collect <- list()
    start_cnt <- 0
    for (name_initialsi in init_argument) {
      if (grepl("_cov", name_initialsi))
        eit_cov <- TRUE
      else
        eit_cov <- FALSE
      if (eit_cov)
        addcovname <- 'cov'
      else
        addcovname <- NULL
      
      start_cnt <- start_cnt + 1
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            evaluated_init <- rep(0, nrep_of_parms)
          } else if (ept(name_initialsi) == 'ysd') {
            eit <- gsub("ysd", paste0("ysd", resp_), ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ymad') {
            eit <- gsub("ymad",
                        paste0("ymad", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'lme_sd_a') {
            eit <-
              gsub("lme_sd_a",
                   paste0("lme_sd_a", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
            
            
          } else if (ept(name_initialsi) == 'ysdxmin') {
            eit <-
              gsub("ysdxmin",
                   paste0("ysdxmin", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysdxmax') {
            eit <-
              gsub("ysdxmax",
                   paste0("ysdxmax", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysdxmid') {
            eit <-
              gsub("ysdxmid",
                   paste0("ysdxmid", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysdxmidxmaxdiff') {
            eit <-
              gsub("ysdxmidxmaxdiff",
                   paste0("ysdxmidxmaxdiff", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'prior') {
            stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
            stanvname_cnt <- 0
            for (stanvname_i in stanvname_) {
              stanvname_cnt <- stanvname_cnt + 1
              stanvname_cnt_name <-
                stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
              assign(stanvname_cnt_name,
                     stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              length_args <-
                length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
            }
            evaluated_init <-
              eval_prior_based_init(
                dist = dist,
                class = class,
                lowerbound = lowerbound,
                upperbound = upperbound,
                length_args = length_args
              )
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = 'dis',
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              evaluated_init <- ept(name_initialsi)
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        }
      } else if (ept(name_initialsi) == 'NULL') {
        evaluated_init <- NULL
      }
      list_collect[[name_initialsi]] <- evaluated_init
    }
    
    if (!is.null(list_collect[[name_initialsi]])) {
      tempv <- list_collect %>% unlist() %>% unname()
      tempv <- as.numeric(tempv)
      if (ept(check_form_0_gr) & !is.null(ept(nparcov_gr))) {
        if (length(tempv) < length(ept(nparcov_gr)) + 1)
          tempv <- rep(tempv, length(ept(nparcov_gr)) + 1)
      }
      for (tempvi in 1:length(tempv)) {
        if (tempv[tempvi] == 0)
          tempv[tempvi] <- 1
      }
      if (dist == 'exponential')
        tempv <- 1 / tempv
      attr(tempv, 'names') <- paste0(parm, addcovname, ii)
      out_list[[name_parm]] <- tempv
    } else {
      out_list[[name_parm]] <- NULL
    }
  }
  
  
  
  
  
  ################################################
  # class cor - cor
  if (class == 'cor' & suffix == 'cor') {
    name_parm <- paste0('L', "_", ii)
    suffix <- 'cor'
    
    allowed_init_options <- NULL
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    if (!is.null(c_t_nabcri)) {
      NC_dims <- c_t_nabcri
    } else {
      NC_dims <- nabcrei
    }
    
    NC_cor_elements <- (NC_dims * (NC_dims - 1)) / 2
    
    lowerbound <- lowerbound
    upperbound <- upperbound
    
    out_list <- list_collect <- list()
    start_cnt <- 0
    for (name_initialsi in init_argument) {
      if (grepl("_cov", name_initialsi))
        eit_cov <- TRUE
      else
        eit_cov <- FALSE
      if (eit_cov)
        addcovname <- 'cov'
      else
        addcovname <- NULL
      start_cnt <- start_cnt + 1
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            L_elements     <- rep(0, NC_cor_elements)
            evaluated_init <- create_cor_mat(NC_dims, L_elements)
          } else if (ept(name_initialsi) == 'prior') {
            stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
            stanvname_cnt <- 0
            for (stanvname_i in stanvname_) {
              stanvname_cnt <- stanvname_cnt + 1
              stanvname_cnt_name <-
                stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
              assign(stanvname_cnt_name,
                     stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              length_args <-
                length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
            }
            evaluated_init <-
              eval_prior_based_init(
                dist = dist,
                class = class,
                lowerbound = lowerbound,
                upperbound = upperbound,
                length_args = NC_dims
              )
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = 'dis',
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              if (length(ept(name_initialsi)) == 1) {
                L_elements <- rep(ept(name_initialsi), NC_cor_elements)
              } else if (length(ept(name_initialsi)) == NC_cor_elements) {
                L_elements <- ept(name_initialsi)
              } else {
                stop(
                  "length of correlation vector must be ",
                  NC_cor_elements,
                  ", but found ",
                  length(ept(name_initialsi)),
                  "\n ",
                  " Please check the following init arg :",
                  'gr_init_cor'
                )
              }
              evaluated_init <- create_cor_mat(NC_dims, L_elements)
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        }
      } else if (ept(name_initialsi) == 'NULL') {
        evaluated_init <- NULL
      }
      list_collect[[name_initialsi]] <- evaluated_init
    }
    
    if (!is.null(list_collect[[name_initialsi]])) {
      tempv <- list_collect %>% unname()
      tempv <- tempv[[1]]
      colnames(tempv) <-
        rownames(tempv) <- paste0(abcrandomelements, addcovname, ii)
      out_list[[name_parm]] <- tempv
    } else {
      out_list[[name_parm]] <- NULL
    }
    
    ####### add "r_init_z"
    set_r_init_z <- TRUE
    
    if (!is.null(c_t_nabcri)) {
      nabcrei_z <- c_t_nabcri
    } else {
      nabcrei_z <- nabcrei
    }
    
    if (set_r_init_z) {
      name_parm <- paste0('z', "_", addcovname, ii)
      name_initialsi <- init_argument_z
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            evaluated_init <- matrix(0, nabcrei_z, N_J_all)
          } else if (ept(name_initialsi) == 'prior') {
            evaluated_init <-
              matrix(rnorm(nabcrei_z * N_J_all, 0, 1), nabcrei_z, N_J_all)
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = 'dis',
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                is.matrix(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              if (is.numeric(ept(name_initialsi)) &
                  length(ept(name_initialsi)) == 1) {
                z_std <- ept(name_initialsi) %>% as.numeric()
                if (z_std > 1 |
                    z_std < 0)
                  stop("sd for standardized matrix must be between 0 and 1")
                evaluated_init <-
                  matrix(rnorm(nabcrei_z * N_J_all, 0, z_std),
                         nabcrei_z,
                         N_J_all)
              } else if (is.matrix(ept(name_initialsi))) {
                evaluated_init <- ept(name_initialsi) %>% as.numeric()
                if (nrow(evaluated_init) != nabcrei_z &
                    nrow(evaluated_init) != N_J_all) {
                  stop(
                    "standardized matrix must have ",
                    nabcrei_z,
                    " rows and ",
                    N_J_all,
                    " columns"
                  )
                }
              } else {
                stop(
                  "initails for standardized matrix must be a single",
                  "\n",
                  " value of sd between 0 and 1",
                  "\n ",
                  " standardized matrix must have ",
                  nabcrei_z,
                  " rows and ",
                  N_J_all,
                  "columns",
                  "\n ",
                  " Please check the following init arg :",
                  'r_init_z'
                )
              }
              evaluated_init <- evaluated_init
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        } else if (ept(name_initialsi) == 'NULL') {
          evaluated_init <- NULL
        }
      }
      if (!is.null(evaluated_init)) {
        out_list[[name_parm]] <- evaluated_init
      }
    }
  }
  
  
  
  
  ################################################
  # class Lrescor - mvr rescor
  if (class == 'rescor' & suffix == 'rescor') {
    name_parm <- paste0('Lrescor', "_", ii)
    suffix <- 'rescor'
    
    allowed_init_options <- NULL
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    NC_dims         <- ept(nys) %>% as.numeric()
    NC_cor_elements <- (NC_dims * (NC_dims - 1)) / 2
    
    lowerbound <- lowerbound
    upperbound <- upperbound
    
    out_list <- list_collect <- list()
    start_cnt <- 0
    for (name_initialsi in init_argument) {
      if (grepl("_cov", name_initialsi))
        eit_cov <- TRUE
      else
        eit_cov <- FALSE
      start_cnt <- start_cnt + 1
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            L_elements     <- rep(0, NC_cor_elements)
            evaluated_init <- create_cor_mat(NC_dims, L_elements)
          } else if (ept(name_initialsi) == 'prior') {
            stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
            stanvname_cnt <- 0
            for (stanvname_i in stanvname_) {
              stanvname_cnt <- stanvname_cnt + 1
              stanvname_cnt_name <-
                stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
              assign(stanvname_cnt_name,
                     stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              length_args <-
                length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
            }
            evaluated_init <-
              eval_prior_based_init(
                dist = dist,
                class = class,
                lowerbound = lowerbound,
                upperbound = upperbound,
                length_args = NC_dims
              )
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = 'dis',
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              if (length(ept(name_initialsi)) == 1) {
                L_elements <-
                  rep(ept(name_initialsi), NC_cor_elements) %>% as.numeric()
              } else if (length(ept(name_initialsi)) == NC_cor_elements) {
                L_elements <- ept(name_initialsi) %>% as.numeric()
              } else {
                stop(
                  "length of correlation vector must be ",
                  NC_cor_elements,
                  ", but found ",
                  length(ept(name_initialsi)),
                  "\n ",
                  " Please check the following init arg :",
                  'gr_init_cor'
                )
              }
              evaluated_init <- create_cor_mat(NC_dims, L_elements)
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        }
      } else if (ept(name_initialsi) == 'NULL') {
        evaluated_init <- NULL
      }
      list_collect[[name_initialsi]] <- evaluated_init
    }
    
    if (!is.null(list_collect[[name_initialsi]])) {
      tempv <- list_collect %>% unname()
      tempv <- tempv[[1]]
      out_list[[name_parm]] <- tempv
    } else {
      out_list[[name_parm]] <- NULL
    }
  }
  
  ################################################
  # class sigma
  if (class == 'sigma' & suffix == 'sigma') {
    name_parm <- paste0('sigma', resp_)
    suffix <- 'sigma'
    
    if (dist == 'normal' |
        dist == 'cauchy' |
        dist == 'student_t' | dist == 'student_nu') {
      allowed_init_options <- allowed_init_options_sd
    } else if (dist == 'exponential') {
      allowed_init_options <- allowed_init_options_rate
    } else if (dist == 'gamma') {
      allowed_init_options <- allowed_init_options_shape
    } else {
      allowed_init_options <- NULL
    }
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    lowerbound <- 0
    upperbound <- upperbound
    
    out_list <- list_collect <- list()
    start_cnt <- 0
    for (name_initialsi in init_argument) {
      if (grepl("_cov", name_initialsi))
        eit_cov <- TRUE
      else
        eit_cov <- FALSE
      start_cnt <- start_cnt + 1
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            evaluated_init <- ept(name_initialsi) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysd') {
            eit <- gsub("ysd", paste0("ysd", resp_), ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ymad') {
            eit <- gsub("ymad",
                        paste0("ymad", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'lme_rsd') {
            eit <-
              gsub("lme_rsd",
                   paste0("lme_rsd", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'lm_rsd') {
            eit <- gsub("lm_rsd",
                        paste0("lm_rsd", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysdxmin') {
            eit <-
              gsub("ysdxmin",
                   paste0("ysdxmin", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysdxmax') {
            eit <-
              gsub("ysdxmax",
                   paste0("ysdxmax", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysdxmid') {
            eit <-
              gsub("ysdxmid",
                   paste0("ysdxmid", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysdxmidxmaxdiff') {
            eit <-
              gsub("ysdxmidxmaxdiff",
                   paste0("ysdxmidxmaxdiff", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'prior') {
            stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
            stanvname_cnt <- 0
            for (stanvname_i in stanvname_) {
              stanvname_cnt <- stanvname_cnt + 1
              stanvname_cnt_name <-
                stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
              assign(stanvname_cnt_name,
                     stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              length_args <-
                length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
            }
            evaluated_init <-
              eval_prior_based_init(
                dist = dist,
                class = class,
                lowerbound = lowerbound,
                upperbound = upperbound,
                length_args = length_args
              )
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = 'dis',
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              evaluated_init <- ept(name_initialsi) %>% as.numeric()
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        }
      } else if (ept(name_initialsi) == 'NULL') {
        evaluated_init <- NULL
      }
      list_collect[[name_initialsi]] <- evaluated_init
    }
    
    if (!is.null(list_collect[[name_initialsi]])) {
      tempv <- list_collect %>% unlist() %>% unname()
      tempv <- as.numeric(tempv)
      for (tempvi in 1:length(tempv)) {
        if (tempv[tempvi] == 0)
          tempv[tempvi] <- 1
      }
      if (dist == 'exponential')
        tempv <- 1 / tempv
      out_list[[name_parm]] <- tempv
    } else {
      out_list[[name_parm]] <- NULL
    }
    
  }
  
  
  # sigma modeling via dpar arguments that allows for inclusion of covariates
  
  if (class == '' &
      dpar != "" & suffix == 'sigma' & !ept(check_sigma_form_0)) {
    if (grepl("dpar_init_sigma", init_argument)) {
      name_parm <- paste0('Intercept_sigma', resp_)
      lowerbound <- 0
      upperbound <- upperbound
    }
    
    if (grepl("dpar_cov_init_sigma", init_argument)) {
      name_parm <- paste0('b_sigma', resp_)
      lowerbound <- lowerbound
      upperbound <- upperbound
    }
    
    name_parm <- paste0(name_parm, resp_)
    
    suffix <- 'sigma'
    
    
    if (grepl('Intercept_sigma', name_parm)) {
      if (dist == 'normal' |
          dist == 'cauchy' |
          dist == 'student_t' | dist == 'student_nu') {
        allowed_init_options <- allowed_init_options_sd
      } else if (dist == 'exponential') {
        allowed_init_options <- allowed_init_options_rate
      } else if (dist == 'gamma') {
        allowed_init_options <- allowed_init_options_shape
      } else {
        allowed_init_options <- NULL
      }
    } else if (!grepl('Intercept_sigma', name_parm)) {
      allowed_init_options <- NULL
    }
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    
    out_list <- list_collect <- list()
    start_cnt <- 0
    for (name_initialsi in init_argument) {
      if (grepl("_cov", name_initialsi))
        eit_cov <- TRUE
      else
        eit_cov <- FALSE
      start_cnt <- start_cnt + 1
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            evaluated_init <- ept(name_initialsi) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ysd') {
            eit <- gsub("ysd", paste0("ysd", resp_), ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'ymad') {
            eit <- gsub("ymad",
                        paste0("ymad", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'lme_rsd') {
            eit <-
              gsub("lme_rsd",
                   paste0("lme_rsd", resp_),
                   ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'lm_rsd') {
            eit <- gsub("lm_rsd",
                        paste0("lm_rsd", resp_),
                        ept(name_initialsi))
            evaluated_init <- ept(eit) %>% as.numeric()
          } else if (ept(name_initialsi) == 'prior') {
            stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
            stanvname_cnt <- 0
            for (stanvname_i in stanvname_) {
              stanvname_cnt <- stanvname_cnt + 1
              stanvname_cnt_name <-
                stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
              assign(stanvname_cnt_name,
                     stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              length_args <-
                length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
            }
            evaluated_init <-
              eval_prior_based_init(
                dist = dist,
                class = class,
                lowerbound = lowerbound,
                upperbound = upperbound,
                length_args = length_args
              )
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = 'dis',
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              evaluated_init <- ept(name_initialsi) %>% as.numeric()
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        }
      } else if (ept(name_initialsi) == 'NULL') {
        evaluated_init <- NULL
      }
      list_collect[[name_initialsi]] <- evaluated_init
    }
    
    if (!is.null(list_collect[[name_initialsi]])) {
      tempv <- list_collect %>% unlist() %>% unname()
      tempv <- as.numeric(tempv)
      if (grepl('Intercept_sigma', name_parm) &
          tempv == 0)
        tempv <- 1
      if (dist == 'exponential')
        tempv <- 1 / tempv
      out_list[[name_parm]] <- tempv
    } else {
      out_list[[name_parm]] <- NULL
    }
    
  }
  
  
  
  if (class == '' &
      dpar != "" & suffix == 'sigma' & ept(check_sigma_form_0)) {
    name_parm <- paste0('b_sigma', resp_)
    
    suffix <- 'sigma'
    
    covplus1 <- ndparcov + 1
    
    allowed_init_options <- NULL
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    lowerbound <- lowerbound
    upperbound <- upperbound
    
    out_list <- list_collect <- list()
    start_cnt <- 0
    for (name_initialsi in init_argument) {
      if (grepl("_cov", name_initialsi))
        eit_cov <- TRUE
      else
        eit_cov <- FALSE
      start_cnt <- start_cnt + 1
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            evaluated_init <- rep(0, nrep_of_parms)
          } else if (ept(name_initialsi) == 'lm_bsigma') {
            if (ept(check_form_0)) {
              lm_gsubby <- paste0("lm_bsigma", "_", nlpar, "_", "all", resp_)
            }
            if (!ept(check_form_0)) {
              if (!eit_cov)
                lm_gsubby <-
                  paste0("lm_bsigma", "_", nlpar, "", "", resp_)
              if (eit_cov)
                lm_gsubby <-
                  paste0("lm_bsigma", "_", nlpar, "_", "cov", resp_)
            }
            evaluated_init <- ept(lm_gsubby) %>% unname()
          } else if (ept(name_initialsi) == 'coef_bsigma') {
            evaluated_init <- ept('coef_bsigma')
          } else if (ept(name_initialsi) == 'prior') {
            stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
            stanvname_cnt <- 0
            for (stanvname_i in stanvname_) {
              stanvname_cnt <- stanvname_cnt + 1
              stanvname_cnt_name <-
                stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
              assign(stanvname_cnt_name,
                     stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              length_args <-
                length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
            }
            evaluated_init <-
              eval_prior_based_init(
                dist = dist,
                class = class,
                lowerbound = lowerbound,
                upperbound = upperbound,
                length_args = length_args
              )
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = 'dis',
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              evaluated_init <- ept(name_initialsi) %>% as.numeric()
              if (length(evaluated_init) < covplus1)
                evaluated_init <- rep(evaluated_init, covplus1)
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        }
      } else if (ept(name_initialsi) == 'NULL') {
        evaluated_init <- NULL
      }
      list_collect[[name_initialsi]] <- evaluated_init
    }
    
    if (!is.null(list_collect[[name_initialsi]])) {
      tempv <- list_collect %>% unlist() %>% unname()
      tempv <- as.numeric(tempv)
      out_list[[name_parm]] <- tempv
    } else {
      out_list[[name_parm]] <- NULL
    }
    
  }
  
  
  # autocorrelation parameters (ar, ma, and arms)
  if (class == 'ar' | class == 'ma' & suffix == 'acor') {
    if (acorclass == 'arma') {
      init_argument <- rep(init_argument, 2)
      name_parm_s <- c("ar", "ma")
    } else {
      name_parm_s <- class
    }
    
    name_parm_s <- paste0(name_parm_s, resp_)
    
    allowed_init_options <- NULL
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    lowerbound <- lowerbound
    upperbound <- upperbound
    
    suffix <- 'acor'
    
    out_list <- list()
    for (name_parm in name_parm_s) {
      
      if(grepl("ar", name_parm)) nrep_of_parms <- nrep_of_parms_p
      if(grepl("ma", name_parm)) nrep_of_parms <- nrep_of_parms_q
      
      list_collect <- list()
      start_cnt <- 0
      for (name_initialsi in init_argument) {
        if (grepl("_cov", name_initialsi))
          eit_cov <- TRUE
        else
          eit_cov <- FALSE
        start_cnt <- start_cnt + 1
        if (ept(name_initialsi) != 'NULL') {
          if (ept(name_initialsi) != 'random') {
            if (ept(name_initialsi) == '0') {
              evaluated_init <- rep(0, nrep_of_parms)
            } else if (ept(name_initialsi) == 'prior') {
              stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
              stanvname_cnt <- 0
              for (stanvname_i in stanvname_) {
                stanvname_cnt <- stanvname_cnt + 1
                stanvname_cnt_name <-
                  stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
                assign(stanvname_cnt_name,
                       stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
                length_args <-
                  length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              }
              evaluated_init <-
                eval_prior_based_init(
                  dist = dist,
                  class = class,
                  lowerbound = lowerbound,
                  upperbound = upperbound,
                  length_args = length_args
                )
            } else {
              check_evalation_of_numeric_init_obj(
                ept(name_initialsi),
                check = 'args',
                x = name_initialsi,
                pname_ = 'xxx',
                dist = 'dis',
                nlpar = parm,
                class = class,
                allowed_init_options = allowed_init_options,
                splitmvar_w2 = splitmvar_w2
              )
              name_initialsi <- ept(name_initialsi)
              if (is.numeric(ept(name_initialsi)) |
                  !is.null(ept(name_initialsi))) {
                evaluated_init <- ept(name_initialsi) %>% as.numeric()
                if (evaluated_init > 1 |
                    evaluated_init < -1)
                  stop("initials for autocorrelation must be between -1 and 1")
              }
            }
          } else if (ept(name_initialsi) == 'random') {
            evaluated_init <- NULL
          }
        } else if (ept(name_initialsi) == 'NULL') {
          evaluated_init <- NULL
        }
        list_collect[[name_initialsi]] <- evaluated_init
      }
      if (!is.null(list_collect[[name_initialsi]])) {
        tempv <- list_collect %>% unlist() %>% unname()
        tempv <- as.numeric(tempv)
        if(length(tempv) == 1 & nrep_of_parms == 1) {
          tempv <- tempv
        } else if(length(tempv) == 1 & nrep_of_parms > 1) {
          tempv <- rep(tempv, nrep_of_parms)
        } else if(length(tempv) != nrep_of_parms) {
          stop("Length of initials should be either 1 of equal to dims")
        }
        out_list[[name_parm]] <- tempv
      } else {
        out_list[[name_parm]] <- NULL
      }
    }
  }
  
  
  # autocorrelation parameters (unst) -> similar to cor / lrescor
  
  if (class == 'Lcortime' & suffix == 'acor') {
    
    # name_parm <- paste0('Cortime', "_", ii)
    
    name_parm_s <- class
    name_parm_s <- paste0(name_parm_s, resp_)
    
    name_parm <- name_parm_s # paste0('Lcortime', "_", ii)
    
    suffix <- 'acor'
    
    allowed_init_options <- NULL
    
    if(!exists('allowed_init_options')) allowed_init_options <- NULL
    
    NC_dims         <- ept(cortimeNlags) %>% as.numeric()
    NC_cor_elements <- (NC_dims * (NC_dims - 1)) / 2
    
    lowerbound <- lowerbound
    upperbound <- upperbound
    
    out_list <- list_collect <- list()
    start_cnt <- 0
    for (name_initialsi in init_argument) {
      if (grepl("_cov", name_initialsi))
        eit_cov <- TRUE
      else
        eit_cov <- FALSE
      start_cnt <- start_cnt + 1
      if (ept(name_initialsi) != 'NULL') {
        if (ept(name_initialsi) != 'random') {
          if (ept(name_initialsi) == '0') {
            L_elements     <- rep(0, NC_cor_elements)
            evaluated_init <- create_cor_mat(NC_dims, L_elements)
          } else if (ept(name_initialsi) == 'prior') {
            stanvname_ <- unique(unlist(lapply(stanvars_datazz, names)))
            stanvname_cnt <- 0
            for (stanvname_i in stanvname_) {
              stanvname_cnt <- stanvname_cnt + 1
              stanvname_cnt_name <-
                stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$name
              assign(stanvname_cnt_name,
                     stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
              length_args <-
                length(stanvars_datazz[[stanvname_cnt]][[stanvname_i]]$sdata)
            }
            evaluated_init <-
              eval_prior_based_init(
                dist = dist,
                class = class,
                lowerbound = lowerbound,
                upperbound = upperbound,
                length_args = NC_dims
              )
          } else {
            check_evalation_of_numeric_init_obj(
              ept(name_initialsi),
              check = 'args',
              x = name_initialsi,
              pname_ = 'xxx',
              dist = 'dis',
              nlpar = parm,
              class = class,
              allowed_init_options = allowed_init_options,
              splitmvar_w2 = splitmvar_w2
            )
            name_initialsi <- ept(name_initialsi)
            if (is.numeric(ept(name_initialsi)) |
                !is.null(ept(name_initialsi))) {
              if (length(ept(name_initialsi)) == 1) {
                L_elements <-
                  rep(ept(name_initialsi), NC_cor_elements) %>% as.numeric()
              } else if (length(ept(name_initialsi)) == NC_cor_elements) {
                L_elements <- ept(name_initialsi) %>% as.numeric()
              } else {
                stop(
                  "length of correlation vector must be ",
                  NC_cor_elements,
                  ", but found ",
                  length(ept(name_initialsi)),
                  "\n ",
                  " Please check the following init arg :",
                  'gr_init_cor'
                )
              }
              evaluated_init <- create_cor_mat(NC_dims, L_elements)
            }
          }
        } else if (ept(name_initialsi) == 'random') {
          evaluated_init <- NULL
        }
      } else if (ept(name_initialsi) == 'NULL') {
        evaluated_init <- NULL
      }
      list_collect[[name_initialsi]] <- evaluated_init
    }
    
    if (!is.null(list_collect[[name_initialsi]])) {
      tempv <- list_collect %>% unname()
      tempv <- tempv[[1]]
      out_list[[name_parm]] <- tempv
    } else {
      out_list[[name_parm]] <- NULL
    }
  }
  
  
  
  # standardized group level effects (part of centred parametrisation)
  
  if (any(grepl("_init_sd", init_argument)) &
      gr_init_cor == "NULL") {
    nabcrei_z <- 1
    ####### add "r_init_z"
    set_r_init_z <- TRUE
    if (set_r_init_z) {
      if (!is.null(ept(nparcov_gr))) {
        # ept(check_form_0_gr) &
        nabcrei_z <- nabcrei_z * ept(nparcov_gr)
      }
      for (name_initialsi in init_argument_z) {
        name_parm <-
          paste0('z',
                 "_",
                 strsplit(init_argument, "_")[[1]][1],
                 "",
                 addcovname,
                 ii)
        if (ept(name_initialsi) != 'NULL') {
          if (ept(name_initialsi) != 'random') {
            if (ept(name_initialsi) == '0') {
              evaluated_init <- matrix(0, nabcrei_z, N_J_all)
            } else if (ept(name_initialsi) == 'prior') {
              evaluated_init <-
                matrix(rnorm(nabcrei_z * N_J_all, 0, 1),
                       nabcrei_z,
                       N_J_all)
            } else {
              check_evalation_of_numeric_init_obj(
                ept(name_initialsi),
                check = 'args',
                x = name_initialsi,
                pname_ = 'xxx',
                dist = 'dis',
                nlpar = parm,
                class = class,
                allowed_init_options = allowed_init_options,
                splitmvar_w2 = splitmvar_w2
              )
              name_initialsi <- ept(name_initialsi)
              if (is.numeric(ept(name_initialsi)) |
                  is.matrix(ept(name_initialsi)) |
                  !is.null(ept(name_initialsi))) {
                if (is.numeric(ept(name_initialsi)) &
                    length(ept(name_initialsi)) == 1) {
                  z_std <- ept(name_initialsi) %>% as.numeric()
                  if (z_std > 1 |
                      z_std < 0)
                    stop("sd for standardized matrix must be between 0 and 1")
                  evaluated_init <-
                    matrix(rnorm(nabcrei_z * N_J_all, 0, z_std),
                           nabcrei_z,
                           N_J_all)
                } else if (is.matrix(ept(name_initialsi))) {
                  evaluated_init <- ept(name_initialsi) %>% as.numeric()
                  if (nrow(evaluated_init) != nabcrei_z &
                      nrow(evaluated_init) != N_J_all) {
                    stop(
                      "standardized matrix must have ",
                      nabcrei_z,
                      " rows and ",
                      N_J_all,
                      " columns"
                    )
                  }
                } else {
                  stop(
                    "initails for standardized matrix must be a",
                    "\n",
                    " single value of sd between 0 and 1",
                    "\n ",
                    " standardized matrix must have ",
                    nabcrei_z,
                    " rows and ",
                    N_J_all,
                    "columns",
                    "\n ",
                    " Please check the following init arg :",
                    'r_init_z'
                  )
                }
                evaluated_init <- evaluated_init
              }
            }
          } else if (ept(name_initialsi) == 'random') {
            evaluated_init <- NULL
          } else if (ept(name_initialsi) == 'NULL') {
            evaluated_init <- NULL
          }
        }
        if (!is.null(evaluated_init)) {
          out_list[[name_parm]] <- evaluated_init
        }
      }
    }
  }
  
  if (!exists('out_list'))
    out_list <- NULL
  
  out_list
}






# optimize_model.bgmfit here because it caused error numeric in the latest
# versions. Following dropped (full optimize_model.bgmfit saved i folder 4)

# @export optimize_model.bgmfit
# @export
# @examples
# @rdname optimize_model.bgmfit
# @export
# optimize_model <- function(model, ...) {
#  UseMethod("optimize_model")
# }




#'Optimize model
#'
#'@param model An object of class \code{bgmfit}.
#'
#'@param newdata An optional \code{data.frame} to be used when optimizing the
#'  model. If \code{NULL} (default), the same same data used for original model
#'  fit is used. Note that data-dependent default priors will not be updated
#'  automatically.
#'
#'@param optimize_df A vector of integers specifying the degree of freedom
#'  (\code{df}) values to be updated. If \code{NULL} (default), the \code{df} is
#'  taken from the original model. For \code{univariate-by-sungroup} and
#'  \code{multivariate} models (see [bsitar::bsitar()] for details),
#'  \code{optimize_df} can be a single integer (e.g., \code{optimize_df = 4}) or
#'  a list (e.g., \code{optimize_df = list(4,5)}). For optimization over
#'  different \code{df}, say for example \code{df} 4 and \code{df} 5 for
#'  univariate model, the corresponding code is \code{optimize_df = list(4,5)}.
#'  For a multivariate model fit to two outcomes with different \code{df}, the
#'  optimization over \code{df} 4 and \code{df} 5 for the first sub model, and
#'  \code{df} 5 and \code{df} 6 for the second
#'  sub model, the corresponding \code{optimize_df} code is \code{optimize_df =
#'  list(list(4,5), list(5,6))} i.e, a list of lists.
#'
#'@param optimize_x A vector specifying the transformations of predictor
#'  variable (i.e., \code{xvar} which is typically \code{age}). The option are
#'  \code{NULL}, \code{log}, \code{sqrt}, and their combinations. Note that user
#'  need not to enclose these options in a single or double quotes as they are
#'  take care of internally. The default setting is to explore all possible
#'  combination i.e., \code{optimize_x = list(NULL, log,  sqrt)}. Similar to the
#'  \code{optimize_df}, user can specify different \code{optimize_x} for
#'  \code{univariate-by-sungroup} and \code{multivariate} sub models.
#'
#'@param optimize_y A vector specifying the transformations of the response
#'  variable (i.e., \code{yvar} which is typically repeated height
#'  measurements). The approach and options available for \code{optimize_y} are
#'  identical to those described above for the \code{optimize_x}.
#'
#'@param exclude_default_funs A logical to indicate whether transformations for
#'  (\code{xvar} and \code{yvar}) used in the original model fit should be
#'  excluded. If \code{TRUE} (default), the \code{xvar} and \code{yvar}
#'  transformations specified for the original model fit are excluded from the
#'  \code{optimize_x} and \code{optimize_y}. From example, if original model is
#'  fit with \code{xvar = log} and \code{yvar = NULL}, then \code{optimize_x} is
#'  translated into \code{optimize_x = list(NULL, sqrt)} and \code{optimize_y}
#'  as \code{optimize_y = list(log, sqrt)}.
#'
#'@param add_fit_criteria An optional (default \code{NULL}) indicator to add fit
#'  criteria to returned model fit. Options are \code{loo} and \code{waic}. 
#'  Please see [brms::add_criterion()] for details.
#'
#'@param add_bayes_R An optional (default \code{NULL}) to add Bayesian R
#'  square. To estimate \code{bayes_R2}, please use 
#'  \code{add_bayes_R = bayes_R2}
#'
#'@param byresp A logical (default \code{FALSE}) to indicate if response wise
#'  fit criteria to be calculated. This argument is evaluated only for the
#'  \code{multivariate} model for which user can select whether to get joint
#'  calculation of point wise log likelihood or response specific. For,
#'  \code{univariate-by-subgroup} model, the only option available is to
#'  calculate separate point wise log likelihood for each sub-model.
#'
#'@param digits An integer to set the number of decimal places.
#'
#'@param cores The number of cores to used in processing (default \code{1}).
#' This passed to the [brms::add_criterion()].
#' 
#' @param verbose A logical (default \code{FALSE}) to indicate whether to print
#' information.
#'
#'@param envir A environment. Default is \code{glonalenv()}.
#'
#'@param ... Other arguments passed to \code{\link{update_model}}.
#'
#'@return A list containing the optimized models of class \code{bgmfit}, and the
#'  the summary statistics if \code{add_fit_criteria} and/or
#'  \code{add_bayes_R} are specified.
#'  
#'@seealso [brms::add_criterion()]
#'
#'@author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#'@importFrom loo pareto_k_table
#' 
#' @keywords internal
#' @noRd
#' 
optimize_model.bgmfit <- function(model,
                                  newdata = NULL,
                                  optimize_df = NULL,
                                  optimize_x = list(NULL, log,  sqrt),
                                  optimize_y = list(NULL, log,  sqrt),
                                  exclude_default_funs = TRUE,
                                  add_fit_criteria = c('waic'),
                                  add_bayes_R = c('bayes_R2'),
                                  byresp = FALSE,
                                  digits = 2,
                                  cores = 1,
                                  verbose = FALSE,
                                  envir = globalenv(),
                                  ...) {
  
  
  ##############################################
  # Initiate non formalArgs()
  ##############################################
  outcome <- NULL;
  xfun <- NULL; 
  yfun <- NULL;
  Parameter <- NULL;
  Estimate <- NULL;
  . <- NULL;
  
  if (is.null(newdata)) {
    newdata <- model$model_info$bgmfit.data
  } else {
    newdata <- newdata
  }
  
  if(!is.null(optimize_x)) {
    if(!is.list(optimize_x)) stop("argument 'optimize_x' must be a list")
  }
  
  if(!is.null(optimize_y)) {
    if(!is.list(optimize_y)) stop("argument 'optimize_y' must be a list")
  }
  
  o <-
    post_processing_checks(model = model,
                           xcall = match.call(),
                           resp = NULL,
                           envir = envir,
                           deriv = 0)
  
  call_o <- match.call()
  call_o_args <- as.list(call_o)[-1]
  
  args_o <- as.list(model$model_info$call.full.bgmfit)[-1]
  args_o_dots_ <- list(...)
  if (length(args_o_dots_) > 0) {
    for (i in names(args_o_dots_)) {
      args_o[[i]] <- args_o_dots_[[i]]
    }
  }
  
  
  
  # This to evaluate T/F to TRUE/FALSE
  for (i in names(args_o)) {
    if (is.symbol(args_o[[i]])) {
      if (args_o[[i]] == "T")
        args_o[[i]] <- eval(args_o[[i]])
      if (args_o[[i]] == "F")
        args_o[[i]] <- eval(args_o[[i]])
    }
  }
  
  for (add_fit_criteriai in add_fit_criteria) {
    if (!add_fit_criteriai %in% c("loo", "waic")) {
      stop("only loo and waic criteria are supported")
    }
  }
  
  for (bayes_Ri in add_bayes_R) {
    if (!bayes_Ri %in% c("bayes_R2")) {
      stop("only bayes_R2 as R square measure is supported")
    }
  }
  
  
  # Not must for expose_function to be true
  if (!args_o$expose_function) {
    if (!is.null(add_fit_criteria) | !is.null(add_bayes_R)) {
      # stop(
      #   "Argument expose_function must be set to TRUE when ",
      #   "\n ",
      #   " adding fit criteria and bayes_R2"
      # )
    }
  }
  
  
  
  
  get_args_opt <- function(xo) {
    get_within_fist_last_paranthesese <- function(x__) {
      x__ <- sub('\\(', '[', x__)
      x__ <- sub("\\)([^)]*)$", "]\\1", x__)
      x__ <-
        gsub("[\\[\\]]", "", regmatches(x__, gregexpr("\\[.*?\\]", x__))[[1]])
      x__ <- gsub("\\[|\\]", "", x__)
      x__
    }
    gsub_comma_within_paranthesese <-
      function(x__, replace_comma_by) {
        tt <-
          gsub("[\\(\\)]", "", regmatches(x__, gregexpr("\\(.*?\\)", x__))[[1]])
        tt2 <- gsub(",", replace_comma_by, tt, fixed = T)
        j <- 0
        for (i in tt) {
          j <- j + 1
          x__ <- gsub(tt[j], tt2[j], x__, fixed = T)
        }
        x__
      }
    
    xxo <- gsub("[[:space:]]", "", xo)
    
    xxo_g <- gsub('\"', "", xxo)
    xxo_g2 <- 
      grepl(
        "[-]?[0-9]+[.]?[0-9]*|[-]?[0-9]+[L]?|[-]?[0-9]+[.]?[0-9]*[eE][0-9]+", 
        xxo_g)
    
    if(any(xxo_g2)) xxo_g3 <- TRUE else xxo_g3 <- FALSE
    
    numeric_dx <- xxo_g3
    
    if (xxo != "NULL" & xxo != "\"NULL\"" & !numeric_dx) {
      xxo <- get_within_fist_last_paranthesese(xxo)
      xxo <- gsub_comma_within_paranthesese(xxo, "_comma_")
      xxo <- strsplit(xxo, ",")[[1]]
      xxo <- gsub("_comma_" , ",", xxo)
      xxo <- gsub('\"', "", xxo)
    } else {
      xxo <- xxo
      xxo <- gsub('\"', "", xxo)
    }
    xxo
  }
  
  optimize_df <- get_args_opt(deparse(substitute(optimize_df)))
  optimize_x  <- get_args_opt(deparse(substitute(optimize_x)))
  optimize_y  <- get_args_opt(deparse(substitute(optimize_y)))
  
  if (exclude_default_funs) {
    optimize_x <- optimize_x[!optimize_x %in% model$model_info$xfuns]
    optimize_y <-
      optimize_y[!optimize_y %in% model$model_info$xfuns]
    if (identical(optimize_x, character(0)))
      optimize_x <- "NULL"
    if (identical(optimize_y, character(0)))
      optimize_y <- "NULL"
  }
  
  optimize_df_x_y <-
    expand.grid(optimize_df, optimize_x, optimize_y)
  
  colnames(optimize_df_x_y) <- c("df", "xfun", "yfun")
  
  add_summary_waic <- NULL
  Count <- Est.Error <- Inference <- Min..n_eff <- where <- NULL
  Min.n_eff <- Percent <- Proportion <- Range <- SE <- NULL
  
  
  combine_summaries <- function(model_list, summary_obj) {
    ic = 0
    list_c <- list()
    for (model_listi in 1:length(model_list)) {
      if (!is.null(model_list[[model_listi]][[summary_obj]])) {
        ic <- ic + 1
        list_c[[ic]] <- model_list[[model_listi]][[summary_obj]]
      }
      summary_of_obj <-
        list_c %>% do.call(rbind, .) %>% data.frame()
    }
    if (nrow(summary_of_obj) < 1)
      summary_of_obj <- NULL
    summary_of_obj
  }
  
  # resp = NULL is only used as a placeholder that too only for multivariate
  # if NULL, then combined log likelihood used for multivariate model
  # if anything else e.g., resp = 'NULL' or anything, '
  # then separate likelihood for responses
  
  add_citeria_fun <- function(fit,
                              add_fit_criteria = NULL,
                              add_bayes_R = NULL,
                              resp = NULL,
                              digits = 2,
                              df,
                              xfun_print,
                              yfun_print,
                              ...) {
    
    if (!fit$model_info$call.full.bgmfit$expose_function) {
      assign(o[[1]], fit$model_info[['exefuns']][[o[[1]]]], envir = envir)
    }
    
    if (!is.null(add_fit_criteria)) {
      what_ <- paste(add_fit_criteria, collapse = ", ")
      message(" Adding", " ", what_, " ", "...")
      if(verbose) cat("\n")
      if (is.na(fit$model_info$univariate_by) |
          !fit$model_info$multivariate) {
        if (!fit$model_info$multivariate) {
          suppressWarnings(fit <- brms::add_criterion(fit,
                                                      add_fit_criteria, 
                                                      cores = cores))
        }
        if (fit$model_info$multivariate) {
          if (is.null(resp)) {
            suppressWarnings(fit <- brms::add_criterion(fit,
                                                        add_fit_criteria, 
                                                        cores = cores))
          }
          if (!is.null(resp)) {
            for (aci in fit$model_info$ys) {
              suppressWarnings(fit <- brms::add_criterion(
                fit,
                add_fit_criteria,
                resp = aci,
                cores = cores
              ))
              aci_names <- paste0(names(fit$criteria), aci)
              names(fit$criteria) <- aci_names
            }
            aci_names <- c()
            for (aci in fit$model_info$ys) {
              aci_names <- c(aci_names, paste0(add_fit_criteria, aci))
            }
            names(fit$criteria) <- aci_names
          }
        }
      }
      
      if (!is.na(fit$model_info$univariate_by)) {
        for (aci in fit$model_info$ys) {
          suppressWarnings(fit <- brms::add_criterion(
            fit,
            add_fit_criteria,
            resp = aci,
            cores = cores
          ))
          aci_names <- paste0(names(fit$criteria), aci)
          names(fit$criteria) <- aci_names
        }
        aci_names <- c()
        for (aci in fit$model_info$ys) {
          aci_names <- c(aci_names, paste0(add_fit_criteria, aci))
        }
        names(fit$criteria) <- aci_names
      }
    } # if (!is.null(add_fit_criteria))
    
    
    if (!is.null(add_bayes_R)) {
      what_ <- paste(add_bayes_R, collapse = ", ")
      if(verbose) message(" Adding", " ", what_, " ", "...")
      if(verbose) cat("\n")
      if (is.na(fit$model_info$univariate_by)) {
        if (!fit$model_info$multivariate) {
          aci_names <- paste0(add_bayes_R, '')
          suppressWarnings(fit$criteria[[aci_names]] <-
                             brms::bayes_R2(fit, cores = cores))
          fit$criteria[[aci_names]] <-
            fit$criteria[[aci_names]] %>%
            data.frame() %>% dplyr::mutate(Parameter = rownames(.)) %>%
            dplyr::relocate(Parameter)
          rownames(fit$criteria[[aci_names]]) <- NULL
        }
        if (fit$model_info$multivariate) {
          if (is.null(resp)) {
            aci_names <- paste0(add_bayes_R, '')
            suppressWarnings(fit$criteria[[aci_names]] <-
                               brms::bayes_R2(fit,
                                              cores = cores))
            fit$criteria[[aci_names]] <-
              fit$criteria[[aci_names]] %>%
              data.frame() %>% dplyr::mutate(Parameter = rownames(.)) %>%
              dplyr::relocate(Parameter)
            rownames(fit$criteria[[aci_names]]) <- NULL
          }
          if (!is.null(resp)) {
            for (aci in fit$model_info$ys) {
              aci_names <- paste0(add_bayes_R, aci)
              suppressWarnings(fit$criteria[[aci_names]] <-
                                 brms::bayes_R2(fit,
                                                resp = aci,
                                                cores = cores))
              fit$criteria[[aci_names]] <-
                fit$criteria[[aci_names]] %>%
                data.frame() %>% dplyr::mutate(Parameter = rownames(.)) %>%
                dplyr::relocate(Parameter)
              rownames(fit$criteria[[aci_names]]) <- NULL
            }
          }
        }
      }
      
      
      
      
      if (!is.na(fit$model_info$univariate_by)) {
        for (aci in fit$model_info$ys) {
          aci_names <- paste0(add_bayes_R, aci)
          suppressWarnings(fit$criteria[[aci_names]] <-
                             brms::bayes_R2(fit,
                                            resp = aci,
                                            cores = cores))
          fit$criteria[[aci_names]] <-
            fit$criteria[[aci_names]] %>%
            data.frame() %>% dplyr::mutate(Parameter = rownames(.)) %>%
            dplyr::relocate(Parameter)
          rownames(fit$criteria[[aci_names]]) <- NULL
        }
      }
    } # if (!is.null(add_bayes_R)) {
    
    
    
    ################
    add_summary_waic <- function(x, digits = 1) {
      summary_waic <- x
      summary_waic$pointwise <- NULL
      summary_waic <- summary_waic$estimates
      summary_waic <- summary_waic %>% data.frame()
      summary_waic <- summary_waic %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                    ~ round(., digits = digits)))
      summary_waic$Parameter <- row.names(summary_waic)
      row.names(summary_waic) <- NULL
      summary_waic <-
        summary_waic %>% dplyr::relocate(Parameter, Estimate, SE)
      summary_waic
    }
    
    add_summary_bayes_R2 <- function(x, digits = 2) {
      summary_bayes_R <- x
      summary_bayes_R <- summary_bayes_R %>% data.frame()
      summary_bayes_R <- summary_bayes_R %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                    ~ round(., digits = digits)))
      row.names(summary_bayes_R) <- NULL
      summary_bayes_R$Parameter <- 'bayes_R2'
      summary_bayes_R$SE <- summary_bayes_R$Est.Error
      summary_bayes_R <-
        summary_bayes_R %>% dplyr::select(-c(Est.Error))
      summary_bayes_R <- summary_bayes_R %>%
        dplyr::relocate(Parameter, Estimate, SE)
      summary_bayes_R
    }
    
    
    
    add_summary_loo <- function(x, digits = 1) {
      summary_loo <- x
      summary_loo$pointwise <- NULL
      summary_loo <- summary_loo$estimates
      summary_loo <- summary_loo %>% data.frame()
      summary_loo <- summary_loo %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                    ~ round(., digits = digits)))
      summary_loo$Parameter <- row.names(summary_loo)
      row.names(summary_loo) <- NULL
      summary_loo <-
        summary_loo %>% dplyr::relocate(Parameter, Estimate, SE)
      summary_loo
    }
    
    add_diagnostic_loo <- function(x, digits = 1) {
      summary_loo_diagnostic <- loo::pareto_k_table(x) %>% data.frame()
      row.names(summary_loo_diagnostic) <- NULL
      summary_loo_diagnostic$Range <- attr(loo::pareto_k_table(x),
                                           "dimnames")[[1]]
      summary_loo_diagnostic$Inference <-
        c('Good', "Ok", "Bad", "Very bad")
      summary_loo_diagnostic$Percent <-
        round(summary_loo_diagnostic$Proportion * 100, digits)
      summary_loo_diagnostic$Min.n_eff  <-
        summary_loo_diagnostic$Min..n_eff
      summary_loo_diagnostic$Min.n_eff <-
        round(summary_loo_diagnostic$Min.n_eff)
      summary_loo_diagnostic <- summary_loo_diagnostic %>%
        dplyr::select(-c(Proportion, Min..n_eff))
      summary_loo_diagnostic <- summary_loo_diagnostic %>%
        dplyr::relocate(Range, Inference, Count, Percent, Min.n_eff)
      summary_loo_diagnostic
    }
    
    
    if ('waic' %in% add_fit_criteria) {
      err. <- FALSE
      tryCatch(
        expr = {
          if (!is.na(fit$model_info$univariate_by)) {
            list_c_ <- list()
            for (aci in fit$model_info$ys) {
              getit_ <- paste0('waic', aci)
              list_c_[[aci]] <-
                add_summary_waic(fit$criteria[[getit_]], 
                                 digits = digits) %>%
                dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
            }
            summary_waic <-
              list_c_ %>%  do.call(rbind, .) %>% data.frame()
          } else if (fit$model_info$multivariate & !is.null(resp)) {
            list_c_ <- list()
            for (aci in fit$model_info$ys) {
              getit_ <- paste0('waic', aci)
              list_c_[[aci]] <-
                add_summary_waic(fit$criteria[[getit_]], 
                                 digits = digits) %>%
                dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
            }
            summary_waic <-
              list_c_ %>%  do.call(rbind, .) %>% data.frame()
          } else if (fit$model_info$multivariate & is.null(resp)) {
            getit_ <- paste0('waic', '')
            summary_waic <-
              add_summary_waic(fit$criteria[[getit_]], digits = digits)
          } else if (is.na(fit$model_info$univariate_by) &
                     !fit$model_info$multivariate) {
            getit_ <- paste0('waic', '')
            summary_waic <-
              add_summary_waic(fit$criteria[[getit_]], digits = digits)
          }
          summary_waic$df <- df
          summary_waic$xfun <- xfun_print
          summary_waic$yfun <- yfun_print
          summary_waic <-
            summary_waic %>% dplyr::relocate(df, xfun, yfun)
          rownames(summary_waic) <- NULL
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (err.) {
        summary_waic <- NULL
      } else {
        summary_waic <- summary_waic
      }
      fit$summary_waic <- summary_waic
    }
    
    
    
    
    
    if ('bayes_R2' %in% add_bayes_R) {
      err. <- FALSE
      tryCatch(
        expr = {
          if (!is.na(fit$model_info$univariate_by)) {
            list_c_ <- list()
            for (aci in fit$model_info$ys) {
              getit_ <- paste0(add_bayes_R, aci)
              list_c_[[aci]] <-
                add_summary_bayes_R2(fit$criteria[[getit_]], 
                                     digits = digits) %>%
                dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
            }
            summary_bayes_R2 <-
              list_c_ %>%  do.call(rbind, .) %>% data.frame()
          } else if (fit$model_info$multivariate & !is.null(resp)) {
            list_c_ <- list()
            for (aci in fit$model_info$ys) {
              getit_ <- paste0(add_bayes_R, aci)
              list_c_[[aci]] <-
                add_summary_bayes_R2(fit$criteria[[getit_]], 
                                     digits = digits) %>%
                dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
            }
            summary_bayes_R2 <-
              list_c_ %>%  do.call(rbind, .) %>% data.frame()
          } else if (fit$model_info$multivariate & is.null(resp)) {
            getit_ <- paste0('bayes_R2', '')
            summary_bayes_R2 <-
              add_summary_bayes_R2(fit$criteria[[getit_]], 
                                   digits = digits)
          } else if (is.na(fit$model_info$univariate_by) &
                     !fit$model_info$multivariate) {
            getit_ <- paste0('bayes_R2', '')
            summary_bayes_R2 <-
              add_summary_bayes_R2(fit$criteria[[getit_]], 
                                   digits = digits)
          }
          summary_bayes_R2$df <- df
          summary_bayes_R2$xfun <- xfun_print
          summary_bayes_R2$yfun <- yfun_print
          summary_bayes_R2 <-
            summary_bayes_R2 %>% dplyr::relocate(df, xfun, yfun)
          rownames(summary_bayes_R2) <- NULL
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (err.) {
        summary_bayes_R2 <- NULL
      } else {
        summary_bayes_R2 <- summary_bayes_R2
      }
      fit$summary_bayes_R2 <-
        summary_bayes_R2 %>% dplyr::select(-Parameter)
    }
    
    
    
    if ('loo' %in% add_fit_criteria) {
      if ('loo' %in% add_fit_criteria) {
        err. <- FALSE
        tryCatch(
          expr = {
            if (!is.na(fit$model_info$univariate_by)) {
              list_c_ <- list()
              for (aci in fit$model_info$ys) {
                getit_ <- paste0('loo', aci)
                list_c_[[aci]] <-
                  add_summary_loo(fit$criteria[[getit_]], 
                                  digits = digits) %>%
                  dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
              }
              summary_loo <-
                list_c_ %>%  do.call(rbind, .) %>% data.frame()
            } else if (fit$model_info$multivariate &
                       !is.null(resp)) {
              list_c_ <- list()
              for (aci in fit$model_info$ys) {
                getit_ <- paste0('loo', aci)
                list_c_[[aci]] <-
                  add_summary_loo(fit$criteria[[getit_]], 
                                  digits = digits) %>%
                  dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
              }
              summary_loo <-
                list_c_ %>%  do.call(rbind, .) %>% data.frame()
            } else if (fit$model_info$multivariate &
                       is.null(resp)) {
              getit_ <- paste0('loo', '')
              
              summary_loo <-
                add_summary_loo(fit$criteria[[getit_]], digits = digits)
            } else if (is.na(fit$model_info$univariate_by) &
                       !fit$model_info$multivariate) {
              getit_ <- paste0('loo', '')
              summary_loo <-
                add_summary_loo(fit$criteria[[getit_]], digits = digits)
            }
            summary_loo$df <- df
            summary_loo$xfun <- xfun_print
            summary_loo$yfun <- yfun_print
            summary_loo <-
              summary_loo %>% dplyr::relocate(df, xfun, yfun)
            rownames(summary_loo) <- NULL
          },
          error = function(e) {
            err. <<- TRUE
          }
        )
        if (err.) {
          summary_loo <- NULL
        } else {
          summary_loo <- summary_loo
        }
        fit$summary_loo <- summary_loo
      }
      
      if ('loo' %in% add_fit_criteria) {
        err. <- FALSE
        tryCatch(
          expr = {
            if (!is.na(fit$model_info$univariate_by)) {
              list_c_ <- list()
              for (aci in fit$model_info$ys) {
                getit_ <- paste0('loo', aci)
                list_c_[[aci]] <-
                  add_diagnostic_loo(fit$criteria[[getit_]], 
                                     digits = digits) %>%
                  dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
              }
              diagnostic_loo <-
                list_c_ %>%  do.call(rbind, .) %>% data.frame()
            } else if (fit$model_info$multivariate &
                       !is.null(resp)) {
              list_c_ <- list()
              for (aci in fit$model_info$ys) {
                getit_ <- paste0('loo', aci)
                list_c_[[aci]] <-
                  add_diagnostic_loo(fit$criteria[[getit_]], 
                                     digits = digits) %>%
                  dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
              }
              diagnostic_loo <-
                list_c_ %>%  do.call(rbind, .) %>% data.frame()
            } else if (fit$model_info$multivariate &
                       is.null(resp)) {
              getit_ <- paste0('loo', '')
              diagnostic_loo <-
                add_diagnostic_loo(fit$criteria[[getit_]], 
                                   digits = digits)
            } else if (is.na(fit$model_info$univariate_by) &
                       !fit$model_info$multivariate) {
              diagnostic_loo <-
                add_diagnostic_loo(fit$criteria[[getit_]], 
                                   digits = digits)
            }
            diagnostic_loo$df <- df
            diagnostic_loo$xfun <- xfun_print
            diagnostic_loo$yfun <- yfun_print
            diagnostic_loo <-
              diagnostic_loo %>% dplyr::relocate(df, xfun, yfun)
            rownames(diagnostic_loo) <- NULL
          },
          error = function(e) {
            err. <<- TRUE
          }
        )
        if (err.) {
          diagnostic_loo <- NULL
        } else {
          diagnostic_loo <- diagnostic_loo
        }
        fit$diagnostic_loo <- diagnostic_loo
      }
    } # if('loo' %in% add_fit_criteria) {
    
    return(fit)
  } # add_citeria_fun
  
  
  
  optimize_fun <- function(.x, model) {
    message("\nOptimizing model no. ",
            .x,
            " (total ",
            nrow(optimize_df_x_y),
            " models)")
    exe_row <- optimize_df_x_y[.x, ]
    df <- levels(droplevels(exe_row$df))
    xfun <- levels(droplevels(exe_row$xfun))
    yfun <- levels(droplevels(exe_row$yfun))
    if (df == 'NULL')
      df <-
      paste0("list(", paste(model$model_info$dfs, collapse = ","), ")")
    else
      df <- df
    if (xfun == 'NULL')
      xfun <- NULL
    else
      xfun <- xfun
    if (yfun == 'NULL')
      yfun <- NULL
    else
      yfun <- yfun
    
    if (is.null(xfun))
      xfun_print <- deparse(xfun)
    else
      xfun_print <- xfun
    if (is.null(yfun))
      yfun_print <- deparse(yfun)
    else
      yfun_print <- yfun
    
    if(verbose) {
      cat("\n")
      cat(paste0("df = ", df, "; xfun = ", xfun_print, "; yfun = ", yfun_print),
          "\n")
    }
    
    
    optimization_info <-
      paste0("df = ", df, "; xfun = ", xfun_print, "; yfun = ", yfun_print)
    
    
    args_o$model <- model
    args_o$df    <- eval(parse(text = df))
    args_o$xfun  <- xfun
    args_o$yfun  <- yfun
    args_o$data  <- newdata %>% data.frame()
    
    
    
    args_o$model  <- NULL
    
    args_o_new <- args_o
    calling    <- model$model_info$call.full.bgmfit
    
    args_o_org <- calling
    args_o_org[[1]] <- NULL
    
    args_o_new$data <- NULL
    args_o_org$data <- NULL
    
    if(is.na(model$model_info$univariate_by) &
       !model$model_info$multivariate) {
      if(length(args_o_new$df) == 1)   args_o_new$df   <- args_o_new$df[[1]]
      if(length(args_o_new$xfun) == 1) args_o_new$xfun <- args_o_new$xfun[[1]]
      if(length(args_o_new$yfun) == 1) args_o_new$yfun <- args_o_new$yfun[[1]]
    }
    
    
    all_same_args_c <- all_same_args <- c()
    # args_o_org_updated <- list()
    for (args_oi in names(args_o_new)) {
      all_same_args_c <- c(all_same_args_c, identical(args_o_org[[args_oi]],
                                                      args_o_new[[args_oi]]) 
      )
    }
    
    
    all_same_args_c_diffs <- args_o_new[!all_same_args_c]
    
    if(length(all_same_args_c_diffs) > 0) {
      all_same_args <- FALSE 
    } else {
      all_same_args <- TRUE
    }
    
    mandatory_opts <- c('df', 'xfun', 'yfun')
    
    
    if(all_same_args) {
      if(verbose) {
        cat("\n")
        cat("The arguemnets supplied for optimization are identical to the", 
            "\n ",
            "original model fit. Therefore, returning the original model fit")
        cat("\n")
      }
      return(model)
    } else if(!all_same_args) {
      user_call   <- calling
      user_call   <- rlang::call_match(user_call, bsitar::bsitar)
      newargs     <- all_same_args_c_diffs
      for (newargsi in names(newargs)) {
        user_call[[newargsi]] <- NULL
      }
      user_call_data_name <- user_call$data
      assign(deparse(user_call_data_name), newdata)
      user_call <- rlang::call_modify(user_call, !!!newargs)
      fit <- eval(user_call)
    }
    
    
    
    fit$model_info$optimization_info <- optimization_info
    fit$model_info$optimize_df <- df
    fit$model_info$optimize_x <- xfun_print
    fit$model_info$optimize_y <- yfun_print
    
    # Add fit_criteria and bares_R to the fit
    # Add summary data frames for criteria and R square
    
    # setresp to anything so that even multivariate will be response wise
    # if desired, this behavious
    # if(length(fit$model_info$ys) == 1) setresp <- NULL
    # if(length(fit$model_info$ys) > 1) setresp <- 'TRUE'
    
    if (fit$model_info$multivariate) {
      if (byresp) {
        setresp <- 'TRUE'
      } else if (!byresp) {
        setresp <- NULL
      }
    } else if (!fit$model_info$multivariate) {
      setresp <- NULL
    }
    
    if (!is.null(add_fit_criteria)) {
      fit <- add_citeria_fun(
        fit,
        add_fit_criteria = add_fit_criteria,
        add_bayes_R =  NULL,
        resp = setresp,
        digits = digits,
        df = df,
        xfun_print = xfun_print,
        yfun_print = yfun_print
      )
    }
    
    if (!is.null(add_bayes_R)) {
      fit <- add_citeria_fun(
        fit,
        add_fit_criteria = NULL,
        add_bayes_R =  add_bayes_R,
        resp = setresp,
        digits = digits,
        df = df,
        xfun_print = xfun_print,
        yfun_print = yfun_print
      )
    }
    return(fit)
  }
  
  optimize_list <- lapply(1:nrow(optimize_df_x_y), function(.x)
    optimize_fun(.x, model))
  
  loo_fit             <- combine_summaries(optimize_list, 'summary_loo')
  loo_diagnostic_fit  <-
    combine_summaries(optimize_list, 'diagnostic_loo')
  waic_fit            <-
    combine_summaries(optimize_list, 'summary_waic')
  bayes_R2_fit        <-
    combine_summaries(optimize_list, 'summary_bayes_R2')
  
  # Parameter column is not created earlier for the bayes_R2_fit
  bayes_R2_fit <- bayes_R2_fit %>% 
    dplyr::mutate(Parameter = 'bayes_R2') %>% 
    dplyr::relocate(df, xfun, yfun, Parameter)
  
  
  attributes(optimize_list) <- NULL
  
  optimize_summary <- data.frame()
  
  if(exists('loo_fit')) optimize_summary <- optimize_summary %>% 
    dplyr::bind_rows(., loo_fit)
  
  if(exists('loo_diagnostic_fit')) optimize_summary <- optimize_summary %>% 
    dplyr::bind_rows(., loo_diagnostic_fit)
  
  if(exists('waic_fit')) optimize_summary <- optimize_summary %>% 
    dplyr::bind_rows(., waic_fit)
  
  if(exists('bayes_R2_fit')) optimize_summary <- optimize_summary %>% 
    dplyr::bind_rows(., bayes_R2_fit)
  
  out <- list(models = optimize_list, optimize_summary = optimize_summary)
  
  return(out)
}




