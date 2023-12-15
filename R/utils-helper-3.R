

#' An internal function to set up data for \code{bsitar} model
#'
#' @param data An input data frame.
#'
#' @param x The predictor (typically age) variables.
#'
#' @param y The outcome variables. Must be a single name except when fitting a
#'   multivariate model.
#'
#' @param id The group identifier.
#'
#' @param uvarby An optional (default \code{NA}) to specify the indicator
#'   variable for fitting univariate-by-subgroup model. See \code{univariate_by}
#'   argument in the [bsitar::bsitar()] function. If not \code{NA}, then it should
#'   be a valid factor variable present in the \code{data}.
#'
#' @param mvar A logical (default \code{FALSE}) to specify the the multivariate
#'   model. See \code{multivariate} argument in the [bsitar::bsitar()] function.
#'
#' @param xfuns Optional name(s) of the transformation function(s) applied to
#'   the predictor variable (typically age). Default \code{NULL}.
#'
#' @param yfuns Optional name(s) of the transformation function(s) applied to
#'   the outcome variable. Default \code{NULL}.
#'
#' @param outliers An optional (default \code{NULL}) to remove velocity
#'   outliers. The argument should be a named list to pass options to the
#'   [bsitar::outliers()] function. See [bsitar::outliers()] for details.
#'   
#' @param subset A logical (default \code{TRUE}) to indicate whether to create
#'   data for each level of the \code{univariate_by} variable, or only for a
#'   subset of levels. The \code{subset = TRUE} is typically used during model
#'   fit and \code{subset = FALSE} during post processing of each sub model. The
#'   argument \code{subset} is ignored when \code{univariate_by} is \code{NA} or
#'   \code{NULL}.
#' 
#' @return A data frame with necessary information added a attributes.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#' 
prepare_data <- function(data,
                         x,
                         y,
                         id,
                         uvarby = NA,
                         mvar = FALSE,
                         xfuns = NULL,
                         yfuns = NULL,
                         outliers = NULL,
                         subset = TRUE) {
  
  data <- data %>% droplevels()
  
  if (!is.null(outliers)) {
    remove_ <- outliers$remove
    icode_ <- outliers$icode
    icode_ <- deparse(substitute(icode_))
    limit_ <- outliers$limit
    velpower_ <- outliers$velpower
    lag_ <- outliers$lag
    linearise_ <- outliers$linearise
    verbose_ <- outliers$verbose
    
    for (yi in 1:length(y)) {
      if (!y[yi] %in% colnames(data)) {
        stop(
          "When model is fit with argument outliers (i.e., outliers not NULL), ",
          "\n",
          "  then outcome variable should be part of the newdata specified.",
          "\n",
          "  please check the missing outcome varibale: ",
          y[yi]
        )
      }
      if (!x[yi] %in% colnames(data)) {
        stop(
          "When model is fit with argument outliers (i.e., outliers not NULL),",
          " \n ",
          "  then predictor variable should be part of the newdata specified.",
          "\n",
          "  please check the missing predictor varibale: ",
          x[yi]
        )
      }
      if (!id[yi] %in% colnames(data)) {
        stop(
          "When model is fit with argument outliers 
          (i.e., outliers not NULL), ",
          "\n",
          "  then group identifier variable should be 
          part of the newdata specified.",
          "\n",
          "  please check the missing group identifier varibale: ",
          id[yi]
        )
      }
      data <-
        outliers(
          x = x[yi],
          y =  y[yi],
          id = id[yi],
          data = data,
          icode = icode_,
          lag = lag_,
          velpower = velpower_,
          limit = limit_,
          linearise = linearise_,
          remove = remove_,
          verbose = verbose_
        )
      
    }
  } # if(!is.null(outliers)) {
  
  org.data <- data
  
  # Note that x tarnsformation is done within the prepare_function
  transform_y <- function(y, yfuns) {
    for (myfunsi in 1:length(y)) {
      mysi <- y[[myfunsi]]
      myfunsi <- yfuns[[myfunsi]]
      if (grepl('.Primitive', myfunsi, fixed = T) &
          grepl('log', myfunsi, fixed = T)) {
        myfunsi <- 'log'
      }
      if (grepl('.Primitive', myfunsi, fixed = T) &
          grepl('sqrt', myfunsi, fixed = T)) {
        myfunsi <- 'sqrt'
      }
      if (myfunsi == 'log')
        if(!is.null(data[[mysi]])) data[[mysi]] <- log(data[[mysi]])
      if (myfunsi == 'sqrt')
        if(!is.null(data[[mysi]])) data[[mysi]] <- sqrt(data[[mysi]])
    }
    return(data)
  }
  
  if (!(is.na(uvarby) | uvarby == "NA")) {
    if (!uvarby %in% colnames(data)) {
      stop(paste(
        "\nvariable",
        uvarby,
        "used for setting univariate submodels is missing"
      ))
    }
    if (!is.factor(data[[uvarby]])) {
      stop("subset by variable '",
           uvarby,
           "' should be a factor variable")
    }
    for (l in levels(data[[uvarby]])) {
      if(!subset) data[[l]] <- data[[y[1]]]
      if(subset) data[[l]]  <- data[[l]]
    }
    unibyimat <-
      model.matrix(~ 0 + eval(parse(text = uvarby)), data)
    subindicators <- paste0(uvarby, levels(data[[uvarby]]))
    colnames(unibyimat) <- subindicators
    y <- levels(data[[uvarby]])
    data <- as.data.frame(cbind(data, unibyimat))
    data <- transform_y(y, yfuns)
    attr(data, "ys") <- y
    attr(data, "multivariate") <- FALSE
    attr(data, "uvarby") <- uvarby
    attr(data, "subindicators") <- subindicators
    # data_out <- data
  } else if (mvar) {
    data <- org.data
    data <- transform_y(y, yfuns)
    attr(data, "ys") <- y
    attr(data, "multivariate") <- TRUE
    attr(data, "uvarby") <- NULL
    attr(data, "subindicators") <- NULL
  } else {
    # data_out <- org.data
    data <- org.data
    data <- transform_y(y, yfuns)
    attr(data, "ys") <- y
    attr(data, "multivariate") <- FALSE
    attr(data, "uvarby") <- NULL
    attr(data, "subindicators") <- NULL
  }
  return(data)
}











#' An internal function to prepare formula for Bayesian SITAR growth curve model
#' 
#' The \code{prepare_formula}) prepares \code{brms::brmsformual} which is  
#' passed on to the [bsitar::bsitar()] function. For univariate-by-
#' subgroup model (specified by using the \code{univariate_by}) and 
#' multivariate model (specified by using the \code{multivariate}),
#' the \code{x}, \code{y}, \code{id}, \code{knots}, \code{nknots}, are 
#' automatically set to match the sub-model(s). See \code{brms::brmsformual} 
#' for details. 
#'
#' @param x vector of predictor (typically age in years).
#' 
#' @param y vector of outcome (i.e., repeated growth measurements). 
#' 
#' @param id a factor variable identifying the groups (typically individuals).
#' 
#' @param knots vector of values for knots.
#' 
#' @param nknots an integer specifying the number of knots.
#' 
#' @param data data frame containing variables \code{x}, \code{y} and \code{id}.
#' 
#' @param internal_formula_args Other internal arguments passed from the 
#' [bsitar::bsitar()] to the \code{prepare_formula}).
#'
#' @return An object of class \code{brmsformula}, which is a \code{list} 
#'   containing formulas.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'   
#' @keywords internal
#' @noRd
#' 
prepare_formula <- function(x,
                            y,
                            id,
                            knots,
                            nknots,
                            data,
                            internal_formula_args) {
  
  
  # Initiate non formalArgs()
  randomsi <- NULL;
  sigma_formula_gr_strsi <- NULL;
  fixedsi <- NULL;
  select_model <- NULL;
  match_sitar_d_form <- NULL;
  spfncname <- NULL;
  terms_rhssi <- NULL;
  subindicatorsi <- NULL;
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
  sigma_formulasi <- NULL;
  sigma_formula_grsi <- NULL;
  set_higher_levels <- NULL;
  a_formula_gr_strsi_present <- NULL;
  a_formula_gr_strsi <- NULL;
  b_formula_gr_strsi_present <- NULL;
  b_formula_gr_strsi <- NULL;
  c_formula_gr_strsi_present <- NULL;
  c_formula_gr_strsi <- NULL;
  d_formula_gr_strsi_present <- NULL;
  d_formula_gr_strsi <- NULL;
  e_formula_gr_strsi_present <- NULL;
  e_formula_gr_strsi <- NULL;
  f_formula_gr_strsi_present <- NULL;
  f_formula_gr_strsi <- NULL;
  g_formula_gr_strsi_present <- NULL;
  g_formula_gr_strsi <- NULL;
  h_formula_gr_strsi_present <- NULL;
  h_formula_gr_strsi <- NULL;
  i_formula_gr_strsi_present <- NULL;
  i_formula_gr_strsi <- NULL;
  s_formula_gr_strsi_present <- NULL;
  s_formula_gr_strsi <- NULL;
  sigma_set_higher_levels <- NULL;
  sigma_formula_gr_strsi_present <- NULL;
  autocor_formi <- NULL;
  unusedsi <- NULL;
  familysi <- NULL;
  mat_s <- NULL;
  ancov_gr_str <- NULL;
  bncov_gr_str <- NULL;
  cncov_gr_str <- NULL;
  dncov_gr_str <- NULL;
  encov_gr_str <- NULL;
  fncov_gr_str <- NULL;
  gncov_gr_str <- NULL;
  hncov_gr_str <- NULL;
  incov_gr_str <- NULL;
  sncov_gr_str <- NULL;
  acovcoefnames_gr_str_id <- NULL;
  bcovcoefnames_gr_str_id <- NULL;
  ccovcoefnames_gr_str_id <- NULL;
  dcovcoefnames_gr_str_id <- NULL;
  ecovcoefnames_gr_str_id <- NULL;
  fcovcoefnames_gr_str_id <- NULL;
  gcovcoefnames_gr_str_id <- NULL;
  hcovcoefnames_gr_str_id <- NULL;
  icovcoefnames_gr_str_id <- NULL;
  scovcoefnames_gr_str_id <- NULL;
  acovcoefnames_gr_str_form <- NULL;
  bcovcoefnames_gr_str_form <- NULL;
  ccovcoefnames_gr_str_form <- NULL;
  dcovcoefnames_gr_str_form <- NULL;
  ecovcoefnames_gr_str_form <- NULL;
  fcovcoefnames_gr_str_form <- NULL;
  gcovcoefnames_gr_str_form <- NULL;
  hcovcoefnames_gr_str_form <- NULL;
  icovcoefnames_gr_str_form <- NULL;
  scovcoefnames_gr_str_form <- NULL;
  sigmancov_gr_str <- NULL;
  sigmacovcoefnames_gr_str <- NULL;
  sigmacovcoefnames_gr_str_id <- NULL;
  sigmacovcoefnames_gr_str_form <- NULL;
  brms_arguments <- NULL;
  
  d_adjustedsi <- NULL;
  
  
  
  if (!is.null(internal_formula_args)) {
    eout <- list2env(internal_formula_args)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  fixedelements <- strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
  for (cfi in letters[1:26]) {
    if(!cfi %in% fixedelements) {
      assign(paste0(cfi, '_', 'formulasi'), NULL)
      assign(paste0(cfi, '', 'form'), NULL)
    }
  }
  
  randomelements <- strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
  for (cfi in letters[1:26]) {
    if(!cfi %in% randomelements) {
      assign(paste0(cfi, '_', 'formula_grsi'), NULL)
    }
  }
  
  
  # Need separate true false for random effects and sigma
  
  # First assign all FALSE
  # This to ignore _str formulae even if present when parmare not random
  for (set_randomsi_higher_levsli in c(letters[1:26], 'sigma')) {
    set_nlpar_what <- set_randomsi_higher_levsli
    assign(paste0(set_nlpar_what, '_formula_gr_strsi_present'), FALSE)
  }
  
  set_randomsi_higher_levsl_nosigma <- strsplit(gsub("\\+", " ", 
                                                     randomsi), " ")[[1]]
  set_randomsi_higher_levsl <- strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
  set_randomsi_higher_levsl <- c(set_randomsi_higher_levsl, 'sigma')
  for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
    set_nlpar_what <- set_randomsi_higher_levsli
    in_gr_strsi <- paste0(set_nlpar_what, '_formula_gr_strsi')
    if(!is.null(ept(in_gr_strsi)[[1]]) & !is.null(ept(in_gr_strsi))) {
      assign(paste0(set_nlpar_what, '_formula_gr_strsi_present'), TRUE)
    } else {
      assign(paste0(set_nlpar_what, '_formula_gr_strsi_present'), FALSE)
    }
  }
  
  
  # Now true false for sigma 
  if(!is.null(sigma_formula_gr_strsi[[1]]) & !is.null(sigma_formula_gr_strsi)) {
    assign(paste0('sigma', '_formula_gr_strsi_present'), TRUE)
  } else {
    assign(paste0('sigma', '_formula_gr_strsi_present'), FALSE)
  }
  
  
  gr_str_id_all   <- sigma_str_id_all   <- c()
  gr_str_form_all <- sigma_str_form_all <- c()
  gr_str_coef_all <- sigma_str_coef_all <- c()
  gr_str_ncov_all  <- sigma_str_ncov_all   <- c()
  gr_str_corr_all  <- sigma_str_corr_all   <- c()
  gr_str_corr_tf_all  <- sigma_str_corr_tf_all   <- c()
  gr_str_id_all_list <- gr_str_corr_all_list <- gr_str_corr_tf_all_list <-list()
  sigma_str_id_all_list <- sigma_str_corr_all_list <- list()
  sigma_str_corr_tf_all_list <- list()
  
  # First assign all FALSE
  # This to ignore _str formulae even if present when parmare not random
  for (set_randomsi_higher_levsli in c(letters[1:26], 'sigma')) {
    set_nlpar_what <- set_randomsi_higher_levsli
    assign(paste0(set_nlpar_what, 'covcoefnames_gr_str'), NULL)
    assign(paste0(set_nlpar_what, 'covcoefnames_gr_str_id'), NULL)
    assign(paste0(set_nlpar_what, 'covcoefnames_gr_str_form'), NULL)
    assign(paste0(set_nlpar_what, 'ncov_gr_str'), NULL)
  }
  
  for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
    set_nlpar_what <- set_randomsi_higher_levsli
    if(ept(paste0(set_nlpar_what, '_formula_gr_strsi_present'))) {
      in_gr_strsi <- paste0(set_nlpar_what, '_formula_gr_strsi')
      if(!grepl("|" , ept(in_gr_strsi), fixed = T)) {
        stop("For '_str' approach of setting the random effects,",
             "\n ", 
             " only the vertical bar '|' approach is aloowed.",
             "\n ", 
             " Please check and correct the '", in_gr_strsi, "' argument ",
             "\n ", 
             " which currently specified as ", ept(in_gr_strsi)
        )
      }
      get_gr_str_coef_id_it   <- get_gr_str_coef_id(ept(in_gr_strsi),
                                                    data = data)
      assign(paste0(set_nlpar_what, 'covcoefnames_gr_str'),
             get_gr_str_coef_id_it[['tsx_c_coef']])
      assign(paste0(set_nlpar_what, 'covcoefnames_gr_str_id'),
             get_gr_str_coef_id_it[['tsx_c_id']])
      assign(paste0(set_nlpar_what, 'covcoefnames_gr_str_form'),
             get_gr_str_coef_id_it[['set_form_gr_it']])
      assign(paste0(set_nlpar_what, 'ncov_gr_str'),
             get_gr_str_coef_id_it[['set_ncov_it']])
      if(set_randomsi_higher_levsli != 'sigma') {
        gr_str_id_all   <- c(gr_str_id_all,    
                             get_gr_str_coef_id_it[['tsx_c_id']] %>% 
                               unlist())
        gr_str_form_all <- c(gr_str_form_all,  
                             get_gr_str_coef_id_it[['set_form_gr_it']] %>% 
                               unlist())
        gr_str_coef_all <- c(gr_str_coef_all,  
                             get_gr_str_coef_id_it[['tsx_c_coef']] %>% 
                               unlist())
        gr_str_ncov_all <- c(gr_str_ncov_all,  
                             get_gr_str_coef_id_it[['set_ncov_it']] %>% 
                               unlist())
        gr_str_corr_all <- c(gr_str_corr_all,  
                             get_gr_str_coef_id_it[['set_corr_it']] %>% 
                               unlist() )
        gr_str_corr_tf_all  <- 
          c(gr_str_corr_tf_all, 
            get_gr_str_coef_id_it[['set_corr_true_false']] %>% 
              unlist())
        gr_str_id_all_list[[set_randomsi_higher_levsli]] <- 
          get_gr_str_coef_id_it[['tsx_c_id']]
        gr_str_corr_all_list[[set_randomsi_higher_levsli]] <- 
          get_gr_str_coef_id_it[['set_corr_it']]
        gr_str_corr_tf_all_list[[set_randomsi_higher_levsli]] <- 
          get_gr_str_coef_id_it[['set_corr_true_false']]
      }
      if(set_randomsi_higher_levsli == 'sigma') {
        sigma_str_id_all   <- c(sigma_str_id_all, 
                                get_gr_str_coef_id_it[['tsx_c_id']] %>% 
                                  unlist())
        sigma_str_form_all <- c(sigma_str_form_all, 
                                get_gr_str_coef_id_it[['set_form_gr_it']] %>% 
                                  unlist())
        sigma_str_coef_all <- c(sigma_str_coef_all, 
                                get_gr_str_coef_id_it[['tsx_c_coef']] %>% 
                                  unlist())
        sigma_str_ncov_all <- c(sigma_str_ncov_all, 
                                get_gr_str_coef_id_it[['set_ncov_it']] %>% 
                                  unlist())
        sigma_str_corr_all <- c(sigma_str_corr_all, 
                                get_gr_str_coef_id_it[['set_corr_it']] %>% 
                                  unlist() )
        sigma_str_corr_tf_all  <- 
          c(sigma_str_corr_tf_all, 
            get_gr_str_coef_id_it[['set_corr_true_false']] %>% 
              unlist())
        sigma_str_id_all_list[[set_randomsi_higher_levsli]] <- 
          get_gr_str_coef_id_it[['tsx_c_id']]
        sigma_str_corr_all_list[[set_randomsi_higher_levsli]] <- 
          get_gr_str_coef_id_it[['set_corr_it']]
        sigma_str_corr_tf_all_list[[set_randomsi_higher_levsli]] <- 
          get_gr_str_coef_id_it[['set_corr_true_false']]
      }
    } # if(ept(paste0(set_nlpar_what, '_formula_gr_strsi_present'))) {
    if(!ept(paste0(set_nlpar_what, '_formula_gr_strsi_present'))) {
      asin_all_null <- NULL
      assign(paste0(set_nlpar_what, 'covcoefnames_gr_str'),
             asin_all_null)
      assign(paste0(set_nlpar_what, 'covcoefnames_gr_str_id'),
             asin_all_null)
      assign(paste0(set_nlpar_what, 'covcoefnames_gr_str_form'),
             asin_all_null)
      assign(paste0(set_nlpar_what, 'ncov_gr_str'),
             asin_all_null)
    } # if(!ept(paste0(set_nlpar_what, '_formula_gr_strsi_present'))) {
  } # for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
  
  
  
  gr_gsctfnb <- get_str_corr_tf_function_new_better(gr_str_id_all_list, 
                                                    gr_str_corr_all_list, 
                                                    gr_str_corr_tf_all_list)
  
  
  gr_str_corr_tf <- gr_gsctfnb[['str_corr_tf']]
  gr_str_unique_id <- gr_gsctfnb[['group_id_unique']]
  
  
  
  sigma_gsctfnb <- 
    get_str_corr_tf_function_new_better(sigma_str_id_all_list, 
                                        sigma_str_corr_all_list, 
                                        sigma_str_corr_tf_all_list)
  
  
  sigma_str_corr_tf <- sigma_gsctfnb[['str_corr_tf']]
  sigma_str_unique_id <- sigma_gsctfnb[['group_id_unique']]
  
  
  
  # For first str i.e., level 2 id, group should be unique
  
  if(!is.null(gr_str_id_all[[1]])) {
    group_arg$groupvar <- gr_str_id_all[[1]] # acovcoefnames_gr_str_id[[1]]
  } else {
    group_arg$groupvar <- group_arg$groupvar
  }
  
  
  # This to overide default 
  
  if(!is.null(sigma_str_id_all[[1]])) {
    sigma_group_arg$groupvar <- sigma_str_id_all[[1]] 
  } else {
    sigma_group_arg$groupvar <- sigma_group_arg$groupvar
  }
  
  
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
    #    abccorr <- FALSE
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
  
  abcnames <-
    paste0(strsplit(gsub("\\+", " ", fixedsi), " ")[[1]], sep = ",")
  
  snames <- c()
  for (i in 1:(nknots - 1)) {
    if (i < (nknots - 1)) {
      name1 <- paste0("s", i, sep = ",")
    }
    else {
      name1 <- paste0("s", i, sep = "")
    }
    snames[i] <- name1
  }
  
  
  if(select_model == "sitar") {
    ###########
    # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
    # In fact for df > 1, it forces 'd' to be random parameter only
    if (match_sitar_d_form) {
      if (!grepl("d", fixedsi, fixed = F) &
          grepl("d", randomsi, fixed = F)) {
        abcnames <- c(abcnames, "d,")
      }
    }
    ###########
  }
  
  
  if(select_model == "sitar" | select_model == "rcs") {
    if(select_model == "sitar") {
      if(any(grepl("s", abcnames))) abcnames <- abcnames[-length(abcnames)]
      if(match_sitar_d_form) abcnames <- gsub('s', 'd', abcnames, fixed = T)
      fullabcsnames <- c(abcnames, snames)
    }
    if(select_model == "rcs") {
      fullabcsnames <- c('a', ",", snames)
    }
  }
  
  
  
  if(select_model != "sitar" & select_model != "rcs") {
    fullabcsnames <- abcnames
  }
  
  
  
  abcselements  <- paste0(fullabcsnames, collapse = "")
  abcselements  <- paste0(spfncname, "(", x, ",", abcselements, ")")
  # abcsformfit   <- paste0(y, " ~ ", abcselements) # as.formula
  
  abcselements <- gsub(",)" , ")" , abcselements, fixed = TRUE)
  
  if (terms_rhssi == "NULL") {
    abcsformfit   <- paste0(y, " ~ ", abcselements)
  }
  
  if (terms_rhssi != "NULL") {
    abcsformfit   <- paste0(y, "|", terms_rhssi, " ~ ", abcselements)
  }
  
 
  
  
  if (!(is.na(univariate_by$by) |
        univariate_by$by == "NA") &
      !is.null(subindicatorsi)) {
    if (terms_rhssi == "NULL") {
      abcsformfit <- (paste0(y,  "|", "subset(", subindicatorsi, ")",
                             " ~ ", abcselements))
    }
    if (terms_rhssi != "NULL") {
      abcsformfit <- (paste0(y,  "|", terms_rhssi, "+", "subset(", 
                             subindicatorsi, ")",
                             " ~ ", abcselements))
    }
  }
  
  
  
  
  
  if ((is.na(univariate_by$by) |
       univariate_by$by == "NA") &
      !multivariate$mvar &  abccorr) {
    coridv <- 1 # Don't use "C" coz marginalefefcts search for it as variable
  }
  
  if ((is.na(univariate_by$by) |
       univariate_by$by == "NA") &
      !multivariate$mvar & !abccorr) {
    coridv <- ""
  }
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (uvarabccorr) {
      coridv <- y
    } else {
      coridv <- ""
    }
  }
  
  if (multivariate$mvar && mvarccorr == "none")
    coridv <- ""
  if (multivariate$mvar && mvarccorr == "W")
    coridv <- y
  if (multivariate$mvar && mvarccorr == "WB")
    coridv <- "MVC"
  
  dmatnames <- NULL
  
  if (grepl("a", fixedsi, fixed = T)) {
    afixed <- a_formulasi
  } else {
    afixed <- NULL
  }
  if (grepl("b", fixedsi, fixed = T)) {
    bfixed <- b_formulasi
  } else {
    bfixed <- NULL
  }
  if (grepl("c", fixedsi, fixed = T)) {
    cfixed <- c_formulasi
  } else {
    cfixed <- NULL
  }
  if (grepl("d", fixedsi, fixed = T)) {
    dfixed <- d_formulasi
  } else {
    dfixed <- NULL
  }
  
  if (grepl("e", fixedsi, fixed = T)) {
    efixed <- e_formulasi
  } else {
    efixed <- NULL
  }
  
  
  if (grepl("f", fixedsi, fixed = T)) {
    ffixed <- f_formulasi
  } else {
    ffixed <- NULL
  }
  
  
  if (grepl("g", fixedsi, fixed = T)) {
    gfixed <- g_formulasi
  } else {
    gfixed <- NULL
  }
  
  if (grepl("h", fixedsi, fixed = T)) {
    hfixed <- h_formulasi
  } else {
    hfixed <- NULL
  }
  
  if (grepl("i", fixedsi, fixed = T)) {
    ifixed <- i_formulasi
  } else {
    ifixed <- NULL
  }
  
  if (grepl("s", fixedsi, fixed = T)) {
    sfixed <- s_formulasi
  } else {
    sfixed <- NULL
  }
  
  # sfixed <- s_formulasi
  
  
  if (grepl("a", randomsi, fixed = T)) {
    arandom <- a_formula_grsi
  } else {
    arandom <- NULL
  }
  if (grepl("b", randomsi, fixed = T)) {
    brandom <- b_formula_grsi
  } else {
    brandom <- NULL
  }
  if (grepl("c", randomsi, fixed = T)) {
    crandom <- c_formula_grsi
  } else {
    crandom <- NULL
  }
  if (grepl("d", randomsi, fixed = T)) {
    drandom <- d_formula_grsi
  } else {
    drandom <- NULL
  }
  
  if (grepl("e", randomsi, fixed = T)) {
    erandom <- e_formula_grsi
  } else {
    erandom <- NULL
  }
  
  
  if (grepl("f", randomsi, fixed = T)) {
    frandom <- f_formula_grsi
  } else {
    frandom <- NULL
  }
  
  if (grepl("g", randomsi, fixed = T)) {
    grandom <- g_formula_grsi
  } else {
    grandom <- NULL
  }
  
  if (grepl("h", randomsi, fixed = T)) {
    hrandom <- h_formula_grsi
  } else {
    hrandom <- NULL
  }
  
  
  if (grepl("i", randomsi, fixed = T)) {
    irandom <- i_formula_grsi
  } else {
    irandom <- NULL
  }
  
  if (grepl("s", randomsi, fixed = T)) {
    srandom <- s_formula_grsi
  } else {
    srandom <- NULL
  }
  
  
  if (sigma_formulasi != "NULL") {
    sigmafixed <- sigma_formulasi
  } else {
    sigmafixed <- NULL
  }
  
  
  
  
  sigmarandom <- sigma_formula_grsi
  
  
  
  if(is.null(sigma_formulasi[[1]]) | sigma_formulasi == 'NULL') {
    display_message <-  paste("Please specify the fixed effect structure 
                              for sigma parameter ",
                              "\n ", 
                              " (by using the sigma_formula argument) 
                              since you have specified the ",
                              "\n ", 
                              " random effect structure via 
                              the sigma_formula_gr argument")
    if(!is.null(sigma_formula_grsi)) {
      stop(display_message)
    }  
  }
  
  
  
  
  
  if(!is.null(sigma_formulasi[[1]]) & !is.null(sigma_formula_grsi)) {
    if(!identical(substr(strsplit(sigma_formulasi, "~", 
                                  fixed = T)[[1]][2], 1, 1),
                  substr(strsplit(sigma_formula_grsi, "~", 
                                  fixed = T)[[1]][2], 1, 1))) {
      stop("Formulae for sigma_formula and sigma_formula_gr should be same",
           "\n ", 
           " in terms of Intercept i.e, both should be either ~ 0 or ~ 1")
    }
  }
  
  
  
  if (!is.null(dpar_formulasi) & sigma_formulasi != 'NULL') {
    stop("Either use dpar_formula or sigma_formula")
  }
  
  
  
  
  if (sigma_formulasi == 'NULL' ) {
    if (!is.null(sigma_formula_grsi) & !is.null(sigma_formula_gr_strsi)) {
      stop("You have set 'sigma_formula_gr' but 'sigma_formula' is NULL ")
    }
  }
  
  
  
  
  
  arandom_wb <- NULL
  arandom_wb_ <- FALSE
  brandom_wb <- crandom_wb <- drandom_wb <- erandom_wb <- arandom_wb
  brandom_wb_ <- crandom_wb_ <- drandom_wb_ <- erandom_wb_ <- arandom_wb_
  
  frandom_wb  <- arandom_wb
  frandom_wb_ <- arandom_wb_
  
  grandom_wb  <- arandom_wb
  grandom_wb_ <- arandom_wb_
  
  hrandom_wb  <- arandom_wb
  hrandom_wb_ <- arandom_wb_
  
  irandom_wb  <- arandom_wb
  irandom_wb_ <- arandom_wb_
  
  srandom_wb  <- arandom_wb
  srandom_wb_ <- arandom_wb_
  
  
  sigma_random_wb <- arandom_wb
  sigma_random_wb_ <- arandom_wb_
  
  get_x_random <- function(x) {
    x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
    x <- gsub("~1+1", "~1", x, fixed = T)
    if(!grepl("~", x)) x <- paste0("~", x)
    x <- sub("\\|.*", "", x)
    x
  }
  
  if(!is.null(arandom)) {
    if(grepl("|", arandom, fixed = TRUE)) {
      arandom_wb <- gsub("~", "", arandom) # with bar
      arandom_wb_ <- TRUE
      arandom <- get_x_random(arandom)
    }
  }
  if(!is.null(brandom)) {
    if(grepl("|", brandom, fixed = TRUE)) {
      brandom_wb <- gsub("~", "", brandom)
      brandom_wb_ <- TRUE
      brandom <- get_x_random(brandom)
    }
  }
  if(!is.null(crandom)) {
    if(grepl("|", crandom, fixed = TRUE)) {
      crandom_wb <- gsub("~", "", crandom)
      crandom_wb_ <- TRUE
      crandom <- get_x_random(crandom)
    }
  }
  
  if(!is.null(drandom)) {
    if(grepl("|", drandom, fixed = TRUE)) {
      drandom_wb <- gsub("~", "", drandom)
      drandom_wb_ <- TRUE
      drandom <- get_x_random(drandom)
    }
  }
  
  if(!is.null(erandom)) {
    if(grepl("|", erandom, fixed = TRUE)) {
      erandom_wb <- gsub("~", "", erandom)
      erandom_wb_ <- TRUE
      erandom <- get_x_random(erandom)
    }
  }
  
  
  if(!is.null(frandom)) {
    if(grepl("|", frandom, fixed = TRUE)) {
      frandom_wb <- gsub("~", "", frandom)
      frandom_wb_ <- TRUE
      frandom <- get_x_random(frandom)
    }
  }
  
  
  if(!is.null(grandom)) {
    if(grepl("|", grandom, fixed = TRUE)) {
      grandom_wb <- gsub("~", "", grandom)
      grandom_wb_ <- TRUE
      grandom <- get_x_random(grandom)
    }
  }
  
  
  if(!is.null(hrandom)) {
    if(grepl("|", hrandom, fixed = TRUE)) {
      hrandom_wb <- gsub("~", "", hrandom)
      hrandom_wb_ <- TRUE
      hrandom <- get_x_random(hrandom)
    }
  }
  
  
  if(!is.null(irandom)) {
    if(grepl("|", irandom, fixed = TRUE)) {
      irandom_wb <- gsub("~", "", irandom)
      irandom_wb_ <- TRUE
      irandom <- get_x_random(irandom)
    }
  }
  
  
  if(!is.null(srandom)) {
    if(grepl("|", srandom, fixed = TRUE)) {
      srandom_wb <- gsub("~", "", srandom)
      srandom_wb_ <- TRUE
      srandom <- get_x_random(srandom)
    }
  }
  
  
  if(!is.null(sigmarandom)) {
    if(grepl("|", sigmarandom, fixed = TRUE)) {
      sigma_random_wb <- gsub("~", "", sigmarandom) # with bar
      sigma_random_wb_ <- TRUE
      sigmarandom <- get_x_random(sigmarandom)
    }
  }
  
  arandom_wb <- gsub("1+1", "1", arandom_wb, fixed = T)
  brandom_wb <- gsub("1+1", "1", brandom_wb, fixed = T)
  crandom_wb <- gsub("1+1", "1", crandom_wb, fixed = T)
  drandom_wb <- gsub("1+1", "1", drandom_wb, fixed = T)
  erandom_wb <- gsub("1+1", "1", erandom_wb, fixed = T)
  frandom_wb <- gsub("1+1", "1", frandom_wb, fixed = T)
  grandom_wb <- gsub("1+1", "1", grandom_wb, fixed = T)
  hrandom_wb <- gsub("1+1", "1", hrandom_wb, fixed = T)
  irandom_wb <- gsub("1+1", "1", irandom_wb, fixed = T)
  srandom_wb <- gsub("1+1", "1", srandom_wb, fixed = T)
  
  sigma_random_wb <- gsub("1+1", "1", sigma_random_wb, fixed = T)
  
  
  if (!is.null(afixed)) {
    acovmat <- eval(parse(text = paste0(
      "model.matrix(",
      afixed, ",data = data)"
    )))
    if (ncol(acovmat) == 1)
      ancov <- NULL
    else
      ancov <- ncol(acovmat) - 1
    acovcoefnames <- colnames(acovmat)
    acovcoefnames <- gsub("\\(|)", "", acovcoefnames)
  } else if (is.null(afixed)) {
    ancov <- acovcoefnames <- NULL
  }
  
  if (!is.null(bfixed)) {
    bcovmat <- eval(parse(text = paste0(
      "model.matrix(",
      bfixed, ",data = data)"
    )))
    if (ncol(bcovmat) == 1)
      bncov <- NULL
    else
      bncov <- ncol(bcovmat) - 1
    bcovcoefnames <- colnames(bcovmat)
    bcovcoefnames <- gsub("\\(|)", "", bcovcoefnames)
  } else if (is.null(bfixed)) {
    bncov <- bcovcoefnames <- NULL
  }
  
  if (!is.null(cfixed)) {
    ccovmat <- eval(parse(text = paste0(
      "model.matrix(",
      cfixed, ",data = data)"
    )))
    if (ncol(ccovmat) == 1)
      cncov <- NULL
    else
      cncov <- ncol(ccovmat) - 1
    ccovcoefnames <- colnames(ccovmat)
    ccovcoefnames <- gsub("\\(|)", "", ccovcoefnames)
  } else if (is.null(cfixed)) {
    cncov <- ccovcoefnames <- NULL
  }
  
  if (!is.null(dfixed)) {
    dcovmat <- eval(parse(text = paste0(
      "model.matrix(",
      dfixed, ",data = data)"
    )))
    if (ncol(dcovmat) == 1)
      dncov <- NULL
    else
      dncov <- ncol(dcovmat) - 1
    dcovcoefnames <- colnames(dcovmat)
    dcovcoefnames <- gsub("\\(|)", "", dcovcoefnames)
  } else if (is.null(dfixed)) {
    dncov <- dcovcoefnames <- NULL
  }
  
  
  
  if (!is.null(efixed)) {
    ecovmat <- eval(parse(text = paste0(
      "model.matrix(",
      efixed, ",data = data)"
    )))
    if (ncol(ecovmat) == 1)
      encov <- NULL
    else
      encov <- ncol(ecovmat) - 1
    ecovcoefnames <- colnames(ecovmat)
    ecovcoefnames <- gsub("\\(|)", "", ecovcoefnames)
  } else if (is.null(efixed)) {
    encov <- ecovcoefnames <- NULL
  }
  
  
  
  if (!is.null(ffixed)) {
    fcovmat <- eval(parse(text = paste0(
      "model.matrix(",
      ffixed, ",data = data)"
    )))
    if (ncol(fcovmat) == 1)
      fncov <- NULL
    else
      fncov <- ncol(fcovmat) - 1
    fcovcoefnames <- colnames(fcovmat)
    fcovcoefnames <- gsub("\\(|)", "", fcovcoefnames)
  } else if (is.null(ffixed)) {
    fncov <- fcovcoefnames <- NULL
  }
  
  
  if (!is.null(gfixed)) {
    gcovmat <- eval(parse(text = paste0(
      "model.matrix(",
      gfixed, ",data = data)"
    )))
    if (ncol(gcovmat) == 1)
      gncov <- NULL
    else
      gncov <- ncol(gcovmat) - 1
    gcovcoefnames <- colnames(gcovmat)
    gcovcoefnames <- gsub("\\(|)", "", gcovcoefnames)
  } else if (is.null(gfixed)) {
    gncov <- gcovcoefnames <- NULL
  }
  
  
  if (!is.null(hfixed)) {
    hcovmat <- eval(parse(text = paste0(
      "model.matrix(",
      hfixed, ",data = data)"
    )))
    if (ncol(hcovmat) == 1)
      hncov <- NULL
    else
      hncov <- ncol(hcovmat) - 1
    hcovcoefnames <- colnames(hcovmat)
    hcovcoefnames <- gsub("\\(|)", "", hcovcoefnames)
  } else if (is.null(hfixed)) {
    hncov <- hcovcoefnames <- NULL
  }
  
  
  if (!is.null(ifixed)) {
    icovmat <- eval(parse(text = paste0(
      "model.matrix(",
      ifixed, ",data = data)"
    )))
    if (ncol(icovmat) == 1)
      incov <- NULL
    else
      incov <- ncol(icovmat) - 1
    icovcoefnames <- colnames(icovmat)
    icovcoefnames <- gsub("\\(|)", "", icovcoefnames)
  } else if (is.null(ifixed)) {
    incov <- icovcoefnames <- NULL
  }
  
  
  if (!is.null(sfixed)) {
    scovmat <- eval(parse(text = paste0(
      "model.matrix(",
      sfixed, ",data = data)"
    )))
    if (ncol(scovmat) == 1) {
      sncov <- NULL
    } else {
      sncov <- ncol(scovmat) - 1
    }
    scovcoefnames <- colnames(scovmat)
    scovcoefnames <- gsub("\\(|)", "", scovcoefnames)
  } else if (is.null(sfixed)) {
    sncov <- scovcoefnames <- NULL
  }
  
  
  if (!is.null(sigmafixed)) {
    sigma_covmat <- eval(parse(text = paste0(
      "model.matrix(",
      sigmafixed, ",data = data)"
    )))
    if (ncol(sigma_covmat) == 1)
      sigmancov <- NULL
    else
      sigmancov <- ncol(sigma_covmat) - 1
    sigmacovcoefnames <- colnames(sigma_covmat)
    sigmacovcoefnames <- gsub("\\(|)", "", sigmacovcoefnames)
  } else if (is.null(sigmafixed)) {
    sigmancov <- sigmacovcoefnames <- NULL
  }
  
  
  if (!is.null(arandom)) {
    acovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      arandom, ",data = data)"
    )))
    if (ncol(acovmat_gr) == 1) {
      ancov_gr <- NULL
    } else {
      ancov_gr <- ncol(acovmat_gr) - 1
    }
    acovcoefnames_gr <- colnames(acovmat_gr)
    acovcoefnames_gr <- gsub("\\(|)", "", acovcoefnames_gr)
  } else if (is.null(arandom)) {
    ancov_gr <- acovcoefnames_gr <- NULL
  }
  
  if (!is.null(brandom)) {
    bcovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      brandom, ",data = data)"
    )))
    if (ncol(bcovmat_gr) == 1) {
      bncov_gr <- NULL
    } else {
      bncov_gr <- ncol(bcovmat_gr) - 1
    }
    bcovcoefnames_gr <- colnames(bcovmat_gr)
    bcovcoefnames_gr <- gsub("\\(|)", "", bcovcoefnames_gr)
  } else if (is.null(brandom)) {
    bncov_gr <- bcovcoefnames_gr <- NULL
  }
  
  if (!is.null(crandom)) {
    ccovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      crandom, ",data = data)"
    )))
    if (ncol(ccovmat_gr) == 1) {
      cncov_gr <- NULL
    } else {
      cncov_gr <- ncol(ccovmat_gr) - 1
    }
    ccovcoefnames_gr <- colnames(ccovmat_gr)
    ccovcoefnames_gr <- gsub("\\(|)", "", ccovcoefnames_gr)
  } else if (is.null(crandom)) {
    cncov_gr <- ccovcoefnames_gr <- NULL
  }
  
  if (!is.null(drandom)) {
    dcovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      drandom, ",data = data)"
    )))
    if (ncol(dcovmat_gr) == 1) {
      dncov_gr <- NULL
    } else {
      dncov_gr <- ncol(dcovmat_gr) - 1
    }
    dcovcoefnames_gr <- colnames(dcovmat_gr)
    dcovcoefnames_gr <- gsub("\\(|)", "", dcovcoefnames_gr)
  } else if (is.null(drandom)) {
    dncov_gr <- dcovcoefnames_gr <- NULL
  }
  
  if (!is.null(erandom)) {
    ecovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      erandom, ",data = data)"
    )))
    if (ncol(ecovmat_gr) == 1) {
      encov_gr <- NULL
    } else {
      encov_gr <- ncol(ecovmat_gr) - 1
    }
    ecovcoefnames_gr <- colnames(ecovmat_gr)
    ecovcoefnames_gr <- gsub("\\(|)", "", ecovcoefnames_gr)
  } else if (is.null(erandom)) {
    encov_gr <- ecovcoefnames_gr <- NULL
  }
  
  
  if (!is.null(frandom)) {
    fcovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      frandom, ",data = data)"
    )))
    if (ncol(fcovmat_gr) == 1) {
      fncov_gr <- NULL
    } else {
      fncov_gr <- ncol(fcovmat_gr) - 1
    }
    fcovcoefnames_gr <- colnames(fcovmat_gr)
    fcovcoefnames_gr <- gsub("\\(|)", "", fcovcoefnames_gr)
  } else if (is.null(frandom)) {
    fncov_gr <- fcovcoefnames_gr <- NULL
  }
  
  if (!is.null(grandom)) {
    gcovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      grandom, ",data = data)"
    )))
    if (ncol(gcovmat_gr) == 1) {
      gncov_gr <- NULL
    } else {
      gncov_gr <- ncol(gcovmat_gr) - 1
    }
    gcovcoefnames_gr <- colnames(gcovmat_gr)
    gcovcoefnames_gr <- gsub("\\(|)", "", gcovcoefnames_gr)
  } else if (is.null(grandom)) {
    gncov_gr <- gcovcoefnames_gr <- NULL
  }
  
  
  if (!is.null(hrandom)) {
    hcovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      hrandom, ",data = data)"
    )))
    if (ncol(hcovmat_gr) == 1) {
      hncov_gr <- NULL
    } else {
      hncov_gr <- ncol(hcovmat_gr) - 1
    }
    hcovcoefnames_gr <- colnames(hcovmat_gr)
    hcovcoefnames_gr <- gsub("\\(|)", "", hcovcoefnames_gr)
  } else if (is.null(hrandom)) {
    hncov_gr <- hcovcoefnames_gr <- NULL
  }
  
  
  if (!is.null(irandom)) {
    icovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      irandom, ",data = data)"
    )))
    if (ncol(icovmat_gr) == 1) {
      incov_gr <- NULL
    } else {
      incov_gr <- ncol(icovmat_gr) - 1
    }
    icovcoefnames_gr <- colnames(icovmat_gr)
    icovcoefnames_gr <- gsub("\\(|)", "", icovcoefnames_gr)
  } else if (is.null(irandom)) {
    incov_gr <- icovcoefnames_gr <- NULL
  }
  
  
  
  if (!is.null(srandom)) {
    scovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      srandom, ",data = data)"
    )))
    if (ncol(scovmat_gr) == 1) {
      sncov_gr <- NULL
    } else {
      sncov_gr <- ncol(scovmat_gr) - 1
    }
    scovcoefnames_gr <- colnames(scovmat_gr)
    scovcoefnames_gr <- gsub("\\(|)", "", scovcoefnames_gr)
  } else if (is.null(srandom)) {
    sncov_gr <- scovcoefnames_gr <- NULL
  }
  
  
  
  
  if (!is.null(sigmarandom)) {
    sigma_covmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      sigmarandom, ",data = data)"
    )))
    if (ncol(sigma_covmat_gr) == 1) {
      sigmancov_gr <- NULL
    } else {
      sigmancov_gr <- ncol(sigma_covmat_gr) - 1
    }
    sigmacovcoefnames_gr <- colnames(sigma_covmat_gr)
    sigmacovcoefnames_gr <- gsub("\\(|)", "", sigmacovcoefnames_gr)
  } else if (is.null(sigmarandom)) {
    sigmancov_gr <- sigmacovcoefnames_gr <- NULL
  }
  
  if (!is.null(afixed)) {
    aform <- paste0("a", afixed)
  } else {
    aform <- NULL
  }
  if (!is.null(bfixed)) {
    bform <- paste0("b", bfixed)
  } else {
    bform <- NULL
  }
  if (!is.null(cfixed)) {
    cform <- paste0("c", cfixed)
  } else {
    cform <- NULL
  }
  if (!is.null(dfixed)) {
    dform <- paste0("d", dfixed)
  } else {
    dform <- NULL
  }
  
  if (!is.null(efixed)) {
    eform <- paste0("e", efixed)
  } else {
    eform <- NULL
  }
  
  
  if (!is.null(ffixed)) {
    fform <- paste0("f", ffixed)
  } else {
    fform <- NULL
  }
  
  if (!is.null(gfixed)) {
    gform <- paste0("g", gfixed)
  } else {
    gform <- NULL
  }
  
  
  if (!is.null(hfixed)) {
    hform <- paste0("h", hfixed)
  } else {
    hform <- NULL
  }
  
  
  if (!is.null(ifixed)) {
    iform <- paste0("i", ifixed)
  } else {
    iform <- NULL
  }
  
  
  if (!is.null(sfixed)) {
    sform <- paste0("s", ifixed)
  } else {
    sform <- NULL
  }
  
  
  if (!is.null(sigmafixed)) {
    sigmaform <- paste0("sigma", sigmafixed)
  } else {
    sigmaform <- NULL
  }
  
  
  sinterceptelements <- paste0(snames, collapse = "+")
  sinterceptelements <- gsub("," , "" , sinterceptelements)
  sform     <- paste0(sinterceptelements, sfixed)
  
  # sform <- NULL
  
  if (!is.null(arandom)) {
    aform_gr <- gsub("^~", "", arandom, fixed = F)
  } else {
    aform_gr <- NULL
  }
  if (!is.null(brandom)) {
    bform_gr <- gsub("^~", "", brandom, fixed = F)
  } else {
    bform_gr <- NULL
  }
  if (!is.null(crandom)) {
    cform_gr <- gsub("^~", "", crandom, fixed = F)
  } else {
    cform_gr <- NULL
  }
  if (!is.null(drandom)) {
    dform_gr <- gsub("^~", "", drandom, fixed = F)
  } else {
    dform_gr <- NULL
  }
  
  if (!is.null(erandom)) {
    eform_gr <- gsub("^~", "", erandom, fixed = F)
  } else {
    eform_gr <- NULL
  }
  
  
  if (!is.null(frandom)) {
    fform_gr <- gsub("^~", "", frandom, fixed = F)
  } else {
    fform_gr <- NULL
  }
  
  
  if (!is.null(grandom)) {
    gform_gr <- gsub("^~", "", grandom, fixed = F)
  } else {
    gform_gr <- NULL
  }
  
  if (!is.null(hrandom)) {
    hform_gr <- gsub("^~", "", hrandom, fixed = F)
  } else {
    hform_gr <- NULL
  }
  
  
  if (!is.null(irandom)) {
    iform_gr <- gsub("^~", "", irandom, fixed = F)
  } else {
    iform_gr <- NULL
  }
  
  if (!is.null(srandom)) {
    sform_gr <- gsub("^~", "", srandom, fixed = F)
  } else {
    sform_gr <- NULL
  }
  
  
  if (!is.null(sigmarandom)) {
    sigma_form_gr <- gsub("^~", "", sigmarandom, fixed = F)
  } else {
    sigma_form_gr <- NULL
  }
  
  
  
  if (!is.null(group_arg)) {
    gr_prefixss <- "gr"
    gr_varss <- group_arg$groupvar
    gr_by = group_arg$by
    gr_cov = group_arg$cov
    gr_dist = group_arg$dist
  } else {
    gr_prefixss <- NULL
    gr_varss    <- NULL
  }
  
  
  
  if (!is.null(randomsi)) {
    # these two arguments are set automaticaly if missing
    if (is.null(gr_prefixss))
      gr_prefixss <- 'gr'
    if (is.null(gr_varss))
      gr_varss <- id
    if (!is.null(group_arg)) {
      if (!is.null(gr_by))  {
        gr_byss2  <- paste0(",by=", gr_by)
      } else {
        gr_byss2 <- NULL
      }
      if (!is.null(gr_cov)) {
        gr_covss2  <- paste0(",cov=", gr_cov)
      } else {
        gr_covss2 <- NULL
      }
      if (!is.null(gr_dist)) {
        gr_distss2 <- paste0(",dist=", paste0("'", gr_dist, "'"))
      } else {
        gr_distss2 <- ''
      }
      gr__args <- paste0(gr_prefixss,
                         "(",
                         gr_varss,
                         gr_byss2,
                         gr_covss2,
                         gr_distss2,
                         ")")
      
      
      
      if (!is.null(aform) & grepl("a", randomsi, fixed = T)) {
        if(!arandom_wb_) {
          aform <- paste0(aform,
                          " + (", aform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(arandom_wb_) {
          aform <- paste0(aform,
                          " + (", arandom_wb, ")")
        }
      }
      
      
      if (!is.null(bform) & grepl("b", randomsi, fixed = T)) {
        if(!brandom_wb_) {
          bform <- paste0(bform,
                          " + (", bform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(brandom_wb_) {
          bform <- paste0(bform,
                          " + (", brandom_wb, ")")
        }
      }
      
      if (!is.null(cform) & grepl("c", randomsi, fixed = T)) {
        if(!crandom_wb_) {
          cform <- paste0(cform,
                          " + (", cform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(crandom_wb_) {
          cform <- paste0(cform,
                          " + (", crandom_wb, ")")
        }
      }
      
      if (!is.null(dform) & grepl("d", randomsi, fixed = T)) {
        if(!drandom_wb_) {
          dform <- paste0(dform,
                          " + (", dform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(drandom_wb_) {
          dform <- paste0(dform,
                          " + (", drandom_wb, ")")
        }
      }
      
      
      
      if (!is.null(eform) & grepl("e", randomsi, fixed = T)) {
        if(!erandom_wb_) {
          eform <- paste0(eform,
                          " + (", eform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(erandom_wb_) {
          eform <- paste0(eform,
                          " + (", erandom_wb, ")")
        }
      }
      
      
      
      if (!is.null(fform) & grepl("f", randomsi, fixed = T)) {
        if(!frandom_wb_) {
          fform <- paste0(fform,
                          " + (", fform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(frandom_wb_) {
          fform <- paste0(fform,
                          " + (", frandom_wb, ")")
        }
      }
      
      
      
      if (!is.null(gform) & grepl("g", randomsi, fixed = T)) {
        if(!grandom_wb_) {
          gform <- paste0(gform,
                          " + (", gform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(grandom_wb_) {
          gform <- paste0(gform,
                          " + (", grandom_wb, ")")
        }
      }
      
      
      if (!is.null(hform) & grepl("h", randomsi, fixed = T)) {
        if(!hrandom_wb_) {
          hform <- paste0(hform,
                          " + (", hform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(hrandom_wb_) {
          hform <- paste0(hform,
                          " + (", hrandom_wb, ")")
        }
      }
      
      
      if (!is.null(iform) & grepl("i", randomsi, fixed = T)) {
        if(!irandom_wb_) {
          iform <- paste0(iform,
                          " + (", iform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(irandom_wb_) {
          iform <- paste0(iform,
                          " + (", irandom_wb, ")")
        }
      }
      
      
      if (!is.null(sform) & grepl("s", randomsi, fixed = T)) {
        if(!srandom_wb_) {
          sform <- paste0(sform,
                          " + (", sform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(srandom_wb_) {
          sform <- paste0(sform,
                          " + (", srandom_wb, ")")
        }
      }
      
      
      
      ###########
      # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
      # In fact for df > 1, it forces 'd' to be random parameter only
      # allow random effect for 'd' even corresponding fixed effect is missing
      # only allowing for 'd'
      if(select_model == "sitar") {
        if (match_sitar_d_form) {
          if (is.null(dform) & grepl("d", randomsi, fixed = T)) {
            if(!drandom_wb_) {
              dform <- paste0(dform, "d ~ 0 + (1|", coridv , "|" , 
                              gr__args, ")")
            } else if(drandom_wb_) {
              dform <- paste0(dform, "d ~ 0 + (1|", coridv , "|" , 
                              gr__args, ")")
            }
          }
        }
      }
      ###########
      
    }
    
    
    
    
    
    if (!is.null(sigma_group_arg)) {
      sigma_gr_prefixss <- "gr"
      sigma_gr_varss <- sigma_group_arg$groupvar
      sigma_gr_by = sigma_group_arg$by
      sigma_gr_cov = sigma_group_arg$cov
      sigma_gr_dist = sigma_group_arg$dist
    } else {
      sigma_gr_prefixss <- NULL
      sigma_gr_varss    <- NULL
    }
    
    if (!is.null(sigma_formula_grsi)) {
      # these two arguments are set automaticaly if missing
      if (is.null(sigma_gr_prefixss))
        sigma_gr_prefixss <- 'gr'
      if (is.null(sigma_gr_varss))
        sigma_gr_varss <- id
      if (!is.null(sigma_group_arg)) {
        if (!is.null(sigma_gr_by))  {
          sigma_gr_byss2  <- paste0(",by=", sigma_gr_by)
        } else {
          sigma_gr_byss2 <- NULL
        }
        if (!is.null(sigma_gr_cov)) {
          sigma_gr_covss2  <- paste0(",cov=", sigma_gr_cov)
        } else {
          sigma_gr_covss2 <- NULL
        }
        if (!is.null(sigma_gr_dist)) {
          sigma_gr_distss2 <- paste0(",dist=", paste0("'", sigma_gr_dist, "'"))
        } else {
          sigma_gr_distss2 <- ''
        }
        sigma_gr__args <- paste0(sigma_gr_prefixss,
                                 "(",
                                 sigma_gr_varss,
                                 sigma_gr_byss2,
                                 sigma_gr_covss2,
                                 sigma_gr_distss2,
                                 ")")
        
        
        
        if (!is.null(sigmaform) & !is.null(sigma_formula_grsi)) {
          if(!sigma_random_wb_) {
            sigmaform <- paste0(sigmaform,
                                " + (", sigma_form_gr, "|", coridv, "|" , 
                                sigma_gr__args, ")")
          } else if(sigma_random_wb_) {
            sigmaform <- paste0(sigmaform,
                                " + (", sigma_random_wb, ")")
          }
        }
        ##
      }
    }
    
    if (is.null(sigma_group_arg)) {
      if (!is.null(sigmaform) & !is.null(sigma_formula_grsi)) {
        if(!sigma_random_wb_) {
          sigmaform <- paste0(sigmaform, " + (", sigma_form_gr, "|", 
                              coridv, "|" , id, ")")
        } else if(arandom_wb_) {
          sigmaform <- paste0(sigmaform, " + (", sigma_random_wb, ")")
        }
      }
    }
    
    
    
    
    if (is.null(group_arg)) {
      if (!is.null(aform) & grepl("a", randomsi, fixed = T)) {
        if(!arandom_wb_) { 
          aform <- paste0(aform, " + (", aform_gr, "|", coridv, "|" , id, ")")
        } else if(arandom_wb_) {
          aform <- paste0(aform, " + (", arandom_wb, ")")
        }
      }
      
      if (!is.null(bform) & grepl("b", randomsi, fixed = T)) {
        if(!brandom_wb_) { 
          bform <- paste0(bform, " + (", bform_gr, "|", coridv, "|" , id, ")")
        } else if(brandom_wb_) {
          bform <- paste0(bform, " + (", brandom_wb, ")")
        }
      }
      
      if (!is.null(cform) & grepl("c", randomsi, fixed = T)) {
        if(!crandom_wb_) { 
          cform <- paste0(cform, " + (", cform_gr, "|", coridv, "|" , id, ")")
        } else if(crandom_wb_) {
          cform <- paste0(cform, " + (", crandom_wb, ")")
        }
      }
      
      if (!is.null(dform) & grepl("d", randomsi, fixed = T)) {
        if(!drandom_wb_) { 
          dform <- paste0(dform, " + (", dform_gr, "|", coridv, "|" , id, ")")
        } else if(drandom_wb_) {
          dform <- paste0(dform, " + (", drandom_wb, ")")
        }
      }
      
      
      
      if (!is.null(eform) & grepl("e", randomsi, fixed = T)) {
        if(!erandom_wb_) { 
          eform <- paste0(eform, " + (", eform_gr, "|", coridv, "|" , id, ")")
        } else if(erandom_wb_) {
          eform <- paste0(eform, " + (", erandom_wb, ")")
        }
      }
      
      
      if (!is.null(fform) & grepl("f", randomsi, fixed = T)) {
        if(!frandom_wb_) { 
          fform <- paste0(fform, " + (", fform_gr, "|", coridv, "|" , id, ")")
        } else if(frandom_wb_) {
          fform <- paste0(fform, " + (", frandom_wb, ")")
        }
      }
      
      
      
      if (!is.null(gform) & grepl("g", randomsi, fixed = T)) {
        if(!grandom_wb_) { 
          gform <- paste0(gform, " + (", gform_gr, "|", coridv, "|" , id, ")")
        } else if(grandom_wb_) {
          gform <- paste0(gform, " + (", grandom_wb, ")")
        }
      }
      
      
      if (!is.null(hform) & grepl("h", randomsi, fixed = T)) {
        if(!hrandom_wb_) { 
          hform <- paste0(hform, " + (", hform_gr, "|", coridv, "|" , id, ")")
        } else if(hrandom_wb_) {
          hform <- paste0(hform, " + (", hrandom_wb, ")")
        }
      }
      
      
      if (!is.null(iform) & grepl("i", randomsi, fixed = T)) {
        if(!irandom_wb_) { 
          iform <- paste0(iform, " + (", iform_gr, "|", coridv, "|" , id, ")")
        } else if(irandom_wb_) {
          iform <- paste0(iform, " + (", irandom_wb, ")")
        }
      }
      
      if (!is.null(sform) & grepl("s", randomsi, fixed = T)) {
        if(!srandom_wb_) { 
          sform <- paste0(sform, " + (", sform_gr, "|", coridv, "|" , id, ")")
        } else if(srandom_wb_) {
          sform <- paste0(sform, " + (", srandom_wb, ")")
        }
      }
      
      
      ###########
      # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
      # In fact for df > 1, it forces 'd' to be random parameter only
      # allow random effect for 'd' even corresponding fixed effect is missing
      # only allowing for 'd'
      if(select_model == "sitar") {
        if (match_sitar_d_form) {
          if (is.null(dform) & grepl("d", randomsi, fixed = T)) {
            if(!drandom_wb_) {
              dform <- paste0(dform, "d ~ 0 + (1|", coridv , "|" , id, ")")
            } else if(drandom_wb_) {
              dform <- paste0(dform, "d ~ 0 + (1|", coridv , "|" , id, ")")
            }
          }
        }
      }
      ###########
      
    }
  }
  
  
  add_higher_level_str <- function(form, str) {
    if(!is.null(str[[1]])) {
      get_n_str <- strsplit(str, ")+(", fixed = T)[[1]][-1]
      get_n_str_length <- length(get_n_str)
    } else {
      get_n_str_length <- 0
    }
    if(get_n_str_length != 0) {
      str <- gsub(")+(", ")_(", str, fixed = T)
      str_ <- sub("^[^_]*_", "", str)
      str_ <- gsub(")_(", ")+(", str_, fixed = T)
      form <- paste0(form, "+", str_)
    } else {
      form <- form
    }
    form <- gsub("[[:space:]]", "", form)
    form
  }
  
  
  add_higher_level_str_id <- function(x) {
    extract_xx <- x
    extract_xx <- gsub("[[:space:]]", "", extract_xx)
    extract_xx <- strsplit(extract_xx, ")+(", fixed = T)[[1]][1]
    extract_xx <- gsub("\\)", "", extract_xx)
    gr_varss <- sub(".*\\|", "", extract_xx) 
    gr_varss
  }
  
  
  
  if(set_higher_levels) {
    if(a_formula_gr_strsi_present) {
      aform <- add_higher_level_str(aform, a_formula_gr_strsi)
      aform <- restore_paranthese_grgr_str_form(aform)
      agr_varss             <- add_higher_level_str_id(a_formula_gr_strsi)
      aform_gr_names        <- lapply(aform, get_x_random2)[[1]]
      aform_gr_names_asitis <- lapply(aform, get_x_random2_asitis)[[1]]
    } else {
      agr_varss <- aform_gr_names <- aform_gr_names_asitis <- NULL
    }
    
    if(b_formula_gr_strsi_present) {
      bform <- add_higher_level_str(bform, b_formula_gr_strsi)
      bform <- restore_paranthese_grgr_str_form(bform)
      bgr_varss             <- add_higher_level_str_id(b_formula_gr_strsi)
      bform_gr_names        <- lapply(bform, get_x_random2)[[1]]
      bform_gr_names_asitis <- lapply(bform, get_x_random2_asitis)[[1]]
    } else {
      bgr_varss <- bform_gr_names <- bform_gr_names_asitis <- NULL
    }
    
    if(c_formula_gr_strsi_present) {
      cform <- add_higher_level_str(cform, c_formula_gr_strsi)
      cform <- restore_paranthese_grgr_str_form(cform)
      cgr_varss             <- add_higher_level_str_id(c_formula_gr_strsi)
      cform_gr_names        <- lapply(cform, get_x_random2)[[1]]
      cform_gr_names_asitis <- lapply(cform, get_x_random2_asitis)[[1]]
    } else {
      cgr_varss <- cform_gr_names <- cform_gr_names_asitis <- NULL
    }
    
    if(d_formula_gr_strsi_present) {
      dform <- add_higher_level_str(dform, d_formula_gr_strsi)
      dform <- restore_paranthese_grgr_str_form(dform)
      dgr_varss             <- add_higher_level_str_id(d_formula_gr_strsi)
      dform_gr_names        <- lapply(dform, get_x_random2)[[1]]
      dform_gr_names_asitis <- lapply(dform, get_x_random2_asitis)[[1]]
    } else {
      dgr_varss <- dform_gr_names <- dform_gr_names_asitis <- NULL
    }
    
    if(e_formula_gr_strsi_present) {
      eform <- add_higher_level_str(eform, e_formula_gr_strsi)
      eform <- restore_paranthese_grgr_str_form(eform)
      egr_varss             <- add_higher_level_str_id(e_formula_gr_strsi)
      eform_gr_names        <- lapply(eform, get_x_random2)[[1]]
      eform_gr_names_asitis <- lapply(eform, get_x_random2_asitis)[[1]]
    } else {
      egr_varss <- eform_gr_names <- eform_gr_names_asitis <- NULL
    }
    
    if(f_formula_gr_strsi_present) {
      fform <- add_higher_level_str(fform, f_formula_gr_strsi)
      fform <- restore_paranthese_grgr_str_form(fform)
      fgr_varss             <- add_higher_level_str_id(f_formula_gr_strsi)
      fform_gr_names        <- lapply(fform, get_x_random2)[[1]]
      fform_gr_names_asitis <- lapply(fform, get_x_random2_asitis)[[1]]
    } else {
      fgr_varss <- fform_gr_names <- fform_gr_names_asitis <- NULL
    }
    
    
    if(g_formula_gr_strsi_present) {
      gform <- add_higher_level_str(gform, g_formula_gr_strsi)
      gform <- restore_paranthese_grgr_str_form(gform)
      ggr_varss             <- add_higher_level_str_id(g_formula_gr_strsi)
      gform_gr_names        <- lapply(gform, get_x_random2)[[1]]
      gform_gr_names_asitis <- lapply(gform, get_x_random2_asitis)[[1]]
    } else {
      ggr_varss <- gform_gr_names <- gform_gr_names_asitis <- NULL
    }
    
    if(h_formula_gr_strsi_present) {
      hform <- add_higher_level_str(hform, h_formula_gr_strsi)
      hform <- restore_paranthese_grgr_str_form(hform)
      hgr_varss             <- add_higher_level_str_id(h_formula_gr_strsi)
      hform_gr_names        <- lapply(hform, get_x_random2)[[1]]
      hform_gr_names_asitis <- lapply(hform, get_x_random2_asitis)[[1]]
    } else {
      hgr_varss <- hform_gr_names <- hform_gr_names_asitis <- NULL
    }
    
    if(i_formula_gr_strsi_present) {
      iform <- add_higher_level_str(iform, i_formula_gr_strsi)
      iform <- restore_paranthese_grgr_str_form(iform)
      igr_varss             <- add_higher_level_str_id(i_formula_gr_strsi)
      iform_gr_names        <- lapply(iform, get_x_random2)[[1]]
      iform_gr_names_asitis <- lapply(iform, get_x_random2_asitis)[[1]]
    } else {
      igr_varss <- iform_gr_names <- iform_gr_names_asitis <- NULL
    }
    
    if(s_formula_gr_strsi_present) {
      sform <- add_higher_level_str(sform, s_formula_gr_strsi)
      sform <- restore_paranthese_grgr_str_form(sform)
      sgr_varss             <- add_higher_level_str_id(s_formula_gr_strsi)
      sform_gr_names        <- lapply(sform, get_x_random2)[[1]]
      sform_gr_names_asitis <- lapply(sform, get_x_random2_asitis)[[1]]
    } else {
      sgr_varss <- sform_gr_names <- sform_gr_names_asitis <- NULL
    }
    
    gr_varss <- c(agr_varss, bgr_varss, cgr_varss, dgr_varss, 
                  egr_varss, fgr_varss)
    
    gr_varss <- unique(gr_varss)
    
    hierarchical_gr_names <- c(aform_gr_names, bform_gr_names, 
                               cform_gr_names, dform_gr_names, 
                               eform_gr_names, fform_gr_names,
                               gform_gr_names, hform_gr_names,
                               iform_gr_names, sform_gr_names)
    
    hierarchical_gr_names_asitis <- 
      c(aform_gr_names_asitis, bform_gr_names_asitis, 
        cform_gr_names_asitis, dform_gr_names_asitis,
        eform_gr_names_asitis, fform_gr_names_asitis,
        gform_gr_names_asitis, hform_gr_names_asitis,
        iform_gr_names_asitis, sform_gr_names_asitis)
    
    
    hierarchical_gr_names <- unique(hierarchical_gr_names)
    
    if(!is.null(gr_varss)) {
      if(length(unique(gr_varss)) != 1) 
        stop("Duplicate group identifiers ", gr_varss)
    }
    
    
    hierarchical_gr_names_asitis <- unique(hierarchical_gr_names_asitis)
    gr_varss <- hierarchical_gr_names_asitis
  } # if(set_higher_levels) {
  
  
  
  if(sigma_set_higher_levels) {
    
    if(sigma_formula_gr_strsi_present) {
      sigmaform <- add_higher_level_str(sigmaform, sigma_formula_gr_strsi)
      sigmaform <- restore_paranthese_grgr_str_form(sigmaform)
      sigmagr_varss             <- 
        add_higher_level_str_id(sigma_formula_gr_strsi)
      sigmaform_gr_names        <- lapply(sigmaform, get_x_random2)[[1]]
      sigmaform_gr_names_asitis <- lapply(sigmaform, get_x_random2_asitis)[[1]]
    } else {
      sigmaform_gr_names <- sigmaform_gr_names_asitis <- NULL
    }
    
    sigma_hierarchical_gr_names <- c(sigmaform_gr_names)
    sigma_hierarchical_gr_names <- unique(sigma_hierarchical_gr_names)
    
    sigma_hierarchical_gr_names_asitis <- c(sigmaform_gr_names_asitis)
    sigma_hierarchical_gr_names_asitis <- 
      unique(sigma_hierarchical_gr_names_asitis)
    
    sigma_gr_varss <- sigma_gr_varss 
    sigma_gr_varss_asitis <- sigmaform_gr_names_asitis
  } # if(sigma_set_higher_levels) {
  
  
  
  
  if(!set_higher_levels) hierarchical_gr_names <- NULL
  if(!set_higher_levels) hierarchical_gr_names_asitis <- NULL
  
  if(!sigma_set_higher_levels) sigma_hierarchical_gr_names <- NULL
  if(!sigma_set_higher_levels) sigma_hierarchical_gr_names_asitis <- NULL
  
  
  
  if(!is.null(sigma_formula_grsi)) {
    if(!sigma_set_higher_levels & grepl("|", 
                                        sigma_formula_grsi, fixed = TRUE)) {
      extract_xx <- sigma_formula_grsi
      extract_xx <- gsub("[[:space:]]", "", extract_xx)
      extract_xx <- strsplit(extract_xx, ")+(", fixed = T)[[1]][1]
      extract_xx <- gsub("\\)", "", extract_xx)
      sigma_gr_varss <- sub(".*\\|", "", extract_xx) 
      sigma_hierarchical_gr_names <- c(sigma_gr_varss)
    }
  } else if(!sigma_set_higher_levels) { 
    sigma_hierarchical_gr_names <- NULL
    sigma_hierarchical_gr_names_asitis <- NULL
  } else if(!sigma_set_higher_levels) {
    sigma_hierarchical_gr_names <- NULL
    sigma_hierarchical_gr_names_asitis <- NULL
  }
  
  
  
  
  if (!is.null(dpar_formulasi)) {
    if (!grepl("lf\\(", dpar_formulasi) |
        !grepl("nlf\\(", dpar_formulasi)) {
      dpar_covi_mat_form <- get_o_paranthesis(dpar_formulasi)
      dpar_covi_mat_form <- paste0("~", dpar_covi_mat_form)
    } else {
      dpar_covi_mat_form <- get_o_paranthesis2(dpar_formulasi)
      dpar_covi_mat_form <- paste0("~", dpar_covi_mat_form)
    }
  }
  
  
  
  if (!is.null(dpar_formulasi)) {
    dparcovmat <- eval(parse(
      text =
        paste0("model.matrix(",
               dpar_covi_mat_form, ",data = data)")
    ))
    
    
    if (ncol(dparcovmat) == 1) {
      dparncov <- NULL
    } else {
      dparncov <- ncol(dparcovmat) - 1
    }
    dparcovcoefnames <- colnames(dparcovmat)
    dparcovcoefnames <- gsub("\\(|)", "", dparcovcoefnames)
  }
  
  
  
  
  
  if (is.null(dpar_formulasi)) {
    dparncov <- NULL
    dparcovcoefnames <- NULL
  }
  
  # The brms replace equal sign = with EQ and removes comma from names
  
  if(!is.null(dparcovcoefnames)) {
    if(!is.null(dparcovcoefnames)) {
      dparcovcoefnames_c2 <- c()
      for (dparcovcoefnamesi in dparcovcoefnames) {
        dparcovcoefnamesi_c <- gsub("[[:space:]]", "", dparcovcoefnamesi)
        if(grepl("rcspline.eval", dparcovcoefnamesi_c)) {
          dparcovcoefnamesi_c <- gsub(",", "", dparcovcoefnamesi_c)
          dparcovcoefnamesi_c <- gsub("=", "EQ", dparcovcoefnamesi_c)
        } else {
          dparcovcoefnamesi_c <- dparcovcoefnamesi_c
        }
        dparcovcoefnames_c2 <- c(dparcovcoefnames_c2, dparcovcoefnamesi_c)
      }
      dparcovcoefnames <- dparcovcoefnames_c2
    }
  }
  
  
  if(select_model != "sitar" & select_model != "rcs" ) {
    abcform <-
      paste(cbind(aform, bform, cform, dform, eform, fform, 
                  gform, hform, iform), 
            collapse = ",")
  }
  
  
  if(select_model == "sitar" | select_model == "rcs" ) {
    abcform <-
      paste(cbind(aform, bform, cform, dform, eform, fform, 
                  gform, hform, iform), 
            collapse = ",")
  }
  
  
  
  
  
  # Imp that beyond this point the bform will be changed to combine
  
  if (!is.null(dpar_formulasi)) {
    if (!grepl("lf\\(", dpar_formulasi) &
        !grepl("nlf\\(", dpar_formulasi)) {
      if(!grepl('sigma~', dpar_formulasi, fixed = T)) {
        sform <- paste0(sform, ",", paste0("sigma", dpar_formulasi))
      } else if(grepl('sigma~', dpar_formulasi, fixed = T)) {
        sform <- paste0(sform, ",", paste0("", dpar_formulasi))
      }
    } else {
      sform <- sform
    }
  }
  
  
  # sform <- NULL
  
  if(!is.null(sigmaform)) {
    if(select_model == "sitar" | select_model == "rcs") {
      sform <- paste0(sform, ",", paste0("", sigmaform))
      # sform <- sigmaform
    }
    if(select_model != "sitar" & select_model != "rcs") {
      abcform <- paste0(abcform, ",", paste0("", sigmaform))
    }
  }
  
  
  if (!is.null(autocor_formi)) {
    if(select_model == "sitar" | select_model == "rcs") {
      sform <- paste0(sform,
                      ",",
                      paste0("autocor=", autocor_formi))
    }
    if(select_model != "sitar" & select_model != "rcs") {
      abcform <- paste0(abcform,
                        ",",
                        paste0("autocor=", autocor_formi))
    }
  }
  
  
  
  if(select_model == "sitar" | select_model == "rcs") {
    setbformula <- paste0("brms::bf(",
                          abcsformfit,
                          ", " ,
                          abcform,
                          ", " ,
                          sform,
                          ", " ,
                          "nl=TRUE,loop=FALSE)")
  }
  
  
  if(select_model != "sitar" & select_model != "rcs") {
    setbformula <- paste0("brms::bf(",
                          abcsformfit,
                          ", " ,
                          abcform,
                          ", " ,
                          # sform,
                          # ", " ,
                          "nl=TRUE,loop=FALSE)")
  }
  
  
  if (!is.null(unusedsi[[1]][1]) & unusedsi != "NULL") {
    setbformula <- paste0(bform, "unused=", unusedsi)
    setbformula <- gsub(")unused=", ",unused=", setbformula, fixed = T)
    setbformula <- paste0(setbformula, ")")
  }
  
  
  
  
  setbformula <- gsub("\\s", "", setbformula)
  
  
  
  if (!is.null(dpar_formulasi)) {
    if (grepl("lf\\(", dpar_formulasi) |
        grepl("nlf\\(", dpar_formulasi)) {
      if (!grepl("\\(~", dpar_formulasi) &
          strsplit(strsplit(dpar_formulasi, "~")[[1]][1],
                   "\\(")[[1]][2] != "sigma") {
        stop(
          "The distributional parameter name on the left hand side of ",
          "\n ",
          " lf/nlf formula (i.e., before ~) should be 'sigma'"
        )
      } else if (grepl("sigma~", dpar_formulasi, fixed = T)) {
        dpar_formulasi <- dpar_formulasi
      } else  {
        dpar_formulasi <- gsub("~", "sigma~", dpar_formulasi)
      }
      setbformula <- paste0(setbformula, "+", dpar_formulasi)
    } else {
      setbformula <- setbformula
    }
  }
  
  
  if (!is.null(familysi)) {
    setbformula <- paste0(setbformula, "+", familysi)
  }
  
  
  
  group_arg_groupvar <- gr_varss
  multivariate_rescor <-  multivariate$rescor
  univariate_by_by <- univariate_by$by
  
  
  if(sigma_set_higher_levels) { 
    sigma_arg_groupvar <- sigma_gr_varss_asitis
  } else if(!is.null(sigma_formula_grsi)) {
    if(!sigma_set_higher_levels & !grepl("|", 
                                         sigma_formula_grsi, fixed = TRUE)) {
      sigma_arg_groupvar <- sigma_gr_varss
    }
  } else if(!is.null(sigma_group_arg)) {
    sigma_arg_groupvar <- sigma_group_arg$groupvar
  } else {
    sigma_arg_groupvar <- NULL
  }
  
  
  if(!is.null(sigma_formula_grsi)) {
    if(!sigma_set_higher_levels & 
       grepl("|", sigma_formula_grsi, fixed = TRUE)) {
      sigma_arg_groupvar <- sigma_gr_varss
    } 
    if(!sigma_set_higher_levels & 
       !grepl("|", sigma_formula_grsi, fixed = TRUE)) {
      sigma_arg_groupvar <- gr_varss
    }  
  }
  
  # New add on 3 06 2023
  if(!sigma_set_higher_levels & !identical(sigma_gr_varss, gr_varss)) {
    sigma_arg_groupvar <- sigma_group_arg$groupvar
  }
  
  
  
  
  # Depending on select_model, assign null values to which not part of the model
  for (set_randomsi_higher_levsli in c(letters[1:26])) {
    set_nlpar_what <- set_randomsi_higher_levsli
    if(!exists(paste0(set_randomsi_higher_levsli, 'form'))) {
      assign(paste0(set_nlpar_what, 'covcoefnames'), NULL)
      assign(paste0(set_nlpar_what, 'covcoefnames_gr'), NULL)
      assign(paste0(set_nlpar_what, 'ncov'), NULL)
      assign(paste0(set_nlpar_what, 'ncov_gr'), NULL)
    } else if(is.null(ept(paste0(set_randomsi_higher_levsli, 'form')))) {
      assign(paste0(set_nlpar_what, 'covcoefnames'), NULL)
      assign(paste0(set_nlpar_what, 'covcoefnames_gr'), NULL)
      assign(paste0(set_nlpar_what, 'ncov'), NULL)
      assign(paste0(set_nlpar_what, 'ncov_gr'), NULL)
    } 
  } 
  
  
  
  
  # fit lm model
  a_covariate <- getcovlist(a_formulasi)
  b_covariate <- getcovlist(b_formulasi)
  c_covariate <- getcovlist(c_formulasi)
  d_covariate <- getcovlist(d_formulasi)
  e_covariate <- getcovlist(e_formulasi)
  f_covariate <- getcovlist(f_formulasi)
  g_covariate <- getcovlist(g_formulasi)
  h_covariate <- getcovlist(h_formulasi)
  i_covariate <- getcovlist(i_formulasi)
  s_covariate <- getcovlist(s_formulasi)
  
  sigma_covariate <- getcovlist(sigma_formulasi)
  
  covariates <- c(a_covariate, b_covariate, c_covariate, d_covariate, 
                  e_covariate, f_covariate, 
                  s_covariate)
  
  covariates_ <- unique(covariates)
  
  covariates_sigma_ <- unique(sigma_covariate)
  
  
  
  
  
  
  if(select_model == "sitar" |
     select_model == "rcs") {
    a_covariate_i <- c()
    if (grepl("\\*", a_formulasi)) {
      for (a_covariatei in a_covariate) {
        if (grepl("\\*", a_covariatei)) {
          t <- strsplit(a_covariatei, "\\*")[[1]]
          t <- c(paste0(t, collapse = "+"), paste0(t, collapse = ":"))
        } else {
          t <- a_covariatei
        }
        a_covariate_i <- c(a_covariate_i, t)
      }
    } else {
      a_covariate_i <- a_covariate
    }
    
    a_covariate <- a_covariate_i
    
    
    
    
    s_covariate_i <- c()
    if (grepl("\\*", s_formulasi)) {
      for (s_covariatei in s_covariate) {
        if (grepl("\\*", s_covariatei)) {
          t <- strsplit(s_covariatei, "\\*")[[1]]
          t <- c(paste0(t, collapse = "+"), paste0(t, collapse = ":"))
        } else {
          t <- s_covariatei
        }
        s_covariate_i <- c(s_covariate_i, t)
      }
    } else {
      s_covariate_i <- s_covariate
    }
    
    s_covariate <- s_covariate_i
    
    
    
    if (!grepl("^~1$", s_formulasi)) {
      if (!identical(strsplit(a_formulasi, "+")[[1]][1:2],
                     strsplit(s_formulasi, "+")[[1]][1:2])) {
        stop(
          "a_formula and s_formula should have the identical ",
          "\n ",
          " intercept structure i.e., ~ 1 or ~ 0"
        )
      }
      
      if (length(a_covariate) != length(s_covariate)) {
        stop(
          "s_formula formula should be intercept only or must have ",
          "\n ",
          " the same number of covariates as a_formula"
        )
      }
      
      if (length(a_covariate) == 1) {
        if (grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~0+",  a_covariate, "*", "mat_s"))
        } else if (!grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~1+",  a_covariate, "*", "mat_s"))
        }
      } else if (length(a_covariate) > 1) {
        main_cov <- a_covariate
        main_cov <-
          paste(unlist(strsplit(main_cov, "+", fixed = T)), sep = " ")
        inte_cov <- paste0(main_cov, ":", "mat_s")
        main_inte_cov <- c(main_cov, "mat_s", inte_cov)
        main_inte_cov <- paste(main_inte_cov, collapse = "+")
        if (grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~0+",  main_inte_cov))
        } else if (!grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~1+",  main_inte_cov))
        }
      }
    }
    
    if (grepl("^~1$", s_formulasi)) {
      if (length(a_covariate) == 1) {
        if (grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~0+",  a_covariate, "+", "mat_s"))
        } else if (!grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~1+",  a_covariate, "+", "mat_s"))
        }
      } else if (length(a_covariate) > 1) {
        main_cov <- a_covariate
        main_inte_cov <- c(main_cov, "mat_s")
        main_inte_cov <- paste(main_inte_cov, collapse = "+")
        if (grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~0+",  main_inte_cov))
        } else if (!grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~1+",  main_inte_cov))
        }
      }
    }
    
    
    if (grepl("^~1$", s_formulasi))
      s_covariate <- NULL
    
    
    if (grepl("^~1$", a_formulasi)) {
      if (grepl("~0", a_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~0+", "mat_s"))
      } else if (!grepl("~0", a_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~1+", "mat_s"))
      }
    }
    
    
    
    lm_fit  <- lm(lmform, data = data)
    lm_coef <- coef(lm_fit)
    
    lm_rsd  <- summary(lm_fit)$sigma
    
    
    # lme
    enverr. <- parent.frame()
    # err. <- FALSE
    assign('err.', FALSE, envir = enverr.)
    tryCatch(
      expr = {
        datalme <- data
        datalme[['mat_s']] <- eval(parse(text = 'mat_s'))
        randomlmer <- "Intercept"
        if (randomlmer == "slope")
          randomform <- paste0("~ 1 + ", x , " | ", id)
        if (randomlmer == "Intercept")
          randomform <- paste0("~ 1 ", " | ", id)
        randomform <- as.formula(randomform)
        lme_fit <-
          nlme::lme(fixed = lmform,
                    random = randomform,
                    data = datalme)
      },
      error = function(e) {
        # err. <<- TRUE
        assign('err.', TRUE, envir = enverr.)
      }
    )
    err. <- get('err.', envir = enverr.)
    if (err.) {
      lme_coef <- lm_coef
      lme_sd_a <- sd(predict(lm_fit))
      lme_rsd <- lm_rsd
    } else if (!err.) {
      lme_coef <- unname(nlme::fixed.effects(lme_fit))
      VarCorrnumeric <- nlme::VarCorr(lme_fit)[, 2] %>% as.numeric()
      lme_sd_a <- VarCorrnumeric[1]
      if (randomlmer == "Intercept") {
        lme_rsd <- VarCorrnumeric[2]
      } else if (randomlmer == "slope") {
        lme_rsd <- VarCorrnumeric[3]
      }
    }
    
    
    
    if (grepl("\\*", a_formulasi) & grepl("^~1$", s_formulasi)) {
      intercept_ <- lm_coef[!grepl("^mat_s", names(lm_coef))]
      spls_ <- lm_coef[grepl("^mat_s", names(lm_coef))]
      lm_coef <- c(intercept_, spls_)
    }
    
    
    if (grepl("\\*", a_formulasi) & !grepl("^~1$", s_formulasi)) {
      intercept_ <- lm_coef[!grepl("mat_s", names(lm_coef))]
      spls_ <- lm_coef[grepl("mat_s", names(lm_coef))]
      lm_coef <- c(intercept_, spls_)
    }
    
    
    
    
    
    
    lm_a_all <- lm_coef[1:ncol(acovmat)]
    if (!is.null(bfixed)) {
      lm_b_all <- rep(0, ncol(bcovmat))
    } else {
      lm_b_all <- NULL
    }
    if (!is.null(cfixed)) {
      lm_c_all <- rep(0, ncol(ccovmat))
    } else {
      lm_c_all <- NULL
    }
    if (!is.null(dfixed)) {
      lm_d_all <- rep(0, ncol(dcovmat))
    } else {
      lm_d_all <- NULL
    }
    
    lm_s_all <- lm_coef[(ncol(acovmat) + 1):length(lm_coef)]
    
    
    
    if (grepl("~1", a_formulasi, fixed = T)) {
      lm_a_all[1] <- lm_a_all[1] + lm_s_all[1] * min(knots)
    }
    
    
    names(lm_a_all) <- acovcoefnames
    names(lm_b_all) <- bcovcoefnames
    names(lm_c_all) <- ccovcoefnames
    names(lm_d_all) <- dcovcoefnames
    
    lm_a <- lm_a_all[1]
    lm_b <- lm_b_all[1]
    lm_c <- lm_c_all[1]
    lm_d <- lm_d_all[1]
    
    if (!is.null(ancov)) {
      lm_a_cov <- lm_a_all[2:length(lm_a_all)]
    } else {
      lm_a_cov <- NULL
    }
    if (!is.null(bncov)) {
      lm_b_cov <- lm_b_all[2:length(lm_b_all)]
    } else {
      lm_b_cov <- NULL
    }
    if (!is.null(cncov)) {
      lm_c_cov <- lm_c_all[2:length(lm_c_all)]
    } else {
      lm_c_cov <- NULL
    }
    if (!is.null(dncov)) {
      lm_d_cov <- lm_d_all[2:length(lm_d_all)]
    } else {
      lm_d_cov <- NULL
    }
    
    
    
    
    xnames <- names(lm_s_all)
    names_mat_s_scovmat <- c()
    for (x in xnames) {
      x.a <- rev(strsplit(x, "_")[[1]])
      x.a <- gsub("^mat", "Intercept", x.a)
      x.a <- gsub(":mat", "", x.a)
      x.a <- paste(x.a, collapse = "_")
      names_mat_s_scovmat <- c(names_mat_s_scovmat, x.a)
    }
    
    if (grepl("^~1$", s_formulasi)) {
      mat_s_scovmat <- mat_s
    } else if (!grepl("^~1$", s_formulasi)) {
      mat_s_scovmat <- model.matrix(lm_fit)
      mat_s_scovmat <-
        mat_s_scovmat[, (ncol(acovmat) + 1):ncol(mat_s_scovmat)]
    }
    sds_X <- rep(0, ncol(mat_s_scovmat))
    
    for (i in 1:ncol(mat_s_scovmat)) {
      sds_X[i] = sd(mat_s_scovmat[, i])
    }
    lm_sdx_all <- sd(data[[y]]) / sds_X
    
    
    names(lm_s_all) <- names_mat_s_scovmat
    names(lm_sdx_all) <- names_mat_s_scovmat
    
    
    lm_s   <- lm_s_all[1:(nknots - 1)]
    lm_sdx <- lm_sdx_all[1:(nknots - 1)]
    
    
    if (!is.null(s_covariate) & length(s_covariate) > 1) {
      lm_s_cov <- lm_s_all[nknots:length(lm_s_all)]
      lm_sdx_cov <- lm_sdx_all[nknots:length(lm_sdx_all)]
      
      inname_c_all <- c()
      for (inname in paste0("s", 1:df)) {
        t <- names(lm_a_all)[2:length(names(lm_a_all))]
        inname_c_all <- c(inname_c_all, paste0(inname, "_", t))
      }
      lm_s_cov   <- lm_s_cov[order(factor(names(lm_s_cov),
                                          levels = inname_c_all))]
      lm_sdx_cov <- lm_sdx_cov[order(factor(names(lm_sdx_cov),
                                            levels = inname_c_all))]
      
      lm_s_all <- lm_sdx_all <- c()
      for (idfi in 1:df) {
        lm_s_all <- c(lm_s_all, c(lm_s[idfi],
                                  lm_s_cov[grep(paste0("s", idfi, "_"),
                                                names(lm_s_cov))]))
        lm_sdx_all <- c(lm_sdx_all, c(lm_sdx[idfi],
                                      lm_sdx_cov[grep(paste0("s", idfi, "_"),
                                                      names(lm_sdx_cov))]))
      }
    } else if ((!is.null(s_covariate) &
                length(s_covariate) == 1) |
               is.null(s_covariate)) {
      if (!grepl("~0", s_formulasi, fixed = T)) {
        lm_s_all <- lm_coef[(ncol(acovmat) + 1):length(lm_coef)]
        names(lm_s_all) <- names_mat_s_scovmat
        names(lm_sdx_all) <- names_mat_s_scovmat
        lm_s <- lm_s_all[1:(nknots - 1)]
        lm_sdx <- lm_sdx_all[1:(nknots - 1)]
        if (length(lm_s_all) > (nknots - 1)) {
          lm_s_cov <- lm_s_all[nknots:length(lm_s_all)]
          lm_sdx_cov <- lm_sdx_all[nknots:length(lm_sdx_all)]
          tnames_s <-
            names_mat_s_scovmat[nknots:length(names_mat_s_scovmat)]
          names(lm_s_cov) <- tnames_s
          tnames_sdx <-
            names_mat_s_scovmat[nknots:length(names_mat_s_scovmat)]
          names(lm_sdx_cov) <- tnames_sdx
        } else {
          lm_s_cov <- NULL
          lm_sdx_cov <- NULL
        }
      }
      
      
      if (grepl("~0", s_formulasi, fixed = T)) {
        lm_s_all <- lm_coef[(ncol(acovmat) + 1):length(lm_coef)]
        names(lm_s_all) <- names_mat_s_scovmat
        names(lm_sdx_all) <- names_mat_s_scovmat
        inname_c_all <- c()
        for (inname in paste0("s", 1:df)) {
          t <- names(lm_a_all)[1:length(names(lm_a_all))]
          inname_c_all <- c(inname_c_all, paste0(inname, "_", t))
        }
        lm_s <- NULL
        lm_s_cov <- NULL
        lm_sdx <- NULL
        lm_sdx_cov <- NULL
      }
    }
    
    if (any(is.na(lm_coef))) {
      stop(
        "Inclusion of covariates resulted in rank-deficient design  matrix",
        "\n ",
        "(with some NA coefficients for the 'lm' model fit)",
        "\n ",
        "Please simplyfy the model"
      )
    }
  } # if(select_model == 'sitar' | select_model == "rcs")) {
  
  
  
  
  
  
  
  
  
  if(select_model != 'sitar' & select_model != "rcs") {
    lm_data_at_max_x <- data
    
    if (grepl("~0", a_formulasi, fixed = T)) {
      lmform  <- as.formula(paste0(y, "~0+", 'mat_s'))
    } else if (!grepl("~0", a_formulasi, fixed = T)) {
      lmform  <- as.formula(paste0(y, "~1+", 'mat_s'))
    }
    
    lm_data_at_max_x[['mat_s']] <- 
      lm_data_at_max_x[[x]] - max(lm_data_at_max_x[[x]])
    
    lm_fit  <- lm(lmform, data = lm_data_at_max_x)
    lm_coef <- coef(lm_fit)
    lm_rsd  <- summary(lm_fit)$sigma
    
    if (any(is.na(lm_coef))) {
      stop(
        "Inclusion of covariates resulted in rank-deficient design  matrix",
        "\n ",
        "(with some NA coefficients for the 'lm' model fit)",
        "\n ",
        "Please simplyfy the model"
      )
    }
    
    # lme
    enverr. <- parent.frame()
    # err. <- FALSE
    assign('err.', FALSE, envir = enverr.)
    tryCatch(
      expr = {
        datalme <- data
        datalme[['mat_s']] <- eval(parse(text = 'mat_s'))
        randomlmer <- "Intercept"
        if (randomlmer == "slope")
          randomform <- paste0("~ 1 + ", x , " | ", id)
        if (randomlmer == "Intercept")
          randomform <- paste0("~ 1 ", " | ", id)
        randomform <- as.formula(randomform)
        lme_fit <-
          nlme::lme(fixed = lmform,
                    random = randomform,
                    data = datalme)
      },
      error = function(e) {
        # err. <<- TRUE
        assign('err.', TRUE, envir = enverr.)
      }
    )
    err. <- get('err.', envir = enverr.)
    if (err.) {
      lme_coef <- lm_coef
      lme_sd_a <- sd(predict(lm_fit))
      lme_rsd <- lm_rsd
    } else if (!err.) {
      lme_coef <- unname(nlme::fixed.effects(lme_fit))
      VarCorrnumeric <- nlme::VarCorr(lme_fit)[, 2] %>% as.numeric()
      lme_sd_a <- VarCorrnumeric[1]
      if (randomlmer == "Intercept") {
        lme_rsd <- VarCorrnumeric[2]
      } else if (randomlmer == "slope") {
        lme_rsd <- VarCorrnumeric[3]
      }
    }
    
    if (grepl("\\*", a_formulasi) & grepl("^~1$", a_formulasi)) {
      intercept_ <- lm_coef[!grepl("^mat_s", names(lm_coef))]
      spls_ <- lm_coef[grepl("^mat_s", names(lm_coef))]
      lm_coef <- c(intercept_, spls_)
    }
    
    if (grepl("\\*", a_formulasi) & !grepl("^~1$", a_formulasi)) {
      intercept_ <- lm_coef[!grepl("mat_s", names(lm_coef))]
      spls_ <- lm_coef[grepl("mat_s", names(lm_coef))]
      lm_coef <- c(intercept_, spls_)
    }
    
    lm_a_all <- lm_coef[1:ncol(acovmat)]
    if (!is.null(bfixed)) {
      if(grepl("^pb", select_model)) lm_b_all <- lm_a_all * 0.9
      if(!grepl("^pb", select_model)) lm_b_all <- lm_a_all * 0.9
    } else {
      lm_b_all <- NULL
    }
    if (!is.null(cfixed)) {
      lm_c_all <- rep(0, ncol(ccovmat))
    } else {
      lm_c_all <- NULL
    }
    if (!is.null(dfixed)) {
      lm_d_all <- rep(0, ncol(dcovmat))
    } else {
      lm_d_all <- NULL
    }
    
    if (!is.null(efixed)) {
      lm_e_all <- rep(0, ncol(ecovmat))
    } else {
      lm_e_all <- NULL
    }
    
    if (!is.null(ffixed)) {
      lm_f_all <- rep(0, ncol(fcovmat))
    } else {
      lm_f_all <- NULL
    }
    
    if (!is.null(gfixed)) {
      lm_g_all <- rep(0, ncol(gcovmat))
    } else {
      lm_g_all <- NULL
    }
    
    if (!is.null(hfixed)) {
      lm_h_all <- rep(0, ncol(hcovmat))
    } else {
      lm_h_all <- NULL
    }
    
    if (!is.null(ifixed)) {
      lm_i_all <- rep(0, ncol(icovmat))
    } else {
      lm_i_all <- NULL
    }
    
    names(lm_a_all) <- acovcoefnames
    names(lm_b_all) <- bcovcoefnames
    names(lm_c_all) <- ccovcoefnames
    names(lm_d_all) <- dcovcoefnames
    names(lm_e_all) <- ecovcoefnames
    names(lm_f_all) <- fcovcoefnames
    names(lm_g_all) <- hcovcoefnames
    names(lm_h_all) <- hcovcoefnames
    names(lm_i_all) <- icovcoefnames
    
    lm_a <- lm_a_all[1]
    lm_b <- lm_b_all[1]
    lm_c <- lm_c_all[1]
    lm_d <- lm_d_all[1]
    lm_e <- lm_e_all[1]
    lm_f <- lm_f_all[1]
    lm_g <- lm_h_all[1]
    lm_h <- lm_h_all[1]
    lm_i <- lm_i_all[1]
    
    if (!is.null(ancov)) {
      lm_a_cov <- lm_a_all[2:length(lm_a_all)]
    } else {
      lm_a_cov <- NULL
    }
    if (!is.null(bncov)) {
      lm_b_cov <- lm_b_all[2:length(lm_b_all)]
    } else {
      lm_b_cov <- NULL
    }
    if (!is.null(cncov)) {
      lm_c_cov <- lm_c_all[2:length(lm_c_all)]
    } else {
      lm_c_cov <- NULL
    }
    if (!is.null(dncov)) {
      lm_d_cov <- lm_d_all[2:length(lm_d_all)]
    } else {
      lm_d_cov <- NULL
    }
    if (!is.null(encov)) {
      lm_e_cov <- lm_e_all[2:length(lm_e_all)]
    } else {
      lm_e_cov <- NULL
    }
    if (!is.null(fncov)) {
      lm_f_cov <- lm_f_all[2:length(lm_f_all)]
    } else {
      lm_f_cov <- NULL
    }
    if (!is.null(gncov)) {
      lm_g_cov <- lm_g_all[2:length(lm_g_all)]
    } else {
      lm_g_cov <- NULL
    }
    if (!is.null(hncov)) {
      lm_h_cov <- lm_h_all[2:length(lm_h_all)]
    } else {
      lm_h_cov <- NULL
    }
    if (!is.null(incov)) {
      lm_i_cov <- lm_i_all[2:length(lm_i_all)]
    } else {
      lm_i_cov <- NULL
    }
    
  } # if(select_model != 'sitar' & select_model != "rcs") {
  
  
  
  
  
  #######################################################
  
  for (set_randomsi_higher_levsli in c(letters[1:26])) {
    set_nlpar_what <- set_randomsi_higher_levsli
    if(!exists(paste0('lm_', set_randomsi_higher_levsli, '_all'))) {
      assign(paste0('lm_', set_randomsi_higher_levsli, '_all'), NULL)
      assign(paste0('lm_', set_randomsi_higher_levsli, '_cov'), NULL)
      assign(paste0('lm_', set_randomsi_higher_levsli, ''), NULL)
    } else if(exists(paste0('lm_', set_randomsi_higher_levsli, '_all'))) {
      if(is.null(ept(paste0('lm_', set_randomsi_higher_levsli, '_all')))) {
        assign(paste0('lm_', set_randomsi_higher_levsli, '_all'), NULL)
        assign(paste0('lm_', set_randomsi_higher_levsli, '_cov'), NULL)
        assign(paste0('lm_', set_randomsi_higher_levsli, ''), NULL)
      } # if(is.null(ept(paste0('f', 'form')))) {
    }
  } # for (set_randomsi_higher_levsli in c(letters[1:26])) {
  
  if(is.null(lm_s_all)) {
    assign(paste0('lm_', 'sdx'), NULL)
    assign(paste0('lm_', 'sdx', '_all'), NULL)
    assign(paste0('lm_', 'sdx', '_cov'), NULL)
  }
  
  # brms removes white spaces from the coefficient names
  # mimicking that behaviors but keeping it separate here as 
  # brms may later change it to underscore or something else
  
  gsubitbt <- ""
  acovcoefnames <- gsub("[[:space:]]", gsubitbt, acovcoefnames)
  bcovcoefnames <- gsub("[[:space:]]", gsubitbt, bcovcoefnames)
  ccovcoefnames <- gsub("[[:space:]]", gsubitbt, ccovcoefnames)
  dcovcoefnames <- gsub("[[:space:]]", gsubitbt, dcovcoefnames)
  ecovcoefnames <- gsub("[[:space:]]", gsubitbt, ecovcoefnames)
  fcovcoefnames <- gsub("[[:space:]]", gsubitbt, fcovcoefnames)
  gcovcoefnames <- gsub("[[:space:]]", gsubitbt, gcovcoefnames)
  hcovcoefnames <- gsub("[[:space:]]", gsubitbt, hcovcoefnames)
  icovcoefnames <- gsub("[[:space:]]", gsubitbt, icovcoefnames)
  
  scovcoefnames <- gsub("[[:space:]]", gsubitbt, scovcoefnames)
  
  acovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, acovcoefnames_gr)
  bcovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, bcovcoefnames_gr)
  ccovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, ccovcoefnames_gr)
  dcovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, dcovcoefnames_gr)
  ecovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, ecovcoefnames_gr)
  fcovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, fcovcoefnames_gr)
  gcovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, gcovcoefnames_gr)
  hcovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, hcovcoefnames_gr)
  icovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, icovcoefnames_gr)
  scovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, scovcoefnames_gr)
  
  dparcovcoefnames <- gsub("[[:space:]]", gsubitbt, dparcovcoefnames)
  
  sigmacovcoefnames <- gsub("[[:space:]]", gsubitbt, sigmacovcoefnames)
  sigmacovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, sigmacovcoefnames_gr)
  
  acovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, acovcoefnames_gr_str)
  bcovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, bcovcoefnames_gr_str)
  ccovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, ccovcoefnames_gr_str)
  dcovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, dcovcoefnames_gr_str)
  ecovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, ecovcoefnames_gr_str)
  fcovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, fcovcoefnames_gr_str)
  gcovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, gcovcoefnames_gr_str)
  hcovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, hcovcoefnames_gr_str)
  icovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, icovcoefnames_gr_str)
  scovcoefnames_gr_str <- gsub("[[:space:]]", gsubitbt, scovcoefnames_gr_str)
  
  
  
  list_out <- list(
    aform = aform, 
    bform = bform, 
    cform = cform, 
    dform = dform, 
    eform = eform, 
    fform = fform,
    gform = gform,
    hform = hform,
    iform = iform,
    sform = sform, 
    
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
    
    acovcoefnames = acovcoefnames,
    bcovcoefnames = bcovcoefnames,
    ccovcoefnames = ccovcoefnames,
    dcovcoefnames = dcovcoefnames,
    ecovcoefnames = ecovcoefnames,
    fcovcoefnames = fcovcoefnames,
    gcovcoefnames = gcovcoefnames,
    hcovcoefnames = hcovcoefnames,
    icovcoefnames = icovcoefnames,
    scovcoefnames = scovcoefnames,
    
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
    
    acovcoefnames_gr = acovcoefnames_gr,
    bcovcoefnames_gr = bcovcoefnames_gr,
    ccovcoefnames_gr = ccovcoefnames_gr,
    dcovcoefnames_gr = dcovcoefnames_gr,
    ecovcoefnames_gr = ecovcoefnames_gr,
    fcovcoefnames_gr = fcovcoefnames_gr,
    gcovcoefnames_gr = gcovcoefnames_gr,
    hcovcoefnames_gr = hcovcoefnames_gr,
    icovcoefnames_gr = icovcoefnames_gr,
    scovcoefnames_gr = scovcoefnames_gr,
    
    ancov_gr_str = ancov_gr_str,
    bncov_gr_str = bncov_gr_str,
    cncov_gr_str = cncov_gr_str,
    dncov_gr_str = dncov_gr_str,
    encov_gr_str = encov_gr_str,
    fncov_gr_str = fncov_gr_str,
    gncov_gr_str = gncov_gr_str,
    hncov_gr_str = hncov_gr_str,
    incov_gr_str = incov_gr_str,
    sncov_gr_str = sncov_gr_str,
    
    acovcoefnames_gr_str    = acovcoefnames_gr_str,
    bcovcoefnames_gr_str    = bcovcoefnames_gr_str,
    ccovcoefnames_gr_str    = ccovcoefnames_gr_str,
    dcovcoefnames_gr_str    = dcovcoefnames_gr_str,
    ecovcoefnames_gr_str    = ecovcoefnames_gr_str,
    fcovcoefnames_gr_str    = fcovcoefnames_gr_str,
    gcovcoefnames_gr_str    = gcovcoefnames_gr_str,
    hcovcoefnames_gr_str    = hcovcoefnames_gr_str,
    icovcoefnames_gr_str    = icovcoefnames_gr_str,
    scovcoefnames_gr_str    = scovcoefnames_gr_str,
    
    acovcoefnames_gr_str_id = acovcoefnames_gr_str_id,
    bcovcoefnames_gr_str_id = bcovcoefnames_gr_str_id,
    ccovcoefnames_gr_str_id = ccovcoefnames_gr_str_id,
    dcovcoefnames_gr_str_id = dcovcoefnames_gr_str_id,
    ecovcoefnames_gr_str_id = ecovcoefnames_gr_str_id,
    fcovcoefnames_gr_str_id = fcovcoefnames_gr_str_id,
    gcovcoefnames_gr_str_id = gcovcoefnames_gr_str_id,
    hcovcoefnames_gr_str_id = hcovcoefnames_gr_str_id,
    icovcoefnames_gr_str_id = icovcoefnames_gr_str_id,
    scovcoefnames_gr_str_id = scovcoefnames_gr_str_id,
    
    acovcoefnames_gr_str_form = acovcoefnames_gr_str_form,
    bcovcoefnames_gr_str_form = bcovcoefnames_gr_str_form,
    ccovcoefnames_gr_str_form = ccovcoefnames_gr_str_form,
    dcovcoefnames_gr_str_form = dcovcoefnames_gr_str_form,
    ecovcoefnames_gr_str_form = ecovcoefnames_gr_str_form,
    fcovcoefnames_gr_str_form = fcovcoefnames_gr_str_form,
    gcovcoefnames_gr_str_form = gcovcoefnames_gr_str_form,
    hcovcoefnames_gr_str_form = hcovcoefnames_gr_str_form,
    icovcoefnames_gr_str_form = icovcoefnames_gr_str_form,
    scovcoefnames_gr_str_form = scovcoefnames_gr_str_form,
    
    dparncov             = dparncov,
    sigmancov            = sigmancov,
    sigmacovcoefnames    = sigmacovcoefnames,
    sigmancov_gr         = sigmancov_gr,
    sigmacovcoefnames_gr = sigmacovcoefnames_gr,
    
    sigmancov_gr_str            = sigmancov_gr_str,
    sigmacovcoefnames_gr_str    = sigmacovcoefnames_gr_str,
    sigmacovcoefnames_gr_str_id = sigmacovcoefnames_gr_str_id,
    sigmacovcoefnames_gr_str_form = sigmacovcoefnames_gr_str_form,
    
    
    gr_str_unique_id = gr_str_unique_id,
    gr_str_id_all    = gr_str_id_all,
    gr_str_form_all  = gr_str_form_all, 
    gr_str_coef_all  = gr_str_coef_all,
    gr_str_ncov_all  = gr_str_ncov_all,
    gr_str_corr_all  = gr_str_corr_all,
    gr_str_corr_tf   = gr_str_corr_tf,
    
    sigma_str_unique_id = sigma_str_unique_id,
    sigma_str_id_all    = sigma_str_id_all,
    sigma_str_form_all  = sigma_str_form_all, 
    sigma_str_coef_all  = sigma_str_coef_all,
    sigma_str_ncov_all  = sigma_str_ncov_all,
    sigma_str_corr_all  = sigma_str_corr_all,
    sigma_str_corr_tf   = sigma_str_corr_tf,
    
    
    dparcovcoefnames = dparcovcoefnames,
    group_arg_groupvar = group_arg_groupvar,
    sigma_arg_groupvar = sigma_arg_groupvar,
    multivariate_rescor = multivariate_rescor,
    univariate_by_by = univariate_by_by,
    set_higher_levels = set_higher_levels,
    hierarchical_gr_names = hierarchical_gr_names,
    hierarchical_gr_names_asitis = hierarchical_gr_names_asitis,
    sigma_hierarchical_gr_names = sigma_hierarchical_gr_names,
    sigma_hierarchical_gr_names_asitis = sigma_hierarchical_gr_names_asitis,
    covariates_ = covariates_,
    covariates_sigma_ = covariates_sigma_,
    lm_a_all = lm_a_all,
    lm_b_all = lm_b_all,
    lm_c_all = lm_c_all,
    lm_d_all = lm_d_all,
    lm_e_all = lm_e_all,
    lm_f_all = lm_f_all,
    lm_g_all = lm_g_all,
    lm_h_all = lm_h_all,
    lm_i_all = lm_i_all,
    
    lm_a = lm_a,
    lm_b = lm_b,
    lm_c = lm_c,
    lm_d = lm_d,
    lm_e = lm_e,
    lm_f = lm_f,
    lm_g = lm_g,
    lm_h = lm_h,
    lm_i = lm_i,
    
    lm_a_cov = lm_a_cov,
    lm_b_cov = lm_b_cov,
    lm_c_cov = lm_c_cov,
    lm_d_cov = lm_d_cov,
    lm_e_cov = lm_e_cov,
    lm_f_cov = lm_f_cov,
    lm_h_cov = lm_g_cov,
    lm_h_cov = lm_h_cov,
    lm_i_cov = lm_i_cov,
    
    lm_s = lm_s,
    lm_s_cov = lm_s_cov,
    lm_s_all = lm_s_all,
    lm_sdx = lm_sdx,
    lm_sdx_cov = lm_sdx_cov,
    lm_sdx_all = lm_sdx_all,
    lm_rsd = lm_rsd,
    lme_sd_a = lme_sd_a,
    lme_rsd = lme_rsd
  )
  
  attr(setbformula, "list_out") <- as.list(list_out)
  
  return(setbformula)
}









#' An internal function to prepare Stan function for Bayesian SITAR growth
#' curve model
#'
#' The \code{prepare_function}) constructs custom Stan function  which is passed
#' on to the [bsitar::bsitar()] function. For univariate-by- subgroup model
#' (\code{univariate_by}) and multivariate (\code{multivariate}) models (see
#' [bsitar::bsitar()]), the \code{x}, \code{y}, \code{id}, \code{knots},
#' \code{nknots}, are automatically matched with the sub-models.
#'
#' @param x Predictor variable in the data. See [bsitar::bsitar()] for details.
#'
#' @param y Response variable in the data. See [bsitar::bsitar()] for details.
#'
#' @param id A vector specifying a unique group identifier for each individual.
#' See [bsitar::bsitar()] for details.
#'
#' @param knots A vector of knots used for constructing the spline design
#' matrix. See [bsitar::bsitar()] for details.
#'
#' @param nknots An integer specifying the number of knots.
#'
#' @param data Data frame containing variables \code{x}, \code{y} and \code{id}.
#'
#' @param internal_function_args Internal arguments passed from the
#'   [bsitar::bsitar()] to the \code{prepare_formula}).
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
prepare_function <- function(x,
                             y,
                             id,
                             knots,
                             nknots,
                             data,
                             internal_function_args) {

  # Initiate non formalArgs()
  brms_arguments <- NULL;
  xfunsi <- NULL;
  Var1 <- NULL;
  Var2 <- NULL;
  select_model <- NULL;
  fixedsi <- NULL;
  match_sitar_d_form <- NULL;
  d_adjustedsi <- NULL;
  randomsi <- NULL;
  getxname <- NULL;
  getknotsname <- NULL;
  spfncname <- NULL;
  xoffset <- NULL;
  yfunsi <- NULL;
  all_raw_str <- NULL;
  all_raw_str <- NULL;
  decomp <- NULL;
  nys <- NULL;
  gsub_out_unscaled <- NULL;
  checkscovsi <- NULL;
  add_rcsfunmatqrinv_genquant <- NULL;
  add_b_Qr_genquan_s_coef <- NULL;
  
  
  
  if (!is.null(internal_function_args)) {
    eout <- list2env(internal_function_args)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  
  
  
  backend <- eval(brms_arguments$backend)

  vector_X_name <- "Xp"
  
  if (nys == 1)
    resp_ <- ""
  if (nys > 1)
    resp_ <- paste0("_", y)
  
  XR_inv_name <- 'XR_inv'
  XR_inv_name_resp <- paste0(XR_inv_name, resp_)
  
  b_sx_name <- 'b_sx'
  b_sx_name_resp <- paste0('b', resp_, "_", 'sx')
  
  
  iysxi <- 's'
  setp <- paste0(iysxi, '')
  setnlp <- paste0('nlp', resp_, '_', iysxi)
  setnlp_vector <- paste0('b', resp_)
  
  
  if (resp_ == "") {
    szx_name <- paste0('s', 1:(nknots - 1))
    szx_name_resp <- paste0('b', "_", 's', 1:(nknots - 1))
  } else {
    szx_name <- paste0('s', 1:(nknots - 1))
    szx_name_resp <- paste0('b', resp_, "_", szx_name)
  }
  
  
  b_s_name <- szx_name_resp
  
  szxbq_vector <- paste0('v', resp_, "_", 's', 'x')
  
  #########
  
  if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
    if (xfunsi == "log") {
      tranform_x_int <- 1
    } else if (xfunsi == "sqrt") {
      tranform_x_int <- 2
    } else if (xfunsi != "log" | xfunsi == "sqrt") {
      tranform_x_int <- 0
    }
  } else if (is.null(xfunsi[[1]][1]) | xfunsi == "NULL") {
    tranform_x_int <- 0
  }
  
  
  
  
  
  #########
  set_x_y_scale_factror <- function(xfunsi = NULL,
                                    yfunsi = NULL,
                                    tranformations = "identity") {
    scale_set_comb <- tranformations
    scale_set_comb1 <-
      paste(scale_set_comb, scale_set_comb, sep = "_")
    scale_set_comb2 <-
      with(subset(expand.grid(scale_set_comb, scale_set_comb),
                  Var1 != Var2), paste0(Var1, '_', Var2))
    scale_set_comb <- c(scale_set_comb1, scale_set_comb2)
    
    if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
      if (xfunsi == "log") {
        xscale_set <- "log"
      } else if (xfunsi == "sqrt") {
        xscale_set <- "sqrt"
      } else if (xfunsi != "log" | xfunsi == "sqrt") {
        xscale_set <- "identity"
      }
    } else if (is.null(xfunsi[[1]][1]) | xfunsi == "NULL") {
      xscale_set <- "identity"
    }
    
    if (!is.null(yfunsi[[1]][1]) & yfunsi != "NULL") {
      if (yfunsi == "log") {
        yscale_set <- "log"
      } else if (yfunsi == "sqrt") {
        yscale_set <- "sqrt"
      } else if (yfunsi != "log" | yfunsi == "sqrt") {
        yscale_set <- "identity"
      }
    } else if (is.null(yfunsi[[1]][1]) | yfunsi == "NULL") {
      yscale_set <- "identity"
    }
    
    
    
    
    if (xscale_set == "identity" & yscale_set == "identity") {
      xscale_factor_str_d1 <- "rep_vector(1, N);"
      xscale_factor_str_d2 <- "rep_vector(1, N);"
      yscale_factor_str_d1 <- "rep_vector(1, N);"
      yscale_factor_str_d2 <- "rep_vector(1, N);"
    } else if (xscale_set == "log" & yscale_set == "log") {
      xscale_factor_str_d1 <- "exp(Xm + xoffset);"
      xscale_factor_str_d2 <- "exp(Xm + xoffset);"
      yscale_factor_str_d1 <- "(pred_d0);"
      yscale_factor_str_d2 <- "(pred_d0);"
    } else if (xscale_set == "sqrt" & yscale_set == "sqrt") {
      xscale_factor_str_d1 <- "(Xm + xoffset);"
      xscale_factor_str_d2 <- "(Xm + xoffset);"
      yscale_factor_str_d1 <- "(sqrt(pred_d0));"
      yscale_factor_str_d2 <- "(sqrt(pred_d0));"
    } else if (xscale_set == "log" & yscale_set == "identity") {
      xscale_factor_str_d1 <- "exp(Xm + xoffset);"
      xscale_factor_str_d2 <- "exp(Xm + xoffset);"
      yscale_factor_str_d1 <- "rep_vector(1, N);"
      yscale_factor_str_d2 <- "rep_vector(1, N);"
    } else if (xscale_set == "sqrt" & yscale_set == "identity") {
      xscale_factor_str_d1 <- "(Xm + xoffset);"
      xscale_factor_str_d2 <- "(Xm + xoffset);"
      yscale_factor_str_d1 <- "rep_vector(0.5, N);"
      yscale_factor_str_d2 <- "rep_vector(0.5, N);"
    } else if (xscale_set == "identity" & yscale_set == "log") {
      xscale_factor_str_d1 <- "rep_vector(1, N);"
      xscale_factor_str_d2 <- "rep_vector(1, N);"
      yscale_factor_str_d1 <- "(pred_d0);"
      yscale_factor_str_d2 <- "(pred_d0);"
    } else if (xscale_set == "sqrt" & yscale_set == "log") {
      xscale_factor_str_d1 <- "(Xm + xoffset);"
      xscale_factor_str_d2 <- "(Xm + xoffset);"
      yscale_factor_str_d1 <- "(rep_vector(0.5, N) .* (pred_d0));"
      yscale_factor_str_d2 <- "(rep_vector(0.5, N) .* (pred_d0));"
    } else if (xscale_set == "identity" & yscale_set == "sqrt") {
      xscale_factor_str_d1 <- "rep_vector(1, N);"
      xscale_factor_str_d2 <- "rep_vector(1, N);"
      yscale_factor_str_d1 <-
        "(rep_vector(2.0, N) .* sqrt(pred_d0));"
      yscale_factor_str_d2 <-
        "(rep_vector(2.0, N) .* sqrt(pred_d0));"
    } else if (xscale_set == "log" & yscale_set == "sqrt") {
      xscale_factor_str_d1 <- "exp(Xm + xoffset);"
      xscale_factor_str_d2 <- "exp(Xm + xoffset);"
      yscale_factor_str_d1 <-
        "(rep_vector(2.0, N) .* sqrt(pred_d0));"
      yscale_factor_str_d2 <-
        "(rep_vector(2.0, N) .* sqrt(pred_d0));"
    }
    
    list(
      xscale_factor_str_d1 = xscale_factor_str_d1,
      xscale_factor_str_d2 = xscale_factor_str_d2,
      yscale_factor_str_d1 = yscale_factor_str_d1,
      yscale_factor_str_d2 = yscale_factor_str_d2
    )
    
  } # end set_x_y_scale_factror
  
  
  setxoffset <- paste0("real xoffset = ", xoffset, ";")
  setxoffset_plane <- paste0("real xoffset = ", xoffset, ";")
  
  
  
  
  # https://mc-stan.org/users/documentation/case-studies/qr_regression.html
  
  decomp_code_qr <-
    "
      int QK = nknots - 1;
      matrix[N, QK] Qc = Spl;
      matrix[N, QK] XQ;
      matrix[QK, QK] XR;
      matrix[QK, QK] XR_inv;
      XQ = qr_thin_Q(Qc) * sqrt(N - 1);
      XR = qr_thin_R(Qc) / sqrt(N - 1);
      XR_inv = inverse(XR);
      "
  
  decomp_code_qr <-
    gsub("XR_inv", XR_inv_name, decomp_code_qr, fixed = T)
  
  
  ##########
  
  
  add_context_getx_fun <-
    "/* Transform x variable
 * Args:
 * Xp: x variable
 * Transformation code (tranform_x, 0 to 2)
 * 0, no transformation, 1 log, 2 square rooot
 * Note that the xoffset  is already transformed
 * Returns:
 * x variable with log/sqrt transformation
 */"
  
  add_context_getknots_fun <-
    "/* Knots
 * xoffset and Knots already transformed:
 * Returns:
 * Knots
 */"
  
  ##########
  
  create_internal_function <-
    function(y,
             function_str,
             fname,
             fnameout,
             spl,
             splout,
             xfunsi,
             yfunsi,
             setxoffset,
             gsub_out_unscaled,
             body,
             vectorA,
             decomp,
             fixedsi) {
      split1 <- strsplit(function_str, gsub("\\[", "\\\\[", spl))[[1]][-1]
      split2 <- strsplit(split1, "return")[[1]][-2]
      out <- gsub(split2, body, function_str, fixed = T)
      out <- gsub(spl, splout, out, fixed = T)
      out <- gsub(fname, fnameout, out, fixed = T)
      
      if (grepl("d0", fnameout)) {
        out <- out
      } else if (grepl("d1", fnameout) | grepl("d2", fnameout)) {
        out <- gsub("return(A+", "return(0+", out, fixed = T)
      }
      
      if (grepl("c", fixedsi, fixed = T)) {
        if (grepl("d0", fnameout)) {
          out <- gsub("]));", "]));", out, fixed = T)
          out <-
            gsub(
              "end of spline function",
              paste0("end of spline function", "_", y, "d", 0),
              out,
              fixed = T
            )
        } else if (grepl("d1", fnameout)) {
          out <- gsub("]));", "]).*exp(c));", out, fixed = T)
          out <-
            gsub(
              "end of spline function",
              paste0("end of spline function", "_", y, "d", 1),
              out,
              fixed = T
            )
        } else if (grepl("d2", fnameout)) {
          out <- gsub("]));", "]).*exp(c)^2);", out, fixed = T)
          out <-
            gsub(
              "end of spline function",
              paste0("end of spline function", "_", y, "d", 2),
              out,
              fixed = T
            )
        } else {
          out <- out
        }
      }
      
      ####
      if (grepl("d0", fnameout)) {
        pattern <- "return\\(\\s*(.*?)\\s*\\);"
        result <- regmatches(out, regexec(pattern, out))
        out_unscaled <-
          paste0("vector[N] out_unscaled=", result[[1]][2], ";")
        
        
        if (!is.null(gsub_out_unscaled)) {
          if (length(gsub_out_unscaled) != 2)
            stop('Length of gsub_out_unscaled should be 2')
          out_unscaled <-
            gsub(gsub_out_unscaled[1],
                 gsub_out_unscaled[2],
                 out_unscaled,
                 fixed = T)
        }
        
        
        if (yfunsi == "log") {
          out_scaled <-
            paste0("    vector[N] out_scaled=",
                   "exp",
                   "(",
                   "out_unscaled",
                   ")",
                   ";")
        } else if (yfunsi == "sqrt") {
          if ((backend == "rstan" &
               utils::packageVersion("rstan") >= "2.26.1") |
              backend == "mock" |
              backend == "cmdstanr") {
            out_scaled <-
              paste0("    vector[N] out_scaled=",
                     "",
                     "(",
                     "out_unscaled",
                     ")^2",
                     ";")
          }
          if ((backend == "rstan" &
               utils::packageVersion("rstan") < "2.26.1") | # &
              backend == "mock" &
              backend != "cmdstanr") {
            out_scaled <-
              paste0("    vector[N] out_scaled=",
                     "",
                     "(",
                     "pow(",
                     "out_unscaled",
                     ", 2)" ,
                     ")",
                     ";")
          }
        } else if (yfunsi != "log" & yfunsi != "sqrt") {
          out_scaled <-
            paste0("    vector[N] out_scaled=",
                   "",
                   "(",
                   "out_unscaled",
                   ")",
                   ";")
        }
        
        out <- gsub(result[[1]][2], "out_scaled", out, fixed = T)
        out_return <- paste0(out_unscaled, "\n", out_scaled)
        
        # setxoffset <- paste0("real xoffset = ", xoffset, ";")
        
        if (is.null(decomp))
          setxoffset <- paste0(setxoffset, vectorA)
        
        out_return <- paste0(setxoffset,
                             "\n    ",
                             out_return)
        out_return_p <- paste0(out_return, "\n", "    return")
        out <- gsub("return", out_return_p, out, fixed = T)
        
      } else if (grepl("d1", fnameout) | grepl("d2", fnameout)) {
        pattern <- "return\\(\\s*(.*?)\\s*\\);"
        result <- regmatches(out, regexec(pattern, out))
        
        set_x_y_scale <- set_x_y_scale_factror(
          xfunsi = xfunsi,
          yfunsi = yfunsi,
          tranformations =
            c("identity", "log", "sqrt")
        )
        
        if (grepl("d1", fnameout)) {
          xscale_factor <- set_x_y_scale[['xscale_factor_str_d1']]
          yscale_factor <- set_x_y_scale[['yscale_factor_str_d1']]
        } else if (grepl("d2", fnameout)) {
          xscale_factor <- set_x_y_scale[['xscale_factor_str_d2']]
          yscale_factor <- set_x_y_scale[['yscale_factor_str_d2']]
        }
        
        xscale_factor <- gsub(";", "", xscale_factor)
        yscale_factor <- gsub(";", "", yscale_factor)
        out_unscaled <-
          paste0("vector[N] out_unscaled=", result[[1]][2], ";")
        out_scaled <- paste0(
          "    vector[N] out_scaled=",
          "(",
          "(",
          yscale_factor,
          ")",
          " .* ",
          "(",
          'out_unscaled',
          ")",
          ")",
          " ./ ",
          "(",
          xscale_factor,
          ")",
          ";"
        )
        out <- gsub(result[[1]][2], "out_scaled", out, fixed = T)
        out_return <- paste0(out_unscaled, "\n", out_scaled)
        addpdo <- paste0("vector[N] pred_d0=", spl_fun_ford, ";")
        # setxoffset <- paste0("real xoffset = ", xoffset, ";")
        out_return <- paste0(addpdo,
                             "\n    ",
                             setxoffset,
                             "\n    ",
                             out_return)
        out_return_p <- paste0(out_return, "\n", "    return")
        
        if (!is.null(decomp)) {
          if (decomp == 'QR') {
            if (grepl("d1", fnameout)) {
              out_return_p <- gsub(
                vectorA,
                paste0(vectorA, "\n", "XQ[,1] = rep_vector(1, N);"),
                out_return_p,
                fixed = T
              )
              out_return_p <-
                gsub(vectorA, "", out_return_p, fixed = T)
            }
            if (grepl("d2", fnameout)) {
              out_return_p <- gsub(
                vectorA,
                paste0(vectorA, "\n", "XQ[,1] = rep_vector(0, N);"),
                out_return_p,
                fixed = T
              )
              out_return_p <-
                gsub(vectorA, "", out_return_p, fixed = T)
            }
          }
        }
        out <- gsub("return", out_return_p, out, fixed = T)
      }
      
      return(out)
    }
  
  
  
  
  
  ##########
  
  create_internal_function_nonsitar <-
    function(y,
             function_str,
             fname,
             fnameout,
             returnmu,
             xfunsi,
             yfunsi,
             setxoffset,
             gsub_out_unscaled = NULL,
             spl_fun_ford,
             body,
             decomp,
             fixedsi) {
      out <- function_str
      for_out <- gsub(fname, fnameout, out)
      
      ####
      if (grepl("d0", fnameout)) {
        out_unscaled <- paste0("vector[N] out_unscaled=", body, ";")
        
        if (!is.null(gsub_out_unscaled)) {
          if (length(gsub_out_unscaled) != 2)
            stop('Length of gsub_out_unscaled should be 2')
          out_unscaled <-
            gsub(gsub_out_unscaled[1],
                 gsub_out_unscaled[2],
                 out_unscaled,
                 fixed = T)
        }
        
        
        if (yfunsi == "log") {
          out_scaled <-
            paste0("    vector[N] out_scaled=",
                   "exp",
                   "(",
                   "out_unscaled",
                   ")",
                   ";")
        } else if (yfunsi == "sqrt") {
          if ((backend == "rstan" &
               utils::packageVersion("rstan") >= "2.26.1") |
              backend == "mock" |
              backend == "cmdstanr") {
            out_scaled <-
              paste0("    vector[N] out_scaled=",
                     "",
                     "(",
                     "out_unscaled",
                     ")^2.0",
                     ";")
          }
          if ((backend == "rstan" &
               utils::packageVersion("rstan") < "2.26.1") | # &
              backend == "mock" &
              backend != "cmdstanr") {
            out_scaled <-
              paste0("    vector[N] out_scaled=",
                     "",
                     "(",
                     "pow(",
                     "out_unscaled",
                     ", 2)" ,
                     ")",
                     ";")
          }
        } else if (yfunsi != "log" & yfunsi != "sqrt") {
          out_scaled <-
            paste0("    vector[N] out_scaled=",
                   "",
                   "(",
                   "out_unscaled",
                   ")",
                   ";")
        }
        
        out_return <- paste0(out_unscaled, "\n", out_scaled)
        # setxoffset <- paste0("real xoffset = ", xoffset, ";")
        out_return <- paste0(setxoffset,
                             "\n    ",
                             out_return)
        out_return_p <- paste0(out_return, "\n", "    return")
        out_scaled_with_parentehsis <-
          paste0("(", 'out_scaled', ")")
        out <- paste(out_return_p, out_scaled_with_parentehsis, ";")
        out <- paste0(gsub("return.*", "", for_out),
                      out,
                      "\n}")
        
      } else if (grepl("d1", fnameout) | grepl("d2", fnameout)) {
        set_x_y_scale <- set_x_y_scale_factror(
          xfunsi = xfunsi,
          yfunsi = yfunsi,
          tranformations =
            c("identity", "log", "sqrt")
        )
        
        if (grepl("d1", fnameout)) {
          xscale_factor <- set_x_y_scale[['xscale_factor_str_d1']]
          yscale_factor <- set_x_y_scale[['yscale_factor_str_d1']]
        } else if (grepl("d2", fnameout)) {
          xscale_factor <- set_x_y_scale[['xscale_factor_str_d2']]
          yscale_factor <- set_x_y_scale[['yscale_factor_str_d2']]
        }
        
        xscale_factor <- gsub(";", "", xscale_factor)
        yscale_factor <- gsub(";", "", yscale_factor)
        out_unscaled <- paste0("vector[N] out_unscaled=", body, ";")
        out_scaled <-
          paste0(
            "    vector[N] out_scaled=",
            "(",
            "(",
            yscale_factor,
            ")",
            " .* ",
            "(",
            'out_unscaled',
            ")",
            ")",
            " ./ ",
            "(",
            xscale_factor,
            ")",
            ";"
          )
        addpdo <- paste0("vector[N] pred_d0=", spl_fun_ford, ";")
        out_return <- paste0(out_unscaled, "\n", out_scaled)
        # setxoffset <- paste0("real xoffset = ", xoffset, ";")
        out_return <- paste0(addpdo,
                             "\n    ",
                             setxoffset,
                             "\n    ",
                             out_return)
        out_return_p <- paste0(out_return, "\n", "    return")
        out_scaled_with_parentehsis <-
          paste0("(", 'out_scaled', ")")
        out <- paste(out_return_p, out_scaled_with_parentehsis, ";")
        out <- paste0(gsub("return.*", "", for_out),
                      out,
                      "\n}")
      }
      return(out)
    }
  
  ##########
  
  
  
  if (select_model == 'sitar' | select_model == 'rcs') {
    abcnames <-
      paste0(strsplit(gsub("\\+", " ", fixedsi), " ")[[1]], sep = ",")
    
    snames <- c()
    for (i in 1:(nknots - 1)) {
      if (i < (nknots - 1)) {
        name1 <- paste0("s", i, sep = ",")
      }
      else {
        name1 <- paste0("s", i, sep = "")
      }
      snames[i] <- name1
    }
    
    # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
    # In fact for df > 1, it forces 'd' to be random parameter only
    if (match_sitar_d_form) {
      if (!grepl("d", fixedsi, fixed = T) &
          grepl("d", randomsi, fixed = T)) {
        abcnames <- c(abcnames, "d,")
      }
    }
    
    
    if (select_model == 'sitar' | select_model == 'rcs') {
      if (any(grepl("s", abcnames)))
        abcnames <- abcnames[-length(abcnames)]
      if (match_sitar_d_form)
        abcnames <- gsub('s', 'd', abcnames, fixed = T)
    }
    
    fullabcsnames <- c(abcnames, snames)
    fullabcsnames_v <-
      paste("vector", fullabcsnames, collapse = " ")
    
    if (select_model == 'sitar') {
      fullabcsnames_for_mat <- abcnames
      fullabcsnames_for_mat <-
        gsub('.{1}$', '', fullabcsnames_for_mat)
      fullabcsnames_v_for_mat <-
        paste("vector", fullabcsnames_for_mat, collapse = ", ")
    }
    
    if (select_model == 'rcs') {
      fullabcsnames_for_mat <- 'a'
      fullabcsnames_v_for_mat <-
        paste("vector", fullabcsnames_for_mat, collapse = " ")
    }

    if (grepl("b", fixedsi, fixed = T) &
        grepl("c", fixedsi, fixed = T)) {
      defineEx <- paste0("(Xm-b).*exp(c)")
    }
    if (grepl("b", fixedsi, fixed = T) &
        !grepl("c", fixedsi, fixed = T)) {
      defineEx <- paste0("(Xm-b)")
    }
    if (!grepl("b", fixedsi, fixed = T) &
        grepl("c", fixedsi, fixed = T)) {
      defineEx <- paste0("(Xm).*exp(c)")
    }
    if (!grepl("b", fixedsi, fixed = T) &
        !grepl("c", fixedsi, fixed = T)) {
      defineEx <- paste0("(Xm)")
    }
    
    add_knotinfo <- paste0(
      "\n  int N=num_elements(",
      vector_X_name,
      ");",
      paste0(
        "\n  vector[N] Xm=",
        paste0(getxname,
               "(", vector_X_name, ")"),
        ";"
      ),
      paste0("\n  vector[N] X=", defineEx, ";"),
      paste0("\n  int nknots=", eval(parse(text = nknots)), ";"),
      paste0(
        "\n  vector[nknots] knots=",
        paste0(getknotsname, "(", '', ")"),
        ";"
      )
    )
    
    
    if (select_model == 'sitar') {
      vectorA <- "\n  vector[N] A=a-(s1*min(knots));"
      if (!is.null(decomp)) {
        if (decomp == 'QR') {
          vectorA <- "\n  vector[N] A=a;"
          # vectorA <- "\n  vector[N] A=a-((XQ[,1].*s1)*min(knots));"
        }
      }
    }
    
    
    
    if (select_model == 'rcs') {
      vectorA <- "\n  vector[N] A=a;"
    }
    
    
    
    # add_knotinfo <- paste0(add_knotinfo, vectorA)
    
    if ((backend == "rstan" &
         utils::packageVersion("rstan") >= "2.26.1") |
        backend == "mock" |
        backend == "cmdstanr") {
      fun_body <- "
    matrix[N, nknots-1] Spl;
    matrix[nknots-1, N] rcs;
    matrix[N, nknots] Xx;
    int km1 = nknots - 1;
    int jp1;
    int j=1;
    for(ia in 1:N) {
     for(ja in 1:nknots) {
      Xx[ia,ja] = (X[ia] - knots[ja] > 0 ? X[ia] - knots[ja] : 0);
        }
    }
     Spl[,1]=X;
     while (j <= nknots - 2) {
      jp1 = j + 1;
       Spl[,jp1] = (Xx[,j]^3-(Xx[,km1]^3)*(knots[nknots]-knots[j])/
       (knots[nknots]-knots[km1]) + (Xx[,nknots]^3)*(knots[km1]-knots[j])/
       (knots[nknots]-knots[km1])) / (knots[nknots]-knots[1])^2;
        j = j + 1;
      }"
    }
    
    if ((backend == "rstan" &
         utils::packageVersion("rstan") < "2.26.1") | # &
        backend == "mock" &
        backend != "cmdstanr") {
      fun_body <- "
    matrix[N, nknots-1] Spl;
    matrix[nknots-1, N] rcs;
    matrix[N, nknots] Xx;
    int km1 = nknots - 1;
    int jp1;
    int j=1;
    for(ia in 1:N) {
     for(ja in 1:nknots) {
      Xx[ia,ja] = (X[ia] - knots[ja] > 0 ? X[ia] - knots[ja] : 0);
        }
    }
     Spl[,1]=X;
     while (j <= nknots - 2) {
     for(i in 1:N) {
      jp1 = j + 1;
       Spl[i,jp1] = (pow(Xx[i,j],3)-(pow(Xx[i,km1],3))*(knots[nknots]-knots[j])/
       (knots[nknots]-knots[km1]) +
       (pow(Xx[i,nknots],3))*(knots[km1]-knots[j])/(knots[nknots]-knots[km1]))/
       pow((knots[nknots]-knots[1]),2);
      }
      j = j + 1;
    }
    "
    }
    
    name4 <- c()
    for (i in 1:(nknots - 1)) {
      name1 <- paste0("", "s", i, sep = "")
      if (i < (nknots - 1)) {
        # name2 <- paste0(' .* to_vector(Spl[,',i,"]') +")
        name2 <- paste0(' .* Spl[,', i, "] +")
      }
      else {
        # name2 <- paste0(' .* to_vector(Spl[,',i,"]') ;\n")
        name2 <- paste0(' .* Spl[,', i, "]")
      }
      name3 <- paste0(name1, name2, sep = "")
      name4[i] <- name3
    }
    name50 <- paste("", name4, collapse = " ")
    
    nameadja <- "A"
    
    ###########
    # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
    # In fact for df > 1, it forces 'd' to be random parameter only
    
    if (match_sitar_d_form) {
      if (grepl("d", randomsi, fixed = T)) {
        if( ept(d_adjustedsi)) nameadja <- "A+(d . * Spl[,1])"
        if(!ept(d_adjustedsi)) nameadja <- "A+(d . * Xm)"
      }
    }
    
    if (!match_sitar_d_form) {
      if (grepl("d", fixedsi, fixed = T)) {
        if( ept(d_adjustedsi)) nameadja <- "A+(d . * Spl[,1])"
        if(!ept(d_adjustedsi)) nameadja <- "A+(d . * Xm)"
      }
    }
    
    
    name5 <- paste(" (", name50, ");\n")
    
    if (grepl("c", fixedsi, fixed = T)) {
      name51 <- paste(" (", name50, ") .* exp(c) ;\n")
      name52 <- paste(" (", name50, ") .* exp(c)^2 ;\n")
    }
    if (!grepl("c", fixedsi, fixed = T)) {
      name51 <- paste(" (", name50, ") ;\n")
      name52 <- paste(" (", name50, ") ;\n")
    }
    
    returnmu <-
      paste0("return(",   paste0(nameadja, "+",
                                 gsub(";", "", name5))     , ");")
    
    # need spaces otherwise rstan 2.21 throws error: variable s1. not found
    returnmu <- gsub("\\s", "", returnmu)
    returnmu <- gsub("\\." , " \\." , returnmu, fixed = FALSE)
    returnmu <- gsub("\\*" , "\\* " , returnmu, fixed = FALSE)
    # don't create space for +
    # returnmu <- gsub("+" , " + " , returnmu, fixed = TRUE)
    
    
    setxoffset_d0_noqr <- paste0(setxoffset,  vectorA)
    returnmu_d0_noqr <- paste0(setxoffset,  vectorA)
    
    if (!is.null(decomp)) {
      if (decomp == 'QR') {
        returnmu <- gsub('Spl', 'XQ', returnmu, fixed = T)
        setxoffset <- paste0(setxoffset,  decomp_code_qr, vectorA)
        returnmu <- gsub('Spl', 'XQ', returnmu, fixed = T)
      }
    }
    
    if (is.null(decomp)) {
      fun_body <- paste0(fun_body, "\n", vectorA)
    }
    
    
    
    endof_fun <-
      paste0("\n    ", returnmu, "\n  } // end of spline function", sep = " ")
    
    
    
    
    start_fun <-
      paste0(
        "\nvector ",
        spfncname,
        "(vector ",
        vector_X_name,
        ", ",
        fullabcsnames_v,
        ") {" ,
        collapse = " "
      )
    
    
    
    
    
    rcsfun <-
      paste(start_fun,
            add_knotinfo,
            fun_body,
            "\n",
            setxoffset,
            endof_fun)
    rcsfun_raw <- rcsfun
    
    
    
    rcsfunmultadd <- NULL
    spfncname_multadd <- paste0(spfncname, "X")
    
    start_fun_multadd <-
      paste0(
        "\nvector ",
        spfncname_multadd,
        "(matrix ",
        vector_X_name,
        ", ",
        fullabcsnames_v,
        ") {" ,
        "\n",
        "  int N=num_elements(",
        vector_X_name,
        "[,1]);",
        "\n",
        "  vector[N] Xm=(",
        vector_X_name,
        "[,1]);",
        # getX
        # "\n",
        # insert_getX_name,
        collapse = " "
      )
    
    
    add_knotinfo_multadd <- paste0(
      "\n  int mcolsmat=cols(",
      vector_X_name,
      ");",
      paste0(
        "\n  vector[mcolsmat+1] knots=",
        paste0(getknotsname, "(", '', ")"),
        ";"
      )
    )
    returnmu_multadd <- returnmu
    returnmu_multadd <-
      gsub("XQ",  vector_X_name, returnmu_multadd, fixed = T)
    returnmu_multadd <-
      gsub("Spl", vector_X_name, returnmu_multadd, fixed = T)
    
    start_fun_multadd <-
      gsub(",)" , ")" , start_fun_multadd, fixed = TRUE)
    endof_fun <- paste0("\n    ",
                        returnmu_multadd,
                        ";",
                        "\n  } // end of spline mat function",
                        sep = " ")
    
    endof_fun <- gsub(";;", ";", endof_fun, fixed = T)
    
    if (grepl('s1', vectorA)) {
      rcsfunmultadd <-
        paste(start_fun_multadd,
              add_knotinfo_multadd,
              vectorA,
              endof_fun)
    } else {
      rcsfunmultadd <- paste(start_fun_multadd, vectorA, endof_fun)
    }
    
    
    
    add_rcsfunmat <- TRUE
    add_rcsfunmatqr <- TRUE
    add_rcsfunmatqrinv <- TRUE
    
    
    funmats <- paste0('', '')
    
    if (add_rcsfunmat) {
      rcsfunmat_name <- paste0(spfncname, 'smat', '')
      start_funmat <-
        paste0(
          "\nmatrix ",
          rcsfunmat_name,
          "(vector ",
          vector_X_name,
          ", ",
          fullabcsnames_v_for_mat,
          ") {" ,
          collapse = " "
        )
      setxoffset_format <- paste0(setxoffset_plane, vectorA)
      rcsfunmat <-
        paste(start_funmat,
              add_knotinfo,
              fun_body,
              "\n",
              setxoffset_format)
      rcsfunmat <- gsub(vectorA, "", rcsfunmat, fixed = T)
      rcsfunmat <- paste0(rcsfunmat, "\n", 'return Spl;')
      rcsfunmat <- paste0(rcsfunmat, '\n}')
      funmats <- paste0(funmats, "\n", rcsfunmat)
    }
    
    if (add_rcsfunmatqr) {
      rcsfunmatgr_name <- paste0(spfncname, 'QRsmat', '')
      
      start_funmat <-
        paste0(
          "\nmatrix ",
          rcsfunmatgr_name,
          "(vector ",
          vector_X_name,
          ", ",
          fullabcsnames_v_for_mat,
          ") {" ,
          collapse = " "
        )
      
      rcsfunmatqr <-
        paste(start_funmat, add_knotinfo, fun_body, "\n", setxoffset)
      rcsfunmatqr <- gsub(vectorA, "", rcsfunmatqr, fixed = T)
      
      
      
      szx <- paste0('s', 1:(nknots - 1))
      cnt <- 0
      cn_c <- c()
      for (vi in szx) {
        cnt <- cnt + 1
        tmx <-
          paste0(
            paste0('s', 'x', '[,', cnt, "]", " = "),
            XR_inv_name_resp,
            '[',
            cnt,
            ",",
            cnt,
            "] * ",
            vi,
            ";"
          )
        cn_c <- c(cn_c, tmx)
      }
      rcsfunmatqr <- paste0(rcsfunmatqr, "\n", 'return XQ;')
      rcsfunmatqr <- paste0(rcsfunmatqr, '\n}')
      funmats <- paste0(funmats, "\n", rcsfunmatqr)
    } # if(add_rcsfunmatqr) {
    
    
    
    if (add_rcsfunmatqrinv) {
      rcsfunmatgrinv_name <- paste0(spfncname, 'QRsmat', 'inv')
      
      start_funmat <-
        paste0(
          "\nmatrix ",
          rcsfunmatgrinv_name,
          "(vector ",
          vector_X_name,
          ", ",
          fullabcsnames_v_for_mat,
          ") {" ,
          collapse = " "
        )
      
      rcsfunmatgrinv <-
        paste(start_funmat, add_knotinfo, fun_body, "\n", setxoffset)
      rcsfunmatgrinv <- gsub(vectorA, "", rcsfunmatgrinv, fixed = T)
      
      szx <- paste0('b_s', 1:(nknots - 1))
      cnt <- 0
      cn_c <- c()
      for (vi in szx) {
        cnt <- cnt + 1
        tmx <-
          paste0(
            paste0('b_s', 'q', '[,', cnt, "]", " = "),
            XR_inv_name_resp,
            '_mat',
            '[',
            cnt,
            ",",
            cnt,
            "] * ",
            vi,
            ";"
          )
        cn_c <- c(cn_c, tmx)
      }
      rcsfunmatgrinv <-
        paste0(rcsfunmatgrinv, "\n", 'return ', XR_inv_name, ';')
      rcsfunmatgrinv <- paste0(rcsfunmatgrinv, '\n}')
      funmats <- paste0(funmats, "\n", rcsfunmatgrinv)
    } # if(add_rcsfunmatqrinv) {
    
    
    
    
    funmats_genquant <- ""
    if (add_rcsfunmatqrinv_genquant) {
      rcsfunmatgrinv_name <- paste0(spfncname, 'QRsmat', 'inv')
      
      setnqq <- nknots - 1
      start_funmat <-
        paste0(
          "\nmatrix ",
          rcsfunmatgrinv_name,
          "(vector ",
          vector_X_name,
          ", ",
          fullabcsnames_v_for_mat,
          ") {" ,
          collapse = " "
        )
      
      rcsfunmatqrinv_genquant <- paste(start_funmat, '')
      rcsfunmatqrinv_genquant <-
        gsub(vectorA, "", rcsfunmatqrinv_genquant, fixed = T)
      
      szx <- paste0(setnlp, 1:(nknots - 1))
      szxbq <- b_sx_name_resp
      
      szx_vector <- paste0(setnlp_vector, 1:(nknots - 1))
      
      
      # setnqq <- 4
      tems <- paste0('placeholder', '_', b_sx_name_resp)
      # b_s_name <- paste0('b_s', '')
      b_s_dims1 <- '[1]'
      svector_ <-
        paste(paste0(b_s_name,  b_s_dims1), collapse = ", ")
      svector_ <-
        paste0('vector[',
               setnqq,
               ']',
               " ",
               tems,
               ' = ',
               '[ ',
               svector_,
               ' ]' ,
               "'" ,
               ";")
      
      tems2 <-  paste0('temp_', b_sx_name_resp)
      
      tems3 <- paste0('vector[', setnqq, ']', " ", tems2, ' = ')
      qs_vector <-
        paste0(tems3, XR_inv_name_resp , ' * ', tems, ";")
      
      cn_c <- c()
      cn_c2 <- c()
      cnt <- 0
      for (vi in szx) {
        cnt <- cnt + 1
        tmx <-
          paste0(
            'vector[N] ',
            paste0(szxbq_vector, cnt, " = "),
            XR_inv_name_resp,
            '[',
            cnt,
            ",",
            cnt,
            "] * ",
            vi,
            ";"
          )
        cn_c <- c(cn_c, tmx)
      }
      cnt <- 0
      for (vi in szx_vector) {
        cnt <- cnt + 1
        # '[', cnt, ',', cnt, "]"
        tmx2 <-
          paste0(paste0(szxbq, '[', cnt, ']', '', " = "),
                 tems2,
                 '[',
                 cnt,
                 "]",
                 ";")
        
        cn_c2 <- c(cn_c2, tmx2)
      }
      
      # if no covariate. then only add re-scaled s betas
      if (add_b_Qr_genquan_s_coef) {
        cn_c <- paste(cn_c, collapse = "\n")
        cn_c2 <- paste(cn_c2, collapse = "\n")
        cn_c <- paste0(cn_c2, "\n", cn_c)
        addcn_c2 <-
          paste0('vector[', "", setnqq, ']', " ", szxbq, ';')
        addcn_c2 <- paste0(addcn_c2, "\n", svector_)
        addcn_c2 <- paste0(addcn_c2, "\n", qs_vector)
        addcn_c2 <- paste0(addcn_c2, "\n", cn_c)
      } else if (!add_b_Qr_genquan_s_coef) {
        cn_c <- paste(cn_c, collapse = "\n")
        addcn_c2 <- cn_c
      }
      replacematrixby <-
        paste0('matrix[', setnqq, ", ", setnqq, ']', XR_inv_name_resp , ' =')
      rcsfunmatqrinv_genquant <-
        gsub('matrix',
             replacematrixby,
             rcsfunmatqrinv_genquant,
             fixed = T)
      rcsfunmatqrinv_genquant <-
        gsub('{', ';', rcsfunmatqrinv_genquant, fixed = T)
      rcsfunmatqrinv_genquant <-
        paste0(rcsfunmatqrinv_genquant, '\n', addcn_c2)
    } # if(add_rcsfunmatqrinv_genquant) {
    
   
    if (funmats == "") {
      add_funmats <- FALSE
    } else {
      add_funmats <- TRUE
    }
    
   
    
    
    getknots_fun_raw <-
      paste0(
        "vector ",
        getknotsname,
        "() {" ,
        paste0("\n  int nknots=", eval(parse(text = nknots)), ";"),
        paste0(
          "\n  vector[nknots] knots=",
          "[",
          paste(knots, collapse = ","),
          "]';"
        ),
        "\n  ",
        "return(knots);",
        "\n}  "
        ,
        paste0("// end of ", getknotsname),
        collapse = " "
      )
    
    # add_context_getx_fun
    getx_fun_raw <-
      paste0(
        "vector ",
        getxname,
        paste0(" (vector ", vector_X_name, ") {") ,
        "\n  ",
        paste0("int N=num_elements(", vector_X_name, ");"),
        "\n  ",
        paste0("real xoffset = ", xoffset, ";"),
        paste0("\n  int tranform_x = ",
               eval(parse(text = tranform_x_int)), ";"),
        paste0("\n  vector[N] x;"),
        "\n",
        paste0(
          "  if(tranform_x == 0 ) {",
          "\n   ",
          "x = ",
          vector_X_name,
          " - xoffset",
          ";",
          "\n  }",
          "\n  ",
          "if(tranform_x == 1 ) {",
          "\n   ",
          "x = log(",
          vector_X_name,
          ") - xoffset;",
          "\n  }",
          "\n  ",
          "if(tranform_x == 2 ) {",
          "\n    ",
          "x = sqrt(",
          vector_X_name,
          ") - xoffset;",
          "\n  }"
        ),
        "\n  ",
        "return(x);",
        "\n}  "
        ,
        paste0("// end of ", getxname),
        collapse = " "
      )
    
    
    
    
    getx_fun     <- paste0(add_context_getx_fun, "\n", getx_fun_raw)
    getknots_fun <-
      paste0(add_context_getknots_fun, "\n", getknots_fun_raw)
    
    getx_knots_fun <- paste0(getx_fun,
                             "\n",
                             getknots_fun)
    
    
    ##########
    
    # Create function d0
    fnameout <- paste0(spfncname, "_", "d0")
    spl <- "Spl[,1]=X;"
    splout <- spl
    spl_fun_ford <- paste0(fnameout,
                           "(vector ",
                           vector_X_name,
                           ", ",
                           fullabcsnames_v,
                           ")")
    spl_fun_ford <- gsub("vector", "", spl_fun_ford, fixed = T)
    
    if ((backend == "rstan" &
         utils::packageVersion("rstan") >= "2.26.1") |
        backend == "mock" |
        backend == "cmdstanr") {
      body <- "
     while (j <= nknots - 2) {
      jp1 = j + 1;
      Spl[,jp1] =
        (1*Xx[,j]^3) * (1/((knots[nknots]-knots[1])^2))  -
        (1*Xx[,km1]^3) * (knots[nknots]-knots[j]) /
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) +
        (1*Xx[,nknots]^3) * (knots[km1]-knots[j]) /
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) ;
      j = j + 1;
    }
    "
    }
    
    if ((backend == "rstan" &
         utils::packageVersion("rstan") < "2.26.1") | # &
        backend == "mock" &
        backend != "cmdstanr") {
      body <- "
     while (j <= nknots - 2) {
      for(i in 1:N) {
          jp1 = j + 1;
          Spl[i,jp1] = (1*pow(Xx[i,j],3) -
          (1*pow(Xx[i,km1],3))*(knots[nknots]-knots[j]) /
          (knots[nknots]-knots[km1]) +
          (1*pow(Xx[i,nknots],3))*(knots[km1]-knots[j]) /
          (knots[nknots]-knots[km1])) /
          pow((knots[nknots]-knots[1]),2);
      }
      j = j + 1;
    }
    "
    }
    
    
    spl_d0 <- create_internal_function(
      y = y,
      function_str = rcsfun,
      fname = spfncname,
      fnameout = fnameout,
      spl = spl,
      splout = splout,
      xfunsi = xfunsi,
      yfunsi = yfunsi,
      setxoffset = setxoffset,
      # setxoffset setxoffset_d0_noqr
      gsub_out_unscaled = NULL,
      # gsub_out_unscaled = c('QR', 'Spl')
      body = body,
      vectorA = vectorA,
      decomp = decomp,
      fixedsi = fixedsi
    )
    
    
    
    # Create function d1
    fnameout <- paste0(spfncname, "_", "d1")
    spl <- "Spl[,1]=X;"
    splout <- "Spl[,1]=rep_vector(1, N);"
    
    if ((backend == "rstan" &
         utils::packageVersion("rstan") >= "2.26.1") |
        backend == "mock" |
        backend == "cmdstanr") {
      body <- "
     while (j <= nknots - 2) {
      jp1 = j + 1;
      Spl[,jp1] =
        (3*Xx[,j]^2) * (1/((knots[nknots]-knots[1])^2))  -
        (3*Xx[,km1]^2) * (knots[nknots]-knots[j]) /
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) +
        (3*Xx[,nknots]^2) * (knots[km1]-knots[j]) /
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) ;
      j = j + 1;
    }
    "
    }
    
    if ((backend == "rstan" &
         utils::packageVersion("rstan") < "2.26.1") | # &
        backend == "mock" &
        backend != "cmdstanr") {
      body <- "
     while (j <= nknots - 2) {
      for(i in 1:N) {
          jp1 = j + 1;
          Spl[i,jp1] = (3*pow(Xx[i,j],2) -
          (3*pow(Xx[i,km1],2))*(knots[nknots]-knots[j]) /
          (knots[nknots]-knots[km1]) +
          (3*pow(Xx[i,nknots],2))*(knots[km1]-knots[j]) /
          (knots[nknots]-knots[km1])) /
          pow((knots[nknots]-knots[1]),2);
      }
      j = j + 1;
    }
    "
    }
    
    
    spl_d1 <- create_internal_function(
      y = y,
      function_str = rcsfun,
      fname = spfncname,
      fnameout = fnameout,
      spl = spl,
      splout = splout,
      xfunsi = xfunsi,
      yfunsi = yfunsi,
      setxoffset = setxoffset,
      # setxoffset setxoffset_d0_noqr
      gsub_out_unscaled = NULL,
      # gsub_out_unscaled = c('QR', 'Spl')
      body = body,
      vectorA = vectorA,
      decomp = decomp,
      fixedsi = fixedsi
    )
    
    
    # Create function d2
    fnameout <- paste0(spfncname, "_", "d2")
    spl <- "Spl[,1]=X;"
    splout <- "Spl[,1]=rep_vector(0, N);"
    
    if ((backend == "rstan" &
         utils::packageVersion("rstan") >= "2.26.1") |
        backend == "mock" |
        backend == "cmdstanr") {
      body <- "
     while (j <= nknots - 2) {
      jp1 = j + 1;
      Spl[,jp1] =
        (6*Xx[,j]^1) * (1/((knots[nknots]-knots[1])^2))  -
        (6*Xx[,km1]^1) * (knots[nknots]-knots[j]) /
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) +
        (6*Xx[,nknots]^1) * (knots[km1]-knots[j]) /
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) ;
      j = j + 1;
   }
    "
    }
    
    if ((backend == "rstan" &
         utils::packageVersion("rstan") < "2.26.1") | # &
        backend == "mock" &
        backend != "cmdstanr") {
      body <- "
     while (j <= nknots - 2) {
      for(i in 1:N) {
          jp1 = j + 1;
          Spl[i,jp1] = (6*pow(Xx[i,j],1) -
          (6*pow(Xx[i,km1],1))*(knots[nknots]-knots[j]) /
          (knots[nknots]-knots[km1]) +
          (6*pow(Xx[i,nknots],1))*(knots[km1]-knots[j]) /
          (knots[nknots]-knots[km1])) /
          pow((knots[nknots]-knots[1]),2);
      }
      j = j + 1;
    }
    "
    }
    
    
    
    spl_d2 <- create_internal_function(
      y = y,
      function_str = rcsfun,
      fname = spfncname,
      fnameout = fnameout,
      spl = spl,
      splout = splout,
      xfunsi = xfunsi,
      yfunsi = yfunsi,
      setxoffset = setxoffset,
      # setxoffset setxoffset_d0_noqr
      gsub_out_unscaled = NULL,
      # gsub_out_unscaled = c('QR', 'Spl')
      body = body,
      vectorA = vectorA,
      decomp = decomp,
      fixedsi = fixedsi
    )
    
    
    if (utils::packageVersion('rstan') < 2.26) {
      rcsfun <- paste(getx_knots_fun, rcsfun)
    }
    
    if (utils::packageVersion('rstan') > 2.26 & is.null(decomp)) {
      rcsfun <- paste0(getx_knots_fun,
                       rcsfun,
                       rcsfunmultadd,
                       spl_d0,
                       spl_d1,
                       spl_d2,
                       sep = "\n")
    }
    
    if (utils::packageVersion('rstan') > 2.26 & !is.null(decomp)) {
      if (decomp == 'QR') {
        if (add_funmats) {
          rcsfun <- paste0(
            getx_knots_fun,
            funmats,
            rcsfun,
            rcsfunmultadd,
            spl_d0,
            spl_d1,
            spl_d2,
            sep = "\n"
          )
        } else if (!add_funmats) {
          rcsfun <- paste0(getx_knots_fun,
                           rcsfun,
                           rcsfunmultadd,
                           spl_d0,
                           spl_d1,
                           spl_d2,
                           sep = "\n")
        }
      }
    }
    
    
  } # if(select_model == 'sitar') {
  
  
  
  
  
  if (grepl("^pb", select_model) |
      grepl("^logistic", select_model)) {
    abcnames <- paste0(strsplit(gsub("\\+", " ",
                                     fixedsi), " ")[[1]], sep = ",")
    fullabcsnames <- abcnames
    fullabcsnames_v <-
      paste("vector", fullabcsnames, collapse = " ")
    defineEx <- paste0("(Xm)")
    getx_fun_raw <-
      paste0(
        "vector ",
        getxname,
        paste0(" (vector ", vector_X_name, ") {") ,
        "\n  ",
        paste0("int N=num_elements(", vector_X_name, ");"),
        "\n  ",
        setxoffset,
        paste0("\n  int tranform_x = ",
               eval(parse(text = tranform_x_int)), ";"),
        paste0("\n  vector[N] x;"),
        "\n",
        paste0(
          "  if(tranform_x == 0 ) {",
          "\n   ",
          "x = ",
          vector_X_name,
          " - xoffset",
          ";",
          "\n  }",
          "\n  ",
          "if(tranform_x == 1 ) {",
          "\n   ",
          "x = log(",
          vector_X_name,
          ") - xoffset;",
          "\n  }",
          "\n  ",
          "if(tranform_x == 2 ) {",
          "\n    ",
          "x = sqrt(",
          vector_X_name,
          ") - xoffset;",
          "\n  }"
        ),
        "\n  ",
        "return(x);",
        "\n}  "
        ,
        paste0("// end of ", getxname),
        collapse = " "
      )
    
    
    getx_fun     <- paste0(add_context_getx_fun, "\n", getx_fun_raw)
    
    # a - asymtote
    # b - size at theta
    # c - s0
    # d - s1
    # e - time (theta)
    
    if (select_model == 'pb1') {
      funstring <- "a-2.0*(a-b)./(exp(c.*(Xm-e))+exp(d.*(Xm-e)))"
      if (utils::packageVersion('rstan') < 2.26)
        funstring <-
          gsub(".*", " .* ", funstring, fixed = T)
      returnmu    <- paste0("return ", "(",  funstring, ")")
      returnmu_d0 <- funstring
      returnmu_d1 <- "rep_vector(2.0,N).*(a-b).*(c.*exp(c.*(Xm-e))+
      d.*exp(d.*(Xm-e)))./(exp(c.*(Xm-e))+exp(d.*(Xm-e)))^2.0"
      
      returnmu_d2 <- "-rep_vector(4.0,N).*(a-b).*(c.*exp(c.*(Xm-e))+
      d.*exp(d.*(Xm-e)))^2.0./(exp(c.*(Xm-e))+exp(d.*(Xm-e)))^3.0+
      rep_vector(2.0,N).*(a-b).*(c^2.0.*exp(c.*(Xm-e))+
      d^2.0.*exp(d.*(Xm-e)))./(exp(c.*(Xm-e))+exp(d.*(Xm-e)))^2.0"
      
      returnmu_d3 <- "rep_vector(12.0,N).*(a-b).*(c.*exp(c.*(Xm-e))+
      d.*exp(d.*(Xm-e)))^3.0./(exp(c.*(Xm-e))+
      exp(d.*(Xm-e)))^4.0-rep_vector(12.0,N).*(a-b).*(c.*exp(c.*(Xm-e))+
      d.*exp(d.*(Xm-e))).*(c^2.0.*exp(c.*(Xm-e))+
      d^2.0.*exp(d.*(Xm-e)))./(exp(c.*(Xm-e))+
      exp(d.*(Xm-e)))^3.0+rep_vector(2.0,N).*(a-b).*(c^3.0.*exp(c.*(Xm-e))+
      d^3.0.*exp(d.*(Xm-e)))./(exp(c.*(Xm-e))+exp(d.*(Xm-e)))^2.0"
      
      returnmu_d1 <- paste0(returnmu_d1, ";")
      returnmu_d2 <- paste0(returnmu_d2, ";")
      returnmu_d3 <- paste0(returnmu_d3, ";")
    } # if(select_model == 'pb') {
    
    # a - asymtote
    # b - size at theta
    # c - s0
    # d - s1
    # e - time (theta)
    # f - gamma
    
    if (select_model == 'pb2') {
      funstring <- "a-((a-b)./(((0.5*exp((f.*c).*(Xm-e)))+
      (0.5*exp((f.*d).*(Xm-e))))^(1.0./f)))"
      if (utils::packageVersion('rstan') < 2.26)
        funstring <-
          gsub(".*", " .* ", funstring, fixed = T)
      returnmu    <- paste0("return ", "(",  funstring, ")")
      returnmu_d0 <- funstring
      returnmu_d1 <-
        "(a-b).*(rep_vector(0.5,N).*f.*c.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f.*d.*exp(f.*d.*(Xm-e)))./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e))))"
      
      returnmu_d2 <-
        "-(a-b).*(rep_vector(0.5,N).*f.*c.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f.*d.*exp(f.*d.*(Xm-e)))^2.0./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f^2.0.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^2.0)+
      (a-b).*(rep_vector(0.5,N).*f^2.0.*c^2.0.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f^2.0.*d^2.0.*exp(f.*d.*(Xm-e)))./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e))))-
      (a-b).*(rep_vector(0.5,N).*f.*c.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f.*d.*exp(f.*d.*(Xm-e)))^2.0./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^2.0)"
      
      returnmu_d3 <-
        "(a-b).*(rep_vector(0.5,N).*f.*c.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f.*d.*exp(f.*d.*(Xm-e)))^
      3.0./((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f^3.0.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^3.0)-
      rep_vector(3.0,N).*(a-b).*(rep_vector(0.5,N).*f.*c.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f.*d.*exp(f.*d.*(Xm-e))).*(rep_vector(0.5,N).*f^
      2.0.*c^2.0.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f^2.0.*d^2.0.*exp(f.*d.*(Xm-e)))./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f^2.0.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^2.0)+
      rep_vector(3.0,N).*(a-b).*(rep_vector(0.5,N).*f.*c.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f.*d.*exp(f.*d.*(Xm-e)))^3.0./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^(rep_vector(1.0,N)./f).*f^
      2.0.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^3.0)+
      (a-b).*(rep_vector(0.5,N).*f^3.0.*c^
      3.0.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f^3.0.*d^3.0.*exp(f.*d.*(Xm-e)))./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e))))-
      rep_vector(3.0,N).*(a-b).*(rep_vector(0.5,N).*f^
      2.0.*c^2.0.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f^2.0.*d^
      2.0.*exp(f.*d.*(Xm-e))).*(rep_vector(0.5,N).*f.*c.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f.*d.*exp(f.*d.*(Xm-e)))./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^2.0)+
      rep_vector(2.0,N).*(a-b).*(rep_vector(0.5,N).*f.*c.*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*f.*d.*exp(f.*d.*(Xm-e)))^3.0./
      ((rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^
      (rep_vector(1.0,N)./f).*f.*(rep_vector(0.5,N).*exp(f.*c.*(Xm-e))+
      rep_vector(0.5,N).*exp(f.*d.*(Xm-e)))^3.0)"
      
      returnmu_d1 <- paste0(returnmu_d1, ";")
      returnmu_d2 <- paste0(returnmu_d2, ";")
      returnmu_d3 <- paste0(returnmu_d1, ";")
    } # if(select_model == 'pb2') {
    
    
    # a - asymtote
    # b - size at theta
    # c - s0
    # d - s1
    # e - time (theta)
    # f - gamma
    
    if (select_model == 'pb3') {
      funstring <- "a-((4.0*(a-b))./((exp(f.*(Xm-e))+
      exp(c.*(Xm-e))).*(1.0+exp(d.*(Xm-e)))))"
      if (utils::packageVersion('rstan') < 2.26)
        funstring <- gsub(".*", " .* ",
                          funstring,
                          fixed = T)
      returnmu    <- paste0("return ", "(",  funstring, ")")
      returnmu_d0 <- funstring
      returnmu_d1 <-
        "rep_vector(4.0,N).*(a-b).*(f.*exp(f.*(Xm-e))+
      c.*exp(c.*(Xm-e)))./((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^
      2.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e))))+
      rep_vector(4.0,N).*(a-b).*d.*exp(d.*(Xm-e))./
      ((exp(f.*(Xm-e))+exp(c.*(Xm-e))).*(rep_vector(1.0,N)+
      exp(d.*(Xm-e)))^2.0)"
      
      returnmu_d2 <-
        "-rep_vector(8.0,N).*(a-b).*(f.*exp(f.*(Xm-e))+
      c.*exp(c.*(Xm-e)))^2.0./((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^
      3.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e))))-
      rep_vector(8.0,N).*(a-b).*(f.*exp(f.*(Xm-e))+
      c.*exp(c.*(Xm-e))).*d.*exp(d.*(Xm-e))./((exp(f.*(Xm-e))+
      exp(c.*(Xm-e)))^2.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^2.0)+
      rep_vector(4.0,N).*(a-b).*(f^2.0.*exp(f.*(Xm-e))+c^2.0.*exp(c.*(Xm-e)))./
      ((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^2.0.*(rep_vector(1.0,N)+
      exp(d.*(Xm-e))))-rep_vector(8.0,N).*(a-b).*d^
      2.0.*exp(d.*(Xm-e))^2.0./((exp(f.*(Xm-e))+
      exp(c.*(Xm-e))).*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^3.0)+
      rep_vector(4.0,N).*(a-b).*d^2.0.*exp(d.*(Xm-e))./((exp(f.*(Xm-e))+
      exp(c.*(Xm-e))).*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^2.0)"
      
      returnmu_d3 <-
        "rep_vector(24.0,N).*(a-b).*(f.*exp(f.*(Xm-e))+
      c.*exp(c.*(Xm-e)))^3.0./((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^
      4.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e))))+
      rep_vector(24.0,N).*(a-b).*(f.*exp(f.*(Xm-e))+c.*exp(c.*(Xm-e)))^
      2.0.*d.*exp(d.*(Xm-e))./((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^
      3.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^2.0)-
      rep_vector(24.0,N).*(a-b).*(f.*exp(f.*(Xm-e))+
      c.*exp(c.*(Xm-e))).*(f^2.0.*exp(f.*(Xm-e))+c^2.0.*exp(c.*(Xm-e)))./
      ((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^
      3.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e))))+
      rep_vector(24.0,N).*(a-b).*(f.*exp(f.*(Xm-e))+
      c.*exp(c.*(Xm-e))).*d^2.0.*exp(d.*(Xm-e))^
      2.0./((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^
      2.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^3.0)-
      rep_vector(12.0,N).*(a-b).*(f^
      2.0.*exp(f.*(Xm-e))+c^2.0.*exp(c.*(Xm-e))).*d.*exp(d.*(Xm-e))./
      ((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^
      2.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^2.0)-
      rep_vector(12.0,N).*(a-b).*(f.*exp(f.*(Xm-e))+
      c.*exp(c.*(Xm-e))).*d^2.0.*exp(d.*(Xm-e))./((exp(f.*(Xm-e))+
      exp(c.*(Xm-e)))^2.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^2.0)+
      rep_vector(4.0,N).*(a-b).*(f^3.0.*exp(f.*(Xm-e))+c^
      3.0.*exp(c.*(Xm-e)))./((exp(f.*(Xm-e))+exp(c.*(Xm-e)))^
      2.0.*(rep_vector(1.0,N)+exp(d.*(Xm-e))))+
      rep_vector(24.0,N).*(a-b).*d^3.0.*exp(d.*(Xm-e))^3.0./
      ((exp(f.*(Xm-e))+exp(c.*(Xm-e))).*(rep_vector(1.0,N)+
      exp(d.*(Xm-e)))^4.0)-rep_vector(24.0,N).*(a-b).*d^3.0.*exp(d.*(Xm-e))^
      2.0./((exp(f.*(Xm-e))+
      exp(c.*(Xm-e))).*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^3.0)+
      rep_vector(4.0,N).*(a-b).*d^3.0.*exp(d.*(Xm-e))./
      ((exp(f.*(Xm-e))+exp(c.*(Xm-e))).*(rep_vector(1.0,N)+exp(d.*(Xm-e)))^2.0)"
      
      returnmu_d1 <- paste0(returnmu_d1, ";")
      returnmu_d2 <- paste0(returnmu_d2, ";")
      returnmu_d3 <- paste0(returnmu_d3, ";")
    } # if(select_model == 'pb3') {
    
    
    
    # a - asymptote
    # b - rate constant
    # c - time at midpoint
    
    if (select_model == 'logistic1') {
      funstring <- "a./(1+exp(-b.*(Xm-c)))"
      if (utils::packageVersion('rstan') < 2.26)
        funstring <-
          gsub(".*", " .* ", funstring, fixed = T)
      returnmu    <- paste0("return ", "(",  funstring, ")")
      returnmu_d0 <- funstring
      returnmu_d1 <-
        "a.*b.*exp(-b.*(Xm - c))./(1 + exp(-b.*(Xm - c)))^2.0"
      
      returnmu_d2 <-
        "rep_vector(2.0,N).*a.*b^2.0.*exp(-b.*(Xm - c))^2.0./
        (1 + exp(-b.*(Xm - c)))^3.0 -
        a.*b^2.0.*exp(-b.*(Xm - c))./
        (1 + exp(-b.*(Xm - c)))^2.0"
      
      returnmu_d3 <-
        "rep_vector(6.0,N).*a.*b^3.0.*exp(-b.*(Xm - c))^3.0./
        (1 + exp(-b.*(Xm - c)))^4.0 -
        rep_vector(6.0,N).*a.*b^3.0.*exp(-b.*(Xm - c))^2.0./
        (1 + exp(-b.*(Xm - c)))^3.0 +
        a.*b^3.0.*exp(-b.*(Xm - c))./
        (1 + exp(-b.*(Xm - c)))^2.0"
      
      returnmu_d1 <- paste0(returnmu_d1, ";")
      returnmu_d2 <- paste0(returnmu_d2, ";")
      returnmu_d3 <- paste0(returnmu_d3, ";")
    } # if(select_model == 'logistic1') {
    
    
    
    # a - asymtote
    # b - size at theta
    # c - s0
    # d - theta1
    # e - s1
    # f - theta2
    
    # wolfram suggests -> (b/(1+exp(-c*(x-d)))) + ((a-b)/(1+exp(-e*(x-f))))
    # maple  -> ((a-b)/(1+exp(-e*(x-f)))) + (b/(1+exp(-c*(x-d))))
    
    # wolfram suggests
    # (a - b)/(exp(-e (x - f)) + 1) + b/(exp(-c (x - d)) + 1)
    
    if (select_model == 'logistic2') {
      funstring <-
        "((a-b)./(1+exp(-e.*(Xm-f)))) + (b./(1+exp(-c.*(Xm-d))))"
      if (utils::packageVersion('rstan') < 2.26)
        funstring <-
          gsub(".*", " .* ", funstring, fixed = T)
      returnmu    <- paste0("return ", "(",  funstring, ")")
      returnmu_d0 <- funstring
      returnmu_d1 <-
        "(e .* (a - b) .* exp(e .*(f + Xm)))./(exp(e .*f) + 
      exp(e .*Xm))^2.0 + (b .*c .* exp(c.* (d + Xm)))./(exp(c .*d) + 
      exp(c.* Xm))^2.0"
      
      
      returnmu_d2 <-
        "(a - b) .* ((2.0 .* e^2.0 .* exp(-2.0 .* e .* (Xm - f))) ./
      (exp(e .*(-(Xm - f))) + 1.0)^3.0- (e^2.0 .* exp(e .* (-(Xm - f))))./
      (exp(e .* (-(Xm - f))) + 1.0)^2.0) + b .* ((2.0 .* c^2.0 .* 
      exp(-2.0 .* c .* (Xm - d)))./(e^(-c .* (Xm - d)) + 1.0)^3.0 - 
      (c^2.0 .* exp(-c .*(Xm - d)))./(exp(-c .* (Xm - d)) + 1.0)^2.0)"
      
      
      returnmu_d3 <-
        "(a - b) .* ((e^3.0 .* exp(e .*(-(Xm - f))))./
      (exp(e .* (-(Xm - f))) + 1.0)^2.0 - (6.0 .* e^3.0 .* 
      exp(-2.0 .* e .* (Xm - f)))./(exp(e .* (-(Xm - f))) + 1.0)^3.0 + 
      (6.0 .* e^3.0 .* exp(-3.0 .* e .* (Xm - f)))./
      (exp(e .* (-(Xm - f))) + 1.0)^4.0.0) + b .* 
      ((c^3.0 .* exp(-c .* (Xm - d)))./(exp(-c .* (Xm - d)) + 1.0)^2.0 - 
      (6.0 .* c^3.0 .* exp(-2.0 .* c .* (Xm - d)))./
      (exp(-c .* (Xm - d)) + 1.0)^3.0 + (6.0 .* c^3.0 .* 
      exp(-3.0 .* c (Xm - d)))./(exp(-c .* (Xm - d)) + 1.0)^4.0)"
      
      
      returnmu_d1 <- paste0(returnmu_d1, ";")
      returnmu_d2 <- paste0(returnmu_d2, ";")
      returnmu_d3 <- paste0(returnmu_d3, ";")
    } # if(select_model == 'logistic2') {
    
    
    
    # a - size at infancy
    # b - rate at infancy
    # c - time at infancy
    # d - size at preadolescence
    # e - rate at preadolescence
    # f - time at preadolescence
    # g - size at adolescence
    # h - rate at adolescence
    # i - time at adolescence
    
    if (select_model == 'logistic3') {
      funstring <-
        "(a ./ (1 + exp(-b .* (Xm - c)))) +
        (d ./ (1 + exp(-e .* (Xm - f)))) +
        (g ./ (1 + exp(-h .* (Xm - i))))"
      if (utils::packageVersion('rstan') < 2.26)
        funstring <-
          gsub(".*", " .* ", funstring, fixed = T)
      returnmu    <- paste0("return ", "(",  funstring, ")")
      returnmu_d0 <- funstring
      
      
      returnmu_d1 <-
        "(a .* b .* exp(-b .* (Xm - c)))./(exp(-b .* (Xm - c)) + 1)^2.0 +
        (d .* e .* exp(-e .* (Xm - f)))./(exp(-e .* (Xm - f)) + 1)^2.0 +
        (g .* h .* exp(-h .* (Xm - i)))./(exp(-h .* (Xm - i)) + 1)^2.0"
      
      returnmu_d2 <-
        "(a .* ((rep_vector(2.0,N) .* b^2.0 .*
                  exp(-rep_vector(2.0,N) .* b .* (Xm - c))) ./
                 (exp(-b .* (Xm - c)) + 1)^3.0 -
                 (b^2.0 .* exp(-b .* (Xm - c))) ./
                 (exp(-b .* (Xm - c)) + 1.0)^2.0)) +
        (d .* ((rep_vector(2.0,N) .* e^2.0 .*
                  exp(-rep_vector(2.0,N) .* e .* (Xm - f))) ./
                 (exp(-e .* (Xm - f)) + 1)^3.0 -
                 (e^2.0 .* exp(-e .* (Xm - f))) ./
                 (exp(-e .* (Xm - f)) + 1.0)^2.0)) +
        (g .* ((rep_vector(2.0,N) .* h^2.0 .*
                  exp(-rep_vector(2.0,N) .* h .* (Xm - i))) ./
                 (exp(-h .* (Xm - i)) + 1)^3.0 -
                 (h^2.0 .* exp(-h .* (Xm - i))) ./
                 (exp(-h .* (Xm - i)) + 1.0)^2.0)) "
      
      returnmu_d3 <- returnmu_d2
      
      
      
      
      returnmu_d1 <- paste0(returnmu_d1, ";")
      returnmu_d2 <- paste0(returnmu_d2, ";")
      returnmu_d3 <- paste0(returnmu_d3, ";")
    } # if(select_model == 'logistic3') {
    
    
    
   
    
    
    insert_getX_name <-
      paste0("  vector[N] Xm = ", getxname, "(Xp);")
    start_fun <-
      paste0(
        "\nvector ",
        spfncname,
        "(vector ",
        vector_X_name,
        ", ",
        fullabcsnames_v,
        ") {" ,
        "\n",
        "  int N = num_elements(Xp);",
        "\n",
        insert_getX_name,
        collapse = " "
      )
    
    start_fun <- gsub(",)" , ")" , start_fun, fixed = TRUE)
    endof_fun <- paste0("\n    ", returnmu,
                        ";", "\n  } // end of spline function", sep = " ")
    
    
    rcsfun <- paste(start_fun, endof_fun)
    rcsfun_raw <- rcsfun
    
    # Create function d0
    fnameout <- paste0(spfncname, "_", "d0")
    
    spl_fun_ford <-
      paste0(fnameout,
             "(vector ",
             vector_X_name,
             ", ",
             fullabcsnames_v,
             ")")
    spl_fun_ford <- gsub("vector", "", spl_fun_ford, fixed = T)
    spl_fun_ford <- gsub("[[:space:]]", "", spl_fun_ford)
    spl_fun_ford <- gsub(",)", ")", spl_fun_ford, fixed = T)
    
    
    spl_d0 <- create_internal_function_nonsitar(
      y = y,
      function_str = rcsfun,
      fname = spfncname,
      fnameout = fnameout,
      returnmu = returnmu,
      xfunsi = xfunsi,
      yfunsi = yfunsi,
      setxoffset = setxoffset,
      gsub_out_unscaled = NULL,
      spl_fun_ford = spl_fun_ford,
      body = returnmu_d0,
      decomp = decomp,
      fixedsi = fixedsi
    )
    
    # Create function d1
    fnameout <- paste0(spfncname, "_", "d1")
    spl_d1 <- create_internal_function_nonsitar(
      y = y,
      function_str = rcsfun,
      fname = spfncname,
      fnameout = fnameout,
      returnmu = returnmu,
      xfunsi = xfunsi,
      yfunsi = yfunsi,
      setxoffset = setxoffset,
      # setxoffset setxoffset_d0_noqr
      gsub_out_unscaled = NULL,
      # gsub_out_unscaled = c('QR', 'Spl')
      spl_fun_ford = spl_fun_ford,
      body = returnmu_d1,
      decomp = decomp,
      fixedsi = fixedsi
    )
    
    # Create function d2
    fnameout <- paste0(spfncname, "_", "d2")
    spl_d2 <- create_internal_function_nonsitar(
      y = y,
      function_str = rcsfun,
      fname = spfncname,
      fnameout = fnameout,
      returnmu = returnmu,
      xfunsi = xfunsi,
      yfunsi = yfunsi,
      setxoffset = setxoffset,
      # setxoffset setxoffset_d0_noqr
      gsub_out_unscaled = NULL,
      # gsub_out_unscaled = c('QR', 'Spl')
      spl_fun_ford = spl_fun_ford,
      body = returnmu_d2,
      decomp = decomp,
      fixedsi = fixedsi
    )
    
    
    rcsfunmultadd <- NULL
    
    if (utils::packageVersion('rstan') > 2.26 & is.null(decomp)) {
      rcsfun <- paste0(getx_fun,
                       rcsfun,
                       rcsfunmultadd,
                       spl_d0,
                       spl_d1,
                       spl_d2,
                       sep = "\n")
    }
    
    if (utils::packageVersion('rstan') > 2.26 & !is.null(decomp)) {
      if (decomp == 'QR') {
        if (add_funmats) {
          rcsfun <- paste0(getx_fun,
                           funmats,
                           rcsfun,
                           rcsfunmultadd,
                           spl_d0,
                           spl_d1,
                           spl_d2,
                           sep = "\n")
        } else if (!add_funmats) {
          rcsfun <- paste0(getx_fun,
                           rcsfun,
                           rcsfunmultadd,
                           spl_d0,
                           spl_d1,
                           spl_d2,
                           sep = "\n")
        }
      }
    }
    
    
    
  } # if(select_model != 'sitar') { # pb models
  
  
  
  
  
  #################
  extract_r_fun_from_scode <-
    function(xstaring, what = NULL, decomp, spfncname) {
      xstaring <- gsub("[[:space:]]" , "", xstaring)
      xstaring <- gsub(";" , ";\n", xstaring)
      xstaring <- gsub("\\{" , "{\n", xstaring)
      xstaring <- gsub("}" , "}\n", xstaring)
      xstaring <- gsub("vector[N]" , "", xstaring, fixed = T)
      xstaring <- gsub("vector" , "", xstaring, fixed = T)
      xstaring <- gsub("int" , "", xstaring, fixed = T)
      xstaring <- gsub("real" , "", xstaring, fixed = T)
      xstaring <- gsub(paste0("jp1;", "\n"), "", xstaring, fixed = T)
      xstaring <- gsub("rep_vector" , "rep", xstaring, fixed = T)
      xstaring <- gsub("rep_" , "rep", xstaring, fixed = T)
      xstaring <-
        gsub(
          "Xx[ia,ja]=(X[ia]-knots[ja]>0?X[ia]-knots[ja]:0);" ,
          "Xx[ia,ja]=ifelse(X[ia]-knots[ja]>0,X[ia]-knots[ja],0);",
          xstaring,
          fixed = T
        )
      xstaring <- gsub("num_elements" , "length", xstaring, fixed = T)
      xstaring <-
        gsub("matrix[N,nknots-1]Spl" ,
             "Spl=matrix(0,N,nknots-1)",
             xstaring,
             fixed = T)
      xstaring <-
        gsub("matrix[nknots-1,N]rcs" ,
             "rcs=matrix(0,nknots-1,N)",
             xstaring,
             fixed = T)
      xstaring <-
        gsub("matrix[N,nknots]Xx" ,
             "Xx=matrix(0,N, nknots)",
             xstaring,
             fixed = T)
      xstaring <-
        gsub("for(iain1:N)" , "for(ia in 1:N)", xstaring, fixed = T)
      xstaring <- gsub("for(jain1:nknots)" ,
                       "for(ja in 1:nknots)",
                       xstaring,
                       fixed = T)
      xstaring <- gsub(".*" , "*", xstaring, fixed = T)
      xstaring <- gsub("./" , "/", xstaring, fixed = T)
      funame__ <- strsplit(xstaring, "\\(")[[1]][1]
      xstaring <- gsub(funame__ , paste0(funame__, "<-function"),
                       xstaring, fixed = T)
      xstaring <- sub("//[^//]+$", "", xstaring)
      # To remove stanadlon ";
      xstaring <-
        gsub(paste0(";\n;\n", ""), ";\n", xstaring, fixed = T)
      xstaring <- gsub("[nknots]knots" , "knots", xstaring, fixed = T)
      
      if (!is.null(what)) {
        if (what == 'getX') {
          xstaring <- gsub(paste0("x;", "\n"), "", xstaring)
        }
        if (what == 'getKnots') {
          xstaring <- gsub("\\[", "c\\(", xstaring)
          xstaring <- gsub("\\]'", "\\)", xstaring)
        }
      }
      
      if (!is.null(decomp)) {
        if (decomp == 'QR') {
          xstaring <-
            gsub(paste0("matrix[N,QK]XQ;", "\n") , "", xstaring, fixed = T)
          xstaring <- gsub("matrix[N,QK]" , "", xstaring, fixed = T)
          xstaring <-
            gsub(paste0("matrix[QK,QK]", XR_inv_name, ";", "\n") ,
                 "",
                 xstaring,
                 fixed = T)
          xstaring <-
            gsub(paste0("matrix[QK,QK]XR;", "\n"), "", xstaring, fixed = T)
          xstaring <- gsub("qr_thin_Q" , "qr", xstaring, fixed = T)
          xstaring <- gsub("qr_thin_R" , "qr.R", xstaring, fixed = T)
          xstaring <-
            gsub(XR_inv_name,
                 "=inverse" ,
                 XR_inv_name,
                 "=chol2inv",
                 xstaring,
                 fixed = T)
        }
      }
      
      xstaring
    } # extract_r_fun_from_scode
  
  
  rcsfun_raw_str   <- extract_r_fun_from_scode(rcsfun_raw,
                                               what = NULL,
                                               decomp = decomp,
                                               spfncname = spfncname)
  spl_d0_str   <- extract_r_fun_from_scode(spl_d0,
                                           what = NULL,
                                           decomp = decomp,
                                           spfncname = spfncname)
  spl_d1_str   <- extract_r_fun_from_scode(spl_d1,
                                           what = NULL,
                                           decomp = decomp,
                                           spfncname = spfncname)
  spl_d2_str   <- extract_r_fun_from_scode(spl_d2,
                                           what = NULL,
                                           decomp = decomp,
                                           spfncname = spfncname)
  getX_str     <- extract_r_fun_from_scode(
    getx_fun_raw,
    what = 'getX',
    decomp = decomp,
    spfncname = spfncname
  )
  getknots_str <- NULL
  if (select_model == 'sitar' | select_model == 'rcs') {
    getknots_str <- extract_r_fun_from_scode(
      getknots_fun_raw,
      what = 'getKnots',
      decomp = decomp,
      spfncname = spfncname
    )
  }
  
  all_raw_str <- c(rcsfun_raw_str,
                   spl_d0_str,
                   spl_d1_str,
                   spl_d2_str,
                   getX_str,
                   getknots_str)
  
  
  if (!add_rcsfunmatqrinv_genquant) {
    out <- list(rcsfun = rcsfun, r_funs = all_raw_str)
  } else if (add_rcsfunmatqrinv_genquant) {
    out <- list(rcsfun = rcsfun,
                r_funs = all_raw_str,
                gq_funs = rcsfunmatqrinv_genquant)
  }
  
  out
}






