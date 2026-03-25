

#' An internal function to prepare formula
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
  match_sitar_a_form <- NULL;
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
  sigma_formula_manualsi <- NULL;
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
  acovcoefnames_gr_str <- NULL;
  bcovcoefnames_gr_str <- NULL;
  ccovcoefnames_gr_str <- NULL;
  dcovcoefnames_gr_str <- NULL;
  ecovcoefnames_gr_str <- NULL;
  fcovcoefnames_gr_str <- NULL;
  gcovcoefnames_gr_str <- NULL;
  hcovcoefnames_gr_str <- NULL;
  icovcoefnames_gr_str <- NULL;
  scovcoefnames_gr_str <- NULL;
  smat <- NULL;
  smat_intercept <- NULL;
  set_model_sigma_by_ba <- NULL;
  set_model_sigma_by_no <- NULL;
  set_model_sigma_by_ls <- NULL;
  sigma_formula_manual_prior_via_sigma_formula <- NULL;
  
  nys <- NULL;
  ysi <- NULL;
  SbasisN <- NULL;
  sigma_formula_manualsi_set <- NULL;
  d_as_random_only <- NULL;
  
  
  # SbasisN = nknots - 1 for nsp nsk and rcs
  
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
  
  
  # Need separate TRUE FALSE for random effects and sigma
  
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
        stop2c("For '_str' approach of setting the random effects,",
             "\n ", 
             " only the vertical bar '||' approach is allowed.",
             "\n ", 
             " Please check and correct the '", in_gr_strsi, "' argument ",
             "\n ", 
             " which is currently specified as ", ept(in_gr_strsi)
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
  for (i in 1:(SbasisN)) {
    if (i < (SbasisN)) {
      name1 <- paste0("s", i, sep = ",")
    }
    else {
      name1 <- paste0("s", i, sep = "")
    }
    snames[i] <- name1
  }
  
  
  ##########
  # covers all sitar models
  # if(grepl("^sitar", select_model)) {
  #   select_model <- 'sitar'
  # }
  
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

  abcselements  <- gsub(",)" , ")" , abcselements, fixed = TRUE)
  
  if (terms_rhssi == "NULL") {
    abcsformfit   <- paste0(y, " ~ ", abcselements)
  }
  
  if (terms_rhssi != "NULL") {
    abcsformfit   <- paste0(y, "|", terms_rhssi, " ~ ", abcselements)
  }
  
 

  if (!(is.na(univariate_by$by) |
        univariate_by$by == "NA") &
      !is.null(subindicatorsi)) {
    termsvia <- univariate_by$terms
    if (terms_rhssi == "NULL") {
      abcsformfit <- (paste0(y,  "|", termsvia, "(", subindicatorsi, ")",
                             " ~ ", abcselements))
    }
    if (terms_rhssi != "NULL") {
      abcsformfit <- (paste0(y,  "|", terms_rhssi, "+", termsvia, "(", 
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
  
  
  
  # add_sigma_by_mu
  # This below will allow setting prior using 
  # sigma_formula and sigma_formula_gr
  # Also, dpar_formulasi must be NULL when using sigma_formula_manualsi
  
  # when sigma_formula_manualsi = TRUE, then dpar_formulasi is always NULL
  # But if sigma_formulasi and sigma_formula_gr_strsi are used to set prior
  # then they must not be set as NULL

  # set_model_sigma_by_no
  if(sigma_formula_manualsi_set) {
    dpar_formulasi <- NULL
    if(sigma_formula_manual_prior_via_sigma_formula) {
      sigmafixed <- sigmafixed
      sigmarandom <- sigmarandom
    } else {
      sigmafixed <- NULL
      sigmarandom <- NULL
    }
  }
  
  
  
  
  if(is.null(sigma_formulasi[[1]]) | sigma_formulasi == 'NULL') {
    display_message <-  paste("Please specify the fixed effect structure 
                              for 'sigma' parameter ",
                              "\n ", 
                              " (by using the 'sigma_formula' argument) 
                              since you have specified the ",
                              "\n ", 
                              " random effect structure via 
                              the 'sigma_formula_gr' argument")
    display_message <- clean_text_spaces(display_message)
    if(!is.null(sigma_formula_grsi)) {
      stop(display_message)
    }  
  }
  
  
  

  
  if(!is.null(sigma_formulasi[[1]]) & !is.null(sigma_formula_grsi)) {
    if(!identical(substr(strsplit(sigma_formulasi, "~", 
                                  fixed = T)[[1]][2], 1, 1),
                  substr(strsplit(sigma_formula_grsi, "~", 
                                  fixed = T)[[1]][2], 1, 1))) {
      # commented out on 07 03 2024
      # stop("Formula for 'sigma_formula' and 'sigma_formula_gr' should be same",
      #      "\n ", 
      #      " in terms of Intercept i.e, both should be either ~ 0 or ~ 1")
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
    strpartstrx <- strsplit(x, "|", fixed = T)[[1]][1]
    strpartstrx_form <- strpartstrx[1]
    if(length(strpartstrx) >1 ) {
      strpartstrx_grpa <- strpartstrx[2:length(strpartstrx)]
      strpartstrx_grpa <- gsub("[()]", "", strpartstrx_grpa)
      xtemp <- paste(strpartstrx_form, strpartstrx_grpa, sep = "|")
    } else {
      xtemp <- strpartstrx_form
    }
    x <- gsub("[[:space:]]", "", xtemp)
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
  # if (!is.null(dfixed)) {
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
  

  # New
  if(d_as_random_only) {
    dcovmat <- NULL
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
 
  
  # New
  if(!is.null(drandom)) {
    if(d_as_random_only) {
      dform <- "d~0"
    }
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
  
  # 28 01 2024
  gr_varss <- id
 
  
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
      # 24.08.204 was ok but xxxxxx conflicted with inclx of Hmisc::rcspline.eval
      str  <- gsub(")+(", ")ZZZZZZZZZZ(", str, fixed = T)
      str_ <- sub("^[^ZZZZZZZZZZ]*ZZZZZZZZZZ", "", str)
      str_ <- gsub(")ZZZZZZZZZZ(", ")+(", str_, fixed = T)
      # 24.08.204
      # str <- gsub(")+(", ")xxxxxx(", str, fixed = T)
      # str_ <- sub("^[^xxxxxx]*xxxxxx", "", str)
      # str_ <- gsub(")xxxxxx(", ")+(", str_, fixed = T)
      form <- paste0(form, "+", str_)
    } else {
      form <- form
    }
    form <- gsub("[[:space:]]", "", form)
    return(form)
  }
  
  
  add_higher_level_str_id <- function(x) {
    extract_xx <- x
    extract_xx <- gsub("[[:space:]]", "", extract_xx)
    extract_xx <- strsplit(extract_xx, ")+(", fixed = T)[[1]][1]
    extract_xx <- gsub("\\)", "", extract_xx)
    gr_varss <- sub(".*\\|", "", extract_xx) 
    gr_varss
  }
  
  
  ##########################
  
  getforvasasits <-  function(xform, x_formula_gr_strsi, 
                              setparantopen = NULL,
                              setparantclose = NULL) {
    xform <- xform %>% gsub_space()
    xform_temp <- xform
    if(is.null(setparantopen)) setparantopen <- "("
    if(is.null(setparantclose)) setparantclose <- ")"
    xform_temp <- gsub("(", setparantopen, xform_temp, fixed = T)
    xform_temp <- gsub(")", setparantclose, xform_temp, fixed = T)
    xform_temp <- add_higher_level_str(xform_temp, x_formula_gr_strsi)
    xform_temp <- restore_paranthese_grgr_str_form(xform_temp, 
                                                   exclude_first = NULL)
    xform_gr_names <- lapply(xform_temp, 
                             get_x_random2_new, 
                             gsubit = setparantclose)[[1]]
    xform_gr_names_asitis <- lapply(xform_temp, get_x_random2_asitis)[[1]]
    xform <- add_higher_level_str(xform, x_formula_gr_strsi)
    xform <- restore_paranthese_grgr_str_form(xform, exclude_first = NULL)
    xform <- xform %>% gsub_space()
    xgr_varss <- add_higher_level_str_id(x_formula_gr_strsi)
    list(varss = xgr_varss, gr_names = xform_gr_names,
         gr_names_asitis = xform_gr_names_asitis)
  }
  
  #########################
 
  
  if(set_higher_levels) {
    if(a_formula_gr_strsi_present) {
      aform <- add_higher_level_str(aform, a_formula_gr_strsi)
      aform <- restore_paranthese_grgr_str_form(aform, exclude_first = NULL)
      getvassetc <- getforvasasits(aform, a_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
                                   )
      agr_varss             <- getvassetc[['varss']]
      aform_gr_names        <- getvassetc[['gr_names']]
      aform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      agr_varss <- aform_gr_names <- aform_gr_names_asitis <- NULL
    }
    
    if(b_formula_gr_strsi_present) {
      bform <- add_higher_level_str(bform, b_formula_gr_strsi)
      bform <- restore_paranthese_grgr_str_form(bform, exclude_first = NULL)
      getvassetc <- getforvasasits(bform, b_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      bgr_varss             <- getvassetc[['varss']]
      bform_gr_names        <- getvassetc[['gr_names']]
      bform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      bgr_varss <- bform_gr_names <- bform_gr_names_asitis <- NULL
    }
    
    if(c_formula_gr_strsi_present) {
      cform <- add_higher_level_str(cform, c_formula_gr_strsi)
      cform <- restore_paranthese_grgr_str_form(cform, exclude_first = NULL)
      getvassetc <- getforvasasits(cform, c_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      cgr_varss             <- getvassetc[['varss']]
      cform_gr_names        <- getvassetc[['gr_names']]
      cform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      cgr_varss <- cform_gr_names <- cform_gr_names_asitis <- NULL
    }
    
    if(d_formula_gr_strsi_present) {
      dform <- add_higher_level_str(dform, d_formula_gr_strsi)
      dform <- restore_paranthese_grgr_str_form(dform, exclude_first = NULL)
      getvassetc <- getforvasasits(dform, d_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      dgr_varss             <- getvassetc[['varss']]
      dform_gr_names        <- getvassetc[['gr_names']]
      dform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      dgr_varss <- dform_gr_names <- dform_gr_names_asitis <- NULL
    }
    
    if(e_formula_gr_strsi_present) {
      eform <- add_higher_level_str(eform, e_formula_gr_strsi)
      eform <- restore_paranthese_grgr_str_form(eform, exclude_first = NULL)
      getvassetc <- getforvasasits(eform, e_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      egr_varss             <- getvassetc[['varss']]
      eform_gr_names        <- getvassetc[['gr_names']]
      eform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      egr_varss <- eform_gr_names <- eform_gr_names_asitis <- NULL
    }
    
    if(f_formula_gr_strsi_present) {
      fform <- add_higher_level_str(fform, f_formula_gr_strsi)
      fform <- restore_paranthese_grgr_str_form(fform, exclude_first = NULL)
      getvassetc <- getforvasasits(fform, f_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      fgr_varss             <- getvassetc[['varss']]
      fform_gr_names        <- getvassetc[['gr_names']]
      fform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      fgr_varss <- fform_gr_names <- fform_gr_names_asitis <- NULL
    }
    
    
    if(g_formula_gr_strsi_present) {
      gform <- add_higher_level_str(gform, g_formula_gr_strsi)
      gform <- restore_paranthese_grgr_str_form(gform, exclude_first = NULL)
      getvassetc <- getforvasasits(gform, g_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      ggr_varss             <- getvassetc[['varss']]
      gform_gr_names        <- getvassetc[['gr_names']]
      gform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      ggr_varss <- gform_gr_names <- gform_gr_names_asitis <- NULL
    }
    
    if(h_formula_gr_strsi_present) {
      hform <- add_higher_level_str(hform, h_formula_gr_strsi)
      hform <- restore_paranthese_grgr_str_form(hform, exclude_first = NULL)
      getvassetc <- getforvasasits(hform, h_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      hgr_varss             <- getvassetc[['varss']]
      hform_gr_names        <- getvassetc[['gr_names']]
      hform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      hgr_varss <- hform_gr_names <- hform_gr_names_asitis <- NULL
    }
    
    if(i_formula_gr_strsi_present) {
      iform <- add_higher_level_str(iform, i_formula_gr_strsi)
      iform <- restore_paranthese_grgr_str_form(iform, exclude_first = NULL)
      getvassetc <- getforvasasits(iform, i_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      igr_varss             <- getvassetc[['varss']]
      iform_gr_names        <- getvassetc[['gr_names']]
      iform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      igr_varss <- iform_gr_names <- iform_gr_names_asitis <- NULL
    }
    
    if(s_formula_gr_strsi_present) {
      sform <- add_higher_level_str(sform, s_formula_gr_strsi)
      sform <- restore_paranthese_grgr_str_form(sform, exclude_first = NULL)
      getvassetc <- getforvasasits(sform, s_formula_gr_strsi,
                                   setparantopen = "(", 
                                   setparantclose = ")"
      )
      sgr_varss             <- getvassetc[['varss']]
      sform_gr_names        <- getvassetc[['gr_names']]
      sform_gr_names_asitis <- getvassetc[['gr_names_asitis']]
    } else {
      sgr_varss <- sform_gr_names <- sform_gr_names_asitis <- NULL
    }
    
    # 24.08.2024
    # now get the .gr_varss from the first element of the .form_gr_names
    # this is because now bsitar higher loop is from 1:.. and not 2:..
    agr_varss <- aform_gr_names[1]
    bgr_varss <- bform_gr_names[1]
    cgr_varss <- cform_gr_names[1]
    dgr_varss <- dform_gr_names[1]
    egr_varss <- eform_gr_names[1]
    fgr_varss <- fform_gr_names[1]
    ggr_varss <- gform_gr_names[1]
    hgr_varss <- hform_gr_names[1]
    igr_varss <- iform_gr_names[1]
    
    gr_varss <- c(agr_varss, bgr_varss, cgr_varss, dgr_varss, 
                  egr_varss, fgr_varss, ggr_varss,
                  hgr_varss, igr_varss)
    
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
      if(length(unique(gr_varss)) != 1) {
        gr_varssmsg <- paste0("Duplicate group identifiers ", 
                              paste(gr_varss, collapse = ", "), 
                              "\n",
                              "  This is because your gr(...) argument",
                              "\n",
                              "  is not same across a b c ... 'gr_str' formulae")
        
        stop(gr_varssmsg)
      }
    }
    
    hierarchical_gr_names_asitis <- unique(hierarchical_gr_names_asitis)
    gr_varss <- gr_varss
  } # if(set_higher_levels) {
  
 
  
  if(sigma_set_higher_levels) {
    if(sigma_formula_gr_strsi_present) {
      # 24.08.2024
      # sigmaform_gr_names and sigmaform_gr_names_asitis is collected as follow
      # The "_________" does not allow split at ( or )
      sigmaform <- sigmaform %>% gsub_space()
      sigmaform_temp <- sigmaform
      setparantopen <- "######"
      setparantclose <- "########"
      sigmaform_temp <- gsub("(", setparantopen, sigmaform_temp, fixed = T)
      sigmaform_temp <- gsub(")", setparantclose, sigmaform_temp, fixed = T)
      sigmaform_temp <- add_higher_level_str(sigmaform_temp, 
                                             sigma_formula_gr_strsi)
      sigmaform_temp <- restore_paranthese_grgr_str_form(sigmaform_temp, 
                                                         exclude_first = FALSE)
      
      sigmaform_gr_names        <- lapply(sigmaform_temp,
                                          get_x_random2_new,
                                          gsubit = setparantclose)[[1]]
      
      sigmaform_gr_names_asitis <- lapply(sigmaform_temp, 
                                          get_x_random2_asitis)[[1]]
      sigmaform <- add_higher_level_str(sigmaform, sigma_formula_gr_strsi)
      sigmaform <- restore_paranthese_grgr_str_form(sigmaform, 
                                                    exclude_first = TRUE)
      sigmaform <- sigmaform %>% gsub_space()
      sigmagr_varss <- add_higher_level_str_id(sigma_formula_gr_strsi)
    } else {
      sigmaform_gr_names <- sigmaform_gr_names_asitis <- NULL
    }
    
    sigma_hierarchical_gr_names <- c(sigmaform_gr_names)
    sigma_hierarchical_gr_names <- unique(sigma_hierarchical_gr_names)
    
    sigma_hierarchical_gr_names_asitis <- c(sigmaform_gr_names_asitis)
    sigma_hierarchical_gr_names_asitis <- 
      unique(sigma_hierarchical_gr_names_asitis)
     # 24.08.2024
    if(length(sigmaform_gr_names) > 1) {
      sigma_gr_varss <- sigmaform_gr_names[1]
    } else {
      sigma_gr_varss <- sigma_gr_varss
    }
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
  
  
  # If no custom_family, then only add family to the call
  if(is.null(custom_family)) {
    if (!is.null(familysi)) {
      setbformula <- paste0(setbformula, "+", familysi)
    }
  }
  
  
  if(sigma_formula_manualsi != "NULL") {
    setbformula <- paste0(setbformula, "+", sigma_formula_manualsi)
  }
  
  # if (!is.null(familysi)) {
  #   setbformula <- paste0(setbformula, "+", familysi)
  # }
  
 
  
  group_arg_groupvar <- gr_varss
  multivariate_rescor <-  multivariate$rescor
  univariate_by_by <- univariate_by$by
  
  
 
  
  
  if(sigma_set_higher_levels) { 
    # 24.08.2024
    # now hierarchical_gr_names_asitis is same as gr_varss
    # this is because now bsitar higher loop is from 1:.. and not 2:..
    # sigma_arg_groupvar <- sigma_gr_varss_asitis
    
    sigma_arg_groupvar <- sigma_gr_varss
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
  
 
  # 24.08.2024 -> !sigma_set_higher_levels to sigma_set_higher_levels
  # if(!is.null(sigma_formula_grsi)) {
  #   if(!sigma_set_higher_levels & 
  #      grepl("|", sigma_formula_grsi, fixed = TRUE)) {
  #     sigma_arg_groupvar <- sigma_gr_varss
  #   } 
  #   if(!sigma_set_higher_levels & 
  #      !grepl("|", sigma_formula_grsi, fixed = TRUE)) {
  #     sigma_arg_groupvar <- gr_varss
  #   }  
  # }
  
  
  if(!is.null(sigma_formula_grsi)) {
    if(sigma_set_higher_levels & 
       grepl("|", sigma_formula_grsi, fixed = TRUE)) {
      sigma_arg_groupvar <- sigma_gr_varss
    } 
    if(sigma_set_higher_levels & 
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
  
  covariates <- c(a_covariate, b_covariate, c_covariate, d_covariate, 
                  e_covariate, f_covariate, g_covariate, 
                  h_covariate, i_covariate, 
                  s_covariate)
  
  covariates_ <- unique(covariates)
  
  # sigmacovariates  <- getcovlist(sigma_formulasi)
  # sigmacovariates_ <- unique(sigmacovariates)
  
  
  if(sigma_formula_manualsi_set) {
    sigmacovariates <- add_default_args_to_nlf_lf(sigma_formula_manualsi,
                                                  nys = nys,
                                                  ysi = ysi,
                                                  extract_covar = TRUE,
                                                  extract_nlpar = FALSE,
                                                  data_varnames = NULL,
                                                  verbose = FALSE)
  } else {
    sigmacovariates  <- getcovlist(sigma_formulasi)
  }
  sigmacovariates_ <- unique(sigmacovariates)
  
  
  # if(sigma_formula_manualsi_set) {
  #   if(set_model_sigma_by_ls) {
  #     sigmacovariates <- add_default_args_to_nlf_lf(sigma_formula_manualsi,
  #                                                   nys = nys,
  #                                                   ysi = ysi,
  #                                                   extract_covar = TRUE,
  #                                                   extract_nlpar = FALSE,
  #                                                   data_varnames = NULL,
  #                                                   verbose = FALSE)
  #   } else {
  #     sigmacovariates  <- getcovlist(sigma_formulasi)
  #   }
  # } else {
  #   sigmacovariates  <- getcovlist(sigma_formulasi)
  # }
  # sigmacovariates_ <- unique(sigmacovariates)
  
  
  
  
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
        # inte_cov <- paste0(main_cov, ":", "mat_s")
        # Even when ~0 or ~1, interaction is excluded for the first covariate
        inte_cov <- paste0(main_cov[2:length(main_cov)], ":", "mat_s")
        main_inte_cov <- c(main_cov, "mat_s", inte_cov)
        main_inte_cov <- paste(main_inte_cov, collapse = "+")
        # if (grepl("~0", s_formulasi, fixed = T)) {
        # here and below in else if added a_formulasi
        if (grepl("~0", s_formulasi, fixed = T) |
            grepl("~0", a_formulasi, fixed = T)
            ) {
          lmform  <- as.formula(paste0(y, "~0+",  main_inte_cov))
        } else if (!grepl("~0", s_formulasi, fixed = T) &
                   !grepl("~0", a_formulasi, fixed = T)
                   ) {
          lmform  <- as.formula(paste0(y, "~1+",  main_inte_cov))
        }
      }
    }
    
    
    if (grepl("^~1$", s_formulasi)) {
      if (length(a_covariate) == 1) {
        # changed 's_formulasi' to 'a_formulasi'
        if (grepl("~0", a_formulasi, fixed = T)) { 
          lmform  <- as.formula(paste0(y, "~0+",  a_covariate, "+", "mat_s"))
          # changed 's_formulasi' to 'a_formulasi'
        } else if (!grepl("~0", a_formulasi, fixed = T)) { 
          lmform  <- as.formula(paste0(y, "~1+",  a_covariate, "+", "mat_s"))
        }
      } else if (length(a_covariate) > 1) {
        main_cov      <- a_covariate
        main_inte_cov <- c(main_cov, "mat_s")
        main_inte_cov <- paste(main_inte_cov, collapse = "+")
        # here and below in else if added a_formulasi
        if (grepl("~0", s_formulasi, fixed = T) |
            grepl("~0", a_formulasi, fixed = T)
        ) {
          lmform  <- as.formula(paste0(y, "~0+",  main_inte_cov))
        } else if (!grepl("~0", s_formulasi, fixed = T) &
                   !grepl("~0", a_formulasi, fixed = T)
        ) {
          lmform  <- as.formula(paste0(y, "~1+",  main_inte_cov))
        }
      }
    }
    
    
    
    if (grepl("^~1$", s_formulasi)) {
      s_covariate <- NULL
    }
      
    
    
    if (grepl("^~1$", a_formulasi)) {
      if (grepl("~0", a_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~0+", "mat_s"))
      } else if (!grepl("~0", a_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~1+", "mat_s"))
      }
    }
    
    # check and adjust intercept for smat nsp nk
    if(smat == 'nsk' | smat == 'nsp') {
      if(smat_intercept) {
        if (grepl("~0", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~0+", "cons:mat_s"))
        } 
        if (grepl("~1", s_formulasi, fixed = T)) {
          lmform  <- as.formula(paste0(y, "~0+", "mat_s"))
        }
      } # if(smat_intercept) {
    } # if(smat == 'nsk'smat == 'nsp') {
    
    
    
    lm_fit    <- lm(lmform, data = data)
    lm_coef   <- coef(lm_fit)
    
    lm_rsd    <- summary(lm_fit)$sigma
    
    lm_coefse <- summary(lm_fit)$coefficients[, 2]
    
    
    # library(ggplot2)
    # dataset_in %>%
    #   ggplot(., aes(x = age)) + geom_line(aes(y = ppp, group = class))
    

    # lme
    enverr. <- parent.frame()
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
    


    if (!is.null(dfixed)) {
      lmform_dparm  <- as.formula(paste0(y, "~0+", "dcovmat", ":", x))
      # lmform_dparm  <- as.formula(paste0(y, "~0+", "", x))
      lm_fit_dparm  <- lm(lmform_dparm, data = data)
      lm_dparm_coef <- coef(lm_fit_dparm)
      # belwo if as.formula(paste0(y, "~0+", "dcovmat", "+", x))
      # lm_d_all      <- lm_dparm_coef[2:(ncol(dcovmat)+1)] # note 2: ...+1 
      # belwo if as.formula(paste0(y, "~0+", "dcovmat", ":", x))
      lm_d_all      <- lm_dparm_coef[1:ncol(dcovmat)]

    } # if (!is.null(dfixed)) {
    
    
    if (is.null(dfixed)) {
      lm_d_all <- NULL
    }
    
    
   
    if(d_as_random_only) {
      lm_d_all <- NULL
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
    
    if (!is.null(dfixed) & !is.null(lm_d_all)) {
      lm_d_all <- lm_d_all # rep(0, ncol(dcovmat))
    } else {
      lm_d_all <- NULL
    }
    

    lm_s_all <- lm_coef[(ncol(acovmat) + 1):length(lm_coef)]
    
    lm_s_se_all <- lm_coefse[(ncol(acovmat) + 1):length(lm_coefse)]
    
    
    
    if (grepl("~1", a_formulasi, fixed = T)) {
      if(match_sitar_a_form) lm_a_all[1] <- 
          lm_a_all[1] + lm_s_all[1] * min(knots)
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
    lm_sdx_all    <- sd(data[[y]]) / sds_X
    
    lm_sdx_se_all <- lm_s_se_all / sds_X
    
    
    names(lm_s_all)   <- names_mat_s_scovmat
    names(lm_sdx_all) <- names_mat_s_scovmat
    
    names(lm_s_se_all) <- names_mat_s_scovmat
    names(lm_sdx_se_all) <- names_mat_s_scovmat

    
    
    lm_s   <- lm_s_all[1:(SbasisN)]
    lm_sdx <- lm_sdx_all[1:(SbasisN)]
    
    lm_s_se   <- lm_s_se_all[1:(SbasisN)]
    lm_sdx_se <- lm_sdx_se_all[1:(SbasisN)]
    
    #########################################
    
    
    
    
    #########################################
    
    # zzzzzz
    
    if (!is.null(s_covariate) & length(s_covariate) > 1) {
      lm_s_cov <- lm_s_all[SbasisN:length(lm_s_all)]
      lm_sdx_cov <- lm_sdx_all[SbasisN:length(lm_sdx_all)]
      
      # lm_s_cov <- lm_s_all[nknots:length(lm_s_all)]
      # lm_sdx_cov <- lm_sdx_all[nknots:length(lm_sdx_all)]
     
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
        lm_s <- lm_s_all[1:(SbasisN)]
        lm_sdx <- lm_sdx_all[1:(SbasisN)]
        if (length(lm_s_all) > (SbasisN)) {
          lm_s_cov <- lm_s_all[SbasisN:length(lm_s_all)]
          lm_sdx_cov <- lm_sdx_all[SbasisN:length(lm_sdx_all)]
          tnames_s <-
            names_mat_s_scovmat[SbasisN:length(names_mat_s_scovmat)]
          names(lm_s_cov) <- tnames_s
          tnames_sdx <-
            names_mat_s_scovmat[SbasisN:length(names_mat_s_scovmat)]
          names(lm_sdx_cov) <- tnames_sdx
          
          # lm_s_cov <- lm_s_all[nknots:length(lm_s_all)]
          # lm_sdx_cov <- lm_sdx_all[nknots:length(lm_sdx_all)]
          # tnames_s <-
          #   names_mat_s_scovmat[nknots:length(names_mat_s_scovmat)]
          # names(lm_s_cov) <- tnames_s
          # tnames_sdx <-
          #   names_mat_s_scovmat[nknots:length(names_mat_s_scovmat)]
          # names(lm_sdx_cov) <- tnames_sdx
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
    
    ############################
    
    if (!is.null(s_covariate) & length(s_covariate) > 1) {
      lm_s_se_cov <- lm_s_se_all[SbasisN:length(lm_s_se_all)]
      lm_sdx_se_cov <- lm_sdx_se_all[SbasisN:length(lm_sdx_se_all)]
      
      inname_c_se_all <- c()
      for (inname in paste0("s", 1:df)) {
        t <- names(lm_a_all)[2:length(names(lm_a_all))]
        inname_c_se_all <- c(inname_c_se_all, paste0(inname, "_", t))
      }
      lm_s_se_cov   <- lm_s_se_cov[order(factor(names(lm_s_se_cov),
                                                levels = inname_c_se_all))]
      lm_sdx_se_cov <- lm_sdx_se_cov[order(factor(names(lm_sdx_se_cov),
                                                  levels = inname_c_se_all))]
      
      lm_s_se_all <- lm_sdx_se_all <- c()
      for (idfi in 1:df) {
        lm_s_se_all <- c(lm_s_se_all, c(lm_s_se[idfi],
                                        lm_s_se_cov[grep(paste0("s", idfi, "_"),
                                                         names(lm_s_se_cov))]))
        lm_sdx_se_all <- c(lm_sdx_se_all, c(lm_sdx_se[idfi],
                                            lm_sdx_se_cov[grep(paste0("s", idfi, "_"),
                                                               names(lm_sdx_se_cov))]))
      }
    } else if ((!is.null(s_covariate) &
                length(s_covariate) == 1) |
               is.null(s_covariate)) {
      if (!grepl("~0", s_formulasi, fixed = T)) {
        lm_s_se_all <- lm_coefse[(ncol(acovmat) + 1):length(lm_coefse)]
        names(lm_s_se_all) <- names_mat_s_scovmat
        names(lm_sdx_se_all) <- names_mat_s_scovmat
        lm_s_se <- lm_s_se_all[1:(SbasisN)]
        lm_sdx_se <- lm_sdx_se_all[1:(SbasisN)]
        if (length(lm_s_se_all) > (SbasisN)) {
          lm_s_se_cov <- lm_s_se_all[SbasisN:length(lm_s_se_all)]
          lm_sdx_se_cov <- lm_sdx_se_all[SbasisN:length(lm_sdx_se_all)]
          tnames_s <-
            names_mat_s_scovmat[SbasisN:length(names_mat_s_scovmat)]
          names(lm_s_se_cov) <- tnames_s
          tnames_sdx <-
            names_mat_s_scovmat[SbasisN:length(names_mat_s_scovmat)]
          names(lm_sdx_se_cov) <- tnames_sdx
        } else {
          lm_s_se_cov <- NULL
          lm_sdx_se_cov <- NULL
        }
      }
      
      if (grepl("~0", s_formulasi, fixed = T)) {
        lm_s_se_all <- lm_coefse[(ncol(acovmat) + 1):length(lm_coefse)]
        names(lm_s_se_all) <- names_mat_s_scovmat
        names(lm_sdx_se_all) <- names_mat_s_scovmat
        inname_c_se_all <- c()
        for (inname in paste0("s", 1:df)) {
          t <- names(lm_a_all)[1:length(names(lm_a_all))]
          inname_c_se_all <- c(inname_c_se_all, paste0(inname, "_", t))
        }
        lm_s_se <- NULL
        lm_s_se_cov <- NULL
        lm_sdx_se <- NULL
        lm_sdx_se_cov <- NULL
      }
    }
    ############################
    
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
  
  
  if(is.null(lm_s_se_all)) {
    assign(paste0('lm_', 'sdx_se'), NULL)
    assign(paste0('lm_', 'sdx_se', '_all'), NULL)
    assign(paste0('lm_', 'sdx_se', '_cov'), NULL)
  }
  
  # brms removes white spaces from the coefficient names
  # mimicking that behaviors but keeping it separate here as 
  # brms may later change it to underscore or something else
  
  
  gsubitbt <- ""
  gsubitbt_fun <- function(x, gsubitbt) {
    gsubitbtx <- function(x, gsubitbt) gsub("[[:space:]]", gsubitbt, x)
    if(is.list(x)) out <- lapply(x, gsubitbtx, gsubitbt)
    if(!is.list(x)) out <- lapply(x, gsubitbtx, gsubitbt) %>% unlist()
    out
  }
  
  
 
  
  for (set_randomsi_higher_levsli in c(letters[1:26], 'sigma', 'dpar')) {
    set_nlpar_what <- set_randomsi_higher_levsli
    
    if(exists(paste0(set_nlpar_what, 'covcoefnames'))) {
      if(!is.null(ept(paste0(set_nlpar_what, 'covcoefnames'))) | 
         ept(paste0(set_nlpar_what, 'covcoefnames') != "" )) {
        assign(paste0(set_nlpar_what, 'covcoefnames'), 
               gsubitbt_fun(ept(paste0(set_nlpar_what, 'covcoefnames')), 
                            gsubitbt = gsubitbt))
      } # if(!is.null(
    } # if(exists(... covcoefnames
    
    if(exists(paste0(set_nlpar_what, 'covcoefnames_gr'))) {
      if(!is.null(ept(paste0(set_nlpar_what, 'covcoefnames_gr'))) | 
         ept(paste0(set_nlpar_what, 'covcoefnames_gr') != "" )) {
        assign(paste0(set_nlpar_what, 'covcoefnames_gr'), 
               gsubitbt_fun(ept(paste0(set_nlpar_what, 'covcoefnames_gr')), 
                            gsubitbt = gsubitbt))
      } # if(!is.null(
    } # if(exists(... covcoefnames_gr
    
    if(exists(paste0(set_nlpar_what, 'covcoefnames_gr_str'))) {
      if(!is.null(ept(paste0(set_nlpar_what, 'covcoefnames_gr_str'))) | 
         ept(paste0(set_nlpar_what, 'covcoefnames_gr_str') != "" )) {
        assign(paste0(set_nlpar_what, 'covcoefnames_gr_str'), 
               gsubitbt_fun(ept(paste0(set_nlpar_what, 'covcoefnames_gr_str')), 
                            gsubitbt = gsubitbt))
      } # if(!is.null(
    } # if(exists(... covcoefnames_gr_str
    
    if(exists(paste0(set_nlpar_what, 'covcoefnames_gr_str'))) {
      if(!is.null(ept(paste0(set_nlpar_what, 'covcoefnames_gr_str'))) | 
         ept(paste0(set_nlpar_what, 'covcoefnames_gr_str') != "" )) {
        assign(paste0(set_nlpar_what, 'covcoefnames_gr_str'), 
               gsubitbt_fun(ept(paste0(set_nlpar_what, 'covcoefnames_gr_str')), 
                            gsubitbt = gsubitbt))
      } # if(!is.null(
    } # if(exists(... covcoefnames_gr_str
  } # for (set_randomsi_higher_levsli
  
  
  
  
  # brms replace "=" with "EQ" and removes ","
  # When using function such as Hmisc::rcspline.eval(age, nk=4,inclx = FALSE)
  # Below we are doing gsub fro sigma only, but all names can be part of the 
  # set_vector_for_gsub_names
  
  # Also brms replaces - sign with M e.g., age-mean(age), also + with P
  
  # for (i in 1:length(sigmacovcoefnames)) {
  #   sigmacovcoefnames[i] <- gsub("=", "EQ", sigmacovcoefnames[i])
  #   sigmacovcoefnames[i] <- gsub(",", "", sigmacovcoefnames[i])
  # }
  
  set_vector_for_gsub_names <- c('sigmacovcoefnames',
                                 'sigmacovcoefnames_gr',
                                 'sigmacovcoefnames_gr_str')
  
  for (set_vector_for_gsub_names_i in 1:length(set_vector_for_gsub_names)) {
    setnamesxx <- get(set_vector_for_gsub_names[set_vector_for_gsub_names_i])
    if(!is.null(setnamesxx)) {
      if(is.list(setnamesxx)) {
        for (ij in 1:length(setnamesxx)) {
          for (i in 1:length(setnamesxx [[ij]] )) {
            setnamesxx[[ij]] [i] <- gsub("=", "EQ", setnamesxx[[ij]] [i])
            setnamesxx[[ij]] [i] <- gsub(",", "",   setnamesxx[[ij]] [i])
            setnamesxx[[ij]] [i] <- gsub("[<>-]", "M",   setnamesxx[[ij]] [i])
            setnamesxx[[ij]] [i] <- gsub("[<>+]", "P",   setnamesxx[[ij]] [i])
          }
        }
      } else if(!is.list(setnamesxx)) {
        setnamesxx <- setnamesxx
        for (i in 1:length(setnamesxx)) {
          setnamesxx[i] <- gsub("=", "EQ", setnamesxx[i])
          setnamesxx[i] <- gsub(",", "",   setnamesxx[i])
          setnamesxx[i] <- gsub("[<>-]", "M",  setnamesxx[i])
          setnamesxx[i] <- gsub("[<>+]", "P",  setnamesxx[i])
        }
        setnamesxx <- setnamesxx
      }
    }
    assign(set_vector_for_gsub_names[set_vector_for_gsub_names_i], setnamesxx)
  }
  
  
#  mysr <- "(1|11|gr(id,by=classid))+(1|33|gr(study,by=NULL))"
  
  
  
  #####################################
  
  # get_all_grby_vars_names <- function(elements = NULL, envir = NULL) {
  #   if(is.null(elements)) {
  #     elements <- letters[1:12]
  #     elements <- c(elements, 'sigma')
  #   }
  #   if(is.null(envir)) {
  #     envir <- parent.frame()
  #   }
  #   
  #   abc_grby_vars_grsi_c <- abc_grby_vars_gr_strsi_c <- c()
  #   sigma_grby_vars_grsi_c <- sigma_grby_vars_gr_strsi_c <- c()
  #   
  #   for (i in elements) {
  #     
  #     if(i != 'sigma') {
  #       
  #       if(exists(paste0(i, "_formula_grsi"), envir = envir)) {
  #         mysr <- get(paste0(i, "_formula_grsi"), envir = envir)
  #         if(is.null(mysr)) {
  #           # 
  #         } else if(!is.null(mysr)) {
  #           mysri <- regmatches(mysr, gregexpr("(?<=by=)[^,)]+", 
  #                                              mysr, perl = TRUE))[[1]]
  #           mysri <- mysri[mysri != "NULL"]
  #           abc_grby_vars_grsi_c <- c(abc_grby_vars_grsi_c, mysri)
  #         } # if(is.null(mysr)) {
  #       } # if(exists(paste0(i, "_formula_grsi"))) {
  #       
  #       if(exists(paste0(i, "_formula_gr_strsi"), envir = envir)) {
  #         mysr <- get(paste0(i, "_formula_gr_strsi"), envir = envir)
  #         if(is.null(mysr)) {
  #           # 
  #         } else if(!is.null(mysr)) {
  #           mysri <- regmatches(mysr, gregexpr("(?<=by=)[^,)]+", 
  #                                              mysr, perl = TRUE))[[1]]
  #           mysri <- mysri[mysri != "NULL"]
  #           abc_grby_vars_gr_strsi_c <- c(abc_grby_vars_gr_strsi_c, mysri)
  #         } # if(is.null(mysr)) {
  #       } # if(exists(paste0(i, "_formula_gr_strsi"))) {
  #       
  #     } else if(i == 'sigma') {
  #       if(exists(paste0(i, "_formula_grsi"), envir = envir)) {
  #         mysr <- get(paste0(i, "_formula_grsi"), envir = envir)
  #         if(is.null(mysr)) {
  #           # 
  #         } else if(!is.null(mysr)) {
  #           mysri <- regmatches(mysr, gregexpr("(?<=by=)[^,)]+", 
  #                                              mysr, perl = TRUE))[[1]]
  #           mysri <- mysri[mysri != "NULL"]
  #           sigma_grby_vars_grsi_c <- c(sigma_grby_vars_grsi_c, mysri)
  #         } # if(is.null(mysr)) {
  #       } # if(exists(paste0(i, "_formula_grsi"))) {
  #       
  #       if(exists(paste0(i, "_formula_gr_strsi"), envir = envir)) {
  #         mysr <- get(paste0(i, "_formula_gr_strsi"), envir = envir)
  #         if(is.null(mysr)) {
  #           # 
  #         } else if(!is.null(mysr)) {
  #           mysri <- regmatches(mysr, gregexpr("(?<=by=)[^,)]+", 
  #                                              mysr, perl = TRUE))[[1]]
  #           mysri <- mysri[mysri != "NULL"]
  #           sigma_grby_vars_gr_strsi_c <- c(sigma_grby_vars_gr_strsi_c, mysri)
  #         } # if(is.null(mysr)) {
  #       } # if(exists(paste0(i, "_formula_gr_strsi"))) {
  #       
  #     } # if(i != 'sigma') { else if(i == 'sigma') {
  #     
  #     
  #   } # for (i in elements) {
  #   
  #   
  #   out <- list()
  #   # out[['abc_grsi']] <- abc_grby_vars_grsi_c
  #   # out[['abc_gr_strsi']] <- abc_grby_vars_gr_strsi_c
  #   # out[['sigma_grsi']] <- sigma_grby_vars_grsi_c
  #   # out[['sigma_gr_strsi']] <- sigma_grby_vars_gr_strsi_c
  #   
  #   abc_grby <- unique(abc_grby_vars_grsi_c, abc_grby_vars_gr_strsi_c)
  #   sigma_grby <- unique(sigma_grby_vars_grsi_c, sigma_grby_vars_gr_strsi_c)
  #   out[['abc_grby']] <- abc_grby_vars_grsi_c
  #   out[['sigma_grby']] <- abc_grby_vars_gr_strsi_c
  #   return(out)
  # }
  
  
  abc_grby_sigma_grby <- get_all_grby_vars_names(elements = NULL, envir = NULL)
  abc_grby   <- abc_grby_sigma_grby[['abc_grby']]
  sigma_grby <- abc_grby_sigma_grby[['sigma_grby']]

 
  
  #####################################

  
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
    sigmacovariates_ = sigmacovariates_,
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
    lme_rsd = lme_rsd,
    
    lm_s_se = lm_s_se,
    lm_s_se_cov = lm_s_se_cov,
    lm_s_se_all = lm_s_se_all,
    lm_sdx_se = lm_sdx_se,
    lm_sdx_se_cov = lm_sdx_se_cov,
    lm_sdx_se_all = lm_sdx_se_all,
    
    abc_grby = abc_grby,
    sigma_grby = sigma_grby
    
  ) # End list_out <- list(
  
 
  # add_sigma_by_mu
  # sigmaform which is actuallu lf(sigmatau ~) that was used for prior setting
  if(sigma_formula_manual_prior_via_sigma_formula) {
    sigmaform_rm <- gsub_space(sigmaform) 
    sigmaform_rm <- paste0(sigmaform_rm, ",")
    setbformula <- gsub(sigmaform_rm, "", setbformula, fixed = TRUE)
  }
 

  attr(setbformula, "list_out") <- as.list(list_out)
  return(setbformula)
}


