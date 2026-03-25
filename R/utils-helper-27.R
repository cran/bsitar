

##############################################################################
# get_form
##############################################################################

#' Title
#'
#' @param x Internal
#' @param method Internal
#' @param smat Internal
#' @param method Internal
#' @param knots Internal
#' @param bknots Internal
#' @param df Internal
#' @param degree Internal
#' @param intercept Internal
#' @param derivs Internal
#' @param centerval Internal
#' @param normalize Internal
#' @param preH Internal
#' @param sfirst Internal
#' @param sparse Internal
#' @param verbose Internal
#'
#' @returns A character string
#' 
#' @keywords internal
#' @noRd
#'
get_form <- function(x,
                     method = 'bs',
                     smat, 
                     knots, 
                     bknots, 
                     df = NULL,
                     degree = 3,
                     intercept = FALSE, 
                     derivs = 0,
                     centerval = FALSE, 
                     normalize = FALSE, 
                     preH = FALSE, 
                     sfirst = FALSE, 
                     sparse = FALSE,
                     verbose = FALSE) {
  
  x <- str2lang(x)
  
  if(is.symbol(smat)) {
    smat <- deparse(smat)
  }
  
  if (length(knots) > 0) {
    knots <- paste0("c(", paste(knots, collapse = ", "), ")")
    df <- NULL
  } else {
    if(is.null(df)) {
      if(smat == "rcs") {
        df <- 2 # For rcs, nk = df + 1, so minimum 3 nk
        knots <- ""
        bknots <- ""
      } else {
        df <- 1
        knots <- ""
      }
    } # if(is.null(df)) {
  } # if (length(knots) > 0) { else {
  
  
  knots <- gsub_space(knots)
  bknots <- gsub_space(bknots)
  knots <- sub('.*=', '', knots)
  bknots <- sub('.*=', '', bknots)
  
  if(is_emptyx(df))     df     <- NULL
  if(is_emptyx(knots))  knots  <- NULL
  if(is_emptyx(bknots)) bknots <- NULL
 
  if(is.character(df))     df     <- str2lang(df)
  if(is.character(knots)) knots <- str2lang(knots)
  if(is.character(bknots)) bknots <- str2lang(bknots)
  
  degree = degree
  intercept = intercept
  derivs = derivs
  centerval = centerval
  normalize = normalize
  preH = preH
  sfirst = sfirst
  sparse = sparse

  if(smat == 'ns') {
    SplineCall <- substitute(TEMPNAME(x = x,
                                      # df = df,
                                      knots = knots,
                                      Boundary.knots = bknots,
                                      intercept = intercept))
  } else {
    SplineCall <- substitute(TEMPNAME(x = x,
                                      knots = knots,
                                      bknots = bknots,
                                      degree = degree,
                                      intercept = intercept,
                                      derivs = derivs,
                                      centerval = centerval,
                                      normalize = normalize,
                                      preH = preH,
                                      sfirst = sfirst,
                                      sparse = sparse))
  }

  SplineCall$df <- df
  
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
  } else if(smat == 'ns') {
    SplineCall[[1]] <- quote(splines::ns)
  }
  
  environment(SplineCall) <- as.environment(getNamespace('bsitar'))
  
  SplineCall <- paste0(gsub_space(deparse(SplineCall)), collapse = "")
 
  return(SplineCall)
}

##############################################################################
# check_set_criteria
##############################################################################

#' Title
#'
#' @param icr_fn Internal
#'
#' @returns A character string a function
#' 
#' @keywords internal
#' @noRd
#'
check_set_criteria <- function(icr_fn, add_attr = FALSE, verbose = FALSE) {
  icr_fn_choices <- c("AIC", "BIC", "CV")
  if (missing(icr_fn)) {
    icr_fn <- stats::AIC
  }
  if(is.character(icr_fn)) {
    icr_fn <- gsub("\"", "", icr_fn)
    score_type <- icr_fn
    if(icr_fn == "AIC" | icr_fn == "stats::AIC") {
      icr_fn <- stats::AIC
    } else if(icr_fn == "BIC" | icr_fn == "stats::BIC") {
      icr_fn <- stats::BIC
    } else if(icr_fn == "CV") {
      icr_fn  <- "CV"
    } else {
      stop2c("For argument 'knots_selection', the option 'criteria' must be 
         one of the following: ", 
             collapse_comma(icr_fn_choices), 
             ". The current choice ", 
             collapse_comma(icr_fn),
             " is invalid")
    }
  } else if(is.function(icr_fn)) {
    icr_fn <- icr_fn
    score_type_temp <- deparse(substitute(icr_fn))
    score_type_temp_c <- paste0(score_type_temp, collapse = "")
    if(grepl("AIC", score_type_temp_c)) {
      score_type <- 'AIC'
    } else if(grepl("BIC", score_type_temp_c)) {
      score_type <- 'BIC'
    } else {
      score_type <- score_type_temp
    }
  } else {
    stop2c("Argument 'icr_fn' must be either a function, or a character 
         string from the following: ", collapse_comma(icr_fn_choices))
  }
  if(add_attr) {
    attr(icr_fn, 'score_type') <- score_type
  }
  return(icr_fn)
}


##############################################################################
# eval_icr_fn_fun
##############################################################################

#' Title
#'
#' @param fun Internal
#' @param model Internal
#' @param dataset Internal
#' @param cvk Internal
#' @param cviter Internal
#' @param forms Internal
#' @param cost Internal
#' @param verbose Internal
#'
#' @returns A character string a function
#' 
#' @keywords internal
#' @noRd
#'
eval_icr_fn_fun <- function(fun, model, dataset, cvk, cviter, forms = NULL, 
                            cost=function(y,yhat) mean((y-yhat)^2),
                            verbose = FALSE) {
  tmpgetstats <- NULL
  if(is.null(forms)) {
    nforms <- 1 
  } else {
    nforms <- length(forms)
  }
  if(is.character(fun)) {
    if(fun == "CV"){
      tmpgetstats <- NULL
      for(j in 1:cviter){
        mods <- list()
        for(i in 1:nforms){ # for(i in 1:length(forms)){
          tmpdat    <- dataset[rownames(stats::model.frame(model)), ]
          mods[[i]] <- boot::cv.glm(data = tmpdat, glmfit = model,
                                    cost=cost,
                                    K = cvk)
        }
        tmpgetstats <- rbind(tmpgetstats, sapply(mods, function(x)x$delta[1]))
      }
      icr_score <- colMeans(tmpgetstats)
    } else if(fun == "AIC"){
      icr_score <- stats::AIC(model)
    } else if(fun == "BIC"){
      icr_score <- stats::BIC(model)
    }
  } else if(is.function(fun)) {
    icr_score <- fun(model)
  }
  return(icr_score)
}



##############################################################################
# get_choose_model
# importFrom mfp fp
##############################################################################

#' Function assessing different types of regression models for a range of
#' degrees of freedom (knots in the case of splines), returning the one yielding
#' the best results according to an information criterion. The function uses
#' fractional polynomials, restricted cubic splines, and restricted cubic
#' splines with non-uniform knot placements.
#'
#' @param dataset The data frame
#'
#' @param dependent The dependent variable in the formula
#'
#' @param independents The independent variables in the formula
#'
#' @param icr_fn The information criterion function for comparing different
#'   models with different degress for freedom or knots (default AIC)
#'
#' @param cost_fn The criterion used to choose which knots to remove, passed to
#'   the function get_choose_removal. Defaults to AIC
#'   
#' @param search_df A logical (default \code{TRUE}) to indicate whether to
#'   first search for the optimal number of knots \code{df} by first running the
#'   \code{get_suggest_knotcount} function. When calling the
#'   \code{get_choose_model} from within the \pkg{bsitar}, the
#'   \code{search_df} is set to \code{FALSE}.
#'   
#' @param return_df A logical (default \code{FALSE}). Returns df which is the 
#' total number of knots (internal + boundary) minus one.
#'
#' @param max_nknots The max number of inner knots for restricted cubic splines
#'   (default 4)
#'
#'
#' @param verbose Verbose output, default FALSE
#'
#' @param ... Internal
#' 
#' @param method Internal
#' @param smat Internal
#' @param method Internal
#' @param df Internal
#' @param knots Internal
#' @param bknots Internal
#' @param degree Internal 
#' @param intercept Internal
#' @param derivs Internal
#' @param centerval Internal
#' @param normalize Internal
#' @param preH Internal
#' @param sfirst Internal
#' @param sparse Internal
#' @param nk Internal
#' @param inclx Internal
#' @param knots.only Internal
#' @param type Internal
#' @param norm Internal
#' @param rpm Internal
#' @param pc Internal
#' @param fractied Internal
#' @param bkrange Internal
#' @param fix_bknots Internal
#' @param bound Internal
#' @param userdata Internal
#'
#' @return A list with named elements, such as 'model', 'type', 'score'. The
#'   function returns a list with named elements and sublists, see examples for
#'   full overview of the returned values.
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' my_model <- get_choose_model(d, y, x)$model
#' result <- get_choose_model(d, y, x, icr_fn = AIC, verbose = FALSE)
#'
#' ret <- get_choose_model(d, y, x)
#'
#' ret$labels[[ret$type]] # Human readable name of chosen model type
#' ret[[ret$type]]        # Gives more values for the chosen model if available
#'
#' ret$model      # The chosen model
#' ret$score      # The chosen model's score
#' ret$type       # The type of model chosen as a string
#' ret$score_name # The type of score used as a string
#' ret$score_fn   # The function used for scores
#'
#' ret$labels     # Description string for the types of models
#'
#' ret$mfp        # Multivariate fractional polynomial
#' ret$mfp$model  # The model
#' ret$mfp$score  # The score
#'
#' ret$ns         # Restricted cubic splines from regular sequence of quantiles
#' ret$ns$model   # The model
#' ret$ns$score   # The score
#' ret$ns$knot_cnt_arg      # The number of inner knots as input argument
#' ret$ns$knot_cnt_distinct # The number of distinct placements in the result
#' ret$ns$knot_placements   # Knots and boundary knots
#' ret$ns$knot_placements$knots           # The knot placements as a list
#' ret$ns$knot_placements$Boundary.knots  # The boundary knots as a list
#'
#' ret$ns_nu         # Restricted cubic splines w/non-uniform placements
#' ret$ns_nu$model   # The model
#' ret$ns_nu$score   # The score
#' ret$ns_nu$knot_cnt_distinct # The number of distinct placements in the result
#' ret$ns_nu$knot_placements   # Knots and boundary knots
#' ret$ns_nu$knot_placements$knots          # The knot placements as a list
#' ret$ns_nu$knot_placements$Boundary.knots # The boundary knots as a list
#' 
get_choose_model <- function(dataset,
                         dependent,
                         independents,
                         ...,
                         icr_fn = stats::AIC,
                         cost_fn = stats::AIC,
                         cvk = 10, 
                         cviter = 10,
                         search_df = TRUE,
                         return_df = FALSE,
                         max_nknots = 4,
                         all_scores = FALSE,
                         plot_all_scores = FALSE,
                         method = 'bs',
                         smat = 'ns',
                         df = NULL, 
                         knots = NULL, 
                         bknots = NA, 
                         degree = 3, 
                         intercept = FALSE, 
                         derivs = 0, 
                         centerval = FALSE, 
                         normalize = FALSE, 
                         preH = FALSE, 
                         sfirst = FALSE, 
                         sparse = FALSE,
                         nk = NULL,
                         inclx = TRUE,
                         knots.only = TRUE,
                         type = "ordinary",
                         norm = 2,
                         rpm = NULL,
                         pc = FALSE,
                         fractied = 0.05,
                         bkrange = FALSE,
                         fix_bknots = TRUE,
                         bound = NULL,
                         userdata = NULL,
                         verbose = TRUE) {
  
  if(is.character(dependent))    dependent    <- str2lang(dependent)
  if(is.character(independents)) independents <- str2lang(independents)
  
  dependent <- rlang::enquo(dependent)
  independents <- rlang::enquo(independents)
  score_type <- NULL

  icr_fn     <- check_set_criteria(icr_fn, add_attr = TRUE, verbose = verbose)
  score_type <- attr(icr_fn, 'score_type')
  attr(icr_fn, 'score_type') <- NULL
  
  cost_fn   <- check_set_criteria(cost_fn, add_attr = FALSE, verbose = verbose)
 
  if (missing(max_nknots)) max_nknots <- 4
  if (missing(verbose)) verbose <- TRUE
  if (missing(bknots)) bknots <- NA
  
  ret_desc <- list(
    "ns_nu" = "Restricted cubic splines with freely placed knots",
    "ns" = "Restricted cubic splines with knots placed at quantiles")
  
  ret <- list(labels = ret_desc, score_fn = icr_fn, score_name = score_type)
  
  suppressWarnings({
    knotcnt_suggestion <-
      get_suggest_knotcount(dataset = dataset, 
                            dependent = dependent, 
                            independents = independents, 
                            max_nknots  = max_nknots,
                            icr_fn = icr_fn, 
                            all_scores = all_scores,
                            plot_all_scores = plot_all_scores,
                            method = method,
                            smat = smat,
                            df = df, 
                            knots = knots, 
                            bknots = bknots, 
                            degree = degree, 
                            intercept = intercept, 
                            derivs = derivs, 
                            centerval = centerval, 
                            normalize = normalize, 
                            preH = preH, 
                            sfirst = sfirst, 
                            sparse = sparse,
                            nk = nk,
                            inclx = inclx,
                            knots.only = knots.only,
                            type = type,
                            norm = norm,
                            rpm = rpm,
                            pc = pc,
                            fractied = fractied,
                            bkrange = bkrange,
                            fix_bknots = fix_bknots,
                            bound = bound,
                            userdata = userdata,
                            verbose = verbose)
    
    if(!search_df) {
      knotcnt_suggestion$nknots <- max_nknots
    }
    if(return_df) {
      return(knotcnt_suggestion$df)
    }

    if(knotcnt_suggestion$nknots == 0) {
      stop2c("Search for knotcnt_suggestion using 'get_suggest_knotcount()' 
      resulted in '0' internal knots which must be at least 1.  
             Please change the criterion (e.g., from BIC to AIC) or 
             the model type and retry")
    }
    
    ns_mod <- get_model_by_count(dataset = dataset, 
                             dependent = dependent, 
                             independents = independents,
                             nknots = knotcnt_suggestion$nknots, 
                             method = method,
                             smat = smat,
                             df = df, 
                             knots = knots, 
                             bknots = bknots, 
                             degree = degree, 
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize, 
                             preH = preH, 
                             sfirst = sfirst, 
                             sparse = sparse,
                             nk = nk,
                             inclx = inclx,
                             knots.only = knots.only,
                             type = type,
                             norm = norm,
                             rpm = rpm,
                             pc = pc,
                             fractied = fractied,
                             bkrange = bkrange,
                             fix_bknots = fix_bknots,
                             bound = bound,
                             userdata = userdata)
    
    ns_score <- eval_icr_fn_fun(fun = icr_fn, model = ns_mod, 
                                 dataset = dataset,
                                 cvk = cvk, cviter = cviter, forms = NULL,
                                 cost=function(y,yhat) mean((y-yhat)^2),
                                 verbose = verbose)

    extracted_knots <- get_extract_knots(ns_mod)
    ret <-
      append(ret, list(ns =
                         list(model = ns_mod,
                              score = ns_score,
                              knot_cnt_arg = knotcnt_suggestion$nknots,
                              knot_cnt_distinct = length(extracted_knots$knots),
                              knot_placements = extracted_knots)))
    
    if (verbose) {
      R.utils::printf("%s\n%s: %f\n", ret_desc[["ns"]], score_type, ns_score)
      R.utils::printf("Suggested knot count: %d\n", knotcnt_suggestion$nknots)
      get_print_knots(get_extract_knots(ns_mod))
      R.utils::printf("\n")
    }
    
    # Restricted cubic splines with freely placed knots between the boundaries
    knutar_res <- get_choose_splines(dataset = dataset, 
                                     dependent = dependent, 
                                     independents = independents,
                                     max_nknots = max_nknots, 
                                     icr_fn = icr_fn, 
                                     cost_fn = cost_fn,
                                     method = method,
                                     smat = smat,
                                     df = df, 
                                     knots = knots, 
                                     bknots = bknots, 
                                     degree = degree, 
                                     intercept = intercept, 
                                     derivs = derivs, 
                                     centerval = centerval, 
                                     normalize = normalize, 
                                     preH = preH, 
                                     sfirst = sfirst, 
                                     sparse = sparse,
                                     nk = nk,
                                     inclx = inclx,
                                     knots.only = knots.only,
                                     type = type,
                                     norm = norm,
                                     rpm = rpm,
                                     pc = pc,
                                     fractied = fractied,
                                     bkrange = bkrange,
                                     fix_bknots = fix_bknots,
                                     bound = bound,
                                     userdata = userdata)
    
    ret <-
      append(ret, list(ns_nu =
                         list(model = knutar_res$model,
                              score = knutar_res$score,
                              knot_cnt_distinct= length(knutar_res$knots$knots),
                              knot_placements = knutar_res$knots)))
    
    if (verbose) {
      R.utils::printf("%s\n%s: %f\n",
                      ret_desc[["ns_nu"]], score_type, knutar_res$score)
      get_print_knots(knutar_res$knots)
      R.utils::printf("\n")
    }
    
    
    if (ns_score <= knutar_res$score) {
      ret <- append(ret, list(model = ns_mod, type = "ns", score = ns_score))
    } else {
      ret <- append(ret, list(model = knutar_res$model, type = "ns_nu",
                              score = knutar_res$score))
    }
    
    if (verbose) {
      R.utils::printf("Chosen model type:\n%s\n", ret_desc[[ret$type]])
    }
    
  }) # end of suppressWarnings({
  
  return(ret)
}


##############################################################################
# get_choose_removal 
##############################################################################

#' Finds the best knot to remove from the given model
#'
#' @details The function searches through the knots of the given model to find
#' the inner knot yielding the best resulting score if one knot has to be
#' removed. The selection criterion function can be seen as a cost function,
#' where lower scores are better.
#'
#' @param dataset The data frame
#'
#' @param dependent The dependent variable in the formula
#'
#' @param independents The independent variables in the formula
#'
#' @param cost_fn The function for the selection criterion score (AIC default)
#'
#' @inheritParams get_choose_model
#'
#' @return A named list with 'index' of the chosen knot, 'model', and 'score'
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' removal_index <- get_choose_removal(my_data, y, x, knots, bknots)$index
#' removal_index <- get_choose_removal(my_data, y, x, knots, bknots,
#'  cost_fn = function(model) { return(-summary(model)$fstatistic[[1]]) })$index
#'  
get_choose_removal <- function(dataset,
                           dependent,
                           independents,
                           cost_fn = stats::AIC,
                           cvk = 10, 
                           cviter = 10,
                           method = 'bs',
                           smat = 'ns',
                           df = NULL, 
                           knots = NULL, 
                           bknots = NA, 
                           degree = 3, 
                           intercept = FALSE, 
                           derivs = 0, 
                           centerval = FALSE, 
                           normalize = FALSE, 
                           preH = FALSE, 
                           sfirst = FALSE, 
                           sparse = FALSE,
                           nk = NULL,
                           inclx = TRUE,
                           knots.only = TRUE,
                           type = "ordinary",
                           norm = 2,
                           rpm = NULL,
                           pc = FALSE,
                           fractied = 0.05,
                           bkrange = FALSE,
                           fix_bknots = TRUE,
                           bound = NULL,
                           userdata = NULL,
                           verbose = FALSE) {
  
  if(is.character(dependent))    dependent    <- str2lang(dependent)
  if(is.character(independents)) independents <- str2lang(independents)
  
  independents <- rlang::enquo(independents)
  dependent <- rlang::enquo(dependent)
  
  cost_fn <- check_set_criteria(cost_fn, add_attr = FALSE, verbose = verbose)
  
  
  model_scores <- lapply(seq_along(knots), function(i) {
    mod <- get_model_by_knots(dataset = dataset, 
                          dependent = dependent, 
                          independents = independents,
                          method = method,
                          smat = smat,
                          df = df, 
                          knots = knots[-i],  
                          bknots = bknots, 
                          degree = degree, 
                          intercept = intercept, 
                          derivs = derivs, 
                          centerval = centerval, 
                          normalize = normalize, 
                          preH = preH, 
                          sfirst = sfirst, 
                          sparse = sparse,
                          nk = nk,
                          inclx = inclx,
                          knots.only = knots.only,
                          type = type,
                          norm = norm,
                          rpm = rpm,
                          pc = pc,
                          fractied = fractied,
                          bkrange = bkrange,
                          fix_bknots = fix_bknots,
                          bound = bound,
                          userdata = userdata)
    
    mod_score <- eval_icr_fn_fun(fun = cost_fn, model = mod, 
                                 dataset = dataset,
                                 cvk = cvk, cviter = cviter, forms = NULL,
                                 cost=function(y,yhat) mean((y-yhat)^2),
                                 verbose = verbose)
    
    return(list(model = mod, score = mod_score))
  }) # end model_scores <- lapply(seq_along(knots), function(i) {
  
  
  scores <- unlist(lapply(model_scores, "[[", "score"))
  index <- which.min(scores)
  min_score <- scores[[index]]
  
  return(list(model = model_scores[[index]][["model"]],
              score = min_score, index = index))
}


##############################################################################
# get_choose_splines 
##############################################################################

#' Chooses the best from a set of restricted cubic spline (RCS) models with an
#' inner knot count lower than or equal to a specified maximum number
#'
#' @details The maximum number of inner knots is given as an input argument.
#'
#' @param dataset The data frame
#'
#' @param dependent The dependent variable in the formula
#'
#' @param independents The independent variables in the formula
#'
#' @param max_nknots The maximum number of inner knots wanted
#'
#' @param icr_fn The information criterion function comparing models with
#'   different knot counts (AIC default)
#'
#' @param cost_fn The criterion used to choose which inner knot to remove, used
#'   by the function 'get_choose_removal'. Default is AIC
#'
#' @param initial_nknots The initial high number inner of knots for the
#'   algorithm (default is the value from the 'get_suggest_knotcount'-function)
#'
#' @param diff_better How much lower must the score be to make a higher knot
#'   model be deemed a better model than an alternative lower knot model?
#'
#' @param all_models If TRUE, the function will include all intermediate models
#'   in the results as 'all_models'. Default is FALSE.
#' @param ... Internal
#' 
#' @inheritParams get_choose_model
#'
#' @return The chosen 'model', 'score', 'knots', and 'all_models'
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' my_model <- get_choose_splines(d, y, x, 7)
#' my_model <- get_choose_splines(d, y, x, 7, AIC)
#' 
get_choose_splines <- function(dataset,
                           dependent,
                           independents,
                           max_nknots = 10,
                           ...,
                           icr_fn = stats::AIC,
                           cost_fn = stats::AIC,
                           cvk = 10, 
                           cviter = 10,
                           initial_nknots = -1,
                           diff_better = 0,
                           all_models = FALSE,
                           all_scores = FALSE,
                           plot_all_scores = FALSE,
                           method = 'bs',
                           smat = 'ns',
                           df = NULL, 
                           knots = NULL, 
                           bknots = NA, 
                           degree = 3, 
                           intercept = FALSE, 
                           derivs = 0, 
                           centerval = FALSE, 
                           normalize = FALSE, 
                           preH = FALSE, 
                           sfirst = FALSE, 
                           sparse = FALSE,
                           nk = NULL,
                           inclx = TRUE,
                           knots.only = TRUE,
                           type = "ordinary",
                           norm = 2,
                           rpm = NULL,
                           pc = FALSE,
                           fractied = 0.05,
                           bkrange = FALSE,
                           fix_bknots = TRUE,
                           bound = NULL,
                           userdata = NULL,
                           verbose = FALSE) {
  
  if(is.character(dependent))    dependent    <- str2lang(dependent)
  if(is.character(independents)) independents <- str2lang(independents)
  
  independents <- rlang::enquo(independents)
  dependent    <- rlang::enquo(dependent)
  
  
  if (missing(max_nknots)) max_nknots <- 10

  icr_fn     <- check_set_criteria(icr_fn, add_attr = TRUE, verbose = verbose)
  score_type <- attr(icr_fn, 'score_type')
  attr(icr_fn, 'score_type') <- NULL
  
  cost_fn   <- check_set_criteria(cost_fn, add_attr = FALSE, verbose = verbose)
  
  if (missing(initial_nknots)) initial_nknots <- -1
  if (missing(diff_better)) diff_better <- 0
  if (missing(all_models)) all_models <- FALSE
  if (missing(bknots)) bknots <- NA
  
  if (initial_nknots == -1) {
    initial_nknots <-
      get_suggest_knotcount(dataset = dataset, 
                            dependent = dependent, 
                            independents = independents, 
                            max_nknots  = max_nknots,
                            target_nknots  = max_nknots,
                            icr_fn = icr_fn, 
                            all_scores = all_scores,
                            plot_all_scores = plot_all_scores,
                            method = method,
                            smat = smat,
                            df = df, 
                            knots = knots, 
                            bknots = bknots, 
                            degree = degree, 
                            intercept = intercept, 
                            derivs = derivs, 
                            centerval = centerval, 
                            normalize = normalize, 
                            preH = preH, 
                            sfirst = sfirst, 
                            sparse = sparse,
                            nk = nk,
                            inclx = inclx,
                            knots.only = knots.only,
                            type = type,
                            norm = norm,
                            rpm = rpm,
                            pc = pc,
                            fractied = fractied,
                            bkrange = bkrange,
                            fix_bknots = fix_bknots,
                            bound = bound,
                            userdata = userdata,
                            verbose = verbose)$nknots
  }
 
  upper_model <- get_suggest_splines(dataset = dataset, 
                                     dependent = dependent, 
                                     independents = independents,
                                     max_nknots = max_nknots, 
                                     # target_nknots is same as 
                                     target_nknots = max_nknots,
                                     initial_nknots = initial_nknots,
                                     cost_fn = cost_fn,
                                     all_scores = all_scores,
                                     plot_all_scores = plot_all_scores,
                                     method = method,
                                     smat = smat,
                                     df = df, 
                                     knots = knots, 
                                     bknots = bknots, 
                                     degree = degree, 
                                     intercept = intercept, 
                                     derivs = derivs, 
                                     centerval = centerval, 
                                     normalize = normalize, 
                                     preH = preH, 
                                     sfirst = sfirst, 
                                     sparse = sparse,
                                     nk = nk,
                                     inclx = inclx,
                                     knots.only = knots.only,
                                     type = type,
                                     norm = norm,
                                     rpm = rpm,
                                     pc = pc,
                                     fractied = fractied,
                                     bkrange = bkrange,
                                     fix_bknots = fix_bknots,
                                     bound = bound,
                                     userdata = userdata,
                                     verbose = verbose)
  
  cur_model <- upper_model
  best_model <- cur_model
  
  cur_score <- eval_icr_fn_fun(fun = icr_fn, model = best_model, 
                              dataset = dataset,
                              cvk = cvk, cviter = cviter, forms = NULL,
                              cost=function(y,yhat) mean((y-yhat)^2),
                              verbose = verbose)

  best_score <- cur_score
  best_knots <- get_extract_knots(best_model)
  cur_nknots <- length(best_knots$knots)
  intermediate_models <- list()
  
  
  while (cur_nknots > 0) {
    these_knots <- get_extract_knots(cur_model)
    chosen <- get_choose_removal(dataset = dataset, 
                             dependent = dependent, 
                             independents = independents,
                             cost_fn = cost_fn,
                             method = method,
                             smat = smat,
                             knots = these_knots$knots,
                             bknots = these_knots$Boundary.knots, 
                             degree = degree, 
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize, 
                             preH = preH, 
                             sfirst = sfirst, 
                             sparse = sparse,
                             nk = nk,
                             inclx = inclx,
                             knots.only = knots.only,
                             type = type,
                             norm = norm,
                             rpm = rpm,
                             pc = pc,
                             fractied = fractied,
                             bkrange = bkrange,
                             fix_bknots = fix_bknots,
                             bound = bound,
                             userdata = userdata)
    
    cur_score <- eval_icr_fn_fun(fun = icr_fn, model = chosen$model, 
                                 dataset = dataset,
                                 cvk = cvk, cviter = cviter, forms = NULL,
                                 cost=function(y,yhat) mean((y-yhat)^2),
                                 verbose = verbose)

    if (cur_score <= (best_score + diff_better)) {
      best_model <- chosen$model
      best_score <- cur_score
      best_knots <- get_extract_knots(best_model)
    }
    cur_model <- chosen$model
    cur_nknots <- length(get_extract_knots(cur_model)$knots)
    if (all_models) {
      intermediate_models <- append(intermediate_models, list(cur_model))
    }
  }
  
  return(list(model = best_model, 
              score = best_score, 
              knots = best_knots,
              all_models = intermediate_models))
}


##############################################################################
# get_extract_knots 
##############################################################################

#' Extracts the distinct knot placements of natural spline regression models
#'
#' @param ns_model The restricted cubic spline regression model
#'
#' @return A list with named elements 'knots' (inner) and 'Boundary.knots'
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' get_extract_knots(my_model)
#' 
get_extract_knots <- function(ns_model) {
  knots <- attr(ns_model$model[[2]], "knots")
  knots <- unique(knots)
  bknots <- attr(ns_model$model[[2]], "Boundary.knots")
  fullknots <- c(bknots[1], knots, bknots[2])
  return(list(knots = knots, Boundary.knots = bknots, fullknots = fullknots))
}


##############################################################################
# get_generate_data 
##############################################################################

#' Generates synthetic data
#'
#' Se the examples section for information about the returned column names
#'
#' @param n The number of units in the sample
#'
#' @param x_accr The accuracy (number of decimals) of independent variable
#'
#' @param y_accr The accuracy (number of decimals) of dependent variable
#'
#' @param f_x_dist Function for the distribution of the independent variable
#'
#' @param f_signal The true function for the population means. Takes the x
#'   values as input.
#'
#' @param f_noise The noise function, i.e., variance distribution(s) around the
#'   population means for y, which can be heteroscedastic.
#'
#' @return A dataframe with generated data
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' The returned data frame has the following columns:
#' ID               # IDs
#' Independent      # The measured X values (rounded according to accuracy)
#' Dependent        # The measured Y values (rounded according to accuracy)
#' IndependentRaw   # The raw X values (not rounded)
#' DependentRaw     # The raw Y values (not rounded)
#' SignalRaw        # The raw signal part of the population mean curve for Y
#' Noise            # The raw noise (variance) for the DependentRaw values
#' SignalMeasured   # The signal rounded according to accuracy
#' 
get_generate_data <- function(n, x_accr, y_accr, f_x_dist,
                          f_signal, f_noise) {
  ids <- 1:n
  xs_raw <- f_x_dist(n)
  ys_signal <- f_signal(xs_raw)
  ys_noise <- f_noise(xs_raw)
  ys_raw <- ys_signal + ys_noise
  
  xs_measured <- xs_raw
  if (!missing(x_accr) && !is.null(x_accr) && x_accr > -1) {
    xs_measured <- round(xs_raw, x_accr)
  }
  
  ys_measured <- ys_raw
  if (!missing(y_accr) && !is.null(y_accr) && y_accr > -1) {
    ys_measured <- round(ys_raw, y_accr)
  }
  
  xs_measured_signal <- f_signal(xs_measured)
  
  return(data.frame(
    ID = ids,
    Independent = xs_measured,
    Dependent = ys_measured,
    IndependentRaw = xs_raw,
    DependentRaw = ys_raw,
    SignalRaw = ys_signal,
    Noise = ys_noise,
    SignalMeasured = xs_measured_signal
  ))
  
}


##############################################################################
# get_plot_model 
##############################################################################

#' Plots the (glm) model with data, curve, confidence bands, etc.
#'
#' @param dataset The data frame
#'
#' @param mod The model
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' get_plot_model(d, y, x, my_mod)
#' 
get_plot_model <- function(dataset, 
                           model,
                           x,
                           y,
                           fullknots = NULL,
                           title = NULL,
                           subtitle = NULL,
                           verbose = FALSE) {
  fit_link <- NULL;
  se_link  <- NULL;
  lwr      <- NULL;
  upr      <- NULL;
  pred     <- NULL;
  if(is.symbol(x)) x <- deparse(x)
  if(is.symbol(y)) y <- deparse(y)
  dependent <- y
  independent <- x
  if(is.character(dependent))    dependent    <- str2lang(dependent)
  if(is.character(independent)) independent <- str2lang(independent)
  
  
  if(is.null(fullknots)) {
    knots <- get_extract_knots(model)
  } else if(!is.null(fullknots)) {
    knots <- list()
    knots$knots          <- checkgetiknotsbknots(fullknots, 'iknots')
    knots$Boundary.knots <- checkgetiknotsbknots(fullknots, 'bknots')
  }
  
  # The confidence bands has to be in the scale of the link function,
  # and we need to get it's inverse for plotting
  
  ilink <- stats::family(model)$linkinv
  
  d <- dataset
  d <- d %>% data.frame() %>% 
    dplyr::mutate(pred = stats::predict(model)) %>%
    dplyr::bind_cols(stats::setNames(dplyr::as_tibble(
      stats::predict(model, newdata = d, se.fit = TRUE)[1:2]),
      c("fit_link", "se_link"))) %>%
    dplyr::mutate(fit_resp = ilink(fit_link),
                  upr = ilink(fit_link + (2 * se_link)),
                  lwr = ilink(fit_link - (2 * se_link)))
  
  fig <- ggplot2::ggplot(d,
                         ggplot2::aes(x = {{ independent }}, 
                                      y = {{ dependent }})) +
    ggplot2::geom_point() +
    ggplot2::geom_ribbon(data = d,
                         ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = pred),
                       color = "blue", linewidth = 1)
  
  
  fig <- fig + 
    ggplot2::geom_vline(xintercept = knots$knots, linetype = "dashed")
  fig <- fig + 
    ggplot2::geom_vline(xintercept = knots$Boundary.knots, linetype = "solid")
  
  fig <- fig + ggplot2::labs(x = x, y = y)
  
  if(!is.null(title)) {
    fig <- fig + ggplot2::labs(title = title)
  }
  
  if(!is.null(subtitle)) {
    fig <- fig + ggplot2::labs(subtitle = subtitle)
  }
  
  return(fig)
}



##############################################################################
# get_plot_rstandard 
##############################################################################

#' Plots the standardized residuals for the model
#'
#' @param dataset The data frame
#'
#' @param mod The model
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' get_plot_rstandard(d, my_mod)
#' 
get_plot_rstandard <- function(dataset, 
                           model,
                           x = NULL,
                           y = NULL,
                           fullknots = NULL,
                           title = NULL,
                           subtitle = NULL,
                           verbose = FALSE) {
  
  
  if(is.null(fullknots)) {
    knots <- get_extract_knots(model)
  } else if(!is.null(fullknots)) {
    knots <- list()
    knots$knots          <- checkgetiknotsbknots(fullknots, 'iknots')
    knots$Boundary.knots <- checkgetiknotsbknots(fullknots, 'bknots')
  }
  
  fig <- ggplot2::ggplot(dataset, 
                         ggplot2::aes(x = stats::predict(model, dataset),
                                        y = stats::rstandard(model))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x)
  
  if(!is.null(x)) {
    fig <- fig + ggplot2::labs(x = x)
  }
  
  if(!is.null(y)) {
    fig <- fig + ggplot2::labs(y = y)
  }
  
  if(!is.null(title)) {
    fig <- fig + ggplot2::labs(title = title)
  }
  
  if(!is.null(subtitle)) {
    fig <- fig + ggplot2::labs(subtitle = subtitle)
  }
  
  return(fig)
}



##############################################################################
# get_print_knots 
##############################################################################

#' Print knot placements
#'
#' @param knot_placements A list with two named lists, 'knots' and
#'   'Boundary.knots'
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' get_print_knots(knot_placements)
#' 
get_print_knots <- function(knot_placements) {
  R.utils::printf("Inner knots count: %d\n", length(knot_placements$knots))
  knots_str <- paste0("[", paste0(knot_placements$knots, collapse = ", "), "]")
  boundary_str <-
    paste0("[", paste0(knot_placements$Boundary.knots, collapse = ", "), "]")
  R.utils::printf("Inner knots: %s\nBoundary knots: %s\n", knots_str,
                  boundary_str)
}


##############################################################################
# restricted_cubic_splines
##############################################################################

#' Utility function for restricted cubic splines that wraps the 'glm' and 'ns'
#' functions
#'
#' @details Creates a restricted cubic splines regression model given the wanted
#' number of inner knots
#'
#' @param dataset The data frame
#'
#' @param dependent The dependent variable in the formula
#'
#' @param independents The independent variables in the formula
#'
#' @param nknots The requested number of inner knots, excluding the boundary
#'   knots
#' 
#' @inheritParams get_choose_model
#'
#' @return The regression model
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' my_model <- create_model(my_data, y, x, 7)
#' 
get_model_by_count <- function(dataset, 
                           dependent, 
                           independents, 
                           nknots,
                           method = 'bs',
                           smat = 'ns',
                           df = NULL, 
                           knots = NULL, 
                           bknots = NA, 
                           degree = 3, 
                           intercept = FALSE, 
                           derivs = 0, 
                           centerval = FALSE, 
                           normalize = FALSE, 
                           preH = FALSE, 
                           sfirst = FALSE, 
                           sparse = FALSE,
                           nk = NULL,
                           inclx = TRUE,
                           knots.only = TRUE,
                           type = "ordinary",
                           norm = 2,
                           rpm = NULL,
                           pc = FALSE,
                           fractied = 0.05,
                           bkrange = FALSE,
                           fix_bknots = TRUE,
                           bound = NULL,
                           userdata = NULL,
                           verbose = FALSE) {
  
  if(is.character(dependent))    dependent    <- str2lang(dependent)
  if(is.character(independents)) independents <- str2lang(independents)
  
  independents <- rlang::enquo(independents)
  dependent <- rlang::enquo(dependent)
  
  independents_str <- sub("~", "", deparse(independents))
  
  if (missing(bknots)) {
    bknots <- NA
  } else if(is.null(bknots)) {
    bknots <- NA
  }
  
  
  if(!is_emptyx(bknots)) {
    if (length(bknots) != 2) {
      stop("'bknots' must be a length of 2") 
    }
  }
  
  if(bkrange) {
    bknots <- range(dataset[[independents_str]])
  }
  
  fullknots <- get_knost_from_df(x = dataset[[independents_str]], 
                                 df = nknots + 1, 
                                 method = method,
                                 smat = smat,
                                 knots = knots, 
                                 bknots = bknots, 
                                 degree = degree, 
                                 intercept = intercept, 
                                 derivs = derivs, 
                                 centerval = centerval, 
                                 normalize = normalize, 
                                 preH = preH, 
                                 sfirst = sfirst, 
                                 sparse = sparse,
                                 nk = nk,
                                 inclx = inclx,
                                 knots.only = knots.only,
                                 type = type,
                                 norm = norm,
                                 rpm = rpm,
                                 pc = pc,
                                 fractied = fractied,
                                 bkrange = bkrange,
                                 fix_bknots = fix_bknots,
                                 bound = bound,
                                 userdata = userdata) 
  
  
  knots <- checkgetiknotsbknots(fullknots, 'iknots')
  bknots <- checkgetiknotsbknots(fullknots, 'bknots')
  
  bknots_str <- paste0("c(", paste0(bknots, collapse = ", "), ")")
  
  myform <- get_form(x = independents_str,
                     method = method,
                          smat = smat, 
                          knots = knots,
                          bknots = bknots_str,
                          df = df,
                          degree = degree,
                          intercept = intercept, 
                          derivs = derivs,
                          centerval = centerval, 
                          normalize = normalize, 
                          preH = preH, 
                          sfirst = sfirst, 
                          sparse = sparse,
                          verbose = verbose)
  
  model_formula <- paste0(rlang::as_label(dependent), " ~ ", myform)
  
  model_formula_str <- model_formula
  
  model_formula <- ept(model_formula_str, envir = environment(myform))
  
  # Somehow, direct model fails for GS_bsp etc
  if(is.null(environment(myform))) {
    model_formula_str <- gsub("GS_", "bsitar:::GS_", model_formula_str)
    ns_model_str <- paste0("stats::glm", "(", 
                           model_formula_str, ",", 
                           "data = dataset", 
                           ")" 
    )
    ns_model <- ept(ns_model_str)
  } else {
    ns_model <- stats::glm(model_formula, data = dataset)
  }

  attr(ns_model$model[[2]], 'knots')          <- knots
  attr(ns_model$model[[2]], 'Boundary.knots') <- bknots
  
  return(ns_model)
}


#' Utility function for restricted cubic splines that wraps the 'glm' and 'ns'
#' functions
#'
#' Creates a restricted cubic spline regression model from the given placements
#' of the inner knots and boundary knots
#'
#' @param dataset The data frame
#' @param dependent The dependent variable in the formula
#' @param independents The independent variables in the formula
#' 
#' @inheritParams get_choose_model
#'
#' @return The regression model
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' my_model <- get_model_by_knots(my_data, y, x, c(0.1, 0.2), c(0.0, 0.3))
#' 
get_model_by_knots <- function(dataset,
                           dependent,
                           independents,
                           method = 'bs',
                           smat = 'ns',
                           df = NULL, 
                           knots = NULL, 
                           bknots = NA, 
                           degree = 3, 
                           intercept = FALSE, 
                           derivs = 0, 
                           centerval = FALSE, 
                           normalize = FALSE, 
                           preH = FALSE, 
                           sfirst = FALSE, 
                           sparse = FALSE,
                           nk = NULL,
                           inclx = TRUE,
                           knots.only = TRUE,
                           type = "ordinary",
                           norm = 2,
                           rpm = NULL,
                           pc = FALSE,
                           fractied = 0.05,
                           bkrange = FALSE,
                           fix_bknots = TRUE,
                           bound = NULL,
                           userdata = NULL,
                           verbose = FALSE) {
  
  if(is.character(dependent))    dependent    <- str2lang(dependent)
  if(is.character(independents)) independents <- str2lang(independents)
  
  independents <- rlang::enquo(independents)
  dependent <- rlang::enquo(dependent)
 
  independents_str <- sub("~", "", deparse(independents))
  
  knots_str <- paste0(
    "c(", paste0(knots, collapse = ", "), ")")
  bknots_str <- paste0(
    "c(", paste0(bknots, collapse = ", "), ")")
  
  
  myform <- get_form(x = independents_str,
                     method = method,
                          smat = smat, 
                          df = NULL,
                          knots = knots,
                          bknots = bknots_str)
  
  model_formula <- paste0(rlang::as_label(dependent), " ~ ", myform)
  
  model_formula_str <- model_formula
  
  model_formula <- ept(model_formula_str, envir = environment(myform))

  # Somehow, direct model fails for GS_bsp etc
  if(is.null(environment(myform))) {
    model_formula_str <- gsub("GS_", "bsitar:::GS_", model_formula_str)
    ns_model_str <- paste0("stats::glm", "(", 
                           model_formula_str, ",", 
                           "data = dataset", 
                           ")" 
    )
    ns_model <- ept(ns_model_str)
  } else {
    ns_model <- stats::glm(model_formula, data = dataset)
  }
  
  
  attr(ns_model$model[[2]], 'knots')          <- knots
  attr(ns_model$model[[2]], 'Boundary.knots') <- bknots
  
  return(ns_model)
}



##############################################################################
# get_suggest_knotcount 
##############################################################################

#' Finds the number of restricted cubic spline inner knots that gives the lowest
#' score for an information criterion and a data set.
#'
#' @param dataset The data frame
#' @param dependent The dependent variable in the formula
#' @param independents The independent variable(s) in the formula
#' @param icr_fn The information criterion function. Defaults to AIC
#' @param max_nknots maximum number of knots
#' @param ... Internal
#' @param all_scores If TRUE, all scores are returned in a list 'all_scores'
#' 
#' @inheritParams get_choose_model
#'
#' @return A list with named elements 'nknots', 'score', and 'all_scores'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' get_suggest_knotcount(d, nwsize, age_dec)
#' 
get_suggest_knotcount <- function(dataset,
                              dependent,
                              independents,
                              max_nknots = -1,
                              initial_nknots = -1,
                              ...,
                              icr_fn = stats::AIC,
                              cvk = 10, 
                              cviter = 10,
                              all_scores = FALSE,
                              plot_all_scores = FALSE,
                              method = 'bs',
                              smat = 'ns',
                              df = NULL, 
                              knots = NULL, 
                              bknots = NA, 
                              degree = 3, 
                              intercept = FALSE, 
                              derivs = 0, 
                              centerval = FALSE, 
                              normalize = FALSE, 
                              preH = FALSE, 
                              sfirst = FALSE, 
                              sparse = FALSE,
                              nk = NULL,
                              inclx = TRUE,
                              knots.only = TRUE,
                              type = "ordinary",
                              norm = 2,
                              rpm = NULL,
                              pc = FALSE,
                              fractied = 0.05,
                              bkrange = FALSE,
                              fix_bknots = TRUE,
                              bound = NULL,
                              userdata = NULL,
                              verbose = FALSE) {
  
 
  if(is.character(dependent))    dependent    <- str2lang(dependent)
  if(is.character(independents)) independents <- str2lang(independents)
 
  dependent    <- rlang::enquo(dependent)
  independents <- rlang::enquo(independents)
  
  if (missing(max_nknots) || max_nknots == -1) {
    max_nknots <- min(50, nrow(dataset) %/% 2)
  }
  
  icr_fn     <- check_set_criteria(icr_fn, add_attr = TRUE, verbose = verbose)
  score_type <- attr(icr_fn, 'score_type')
  attr(icr_fn, 'score_type') <- NULL
  
  
  if (missing(all_scores)) all_scores <- FALSE
  if (missing(bknots)) bknots <- NA
  
  if(!all_scores) {
    plot_all_scores <- FALSE
  }
  
  min_icr <- Inf
  min_ndf <- Inf
  
  n_knots <- list()
  scores <- list()
  
  independents_str <- sub("~", "", deparse(independents))
  
  dependent_str <- sub("~", "", deparse(dependent))
  
  
  
  if(!is_emptyx(bknots)) {
    if (length(bknots) != 2) {
      stop("'bknots' must be a length of 2") 
    }
  }
  
  if(bkrange) {
    bknots <- range(dataset[[independents_str]])
  }
  
  
  bknots_str <- paste0(
    ", Boundary.knots = c(", paste0(bknots, collapse = ", "), ")"
  )

 
  
  consecutive_non_convergance <- 0
  
  for (i in 1:(max_nknots + 1)) {
    # Because of a change between R 4.2.0 and 4.3.0, quantiles coinciding with
    # boundary knots must now be removed.
    #
    # This caused a bit of rewriting...
    #
    # See: https://stat.ethz.ch/pipermail/r-announce/2023/000691.html
    # "bs() and ns() in the (typical) case of automatic knot construction, when
    # some of the supposedly inner knots coincide
    # with boundary knots, now moves them inside (with a warning),
    # building on PR#18442 by Ben Bolker."
    
    
    subset_data <- dataset %>%
      dplyr::filter(
        !!independents >= bknots[[1]],
        !!independents <= bknots[[2]]
      )
    

    n <- i - 1
    knots <- c()
  
    if (n > 0) {
      quantiles <- unique(quantile(subset_data[[rlang::as_label(independents)]],
                                   probs = c(0, seq(1 / n, 1, by = 1 / n))
      ))
      
      if (quantiles[[1]] == bknots[[1]]) {
        quantiles <- quantiles[-1]
      }
      
      if (!length(quantiles) == 0 &&
          quantiles[[length(quantiles)]] == bknots[[2]]) {
        quantiles <- quantiles[-length(quantiles)]
      }
      
      knots <- quantiles
    }
   
    model_formula_str <- NULL
   
    myform <- get_form(x = independents_str,
                       method = method,
                       smat = smat, 
                       df = NULL,
                       knots = knots, 
                       bknots = bknots_str)
    
    model_formula <- paste0(rlang::as_label(dependent), " ~ ", myform)
   
    model_formula_str <- model_formula
    
    model_formula <- ept(model_formula_str, envir = environment(myform))
    
    mod_spline <- NULL
    
    # Somehow, direct model fails for GS_bsp etc
    if(is.null(environment(myform))) {
      model_formula_str <- gsub("GS_", "bsitar:::GS_", model_formula_str)
      ns_model_str <- paste0("stats::glm", "(", 
                             model_formula_str, ",", 
                             "data = dataset", 
                             ")" 
      )
      try(mod_spline <-  ept(ns_model_str))
    } else {
      try(mod_spline <- stats::glm(model_formula, data = dataset))
    }
    

    if (!is.null(mod_spline) && mod_spline$converged) {
      consecutive_non_convergance <- 0
    } else {
      consecutive_non_convergance <- consecutive_non_convergance + 1
    }
    
    
   
    if (consecutive_non_convergance == 0) {
      icr_score <- eval_icr_fn_fun(fun = icr_fn, model = mod_spline, 
                                   dataset = dataset,
                                   cvk = cvk, cviter = cviter, forms = NULL,
                                   cost=function(y,yhat) mean((y-yhat)^2),
                                   verbose = verbose)
      
      if (all_scores) {
        scores <- append(scores, icr_score)
        n_knots <- append(n_knots, i - 1)
      }
      
      if (icr_score < min_icr) {
        min_icr <- icr_score
        min_ndf <- i 
      }
    } else if (consecutive_non_convergance >= 3) {
      warning(paste(
        "Models failed to converge three consecutive times,",
        "will not assess any higher knot counts."
      ))
      break
    }
  }
  
  
  if (all_scores) {
    all_scores_df <- cbind(n_knots %>% unlist(), scores %>% unlist()) %>% 
      data.frame() %>% setNames(c('knot', 'score'))
  } else {
    all_scores_df <- NULL
  }
  
  
  if (!all_scores) {
    if(plot_all_scores) {
      plot_all_scores <- FALSE
      if(verbose) message2c("setting 'plot_all_scores = FALSE' because of 
                          'all_scores =  FALSE'")
    }
  }
  
  
  if(plot_all_scores) {
    ylab_str <- paste0(score_type, " ", "(", dependent_str, ")")
    xlab_str <- paste0("Degrees of Freedom (full knots -1)",
                       "\n",
                       "full knots =  knots + boundary knots")
    getallscores  <- all_scores_df$score
    kgetallscores <- 1:length(getallscores)
    plot(kgetallscores, getallscores, type="o", pch=16, col="black", 
         xlab=xlab_str, 
         main = "",
         sub = "",
         ylab = ylab_str)
    graphics::points(kgetallscores[which.min(getallscores)], 
                     min(getallscores), pch=16, col="red")
  }
  
  # nknots = internal knots
  return(list(
    df = min_ndf,
    nknots = min_ndf - 1, score = min_icr, 
    method = method, smat = smat,
    all_scores_df = all_scores_df,
    all_scores = list(scores = scores, n_knots = n_knots)
  ))
  
}


##############################################################################
# get_suggest_splines 
##############################################################################

#' Suggests a restricted cubic splines regression model with inner knot
#' placements that can be non-uniform with respect to quantiles and widths
#'
#' The target number of inner knots for the model is given as an input argument.
#' The algorithm starts with a model with a high number of knots and
#' systematically removes inner knots until the target number of inner knots is
#' reached. The initial number of inner knots, before starting to remove knots,
#' can be given as an argument as well, but  defaults to the suggested number
#' obtained from the function 'get_suggest_knotcount'.
#'
#' @param dataset The data frame
#' @param dependent The dependent variable in the formula
#' @param independents The independent variables in the formula
#' @param target_nknots The target and maximum number of inner knots for the
#'   model
#' @param initial_nknots The number of inner knots initially, defaults to the
#'   result from the function 'get_suggest_knotcount'
#' @param cost_fn The function for the selection criterion score (AIC default)
#'   used to compare which inner knot should be removed, passed to
#'   get_choose_removal
#' @param all_knots If TRUE, then knots for all intermediate models will be
#' @param ... Internal
#' 
#' @inheritParams get_choose_model
#'
#' @return Returns the suggested natural splines model, or if the 'all_knots'
#'   argument was TRUE, then a list with named elements 'model', 'all_knots',
#'   and 'Boundary.knots' is returned.
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' my_model <- get_suggest_splines(d, y, x, 4)
#' my_model <- get_suggest_splines(d, y, x, 4,initial_nknots=100, cost_fn = AIC)
#' 
get_suggest_splines <- function(dataset,
                            dependent,
                            independents,
                            target_nknots,
                            ...,
                            initial_nknots = -1,
                            cost_fn = stats::AIC,
                            icr_fn = stats::AIC,
                            cvk = 10, 
                            cviter = 10,
                            all_knots = FALSE,
                            all_scores = FALSE,
                            plot_all_scores = FALSE,
                            method = 'bs',
                            smat = 'ns',
                            df = NULL, 
                            knots = NULL, 
                            bknots = NA, 
                            degree = 3, 
                            intercept = FALSE, 
                            derivs = 0, 
                            centerval = FALSE, 
                            normalize = FALSE, 
                            preH = FALSE, 
                            sfirst = FALSE, 
                            sparse = FALSE,
                            nk = NULL,
                            inclx = TRUE,
                            knots.only = TRUE,
                            type = "ordinary",
                            norm = 2,
                            rpm = NULL,
                            pc = FALSE,
                            fractied = 0.05,
                            bkrange = FALSE,
                            fix_bknots = TRUE,
                            bound = NULL,
                            userdata = NULL,
                            verbose = FALSE) {
  
  if(is.character(dependent))    dependent    <- str2lang(dependent)
  if(is.character(independents)) independents <- str2lang(independents)
  
  independents <- rlang::enquo(independents)
  dependent    <- rlang::enquo(dependent)
  
  dependent_str <- sub("~", "", deparse(dependent))
  independents_str <- sub("~", "", deparse(independents))
  
  
  if (missing(bknots)) bknots <- NA
  

  if(!is_emptyx(bknots)) {
    if (length(bknots) != 2) {
     stop("'bknots' must be a length of 2") 
    }
  }
  
  if(bkrange) {
    bknots <- range(dataset[[independents_str]])
  }
  
  
  bknots_org <- bknots
  
  if(fix_bknots) {
    if(any(is.na(bknots))) {
      stop2c("you have set 'fix_bknots = TRUE' but the argument 
            'bknots' is missing")
    }
    bknots <- bknots_org
  }
  
  if(smat == 'rcs') {
    if (initial_nknots == -1) {
     # initial_nknots <- 1
    }
  }
  
  
  ###########################################################
  
  if(method == 'rs') {
    if(all_scores | all_knots) {
      stop2c("arguments all_scores and all_knots are 
           not supported for method = rs")
    }
    bknots_str <- paste0(
      ", Boundary.knots = c(", paste0(bknots, collapse = ", "), ")")
    myform <- get_form(x = independents_str,
                       method = method,
                       smat = smat, 
                       df = NULL,
                       knots = knots, 
                       bknots = bknots_str)
    
    model_formula <- paste0(rlang::as_label(dependent), " ~ ", myform)
    
    model_formula_str <- model_formula
    model_formula     <- ept(model_formula_str, envir = environment(myform))
    
    # Somehow, direct model fails for GS_bsp etc
    if(is.null(environment(myform))) {
      model_formula_str <- gsub("GS_", "bsitar:::GS_", model_formula_str)
      ns_model_str <- paste0("stats::glm", "(", 
                             model_formula_str, ",", 
                             "data = dataset", 
                             ")" )
    } else {
      model_formula_str <- model_formula_str
      ns_model_str <- paste0("stats::glm", "(", 
                             model_formula_str, ",", 
                             "data = dataset", 
                             ")" )
    }
    
    result <- design_sigmoid_knots(dataset = dataset,
                                   x = independents_str, 
                                   y = dependent_str, 
                                   max_knots = initial_nknots, 
                                   min_knots = target_nknots + 1, 
                                   criterion = cost_fn,
                                   smat = smat, 
                                   ns_model_str = ns_model_str,
                                   model_formula = model_formula,
                                   fix_bknots = fix_bknots,
                                   bknots = bknots)
    final_mod <- result$model
    attr(final_mod$model[[2]], 'knots')          <- result$knots
    attr(final_mod$model[[2]], 'Boundary.knots') <- result$Boundary.knots
    return(final_mod)
  }
  
  
  ###########################################################

  if (initial_nknots == -1) {
    initial_nknots <-
      get_suggest_knotcount(dataset = dataset, 
                            dependent = dependent, 
                            independents = independents, 
                            # max_nknots  = max_nknots,
                            icr_fn = icr_fn, 
                            all_scores = all_scores,
                            plot_all_scores = plot_all_scores,
                            method = method,
                            smat = smat,
                            df = df, 
                            knots = knots, 
                            bknots = bknots, 
                            degree = degree, 
                            intercept = intercept, 
                            derivs = derivs, 
                            centerval = centerval, 
                            normalize = normalize, 
                            preH = preH, 
                            sfirst = sfirst, 
                            sparse = sparse,
                            nk = nk,
                            inclx = inclx,
                            knots.only = knots.only,
                            type = type,
                            norm = norm,
                            rpm = rpm,
                            pc = pc,
                            fractied = fractied,
                            bkrange = bkrange,
                            fix_bknots = fix_bknots,
                            bound = bound,
                            userdata = userdata,
                            verbose = verbose)$nknots
  }
  
  
  if(initial_nknots == 0) {
    stop2c("Search for initial_nknots using 'get_suggest_knotcount()' 
      resulted in '0' internal knots which must be at least 1.  
             Please change the criterion (e.g., from BIC to AIC) or 
             the model type and retry")
  }
  
  cost_fn   <- check_set_criteria(cost_fn, add_attr = FALSE, verbose = verbose)
  
  
  if (missing(all_knots)) all_knots <- FALSE

  ns_model <-
    get_model_by_count(dataset = dataset, 
                   dependent = dependent, 
                   independents = independents, 
                   nknots = initial_nknots,
                   method = method,
                   smat = smat,
                   df = df, 
                   knots = knots, 
                   bknots = bknots, 
                   degree = degree, 
                   intercept = intercept, 
                   derivs = derivs, 
                   centerval = centerval, 
                   normalize = normalize, 
                   preH = preH, 
                   sfirst = sfirst, 
                   sparse = sparse,
                   nk = nk,
                   inclx = inclx,
                   knots.only = knots.only,
                   type = type,
                   norm = norm,
                   rpm = rpm,
                   pc = pc,
                   fractied = fractied,
                   bkrange = bkrange,
                   fix_bknots = fix_bknots,
                   bound = bound,
                   userdata = userdata)
  
  knots <- get_extract_knots(ns_model)
  
  intermediate_knots <- list()
  
  # Initialize the variables that will hold the final knot placements
  final_knots    <- knots$knots
  bknots <- knots$Boundary.knots
  
  if(fix_bknots) {
    bknots <- bknots_org
  }
 
  if (all_knots) {
    intermediate_knots <- append(intermediate_knots, list(final_knots))
  }
  
  
  if (length(knots$knots) > target_nknots) {
    for (i in 1:(length(knots$knots) - target_nknots)) {
      rm_index <- get_choose_removal(dataset = dataset, 
                                 dependent = dependent, 
                                 independents = independents,
                                 cost_fn = cost_fn,
                                 method = method,
                                 smat = smat,
                                 knots = final_knots, 
                                 bknots = bknots, 
                                 degree = degree, 
                                 intercept = intercept, 
                                 derivs = derivs, 
                                 centerval = centerval, 
                                 normalize = normalize, 
                                 preH = preH, 
                                 sfirst = sfirst, 
                                 sparse = sparse,
                                 nk = nk,
                                 inclx = inclx,
                                 knots.only = knots.only,
                                 type = type,
                                 norm = norm,
                                 rpm = rpm,
                                 pc = pc,
                                 fractied = fractied,
                                 bkrange = bkrange,
                                 fix_bknots = fix_bknots,
                                 bound = bound,
                                 userdata = userdata)$index
      
      final_knots <- final_knots[-rm_index]
      
      if (all_knots) {
        intermediate_knots <- append(intermediate_knots, list(final_knots))
      }
    }
  }
 
  
  final_mod <- get_model_by_knots(dataset = dataset, 
                              dependent = dependent, 
                              independents = independents,
                              method = method,
                              smat = smat,
                              df = df, 
                              knots = final_knots, 
                              bknots = bknots, 
                              degree = degree, 
                              intercept = intercept, 
                              derivs = derivs, 
                              centerval = centerval, 
                              normalize = normalize, 
                              preH = preH, 
                              sfirst = sfirst, 
                              sparse = sparse,
                              nk = nk,
                              inclx = inclx,
                              knots.only = knots.only,
                              type = type,
                              norm = norm,
                              rpm = rpm,
                              pc = pc,
                              fractied = fractied,
                              bkrange = bkrange,
                              fix_bknots = fix_bknots,
                              bound = bound,
                              userdata = userdata)
  
  if (all_knots) {
    return(list(model = final_mod, all_knots = intermediate_knots,
                Boundary.knots = bknots))
  } else {
    return(final_mod)
  }
}



##############################################################################
# get_create_figure
##############################################################################

# https://github.com/jo-inge-arnes/knutar-experiments

#' Title
#'
#' @param d Internal
#' @param mod Internal
#' @param title Internal
#' @param x Internal
#' @param y Internal
#'
#' @returns A \code{ggplot} object
#' @keywords internal
#' @noRd
#'
get_create_figure <- function(dataset, 
                              model, 
                              x,
                              y,
                              fullknots = NULL,
                              title = NULL, 
                              subtitle = NULL,
                              verbose = FALSE) {
  
  d   <- dataset
  mod <- model
  if(is.symbol(x)) x <- deparse(x)
  if(is.symbol(y)) y <- deparse(y)
  if(is.character(x)) x <- str2lang(x)
  if(is.character(y)) y <- str2lang(y)
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  x <- sub("~", "", deparse(x))
  y <- sub("~", "", deparse(y))
  fig <- ggplot2::ggplot()
  fig <- fig + ggplot2::theme_bw()
  fig <- fig + ggplot2::geom_point(data = d, 
                                   ggplot2::aes(.data[[x]], .data[[y]]), 
                                   shape = 1,
                                   color = "gray50")
 
  
  fig <- fig + ggplot2::xlab(scales::parse_format()("'Predictor'~X"))
  fig <- fig + ggplot2::ylab(scales::parse_format()("'Response'~Y"))
  
  
  if(is.null(fullknots)) {
    knots <- get_extract_knots(mod)
  } else if(!is.null(fullknots)) {
    knots <- list()
    knots$knots          <- checkgetiknotsbknots(fullknots, 'iknots')
    knots$Boundary.knots <- checkgetiknotsbknots(fullknots, 'bknots')
  }
  
  fig <- fig + 
    ggplot2::geom_vline(xintercept = knots$knots, linetype = "dashed")
  fig <- fig + 
    ggplot2::geom_vline(xintercept = knots$Boundary.knots, linetype = "solid")
  
  fig <- fig + 
    ggplot2::geom_line(ggplot2::aes(x = mod$data[[x]], y = mod$fitted.values),
                         linetype = "solid", color = "black", linewidth = 0.5)
  
  # fig <- fig + force_panelsizes(rows = unit(3.5, "in"), cols = unit(7, "in"))
  
  if(!is.null(x)) {
    fig <- fig + ggplot2::labs(x = x)
  }
  
  if(!is.null(y)) {
    fig <- fig + ggplot2::labs(y = y)
  }
  
  if(!is.null(title)) {
    fig <- fig + ggplot2::labs(title = title)
  }
  
  if(!is.null(subtitle)) {
    fig <- fig + ggplot2::labs(subtitle = subtitle)
  }
  
  return(fig)
}



##############################################################################
# get_print_return_obj
##############################################################################

# https://github.com/jo-inge-arnes/knutar-experiments

#' Title
#'
#' @param knots_old Internal
#' @param model_new Internal
#' @param knots_new Internal
#' @param what Internal
#' @param print Internal
#' @param return Internal
#' @param verbose Internal
#'
#' @returns A \code{ggplot} object
#' @keywords internal
#' @noRd
#'
get_model_by_knots_wrapper <- function(list_arg, knots) {
  list_arg[['knots']]  <-  checkgetiknotsbknots(knots, 'iknots')
  list_arg[['bknots']] <-  checkgetiknotsbknots(knots, 'bknots')
  
  remove_these_args <- c('x', 
                         'max_nknots', 
                         'target_nknots', 
                         'initial_nknots',
                         'print', 
                         'return', 
                         'select', 
                         'return_df', 
                         'search_df', 
                         'cost_fn', 
                         'icr_fn', 
                         'kspace',
                         'cvk',
                         'cviter',
                         'plot_all_scores',
                         'all_scores',
                         'all_knots')
  for (i in remove_these_args) {
    list_arg[i] <- NULL
  }
  
  model <- do.call(get_model_by_knots, list_arg)
  return(model)
}


##############################################################################
# get_print_return_obj
##############################################################################

# https://github.com/jo-inge-arnes/knutar-experiments

#' Title
#'
#' @param knots_old Internal
#' @param model_new Internal
#' @param knots_new Internal
#' @param select Internal
#' @param kspace Internal
#' @param what Internal
#' @param print Internal
#' @param return Internal
#' @param verbose Internal
#'
#' @returns A \code{ggplot} object
#' @keywords internal
#' @noRd
#'
get_print_return_obj <- function(knots = NULL, 
                                 model = NULL,  
                                 model_new = NULL,  
                                 knots_new = NULL, 
                                 nys,
                                 xsi, 
                                 ysi,
                                 select,
                                 kspace,
                                 what,
                                 print = FALSE, 
                                 return = FALSE,
                                 list_arg = NULL,
                                 verbose = FALSE) {
  

  knots_old              <- knots
  knots_selection_what   <- what 
  knots_selection_print  <- print 
  knots_selection_return <- return
  
  model_old     <- get_model_by_knots_wrapper(list_arg, knots_old) 
  knots_old_str <-  paste(knots_old, collapse = " ")
  
  # y -axis is labelled by x, so no need of substile
  if(nys > 1) {
    # setsubtitle <- ysi
    setsubtitle <- NULL
  } else {
    setsubtitle <- NULL
  }
  
  
  ###########################################################################
  ###########################################################################
  if(is.null(model_new)) {
    # if(is_emptyx(what)) what <- 'plot1'
    knots_selection_what_choices <- c('knots', 
                                      'plot1', 
                                      'plot2', 
                                      'plot3')
    
    if(!what %in% knots_selection_what_choices) {
      stop("knots_selection_what must be one of the following: ",
           "\n ",
           collapse_comma(knots_selection_what_choices))
    }
    
    if(what == 'plot1') {
      plot_object <- get_create_figure(dataset = list_arg[['dataset']], 
                                   model = model_old, 
                                   x = xsi, y = ysi,
                                   fullknots = knots_old,
                                   title = paste0("knots: ", knots_old_str),
                                   subtitle = setsubtitle,
                                   verbose = verbose)
    }
    
    if(what == 'plot2') {
      plot_object <- get_plot_model(dataset = list_arg[['dataset']], 
                                     model = model_old, 
                                     x = xsi, y = ysi,
                                     fullknots = knots_old,
                                     title = paste0("knots: ", knots_old_str),
                                    subtitle = setsubtitle,
                                     verbose = verbose)
    }
    
    if(what == 'plot3') {
      plot_object <- get_plot_rstandard(dataset = list_arg[['dataset']], 
                                        model = model_old, 
                                        x = xsi, y = ysi,
                                        fullknots = knots_old,
                                        title = 
                                          paste0("knots: ", knots_old_str),
                                        subtitle = setsubtitle,
                                        verbose = verbose)
    }
    return(plot_object)
  } # if(is.null(model_new)) {
  
  ###########################################################################
  ###########################################################################
  
  # if(!is.null(model_new)) ... continue below
  
  if(select == 'knots') {
    if(select == 'un') {
      if(length(knots_old) > length(knots_new)) {
        stop2c("The length of old knots is greater than the new knots. ",
               "The difference is ", 
               length(knots_old), ' - ', length(knots_new),
               " = ", length(knots_old) - length(knots_new),
               ". Either decrease the 'df', or increase the 'nsearch' ",
               "in 'knots_selection' argument by ", 
               length(knots_old) - length(knots_new))
      } # if(length(knots_old) > length(knots_new)) {
    } # if(select == 'un') {
  } # if(select == 'knots') {
  
  
  model_new <- get_model_by_knots_wrapper(list_arg, knots_new) 
  
  ######
  knots_new_str <-  paste(knots_new, collapse = " ")
  
  knots_selection_what_choices <- c('knots', 
                                    'plot1',
                                    'plot2', 
                                    'plot3', 
                                    'plot4',
                                    'plot5',
                                    'plot6',
                                    'plot7',
                                    'plot8',
                                    'plot9')
  
  if(!knots_selection_what %in% knots_selection_what_choices) {
    stop("knots_selection_what must be one of the following: ",
         "\n ",
         collapse_comma(knots_selection_what_choices))
  }
  
  
  if(knots_selection_what == 'plot1' | 
     knots_selection_what == 'plot2' | 
     knots_selection_what == 'plot3') {
    fig_old <- get_create_figure(dataset = 
                                   list_arg[['dataset']], 
                                 model = model_old, 
                                 x = xsi, y = ysi,
                                 fullknots = knots_old,
                                 title = 
                                   paste0("Old knots: ", knots_old_str),
                                 subtitle = setsubtitle,
                                 verbose = verbose)
    fig_new <- get_create_figure(dataset = 
                                   list_arg[['dataset']], 
                                 model = model_new, 
                                 x = xsi, y = ysi,
                                 fullknots = knots_new,
                                 title = 
                                   paste0("New knots: ", knots_new_str),
                                 subtitle = setsubtitle,
                                 verbose = verbose)
    fig_old_new <- 
      patchwork::wrap_plots(fig_old / fig_new)
  }
  
  if(knots_selection_what == 'plot4' | 
     knots_selection_what == 'plot5' | 
     knots_selection_what == 'plot6') {
    fig_old_plot <- get_plot_model(dataset = 
                                     list_arg[['dataset']], 
                                   model = model_old, 
                                   x = xsi, y = ysi,
                                   fullknots = knots_old,
                                   title = 
                                     paste0("Old knots: ", knots_old_str),
                                   subtitle = setsubtitle,
                                   verbose = verbose)
    fig_new_plot <- get_plot_model(dataset = 
                                     list_arg[['dataset']], 
                                   model = model_new, 
                                   x = xsi, y = ysi,
                                   fullknots = knots_new,
                                   title = 
                                     paste0("New knots: ", knots_new_str),
                                   subtitle = setsubtitle,
                                   verbose = verbose)
    fig_old_new_plot <- 
      patchwork::wrap_plots(fig_old_plot / fig_new_plot)
  }
  
  if(knots_selection_what == 'plot7' | 
     knots_selection_what == 'plot8' | 
     knots_selection_what == 'plot9') {
    fig_old_resid <- 
      get_plot_rstandard(dataset = 
                           list_arg[['dataset']], 
                         model = model_old, 
                         x = xsi, y = ysi,
                         fullknots = knots_old,
                         title = 
                           paste0("Old knots: ", knots_old_str),
                         subtitle = setsubtitle,
                         verbose = verbose)
    fig_new_resid <- 
      get_plot_rstandard(dataset = 
                           list_arg[['dataset']], 
                         model = model_new, 
                         x = xsi, y = ysi,
                         fullknots = knots_new,
                         title = 
                           paste0("New knots: ", knots_new_str),
                         subtitle = setsubtitle,
                         verbose = verbose)
    fig_old_new_resid <- 
      patchwork::wrap_plots(fig_old_resid / fig_new_resid)
  }
  
 
  if(knots_selection_what == 'knots') {
    knots_old_knots_new_msg <- paste0("Default knots ", 
                                      deparse(knots_old), " ",
                                      "replaced by new knots ", 
                                      deparse(knots_new))
    plot_object <- knots_old_knots_new_msg
  } else if(knots_selection_what == 'plot1') {
    plot_object <- fig_old
  } else if(knots_selection_what == 'plot2') {
    plot_object <- fig_new
  } else if(knots_selection_what == 'plot3') {
    plot_object <- fig_old_new
  } else if(knots_selection_what == 'plot4') {
    plot_object <- fig_old_plot
  } else if(knots_selection_what == 'plot5') {
    plot_object <- fig_new_plot
  } else if(knots_selection_what == 'plot6') {
    plot_object <- fig_old_new_plot
  } else if(knots_selection_what == 'plot7') {
    plot_object <- fig_old_resid
  } else if(knots_selection_what == 'plot8') {
    plot_object <- fig_new_resid
  } else if(knots_selection_what == 'plot9') {
    plot_object <- fig_old_new_resid
  } else if(knots_selection_what == 'plot10') {
    plot_object <- NULL
  } else if(knots_selection_what == 'plot11') {
    plot_object <- NULL
  } 

  return(plot_object)
}




##############################################################################
##############################################################################

# https://github.com/davidaarmstrong/damisc


##############################################################################
##############################################################################




##############################################################################
# get_NKnots 
##############################################################################

# https://github.com/davidaarmstrong/damisc/blob/master/R/DAMisc_functions.R

#' AIC and BIC selection of number of spline knots
#'
#' Calculates AIC and BIC for the selection of knots in a spline over values
#' (potentially including polynomials) up to a user-defined maximum.
#'
#'
#' @param form A formula detailing the model for which smoothing is to be
#'   evaluated.
#' @param var A character string identifying the variable for which smoothing is
#'   to be evaluated.
#' @param data Data frame providing values of all variables in \code{form}.
#' @param degree Degree of polynomial in B-spline basis functions.
#' @param min.knots Minimum number of internal B-spline knots to be evaluated.
#' @param max.knots Maximum number of internal B-spline knots to be evaluated.
#' @param includePoly Include linear and polynomial models up to, and including
#'   \code{degree}-th order polynomials.
#' @param plot Logical indicating whether a plot should be returned.
#' @param criterion Statistical criterion to minimize in order to find the best
#'   number of knots - AIC, BIC or Cross-validation. Ignored here
#' @param cvk Number of groups for cross-validation. cross-validation ignored
#'   here.
#' @param cviter Number of iterations of cross-validation to average over.
#'   10 is the default but in real-world applications, this should be
#'   somewhere around 200.
#' @returns A plot, if \code{plot=TRUE}, otherwise a data frame with the degrees
#'   of freedom and corresponding fit measure.
#' @author Dave Armstrong
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' data(Prestige, package="carData")
#' get_NKnots(prestige ~ education + type, var="income", 
#' data=na.omit(Prestige), plot=FALSE)
#' 
get_NKnots <- function(form, var, data, degree=3, min.knots=1,
                   max.knots=10, includePoly = FALSE, plot=FALSE, 
                   criterion=c("AIC", "BIC", "CV"),
                   cvk=10, cviter=10, smat = 'nsk'){
  crit <- match.arg(criterion)
  k <- seq(min.knots, max.knots, by=1)
  forms <- vector("list", ifelse(includePoly, length(k)+3, length(k)))
  df_poly <- NULL
  m <- 1
  if(includePoly){
    df_poly <- 1:3
    forms[[1]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], " + ", var, sep=""))
    forms[[2]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], 
                                  "+ poly(", var,  ", 2)", sep=""))
    forms[[3]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], 
                                  "+ poly(", var,  ", 3)", sep=""))
    m <- 4
  }
  
  df_spline <- NULL
  for(i in 1:length(k)){
    df_spline <- c(df_spline, degree+i)
    forms[[m]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], 
                                  "+ splines2::bsp(", var, ", df=", degree+k[i],
                                  ", Boundary.knots=c(", 
                                  min(data[[var]], na.rm=TRUE),", ", 
                                  max(data[[var]], 
                                      na.rm=TRUE), "))", 
                                  sep=""))
    m <- m+1
  }
  
  if(crit %in% c("AIC", "BIC")){
    mods <- lapply(forms, function(x)lm(x, data=data))
    getallscores <- sapply(mods, function(x)do.call(crit, list(object=x)))
  }
  if(crit == "CV"){
    stop("criterion 'CV' not supported")
  }
  
  if(plot){
    kgetallscores <- k+3
    if(includePoly){kgetallscores <- c(1:3, k)}
    plot(kgetallscores, getallscores, type="o", pch=16, 
         col="black", xlab="# Degrees of Freedom", ylab = crit)
    graphics::points(kgetallscores[which.min(getallscores)], 
                     min(getallscores), pch=16, col="red")
  }else{
    return(data.frame(df = c(df_poly, df_spline), stat=getallscores))
  }
}



##############################################################################
# get_NKnotsTest 
##############################################################################

#' Test of functional form assumption using B-splines
#'
#' Estimate hypothesis test of lower- and higher-order non-linear relationships
#' against an assumed target relationship.
#'
#'
#' @param form A formula detailing the model for which smoothing is to be
#'   evaluated.
#' @param var A character string identifying the variable for which smoothing is
#'   to be evaluated.
#' @param data Data frame providing values of all variables in \code{form}.
#' @param targetdf The assumed degrees of freedom against which the tests will
#'   be conducted.
#' @param degree Degree of polynomial in B-spline basis functions.
#' @param min.knots Minimum number of internal B-spline knots to be evaluated.
#' @param max.knots Maximum number of internal B-spline knots to be evaluated.
#' @param adjust Method by which p-values will be adjusted (see
#'   \code{\link{p.adjust}})
#' @returns A matrix with the following columns: \item{F}{F statistics of test
#'   of candidate models against target model} \item{DF1}{Numerator DF from
#'   F-test} \item{DF2}{Denominator DF from F-test} \item{p(F)}{p-value from the
#'   F-test} \item{Clarke}{Test statistic from the Clarke test}
#' \item{Pr(Better)}{The Clarke statistic divided by the number of
#' observations} \item{p(Clarke)}{p-value from the Clarke test.  (T) means that
#' the significant p-value is in favor of the Target model and (C) means the
#' significant p-value is in favor of the candidate (alternative) model.}
#'   \item{Delta_AIC}{AIC(candidate model) - AIC(target model)}
#'   \item{Delta_AICc}{AICc(candidate model) - AICc(target model)}
#'   \item{Delta_BIC}{BIC(candidate model) - BIC(target model)}
#' @author Dave Armstrong
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#'
#' data(Prestige, package="carData")
#' get_NKnotsTest(prestige ~ education + type, var="income", data=na.omit(Prestige), targetdf=3)
#' 
get_NKnotsTest <- function(form, var, data, targetdf = 1, degree=3, min.knots=1,
                       max.knots=10, adjust="none"){
  k <- seq(min.knots, max.knots, by=1)
  forms <- vector("list", length(k)+3)
  m <- 1
  forms[[1]]<- as.formula(paste(as.character(form)[2], "~",
                                as.character(form)[3], " + ", var, sep=""))
  forms[[2]]<- as.formula(paste(as.character(form)[2], "~",
                                as.character(form)[3], 
                                "+ poly(", var,  ", 2)", sep=""))
  forms[[3]]<- as.formula(paste(as.character(form)[2], "~",
                                as.character(form)[3], 
                                "+ poly(", var,  ", 3)", sep=""))
  m <- 4
  for(i in 1:length(k)){
    forms[[m]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], 
                                  "+ bs(", var, ", df=", degree+k[i], ")", 
                                  sep=""))
    m <- m+1
  }
  mods <- lapply(forms, function(x)lm(x, data=data))
  mods.df <- c(1:3, k+3)
  
  target.mod <- mods[[which(mods.df == targetdf)]]
  cand.mods <- mods
  cand.mods[[which(mods.df == targetdf)]] <- NULL
  tests <- lapply(cand.mods, function(x)as.matrix(stats::anova(target.mod, x)))
  # tests2 <- lapply(cand.mods, function(x)clarke_test(target.mod, x))
  num.df <- sapply(tests, function(x)abs(diff(x[,1])))
  denom.df <- sapply(tests, function(x)min(x[,1]))
  Fstats <- sapply(tests, function(x)x[2,5])
  pval <- stats::p.adjust(sapply(tests, function(x)x[2,6]), method=adjust)
# cstats <- sapply(tests2, function(x)x$stat)
# cprobs <- sapply(tests2, function(x)x$stat/x$nobs)
# cminstat <- sapply(tests2, function(x)min(x$stat, x$nobs - x$stat))
# cbetter <-cp<-2*stats::pbinom(cminstat, sapply(tests2,function(x)x$nobs), 0.5)
# pref<-sapply(tests2,function(x)ifelse(x$stat > x$nobs-x$stat, "(T)", "(C)"))
# pref <- ifelse(cp > .05, "", pref)
  cstats <- NULL
  cprobs <- NULL
  pref   <- NULL
  cp     <- NULL
  
  delta.aic <- sapply(cand.mods, stats::AIC) - stats::AIC(target.mod)
  delta.bic <- sapply(cand.mods, stats::BIC) - stats::BIC(target.mod)
# delta.aicc<-sapply(cand.mods, AICcmodavg::AICc) - AICcmodavg::AICc(target.mod)
  delta.aicc <- NULL
  res <- cbind(Fstats, num.df, denom.df, pval, cstats, cprobs, 
               cp, delta.aic, delta.aicc, delta.bic)
  sigchar <- ifelse(res[,4] < .05, "*", " ")
  sigchar2 <- ifelse(res[,7] < .05, "*", " ")
  strres <- NULL
  digs <- c(3,0,0,3, 0, 3, 3, 3, 3, 3)
  for(i in 1:10){
    tmp <- sprintf(paste("%.", digs[i], "f", sep=""), res[,i])
    if(i == 1){
      tmp <- paste(tmp, sigchar, sep="")
    }
    if(i == 5){
      tmp <- paste(tmp, sigchar2, sep="")
    }
    if(i == 7){
      tmp <- paste(tmp, pref,  sep=" ")
    }
    
    strres <- cbind(strres,tmp )
  }
  colnames(strres) <- c("F", "DF1", "DF2", "p(F)", 
                        "Clarke", "Pr(Better)", 
                        "p(Clarke)", "Delta_AIC", 
                        "Delta_AICc", "Delta_BIC")
  rownames(strres) <- paste("DF=", targetdf, " vs. DF=", 
                            mods.df[-targetdf], sep="")
  if(targetdf > 1){
    below <- strres[1:(targetdf-1), , drop=F]
    above <- strres[targetdf:nrow(strres),, drop=F]
    strres <- rbind(below, rep("", 10), above)
    rownames(strres)[targetdf] <- "   Target"
  }
  print(strres, quote=FALSE)
}


##############################################################################
# design_sigmoid_knots
##############################################################################


design_sigmoid_knots <- function(dataset,
                                 x, 
                                 y, 
                                 max_knots = 8, 
                                 min_knots = 3, 
                                 criterion = stats::AIC,
                                 smat = 'ns', 
                                 ns_model_str = NULL,
                                 model_formula = NULL,
                                 fix_bknots = TRUE,
                                 bknots = NULL) {
  
  if(is.null(ns_model_str) & is.null(model_formula)) {
    stop("specify one of the 'ns_model_str' or 'model_formula'")
  }
           
  x <- dataset[[x]]
  y <- dataset[[y]]
  
  stopifnot(length(x) == length(y))
  
  min_knots <- min_knots - 1
  
  
  if(smat == 'rcs') {
    min_knots_toadd <- min_knots + 3
  } else {
    min_knots_toadd <- min_knots + 1
  }
                                  
  knots <- 
    stats::quantile (x, 
                     probs = 
                       seq(0, 1, 
                           length.out = min_knots_toadd))[-c(1, 
                                                             min_knots_toadd)]  
 
  
  knots  <- knots %>% unname()
  
  best_score <- Inf
  best_model <- NULL
  best_knots <- knots
  
  
  
  if(is.null(bknots)) bknots <- range(x)
  
  bknots <- bknots %>% unname()
  
  old_bknots <- bknots
  
  if(smat == 'rcs') {
    fullknots <- c(bknots[1], knots, bknots[2])
    knots     <- checkgetiknotsbknots(fullknots, 'iknots')
    bknots     <- checkgetiknotsbknots(fullknots, 'bknots')
  } else {
    knots  <- knots
    bknots <- bknots
  }
  
  
  knots <- unique(knots)
  
  
 
  
  set_min_knots_k <- min_knots
  set_max_knots_k <- max_knots
  
  if(smat == 'rcs') {
    set_min_knots_k <- min_knots + 0
    set_max_knots_k <- max_knots - 0
  }
  
  for (k in set_min_knots_k:set_max_knots_k) {
    old_knots <- knots
    if(fix_bknots) old_bknots <- bknots
    repeat {
      # Fit current spline
      model_formula[[3]][['knots']]  <- knots
      model_formula[[3]][['bknots']] <- bknots
      lm_call <- bquote(stats::glm(formula = .(model_formula), data = dataset))
      fit    <- eval(lm_call)
      preds  <- predict(fit)
      resids <- abs(y - preds)
      
      # Candidate for improved knot placement: max residual
      candidate_x <- x[which.max(resids)]
      knots <- sort(unique(c(knots, candidate_x)))
      
      if(smat == 'rcs') {
        # if(min_knots < 2) break
        fullknots <- c(bknots[1], knots, bknots[2])
        knots     <- checkgetiknotsbknots(fullknots, 'iknots')
        bknots     <- checkgetiknotsbknots(fullknots, 'bknots')
      } else {
        knots  <- knots
        bknots <- bknots
      }
     
     
      if (length(knots) > k) break
      # Evaluate model with updated knots
      model_formula[[3]][['knots']]  <- knots 
      model_formula[[3]][['bknots']] <- bknots 
      lm_call <- bquote(stats::glm(formula = .(model_formula), data = dataset))
      fit <- eval(lm_call)
      # fit <- lm(y ~ ns(x, knots = knots, Boundary.knots = bknots))
      score <- criterion(fit)
      if (score < best_score) {
        best_score <- score
        best_model <- fit
        best_knots <- knots
      } else {
        break
      }
    }
    knots <- old_knots
    if(fix_bknots) bknots <- old_bknots
  }
  
  
  best_knots_x <- best_knots
  if(smat == 'rcs') {
    best_knots     <- checkgetiknotsbknots(best_knots_x, 'iknots')
    Boundary.knots <- checkgetiknotsbknots(best_knots_x, 'bknots')
  } else {
    best_knots     <- best_knots_x
    Boundary.knots <- bknots
  }
  
  if(fix_bknots) Boundary.knots <- old_bknots
  
  # returned is the best_knots
  list(model = best_model, knots = best_knots, Boundary.knots = Boundary.knots)
}



