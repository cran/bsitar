


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
        df <- 2 
        knots <- ""
        bknots <- ""
      } else {
        df <- 1
        knots <- ""
      }
    } 
  } 
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
        for(i in 1:nforms){ 
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
  }) 
  return(ret)
}





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
  }) 
  scores <- unlist(lapply(model_scores, "[[", "score"))
  index <- which.min(scores)
  min_score <- scores[[index]]
  return(list(model = model_scores[[index]][["model"]],
              score = min_score, index = index))
}





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




get_extract_knots <- function(ns_model) {
  knots <- attr(ns_model$model[[2]], "knots")
  knots <- unique(knots)
  bknots <- attr(ns_model$model[[2]], "Boundary.knots")
  fullknots <- c(bknots[1], knots, bknots[2])
  return(list(knots = knots, Boundary.knots = bknots, fullknots = fullknots))
}




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





get_print_knots <- function(knot_placements) {
  R.utils::printf("Inner knots count: %d\n", length(knot_placements$knots))
  knots_str <- paste0("[", paste0(knot_placements$knots, collapse = ", "), "]")
  boundary_str <-
    paste0("[", paste0(knot_placements$Boundary.knots, collapse = ", "), "]")
  R.utils::printf("Inner knots: %s\nBoundary knots: %s\n", knots_str,
                  boundary_str)
}






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
  return(list(
    df = min_ndf,
    nknots = min_ndf - 1, score = min_icr, 
    method = method, smat = smat,
    all_scores_df = all_scores_df,
    all_scores = list(scores = scores, n_knots = n_knots)
  ))
}





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

  if (initial_nknots == -1) {
    initial_nknots <-
      get_suggest_knotcount(dataset = dataset, 
                            dependent = dependent, 
                            independents = independents, 
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
  if(nys > 1) {
    setsubtitle <- NULL
  } else {
    setsubtitle <- NULL
  }
  if(knots_selection_what == "none") {
    return(knots_new)
  }
  if(is.null(model_new)) {
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
  } 
  
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
      } 
    } 
  } 
  
  model_new <- get_model_by_knots_wrapper(list_arg, knots_new) 
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
  num.df <- sapply(tests, function(x)abs(diff(x[,1])))
  denom.df <- sapply(tests, function(x)min(x[,1]))
  Fstats <- sapply(tests, function(x)x[2,5])
  pval <- stats::p.adjust(sapply(tests, function(x)x[2,6]), method=adjust)
  cstats <- NULL
  cprobs <- NULL
  pref   <- NULL
  cp     <- NULL
  delta.aic <- sapply(cand.mods, stats::AIC) - stats::AIC(target.mod)
  delta.bic <- sapply(cand.mods, stats::BIC) - stats::BIC(target.mod)
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
      model_formula[[3]][['knots']]  <- knots
      model_formula[[3]][['bknots']] <- bknots
      lm_call <- bquote(stats::glm(formula = .(model_formula), data = dataset))
      fit    <- eval(lm_call)
      preds  <- predict(fit)
      resids <- abs(y - preds)
      candidate_x <- x[which.max(resids)]
      knots <- sort(unique(c(knots, candidate_x)))
      if(smat == 'rcs') {
        fullknots <- c(bknots[1], knots, bknots[2])
        knots     <- checkgetiknotsbknots(fullknots, 'iknots')
        bknots     <- checkgetiknotsbknots(fullknots, 'bknots')
      } else {
        knots  <- knots
        bknots <- bknots
      }
      if (length(knots) > k) break
      model_formula[[3]][['knots']]  <- knots 
      model_formula[[3]][['bknots']] <- bknots 
      lm_call <- bquote(stats::glm(formula = .(model_formula), data = dataset))
      fit <- eval(lm_call)
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
  list(model = best_model, knots = best_knots, Boundary.knots = Boundary.knots)
}





