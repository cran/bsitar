

#' An internal function to get adjusted curves (for random effects)
#' 
#' @exportS3Method xyadj_curves bgmfit
#' @noRd
#'
xyadj_curves.bgmfit <-
  function (model,
            x = NULL,
            y = NULL,
            id = NULL,
            v = NULL,
            resp = NULL,
            newdata = NULL,
            ndraws = NULL,
            draw_ids = NULL,
            levels_id = NULL,
            abc = NULL,
            summary = FALSE,
            conf = 0.95,
            robust = FALSE,
            tomean = TRUE,
            ipts = NULL,
            xrange = NULL,
            aux_variables = NULL,
            numeric_cov_at = NULL,
            idata_method = NULL,
            verbose = FALSE,
            model_deriv = NULL,
            deriv = NULL, 
            envir = NULL,
            ...) {
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- envir
    }
    if(!is.na(model$model_info$univariate_by$by)) {
      stop("option = 'a' is not yet available for 'univariate_by' model")
    } else if(model$model_info$multivariate$mvar) {
      stop("option = 'a' is not yet available for 'multivariate' model")
    } 
    if(!is.null(ipts)) 
      stop("It does not a make sense to interploate data when estimating",
           "\n ",
           " adjusted curves. Please set ipts = NULL")
    xvar <- NULL;
    yvar <- NULL;
    idvar <- NULL;
    cov_vars <- NULL;
    cov_factor_vars <- NULL;
    cov_numeric_vars <- NULL;
    groupby_fstr <- NULL;
    groupby_fistr <- NULL;
    uvarby <- NULL;
    subindicatorsi <- NULL;
    . <- NULL;
    setxcall_ <- match.call()
    post_processing_checks_args <- list()
    post_processing_checks_args[['model']]    <- model
    post_processing_checks_args[['xcall']]    <- setxcall_
    post_processing_checks_args[['resp']]     <- resp
    post_processing_checks_args[['envir']]    <- envir
    post_processing_checks_args[['deriv']]    <- ''
    post_processing_checks_args[['all']]      <- FALSE
    post_processing_checks_args[['verbose']]  <- verbose
    post_processing_checks_args[['check_d0']] <- FALSE
    post_processing_checks_args[['check_d1']] <- TRUE
    post_processing_checks_args[['check_d2']] <- FALSE
     o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
    if (is.null(ndraws))
      ndraws  <- brms::ndraws(model)
    else
      ndraws <- ndraws
    if (is.null(idata_method)) {
      idata_method <- 'm2'
    }
    if (is.null(resp)) {
      resp_rev_ <- resp
    } else if (!is.null(resp)) {
      resp_rev_ <- paste0("_", resp)
    }
    if (is.null(newdata)) {
      stop("Please specify newdata")
    }
    newdata_dummy <- get.newdata(model,
                           newdata = newdata,
                           resp = resp,
                           numeric_cov_at = numeric_cov_at,
                           aux_variables = aux_variables,
                           levels_id = levels_id,
                           ipts = ipts,
                           xrange = xrange,
                           idata_method = idata_method,
                           verbose = verbose)
    list_c <- attr(newdata_dummy, 'list_c')
    for (list_ci in names(list_c)) {
      assign(list_ci, list_c[[list_ci]])
    }
    rm('newdata_dummy')
    check__ <- c('xvar', 'yvar', 'idvar', 'cov_vars', 'cov_factor_vars', 
                 'cov_numeric_vars', 'groupby_fstr', 'groupby_fistr', 
                 'uvarby', 'subindicatorsi')
    for (check___ in check__) {
      if(!exists(check___)) assign(check___, NULL)
    }
    if(is.null(uvarby)) uvarby <- NA
    Xx <- xvar
    Yy <- yvar

    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', 'Est.Error', probtitles)
    
    xoffset <- 0
    d_adjustedXnames <- 'd_adjusted'
    if (!is.null(resp)) d_adjustedXnames <- paste0(d_adjustedXnames, 
                                                   resp_rev_)
    d_adjustedXnames <- model$model_info[[d_adjustedXnames]]
    d_adjusted <- d_adjustedXnames
    if(is.null(x)) {
      x <- newdata[[Xx]]
    } else {
      if(!is.numeric(x)) {
        stop("'x' must be numeric")
      }
      x <- x
    }
    if(is.null(y)) {
      y <- newdata[[Yy]]
    } else {
      if(!is.numeric(y)) {
        stop("'y' must be numeric")
      }
      y <- y
    }
    if(is.null(x)) {
      x <- newdata[[Xx]]
    } else {
      if(!is.numeric(x)) {
        stop("'x' must be numeric")
      }
      x <- x
    }
    if(is.null(v)) {
      v <- 0
    } else {
      if(!is.numeric(v)) {
        stop("'v' must be numeric")
      }
      v <- v
    }
    if(is.null(id)) {
      idvar <- model$model_info$idvars
      idvar <- idvar[1]
      id <- newdata[[idvar]][1]
    }
    if(!is.null(ipts)) {
      add_outcome <- model$data %>%
        dplyr::select(dplyr::all_of(c(Yy, idvar)))
      newdata <- newdata %>% 
        dplyr::left_join(., add_outcome, by = c(idvar))
      x <- newdata[[Xx]]
      y <- newdata[[Yy]]
      id <- newdata[[idvar]][1]
    }
    nrowdatadims <- nrow(newdata)
    predprep <- brms::prepare_predictions(model, resp = resp, 
                                          newdata = newdata)
    rparnames <- names(predprep$nlpars)
    respstr <- "" 
    septsr  <- ""
    if(any(grepl(paste0(respstr, septsr, "a"), rparnames)) |
       any(grepl(paste0(respstr, "a", septsr), rparnames))) {
      a_r <- TRUE
    } else {
      a_r <- FALSE
    }
    if(any(grepl(paste0(respstr, septsr, "b"), rparnames)) |
       any(grepl(paste0(respstr, "b", septsr), rparnames))) {
      b_r <- TRUE
    } else {
      b_r <- FALSE
    }
    if(any(grepl(paste0(respstr, septsr, "c"), rparnames)) |
       any(grepl(paste0(respstr, "c", septsr), rparnames))) {
      c_r <- TRUE
    } else {
      c_r <- FALSE
    }
    if(any(grepl(paste0(respstr, septsr, "d"), rparnames)) |
       any(grepl(paste0(respstr, "d", septsr), rparnames))) {
      d_r <- TRUE
    } else {
      d_r <- FALSE
    }
    
    if(!is.null(draw_ids)) {
      setloopdins <- length(eval(draw_ids))
    } else if(!is.null(ndraws)) {
      setloopdins <- ndraws
    } 
    
    if(a_r) {
      null_a <- fitted(model, resp = resp, newdata = newdata, 
                       nlpar="a", ndraws = ndraws,  draw_ids = draw_ids,
                       re_formula = NULL, summary = summary,
                       fullframe = NULL, itransform = "")
      naaa_a <- fitted(model, resp = resp, newdata = newdata, 
                       nlpar="a", ndraws = ndraws,  draw_ids = draw_ids,
                       re_formula = NA, summary = summary,
                       fullframe = NULL, itransform = "")
    } else {
      if(summary) {
        null_a <- matrix(0, nrowdatadims, 1)
        naaa_a <- matrix(0, nrowdatadims, 1)
      } else if(!summary) {
        null_a <- matrix(0, setloopdins, nrowdatadims)
        naaa_a <- matrix(0, setloopdins, nrowdatadims)
      }
    }
    if(b_r) {
      null_b <- fitted(model, resp = resp, newdata = newdata, 
                       nlpar="b", ndraws = ndraws,  draw_ids = draw_ids,
                       re_formula = NULL, summary = summary,
                       fullframe = NULL, itransform = "")
      naaa_b <- fitted(model, resp = resp, newdata = newdata, 
                       nlpar="b", ndraws = ndraws,  draw_ids = draw_ids,
                       re_formula = NA, summary = summary,
                       fullframe = NULL, itransform = "")
    } else {
      if(summary) {
        null_b <- matrix(0, nrowdatadims, 1)
        naaa_b <- matrix(0, nrowdatadims, 1)
      } else if(!summary) {
        null_b <- matrix(0, setloopdins, nrowdatadims)
        naaa_b <- matrix(0, setloopdins, nrowdatadims)
      }
    }
    if(c_r) {
      null_c <- fitted(model, resp = resp, newdata = newdata, 
                       nlpar="c", ndraws = ndraws,  draw_ids = draw_ids,
                       re_formula = NULL, summary = summary,
                       fullframe = NULL, itransform = "")
      naaa_c <- fitted(model, resp = resp, newdata = newdata, 
                       nlpar="c", ndraws = ndraws,  draw_ids = draw_ids,
                       re_formula = NA, summary = summary,
                       fullframe = NULL, itransform = "")
    } else {
      if(summary) {
        null_c <- matrix(0, nrowdatadims, 1)
        naaa_c <- matrix(0, nrowdatadims, 1)
      } else if(!summary) {
        null_c <- matrix(0, setloopdins, nrowdatadims)
        naaa_c <- matrix(0, setloopdins, nrowdatadims)
      }
    }
    if(d_r) {
      null_d <- fitted(model, resp = resp, newdata = newdata, 
                       nlpar="d", ndraws = ndraws,  draw_ids = draw_ids,
                       re_formula = NULL, summary = summary,
                       fullframe = NULL, itransform = "")
      naaa_d <- fitted(model, resp = resp, newdata = newdata, 
                       nlpar="d", ndraws = ndraws,  draw_ids = draw_ids,
                       re_formula = NA, summary = summary,
                       fullframe = NULL, itransform = "")
    } else {
      if(summary) {
        null_d <- matrix(0, nrowdatadims, 1)
        naaa_d <- matrix(0, nrowdatadims, 1)
      } else if(!summary) {
        null_d <- matrix(0, setloopdins, nrowdatadims)
        naaa_d <- matrix(0, setloopdins, nrowdatadims)
      }
    }
    
    dots <- list(...)
    set_get_dv <- FALSE
    if(!is.null(dots$get_dv)) {
      if(dots$get_dv) {
        if(verbose) message("executing 'get_dv'!")
        set_get_dv <- TRUE
      }
    }
    
    if(!summary) {
      xadj_tmt <- yadj_tmt <- vadj_tmt <- list()
      xadj_tmf <- yadj_tmf <- vadj_tmf <- list()
      for (i in 1:setloopdins) {
        r_a <- null_a[ i, ]
        r_b <- null_b[ i, ]
        r_c <- null_c[ i, ]
        r_d <- null_d[ i, ]
        na_a <- naaa_a[ i, ]
        na_b <- naaa_b[ i, ]
        na_c <- naaa_c[ i, ]
        na_d <- naaa_d[ i, ]
        rz_a <- r_a - na_a
        rz_b <- r_b - na_b
        rz_c <- r_c - na_c
        rz_d <- r_d - na_d
        r_data_ <- cbind(rz_a, rz_b, rz_c, rz_d) %>% data.frame()
        colnames(r_data_) <- letters[1:4]
        r_data_ <- r_data_ %>% 
          dplyr::mutate(x = x) %>% 
          dplyr::mutate(d.adjusted = d_adjusted %||% FALSE)
        adj_tmt <- r_data_ %>%
          dplyr::mutate(x.adj = (x - .data$b) * exp(.data$c) + xoffset,
                        y.adj = y - .data$a - 
                          .data$d * dplyr::if_else(.data$d.adjusted,
                                                   .data$x.adj - xoffset,
                                                   x),
                        v.adj = dplyr::if_else(.data$d.adjusted,
                                               v / exp(.data$c) - .data$d,
                                               (v - .data$d) / exp(.data$c)))
        adj_tmf <- r_data_ %>%
          dplyr::mutate(x.adj = x / exp(.data$c) + .data$b + xoffset,
                        y.adj = y + .data$a + 
                          .data$d * dplyr::if_else(.data$d.adjusted,
                                                   .data$x.adj - xoffset,
                                                   x),
                        v.adj = dplyr::if_else(.data$d.adjusted,
                                               (v + .data$d) * exp(.data$c),
                                               v * exp(.data$c) + .data$d))
        adj_tmt <- adj_tmt %>% data.frame()
        adj_tmf <- adj_tmf %>% data.frame()
        xadj_tmt[[i]] <- adj_tmt %>% dplyr::select(dplyr::all_of('x.adj')) %>% 
          unlist() %>% as.numeric()
        yadj_tmt[[i]] <- adj_tmt %>% dplyr::select(dplyr::all_of('y.adj')) %>% 
          unlist() %>% as.numeric()
        vadj_tmt[[i]] <- adj_tmt %>% dplyr::select(dplyr::all_of('v.adj')) %>% 
          unlist() %>% as.numeric()
        xadj_tmf[[i]] <- adj_tmf %>% dplyr::select(dplyr::all_of('x.adj')) %>% 
          unlist() %>% as.numeric()
        yadj_tmf[[i]] <- adj_tmf %>% dplyr::select(dplyr::all_of('y.adj')) %>% 
          unlist() %>% as.numeric()
        vadj_tmf[[i]] <- adj_tmf %>% dplyr::select(dplyr::all_of('v.adj')) %>% 
          unlist() %>% as.numeric()
      } 
      xadj_tmt <- array(unlist(xadj_tmt), 
                        dim=c(length(xadj_tmt[[1]]), length(xadj_tmt)  ))
      xadj_tmt <- t(xadj_tmt)
      yadj_tmt <- array(unlist(yadj_tmt), 
                        dim=c(length(yadj_tmt[[1]]), length(yadj_tmt)  ))
      yadj_tmt <- t(yadj_tmt)
      vadj_tmt <- array(unlist(vadj_tmt), 
                        dim=c(length(vadj_tmt[[1]]), length(vadj_tmt)  ))
      vadj_tmt <- t(vadj_tmt)
      xadj_tmf <- array(unlist(xadj_tmf), 
                        dim=c(length(xadj_tmf[[1]]), length(xadj_tmf)  ))
      xadj_tmf <- t(xadj_tmf)
      yadj_tmf <- array(unlist(yadj_tmf), 
                        dim=c(length(yadj_tmf[[1]]), length(yadj_tmf)  ))
      yadj_tmf <- t(yadj_tmf)
      vadj_tmf <- array(unlist(vadj_tmf), 
                        dim=c(length(vadj_tmf[[1]]), length(vadj_tmf)  ))
      vadj_tmf <- t(vadj_tmf)
      if(set_get_dv) {
        if(tomean) stop("'tomean' must be FALSE when 'get_dv = TRUE'")
        if(!tomean) return(xadj_tmf)
        if(tomean)  return(xadj_tmt)
      }
      if(!is.null(dots$xadj_tmt)) {
        if(dots$xadj_tmt) {
          if(verbose) message("returning 'xadj' tomean = TRUE")
          return(xadj_tmt)
        }
      }
      if(!is.null(dots$xadj_tmf)) {
        if(dots$xadj_tmf) {
          if(verbose) message("returning xadj tomean = FALSE")
          return(xadj_tmf)
        }
      }
      xadj_tmt <- brms::posterior_summary(xadj_tmt, probs = probs, 
                                          robust = robust) 
      yadj_tmt <- brms::posterior_summary(yadj_tmt, probs = probs, 
                                          robust = robust)
      vadj_tmt <- brms::posterior_summary(vadj_tmt, probs = probs, 
                                          robust = robust)
      xadj_tmf <- brms::posterior_summary(xadj_tmf, probs = probs, 
                                          robust = robust)
      yadj_tmf <- brms::posterior_summary(yadj_tmf, probs = probs, 
                                          robust = robust)
      vadj_tmf <- brms::posterior_summary(vadj_tmf, probs = probs, 
                                          robust = robust)
      if (tomean) {
        x.adj <- xadj_tmt
        y.adj <- yadj_tmt
        v.adj <- vadj_tmt
      }
      else {
        x.adj <- xadj_tmf
        y.adj <- yadj_tmf
        v.adj <- vadj_tmf
      }
      out <- newdata
      out[[Xx]] <- x.adj[, 1]
      out[[Yy]] <- y.adj[, 1]
      out <- out %>% dplyr::relocate(dplyr::all_of(c(Xx, Yy, idvar))) 
      out <- cbind(out, y.adj)
    } 
    
    if(summary) {
      r_a <- null_a[ , 1]
      r_b <- null_b[ , 1]
      r_c <- null_c[ , 1]
      r_d <- null_d[ , 1]
      na_a <- naaa_a[ , 1]
      na_b <- naaa_b[ , 1]
      na_c <- naaa_c[ , 1]
      na_d <- naaa_d[ , 1]
      rz_a <- r_a - na_a
      rz_b <- r_b - na_b
      rz_c <- r_c - na_c
      rz_d <- r_d - na_d
      r_data_ <- cbind(rz_a, rz_b, rz_c, rz_d) %>% data.frame()
      colnames(r_data_) <- letters[1:4]
      r_data_ <- r_data_ %>% 
        dplyr::mutate(x = x) %>% 
        dplyr::mutate(d.adjusted = d_adjusted %||% FALSE)
      adj_tmt <- r_data_ %>%
        dplyr::mutate(x.adj = (x - .data$b) * exp(.data$c) + xoffset,
                      y.adj = y - .data$a - 
                        .data$d * dplyr::if_else(.data$d.adjusted,
                                                 .data$x.adj - xoffset,
                                                 x),
                      v.adj = dplyr::if_else(.data$d.adjusted,
                                             v / exp(.data$c) - .data$d,
                                             (v - .data$d) / exp(.data$c)))
      adj_tmf <- r_data_ %>%
        dplyr::mutate(x.adj = x / exp(.data$c) + .data$b + xoffset,
                      y.adj = y + .data$a + 
                        .data$d * dplyr::if_else(.data$d.adjusted,
                                                 .data$x.adj - xoffset,
                                                 x),
                      v.adj = dplyr::if_else(.data$d.adjusted,
                                             (v + .data$d) * exp(.data$c),
                                             v * exp(.data$c) + .data$d))
      adj_tmt <- adj_tmt %>% data.frame()
      adj_tmf <- adj_tmf %>% data.frame()
      xadj_tmt <- adj_tmt %>% dplyr::select(dplyr::all_of('x.adj')) %>% 
        unlist() %>% as.numeric()
      yadj_tmt <- adj_tmt %>% dplyr::select(dplyr::all_of('y.adj')) %>% 
        unlist() %>% as.numeric()
      vadj_tmt <- adj_tmt %>% dplyr::select(dplyr::all_of('v.adj')) %>% 
        unlist() %>% as.numeric()
      xadj_tmf <- adj_tmf %>% dplyr::select(dplyr::all_of('x.adj')) %>% 
        unlist() %>% as.numeric()
      yadj_tmf <- adj_tmf %>% dplyr::select(dplyr::all_of('y.adj')) %>% 
        unlist() %>% as.numeric()
      vadj_tmf <- adj_tmf %>% dplyr::select(dplyr::all_of('v.adj')) %>% 
        unlist() %>% as.numeric()
      if (tomean) {
        x.adj <- xadj_tmt
        y.adj <- yadj_tmt
        v.adj <- vadj_tmt
      }
      else {
        x.adj <- xadj_tmf
        y.adj <- yadj_tmf
        v.adj <- vadj_tmf
      }
      out <- cbind(x.adj, y.adj)
      setadnamex <- paste0("adj", "_", Xx)
      setadnamey <- 'Estimate'
      colnames(out) <- c(setadnamex, setadnamey)
      out <- cbind(newdata, out)
      out <- newdata
      out[[Xx]] <- x.adj
      out[[Yy]] <- y.adj
      out <- out %>% dplyr::relocate(dplyr::all_of(c(Xx, Yy, idvar)))
    } 
    out
  } 


#' @noRd
#' @exportS3Method xyadj_curves bgmfit
xyadj_curves <- function(model, ...) {
  UseMethod("xyadj_curves")
}





#' Title An internal function to get unadjusted curves
#'
#' @exportS3Method xyunadj_curves bgmfit
#' @noRd
#'
xyunadj_curves.bgmfit <- function (model,
                                   x = NULL,
                                   y = NULL,
                                   id = NULL,
                                   newdata = NULL,
                                   ndraws = NULL,
                                   draw_ids = NULL,
                                   resp = NULL,
                                   verbose = FALSE,
                                   model_deriv = NULL,
                                   deriv = NULL, 
                                   envir = NULL,
                                   ...) {
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- envir # parent.frame()
  }
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  xvar <- NULL;
  yvar <- NULL;
  idvar <- NULL;
  cov_vars <- NULL;
  cov_factor_vars <- NULL;
  cov_numeric_vars <- NULL;
  groupby_fstr <- NULL;
  groupby_fistr <- NULL;
  uvarby <- NULL;
  subindicatorsi <- NULL;
  setxcall_ <- match.call()
  post_processing_checks_args <- list()
  post_processing_checks_args[['model']]    <- model
  post_processing_checks_args[['xcall']]    <- setxcall_
  post_processing_checks_args[['resp']]     <- resp
  post_processing_checks_args[['envir']]    <- envir
  post_processing_checks_args[['deriv']]    <- ''
  post_processing_checks_args[['all']]      <- FALSE
  post_processing_checks_args[['verbose']]  <- verbose
  post_processing_checks_args[['check_d0']] <- FALSE
  post_processing_checks_args[['check_d1']] <- TRUE
  post_processing_checks_args[['check_d2']] <- FALSE
  o <- CustomDoCall(post_processing_checks, post_processing_checks_args)
  newdata <- get.newdata(model, 
                         newdata = newdata, 
                         resp = resp, 
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
  if(is.null(uvarby)) uvarby <- NA
  Xx <- xvar
  Yy <- yvar
  if(!is.na(uvarby)) {
    newdata <- newdata %>%
      dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
      droplevels()
  }
  if(is.null(x)) {
    x <- newdata[[Xx]]
  } else {
    if(!is.numeric(x)) {
      stop("'x' must be numeric")
    }
    x <- x
  }
  if(is.null(y)) {
    y <- newdata[[Yy]]
  } else {
    if(!is.numeric(y)) {
      stop("'y' must be numeric")
    }
    y <- y
  }
  if(is.null(id)) {
    idvar <- model$model_info$idvars
    idvar <- idvar[1]
    id <- newdata[[idvar]][1]
  }
  out <- as.data.frame(as.factor(newdata[[idvar]]))
  out <- cbind(x, y, out)
  colnames(out) <- c(Xx, Yy, idvar)
  if(!is.na(uvarby)) {
    out[[uvarby]] <- resp
  }
  out
} 


#' @noRd
#' @exportS3Method xyunadj_curves bgmfit
xyunadj_curves <- function(model, ...) {
  UseMethod("xyunadj_curves")
}



#' Title An internal function to trim growth curves
#' 
#' @exportS3Method trimlines_curves bgmfit
#' @noRd
#' 
trimlines_curves.bgmfit <-
  function(model,
           x = NULL,
           y = NULL,
           id = NULL,
           newdata = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           resp = NULL,
           level = 0,
           trim = 0,
           estimation_method = 'fitted',
           verbose = FALSE,
           model_deriv = NULL,
           deriv = NULL, 
           envir = NULL,
           ...) {
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- envir # parent.frame()
    }
    if (is.null(ndraws))
      ndraws  <- brms::ndraws(model)
    else
      ndraws <- ndraws
    if (is.null(resp)) {
      resp_rev_ <- resp
    } else if (!is.null(resp)) {
      resp_rev_ <- paste0("_", resp)
    }
    xvar <- NULL;
    yvar <- NULL;
    idvar <- NULL;
    cov_vars <- NULL;
    cov_factor_vars <- NULL;
    cov_numeric_vars <- NULL;
    groupby_fstr <- NULL;
    groupby_fistr <- NULL;
    uvarby <- NULL;
    subindicatorsi <- NULL;
    Xx <- NULL;
    Yy <- NULL;
    dy <- NULL;
    setxcall_ <- match.call()
    post_processing_checks_args <- list()
    post_processing_checks_args[['model']]    <- model
    post_processing_checks_args[['xcall']]    <- setxcall_
    post_processing_checks_args[['resp']]     <- resp
    post_processing_checks_args[['envir']]    <- envir
    post_processing_checks_args[['deriv']]    <- ''
    post_processing_checks_args[['check_d0']] <- FALSE
    post_processing_checks_args[['check_d1']] <- TRUE
    post_processing_checks_args[['check_d2']] <- FALSE
    o <- CustomDoCall(post_processing_checks, post_processing_checks_args)
    uvarby <- model$model_info$univariate_by$by  
    newdata.o <- newdata
    if (trim == 0) {
      return(newdata)
    }
    if(is.null(x)) {
      .x <- Xx
    } else {
      if(is.symbol(x)) {
        .x <- deparse(x)
      } else if(is.character(x)) {
        .x <- x
      } else {
        stop("'x' must be NULL, a symbol or a character string")
      }
    }
    if(is.null(y)) {
      .y <- Yy
    } else {
      if(is.symbol(y)) {
        .y <- deparse(y)
      } else if(is.character(y)) {
        .y <- y
      } else {
        stop("'y' must be NULL, a symbol or a character string")
      }
    }
    if(is.null(id)) {
      idvar <- model$model_info$idvars
      idvar <- idvar[1]
      .id <- idvar
    } else {
      if(is.symbol(id)) {
        .id <- deparse(id)
      } else if(is.character(id)) {
        .id <- id
      } else {
        stop("'id' must be NULL, a symbol or a character string")
      }
    }
    newdata <- with(newdata, newdata[order(newdata[[.id]], newdata[[.x]]), ])
    extra <- dplyr::as_tibble(diff(as.matrix(newdata[, 1:2])))
    extra[[.id]] <- newdata[[.id]][-1]
    did <- diff(as.integer(newdata[[.id]]))
    extra$dx <- extra[[.x]]
    extra[, 1:2] <- newdata[-1, 1:2] - extra[, 1:2] / 2
    extra <- extra[!did, ]
    if(!is.na(uvarby)) {
      extra[[subindicatorsi]] <- 1
    }
    if (level == 0) {
      re_formula <- NA
    } else if (level == 1) {
      re_formula <- NULL
    }
    estimation_method_args <- list()
    estimation_method_args[['model']]      <- model
    estimation_method_args[['resp']]       <- resp
    estimation_method_args[['newdata']]    <- extra
    estimation_method_args[['ndraws']]     <- ndraws
    estimation_method_args[['re_formula']] <- re_formula
    estimation_method_args[['summary']]    <- TRUE
    estimation_method_args[['fullframe']]  <- NULL
    estimation_method_args[['itransform']] <- ""
    estimation_method_args[['envir']]      <- envir
    if (estimation_method == 'fitted') {
       extra$ey <- CustomDoCall(fitted_draws, estimation_method_args)
    } else if (estimation_method == 'predict') {
       extra$ey <- CustomDoCall(predict_draws, estimation_method_args)
    }
    extra$ey <- extra$ey[, 1]
    extra <- extra %>%
      dplyr::mutate(dy = abs(extra[[.y]] - extra$ey),
                    xy = extra$dx / mad(extra$dx) + dy / mad(dy))
    outliers <- order(extra$xy, decreasing = TRUE)[1:trim]
    extra <- extra[outliers, 1:3]
    extra[[.y]] <- NA
    if(!is.na(uvarby)) {
      newdata_tt <- newdata
      common_colsnms <- intersect(colnames(newdata) , colnames(extra))
      newdata <-newdata %>% dplyr::select(dplyr::all_of(common_colsnms))
    }
    newdata <- newdata %>% dplyr::select(dplyr::all_of(colnames(extra)))
    newdata <- rbind(newdata, extra)
    newdata <- with(newdata, newdata[order(newdata[[.id]], newdata[[.x]]), ])
    if(!is.na(uvarby)) {
      tempotnames <- c(idvar, Xx, Yy)
      tempot <- newdata_tt %>%  dplyr::select(-dplyr::all_of(tempotnames))
      newdata <- cbind(newdata[-1, ], tempot) %>% data.frame()
    }
    return(newdata)
  }



#' @noRd
#' @exportS3Method trimlines_curves bgmfit
trimlines_curves <- function(model, ...) {
  UseMethod("trimlines_curves")
}





set_lines_colors <- function(plot, ngroups, 
                             linetype.groupby,
                             color.groupby) {
  nrepvals <- ngroups
  if(is.null(linetype.groupby)) {
    linetype.groupby <- deparse(linetype.groupby)
  } else  if(is.na(linetype.groupby)) {
    linetype.groupby <- deparse(linetype.groupby)
  } else {
    linetype.groupby <- linetype.groupby
  }
  if(is.null(color.groupby)) {
    color.groupby <- deparse(color.groupby)
  } else  if(is.na(color.groupby)) {
    color.groupby <- deparse(color.groupby)
  } else {
    color.groupby <- color.groupby
  }
  ggplotColors <- function(g){
    g <- g - 1
    d <- 360/g
    h <- cumsum(c(15, rep(d,g - 1)))
    O <- grDevices::hcl(h = h, c = 100, l = 65)
    O <- c('black', O)
    O
  }
  ggplotlines <- function(g){
    lineTypes1 <- c("solid", "22", "42", "44", "13", "1343", "73", "2262")
    lineTypes2 <- apply(expand.grid(1:3, 1:3, 1:3, 1:3), 1, 
                        paste0, collapse="")
    lineTypes3 <- apply(expand.grid(1:2, 1:2, 1:2, 1:2), 1, 
                        paste0, collapse="")
    lineTypes <- c(lineTypes1, lineTypes2, lineTypes3)
    lineTypes[1:g]
  }
  default.set.line.groupby <- 'solid'
  default.set.color.groupby <- 'black'
  line.guide <- "none"
  color.guide <- "none"
  if(linetype.groupby == 'NA' & color.groupby == 'NA') {
    if(nrepvals == 1) {
      set.line.groupby <- default.set.line.groupby
      set.color.groupby <- default.set.color.groupby
    }
    if(nrepvals > 1) {
      set.line.groupby <- rep(default.set.line.groupby, nrepvals)
      set.color.groupby <- rep(default.set.color.groupby, nrepvals)
      line.guide <- "none"
      color.guide <- "legend"
    }
  } 
  if(linetype.groupby == 'NA' & color.groupby != 'NA') {
    set.line.groupby <- rep(default.set.line.groupby, nrepvals)
    if(nrepvals == 1) {
      if(color.groupby == 'NULL') {
        set.color.groupby <- default.set.color.groupby
      } else if(color.groupby != 'NULL') {
        set.color.groupby <- color.groupby[1]  
      }
    }
    if(nrepvals > 1) {
      set.line.groupby <- rep(default.set.line.groupby, nrepvals)
      if(color.groupby == 'NULL') {
        set.color.groupby <- ggplotColors(nrepvals)
      }
      if(color.groupby != 'NULL') {
        if(length(color.groupby) == nrepvals) {
          set.color.groupby <- color.groupby
        } else if(length(color.groupby) != nrepvals) {
          set.color.groupby <- rep(color.groupby, nrepvals)
        }
      }
      line.guide <- "none"
      color.guide <- "legend"
    }
  } 
  if(linetype.groupby != 'NA' & color.groupby == 'NA') {
    set.color.groupby <- rep(default.set.color.groupby, nrepvals)
    if(nrepvals == 1) {
      if(linetype.groupby == 'NULL') {
        set.line.groupby <- default.set.line.groupby
      } else if(linetype.groupby != 'NULL') {
        set.line.groupby <- linetype.groupby[1]  
      }
    }
    if(nrepvals > 1) {
      if(linetype.groupby == 'NULL') {
        set.line.groupby <- ggplotlines(nrepvals)
        if(length(set.line.groupby) < nrepvals) {
          set.line.groupby <- rep(set.line.groupby, nrepvals)
        }
      }
      if(linetype.groupby != 'NULL') {
        if(length(linetype.groupby) == nrepvals) {
          set.line.groupby <- linetype.groupby
        } else if(length(color.groupby) != nrepvals) {
          set.line.groupby <- rep(linetype.groupby, nrepvals)
        }
      }
      line.guide <- "none" 
      color.guide <- "legend"  
    }
  } 
  if(linetype.groupby != 'NA' & color.groupby != 'NA') {
    if(nrepvals == 1) {
      if(color.groupby == 'NULL') {
        set.color.groupby <- 'black'
      } else if(color.groupby != 'NULL') {
        set.color.groupby <- color.groupby[1]   
      }
      if(linetype.groupby == 'NULL') {
        set.line.groupby <- 'solid'
      } else if(linetype.groupby != 'NULL') {
        set.line.groupby <- linetype.groupby[1]   
      }
    }
    if(nrepvals > 1) {
      if(color.groupby == 'NULL') {
        set.color.groupby <- ggplotColors(nrepvals)
        if(length(set.color.groupby) < nrepvals) {
          set.color.groupby <- rep(set.color.groupby, nrepvals)
        }
      }
      if(linetype.groupby == 'NULL') {
        set.line.groupby <- ggplotlines(nrepvals)
        if(length(set.line.groupby) < nrepvals) {
          set.line.groupby <- rep(set.line.groupby, nrepvals)
        }
      }
      if(color.groupby != 'NULL') {
        if(length(color.groupby) == nrepvals) {
          set.color.groupby <- color.groupby
        } else if(length(color.groupby) != nrepvals) {
          set.color.groupby <- rep(color.groupby, nrepvals)
        }
      }
      if(linetype.groupby != 'NULL') {
        if(length(linetype.groupby) == nrepvals) {
          set.line.groupby <- linetype.groupby
        } else if(length(linetype.groupby) != nrepvals) {
          set.line.groupby <- rep(linetype.groupby, nrepvals)
        }
      }
      line.guide <- "none"
      color.guide <- "legend"
    }
  } 
  suppressMessages({
    plot <- plot + 
      ggplot2::scale_linetype_manual(values=set.line.groupby, 
                                     guide = line.guide) +
      ggplot2::scale_color_manual(values=set.color.groupby, 
                                  guide = color.guide)
  })
  return(plot)
} 







set_lines_colors_ribbon <- function(plot, guideby = NULL) {
  getbuiltingg <- ggplot2::ggplot_build(plot)
  get_line_  <- getbuiltingg$data[[1]]["linetype"]
  get_color_ <- getbuiltingg$data[[1]]["colour"]
  get_fill_  <- getbuiltingg$data[[1]]["colour"]
  ngrpanels  <- getbuiltingg$data[[1]]["group"]
  get_line_  <- unique(unlist(get_line_))
  get_color_ <- unique(unlist(get_color_))
  get_fill_  <- unique(unlist(get_fill_))
  ngrpanels <- length(unique(unlist(ngrpanels)))
  if(length(get_line_) != ngrpanels) get_line_ <- 
    rep(get_line_, ngrpanels)
  if(length(get_color_) != ngrpanels) get_color_ <- 
    rep(get_color_, ngrpanels)
  if(length(get_fill_) != ngrpanels) get_fill_ <- 
    rep(get_fill_, ngrpanels)
  setguide_line <- setguide_color <- setguide_fill <- 'none'
  if(is.null(guideby)) {
    setguide_line <- setguide_color <- setguide_fill <- 'none'
  } else if(guideby == 'line') {
    setguide_line <- 'legend'
  } else if(guideby == 'color') {
    setguide_color <- 'legend'
  } else if(guideby == 'fill') {
    setguide_fill <- 'legend'
  }
  suppressMessages({
    plot <- plot +
      ggplot2::scale_linetype_manual(values=get_line_, 
                                     guide = setguide_line) +
      ggplot2::scale_color_manual(values=get_color_, 
                                  guide = setguide_color) +
      ggplot2::scale_fill_manual(values=get_fill_, 
                                 guide = setguide_fill)
  })
  plot
}





#' An internal function to transform y axis when plotting with dual y axis
#' 
#' @noRd
#'
transform_sec_axis <- function(primary,
                               secondary,
                               na.rm = TRUE) {
  from <- range(secondary, na.rm = na.rm)
  to   <- range(primary, na.rm = na.rm)
  zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
    if (length(x) == 1) {
      return(TRUE)
    }
    if (length(x) != 2)
      stop("x must be length 1 or 2")
    if (any(is.na(x))) {
      return(NA)
    }
    if (x[1] == x[2]) {
      return(TRUE)
    }
    if (all(is.infinite(x))) {
      return(FALSE)
    }
    m <- min(abs(x))
    if (m == 0) {
      return(FALSE)
    }
    abs((x[1] - x[2]) / m) < tol
  }
  rescale.numeric_ <-
    function(x,
             to = c(0, 1),
             from = range(x, na.rm = TRUE, finite = TRUE),
             ...) {
      if (zero_range(from) || zero_range(to)) {
        return(ifelse(is.na(x), NA, mean(to)))
      }
      (x - from[1]) / diff(from) * diff(to) + to[1]
    }
  forward <- function(x) {
    rescale.numeric_(x, from = from, to = to)
  }
  reverse <- function(x) {
    rescale.numeric_(x, from = to, to = from)
  }
  list(fwd = forward, rev = reverse)
}




add_global_label <-
  function(pwobj,
           Xlab = NULL,
           Ylab = NULL,
           Xgap = 0.08,
           Ygap = 0.03,
           ...) {
    ylabgrob <- patchwork::plot_spacer()
    if (!is.null(Ylab)) {
      ylabgrob <- ggplot2::ggplot() +
        ggplot2::geom_text(ggplot2::aes(x = .5, y = .5),
                           label = Ylab,
                           angle = 90,
                           ...) +
        ggplot2::theme_void()
    }
    if (!is.null(Xlab)) {
      xlabgrob <- ggplot2::ggplot() +
        ggplot2::geom_text(ggplot2::aes(x = .5, y = .5), label = Xlab, ...) +
        ggplot2::theme_void()
    }
    if (!is.null(Ylab) & is.null(Xlab)) {
      return((ylabgrob + patchwork::patchworkGrob(pwobj)) +
               patchwork::plot_layout(widths = 100 * c(Ygap, 1 - Ygap))
      )
    }
    if (is.null(Ylab) & !is.null(Xlab)) {
      return((ylabgrob + pwobj) +
               (xlabgrob) +
               patchwork::plot_layout(
                 heights = 100 * c(1 - Xgap, Xgap),
                 widths = c(0, 100),
                 design = "
                                   AB
                                   CC
                                   "
               )
      )
    }
    if (!is.null(Ylab) & !is.null(Xlab)) {
      return((ylabgrob + pwobj) +
               (xlabgrob) +
               patchwork::plot_layout(
                 heights = 100 * c(1 - Xgap, Xgap),
                 widths = 100 * c(Ygap, 1 - Ygap),
                 design = "
                                   AB
                                   CC
                                   ")
             )
    }
    return(pwobj)
  }






line_color_key <- function(data, colorkey = NULL) {
  if(!is.null(colorkey)) {
    if(length(colorkey) == 1) {
      if(is.na(colorkey)) colorkey <- NULL
    }
  }
  if (is.null(colorkey)) {
    out <- NULL
  } else if (length(colorkey) == 1) {
    out <- data[[colorkey]]
  } else {
    cols_df <- data.frame(lapply(colorkey, function(name) data[[name]]))
    names(cols_df) <- colorkey
    out <- interaction(cols_df, drop = TRUE)
  }
  return(out)
}





max_decimals <- function(x) {
  parts <- strsplit(as.character(x), "\\.")
  decimals <- sapply(parts, function(p) if (length(p) == 2) nchar(p[2]) else 0L)
  max(decimals)
}



scale_x_continuous_fun <- function(x_minimum, x_maximum, by = 1, fun) {
  xseq <- seq(x_minimum, x_maximum, by)
  dec <- max_decimals(xseq)
  accuracy <- 10^-dec
  label_fun <- scales::label_number(accuracy = accuracy)
  body(fun) <- call("(", body(fun))
  ggplot2::scale_x_continuous(breaks = xseq, labels = fun(xseq) %>% 
                                label_fun())
}




get_fun_form <- function(xaxis_fun = NULL) {
  if (base::is.null(xaxis_fun)) {
    return(base::identity)
  }
  
  if (base::is.character(xaxis_fun)) {
    return(base::switch(
      xaxis_fun,
      identity = base::identity,
      log      = base::log,
      log10    = base::log10,
      exp      = base::exp,
      sqrt     = base::sqrt,
      square   = function(x) x^2,
      base::stop("Unknown xaxis_fun")
    ))
  }
  
  if (base::is.function(xaxis_fun)) {
    return(xaxis_fun)
  }
  
  base::stop("'xaxis_fun' must be NULL, character, or function.")
}

build_scale_x <- function(x_min = NULL,
                          x_max = NULL,
                          n = 5,
                          by = NULL,
                          xaxis_fun = NULL,
                          accuracy = 0.1,
                          expand_mult = 0,
                          use_coord_limits = TRUE,
                          coord_pad_mult = 0.05,
                          coord_pad_add = 0.0) {
  
  axis_fun <- get_fun_form(xaxis_fun)
  
  breaks_fun <- function(lims) {
    xmin <- if (base::is.null(x_min)) lims[1] else x_min
    xmax <- if (base::is.null(x_max)) lims[2] else x_max
    
    if (!base::is.finite(xmin) || !base::is.finite(xmax)) {
      return(base::numeric(0))
    }
    
    if (xmin == xmax) {
      return(xmin)
    }
    
    br <- if (!base::is.null(by)) {
      base::seq(from = xmin, to = xmax, by = by)
    } else {
      base::seq(from = xmin, to = xmax, length.out = n)
    }
    
    br[base::is.finite(br)]
  }
  
  label_fun <- function(x) {
    scales::label_number(accuracy = accuracy)(axis_fun(x))
  }
  
  scale_obj <- ggplot2::scale_x_continuous(
    limits = NULL,
    breaks = breaks_fun,
    labels = label_fun,
    expand = ggplot2::expansion(mult = expand_mult)
  )
  
  if (use_coord_limits && !(base::is.null(x_min) && base::is.null(x_max))) {
    xrng <- x_max - x_min
    pad <- xrng * coord_pad_mult + coord_pad_add
    
    base::list(
      scale_obj,
      ggplot2::coord_cartesian(
        xlim = base::c(x_min - pad, x_max + pad),
        expand = FALSE
      )
    )
  } else {
    scale_obj
  }
}



# build_scale_x <- function(x_min = NULL,
#                           x_max = NULL,
#                           n = 5,
#                           by = NULL,
#                           xaxis_fun = NULL,
#                           accuracy = 0.1,
#                           expand_mult = 0,
#                           use_coord_limits = TRUE,
#                           coord_pad_mult = 0.05,
#                           coord_pad_add = 0) {
#   
#   axis_fun <- get_fun_form(xaxis_fun)
#   
#   breaks_fun <- function(lims) {
#     xmin <- if (base::is.null(x_min)) lims[1] else x_min
#     xmax <- if (base::is.null(x_max)) lims[2] else x_max
#     
#     if (!base::is.finite(xmin) || !base::is.finite(xmax)) {
#       return(base::numeric(0))
#     }
#     
#     if (xmin == xmax) {
#       return(xmin)
#     }
#     
#     br <- if (!base::is.null(by)) {
#       base::seq(from = xmin, to = xmax, by = by)
#     } else {
#       base::seq(from = xmin, to = xmax, length.out = n)
#     }
#     
#     br[base::is.finite(br)]
#   }
#   
#   label_fun <- function(x) {
#     scales::label_number(accuracy = accuracy)(axis_fun(x))
#   }
#   
#   scale_obj <- ggplot2::scale_x_continuous(
#     limits = NULL,
#     breaks = breaks_fun,
#     labels = label_fun,
#     expand = ggplot2::expansion(mult = expand_mult)
#   )
#   
#   if (use_coord_limits && !(base::is.null(x_min) && base::is.null(x_max))) {
#     xrng <- x_max - x_min
#     pad <- xrng * coord_pad_mult + coord_pad_add
#     
#     base::list(
#       scale_obj,
#       ggplot2::coord_cartesian(
#         xlim = base::c(x_min - pad, x_max + pad),
#         expand = FALSE
#       )
#     )
#   } else {
#     scale_obj
#   }
# }


#######################################

check_unique_cap_opt <- function(opt) {
  chars <- strsplit(opt, "")[[1]]
  letters_lower <- tolower(chars)
  unique_letters <- unique(letters_lower)
  for (letter in unique_letters) {
    has_lower <- letter %in% chars
    has_upper <- toupper(letter) %in% chars
    if (has_lower && has_upper) {
      return(FALSE)  # Both cases found - invalid
    }
  }
  return(TRUE)
}

get_unique_opt_bands <- function(opt, bands, upper = FALSE) {
  if(is.null(bands)) bands <- ""
  opt_chars   <- strsplit(opt, "")  [[1]]  
  bands_chars <- strsplit(bands, "")[[1]] 
  if(upper) {
    opt_chars   <- opt_chars[opt_chars %in% LETTERS]
    bands_chars <- bands_chars[bands_chars %in% LETTERS]
  }
  if(is_emptyx(bands_chars)) bands_chars <- ""
  opt_chars <- unique(opt_chars)
  bands_chars <- unique(bands_chars)
  opt_chars <- paste0(opt_chars, collapse = "")
  bands_chars <- paste0(bands_chars, collapse = "")
  out <- list(opt = opt_chars, bands = bands_chars)
  return(out)
}


get_opt_bands <- function(opt, bands, upper = TRUE) {
  opt_chars   <- strsplit(opt, "")  [[1]] 
  bands_chars <- strsplit(bands, "")[[1]] 
  if(!is.null(upper)) {
    if(upper) {
      opt_chars   <- opt_chars[opt_chars %in% LETTERS]
      bands_chars <- bands_chars[bands_chars %in% LETTERS]
    } else {
      opt_chars   <- opt_chars[opt_chars %in% letters]
      bands_chars <- bands_chars[bands_chars %in% letters]
    }
  }
  opt_chars <- unique(opt_chars)
  bands_chars <- unique(bands_chars)
  bands_match_opt <- sapply(opt_chars, function(x) {
    idx <- match(x, bands_chars)
    if (is.na(idx)) "" else bands_chars[idx]
  })
  out <- list(opt = opt_chars, bands = unname(bands_match_opt))
  return(out)
}




loop_opt_bands <- function(opti, 
                           bandsi, 
                           arguments,
                           internal_formula_args,
                           ...) {
  yvar <- NULL;
  groupbytest <- NULL;
  subindicatorsi <- NULL;
  dy <- NULL;
  yvar <- NULL;
  subindicatorsi <- NULL;
  Estimate <- NULL;
  groupby <- NULL;
  groupby_line <- NULL;
  groupby_color <- NULL;
  Parameter <- NULL;
  cov_factor_vars <- NULL;
  Estimate.x <- NULL;
  groupby.x <- NULL;
  groupby_line.x <- NULL;
  groupby_color.x <- NULL;
  Estimate.y <- NULL;
  groupby.y <- NULL;
  groupby_line.y <- NULL;
  groupby_color.y <- NULL;
  groupby_fistr <- NULL;
  cov_vars <- NULL;
  ':=' <- NULL;
  . <- NULL;
  uvarby <- NULL;
  returndata <- NULL;
  difx <- NULL;
  dpar <- NULL;
  resp <- NULL;
  verbose <- NULL;
  itransform <- NULL;
  grid_add <- NULL;
  idvar <- NULL;
  numeric_cov_at <- NULL;
  aux_variables <- NULL;
  levels_id <- NULL;
  xrange <- NULL;
  idata_method <- NULL;
  newdata_fixed <- NULL;
  draw_ids <- NULL;
  groupby <- NULL;
  show_age_takeoff <- NULL;
  show_age_peak <- NULL;
  show_age_cessation <- NULL;
  show_vel_takeoff <- NULL;
  show_vel_peak <- NULL;
  show_vel_cessation <- NULL;
  linecolor <- NULL;
  linecolor1 <- NULL;
  linecolor2 <- NULL;
  label.y <- NULL;
  color.groupby <- NULL;
  xaxis_fun <- NULL;
  linetype.groupby <- NULL;
  band.legends <- NULL;
  conf <- NULL;
  robust <- NULL;
  envir <- NULL;
  trim <- NULL;
  estimation_method <- NULL;
  
  transform_xaxis <- NULL;
  addylab_dv <- NULL;
  estimation_method <- NULL;
  estimation_method <- NULL;
  estimation_method <- NULL;
  
  if (!is.null(internal_formula_args)) {
    eout <- list2env(internal_formula_args)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  opt <- opt.org <- opti
  bands <- tolower(bandsi)
  if (opt == 'd' | opt == 'D') {
    only_distance_curve <- TRUE
  } else {
    only_distance_curve <- FALSE
  }
  if (grepl("v", opt, ignore.case = F) |
      grepl("V", opt, ignore.case = F)) {
    need_velocity_curve <- TRUE
  } else {
    need_velocity_curve <- FALSE
  }
  if(only_distance_curve) {
    need_velocity_curve <- FALSE
  }
  if(need_velocity_curve) {
    need_xvar_must <- TRUE
  } else {
    need_xvar_must <- FALSE
  }
  if(returndata) {
    need_xvar_must <- need_xvar_must
  } else {
    need_xvar_must <- TRUE
  }
  arguments$opt <- opt
  arguments$ndraws <- ndraws
  arguments$draw_ids <- draw_ids

  if(is.null(difx)) difx <- xvar
  
  arguments$model$model_info[['difx']] <- arguments[['difx']] <- difx
  if(dpar == "sigma") {
    sigma_model <- get_sigmamodel_info(model = model,
                                       newdata = newdata,
                                       dpar = dpar, 
                                       resp = resp, 
                                       what = 'model',
                                       cov = cov, 
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
    
    if(is.null(grid_call)) {
      grid_call <- FALSE
    }
    
    if(is.null(grid_by)) grid_by <- cov
    
    if(sigma_model != "ls" && need_velocity_curve) {
      xvar <- check_set_xvar_sigma(model = model, 
                                   dpar = dpar, 
                                   xvar = xvar, 
                                   resp = resp, 
                                   auto = TRUE,
                                   verbose = verbose)
      
      if(grid_call)
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
                                       by = grid_by,
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
    } 
  } 
  
  assign_function_to_environment(transform_draws, 'transform_draws',
                                 envir = NULL)
  
  arguments$model$model_info[['transform_draws']] <-
    model$model_info[['transform_draws']] <- transform_draws
  get.newdata_args <- list()
  get.newdata_args[['model']]          <- model
  get.newdata_args[['newdata']]        <- newdata
  get.newdata_args[['xvar']]           <- xvar
  get.newdata_args[['idvar']]          <- idvar
  get.newdata_args[['resp']]           <- resp
  get.newdata_args[['numeric_cov_at']] <- numeric_cov_at
  get.newdata_args[['aux_variables']]  <- aux_variables
  get.newdata_args[['levels_id']]      <- levels_id
  get.newdata_args[['xrange']]         <- xrange
  get.newdata_args[['idata_method']]   <- idata_method
  get.newdata_args[['newdata_fixed']]  <- newdata_fixed
  get.newdata_args[['verbose']]        <- verbose
  get.newdata_args[['ipts']]           <- NULL
  get.newdata_args[['dpar']]           <- dpar
  get.newdata_args[['cov']]            <- cov
  newdata.xyadj              <- CustomDoCall(get.newdata, get.newdata_args)
  
  if(dpar == 'sigma') {
    if(!is.null(ipts)) {
      if(is.logical(ipts)) {
        if(ipts) {
          get.newdata_args[['ipts']]     <- ipts
          get.newdata_args[['newdata']]  <- newdata
          newdata <- CustomDoCall(get.newdata, get.newdata_args)
          attr(newdata, 'ipts_nocall') <- TRUE
        } else if(!ipts) {
          newdata <- newdata
          attr(newdata, 'ipts_nocall') <- FALSE
        }
      }
    } else {
      get.newdata_args[['ipts']]     <- ipts
      get.newdata_args[['newdata']]  <- newdata
      newdata <- CustomDoCall(get.newdata, get.newdata_args)
      attr(newdata, 'ipts_nocall') <- TRUE
    }
  } else if(dpar == 'mu') {
    get.newdata_args[['ipts']]     <- ipts
    get.newdata_args[['newdata']]  <- newdata
    newdata                    <- CustomDoCall(get.newdata, get.newdata_args)
    attr(newdata, 'ipts_nocall') <- TRUE
  }
  
  arguments$newdata          <- newdata
  
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
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  if (is.null(bands)) {
    bands <- ''
  }
  
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
    if(is.null(cov)) cov <- model$model_info[[cov_]] else cov <- cov
  } else if(dpar == "sigma") {
    if(!is.na(model$model_info[[sigmaxvar_]])) {
      xvar   <- model$model_info[[sigmaxvar_]]
    } else if(is.na(model$model_info[[sigmaxvar_]]) & 
              !is.null(model$model_info[[xvar_]])) {
      if(is.null(xvar)) {
        xvar   <- model$model_info[[xvar_]]
      }
    }
    if(is.null(cov)) cov <- model$model_info[[sigmacov_]] else cov <- cov
  } 
  
  groupvar_     <- paste0('groupvar', resp_rev_)
  yvar_         <- paste0('yvar', resp_rev_)
  yvar          <- model$model_info[[yvar_]]
  hierarchical_ <- paste0('hierarchical', resp_rev_)
  if(is.null(levels_id) & is.null(idvar)) {
    idvar <- model$model_info[[groupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      idvar <- model$model_info[[hierarchical_]]
    }
    model$model_info[[groupvar_]] <- idvar # idvar[1]
  } else if (!is.null(levels_id)) {
    idvar <- levels_id
  } else if (!is.null(idvar)) {
    idvar <- idvar
  }
  if(is.null(idvar)) {
    if(is.null(idvar)) {
      if(!is.null(model$model_info[['idvars']])) {
        idvar <- model$model_info[['idvars']]
      }
    }
  }
  
  Xx <- xvar
  Yy <- yvar
  
  if (grepl("p", bands, ignore.case = T) & summary) {
    stop2c(
      "To construct bands (e.g., 95%) around the parameter estimates",
      "\n ",
      " (such as APGV, PGV), they are first calculated for each",
      "\n ",
      " posterior draw and then summarised for the give conf limit.",
      "\n ",
      " Therefore,summary option must be set to FALSE"
    )
  }
  if (grepl("a", bands, ignore.case = T) & summary) {
    stop2c(
      "To construct bands (e.g., 95%) around the adjusted curve estimates, ",
      "\n ",
      " the summary option must be set to FALSE"
    )
  }
  if (grepl("a", opt, ignore.case = F) |
      grepl("u", opt, ignore.case = F)) {
    ipts <- NULL
    if(verbose) {
      message("The ipts has been set to NULL i.e., ipts = NULL",
              "\n ",
              "This because it does't a make sense to interploate data when",
              "\n ",
              " estimating adjusted/unadjusted curves")
    }
    testdata1 <- model$data %>% dplyr::select(dplyr::all_of(idvar)) %>% 
      droplevels() %>% 
      dplyr::mutate(
        groupbytest = interaction(dplyr::across(dplyr::all_of(idvar)))
      ) %>% 
      dplyr::select(dplyr::all_of(groupbytest)) %>% dplyr::ungroup()
    testdata2 <- newdata %>% dplyr::select(dplyr::all_of(idvar)) %>% 
      droplevels() %>% 
      dplyr::mutate(groupbytest = 
                      interaction(dplyr::across(dplyr::all_of(idvar)))) %>% 
      dplyr::select(dplyr::all_of(groupbytest)) %>% dplyr::ungroup()
  }
  pv <- FALSE
  if (returndata & nchar(opt) > 1) {
    stop2c(
      "For returndata, please specify only one option at a time",
      "\n ",
      " (out of the total six optiona available, i.e., dvDVau)",
      "\n ",
      " For example, opt = 'd'"
    )
  }
  if (!grepl("v", opt, ignore.case = F) &
      !grepl("V", opt, ignore.case = F)) {
    apv <- arguments$apv <- FALSE
    pv <- arguments$pv <- FALSE
  }
  if(is.null(arguments$pv)) arguments$pv <- FALSE
  if(length(list(...)) != 0) arguments <- c(arguments, list(...))
  arguments$draw_ids <- draw_ids
  for (i in names(arguments)) {
    if(is.symbol(arguments[[i]])) {
      if(deparse(arguments[[i]]) == "") {
        arguments[[i]] <- NULL
      }
    }
  }
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
  d. <- CustomDoCall(growthparameters, arguments)
  if(is.null(d.)) return(invisible(NULL))
  p.                    <- d.[['parameters']]
  probtitles            <- d.[['probtitles']]
  groupby_str_d         <- d.[['groupby_str_d']]
  groupby_str_v         <- d.[['groupby_str_v']]
  p.as.d.out_attr       <- p.
  d.[['parameters']]    <- NULL
  d.[['probtitles']]    <- NULL
  d.[['groupby_str_d']] <- NULL
  d.[['groupby_str_v']] <- NULL
  groupby_str_d <- unique(c(groupby_str_d, groupby))
  groupby_str_v <- unique(c(groupby_str_v, groupby))
  d. <- d. %>% CustomDoCall(rbind, .) %>% data.frame()
  row.names(d.) <- NULL
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
    d. <- prepare_transformations(data = d., model = model,
                                  itransform = itransform_set)
    newdata <- prepare_transformations(data = newdata, model = model,
                                       itransform = itransform_set)
  }
  
  curve.d <- 'distance'
  curve.v <- 'velocity'
  name.apv <- "APGV"
  name.pv <- "PGV"
  name.atv <- "ATGV"
  name.tv <- "TGV"
  name.acv <- "ACGV"
  name.cv <- "CGV"
  name.vline <- c()
  if(show_age_takeoff)   name.vline <- c(name.vline, name.atv)
  if(show_age_peak)      name.vline <- c(name.vline, name.apv)
  if(show_age_cessation) name.vline <- c(name.vline, name.acv)
  name.hline <- c()
  if(show_vel_takeoff)   name.hline <- c(name.hline, name.tv)
  if(show_vel_peak)      name.hline <- c(name.hline, name.pv)
  if(show_vel_cessation) name.hline <- c(name.hline, name.cv)
  name.hline <- c()
  x_minimum <- min(newdata[[Xx]])
  x_maximum <- max(newdata[[Xx]])
  x_minimum <- floor(x_minimum)
  x_maximum <- ceiling(x_maximum)
  single_plot_pair_color_dv_au <- c('black', 'red')
  if (nchar(opt) > 2) {
    layout <- 'facet'
  }
  if (is.null(linecolor)) {
    color_single <- "grey50"
  }
  if (is.null(linecolor1)) {
    color.d <- "orange2"
    color.adj <- "orange2"
  }
  if (is.null(linecolor2)) {
    color.v <- "green4"
    color.unadj <- "green4"
  }
  if (is.null(label.y)) {
    label.d     <- firstup(curve.d)
    label.v     <- firstup(curve.v)
    label.adj   <- firstup('adjusted')
    label.unadj <- firstup('unadjusted')
  } else {
    label.d     <- label.v     <- label.y
    label.adj   <- label.unadj <- label.y
  }
  # Note that labs y = "", change them to addylab_d / addylab_v if needed
  y_lab_d <- "" # addylab_d
  y_lab_v <- "" # addylab_v
  y_lab_o <- "" # "Observed Trajectories"
  
  if (is.null(label.x)) {
    label.x     <- paste0(Xx, "")
  }
  if (is.null(legendpos)) {
    if (grepl("D", opt, ignore.case = F) |
        grepl("V", opt, ignore.case = F)) {
      legendpos <- "none"
    } else {
      legendpos <- "bottom"
      legendpos.adj.unadj <- "topleft"
    }
  } else if (!is.null(legendpos)) {
    legendpos <- legendpos.adj.unadj <- legendpos
  }
  if (is.null(linetype.apv)) {
    linetype.apv <- 'dotted'
    linetype.pv <- 'dotted'
  } else {
    linetype.apv <- linetype.pv <- linetype.apv
  }
  if (is.null(linewidth.main)) {
    linewidth.main <- 0.5
  }
  if (is.null(linewidth.apv)) {
    linewidth.apv <- 0.5
    linewidth.pv <- 0.5
  } else {
    linewidth.apv <- linewidth.pv <- linewidth.apv
  }
  if (is.null(band.alpha)) {
    band.alpha <- 0.25
  }
  if(is.null(color.groupby)) {
    if(dpar == "mu")    color.groupby <- cov 
    if(dpar == "sigma") color.groupby <- cov 
  }
  
  if(is.null(fill.groupby)) {
    fill.groupby <- color.groupby
  } else if(!is.null(fill.groupby)) {
    if(isFALSE(fill.groupby)) fill.groupby <- NA
  }
  
  
  build_scale_x_args <- list()
  for (i in names(formals(build_scale_x))) {
    build_scale_x_args[[i]] <- transform_xaxis[[i]]
  }
  if(is.null(build_scale_x_args[["x_min"]])) {
    build_scale_x_args[["x_min"]] <- x_minimum
  }
  if(is.null(build_scale_x_args[["x_max"]])) {
    build_scale_x_args[["x_max"]] <- x_maximum
  }
  
  
  defaults <- base::as.list(base::formals(build_scale_x))
  build_scale_x_args <- utils::modifyList(defaults, build_scale_x_args)
  add_build_scale_x <- base::do.call(build_scale_x, build_scale_x_args)
  
  
  if (grepl("O", opt.org, ignore.case = T) |
      grepl("O", opt.org, ignore.case = T)) {
    groupby_str_d.o <- idvar
    curves <- 'xxxxx' 
    if (length(curves) == 1) {
      layout <- 'facet'
    } else {
      layout <- layout
    }
    d.o <- model[['data']]
    if (layout == 'facet')
      color.d <- color.v <- color_single
    
    if (grepl("O", opt.org, ignore.case = T)) {
      d.o <- d.o
      index_opt <- gregexpr("O", opt.org, ignore.case = T)[[1]]
      dist.. <- substr(opt.org, index_opt, index_opt)
      if (grepl("^[[:upper:]]+$", dist..)) {
        d.o <-
          d.o %>%
          dplyr::mutate(
            groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d.o)))
          )
      } else if (!grepl("^[[:upper:]]+$", dist..)) {
        if (is.null(groupby_str_d.o))
          d.o <- d.o %>% dplyr::mutate(groupby = NA)
        if (!is.null(groupby_str_d.o))
          d.o <-
            d.o %>%
            dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d.o)))
            )
      }
      if(is.na(d.o[['groupby']][1])) {
        d.o$groupby_line  <- 'solid'
        d.o$groupby_color <- 'black'
      } else {
        d.o$groupby_line  <- d.o$groupby
        d.o$groupby_color <- d.o$groupby
      }
      d.o[[Xx]] <- ifunx_(d.o[[Xx]])
      plot.o.O <- d.o %>% 
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = .data[[Yy]],
            group = groupby,
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, color.groupby)
          ),
          linewidth = linewidth.main
        ) +
        add_build_scale_x + 
        ggplot2::labs(x = "", y = y_lab_o, title = "Observed") +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
  }
  
  if (grepl("d", opt, ignore.case = T) |
      grepl("v", opt, ignore.case = T)) {
    curves <- unique(d.$curve)
    if (length(curves) == 1) {
      layout <- 'facet'
    } else {
      layout <- layout
    }
    curves <- 'xxxxx' # placeholder for setting curves == 1
    if (layout == 'facet')
      color.d <- color.v <- color_single
    if (grepl("d", opt, ignore.case = T)) {
      d.o <- d.
      index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
      dist.. <- substr(opt, index_opt, index_opt)
      if (grepl("^[[:upper:]]+$", dist..)) {
        d. <-
          d. %>% 
          dplyr::mutate(
            groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
          )
      } else if (!grepl("^[[:upper:]]+$", dist..)) {
        if (is.null(groupby_str_d))
          d. <- d. %>% dplyr::mutate(groupby = NA)
        if (!is.null(groupby_str_d))
          d. <-
            d. %>% 
            dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
            )
      }
      if(is.na(d.[['groupby']][1])) {
        d.$groupby_line  <- 'solid'
        d.$groupby_color <- 'black'
      } else {
        d.$groupby_line  <- d.$groupby
        d.$groupby_color <- d.$groupby
      }
      plot.o.d <- d. %>% dplyr::filter(curve == curve.d) %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate,
            group = groupby,
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, color.groupby)
          ),
          linewidth = linewidth.main
        ) +
        add_build_scale_x + 
        ggplot2::labs(x = "", y = y_lab_d, title = label.d) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      if (grepl("d", bands, ignore.case = T)) {
        plot.o.d <- plot.o.d +
          ggplot2::geom_ribbon(
            data = d. %>% dplyr::filter(curve == curve.d),
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '')]],
              ymax = .data[[paste0(probtitles[2], '')]],
              group = groupby,
              linetype = line_color_key(.data, linetype.groupby),
              color = line_color_key(.data, color.groupby),
              fill = line_color_key(.data, fill.groupby)
            ),
            alpha = band.alpha, show.legend = band.legends
          )
        if(!band.legends) plot.o.d <- plot.o.d + ggplot2::guides(fill = "none")
      }
      d. <- d.o
      if ('curve' %in% names(d.)) {
        d.out <- d. %>% dplyr::select(-dplyr::all_of('curve')) # curve to 'curve'
      } else {
        d.out <- d.
      }
    }
    if (!grepl("d", opt, ignore.case = T)) {
      plot.o.d <- NULL
    }
    if (grepl("v", opt, ignore.case = T)) {
      d.o <- d.
      index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
      velc.. <- substr(opt, index_opt, index_opt)
      if (grepl("^[[:upper:]]+$", velc..)) {
        d. <-
          d. %>% 
          dplyr::mutate(
            groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
          )
      } else if (!grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_v))
          d. <- d. %>% dplyr::mutate(groupby = NA)
        if (!is.null(groupby_str_v))
          d. <-
            d. %>% dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
            )
      }
      if(is.na(d.[['groupby']][1])) {
        d.$groupby_line <- 'solid'
        d.$groupby_color <- 'black'
      } else {
        d.$groupby_line <- d.$groupby
        d.$groupby_color <- d.$groupby
      }
      plot.o.v <- d. %>% dplyr::filter(curve == curve.v) %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate,
            group = groupby,
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, color.groupby)
          ),
          linewidth = linewidth.main
        ) +
        add_build_scale_x + 
        ggplot2::labs(x = "", y = y_lab_v, title = label.v) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      if (grepl("v", bands, ignore.case = T)) {
        plot.o.v <- plot.o.v +
          ggplot2::geom_ribbon(
            data = d. %>% dplyr::filter(curve == curve.v),
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '')]],
              ymax = .data[[paste0(probtitles[2], '')]],
              group = groupby,
              linetype = line_color_key(.data, linetype.groupby),
              color = line_color_key(.data, color.groupby),
              fill = line_color_key(.data, fill.groupby)
            ),
            alpha = band.alpha, show.legend = band.legends
          )
        if(!band.legends) plot.o.v <- plot.o.v + ggplot2::guides(fill = "none")
      }
      if (pv) {
        data_hline <- p. %>% dplyr::filter(Parameter == name.pv)
        plot.o.v <- plot.o.v +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = .data[['Estimate']]),
            linewidth = linewidth.pv,
            linetype = linetype.pv
          )
        if (grepl("p", bands, ignore.case = T)) {
          plot.o.v <- plot.o.v +
            ggplot2::annotate(
              "rect",
              xmin = -Inf,
              xmax = Inf,
              ymin = (data_hline[[paste0(probtitles[1], '')]]),
              ymax = (data_hline[[paste0(probtitles[2], '')]]),
              alpha = band.alpha
            )
        }
      }
      if (!is.null(name.vline) & !is.null(p.)  ) {
        data_vline <- 
          dplyr::filter(p., grepl(paste(name.vline, collapse  = "|"),
                                  Parameter))
        plot.o.v <- plot.o.v +
          ggplot2::geom_vline(
            data = data_vline,
            ggplot2::aes(xintercept = .data[['Estimate']]),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        if (grepl("p", bands, ignore.case = T)) {
          plot.o.v <- plot.o.v +
            ggplot2::annotate(
              "rect",
              xmin = (data_vline[[paste0(probtitles[1], '')]]),
              xmax = (data_vline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      } 
      if (!is.null(name.hline) & !is.null(p.) ) {
        data_hline <- 
          dplyr::filter(p., grepl(paste(name.hline, collapse  = "|"), 
                                  Parameter))
        plot.o.v <- plot.o.v +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = t.s.axis$fwd(.data[['Estimate']]) ),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        if (grepl("p", bands, ignore.case = T)) {
          plot.o.v <- plot.o.v +
            ggplot2::annotate(
              "rect",
              xmin = (data_vline[[paste0(probtitles[1], '')]]),
              xmax = (data_vline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      } 
      d. <- d.o
      if ('curve' %in% names(d.)) {
        d.out <- d. %>% dplyr::select(-dplyr::all_of('curve'))
      } else {
        d.out <- d.
      }
    }
    if (!grepl("v", opt, ignore.case = T)) {
      plot.o.v <- NULL
    }
    if (length(curves) > 1 & layout == 'facet') {
      plot.o <- patchwork::wrap_plots(plot.o.d + plot.o.v) %>%
        add_global_label(
          Xlab = label.x,
          Ylab = "",
          size = 5,
          Xgap = 0.08,
          Ygap = 0.04
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = legendpos,
                       legend.direction = 'horizontal')
      
    } else if (length(curves) == 1 & layout == 'facet') {
      if (!is.null(plot.o.d)) {
        plot.o <- plot.o.d +
          ggplot2::labs(x = label.x, y = label.d) +
          ggplot2::theme(plot.title = ggplot2::element_blank()) +
          ggplot2::theme(legend.position = legendpos,
                         legend.direction = 'horizontal')
      } else if (!is.null(plot.o.v)) {
        plot.o <- plot.o.v +
          ggplot2::labs(x = label.x, y = label.v) +
          ggplot2::theme(plot.title = ggplot2::element_blank()) +
          ggplot2::theme(legend.position = legendpos,
                         legend.direction = 'horizontal')
      }
    }
    if (length(curves) > 1 & layout == 'single') {
      data_d <- subset(d., curve == "distance")
      data_v <- subset(d., curve == "velocity")
      by_join_ <- c(idvar, Xx, groupby_str_d)
      by_join_ <- unique(by_join_)
      data_dv <- dplyr::left_join(data_d, data_v, by = by_join_)
      data_dv.o <- data_dv
      if (grepl("d", opt, ignore.case = T)) {
        index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
        dist.. <- substr(opt, index_opt, index_opt)
        if (grepl("^[[:upper:]]+$", dist..)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
            ) %>%
            dplyr::mutate(
              groupby.x = 
                interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
            )
        } else if (!grepl("^[[:upper:]]+$", dist..)) {
          if (is.null(groupby_str_d)) {
            data_dv <- data_dv %>% dplyr::mutate(groupby = NA) %>%
              dplyr::mutate(groupby.x = NA)
          } else if (!is.null(groupby_str_d)) {
            data_dv <- data_dv %>%
              dplyr::mutate(
                groupby = 
                  interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
              ) %>%
              dplyr::mutate(
                groupby.x = 
                  interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
              )
          }
        }
      }
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
        if (grepl("^[[:upper:]]+$", velc..)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
            ) %>%
            dplyr::mutate(groupby.y = groupby)
        } else if (!grepl("^[[:upper:]]+$", velc..)) {
          if (is.null(groupby_str_v)) {
            data_dv <- data_dv %>% dplyr::mutate(groupby = NA) %>%
              dplyr::mutate(groupby.y = NA)
          } else if (!is.null(groupby_str_v)) {
            data_dv <- data_dv %>%
              dplyr::mutate(
                groupby = 
                  interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
              ) %>%
              dplyr::mutate(
                groupby.y = 
                  interaction(dplyr::across(dplyr::all_of(groupby_str_v)))
              )
          }
        }
      }
      if (grepl("^[[:upper:]]+$", dist..) &
          !grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_v)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby.x = 
                interaction(dplyr::across(dplyr::all_of(groupby_str_d))),
              groupby.y = NA)
        } else if (!is.null(groupby_str_v)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby.x = 
                interaction(dplyr::across(dplyr::all_of(groupby_str_d))),
              groupby.y =
                interaction(dplyr::across(dplyr::all_of(groupby_str_v))))
        }
      }
      if (!grepl("^[[:upper:]]+$", dist..) &
          grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_d)) {
          data_dv <-
            data_dv %>% dplyr::mutate(
              groupby.x = NA,
              groupby.y =
                interaction(dplyr::across(dplyr::all_of(groupby_str_v))))
        } else if (!is.null(groupby_str_d)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(
              groupby.x = 
                interaction(dplyr::across(dplyr::all_of(groupby_str_d))),
              groupby.y =
                interaction(dplyr::across(dplyr::all_of(groupby_str_v))))
        }
      }
      t.s.axis <- with(data_dv, transform_sec_axis(Estimate.x, Estimate.y))
      if(is.na(uvarby)) {
        if(is.na(data_dv[['groupby.x']][1])) {
          legendlabs_mult_singel <- addylab_dv # c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          data_dv$groupby_line.x <- legendlabs_mult_singel[1]
          data_dv$groupby_color.x <- legendlabs_mult_singel[1]
          data_dv$groupby_line.y <- legendlabs_mult_singel[2]
          data_dv$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          data_dv$groupby_line.x <- data_dv$groupby.x
          data_dv$groupby_color.x <- data_dv$groupby.x
          data_dv$groupby_line.y <- data_dv$groupby.y
          data_dv$groupby_color.y <- data_dv$groupby.y
          legendlabs_mult_mult <- unique(data_dv[['groupby.x']])
        }
      }
      if(!is.na(uvarby)) {
        if(is.null(cov_factor_vars)) {
          legendlabs_mult_singel <- addylab_dv # c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          data_dv$groupby_line.x <- legendlabs_mult_singel[1]
          data_dv$groupby_color.x <- legendlabs_mult_singel[1]
          data_dv$groupby_line.y <- legendlabs_mult_singel[2]
          data_dv$groupby_color.y <- legendlabs_mult_singel[2]
          legendlabs_mult_mult <- NULL # unique(data_dv[['groupby.x']])
        } else {
          data_dv$groupby_line.x <- data_dv$groupby.x
          data_dv$groupby_color.x <- data_dv$groupby.x
          data_dv$groupby_line.y <- data_dv$groupby.y
          data_dv$groupby_color.y <- data_dv$groupby.y
          legendlabs_mult_mult <- unique(data_dv[['groupby.x']])
        }
      }
      if(length( unique(round(data_dv$Estimate.y, 6))) == 1) {
        stop2c("The velocity estimates are identical over the entire range of x",
               "\n  ", 
               "Therefore, can't draw distance and velocity curves together")
      }
      plot.o <- data_dv %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate.x,
            group = groupby.x,
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, color.groupby)
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = t.s.axis$fwd(Estimate.y),
            group = groupby.y,
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, color.groupby)
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::scale_y_continuous(sec.axis =
                                      ggplot2::sec_axis(~ t.s.axis$rev(.),
                                                        name = label.v)) +
        ggplot2::labs(x = label.x, y = label.d, title  = "") +
        add_build_scale_x + 
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90))
      getbuiltingg <- ggplot2::ggplot_build(plot.o)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- 
        rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- 
        rep(get_color_, ngrpanels)
      if(!exists('legendlabs_mult_line')) legendlabs_mult_line <- 'solid'
      if(!exists('legendlabs_mult_color')) legendlabs_mult_color <- 'black'
      if(!exists('legendlabs_mult_singel')) legendlabs_mult_singel <- 'solid'
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- legendlabs_mult_line
        get_color_ <- legendlabs_mult_color
        legendlabs_ <- legendlabs_mult_singel
      }
      plot.o <- plot.o + 
        ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
        ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
      if (grepl("d", bands, ignore.case = T)) {
        plot.o <- plot.o +
          ggplot2::geom_ribbon(
            data = data_dv,
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '.x')]],
              ymax = .data[[paste0(probtitles[2], '.x')]],
              group = groupby.x,
              linetype = line_color_key(.data, linetype.groupby),
              color = line_color_key(.data, color.groupby),
              fill = line_color_key(.data, fill.groupby)
            ),
            alpha = band.alpha, show.legend = band.legends
          )
        if(!band.legends) plot.o <- plot.o + ggplot2::guides(fill = "none")
      }
      if (grepl("v", bands, ignore.case = T)) {
        plot.o <- plot.o +
          ggplot2::geom_ribbon(
            data = data_dv,
            ggplot2::aes(
              ymin = t.s.axis$fwd(.data[[paste0(probtitles[1], '.y')]]),
              ymax = t.s.axis$fwd(.data[[paste0(probtitles[2], '.y')]]),
              group = groupby.y,
              linetype = line_color_key(.data, linetype.groupby),
              color = line_color_key(.data, color.groupby),
              fill = line_color_key(.data, fill.groupby)
            ),
            alpha = band.alpha, show.legend = band.legends
          )
        if(!band.legends) plot.o <- plot.o + ggplot2::guides(fill = "none")
      }
      if((grepl("d", bands, ignore.case = T) & 
          !grepl("v", bands, ignore.case = T)) |
         !grepl("d", bands, ignore.case = T) & 
         grepl("v", bands, ignore.case = T)
      ) {
        one_band <- TRUE
      } else {
        one_band <- FALSE
      }
      if(one_band & ngrpanels == 1) {
        if(grepl("d", bands, ignore.case = T) & 
           !grepl("v", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::scale_fill_manual(values=legendlabs_mult_color[1], 
                                       guide = 'none')
        }
        if(!grepl("d", bands, ignore.case = T) & 
           grepl("v", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::scale_fill_manual(values=legendlabs_mult_color[2], 
                                       guide = 'none')
        }
      }
      if (pv) {
        data_hline <- p. %>% dplyr::filter(Parameter == name.pv)
        plot.o <- plot.o +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = t.s.axis$fwd(.data[['Estimate']])),
            linewidth = linewidth.pv,
            linetype = linetype.pv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::annotate(
              "rect",
              xmin = -Inf,
              xmax = Inf,
              ymin = t.s.axis$fwd(data_hline[[paste0(probtitles[1], '')]]),
              ymax = t.s.axis$fwd(data_hline[[paste0(probtitles[2], '')]]),
              alpha = band.alpha
            )
        }
      }
      if (!is.null(name.vline) & !is.null(p.)) {
        data_vline <- dplyr::filter(p., grepl(paste(name.vline, 
                                                    collapse  = "|"), 
                                              Parameter))
        plot.o <- plot.o +
          ggplot2::geom_vline(
            data = data_vline,
            ggplot2::aes(xintercept = .data[['Estimate']]),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        if (grepl("p", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::annotate(
              "rect",
              xmin = (data_vline[[paste0(probtitles[1], '')]]),
              xmax = (data_vline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      }
      if (!is.null(name.hline) & !is.null(p.) ) {
        data_hline <- dplyr::filter(p., grepl(paste(name.hline, 
                                                    collapse  = "|"), 
                                              Parameter))
        plot.o <- plot.o +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = t.s.axis$fwd(.data[['Estimate']]) ),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        if (grepl("p", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::annotate(
              "rect",
              xmin = (data_hline[[paste0(probtitles[1], '')]]),
              xmax = (data_hline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      } 
      data_dv <- data_dv.o
      if ('curve' %in% names(data_dv)) {
        d.out <- data_dv %>% dplyr::select(-dplyr::all_of('curve'))
      } else {
        d.out <- data_dv
      }
    }
  }
  groupby_str_au <- groupby_fistr
  if (grepl("a", opt, ignore.case = T) |
      grepl("u", opt, ignore.case = T)) {
    if (grepl("a", opt, ignore.case = T)) {
      xyadj_ed <- xyadj_curves(model, 
                               x = NULL,
                               y = NULL,
                               id = NULL,
                               v = NULL,
                               newdata = newdata.xyadj, 
                               ndraws = ndraws,
                               draw_ids = draw_ids,
                               resp = resp, 
                               tomean = TRUE,
                               conf = conf, 
                               robust = robust,
                               summary = summary, 
                               numeric_cov_at = numeric_cov_at,
                               aux_variables = aux_variables,
                               levels_id = levels_id,
                               ipts = ipts,
                               xrange = xrange, 
                               idata_method = idata_method,
                               verbose = verbose,
                               model_deriv = NULL,
                               deriv = NULL, 
                               envir = envir,
                               ...) 
      out_a_ <- trimlines_curves(model, 
                                 x = Xx,
                                 y = Yy,
                                 id = idvar,
                                 newdata = xyadj_ed, 
                                 ndraws = ndraws,
                                 draw_ids = draw_ids,
                                 resp = resp, 
                                 level = 0,
                                 trim = trim, 
                                 estimation_method = estimation_method,
                                 verbose = verbose,
                                 model_deriv = NULL,
                                 deriv = NULL, 
                                 envir = envir,
                                 ...)
      if(any(itransform_set != "")) {
        out_a_ <- prepare_transformations(data = out_a_, model = model,
                                          itransform = itransform_set)
      }
      d.out <- out_a_
      dots <- list(...)
      set_get_dv <- FALSE
      if(!is.null(dots$get_dv)) {
        if(dots$get_dv) {
          if(verbose) message("executing 'get_dv'!")
          set_get_dv <- TRUE
        }
      }
      if(set_get_dv) {
        return(out_a_)
      }
      if(!is.null(dots$xadj_tmt)) {
        if(dots$xadj_tmt) {
          return(out_a_)
        }
      }
      if(!is.null(dots$xadj_tmf)) {
        if(dots$xadj_tmf) {
          return(out_a_)
        }
      }
      out_a_ <-
        out_a_ %>%
        dplyr::mutate(
          groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_au)))
        )
      

      build_scale_x_args <- list()
      for (i in names(formals(build_scale_x))) {
        build_scale_x_args[[i]] <- transform_xaxis[[i]]
      }
      if(is.null(build_scale_x_args[["x_min"]])) {
        build_scale_x_args[["x_min"]] <- x_minimum
      }
      if(is.null(build_scale_x_args[["x_max"]])) {
        build_scale_x_args[["x_max"]] <- x_maximum
      }
      
      
      defaults <- base::as.list(base::formals(build_scale_x))
      build_scale_x_args <- utils::modifyList(defaults, build_scale_x_args)
      add_build_scale_x <- base::do.call(build_scale_x, build_scale_x_args)
      
      
      out_a_ <- out_a_[out_a_[[Xx]] >= x_minimum & 
                         out_a_[[Xx]] <= x_maximum, ]
      out_a_ <- out_a_ %>% dplyr::mutate(groupby.x = groupby, 
                                         groupby.y = groupby.x)
      if(is.na(uvarby)) {
        if(is.na(out_a_[['groupby']][1])) {
          legendlabs_mult_singel <- addylab_dv # c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_a_$groupby_line.x <- legendlabs_mult_singel[1]
          out_a_$groupby_color.x <- legendlabs_mult_singel[1]
          out_a_$groupby_line.y <- legendlabs_mult_singel[2]
          out_a_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_a_$groupby_line.x <- out_a_$groupby.x
          out_a_$groupby_color.x <- out_a_$groupby.x
          out_a_$groupby_line.y <- out_a_$groupby.y
          out_a_$groupby_color.y <- out_a_$groupby.y
          legendlabs_mult_mult <- unique(out_a_[['groupby']])
        }
      }
      if(!is.na(uvarby)) {
        if(is.null(cov_factor_vars)) {
          legendlabs_mult_singel <- addylab_dv # c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_a_$groupby_line.x <- legendlabs_mult_singel[1]
          out_a_$groupby_color.x <- legendlabs_mult_singel[1]
          out_a_$groupby_line.y <- legendlabs_mult_singel[2]
          out_a_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_a_$groupby_line.x <- out_a_$groupby.x
          out_a_$groupby_color.x <- out_a_$groupby.x
          out_a_$groupby_line.y <- out_a_$groupby.y
          out_a_$groupby_color.y <- out_a_$groupby.y
          legendlabs_mult_mult <- unique(out_a_[['groupby']])
        }
      }
      plot.o.a <- out_a_ %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = !!as.name(Yy),
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, color.groupby),
            group = groupby.x # ,
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::labs(x = label.x, y = label.d, title  = "") +
        add_build_scale_x + 
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(y = paste0("Adjusted ", "Individual Curves")) +
        ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90)) +
        ggplot2::labs(title = label.adj) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      getbuiltingg <- ggplot2::ggplot_build(plot.o.a)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- 
        rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- 
        rep(get_color_, ngrpanels)
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- 'solid' # legendlabs_mult_line
        get_color_ <- 'black' # legendlabs_mult_color
        legendlabs_ <- NULL # legendlabs_mult_singel
      }
      if (grepl("a", bands, ignore.case = T)) {
        plot.o.a <- plot.o.a +
          ggplot2::geom_ribbon(
            data = out_a_,
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '')]],
              ymax = .data[[paste0(probtitles[2], '')]],
              group = groupby.x,
              linetype = line_color_key(.data, linetype.groupby),
              color = line_color_key(.data, color.groupby),
              fill = line_color_key(.data, fill.groupby)
            ),
            alpha = band.alpha, show.legend = band.legends
          )
        if(!band.legends) plot.o.a <- plot.o.a + ggplot2::guides(fill = "none")
      }
      if (nchar(opt) == 1) {
        plot.o.a <- plot.o.a +
          ggplot2::theme(plot.title = ggplot2::element_blank())
      } else if (nchar(opt) > 1) {
        plot.o.a <- plot.o.a +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
      if (!grepl("u", opt, ignore.case = T)) {
        suppressMessages({
        })
      }
    } else if (!grepl("a", opt, ignore.case = T)) {
      plot.o.a <- NULL
    }
    if (grepl("u", opt, ignore.case = T)) {
      xyunadj_ed <- xyunadj_curves(model, 
                                   x = NULL,
                                   y = NULL,
                                   id = NULL,
                                   newdata = NULL,
                                   ndraws = ndraws,
                                   draw_ids = draw_ids,
                                   resp = resp, 
                                   verbose = verbose,
                                   model_deriv = NULL,
                                   deriv = NULL, 
                                   envir = envir,
                                   ...)
      out_u_ <- trimlines_curves(model, 
                                 x = Xx,
                                 y = Yy,
                                 id = idvar,
                                 newdata = xyunadj_ed, 
                                 ndraws = ndraws,
                                 draw_ids = draw_ids,
                                 resp = resp, 
                                 level = 0,
                                 trim = trim, 
                                 estimation_method = estimation_method,
                                 verbose = verbose,
                                 model_deriv = NULL,
                                 deriv = NULL, 
                                 envir = envir,
                                 ...)
      if(any(itransform_set != "")) {
        out_u_ <- prepare_transformations(data = out_u_, model = model,
                                          itransform = itransform_set)
      }
      d.out <- out_u_
      groupby_cov <- setdiff(groupby_str_au, idvar)
      out_u_x_unique <- out_u_ %>% 
        dplyr:: distinct(!!as.name(idvar), .keep_all = F)
      zzz_unique     <- model$data %>% 
        dplyr::select(dplyr::all_of(c(idvar, groupby_cov))) %>% 
        dplyr:: distinct(!!as.name(idvar), .keep_all = TRUE)
      out_u_x_unique <- out_u_x_unique %>% 
        dplyr:: left_join(zzz_unique, by = idvar, keep = F)
      out_u_ <- dplyr::left_join(out_u_, out_u_x_unique %>% 
                                   dplyr::select(dplyr::all_of(c(idvar, groupby_cov))), by = idvar, keep = F)
      out_u_ <-
        out_u_ %>%
        dplyr::mutate(
          groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_au)))
        )
      out_u_ <- out_u_ %>% dplyr::mutate(groupby.x = groupby, 
                                         groupby.y = groupby.x)
      

      
      build_scale_x_args <- list()
      for (i in names(formals(build_scale_x))) {
        build_scale_x_args[[i]] <- transform_xaxis[[i]]
      }
      if(is.null(build_scale_x_args[["x_min"]])) {
        build_scale_x_args[["x_min"]] <- x_minimum
      }
      if(is.null(build_scale_x_args[["x_max"]])) {
        build_scale_x_args[["x_max"]] <- x_maximum
      }
     
      
      defaults <- base::as.list(base::formals(build_scale_x))
      build_scale_x_args <- utils::modifyList(defaults, build_scale_x_args)
      add_build_scale_x <- base::do.call(build_scale_x, build_scale_x_args)
      
      if(is.na(uvarby)) { 
        if(is.na(out_u_[['groupby']][1])) {
          legendlabs_mult_singel <- addylab_dv # c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_u_$groupby_line.x <- legendlabs_mult_singel[1]
          out_u_$groupby_color.x <- legendlabs_mult_singel[1]
          out_u_$groupby_line.y <- legendlabs_mult_singel[2]
          out_u_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_u_$groupby_line.x <- out_u_$groupby.x
          out_u_$groupby_color.x <- out_u_$groupby.x
          out_u_$groupby_line.y <- out_u_$groupby.y
          out_u_$groupby_color.y <- out_u_$groupby.y
          legendlabs_mult_mult <- unique(out_u_[['groupby']])
        }
      }
      if(!is.na(uvarby)) { 
        if(is.null(cov_factor_vars)) {
          legendlabs_mult_singel <- addylab_dv # c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_u_$groupby_line.x <- legendlabs_mult_singel[1]
          out_u_$groupby_color.x <- legendlabs_mult_singel[1]
          out_u_$groupby_line.y <- legendlabs_mult_singel[2]
          out_u_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_u_$groupby_line.x <- out_u_$groupby.x
          out_u_$groupby_color.x <- out_u_$groupby.x
          out_u_$groupby_line.y <- out_u_$groupby.y
          out_u_$groupby_color.y <- out_u_$groupby.y
          legendlabs_mult_mult <- unique(out_u_[['groupby']])
        }
      }
      plot.o.u <- out_u_ %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = !!as.name(Yy),
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, color.groupby),
            group = groupby.y
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::labs(x = label.x, y = label.d, title  = "") +
        add_build_scale_x + 
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(y = paste0("Unadjusted ", "Individual Curves")) +
        ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90)) +
        ggplot2::labs(title = label.unadj) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      getbuiltingg <- ggplot2::ggplot_build(plot.o.u)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- 
        rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- 
        rep(get_color_, ngrpanels)
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- 'solid' # legendlabs_mult_line
        get_color_ <- 'black' # legendlabs_mult_color
        legendlabs_ <- NULL # legendlabs_mult_singel
      }
      if (!grepl("a", opt, ignore.case = T)) {
        suppressMessages({
        })
      }
      if (nchar(opt) == 1) {
        plot.o.u <- plot.o.u +
          ggplot2::theme(plot.title = ggplot2::element_blank())
      } else if (nchar(opt) > 1) {
        plot.o.u <- plot.o.u +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
    } else if (!grepl("u", opt, ignore.case = T)) {
      plot.o.u <- NULL
    }
    if (grepl("a", opt, ignore.case = T) &
        !grepl("u", opt, ignore.case = T)) {
      plot.o <- plot.o.a
    } else if (!grepl("a", opt, ignore.case = T) &
               grepl("u", opt, ignore.case = T)) {
      plot.o <- plot.o.u
    } else if (grepl("a", opt, ignore.case = T) &
               grepl("u", opt, ignore.case = T)) {
      if (layout == 'facet') {
        out_a_u_ <-
          d.out <- out_a_ %>% dplyr::mutate(curve = 'Adjusted') %>%
          dplyr::bind_rows(., out_u_ %>%
                             dplyr::mutate(curve = 'Unadjusted')) %>%
          data.frame()
        

        build_scale_x_args <- list()
        for (i in names(formals(build_scale_x))) {
          build_scale_x_args[[i]] <- transform_xaxis[[i]]
        }
        if(is.null(build_scale_x_args[["x_min"]])) {
          build_scale_x_args[["x_min"]] <- x_minimum
        }
        if(is.null(build_scale_x_args[["x_max"]])) {
          build_scale_x_args[["x_max"]] <- x_maximum
        }
        
        defaults <- base::as.list(base::formals(build_scale_x))
        build_scale_x_args <- utils::modifyList(defaults, build_scale_x_args)
        add_build_scale_x <- base::do.call(build_scale_x, build_scale_x_args)
        
        
        out_a_u_ <- out_a_u_[out_a_u_[[Xx]] >= x_minimum & 
                               out_a_u_[[Xx]] <= x_maximum, ]
        out_a_u_ <- out_a_u_ %>% dplyr::mutate(groupby.x = groupby, 
                                               groupby.y = groupby.x)
        if(is.na(uvarby)) {
          if(is.na(out_a_u_[['groupby']][1])) {
            legendlabs_mult_singel <- addylab_dv # c('Distance', 'Velocity')
            legendlabs_mult_color <- single_plot_pair_color_dv_au
            legendlabs_mult_line <- c('solid', 'solid')
            out_a_u_$groupby_line.x <- legendlabs_mult_singel[1]
            out_a_u_$groupby_color.x <- legendlabs_mult_singel[1]
            out_a_u_$groupby_line.y <- legendlabs_mult_singel[2]
            out_a_u_$groupby_color.y <- legendlabs_mult_singel[2]
          } else {
            out_a_u_$groupby_line.x <- out_a_u_$groupby.x
            out_a_u_$groupby_color.x <- out_a_u_$groupby.x
            out_a_u_$groupby_line.y <- out_a_u_$groupby.y
            out_a_u_$groupby_color.y <- out_a_u_$groupby.y
            legendlabs_mult_mult <- unique(out_a_u_[['groupby']])
          }
        }
        if(!is.na(uvarby)) {
          if(is.null(cov_factor_vars)) {
            legendlabs_mult_singel <- addylab_dv # c('Distance', 'Velocity')
            legendlabs_mult_color <- single_plot_pair_color_dv_au
            legendlabs_mult_line <- c('solid', 'solid')
            out_a_u_$groupby_line.x <- legendlabs_mult_singel[1]
            out_a_u_$groupby_color.x <- legendlabs_mult_singel[1]
            out_a_u_$groupby_line.y <- legendlabs_mult_singel[2]
            out_a_u_$groupby_color.y <- legendlabs_mult_singel[2]
          } else {
            out_a_u_$groupby_line.x <- out_a_u_$groupby.x
            out_a_u_$groupby_color.x <- out_a_u_$groupby.x
            out_a_u_$groupby_line.y <- out_a_u_$groupby.y
            out_a_u_$groupby_color.y <- out_a_u_$groupby.y
            legendlabs_mult_mult <- unique(out_a_u_[['groupby']])
          }
        }
        plot.o <- out_a_u_ %>%
          ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
          ggplot2::geom_line(
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby.x,
              linetype = line_color_key(.data, linetype.groupby),
              color = line_color_key(.data, color.groupby)
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::labs(x = label.x,
                        y = label.d,
                        title  = "") +
          add_build_scale_x + 
          jtools::theme_apa(legend.pos = legendpos) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::theme(axis.title.y.right =
                           ggplot2::element_text(angle = 90)) +
          ggplot2::theme(legend.position = "none") +
          ggplot2::labs(y = paste0("Individual curves")) +
          ggplot2::facet_wrap(~ curve, scales = 'free_x')
        
        getbuiltingg <- ggplot2::ggplot_build(plot.o)
        get_line_  <- getbuiltingg$data[[1]]["linetype"]
        get_color_ <- getbuiltingg$data[[1]]["colour"]
        get_fill_  <- getbuiltingg$data[[1]]["colour"]
        ngrpanels  <- getbuiltingg$data[[1]]["group"]
        ngrpanels <- length(unique(unlist(ngrpanels)))
        get_line_ <- unique(unlist(get_line_))
        get_color_ <- unique(unlist(get_color_))
        if(length(get_line_) != ngrpanels) get_line_ <- 
          rep(get_line_, ngrpanels)
        if(length(get_color_) != ngrpanels) get_color_ <- 
          rep(get_color_, ngrpanels)
        if(ngrpanels > 1) {
          get_line_ <- get_line_
          get_color_ <- get_color_
          legendlabs_ <- legendlabs_mult_mult
        } else if(ngrpanels == 1) {
          get_line_ <- 'solid'  # legendlabs_mult_line
          get_color_ <- 'black' # legendlabs_mult_color
          legendlabs_ <- NULL   # legendlabs_mult_singel
        }
        plot.o <- plot.o +
          ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
          ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
        
      }
      if (layout == 'single') {
        plot.o <- out_a_ %>%
          ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
          ggplot2::geom_line(
            data = out_u_,
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby,
              colour = label.unadj
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::geom_line(
            data = out_a_,
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby,
              colour = label.adj
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::labs(x = label.x,
                        y = label.d,
                        title  = "") +
          ggplot2::scale_color_manual(values = single_plot_pair_color_dv_au) +
          add_build_scale_x + 
          jtools::theme_apa(legend.pos = legendpos.adj.unadj) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::labs(y = paste0("Individual curves")) +
          ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90))
      }
    }
  }

  if(opt == 'd' | opt == 'D') return(plot.o.d)
  if(opt == 'v' | opt == 'V') return(plot.o.v)
  if(opt == 'a' | opt == 'a') return(plot.o.a)
  if(opt == 'u' | opt == 'u') return(plot.o.u)
  if(opt.org == 'o' | opt.org == 'O') return(plot.o.O)
} 





rename_plot_list <- function(plot.list, unique_opt_sort, opt) {
  # print(opt)
  # print(unique_opt_sort)
  
  # "" and character(0) - example opt = "DO",
  if(is_emptyx(unique_opt_sort) & is_emptyx(opt)) {
    return(plot.list)
  }
  
  if(unique_opt_sort != opt) {
    if(grepl("D", unique_opt_sort) & grepl("D", opt)) {
      if("d" %in% names(plot.list)) {
        plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
      }
    }
    if(grepl("V", unique_opt_sort) & grepl("V", opt)) {
      if("d" %in% names(plot.list)) {
        plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
      }
    }  
  } else if(unique_opt_sort == opt & opt == "DV") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "dV") {
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Dv") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Da") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Du") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Va") {
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Vu") {
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "dVa") {
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Dva") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
  } else if(unique_opt_sort == opt & opt == "DVa") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "dVu") {
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Dvu") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
  } else if(unique_opt_sort == opt & opt == "DVu") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Dau") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Vau") {
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "Dvau") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
  } else if(unique_opt_sort == opt & opt == "dVau") {
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } else if(unique_opt_sort == opt & opt == "DVau") {
    plot.list[['D']] <- plot.list[['d']]; plot.list[['d']] <- NULL
    plot.list[['V']] <- plot.list[['v']]; plot.list[['v']] <- NULL
  } 
  
  plot.list <- plot.list[strsplit(unique_opt_sort, "")[[1]]]
  return(plot.list)
}




is_patchwork <- function(x) {
  inherits(x, "patchwork")
}



