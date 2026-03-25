

#' Title An internal function to get adjusted curves (for random effects)
#' 
#' @details
#' Adapted from https://github.com/statist7/sitar/blob/master/R/xyadj.R
#'  v.adj also calculated but not returned 
#'
#' @param model 
#' @param x 
#' @param y 
#' @param id 
#' @param v 
#' @param resp 
#' @param ndraws 
#' @param newdata 
#' @param levels_id 
#' @param abc 
#' @param summary 
#' @param conf 
#' @param robust 
#' @param tomean 
#' @param ipts 
#' @param xrange 
#' @param aux_variables 
#' @param numeric_cov_at 
#' @param idata_method 
#' @param verbose 
#' @param model_deriv 
#' @param deriv 
#' @param envir 
#' @param ... 
#'
#' @returns A data frame
#'
#' @keywords internal
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
      envir <- envir # parent.frame()
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
    
    # o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
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
    
    
    # 6.03.2025
    # Don't evalaute data, error or no new data 
    if (is.null(newdata)) {
      stop("Please specify newdata")
    }
    
    
    
    # This dummy data only to get list_c and its componenrts 
     
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
    
    
    if(!is.null(cov_vars)) {
      # stop("option = 'a' is not yet available for model with covariates")
    }
    
    
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', 'Est.Error', probtitles)
    
    
    
    ######################################################
    # prepare_data2 change 
    # xoffsetXnames <- 'xoffset'
    # randomRnames   <- 'random' 
    # if (!is.null(resp)) xoffsetXnames <- paste0(xoffsetXnames, resp_rev_)
    # xoffsetXnames <- model$model_info[[xoffsetXnames]]
    # xoffset <- xoffsetXnames
    ######################################################
    
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
    

    
    # re_effx <- brms::ranef(berkeley_fit, summary = F)
    # re_effx <- re_effx[['id']]
    # re_effx <- re_effx[match('id', rownames(re_effx)), , drop = FALSE]
    
    
    if(!is.null(ipts)) {
      add_outcome <- model$data %>%
        dplyr::select(dplyr::all_of(c(Yy, idvar)))
      newdata <- newdata %>% 
        dplyr::left_join(., add_outcome, by = c(idvar))
      x <- newdata[[Xx]]
      y <- newdata[[Yy]]
      id <- newdata[[idvar]][1]
    }
    
    ######################################################
    # prepare_data2 change 
    # This x - xoffset is needed for bsitar based computation
    # x <- x - xoffset
    ######################################################
    
    nrowdatadims <- nrow(newdata)
    ### 
    predprep <- brms::prepare_predictions(model, resp = resp, 
                                          newdata = newdata)
    
    rparnames <- names(predprep$nlpars)
    
    respstr <- "" # paste0(resp, "_")
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
      null_a <- matrix(0, nrowdatadims, 1)
      naaa_a <- matrix(0, nrowdatadims, 1)
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
      null_b <- matrix(0, nrowdatadims, 1)
      naaa_b <- matrix(0, nrowdatadims, 1)
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
      null_c <- matrix(0, nrowdatadims, 1)
      naaa_c <- matrix(0, nrowdatadims, 1)
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
      null_d <- matrix(0, nrowdatadims, 1)
      naaa_d <- matrix(0, nrowdatadims, 1)
    }
    
    
    # 6.03.2025
    if(!is.null(draw_ids)) {
      setloopdins <- length(eval(draw_ids))
    } else if(!is.null(ndraws)) {
      setloopdins <- ndraws
    } else {
      # nothing
    }
    
    
    # 6.03.2025
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
        
        # Re create random effects - coef = fixed + random 
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
      } # for (i in 1:ndraws) {
      
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
      
      # 6.03.2025
      ##############################################################
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
      ##############################################################
      
      
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
      
      # This was good but not required
      # out <- cbind(x.adj[, 1], y.adj)
      # setadnamex <- paste0("adj", "_", Xx)
      # setadnamey <- colnames(y.adj)
      # colnames(out) <- c(setadnamex, setadnamey)
      # out <- cbind(newdata, out)
      
      # But for trimline, we need the following order 
      out <- newdata
      out[[Xx]] <- x.adj[, 1]
      out[[Yy]] <- y.adj[, 1]
      out <- out %>% dplyr::relocate(c(Xx, Yy, idvar))
      # now add CI also - Estimate will be same as outcome
      out <- cbind(out, y.adj)
    } # if(!summary) {
    
    
   
    if(summary) {
      r_a <- null_a[ , 1]
      r_b <- null_b[ , 1]
      r_c <- null_c[ , 1]
      r_d <- null_d[ , 1]
      na_a <- naaa_a[ , 1]
      na_b <- naaa_b[ , 1]
      na_c <- naaa_c[ , 1]
      na_d <- naaa_d[ , 1]
      
      # Re create random effects - coef = fixed + random 
      # Thus, random = coef - fixed
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
      
      # This was good
      out <- cbind(x.adj, y.adj)
      setadnamex <- paste0("adj", "_", Xx)
      setadnamey <- 'Estimate'
      colnames(out) <- c(setadnamex, setadnamey)
      out <- cbind(newdata, out)
      
      # But for trimline, we need folowing order 
      out <- newdata
      out[[Xx]] <- x.adj
      out[[Yy]] <- y.adj
      out <- out %>% dplyr::relocate(dplyr::all_of(c(Xx, Yy, idvar)))
    } # if(summary) {
    out
  } 


#' @noRd
#' @exportS3Method xyadj_curves bgmfit
xyadj_curves <- function(model, ...) {
  UseMethod("xyadj_curves")
}

#########################################################################
#########################################################################


#' Title An internal function to get unadjusted curves
#'
#' @param model 
#' @param x 
#' @param y 
#' @param id 
#' @param resp 
#' @param newdata 
#' @param verbose 
#' @param model_deriv 
#' @param deriv 
#' @param envir 
#' @param ndraws 
#' @param draw_ids 
#' @param ...
#'
#' @returns A data frame
#' 
#' @keywords internal
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
  
  
  if(!is.null(cov_vars)) {
    # stop("option = 'a' is not yet available for model with covariates")
  }
  
  
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

#########################################################################
#########################################################################



#' Title An internal function to trim growth curves
#'
#' @param model 
#' @param x 
#' @param y 
#' @param id 
#' @param newdata 
#' @param resp 
#' @param ndraws 
#' @param level 
#' @param trim 
#' @param envir 
#' @param verbose 
#' @param model_deriv 
#' @param deriv 
#' @param draw_ids 
#' @param estimation_method 
#' @param ... 
#' 
#' @returns A data frame
#' 
#' @keywords internal
#' @noRd
#'
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
    # post_processing_checks_args[['all']]      <- FALSE
    # post_processing_checks_args[['verbose']]  <- verbose
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
    
    
    # if (missing(x)) {
    #   .x <- Xx
    # } else {
    #   .x <- x
    # }
    # 
    # if (missing(y)) {
    #   .y <- Yy
    # } else {
    #   .y <- y
    # }
    # 
    # if (missing(id)) {
    #   .id <- idvar
    # } else {
    #   .id <- id
    # }
    
   
    
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
    # 6.03.2025
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


#########################################################################
#########################################################################

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
  
  # https://data.library.virginia.edu/setting-up-color-palettes-in-r/
  ggplotColors <- function(g){
    g <- g - 1
    d <- 360/g
    h <- cumsum(c(15, rep(d,g - 1)))
    O <- grDevices::hcl(h = h, c = 100, l = 65)
    O <- c('black', O)
    O
  }
  
  # https://groups.google.com/g/ggplot2/c/XIcXU3KlxW0
  ggplotlines <- function(g){
    lineTypes1 <- c("solid", "22", "42", "44", "13", "1343", "73", "2262")
    # lineTypes1 <- c("solid", "solid", "solid", "13", "1343", "73", "2262")
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
  } # if(is.na(linetype.groupby) & is.na(color.groupby)) {
  
  
  
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
  } # if(is.na(linetype.groupby) & !is.na(color.groupby)) {
  
  
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
      line.guide <- "none" # "legend"
      color.guide <- "legend"  # "none"
    }
  } # if(!is.na(linetype.groupby) & is.na(color.groupby)) {
  
  
  
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
  } # if(!is.na(linetype.groupby) & !is.na(set.color.groupby)) {
  
  
  suppressMessages({
    plot <- plot + 
      ggplot2::scale_linetype_manual(values=set.line.groupby, 
                                     guide = line.guide) +
      ggplot2::scale_color_manual(values=set.color.groupby, 
                                  guide = color.guide)
  })
  
  plot
} # set_lines_colors


#########################################################################
#########################################################################


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


#########################################################################
#########################################################################

#' An internal function to transform y axis when plotting with dual y axis
#'
#' @param primary Primary y axis.
#' @param secondary secondary y axis.
#' @keywords internal
#' @return A plot object.
#' @noRd
#'
transform.sec.axis <- function(primary,
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



#########################################################################
#########################################################################

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
                                   "
               )
      )
    }
    return(pwobj)
  }

