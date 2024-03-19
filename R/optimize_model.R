

#' @title Optimize SITAR model
#' 
#' @description Select the best fitting SITAR model that involves choosing the
#'   optimum degrees of freedom (\code{df}) for the natural cubic-spline curve
#'   and the appropriate transformations of the predictor \code{x} and response
#'   \code{y} variables.
#'
#' @param optimize_df A list of integers specifying the degree of freedom
#'   (\code{df}) values to be optimized. If \code{NULL} (default), the \code{df}
#'   is taken from the original model. For optimization over different
#'   \code{df}, say for example \code{df} 4 and \code{df} 5, the corresponding
#'   code is \code{optimize_df = list(4,5)}. For \code{univariate_by} and
#'   \code{multivariate} models, \code{optimize_df} can be a single integer
#'   (e.g., \code{optimize_df = 4}) or a list (e.g., \code{optimize_df =
#'   list(4,5)}), or a a list of lists. As an example, consider optimization
#'   over \code{df} 4 and \code{df} 5 for the first sub model, and \code{df} 5
#'   and \code{df} 6 for the second sub model, the corresponding code is
#'   \code{optimize_df = list(list(4,5), list(5,6))}.
#'
#' @param optimize_x A vector specifying the transformations for the predictor
#'   variable (i.e., \code{x}). The options available are \code{NULL},
#'   \code{'log'}, \code{'sqrt'}, or their combinations. Note that user need not
#'   to enclose these options in a single or double quotes as they are take care
#'   of internally. The default setting is to explore all possible combination
#'   i.e., \code{optimize_x = list(NULL, log,  sqrt)}. Similar to the
#'   \code{optimize_df}, user can specify different \code{optimize_x} for
#'   \code{univariate_by} and \code{multivariate} sub models.
#'
#' @param optimize_y A vector specifying the transformations of the the response
#'   variable (i.e., \code{y}). The approach and options available for
#'   \code{optimize_y} are same as described above for the \code{optimize_x}.
#'   
#' @param transform_prior_class A character vector (default \code{NULL})
#'   specifying the transformations of location-scale based priors such as
#'   \code{normal()} when response variable (i.e., \code{y}) is \code{'log'} or
#'   \code{'sqrt'} transformed. The prior type that could be transformed are
#'   \code{'beta'}, \code{'sd'}, \code{'rsd'}, \code{'sigma'} and \code{'dpar'}.
#'   Currently it is available only for \code{'log'} transformed \code{y}. Each
#'   prior type (i.e., \code{'beta', 'sd', 'rsd', 'sigma', 'dpar'}) specified
#'   via \code{transform_prior_class} is log transformed as follows: \cr
#'  \code{log_location = log(location / sqrt(scale^2 / location^2 + 1))}, \cr 
#'  \code{log_scale = sqrt(log(scale^2 / location^2 + 1))}, \cr 
#'  where location and scale are the original parameters supplied by the user
#'  and the log_location and log_scale are the equivalent parameters on the log
#'  scale. For more details, see \code{a_prior_beta} argument in [bsitar()]
#'  function. Note that \code{transform_prior_class} is used as an experiment
#'  and therefore results may not be what user intended. Thus we recommend to
#'  explicitly set the desired prior and not to use
#'  \code{transform_prior_class}.
#'  
#' @param transform_beta_coef A character vector (default \code{NULL})
#'   specifying the transformations of location-scale based priors for specific
#'   regression coefficient(s) when response variable (i.e., \code{y}) is
#'   \code{'log'} or \code{'sqrt'} transformed. The coefficient that could be
#'   transformed are \code{'a'}, \code{'b'}, \code{'c'}, \code{'d'} and
#'   \code{'s'}. The default is \code{transform_beta_coef = c('b',' b', 'd')}
#'   which implies that parameters \code{'a'}, \code{'a'} and \code{'a'} will be
#'   transformed whereas parameter \code{'a'} will be left unchanged because
#'   default prior for parameter \code{'a'} is based on outcome  \code{y} itself
#'   (e.g., \code{a_prior_beta = normal(ymean, ysd)}) which has be transformed.
#'   However, we strongly suggest that user explicitly set the desired prior and
#'   not to rely on \code{transform_beta_coef} because it is included on
#'   experimental basis. See \code{transform_prior_class} for details.
#' 
#' @param transform_sd_coef A character vector (default \code{NULL}) specifying
#'   the transformations of location-scale based priors for specific group level
#'   coefficient(s) when response variable (i.e., \code{y}) is \code{'log'} or
#'   \code{'sqrt'} transformed. The coefficient that could be transformed are
#'   \code{'a'}, \code{'b'}, \code{'c'}, \code{'d'} and \code{'s'}. The default
#'   is \code{transform_beta_coef = c('b',' b', 'd')}. See
#'   \code{transform_prior_class} and \code{transform_beta_coef}  for details.
#'  
#' @param exclude_default_funs A logical to indicate whether transformations for
#'   (\code{x} and \code{y}) variables used in the original model fit should be
#'   excluded. If \code{TRUE} (default), the transformations specified for the
#'   \code{x} and \code{y} variables in the original model fit are excluded from
#'   the \code{optimize_x} and \code{optimize_y}. From example, if original
#'   model is fit with \code{xvar = log} and \code{yvar = NULL}, then
#'   \code{optimize_x} is translated into \code{optimize_x = list(NULL, sqrt)},
#'   and similarly \code{optimize_y} is reset as \code{optimize_y = list(log,
#'   sqrt)}.
#'
#' @param add_fit_criteria An optional argument (default \code{NULL}) to
#'   indicate whether to add fit criteria to the returned model fit. Options
#'   available are \code{'loo'} and \code{'waic'}. Please see
#'   [brms::add_criterion()] for details.
#'
#' @param add_bayes_R An optional argument (default \code{NULL}) to indicate
#'   whether to add Bayesian R square to the returned model fit. To estimate and
#'   add \code{bayes_R2} to the model fit, the argument \code{add_bayes_R} is
#'   set as \code{add_bayes_R = 'bayes_R2'}.
#'
#' @param byresp A logical (default \code{FALSE}) to indicate if response wise
#'   fit criteria to be calculated. This argument is evaluated only for the
#'   \code{multivariate} model in which user can select whether to get joint
#'   calculation of point wise log likelihood (\code{byresp = FALSE}) or
#'   response specific (\code{byresp = TRUE}). For, \code{univariate_by} model,
#'   the only option available is to calculate separate point wise log
#'   likelihood for each sub-model, i.e., \code{byresp = TRUE}.
#'
#' @param cores The number of cores to used in parallel processing (default
#'   \code{1}). The argument \code{cores} is passed to the
#'   [brms::add_criterion()].
#' 
#' @param ... Other arguments passed to \code{\link{update_model}}.
#' 
#' @inheritParams growthparameters.bgmfit
#'
#' @return A list containing the optimized models of class \code{bgmfit}, and
#'   the the summary statistics if \code{add_fit_criteria} and/or
#'   \code{add_bayes_R} are specified.
#'  
#' @export optimize_model.bgmfit
#' @export
#' 
#' @importFrom loo pareto_k_table
#'  
#' @seealso [brms::add_criterion()]
#'
#' @inherit berkeley author
#'
#' @examples
#' 
#' \donttest{
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Below example shows dummy call to optimization to save time. 
#' # Note that in case degree of freedom and both  optimize_x and optimize_y are
#' # NULL (i.e., nothing to optimize), the original model object is returned.   
#' # To explicitly get this information whether model is being optimized or not, 
#' # user can set verbose = TRUE. The verbose = TRUE also useful in getting the
#' # information regarding what all arguments have been changed as compared to
#' # the original model.
#' 
#' model2 <- optimize_model(model, 
#'   optimize_df = NULL, 
#'   optimize_x = NULL, 
#'   optimize_y = NULL,
#'   verbose = TRUE)
#' 
#' }
#' 
optimize_model.bgmfit <- function(model,
                                  newdata = NULL,
                                  optimize_df = NULL,
                                  optimize_x = list(NULL, log,  sqrt),
                                  optimize_y = list(NULL, log,  sqrt),
                                  transform_prior_class = c('beta', 'sd', 
                                                       'rsd', 'sigma', 'dpar'),
                                  transform_beta_coef = c('b', 'c', 'd'),
                                  transform_sd_coef = c('b', 'c', 'd'),
                                  exclude_default_funs = TRUE,
                                  add_fit_criteria = NULL,
                                  add_bayes_R = NULL,
                                  byresp = FALSE,
                                  digits = 2,
                                  cores = 1,
                                  verbose = FALSE,
                                  expose_function = NULL,
                                  usesavedfuns = FALSE,
                                  clearenvfuns = NULL,
                                  envir = NULL,
                                  ...) {
  
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- parent.frame()
  }
  

  check_if_package_installed(model, xcall = NULL)
  
  # Initiate non formalArgs()
  outcome <- NULL;
  xfun <- NULL; 
  yfun <- NULL;
  Parameter <- NULL;
  Estimate <- NULL;
  . <- NULL;
  Criterion <- NULL;
  
  
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
                           deriv = 0,
                           all = FALSE)
  
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
  
  
  need_exposed_function <- FALSE
  if(!is.null(add_fit_criteria)) {
    need_exposed_function <- TRUE
  } else if(is.list(add_fit_criteria)) {
    if(!any(is.null(add_fit_criteria[[1]]))) need_exposed_function <- TRUE
  } else if(!is.null(add_bayes_R)) {
    need_exposed_function <- TRUE
  } else if(is.list(add_bayes_R)) {
    if(!any(is.null(add_bayes_R[[1]]))) need_exposed_function <- TRUE
  }

  
  
  
  # The 'expose_function' must TRUE when adding fit criteria or bayes R2
  if (need_exposed_function) {
    if(is.null(expose_function)) {
      args_o$expose_function <- expose_function <- TRUE
      if(verbose) 
        message("Argument 'expose_function' set to TRUE for fit criteria")
    } else if(!is.null(expose_function)) {
      if (!args_o$expose_function) {
        stop(
          "Argument 'expose_function' must be set to TRUE ",
          "\n ",
          " when adding 'fit criteria' or 'bayes_R'"
        )
      }
    } # if(is.null(expose_function)) {
  } # if (need_exposed_function) {
  
  
  
  if(!is.null(call_o_args$expose_function)) {
    args_o$expose_function <- call_o_args$expose_function
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
  
  if(!is.null(optimize_df)) {
    if(is.list(optimize_df)) {
      optimize_df <- unlist(optimize_df)
    } else {
      optimize_df <- optimize_df
    }
  } else if(is.null(optimize_df)) {
    optimize_df <- model$model_info$dfs
  }
  
  
  optimize_df <- as.factor(optimize_df)
  
  # optimize_df <- get_args_opt(deparse(substitute(optimize_df)))
  
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
  # If NULL, then combined log likelihood used for multivariate model
  # otherwise separate log likelihood  for each response

  add_citeria_fun <- function(fit,
                              add_fit_criteria = NULL,
                              add_bayes_R = NULL,
                              resp = NULL,
                              digits = 2,
                              df,
                              xfun_print,
                              yfun_print,
                              usesavedfuns,
                              clearenvfuns,
                              envir,
                              ...) {
    
    assign(o[[1]], fit$model_info[['exefuns']][[o[[1]]]], envir = envir)
   
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
            dplyr::relocate(dplyr::all_of('Parameter'))
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
              dplyr::relocate(dplyr::all_of('Parameter'))
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
                dplyr::relocate(dplyr::all_of('Parameter'))
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
            dplyr::relocate(dplyr::all_of('Parameter'))
          rownames(fit$criteria[[aci_names]]) <- NULL
        }
      }
    } # if (!is.null(add_bayes_R)) {
    
    
    
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
        c('Good', "Bad", "Very bad")
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
      # enverr. <- parent.frame()
      enverr. <- environment()
      assign('err.', FALSE, envir = enverr.)
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
          assign('err.', TRUE, envir = enverr.)
        }
      )
      err. <- get('err.', envir = enverr.)
      if (err.) {
        summary_waic <- NULL
      } else {
        summary_waic <- summary_waic
      }
      fit$summary_waic <- summary_waic
    }
    
    
    
    
    
    if ('bayes_R2' %in% add_bayes_R) {
      # enverr. <- parent.frame()
      enverr. <- environment()
      assign('err.', FALSE, envir = enverr.)
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
          assign('err.', TRUE, envir = enverr.)
        }
      )
      err. <- get('err.', envir = enverr.)
      if (err.) {
        summary_bayes_R2 <- NULL
      } else {
        summary_bayes_R2 <- summary_bayes_R2
      }
      fit$summary_bayes_R2 <-
        summary_bayes_R2 %>% dplyr::select(-dplyr::all_of('Parameter'))
    }
    
    
    
    if ('loo' %in% add_fit_criteria) {
      if ('loo' %in% add_fit_criteria) {
        # enverr. <- parent.frame()
        enverr. <- environment()
        assign('err.', FALSE, envir = enverr.)
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
            assign('err.', TRUE, envir = enverr.)
          }
        )
        err. <- get('err.', envir = enverr.)
        if (err.) {
          summary_loo <- NULL
        } else {
          summary_loo <- summary_loo
        }
        fit$summary_loo <- summary_loo
      }
      
      if ('loo' %in% add_fit_criteria) {
        # enverr. <- parent.frame()
        enverr. <- environment()
        assign('err.', FALSE, envir = enverr.)
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
            assign('err.', TRUE, envir = enverr.)
          }
        )
        err. <- get('err.', envir = enverr.)
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
  
  
  
  optimize_fun <- function(.x, model, exe_model_fit) {
    message("\nOptimizing model no. ",
            .x,
            " (total ",
            nrow(optimize_df_x_y),
            " models)")
    exe_row <- optimize_df_x_y[.x, ]
    df <- levels(droplevels(exe_row$df))
    xfun <- levels(droplevels(exe_row$xfun))
    yfun <- levels(droplevels(exe_row$yfun))
    
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
    
    
    df_print <- deparse(df)
    
    if(grepl("\\(", df_print)) {
      df_print <- 
        regmatches(df_print, gregexpr("(?<=\\().*?(?=\\))", 
                                      df_print, perl=T))[[1]]
    }
    
    df_print <- eval(parse(text = df_print))
    
    
    if(verbose) {
      cat("\n")
      cat(paste0("df = ", df_print, "; xfun = ", 
                 xfun_print, "; yfun = ", yfun_print),
          "\n")
    }
    
    
    optimization_info <-
      paste0("df = ", df_print, "; xfun = ", 
             xfun_print, "; yfun = ", yfun_print)
    
    
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
        message("Arguemnets supplied for optimize_model()' call are same as ",
                "the original model fit.", 
                "\n ",
                "Therefore, returning the original model fit")
        cat("\n")
      }
      fit <- NULL
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
      # Setting it to FALSE because we are exposing it anyways below
      user_call$expose_function <- FALSE
      ####
      # Modify priors for log transformed outcome y
      transform_allowed_dist <- 
        c('normal', 'student_t', 'student_nu', 'cauchy', 'lognormal')
      
      if(is.null(transform_beta_coef)) {
        transform_beta_a <- transform_beta_b <- 
          transform_beta_c <- transform_beta_d <- transform_beta_s <- FALSE
      } else if(!is.null(transform_beta_coef)) {
        if('a' %in% transform_beta_coef) {
          transform_beta_a <- TRUE 
        } else {
          transform_beta_a <- FALSE
        }
        if('b' %in% transform_beta_coef) {
          transform_beta_b <- TRUE 
        } else {
          transform_beta_b <- FALSE
        }
        if('c' %in% transform_beta_coef) {
          transform_beta_c <- TRUE 
        } else {
          transform_beta_c <- FALSE
        }
        if('d' %in% transform_beta_coef) {
          transform_beta_d <- TRUE 
        } else {
          transform_beta_d <- FALSE
        }
        if('s' %in% transform_beta_coef) {
          transform_beta_s <- TRUE 
        } else {
          transform_beta_s <- FALSE
        }
      }
      
      
      if(is.null(transform_sd_coef)) {
        transform_sd_a <- transform_sd_b <- 
          transform_sd_c <- transform_sd_d <- transform_sd_s <- FALSE
      } else if(!is.null(transform_sd_coef)) {
        if('a' %in% transform_sd_coef) {
          transform_sd_a <- TRUE 
        } else {
          transform_sd_a <- FALSE
        }
        if('b' %in% transform_sd_coef) {
          transform_sd_b <- TRUE 
        } else {
          transform_sd_b <- FALSE
        }
        if('c' %in% transform_sd_coef) {
          transform_sd_c <- TRUE 
        } else {
          transform_sd_c <- FALSE
        }
        if('d' %in% transform_sd_coef) {
          transform_sd_d <- TRUE 
        } else {
          transform_sd_d <- FALSE
        }
        if('s' %in% transform_sd_coef) {
          transform_sd_s <- TRUE 
        } else {
          transform_sd_s <- FALSE
        }
      }
      
      if(is.null(transform_prior_class)) {
        transform_class_beta <- transform_class_sd <- 
          transform_class_rsd <- transform_class_sigma <- 
          transform_class_dpar <- FALSE
      } else if(!is.null(transform_prior_class)) {
        if('beta' %in% transform_prior_class) {
          transform_class_beta <- TRUE 
        } else {
          transform_class_beta <- FALSE
        }
        if('sd' %in% transform_prior_class) {
          transform_class_sd <- TRUE 
        } else {
          transform_class_sd <- FALSE
        }
        if('rsd' %in% transform_prior_class) {
          transform_class_rsd <- TRUE 
        } else {
          transform_class_rsd <- FALSE
        }
        if('sigma' %in% transform_prior_class) {
          transform_class_sigma <- TRUE 
        } else {
          transform_class_sigma <- FALSE
        }
        if('dpar' %in% transform_prior_class) {
          transform_class_dpar <- TRUE 
        } else {
          transform_class_dpar <- FALSE
        }
      }
      
     
      
      if(!is.null(args_o$yfun)) {
        if(args_o$yfun == "log") {
          set_fxls <- 'log'
          for (user_calli in names(user_call)) {
            if(grepl("_prior", user_calli)) {
              if(grepl("_beta", user_calli)) {
                transform_beta_coef_tf <- FALSE
                if(grepl("a_prior", user_calli)) {
                  if(transform_beta_a & transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(transform_beta_a & !transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(!transform_beta_a & !transform_class_beta) 
                    transform_beta_coef_tf <- FALSE
                } 
                if(grepl("b_prior", user_calli)) {
                  if(transform_beta_b & transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(transform_beta_b & !transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(!transform_beta_b & !transform_class_beta) 
                    transform_beta_coef_tf <- FALSE
                } 
                if(grepl("c_prior", user_calli)) {
                  if(transform_beta_c & transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(transform_beta_c & !transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(!transform_beta_c & !transform_class_beta) 
                    transform_beta_coef_tf <- FALSE
                } 
                if(grepl("d_prior", user_calli)) {
                  if(transform_beta_d & transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(transform_beta_d & !transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(!transform_beta_d & !transform_class_beta) 
                    transform_beta_coef_tf <- FALSE
                } 
                if(grepl("s_prior", user_calli)) {
                  if(transform_beta_s & transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(transform_beta_s & !transform_class_beta) 
                    transform_beta_coef_tf <- TRUE
                  if(!transform_beta_s & !transform_class_beta) 
                    transform_beta_coef_tf <- FALSE
                } 
                
                if(!is.null(user_call[[user_calli]])) {
                  if(transform_beta_coef_tf) {
                    tem_dist <- deparse(user_call[[user_calli]][[1]])
                    tem_as_str <- FALSE
                    if(grepl("\"",  tem_dist)) {
                      tem_as_str <- TRUE
                      tem_dist <- gsub("\"", "", tem_dist)
                      tem_dist <- strsplit(tem_dist, "\\(")[[1]][1]
                    }
                    tem_prior <- deparse(user_call[[user_calli]])
                    if(!tem_dist %in% transform_allowed_dist) {
                      stop("Tranformation of '", tem_dist, "; distribution ", 
                           "for ", user_calli, " is not allowed.",
                           "\n", 
                           "  Please adjust 'transform_prior_class' argument",
                           " which at present is specified as:",
                           "\n", 
                           "  ",  paste(transform_prior_class, collapse = ", "), 
                           "\n", 
                           "  Else, change prior '", user_calli, "' from its",
                           " current formulation ", tem_prior, 
                           "\n",
                           " ", " to one that",
                           " is based on one of the following distributions",
                           "\n", 
                           "  ",paste(transform_allowed_dist, collapse = ", ")
                      )
                    }
                    tem <- user_call[[user_calli]]
                    if(tem_as_str) tem <- parse(text = tem)[[1]]
                    tem <- rlang::call_modify(tem, fxls = set_fxls)
                    user_call[[user_calli]] <- tem
                  }
                }
              }
              if(grepl("_sd", user_calli) & !grepl("$sigma", user_calli)) {
                transform_sd_coef_tf <- FALSE
                if(grepl("a_prior", user_calli)) {
                  if(transform_sd_a & transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(transform_sd_a & !transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(!transform_sd_a & !transform_class_beta) 
                    transform_sd_coef_tf <- FALSE
                } 
                if(grepl("b_prior", user_calli)) {
                  if(transform_sd_b & transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(transform_sd_b & !transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(!transform_sd_b & !transform_class_beta) 
                    transform_sd_coef_tf <- FALSE
                } 
                if(grepl("c_prior", user_calli)) {
                  if(transform_sd_c & transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(transform_sd_c & !transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(!transform_sd_c & !transform_class_beta) 
                    transform_sd_coef_tf <- FALSE
                } 
                if(grepl("d_prior", user_calli)) {
                  if(transform_sd_d & transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(transform_sd_d & !transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(!transform_sd_d & !transform_class_beta) 
                    transform_sd_coef_tf <- FALSE
                } 
                if(grepl("s_prior", user_calli)) {
                  if(transform_sd_s & transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(transform_sd_s & !transform_class_beta) 
                    transform_sd_coef_tf <- TRUE
                  if(!transform_sd_s & !transform_class_beta) 
                    transform_sd_coef_tf <- FALSE
                } 
                if(!is.null(user_call[[user_calli]])) {
                  if(transform_sd_coef_tf) {
                    tem_dist <- deparse(user_call[[user_calli]][[1]])
                    tem_as_str <- FALSE
                    if(grepl("\"",  tem_dist)) {
                      tem_as_str <- TRUE
                      tem_dist <- gsub("\"", "", tem_dist)
                      tem_dist <- strsplit(tem_dist, "\\(")[[1]][1]
                    }
                    tem_prior <- deparse(user_call[[user_calli]])
                    if(!tem_dist %in% transform_allowed_dist) {
                      stop("Tranformation of '", tem_dist, "; distribution ", 
                           "for ", user_calli, " is not allowed.",
                           "\n", 
                           "  Please adjust 'transform_prior_class' argument",
                           " which at present is specified as:",
                           "\n", 
                           "  ",  paste(transform_prior_class, collapse = ", "), 
                           "\n", 
                           "  Else, change prior '", user_calli, "' from its",
                           " current formulation ", tem_prior, 
                           "\n",
                           " ", " to one that",
                           " is based on one of the following distributions",
                           "\n", 
                           "  ",paste(transform_allowed_dist, collapse = ", ")
                      )
                    }
                    tem <- user_call[[user_calli]]
                    if(tem_as_str) tem <- parse(text = tem)[[1]]
                    tem <- rlang::call_modify(tem, fxls = set_fxls)
                    user_call[[user_calli]] <- tem
                  }
                }
              }
              if(grepl("rsd_", user_calli) & grepl("sigma", user_calli)) {
                if(!is.null(user_call[[user_calli]])) {
                  if(transform_class_rsd) {
                    tem_dist <- deparse(user_call[[user_calli]][[1]])
                    tem_as_str <- FALSE
                    if(grepl("\"",  tem_dist)) {
                      tem_as_str <- TRUE
                      tem_dist <- gsub("\"", "", tem_dist)
                      tem_dist <- strsplit(tem_dist, "\\(")[[1]][1]
                    }
                    tem_prior <- deparse(user_call[[user_calli]])
                    if(!tem_dist %in% transform_allowed_dist) {
                      stop("Tranformation of '", tem_dist, "; distribution ", 
                           "for ", user_calli, " is not allowed.",
                           "\n", 
                           "  Please adjust 'transform_prior_class' argument",
                           " which at present is specified as:",
                           "\n", 
                           "  ",  paste(transform_prior_class, collapse = ", "), 
                           "\n", 
                           "  Else, change prior '", user_calli, "' from its",
                           " current formulation ", tem_prior, 
                           "\n",
                           " ", " to one that",
                           " is based on one of the following distributions",
                           "\n", 
                           "  ",paste(transform_allowed_dist, collapse = ", ")
                      )
                    }
                    tem <- user_call[[user_calli]]
                    if(tem_as_str) tem <- parse(text = tem)[[1]]
                    tem <- rlang::call_modify(tem, fxls = set_fxls)
                    user_call[[user_calli]] <- tem
                  }
                }
              }
              if(grepl("dpar_", user_calli) & grepl("sigma", user_calli)) {
                if(!is.null(user_call[[user_calli]])) {
                  if(transform_class_dpar) {
                    tem_dist <- deparse(user_call[[user_calli]][[1]])
                    tem_as_str <- FALSE
                    if(grepl("\"",  tem_dist)) {
                      tem_as_str <- TRUE
                      tem_dist <- gsub("\"", "", tem_dist)
                      tem_dist <- strsplit(tem_dist, "\\(")[[1]][1]
                    }
                    tem_prior <- deparse(user_call[[user_calli]])
                    if(!tem_dist %in% transform_allowed_dist) {
                      stop("Tranformation of '", tem_dist, "; distribution ", 
                           "for ", user_calli, " is not allowed.",
                           "\n", 
                           "  Please adjust 'transform_prior_class' argument",
                           " which at present is specified as:",
                           "\n", 
                           "  ",  paste(transform_prior_class, collapse = ", "), 
                           "\n", 
                           "  Else, change prior '", user_calli, "' from its",
                           " current formulation ", tem_prior, 
                           "\n",
                           " ", " to one that",
                           " is based on one of the following distributions",
                           "\n", 
                           "  ",paste(transform_allowed_dist, collapse = ", ")
                      )
                    }
                    tem <- user_call[[user_calli]]
                    if(tem_as_str) tem <- parse(text = tem)[[1]]
                    tem <- rlang::call_modify(tem, fxls = set_fxls)
                    user_call[[user_calli]] <- tem
                  }
                }
              }
              if(grepl("sigma", user_calli) & grepl("_sd", user_calli)) {
                if(!is.null(user_call[[user_calli]])) {
                  if(transform_class_sigma) {
                    tem_dist <- deparse(user_call[[user_calli]][[1]])
                    tem_as_str <- FALSE
                    if(grepl("\"",  tem_dist)) {
                      tem_as_str <- TRUE
                      tem_dist <- gsub("\"", "", tem_dist)
                      tem_dist <- strsplit(tem_dist, "\\(")[[1]][1]
                    }
                    tem_prior <- deparse(user_call[[user_calli]])
                    if(!tem_dist %in% transform_allowed_dist) {
                      stop("Tranformation of '", tem_dist, "; distribution ", 
                           "for ", user_calli, " is not allowed.",
                           "\n", 
                           "  Please adjust 'transform_prior_class' argument",
                           " which at present is specified as:",
                           "\n", 
                           "  ",  paste(transform_prior_class, collapse = ", "), 
                           "\n", 
                           "  Else, change prior '", user_calli, "' from its",
                           " current formulation ", tem_prior, 
                           "\n",
                           " ", " to one that",
                           " is based on one of the following distributions",
                           "\n", 
                           "  ",paste(transform_allowed_dist, collapse = ", ")
                      )
                    }
                    tem <- user_call[[user_calli]]
                    if(tem_as_str) tem <- parse(text = tem)[[1]]
                    tem <- rlang::call_modify(tem, fxls = set_fxls)
                    user_call[[user_calli]] <- tem
                  }
                }
              }
              
            } # if(grepl("_prior", user_calli)) {
          } # for (user_calli in names(user_call)) {
        } # if(args_o$yfun == "log") {
      } # if(!is.null(args_o$yfun == "log")) {
      
      
      
  
      ###
      fit <- eval(user_call)
  
      if(!exe_model_fit) {
        return(fit)
      }
      
      # if(!exe_model_fit) {
      #   if(get_priors) {
      #     return(do.call(brms::get_prior, brm_args))
      #   } else if(get_standata) {
      #     return(do.call(brms::make_standata, brm_args))
      #   } else if(get_stancode) {
      #     return(scode_final)
      #   } else if(get_priors_eval) {
      #     return(get_priors_eval_out)
      #   } else if(validate_priors) {
      #     return(do.call(brms::validate_prior, brm_args))
      #   } else if(get_init_eval) {
      #     return(brm_args$init)
      #   } else if(get_formula) {
      #     return(brm_args$formula)
      #   } else if(get_stanvars) {
      #     return(brm_args$stanvars)
      #   }
      # } 
      
      if("brmsfit" %in% class(fit) | "bgmfit" %in% class(fit)) {
        class_fit <- TRUE
      } else {
        class_fit <- FALSE
      }
      
      
      
      if(args_o$expose_function) {
        if(is.null(envir)) {
          if(!is.null(fit$model_info$exefuns[[1]])) {
            envir <- environment(fit$model_info$exefuns[[1]])
          } else {
            envir <- parent.frame()
          }
        }
        
        if(is.null(usesavedfuns)) {
          if(!is.null(fit$model_info$exefuns[[1]])) {
            usesavedfuns <- TRUE
          } else if(is.null(fit$model_info$exefuns[[1]])) {
            if(expose_function) {
              fit <- 
                expose_model_functions(fit, envir = envir, verbose = verbose)
              usesavedfuns <- TRUE
            } else if(!expose_function) {
              usesavedfuns <- FALSE
            }
          }
        } else { # if(!is.null(usesavedfuns)) {
          if(!usesavedfuns) {
            if(expose_function) {
              fit <- 
                expose_model_functions(fit, envir = envir, verbose = verbose)
              usesavedfuns <- TRUE
            }
          } else if(usesavedfuns) {
            check_if_functions_exists(fit, checks = TRUE, 
                                      usesavedfuns = usesavedfuns)
          }
        } # if(is.null(envir)) {...
      } # if(args_o$expose_function) {
    } # else if(!all_same_args) {
    
    
    if(!is.null(fit)) {
      fit$model_info$optimization_info <- optimization_info
      fit$model_info$optimize_df <- df_print
      fit$model_info$optimize_x <- xfun_print
      fit$model_info$optimize_y <- yfun_print
      
      # Add fit_criteria and bares_R to the fit
      # Also, add summary data frames for criteria and R square
      
      # 'setresp' to anything so that even multivariate will be response wise
      # if desired, this behavior
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
        enverr. <- environment()
        assign('err.', FALSE, envir = enverr.)
        tryCatch(
          expr = {
            fit_ac <- add_citeria_fun(
              fit,
              add_fit_criteria = add_fit_criteria,
              add_bayes_R =  NULL,
              resp = setresp,
              digits = digits,
              df = df,
              xfun_print = xfun_print,
              yfun_print = yfun_print,
              usesavedfuns = usesavedfuns,
              clearenvfuns = clearenvfuns,
              envir = envir
            )
          },
          error = function(e) {
            assign('err.', TRUE, envir = enverr.)
          }
        )
        err. <- get('err.', envir = enverr.)
        if (err.) {
          fit <- fit
        } else {
          fit <- fit_ac
        } # tryCatch
        # fit <- add_citeria_fun(
        #   fit,
        #   add_fit_criteria = add_fit_criteria,
        #   add_bayes_R =  NULL,
        #   resp = setresp,
        #   digits = digits,
        #   df = df,
        #   xfun_print = xfun_print,
        #   yfun_print = yfun_print,
        #   usesavedfuns = usesavedfuns,
        #   clearenvfuns = clearenvfuns,
        #   envir = envir
        # )
      } # if (!is.null(add_fit_criteria)) {
      
      if (!is.null(add_bayes_R)) {
        enverr. <- environment()
        assign('err.', FALSE, envir = enverr.)
        tryCatch(
          expr = {
            fit_rs <- add_citeria_fun(
              fit,
              add_fit_criteria = NULL,
              add_bayes_R =  add_bayes_R,
              resp = setresp,
              digits = digits,
              df = df,
              xfun_print = xfun_print,
              yfun_print = yfun_print,
              usesavedfuns = usesavedfuns,
              clearenvfuns = clearenvfuns,
              envir = envir
            )
          },
          error = function(e) {
            assign('err.', TRUE, envir = enverr.)
          }
        )
        err. <- get('err.', envir = enverr.)
        if (err.) {
          fit <- fit
        } else {
          fit <- fit_rs
        } # tryCatch
        # fit <- add_citeria_fun(
        #   fit,
        #   add_fit_criteria = NULL,
        #   add_bayes_R =  add_bayes_R,
        #   resp = setresp,
        #   digits = digits,
        #   df = df,
        #   xfun_print = xfun_print,
        #   yfun_print = yfun_print,
        #   usesavedfuns = usesavedfuns,
        #   clearenvfuns = clearenvfuns,
        #   envir = envir
        # )
      } # if (!is.null(add_bayes_R)) {
    } # if(!is.null(fit)) {
    
    return(fit)
  }
  
  
  exe_model_fit <- TRUE
  if(args_o$get_stancode |
     args_o$get_standata |
     args_o$get_formula |
     args_o$get_stanvars |
     args_o$get_priors |
     args_o$get_priors_eval |
     args_o$validate_priors |
     args_o$get_init_eval) {
    exe_model_fit <- FALSE
  }
  
  optimize_list <- lapply(1:nrow(optimize_df_x_y), function(.x)
    optimize_fun(.x, model, exe_model_fit))
  
  
  if(!exe_model_fit) {
    return(optimize_list)
  }
  
  
  if(!is.null(optimize_list[[1]])) {
    loo_fit             <- combine_summaries(optimize_list, 'summary_loo')
    loo_diagnostic_fit  <-
      combine_summaries(optimize_list, 'diagnostic_loo')
    waic_fit            <-
      combine_summaries(optimize_list, 'summary_waic')
    bayes_R2_fit        <-
      combine_summaries(optimize_list, 'summary_bayes_R2')
    
    # Parameter column is not created earlier for the 'bayes_R2_fit'
    if(!is.null(bayes_R2_fit)) {
      bayes_R2_fit <- bayes_R2_fit %>% 
        dplyr::mutate(Parameter = 'bayes_R2') %>% 
        dplyr::relocate(df, xfun, yfun, Parameter)
    }
    
    
    ########################
    
    
    
    
    
    ########################
    
    
    
    # loo_fitx <<- loo_fit
    # loo_diagnostic_fitx <<- loo_diagnostic_fit
    # waic_fitx <<- waic_fit
    # bayes_R2_fitx <<- bayes_R2_fit
    
    
    attributes(optimize_list) <- NULL
    
    optimize_summary <- data.frame()
    
    if(exists('loo_fit')) {
      if(!is.null(loo_fit)) {
        loo_fit <- loo_fit %>% 
          dplyr::mutate(Criterion = 'loo') %>% 
          dplyr::relocate(Criterion, .before = dplyr::all_of('Parameter'))
      }
      
      optimize_summary <- optimize_summary %>% 
        dplyr::bind_rows(., loo_fit)
    }
    
    if(exists('loo_diagnostic_fit')) {
      if(!is.null(loo_diagnostic_fit)) {
        loo_diagnostic_fit <- loo_diagnostic_fit %>% 
          dplyr::rename(Parameter = Inference) %>% 
          dplyr::mutate(Criterion = 'loo') %>% 
          dplyr::relocate(Criterion, .before = dplyr::all_of('Parameter'))
        
        loo_diagnostic_fit <- loo_diagnostic_fit %>% 
          dplyr::mutate(Parameter = paste0("Pareto k: ", Parameter)) %>% 
          dplyr::mutate(Estimate = NA, SE = NA) %>% 
          dplyr::relocate(Parameter, .before = dplyr::all_of('Range')) %>% 
          dplyr::relocate(Estimate, .before = dplyr::all_of('Range'))
      }
      
        
      optimize_summary <- optimize_summary %>% 
        dplyr::bind_rows(., loo_diagnostic_fit)
    }
    
    if(exists('waic_fit')) {
      if(!is.null(waic_fit)) {
        waic_fit <- waic_fit %>% 
          dplyr::mutate(Criterion = 'waic') %>% 
          dplyr::relocate(Criterion, .before = dplyr::all_of('Parameter'))
      }
     
      optimize_summary <- optimize_summary %>% 
        dplyr::bind_rows(., waic_fit)
    }
    
    if(exists('bayes_R2_fit')) {
      # Somehow bayes_R2_fit throw error: mutate applied to 'NULL'
      if(!is.null(bayes_R2_fit)) {
      bayes_R2_fit <- bayes_R2_fit %>% 
        dplyr::mutate(Criterion = 'Bayes_R2') %>% 
        dplyr::relocate(Criterion, .before = dplyr::all_of('Parameter'))
      }
   
      optimize_summary <- optimize_summary %>% 
        dplyr::bind_rows(., bayes_R2_fit)
    }
    
    
    # optimize_summary <- optimize_summary %>%
    #   dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
    #                               ~ round(., digits = digits)))
    
    
    if(nrow(optimize_summary) > 0) {
      comment_optimize_summary <- 
        paste("Notes:",
              "\n",
              "Columns 'Range' 'Count' 'Percent' and 'Min.n_eff'",
              " are relevant only for the 'loo Inference'.",
              "\n",
              "Columns CI limits are applicable only for the 'Bayes_R2'",
              "\n",
              "Columns Estimate and SE are applicable for",
              "'loo Estimate', 'waic' and 'Bayes_R2'"
        )
      attr(optimize_summary, 'comment') <- comment_optimize_summary
    } else {
      optimize_summary <- NULL
    }
    
    out <- list(models = optimize_list, optimize_summary = optimize_summary)
    
    return(out)
  } # if(!is.null(optimize_list[[1]])) {
 
  
}



#' @rdname optimize_model.bgmfit
#' @export
optimize_model <- function(model, ...) {
  UseMethod("optimize_model")
}
