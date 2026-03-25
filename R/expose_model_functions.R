

#' @title Expose user-defined Stan functions for the Bayesian SITAR model
#'
#' @description The \strong{expose_model_functions()} function is a wrapper
#'   around [rstan::expose_stan_functions()] that exposes user-defined Stan
#'   function(s). These functions are necessary for post-processing the
#'   posterior draws.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @param scode A character string containing the user-defined Stan function(s)
#'   in \code{Stan} code. If \code{NULL} (the default), the \code{scode} will be
#'   retrieved from the \code{model}.
#'
#' @param expose A logical (default \code{TRUE}) to indicate whether to expose
#'   the functions and add them as an attribute to the \code{model}.
#'
#' @param select_model A character string (default \code{NULL}) to specify the
#'   model name. This parameter is for internal use only.
#'
#' @param returnobj A logical (default \code{TRUE}) to specify whether to return
#'   the model object. If \code{expose = TRUE}, it is advisable to set
#'   \code{returnobj = TRUE}.
#'
#' @param vectorize A logical (default \code{FALSE}) to indicate whether the
#'   exposed functions should be vectorized using [base::Vectorize()]. Note that
#'   currently, \code{vectorize} should be set to \code{FALSE}, as setting it to
#'   \code{TRUE} may not work as expected.
#'   
#' @param sigmafun A logical (default \code{FALSE}) to indicate whether to return
#'   the sigma functions. This parameter is for internal use only. Ignored
#'   
#' @param backend A character string (default \code{NULL}) to set up the
#'   the \code{backend} method (\code{'rstan'} or \code{'cmdstanr'}) for
#'   compiling and exposing \code{stan} functions. If \code{NULL}, then the
#'   \code{backend} is same as the \code{backend} method used for model
#'   fitting. Note that when \code{backend = FALSE}, then default
#'   \code{backend}  will be set as \code{'rstan'}. This is particularly useful
#'   because \code{backend = 'cmdstanr'} will fail on \code{WSL} system.
#'   
#' @param path A character string to set up the path of installed
#'   \code{CmdStan}. If \code{NULL} (default)
#'
#' @inherit growthparameters.bgmfit params
#'
#' @param ... Additional arguments passed to the
#'   [rstan::expose_stan_functions()] function. The "..." can be used to set the
#'   compiler, which can be either [rstan::stanc()] or [rstan::stan_model()].
#'   You can also pass other compiler-specific arguments such as \code{save_dso}
#'   for [rstan::stan_model()]. Note that while both [rstan::stanc()] and
#'   [rstan::stan_model()] can be used as compilers before calling
#'   [rstan::expose_stan_functions()], it is important to note that the
#'   execution time for [rstan::stan_model()] is approximately twice as long as
#'   [rstan::stanc()].
#'
#' @return An object of class \code{bgmfit} if \code{returnobj = TRUE};
#'   otherwise, it returns \code{NULL} invisibly.
#'
#' @rdname expose_model_functions
#' @export
#'
#' @seealso [rstan::expose_stan_functions()]
#'
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether the model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # To save time, argument expose is set as FALSE, which runs a dummy test 
#' # and avoids model compilation that often takes time.
#' 
#' expose_model_functions(model, expose = FALSE)
#' }
#' 
expose_model_functions.bgmfit <- function(model, 
                                 scode = NULL, 
                                 expose = TRUE, 
                                 select_model = NULL, 
                                 returnobj = TRUE,
                                 vectorize = FALSE,
                                 verbose = FALSE,
                                 sigmafun = FALSE,
                                 backend = NULL,
                                 path = NULL,
                                 envir = NULL,
                                 ...) {
  
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- envir
  }
  
  fun_env <- envir
  
  if(is.null(select_model)) {
    select_model <- model$model_info$select_model
  }
  
  if(expose) {
    if (verbose) {
      setmsgtxt <-
        paste0("\n Exposing Stan functions for post-processing\n")
      message(setmsgtxt)
    }
  }
  
  sigma_model_all_      <- paste0('sigmamodel_all', "")
  sigma_model_all       <- model$model_info[[sigma_model_all_]]
  
  expose_sigma_ls_model_fun <- FALSE
  expose_sigma_var_model_fun <- FALSE
  expose_sigma_basic_model_fun <- FALSE
  if(!is.null(sigma_model_all)) {
    if(all(grepl("ls", sigma_model_all))) {
      expose_sigma_ls_model_fun <- TRUE
    } else if(all(grepl("basic", sigma_model_all))) {
      expose_sigma_basic_model_fun <- TRUE
    } else if(any(grepl("ls", sigma_model_all)) |
              any(grepl("basic", sigma_model_all))) {
      expose_sigma_ls_model_fun <- TRUE
      expose_sigma_var_model_fun <- TRUE
    } else {
      expose_sigma_ls_model_fun <- FALSE
      expose_sigma_var_model_fun <- FALSE
    }
  } else if(is.null(sigma_model_all)) {
    expose_sigma_ls_model_fun <- FALSE
    expose_sigma_var_model_fun <- FALSE
    expose_sigma_basic_model_fun <- FALSE
  }
  
  
  if(!expose) {
    if (is.null(model$model_info$decomp))  expose_r_from_stan <- TRUE
    if (!is.null(model$model_info$decomp)) expose_r_from_stan <- TRUE # FALSE
  } else {
    expose_r_from_stan <- FALSE
  }
  

  if(expose) {
    if (is.null(scode)) {
      if(model$model_info[['fit_edited_scode']]) {
        if(!is.null(model$model_info$fun_scode)) {
          exposecode <- model$model_info$fun_scode 
        } else if(is.null(model$model_info$fun_scode)) {
          exposecode <- model$model_info$emodel
        }
      } else {
        if(!is.null(model$model_info$fun_scode)) {
          exposecode <- model$model_info$fun_scode 
        } else if(is.null(model$model_info$fun_scode)) {
          exposecode <- brms::stancode(model)
        }
      }
    } else if (!is.null(scode)) {
      exposecode <- scode
    }
  }
  

  if(is.null(backend)) {
    backend <- model$backend
  } else if(!is.null(backend)) {
    if(isFALSE(backend)) {
      backend <- 'rstan'
    } else {
      backend <- backend
    }
  }
  
  
  if(!expose) {
    backend <- "rstan"
  }
  
  
  if(backend == "cmdstanr") {
    if(isTRUE(check_if_cmdstanr_available())) {
      write_stan_file <- 
        utils::getFromNamespace("write_stan_file", "cmdstanr")
      cmdstan_model   <- 
        utils::getFromNamespace("cmdstan_model", "cmdstanr")
      expose_model   <- 
        utils::getFromNamespace("expose_functions", "cmdstanr")
      get_path   <- 
        utils::getFromNamespace("cmdstan_path", "cmdstanr")
      get_cmdstan_default_path   <- 
        utils::getFromNamespace("cmdstan_default_path", "cmdstanr")
      get_cmdstan_default_install_path   <- 
        utils::getFromNamespace("cmdstan_default_install_path", "cmdstanr")
      set_path   <- 
        utils::getFromNamespace("set_cmdstan_path", "cmdstanr")
      
      restore_path <- get_path()
      if(is.null(path)) {
        if(grepl("wsl", restore_path, fixed = TRUE)) {
          if(!grepl("wsl", get_cmdstan_default_path(), fixed = TRUE)) {
            path <- get_cmdstan_default_path()
          } else {
            stop("cmdstanr 'expose_functions' does not work for 'WSL'")
          }
        } else if(!grepl("wsl", restore_path, fixed = TRUE)) {
          path <- restore_path
        }
      } else if(!is.null(path)) {
        path <- path
      }
      suppressWarnings(suppressMessages({set_path(path)}))
    } # if(isTRUE(check_if_cmdstanr_available())) {
  } # if(backend == "cmdstanr") {
  
  
  if(expose &  backend == "cmdstanr") {
    if(is.null(model$functions$fun_names)) {
      suppressWarnings(suppressMessages({
        c_model <- cmdstan_model(write_stan_file(exposecode),
                                 quiet = TRUE,
                                 cpp_options = attr(model$fit, 'cpp_options'),
                                 # stanc_options = attr(model$fit, 'stanc_options'),
                                 dir = NULL,
                                 pedantic = FALSE,
                                 include_paths = attr(model$fit, 'include_paths'),
                                 user_header = NULL,
                                 compile_model_methods = FALSE,
                                 compile_hessian_method = FALSE,
                                 force_recompile = TRUE,
                                 compile_standalone = TRUE)
        # c_model <- attr(model$fit, 'CmdStanModel')
        c_model$expose_model
        set_path(restore_path)
      }))
    } # if(is.null(model$functions$fun_names)) {
    if(!is.null(model$functions$fun_names)) {
      c_model <- attr(model$fit, 'CmdStanModel')
    }
  } # if(expose &  backend == "cmdstanr") {
  

  if(expose & backend == "rstan") {
    stanc_arguments <- list()
    stanc_arguments[['model_code']] <- exposecode
    stan_model_arguments <- stanc_arguments
    dots_args <- list(...)
    if(!is.null(dots_args)) {
      stan_model_arguments <- c(stan_model_arguments, dots_args)
    }
  
    if(is.null(dots_args[['setcompiler']])) {
      compiled_code_via <- 'stanc'
    } else if(!is.null(dots_args[['setcompiler']])) {
      if(dots_args[['setcompiler']] == 'stanc') {
        compiled_code_via <- 'stanc'
      } else if(dots_args[['setcompiler']] == 'stan_model') {
        compiled_code_via <- 'stan_model'
      } else {
        stop(paste("'setcompiler' argument should be 'stanc' or 'stan_model'", 
             "or else NULL", collapse =","))
      }
    }
    
    stanc_arguments[['setcompiler']] <- NULL
    stan_model_arguments[['setcompiler']] <- NULL
    
    if(compiled_code_via == 'stanc') {
      # compiled_code <- do.call(rstan::stanc, stanc_arguments)
      compiled_code <- CustomDoCall(rstan::stanc, stanc_arguments)
    }
    if(compiled_code_via == 'stan_model') {
      # compiled_code <- do.call(rstan::stan_model, stan_model_arguments)
      compiled_code <- CustomDoCall(rstan::stan_model, stan_model_arguments)
    }
    
    compiled_code_args <- list()
    compiled_code_args[['env']] <- fun_env
    compiled_code_args[['stanmodel']] <- compiled_code
    compiled_code_args[['includes']] <- NULL
    compiled_code_args[['show_compiler_warnings']] <- FALSE
    # do.call(rstan::expose_stan_functions, compiled_code_args)
    CustomDoCall(rstan::expose_stan_functions, compiled_code_args)
  }
  
  
  if(expose_r_from_stan) {
    for (funi in 1:length(model$model_info$funlist_r)) {
      assign(gsub("<-.*$", "", model$model_info$funlist_r[funi]),
             ept(model$model_info$funlist_r[funi]), envir = envir)
    }
    
    if(expose_sigma_ls_model_fun) {
      for (funi in 1:length(model$model_info$sigmafunlist_r)) {
        assign(gsub("<-.*$", "", model$model_info$sigmafunlist_r[funi]),
               ept(model$model_info$sigmafunlist_r[funi]), envir = envir)
      }
    } # expose_sigma_ls_model_fun
    
    if(expose_sigma_var_model_fun) {
      for (funi in 1:length(model$model_info$sigmavarfunlist_r)) {
        assign(gsub("<-.*$", "", model$model_info$sigmavarfunlist_r[funi]),
               ept(model$model_info$sigmavarfunlist_r[funi]), envir = envir)
      }
    } # expose_sigma_var_model_fun
    
    if(expose_sigma_basic_model_fun) {
      for (funi in 1:length(model$model_info$sigmabasicfunlist_r)) {
        # no use of attaching explicit namespace -> not possible also
        if(!grepl("::", model$model_info$sigmabasicfunlist_r[funi]) &
           !grepl(":::", model$model_info$sigmabasicfunlist_r[funi])) {
          assign(gsub("<-.*$", "", model$model_info$sigmabasicfunlist_r[funi]),
                 ept(model$model_info$sigmabasicfunlist_r[funi]), envir = envir)
        } # if(!grepl("::",
      } # for (funi in 1:l
    } # expose_sigma_basic_model_fun
    
    sigmavarSplineFun_name_val <- model$model_info[['sigmavarfunlist_r']]
    if(!is_emptyx(sigmavarSplineFun_name_val)) {
      model$model_info$funlist_r <- c( model$model_info$funlist_r,
                                       sigmavarSplineFun_name_val)
      
      for (funi in 1:length(sigmavarSplineFun_name_val)) {
        assign(gsub("<-.*$", "", sigmavarSplineFun_name_val[funi]),
               ept(sigmavarSplineFun_name_val[funi]), envir = envir)
      }
    } # if(!is_emptyx(sigmavarSplineFun_name_val) &
  } # if(expose_r_from_stan) {
  
  
  SplineFun_name      <- model$model_info[['StanFun_name']]
  sigmaSplineFun_name <- model$model_info[['sigmaStanFun_name']]
  spfun_collect       <- model$model_info$include_fun_names

  
  if(expose & backend == "rstan") {
    Spl_funs <- list()
    spfun_collectic <- -1
    for (spfun_collecti in spfun_collect) {
      spfun_collectic <- spfun_collectic + 1
      spfun_collecti_name <- spfun_collecti
      spfun_collecti_name_org <- spfun_collecti_name
      spfun_collecti_name <- gsub("_d0", "0", spfun_collecti_name)
      spfun_collecti_name <- gsub("_d1", "1", spfun_collecti_name)
      spfun_collecti_name <- gsub("_d2", "2", spfun_collecti_name)
      getfun_ <- spfun_collecti
      getfun_ <- eval(parse(text = getfun_), envir = envir)
      if(vectorize) getfun_ <- Vectorize(getfun_, SIMPLIFY = TRUE)
      assign(spfun_collecti_name, getfun_, envir = envir)
      Spl_funs[[paste0(spfun_collecti_name, "")]] <- getfun_
      if(grepl("_d", spfun_collecti_name_org)) {
        if(exists(spfun_collecti_name_org, envir = envir )) {
          remove(list=spfun_collecti_name_org, envir = envir)
        }
      }
    }
  } 
  
  
  if(expose & backend == "cmdstanr") {
    spfun_collect      <- c_model$functions
    Spl_funs <- list()
    for (spfun_collecti in spfun_collect$fun_names) {
      spfun_collecti_name <- spfun_collecti
      spfun_collecti_name_org <- spfun_collecti_name
      spfun_collecti_name <- gsub("_d0", "0", spfun_collecti_name)
      spfun_collecti_name <- gsub("_d1", "1", spfun_collecti_name)
      spfun_collecti_name <- gsub("_d2", "2", spfun_collecti_name)
      getfun_ <- spfun_collect[[spfun_collecti]]
      assign(spfun_collecti_name, getfun_, envir = envir)
      Spl_funs[[paste0(spfun_collecti_name, "")]] <- getfun_
      if(grepl("_d", spfun_collecti_name_org)) {
        if(exists(spfun_collecti_name_org, envir = envir )) {
          remove(list=spfun_collecti_name_org, envir = envir)
        }
      }
    }
  }
  
  
  if(expose_r_from_stan) {
    Spl_funs <- list()
    spfun_collectic <- -1
    for (spfun_collecti in spfun_collect) {
      spfun_collectic <- spfun_collectic + 1
      spfun_collecti_name <- spfun_collecti
      spfun_collecti_name <- gsub("_d0", "0", spfun_collecti_name)
      spfun_collecti_name <- gsub("_d1", "1", spfun_collecti_name)
      spfun_collecti_name <- gsub("_d2", "2", spfun_collecti_name)
      getfun_ <- spfun_collecti
      getfun__ <- deparse(ept(getfun_))
      gsub_it <- '_d0'
      gsub_by <- "0"
      getfun__ <- gsub(gsub_it, gsub_by, getfun__, fixed = T)
      getfun__ <- paste0(getfun__, collapse =  "\n")
      getfun__ <- eval(parse(text = getfun__), envir = envir)
      if(vectorize) getfun__ <- Vectorize(getfun__, SIMPLIFY = TRUE)
      Spl_funs[[paste0(spfun_collecti_name, "")]] <- getfun__
    }
  } 
  

  if(!expose & !expose_r_from_stan) {
    Spl_funs <- NULL
  }
  
  
  allSplineFun_name <- SplineFun_name
  
  if(expose_sigma_ls_model_fun) {
    allSplineFun_name <- c(allSplineFun_name, sigmaSplineFun_name)
  } else if(expose_sigma_var_model_fun) {
    # allSplineFun_name <- c(allSplineFun_name, sigmavarSplineFun_name)
  } else if(expose_sigma_basic_model_fun) {
    sigmabasicSplineFun_name <- model$model_info[['sigmabasicSplineFun_name']]
    allSplineFun_name <- c(allSplineFun_name, sigmabasicSplineFun_name)
  } 
  
  
  # The 'sigma_basic_model_fun' are not part of the stan code but only R
  # Therefore, we need to attach them back to the exposed stan functions 
  if(expose) {
    if(expose_sigma_basic_model_fun) {
      for (funi in 1:length(model$model_info$sigmabasicfunlist_r)) {
        name_ix <- gsub("<-.*$", "", model$model_info$sigmabasicfunlist_r[funi])
        assign(name_ix,
               ept(model$model_info$sigmabasicfunlist_r[funi]), envir = envir)
        Spl_funs[[name_ix]] <- ept(model$model_info$sigmabasicfunlist_r[funi])
        environment(Spl_funs[[funi]]) <- envir
      }
    } # expose_sigma_basic_model_fun
  } # if(expose) {
  
  
  model$model_info[['namesexefuns']] <- allSplineFun_name
  model$model_info[['exefuns']]      <- Spl_funs
  
 
  scode_include <- brms::stancode(model)
  model$bmodel  <- scode_include
  
  if(expose_r_from_stan) {
    if(is.null(model$model_info[['expose_method']])) {
      model$model_info[['expose_method']] <- 'R'
    }
  } else {
    if(is.null(model$model_info[['expose_method']])) {
      model$model_info[['expose_method']] <- 'S'
    }
  }
  
  if(returnobj) {
    model$model <- model$bmodel
    return(invisible(model))
  } else {
    return(invisible(NULL))
  }
  
}



#' @rdname expose_model_functions
#' @export
expose_model_functions <- function(model, ...) {
  UseMethod("expose_model_functions")
}




