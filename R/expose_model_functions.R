


#' @title Expose User-Defined Stan Functions for Post-Processing
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
#' @export expose_model_functions.bgmfit
#' @export
#'
#' @seealso [rstan::expose_stan_functions()]
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
                                 envir = NULL,
                                 ...) {
  
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- parent.frame()
  }
  
  sigmanlflf <- FALSE
  if(!is.null(model$model_info[['sigmaStanFun_name']])) {
    sigmanlflf <- TRUE
  }
  
  fun_env <- envir
  
  if(is.null(select_model)) select_model <- model$model_info$select_model
  
  if(!expose) {
    if (is.null(model$model_info$decomp))  expose_r_from_stan <- TRUE
    if (!is.null(model$model_info$decomp)) expose_r_from_stan <- FALSE
  } else {
    expose_r_from_stan <- FALSE
  }
  
  
  # mcall <- match.call()
  # arguments <- as.list(mcall)[-1]
  
  #dots <- list(...)
  #arguments <- c(arguments, dots)

  
  if(expose) {
    if (verbose) {
      setmsgtxt <-
        paste0("\n Exposing Stan functions for post-processing\n")
      message(setmsgtxt)
    }
  }
  
  if(expose) {
    if (is.null(scode)) {
      exposecode <- brms::stancode(model)
    } else if (!is.null(scode)) {
      exposecode <- scode
    }
    # rstan::expose_stan_functions(rstan::stanc(model_code = exposecode), 
    #                              env = fun_env)
    
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
        stop(paste("The 'setcompiler' argument should be 'stanc' or 'stan_model'", 
             "or else NULL", collapse =","))
      }
    }
    
    stanc_arguments[['setcompiler']] <- NULL
    stan_model_arguments[['setcompiler']] <- NULL
    
    if(compiled_code_via == 'stanc') {
      compiled_code <- do.call(rstan::stanc, stanc_arguments)
    }
    if(compiled_code_via == 'stan_model') {
      compiled_code <- do.call(rstan::stan_model, stan_model_arguments)
    }
    
    compiled_code_args <- list()
    compiled_code_args[['env']] <- fun_env
    compiled_code_args[['stanmodel']] <- compiled_code
    compiled_code_args[['includes']] <- NULL
    compiled_code_args[['show_compiler_warnings']] <- FALSE
    
    do.call(rstan::expose_stan_functions, compiled_code_args)
  }
  
  
  

  if(expose_r_from_stan) {
    for (funi in 1:length(model$model_info$funlist_r)) {
      assign(gsub("<-.*$", "", model$model_info$funlist_r[funi]),
             ept(model$model_info$funlist_r[funi]), envir = envir)
    }
    
    if(sigmanlflf) {
      for (funi in 1:length(model$model_info$sigmafunlist_r)) {
        assign(gsub("<-.*$", "", model$model_info$sigmafunlist_r[funi]),
               ept(model$model_info$sigmafunlist_r[funi]), envir = envir)
      }
    }
  } # if(expose_r_from_stan) {
  
  
  SplineFun_name <- model$model_info[['StanFun_name']]
  # spfun_collect <- c(SplineFun_name,
  #                    paste0(SplineFun_name, "_", 
  #                           c("d0", 
  #                             "d1",
  #                             "d2")))
  # print(spfun_collect)
  # stop()
  
  spfun_collect <- model$model_info$include_fun_names

  
  if(sigmanlflf) {
    sigmaSplineFun_name <- model$model_info[['sigmaStanFun_name']]
    sigmaspfun_collect <- c(sigmaSplineFun_name,
                       paste0(sigmaSplineFun_name, "_", 
                              c("d0", 
                                "d1",
                                "d2")))
    
    spfun_collect <- c(spfun_collect, sigmaspfun_collect)
  }
  

  if(expose) {
    additionlsfuns <- c()
    if(sigmanlflf) {
      sigmaadditionlsfuns <- c('sigmagetX')
      additionlsfuns <- c(additionlsfuns, sigmaadditionlsfuns)
    }
    if(model$model_info[['select_model']] == 'sitar' |
       model$model_info[['select_model']] == 'rcs') {
      if(sigmanlflf) {
        sigmaadditionlsfuns <- c('sigmagetKnots')
        additionlsfuns <- c(additionlsfuns, sigmaadditionlsfuns)
      }
    }
    spfun_collect <- c(spfun_collect, additionlsfuns)
  }


  if(expose_r_from_stan) {
    spfun_collect <- c(spfun_collect)
    if(sigmanlflf) {
      spfun_collect <- c(spfun_collect, 'sigmagetX')
    }
    if(select_model == 'sitar' | select_model == 'rcs') {
      if(sigmanlflf) {
        spfun_collect <- c(spfun_collect, 'sigmagetKnots')
      }
    }
  }
  
  
  
  
  # if(expose) {
  #   additionlsfuns <- c('getX')
  #   if(sigmanlflf) {
  #     sigmaadditionlsfuns <- c('sigmagetX')
  #     additionlsfuns <- c(additionlsfuns, sigmaadditionlsfuns)
  #   }
  #   if(model$model_info[['select_model']] == 'sitar' |
  #      model$model_info[['select_model']] == 'rcs') {
  #     additionlsfuns <- c(additionlsfuns, 'getKnots')
  #     if(sigmanlflf) {
  #       sigmaadditionlsfuns <- c('sigmagetKnots')
  #       additionlsfuns <- c(additionlsfuns, sigmaadditionlsfuns)
  #     }
  #   }
  #   spfun_collect <- c(spfun_collect, additionlsfuns)
  # }
  # 
  # 
  # if(expose_r_from_stan) {
  #   spfun_collect <- c(spfun_collect, 'getX')
  #   if(sigmanlflf) {
  #     spfun_collect <- c(spfun_collect, 'sigmagetX')
  #   }
  #   if(select_model == 'sitar' | select_model == 'rcs') {
  #     spfun_collect <- c(spfun_collect, 'getKnots')
  #     if(sigmanlflf) {
  #       spfun_collect <- c(spfun_collect, 'sigmagetKnots')
  #     }
  #   }
  # }
  
  
  
  
  
  
  
  nys <- model$model_info$nys
  ys <- model$model_info$ys
  if(nys > 1) {
    spfun_collect2 <- c()
    for (ysii in ys) {
      tempysi <- paste0(ysii, "_", spfun_collect)
      spfun_collect2 <- c(spfun_collect2, tempysi)
    }
    spfun_collect <- spfun_collect2
  }
  
  
  
  if(expose) {
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
  
  
  
  if(!expose & !expose_r_from_stan) Spl_funs <- NULL
  
  allSplineFun_name <- SplineFun_name
  if(sigmanlflf) {
    allSplineFun_name <- c(allSplineFun_name, sigmaSplineFun_name)
  }
  
  model$model_info[['namesexefuns']] <- allSplineFun_name
  # model$model_info[['namesexefuns']] <- SplineFun_name
  model$model_info[['exefuns']]      <- Spl_funs
 
  
  scode_include <- brms::stancode(model)
  model$bmodel <- scode_include
  if (nys == 1 | nys > 1) {
    for (nys__i in 1:nys) {
      cont_ <- 0
      for (cont_i in 0:2) {
        cont_ <- cont_ + 1
        if (nys == 1) {
          gsubit <- paste0(
            "vector",
            " ",
            paste0("", "", SplineFun_name),
            "_",
            "d",
            cont_i,
            paste0(".*end of spline function", "_", ys[nys__i],
                   "d", cont_i, "")
          )
        } else if (nys > 1) {
          gsubit <-
            paste0(
              "vector",
              " ",
              paste0(ys[nys__i], "_", SplineFun_name),
              "_",
              "d",
              cont_i,
              paste0(".*end of spline function", "_", ys[nys__i],
                     "d", cont_i, "")
            )
        }
        scode_include <-
          gsub(gsubit, "", scode_include, fixed = F)
      }
    }
  }
  
  if(returnobj) {
    model$model <- model$bmodel
    return(invisible(model))
  } else {
    return(invisible(NULL))
  }
  
}



#' @rdname expose_model_functions.bgmfit
#' @export
expose_model_functions <- function(model, ...) {
  UseMethod("expose_model_functions")
}

