


#' Expose user defined Stan function for post-processing
#' 
#' @description The \strong{expose_model_functions()} is a wrapper around the
#'   [rstan::expose_stan_functions()] to expose user defined
#'   \code{Stan} function(s). These exposed functions are needed during the
#'   post-processing of the posterior draws.
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @param scode A character string (\code{Stan code}) containing the
#'   user-defined Stan function(s). If \code{NULL} (default), the \code{scode}
#'   is retrieved from the \code{model}.
#' 
#' @param expose A logical (default \code{TRUE}) to indicate whether to expose
#' functions and add them to the \code{model} as an attribute.
#' 
#' @param select_model A character string (default \code{NULL}) to indicate the
#'   model name. This is for internal use only.
#' 
#' @param returnobj A logical (default \code{TRUE}) to indicate whether to
#'   return the model object. When \code{expose = TRUE}, then it is advisable to
#'  set \code{returnobj = TRUE} too.
#'
#' @param vectorize A logical (default \code{FALSE}) to indicate whether the
#'   exposed functions should be vectorized via [base::Vectorize()]. Note that
#'   currently \code{vectorize} should be set to \code{FALSE} because setting it
#'   \code{TRUE} may not work as expected.
#' 
#' @inherit growthparameters.bgmfit params
#' 
#' @param ... Additional arguments passed to the
#'   [rstan::expose_stan_functions()] function.
#'
#' @return An object of class \code{bgmfit} if \code{returnobj=TRUE}, otherwise
#'   invisible \code{NULL}.
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
#' # To save time, argument expose is set as FALSE which runs a dummy test 
#' # and avoid model compilation which often takes time
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
  
  fun_env <- envir
  
  if(is.null(select_model)) select_model <- model$model_info$select_model
  
  if(!expose) {
    if (is.null(model$model_info$decomp))  expose_r_from_stan <- TRUE
    if (!is.null(model$model_info$decomp)) expose_r_from_stan <- FALSE
  } else {
    expose_r_from_stan <- FALSE
  }
  
  
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
    rstan::expose_stan_functions(rstan::stanc(model_code = exposecode), 
                                 env = fun_env)
  }
  
  
  if(expose_r_from_stan) {
    for (funi in 1:length(model$model_info$funlist_r)) {
      assign(gsub("<-.*$", "", model$model_info$funlist_r[funi]),
             ept(model$model_info$funlist_r[funi]), envir = envir)
    }
  }
  
  
  SplineFun_name <- model$model_info[['StanFun_name']]
  spfun_collect <- c(SplineFun_name,
                     paste0(SplineFun_name, "_", 
                            c("d0", 
                              "d1",
                              "d2")))
  
  
  if(expose) {
    additionlsfuns <- c('getX')
    if(model$model_info[['select_model']] == 'sitar' |
       model$model_info[['select_model']] == 'rcs') {
      additionlsfuns <- c(additionlsfuns, 'getKnots')
    }
    spfun_collect <- c(spfun_collect, additionlsfuns)
  }
  
  
  if(expose_r_from_stan) {
    spfun_collect <- c(spfun_collect, 'getX')
    if(select_model == 'sitar' | select_model == 'rcs') {
      spfun_collect <- c(spfun_collect, 'getKnots')
    }
  }
  
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
  
  
  model$model_info[['namesexefuns']] <- SplineFun_name
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

