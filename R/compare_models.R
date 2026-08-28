

#' Model comparison with the loo package
#' 
#' @description 
#' A wrapper around [brms::loo_compare()] to compare models using fit criteria.
#' 
#' @details
#' All models should have pre computed fit criterion See [add_model_criterion()]
#' for more help.
#' 
#'
#' @param model An \code{R} object of class \code{bsitar} or a list of models.
#' 
#' @param check_criterion A logical (default \code{FALSE}) indicating whether to
#'   check and add fit criterion if not already added to model objects
#'   
#' @param add_criterion_args A named list to set up the arguments for the
#'   [add_model_criterion()] function. Ignored when \code{check_criterion =
#'   FALSE} or when each model objects has pre computed fit criteria.
#'   
#' @param expose_model_functions_args A named list to set up the arguments for the
#'   [expose_model_functions()] function. Ignored when \code{check_criterion =
#'   FALSE} or when each model objects has pre computed fit criteria.
#'   
#' @param verbose A logical (default \code{FALSE}) to print some useful
#'   information.
#'   
#' 
#' @inheritParams brms::loo_compare
#' @inherit brms::loo_compare return
#' @inheritParams add_model_criterion
#'
#'
#' @rdname compare_models
#' @export
#' 
#' @seealso [brms::loo_compare()]
#' 
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid model estimation which can take time, the Bayesian SITAR model fit
#' # to the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' model <- getNsObject(berkeley_exfit)
#' 
#' # For illustration purposes, we compare model with itself
#' # In the example below, compare_models() should indicate no difference 
#' # between model_1 and model_2 as both these models are exactly identical
#' 
#' # Add model fit criteria (e.g., WAIC). 
#' model_1 <- add_model_criterion(model, criterion = c("waic"))
#' model_2 <- add_model_criterion(model, criterion = c("waic"))
#' 
#' # compare models model_1 and model_2
#' compare_model_12 <- compare_models(model_1, model_2, criterion = c("waic"))
#' 
#' # compare models model_1 and model_2 - model names: "mods[[1L]]" "mods[[2L]]"
#' mods <- list(model_1, model_2)
#' compare_model_12 <- compare_models(mods)
#' 
#' # compare models model_1 and model_2 - model names extracted as such
#' compare_model_12 <- compare_models(list(model_1, model_2))
#' compare_model_12 <- compare_models("list(model_1, model_2)")
#' 
#' }
#' 
compare_models.bgmfit <- function(model, 
                                 ..., 
                                 criterion = "loo", 
                                 model_name = NULL,
                                 check_criterion = TRUE,
                                 add_criterion_args = list(),
                                 expose_model_functions_args = list(),
                                 expose_function = FALSE,
                                 verbose = FALSE) {

  only_object <- NULL
  if (is.character(model) && length(model) == 1 && 
      grepl("^\\s*list\\s*\\(", model)) {
    only_object <- FALSE
  } 
  if(is.list(model)) {
    only_object <- FALSE
  } 
  if(is.bgmfit(model)) {
    only_object <- TRUE
  }
  if(length(c(list(model), list(...))) > 1) {
    only_object <- FALSE
  }
  
  if (is.character(model) && length(model) == 1 && 
      grepl("^\\s*list\\s*\\(", model)) {
    model_list_str <-  model
    model <- eval(parse(text = model), envir = parent.frame())
  } else {
    model_list_str <-  deparse(substitute(model))
  }
  gsub_namespace <- FALSE
  model_list_names <- NULL
  if(grepl("list\\(", model_list_str)) {
    model_list_str <- gsub("list", "", model_list_str)
    model_list_str <- strsplit(model_list_str, ",")[[1]]
    # replace :: with _ , not when :::
    model_list_str <- gsub("[^A-Za-z0-9_:]", "", model_list_str)
    if(gsub_namespace) {
      tmp <- gsub(":::", "@@@COLONPAIR@@@", model_list_str)
      tmp <- gsub("::", "_", tmp)
      model_list_str <- gsub("@@@COLONPAIR@@@", ":::", tmp)
    }
    model_list_names <- model_list_str
  }
  
  if (is.null(model_name)) {
    if(!is.null(model_list_names)) model_names <- model_list_names
    if( is.null(model_list_names)) model_names <- NULL
  } else {
    model_names <- model_name
  }
  
  if(!is.list(add_criterion_args)) {
    stop("Argument 'add_criterion_args' must be a named list")
  }
 
  exprs <- as.list(substitute(list(model, ...)))[-1]
  vals  <- c(list(model), list(...))
  
  models <- unlist(lapply(vals, flatten_models), recursive = FALSE)
  
  if (is.null(model_names)) {
    model_names <- unlist(Map(flatten_exprs, exprs, vals), use.names = FALSE)
  }
  
  if (length(models) < 2) {
    stop("Need at least two models for loo_compare().")
  }
  
  if(!check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion(fit, criterion)) {
        stop("No precomputed criterion availabel for one or more models.",
             " Either add criterion before hand using 'add_model_criterion()'",
             " or else set check_criterion = TRUE. Note that arguments to",
             " 'add_model_criterion()' function can be set by using",
             " 'add_criterion_args' which must be a named list")
      }
    })
  } 
  
  if(check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion(fit, criterion)) {
        fit <- do.call(
          expose_model_functions,
          c(list(model = fit, expose = expose_function), 
            expose_model_functions_args))
      }
      fit
    })
  } 
  
  if(check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion(fit, criterion)) {
        fit <- do.call(
          expose_model_functions,
          c(list(model = fit, expose = expose_function), 
            expose_model_functions_args))
        suppressWarnings({
        fit <- do.call(
          brms::add_criterion,
          c(list(x = fit, criterion = criterion), add_criterion_args))
        })
      }
      fit
    })
  }
  
  if (length(model_names) != length(models)) {
    nnames <- length(model_names)
    model_names_all <- paste0("model", seq_along(models)-length(model_names))
    model_names_all[1:nnames] <- model_names
    model_names <- model_names_all
    message2c("The number of model names is not same as the number of models.  
              The remaining models are named sequentially as model1,...")
  }
  
  do.call(brms::loo_compare, c(models, 
                               list(criterion = criterion, 
                                    model_names = model_names)))

}



#' @rdname compare_models
#' @export
compare_models <- function(model, ...) {
  UseMethod("compare_models")
}


#' @rdname compare_models
#' @export
compare_models.list <- function(model, ...) {
  compare_models.bgmfit(model, ...)
}


#' @rdname compare_models
#' @export
compare_models.character <- function(model, ...) {
  compare_models.bgmfit(model, ...)
}


#' @rdname compare_models
#' @export
compare_models.default <- function(model, ...) {
  if(!inherits(model, 'bgmfit') & 
     !inherits(model, 'list') &
     !inherits(model, 'character'))
    stop(
      "`model` must be an object of class 'bgmfit', a list, or a string",
      call. = FALSE
    )
}
  

#' @rdname compare_models
#' @export
compare_model <- function(model, ...) {
  UseMethod("compare_models")
}



