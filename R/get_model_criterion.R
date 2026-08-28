

#' @title Return Model-Fit Criteria for a Bayesian SITAR Model
#'
#' @description
#' \code{get_model_criterion()} is a wrapper around [add_model_criterion()] that
#' computes and returns model-fit criteria. See [add_model_criterion()] for
#' details and available arguments.
#'
#' In addition to the criteria supported by [add_model_criterion()]
#' (\code{'loo', 'waic', 'kfold', 'loo_subsample', 'bayes_R2', 'loo_R2',
#' 'marglik'}), \code{get_model_criterion()} also computes conditional and
#' marginal versions of:
#' \itemize{
#'   \item Bayesian \code{R2} (\code{bayes_R2_conditional} and
#'         \code{bayes_R2_marginal}) via [performance::r2_bayes()], and
#'   \item LOO-adjusted Bayesian \code{R2} (\code{loo_R2_conditional} and
#'         \code{loo_R2_marginal}) via [performance::r2_loo()].
#' }
#'
#' The LOO-adjusted \code{R2} is conceptually analogous to the adjusted
#' \code{R2} in classical regression. The marginal \code{R2} reflects the
#' variance explained by the fixed effects alone, whereas the conditional
#' \code{R2} reflects the variance explained by both the fixed and random
#' effects. See [performance::r2_bayes()] for further details.
#' 
#' @param reformat Logical indicating whether to round numeric variables in the
#'   output \code{data.frame} using [base::round()]. The default is \code{NULL},
#'   which is treated as \code{TRUE}. When \code{TRUE}, numeric variables are
#'   rounded to the number of decimal places specified by \code{digits}.
#'   
#' @param tibble_table Logical indicating whether to return the table as a
#'   \code{data.frame} (\code{FALSE}) or as a \code{tibble} (\code{TRUE}).
#'   Default is \code{FALSE}.
#'
#' @param print_table Logical indicating whether to print the table using
#'   [knitr::kable()] (\code{TRUE}) or return it as an object (\code{FALSE}).
#'   When \code{print_table = TRUE}, the table is printed and the function
#'   returns \code{invisible(NULL)}. Default is \code{FALSE}.
#'
#' @param add_attr Logical indicating whether to attach complex list elements of
#'   the criteria as attributes to the returned \code{data.frame}. Default is
#'   \code{FALSE}.
#'
#' @param ... Additional arguments passed to [add_model_criterion()].
#' 
#' @inheritParams compare_models.bgmfit
#' @inheritParams add_model_criterion.bgmfit
#' @inheritParams growthparameters.bgmfit
#' @inheritParams brms::bayes_R2.brmsfit
#' @inheritParams brms::add_criterion.brmsfit
#' @inheritParams brms::waic.brmsfit
#' @inheritParams fitted_draws.bgmfit
#' 
#' @return A \code{data.frame}. If \code{add_attr = TRUE}, an additional
#'   attribute named \code{"attr_object"} is attached to the returned data
#'   frame. This attribute contains the complex, nested list components from
#'   \code{add_model_criterion()} that cannot be represented as regular columns.
#' 
#' @rdname get_model_criterion
#' @export
#' 
#' @seealso [brms::add_loo], [brms::add_ic()], [brms::add_waic()],
#'   [brms::bayes_R2()]
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
#' # For illustration purposes, we make a copy of model with itself
#' # In the example below, get_model_criterion() should indicate no difference 
#' # between model_1 and model_2 as both these models are exactly identical
#' model_1 <- model
#' model_2 <- model
#' 
#' # Add model fit criteria (e.g., WAIC). 
#' out_1 <- get_model_criterion(model_1, criterion = c("waic"))
#' out_2 <- get_model_criterion(model_2, criterion = c("loo"))
#' 
#' # compare models model_1 and model_2
#' out_12 <- get_model_criterion(model_1, model_2, criterion = c("waic"))
#' 
#' # compare models model_1 and model_2 - model names: "mods[[1L]]" "mods[[2L]]"
#' mods <- list(model_1, model_2)
#' out_12 <- get_model_criterion(mods)
#' 
#' # compare models model_1 and model_2 - model names extracted as such
#' # Note that list() could be supplied as a string "list()" also.
#' out_12 <- get_model_criterion(list(model_1, model_2))
#' out_12 <- get_model_criterion("list(model_1, model_2)")
#' 
#' # compare different versions of R2, from brms and performance packages
#' set_criterion <- c("bayes_R2", "loo_R2", 
#' "bayes_R2_conditional", "bayes_R2_marginal", 
#' "loo_R2_conditional", "loo_R2_marginal")
#' out_R2 <- get_model_criterion(model_2, criterion = set_criterion)
#' 
#' }
#' 
get_model_criterion.bgmfit <- function(model,
                                       ...,
                                       criterion = "loo",
                                       ndraws = NULL,
                                       draw_ids = NULL,
                                       pointwise = FALSE,
                                       model_name = NULL,
                                       summary = TRUE,
                                       robust = FALSE,
                                       probs = c(0.025, 0.975),
                                       newdata = NULL,
                                       resp = NULL,
                                       cores = 1,
                                       expose_function = FALSE, 
                                       verbose = FALSE,
                                       reformat = NULL,
                                       tibble_table = FALSE,
                                       print_table = FALSE,
                                       digits = 3,
                                       add_attr = FALSE) {

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
  
  conf_level <- probs[2] - probs[1]
    
  conf <- conf_level
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', "Est.Error", probtitles)
  
  
  add_args <- as.list(match.call(expand.dots = FALSE))
  defaults_it <- base::as.list(base::formals(get_model_criterion))
  for (i in names(defaults_it)) {
    if(is.null(add_args[[i]])) add_args[[i]] <- defaults_it[[i]]
  }
  
  defaults <- base::as.list(base::formals(get_model_criterion.bgmfit))
  defaults[['model']] <- NULL
  build_args <- utils::modifyList(defaults, add_args)

  check_criterion <- TRUE
  add_criterion_args                      <- build_args
  add_criterion_args[["model_name"]]      <- NULL
  add_criterion_args[["check_criterion"]] <- FALSE
  add_criterion_args[["compare"]]         <- FALSE
  add_criterion_args[["return_criteria"]] <- TRUE
  add_criterion_args[["return_model"]]    <- FALSE
  add_criterion_args[["add_attr"]]         <- NULL
  
  add_criterion_args <- move_to_front_list(add_criterion_args,c("model", "..."))
  
  add_criterion_args[['model']] <- NULL
  add_criterion_args[['...']] <- NULL
  add_criterion_args[['compare']] <- FALSE
  
  if (!is.list(add_criterion_args)) {
    stop("Argument 'add_criterion_args' must be a named list")
  }
  
  exprs <- as.list(substitute(list(model, ...)))[-1]
  vals <- c(list(model), list(...))
  models <- unlist(lapply(vals, flatten_models), recursive = FALSE)
  if (is.null(model_names)) {
    model_names <- unlist(Map(flatten_exprs, exprs, vals), 
                          use.names = FALSE)
  }
  
  only_one_model <- FALSE
  if (length(models) < 2) {
    only_one_model <- TRUE
  }
  
  only_one_criterion <- FALSE
  if (length(criterion) < 2) {
    only_one_criterion <- TRUE
  }
  
  if (!check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion_multiple(fit, criterion)) {
        stop("No precomputed criterion availabel for one or more models.", 
             " Either add criterion before hand using 'add_model_criterion()'", 
             " or else set check_criterion = TRUE. Note that arguments to", 
             " 'add_model_criterion()' function can be set by using", 
             " 'add_criterion_args' which must be a named list")
      }
    })
  }
  
  add_model_criterion_criterion <- c('loo', 'waic', 'kfold', 
                                     'loo_subsample', 'bayes_R2', 'loo_R2', 
                                     'marglik')
  
  performance_criterion  <- c("bayes_R2_conditional", "bayes_R2_marginal",
                              "loo_R2_conditional", "loo_R2_marginal")
  
  add_criterion_args_ele <- eval(add_criterion_args[['criterion']])
  
  all_criterion <- c(add_model_criterion_criterion, performance_criterion)
  
  add_criterion_args_ele_c <- c()
  for (add_criterion_args_elei in add_criterion_args_ele) {
    if(!add_criterion_args_elei %in% all_criterion) {
      add_criterion_args_ele_c <- c(add_criterion_args_ele_c, 
                                    add_criterion_args_elei)
    }
  }
  
  if(!is_emptyx(add_criterion_args_ele_c)) {
    stop2c("Following criterion are invalid: ", 
           collapse_comma(add_criterion_args_ele_c),
           ". Allowed criterion are: ",
           collapse_comma(all_criterion))
  }
  
  add_criterion_args[['criterion']] <- 
    add_criterion_args_ele[!add_criterion_args_ele %in% performance_criterion]
  rm('add_criterion_args_ele')
  
  call_add_criterion <- TRUE
  if(is_emptyx(add_criterion_args[['criterion']])) {
    call_add_criterion <- FALSE
  }
 
  performance_criterion_args <- list()
  performance_criterion_args[["robust"]]  <- robust
  performance_criterion_args[["conf"]]    <- conf
  performance_criterion_args[["verbose"]] <- verbose
  
  call_r2_bayes <- call_r2_loo <- FALSE
  add_bayes_R2_conditional <- add_bayes_R2_marginal <- FALSE
  add_bayes_R2_both <- add_loo_R2_conditional <- FALSE
  add_loo_R2_marginal <- add_loo_R2_both <- FALSE
  
  if("bayes_R2_conditional" %in% criterion | 
     "bayes_R2_marginal" %in% criterion) {
    call_r2_bayes <- TRUE
  }
  if("bayes_R2_conditional" %in% criterion) {
    add_bayes_R2_conditional <- TRUE
  }
  if("bayes_R2_marginal" %in% criterion ) {
    add_bayes_R2_marginal <- TRUE
  }
  if(add_bayes_R2_conditional & add_bayes_R2_marginal) {
    add_bayes_R2_both <- TRUE
  }
  
  if("loo_R2_conditional" %in% criterion | 
     "loo_R2_marginal" %in% criterion) {
    call_r2_loo <- TRUE
  }
  if("loo_R2_conditional" %in% criterion) {
    add_loo_R2_conditional <- TRUE
  }
  if("loo_R2_marginal" %in% criterion ) {
    add_loo_R2_marginal <- TRUE
  }
  if(add_loo_R2_conditional & add_loo_R2_marginal) {
    add_loo_R2_both <- TRUE
  }
  
  if(!call_add_criterion & !call_r2_bayes & !call_r2_loo) {
    stop2c("No valid criterion specified")
  }
  
  if(call_r2_bayes | call_r2_loo) {
    insight::check_if_installed("performance")
  }
  
  make_r2_bayes_loo_out <- function(df, criterion, set_names_) {
    Component <- NULL;
    df_out <- as.data.frame(df )
    select_vars <- c("R2", "SD", "CI_low", "CI_high", "Component")
    df_out <- df_out %>% dplyr::select(dplyr::all_of(select_vars))
    df_out <- df_out %>%
      dplyr::rename(!!as.symbol(set_names_[1]) := 
                      dplyr::all_of('R2')) %>% 
      dplyr::rename(!!as.symbol(set_names_[2]) := 
                      dplyr::all_of('SD')) %>% 
      dplyr::rename(!!as.symbol(set_names_[3]) := 
                      dplyr::all_of('CI_low')) %>% 
      dplyr::rename(!!as.symbol(set_names_[4]) := 
                      dplyr::all_of('CI_high')) 
    df_out_conditional <- df_out %>%
      dplyr::filter(Component == "conditional") %>% 
      dplyr::select(-dplyr::all_of("Component"))
    df_out_marginal <- df_out %>% 
      dplyr::filter(Component == "marginal") %>% 
      dplyr::select(-dplyr::all_of("Component"))
    out <- list()
    out[[paste0(criterion, "_", 'conditional')]] <- df_out_conditional
    out[[paste0(criterion, "_",  'marginal')]] <- df_out_marginal
    return(out)
  }
  
  build_out_list <- function(out, df_out_list, criterion_name) {
    out <- c(out, 
             setNames(list(df_out_list[[criterion_name]]), criterion_name)
             )
    return(out)
  }
  
  if (check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion_multiple(fit, criterion)) {
        fit <- do.call(expose_model_functions, 
                       c(list(model = fit, expose = expose_function) ))
        
        if(call_add_criterion) {
          add_criterion_args[['clearenvfuns']] <- FALSE
          suppressWarnings({
            out <- do.call(add_model_criterion, c(list(model = fit), 
                                                  add_criterion_args))
          })
        } else {
          out <- list()
        }
        
        if(call_r2_bayes) {
          df_out <- do.call(performance::r2_bayes, 
                            c(list(model = fit), performance_criterion_args))
          df_out_list <- make_r2_bayes_loo_out(df_out, "bayes_R2", set_names_)
          if(add_bayes_R2_both) {
            out <- c(out, df_out_list)
          } else if(add_bayes_R2_conditional) {
            out <- build_out_list(out, df_out_list, 'bayes_R2_conditional')
          } else if(add_bayes_R2_marginal) {
            out <- build_out_list(out, df_out_list, 'bayes_R2_marginall')
          }
        }
        
        if(call_r2_loo) {
          suppressWarnings({
            df_out <- do.call(performance::r2_loo, 
                              c(list(model = fit), performance_criterion_args))
          })
          df_out_list <- make_r2_bayes_loo_out(df_out, "loo_R2", set_names_)
          if(add_loo_R2_both) {
            out <- c(out, df_out_list)
          } else if(add_loo_R2_conditional) {
            out <- build_out_list(out, df_out_list, 'loo_R2_conditional')
          } else if(add_loo_R2_marginal) {
            out <- build_out_list(out, df_out_list, 'loo_R2_marginal')
          }
        }
        
      }
      
      out
    })
  }
  

  if (length(model_names) != length(models)) {
    nnames <- length(model_names)
    model_names_all <- paste0("model", seq_along(models) - 
                                length(model_names))
    model_names_all[1:nnames] <- model_names
    model_names <- model_names_all
    message2c("The number of model names is not same as the number of models.\n              The remaining models are named sequentially as model1,...")
  }

  out <- nested_to_df(models, model_names = model_names, add_attr = add_attr,
                      summary = summary, robust = robust, probs = probs,
                      verbose = verbose)
  
  if(is.null(reformat)) {
    reformat <- TRUE
  }
  
  if(reformat) {
    out <- out %>% 
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), 
                                  ~ round(.x, digits = digits)))
  }
  
  if(tibble_table) out <- out %>% tibble::as_tibble()
  
  if(print_table) {
    print(knitr::kable(out))
    return(invisible(NULL))
  } else {
    return(out)
  }
  
  return(invisible(NULL))
}




#' @rdname get_model_criterion
#' @export
get_model_criterion <- function(model, ...) {
  UseMethod("get_model_criterion")
}


#' @rdname get_model_criterion
#' @export
get_model_criterion.list <- function(model, ...) {
  get_model_criterion.bgmfit(model, ...)
}


#' @rdname get_model_criterion
#' @export
get_model_criterion.character <- function(model, ...) {
  get_model_criterion.bgmfit(model, ...)
}


#' @rdname get_model_criterion
#' @export
get_model_criterion.default <- function(model, ...) {
  if(!inherits(model, 'bgmfit') & 
     !inherits(model, 'list') &
     !inherits(model, 'character'))
  stop(
    "`model` must be an object of class 'bgmfit', a list, or a string",
    call. = FALSE
  )
}



