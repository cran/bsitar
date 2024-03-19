


#' An internal function to get arguments from the call
#'
#' @param arguments A list of default function arguments.
#' @param xcall A character string specifying the name of the calling function.
#' @keywords internal
#' @return A list comprised of function arguments.
#' @keywords internal
#' @noRd
#'
get_args_ <- function(arguments, xcall) {
  `%!in%` <- Negate(`%in%`)
  f_funx_arg <- formals(paste0(xcall, ".", 'bgmfit'))
  nf_funx_arg_names <-
    intersect(names(arguments), names(f_funx_arg))
  arguments <-
    c(arguments, f_funx_arg[names(f_funx_arg) %!in% nf_funx_arg_names])
  arguments
}



#' An internal function to deparse a symbol argument and remove spaces
#'
#' @param deparseobj A symbol
#' @keywords internal
#' @return A character string.
#' @noRd
#'
deparse_0 <- function(deparseobj) {
  deparseobj <- paste(deparse(deparseobj), collapse = "")
  deparseobj <- gsub("[[:space:]]", "", deparseobj)
  deparseobj
}


#' An internal function to substitute and deparse a symbol argument
#'
#' @param deparseobj A symbol
#' @keywords internal
#' @return A character string.
#' @noRd
#'
deparse_0s <- function(deparseobj) {
  deparseobj <- paste(deparse(substitute(deparseobj)), collapse = "")
  deparseobj
}


#' An internal function to remove spaces from the string
#'
#' @param deparseobj A character string
#' @keywords internal
#' @return A character string.
#' @noRd
#'

gsub_space <- function(deparseobj) {
  deparseobj <- gsub("[[:space:]]", "", deparseobj)
  deparseobj
}



#' An internal function to expose function after optimization
#'
#' @param model An object of class \code{bgmfit}.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'

expose_optimize_fit <- function(model,
                                subset_list = NULL,
                                expose_function = T) {
  optimize_fit_models <-  model

  if(!is.null(subset_list)) {
    if(!is.numeric(subset_list)) stop("models must a numeric vector")
    optimize_fit_models <- optimize_fit_models[subset_list]
  } else {
    optimize_fit_models <- optimize_fit_models
  }

  m_list <- list()
  for (il in 1:length(optimize_fit_models)) {
    if(is.null(expose_function)) {
      m_list[[il]] <- optimize_fit_models[[il]]
    } else if(!is.null(expose_function)) {
      if(expose_function) {
        message("Exposing Stan function...", " (for model no. ", il, ")")
        m_list[[il]] <-
          expose_model_functions(optimize_fit_models[[il]],
                                optimize_fit_models[[il]]$bmodel,
                                expose = TRUE,
                                select_model = NULL,
                                returnobj = TRUE,
                                envir = NULL)
      } else if(!expose_function) {
        m_list[[il]] <- optimize_fit_models[[il]]
      }
    }
  }
  m_list <- m_list[!sapply(m_list, is.null)]
  m_list
}



#' An internal function to models after optimization
#'
#' @param model An object of class \code{bgmfit}.
#' @keywords internal
#' @return A list comprised of plot objects.
#' @noRd
#'

plot_optimize_fit <- function(model,
                              subset_list = NULL,
                              what = "plot",
                              expose_function = FALSE,
                              print= TRUE, ...) {

  optimize_fit_models <-  model
  dots <- list(...)
  for (i in names(dots)) {
    if(!i %in% formalArgs(plot_curves.bgmfit))
      stop("arguments must be be one of the following",
           "\n ",
           formalArgs(plot_curves.bgmfit))
  }

  if(!is.null(subset_list)) {
    if(!is.numeric(subset_list)) stop("models must a numeric vector")
    optimize_fit_models <- optimize_fit_models[subset_list]
  } else {
    optimize_fit_models <- optimize_fit_models
  }

  m_list <- list()
  for (il in 1:length(optimize_fit_models)) {
    if(is.null(expose_function)) {
      m_list[[il]] <- optimize_fit_models[[il]]
    } else if(!is.null(expose_function)) {
      if(expose_function) {
        message("Exposing Stan function...", " (for model no. ", il, ")")
        m_list[[il]] <-
          expose_model_functions(model = optimize_fit_models[[il]],
                                scode = optimize_fit_models[[il]]$bmodel,
                                expose = TRUE,
                                select_model = NULL,
                                returnobj = TRUE,
                                envir = NULL)
      } else if(!expose_function) {
        m_list[[il]] <- optimize_fit_models[[il]]
      }
    }
  }

  m_list <- m_list[!sapply(m_list, is.null)]
  nx <- function(.x, bx, args_) {
    message("Working on model no. ", .x)
    dots$model <- bx[[.x]]
    dots$... <- NULL
    if(is.null(what)) what <- 'plot'
    if(what == "plot") {
      out_ <- do.call(plot_curves, dots)
      title_ <- bx[[.x]]$model_info$optimization_info
      out_ <- out_ + ggplot2::labs(title = title_)
    }

    if(what == "growthparameters") {
      out_ <- do.call(growthparameters, dots)
    }

    if(!is.null(print)) {
      if(print) {
        print(out_)
      }
    }
    return(out_)
  }
  out <- purrr::map(1:length(m_list), ~nx(.x, m_list, args_))
  return(out)
}


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





#' An internal function to evaluate arguments ending with _str suffix
#'
#' @param tsx An argument with _str suffix.
#' @param data A data frame.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'

get_gr_str_coef_id <- function(tsx,
                               data) {
  tsx <- strsplit(tsx, "+(", fixed = T)[[1]]
  tsx_id_w_or_wo_gr <- c()
  for (tsx_id_w_or_wo_gri in 1:length(tsx)) {
    tsx_id_w_or_wo_gr_get <- get_x_random2_asitis(tsx[tsx_id_w_or_wo_gri])
    tsx_id_w_or_wo_gr <- c(tsx_id_w_or_wo_gr, tsx_id_w_or_wo_gr_get)
  }
  tsx <- gsub("(", "", tsx, fixed = T)
  tsx <- gsub(")", "", tsx, fixed = T)
  tsx_c_coef  <- tsx_c_id    <- set_form_gr_it      <- list()
  set_ncov_it <- set_corr_it <- set_corr_true_false <- list()
  for (i in 1:length(tsx)) {
    tsx_c <- strsplit(tsx[i], "|", fixed = T)[[1]]
    set_corr_it_get <- tsx_c[2]
    tsx_c1 <- tsx_c[1]
    tsx_c3 <- tsx_c[3]
    if(!grepl("^~", tsx_c1)) tsx_c1 <- paste0("~", tsx_c1)
    if(grepl("^~0", tsx_c1)) set_form_0_gr <- TRUE
    if(grepl("^~1", tsx_c1)) set_form_0_gr <- FALSE
    set_form_gr <- tsx_c1
    tsx_c1_mat <- eval(parse(text = paste0(
      "model.matrix(",
      tsx_c1, ",data = data)"
    )))
    if (ncol(tsx_c1_mat) == 1)
      nlcov <- NULL
    else
      nlcov <- ncol(tsx_c1_mat) - 1
    nlcovoefnames <- colnames(tsx_c1_mat)
    nlcovoefnames <- gsub("\\(|)", "", nlcovoefnames)
    if(length(nlcovoefnames) > 1) tsx_c3 <- rep(tsx_c3, length(nlcovoefnames))
    tsx_c_coef[[i]] <- nlcovoefnames
    tsx_c_id[[i]]   <- tsx_id_w_or_wo_gr[i] # tsx_c3[1]
    set_form_gr_it[[i]]   <- set_form_gr
    if(set_form_0_gr) {
      set_ncov_it_get <- length(nlcovoefnames)
    }
    if(!set_form_0_gr) {
      if(length(nlcovoefnames) == 1) set_ncov_it_get <- NULL
      if(length(nlcovoefnames) > 1) set_ncov_it_get <- length(nlcovoefnames) - 1
    }
    if(is.null(set_ncov_it_get)) {
      set_corr_true_false[[i]] <- FALSE
    } else if(!is.null(set_ncov_it_get)) {
      if(set_corr_it_get == "") set_corr_true_false[[i]] <- FALSE
      if(set_corr_it_get != "") set_corr_true_false[[i]] <- TRUE
    }
    set_corr_it[[i]] <- set_corr_it_get
    set_ncov_it[[i]] <- set_ncov_it_get
  } # for (i in 1:length(tsx)) {

  if(length(tsx_c_coef) != length(tsx_c_id))
    stop("coef and id length should be same")
  list(tsx_c_coef = tsx_c_coef, tsx_c_id = tsx_c_id,
       set_form_gr_it = set_form_gr_it, set_ncov_it = set_ncov_it,
       set_corr_it = set_corr_it, set_corr_true_false = set_corr_true_false)
}




#' An internal function to get corr structure from || syntax for
#'  arguments ending with _str suffix
#'
#' @param str_id_all_list An argument with _str suffix for \code{id}.
#' @param str_corr_all_list An argument with _str suffix for \code{gr_cor}.
#' @param str_corr_tf_all_list An argument with _str suffix for \code{corr}.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'

get_str_corr_tf_function_new_better <- function(str_id_all_list,
                                                str_corr_all_list,
                                                str_corr_tf_all_list) {

  if(length(str_id_all_list) > 0 ) {
    id_corr_tf_bind <- cbind(unlist(str_id_all_list),
                             unlist(str_corr_all_list),
                             unlist(str_corr_tf_all_list))

    checkdi_c <- group_id_unique <- str_corr_tf <- c()

    for (checkdi in 1:length(str_id_all_list)) {
      checkdi_c <- c(checkdi_c,  length(str_corr_all_list[[checkdi]]) )
    }

    for (id_corr_tf_bind_1i in unique(id_corr_tf_bind[ , 1])) {
      temp_mat <- id_corr_tf_bind[which(id_corr_tf_bind == id_corr_tf_bind_1i),]
      if(!is.matrix(temp_mat)) temp_mat <- matrix(temp_mat) %>% t()
      get_check_id_mat      <- temp_mat[ , 1]
      get_check_corr_mat    <- temp_mat[ , 1]
      get_check_corr_tf_mat <- temp_mat[ , 1]
      if(any(grepl("TRUE", get_check_corr_tf_mat)))
        str_corr_tf_single <- TRUE else str_corr_tf_single <- FALSE
      if(max(table(get_check_corr_mat)) > 1) {
        if(!str_corr_tf_single) str_corr_tf_single <- TRUE
      }
      str_corr_tf <- c(str_corr_tf, str_corr_tf_single)
      group_id_unique <- c(group_id_unique, id_corr_tf_bind_1i)
    }
    list(str_corr_tf = str_corr_tf, group_id_unique = group_id_unique)
  } else {
    list(str_corr_tf = NULL, group_id_unique = NULL)
  }
}



#' An internal function to append priors to the bpriors
#'
#' @param tempx A prior object.
#' @keywords internal
#' @return A prior object.
#' @noRd
#'

extract_prior_str_lv <- function(tempx) {
  if(!is.list(tempx) & !is.vector(tempx)) {
    out_prior_str <- tempx
  } else if(is.vector(tempx) & length(tempx) == 1 &
            !grepl("c\\(", tempx) &
            !grepl("list\\(", tempx) ) {
    out_prior_str <- tempx
  } else if(is.vector(tempx) & length(tempx) == 1 &
            grepl("list\\(", tempx) ) {
    tempx2 <- str2lang(tempx)
    tempx2[[1]] <- NULL
    out_prior_str <- c()
    for (tempx2i in 1:length(tempx2)) {
      get_it_ <- tempx2[[tempx2i]] %>% deparse()
      get_it_ <- gsub("\"", "", get_it_)
      out_prior_str <- c(out_prior_str, get_it_)
    }
  } else if(is.list(tempx) | is.vector(tempx)  ) {
    tempx2 <- str2lang(tempx)
    tempx2[[1]] <- NULL
    out_prior_str <- c()
    for (tempx2i in 1:length(tempx2)) {
      get_it_ <- tempx2[[tempx2i]] %>% deparse()
      get_it_ <- gsub("\"", "", get_it_)
      out_prior_str <- c(out_prior_str, get_it_)
    }
  }
  out_prior_str
}




#' An internal function to restore parantheses in formuale objects
#'
#' @param strx A formual object.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'

restore_paranthese_grgr_str_form <- function(strx) {
  restore_paranthese_grgr_str <- function(strx2) {
    if(!grepl("gr", strx2, fixed = T)) {
      strx_ <- strx2
      if(grepl("|", strx2, fixed = T) & !grepl("^\\(", strx2, fixed = F)) {
        strx_ <- paste0("(", strx2, "")
      }
    } else if(grepl("|gr", strx2)) {
        pattern <- "gr\\s*(.*?)\\s*,"
        strx2_check_ <- regmatches(strx2, regexec(pattern, strx2))
        strx2_check_ <- strx2_check_[[1]][2]
        if(grepl("_", strx2_check_)) {
          stop("Underscore '_' is not allowed in the variable name when",
               "\n ",
               " defining the group identifier using the 'gr()' formulation",
               "\n ",
               " please check '", strx2_check_, "' varibale in the random ",
               "formula for '", sub("\\~.*", "", strx), "'"
          )
        }
      if(!grepl("|gr(", strx2, fixed = T)) {
        strx_ <- gsub("|gr" , "|gr(", strx2, fixed = T)
        strx_ <- paste0("(", strx_, ")")
      } else if(grepl("^|gr\\(", strx2)) {
        strx_ <- paste0("(", strx2, "")
      }
    } else if(grepl("|gr", strx2)) {
      strx_ <-  paste0("(", strx2, "")
    }
    strx_
  }
  abx_c <- c()
  strx <- gsub("+(", "_xxxx_", strx, fixed = T)
  abxs <- strsplit(strx, "_xxxx_", fixed = T)[[1]]
  for (abxi in 1:length(abxs)) {
    abx_c[abxi] <- restore_paranthese_grgr_str(abxs[abxi])
    abx_c[abxi]
  }
  abx_c <- paste0(abx_c, collapse = "+")
  abx_c
}


#' An internal function to get random effect formula arguments
#'
#' @param x A character string of random effect formula.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'

get_x_random2 <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  x <- strsplit(x, ")+" )[[1]]
  x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
  if(any(grepl("^|gr", x))) {
    x <- sub(".*gr", "", x)
    x_c <- c()
    for (xi in 1:length(x)) {
      gxi <- strsplit(x[xi], ",")[[1]][1]
      x_c <- c(x_c, gxi)
    }
    x <- x_c
    # x <- strsplit(x, ",")[[1]][1]
  }
  x <- sub(".*\\|", "", x)
  x <- unique(unlist(strsplit(x, ":")) )
  x
}


#' An internal function to get random effect formula arguments with tilde sign
#'
#' @param x A character string of random effect formula.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'

get_x_random2_asitis <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  x <- strsplit(x, ")+" )[[1]]
  x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
  if(any(grepl("^|gr", x))) {
    x <- sub(".*gr", "", x)
    x <- strsplit(x, ",")[[1]][1]
  }
  x <- sub(".*\\|", "", x)
  x
}



#' An internal function to get object enclosed within the parenthesis
#'
#' @param x A character string.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'

get_o_paranthesis <- function(x) {
  if(!grepl("lf\\(", x)) {
    x <- gsub("^lf\\(", "", x)
    x <- gsub(")$", "", x)
  }
  if(!grepl("nlf\\(", x)) {
    x <- gsub("^nlf\\(", "", x)
    x <- gsub(")$", "", x)
  }
  x <- strsplit(x, "~")[[1]][2]
  x
}


#' An internal function to get object enclosed within the parenthesis without
#'  parenthesis.
#'
#' @param x A character string.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'

get_o_paranthesis2 <- function(x) {
  x <- gsub("^\\(", "", x)
  x <- gsub(")$", "", x)
  x <- strsplit(x, "~")[[1]][2]
  x
}


#' An internal function to get covariates from the formula.
#'
#' @param x A character string.
#' @keywords internal
#' @return A vector comprised of character strings.
#' @noRd
#'

getcovlist <- function(x) {
  if (is.character(x))
    x <- x
  else
    x <- deparse(x)
  x <- gsub("~", "", gsub("\\s", "", x))
  x <- strsplit(x, "+", fixed = T)[[1]]
  if (length(x) == 1)
    x <- NULL
  else
    x <- x[-1]
  return(x)
}


#' An internal function to parse and evaluate a character string.
#'
#' @param x A character string.
#' @keywords internal
#' @return An evaluated object.
#' @noRd
#'

ept <- function(x) eval(parse(text = x), envir = parent.frame())


#' An internal function to get parameter names from the stancode.
#'
#' @param code A character string of stancode.
#' @param full A logical (default \code{TRUE}) indicating whether to get full
#' names.
#' @param section A character string specifying the Stan block
#' (default \code{parameters})
#' @param what A character string specifying the name of a particular parameter.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'

get_par_names_from_stancode <- function(code,
                                        full = TRUE,
                                        section =  'parameters',
                                        semicolan = FALSE,
                                        what = '') {
  regex_for_section <- paste(".*(",section,"\\s*\\{.*?\\}).*", sep = '')
  filtered_stan_code <- gsub(code, pattern = regex_for_section,
                             replacement = "\\1")

  zz <- strsplit(filtered_stan_code, "\n")[[1]][-1]
  collect <- c()
  collect_full <- c()
  for (i in 1:length(zz)-1) {
    if(!(identical(zz[i], character(0))))  {
      if(!semicolan) t <- sub(";.*", "", zz[i])
      if( semicolan) t <- sub(";.*", ";", zz[i])
      t_full <- t
      t_full <- gsub("^ *|(?<= ) | *$", "", t_full, perl=T)
      if(what == "") {
        get_t_full <- t_full
        collect_full <- c(collect_full, get_t_full)
      } else if(what != "") {
        get_t_full <- t_full
        get_t_full <- get_t_full[grepl(paste0(what, "_"), get_t_full)]
        collect_full <- c(collect_full, get_t_full)
      }
      t <- tail(strsplit(t,split=" ")[[1]],1)
      collect <- c(collect, t)
    }
  }
  if(!full) {
    out <- collect
  }
  if(full) {
    out <- collect_full
  }
  out
}



#' An internal function to get/set the number of cores
#'
#' @param cores.arg A character string specifying cores argument from the
#' function.
#' @keywords internal
#' @return A list comprised of integers.
#' @noRd
#'

get.cores <- function(cores.arg) {
  cores_ <- eval(cores.arg)
  if (!is.null(cores_)) {
    if (cores_ == "maximise") {
      max.cores <-
        as.numeric(future::availableCores(methods = "system", omit = 0))
      if (max.cores < 1)
        max.cores <- 1
    } else if (cores_ == "optimize") {
      max.cores <-
        as.numeric(future::availableCores(methods = "system", omit = 1))
      if (max.cores < 1)
        max.cores <- 1
    } else {
      max.cores <- eval(cores_)
    }
  } else if (is.null(cores_)) {
    max.cores <- NULL
  }

  if (!is.null(cores_)) {
    if (Sys.info()["sysname"] == "Windows") {
      .cores_ps <- 1
    } else {
      .cores_ps <- max.cores
    }
  } else if (is.null(cores_)) {
    .cores_ps <- 1
  }

  list(max.cores = max.cores, .cores_ps = .cores_ps)
}




#' An internal function to validate the response variable
#'
#' @param model An object of class \code{bgmfit}.
#' @param resp A character string specifying the name of the response variable.
#' Default \code{NULL}.
#' @keywords internal
#' @return An error if evaluation fails.
#' @noRd
#'

validate_response <- function(model,
                              resp = NULL) {
  if (model$model_info$nys == 1 & !is.null(resp)) {
    stop(
      "You have fit a univariate model",
      " but set resp option as ",
      resp,
      ".",
      "\n ",
      " The resp option should be appropriately set to NULL",
      "\n ",
      " (i.e., resp = NULL)"
    )
  }
  if (model$model_info$nys > 1 & is.null(resp)) {
    if (!is.na(model$model_info$univariate_by)) {
      stop(
        "You have fit a univariate-by-subset model for ",
        model$model_info$univariate_by,
        "\n ",
        " but dit not set the the resp options appropriately",
        " (which is NULL at present).",
        "\n ",
        " The response options are ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
    if (model$model_info$multivariate) {
      stop(
        "You have fit a multivariate model ",
        "\n ",
        " but dit not set the the resp options appropriately",
        " (which is NULL at present).",
        "\n ",
        " The response options are ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
  }

  if (!is.null(resp)) {
    if (!resp %in% model$model_info[['ys']]) {
      stop(
        "Response should be one of the following: ",
        paste(model$model_info[['ys']], collapse = " "),
        "\n ",
        " but you have specified: ",
        resp
      )
    }
  }

}





#' An internal function to set up the priors when fitting a model with 3 or
#' more levels of hierarchy.
#'
#' @param new_prior_list A prior object.
#' @keywords internal
#' @return A prior object.
#' @noRd
#'

setup_higher_priors <- function(new_prior_list) {

  ##############################################
  # Initiate non formalArgs()
  ##############################################
  . <- NULL;
  o_l <- list()
  ixi = 0
  group_ <- class_ <- nlpar_ <- resp_ <- cor_check <- sd_check <- NA
  group_ <-
    class_ <- nlpar_ <- resp_ <- cor_check <- sd_check <- c()
  for (new_prior_listi in 1:length(new_prior_list)) {
    ixi <- ixi + 1
    pstr <- new_prior_list[[new_prior_listi]]
    if (is.null(pstr[['prior']]))
      prior_i <- ''
    else
      prior_i <- pstr[['prior']]
    if (is.null(pstr[['group']]))
      group_i <- ''
    else
      group_i <- pstr[['group']]
    if (is.null(pstr[['class']]))
      class_i <- ''
    else
      class_i <- pstr[['class']]
    if (is.null(pstr[['nlpar']]))
      nlpar_i <- ''
    else
      nlpar_i <- pstr[['nlpar']]
    if (is.null(pstr[['coef']]))
      coef_i  <- ''
    else
      coef_i  <- pstr[['coef']]
    if (is.null(pstr[['resp']]))
      resp_i  <- ''
    else
      resp_i  <- pstr[['resp']]
    if (is.null(pstr[['lb']]))
      lb_i    <- ''
    else
      lb_i    <- pstr[['lb']]
    if (is.null(pstr[['ub']]))
      ub_i    <- ''
    else
      ub_i <- pstr[['ub']]
    if (is.null(pstr[['dpar']]))
      dpar_i  <- ''
    else
      dpar_i  <- pstr[['dpar']]

    prior_i <- gsub("[[:space:]]", "", prior_i)
    group_i <- gsub("[[:space:]]", "", group_i)
    class_i <- gsub("[[:space:]]", "", class_i)
    nlpar_i <- gsub("[[:space:]]", "", nlpar_i)
    coef_i  <- gsub("[[:space:]]", "", coef_i)
    resp_i  <- gsub("[[:space:]]", "", resp_i)
    lb_i    <- gsub("[[:space:]]", "", lb_i)
    ub_i    <- gsub("[[:space:]]", "", ub_i)
    dpar_i  <- gsub("[[:space:]]", "", dpar_i)

    group_ <- c(group_, group_i)
    class_ <- c(class_, class_i)
    nlpar_ <- c(nlpar_, nlpar_i)
    resp_ <- c(nlpar_, resp_i)
    if (class_i == 'sd')
      sd_check <- c(sd_check, group_i)
    if (class_i == 'cor')
      cor_check <- c(cor_check, group_i)

    if (lb_i == '' & ub_i == '') {
      o_l[[ixi]] <- prior_string(
        prior_i,
        group = group_i,
        class = class_i,
        nlpar = nlpar_i,
        coef = coef_i,
        resp = resp_i,
        # lb = lb_i,
        # ub = ub_i,
        dpar = dpar_i
      )
    } else if (lb_i != '' | ub_i != '') {
      o_l[[ixi]] <- prior_string(
        prior_i,
        group = group_i,
        class = class_i,
        nlpar = nlpar_i,
        coef = coef_i,
        resp = resp_i,
        # lb = lb_i,
        # ub = ub_i,
        dpar = dpar_i
      )
    } # if(lb_i == '' & ub_i == '' ) {
  } # for (new_prior_listi in 1:length(new_prior_list)) {
  o_l %>%  do.call(rbind, .)
}



#' An internal function to rename patterns in a character vector.
#' This is adapted from the brms package.
#'
#' @param x a character vector to be renamed
#' @param pattern the regular expressions in x to be replaced
#' @param replacement the replacements
#' @param fixed same as for 'gsub'
#' @param check_dup: logical; check for duplications in x after renaming
#' @param ... passed to 'gsub'
#' @keywords internal
#' @return renamed character vector of the same length as x
#' @noRd
#'
rename <- function(x,
                   pattern = NULL,
                   replacement = NULL,
                   fixed = TRUE,
                   check_dup = FALSE, ...) {
  pattern <- as.character(pattern)
  replacement <- as.character(replacement)
  if (!length(pattern) && !length(replacement)) {
    # default renaming to avoid special characters in coeffcient names
    pattern <- c(
      " ", "(", ")", "[", "]", ",", "\"", "'",
      "?", "+", "-", "*", "/", "^", "="
    )
    replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ")
  }
  if (length(replacement) == 1L) {
    replacement <- rep(replacement, length(pattern))
  }
  stopifnot(length(pattern) == length(replacement))
  has_chars <- nzchar(pattern)
  pattern <- pattern[has_chars]
  replacement <- replacement[has_chars]
  out <- x
  for (i in seq_along(pattern)) {
    out <- gsub(pattern[i], replacement[i], out, fixed = fixed, ...)
  }
  dup <- duplicated(out)
  if (check_dup && any(dup)) {
    dup <- x[out %in% out[dup]]
    stop2("Internal renaming led to duplicated names. \n",
          "Occured for: ", collapse_comma(dup))
  }
  out
}




#' An internal function to evaluate priors specified in data block of Stan
#'
#' @param model An object of class \code{bgmfit}
#' @param spriors A prior object. If \code{NULL} (default),
#'   [brms::prior_summary()] is used to \code{spriors} from the  \code{model}
#' @param sdata A Stan data object. If \code{NULL} (default), [brms::standata()]
#'   is used to get \code{sdata} from the  \code{model}.
#' @param prior_name_asit A logical (default \code{FALSE}) to indicate whether
#'   prior names should be returned as it is from the stancode.
#' @param gsub_group A character vector specifying the group identifier that
#'   will be removed from the \code{group} column of the prior object. Default
#'   \code{NULL}.
#' @param sort_response A character vector specifying the order of response
#'   variables that will be used in sorting the \code{resp} column in the prior
#'   object. Default \code{NULL}.
#' @param sort_parameter A character vector specifying the order of parameter
#'   names that will be used in sorting the \code{nlpar} column in the prior
#'   object. Default \code{NULL}.
#' @param sort_coefficient A character vector specifying the order of
#'   coefficient names that will be used in sorting the \code{nlpar} column in
#'   the prior object. Default \code{NULL}.
#' @param sort_class A character vector specifying the order of class names that
#'   will be used in sorting the \code{class} column in the prior object.
#'   Default \code{NULL}.
#' @param digits An integer to set the \code{digits} argument for the
#'   \code{round} function.
#' @param viewer A logical (default \code{FALSE}) to indicate whether to display
#'   the output in R viewer. Currently ignored to avoid dependency on the 'gt'
#'   package.
#' @param raw A logical (default \code{FALSE}) to indicate whether to return
#'   the output in original format.
#'   
#' @keywords internal
#' @return A data frame object.
#' @noRd
#'
priors_to_textdata <- function(model,
                                spriors = NULL,
                                sdata = NULL,
                                prior_name_asit = FALSE,
                                gsub_coef = NULL,
                                gsub_group = NULL,
                                sort_response = NULL,
                                sort_group = NULL,
                                sort_parameter = c(letters[1:26], "sigma"),
                                sort_coefficient = c("Intercept"),
                                sort_class = c("b", "sd", "cor"),
                                digits = 2,
                                viewer = FALSE,
                                raw = FALSE
                               ) {
  arguments <- as.list(match.call())[-1]

  if (missing(model)) {
    model <- NULL
  }

  ##############################################
  # Initiate non formalArgs()
  ##############################################
  nlpar <- NULL;
  coef <- NULL;
  class <- NULL;
  prior <- NULL;
  group <- NULL;
  resp <- NULL;
  dpar <- NULL;
  Response <- NULL;
  Coefficient <- NULL;
  Parameter <- NULL;
  Group <- NULL;
  Class <- NULL;
  . <- NULL;


  if (is.null(model) & is.null(spriors) & is.null(sdata)) {
    stop("Supply either model or spriors and sdata arguments")
  } else if (!is.null(model) &
             !is.null(spriors) & !is.null(sdata)) {
    stop("Supply only model or spriors and sdata arguments")
  } else if (!is.null(model)) {
    spriors <- brms::prior_summary(model)
    sdata <- brms::standata(model)
  } else if (is.null(model)) {
    if (is.null(spriors) & is.null(sdata)) {
      stop("Supply spriors and sdata arguments")
    }
    if (is.null(spriors) & is.null(sdata)) {
      stop("Supply spriors and sdata arguments")
    }
  }

  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }

  if(!raw) spriors <- spriors %>% dplyr::filter(source == 'user')
  if( raw) prior_name_asit <- TRUE

  env_ <- environment()
  list2env(sdata, envir =  env_)

  for (i in 1:nrow(spriors)) {
    getxit <- spriors[i, ]$prior
    if(getxit == "") getxit <- "flat"
    prior_name <- strsplit(getxit, "\\(")[[1]][1]
    if (!prior_name_asit) {
      if (!is.na(prior_name) & prior_name == 'lkj') {
        prior_name_case <- toupper(prior_name)
      } else if (!is.na(prior_name) & prior_name == 'lkj_corr_cholesky') {
        prior_name_case <- 'LKJ'
      } else {
        if(!raw) prior_name_case <- firstup(prior_name)
      }
    }


    if (prior_name_asit) prior_name_case <- prior_name


    getxit_2 <-
      regmatches(getxit, gregexpr("(?<=\\().*?(?=\\))", getxit, perl = T))[[1]]

    # This is for flat priors
    if(identical(getxit_2, character(0))) {
      getxit_7 <- paste0(prior_name_case, '')
    }  else if(!identical(getxit_2, character(0))) {
      getxit_3 <- strsplit(getxit_2, ",")[[1]]
      getxit_4 <- sapply(getxit_3, function(x)
        eval(parse(text = x)))
      getxit_4 <- round(getxit_4, digits = digits)
      getxit_5 <- paste(getxit_4, collapse = ", ")
      getxit_6 <- paste0("(", getxit_5, ")")
      getxit_7 <- paste0(prior_name_case, getxit_6)
    } else {
      getxit_7 <- NULL
    }
    spriors[i, ]$prior <- getxit_7
  }
  
  
  if(raw) {
    return(spriors)
  }

  spriors <-
    spriors %>% data.frame() %>% dplyr::select(-c('lb', 'ub', 'source'))
  spriors <- spriors %>% `rownames<-`(NULL)
  spriors <-
    spriors %>%  dplyr::mutate(class =  dplyr::if_else(class == 'L', 'cor',
                                                       class))


  if (!is.null(gsub_coef)) {
    for (gsub_coefi in gsub_coef) {
      spriors <-
        spriors %>%  dplyr::mutate(coef = gsub(gsub_coefi, "" , coef))
    }
  }

  if (!is.null(gsub_group)) {
    for (gsub_groupi in gsub_group) {
      spriors <-
        spriors %>%  dplyr::mutate(group = gsub(gsub_groupi, "" , group))
    }
  }
  
  
  
  spriors <- spriors %>% dplyr::relocate(nlpar, coef,
                                         class, prior,
                                         group, resp,
                                         dpar)

  # for sigma betas
  spriors <-
    spriors %>%  dplyr::mutate(coef =  dplyr::if_else(coef == '' &
                                                        class == 'Intercept',
                                                      class, coef))

  spriors <-
    spriors %>%  dplyr::mutate(
      class =  dplyr::if_else(
        class == 'Intercept' &
          dpar == 'sigma' &
          class == 'Intercept',
        'b',
        class
      )
    )



  spriors <-
    spriors %>%  dplyr::mutate(nlpar =  dplyr::if_else(nlpar == '' &
                                                         dpar != '',
                                                       dpar, nlpar)) %>%
    dplyr::select(-'dpar')


  spriors <- spriors %>% dplyr::rename(
    Parameter = nlpar,
    Coefficient = coef,
    Class = class,
    Prior = prior,
    Group = group,
    Response = resp
  )


  if(is.null(sort_response)) {
    if (!is.null(model)) {
      if(length(model$model_info$nys) > 1) {
        sort_response <- model$model_info$ys
      }
    }
  }


  spriors <- spriors %>%
    dplyr::arrange(match(Response, sort_response)) %>%
    dplyr::arrange(match(Coefficient, sort_coefficient)) %>%
    dplyr::arrange(match(Parameter, sort_parameter)) %>%
    dplyr::arrange(match(Group, sort_group)) %>%
    dplyr::arrange(match(Class, sort_class))

  if (!is.null(model)) {
    if (is.na(model$model_info$univariate_by) &
        !model$model_info$multivariate) {
      spriors <- spriors %>%  dplyr::select(-'Response')
    }
  }

  return(spriors)
}



#' An internal function to split vector at factor indices
#'
#' @param x A vector.
#' @param pos A vector of indices.
#' @keywords internal
#' @return A vector.
#' @noRd
#'

splitAt2 <- function(x, pos) {
  x <- droplevels(x)
  out <- c()
  pos2 <- c(1, pos, length(x)+1)
  for (i in seq_along(pos2[-1])) {
    out[[i]] <- x[pos2[i]:(pos2[i+1]-1)]
  }
  return(out)
}



#' An internal function to Negate R's in function
#' @param `%in%` R's in function
#' @keywords internal
#' @return An R function.
#' @noRd
#'

`%!in%` <- Negate(`%in%`)


#' An internal function to evaluate NULL and length zero arguments
#' @param x A symbol (argument)
#' @param y A symbol (argument)
#' @keywords internal
#' @return An R function.
#' @noRd
#'

'%||%' <- function(x, y) {
  if (is.null(x)) x <- y
  x
}


#' An internal function to customize R's stop function
#' @param ... An argument
#' @keywords internal
#' @return A string (error message) from R's warning() function.
#' @noRd
#'
stop2 <- function(...) {
  stop(..., call. = FALSE)
}



#' An internal function to customize R's warning function
#' @param ... An argument
#' @keywords internal
#' @return A string (warning message) from R's warning() function.
#' @noRd
#'
warning2 <- function(...) {
  warning(..., call. = FALSE)
}


#' An internal function to collapse elements of vector separated by a comma
#' @param ... An argument
#' @keywords internal
#' @return A character string.
#' @noRd
#'
collapse_comma <- function(...) {
  paste0("'", ..., "'", collapse = ", ")
}



#' An internal function to edit stancode for NCP parametarization
#' @param stancode A string character of stan code
#' @param genq_only A logical (default \code{FALSE}) to indicate whether to
#' return only the generated quantity sub code.
#' @param normalize A logical (default \code{TRUE}) to indicate whether to
#' include the normalizing constant in the prior target density.
#' @keywords internal
#' @return A character string.
#' @noRd
#'
edit_scode_ncp_to_cp <- function(stancode,
                                 genq_only = FALSE,
                                 normalize = TRUE) {

  # Rename transformed parameters and parameters for ease of processing
  true_name_tp  <- 'transformed parameters'
  true_name_p   <- 'parameters'
  tempt_name_tp <- 'transformed_parameters_'
  tempt_name_p  <- 'parameters_'

  clines_tp <- get_par_names_from_stancode(stancode,
                                           section =  true_name_tp,
                                           semicolan = TRUE,
                                           full = TRUE)

  clines_p <- get_par_names_from_stancode(stancode,
                                          section =  true_name_p,
                                          semicolan = TRUE,
                                          full = TRUE)

  clines_m <- get_par_names_from_stancode(stancode,
                                          section =  'model',
                                          semicolan = TRUE,
                                          full = TRUE)


  editedcode    <- stancode
  editedcode    <- gsub(true_name_tp, tempt_name_tp, editedcode, fixed = T)
  editedcode    <- gsub(true_name_p,  tempt_name_p,  editedcode, fixed = T)

  editedcode2 <- editedcode

  clines_tp2 <- c()
  for (il in clines_tp) {
    il <- gsub(pattern = "//", replacement = "//", x = il, fixed = T)
    il <- gsub(pattern = "//[^\\\n]*", replacement = "", x = il)
    if(!grepl('^lprior', gsub_space(il)) &
       !grepl('^reallprior', gsub_space(il)) &
       !grepl('^-', gsub_space(il))) {
      if(!is_emptyx(il)) {
        clines_tp2 <- c(clines_tp2, il)
      }
    }
  }

  clines_tp <- clines_tp2

  how_many_r_1 <- 0
  move_to_p <- move_to_m <- c()
  for (clines_tpi in clines_tp) {
    for (igr in 1:100) {
      if(grepl(paste0("r_", igr), clines_tpi)) {
        if(grepl("^matrix", clines_tpi) |
           grepl("^array", clines_tpi) # This for z_1 when only one sd
        ) {
          how_many_r_1 <- how_many_r_1 + 1
          move_to_p <- c(move_to_p, clines_tpi)
          clines_tpi <- ""
        }
        if(!grepl(paste0("^", "r_", igr, " = "), clines_tpi) &
           !grepl(paste0("^", "r_", igr, "="), clines_tpi) &
           !grepl("//", clines_tpi) &
           clines_tpi != "") {
          move_to_m <- c(move_to_m, clines_tpi)
        }
      }
    }
  }


  prepare_p <- c()
  for (clines_pi in clines_p) {
    if(grepl("z_", clines_pi)) {
      for (igr in 1:100) {
        if(grepl(paste0("z_", igr), clines_pi)) {
          if(grepl("^matrix", clines_pi) |
             grepl("^array", clines_pi) # This for z_1 when only one sd
          ) {
            if(grepl("^array", clines_pi)) {
              what_p <- paste0("// ", clines_pi)
            } else {
              what_p <- paste0("// ", clines_pi)
            }
          }
        }
      }
    } else {
      what_p <- paste0("", clines_pi)
    }
    prepare_p <- c(prepare_p, what_p)
  }




  move_to_p <- paste(move_to_p, collapse = "\n")
  prepare_p <- paste(prepare_p, collapse = "\n")
  prepare_p <- paste0(prepare_p, "\n", move_to_p)

  # Remove duplicate - this happens when some parms are single vector
  # For example, sigma_gr ~ 1
  prepare_p <- strsplit(prepare_p, "\n", fixed = T)[[1]] %>%
    data.frame() %>%
    dplyr::distinct() %>%
    unlist() %>%
    as.vector()

  prepare_p <- paste(prepare_p, collapse = "\n")

  ## Match data from regexpr()
  pattern_r <- pattern_N <- pattern_M <- pattern_sd <- pattern_L <-  c()
  for (rxi in 1:100) {
    pattern     <- paste0('r_', rxi)
    pattern_    <- regexpr(pattern, prepare_p)
    pattern_ri  <- regmatches(prepare_p, pattern_)
    pattern_Ni  <- gsub('r_', 'N_', pattern_ri, fixed = T)
    pattern_Mi  <- gsub('r_', 'M_', pattern_ri, fixed = T)
    pattern_Li  <- gsub('r_', 'L_', pattern_ri, fixed = T)
    pattern_sdi <- gsub('r_', 'sd_', pattern_ri, fixed = T)
    pattern_r   <- c(pattern_r, pattern_ri)
    pattern_N   <- c(pattern_N, pattern_Ni)
    pattern_M   <- c(pattern_M, pattern_Mi)
    pattern_L   <- c(pattern_L, pattern_Li)
    pattern_sd  <- c(pattern_sd, pattern_sdi)
  }


  # Add space to model block elements
  zz_c <- c()
  for (iz in move_to_m) {
    zz_c <- c(zz_c, paste0("  ", iz))
  }
  move_to_m <- paste(zz_c, collapse = '\n')

  # Add space to parameters elements
  zz <- strsplit(prepare_p, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    zz_c <- c(zz_c, paste0("  ", zz[iz]))
  }
  prepare_p <- paste(zz_c, collapse = '\n')

  if(normalize) {
    lprior_target <- "target"
  } else if(!normalize) {
    lprior_target <- "lprior"
  }

  m_n_c_l_c <- c()
  for (h1i in 1:length(pattern_r)) {
    pattern_ri  <- pattern_r[h1i]
    pattern_Ni  <- pattern_N[h1i]
    pattern_Mi  <- pattern_M[h1i]
    pattern_Li  <- pattern_L[h1i]
    pattern_sdi <- pattern_sd[h1i]
    m_n_c_l <-
      paste0("  for(i in 1:", pattern_Ni, ') {\n',
             "    ", lprior_target, " +=  multi_normal_cholesky_lpdf(",
             pattern_ri, "[i, ] |\n",
             "    rep_row_vector(0, ",
             pattern_Mi, "),\n",
             "    diag_pre_multiply(",
             pattern_sdi, ", ", pattern_Li, "));",
             "  \n  }"
      )
    m_n_c_l_c <- c(m_n_c_l_c, m_n_c_l)
  } # for (h1i in 1:length(pattern_r)) {

  m_n_c_l_c <- paste(m_n_c_l_c, collapse = "\n")

  for (il in clines_p) {
    if(!grepl("^array", il)) {
      editedcode2 <- gsub(pattern = "//", replacement = "//",
                          x = editedcode2, fixed = T)
      editedcode2 <- gsub(pattern = "//[^\\\n]*", replacement = "",
                          x = editedcode2)
      editedcode2 <- gsub(paste0(il, ""), "", editedcode2, fixed = T)
    }
  }

  for (il in clines_tp) {
    if(!grepl("^array", il)) {
      editedcode2 <- gsub(pattern = "//", replacement = "//",
                          x = editedcode2, fixed = T)
      editedcode2 <- gsub(pattern = "//[^\\\n]*", replacement = "",
                          x = editedcode2)
      editedcode2 <- gsub(paste0(il, ""), "", editedcode2, fixed = T)
    }
  }

  # Below to flatten code without any empty splace or lines
  # editedcode2 <- gsub("(?m)^\\h*\\R?", "", editedcode2, perl=TRUE)
  # editedcode2 <- gsub("\r", "", editedcode2, fixed=TRUE)
  p_block_syb_by <- paste0("", tempt_name_p, " {")
  p_block_syb_it <- paste0(p_block_syb_by, "\n", prepare_p)
  editedcode2 <- gsub(paste0("", p_block_syb_by), p_block_syb_it,
                      editedcode2, fixed=T, perl=F)


  # Remove empty lines
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
      # comment out to_vector(z_
      if(how_many_r_1 > 0) {
        if(grepl("to_vector(z_", zz_in, fixed = T))
          zz_in <- paste0("  //", zz_in)
        if(grepl("scale_r_cor(z_", zz_in, fixed = T))
          zz_in <- paste0("  //", zz_in)
      }
      zz_c <- c(zz_c, zz_in)
    }
  }

  editedcode2 <- paste(zz_c, collapse = '\n')

  add_to_model_block <- paste0(m_n_c_l_c, "\n", move_to_m)
  add_to_genq_block <- paste0( move_to_m)


  # lprior_code <- "real lprior = 0;"
  if(normalize) {
    lprior_code <- "model {"
  } else if(!normalize) {
    lprior_code <- "real lprior = 0;"
  }


  editedcode2 <- gsub(lprior_code, paste0(lprior_code, "\n",
                                          add_to_model_block, "\n"),
                      editedcode2, fixed = T)

  genq_code <- "generated quantities {"
  editedcode2 <- gsub(genq_code, paste0(genq_code, "\n",
                                        add_to_genq_block),
                      editedcode2, fixed = T)

  editedcode2 <- gsub(tempt_name_tp, true_name_tp, editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_p,  true_name_p,  editedcode2, fixed = T)

  # If only one random effects and hence no r_1 etc, then return original code
  if(identical(pattern_r, character(0))) {
    return(stancode)
  } else if(!identical(pattern_r, character(0))) {
    if(genq_only) return(add_to_genq_block)
    return(editedcode2)
  }

}


#' An internal function to get derivatives from distance curve for QR decomp
#' @param model An object of class \code{bgmfit}.
#' @param y0 A matrix comprised of distance curves.
#' @param newdata A data frame. If \code{NULL}, data used in original model
#' fit used.
#' @param deriv An integer (\code{1 or 2}) to specify derivative. Default
#'  \code{deriv = 1} estimates velocity curve whereas \code{deriv = 2} is to
#'  get acceleration curve.
#' @param probs The percentiles to be computed by the quantile function.
#' @param robust If FALSE (the default) the mean is used as the measure of
#' central tendency and the standard deviation as the measure of variability.
#' If TRUE, the median and the median absolute deviation (MAD) are applied
#' instead. Only used if summary is TRUE.
#' @keywords internal
#' @return A character string.
#' @noRd
#'
mapderivqr <- function(model,
                       y0,
                       newdata = NULL,
                       deriv = 1,
                       resp = NULL,
                       probs = c(0.025, 0.975),
                       robust = FALSE) {

  if(is.null(probs)) probs <- c(0.025, 0.975)
  if(is.null(robust)) robust <- FALSE

  if(is.null(newdata)) newdata <- model$data

  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }

  validate_response(model, resp)

  list_c <- list()
  xvar_ <- paste0('xvar', resp_rev_)
  yvar_ <- paste0('yvar', resp_rev_)
  groupvar_ <- paste0('groupvar', resp_rev_)
  xvar <- model$model_info[[xvar_]]
  yvar <- model$model_info[[yvar_]]
  hierarchical_ <- paste0('hierarchical', resp_rev_)

  idvar <- model$model_info[[groupvar_]]
  if(length(idvar) > 1) idvar <- idvar[1]
  yvar  <- 'yvar'


  if(!is.na(model$model_info$univariate_by)) {
    newdata <- newdata %>% data.frame() %>%
      dplyr::filter(!!as.symbol(model$model_info$univariate_by) == 'Male')
  }


  getdydx <- function (x, y, id, data, ndigit = 2) {
    ##############################################
    # Initiate non formalArgs()
    ##############################################
    sorder <- NULL;
    data$sorder <- as.numeric(row.names(data))
    .data <- data %>%
      dplyr::mutate(.x = !!dplyr::sym(x)) %>%
      dplyr::mutate(.y = !!dplyr::sym(y)) %>%
      dplyr::mutate(.id = !!dplyr::sym(id)) %>%
      data.frame()

    .dydx <- function(x, y) {
      n <- length(x); i1 <- 1:2; i2 <- (n - 1):n
      c(diff(y[i1])/diff(x[i1]), (y[-i1] - y[-i2])/(x[-i1] - x[-i2]),
        diff(y[i2])/diff(x[i2]))
    }
    dydx <- lapply(split(.data, as.numeric(.data$.id)),
                   function(x) {x$.v <- .dydx(x$.x, x$.y); x } )
    dydx <- do.call(rbind, dydx) %>% data.frame() %>% dplyr::arrange(sorder)
    return(round(dydx[[".v"]], ndigit))
  }

  mapderiv <- function(.xrow, x = xvar, y = yvar, id = idvar,
                       data = newdata) {
    newdata[[y]] <- .xrow
    getdydx(x = x, y = y, id = id, data = newdata)
  }

  if(deriv == 1) {
    tempx <- apply(y0, 1, mapderiv) %>% t()
  }
  if(deriv == 2) {
    tempx <- apply(y0, 1, mapderiv) %>% t()
    tempx <- apply(tempx , 1, mapderiv) %>% t()
  }
  dout <-  brms::posterior_summary(tempx , probs = probs, robust = robust)
  dout
}



#' Title Check if string, vector, or list is empty
#' @description Adapted from from is_empty.
#' See https://github.com/strengejacke/sjmisc/blob/master/R/is_empty.R
#' @param x String, character vector, list, data.frame or numeric vector or
#' factor.
#' @param first.only Logical, if \code{FALSE} and \code{x} is a character
#' vector, each element of \code{x} will be checked if empty. If
#' \code{TRUE}, only the first element of \code{x} will be checked.
#' @param all.na.empty Logical, if \code{x} is a vector with all \code{NA},
#' \code{is_emptyx} will return \code{FALSE} if \code{all.na.empty = FALSE},
#' and will return \code{TRUE} if \code{all.na.empty = TRUE} (default).
#' @return Logical, \code{TRUE} if \code{x} is a character vector or string
#' and is empty, \code{TRUE} if \code{x} is a vector or list and of length 0,
#' \code{FALSE} otherwise.
#' @keywords internal
#' @noRd
#'
is_emptyx <- function(x, first.only = TRUE, all.na.empty = TRUE) {
  if (!is.null(x)) {
    if (is.character(x)) {
      if (length(x) == 0) return(TRUE)
      zero_len <- nchar(x) == 0
      if (first.only) {
        zero_len <- .is_truex(zero_len[1])
        if (length(x) > 0) x <- x[1]
      } else {
        return(unname(zero_len))
      }
    } else if (is.list(x)) {
      x <- purrr::compact(x)
      zero_len <- length(x) == 0
    } else {
      zero_len <- length(x) == 0
    }
  }
  any(is.null(x) || zero_len || (all.na.empty && all(is.na(x))))
}


#' Title Check if TRUE or False
#'
#' @param x String, character vector
#'
#' @return Logical, \code{TRUE} / \code{FALSE}
#' @keywords internal
#' @noRd
#'
.is_truex <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x) && x
}




# Commenting out for CRAN initial release

#' Fit model via cmdstanr
#'
#' @param scode A character string of model code
#' @param sdata A list of data objects
#' @param brm_args A list of argument passes to the [[brms::brm()]]
#'
#' @return An object of class \code{bgmfit}
#' @keywords internal
#' #noRd
#'

# brms_via_cmdstanr <- function(scode, sdata, brm_args) {
#   if(!is.null(brm_args$threads$threads)) {
#     stan_threads <- TRUE
#   } else {
#     stan_threads <- FALSE
#   }
#
#   if(!is.null(brm_args$opencl)) {
#     stan_opencl <- TRUE
#   } else {
#     stan_opencl <- FALSE
#   }
#
#   cpp_options <- list(stan_threads = stan_threads,
#                       stan_opencl = stan_opencl)
#
#
#   stanc_options <- brm_args$stan_model_args$stanc_options
#
#
#   if(brm_args$silent == 0) {
#     show_messages = TRUE
#     show_exceptions = TRUE
#   }
#   if(brm_args$silent == 1) {
#     show_messages = TRUE
#     show_exceptions = FALSE
#   }
#   if(brm_args$silent == 2) {
#     show_messages = FALSE
#     show_exceptions = FALSE
#   }
#
#
#   c_scode <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(scode),
#                                      quiet = TRUE,
#                                      cpp_options = cpp_options,
#                                      stanc_options = stanc_options,
#                                      dir = NULL,
#                                      pedantic = FALSE,
#                                      include_paths = NULL,
#                                      user_header = NULL,
#                                      compile_model_methods = FALSE,
#                                      compile_hessian_method = FALSE,
#                                      compile_standalone = FALSE)
#
#
#   iter_sampling <- brm_args$iter - brm_args$warmup
#   iter_warmup   <- brm_args$warmup
#
#   cb_fit <- c_scode$sample(
#     data = sdata,
#     seed = brm_args$seed,
#     init = brm_args$init,
#     chains = brm_args$chains,
#     parallel_chains = brm_args$cores,
#     threads_per_chain = brm_args$threads$threads,
#     opencl_ids = brm_args$opencl,
#     iter_sampling = iter_sampling,
#     iter_warmup = iter_warmup,
#     thin = brm_args$thin,
#     max_treedepth = brm_args$control$max_treedepth,
#     adapt_delta = brm_args$control$adapt_delta,
#     adapt_engaged = TRUE,
#     fixed_param = FALSE,
#     show_messages = show_messages,
#     show_exceptions = show_exceptions
#   )
#
#   cb_fit <- rstan::read_stan_csv(cb_fit$output_files())
#   attributes(cb_fit)$CmdStanModel <- c_scode
#
#   brm_args_empty <- brm_args
#   brm_args_empty$empty <- TRUE
#
#   # Create an empty brms object -> Set empty = TRUE
#   bfit <- do.call(brms::brm, brm_args_empty)
#   bfit$fit = cb_fit
#   bfit <- brms::rename_pars(bfit)
#   bfit
# }



#' Fit model via rstan
#'
#' @param scode A character string of model code
#' @param sdata A list of data objects
#' @param brm_args A list of argument passes to the brm
#'
#' @return An object of class \code{bgmfit}
#' @keywords internal
#' @noRd
#'
brms_via_rstan <- function(scode, sdata, brm_args) {
  if(!is.null(brm_args$threads$threads)) {
    stan_threads <- TRUE
  } else {
    stan_threads <- FALSE
  }

  if(stan_threads) {
    rstan::rstan_options(threads_per_chain = brm_args$threads$threads)
  }

  # rstan::rstan_options(auto_write = TRUE)

  algorithm <- "NUTS" # c("NUTS", "HMC", "Fixed_param")

  cpp_options <- list(stan_threads = stan_threads)
  stanc_options <- NULL

  if(brm_args$silent == 0) {
    show_messages = TRUE
    show_exceptions = TRUE
  }
  if(brm_args$silent == 1) {
    show_messages = TRUE
    show_exceptions = FALSE
  }
  if(brm_args$silent == 2) {
    show_messages = FALSE
    show_exceptions = FALSE
  }

  message("Compiling Stan program...")
  c_scode <- rstan::stan_model(
    # file,
    model_name = "anon_model",
    model_code = scode,
    stanc_ret = NULL,
    boost_lib = NULL,
    eigen_lib = NULL,
    save_dso = TRUE,
    verbose = FALSE,
    auto_write = rstan::rstan_options("auto_write"),
    obfuscate_model_name = TRUE,
    allow_undefined = isTRUE(getOption("stanc.allow_undefined", FALSE)),
    allow_optimizations = isTRUE(getOption("stanc.allow_optimizations", FALSE)),
    standalone_functions=isTRUE(getOption("stanc.standalone_functions", FALSE)),
    use_opencl = isTRUE(getOption("stanc.use_opencl", FALSE)),
    warn_pedantic = isTRUE(getOption("stanc.warn_pedantic", FALSE)),
    warn_uninitialized = isTRUE(getOption("stanc.warn_uninitialized", FALSE)),
    includes = NULL,
    isystem = c(if (!missing(file)) dirname(file), getwd())
  )


  message("Start sampling")
  cb_fit <- rstan::sampling(
    object = c_scode,
    data = sdata,
    pars = NA,
    chains = brm_args$chains,
    iter = brm_args$iter,
    warmup = brm_args$warmup,
    thin = brm_args$thin,
    seed = brm_args$seed,
    init = brm_args$init,
    cores = brm_args$cores,
    check_data = TRUE,
    sample_file = NULL,
    diagnostic_file = NULL,
    verbose = FALSE,
    algorithm = algorithm,
    control = brm_args$control,
    include = TRUE,
    open_progress = interactive() && !isatty(stdout()) &&
      !identical(Sys.getenv("RSTUDIO"), "1"),
    show_messages = show_messages
  )

  # Create an empty brms object and populate it with the rsran fit
  brm_args$empty <- TRUE
  bfit      <- do.call(brms::brm, brm_args)
  bfit$fit  <- cb_fit
  bfit      <- brms::rename_pars(bfit)
  bfit
}



#' Transform initial values for parameters with lower bound
#'
#' @param x An initial value on unconstrained parameter space
#' (\code{real value}).
#'
#' @param lb A lower bound on the parameter space (\code{real value}).
#'
#' @return Transformed initial (\code{real value}).
#' @keywords internal
#' @noRd
#'
inits_lb <- function(x, lb = 0) {
  if(x < 1) 1+log(1+x) + lb else log(x) + lb
}


# Adapted from https://rdrr.io/cran/rempsyc/src/R/utils.R

#' @title Check and install package if not already installed
#' @param pkgs Packages to install if not already installed
#' @keywords internal
#' @noRd
#'
check_and_install_if_not_installed <- function(pkgs,
                                               getfun = NULL,
                                               installpkg = TRUE,
                                               verbose = FALSE) {
  successfully_loaded <- vapply(
    pkgs, requireNamespace,
    FUN.VALUE = logical(1L), quietly = TRUE
  )
  required_pkgs <- names(which(successfully_loaded == FALSE))


  if(!is.null(getfun)) {
    if(is.symbol(getfun)) getfun <- deparse(getfun)
    if(verbose) {
      message('Checking required packages for ', getfun, " ",
              "\n ",
              paste(pkgs, collapse = ", "))
    } # if(verbose) {
  } # if(!is.null(getfun)) {

  # Dont install package in function
  # CRAN does not accept it, so comment it out

  if(installpkg) {
    # message('Installing required packages',
    #         paste(required_pkgs, collapse = ", "))
    #
    # utils::install.packages(required_pkgs,
    #                         repos = "http://cran.us.r-project.org")
  }
}




#' Plot tripple logistic model with marked x and y axis
#'
#' @param model An object of class \code{brmsfit}
#' (\code{real value}).
#'
#' @param return_plot A logical (default \code{FALSE}) to indicate whether to
#' return the plot object.
#'
#' @param print_plot A logical (default \code{FALSE}) to indicate whether to
#' print plot along with the output
#'
#' @param digits A integer (default \code{2}) to set the number of decimal
#' places.
#'
#' @return A plot object if (\code{return_plot = TRUE}).
#' @keywords internal
#' @noRd
#'
plot_lositic3 <- function(model,
                          return_plot = FALSE,
                          print_plot = TRUE,
                          digits = 2 ,
                          resp = NULL,
                          envir = NULL,
                          ...) {

  if(is.null(envir)) {
    envir <- parent.frame()
  }

  if (is.null(resp)) {
    resp_ <- resp
  } else if (!is.null(resp)) {
    resp_ <- paste0(resp, "_")
  }

  args <- list(...)
  args$model <- model

  pob    <- do.call(plot_curves, args)
  fixed_ <- brms::fixef(model)

  xintercept_1 <- fixed_[3,1]
  xintercept_2 <- fixed_[6,1] + fixed_[3,1]
  xintercept_3 <- fixed_[9,1]

  Funx0 <- NULL;
  Funx1 <- NULL;

  assign(paste0(resp_,
                model$model_info[['namesexefuns']],
                '0'),
         model$model_info$exefuns[[paste0(resp_,
                                          model$model_info[['namesexefuns']],
                                          '0')]], envir = envir)

  assign('Funx0',
         model$model_info$exefuns[[paste0(resp_,
                                          model$model_info[['namesexefuns']],
                                          '0')]], envir = envir)

  assign('Funx1',
         model$model_info$exefuns[[paste0(resp_,
                                          model$model_info[['namesexefuns']],
                                          '1')]], envir = envir)

  assign('Funx2',
         model$model_info$exefuns[[paste0(resp_,
                                          model$model_info[['namesexefuns']],
                                          '2')]], envir = envir)

  # distance

  yintercept_1 <-
    Funx0(xintercept_1,
          fixed_[1,1], fixed_[2,1], fixed_[3,1],
          fixed_[4,1], fixed_[5,1], fixed_[6,1],
          fixed_[7,1], fixed_[8,1], fixed_[9,1])

  yintercept_2 <-
    Funx0(xintercept_2,
          fixed_[1,1], fixed_[2,1], fixed_[3,1],
          fixed_[4,1], fixed_[5,1], fixed_[6,1],
          fixed_[7,1], fixed_[8,1], fixed_[9,1])

  yintercept_3 <-
    Funx1(xintercept_3,
          fixed_[1,1], fixed_[2,1], fixed_[3,1],
          fixed_[4,1], fixed_[5,1], fixed_[6,1],
          fixed_[7,1], fixed_[8,1], fixed_[9,1])



  # velocity - for secondry axis
  getfb <- transform.sec.axis(pob$data$Estimate.x, pob$data$Estimate.y)

  xyvelocity_1 <-
    Funx1(xintercept_1,
          fixed_[1,1], fixed_[2,1], fixed_[3,1],
          fixed_[4,1], fixed_[5,1], fixed_[6,1],
          fixed_[7,1], fixed_[8,1], fixed_[9,1])

  yintercept_v1 <- getfb$fwd(xyvelocity_1)

  xyvelocity_2 <-
    Funx1(xintercept_2,
          fixed_[1,1], fixed_[2,1], fixed_[3,1],
          fixed_[4,1], fixed_[5,1], fixed_[6,1],
          fixed_[7,1], fixed_[8,1], fixed_[9,1])

  yintercept_v2 <- getfb$fwd(xyvelocity_2)


  xyvelocity_3 <-
    Funx1(xintercept_3,
          fixed_[1,1], fixed_[2,1], fixed_[3,1],
          fixed_[4,1], fixed_[5,1], fixed_[6,1],
          fixed_[7,1], fixed_[8,1], fixed_[9,1])

  yintercept_v3 <- getfb$fwd(xyvelocity_3)


  xintercept_1 <- round(xintercept_1, digits)
  yintercept_1 <- round(yintercept_1, digits)
  xyvelocity_1 <- round(xyvelocity_1, digits)

  xintercept_2 <- round(xintercept_2, digits)
  yintercept_2 <- round(yintercept_2, digits)
  xyvelocity_2 <- round(xyvelocity_2, digits)

  xintercept_3 <- round(xintercept_3, digits)
  yintercept_3 <- round(yintercept_3, digits)
  xyvelocity_3 <- round(xyvelocity_3, digits)



  setprint_1 <-
    paste0("stage 1: ", "\n ",
           "timing = ", xintercept_1, "; velocit = ",
           xyvelocity_1, "; size = ", yintercept_1)

  setprint_2 <-
    paste0("stage 2: ", "\n ",
           "timing = ", xintercept_2, "; velocit = ",
           xyvelocity_2, "; size = ", yintercept_2)

  setprint_3 <-
    paste0("stage 3: ", "\n ",
           "timing = ", xintercept_3, "; velocit = ",
           xyvelocity_3, "; size = ", yintercept_3)
  #


  pob <- pob +
    ggplot2::geom_hline(ggplot2::aes(yintercept = yintercept_1))  +
    ggplot2::geom_hline(ggplot2::aes(yintercept = yintercept_2) ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = yintercept_3))  +

    ggplot2::geom_vline(ggplot2::aes(xintercept = xintercept_1) ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = xintercept_2))  +
    ggplot2::geom_vline(ggplot2::aes(xintercept = xintercept_3) ) +

    ggplot2::geom_hline(ggplot2::aes(yintercept =  yintercept_v1 ) ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept =  yintercept_v2 ) ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept =  yintercept_v3 ) )


  if(print_plot) {
    print(pob)
    setprint_123 <- paste(setprint_1, setprint_2, setprint_3, sep = "\n")
    cat(setprint_123)
  }

  if(return_plot) return(pob)
}



# Adapted from
# https://stackoverflow.com/questions/37149649/randomly-sample-groups

#' @title select a random sample of n groups
#' @param data A data frame
#' @param size The number of groups to be selected
#' @examples
#' # example code
#'  set.seed(1234)
#'  subdata <- berkeley_mdata %>% sample_n_of_groups(size = 2, id)
#' @keywords internal
#' @noRd
#'
sample_n_of_groups <- function(data, size, ...) {
  dots <- rlang::quos(...)

  group_ids <- data %>%
    dplyr::group_by(!!! dots) %>%
    dplyr::group_indices()

  sampled_groups <- sample(unique(group_ids), size)

  data %>%
    dplyr::filter(group_ids %in% sampled_groups) %>%
    droplevels()
}


#' An internal function to check the minimum version of the package
#'
#' @param pkg A character string of package names
#' @param minver A character string of minimum version of the package
#' @param verbose A logical (default \code{FALSE}) to check 
#' @param ... other arguments. Currently ignored.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'

check_pkg_version_exists <- function(pkg, 
                                     minversion = NULL, 
                                     verbose = FALSE,
                                     ...) {
  
  try(zz <- insight::check_if_installed(pkg, 
                                        minimum_version = minversion,
                                        ...))
 
  if(!isTRUE(zz)) {
    if(verbose) {
      message("Please install the latest version of the 'brms' package",
              "\n ",
              "remotes::install_github('paul-buerkner/brms')")
    }
  }
  return(zz)
}


#' An internal function to check for the exposed function
#'
#' @param o An object used as an index for functions
#' @param checks A logical (default \code{FALSE}) to check if funnctions are 
#' attached to the \code{model}.
#' @inherit growthparameters.bgmfit params
#' @param ... other arguments. Currently ignored.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'
check_if_functions_exists <- function(model, 
                                      o = NULL, 
                                      xcall = NULL, 
                                      verbose = TRUE, 
                                      usesavedfuns = FALSE, 
                                      checks = FALSE,...) {
  
  if(!checks) {
    if(is.null(o)) stop("object 'o' must be specified")
  }
  
  check_brms_v <- 
  check_pkg_version_exists('brms', 
                           minversion = get_package_minversion('brms'), 
                           prompt = FALSE,
                           stop = FALSE,
                           verbose = FALSE)
  
  latest_brms_v <- TRUE
  if(!isTRUE(check_brms_v)) {
    latest_brms_v <- FALSE
  }
  
  # if(exists(o[[1]], mode = "function", envir = globalenv())) {
  #   envgtf <- TRUE
  # } else {
  #   envgtf <- FALSE
  # }
  
  if(is.null(xcall)) {
    xcall <- strsplit( deparse(sys.calls()[[sys.nframe()-1]]) , "\\(")[[1]][1]
  }
  
  classname <- attr(model, 'class')[2]
  calname.fun <- xcall # match.call()[1]
  calname.fun <- gsub(paste0(".", classname), "", calname.fun)
  
  msg1 <- paste0(" Please expose user defined Stan function before calling the",
                 "\n ",
                 "'", calname.fun, "()'", " function",
                  "\n ",
                 " (See '?expose_model_functions()' for details).",
                 "\n ",
                 "\n ",
                 "Note that if you have already exposed Stan functions in ",
                 "'bsitar()' call,\n then those saved functions can be used here ",
                 "by setting usesavedfuns = TRUE",
                 "\n ",
                 paste0(calname.fun,
                        "(...,", " usesavedfuns = TRUE"),
                 "\n "              )
  
  
  msg2 <- paste0("Please expose user defined Stan function before calling the",
                 "\n",
                 "'", calname.fun, "()'", " function",
                 # "\n ",
                 " (See '?expose_model_functions()' for details).",
                 "\n ",
                 "\n ",
                 "Note that you can use 'usesavedfuns = TRUE' only if Stan ",
                 "functions have been ",
                 "\n",
                 " exposed and saved within the 'bsitar()' ",
                 "by there using 'expose_functions = TRUE'",
                 "\n "              )
  
  msg3 <- paste0(" Please expose user defined Stan function before calling the ",
                 "'", calname.fun, "()'", " function",
                 "\n ",
                 "(See '?expose_model_functions()' for details).",
                 "\n ",
                 "\n ",
                 "Also, 'envir' should be set as global environment i.e.,",
                 "\n ",
                 paste0(calname.fun, "(...,", " envir = "," .GlobalEnv)"),
                 "\n ",
                 "This is a known issue ",
                 "(https://github.com/paul-buerkner/brms/issues/1577)",
                 "\n ",
                 "\n ",
                 "Note that if you have already exposed the Stan functions in ",
                 "'bsitar()' call,\n then those saved functions can be used ",
                 "here by setting 'usesavedfuns = TRUE'",
                 "\n ",
                 paste0(calname.fun,
                        "(...,", " usesavedfuns = TRUE, envir = "," .GlobalEnv)"),
                 "\n "              )
  
  if(!latest_brms_v) {
    msg3 <- paste0(msg3, 
                   "\n ",
                   "Or else, you can install the lates deveopmental versions ",
                   " of the brms package:",
                   "\n ",
                   "remotes::install_github('paul-buerkner/brms')"
    )
  }
  
  if(checks) {
    if(is.null(model$model_info$exefuns[[1]])) {
      if(!is.null(usesavedfuns)) {
        if(!usesavedfuns & latest_brms_v) message(msg1)
        if(usesavedfuns & latest_brms_v) message(msg2)
        
        if(!usesavedfuns & !latest_brms_v) message(msg3)
        if(usesavedfuns & !latest_brms_v) message(msg3)
        
      }
    }
    return(invisible(NULL))
  }
  
  
  if(exists(o[[1]], mode = "function", envir = globalenv())) {
    envgtf <- TRUE
  } else {
    envgtf <- FALSE
  }

  if(verbose) {
    if(!envgtf) {
      if(verbose) message(msg3)
    }
  }

  if(!envgtf) {
    en <- NULL
  } else if(envgtf) {
    en <- environment(eval(parse(text = o[[1]])))
  }

  return(en)
}




#' An internal function to check required package(s) installed 
#'
#' @param o An object used as an index for functions
#' @param checks A logical (default \code{FALSE}) to check if funnctions are 
#' attached to the \code{model}.
#' @inherit growthparameters.bgmfit params
#' @param ... other arguments. Currently ignored.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'
check_if_package_installed <- function(model, 
                                      xcall = NULL, 
                                      package = NULL, 
                                      reason = "for this function to work",
                                      stop = TRUE,
                                      minimum_version = NULL,
                                      quietly = FALSE,
                                      prompt = FALSE,
                                      verbose = TRUE, 
                                      ...) {
  
  
  
  if(is.null(xcall)) {
    xcall <- strsplit( deparse(sys.calls()[[sys.nframe()-1]]) , "\\(")[[1]][1]
  }
  
  classname <- attr(model, 'class')[2]
  calname.fun <- xcall # match.call()[1]
  calname.fun <- gsub(paste0(".", classname), "", calname.fun)
  
  
  if(is.null(package)) {
    if(calname.fun == "plot_curves") {
      package <- c('ggplot2', 'jtools')
    } else if(calname.fun == "growthparameters_comparison") {
      package <- c('tidyr', 'collapse')
    } else if(calname.fun == "marginal_draws") {
      package <- c('tidyr', 'collapse')
    } else  {
      return(invisible(NULL))
    }
  } # if(is.null(package)) {
  
  
  
  if(!is.null(package)) {
    if(is.null(reason)) {
      reason <- paste0("for ", "'", calname.fun, "()'", " function", " to work")
    } else {
      reason <- reason
    }
    
    insight::check_if_installed(package = package,
                                reason = reason,
                                stop = stop,
                                minimum_version = minimum_version,
                                quietly = quietly,
                                prompt = prompt
    )
    
    return(invisible(NULL))
  } # if(!is.null(package)) {
  
}







#' An internal function to get the environment of an object
#'
#' @param x A symbol or a character string.
#' @param geteval A logical (default \code{TRUE}) to indicate whether to return
#' the object as a character string or as an environment.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'

getEnv <- function(x, geteval = TRUE) {
  if(!is.character(x)) xobj <- deparse(substitute(x)) else xobj <- x
  gobjects <- ls(envir=.GlobalEnv)
  envirs <- gobjects[sapply(gobjects, function(x) is.environment(get(x)))]
  envirs <- c('.GlobalEnv', envirs)
  xin <- sapply(envirs, function(e) xobj %in% ls(envir=get(e)))
  out <- envirs[xin]
  if(geteval) out <- eval(parse(text = out))
  out
}



#' An internal function to get the 'model' name from the arguments
#'
#' @param arguments A list of arguments.
#' @param asstr A logical (default \code{FALSE}) to indicate whether to
#' return the object as a character string.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'

getpipedot <- function(arguments, asstr = FALSE) {
  if(deparse(arguments$model) == ".") {
    first_call <- sys.calls()[[1]] # get the first entry on the call stack
    lhs <- first_call[[2]] # get the second element of this entry
    mymodel <- lhs # rlang::as_name(lhs) # lhs
  } else {
    mymodel <- arguments$model
  }
  if(asstr) mymodel <- deparse(mymodel)
  mymodel
}



#' An internal function to check an argument
#'
#' @param checkarg A list of defined arguments.
#' @param checkcall A list of passed arguments.
#' @param check A aymbol or a character string to be checked.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'

checkifargmiss <- function(checkarg, checkcall, check) {
  defined <- checkarg
  passed <- names(as.list(checkcall)[-1])
  if(!is.character(check)) check <- deparse(substitute(check))
  allargs <- unique(c(defined, passed))
  if (!check %in% allargs) checkresulst <- FALSE else checkresulst <- TRUE
  checkresulst
}



#' An internal function to convert dummy variables to a factor variable
#'
#' @param df A data frame.
#' @param factor.dummy A vector of character strings that will be converted to a
#'   factor variable.
#' @param factor.name A character string to name the newly created factor
#'   variable.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'

convert_dummy_to_factor <- function(df,
                                    factor.dummy = NULL,
                                    factor.level = NULL,
                                    factor.name = NULL) {

  if(!is.data.frame(df)) stop("df should be a data frame")
  all.dfnames  <- colnames(df)
  # factor.dummy <- c("classClassI",  "classClassII" )
  if(is.null(factor.dummy)) {
    dfout <- df
  } else if(!is.null(factor.dummy)) {
    if(!is.null(factor.name)) {
      if(!is.character(factor.name)) {
        stop("factor.name should be a character string")
      }
    } else if(is.null(factor.name)) {
      factor.variable <- "factor.var"
    }
    all_min_factor <- setdiff(all.dfnames, factor.dummy)
    dfout <- cbind(df[, all_min_factor],
                   tempvar = factor(max.col(df[, factor.dummy]),
                                        ordered = TRUE))
    colnames(dfout) <- gsub('tempvar', factor.variable, names(dfout))
    if(!is.null(factor.level)) {
      if(length(factor.dummy) != length(factor.level)) {
        stop("Lengths of factor.dummy and factor.level must be same")
      }
      levels(dfout[[factor.variable]]) <- factor.level
    } else if(is.null(factor.level)) {
      levels(dfout[[factor.variable]]) <- factor.dummy
    }
    dfout <- cbind(dfout, df[, factor.dummy])
  }
  dfout
}




#' An internal function to add growthparameters to the plot_curves data
#'
#' @param data A data frame returned by the \code{plot_curves} function.
#' @param gpdata A data frame with growth parameters. If \code{NULL} (default),
#'   the \code{gpdata} is taken from the data returned by the
#'   \code{plot_curves()} as an \code{attribute}.
#' @param Parametername A character string specifying the name of the Parameter
#'   column in the \code{gpdata}.
#' @param parmcols A character string, or a vector of character strings
#'   specifying the name of growth parameter estimates columns in the
#'   \code{gpdata}. Typically, they are \code{Estimate} and the associated
#'   uncertainty parameters such as \code{Est.Error}, \code{Q2.5}, and
#'   \code{Q97.5}.
#' @param nonparmcols A character string, or a vector of character strings
#'   specifying the name of columns in the \code{gpdata} other than the names
#'   specified as \code{parmcols}.
#' @param byjoincols A character string, or a vector of character strings
#'   specifying the name of columns to be used in joining the \code{data} and
#'   \code{gpdata}. Typically, they are same as \code{nonparmcols}.
#' @param ... Other internal arguments passed to the
#'   \code{add_parms_to_curve_data} function.
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
add_parms_to_curve_data <- function(data,
                                    gpdata = NULL,
                                    Parametername = NULL,
                                    parmcols = NULL,
                                    nonparmcols = NULL,
                                    byjoincols = NULL,
                                    ...) {
  tojoinwith <- data
  if(  is.null(gpdata)) gp <- attr(tojoinwith, "growthparameters")
  if(! is.null(gpdata)) gp <- gpdata

  # Initiate non formalArgs()
  . <- NULL;

  if(is.null(Parametername)) {
    Parametername <- "Parameter"
  }
  if(is.null(parmcols)) {
    parmcols <- c('Estimate', "Est.Error", "Q2.5", "Q97.5")
  }
  if(is.null(nonparmcols))  {
    stop("Please specify the 'nonparmcols'")
  }
  if(is.null(byjoincols)) {
    stop("Please specify the 'byjoincols'")
  }

  parmnames <- gp %>% dplyr::select(dplyr::all_of(Parametername)) %>%
    unique() %>% unlist() %>% as.vector()

  whati_list <- list()
  for (whati in parmnames) {
    addpre <- paste0(whati, ".")
    addsuf <- NULL # paste0(".", whati)

    tojoinit2 <-
      gp %>% dplyr::filter(!!dplyr::sym(Parametername) == whati) %>%
      dplyr::select(dplyr::any_of(parmcols)) %>%
      stats::setNames(paste0(addpre, names(.), addsuf))

    whati_list[[whati]] <- tojoinit2
  }

  tojoinit1 <-
    gp %>%
    dplyr::filter(!!dplyr::sym(Parametername) == names(whati_list)[1]) %>%
    dplyr::select(dplyr::any_of(nonparmcols))

  # Note dplyr::bind_cols instead of cbind. cbind adds again list name as prefix
  tojoinit2all <- whati_list %>% do.call(dplyr::bind_cols, .) %>% data.frame()

  tojoinit12 <- cbind(tojoinit1, tojoinit2all)

  mergebycols <- intersect(nonparmcols, byjoincols)
  setdiffcols <- setdiff(byjoincols, nonparmcols)
  if(length(setdiffcols) != 0) {
    stop("Variable(s) ", "'", paste(setdiffcols, collapse = ", "), "'",
         " missing in nonparmcols" )
  }

  tojoinwith <- tojoinwith %>% dplyr::left_join(., tojoinit12, by = byjoincols)
  return(tojoinwith)
}



#' An internal function to bind rows of unequal lengths (adapted from qpcR:::cbind.na)
#'
#' @param deparse.level An integer to set deparse level. 
#' @param ... A list or name of column vectors.
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
rbind_fill_na1 <- function (..., deparse.level = 1) {
    na <- nargs() - (!missing(deparse.level))
    deparse.level <- as.integer(deparse.level)
    stopifnot(0 <= deparse.level, deparse.level <= 2)
    argl <- list(...)
    while (na > 0 && is.null(argl[[na]])) {
      argl <- argl[-na]
      na <- na - 1
    }
    if (na == 0) 
      return(NULL)
    if (na == 1) {
      if (isS4(..1)) 
        return(methods::rbind2(..1))
      else return(matrix(..., nrow = 1))
    }
    if (deparse.level) {
      symarg <- as.list(sys.call()[-1L])[1L:na]
      Nms <- function(i) {
        if (is.null(r <- names(symarg[i])) || r == "") {
          if (is.symbol(r <- symarg[[i]]) || deparse.level == 
              2) 
            deparse(r)
        }
        else r
      }
    }
    if (na == 0) {
      r <- argl[[2]]
      fix.na <- FALSE
    }
    else {
      nrs <- unname(lapply(argl, ncol))
      iV <- sapply(nrs, is.null)
      fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
      if (deparse.level) {
        if (fix.na) 
          fix.na <- !is.null(Nna <- Nms(na))
        if (!is.null(nmi <- names(argl))) 
          iV <- iV & (nmi == "")
        ii <- if (fix.na) 
          2:(na - 1)
        else 2:na
        if (any(iV[ii])) {
          for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i))) 
            names(argl)[i] <- nmi
        }
      }
      nCol <- as.numeric(sapply(argl, function(x) if (is.null(ncol(x))) length(x) else ncol(x)))
      maxCol <- max(nCol, na.rm = TRUE)
      argl <- lapply(argl, function(x) if (is.null(ncol(x))) 
        c(x, rep(NA, maxCol - length(x)))
        else cbind(x, matrix(, nrow(x), maxCol - ncol(x))))
      namesVEC <- rep(NA, maxCol)
      for (i in 1:length(argl)) {
        CN <- colnames(argl[[i]])
        m <- !(CN %in% namesVEC)
        namesVEC[m] <- CN[m]
      }
      for (j in 1:length(argl)) {
        if (!is.null(ncol(argl[[j]]))) 
          colnames(argl[[j]]) <- namesVEC
      }
      r <- do.call(rbind, c(argl[-1L], list(deparse.level = deparse.level)))
    }
    d2 <- dim(r)
    colnames(r) <- colnames(argl[[1]])
    r <- methods::rbind2(argl[[1]], r)
    if (deparse.level == 0) 
      return(r)
    ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
    ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
    if (ism1 && ism2) 
      return(r)
    Nrow <- function(x) {
      d <- dim(x)
      if (length(d) == 2L) 
        d[1L]
      else as.integer(length(x) > 0L)
    }
    nn1 <- !is.null(N1 <- if ((l1 <- Nrow(..1)) && !ism1) Nms(1))
    nn2 <- !is.null(N2 <- if (na == 2 && Nrow(..2) && !ism2) Nms(2))
    if (nn1 || nn2 || fix.na) {
      if (is.null(rownames(r))) 
        rownames(r) <- rep.int("", nrow(r))
      setN <- function(i, nams) rownames(r)[i] <<- if (is.null(nams)) 
        ""
      else nams
      if (nn1) 
        setN(1, N1)
      if (nn2) 
        setN(1 + l1, N2)
      if (fix.na) 
        setN(nrow(r), Nna)
    }
    r
  }


#' An internal function to bind columns of unequal lengths (adapted from qpcR:::cbind.na)
#'
#' @param deparse.level An integer to set deparse level. 
#' @param ... A list or name of column vectors.
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
cbind_fill_na1 <- function (..., deparse.level = 1) {
  na <- nargs() - (!missing(deparse.level))
  deparse.level <- as.integer(deparse.level)
  stopifnot(0 <= deparse.level, deparse.level <= 2)
  argl <- list(...)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na == 0) 
    return(NULL)
  if (na == 1) {
    if (isS4(..1)) 
      return(methods::cbind2(..1))
    else return(matrix(...))
  }
  if (deparse.level) {
    symarg <- as.list(sys.call()[-1L])[1L:na]
    Nms <- function(i) {
      if (is.null(r <- names(symarg[i])) || r == "") {
        if (is.symbol(r <- symarg[[i]]) || deparse.level == 
            2) 
          deparse(r)
      }
      else r
    }
  }
  if (na == 0) {
    r <- argl[[2]]
    fix.na <- FALSE
  }
  else {
    nrs <- unname(lapply(argl, nrow))
    iV <- sapply(nrs, is.null)
    fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
    if (deparse.level) {
      if (fix.na) 
        fix.na <- !is.null(Nna <- Nms(na))
      if (!is.null(nmi <- names(argl))) 
        iV <- iV & (nmi == "")
      ii <- if (fix.na) 
        2:(na - 1)
      else 2:na
      if (any(iV[ii])) {
        for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i))) 
          names(argl)[i] <- nmi
      }
    }
    nRow <- as.numeric(sapply(argl, function(x) NROW(x)))
    maxRow <- max(nRow, na.rm = TRUE)
    argl <- lapply(argl, function(x) if (is.null(nrow(x))) 
      c(x, rep(NA, maxRow - length(x)))
      else rbind_fill_na1(x, matrix(, maxRow - nrow(x), ncol(x))))
    r <- do.call(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
  }
  d2 <- dim(r)
  r <- methods::cbind2(argl[[1]], r)
  if (deparse.level == 0) 
    return(r)
  ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
  ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
  if (ism1 && ism2) 
    return(r)
  Ncol <- function(x) {
    d <- dim(x)
    if (length(d) == 2L) 
      d[2L]
    else as.integer(length(x) > 0L)
  }
  nn1 <- !is.null(N1 <- if ((l1 <- Ncol(..1)) && !ism1) Nms(1))
  nn2 <- !is.null(N2 <- if (na == 2 && Ncol(..2) && !ism2) Nms(2))
  if (nn1 || nn2 || fix.na) {
    if (is.null(colnames(r))) 
      colnames(r) <- rep.int("", ncol(r))
    setN <- function(i, nams) colnames(r)[i] <<- if (is.null(nams)) 
      ""
    else nams
    if (nn1) 
      setN(1, N1)
    if (nn2) 
      setN(1 + l1, N2)
    if (fix.na) 
      setN(ncol(r), Nna)
  }
  r
}





#' An internal function to bind columns of unequal lengths
#'
#' @param names A vector of character string to name columns. 
#' @param ... A list or names of column vectors.
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
cbind_fill_na2 <- function(..., names = NA) {
  xlist = list(...)
  cbindfill.id <- NULL;
  suppressWarnings({
    y= Reduce(
      function(a,b) {
        if(is.vector(a)) na = length(a)
        if(is.factor(a)) na = levels(a)
        if(is.factor(a)) a = data.frame(a) %>% droplevels()
        if(is.data.frame(a)|is.matrix(a)) na = nrow(a)
        if(is.vector(b)) nb = length(b)
        if(is.factor(b)) nb = levels(b)
        if(is.factor(b)) b = data.frame(b) %>% droplevels()
        if(is.data.frame(b)|is.matrix(b)) nb = nrow(b)
        subset(
          merge(
            cbind(cbindfill.id = 1:na, a),
            cbind(cbindfill.id = 1:nb, b),
            all = TRUE, by = "cbindfill.id"
          ),
          select = -cbindfill.id
        )}
      ,xlist)
    if(is.na(names)) colnames(y) <- names(xlist) else colnames(y) <- names
  })
  return(y)
}



#' An internal function to check and appropriately set brms exported functions
#'
#' @param call A \code{call} object, typically the \code{match.call()}.
#' @param arg A character string or vector of character string of brms exported
#'   functions.
#' @param prefix A character string specifying the namespace i.e, \code{brms::}
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
check_brms_args <- function(call, arg, prefix = NULL) {
  newcall <- call
  if(is.null(prefix)) prefix <- "brms::"
  for (argi in arg) {
    if(!is.null((newcall[[argi]]))) {
      argin <- newcall[[argi]]
      argin <- deparse(substitute(argin))
      if(!grepl(prefix, argin)) {
        newargin <- paste0(argi, " = ", prefix, argin)
        newcall[[argi]] <- NULL
        newcall[[argi]] <- (str2expression(newargin))
      } else {
        newcall[[argi]] <- newcall[[argi]]
      } 
    }
    if(is.null((newcall[[argi]]))) {
      newcall <- newcall
    }
  }
  return(newcall)
}



#' An internal function to check and appropriately set brms exported functions
#' 
#' @details This function is used when no random effects are included and the
#'   groupvar is NULL An an artificial groupvar created to maintain consistency
#'   across various functions.
#' 
#' @param model model An object of class \code{bgmfit}.
#' @param newdata A data frame
#' @param IDvar A character string specifying the group identifier
#' @param resp A character string specifying the response variable (default
#'   \code{NULL})
#' @param verbose A logical to indicate whether to print relevant information.
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
check_newdata_args <- function(model, newdata, IDvar, resp = NULL, verbose = FALSE) {
  # This is when no random effects and this groupvar is NULL
  # Therefore, an artificial group var created
  # see also changes made to the get_idata function lines 17
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  
  if (is.null(model$model_info$groupvar)) {
    name_hypothetical_id <- paste0("id", resp_rev_)
    model$model_info$groupvar <- name_hypothetical_id
    newdata[[name_hypothetical_id]] <- as.factor("tempid")
  } else if (!is.null(model$model_info$groupvar)) {
    if(length(newdata[[model$model_info$groupvar]]) == 0) {
      # name_hypothetical_id <- paste0("hy_id", resp_rev_)
      if(length(IDvar) > 1) {
        name_hypothetical_id <- IDvar[1] 
      } else {
        name_hypothetical_id <- IDvar
      }
      model$model_info$groupvar <- name_hypothetical_id
      newdata[[name_hypothetical_id]] <- as.factor("tempid")
    }
  }
  
  newdata
}





#' An internal function to create interactions within the dplyr framework
#'
#' @param data A data frame.
#' @param vars A character vector specifying the variables included in
#'   interaction.
#' @param varname A character that will be used as a name for the interaction
#'   term created.
#' @param envir An environment for function evaluation.
#' @param full A logical to indicate whether to return the full data frame.
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
vars_to_interaction <- function(data, 
                                vars, 
                                varname, 
                                envir = NULL, 
                                full = FALSE) {
  if(is.null(envir)) envir <- parent.frame()
  `:=` <- NULL;
  data_in <- data
  nested_vars_x <- paste0("interaction(", paste(vars, collapse = ","), ")" )
  data_ou <- data %>% dplyr::mutate(!! varname := eval(parse(text = nested_vars_x),
                                                       envir = envir)) %>%
    dplyr::select(dplyr::all_of(varname)) %>% 
    dplyr::select(dplyr::all_of(varname)) %>% unlist() # %>% as.vector()
  attr(data_ou, "names") <- NULL
  if(full) {
    data_in[[varname]] <- data_ou
    return(data_in)
  } else {
    return(data_ou)
  }
}



#' An internal function to redefine grid with nested variables
#'
#' @param fullgrid A data frame.
#' @param fulldata A data frame.
#' @param all_vars A character vector specifying the variables included in
#'   interaction.
#' @param nested_vars A character vector specifying the variables included in
#'   interaction.
#' @param xvar A character
#' @param yvar A character
#' @param idvar A character
#' @param envir An environment for function evaluation.
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
refine_grid <- function(fullgrid = NULL, 
                        fulldata = NULL, 
                        varsvector = NULL, 
                        all_vars = NULL, 
                        nested_vars = NULL, 
                        xvar = NULL, 
                        yvar = NULL, 
                        idvar = NULL, 
                        envir = NULL) {
  if(is.null(fullgrid)) stop("Please specify fullgrid")
  if(is.null(all_vars)) stop("Please specify all_vars")
  if(is.null(nested_vars)) stop("Please specify nested_vars")
  if(is.null(fulldata) & is.null(varsvector)) 
    stop("Please specify at least fulldata or varsvector")
  if(!is.null(fulldata) & !is.null(varsvector)) 
    stop("Please specify either fulldata or varsvector, not both")
  
  if(!is.null(varsvector)) {
    # if(!is.vector(varsvector)) stop("varsvector must be a vector")
    if(!is.factor(varsvector)) stop("varsvector must be a factor vector")
  }
  
  `.` <- NULL;
  `:=` <- NULL;
  zzz <- NULL;
  nested_vars_name <- 'varname'
  
  if(!is.null(fulldata)) {
    zz <- fulldata %>% dplyr::arrange(!! as.name(all_vars)) %>% droplevels() %>% 
      dplyr::mutate(nested_vars_name = 
                      vars_to_interaction(., nested_vars, nested_vars_name)) %>% 
      dplyr::select(nested_vars_name) %>% unlist() %>% as.vector()
  }
  
  if(!is.null(varsvector)) {
    # envir <- parent.frame()
    # nested_vars_x <- paste0("interaction(", paste(varsvector, collapse = ","), 
    #                         ")")
    # zz <-  eval(parse(text = nested_vars_x), envir = envir)
    zz <- varsvector
  }
  
  zz2 <- fullgrid %>% dplyr::arrange(!! as.name(all_vars)) %>% droplevels() %>% 
    dplyr::mutate(nested_vars_name = 
                    vars_to_interaction(., nested_vars, nested_vars_name)) %>% 
    dplyr::select(nested_vars_name) %>% unlist() %>% as.vector()
  
  zzz3 <- intersect(zz, zz2)
  
  nested_vars_name <- 'zzz'
  out <- fullgrid %>% 
    dplyr::mutate(zzz = 
                    vars_to_interaction(., nested_vars, nested_vars_name)) %>%
    dplyr::filter(zzz %in% zzz3) %>% 
    dplyr::select(-dplyr::all_of('zzz')) %>% 
    dplyr::arrange(!! as.name(all_vars)) %>% droplevels()
  
  out
}




# https://stackoverflow.com/questions/71339547/how-to-add-a-label-to-the-x-y-
# axis-whenever-a-vertical-horizontal-line-is-ad

#' An internal function to extract xintercept label
#'
#' @param plot A \code{ggplot} object
#' @param xval A numeric value
#' @param linewidth Argument to control the width of line drawn
#' @param linetype Argument to control the type of line drawn
#' @param alpha Argument to control the transparency of line drawn
#' @param color_line Argument to control the color of line drawn
#' @param color_text Argument to control the color of value marked on axis
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
mark_value_on_xaxis <- function(plot, xval,
                                linewidth = 1, 
                                linetype = 1, 
                                alpha = 0.7,
                                color_line = 'black', 
                                color_text = 'black'
                                ) {
  
  try(insight::check_if_installed(c("ggplot2", "ggtext"), stop = FALSE, 
                                  prompt = FALSE))
  
  p2 <- ggplot2::ggplot_build(plot)
  breaks <- p2$layout$panel_params[[1]]$x$breaks
  breaks <- breaks[!is.na(breaks)]
  
  color <- c(color_text, rep("black", length(breaks)  ))
  setx <- (c(xval, breaks)) # sort
  labs <- as.character(setx)
  name <- glue::glue("<i style='color:{color}'>{labs}")
  
  plot +
    ggplot2::geom_vline(xintercept = xval, 
                        linewidth = linewidth,
                        linetype = linetype,
                        color = color_line,
                        alpha = alpha) +
    ggplot2::scale_x_continuous(breaks = setx, labels = name) +
    ggplot2::theme(axis.text.x = ggtext::element_markdown())
}

#' An internal function to extract xintercept label
#'
#' @param plot A \code{ggplot} object
#' @param yval A numeric value
#' @param linewidth Argument to control the width of line drawn
#' @param linetype Argument to control the type of line drawn
#' @param alpha Argument to control the transparency of line drawn
#' @param color_line Argument to control the color of line drawn
#' @param color_text Argument to control the color of value marked on axis
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
mark_value_on_yaxis <- function(plot, yval, 
                                linewidth = 1, 
                                linetype = 1, 
                                alpha = 0.7,
                                color_line = 'black', 
                                color_text = 'black'
                                ) {
  
  try(insight::check_if_installed(c("ggplot2", "ggtext"), stop = FALSE, 
                                  prompt = FALSE))
  
  p2 <- ggplot2::ggplot_build(plot)
  breaks <- p2$layout$panel_params[[1]]$y$breaks
  breaks <- breaks[!is.na(breaks)]
  
  color <- c(color_text, rep("black", length(breaks)  ))
  setx <- (c(yval, breaks)) 
  labs <- as.character(setx)
  name <- glue::glue("<i style='color:{color}'>{labs}")
  
  plot +
    ggplot2::geom_hline(yintercept = yval, 
                        linewidth = linewidth,
                        linetype = linetype,
                        color = color_line,
                        alpha = alpha
    ) +
    ggplot2::scale_y_continuous(breaks = setx, labels = name) +
    ggplot2::theme(axis.text.y = ggtext::element_markdown())
}



#' An internal function to extract and add xintercept value label
#'
#' @param p A \code{ggplot} object
#' @param linewidth Argument to control the width of line drawn
#' @param linetype Argument to control the type of line drawn
#' @param alpha Argument to control the transparency of line drawn
#' @param color_line Argument to control the color of line drawn
#' @param color_text Argument to control the color of value marked on axis
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
mark_value_of_xintercept <- function(plot,
                                     linewidth = 1, 
                                     linetype = 1, 
                                     alpha = 0.7,
                                     color_line = 'black', 
                                     color_text = 'black'
                                     ) {
  
  try(insight::check_if_installed(c("ggplot2", "ggtext"), stop = FALSE, 
                                  prompt = FALSE))
  
  p <- plot
  p2 <- ggplot2::ggplot_build(p)
  breaks <- p2$layout$panel_params[[1]]$x$breaks
  breaks <- breaks[!is.na(breaks)]
  
  vals <- unlist(lapply(seq_along(p$layers), function(x) {
    d <- ggplot2::layer_data(p, x)
    if('xintercept' %in% names(d)) d$xintercept else numeric()
  }))
  
  xval <- vals
  
  color <- c(color_text, rep("black", length(breaks)  ))
  setx <- c(xval, breaks)
  labs <- as.character(setx)
  name <- glue::glue("<i style='color:{color}'>{labs}")
  
  plot +
    ggplot2::geom_vline(xintercept = xval, 
                        linewidth = linewidth,
                        linetype = linetype,
                        color = color_line,
                        alpha = alpha) +
    ggplot2::scale_x_continuous(breaks = setx, labels = name) +
    ggplot2::theme(axis.text.x = ggtext::element_markdown())
}



#' An internal function to extract and add xintercept value label
#'
#' @param p A \code{ggplot} object
#' @param linewidth Argument to control the width of line drawn
#' @param linetype Argument to control the type of line drawn
#' @param alpha Argument to control the transparency of line drawn
#' @param color_line Argument to control the color of line drawn
#' @param color_text Argument to control the color of value marked on axis
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
mark_value_of_yintercept <- function(plot,
                                     linewidth = 1, 
                                     linetype = 1, 
                                     alpha = 0.7,
                                     color_line = 'black', 
                                     color_text = 'black'
                                     ) {
  
  try(insight::check_if_installed(c("ggplot2", "ggtext"), stop = FALSE, 
                                  prompt = FALSE))
  
  p <- plot
  p2 <- ggplot2::ggplot_build(p)
  breaks <- p2$layout$panel_params[[1]]$y$breaks
  breaks <- breaks[!is.na(breaks)]
  
  vals <- unlist(lapply(seq_along(p$layers), function(x) {
    d <- ggplot2::layer_data(p, x)
    if('yintercept' %in% names(d)) d$yintercept else numeric()
  }))
  
  yval <- vals
  
  color <- c(color_text, rep("black", length(breaks)  ))
  setx <- c(yval, breaks)
  labs <- as.character(setx)
  name <- glue::glue("<i style='color:{color}'>{labs}")
  
  plot +
    ggplot2::geom_hline(yintercept = yval, 
                        linewidth = linewidth,
                        linetype = linetype,
                        color = color_line,
                        alpha = alpha) +
    ggplot2::scale_y_continuous(breaks = setx, labels = name) +
    ggplot2::theme(axis.text.y = ggtext::element_markdown())
}


#' An internal function to get the minimum version of packahege need
#'
#' @param pkg A character string specifying the package
#' @param version A numeric indicating the version to be returned
#' @param verbose A logical
#' @keywords internal
#' @return A character string.
#' @keywords internal
#' @noRd
#'
get_package_minversion <- function(pkg, version = NULL, verbose = FALSE) {
  if(!is.character(pkg)) stop('pkg must be a character')
  if(pkg == 'brms') {
    if(is.null(version)) {
      out <- '2.20.17' 
    } else {
      if(!is.character(version)) stop('version must be a character')
      out <- version
    }
  }
  if(pkg == 'marginaleffects') {
    if(is.null(version)) {
      out <- '0.18.0.9003'
    } else {
      if(!is.character(version)) stop('version must be a character')
      out <- version
    }
  }
  
  return(out)
}



#' An internal function to sanitize algorithm specific arguments
#'
#' @param args A list of argument to be sanitized
#' @param algorithm A character specifying the algorithm
#' @param verbose A logical
#' @keywords internal
#' @return A named list.
#' @keywords internal
#' @noRd
#'
sanitize_algorithm_args <- function(args, algorithm, verbose = FALSE) {
  if(!is.character(algorithm)) stop('algorithm must be a character')
  
  pathfinderargs <- c('save_latent_dynamics', 'output_dir',
                      'output_basename', 'sig_figs', 
                      'num_threads', 'init_alpha', 'tol_obj',
                      'tol_rel_obj', 'tol_grad', 'tol_rel_grad',
                      'tol_param', 'history_size', 'single_path_draws',
                      'draws', 'num_paths', 'max_lbfgs_iters', 
                      'num_elbo_draws', 'save_single_paths')
  
  laplacerargs <- c('save_latent_dynamics', 'output_dir',
                    'output_basename', 'sig_figs', 
                    'mode', 'opt_args', 'jacobian',
                    'draws')
  
  # if(algorithm == 'sampling') {
  #   return(args)
  # } else if(algorithm == 'meanfield') {
  #   return(args)
  # } else if(algorithm == 'fullrank') {
  #   return(args)
  # } else if(algorithm == 'fixed_param') {
  #   return(args)
  # }  
  
  if(!'pathfinder' %in% algorithm) {
    for (i in pathfinderargs) {
      if(!is.null(args[[i]])) args[[i]] <- NULL
    }
  } else if(!'laplace' %in% algorithm) {
    for (i in laplacerargs) {
      if(!is.null(args[[i]])) args[[i]] <- NULL
    }
  } else {
    args <- args
  }
  
  return(args)
}



