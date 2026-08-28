


#' Range by simulation or quantile method
#'
#' Compute a two-element range vector using either random normal simulation
#' or deterministic quantile-based constructions. The function is designed to
#' feel similar to [base::range()], but adds a `method` argument controlling
#' how values are obtained before the range is calculated.
#'
#' If `x` is `NULL`, the function generates values from a normal distribution.
#' Method `"r"` uses [stats::rnorm()] for random sampling, whereas methods
#' `"q1"`, `"q2"`, and `"q3"` use [stats::qnorm()] to construct deterministic,
#' symmetric values from normal quantiles.
#'
#' If `x` is supplied, then `x` is used directly and the arguments `n`,
#' `mean`, `sd`, and `seed` are ignored. In that case, method `"r"` returns
#' the ordinary range of `x`, while methods `"q1"`, `"q2"`, and `"q3"` return
#' the range of empirical quantiles of `x`, computed with [stats::quantile()].
#'
#' @param x An optional numeric vector. If supplied, the function works on
#'   `x` directly and ignores `n`, `mean`, `sd`, and `seed`.
#' @param ... Additional arguments passed on to [base::range()] when relevant.
#'   In practice this is mainly included for interface similarity with
#'   [base::range()].
#' @param method Character string specifying the method to use. Must be one of
#'   `"r"`, `"q1"`, `"q2"`, or `"q3"`.
#'   \describe{
#'     \item{`"r"`}{Random normal simulation via [stats::rnorm()] when `x` is
#'     `NULL`; otherwise uses `x` directly.}
#'     \item{`"q1"`}{Deterministic quantiles using
#'     `seq(0.001, 0.999, length.out = n)` when `x` is `NULL`, or the analogous
#'     empirical quantiles of `x` when `x` is supplied.}
#'     \item{`"q2"`}{Deterministic quantiles using
#'     `((1:n) - 0.5) / n` when `x` is `NULL`, or the analogous empirical
#'     quantiles of `x` when `x` is supplied.}
#'     \item{`"q3"`}{Explicit lower-tail and upper-tail symmetric quantiles when
#'     `x` is `NULL`, or the analogous empirical quantiles of `x` when `x` is
#'     supplied. Requires an even effective sample size.}
#'   }
#' @param n Integer sample size used only when `x` is `NULL`. Ignored if `x`
#'   is supplied.
#' @param mean Numeric mean of the normal distribution used only when `x` is
#'   `NULL`. Ignored if `x` is supplied.
#' @param sd Numeric standard deviation of the normal distribution used only
#'   when `x` is `NULL`. Ignored if `x` is supplied.
#' @param seed Optional integer seed used only when `x` is `NULL` and
#'   `method = "r"`. If not `NULL`, [base::set.seed()] is called before
#'   simulation. Ignored for quantile methods and ignored if `x` is supplied.
#' @param na.rm Logical; should missing values be removed before computing the
#'   result? This follows the usual meaning of `na.rm` in R.
#' @param finite Logical; should non-finite values be removed before computing
#'   the result?
#'
#' @return A numeric vector of length 2 containing the minimum and maximum of
#'   the constructed values, in the same style as [base::range()].
#'
#' @details
#' The function has two operating modes.
#'
#' \strong{1. `x` is `NULL`:} values are generated from a normal distribution.
#' Method `"r"` generates random values using [stats::rnorm()], so the result is
#' generally not exactly symmetric. Methods `"q1"`, `"q2"`, and `"q3"` use
#' [stats::qnorm()] with symmetric probabilities, so the resulting range is
#' deterministic and symmetric around `mean`.
#'
#' \strong{2. `x` is supplied:} the function does not generate from a normal
#' distribution. Instead, it uses `x` itself. In this case `n`, `mean`, `sd`,
#' and `seed` are ignored. Method `"r"` returns the ordinary range of `x`.
#' Methods `"q1"`, `"q2"`, and `"q3"` compute empirical quantiles of `x` using
#' [stats::quantile()] and then return the range of those quantiles.
#'
#' For method `"q3"`, the effective sample size must be even. If `x` is
#' supplied, this means that the length of `x` after any `na.rm`/`finite`
#' filtering must be even.
#'
#' @examples
#' # Random simulation from N(0, 1)
#' range_method(method = "r", n = 1000, mean = 0, sd = 1, seed = 123)
#'
#' # Deterministic symmetric ranges from normal quantiles
#' range_method(method = "q1", n = 1000, mean = 0, sd = 1)
#' range_method(method = "q2", n = 1000, mean = 0, sd = 1)
#' range_method(method = "q3", n = 1000, mean = 0, sd = 1)
#'
#' # Work directly on observed data
#' zz <- rnorm(1000, 0, 1)
#' range_method(x = zz, method = "r")
#' range_method(x = zz, method = "q1")
#' range_method(x = zz, method = "q2")
#'
#' # When x is supplied, n, mean, sd, and seed are ignored
#' range_method(x = zz, method = "q1", n = 1000, mean = 0, sd = 10, seed = 999)
#'
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#'
range_method <- function(x = NULL, 
                         ...,
                         method = c("r", "q1", "q2", "q3"),
                         n = 1000,
                         mean = 0,
                         sd = 1,
                         seed = 123,
                         na.rm = FALSE,
                         finite = FALSE) {
  
  method <- match.arg(method)
  
  if (!is.null(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    
    if (finite) {
      x <- x[is.finite(x)]
    }
    
    n_x <- length(x)
    
    vals <- switch(
      method,
      
      r = x,
      
      q1 = {
        probs <- seq(0.001, 0.999, length.out = n_x)
        stats::quantile(x, probs = probs, names = FALSE, na.rm = FALSE)
      },
      
      q2 = {
        probs <- ((1:n_x) - 0.5) / n_x
        stats::quantile(x, probs = probs, names = FALSE, na.rm = FALSE)
      },
      
      q3 = {
        if (n_x %% 2 != 0) {
          stop("For method = 'q3', length(x) after filtering must be even.")
        }
        p_low  <- seq(0.001, 0.499, length.out = n_x / 2)
        p_high <- seq(0.501, 0.999, length.out = n_x / 2)
        probs  <- c(p_low, p_high)
        sort(stats::quantile(x, probs = probs, names = FALSE, na.rm = FALSE))
      }
    )
    
    return(base::range(vals, ..., na.rm = FALSE, finite = FALSE))
  }
  
  vals <- switch(
    method,
    
    r = {
      if (!is.null(seed)) {
        base::set.seed(seed)
      }
      stats::rnorm(n = n, mean = mean, sd = sd)
    },
    
    q1 = {
      probs <- seq(0.001, 0.999, length.out = n)
      stats::qnorm(probs, mean = mean, sd = sd)
    },
    
    q2 = {
      probs <- ((1:n) - 0.5) / n
      stats::qnorm(probs, mean = mean, sd = sd)
    },
    
    q3 = {
      if (n %% 2 != 0) {
        stop("For method = 'q3', n must be even.")
      }
      p_low  <- seq(0.001, 0.499, length.out = n / 2)
      p_high <- seq(0.501, 0.999, length.out = n / 2)
      probs  <- c(p_low, p_high)
      sort(stats::qnorm(probs, mean = mean, sd = sd))
    }
  )
  
  base::range(vals, ..., na.rm = na.rm, finite = finite)
}



#' Evaluate priors defined in the Stan data block
#'
#' Extracts and organises prior information for a fitted model by combining
#' prior summaries and the corresponding Stan data into a structured data frame.
#'
#' @param model An object of class \code{bgmfit}.
#' @param spriors A prior object. If \code{NULL} (default),
#'   \code{\link[brms:prior_summary]{brms::prior_summary()}} is used to obtain
#'   \code{spriors} from \code{model}.
#' @param sdata A Stan data object. If \code{NULL} (default),
#'   \code{\link[brms:standata]{brms::standata()}} is used to obtain
#'   \code{sdata} from \code{model}.
#' @param prior_name_asit Logical (default \code{FALSE}) indicating whether
#'   prior names should be returned exactly as they appear in the Stan code.
#' @param gsub_group Character vector specifying group identifiers to remove
#'   from the \code{group} column of the prior object. Default \code{NULL}.
#' @param sort_response Character vector specifying the desired order of
#'   response variables used to sort the \code{resp} column of the prior
#'   object. Default \code{NULL}.
#' @param sort_parameter Character vector specifying the desired order of
#'   parameter names used to sort the \code{nlpar} column of the prior
#'   object. Default \code{NULL}.
#' @param sort_coefficient Character vector specifying the desired order of
#'   coefficient names used to sort the \code{coef} column of the prior
#'   object. Default \code{NULL}.
#' @param sort_class Character vector specifying the desired order of class
#'   names used to sort the \code{class} column of the prior object.
#'   Default \code{NULL}.
#' @param digits Integer giving the number of decimal places to use when
#'   rounding numeric values via \code{round()}.
#' @param viewer Logical (default \code{FALSE}) indicating whether to display
#'   the output in the R Viewer. Currently ignored to avoid a dependency on
#'   the \pkg{gt} package.
#' @param sort_dpar Logical (default \code{FALSE}) indicating whether to sort
#'   rows with \code{dpar == "sigma"} to the end of the prior object.
#' @param raw Logical (default \code{FALSE}) indicating whether to return the
#'   output in its original (unsorted/unrounded) format.
#'
#' @return A data frame containing prior information derived from
#'   \code{model}, \code{spriors}, and \code{sdata}.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
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
                               sort_dpar = TRUE,
                               raw = FALSE) {
  arguments <- as.list(match.call())[-1]
  if (missing(model)) {
    model <- NULL
  }
  if(!is.null(model)) {
    uvarby <- model$model_info$univariate_by$by
    if(is.null(uvarby)) uvarby <- NA 
  }
  
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
    stop2c("Supply either model or spriors and sdata arguments")
  } else if (!is.null(model) &
             !is.null(spriors) & !is.null(sdata)) {
    stop2c("Supply only model or spriors and sdata arguments")
  } else if (!is.null(model)) {
    spriors <- brms::prior_summary(model)
    sdata <- brms::standata(model)
  } else if (is.null(model)) {
    if (is.null(spriors) & is.null(sdata)) {
      stop2c("Supply spriors and sdata arguments")
    }
    if (is.null(spriors) & is.null(sdata)) {
      stop2c("Supply spriors and sdata arguments")
    }
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
    if(sort_dpar) {
      if(!all(remove_empty_string_from_vector(spriors[["dpar"]]) == "")) {
        spriors <- spriors %>% dplyr::arrange(!! as.symbol("dpar"))
      } 
    }
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
        sort_response <- model$model_info$yvars
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
    if (is.na(uvarby) &
        !model$model_info$multivariate$mvar) {
      spriors <- spriors %>%  dplyr::select(-'Response')
    }
  }
  
  if(sort_dpar) {
    if(!all(remove_empty_string_from_vector(spriors[["dpar"]]) == "")) {
      spriors <- spriors %>% dplyr::arrange(!! as.symbol("dpar"))
    } 
  }

  return(spriors)
}








#' Create a prior summary table for the bsitar model
#' @noRd
#' @exportS3Method prior_summary_table bgmfit
prior_summary_table.bgmfit <- function(model,
                                       set_width = c(0.95, 0.9999),
                                       set_digits = 1,
                                       empty = "-",
                                       print = FALSE,
                                       return_table = TRUE,
                                       tibble_table = TRUE,
                                       print_table = FALSE,
                                       return_file = NULL,
                                       flex_table = FALSE,
                                       path = NULL,
                                       title = NULL,
                                       align = "center",
                                       sheet_name = "table",
                                       draw_samples = 100000,
                                       add_range = FALSE,
                                       transform_class = NULL,
                                       transform_parameter = NULL,
                                       transform_fun = NULL,
                                       range_method_arg = NULL,
                                       seed = 123,
                                       verbose = FALSE,
                                       ...) {
  
  .dist_obj <- NULL
  .lower <- NULL
  .row_id  <- NULL
  .upper <- NULL
  .value <- NULL
  .width <- NULL
  .rule_id <- NULL
  ci <- NULL
  coefficient <- NULL
  dpar <- NULL
  draws <- NULL
  group <- NULL
  lb <- NULL
  level_lab <- NULL
  nlpar <- NULL
  parameter <- NULL
  resp <- NULL
  tag <- NULL
  ub <- NULL
  xmax_range <- NULL
  xmin_range <- NULL
  range_vec <- NULL
  zz <- NULL
  dist_key <- NULL
  range <- NULL
  
  insight::check_if_installed("flexlsx", prompt = FALSE)
  insight::check_if_installed("flextable", prompt = FALSE)
  insight::check_if_installed("distributional", prompt = FALSE)
  insight::check_if_installed("ggdist", prompt = FALSE)
  
  ggplot2::theme_set(ggdist::theme_ggdist())
  
  if (is.null(title)) title <- ""
  
  if (is.null(title)) {
    tab_name <- NULL
  } else if (identical(title, "")) {
    tab_name <- "Prior_summary"
  } else {
    tab_name <- title
  }
  
  if (!is.logical(add_range) || length(add_range) != 1 || is.na(add_range)) {
    stop("add_range must be a single TRUE or FALSE value.")
  }
  
  add_set_width <- TRUE
  if (!is.null(set_width)) {
    if(is.logical(set_width)) {
      if(!set_width) {
        add_set_width <- FALSE
      }
    } else if (!is.numeric(set_width)) {
      stop("set_width must be NULL, logical TRUE/FALSE or a numeric vector.")
    }
  }
  
  if(!add_set_width & !add_range) draw_samples <- 1
  
  if(!add_set_width) set_width <- NULL
  
  if (!is.null(set_width)) {
    if (!is.numeric(set_width)) {
      stop("set_width must be NULL or a numeric vector.")
    }
    if (length(set_width) < 1) {
      stop("If set_width is not NULL, it must contain at least one value.")
    }
    if (any(is.na(set_width))) {
      stop("set_width cannot contain NA values.")
    }
    if (any(set_width <= 0 | set_width >= 1)) {
      stop("Each set_width value must be strictly between 0 and 1.")
    }
  }
  
  format_transform_fun <- function(f) {
    if (!is.function(f)) {
      return("<not a function>")
    }
    
    txt <- paste(deparse(f), collapse = " ")
    txt <- gsub("[[:space:]]+", " ", txt)
    txt <- trimws(txt)
    
    if (nchar(txt) > 80) {
      txt <- paste0(substr(txt, 1, 77), "...")
    }
    
    txt
  }
  
  prior_object <-
    priors_to_textdata(
      model,
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
      digits = set_digits,
      viewer = FALSE,
      sort_dpar = TRUE,
      raw = TRUE
    ) %>%
    dplyr::filter(source == "user") %>%
    dplyr::filter(class != "L") %>% 
    dplyr::mutate(
      lb = dplyr::if_else(class == "sd", "0", lb),
      lb = dplyr::if_else(class == "sd" & dpar == "sigma", "0", lb),
      .row_id = dplyr::row_number()
    )
  
  
  
  prior_parsed <-
    prior_object %>%
    ggdist::parse_dist(prior, lb = "lb", ub = "ub") %>%
    dplyr::mutate(
      zz = paste0(class, coef, nlpar, dpar, group, resp)
    )
  
  if ("sigma" %in% transform_parameter) {
    if (any("sigma" %in% prior_parsed[["class"]])) {
      prior_parsed <- prior_parsed %>%
        dplyr::mutate(
          nlpar = dplyr::if_else(class == "sigma", "sigma", nlpar)
        )
      prior_parsed <- prior_parsed %>%
        dplyr::mutate(
          class = dplyr::if_else(nlpar == "sigma", "b", class)
        )
    } else if (any("sigma" %in% prior_parsed[["dpar"]])) {
      prior_parsed <- prior_parsed %>%
        dplyr::mutate(
          nlpar = dplyr::if_else(dpar != "" & nlpar == "", dpar, nlpar)
        )
    }
  }
  
  
  
  
  available_classes <- sort(unique(prior_parsed$class))
  available_parameters <- sort(unique(prior_parsed$nlpar))
  
  has_transform_args <- !is.null(transform_class) ||
    !is.null(transform_parameter) ||
    !is.null(transform_fun)
  
  if (isTRUE(has_transform_args)) {
    if (is.null(transform_class)) transform_class <- "b"
    if('sigma' %in% transform_parameter) {
      if(is.null(transform_fun)) {
        if(model$family$link_sigma == 'log') {
          transform_fun <- function(x)exp(x)
          if(verbose) {
            message2c("The link for the distributional parameter sigma  is 
                      'log', hence the automatic transformation applied is
                      'exp', ")
          }
        } else if(model$family$link_sigma == 'identity') {
          transform_fun <- function(x)(x)
        }
      }
    }
    
    
    if (is.null(transform_class) ||
        is.null(transform_parameter) ||
        is.null(transform_fun)) {
      stop2c("If any transformation arguments are supplied, transform_class,
             transform_parameter, and transform_fun must all be supplied.")
    }
    
    if (!is.character(transform_class)) {
      stop("transform_class must be a character vector.")
    }
    
    if (!is.character(transform_parameter)) {
      stop("transform_parameter must be a character vector.")
    }
    
    transform_class <- unique(transform_class)
    transform_parameter <- unique(transform_parameter)
    
    bad_class <- setdiff(transform_class, available_classes)
    if (length(bad_class) > 0) {
      stop(
        "Unknown transform_class value(s): ",
        paste(bad_class, collapse = ", "),
        ". Available classes are: ",
        paste(available_classes, collapse = ", "),
        "."
      )
    }
    
    bad_parameter <- setdiff(transform_parameter, available_parameters)
    if (length(bad_parameter) > 0) {
      stop(
        "Unknown transform_parameter value(s): ",
        paste(bad_parameter, collapse = ", "),
        ". Available parameters are: ",
        paste(available_parameters, collapse = ", "),
        "."
      )
    }
    
    if (is.function(transform_fun)) {
      transform_fun <- list(transform_fun)
    }
    
    if (!is.list(transform_fun)) {
      stop("transform_fun must be a function or a list of functions.")
    }
    
    if (!all(vapply(transform_fun, is.function, logical(1)))) {
      stop("Every element of transform_fun must be a function.")
    }
    
    n_class <- length(transform_class)
    n_parameter <- length(transform_parameter)
    n_combo <- n_class * n_parameter
    n_fun <- length(transform_fun)
    
    if (n_fun > n_combo) {
      transform_fun <- transform_fun[1:n_combo]
      n_fun <- length(transform_fun)
    }
    
    if (!(n_fun %in% c(1L, n_parameter, n_combo))) {
      stop(
        "transform_fun must have length 1, length(transform_parameter), or ",
        "length(transform_class) * length(transform_parameter)."
      )
    }
    
    transform_rules <- expand.grid(
      class = transform_class,
      parameter = transform_parameter,
      stringsAsFactors = FALSE,
      KEEP.OUT.ATTRS = FALSE
    )
    
    if (n_fun == 1L) {
      transform_rules$fun <- rep(transform_fun, nrow(transform_rules))
    } else if (n_fun == n_parameter) {
      param_fun_map <- transform_fun
      names(param_fun_map) <- transform_parameter
      transform_rules$fun <- unname(param_fun_map[transform_rules$parameter])
    } else {
      transform_rules$fun <- transform_fun
    }
    
    matched_rows <- unique(paste(prior_parsed$class, prior_parsed$nlpar, 
                                 sep = "___"))
    requested_rows <- unique(paste(transform_rules$class, 
                                   transform_rules$parameter, sep = "___"))
    missing_combos <- setdiff(requested_rows, matched_rows)
    
    if (length(missing_combos) > 0) {
      missing_labels <- vapply(
        strsplit(missing_combos, "___", fixed = TRUE),
        function(x) paste0("class = '", x[1], "', parameter = '", x[2], "'"),
        character(1)
      )
      stop(
        "Requested transformation combination(s) not found in priors: ",
        paste(missing_labels, collapse = "; "),
        "."
      )
    }
    
    transform_map_msg <- vapply(
      seq_len(nrow(transform_rules)),
      function(i) {
        paste0(
          "  - class = '", transform_rules$class[i],
          "', parameter = '", transform_rules$parameter[i],
          "' transformed as ", format_transform_fun(transform_rules$fun[[i]])
        )
      },
      character(1)
    )
    
    if (verbose) {
      message2c(
        "Applying transformations to prior draws for the following 
        class/parameter combinations:\n",
        paste(transform_map_msg, collapse = "\n")
      )
    }
  } else {
    transform_rules <- NULL
  }
  
  set.seed(seed)

  sim_tbl0 <-
    prior_parsed %>%
    dplyr::mutate(
      parameter = nlpar,
      dist_key = paste(prior, lb, ub, sep = "||")
    )
  
  # Handle constant priors such as constant(1)
  rconstant <- function (n, rate) rate
  
  draws_by_dist <-
    sim_tbl0 %>%
    dplyr::distinct(dist_key, .dist_obj) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      draws = list(
        as.numeric(distributional::generate(.dist_obj, times = draw_samples))
      )
    ) %>%
    dplyr::ungroup()
  
  sim_tbl <-
    sim_tbl0 %>%
    dplyr::left_join(
      draws_by_dist %>% dplyr::select(dist_key, draws),
      by = "dist_key"
    ) %>%
    dplyr::select(.row_id, class, parameter, .dist_obj, draws) %>%
    dplyr::arrange(.row_id)
  
  if (!is.null(transform_rules)) {
    transform_rules <- transform_rules %>%
      dplyr::mutate(.rule_id = dplyr::row_number())
    
    sim_tbl <- sim_tbl %>%
      dplyr::left_join(
        transform_rules %>% dplyr::select(class, parameter, .rule_id),
        by = c("class", "parameter")
      )
    
    sim_tbl <- sim_tbl %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        draws = list(
          if (is.na(.rule_id)) {
            draws
          } else {
            transform_result <- transform_rules$fun[[.rule_id]](draws)
            
            if (!is.numeric(transform_result)) {
              stop(
                "Each transform_fun must return a numeric vector. ",
                "Problem encountered for class = '", class,
                "', parameter = '", parameter, "'."
              )
            }
            
            if (length(transform_result) != length(draws)) {
              stop2c(
                "Each transform_fun must return a numeric vector of the 
                same length as its input. ",
                "Problem encountered for class = '", class,
                "', parameter = '", parameter, "'."
              )
            }
            
            as.numeric(transform_result)
          }
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.rule_id)
  }
  

  draws_long <-
    sim_tbl %>%
    tidyr::unnest(draws) %>%
    dplyr::rename(.value = draws)
  
  ci_labs <- character(0)
  
  if (is.null(set_width)) {
    ci_tbl_wide <- prior_parsed %>%
      dplyr::select(.row_id) %>%
      dplyr::distinct()
  } else {
    ci_tbl_long <-
      draws_long %>%
      dplyr::group_by(.row_id) %>%
      ggdist::point_interval(
        .value,
        .width = set_width,
        .point = median
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        level_lab = paste0(.width * 100, "% CI"),
        ci = paste0(
          sprintf("%0.2f", .lower),
          ", ",
          sprintf("%0.2f", .upper)
        )
       , ci = paste0("[", ci, "]") # comment out if dont want square brackets[]
      ) %>%
      dplyr::select(.row_id, level_lab, ci)
    
    ci_tbl_wide <-
      ci_tbl_long %>%
      tidyr::pivot_wider(
        id_cols = .row_id,
        names_from = level_lab,
        values_from = ci
      )
    
    ci_labs <- paste0(set_width * 100, "% CI")
  }
 

  if (is.null(range_method_arg)) {
    range_method_arg <- list()
  } else if (!is.list(range_method_arg)) {
    stop2c("range_method_arg must be a named list to pass arguments to 
           the range_method().", " The available arguments are: ",
           collapse_comma(methods::formalArgs(range_method)))
  }
  
  if (is.null(range_method_arg[["method"]])) range_method_arg[["method"]] <- "r"
  if (is.null(range_method_arg[["na.rm"]])) range_method_arg[["na.rm"]] <- TRUE
  if (is.null(range_method_arg[["seed"]])) range_method_arg[["seed"]] <- seed
  
  if (isTRUE(add_range)) {
    range_tbl <-
      draws_long %>%
      dplyr::group_by(.row_id) %>%
      dplyr::summarise(
        range_vec = list(
          do.call(range_method, c(list(x = .value), range_method_arg))
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        xmin_range = vapply(range_vec, `[`, numeric(1), 1),
        xmax_range = vapply(range_vec, `[`, numeric(1), 2),
        range = paste0(
          sprintf("%0.2f", xmin_range),
          ", ",
          sprintf("%0.2f", xmax_range)
        )
        , range = paste0("[", range, "]") # comment out if dont want square brackets[]
      ) %>%
      dplyr::select(.row_id, range)
  } else {
    range_tbl <- NULL
  }
  
  prior_object_range_ci <-
    prior_parsed %>%
    dplyr::left_join(ci_tbl_wide, by = ".row_id")
  
  
  
  if (isTRUE(add_range)) {
    prior_object_range_ci <-
      prior_object_range_ci %>%
      dplyr::left_join(range_tbl, by = ".row_id")
  }
  
  prior_object_range_ci <-
    prior_object_range_ci %>% 
    dplyr::mutate(
      ub = dplyr::if_else(ub == "", "Inf", ub),
      lb = dplyr::if_else(lb == "", "Inf", lb),
      prior = dplyr::if_else(
        ub == "Inf" & lb == "Inf",
        prior,
        # paste0(prior, "{", lb, ", ", ub, "}")
        paste0(prior, " ", "|", lb, ", ", ub, "|")
      )
    ) %>%
    dplyr::select(
      -dplyr::any_of(c(".dist", ".args", ".dist_obj")),
      -lb, -ub, -source
    ) %>%
    dplyr::mutate(
      coef = gsub("class", "", coef, ignore.case = FALSE)
    ) %>%
    dplyr::arrange(.row_id)
  
  range_lab <- if (isTRUE(add_range)) "Range" else NULL
  
  relocate_cols <- ci_labs[ci_labs %in% names(prior_object_range_ci)]
  if (length(relocate_cols) > 0) {
    prior_object_range_ci <-
      prior_object_range_ci %>%
      dplyr::relocate(dplyr::all_of(relocate_cols), .after = nlpar)
  }
  
  if ("range" %in% names(prior_object_range_ci)) {
    if (length(relocate_cols) > 0) {
      prior_object_range_ci <-
        prior_object_range_ci %>%
        dplyr::relocate(range, .after = dplyr::last_col())
    } else {
      prior_object_range_ci <-
        prior_object_range_ci %>%
        dplyr::relocate(range, .after = nlpar)
    }
  }
  

  prior_object_range_ci_out <- prior_object_range_ci
  
  prior_object_range_ci_out <- prior_object_range_ci_out %>% 
    dplyr::mutate(
      class = dplyr::if_else(class == "Intercept" & dpar == "sigma" & coef == "", 
                             "b", class),
      coef = dplyr::if_else(class == "b" & dpar == "sigma" & coef == "", 
                            "Intercept", coef)
    ) %>% 
    dplyr::mutate(
      nlpar = dplyr::if_else(dpar == "sigma" & nlpar == "", 
                            "sigma", nlpar),
      dpar = dplyr::if_else(nlpar == "sigma", 
                            "", dpar)
    )
  
  # when rsd_formual
  prior_object_range_ci_out <- prior_object_range_ci_out %>% 
    dplyr::mutate(
      coef = dplyr::if_else(class == "sigma" & dpar == "" & coef == "", 
                            "Intercept", coef)
    ) 

  prior_object_range_ci_out <- prior_object_range_ci_out %>% 
    dplyr::arrange(.row_id) %>%
    dplyr::mutate(
      coef = dplyr::if_else(class == "sd", paste0(coef, " (", group, ")"), coef)
    ) %>%
    dplyr::select(-resp, -dpar, -group, -zz, -.row_id) %>%
    dplyr::relocate(nlpar, .after = class) %>%
    dplyr::relocate(prior, .after = coef) %>%
    dplyr::rename(
      parameter = nlpar,
      coefficient = coef
    ) 
  
  prior_object_range_ci_out <- prior_object_range_ci_out %>% 
    dplyr::mutate(
      parameter = dplyr::if_else(class == "sigma", class, parameter),
      class = dplyr::if_else(parameter == "sigma", "rsd", class),
      coefficient = gsub("ClassI", "Class I", coefficient)
    )

  if ("tag" %in% names(prior_object_range_ci_out)) {
    tag_all_empty <- all(
      is.na(prior_object_range_ci_out$tag) |
        trimws(as.character(prior_object_range_ci_out$tag)) == ""
    )
    if (isTRUE(tag_all_empty)) {
      prior_object_range_ci_out <- prior_object_range_ci_out %>%
        dplyr::select(-tag)
    }
  }
  
  ci_cols <- ci_labs[ci_labs %in% names(prior_object_range_ci_out)]
  
  header_labs <- list(
    class = "Class",
    parameter = "Parameter",
    coefficient = "Coefficient",
    prior = "Prior distribution"
  )
  
  for (lab in ci_cols) {
    header_labs[[lab]] <- lab
  }
  
  if (isTRUE(add_range) && "range" %in% names(prior_object_range_ci_out)) {
    header_labs[["range"]] <- range_lab
  }
  
  add_header_lines_set <- tab_name
  
  foot_letters <- letters
  class_sym <- foot_letters[1]
  parameter_sym <- foot_letters[2]
  coefficient_sym <- foot_letters[3]
  prior_sym <- foot_letters[4]
  
  ci_syms <- if (length(ci_cols) > 0) {
    foot_letters[seq.int(from = 5, length.out = length(ci_cols))]
  } else {
    character(0)
  }
  
  range_sym <- if (isTRUE(add_range) && 
                   "range" %in% names(prior_object_range_ci_out)) {
    foot_letters[5 + length(ci_cols)]
  } else {
    NULL
  }
  
  transform_note_sym <- if (!is.null(transform_rules)) {
    foot_letters[5 + length(ci_cols) + as.integer(!is.null(range_sym))]
  } else {
    NULL
  }
  
  transform_desc <- NULL
  transform_note_text <- NULL
  
  if (!is.null(transform_rules)) {
    transform_desc <- unique(vapply(
      seq_len(nrow(transform_rules)),
      function(i) {
        paste0(
          "class = '", transform_rules$class[i],
          "', parameter = '", transform_rules$parameter[i], "'"
        )
      },
      character(1)
    ))
    
    transform_note_text <- paste0(
      "Transformations Note: ",
      paste(
        vapply(
          seq_len(nrow(transform_rules)),
          function(i) {
            paste0(
              "class = '", transform_rules$class[i],
              "', parameter = '", transform_rules$parameter[i],
              "' transformed as ",format_transform_fun(transform_rules$fun[[i]])
            )
          },
          character(1)
        ),
        collapse = "; "
      ),
      "."
    )
  }
  
  
  get_lab_map <- function(getnanmesx) {
    lab_map <- c(
      a = "a - size;",
      b = "b - timing;",
      c = "c - intensity;",
      d = "d - post-growth slope;",
      sigma = "sigma - within individual variability"
    )
    ord <- names(lab_map)[names(lab_map) %in% getnanmesx]
    out <- paste0("Parameter: ", paste(lab_map[ord], collapse = " "))
    return(out)
  }
  
  set_lab_map <- get_lab_map(prior_object_range_ci_out$parameter)

  summary_1 <-
    prior_object_range_ci_out %>%
    dplyr::ungroup() %>%
    flextable::flextable() %>%
    flextable::theme_apa() %>%
    flextable::merge_at(i = NULL, j = 1, part = "header") %>%
    flextable::merge_v(j = "class", part = "body") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::set_header_labels(values = header_labs) %>%
    flextable::footnote(
      i = NULL,
      j = 1,
      value = flextable::as_paragraph(
        paste0(
          "Class:",
          " b - regression parameters;",
          " sd - standard deviation for random effects;",
          " rsd - residual standard deviation"
        )
      ),
      ref_symbols = paste0(" ", class_sym, " "),
      part = "header"
    ) %>%
    flextable::footnote(
      i = NULL,
      j = 2,
      value = flextable::as_paragraph(
        set_lab_map
      ),
      ref_symbols = paste0(" ", parameter_sym, " "),
      part = "header"
    ) %>%
    flextable::footnote(
      i = NULL,
      j = 3,
      value = flextable::as_paragraph(
        paste0(
          "Coefficient:",
          " For Class b, the Intercept represents Class I estimate ",
          "whereas Class II denotes the difference between Class I ",
          "and Class II.",
          " The s parameters (s1, s2,...) are spline coefficients;",
          " For Class sd, the Intercept is standard deviation of ",
          "random effects for the group enclosed in the parentheses ",
          ";",
          " For Class rsd, Class I and Class II denote the ",
          "within-individual standard deviation estimates"
        )
      ),
      ref_symbols = paste0(" ", coefficient_sym, " "),
      part = "header"
    ) %>%
    flextable::footnote(
      i = NULL,
      j = 4,
      value = flextable::as_paragraph(
        paste0(
          "Prior distribution:",
          " Each coefficient is assigned a normal distribution ",
          "with mean and standard deviation specified in the ",
          "parentheses.",
          " Lower and/or upper  bounds, if any, are enclosed in the vertical ",
          " bars || next to the coefficient. ",
          " For instance, for half-normal distribution, bounds are |0, Inf| ",
          "."
        )
      ),
      ref_symbols = paste0(" ", prior_sym, " "),
      part = "header"
    )
  
  if (length(ci_cols) > 0) {
    for (i in seq_along(ci_cols)) {
      ci_note <- paste0(ci_cols[i]," credible interval mass for the estimates.")
      
      if (!is.null(transform_desc)) {
        ci_note <- paste0(
          ci_note,
          " Transformations applied for: ",
          paste(transform_desc, collapse = "; "),
          ".",
          " See 'Transformations Note' below for details."
        )
      }
      
      summary_1 <-
        summary_1 %>%
        flextable::footnote(
          i = NULL,
          j = which(names(prior_object_range_ci_out) == ci_cols[i]),
          value = flextable::as_paragraph(ci_note),
          ref_symbols = paste0(" ", ci_syms[i], " "),
          part = "header"
        )
    }
  }
  
  if ("range" %in% names(prior_object_range_ci_out)) {
    range_note <- paste0(
      "Range (min/max) based on base::range() applied to ",
      draw_samples,
      " simulated draws from the parsed prior distribution."
    )
    
    if (!is.null(transform_desc)) {
      range_note <- paste0(
        range_note,
        " Transformations applied for: ",
        paste(transform_desc, collapse = "; "),
        ".",
        " See 'Transformations Note' below for details."
      )
    }
    
    summary_1 <-
      summary_1 %>%
      flextable::footnote(
        i = NULL,
        j = which(names(prior_object_range_ci_out) == "range"),
        value = flextable::as_paragraph(range_note),
        ref_symbols = paste0(" ", range_sym, " "),
        part = "header"
      )
  }
  
  if (!is.null(transform_note_text)) {
    note_col <- if ("range" %in% names(prior_object_range_ci_out)) {
      "range"
    } else if (length(ci_cols) > 0) {
      ci_cols[length(ci_cols)]
    } else {
      "prior"
    }
    
    summary_1 <-
      summary_1 %>%
      flextable::footnote(
        i = NULL,
        j = which(names(prior_object_range_ci_out) == note_col),
        value = flextable::as_paragraph(transform_note_text),
        ref_symbols = paste0(" ", "", " "),
        part = "header",
        inline = FALSE
      )
  }
  
  summary_1 <-
    summary_1 %>%
    flextable::valign(valign = "center", part = "all") %>%
    flextable::add_header_lines(values = add_header_lines_set) %>%
    flextable::align(align = "center") %>%
    flextable::align(i = 1, j = NULL, align = "left", part = "header") %>%
    flextable::align(i = 2, j = NULL, align = "center", part = "header") %>%
    flextable::valign(i = 1, j = NULL, valign = "center", part = "header") %>%
    flextable::valign(i = 2, j = NULL, valign = "center", part = "header") %>%
    flextable::valign(i = NULL, j = 1, valign = "top", part = "body") %>%
    flextable::valign(i = NULL, j = 2, valign = "top", part = "body") %>%
    flextable::valign(i = NULL, j = 3, valign = "top", part = "body") %>%
    flextable::valign(i = NULL, j = 4, valign = "top", part = "body") %>%
    flextable::style(
      pr_t = flextable::fp_text_default(
        color = "black",
        font.size = 10,
        font.family = "Arial"
      ),
      part = "all"
    ) %>%
    flextable::style(
      i = 1,
      j = NULL,
      pr_t = flextable::fp_text_default(
        color = "black",
        font.size = 12,
        font.family = "Times New Roman"
      ),
      part = "header"
    ) %>%
    flextable::padding(padding.top = 0, part = "all") %>%
    flextable::padding(padding.bottom = 0, part = "all") %>%
    flextable::set_table_properties(width = 1, layout = "fixed") %>%
    flextable::width(j = 1, width = 0.6, unit = "in") %>%
    flextable::width(j = 2, width = 0.85, unit = "in") %>%
    flextable::width(j = 3, width = 1.25, unit = "in") %>%
    flextable::width(j = 4, width = 1.5, unit = "in")
  
  ncols_ft <- ncol(prior_object_range_ci_out)
  if (ncols_ft >= 5) {
    for (j in 5:ncols_ft) {
      summary_1 <- flextable::width(summary_1, j = j, width = 1.25, unit = "in")
    }
  }
  
  summary_1 <-
    summary_1 %>%
    flextable::border_remove() %>%
    flextable::border_inner_h(border = NULL, part = "header") %>%
    flextable::hline_top(border = NULL, part = "body") %>%
    flextable::hline_top(border = NULL, part = "footer")
  
  out_flex <- summary_1
  
  if (!is.null(title)) {
    out_flex <- flextable::set_caption(out_flex, caption = title)
  }
  
  if (print) print(out_flex$body$dataset)
  
  out <- export_flextable(
    ft = out_flex,
    return_file = return_file,
    path = path,
    title = title,
    align = align,
    sheet_name = sheet_name
  )
  
  if (is.null(return_file)) {
    #
  } else {
    return_table <- FALSE
  }
  
  if (return_table) {
    if (!flex_table) {
      out <- out$body$dataset
      if(tibble_table) out <- out %>% tibble::as_tibble()
      if(print_table) {
        print(knitr::kable(out))
        return(invisible(NULL))
      } else {
        return(out)
      }
    } else if (flex_table) {
      if (!is.null(title)) {
        out <- flextable::set_caption(out, caption = title)
      }
      return(out)
    }
  }
  
  if (!return_table) {
    export_flextable(
      ft = out_flex,
      return_file = return_file,
      path = path,
      title = title,
      align = align,
      sheet_name = sheet_name
    )
  }
  
  invisible(NULL)
}



#' @noRd
#' @exportS3Method prior_summary_table bgmfit
prior_summary_table <- function(model, ...) {
  UseMethod("prior_summary_table")
}

# Examples: prior_summary_table / prior_table

# prior_table(
#   model = model,
#   print_table = TRUE
# )
# 
# # Apply exp function to 'c' parameter
# prior_table(
#   model = model,
#   transform_class = c("b"),
#   transform_parameter = c("c"),
#   transform_fun = function(x) exp(x),
#   print_table = TRUE
# )
# 
# # Apply functions per selected parameter, recycled across all selected classes
# prior_table(
#   model = model,
#   transform_class = c("b", "sd"),
#   transform_parameter = c("b", "c"),
#   transform_fun = list(
#     function(x) x,
#     function(x) exp(x)
#   ),
#   print_table = TRUE
# )
# 
# # Apply functions per expanded class/parameter combination
# prior_table(
#   model = model,
#   transform_class = c("b", "b"),
#   transform_parameter = c("c", "c"),
#   transform_fun = list(
#     function(x) x,
#     function(x) exp(x),
#     function(x) x,
#     function(x) exp(x)
#   ),
#   print_table = TRUE
# )





#' Check that time is a positive integer in increasing order within each
#' group
#'
#' This function verifies that a given variable (typically time/age) is:
#' \itemize{
#'   \item A positive integer (greater than 0)
#'   \item Increasing by exactly 1 within each group (e.g., 1, 2, 3, ...)
#'   \item Strictly increasing within each group
#' }
#'
#' The checks are performed separately for each unique value of the grouping
#' variable.
#'
#' @param data A data frame containing the variables specified in \code{idvar}
#'   and \code{xvar}.
#' @param idvar A string specifying the name of the grouping variable (e.g.,
#'   individual ID). The checks are performed within each unique value of this
#'   variable.
#' @param xvar A string specifying the name of the variable to check (e.g.,
#'   "time" or "age"). This variable will be converted to numeric for the
#'   checks.
#' @param flag A logical flag indicating whether to stop with an error if the
#'   check fails. If \code{TRUE} (default) and the check fails, an error is
#'   thrown via \code{stop2c()}. If \code{FALSE}, the function returns the
#'   logical result without stopping.
#' 
#' @return A single logical value:
#'   \itemize{
#'     \item \code{TRUE} if all groups satisfy all three conditions (positive
#'     integer, increasing by 1, strictly increasing)
#'     \item \code{FALSE} if any group fails any of the conditions
#'   }
#'   When \code{flag = TRUE} and the result is \code{FALSE}, an error is thrown
#'   instead of returning.
#'
#' @examples
#' # Example 1: Valid data (passes all checks)
#' df_valid <- data.frame(
#'   id = c("A", "A", "A", "B", "B", "B"),
#'   age = c("1", "2", "3", "2", "3", "4")
#' )
#' check_id_xvar(df_valid, "id", "age")  # Returns TRUE
#'
#' # Example 2: Invalid data (fails increasing by 1)
#' df_invalid <- data.frame(
#'   id = c("A", "A", "A", "B", "B", "B"),
#'   age = c("1", "2", "3", "2", "4", "5")  # B jumps from 2 to 4
#' )
#' check_id_xvar(df_invalid, "id", "age")  # Returns FALSE
#'
#' # Example 3: With flag = FALSE (no error thrown)
#' check_id_xvar(df_invalid, "id", "age", flag = FALSE)  # Returns FALSE silently
#'
#' # Example 4: Will throw error with default flag = TRUE
#' # check_id_xvar(df_invalid, "id", "age")  # Throws error via stop2c()
#'
#' @details
#' The function performs three checks within each group:
#' \itemize{
#'   \item \code{positive_integer}: All values are positive integers (> 0)
#'   \item \code{increasing_by_1}: Values increase by exactly 1 (e.g., 1, 2, 3)
#'   \item \code{strictly_increasing}: Values are strictly increasing (any
#'   positive step)
#'  }
#'
#' The function returns \code{TRUE} only if ALL three conditions are satisfied
#' for ALL groups.
#'
#' @note
#' \itemize{
#'   \item The \code{xvar} column can be stored as strings; it will be converted
#'   to numeric internally.
#'   \item Non-numeric strings in \code{xvar} become \code{NA} and cause the
#'   check to fail.
#'   \item The data is sorted by \code{idvar} and \code{xvar} before checking,
#'   so input order doesn't matter.
#'   \item When \code{flag = TRUE} and the check fails, \code{stop2c()} is
#'   called with a descriptive error message.
#' }
#'
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#'
check_id_xvar <- function(data, idvar, xvar, flag = TRUE, append_msg = "") {
  xvar_num <- NULL;
  positive_integer <- NULL;
  increasing_by_1 <- NULL;
  strictly_increasing <- NULL; 
  all_valid <- NULL; 
  variables <- c(idvar, xvar)
  for (j in variables) {
      if(length(data[[j]]) == 0) {
        stop2c("The variable ", collapse_comma(j), 
               " is missing. Check your data")
      }
  }
  out <- data %>%
    dplyr::mutate(xvar_num = as.numeric(get(xvar))) %>%
    dplyr::arrange(get(idvar), xvar_num) %>%
    dplyr::group_by(get(idvar)) %>%
    dplyr::summarise(
      positive_integer = all(
        xvar_num > 0 &
          xvar_num == as.integer(xvar_num) &
          !is.na(xvar_num)
      ),
      increasing_by_1 = all(diff(xvar_num) == 1),
      strictly_increasing = all(diff(xvar_num) > 0),
      .groups = "drop"
    ) %>%
    dplyr::summarise(
      all_valid = all(positive_integer, increasing_by_1, strictly_increasing)
    ) %>%
    dplyr::pull(all_valid)
  
  if(flag) {
    if(!out) {
      msg <- paste0("The time variable ", collapse_comma(xvar),
             " must be a positive integer in increasing order (e.g., 1, 2, 3) ",
             "for each individual identified by the variable ", 
             collapse_comma(idvar), ".")
      if(append_msg != "") msg <- paste0(msg, " ", append_msg)
      stop2c(msg)
    }
  }
  return(out)
}




make_id_xvar <- function(data, idvar, xvar, timevar, resp = NULL, nys = 1) {
  ysi <- NULL;
  oooooooooooo <- NULL;
  if(nys > 1) {
    newtimevar <- paste0(timevar, "_", ysi)
  } else {
    newtimevar <- timevar
  }
  out <- data %>% 
    dplyr::mutate(oooooooooooo = as.factor(dplyr::cur_group_id())) %>% 
    dplyr::arrange(idvar, xvar) %>% 
    dplyr::group_by_at(c(idvar)) %>% 
    dplyr::mutate(!! as.symbol(newtimevar) := 
                    factor(
                      base::seq_along(.data[[xvar]]), 
                      labels = paste0("T", base::seq_along(.data[[xvar]]))
                    )
                  ) %>%
    dplyr::arrange(oooooooooooo) %>% 
    dplyr::select(-dplyr::all_of('oooooooooooo')) %>%
    dplyr::ungroup()
  if(nys > 1) {
    out <- out %>% dplyr::select(-dplyr::all_of(timevar))
  }
  attr(out, 'newtimevar') <- newtimevar
  return(out)
}



#' Checks if object is of class \code{bgmfit}
#'
#' @param x An \R object
#' 
#' @inherit berkeley author
#'
#' @export
is.bgmfit <- function(x) {
  inherits(x, "bgmfit")
}



#' Evaluate Global Arguments in a Call Object
#' 
#' Modifies a call object (typically from [match.call()]) by evaluating and
#' replacing argument expressions that resolve successfully in the specified
#' environment (typically global). Skips arguments listed in \code{exceptions}.
#' 
#' @param mcall A call object, usually obtained from
#'   \code{match.call(expand.dots = FALSE)}. Arguments will be selectively
#'   evaluated and replaced with their values.
#' 
#' @param envir An environment in which to evaluate argument expressions.
#'   Defaults to [globalenv()]. Only expressions that evaluate without error in
#'   this environment (with no parent lookup via [emptyenv()]) are replaced.
#' 
#' @param exceptions [NULL] (default) or a \code{character} vector of argument
#'   names to skip. These arguments retain their original unevaluated
#'   expressions.
#' 
#' @return The modified \code{mcall} object with selected arguments replaced by
#'   their evaluated values (e.g., symbols become actual objects, lists/data
#'   frames are materialized). Suitable for subsequent \code{eval(mcall,
#'   parent.frame())}.
#' 
#' @examples
#' age <- 1:5
#' df <- data.frame(age = 1:10)
#' 
#' myfun <- function(x, data, y) {
#'   mcall <- match.call(expand.dots = FALSE)
#'   mcall <- eval_globals_in_mcall(mcall, exceptions = "data")
#'   eval(mcall, parent.frame())
#' }
#' 
#' myfun(x = age, data = df, y = df$age)
#' # Returns: list(x = 1:5, data = df, y = 1:10)
#' 
#' @seealso [match.call()], [eval()], [globalenv()]
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#'
eval_globals_in_mcall <- function(mcall, envir = globalenv(), 
                                  exceptions = NULL) {
  arg_names <- names(mcall)[-1]
  arg_names <- setdiff(arg_names, exceptions)
  arg_names <- arg_names[arg_names != '...']
  
  arg_names_formula_names <- c()
  for (i in arg_names) {
    if(rlang::is_formula(mcall[[i]])) {
      arg_names_formula_names <- c(arg_names_formula_names, i)
    }
  }
  
  exceptions <- c(exceptions, arg_names_formula_names)
  arg_names <- setdiff(arg_names, exceptions)
  
  for (nm in arg_names) {
    expr <- mcall[[nm]]
    if (is.recursive(expr) && length(expr) == 3 && identical(expr[[1]], 
                                                             quote(`[[`))) {
      obj_expr <- expr[[2]]
      idx_expr <- expr[[3]]
      obj_val <- try(eval(obj_expr, envir = envir, enclos = emptyenv()), 
                     silent = TRUE)
      idx_val <- try(eval(idx_expr, envir = envir, enclos = emptyenv()), 
                     silent = TRUE)
      if (!inherits(obj_val, "try-error") && !inherits(idx_val, "try-error") &&
          (is.list(obj_val) || is.vector(obj_val)) && 
          !is.null(tmp <- obj_val[[idx_val]]) && !is.function(tmp)) {
        if(!is.list(tmp)) {
          if(check_is_numeric_like(tmp)) tmp <- eval(tmp)
        }
        mcall[[nm]] <- tmp  # Single extracted value
        next
      }
    }

    if (is.recursive(expr) && length(expr) == 3 && identical(expr[[1]], 
                                                             quote(`$`))) {
      
      expr <- dollar_to_double_bracket(expr)
      obj_expr <- expr[[2]]
      idx_expr <- expr[[3]]
      obj_val <- try(eval(obj_expr, envir = envir, enclos = emptyenv()), 
                     silent = TRUE)
      idx_val <- try(eval(idx_expr, envir = envir, enclos = emptyenv()), 
                     silent = TRUE)
      if (!inherits(obj_val, "try-error") && !inherits(idx_val, "try-error") &&
          (is.list(obj_val) || is.vector(obj_val)) && 
          !is.null(tmp <- obj_val[[idx_val]]) && !is.function(tmp)) {
        if(check_is_numeric_like(tmp)) tmp <- eval(tmp)
        mcall[[nm]] <- tmp  
        next
      }
    }

    is_wrapped <- FALSE
    wrapper <- NULL
    inner_expr <- expr
    if (is.recursive(expr) && length(expr) >= 2 && 
        (identical(expr[[1]], quote(`list`)) || identical(expr[[1]], 
                                                          quote(`c`)))) {
      wrapper <- expr[[1]]
      inner_expr <- expr[-1]
      is_wrapped <- TRUE
    }

    inner_exprde_e <- deparse_0(inner_expr)
    if(grepl("\\+", inner_exprde_e) & !is_wrapped) {
      inner_exprde_e <- gsub("\"", "", inner_exprde_e)
      inner_expr <- inner_exprde_e
    }
    
    if(grepl("_prior_", nm)) {
      if(!is.list(inner_expr)) {
        inner_expr <- deparse_0(inner_expr)
      }
    }

    if (is.recursive(inner_expr) & !is.null(wrapper)) { 
      evaled_elements <- lapply(inner_expr, function(e) {
        val <- try(eval(e, envir = envir, enclos = emptyenv()), 
                   silent = TRUE)
        if (inherits(val, "try-error") || is.null(val) || is.function(val)) {
          e
        } else {
          de_e <- deparse_0(e)
          de_e <- gsub("\"", "", de_e)
          if(grepl("+", de_e)) {
            val <- de_e
          } else {
            val <- val
          }
          if(check_is_numeric_like(val)) val <- eval(val)
          val
        }
      })

      if (is_wrapped) {
        mcall[[nm]] <- as.call(c(wrapper, evaled_elements))
      } else if (is.vector(inner_expr) && !is.list(inner_expr)) {
        mcall[[nm]] <- unname(as.vector(evaled_elements, mode = typeof(expr)))
      } else {
        names(evaled_elements) <- names(inner_expr)
        mcall[[nm]] <- evaled_elements
      }
    } else if (is.recursive(inner_expr) & is.null(wrapper)) { # added for NULL
      mcall[[nm]] <- NULL # added for NULL
    } else {
      val <- try(eval(expr, envir = envir, enclos = emptyenv()), silent = TRUE)
 
      if (!inherits(val, "try-error") && !is.null(val) && !is.function(val)) {
        if(is.null(expr)) {
          mcall[[nm]] <- expr
        } else if(!is.null(expr)) { 
          if(is.logical(expr)) mcall[[nm]] <- expr
        } else if(!is.null(expr) & rlang::is_bare_numeric(expr)) {
          mcall[[nm]] <- expr
        } else if(!is.null(expr) & !rlang::is_bare_numeric(expr)) {
          de_e <- deparse_0(val)
          de_e <- gsub("\"", "", de_e)
          if(grepl("+", de_e)) {
            val <- de_e
          } else {
            val <- val
          }
          if(check_is_numeric_like(val)) val <- eval(val)
          mcall[[nm]] <- val
        } 
      }
    }
  }
  return(mcall)
}



#' Convert \code{'$'} notation to bracket notation
#'
#' Internal helper to rewrite expressions of the form \code{"list$name"}
#' into bracket notation (for example, \code{list[["name"]]}), typically
#' for safer programmatic evaluation.
#'
#' @param expr A character string or expression, such as \code{"list$name"}.
#'
#' @return A list representing the corresponding function call and arguments.
#' 
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
dollar_to_double_bracket <- function(expr) {
  obj_expr <- expr[[2]]
  idx_expr <- expr[[3]]
  zz <- deparse(idx_expr)
  zx <- deparse(obj_expr)
  expr <- str2lang(paste0(zx, '[[', "'", zz, "'", ']]'))
  return(expr)
}


#' An internal function to remove empty elemnsts
#'
#' @param x A string
#' @return A list comprised of function arguments.
#' 
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#'
remove_empty_string_from_vector <- function(x) {
  if(is.list(x)) {
    x <- names(x)
  }
  x[x!=""]
}


#' Check whether a function was called via \code{do.call}
#'
#' Internal helper to determine whether a function was invoked using
#' \code{do.call}, typically used in metaprogramming contexts.
#'
#' @param x A character string.
#'
#' @return A list representing the corresponding function call and arguments.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
called_via_do_call <- function() {
  calls <- sys.calls()
  any(vapply(calls, function(cl) identical(cl[[1L]], quote(do.call)), 
             logical(1)))
}


#' Check whether a function was called via \code{CustomDoCall}
#'
#' Internal helper to determine whether a function was invoked using
#' \code{CustomDoCall}, typically used in custom metaprogramming contexts.
#'
#' @param x A character string.
#'
#' @return A list representing the corresponding function call and arguments.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
called_via_CustomDoCall <- function() {
  calls <- sys.calls()
  any(vapply(calls, function(cl) identical(cl[[1L]], quote(CustomDoCall)), 
             logical(1)))
}





#' Remove excess spaces from a string
#'
#' Internal helper to collapse multiple consecutive spaces into a single space
#' within a character string.
#'
#' @param x A character string.
#'
#' @return A character string with excess spaces removed.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
clean_text_spaces <- function(x) {
  trimws(gsub("\\s+", " ", x))
}


#' Extract function arguments
#'
#' Internal helper to extract and standardise arguments from a function call.
#'
#' @param arguments A list of default function arguments.
#' @param xcall A character string specifying the name of the calling function.
#'
#' @return A list of function arguments.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
get_args_ <- function(arguments, xcall, xclass = NULL, scallstatus = NULL) {
  `%!in%` <- Negate(`%in%`)
  if(is.null(xclass)) {
    pastexclass <- paste0(".", 'bgmfit')
  } else if(xclass == "") {
    pastexclass <- ""
  } else {
    if( grepl("^\\.", xclass)) pastexclass <- xclass
    if(!grepl("^\\.", xclass)) pastexclass <- paste0(".", xclass)
  }
  if(grepl(pastexclass, xcall)) {
    pastexclass <- ""
  }
  enverr. <- environment()
  assign('err.', FALSE, envir = enverr.)
  tryCatch(
    expr = {
      f_funx_arg <- formals(paste0(xcall, pastexclass))
    },
    error = function(e) {
      assign('err.', TRUE, envir = enverr.)
    }
  )
  err. <- get('err.', envir = enverr.)
  if (err.) {
    if(is.null(scallstatus)) stop2c("set scallstatus via sys.status()")
    f_funx_arg <- get_xcall_byclass(scallstatus, pastexclass)
  } else {
    f_funx_arg <- f_funx_arg
  }
  nf_funx_arg_names <-
    intersect(names(arguments), names(f_funx_arg))
  arguments <-
    c(arguments, f_funx_arg[names(f_funx_arg) %!in% nf_funx_arg_names])
  arguments
}



#' Deparse a symbol and remove spaces
#'
#' Internal helper to convert a symbol to a character string and remove
#' any excess spaces.
#'
#' @param deparseobj A symbol.
#'
#' @return A character string.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
deparse_0 <- function(deparseobj) {
  deparseobj <- paste(deparse(deparseobj), collapse = "")
  deparseobj <- gsub("[[:space:]]", "", deparseobj)
  deparseobj
}




#' Check for a pipe character in a string
#'
#' Internal helper to detect the pipe character (\code{"|"}) in a symbol or
#' character string.
#'
#' @param x A symbol or character string.
#' @param return_name Logical; if \code{TRUE}, return the original string
#'   (when a pipe is found); if \code{FALSE}, return a logical indicating
#'   whether a pipe was found.
#'
#' @return A character string (if \code{return_name = TRUE}) or a logical
#'   (if \code{return_name = FALSE}).
#'   
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
check_forpipe <- function(x, return = 'name') {
  data_name_str   <- x
  if(is.symbol(data_name_str)) {
    data_name_str <- deparse(data_name_str)
  }
  data_name_split <- paste(gsub_space(data_name_str), collapse = "")
  data_name_pipe  <- FALSE
  if(grepl("%>%", data_name_split, fixed = T)) {
    data_name_str_attr <- strsplit(data_name_split, "%>%", fixed = T)[[1]][1]
    data_name_pipe  <- TRUE
  } else if(grepl("(", data_name_split, fixed = T) &&
            grepl(",", data_name_split, fixed = T)) {
    data_name_str_attr <- strsplit(data_name_split, "|>", fixed = T)[[1]][1]
    data_name_pipe  <- TRUE
  } else {
    data_name_pipe  <- FALSE
    data_name_str_attr <- data_name_split
  }
  if(return == 'name') {
    out <- data_name_str
  } else if(return == 'attr') {
    out <- data_name_str_attr
  } else if(return == 'logical') {
    out <- data_name_pipe
  } else {
    stop2c("argument 'return' must be either 'name', 'attr', or 'logical'")
  }
  return(out)
}



#' Substitute and deparse a symbol argument
#'
#' Internal helper to substitute a symbol (typically via \code{substitute()})
#' and then convert it to a character string (via \code{deparse()}).
#'
#' @param deparseobj A symbol.
#'
#' @return A character string.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
deparse_0s <- function(deparseobj) {
  deparseobj <- paste(deparse(substitute(deparseobj)), collapse = "")
  deparseobj
}


#' Remove spaces from a string
#'
#' Internal helper to remove all spaces from a character string.
#'
#' @param deparseobj A character string.
#'
#' @return A character string with spaces removed.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
gsub_space <- function(deparseobj) {
  deparseobj <- gsub("[[:space:]]", "", deparseobj)
  deparseobj
}



#' Remove spaces from a string
#'
#' Internal helper to remove all spaces from a character string.
#'
#' @param deparseobj A character string.
#'
#' @return A character string with spaces removed.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
gsub_quote1 <- function(deparseobj) {
  gsub("\"", "", deparseobj)
}



#' Get arguments from the global environment
#'
#' Internal helper to extract function arguments from the global environment
#' and related environments.
#'
#' @param mcallarg An argument from a \code{mcall()} call.
#' @param envir An environment used for function evaluation.
#' @param search_envir An environment in which to search for objects used as
#'   arguments.
#' @param exceptions A character vector indicating exceptions.
#' @param ... Additional arguments.
#'
#' @return A list of function arguments.
#'
#' @inherit berkeley author
#' 
#' @keywords internal
#' @noRd
#' 
mcall_dictionary <- function(mcallarg, 
                             envir = NULL, 
                             xenvir = NULL, 
                             exceptions = NULL, 
                             ...) {
  mcallx <- mcallarg
  if(is.null(envir)) {
    enverr. <- environment()
  } else {
    enverr. <- envir
  }
  if(is.null(xenvir)) {
    searchenvir. <- .GlobalEnv
  } else {
    searchenvir. <- xenvir
  }
  for(i in names(mcallx)) {
    if(!i %in% exceptions) {
      if(!is.null(mcallx[[i]])) {
        gxz <- mcallx[[i]]
        assign('err.', FALSE, envir = enverr.)
        tryCatch(
          expr = {
            getgxz <- get(gxz, envir = searchenvir.)
          },
          error = function(e) {
            assign('err.', TRUE, envir = enverr.)
          }
        )
        err. <- get('err.', envir = enverr.)
        if (err.) {
          mcallx[[i]] <- gxz
        } else {
          validca <- getgxz
          if(is.symbol(validca)) {
            mcallx[[i]]  <- validca
          } else if(is.character(validca)) {
            mcallx[[i]] <-  validca
          } else if(is.numeric(validca)) {
            mcallx[[i]]  <- validca
          } else if(is.data.frame(validca)) {
            mcallx[[i]] <-  gxz # note gxz and not 
          } else if(tibble::is_tibble(validca)) {
            mcallx[[i]] <-  gxz # note gxz and not 
          } else if(is.language(validca)) {
            mcallx[[i]] <-  validca
          } else if(is.list(validca)) {
            mcallx[[i]] <-  validca 
          } else if(is.vector(validca)) {
            mcallx[[i]] <-  validca 
          }
        }
      } 
    } 
  } 
  return(mcallx)
}



#' Expose function after optimization
#'
#' Internal helper to expose a function after model optimization.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @return A list of exposed functions.
#' 
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
expose_optimize_fit <- function(model,
                                subset_list = NULL,
                                expose_function = T) {
  optimize_fit_models <-  model

  if(!is.null(subset_list)) {
    if(!is.numeric(subset_list)) stop2c("models must a numeric vector")
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
        message2c("Exposing Stan function...", " (for model no. ", il, ")")
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
  return(m_list)
}



#' Process models after optimization
#'
#' Internal helper to process models after optimization.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @return A list of plot objects.
#' 
#' @inherit berkeley author
#'
#' @keywords internal
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
    if(!i %in% methods::formalArgs(plot_curves.bgmfit))
      stop2c("arguments must be be one of the following",
           "\n ",
           methods::formalArgs(plot_curves.bgmfit))
  }
  if(!is.null(subset_list)) {
    if(!is.numeric(subset_list)) stop2c("models must a numeric vector")
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
        message2c("Exposing Stan function...", " (for model no. ", il, ")")
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
    message2c("Working on model no. ", .x)
    dots$model <- bx[[.x]]
    dots$... <- NULL
    if(is.null(what)) what <- 'plot'
    if(what == "plot") {
      out_ <- CustomDoCall(plot_curves, dots)
      title_ <- bx[[.x]]$model_info$optimization_info
      out_ <- out_ + ggplot2::labs(title = title_)
    }
    if(what == "growthparameters") {
      out_ <- CustomDoCall(growthparameters, dots)
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




#' Evaluate arguments ending with the \code{_str} suffix
#'
#' Internal helper to evaluate arguments that end with the \code{_str} suffix
#' and convert them to character strings using the provided data.
#'
#' @param tsx An argument with the \code{_str} suffix.
#' @param data A data frame.
#'
#' @return A list of character strings.
#' 
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get_gr_str_coef_id <- function(tsx,
                               data) {
  
  tsx <- strsplit(tsx, "+(", fixed = T)[[1]]
  tsxi <- tsx
  tsxi_c <- c()
  for (i in 1:length(tsxi)) {
    strpartstrx <- tsxi[i]
    strpartstrx <- strsplit(strpartstrx, "|", fixed = T)[[1]]
    strpartstrx_form <- strpartstrx[1]
    strpartstrx_form <-  gsub("~(", "~",  strpartstrx_form, fixed = T)
    if(length(strpartstrx) > 1 ) {
      strpartstrx_grpa <- strpartstrx[2:length(strpartstrx)]
      strpartstrx_grpa <- gsub("[()]", "", strpartstrx_grpa)
      strpartstrx_grpa2 <- paste0("", strpartstrx_grpa, collapse = "|")
      tsx_t <- paste0(strpartstrx_form, "|", strpartstrx_grpa2)
    } else {
      tsx_t <- strpartstrx_form
    }
    tsx <- c(tsxi_c, tsx_t)
  }
  tsx_id_w_or_wo_gr <- c()
  for (tsx_id_w_or_wo_gri in 1:length(tsx)) {
    tsx_id_w_or_wo_gr_get <- get_x_random2_asitis(tsx[tsx_id_w_or_wo_gri])
    tsx_id_w_or_wo_gr <- c(tsx_id_w_or_wo_gr, tsx_id_w_or_wo_gr_get)
  }
  tsx_c_coef  <- tsx_c_id    <- set_form_gr_it      <- list()
  set_ncov_it <- set_corr_it <- set_corr_true_false <- list()
  for (i in 1:length(tsx)) {
    tsx_c <- strsplit(tsx[i], "|", fixed = T)[[1]]
    set_corr_it_get <- tsx_c[2]
    tsx_c1 <- tsx_c[1]
    tsx_c3 <- tsx_c[3]
    if(grepl("^\\(", tsx_c1)) tsx_c1 <- gsub("^\\(", "", tsx_c1)
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
  } 
  if(length(tsx_c_coef) != length(tsx_c_id))
    stop2c("coef and id length should be same")
  list(tsx_c_coef = tsx_c_coef, tsx_c_id = tsx_c_id,
       set_form_gr_it = set_form_gr_it, set_ncov_it = set_ncov_it,
       set_corr_it = set_corr_it, set_corr_true_false = set_corr_true_false)
}




#' Extract correlation structure from \code{||} syntax for \code{_str} arguments
#'
#' Internal helper to extract correlation structure from the \code{||} syntax
#' for arguments ending with the \code{_str} suffix.
#'
#' @param str_id_all_list An argument with the \code{_str} suffix for \code{id}.
#' @param str_corr_all_list An argument with the \code{_str} suffix for
#'   \code{gr_cor}.
#' @param str_corr_tf_all_list An argument with the \code{_str} suffix for
#'   \code{corr}.
#'
#' @return A list of character strings.
#' 
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
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


#' Append priors to \code{bpriors}
#'
#' Internal helper to append prior information to the \code{bpriors} object.
#'
#' @param tempx A prior object.
#'
#' @return A prior object.
#' 
#' @inherit berkeley author
#'
#' @keywords internal
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




#' Restore parentheses in formula objects
#'
#' Internal helper to restore parentheses in formula objects.
#'
#' @param strx A formula object.
#' @param exclude_first Logical indicating whether to exclude the fixed part
#'   from adding opening and closing parentheses. If \code{NULL}, it is set to
#'   \code{TRUE} for \code{sigma} and \code{FALSE} otherwise. This could be
#'   \code{TRUE} globally, but needs to be checked.
#'
#' @return A character string.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
restore_paranthese_grgr_str_form <- function(strx, exclude_first = NULL) {
  if(is.null(exclude_first)) {
    if(grepl("sigma~", strx)) {
      exclude_first <- TRUE
    } else {
      exclude_first <- FALSE
    }
  } 
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
          stop2c("Underscore '_' is not allowed in the variable name when",
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
    if(exclude_first) {
      if(abxi == 1) {
        set_abx_c <- abxs[abxi]
      } else {
        set_abx_c <- restore_paranthese_grgr_str(abxs[abxi])
      }
    } else if(!exclude_first) {
      set_abx_c <- restore_paranthese_grgr_str(abxs[abxi])
    }
    abx_c[abxi] <- set_abx_c
    # abx_c[abxi] <- restore_paranthese_grgr_str(abxs[abxi])
    # abx_c[abxi]
  }
  abx_c <- paste0(abx_c, collapse = "+")
  return(abx_c)
}


#' Get random effect formula arguments
#'
#' Internal helper to extract random effect formula arguments.
#'
#' @param x A character string of a random effect formula.
#'
#' @return A list of character strings.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get_x_random2 <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  x <- strsplit(x, ")+" )[[1]]
  x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
  if(any(grepl("^|gr", x)) | !any(grepl("^|gr", x))  ) {
    x <- sub(".*gr", "", x)
    x_c <- c()
    for (xi in 1:length(x)) {
      if(!grepl("^\\+", x[xi])) {
        gxi <- strsplit(x[xi], ",")[[1]][1]
      } else {
        gxi <- NULL
      }
      x_c <- c(x_c, gxi)
    }
    x <- x_c
  }
  x <- sub(".*\\|", "", x)
  x <- unique(unlist(strsplit(x, ":")) )
  return(x)
}



#' Get random effect formula arguments
#'
#' Internal helper to extract random effect formula arguments.
#'
#' @param x A character string of a random effect formula.
#' @param gsubit A character string indicating the split location.
#'
#' @return A list of character strings.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get_x_random2_new <- function(x, gsubit = NULL) {
  x <- gsub("[[:space:]]", "", x)
  if(is.null(gsubit)) {
    x <- strsplit(x, ")+" )[[1]]
  } else {
    x <- strsplit(x, gsubit, fixed = T )[[1]]
  }
  x_c <- c()
  for (xi in x) {
    if(grepl("|", xi, fixed = T)) {
      zx <- strsplit(xi, "|", fixed = T)[[1]]
      zx <- zx[length(zx)]
      zx <- sub(".*gr", "", zx)
      zx <- gsub("[()]", "", zx)
      zx <- strsplit(zx, ",")[[1]][1]
    } else {
      zx <- NULL
    }
    x_c <- c(x_c, zx)
  }
  x_c <- sub(".*\\|", "", x_c)
  return(x_c)
}



#' Get random effect formula arguments with tilde sign
#'
#' Internal helper to extract random effect formula arguments that include
#' the tilde sign.
#'
#' @param x A character string of a random effect formula.
#'
#' @return A list of character strings.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get_x_random2_asitis <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
   if(any(grepl("^|gr", x)) | !any(grepl("^|gr", x))  ) {
    x <- sub(".*gr", "", x)
    x <- strsplit(x, ",")[[1]][1]
  }
  x <- sub(".*\\|", "", x)
  return(x)
}



#' Get object enclosed within parentheses
#'
#' Internal helper to extract the object enclosed within parentheses from a
#' character string.
#'
#' @param x A character string.
#'
#' @return A list of character strings.
#'
#' @inherit berkeley author
#'
#' @keywords internal
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
  return(x)
}


#' Get object enclosed within parentheses (without the parentheses)
#'
#' Internal helper to extract the object enclosed within parentheses from a
#' character string, returning the content without the surrounding parentheses.
#'
#' @param x A character string.
#'
#' @return A list of character strings.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get_o_paranthesis2 <- function(x) {
  x <- gsub("^\\(", "", x)
  x <- gsub(")$", "", x)
  x <- strsplit(x, "~")[[1]][2]
  return(x)
}


#' Get covariates from the formula
#'
#' Internal helper to extract covariates from a formula.
#'
#' @param x A character string.
#'
#' @return A character vector of covariate names.
#'
#' @inherit berkeley author
#'
#' @keywords internal
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


#' Parse and evaluate a character string
#'
#' Internal helper to parse and evaluate a character string in a given
#' environment.
#'
#' @param x A character string.
#' @param envir An environment for call evaluation.
#'
#' @return An evaluated object.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
ept <- function(x, envir = NULL) {
  if(is.null(envir)) {
    envir <- parent.frame()
  }
  eval(parse(text = x), envir = envir)
}




#' Get parameter names from the Stan code
#'
#' Internal helper to extract parameter names from a Stan code string.
#'
#' @param code A character string containing Stan code.
#' @param full Logical (default \code{TRUE}) indicating whether to return full
#'   parameter names.
#' @param section A character string specifying the Stan block
#'   (default \code{"parameters"}).
#' @param what A character string specifying the name of a particular parameter.
#'
#' @return A list of character strings.
#'
#' @inherit berkeley author
#'
#' @keywords internal
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



#' Get or set the number of cores
#'
#' Internal helper to get or set the number of cores for parallel computation.
#'
#' @param cores.arg A character string specifying the cores argument from the
#'   function.
#'
#' @return A list of integers.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get.cores <- function(cores.arg) {
  cores_ <- eval(cores.arg, envir = parent.frame())
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
  return(list(max.cores = max.cores, .cores_ps = .cores_ps))
}



#' Set up future arguments
#'
#' Internal helper to set up arguments for the \pkg{future} package.
#'
#' @param future Logical indicating whether to use \pkg{future}.
#' @param future_session A character string specifying the future session.
#' @param oldfutureplan A character string specifying the previous future plan.
#' @param setincores An integer specifying the number of cores.
#' @param verbose Logical indicating whether to print verbose output.
#'
#' @return A list.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get_future_plan_args <- function(future, 
                                 future_session, 
                                 oldfutureplan,
                                 setincores = 1,
                                 verbose = FALSE) {
  if(!future) return(NULL)
  getfutureplan <- oldfutureplan
  if(is.null(future_session)) {
    future_session      <- "sequential"
    future_session_list <- list()
  } else if(is.list(future_session)) {
    future_session_list <- future_session
    future_session      <- future_session_list[['future_session']]
    if(is.null(future_session)) {
      future_session <- "sequential"
    }
  } else if(!is.list(future_session)) {
    if(is.character(future_session)) {
      future_session_list <- list()
      future_session <- future_session
    } else {
      stop2c("'future_session' must be a single character or a named list")
    }
  }
  oldplanin <- attr(getfutureplan, "call")
  if(grepl("mirai_", future_session)) {
    insight::check_if_installed('mirai')
    if(!grepl("future.mirai::", future_session)) {
      future_session <- paste0("future.mirai::", future_session)
    }
  } else if(!grepl("mirai_", future_session)) {
    if(!grepl("future.mirai::", future_session)) {
      future_session <- paste0("future::", future_session)
    }
  }
  if(grepl("sequential", future_session)) {
    setplanis <- "sequential"
  } else if(grepl("multisession", future_session)) {
    setplanis <- "multisession"
  } else if(grepl("multicore", future_session)) {
    setplanis <- "multicore"
  } else if(grepl("cluster", future_session)) {
    setplanis <- "cluster"
  }
  if (inherits(getfutureplan, "sequential")) {
    mirai_daemons_args <- list()
    future_plan_args <- list()
    if(grepl("mirai_", future_session)) {
      if(grepl("mirai_cluster", future_session)) {
        if(is.null(future_session_list[['daemons']])) {
          #
        } else {
          if(!is.list(future_session_list[['daemons']])) {
            stop2c("The daemons in future_session list must be a list")
          }
          mirai_daemons_args <- future_session_list[['daemons']]
        }
        future_plan_args[['strategy']] <- future_session
      } else if(!grepl("mirai_cluster", future_session)) {
        future_plan_args[['workers']]  <- setincores
        future_plan_args[['strategy']] <- future_session
      } 
      mirai_daemons_args[['n']]      <- setincores
    } else if(!grepl("mirai_", future_session)) {
      if(grepl("cluster", future_session)) {
        # workers should ne n1, n2..
        future_plan_args[['workers']]  <- setincores
        future_plan_args[['strategy']] <- future_session
      } else {
        future_plan_args[['workers']]  <- setincores
        future_plan_args[['strategy']] <- future_session
      } 
    }
    if(!is_emptyx(mirai_daemons_args)) {
      do.call(mirai::daemons, mirai_daemons_args)
      if(verbose) {
        message2c("Setting mirai daemons: ", mirai::status())
      }
    }
    if(verbose) {
      message2c("The existing future plan: ", oldplanin, 
                " updated as ", future_session)
    }
  } else if (!inherits(getfutureplan, "sequential")) {
    future_plan_args <- NULL
    if(verbose) {
      message2c("Using the existing future plan: ", oldplanin)
    }
  }
  if(setplanis == "sequential") {
    future_plan_args[['workers']] <- NULL
  }
  out <- list(future_plan_args = future_plan_args, setplanis = setplanis)
  return(out)
} 


#' Validate the response variable
#'
#' Internal helper to validate the response variable in a fitted model.
#'
#' @param model An object of class \code{bgmfit}.
#' @param resp A character string specifying the name of the response variable.
#'   Default \code{NULL}.
#'
#' @return An error if evaluation fails.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
validate_response <- function(model,
                              resp = NULL) {
  uvarby <- model$model_info$univariate_by$by
  if(is.null(uvarby)) uvarby <- NA
  if (model$model_info$nys == 1 & !is.null(resp)) {
    stop_msg <- 
      paste0(
      "You have fit a univariate model",
      " but set resp option as ",
      resp,
      ".",
      "\n ",
      " The resp option should be appropriately set to NULL",
      "\n ",
      " (i.e., resp = NULL)")
    stop2c(clean_text_spaces(stop_msg))
  }
  if (model$model_info$nys > 1 & is.null(resp)) {
    if (!is.na(uvarby)) {
      stop_msg <- 
        paste0(
        "You have fit a univariate-by-subset model for '",
        uvarby, "'",
        "\n ",
        " However, the 'resp' options is not set appropriately",
        " (which is NULL at present).",
        "\n ",
        " The response options are ",
        collapse_comma(model$model_info$yvars))
      stop2c(clean_text_spaces(stop_msg))
    }
    if (model$model_info$multivariate$mvar) {
      stop_msg <- 
        paste0("You have fit a multivariate model. However, the 'resp' option",
        "\n ",
        " is not set appropriately (which is NULL at present).",
        "\n ",
        " The available response options are: ",
        collapse_comma(model$model_info$yvars))
      stop2c(clean_text_spaces(stop_msg))
    }
  }
  if (!is.null(resp)) {
    if (!resp %in% model$model_info[['yvars']]) {
      stop_msg <- 
        paste0(
        "Response should be one of the following: ",
        paste(model$model_info[['yvars']], collapse = " "),
        "\n ",
        " but you have specified: ",
        resp)
      stop2c(clean_text_spaces(stop_msg))
    }
  }
}





#' Set up priors for models with 3 or more hierarchy levels
#'
#' Internal helper to set up priors when fitting a model with three or more
#' levels of hierarchy.
#'
#' @param new_prior_list A prior object.
#'
#' @return A prior object.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
setup_higher_priors <- function(new_prior_list) {
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
        dpar = dpar_i
      )
    } 
  }
  o_l %>%  CustomDoCall(rbind, .)
}



#' Rename patterns in a character vector
#'
#' Internal helper to rename patterns in a character vector. This is adapted
#' from the \pkg{brms} package.
#'
#' @param x A character vector to be renamed.
#' @param pattern The regular expressions in \code{x} to be replaced.
#' @param replacement The replacements.
#' @param fixed Same as for \code{gsub}.
#' @param check_dup Logical indicating whether to check for duplications in
#'   \code{x} after renaming.
#' @param ... Arguments passed to \code{gsub}.
#'
#' @return A renamed character vector of the same length as \code{x}.
#'
#' @inherit berkeley author
#'
#' @keywords internal
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

#' Get call levels
#'
#' Internal helper to extract call levels from a system call status.
#'
#' @param scallstatus A system call object from \code{sys.status()}.
#' @param xclass A character string (default \code{NULL}) indicating the
#'   S3 method class.
#'
#' @return A language object.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get_xcall_byclass <- function(scallstatus, xclass = NULL) {
  if(is.null(xclass)) xclass <- '.bgmfit' 
  for (i in 1:length(scallstatus)) {
    scall <- scallstatus[[i]]
    scall <- gsub_space(paste(deparse(scall), collapse = ""))
    xcall <- NULL
    if(any(grepl(xclass, scall, fixed = F))) {
      xcall <- scall
      break
    } 
    xcall
  }
  return(xcall)
} 


#' Get call levels
#'
#' Internal helper to extract call levels from a system call.
#'
#' @param xcall A character string setting the first calling function.
#' @param scall A system call object from \code{sys.calls()}.
#' @param xstr A character string.
#' @param xclass A character string (default \code{NULL}) indicating the
#'   S3 method class.
#'
#' @return A language object.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
get_xcall__ <- function(xcall, scall, xstr, xclass = NULL) {
  scall <- scall[[length(scall)]]
  if(is.null(xclass)) xclass <- 'bgmfit'
  for (i in xstr) {
    i.xclass <- paste0(i, ".", xclass)
    xstr_x <- NULL
    xstr_FALSE <- FALSE
    if(any(grepl(i, scall, fixed = T)) |
       any(grepl(i.xclass, scall, fixed = T))) {
      xstr_x <- i
      xstr_FALSE <- TRUE
      break
    } 
    if(xstr_FALSE) {
      xcall <- xstr_x
    } else {
      xcall <- xcall
    }
    xcall
  }
  return(xcall)
} 





#' Convert first letter to upper case
#'
#' Internal helper to convert the first letter of a character string to upper
#' case.
#'
#' @param x A character string.
#'
#' @return A character string with the first letter in upper case.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}



#' Split vector at factor indices
#'
#' Internal helper to split a vector at specified factor indices.
#'
#' @param x A vector.
#' @param pos A vector of indices.
#'
#' @return A list of vector segments.
#'
#' @inherit berkeley author
#'
#' @keywords internal
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



#' Negate R's `%in%` operator in a function
#'
#' Internal helper to create a negated version of R's `%in%` operator.
#'
#' @param `%in%` The R `%in%` operator.
#'
#' @return An R function that implements the negated operator.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
`%!in%` <- Negate(`%in%`)


#' Evaluate NULL and length-zero arguments
#'
#' Internal helper to check whether arguments are \code{NULL} or have length
#' zero.
#'
#' @param x A symbol (argument).
#' @param y A symbol (argument).
#'
#' @return An R function.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
'%||%' <- function(x, y) {
  if (is.null(x)) x <- y
  x
}


#' Customize R's \code{stop} function
#'
#' Internal helper to customize R's \code{stop} function for error handling.
#'
#' @param ... Arguments passed to the customized \code{stop} function.
#'
#' @return A character string (error message) from R's \code{warning2c()}
#'   function.
#'
#' @inherit berkeley author
#'
#' @keywords internal
#' @noRd
#' 
stop2 <- function(...) {
  stop(..., call. = FALSE)
}








#' An internal function to customize R's stop function.
#' 
#' @details
#' This a wrapper around the stop function but it cleans the message for
#' excessive white spaces before parsing
#' 
#' @param ... An argument
#' @param call. A logical indicating if the call should become part of the error
#'   message.
#'   @inherit berkeley author
#' @keywords internal
#' @return A string (error message) from R's warning2c() function.
#' @noRd
#'
stop2c <- function(..., 
                   call. = FALSE, 
                   pad_before = NULL,
                   pad_after = NULL) {
  msg <- paste0(list(...), collapse = "")
  msg <- clean_text_spaces(msg)
  msg <- paste0(pad_before, " ", msg, " ", pad_after)
  stop(msg, call. = call.)
}



#' An internal function to customize R's warning function
#' @param ... An argument
#' @keywords internal
#' @return A string (warning message) from R's warning2c() function.
#' @noRd
#'
warning2 <- function(...) {
  warning(..., call. = FALSE)
}

#' An internal function to customize R's warning function.
#' 
#' @details
#' This a wrapper around the warning function but it cleans the message for
#' excessive white spaces before parsing
#' 
#' @param ... An argument
#' @param call. A logical indicating if the call should become part of the error
#'   message.
#' @param immediate. A logical
#' @param noBreaks. A logical
#' @param domain A logical
#' @keywords internal
#' @return A string (error message) from R's warning2c() function.
#' @noRd
#'
warning2c <- function(..., 
                      call = FALSE, 
                      immediate. = FALSE, 
                      noBreaks. = FALSE, 
                      domain = NULL,
                      pad_before = NULL,
                      pad_after = NULL) {
  msg <- paste0(list(...), collapse = "")
  msg <- clean_text_spaces(msg)
  msg <- paste0(pad_before, " ", msg, " ", pad_after)
  warning(msg, 
          call. = call, 
          immediate. = immediate., 
          noBreaks. = noBreaks.,
          domain = domain)
}




#' An internal function to customize R's message function
#' @param ... An argument
#' @keywords internal
#' @return A string (warning message) from R's message() function.
#' @noRd
#'
message2 <- function(...) {
  warning(..., call. = FALSE)
}

#' An internal function to customize R's message function.
#' 
#' @details
#' This a wrapper around the message function but it cleans the message for
#' excessive white spaces before parsing
#' 
#' @param ... An argument
#' @param domain A logical indicating if the call should become part of the
#'   error message.
#' @param appendLF A logical
#' @keywords internal
#' @return A string (error message) from R's warning2c() function.
#' @noRd
#'
message2c <- function(..., 
                      domain = NULL, 
                      appendLF = TRUE,
                      pad_before = NULL,
                      pad_after = NULL) {
  msg <- paste0(list(...), collapse = "")
  msg <- clean_text_spaces(msg)
  msg <- paste0(pad_before, " ", msg, " ", pad_after)
  message(msg, domain = domain, appendLF = appendLF)
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
#' @param cp_via A string character of stan code
#' @keywords internal
#' @return A character string.
#' @noRd
#'
edit_scode_ncp_to_cp_new <- function(stancode,
                                     genq_only = FALSE,
                                     normalize = TRUE, 
                                     cp_via = "multi_normal_cholesky_lpdf") {
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
              # new 03.06.2025
              # what_p <- paste0("// ", clines_pi)
              what_p <- gsub(clines_pi, "", what_p, fixed = T)
            } else {
              # what_p <- paste0("// ", clines_pi)
              what_p <- gsub(clines_pi, "", what_p, fixed = T)
            }
          }
        }
      }
    } else {
      what_p <- paste0("", clines_pi)
    }
    prepare_p <- c(prepare_p, what_p)
  }
  
  # For Nby_
  remove_duplicate_decls <- function(stan_text) {
    # ensure a character vector of lines
    lines <- if (length(stan_text) == 1L) 
      strsplit(stan_text, "\\r?\\n")[[1L]] 
    else 
      stan_text
    lines_trim <- trimws(lines)
    # identify declaration lines matching the pattern like:
    # array[Nby_1] cholesky_factor_corr[M_1] L_1;
    pat <- "^array\\s*\\[.+\\]\\s*cholesky_factor_corr\\s*\\[.+\\]\\s*\\S+\\s*;\\s*$"
    is_decl <- grepl(pat, lines_trim)
    # for declaration lines, detect duplicates based on the full trimmed line
    decl_lines <- lines_trim[is_decl]
    dup_flags <- duplicated(decl_lines)
    # keep first occurrence => mark which original lines to keep
    keep <- rep(TRUE, length(lines))
    # among declaration lines, mark duplicates for removal
    if (length(decl_lines) > 0L) {
      decl_idx <- which(is_decl)
      keep[decl_idx[dup_flags]] <- FALSE
    }
    cleaned_lines <- lines[keep]
    removed_lines <- if (length(decl_lines) > 0L) 
      decl_lines[dup_flags] 
    else 
      character(0)
    # return cleaned text (single string) and removed lines
    list(
      clean = paste(cleaned_lines, collapse = "\n"),
      removed = removed_lines,
      details = list(
        total_lines = length(lines),
        decl_count = length(decl_lines),
        removed_count = length(removed_lines)
      )
    )
  }

  move_to_p <- paste(move_to_p, collapse = "\n")
  prepare_p <- paste(prepare_p, collapse = "\n")
  prepare_p <- paste0(prepare_p, "\n", move_to_p)
  prepare_p <- strsplit(prepare_p, "\n", fixed = T)[[1]] %>%
    data.frame() %>%
    dplyr::distinct() %>%
    unlist() %>%
    as.vector()
  prepare_p <- paste(prepare_p, collapse = "\n")
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
  zz_c <- c()
  for (iz in move_to_m) {
    zz_c <- c(zz_c, paste0("  ", iz))
  }
  move_to_m <- paste(zz_c, collapse = '\n')
  zz <- strsplit(prepare_p, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    zz_c <- c(zz_c, paste0("  ", zz[iz]))
  }
  prepare_p <- paste(zz_c, collapse = '\n')
  lprior_target <- "lprior"
  m_n_c_l_c <- c()
  for (h1i in 1:length(pattern_r)) {
    pattern_ri  <- pattern_r[h1i]
    pattern_Ni  <- pattern_N[h1i]
    pattern_Mi  <- pattern_M[h1i]
    pattern_Li  <- pattern_L[h1i]
    pattern_sdi <- pattern_sd[h1i]
    # For Nby_
    check_Nby_call <- FALSE
    check_Nby_call_str <- paste0("array[Nby_", h1i, "]")
    if(grepl(check_Nby_call_str, prepare_p, fixed = T)) {
      check_Nby_call <- TRUE
    }
    if(check_Nby_call) {
      if(cp_via == "multi_normal_cholesky_lpdf") {
        make_Nby_str <- paste0("Jby_", h1i, "[i]")
        pattern_sdi <- paste0(pattern_sdi, "[,", make_Nby_str, "]")
        pattern_Li <- paste0(pattern_Li, "[", make_Nby_str, "]")
      }
      if(cp_via == "multi_normal_lpdf") {
        make_Nby_str <- paste0("Jby_", h1i, "[i]")
        pattern_sdi <- paste0(pattern_sdi, "[,", make_Nby_str, "]")
        pattern_Li <- paste0(pattern_Li, "[", make_Nby_str, "]")
      }
    }
    
    if(cp_via == "multi_normal_cholesky_lpdf") {
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
    } 
    if(cp_via == "multi_normal_lpdf") {
      m_n_c_l <-
        paste0("  for(i in 1:", pattern_Ni, ') {\n',
               "    ", lprior_target, " +=  multi_normal_lpdf(",
               pattern_ri, "[i, ] |\n",
               "    rep_row_vector(0, ",
               pattern_Mi, "),\n",
               paste0("    quad_form_diag(multiply_lower_tri_self_transpose(", 
                      pattern_Li, "), ", pattern_sdi, "));"),
               "  \n  }"
        )
    } 
    m_n_c_l_c <- c(m_n_c_l_c, m_n_c_l)
  } 
  
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
  p_block_syb_by <- paste0("", tempt_name_p, " {")
  p_block_syb_it <- paste0(p_block_syb_by, "\n", prepare_p)
  editedcode2 <- gsub(paste0("", p_block_syb_by), p_block_syb_it,
                      editedcode2, fixed=T, perl=F)
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
      if(how_many_r_1 > 0) {
        if(grepl("to_vector(z_", zz_in, fixed = T))
          zz_in <- gsub(zz_in, "", zz_in, fixed = T)
        if(grepl("scale_r_cor(z_", zz_in, fixed = T))
          zz_in <- gsub(zz_in, "", zz_in, fixed = T)
      }
      zz_c <- c(zz_c, zz_in)
    }
  }
  editedcode2 <- paste(zz_c, collapse = '\n')
  add_to_model_block <- paste0(m_n_c_l_c, "\n", move_to_m)
  add_to_genq_block <- paste0( move_to_m)
  lprior_code <- "real lprior = 0;"
  editedcode2 <- gsub(lprior_code, paste0(lprior_code, "\n",
                                          add_to_model_block, "\n"),
                      editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_tp, true_name_tp, editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_p,  true_name_p,  editedcode2, fixed = T)
  if(!normalize) {
    editedcode2 <- gsub("_lpdf", "_lupdf", editedcode2, fixed = T)
  }
  
  # For Nby_
  editedcode2 <- remove_duplicate_decls(editedcode2)$clean
  
  if(identical(pattern_r, character(0))) {
    return(stancode)
  } else if(!identical(pattern_r, character(0))) {
    if(genq_only) return(add_to_genq_block)
    return(editedcode2)
  }
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
                                 normalize = TRUE, 
                                 cp_via = "multi_normal_cholesky_lpdf") {
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
  prepare_p <- strsplit(prepare_p, "\n", fixed = T)[[1]] %>%
    data.frame() %>%
    dplyr::distinct() %>%
    unlist() %>%
    as.vector()
  prepare_p <- paste(prepare_p, collapse = "\n")
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
  zz_c <- c()
  for (iz in move_to_m) {
    zz_c <- c(zz_c, paste0("  ", iz))
  }
  move_to_m <- paste(zz_c, collapse = '\n')
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
    if(cp_via == "multi_normal_cholesky_lpdf") {
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
    } 
    if(cp_via == "multi_normal_lpdf") {
      m_n_c_l <-
        paste0("  for(i in 1:", pattern_Ni, ') {\n',
               "    ", lprior_target, " +=  multi_normal_lpdf(",
               pattern_ri, "[i, ] |\n",
               "    rep_row_vector(0, ",
               pattern_Mi, "),\n",
               paste0("    quad_form_diag(multiply_lower_tri_self_transpose(", 
                      pattern_Li, "), ", pattern_sdi, "));"),
               # "    diag_pre_multiply(",
               # pattern_sdi, ", ", pattern_Li, "));",
               "  \n  }"
        )
    } 
    m_n_c_l_c <- c(m_n_c_l_c, m_n_c_l)
  } 
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
  p_block_syb_by <- paste0("", tempt_name_p, " {")
  p_block_syb_it <- paste0(p_block_syb_by, "\n", prepare_p)
  editedcode2 <- gsub(paste0("", p_block_syb_by), p_block_syb_it,
                      editedcode2, fixed=T, perl=F)
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
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
  if(!normalize) {
    editedcode2 <- gsub("_lpdf", "_lupdf", editedcode2, fixed = T)
  }
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
#' @param summary A logical to indicate whether to summarize the posterior
#'   draws.
#' @param robust If FALSE (the default) the mean is used as the measure of
#' central tendency and the standard deviation as the measure of variability.
#' If TRUE, the median and the median absolute deviation (MAD) are applied
#' instead. Only used if summary is TRUE.
#' @param dpar A character string
#' @param verbose Print relevant information 
#' @keywords internal
#' @return A matrix
#' @noRd
#'
mapderivqr <- function(model,
                       y0,
                       xvar = NULL,
                       difx = NULL,
                       idvar = NULL,
                       levels_id = NULL,
                       newdata = NULL,
                       deriv = 1,
                       resp = NULL,
                       probs = c(0.025, 0.975),
                       summary = TRUE,
                       robust = FALSE,
                       dpar = NULL,
                       itransform = NULL,
                       cov = NULL,
                       verbose = FALSE) {

  if(is.null(probs)) {
    probs  <- c(0.025, 0.975)
  }
  if(is.null(robust)) {
    robust <- FALSE
  }

  if(is.null(newdata)) {
    newdata <- model$data
  }
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  if(is.null(dpar)) {
    dpar <- "mu"
  } else {
    dpar <- dpar
  }
  
  if(is.null(difx)) difx <- xvar
  
  validate_response(model, resp)

  list_c <- list()
  xvar_      <- paste0('xvar', resp_rev_)
  sigmaxvar_ <- paste0('sigma', xvar_)
  cov_       <- paste0('cov', resp_rev_)
  sigmacov_  <- paste0('sigma', cov_)
  uvarby     <- model$model_info$univariate_by$by
  
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
      xvar   <- model$model_info[[xvar_]]
    }
    if(is.null(cov)) cov <- model$model_info[[sigmacov_]] else cov <- cov
  } 
  
  yvar_ <- paste0('yvar', resp_rev_)
  yvar  <- model$model_info[[yvar_]]
  groupvar_ <- paste0('groupvar', resp_rev_)
  hierarchical_ <- paste0('hierarchical', resp_rev_)

  if(is.null(difx)) {
    xvar <- xvar
  } else if(!is.null(difx)) {
    xvar <- difx
    if(verbose) {
      message2c("The 'difx' is set as 'xvar' for dpar = ",
              collapse_comma(dpar))
    }
  }
  
  if(is.null(levels_id) & is.null(idvar)) {
    idvar <- model$model_info[[groupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      idvar <- model$model_info[[hierarchical_]]
    }
  } else if (!is.null(levels_id)) {
    idvar <- levels_id
  } else if (!is.null(idvar)) {
    idvar <- idvar
  }
  
  if(length(idvar) > 1) {
    idvar <- idvar[1]
  }

  if(!is.na(uvarby)) {
    newdata <- newdata %>% data.frame() %>%
      dplyr::filter(!!as.symbol(uvarby) == yvar)
  }

  if(!is.factor(newdata[[idvar]])) {
    newdata[[idvar]] <- as.factor(newdata[[idvar]])
    if(verbose) {
      message2c("The ", idvar, 
              " used in 'mapderivqr' has been converted to 'as.factor()'")
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
  newdata[[xvar]] <- ifunx_(newdata[[xvar]])
  getdydx <- function (x, y, id, data, ndigit = 2) {
    sorder <- NULL;
    data$sorder <- as.numeric(row.names(data))
    # data <- data %>% dplyr::mutate(sorder = dplyr::row_number())
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
    dydx <- CustomDoCall(rbind, dydx) %>% data.frame() %>%
      dplyr::arrange(sorder)
    return(round(dydx[[".v"]], ndigit))
  }

  mapderiv <- function(.xrow, x = xvar, y = yvar, id = idvar,
                       data = newdata) {
    newdata[[y]] <- .xrow
    getdydx(x = x, y = y, id = id, data = newdata)
  }
  
  if(is.symbol(y0)) {
    y0 <-  newdata[[deparse(y0)]]
    y0 <- as.matrix(y0) %>% t()
  } else if(is.character(y0)) {
    y0 <- newdata[[y0]]
    y0 <- as.matrix(y0) %>% t()
  } else if(is.vector(y0)) {
    y0 <- as.matrix(y0) %>% t()
  } else if(is.matrix(y0)) {
    y0 <- y0
  }

  if(deriv == 1) {
    tempx <- apply(y0, 1, mapderiv) %>% t()
  }
  if(deriv == 2) {
    tempx <- apply(y0, 1, mapderiv) %>% t()
    tempx <- apply(tempx , 1, mapderiv) %>% t()
  }
  
  if(summary) {
    dout <-  brms::posterior_summary(tempx , probs = probs, robust = robust)
  } else {
    dout <- tempx
  }
  
  if(all(is.infinite(dout))) {
    stop2c("The 'mapderivqr()' resulted in all infinite values.",
         "\n  ", 
         "This could be because of an ncorrect xvar used.",
         "\n  ", 
         "The currect xvar used in mapderivqr() is: ", 
         collapse_comma(xvar))
  }
  return(dout)
}


#' mapderivqr_standalone
#' @noRd
mapderivqr_standalone <- function(xvar, 
                                  yvar, 
                                  idvar, 
                                  newdata, 
                                  deriv = 1,
                                  summary = TRUE,
                                  robust = FALSE,
                                  probs = c(0.025, 0.975),
                                  ndigit = 2) {
  y0 <- yvar
  getdydx <- function(x, y, id, data, ndigit = 2) {
    sorder <- NULL
    # data$sorder <- as.numeric(row.names(data))
    data <- data %>% dplyr::mutate(sorder = dplyr::row_number())
    .data <- data %>% dplyr::mutate(.x = !!dplyr::sym(x)) %>% 
      dplyr::mutate(.y = !!dplyr::sym(y)) %>% dplyr::mutate(.id = !!dplyr::sym(id)) %>% 
      data.frame()
    .dydx <- function(x, y) {
      n <- length(x)
      i1 <- 1:2
      i2 <- (n - 1):n
      c(diff(y[i1])/diff(x[i1]), (y[-i1] - y[-i2])/(x[-i1] - 
                                                      x[-i2]), diff(y[i2])/diff(x[i2]))
    }
    dydx <- lapply(split(.data, as.numeric(.data$.id)), function(x) {
      x$.v <- .dydx(x$.x, x$.y)
      x
    })
    dydx <- CustomDoCall(rbind, dydx) %>% data.frame() %>% 
      dplyr::arrange(sorder)
    return(round(dydx[[".v"]], ndigit))
  }
  mapderiv <- function(.xrow, x = xvar, y = yvar, id = idvar, 
                       data = newdata) {
    newdata[[y]] <- .xrow
    getdydx(x = x, y = y, id = id, data = newdata)
  }
  if (is.symbol(y0)) {
    y0 <- newdata[[deparse(y0)]]
    y0 <- as.matrix(y0) %>% t()
  }
  else if (is.character(y0)) {
    y0 <- newdata[[y0]]
    y0 <- as.matrix(y0) %>% t()
  }
  else if (is.vector(y0)) {
    y0 <- as.matrix(y0) %>% t()
  }
  else if (is.matrix(y0)) {
    y0 <- y0
  }
  if (deriv == 1) {
    tempx <- apply(y0, 1, mapderiv) %>% t()
  }
  if (deriv == 2) {
    tempx <- apply(y0, 1, mapderiv) %>% t()
    tempx <- apply(tempx, 1, mapderiv) %>% t()
  }
  if (summary) {
    dout <- brms::posterior_summary(tempx, probs = probs, 
                                    robust = robust)
  }
  else {
    dout <- tempx
  }
  if (all(is.infinite(dout))) {
    stop2c("The 'mapderivqr()' resulted in all infinite values.", 
           "\n  ", "This could be because of an ncorrect xvar used.", 
           "\n  ", "The currect xvar used in mapderivqr() is: ", 
           collapse_comma(xvar))
  }
  return(dout)
}



#' Compute first derivatives from a distance curve
#'
#' Internal utility to estimate the first derivative (velocity curve) from a
#' distance curve using a range of numerical differentiation methods.
#'
#' @details The function expects \code{y} (response) before \code{x} (predictor)
#' to align with conventions used in downstream workflows (e.g., marginal
#' effects).
#'
#' Six numerical differentiation methods are currently implemented:
#' \itemize{
#'   \item Method 1: Forward finite differences (direct method)
#'   \item Method 2: Smoothing spline followed by analytical derivative
#'   \item Method 3: Centered finite differences
#'   \item Method 4: Higher-order (5-point stencil) finite differences
#'   \item Method 5: Richardson extrapolation
#'   \item Method 6: Hybrid finite differences (similar to \code{.dydx})
#' }
#'
#' @param y Numeric vector representing the distance curve (response variable).
#' @param x Numeric vector representing the predictor variable (e.g., time or
#'   age).
#' @param method Integer specifying the differentiation method (1--6). Default
#'   is 1.
#' @param length.out Integer specifying the number of evaluation points for the
#'   smoothing spline method (\code{method = 2}). Default is \code{NULL}, which
#'   uses 100 points.
#' @param df Numeric specifying the degrees of freedom for the smoothing spline
#'   (\code{method = 2}). Default is \code{NULL}, in which case \code{0.3 *
#'   length(x)} is used.
#' @param verbose Logical; if \code{TRUE}, prints diagnostic information.
#'   Currently unused.
#'
#' @return A data frame with two columns:
#' \itemize{
#'   \item \code{x}: locations at which the derivative is evaluated
#'   \item \code{y}: estimated first derivative values
#' }
#'
#' @examples
#' # Empirical sigmoid-like data
#' x_data <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' y_data <- c(0.1, 0.15, 0.25, 0.45, 0.7, 0.85, 0.92, 0.96, 0.98, 0.99, 0.995)
#'
#' ipts <- 10
#' df <- 10
#'
#' plot(
#'   get_d1_from_d0(y_data, x_data, method = 1)$x,
#'   get_d1_from_d0(y_data, x_data, method = 1)$y,
#'   xlab = "x",
#'   ylab = "dy/dx"
#' )
#'
#' lines(
#'   get_d1_from_d0(y_data, x_data, method = 1, length.out = ipts, df = df)$x,
#'   get_d1_from_d0(y_data, x_data, method = 1, length.out = ipts, df = df)$y
#' )
#' lines(
#'   get_d1_from_d0(y_data, x_data, method = 2, length.out = ipts, df = df)$x,
#'   get_d1_from_d0(y_data, x_data, method = 2, length.out = ipts, df = df)$y
#' )
#' lines(
#'   get_d1_from_d0(y_data, x_data, method = 3, length.out = ipts, df = df)$x,
#'   get_d1_from_d0(y_data, x_data, method = 3, length.out = ipts, df = df)$y
#' )
#' lines(
#'   get_d1_from_d0(y_data, x_data, method = 4, length.out = ipts, df = df)$x,
#'   get_d1_from_d0(y_data, x_data, method = 4, length.out = ipts, df = df)$y
#' )
#' lines(
#'   get_d1_from_d0(y_data, x_data, method = 5, length.out = ipts, df = df)$x,
#'   get_d1_from_d0(y_data, x_data, method = 5, length.out = ipts, df = df)$y
#' )
#' lines(
#'   get_d1_from_d0(y_data, x_data, method = 6, length.out = ipts, df = df)$x,
#'   get_d1_from_d0(y_data, x_data, method = 6, length.out = ipts, df = df)$y
#' )
#'
#' legend(
#'   "topright",
#'   legend = paste("method", 1:6),
#'   col = 1:6,
#'   lty = 1,
#'   bty = "n"
#' )
#' @keywords internal
#' @noRd
#' 
get_d1_from_d0 <- function(y,
                           x,
                           method = 1, 
                           length.out = NULL, 
                           df = NULL, 
                           verbose = FALSE) {
  x_data <- x
  y_data <- y
  out <- list()
  if(method == 1) { # method == 1 -> direct
    direct_diff <- function(x, y) {
      sorted_idx <- order(x)
      x_sorted <- x[sorted_idx]
      y_sorted <- y[sorted_idx]
      dydx <- diff(y_sorted) / diff(x_sorted)
      x_mid <- x_sorted[-length(x_sorted)] + diff(x_sorted)/2
      return(data.frame(x = x_mid, y = dydx))
    }
    temp <- direct_diff(x_data, y_data)
    out[['x']] <- temp$x
    out[['y']] <- temp$y
  } else if(method == 2) { # method == 2 -> Smoothing Spline + Derivative
    smooth_diff <- function(x, y, df, length.out) {
      if (is.null(df))         df <- length(x) * 0.3 
      if (is.null(length.out)) length.out <- 100
      smooth_fit <-  stats::smooth.spline(x, y, df = df)
      x_fine <- seq(min(x), max(x), length.out = length.out)
      deriv_smooth <- stats::predict(smooth_fit, x_fine, deriv = 1)$y
      x_deriv_smooth <- x_fine
      return(data.frame(x = x_deriv_smooth, y = deriv_smooth))
    }
    temp <- smooth_diff(x_data, y_data, df = df, length.out = length.out)
    out[['x']] <- temp$x
    out[['y']] <- temp$y
  } else if(method == 3) { # method == 3 -> Centered difference
    # Centered difference: f'(x) ≈ (f(x+h) - f(x-h)) / 2h
    centered_diff <- function(x, y) {
      n <- length(x)
      deriv <- numeric(n-2)
      x_mid <- numeric(n-2)
      for (i in 2:(n-1)) {
        h1 <- x[i] - x[i-1]
        h2 <- x[i+1] - x[i]
        deriv[i-1] <- (y[i+1] - y[i-1]) / (h1 + h2)
        x_mid[i-1] <- x[i]
      }
      return(data.frame(x = x_mid, y = deriv))
    }
    temp <- centered_diff(x_data, y_data)
    out[['x']] <- temp$x
    out[['y']] <- temp$y
  } else if(method == 4) { #method == 4 -> Higher-Order (4th) Finite Differences
    # 5-point stencil for higher accuracy
    higher_order_diff <- function(x, y) {
      n <- length(x)
      if (n < 5) stop2c("Need at least 5 points")
      h <- mean(diff(x))  # assumes uniform spacing
      deriv <- (-y[5:n] + 8*y[4:(n-1)] - 8*y[2:(n-3)] + y[1:(n-4)]) / (12*h)
      x_mid <- x[3:(n-2)]
      return(data.frame(x = x_mid, y = deriv))
    }
    temp <- higher_order_diff(x_data, y_data)
    out[['x']] <- temp$x
    out[['y']] <- temp$y
  } else if(method == 5) { # method == 5 -> Richardson Extrapolation
    richardson_deriv <- function(x, y) {
      h <- mean(diff(x))
      deriv_h <- diff(y) / diff(x)
      x_h <- x[-length(x)] + diff(x)/2
      x_fine <- seq(min(x), max(x), by = h/2)
      y_interp <- stats::approx(x, y, x_fine)$y
      deriv_h2 <- diff(y_interp) / diff(x_fine)
      x_h2 <- x_fine[-length(x_fine)] + diff(x_fine)/2
      deriv_fine_interp <- stats::approx(x_h2, deriv_h2, x_h)$y
      richardson_est <- (4 * deriv_fine_interp - deriv_h) / 3
      return(data.frame(x = x_h, y = richardson_est))
    }
    temp <- richardson_deriv(x_data, y_data)
    out[['x']] <- temp$x
    out[['y']] <- temp$y
  } else if(method == 6) { # method == 6 -> .dydx
    .dydx_deriv <- function(x, y) {
      sorted_idx <- order(x)
      x_sorted <- x[sorted_idx]
      y_sorted <- y[sorted_idx]
      .dydx <- function(x, y) {
        n <- length(x); i1 <- 1:2; i2 <- (n - 1):n
        c(diff(y[i1])/diff(x[i1]), (y[-i1] - y[-i2])/(x[-i1] - x[-i2]),
          diff(y[i2])/diff(x[i2]))
      }
      deriv_vals <- .dydx(x_sorted, y_sorted)
      return(data.frame(x = x_sorted, y = deriv_vals))
    }
    temp <- .dydx_deriv(x_data, y_data)
    out[['x']] <- temp$x
    out[['y']] <- temp$y
  } else if(method == 7) { # method == 7 -> Functional Data Analysis (FDA)
    stop2c("'method' should be between 1 and 6")
  } else if(method == 8) { # method == 8 -> Smoothing Spline + Derivative
    stop2c("'method' should be between 1 and 7")
  }
  out_d <- data.frame(x = out$x, y = out$y)
  return(out_d)
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




#' An internal function to edit code for rescorr by group fit via cmdstanr/rstan
#'
#' @param stan_code A character string
#' @param threads An integer of \code{NULL} 
#' @param corr_method A character string
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
edit_stancode_for_multivariate_rescor_by <- function(stan_code, 
                                                     threads = NULL,
                                                     normalize = TRUE,
                                                     corr_method = 'lkj') {
  if(is.null(corr_method)) {
    via <- 'cde'
  } else if(corr_method == 'lkj') {
    via <- 'lkj'
  } else if(corr_method == 'cde') {
    via <- 'cde'
  } else {
    stop2c("corr_method must be either 'lkj' or 'cde'")
  }
  brms_code_edited <- stan_code
  gsub_it <- "cholesky_factor_corr[nresp] Lrescor;"
  gsub_by <- "array[Rescor_Nby] cholesky_factor_corr[nresp] Lrescor;"
  if(corr_method == 'lkj') {
    gsub_by <- "array[Rescor_Nby] cholesky_factor_corr[nresp] Lrescor;"
  } else if(corr_method == 'cde') {
    gsub_by <- "array[Rescor_Nby] vector<lower=-1, upper=1>[N_rhos] CoefRescor;"
  }
  brms_code_edited <- gsub(gsub_it, gsub_by, brms_code_edited, fixed = T)
  
  gsub_it <- "corr_matrix[nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor);"
  gsub_by <- "array[Rescor_Nby] corr_matrix[nresp] Rescor;
        for (c in 1:Rescor_Nby) Rescor[c] = multiply_lower_tri_self_transpose(Lrescor[c]);"
  brms_code_edited <- gsub(gsub_it, gsub_by, brms_code_edited, fixed = T)
  
  gsub_it <- "vector<lower=-1,upper=1>[nrescor] rescor;"
  gsub_by <- "array[Rescor_Nby] vector<lower=-1,upper=1>[nrescor] rescor;"
  brms_code_edited <- gsub(gsub_it, gsub_by, brms_code_edited, fixed = T)
  
  gsub_it <- "rescor[choose(k - 1, 2) + j] = Rescor[j, k];"
  gsub_by <- "for (c in 1:Rescor_Nby) rescor[c, choose(k - 1, 2) + j] = Rescor[c, j, k];"
  brms_code_edited <- gsub(gsub_it, gsub_by, brms_code_edited, fixed = T)
  
  gsub_it_start <- "lprior += lkj_corr_cholesky_lpdf(Lrescor"
  gsub_it_end   <- ");"
  if(corr_method == 'lkj') {
    gsub_by <- "for (c in 1:Rescor_Nby) lprior += lkj_corr_cholesky_lpdf(Lrescor[c] | Rescor_prior[c]);"
  } else if(corr_method == 'cde') {
    gsub_by <- "for (c in 1:Rescor_Nby) lprior += uniform_lpdf(CoefRescor[c] | -1, 1);"
  }
  if(!normalize) {
    gsub_it_start <- gsub("_lpdf", "_lupdf", gsub_it_start, fixed = T)
    gsub_by       <- gsub("_lpdf", "_lupdf", gsub_by, fixed = T)
  }
  
  brms_code_edited <- replace_string_part(x = brms_code_edited,
                                          start = gsub_it_start, 
                                          end =  gsub_it_end,
                                          replace = gsub_by)

  if(corr_method == 'lkj') {
    brms_code_edited <- brms_code_edited
  } else if(corr_method == 'cde') {
    gsub_it_start <- "transformed parameters"
    gsub_it_end   <- "{"
    gsub_by <- "transformed parameters {
          // Imp to use local because otherwise interger such as int ind = 0 not allowed in tparameters
          array[Rescor_Nby] matrix[nresp, nresp] Rrescor; 
          array[Rescor_Nby] matrix[nresp, nresp] Lrescor;
          { // start local
              // array[Rescor_Nby] vector[nresp] diag_elements_array = rep_vector(1.0, nresp);
              for (c in 1:Rescor_Nby) {
                matrix[nresp, nresp] current_Rrescor;
                for (i in 1:nresp) {
                   current_Rrescor[i, i] = 1.0;
                  // current_Rrescor[i, i] = diag_elements_array[c, i];
                }
                int ind = 1; 
                for (j in 1:(nresp-1)) { 
                  for (i in (j+1):nresp) { 
                    current_Rrescor[i, j] = CoefRescor[c, ind]; 
                    current_Rrescor[j, i] = CoefRescor[c, ind]; 
                    ind += 1;
                  }
                }
                Rrescor[c] = current_Rrescor;
                Lrescor[c] = cholesky_decompose(Rrescor[c]);
              }
          }  // end local"
    brms_code_edited <- replace_string_part(x = brms_code_edited,
                                            start = gsub_it_start, 
                                            end =  gsub_it_end,
                                            replace = gsub_by)
  } 
  gsub_it_start <- "LSigma[n] = diag_pre_multiply"
  gsub_it_end   <- ");"
  if(is.null(threads)) {
    gsub_by <- "int Index_Rescor = Rescor_by_id[Rescor_gr_id[n]];
        LSigma[n] = diag_pre_multiply(sigma[n], Lrescor[Index_Rescor]);"
  } else if(!is.null(threads)) {
    gsub_by <- "int Index_Rescor = Rescor_by_id[Rescor_gr_id[nn]];
      LSigma[n] = diag_pre_multiply(sigma[n], Lrescor[Index_Rescor]);"
  }
  brms_code_edited <- replace_string_part(x = brms_code_edited,
                                          start = gsub_it_start, 
                                          end =  gsub_it_end,
                                          replace = gsub_by)
  if(is.null(threads)) {
    brms_code_edited <- brms_code_edited
  } else if(!is.null(threads)) {
    gsub_it   <- "matrix Lrescor,"
    gsub_by   <- "array[] matrix Lrescor,"
    brms_code_edited <- gsub(gsub_it, gsub_by, brms_code_edited, fixed = TRUE)
  } 
  if(grepl("sigma[n]", brms_code_edited, fixed = T)) {
    sigma_single_parm <- FALSE
  } else if(!grepl("sigma[n]", brms_code_edited, fixed = T)) {
    sigma_single_parm <- TRUE
  } else {
    stop2c("Something wrong with sigma rescor")
  }
  if(!sigma_single_parm) {
    out_edited_code <- brms_code_edited
  } else if(sigma_single_parm) {
    gsub_it_start <- "Mu[n] = transpose"
    gsub_it_end   <- ");"
    mu_n <- replace_string_part(x = brms_code_edited,
                                start = gsub_it_start, 
                                end =  gsub_it_end,
                                replace = "",
                                extract = T)
    gsub_it_start <- "sigma ="
    gsub_it_end   <- ");"
    sigma_n <- replace_string_part(x = brms_code_edited,
                                   start = gsub_it_start, 
                                   end =  gsub_it_end,
                                   replace = "",
                                   extract = T)
    gsub_it_start <- "vector[nresp] sigma ="
    gsub_it_end   <- ");"
    brms_code_edited <- replace_string_part(x = brms_code_edited,
                                            start = gsub_it_start, 
                                            end =  gsub_it_end,
                                            replace = "array[N] vector[nresp] sigma;",
                                            extract = F)
    gsub_it_start <- "matrix[nresp, nresp] LSigma ="
    gsub_it_end   <- ");"
    brms_code_edited <- replace_string_part(x = brms_code_edited,
                                            start = gsub_it_start, 
                                            end =  gsub_it_end,
                                            replace = "array[N] matrix[nresp, nresp] LSigma;",
                                            extract = F)
    
    sigma_n <- gsub("sigma =", "sigma[n] =", sigma_n, fixed = T)
    if(is.null(threads)) {
      plus_mis <- "int Index_Rescor = Rescor_by_id[Rescor_gr_id[n]];
      LSigma[n] = diag_pre_multiply(sigma[n], Lrescor[Index_Rescor]);"
    } else if(!is.null(threads)) {
      plus_mis <- "int Index_Rescor = Rescor_by_id[Rescor_gr_id[nn]];
      LSigma[n] = diag_pre_multiply(sigma[n], Lrescor[Index_Rescor]);"
    }
    mu_n_sigma_n_plus_mis <- paste0(mu_n, "\n      ", 
                                    sigma_n, "\n      ",
                                    plus_mis)
    
    brms_code_edited <- gsub(mu_n, mu_n_sigma_n_plus_mis, 
                             brms_code_edited, fixed = T)
    if(is.null(threads)) {
      gsub_it_target_n <- "target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);"
      gsub_by_target_n <- "for (n in 1:N) {
    target += multi_normal_cholesky_lpdf(Y[n] | Mu[n], LSigma[n]);
  }"
    } else if(!is.null(threads)) {
      gsub_it_target_n <- "ptarget += multi_normal_cholesky_lpdf(Y[start:end] | Mu, LSigma);"
      gsub_by_target_n <- 
        "for (n in 1:N) {
     int nn = n + start - 1;
     ptarget += multi_normal_cholesky_lpdf(Y[nn] | Mu[n], LSigma[n]);
    }"
    }
    
    if(!normalize) {
      gsub_it_target_n <- gsub("_lpdf", "_lupdf", gsub_it_target_n, fixed = T)
      gsub_by_target_n <- gsub("_lpdf", "_lupdf", gsub_by_target_n, fixed = T)
    }
    
    brms_code_edited <- gsub(gsub_it_target_n, gsub_by_target_n, brms_code_edited, fixed = T)
    
    out_edited_code <- brms_code_edited
  } 
  return(out_edited_code)
}


#' custom_get_data.brmsfit for for \code{insight get_data}
#'
#' @param x A brms objects
#' @param effects levels of group levels - see [insight::get_data()]
#' @param component levels of group levels - see [insight::get_data()]
#' @param source levels of group levels - see [insight::get_data()]
#' @param verbose levels of group levels - see [insight::get_data()]
#' @param ... Other arguments
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
custom_get_data.brmsfit <- function (x, 
                                     effects = "all", 
                                     component = "all", 
                                     source = "environment", 
                                     verbose = FALSE, 
                                     ...) {
  
  .get_data_from_environment <-  .all_elements <-  find_variables <- 
    .is_multi_membership <- .return_combined_data <- 
    .prepare_get_data <- is_multivariate <- 
    .clean_brms_mm <- find_random_slopes <- is_empty_object <- NULL;
  getfrom_ <- c('.get_data_from_environment', 
                '.all_elements', 
                'find_variables', 
                '.is_multi_membership', 
                '.return_combined_data', 
                '.prepare_get_data', 
                'is_multivariate', 
                '.clean_brms_mm',
                'find_random_slopes',
                'is_empty_object')
  
  for (i in getfrom_) {
    assign(i, utils::getFromNamespace(i, 'insight'))
  }
  
  return_option <- 3
  clean.x       <- F
  
  if(!is.null(x$model_info$call.bgmfit$sigma_formula_manual)) {
    return_option <- 4
  }

  if(return_option == 4) {
    out <- x$data
    if(clean.x) {
      # out <- out[, -grep("\\.V\\d+$", names(out))]
      out <- out %>% dplyr::select(!dplyr::matches("\\.\\d+$"))
    }
    return(out)
  }

  data_name <- attr(x$data, "data_name")
  assign(data_name, x$data)
  x$data <- get(data_name)
  model_data <- .get_data_from_environment(x, effects = effects, 
                                           component = component, 
                                           source = source, verbose = verbose, 
                                           data_name = data_name, ...)
  if(return_option == 1) {
    if (!is.null(model_data)) {
      out <- model_data
    } else if (is.null(model_data)) {
      out <- x$data
    }
  }
  effects <- match.arg(effects, choices = c("all", "fixed", 
                                            "random"))
  component <- match.arg(component, choices = c("all", .all_elements()))
  model.terms <- find_variables(x, effects = "all", component = "all", 
                                flatten = FALSE)
  mf <- stats::model.frame(x)
  if (.is_multi_membership(x)) {
    model.terms <- lapply(model.terms, .clean_brms_mm)
    rs <- setdiff(unlist(find_random_slopes(x), use.names = FALSE), 
                  unlist(model.terms, use.names = FALSE))
    if (!is_empty_object(rs)) 
      model.terms$random <- c(rs, model.terms$random)
  }
  if(return_option == 2) {
    outmf <- .return_combined_data(x, 
                                   .prepare_get_data(x, 
                                                     mf, 
                                                     effects = effects, 
                                                     verbose = verbose), 
                                   effects, component, model.terms, 
                                   is_mv = is_multivariate(x), 
                                   verbose = verbose) 
    
    new_cols             <- setdiff(names(outmf), names(model_data))
    model_data[new_cols] <- outmf[new_cols]
    out                  <- model_data
  }
  if(return_option == 3) {
    out <- .prepare_get_data(x, mf, effects = effects, verbose = verbose)
    if(clean.x) {
      # out <- out[, -grep("\\.V\\d+$", names(out))]
      out <- out %>% dplyr::select(!dplyr::matches("\\.V\\d+$"))
    }
  }
  return(out)
}


#' custom_get_predictors.brmsfit for for \code{insight get_predictors}
#'
#' @param x A brms objects
#' @param effects levels of group levels - see [insight::get_data()]
#' @param component levels of group levels - see [insight::get_data()]
#' @param source levels of group levels - see [insight::get_data()]
#' @param verbose levels of group levels - see [insight::get_data()]
#' @param ... Other arguments
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
custom_find_predictors.brmsfit <- function (x, 
                                           effects = "fixed", 
                                           component = "all", 
                                           flatten = FALSE, 
                                           verbose = TRUE, ...) {
  
  
  .get_data_from_environment <-  .all_elements <-  find_variables <- 
    .is_multi_membership <- .return_combined_data <- 
    .prepare_get_data <- is_multivariate <- 
    .clean_brms_mm <- find_random_slopes <- is_empty_object <- NULL;
  
  validate_argument <- find_formula <- .get_elements <- NULL;
  find_auxiliary <- .brms_elements <- .prepare_predictors_brms <- NULL;
  .return_vars <- compact_list <- NULL;
  
  getfrom_ <- c('.get_data_from_environment', 'validate_argument',
                '.all_elements', 'find_formula',
                'find_variables', '.get_elements',
                '.is_multi_membership', 'find_auxiliary',
                '.return_combined_data', '.brms_elements',
                '.prepare_get_data', '.prepare_predictors_brms',
                'is_multivariate', '.return_vars',
                '.clean_brms_mm', "compact_list",
                'find_random_slopes',
                'is_empty_object')
  
  for (i in getfrom_) {
    assign(i, utils::getFromNamespace(i, 'insight'))
  }
  return_option <- 3
  clean.x       <- F
  if(!is.null(x$model_info$call.bgmfit$sigma_formula_manual)) {
    return_option <- 4
  }
  effects <- validate_argument(effects, c("fixed", "random", 
                                          "all"))
  component <- validate_argument(component, c("all", "conditional", 
                                              "zi", "zero_inflated",
                                              "dispersion", "instruments", 
                                              "correlation", 
                                              "smooth_terms", "location", 
                                              "auxiliary", "distributional"))
  f <- find_formula(x, verbose = verbose)
  is_mv <- is_multivariate(f)
  elements <- .get_elements(effects, component, model = x)
  dpars <- find_auxiliary(x)
  elements <- .brms_elements(effects, component, dpars)
  if (is_mv) {
    f <- lapply(f, function(.x) .prepare_predictors_brms(x, 
                                                         .x, elements))
  }
  else {
    f <- .prepare_predictors_brms(x, f, elements)
  }
  if (is_mv) {
    l <- lapply(f, .return_vars, x = x)
  }
  else {
    l <- .return_vars(f, x)
  }
  if (is_empty_object(l) || is_empty_object(compact_list(l))) {
    return(NULL)
  }
  if (any(endsWith(names(l), "random")) && effects == "all") {
    random_slope <- unlist(find_random_slopes(x), use.names = FALSE)
    all_predictors <- unlist(unique(l), use.names = FALSE)
    rs_not_in_pred <- unique(setdiff(random_slope, all_predictors))
    if (length(rs_not_in_pred)) 
      l$random <- c(rs_not_in_pred, l$random)
  }
  if (flatten) {
    out <- unique(unlist(l, use.names = FALSE))
  }
  else {
    out <- compact_list(l)
  }
  return(out)
}



#' custom_get_ci_draws 
#' 
#' @description
#' Useful in \code{marginaleffects method = 'pkg'} \code{get_growthparameters}
#' as NA are not allowed
#'
#' @param x See \code{marginaleffects:::get_ci_draws}
#' @param conf_level See \code{marginaleffects:::get_ci_draws}
#' @param draws See \code{marginaleffects:::get_ci_draws}
#' @param model See \code{marginaleffects:::get_ci_draws}
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
custom_get_ci_draws <- function (x, conf_level, draws, model = NULL) {
  get_eti <- get_hdi <- NULL;
  getfrom_ <- c('get_eti', 'get_hdi')
  for (i in getfrom_) {
    assign(i, utils::getFromNamespace(i, 'marginaleffects'))
  }
  checkmate::check_number(conf_level, lower = 1e-10, upper = 1 - 1e-10)
  critical <- (1 - conf_level)/2
  if (inherits(model, "inferences_simulation")) {
    insight::check_if_installed("collapse", minimum_version = "1.9.0")
    CIs <- collapse::dapply(draws, MARGIN = 1, FUN = collapse::fquantile, 
                            probs = c(critical, 1 - critical))
    x$std.error <- collapse::dapply(draws, MARGIN = 1, FUN = collapse::fsd)
    x$conf.low <- CIs[, 1]
    x$conf.high <- CIs[, 2]
    return(x)
  }
  else if (identical("eti", getOption("marginaleffects_posterior_interval", 
                                      default = "eti")) && 
           identical("median", getOption("marginaleffects_posterior_center", 
                                         default = "median"))) {
    
    insight::check_if_installed("collapse", minimum_version = "1.9.0")
    
    if(!is_emptyx(draws)) {
      if (nrow(draws) > 0) {
        CIs <- collapse::dapply(draws, MARGIN = 1, FUN = collapse::fquantile, 
                                probs = c(critical, 0.5, 1 - critical))
        x$estimate <- CIs[, 2]
        x$conf.low <- CIs[, 1]
        x$conf.high <- CIs[, 3]
      }
    } else {
        x$estimate <- NA
        x$conf.low <- NA
        x$conf.high <- NA
    }
    return(x)
  }
  if (identical("eti", getOption("marginaleffects_posterior_interval", 
                                 default = "eti")) && 
      identical("mean", getOption("marginaleffects_posterior_center", 
                                  default = "median"))) {
    
    insight::check_if_installed("collapse", minimum_version = "1.9.0")
    if(!is_emptyx(draws)) {
      Bs <- collapse::dapply(draws, MARGIN = 1, FUN = collapse::fmean)
      CIs <- collapse::dapply(draws, MARGIN = 1, FUN = collapse::fquantile,
                              probs = c(critical, 1 - critical))
      x$estimate <- Bs
      x$conf.low <- CIs[, 1]
      x$conf.high <- CIs[, 2]
    } else {
      x$estimate <- NA
      x$conf.low <- NA
      x$conf.high <- NA
    }
    return(x)
  }
  
  FUN_INTERVAL <- getOption("marginaleffects_posterior_interval")
  if (is.null(FUN_INTERVAL)) {
    FUN_INTERVAL <- getOption("marginaleffects_credible_interval", 
                              default = "eti")
  }
  checkmate::assert_choice(FUN_INTERVAL, choices = c("eti", 
                                                     "hdi"))
  if (FUN_INTERVAL == "hdi") {
    FUN_INTERVAL <- get_hdi
  }
  else {
    FUN_INTERVAL <- get_eti
  }
  FUN_CENTER <- getOption("marginaleffects_posterior_center", 
                          default = stats::median)
  checkmate::assert(checkmate::check_choice(FUN_CENTER, 
                                            choices = c("mean", "median")), 
                    checkmate::check_function(FUN_CENTER))
  if (identical(FUN_CENTER, "mean")) {
    FUN_CENTER <- mean
  }
  else if (identical(FUN_CENTER, "median")) {
    FUN_CENTER <- stats::median
  }
  colnames(draws) <- row.names(draws) <- NULL
  CIs <- t(apply(draws, 1, FUN_INTERVAL, credMass = conf_level))
  Bs <- apply(draws, 1, FUN_CENTER)
  if (nrow(x) < nrow(CIs)) {
    CIs <- unique(CIs)
    Bs <- unique(Bs)
  }
  x[["estimate"]] <- Bs
  x[["conf.low"]] <- CIs[, "lower"]
  x[["conf.high"]] <- CIs[, "upper"]
  return(x)
}




#' unlock_replace_bind for for \code{insight get_data}
#'
#' @param package A brms objects
#' @param what levels of group levels - see [insight::get_data()]
#' @param replacement levels of group levels - see [insight::get_data()]
#' @param ept_str levels of group levels - see [insight::get_data()]
#' @param verbose levels of group levels - see [insight::get_data()]
#' 
#' @return A function
#' @keywords internal
#' @noRd
#'
unlock_replace_bind <- function(package, 
                                what, 
                                replacement, 
                                ept_str = TRUE, 
                                verbose = FALSE) {
  if(!ept_str) {

  } else if(ept_str) {
    what <- collapse_comma(what)
    package <- collapse_comma(package)
    replacement_str <- collapse_comma('replacement')
    assign(replacement_str, replacement)
    step_1 <- paste0("unlockBinding(", what, ", getNamespace(", package, "))")
    step_2 <- paste0("assign(", what, ",", 
                     ept(replacement_str), ",  envir=getNamespace(", package, "))")
    step_3 <- paste0("lockBinding(", what, ", getNamespace(", package, "))")
    ept(step_1)
    ept(step_2)
    ept(step_3)
  }
}



#' custom_rename_pars for residual corr by group fit via cmdstanr
#'
#' @param x A brms objects
#' @param Rescor_by_levels levels of group levels
#' @param ... Other arguments
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
custom_rename_pars <- function (x, Rescor_by_levels = NULL, ...) {
  if (!length(x$fit@sim)) {
    return(x)
  }
  brmsframe <- rename_predictor <- rename_re <- rename_Xme <- 
    save_old_par_order <- do_renaming <- repair_stanfit <-
    compute_quantities <- reorder_pars <- is.rlist <- NULL;
  getfrom_ <- c('brmsframe', 'rename_predictor', 'rename_re', 'rename_Xme', 
                'save_old_par_order', 'do_renaming', 'repair_stanfit', 
                "compute_quantities", "reorder_pars", "is.rlist")
  
  for (i in getfrom_) {
    assign(i, utils::getFromNamespace(i, 'brms'))
    }
  rename_predictor.mvbrmsterms <- function (x, pars, 
                                            Rescor_by_levels = NULL,
                                            rescor_names = NULL, ...) {
    
    `c<-` <- `lc<-` <- get_cornames <- rlist <- NULL;
    
    getfrom_ <- c("c<-", "lc<-", 'get_cornames', 'rlist')
    for (i in getfrom_) {
      assign(i, utils::getFromNamespace(i, 'brms'))
      }
    out <- list()
    for (i in seq_along(x$terms)) {
      c(out) <- rename_predictor(x$terms[[i]], pars = pars,  ...)
    }
    if (x$rescor) {
      if(is.null(rescor_names)) {
        rescor_names <- get_cornames(x$responses, type = "rescor", 
                                     brackets = FALSE)
        if(!is.null(Rescor_by_levels)) {
          rescor_names_Rescor_by_levels <- c()
          for (i in Rescor_by_levels) {
            tx_i <- gsub("rescor", paste0("rescor", "_", i), rescor_names, 
                         fixed = T )
            rescor_names_Rescor_by_levels <- c(rescor_names_Rescor_by_levels, 
                                               tx_i)
          }
          rescor_names <- rescor_names_Rescor_by_levels
        }
      } 
      lc(out) <- rlist(grepl("^rescor\\[", pars), rescor_names)
    }
    out
  }
  bframe <- brmsframe(x$formula, x$data)
  pars <- variables(x)
  to_rename <- c(rename_predictor(bframe, pars = pars, prior = x$prior, 
                                  Rescor_by_levels = Rescor_by_levels,
                                  rescor_names = NULL),  
                 rename_re(bframe, pars = pars), rename_Xme(bframe, 
                                                            pars = pars))
  x <- save_old_par_order(x)
  x <- do_renaming(x, to_rename)
  x$fit <- repair_stanfit(x$fit)
  x <- compute_quantities(x)
  x <- reorder_pars(x)
  x
}




#' Sanitize pathfinder arguments for fit via cmdstanr
#'
#' @param sdata A list of data objects
#' @param pathfinder_args A list of argument allowed for pathfinder
#' @param brm_args A list of argument passes to the [[brms::brm()]]
#' @param ... Other arguments
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
sanitize_pathfinder_args <- function(sdata, pathfinder_args, brm_args, ...) { 
  pathfinder_args_all <- list(
    data = NULL,
    seed = NULL,
    refresh = NULL,
    init = NULL,
    save_latent_dynamics = FALSE,
    output_dir = getOption("cmdstanr_output_dir"),
    output_basename = NULL,
    sig_figs = NULL,
    opencl_ids = NULL,
    threads = NULL,
    init_alpha = NULL,
    tol_obj = NULL,
    tol_rel_obj = NULL,
    tol_grad = NULL,
    tol_rel_grad = NULL,
    tol_param = NULL,
    history_size = NULL,
    single_path_draws = NULL,
    draws = NULL,
    num_paths = 4,
    max_lbfgs_iters = NULL,
    num_elbo_draws = NULL,
    save_single_paths = NULL,
    psis_resample = NULL,
    calculate_lp = NULL,
    show_messages = TRUE,
    show_exceptions = TRUE,
    save_cmdstan_config = NULL
  )
  pathfinder_args_all_names <- names(pathfinder_args_all)
  pathfinder_args_final <- c(pathfinder_args, brm_args)
  pathfinder_args_final[['data']] <- NULL
  pathfinder_args_final_names <- names(pathfinder_args_final)
  pathfinder_args_final_valid_names <-
    setdiff(pathfinder_args_final_names, pathfinder_args_all_names) 
  for (i in pathfinder_args_final_valid_names) {
    pathfinder_args_final[[i]] <- NULL
  }
  if(!is.null(brm_args$threads$threads)) 
    pathfinder_args_final[['threads']] <- brm_args$threads$threads
  pathfinder_args_final[['data']] <- sdata
  pathfinder_args_final
}



#' check system info
#'
#' 
#' @return A logical
#' @keywords internal
#' @noRd
#'

check_system_info <- function(..., verbose = FALSE) {
  system_info <- Sys.info()
  os_name <- system_info['sysname']
  if (os_name == "Windows") {
    if(verbose) print("The system is running on Windows.")
  } else if (os_name == "Linux") {
    if(verbose) print("The system is running on Linux.")
  } else {
    if(verbose) print(paste("The system is running on:", os_name))
  }
  return(os_name)
}



#' check if cmdstanr available is available
#'
#' 
#' @return A logical
#' @keywords internal
#' @noRd
#'
check_if_cmdstanr_available <- function() {
  minimum_version <- get_package_minversion('cmdstanr')
  try(zz <- insight::check_if_installed(c("cmdstanr"), 
                                        minimum_version = minimum_version, 
                                        prompt = FALSE,
                                        stop = FALSE))
  if(!isTRUE(zz)) {
    message2c("Please install the latest version of the 'cmdstanr' 
              package",
            "\n ",
            paste0("install.packages('cmdstanr', "   ,
                   "repos = c('https://mc-stan.org/r-packages/', "   ,
                   "getOption('repos')))")
    )
    return(invisible(NULL))
  } 
  if(isTRUE(zz)) {
    write_stan_file <- utils::getFromNamespace("write_stan_file", "cmdstanr")
    cmdstan_model   <- utils::getFromNamespace("cmdstan_model", "cmdstanr")
    expose_model   <- utils::getFromNamespace("expose_functions", "cmdstanr")
    get_cmdstan_path   <- utils::getFromNamespace("cmdstan_path", "cmdstanr")
    get_cmdstan_default_path  <- 
      utils::getFromNamespace("cmdstan_default_path", "cmdstanr")
    get_cmdstan_default_install_path  <- 
      utils::getFromNamespace("cmdstan_default_install_path", "cmdstanr")
    set_cmdstan_path  <- utils::getFromNamespace("set_cmdstan_path", "cmdstanr")

  }
  return(zz)
}



# Commenting out for CRAN initial release

#' Fit model via cmdstanr
#'
#' @param scode A character string of model code
#' @param sdata A list of data objects
#' @param brm_args A list of argument passes to the [[brms::brm()]]
#' @param brms_arguments A list of argument passes to the [[brms::brm()]]
#'   especially when passing include_paths
#' @param pathfinder_args A list of argument allowed for pathfinder
#' @param pathfinder_init A logical to indicate whether to use pathfinder
#'   initials.
#'
#' @return An object of class \code{bgmfit}
#' @keywords internal
#' @noRd
#'
brms_via_cmdstanr <- function(scode, 
                              sdata, 
                              brm_args, 
                              brms_arguments = list(),
                              pathfinder_args = NULL,
                              pathfinder_init = FALSE,
                              Rescor_by_levels = NULL, 
                              verbose = FALSE) {
  
  if(isTRUE(check_if_cmdstanr_available())) {
    write_stan_file <- utils::getFromNamespace("write_stan_file", "cmdstanr")
    cmdstan_model   <- utils::getFromNamespace("cmdstan_model", "cmdstanr")
  }
  if(!is.null(brm_args$threads$threads)) {
    stan_threads <- TRUE
  } else {
    stan_threads <- FALSE
  }
  if(!is.null(brm_args$opencl)) {
    stan_opencl <- TRUE
  } else if(is.null(brm_args$opencl)) {
    stan_opencl <- FALSE
  }
  if(!is.null(brms_arguments$stan_model_args$include_paths)) {
    set_includes <- brms_arguments$stan_model_args$include_paths
  } else {
    set_includes <- NULL
  }
  cpp_options <- list(stan_threads = stan_threads)
  stanc_options <- brm_args$stan_model_args$stanc_options
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
  c_scode <- cmdstan_model(write_stan_file(scode),
                            quiet = TRUE,
                            cpp_options = cpp_options,
                            stanc_options = stanc_options,
                            dir = NULL,
                            pedantic = FALSE,
                            include_paths = set_includes,
                            user_header = NULL,
                            compile_model_methods = FALSE,
                            # compile_hessian_method = FALSE,
                            compile_standalone = FALSE)
  # iter_sampling <- brm_args$iter - brm_args$warmup
  # iter_warmup   <- brm_args$warmup
  if(is.null(brm_args$warmup)) {
    if(!is.null(brm_args$iter)) {
      brm_args$warmup <- floor(brm_args$iter / 2)
    }
  }
  if(is.null(brm_args$iter_sampling)) {
    iter_sampling <- brm_args$iter - brm_args$warmup
  } else {
    iter_sampling <- brm_args$iter_sampling
  }
  if(is.null(brm_args$iter_warmup)) {
    iter_warmup <- brm_args$warmup
  } else {
    iter_warmup <- brm_args$iter_warmup
  }
  call_pathfinder_ <- FALSE
  if(pathfinder_init | !is.null(pathfinder_args)) {
    call_pathfinder_ <- TRUE
  }
  if(call_pathfinder_) {
    if(is.null(pathfinder_args)) {
      pathfinder_args_final <- list()
      pathfinder_args_final[['refresh']] <- 0
      pathfinder_args_final[['save_cmdstan_config']] <- TRUE # NULL
      pathfinder_args_final[['show_messages']]       <- FALSE
      pathfinder_args_final[['show_exceptions']]     <- FALSE
      if(!is.null(brm_args$threads$threads)) {
        pathfinder_args_final[['threads']] <- brm_args$threads$threads
      }
      pathfinder_args_final[['data']] <- sdata
      pathfinder_args_final[['init']] <- brm_args$init
      pathfinder_args_final[['history_size']] <- 100
      pathfinder_args_final[['num_paths']] <- brm_args$chains
      
    } else if(!is.null(pathfinder_args)) {
      pathfinder_args_final <- sanitize_pathfinder_args(sdata, 
                                                        pathfinder_args, 
                                                        brm_args)
      pathfinder_args_final[['refresh']] <- 0
      pathfinder_args_final[['save_cmdstan_config']] <- TRUE 
      pathfinder_args_final[['show_messages']]       <- FALSE
      pathfinder_args_final[['show_exceptions']]     <- FALSE
    }
    if(verbose) {
      message2c("Running '$pathfinder()' for initial values")
    }
    enverr. <- environment()
    assign('err.', FALSE, envir = enverr.)
    
    tryCatch(
      expr = {
        cb_pathfinder <- CustomDoCall(c_scode$pathfinder, pathfinder_args_final)
      },
      error = function(e) {
        assign('err.', TRUE, envir = enverr.)
      }
    )
    err. <- get('err.', envir = enverr.)
    if(!err.) {
      if(!exists('cb_pathfinder')) err. <- TRUE
    }
    if (err.) {
      stop2c("Current setting of 'init' argument fails 'pathfinder()'",
           "\n  ",
           "Please try some other configuration of initial values",
           "\n  ",
           "First try 'vcov_init_0 = TRUE', if still fails, then 'init = 0'") 
    } else {
      cb_pathfinder <- cb_pathfinder
    }
    cb_pathfinder_init <- NULL
    if(pathfinder_init) {
      brm_args$init      <-  cb_pathfinder
    } else if(!pathfinder_init) {
      cb_pathfinder <- brms::read_csv_as_stanfit(cb_pathfinder$output_files(), 
                                                 model = c_scode)
      attributes(cb_pathfinder)$CmdStanModel <- c_scode
    }
  } 
  cb_fit <- c_scode$sample(
    data = sdata,
    seed = brm_args$seed,
    init = brm_args$init,
    chains = brm_args$chains,
    parallel_chains = brm_args$cores,
    threads_per_chain = brm_args$threads$threads,
    opencl_ids = brm_args$opencl,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    thin = brm_args$thin,
    max_treedepth = brm_args$control$max_treedepth,
    adapt_delta = brm_args$control$adapt_delta,
    adapt_engaged = TRUE,
    fixed_param = FALSE,
    show_messages = show_messages,
    show_exceptions = show_exceptions)
  
  cb_fit <- brms::read_csv_as_stanfit(cb_fit$output_files(), model = c_scode)
  attributes(cb_fit)$CmdStanModel <- c_scode
  brm_args_empty <- brm_args
  brm_args_empty$empty <- TRUE
  bfit     <- CustomDoCall(brms::brm, brm_args_empty)
  bfit$fit <- cb_fit
  bfit     <- custom_rename_pars(x = bfit, 
                                 Rescor_by_levels = Rescor_by_levels)
  if(!inherits(bfit, 'bgmfit')) {
    attr(bfit, 'class') <- c(attr(bfit, 'class'), 'bgmfit')
  }
  return(bfit)
}



#' Fit model via rstan
#'
#' @param scode A character string of model code
#' @param sdata A list of data objects
#' @param brm_args A list of argument passes to the brm
#' @param brms_arguments A list of argument passes to the [[brms::brm()]]
#'   especially when passing include_paths
#' @return An object of class \code{bgmfit}
#' @keywords internal
#' @noRd
#'
brms_via_rstan <- function(scode, 
                           sdata, 
                           brm_args, 
                           brms_arguments = list(),
                           Rescor_by_levels = NULL,
                           verbose = FALSE) {
  if(is.null(brm_args$warmup)) {
    if(!is.null(brm_args$iter)) {
      brm_args$warmup <- floor(brm_args$iter / 2)
    }
  }
  if(!is.null(brm_args$threads$threads)) {
    stan_threads <- TRUE
  } else {
    stan_threads <- FALSE
  }
  if(stan_threads) {
    rstan::rstan_options(threads_per_chain = brm_args$threads$threads)
  }
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
  if(!is.null(brms_arguments$stan_model_args$include_paths)) {
    set_allow_undefined <- TRUE
    set_includes <- brms_arguments$stan_model_args$include_paths
  } else {
    set_allow_undefined <- isTRUE(getOption("stanc.allow_undefined", FALSE))
    set_includes <- NULL
  }

  message2c("Compiling Stan program...")
 if(is.null(set_includes)) {
   c_scode <- rstan::stan_model(
     # file, 
     model_name = "anon_model",
     model_code = scode, 
     stanc_ret = NULL,
     boost_lib = NULL,
     eigen_lib = NULL,
     # save_dso = TRUE,
     verbose = FALSE,
     auto_write = rstan::rstan_options("auto_write"),
     obfuscate_model_name = TRUE,
     allow_undefined = set_allow_undefined,
     allow_optimizations = isTRUE(getOption("stanc.allow_optimizations", FALSE)),
     standalone_functions=isTRUE(getOption("stanc.standalone_functions", FALSE)),
     use_opencl = isTRUE(getOption("stanc.use_opencl", FALSE)),
     warn_pedantic = isTRUE(getOption("stanc.warn_pedantic", FALSE)),
     warn_uninitialized = isTRUE(getOption("stanc.warn_uninitialized", FALSE)),
     includes = set_includes,
     isystem = c(if (!missing(file)) dirname(file), getwd()))
 } else {
   c_scode <- rstan::stan_model(
     model_code = scode)
 }
  
  message2c("Start sampling")
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

  brm_args$empty <- TRUE
  bfit      <- CustomDoCall(brms::brm, brm_args)
  bfit$fit  <- cb_fit
  bfit      <- custom_rename_pars(x = bfit, 
                                  Rescor_by_levels = Rescor_by_levels) 
  if(!inherits(bfit, 'bgmfit')) {
    attr(bfit, 'class') <- c(attr(bfit, 'class'), 'bgmfit')
  }
  return(bfit)
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
      message2c('Checking required packages for ', getfun, " ",
              "\n ",
              collapse_comma(pkgs))
    } 
  } 
  if(installpkg) {
    # message2c('Installing required packages',
    #         collapse_comma(required_pkgs))
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
  pob    <- CustomDoCall(plot_curves, args)
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
  getfb <- transform_sec_axis(pob$data$Estimate.x, pob$data$Estimate.y)
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
#' @return A list comprised of exposed functions.
#' 
#' 
#' @keywords internal
#' @noRd
#'
check_pkg_version_exists <- function(pkg, 
                                     minimum_version = NULL, 
                                     verbose = FALSE,
                                     ...) {
  try(zz <- insight::check_if_installed(pkg, 
                                        minimum_version = minimum_version,
                                        ...))
  if(!isTRUE(zz)) {
    if(verbose) {
      message2c("Please install the latest version of the 'brms' package",
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
#' @param ... other arguments. Currently ignored.
#' @inherit growthparameters.bgmfit params
#' 
#' @return A list comprised of exposed functions.
#' 
#' @keywords internal
#' @noRd
#'
check_if_functions_exists <- function(model, 
                                      o = NULL, 
                                      xcall = NULL, 
                                      verbose = TRUE, 
                                      usesavedfuns = FALSE, 
                                      checks = FALSE,
                                      envir = NULL, 
                                      ...) {
  if(is.null(envir)) envir <- globalenv()
  if(!checks) {
    if(is.null(o)) stop2c("object 'o' must be specified")
  }
  check_brms_v <- 
    check_pkg_version_exists('brms', 
                             minimum_version = get_package_minversion('brms'), 
                             prompt = FALSE,
                             stop = FALSE,
                             verbose = FALSE)
  latest_brms_v <- TRUE
  if(!isTRUE(check_brms_v)) {
    latest_brms_v <- FALSE
  }
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
        if(!usesavedfuns & latest_brms_v) message2c(msg1)
        if(usesavedfuns & latest_brms_v) message2c(msg2)
        
        if(!usesavedfuns & !latest_brms_v) message2c(msg3)
        if(usesavedfuns & !latest_brms_v) message2c(msg3)
        
      }
    }
    return(invisible(NULL))
  }
  
  if(exists(o[[1]], mode = "function", envir = envir)) {
    envgtf <- TRUE
  } else {
    envgtf <- FALSE
  }
  
  if(verbose) {
    if(!envgtf) {
      if(verbose) message2c(msg3)
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
#' @param ... other arguments. Currently ignored.
#' 
#' @inherit growthparameters.bgmfit params
#' 
#' @keywords internal
#' 
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
    } else if(calname.fun == "get_predictions") {
      package <- c('tidyr', 'collapse')
    } else  {
      return(invisible(NULL))
    }
  } 
  
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
                                prompt = prompt)
    return(invisible(NULL))
  } 
}



#' An internal function to get the environment of an object
#'
#' @param x A symbol or a character string.
#' @param geteval A logical (default \code{TRUE}) to indicate whether to return
#' the object as a character string or as an environment.
#' 
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


#' An internal function to get the 'model' name from the arguments
#'
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'
get_lhs_pipe <- function(){
  calls <- sys.calls()
  call_firsts <- lapply(calls,`[[`,1) 
  pipe_calls <- vapply(call_firsts,identical,logical(1),quote(`%>%`))
  if(all(!pipe_calls)){
    out <- NULL
  } else {
    pipe_calls <- which(pipe_calls)
    pipe_calls <- pipe_calls[length(pipe_calls)]
    #Get the second element of the pipe call
    this_call <- calls[[c(pipe_calls,2)]]
    while(is.call(this_call) && identical(this_call[[1]],quote(`%>%`))){
      this_call <- this_call[[2]]
    }
    out <- this_call
  }
  return(out)
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

  if(!is.data.frame(df)) stop2c("df should be a data frame")
  all.dfnames  <- colnames(df)
  if(is.null(factor.dummy)) {
    dfout <- df
  } else if(!is.null(factor.dummy)) {
    if(!is.null(factor.name)) {
      if(!is.character(factor.name)) {
        stop2c("factor.name should be a character string")
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
        stop2c("Lengths of factor.dummy and factor.level must be same")
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
  . <- NULL;
  if(is.null(Parametername)) {
    Parametername <- "Parameter"
  }
  if(is.null(parmcols)) {
    parmcols <- c('Estimate', "Est.Error", "Q2.5", "Q97.5")
  }
  if(is.null(nonparmcols))  {
    stop2c("Please specify the 'nonparmcols'")
  }
  if(is.null(byjoincols)) {
    stop2c("Please specify the 'byjoincols'")
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
  tojoinit2all <- whati_list %>% CustomDoCall(dplyr::bind_cols, .) %>% data.frame()
  tojoinit12 <- cbind(tojoinit1, tojoinit2all)
  mergebycols <- intersect(nonparmcols, byjoincols)
  setdiffcols <- setdiff(byjoincols, nonparmcols)
  if(length(setdiffcols) != 0) {
    stop2c("Variable(s) ", collapse_comma(setdiffcols),
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
      r <- CustomDoCall(rbind, c(argl[-1L], list(deparse.level = deparse.level)))
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
    r <- CustomDoCall(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
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
      argin <- gsub_space(argin)
      argin <- paste(argin, collapse=",")
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
#' @param idvar A character string specifying the group identifier
#' @param resp A character string specifying the response variable (default
#'   \code{NULL})
#' @param verbose A logical to indicate whether to print relevant information.
#' @keywords internal
#' @return A data frame.
#' @keywords internal
#' @noRd
#'
check_newdata_args <- function(model, 
                               newdata, 
                               idvar, 
                               resp = NULL,
                               verbose = FALSE) {
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  groupvar_ <- paste0('groupvar', resp_rev_)
  if (is.null(model$model_info[[groupvar_]] )) {
    name_hypothetical_id <- paste0("id", resp_rev_)
    model$model_info[[groupvar_]]  <- name_hypothetical_id
    newdata[[name_hypothetical_id]] <- as.factor("tempid")
  } else if (!is.null(model$model_info[[groupvar_]] )) {
    m_m_groupvar_  <- model$model_info[[groupvar_]] 
    m_m_length     <- sapply(m_m_groupvar_, function(x) length(newdata[[x]]))
    max_m_m_length <- max(m_m_length)
    if(max_m_m_length == 0) {
      if(length(idvar) > 1) {
        name_hypothetical_id <- idvar[1] 
      } else {
        name_hypothetical_id <- idvar
      }
      model$model_info[[groupvar_]]   <- name_hypothetical_id
      newdata[[name_hypothetical_id]] <- as.factor("tempid")
    }
  }
  return(newdata)
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
  if(is.null(fullgrid)) stop2c("Please specify fullgrid")
  if(is.null(all_vars)) stop2c("Please specify all_vars")
  if(is.null(nested_vars)) stop2c("Please specify nested_vars")
  if(is.null(fulldata) & is.null(varsvector)) 
    stop2c("Please specify at least fulldata or varsvector")
  if(!is.null(fulldata) & !is.null(varsvector)) 
    stop2c("Please specify either fulldata or varsvector, not both")
  
  if(!is.null(varsvector)) {
    # if(!is.vector(varsvector)) stop2c("varsvector must be a vector")
    if(!is.factor(varsvector)) stop2c("varsvector must be a factor vector")
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



### avoid ggtext - it depends on jpeg whihc fails rmdcheck on ubuntu

#' Title ggtextelementmarkdown
#' 
#' @noRd
#'
ggtextelementmarkdown <- function (family = NULL, 
                                   face = NULL, 
                                   size = NULL,
                                   colour = NULL, 
                                   fill = NULL, 
                                   box.colour = NULL, 
                                   linetype = NULL, 
                                   linewidth = NULL, 
                                   hjust = NULL, 
                                   vjust = NULL, 
                                   halign = NULL, 
                                   valign = NULL, 
                                   angle = NULL, 
                                   lineheight = NULL, 
                                   margin = NULL, 
                                   padding = NULL, 
                                   r = NULL, 
                                   color = NULL, 
                                   box.color = NULL, 
                                   align_widths = NULL, 
                                   align_heights = NULL, 
                                   rotate_margins = NULL, 
                                   debug = FALSE, 
                                   inherit.blank = FALSE) {
    if (!is.null(color)) 
      colour <- color
    if (!is.null(box.color)) 
      box.colour <- box.color
    structure(list(family = family, face = face, size = size, 
                   colour = colour, fill = fill, box.colour = box.colour, 
                   linetype = linetype, linewidth = linewidth, hjust = hjust, 
                   vjust = vjust, halign = halign, valign = valign, 
                   angle = angle, 
                   lineheight = lineheight, margin = margin, padding = padding, 
                   r = r, align_widths = align_widths, align_heights = align_heights, 
                   rotate_margins = rotate_margins, debug = debug, 
                   inherit.blank = inherit.blank), 
              class = c("element_markdown", "element_text", "element"))
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
  try(insight::check_if_installed(c("ggplot2"), stop = FALSE, prompt = FALSE))
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
    ggplot2::theme(axis.text.x = ggtextelementmarkdown())
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
  try(insight::check_if_installed(c("ggplot2"), stop = FALSE, prompt = FALSE))
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
    ggplot2::theme(axis.text.y = ggtextelementmarkdown())
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
  try(insight::check_if_installed(c("ggplot2"), stop = FALSE, prompt = FALSE))
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
    ggplot2::theme(axis.text.x = ggtextelementmarkdown())
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
  try(insight::check_if_installed(c("ggplot2"), stop = FALSE, prompt = FALSE))
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
    ggplot2::theme(axis.text.y = ggtextelementmarkdown())
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
  if(!is.character(pkg)) stop2c('pkg must be a character')
  if(pkg == 'brms') {
    if(is.null(version)) {
      out <- '2.21.0' 
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'marginaleffects') {
    if(is.null(version)) {
      out <- '0.19.0'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'data.table') {
    if(is.null(version)) {
      out <- '1.15.4'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'dtplyr') {
    if(is.null(version)) {
      out <- '1.3.1'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'checkmate') {
    if(is.null(version)) {
      out <- '2.3.1'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'collapse') {
    if(is.null(version)) {
      out <- '2.0.13'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'doParallel') {
    if(is.null(version)) {
      out <- '1.0.17'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'foreach') {
    if(is.null(version)) {
      out <- '1.5.2'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'parallel') {
    if(is.null(version)) {
      out <- '0.0.1' # '4.3.1'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
      out <- version
    }
  }
  if(pkg == 'cmdstanr') {
    if(is.null(version)) {
      out <- '0.7.1'
    } else {
      if(!is.character(version)) stop2c('version must be a character')
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
  if(!is.character(algorithm)) stop2c('algorithm must be a character')
  if(algorithm == "sampling" | 
     algorithm == "meanfield" |
     algorithm == "fullrank" |
     algorithm == "fixed_param") {
    return(args)
  }
  pathfinderargs <- c('save_latent_dynamics', 'output_dir',
                      'output_basename', 'sig_figs', 
                      'threads', 'init_alpha', 'tol_obj',
                      'tol_rel_obj', 'tol_grad', 'tol_rel_grad',
                      'tol_param', 'history_size', 'single_path_draws',
                      'draws', 'num_paths', 'max_lbfgs_iters', 
                      'num_elbo_draws', 'save_single_paths')
  laplacerargs <- c('save_latent_dynamics', 'output_dir',
                    'output_basename', 'sig_figs', 
                    'mode', 'opt_args', 'jacobian',
                    'draws')
  remove_for_brms_pathfinder <- c('adapt_delta', 'max_treedepth', 'control')
  if(!'pathfinder' %in% algorithm) {
    for (i in pathfinderargs) {
      if(!is.null(args[[i]])) args[[i]] <- NULL
    }
  } 
  if(!'laplace' %in% algorithm) {
    for (i in laplacerargs) {
      if(!is.null(args[[i]])) args[[i]] <- NULL
    }
  } 
  if('pathfinder' %in% algorithm) {
    for (i in remove_for_brms_pathfinder) {
      args[[i]] <- NULL
    }
  } 
  return(args)
}


#' An internal function to convert velocity parameter from exp to unit/time and vice versa
#'
#' @param x A numeric value or a vector.
#' @param to A character string to indicate direction of conversion.
#' @keywords internal
#' @return A list comprised of character strings.
#' @noRd
#'
vel_exp_unit_convert <- function(x, to = 'unit') {
  if(to == 'unit') {
    message2c("converted from exp(x) to unit/time")
    out <-  1 - exp(x)^2
  }
  if(to == 'exp') {
    message2c("converted from unit/time to exp(x)")
    out <-  log(sqrt(1-x))
  }
  out
}





#' Title
#' 
#' @details A customized version of 'plot.see_equivalence_test' that allows
#'   plotting equivalence_test without model. In other words, it works for
#'   numeric test also.
#' 
#' @param x A 'equivalence_test' object
#' 
#' @param parms_data A data frame
#' 
#' @param rope_color see bayestestR
#' 
#' @param rope_alpha see bayestestR
#' 
#' @param rope.line.alpha see bayestestR
#' 
#' @param n_columns see bayestestR
#' 
#' @param fill.color see bayestestR
#' 
#' @param legend.title see bayestestR
#'
#' @return A plot object if (\code{return_plot = TRUE}).
#' 
#' @keywords internal
#' @noRd
#'
plot_equivalence_test <-  function(x,
                                   parms_data,
                                   rope_color = "#0171D3",
                                   rope_alpha = 0.5,
                                   rope.line.alpha = NULL,
                                   n_columns = 1,
                                   fill.color = c("#CD423F", 
                                                  "#018F77", 
                                                  "#FCDA3B"),
                                   legend.title = "Decision on H0") {
  predictor <- NULL;
  estimate <- NULL;
  grp <- NULL;
  insight::check_if_installed("ggridges", prompt = FALSE)
  .has_multiple_panels <- function (x) {
    (!"Effects" %in% names(x) || insight::n_unique(x$Effects) <=
       1L) &&
      (!"Component" %in% names(x) || insight::n_unique(x$Component) <=
         1L)
  }
  .clean_parameter_names <- function (params, grid = FALSE)
  {
    params <- unique(params)
    parameter_labels <- params
    params <- gsub("(b_|bs_|bsp_|bcs_)(.*)", "\\2", params, perl = TRUE)
    params <- gsub("^zi_(.*)", "\\1 (Zero-Inflated)", params, perl = TRUE)
    params <- gsub("(.*)_zi$", "\\1 (Zero-Inflated)", params, perl = TRUE)
    params <- gsub("(.*)_disp$", "\\1 (Dispersion)", params, perl = TRUE)
    params <- gsub("r_(.*)\\.(.*)\\.", "(re) \\1", params)
    params <- gsub("b\\[\\(Intercept\\) (.*)\\]", "(re) \\1", params)
    params <- gsub("b\\[(.*) (.*)\\]", "(re) \\2", params)
    params <- gsub("^smooth_sd\\[(.*)\\]", "\\1 (smooth)", params)
    params <- gsub("^sds_", "\\1 (Smooth)", params)
    params <- gsub("(.*)(\\.)(\\d)$", "\\1 \\3", params)
    params <- gsub("(.*)__zi\\s(.*)", "\\1 \\2 (Zero-Inflated)", params, perl = TRUE)
    params <- gsub("\\(re\\)\\s(.*)", "\\1 (Random)", params, perl = TRUE)
    cor_sd <- grepl("(sd_|cor_)(.*)", params)
    if (any(cor_sd)) {
      params[cor_sd] <- paste("SD/Cor: ",
                              gsub("^(sd_|cor_)(.*?)__(.*)", "\\3", params[cor_sd], perl = TRUE))
      cor_only <- !is.na(params[cor_sd]) &
        startsWith(params[cor_sd], "cor_")
      if (any(cor_only)) {
        params[cor_sd][which(cor_sd)[cor_only]] <- sub("__", " ~ ", params[cor_sd][which(cor_sd)[cor_only]], fixed = TRUE)
      }
    }
    cor_sd <- grepl("^Sigma\\[(.*)", params)
    if (any(cor_sd)) {
      parm1 <- gsub("^Sigma\\[(.*):(.*),(.*)\\]", "\\2", params[cor_sd], perl = TRUE)
      parm2 <- gsub("^Sigma\\[(.*):(.*),(.*)\\]", "\\3", params[cor_sd], perl = TRUE)
      params[which(cor_sd)] <- parm1
      rand_cor <- parm1 != parm2
      if (any(rand_cor)) {
        params[which(cor_sd)[rand_cor]] <- paste0(parm1[rand_cor], " ~ ", parm2[rand_cor])
      }
      params[cor_sd] <- paste("SD: ", params[cor_sd])
    }
    if (grid) {
      params <- trimws(gsub("(Zero-Inflated)", "", params, fixed = TRUE))
      params <- trimws(gsub("(Random)", "", params, fixed = TRUE))
      params <- trimws(gsub("(Dispersion)", "", params, fixed = TRUE))
    }
    else {
      params <- gsub("(Zero-Inflated) (Random)",
                     "(Random, Zero-Inflated)",
                     params,
                     fixed = TRUE)
    }
    stats::setNames(params, parameter_labels)
  }
  .fix_facet_names <- function (x)
  {
    if ("Component" %in% names(x)) {
      x$Component <- as.character(x$Component)
      if ("Effects" %in% names(x)) {
        x$Component[x$Component == "conditional"] <- "(Conditional)"
        x$Component[x$Component == "zero_inflated"] <- "(Zero-Inflated)"
        x$Component[x$Component == "dispersion"] <- "(Dispersion)"
        x$Component[x$Component == "simplex"] <- "(Monotonic Effects)"
      }
      else {
        x$Component[x$Component == "conditional"] <- "Conditional"
        x$Component[x$Component == "zero_inflated"] <- "Zero-Inflated"
        x$Component[x$Component == "dispersion"] <- "Dispersion"
        x$Component[x$Component == "simplex"] <- "Monotonic Effects"
      }
    }
    if ("Effects" %in% names(x)) {
      x$Effects <- as.character(x$Effects)
      x$Effects[x$Effects == "fixed"] <- "Fixed Effects"
      x$Effects[x$Effects == "random"] <- "Random Effects"
    }
    x
  }
  .reshape_to_long <- function(x,
                               names_to = "group",
                               values_to = "values",
                               columns = colnames(x),
                               id = "id") {
    if (is.numeric(columns))
      columns <- colnames(x)[columns]
    dat <- stats::reshape(
      as.data.frame(x),
      idvar = id,
      ids = row.names(x),
      times = columns,
      timevar = names_to,
      v.names = values_to,
      varying = list(columns),
      direction = "long"
    )
    if (is.factor(dat[[values_to]])) {
      dat[[values_to]] <- as.character(dat[[values_to]])
    }
    dat[, 1:(ncol(dat) - 1), drop = FALSE]
  }

  x$Effects <- "fixed"
  x$Component <- "conditional"
  attr(x, "Cleaned_Parameter") <- x$Parameter
  attr(x, "object_name") <- "model"
  .rope <- c(x$ROPE_low[1], x$ROPE_high[1])
  tests <- split(x, x$CI)

  result <- lapply(tests, function(i) {
    tmp <- parms_data[, i$Parameter, drop = FALSE]
    tmp2 <- lapply(seq_len(nrow(i)), function(j) {
      p <- i$Parameter[j]
      tmp[[p]][tmp[[p]] < i$HDI_low[j]] <- NA
      tmp[[p]][tmp[[p]] > i$HDI_high[j]] <- NA
      tmp[[p]]
    })
    cnames <- colnames(tmp)
    tmp <- as.data.frame(tmp2)
    colnames(tmp) <- cnames
    tmp <- .reshape_to_long(tmp, names_to = "predictor", values_to = "estimate")
    tmp$grp <- NA
    for (j in seq_len(nrow(i))) {
      tmp$grp[tmp$predictor == i$Parameter[j]] <- i$ROPE_Equivalence[j]
    }
    tmp$predictor <- factor(tmp$predictor)
    tmp$predictor <- factor(tmp$predictor, levels = rev(levels(tmp$predictor)))
    tmp$HDI <- sprintf("%g%% HDI", 100 * i$CI[1])
    tmp
  })
  tmp <- CustomDoCall(rbind, result)
  if (.has_multiple_panels(tmp)) {
    n_columns <- NULL
  }
  labels <- .clean_parameter_names(tmp$predictor, grid = !is.null(n_columns))
  tmp <- .fix_facet_names(tmp)
  if (length(unique(tmp$HDI)) > 1L) {
    x.title <- "Highest Density Region of Posterior Samples"
  } else {
    x.title <- sprintf("%g%% Highest Density Region of Posterior Samples", 100 * x$CI[1])
  }
  
  fill.color <- fill.color[sort(unique(match(
    x$ROPE_Equivalence, c("Accepted", "Rejected", "Undecided")
  )))]
  add.args <- lapply(match.call(expand.dots = FALSE)$`...`, function(x)
    x)
  if ("colors" %in% names(add.args))
    fill.color <- eval(add.args[["colors"]])
  if ("x.title" %in% names(add.args))
    x.title <- eval(add.args[["x.title"]])
  if ("legend.title" %in% names(add.args))
    legend.title <- eval(add.args[["legend.title"]])
  if ("labels" %in% names(add.args))
    labels <- eval(add.args[["labels"]])
  if (is.null(rope.line.alpha)) {
    rope.line.alpha <- 1.25 * rope_alpha
  }
  if (rope.line.alpha > 1)
    rope.line.alpha <- 1
  p <- ggplot2::ggplot(tmp, ggplot2::aes(x = estimate, 
                                         y = predictor, fill = grp)) +
    ggplot2::annotate(
      "rect",
      xmin = .rope[1],
      xmax = .rope[2],
      ymin = 0,
      ymax = Inf,
      fill = rope_color,
      alpha = (rope_alpha / 3),
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = .rope,
      linetype = "dashed",
      colour = rope_color,
      alpha = rope.line.alpha,
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      colour = rope_color,
      linewidth = 0.8,
      alpha = rope.line.alpha,
      na.rm = TRUE
    ) +
    ggridges::geom_density_ridges2(
      rel_min_height = 0.01,
      scale = 2,
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(values = fill.color) +
    ggplot2::labs(x = x.title, y = NULL, fill = legend.title) +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::theme(legend.position = "bottom")
  if (!is.null(n_columns)) {
    if ("Component" %in% names(x) && "Effects" %in% names(x)) {
      if (length(unique(tmp$HDI)) > 1L) {
        p <- p + ggplot2::facet_wrap( ~ Effects + Component + HDI,
                                      scales = "free",
                                      ncol = n_columns)
      } else {
        p <- p + ggplot2::facet_wrap( ~ Effects + Component,
                                      scales = "free",
                                      ncol = n_columns)
      }
    } else if ("Effects" %in% names(x)) {
      if (length(unique(tmp$HDI)) > 1L) {
        p <- p + ggplot2::facet_wrap( ~ Effects + HDI, scales = "free", ncol = n_columns)
      } else {
        p <- p + ggplot2::facet_wrap( ~ Effects, scales = "free", ncol = n_columns)
      }
    } else if ("Component" %in% names(x)) {
      if (length(unique(tmp$HDI)) > 1L) {
        p <- p + ggplot2::facet_wrap( ~ Component + HDI, scales = "free", ncol = n_columns)
      } else {
        p <- p + ggplot2::facet_wrap( ~ Component, scales = "free", ncol = n_columns)
      }
    }
  } else {
    if (length(unique(tmp$HDI)) > 1L) {
      p <- p + ggplot2::facet_wrap( ~ HDI, scales = "free", ncol = n_columns)
    }
  }
  p
}



#' Get information on dpar
#' 
#' @details
#' This checkresp_info is used in getmodel_info and post_processing_checks 
#' functions
#' 
#' @param model An object of class \code{bgmfit} 
#' @param resp A character string
#'
#' @return An object of class \code{bgmfit} 
#' @keywords internal
#' @noRd
#'
checkresp_info <- function(model, resp) {
  uvarby <- model$model_info$univariate_by$by
  if(is.null(uvarby)) uvarby <- NA 
  if(is.null(resp)) check_resp_str <- FALSE else check_resp_str <- TRUE
  error_message <- NULL
  if (model$model_info$nys == 1 & !is.null(resp)) {
    error_message <-  
    paste0("You have fit a univariate model",
      " but set resp argument as: ",
      resp,
      ".",
      "\n ",
      " For univariate model, the resp option should be NULL",
      "\n ",
      " (i.e., resp = NULL)"
    )
  }
  if (model$model_info$nys > 1 & is.null(resp)) {
    if (!is.na(uvarby)) {
      error_message <-  
        paste0("You have fit a univariate_by model for ",
               uvarby,
        "\n ",
        " but did not correctly specified the 'resp' argument",
        " (which is NULL at present).",
        "\n ",
        " The response options are: ",
        collapse_comma(model$model_info$yvars)
      )
    }
    if (model$model_info$multivariate$mvar) {
      error_message <-  
        paste0("You have fit a multivariate model but did not set the ",
        "\n ",
        "'resp' argument correctly (which is 'NULL' at present).",
        "\n  ",
        "The avaialble response options are: ",
        collapse_comma(model$model_info$yvars)
      )
    }
  }
  if(!is.null(error_message)) {
    error_message <- stop2c(clean_text_spaces(error_message))
  }
  
  if(check_resp_str) {
    if (!is.na(uvarby)) {
      model_str <- 'univariate_by'
    } else if (model$model_info$multivariate$mvar) {
      model_str <- 'multivariate'
    } else {
      model_str <- 'univariate'
    }
    if(!resp %in% model$model_info$yvars) {
      error_message <- sprintf(
        "You have fit a '%s' model but did not set the \n '%s' argument 
        correctly (which is '%s' at present).\n  The avaialble 
        response options are: %s",
        model_str,
        "resp",
        resp,
        collapse_comma(model$model_info$yvars)
      )
    } else {
      error_message <- NULL
    }
    if(model_str != 'univariate') {
      if(!is.null(error_message)) {
        error_message <- stop2c(clean_text_spaces(error_message))
      }
    }
  }
  
  return(invisible(NULL))
}




#' allowed_namespace_for_sigma_d1
#'
#' @return A character string
#' @keywords internal
#' @noRd
#'
allowed_namespace_for_sigma_d1 <- function() {
  c('splines2', 'bsitar')
}




#' exclude_global_for_sigma_d1
#'
#' @return A character string
#' @keywords internal
#' @noRd
#'
exclude_global_for_sigma_d1 <- function() {
  c('brms')
}



#' Get information on dpar
#'
#' @param model An object of class \code{bgmfit} 
#' @param dpar A logical or a character string \code{'mu'} or \code{'sigma'}
#' @param resp A character string
#' @param deriv default \code{NULL} 
#' @param remove_cons_cov Logical, calls 'remove_cons_as_numeric_cov'
#' @param strict_1 Logical, passed on to 'remove_cons_as_numeric_cov'
#' @return An object of class \code{bgmfit} 
#' @keywords internal
#' @noRd
#'
getmodel_info <- function(model, 
                          dpar, 
                          resp, 
                          deriv = NULL, 
                          remove_cons_cov = TRUE,
                          strict_1 = TRUE,
                          verbose = FALSE) {
  if(is.null(model$test_mode)) {
    model[['test_mode']] <- FALSE
    if(verbose) {
      message2c("The 'model' must include a test_mode element, accessible as 
                model[['test_mode']], and it must be a logical 
                value (TRUE or FALSE). In the CRAN version of berkeley_exfit,
                model[['test_mode']] is set to TRUE. Setting 
                model[['test_mode']] to FALSE allows access to the full 
                data via insight::get_data(), which is required by 
                the marginaleffects functions. If model[['test_mode']] 
                is NULL, it is set to FALSE.", 
                pad_before = "\n", 
                pad_after = "\n")
    }
  }
  checkresp_info(model, resp)
  if (is.null(resp)) {
    resp_    <- resp
    revresp_ <- ""
  } else if (!is.null(resp)) {
    resp_    <- paste0(resp, "_")
    revresp_ <- paste0("_", resp)
  }
  setsigmaxvars_ <- paste0('setsigmaxvar', revresp_)
  sigma_model_      <- paste0('sigmamodel', revresp_)
  sigma_model_name_ <- paste0('sigmabasicfunname', revresp_)
  sigma_model_attr_ <- paste0('sigmabasicfunattr', revresp_)
  sigma_model       <- model$model_info[[sigma_model_]]
  sigma_model_name  <- model$model_info[[sigma_model_name_]]
  sigma_model_attr  <- model$model_info[[sigma_model_attr_]]
  sigma_model_is_ls <- FALSE
  if(!is.null(sigma_model)) {
    if(sigma_model == "ls") {
      sigma_model_is_ls <- TRUE
    }
  }
  oxx <- model$model_info[['namesexefuns']]
  if(is.null(dpar)) {
    oxx <- oxx[!grepl("sigma", oxx)]
    sigma_fun_mode <- NULL
  } else if(!is.null(dpar)) {
    if(dpar == "mu") {
      oxx <- oxx[!grepl("sigma", oxx)]
      sigma_fun_mode <- NULL
    }
    if(dpar == "sigma") {
      if(model$model_info[[setsigmaxvars_]]) {
        if(is.null(model$model_info[['sigma_fun_mode']])) {
          sigma_fun_mode <- "custom"
        } else {
          sigma_fun_mode <- model$model_info[['sigma_fun_mode']]
        }
      } else if(!model$model_info[[setsigmaxvars_]]) {
        if(is.null(model$model_info[['sigma_fun_mode']])) {
          sigma_fun_mode <- "inline"
        } else {
          sigma_fun_mode <- model$model_info[['sigma_fun_mode']]
        }
      }
      if(sigma_fun_mode == "custom") {
        oxx <- oxx[grepl("sigma", oxx)]
      } else {
        oxx <- oxx[!grepl("sigma", oxx)]
      }
    } 
  } 
  if(sigma_model_is_ls) {
    model$model_info[['namesexefuns']] <- oxx
    model$model_info[['sigma_fun_mode']] <- sigma_fun_mode
  }
  if(remove_cons_cov) {
    cov_       <- paste0('cov', revresp_)
    sigmacov_  <- paste0('sigma', cov_)
    model$model_info[[cov_]] <- 
      remove_cons_as_numeric_cov(model$model_info[[cov_]],
                                 model$model_info[['bgmfit.data']],
                                 strict_1 = strict_1,
                                 verbose = verbose)
    model$model_info[[sigmacov_]] <- 
      remove_cons_as_numeric_cov(model$model_info[[sigmacov_]],
                                 model$model_info[['bgmfit.data']],
                                 strict_1 = strict_1,
                                 verbose = verbose)
  } 
  return(model)
}



#' An internal function to perform checks when calling post-processing functions
#'
#' @description The \code{post_processing_checks} perform essential checks (such
#'   as the validity of model class, response etc.) during post-processing of
#'   posterior draws.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @param xcall The \code{match.call()} from the post-processing function.
#'
#' @param resp Response variable (default \code{NULL}) specified as a character
#'   string. This is needed when processing \code{univariate_by} and
#'   \code{multivariate} models (see \code{bgmfit} function for details).
#'
#' @param deriv An integer value to specify whether to estimate distance curve
#'   (i.e., model estimated curve(s)) or the velocity curve (first derivative of
#'   the model estimated curve(s)). A value \code{0} (default) is for distance
#'   curve and  \code{1} for the velocity curve.
#'   
#' @param all A logical (default \code{NULL}) to specify whether to return all
#'   the exposed functions.
#'
#' @param envir Indicator to set the environment of function evaluation. The
#'   default is \code{parent.frame}.
#'   
#' @param verbose A logical (default \code{FALSE}) to indicate whetehr to print
#'  \code{warnings} and \code{messages} during the function evaluation.
#'
#' @return A string with the error captured or else a list with necessary
#'   information needed when executing the post-processing function
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
post_processing_checks <- function(model, 
                                   xcall, 
                                   resp = NULL, 
                                   deriv = NULL,
                                   all = FALSE,
                                   envir = NULL, 
                                   verbose = FALSE, 
                                   check_d0 = FALSE,
                                   check_d1 = FALSE,
                                   check_d2 = FALSE) {
  if(is.null(envir)) envir <- parent.frame()
  if(is.null(deriv)) deriv <- 0
  if(!'bgmfit' %in% class(model)) {
    stop2c("The class of model object should be 'bgmfit' ")
  }
  excall_ <- c("plot_ppc", "loo_validation")
  xcall_check_it <- paste(deparse(substitute(xcall)), collapse = "")
  xcall_check_it <- gsub_space(xcall_check_it)
  check_it       <- sub(" *\\(.*", "", xcall_check_it)
  check_it_sss <- strsplit(check_it, "\\.")[[1]][1]
  if(!is.null(model$xcall)) {
    if(grepl("get_predictions", model$xcall)) {
      xcall <- "get_predictions"
    } else if(grepl("modelbased_growthparameters", model$xcall)) {
      xcall <- "modelbased_growthparameters"
    } else if(grepl("get_growthparameters", model$xcall)) {
      xcall <- "get_growthparameters"
    }
  } else {
    rlang_trace_back <- rlang::trace_back()
    check_trace_back.bgmfit <- grepl(".bgmfit", rlang_trace_back[[1]])
    if(all(!check_trace_back.bgmfit)) {
      # 
    } else {
      rlang_trace_back.bgmfit_i <- min(which(check_trace_back.bgmfit == TRUE))
      rlang_trace_back.bgmfit <- rlang_trace_back[[1]][[rlang_trace_back.bgmfit_i]]
      rlang_call_name <- rlang::call_name(rlang_trace_back.bgmfit)
      xcall <- rlang_call_name
    }
  }
  call_checkresp_info <- TRUE
  if(xcall == "optimize_model.bgmfit" | xcall == "optimize_model") {
    call_checkresp_info <- FALSE
  }
  if(call_checkresp_info) {
    checkresp_info(model, resp)
  }
  if (is.null(resp)) {
    resp_    <- resp
    revresp_ <- ""
  } else if (!is.null(resp)) {
    resp_    <- paste0(resp, "_")
    revresp_ <- paste0("_", resp)
  }
  if(model$model_info[['expose_method']] == 'R') {
    assign(paste0(resp_, 
                  model$model_info[['namesexefuns']], 
                  '0'), 
           model$model_info$exefuns[[paste0(resp_, 
                                            model$model_info[['namesexefuns']], 
                                            '0')]], envir = envir)
    if(model$model_info[['select_model']] == 'sitar' |
       model$model_info[['select_model']] == 'rcs') {
      assign(paste0(resp_, 'getKnots'), 
             model$model_info$exefuns[[paste0(resp_, 'getKnots')]], 
             envir = envir)
    }
  }
  if(!all) {
    out <-
      list(
        paste0(resp_, model$model_info[['namesexefuns']], ''),
        paste0(resp_, model$model_info[['namesexefuns']], deriv)
      ) 
  }
  if(all) {
    out <- model$model_info[['exefuns']]
  } 
  if(check_d1) {
    available_d0 <- available_d1 <- available_d2 <- FALSE
    for (i in names(model$model_info$exefuns)) {
      check_funds <- ifelse(grepl("\\d$", i), sub(".*?(\\d+)$", "\\1", i), "")
      if(grepl("0", check_funds)) {
        available_d0 <- TRUE
      }
      if(grepl("1", check_funds)) {
        available_d1 <- TRUE
      }
      if(grepl("2", check_funds)) {
        available_d2 <- TRUE
      }
    }
    if(verbose) { 
      if(!available_d0) {
        stop2c("No 'd0' found")
      }
      if(!available_d1) {
        message2c("No 'd1' found, setting 'model_deriv = FALSE', 'deriv = 0'")
      }
      if(!available_d2) {
        #
      }
    }
    setsigmaxvars_ <- paste0('setsigmaxvar', revresp_)
    sigma_model_      <- paste0('sigmamodel', revresp_)
    sigma_model_name_ <- paste0('sigmabasicfunname', revresp_)
    sigma_model_attr_ <- paste0('sigmabasicfunattr', revresp_)
    sigma_model       <- model$model_info[[sigma_model_]]
    sigma_model_name  <- model$model_info[[sigma_model_name_]]
    sigma_model_attr  <- model$model_info[[sigma_model_attr_]]
    sigma_model_is_ba_set_d0_as_d1 <- FALSE
    sigma_model_is_ba_set_d0_as_d1_funs <- list()
    sigma_model_is_ls <- FALSE
    if(!is.null(sigma_model)) {
      if(sigma_model == "ls") {
        sigma_model_is_ls <- TRUE
        available_d1 <- available_d1
        available_d2 <- available_d2
      } else if(sigma_model == "basic") {
        if(all(sigma_model_attr %in% allowed_namespace_for_sigma_d1())) {
          sigma_model_is_ba_set_d0_as_d1 <- TRUE
          available_d1 <- available_d1
          available_d2 <- available_d2
          if(deriv > 0) {
            # stixt <- 0
            for (i in sigma_model_name) {
              # stixt <- stixt + 1
              tempfunx <- model$model_info$exefuns[[i]]
              formals(tempfunx)[['derivs']] <- deriv
              sigma_model_is_ba_set_d0_as_d1_funs[[i]] <- tempfunx
            }
          }
        } else {
          sigma_model_is_ba_set_d0_as_d1 <- FALSE
          available_d1 <- FALSE
          available_d2 <- FALSE
        }
      } else {
        available_d1 <- available_d1
        available_d2 <- available_d2
      }
    } else if(is.null(sigma_model)) {
      available_d1 <- available_d1
      available_d2 <- available_d2
    }
    if(!is.null(model$model_info[['model_deriv']])) {
      if(!model$model_info[['model_deriv']]) {
        available_d1 <- FALSE
        available_d2 <- FALSE
      }
    }
    if(!is.null(sigma_model)) {
      if(sigma_model == "basic") {
        if(sigma_model_is_ba_set_d0_as_d1) {
          out[['sigma_model_is_ba_set_d0_as_d1_funs']] <- sigma_model_is_ba_set_d0_as_d1_funs
          out[['sigma_model_is_ba_set_d0_as_d1_val']] <- 0
        }
        out[['sigma_model_is_ba_set_d0_as_d1']] <- sigma_model_is_ba_set_d0_as_d1
      }
    }
    out[['sigma_model']] <- sigma_model
    out[['available_d0']] <- available_d0
    out[['available_d1']] <- available_d1
    out[['available_d2']] <- available_d2
  }
  return(out)
}




#' An internal function to get arguments from the call
#'
#' @param x A character string
#' @param gsubit A character string specifying the part to be replaced.
#' @param gsubby A character string specifying the replacement
#' @param pasteit A logical indicating whether \code{gsubit} will be pasted
#'   before \code{gsubby}
#' @param fixed A logical indicating whether the fixed option
#' 
#' @return A character string.
#' @keywords internal
#' @noRd
#'
x_gsubit_gsubby <- function(x, 
                            gsubit, 
                            gsubby,
                            pasteit = TRUE, 
                            fixed = FALSE) {
  if(pasteit) gsubby <- paste0(gsubit, gsubby)
  gsub(gsubit, gsubby, x, fixed = fixed)
}





#' An internal function to get arguments from the call
#'
#' @param str A character string
#' @param x A character string specifying the part to be looked for
#' @param allowed_left A character string specifying the allowed valid character
#' on the left side. Typically period and an underscore 
#' @param allowed_right A character string specifying the allowed valid character
#' on the right side. Typically period and an underscore 
#' 
#' @return A character string.
#' @keywords internal
#' @noRd
#'
check_if_varname_exact <- function(str, 
                            x, 
                            allowed_left = "._", 
                            allowed_right = "._") {
  make_check_left  <- paste0("(^|[^[:alnum:]", allowed_left, "])")
  make_check_right <- paste0("($|[^[:alnum:]", allowed_right, "])")
  patxsi_not <- paste0(make_check_left, x, make_check_right)
  if(grepl(patxsi_not, str, fixed = FALSE)) {
    stop2c("Predictor ", collapse_comma(x), "",
         " has already been used for modelling the 'mu'",
         "\n  ", 
         "The same predictor can not be used for both 'mu' and 'sigma'",
         "\n  ", 
         "This is because predictors for 'mu' and 'sigma' often require different",
         "\n  ", 
         "transformations which is not possible with single a common predictor",
         "\n  ", 
         # "Please check and edit the relevent part of the code shown below:",
         # "\n  ", 
         # str, 
         "\n  ",
         "Please define a new predictor which could be the same earlier", 
         "\n  ",
         "predictor ", collapse_comma(x), " but renamed as ", 
         collapse_comma(paste0(x, "2")), 
         ", or any other name"
         )
  }
}



#' An internal function to get arguments from the call
#'
#' @param x A character string
#' @param patterns A character string specifying patterns to be replaced 
#' @param replacements A character string specifying the replacements
#' 
#' @return A character string.
#' @keywords internal
#' @noRd
#'
multi_gsub_exact <- function(x, patterns, replacements) {
  patterns <- paste0("\\b", patterns, "\\b")
  Reduce(
    function(s, i) gsub(patterns[i], replacements[i], s),
    seq_along(patterns),
    init = x
  )
}



#' An internal function to replace T/F with full TRUE/FALSE - string
#' 
#' @description
#' Operates at at formula string level i.e., deparse(formula)
#' 
#' @param x A character string specifying \code{T/F} to be replaced by full.
#' 
#' @return A character string.
#' @keywords internal
#' @noRd
#'
replace_t_f_to_full <- function(x) {
  x <- gsub("T,", "TRUE,"  , x, fixed = TRUE)
  x <- gsub("T)", "TRUE)"  , x, fixed = TRUE)
  x <- gsub("F,", "FALSE," , x, fixed = TRUE)
  x <- gsub("F)", "FALSE)" , x, fixed = TRUE)
  x
}



#' An internal function to replace T/F with full TRUE/FALSE - formula
#' 
#' @description
#' Operates at at formula level i.e., class formula
#' 
#' @param x A formula object
#' 
#' @return A formula object
#' @keywords internal
#' @noRd
#'
replace_T_preserve_env <- function(x) {
  formula_obj <- x
  # Preserve the original environment
  original_env <- environment(formula_obj)
  formula_str <- deparse(formula_obj)
  formula_str <- gsub("\\bT\\b", "TRUE", formula_str)
  formula_str <- gsub("\\bF\\b", "FALSE", formula_str)
  new_formula <- as.formula(paste(formula_str, collapse = ""))
  environment(new_formula) <- original_env
  return(new_formula)
}




#' An internal function to get arguments from the call
#'
#' @param str A character string
#' @param x A character string specifying \code{T/F} to be replaced by full. if
#' \code{NULL}, both \code{T/F} checked 
#' @param what A character string specifying the full form \code{TRUE/FALSE}. if
#' \code{NULL}, both \code{TRUE/FALSE} used as full. Note that length of what 
#' should match the length of x 
#' @param allowed_left A character string specifying the allowed valid character
#' on the left side. Typically it is \code{'='} that is used after argument 
#' @param allowed_right A character string specifying the allowed valid character
#' on the right side. Typically empty string
#' 
#' @return A character string.
#' @keywords internal
#' @noRd
#'
check_and_replace_sort_to_full <- function(str, 
                                   x = NULL, 
                                   what = NULL, 
                                   allowed_left = NULL, 
                                   allowed_right = NULL) {
  if(!is.null(what)) {
    if(is.null(x)) stop2c("'x' must be specified when 'what' is not NULL")
    if(length(what) != length(what)) stop2c("lengths of 'x' and 'what' must be same")
  }
  if(is.null(x)) {
    x <- c("T", "F")
    if(is.null(what)) {
      what <- c("TRUE", "FALSE")
    }
  }
  if(is.null(allowed_left)) {
    allowed_left <- rep("", length(x))
  }
  if(is.null(allowed_right)) {
    allowed_right <- rep("", length(x))
  }
  make_check_left  <- allowed_left
  make_check_right <- allowed_right
  if(length(str) > 1) {
    for (i in 1:length(str)) {
      str_i <- str[i]
      for (i in 1:length(x)) {
        patTF <- paste0(make_check_left, x[i], make_check_right)
        if(x[i] %in% str) {
          str <- gsub(patTF, what[i], str, fixed = F)
        }
      }
    }
  }
  if(length(str) == 1) {
    for (i in 1:length(str)) {
      str_i <- str[i]
      for (i in 1:length(x)) {
        patTF <- paste0(make_check_left, x[i], make_check_right)
        if(grepl(patTF, str_i, fixed = FALSE)) {
          str <- gsub(x[i], what[i], str, fixed = F)
        }
      }
    }
  }
  return(str)
}



#' An internal function to check and set xvar for dpar sigma
#'
#' @param model An object of class bgmfit
#' @param newdata A data frame
#' @param dpar A character string 
#' @param resp A character string 
#' @param all A logical \code{TRUE} for \code{what = 'model'}, will return all 
#' models otherwise resp specific. For for \code{what = 'cov'},  \code{TRUE} 
#' will return all co variates (factor and numeric) otherwise resp specific
#' (factor and numeric)
#' @param cov A character string, \code{'numeric'} will return only numeric
#' co variate whereas \code{'factor'} will return only factor covariates. 
#' \code{'all'} will return both factor and numeric covaraites 
#' @param auto ignored
#' @param verbose print information
#' 
#' @return A character string
#' @keywords internal
#' @noRd
#'
get_sigmamodel_info <- function(model,
                           newdata = NULL,
                           dpar = NULL, 
                           resp = NULL, 
                           what = 'model',
                           cov = NULL, 
                           all = FALSE, 
                           verbose = FALSE) {
  if(is.null(newdata)) {
    newdata <- model$data
  }
  if(is.null(dpar)) {
    dpar <- "mu"
  }
  if(dpar == "mu") {
     return(NULL)
  }
  if(is.null(cov)) {
    cov <- "all"
  }
  if(is.null(newdata)) {
    newdata <- model$data
  }
  if(dpar == "sigma") {
    if(what == 'model') {
      if(all) {
        sigma_model_      <- paste0('sigmamodel', "_", all)
        sigma_model       <- model$model_info[[sigma_model_]]
        out <- sigma_model
      } else if(!all) {
        if (is.null(resp)) {
          resp_    <- resp
          revresp_ <- ""
        } else if (!is.null(resp)) {
          resp_    <- paste0(resp, "_")
          revresp_ <- paste0("_", resp)
        }
        sigma_model_      <- paste0('sigmamodel', revresp_)
        sigma_model       <- model$model_info[[sigma_model_]]
        out <- sigma_model
      }
      if(is.null(out)) {
        out <- 'conventional'
      }
    }
    if(what == 'cov' | what == 'covariate' | what == 'covariates') {
      if(all) {
        sigmacov_sigma_model_ <- paste0('sigmacov', 's')
        sigmacov_cov_vars     <- model$model_info[[sigmacov_sigma_model_]]
        sigmacov_factor_vars  <- names(newdata[sapply(newdata, is.factor)])
        sigmacov_numeric_vars <- names(newdata[sapply(newdata, is.numeric)])
        sigmacov_cov_factor_vars  <- intersect(sigmacov_cov_vars, 
                                               sigmacov_factor_vars)
        sigmacov_cov_numeric_vars <- intersect(sigmacov_cov_vars, 
                                               sigmacov_numeric_vars)
        if(cov == "numeric") {
          out <- sigmacov_cov_numeric_vars
        } else if(cov == "factor") {
          out <- sigmacov_cov_factor_vars
        } else if(cov == "all") {
          out <- c(sigmacov_cov_factor_vars, sigmacov_cov_numeric_vars)
        }
      } else if(!all) {
        if (is.null(resp)) {
          resp_    <- resp
          revresp_ <- ""
        } else if (!is.null(resp)) {
          resp_    <- paste0(resp, "_")
          revresp_ <- paste0("_", resp)
        }
        sigmacov_sigma_model_ <- paste0('sigmacov', revresp_)
        sigmacov_cov_vars     <- model$model_info[[sigmacov_sigma_model_]]
        sigmacov_factor_vars  <- names(newdata[sapply(newdata, is.factor)])
        sigmacov_numeric_vars <- names(newdata[sapply(newdata, is.numeric)])
        sigmacov_cov_factor_vars  <- intersect(sigmacov_cov_vars, 
                                               sigmacov_factor_vars)
        sigmacov_cov_numeric_vars <- intersect(sigmacov_cov_vars, 
                                               sigmacov_numeric_vars)
        if(cov == "numeric") {
          out <- sigmacov_cov_numeric_vars
        } else if(cov == "factor") {
          out <- sigmacov_cov_factor_vars
        } else if(cov == "all") {
          out <- c(sigmacov_cov_factor_vars, sigmacov_cov_numeric_vars)
        }
      }
    } 
  } 
  return(out)
}




#' An internal function to check and set xvar for dpar sigma
#'
#' @param model An object of class bgmfit
#' @param newdata A data frame
#' @param dpar A character string 
#' @param xvar A character string 
#' @param resp A character string 
#' @param dpar A character string 
#' @param auto Set \code{'xvar'} as \code{'sigmacov_cov_numeric_vars'} if it
#' is of length one.
#' @param verbose print information
#' 
#' @return A character string
#' @keywords internal
#' @noRd
#'
check_set_xvar_sigma <- function(model, 
                                 newdata = NULL,
                                 dpar = NULL, 
                                 xvar = NULL, 
                                 resp = NULL, 
                                 auto = TRUE,
                                 verbose = FALSE) {
  if(is.null(newdata)) {
    newdata <- model$data
  }
  if(is.null(dpar)) {
    dpar <- "mu"
  }
  if(dpar == "mu") {
    return(xvar)
  }
  set_sigma_xvar_as_mu_if <- c("varpower", 
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
  if(dpar == "sigma") {
    sigma_model <- get_sigmamodel_info(model,
                                       newdata = newdata,
                                       dpar = dpar, 
                                       resp = resp, 
                                       what = 'model',
                                       cov = NULL, 
                                       all = FALSE, 
                                       verbose = FALSE)
    if(!is.null(sigma_model)) {
      if(sigma_model %in% set_sigma_xvar_as_mu_if) {
        if(is.null(xvar)) {
          xvar <- get_basic_info(model = model, 
                                 dpar = dpar, 
                                 resp = resp, 
                                 auto = TRUE,
                                 what = 'xvar',
                                 component = "mu",
                                 verbose = verbose)
        }
        if(verbose) {
          create_msg_xvar_as_mu <- 
            paste0("For dpar = 'sigma', the 'xvar' was NULL. Now the 'xvar'",
                   "\n  ",
                   "for sigma model ", collapse_comma(sigma_model),
                   " has been set same as the 'xvar' for dpar = 'mu':",
                   "\n ", 
                   collapse_comma(xvar),
                   "\n ")
          message2c(create_msg_xvar_as_mu)
        }
        return(xvar)
      }
    }
    sigmacov_cov_factor_vars <- get_sigmamodel_info(model,
                                                   newdata = newdata,
                                                   dpar = dpar, 
                                                   resp = resp, 
                                                   what = 'cov',
                                                   cov = "factor", 
                                                   all = FALSE, 
                                                   verbose = FALSE)
    sigmacov_cov_numeric_vars <- get_sigmamodel_info(model,
                                                     newdata = newdata,
                                                     dpar = dpar, 
                                                     resp = resp, 
                                                     what = 'cov',
                                                     cov = "numeric", 
                                                     all = FALSE, 
                                                     verbose = FALSE)
    sigmacov_cov_numeric_vars_without_cons <- c()
    for (i in sigmacov_cov_numeric_vars) {
      if(min(newdata[[i]]) != max(newdata[[i]])) {
        sigmacov_cov_numeric_vars_without_cons <- 
          c(sigmacov_cov_numeric_vars_without_cons, i)
      }
    }
    sigmacov_cov_numeric_vars <- sigmacov_cov_numeric_vars_without_cons
    if(length(sigmacov_cov_factor_vars) == 0) {
      sigmacov_cov_factor_vars <- NULL
    }
    if(length(sigmacov_cov_numeric_vars) == 0) {
      sigmacov_cov_numeric_vars <- NULL
    }
    create_msg_xvar_null <- 
      paste0("For dpar = 'sigma', the 'xvar' should be specified",
             "\n  ",
             "The available options are:\n ", 
             collapse_comma(sigmacov_cov_numeric_vars))
    create_msg_xvar_not_null_and_one <- 
      paste0(" The 'xvar' required for plot is automatically set as: ", 
              collapse_comma(sigmacov_cov_numeric_vars),
              "\n ",
              "Note that this 'auto' option works only if there ",
              "\n ",
              "is one unique numeric variable in the sigma formula")
    create_msg_xvar_not_null_and_more_than_one <- 
      paste0(" The 'xvar' required for plot can not be set automatically because", 
             "\n ",
             " of more than one numeric variables used in the sigma formula ",
             collapse_comma(sigmacov_cov_numeric_vars),
             "\n ",
             "Note that the 'auto' option works only if there ",
             "\n ",
             "is one unique numeric variable in the sigma formula")
    create_msg_xvar_used <- 
      paste0("For dpar = 'sigma', you have specified ", 
             collapse_comma(xvar), 
             "\n  ",
             "as 'xvar' which is invalid",
             "\n  ",
             "The available options for 'xvar' are:\n ", 
             collapse_comma(sigmacov_cov_numeric_vars))
    if(is.null(xvar)) {
      if(!is.null(sigma_model)) {
        if(sigma_model != "ls") {
          if(!auto) {
            stop2c(create_msg_xvar_null)
          } else if(auto) {
            if(length(sigmacov_cov_numeric_vars) > 0) {
              if(length(sigmacov_cov_numeric_vars) == 1) {
                xvar <- sigmacov_cov_numeric_vars
                if(verbose) {
                  message2c(create_msg_xvar_not_null_and_one)
                }
              } else {
                stop2c(create_msg_xvar_not_null_and_more_than_one)
              }
            }
          }
        }
      }
    } else if(!is.null(xvar)) {
      for (i in xvar) {
        if(!i %in% sigmacov_cov_numeric_vars) {
          if(sigma_model == "basic") stop2c(create_msg_xvar_used)
        }
      }
    }
  } 
  if(is.null(xvar)) {
    return(invisible(NULL))
  } else {
    return(xvar)
  }
} 


#' An internal function to check and set xvar for dpar sigma
#'
#' @param cov A character vector of all co variates (factor and numeric)
#' @param newdata A data frame
#' @param strict_1 A logical (\code{TRUE}) that indicates to remove only 
#' if all \code{cov == 1}. Therefore \code{strict_1 = TRUE} will only remove
#' variables that were used a manual intercept variable.
#' @param verbose print information
#' 
#' @return A character vector 
#' @keywords internal
#' @noRd
#'
remove_cons_as_numeric_cov <- function(cov, 
                                       newdata, 
                                       strict_1 = TRUE,
                                       verbose = FALSE) {
  cov_vars     <- cov
  factor_vars  <- names(newdata[sapply(newdata, is.factor)])
  numeric_vars <- names(newdata[sapply(newdata, is.numeric)])
  cov_factor_vars  <- intersect(cov_vars, factor_vars)
  cov_numeric_vars <- intersect(cov_vars, numeric_vars)
  cov_numeric_vars_without_cons <- c()
  if(strict_1) {
    for (i in cov_numeric_vars) {
      if(!all(newdata[[i]] == 1)) {
        cov_numeric_vars_without_cons <- 
          c(cov_numeric_vars_without_cons, i)
      }
    }
  } 
  if(!strict_1) {
    for (i in cov_numeric_vars) {
      if(min(newdata[[i]]) != max(newdata[[i]])) {
        cov_numeric_vars_without_cons <- 
          c(cov_numeric_vars_without_cons, i)
      }
    }
  }
  removed_cons_vars <- setdiff(cov_numeric_vars, cov_numeric_vars_without_cons) 
  if(length(removed_cons_vars) > 0) {
    if(verbose) {
      message2c("Following covariate(s) removed because of constant values",
              "\n ",
              collapse_comma(removed_cons_vars))
    }
  }
  cov_numeric_vars <- cov_numeric_vars_without_cons
  if(length(cov_factor_vars) == 0) {
    cov_factor_vars <- NULL
  }
  if(length(cov_numeric_vars) == 0) {
    cov_numeric_vars <- NULL
  }
  cov_out <- c(cov_factor_vars, cov_numeric_vars)
  return(cov_out)
}



#' An internal function to check and set xvar for dpar sigma
#'
#' @param model An object of class bgmfit
#' @param newdata A data frame
#' @param dpar A character string 
#' @param xvar A character string 
#' @param resp A character string 
#' @param dpar A character string 
#' @param auto Set \code{'xvar'} as \code{'sigmacov_cov_numeric_vars'} if it
#' is of length one.
#' @param verbose print information
#' 
#' @return A character string
#' @keywords internal
#' @noRd
#'
check_set_transform_draws_sigma <- function(model, 
                                            newdata = NULL,
                                            dpar = NULL, 
                                            xvar = NULL, 
                                            resp = NULL, 
                                            auto = TRUE,
                                            transform_draws = NULL,
                                            itransform = NULL,
                                            verbose = FALSE) {
  if(is.null(newdata)) {
    newdata <- model$data
  }
  if(is.null(dpar)) {
    dpar <- "mu"
  }
  if(dpar == "mu") {
    return(xvar)
  }
  set_sigma_transform_draws_if <- c("varpower", 
                                    "varconstpower",
                                    "fittedpower", 
                                    "meanpower", 
                                    "residualpower")
  if(dpar == "sigma") {
    sigma_model <- get_sigmamodel_info(model,
                                       newdata = newdata,
                                       dpar = dpar, 
                                       resp = resp, 
                                       what = 'model',
                                       cov = NULL, 
                                       all = FALSE, 
                                       verbose = FALSE)
    set_transform_draws <- NULL
    if(!is.null(sigma_model)) {
      if(is.null(transform_draws)) {
        if(sigma_model %in% set_sigma_transform_draws_if) {
          set_transform_draws <- 'exp'
          if(verbose) {
            create_msg_transform_draws <- 
              paste0("For dpar = 'sigma', the 'transform_draws' was NULL. ",
                     "\n  ",
                     "Now for sigma model ", collapse_comma(sigma_model),
                     ", the 'transform_draws' is set as:",
                     "\n ", 
                     collapse_comma(set_transform_draws),
                     "\n ")
            message2c(create_msg_transform_draws)
            } 
        } 
      } 
    } 
  } 
  if(is.null(set_transform_draws)) {
    fun_eval <- NULL
  } else if(!is.null(set_transform_draws)) {
    if(set_transform_draws == "identity") {
      fun_eval <- function(x)x
    } else if(set_transform_draws == "log") {
      fun_eval <- function(x)log(x)
    } else if(set_transform_draws == "exp") {
      fun_eval <- function(x)exp(x)
    } else if(set_transform_draws == "sqrt") {
      fun_eval <- function(x)sqrt(x)
    }
  }
  if(!is.null(set_transform_draws)) {
    if(is.null(itransform)) {
      if(verbose) {
        message2c("Note that argument 'transform_draws' has bee set as",
                "\n ",
                deparse(fun_eval), " but 'itransform' is NULL. Thegerfore, no",
                "\n ",
                "automatic adhustment applied to the parametrs such as APGV",
                "\n ")
      }
    }
  }
  return(fun_eval)
} 



#' An internal function to check and set ifunx for dpar mu and sigma
#'
#' @param model An object of class bgmfit
#' @param dpar A character string 
#' @param resp A character string 
#' @param itransform A character string ora function
#' @param auto Set \code{'xvar'} as \code{'sigmacov_cov_numeric_vars'} if it
#' is of length one.
#' @param verbose print information
#' 
#' @return A character string
#' @keywords internal
#' @noRd
#'
check_set_fun_transform <- function(model, 
                                    which = NULL,
                                    dpar = NULL, 
                                    resp = NULL, 
                                    transform = NULL,
                                    auto = TRUE, 
                                    verbose = FALSE) {
  if(is.null(which)) {
    which <- "ixfuntransform2"
  }
  if(is.null(dpar)) {
    dpar <- "mu"
  }
  add_prefix_to_fun <- ""
  if(!is.null(dpar)) {
    if(dpar == "mu") {
      add_prefix_to_fun <- ""
    } else if(dpar == "sigma") {
      add_prefix_to_fun <- "sigma"
    }
  } 
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  setfun <- paste0(which, resp_rev_)
  setfun <- setfunname <- paste0(add_prefix_to_fun, setfun)
  setfun <- model$model_info[[setfun]]
  if(is.null(setfun)) {
    was_null <- TRUE
  } else {
    was_null <- FALSE
  }
  allowedstrfun <- c("identity", "log", "exp")
  if(is.null(setfun)) {
    if(!is.null(transform)) {
      if(is.character(transform)) {
        transform <- paste(gsub_space(transform), collapse = "")
        if(transform != "") {
          if(grepl("function(", transform, fixed = T)) {
            setfun <-  ept(transform)
          } else if(transform == "identity") {
            setfun <- function(x)x
          } else if(transform == "log") {
            setfun <- function(x)log(x)
          } else if(transform == "exp") {
            setfun <- function(x)exp(x)
          } else {
            stop2c("The argument 'transform' must be either ", 
                 collapse_comma(allowedstrfun),
                 "\n ",
                 " or else a valid function such as function(x) x")
          }
        }
      } else if(is.function(transform)) {
        setfun <- transform
      }
    } 
  }
  
  setfunas_char <- 'setfun'
  if(is.null(setfun)) {
    if(is.null(transform)) {
      if(auto) {
        setfun <- function(x)x
        if(verbose) {
          message2c(paste0("For dpar = '", dpar, "', the ", 
                         setfunas_char, 
                         " was 'NULL', now set as function(x)x"))
        }
      } 
    } 
    if(!is.null(transform)) {
      setfun <- transform
    } 
  } 
  
  if(is.character(setfun)) {
    assign('setfun', ept(setfun))
  } else {
    assign('setfun', eval(setfun))
  }
  setfun <- as.function(setfun)
  out <- list(setfun = setfun, setfunname = setfunname, was_null = was_null)
  return(out)
}



#' An internal function to extract basic information from the model_info
#'
#' @param model An object of class bgmfit
#' @param dpar A character string 
#' @param resp A character string 
#' @param what A character string indicating what to extrcat. Options are
#'   \code{'xvar'}, \code{'yvar'}, \code{'idvar'}, \code{'cov_vars'},
#'   \code{'uvarby'}.
#' @param component A character string indicating \code{'what'} should be for
#'   \code{'mu'} or  \code{'sigma'}. For example if \code{'mu'}, then
#'   \code{'xvar'} will be extracted and if \code{'sigma'}, then
#'   \code{'sigmaxvar'} will be extracted.
#' 
#' @param auto ignored
#' @param verbose print information
#' 
#' @return A character string
#' @keywords internal
#' @noRd
#'
get_basic_info <- function(model = model, 
                           dpar = dpar, 
                           resp = resp, 
                           what = 'xvar',
                           component = "mu",
                           auto = TRUE,
                           verbose = verbose) {
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  if(is.null(dpar)) {
    dpar <- "mu"
  }
  levels_id <- NULL
  idvar <- NULL
  sigmalevels_id <- NULL
  sigmaidvar <- NULL
  validate_response(model, resp)
  list_c <- list()
  xvar_ <- paste0('xvar', resp_rev_)
  yvar_ <- paste0('yvar', resp_rev_)
  groupvar_ <- paste0('groupvar', resp_rev_)
  hierarchical_ <- paste0('hierarchical', resp_rev_)
  sigmaxvar_ <- paste0('sigma', xvar_)
  sigmayvar_ <- paste0('sigma', yvar_)
  sigmagroupvar_ <- paste0('sigma', groupvar_)
  sigmahierarchical_ <- paste0('sigma', hierarchical_)
  xvar <- model$model_info[[xvar_]]
  yvar <- model$model_info[[yvar_]]
  sigmaxvar <- model$model_info[[sigmaxvar_]]
  sigmayvar <- model$model_info[[sigmayvar_]]
  if(is.null(levels_id) & is.null(idvar)) {
    idvar <- model$model_info[[groupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      idvar <- model$model_info[[hierarchical_]]
    }
  } else if (!is.null(levels_id)) {
    idvar <- levels_id
  } else if (!is.null(idvar)) {
    idvar <- idvar
  }
  if(is.null(idvar)) {
    if(!is.null(model$model_info[['idvars']])) {
      idvar <- model$model_info[['idvars']]
    }
  }
  if(is.null(sigmalevels_id) & is.null(sigmaidvar)) {
    sigmaidvar <- model$model_info[[sigmagroupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      sigmaidvar <- model$model_info[[hierarchical_]]
    }
  } else if (!is.null(sigmalevels_id)) {
    sigmaidvar <- sigmalevels_id
  } else if (!is.null(sigmaidvar)) {
    sigmaidvar <- sigmaidvar
  }
  if(is.null(sigmaidvar)) {
    if(!is.null(model$model_info[['sigmaidvars']])) {
      sigmaidvar <- model$model_info[['sigmaidvars']]
    }
  }
  cov_          <- paste0('cov', resp_rev_)
  sigmacov_     <- paste0('sigma', cov_)
  cov_vars      <-  model$model_info[[cov_]]
  sigmacov_vars <-  model$model_info[[sigmacov_]]
  cov_vars      <- cov_vars[!is.na(cov_vars)]
  sigmacov_vars <- sigmacov_vars[!is.na(sigmacov_vars)]
  if(length(cov_vars) == 0) cov_vars <- NULL
  if(length(sigmacov_vars) == 0) sigmacov_vars <- NULL
  if (!is.null(cov_vars)) {
    cov_vars <- covars_extrcation(cov_vars)
  }
  if (!is.null(sigmacov_vars)) {
    sigmacov_vars <- covars_extrcation(sigmacov_vars)
  }
  uvarby     <- model$model_info$univariate_by$by
  if(is.null(sigmaidvar)) {
    sigmaidvar <- idvar
  }
  out <- list()
  if(component == "mu") {
    if('xvar' %in% what)     out[[what]] <- xvar
    if('yvar' %in% what)     out[[what]] <- yvar
    if('idvar' %in% what)    out[[what]] <- idvar
    if('cov_vars' %in% what) out[[what]] <- cov_vars
    if('uvarby' %in% what)   out[[what]] <- uvarby
  } else if(component == "sigma") {
    if('xvar' %in% what)     out[[what]] <- sigmaxvar
    if('yvar' %in% what)     out[[what]] <- sigmayvar
    if('idvar' %in% what)    out[[what]] <- sigmaidvar
    if('cov_vars' %in% what) out[[what]] <- sigmacov_vars
  }
  if(length(out) == 1) {
    out <- out[[1]]
  }
  return(out)
} 


#' An internal function to extract basic information from the model_info
#' 
#' @details
#' Called in \code{'fitted_draws'}, \code{'growthparameters'}, 
#' \code{'plot_curves'} and \code{'plot_curves'}
#'  check by following. Note that need to call \code{"library(NCmisc)"} before
#' \code{"bsitar:::find_function_used_in_R_files('bsitar', 
#' 'set_sigma_grid_newdata')"}
#' 
#'
#' @param model An object of class bgmfit
#' @param newdata A data frame
#' @param dpar A character string 
#' @param resp A character string 
#' @param idvar A character string 
#' @param xvar A character string
#' @param difx A character string
#' @param difx_asit A character string
#' @param auto ignored
#' @param xrange A character string (\code{'range'} or \code{'minmax'}) or a
#'   numeric vector of length two.
#' @param grid_add A vector of variable names that should be added to the 
#' grid. Note that the first row of these variables is set for the entire 
#' variable.
#' 
#' @param verbose print information
#' 
#' @inheritParams marginaleffects::datagrid
#' 
#' @return A character string
#' @keywords internal
#' @noRd
#'
set_manual_datagrid <- function(model, 
                                newdata = NULL,
                                resp = NULL, 
                                dpar = NULL, 
                                idvar = NULL,
                                xvar = NULL,
                                difx = NULL,
                                difx_asit = FALSE,
                                auto = TRUE,
                                xrange = NULL,
                                length.out = NULL,
                                grid_add = NULL,
                                grid_type = NULL,
                                by = NULL,
                                FUN = NULL,
                                FUN_character = NULL,
                                FUN_factor = NULL,
                                FUN_logical = NULL,
                                FUN_numeric = NULL,
                                FUN_integer = NULL,
                                FUN_binary = NULL,
                                FUN_other = NULL,
                                verbose = FALSE) {
  
  # need to load 'NCmisc'
  ept("library(NCmisc)")
  
  xvar_temp <- idvar_temp <- NULL;
  uvarby <- model$model_info$univariate_by$by
  if(is.null(uvarby)) {
    uvarby <- NA
  }
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  if(is.null(dpar)) {
    dpar <- "sigma"
  }
  if(is.null(newdata)) {
    newdata <- model$model_info$bgmfit.data
  }
  if(is.null(xrange)) {
    if(!is.null(length.out)) {
      stop2c("For 'length.out' to take effect, the 'xrange' should be specified.
             The xrange must be a character string 'range' or 'minmax', or 
             else a numeric vector of lenght two")
    }
  }
  if(is.null(length.out)) {
    length.out <- 10
    if(verbose) message2c("The default 'length.out' is set as 10")
  }
  if(is.null(difx)) {
    difx_range <- NULL
    difx_name  <- NULL
    difx_val   <- NULL
  }
  if(!is.null(difx)) {
    if(is.character(difx)) {
      difx_name <- difx
      difx_val <- newdata[[difx_name]]
      if(is.null(difx_val)) {
        stop2c("The 'difx' variable ", collapse_comma(difx_name),
             " is not available in the dataframe")
      }
      if(length(unique(difx_val)) < 2) {
        stop2c("The 'difx' variable ", collapse_comma(difx_name),
             "\n ",
             "should have at least of a unique length 2",
             "\n ",
             "Currently, the unique elements are:",
             unique(difx_val))
      }
      difx_range <- range(difx_val)
    } else if(!is.character(difx)) {
      difx_name <- 'difx' # any hypotherical names
      if(difx_asit) {
        stop2c("To pass 'difx' variable as it is to the grid", 
             " it must a variable in the dataframe")
      }
      if(!is.numeric(difx)) {
        stop2c("The 'difx' variable", 
             " should be either a character, or a numeric vector of length 2")
      } else if(is.numeric(difx)) {
        if(length(difx) < 2) {
          stop2c("The 'difx' variable", 
               " should be either a character, or a numeric vector of length 2")
        } else if(length(difx) == 2) {
          difx_range <- difx
        } else if(length(difx) > 2) {
          difx_range <- range(difx)
        }
      } 
    } 
  } 
  if(is.null(idvar)) {
    idvar <- get_basic_info(model = model, 
                            dpar = dpar, 
                            resp = resp, 
                            auto = TRUE,
                            what = 'idvar',
                            component = "sigma",
                            verbose = verbose)
    idvar <- idvar[1]
  }
  if(is.null(xvar)) {
    xvar <- get_basic_info(model = model, 
                           dpar = dpar, 
                           resp = resp, 
                           auto = TRUE,
                           what = 'xvar',
                           component = "sigma",
                           verbose = verbose)
  }
  if(is.null(xvar)) {
    xvar <- check_set_xvar_sigma(model = model, 
                                 dpar = dpar, 
                                 xvar = xvar, 
                                 resp = resp, 
                                 auto = TRUE,
                                 verbose = verbose)
  }
  if(dpar == "sigma") {
    if(is.null(xvar) | is.na(xvar)) {
      msg_xvar_not_in_data <-
        paste0("For dpar = 'sigma', the 'xvar' is required but is not",
               "\n  ",
               "specified.",
               "\n  ",
               "Please see the documentation and specify 'xvar' argument")
      stop2c(msg_xvar_not_in_data)
    }
  }
  if(is.null(xrange)) {
    setxvarvec <- newdata[[xvar]]
  } else if(!is.null(xrange)) {
    if(is.character(xrange)) {
      if(xrange == 'range') {
        xrange <- range(newdata[[xvar]])
      } else if(xrange == 'minmax') {
        xrange <- range(newdata[[xvar]])
      }
    } else if(is.numeric(xrange)) {
      if(length(xrange) != 2) {
        stop2c("'xrange' must be a vector of lenght 2")
      }
    } else {
      stop2c("The xrange must be a character string 'range' or 'minmax', or 
             else must be a numeric vector of lenght two")
    }
    setxvarvec <- seq.int(xrange[1], xrange[2], length.out = length.out)
  }
  if(is.null(grid_type)) {
    grid_type <- "mean_or_mode" # grid_type <- "dataframe"
  }
  if(!is.factor(newdata[[idvar]])) {
    newdata[[idvar]] <- as.factor(newdata[[idvar]])
    if(verbose) {
      message2c("The ", idvar, 
              " used in 'mapderivqr' has been converted to 'as.factor()'")
    }
  }
  if(!model$test_mode) {
    unlock_replace_bind(package = "insight", what = "get_data",
                        replacement = custom_get_data.brmsfit, ept_str = T)
    if(verbose) {
      message2c("As model[['test_mode']] = FALSE, the full data are extracted 
                via insight::get_data() using custom_get_data.brmsfit. To 
                override this behavior, set model[['test_mode']] = TRUE.",
              pad_before = "\n", 
              pad_after = "\n")
    }
  } 
  grid_args <- list()
  if(grid_type == 'dataframe') {
    grid_args[[xvar]]          <- NULL
    grid_args[[idvar]]         <- NULL
    if(is.null(FUN)) FUN       <- function(x)x
    if(verbose) {
      message2c("For grid_type = 'dataframe', the 'xvar' and 'idvar' arguments
                for datagrid have been set as 'NULL', and the 'FUN' is defined
                as 'function(x)x'. This ensures that all variables have the 
                same lenght. Note that, user need to set mean/mode etc for each
                varibale in the data passed as newdata argument. For 
                grid_type = 'dataframe', the datagrid will simply combine the
                columns and set the needed attributes for the marginaleffects")
    }
  } else if(grid_type != 'dataframe') {
    grid_args[[xvar]]          <- setxvarvec
    grid_args[[idvar]]         <- levels(newdata[[idvar]])
    FUN                        <- FUN
  }
  grid_args[['model']]         <- model
  grid_args[['newdata']]       <- newdata
  grid_args[['grid_type']]     <- grid_type
  grid_args[['FUN']]           <- FUN 
  grid_args[['FUN_character']] <- FUN_character
  grid_args[['FUN_factor']]    <- FUN_factor
  grid_args[['FUN_logical']]   <- FUN_logical
  grid_args[['FUN_numeric']]   <- FUN_numeric
  grid_args[['FUN_integer']]   <- FUN_integer
  grid_args[['FUN_binary']]    <- FUN_binary
  grid_args[['FUN_other']]     <- FUN_other
  if(!is.na(uvarby)) {
    grid_args[['by']] <- c(uvarby, by)
  } else {
    grid_args[['by']] <- by
  }
  if(!is.null(grid_add)) {
    if(is.list(grid_add)) {
      for (i in grid_add) {
        grid_args <- c(grid_args, i)
      }
      if(verbose) {
        message2c(" Adding following to the grid via grid_add:\n ", 
                collapse_comma(names(grid_add)))
      }
    } else if(!is.list(grid_add)) {
      for (i in grid_add) {
        grid_args[[i]] <- newdata[[i]]
      } 
      if(verbose) {
        message2c("Adding following to the grid via grid_add:\n ", 
                collapse_comma(grid_add))
      }
    } 
  } 
  newdata_all        <- newdata
  newdata_names_all  <- names(newdata_all)
  newdata            <- do.call(marginaleffects::datagrid, grid_args)
  # newdata <- newdata %>% droplevels()
  # print(str(newdata))
  # stop()
  newdata_names_grid <- names(newdata)
  missing_names_grid <- setdiff(newdata_names_all, newdata_names_grid)
  for (i in missing_names_grid) {
    if(!is.null(newdata_all[[i]])) {
      newdata[[i]] <- newdata_all[[i]][1]
    }
  }
  if(length(missing_names_grid) > 0) {
    if(verbose) {
      message2c("Adding following to the grid as first only only:\n ", 
              collapse_comma(missing_names_grid))
    }
  }
  if(!is.null(difx)) {
    if(difx_asit) {
      difx_val <- difx_val
      if(verbose) {
        message2c("The 'difx' variable ", collapse_comma(difx_name),
                " added to the grid") 
      }
    } else if(!difx_asit) {
      if(verbose) {
        message2c("The 'difx' variable ", collapse_comma(difx_name),
                " added to the grid with the following range: ",
                "\n ",
                difx_range) 
      }
    }
  }
  sortxxxxzzz <- NULL;
  newdata <- newdata %>% dplyr::mutate(sortxxxxzzz = dplyr::row_number()) 
  idvar_xvar <- c(idvar, xvar)
  if(!is.null(difx)) {
    newdata <- newdata %>% 
      dplyr::arrange(!!as.name(dplyr::all_of(idvar_xvar))) %>% 
      dplyr::group_by(!!as.name(idvar)) %>%
      dplyr::mutate(
        !! as.name (difx_name) := seq.int(from = difx_range[1], 
                                          to = difx_range[2], 
                                          length.out = dplyr::n())
      ) %>%
      dplyr::ungroup()
  }
  newdata <- newdata %>% dplyr::arrange(sortxxxxzzz) %>% 
    dplyr::select(-sortxxxxzzz)
  if(!is.null(difx)) {
    attr(newdata, 'difx') <- difx_name
  }
  attr(newdata, 'xvar_for_sigma_model_basic')  <- xvar # idvar
  attr(newdata, 'idvar_for_sigma_model_basic') <- idvar # xvar
  return(newdata)
}




#' An internal function to extract co variates
#'
#' @param str A character string
#' 
#' @return A character string
#' @keywords internal
#' @noRd
#'
covars_extrcation <- function(str) {
  str <- gsub("[[:space:]]", "", str)
  for (ci in c("*", "+", ":")) {
    str <- gsub(ci, ' ', str, fixed = T)
  }
  str <- strsplit(str, " ")[[1]]
  str
}



#' An internal function to strip attributes
#'
#' @param vec A data frame
#' @param keep A character string or a vector
#' 
#' @return data frame
#' @keywords internal
#' @noRd
#'
attrstrip <- function(vec, keep = 'class'){
  if(is.null(keep))return(vec)
  att_names <- names(attributes(vec))
  att_names <- att_names[!att_names %in% keep]
  if(length(att_names) == 0)return(vec)
  for(i in att_names){
    attr(vec, i) <- NULL
  }
  vec
}



#' An internal function to check the minimum version of the package
#'
#' @param ipts An integer
#' @param nipts An integer specifying the default ipts
#' @param verbose A logical (default \code{FALSE}) to check 
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'
set_for_check_ipts <- function(ipts, nipts = 50, dpar = 'mu', verbose = FALSE) {
  if(is.logical(ipts)) {
    if(!ipts) {
      out <- NULL
    } else if(ipts) {
      out <- nipts
    }
  } else if(!is.logical(ipts)) {
    if(!is.null(ipts)) {
      out <- ipts
    } else if(is.null(ipts)) {
      if(dpar == "mu") {
        if(verbose) {
          message2c("Argument ipts has been set to ", nipts,
                    pad_before = "\n", 
                    pad_after = "\n")
        }
        out <- nipts
      } else if(dpar == "sigma") {
        out <- ipts
      }
    }
  }
  return(out)
}



#' An internal function to check the minimum version of the package
#'
#' @param ipts An integer
#' @param nipts An integer specifying the default ipts
#' @param check_fun A logical (default \code{FALSE})
#' @param available_d1 A logical (default \code{FALSE})
#' @param xcall A logical (default \code{NULL})
#' @param verbose A logical (default \code{FALSE})
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'
check_ipts <- function(ipts = NULL, 
                       nipts = NULL, 
                       check_fun  = FALSE,
                       available_d1  = FALSE,
                       xcall = NULL,
                       verbose = FALSE) {
  if(is.null(nipts)) {
    niptsavailable_d1_T <- nipts <- 20
    niptsavailable_d1_F <- nipts <- 50
  } else {
    niptsavailable_d1_T <- nipts
    niptsavailable_d1_F <- nipts
  }
  if(is.null(ipts)) {
    if(check_fun) {
      if(available_d1) {
        if(is.null(xcall)) {
          ipts <- niptsavailable_d1_T
        } else if(!is.null(xcall)) {
          if(xcall == 'fitted_draws') ipts <- niptsavailable_d1_T
        }
      } else if(!available_d1) {
        if(is.null(xcall)) {
          ipts <- niptsavailable_d1_F
        } else if(!is.null(xcall)) {
          if(xcall == 'fitted_draws') ipts <- niptsavailable_d1_F
        }
      }
    }
    if(verbose) {
      message2c("Argument ipts has been set to ", nipts,
                pad_before = "\n", 
                pad_after = "\n")
    }
  }
  return(ipts)
}




#' An internal function to find all instances of function used in R files
#'
#' @param package_name A character string
#' @param function_to_find A character string
#' @param print A logical
#' 
#' @return A character vector
#'
#' @return A character vector
#' @noRd
#'
find_function_used_in_R_files <- function(package_name, 
                                          function_to_find, 
                                          print = FALSE) {
  list.functions.in.file <- NULL;
  package_path <- find.package(package_name)
  r_files_path <- file.path(package_path, "R")
  r_files <- list.files(r_files_path, pattern = "\\.R$", full.names = TRUE)
  files_with_function <- character(0)
  for (file in r_files) {
    functions_in_file <- tryCatch({
      list.functions.in.file(file)
    }, error = function(e) {
      message2c(paste("Error processing file:", file, "-", e$message))
      NULL
    })
    if (!is.null(functions_in_file)) {
      for (pkg_fns in functions_in_file) {
        if (function_to_find %in% pkg_fns) {
          files_with_function <- c(files_with_function, basename(file))
          break # Found in this file, move to the next
        }
      }
    }
  }
  if(print) {
    if (length(files_with_function) > 0) {
      cat(paste0("The function '", function_to_find, "' is used in the following R files in the '", package_name, "' package:\n"))
      print(unique(files_with_function))
    } else {
      cat(paste0("The function '", function_to_find, "' was not found in any R files in the '", package_name, "' package.\n"))
    }
  } 
  return(unique(files_with_function))
} 


#' Find R files containing a specific pattern
#'
#' @description
#' Searches for a given text pattern across all `.R` files in a package
#' directory (including tests and subdirectories). Returns file paths of matches
#' and prints basenames for convenience.
#'
#' @param path Character string. Path to package root directory. Defaults to
#'   `"."` (current directory). If `NULL`, also set to `"."`.
#' @param folder Character string or `NULL`. Specific subdirectory to search
#'   (e.g., `"R"`, `"tests"`). If `NULL` (default), searches all `.R` files
#'   recursively across entire path.
#' @param pattern Character string. Exact text string to search for in R files.
#'
#' @return Character vector of full paths to matching `.R` files (invisibly).
#'   Prints basenames of matches or informative message if no matches found.
#'
#' @details 
#' Uses `list.files(recursive = TRUE)` to find all `.R` files, then
#' [readLines()] and `grepl(fixed = TRUE)` for efficient exact string matching.
#' `warn = FALSE` in `readLines()` suppresses warnings for encoding/large files.
#'
#' @examples
#' # Search current package for "abcbd" across all folders
#' find_r_files_with_pattern(pattern = "abcbd")
#'
#' # Search only R/ directory
#' find_r_files_with_pattern(folder = "R", pattern = "abcbd")
#'
#' # Search tests/ in specific package path
#' find_r_files_with_pattern(
#'   path = "/path/to/my-package", 
#'   folder = "tests", 
#'   pattern = "test_that"
#' )
#'
#' @keywords internal
#' @noRd
#' 
find_r_files_with_pattern <- function(path = ".", folder = NULL, pattern) {
  if (is.null(path)) path <- "."
  search_path <- if (is.null(folder)) {
    path
  } else {
    file.path(path, folder)
  }
  r_files <- list.files(
    path = search_path,
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE
  )
  matching_files <- r_files[vapply(r_files, function(file) {
    lines <- readLines(file, warn = FALSE)
    any(grepl(pattern, lines, fixed = TRUE))
  }, logical(1))]
  
  if (length(matching_files) > 0) {
    cat("Files containing '", pattern, "' in ", 
        if (is.null(folder)) "all folders" else folder, ":\n", sep = "")
    print(basename(matching_files))
  } else {
    cat("No .R files contain '", pattern, "'.\n", sep = "")
  }
  invisible(matching_files)
}



#' Find (and optionally replace) patterns in R files
#'
#' Searches for a text pattern across \code{.R} files in a package directory, 
#' optionally replaces matches with new text. Handles symbols, quoted strings 
#' (\code{"abcd"}, \code{'abcd'}).
#'
#' @param path Character. Package root path (default \code{"."}). 
#'   \code{NULL} sets to \code{"."}.
#' @param old_pattern Character. Pattern/symbol to find (e.g., \code{"abcbd"}, 
#'   \code{myfun}).
#' @param new_pattern Character. Replacement text (required if \code{replace =
#'   TRUE}).
#' @param folder Character or \code{NULL}. Subdirectory (e.g., \code{"R"}, 
#'   \code{"tests"}). \code{NULL} searches all folders.
#' @param replace Logical. If \code{TRUE}, replace matches in files (creates 
#'   \code{.bak} backups).
#' @param verbose Logical. If \code{TRUE}, print number of replacements in files
#'
#' @return Character vector of affected full paths (invisibly).
#'
#' @details 
#' - Find mode (\code{replace = FALSE}): Uses \code{grepl(fixed = TRUE)} 
#'   for exact search.
#' - Replace mode: Parses lines with \code{parse(text = ...)} for symbols; 
#'   uses \code{gsub(fixed = TRUE)} for strings. Creates backups.
#' - Backups saved as \code{filename.R.bak}.
#' - Test with \code{replace = FALSE} before replacing.
#'
#' @examples
#' # Find only
#' # get_predictions -> get_predictions
#' # get_comparisons -> get_comparisons
#' # get_growthparameters -> get_growthparameters
#' 
#' find_replace_r_files_with_pattern("abcbd")
#' 
#' # Replace in R/ folder
#' find_replace_r_files_with_pattern(
#'   folder = "R", 
#'   old_pattern = "abcbd", 
#'   new_pattern = "newfun", 
#'   replace = TRUE
#' )
#'
#' @keywords internal
#' @noRd
#' 
find_replace_r_files_with_pattern <- function(path = ".", 
                                              old_pattern,
                                              new_pattern = NULL,
                                              folder = NULL,
                                              replace = FALSE,
                                              verbose = TRUE) {
  if (is.null(path)) path <- "."
  if (replace && is.null(new_pattern)) {
    stop("new_pattern required for replace = TRUE")
  }
  search_path <- if (is.null(folder)) path else file.path(path, folder)
  r_files <- list.files(
    path = search_path,
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE
  )
  matching_files <- character(0)
  for (file in r_files) {
    lines <- readLines(file, warn = FALSE)
    has_match <- any(grepl(old_pattern, lines, fixed = TRUE))
    if (!has_match) next
    matching_files <- c(matching_files, file)
    if (replace) {
      file_bak <- paste0(file, ".bak")
      file.copy(file, file_bak, overwrite = TRUE)
      new_lines <- gsub(old_pattern, new_pattern, lines, fixed = TRUE)
      n_old <- sum(grepl(old_pattern, lines, fixed = TRUE))
      n_new <- sum(grepl(old_pattern, new_lines, fixed = TRUE))
      cat(basename(file), ": replaced ", n_old - n_new,
          " occurrences (may be 0 if patterns overlap).\n", sep = "")
      writeLines(new_lines, file)
    }
  }
  if (length(matching_files) > 0) {
    cat("Files with '", old_pattern, "':\n", sep = "")
    print(basename(matching_files))
  } else {
    cat("No matches for '", old_pattern, "'.\n", sep = "")
  }
  return(invisible(matching_files))
}




#' Function to check if a string contains only letters
#' 
#' @details
#' This will check a word only and flag presence of any paranthesis, comma etc
#' 
#' @param x A string
#' 
#' @return A logical
#'
#' @keywords internal
#' @noRd
#' 
is_only_letters <- function(x) {
  grepl("^[[:alpha:]]+$", x)
}




#' Function to check if a string contains only letters
#' 
#' @details
#' This will get sigma method specified in nlf() via sigma_formula_manualsi
#' such as method=fitted
#' 
#' @param x A string
#' @param search A string to look for custom argument 
#' @param clean A logical to indicate removal of custom argument
#' 
#' @return A string
#'
#' @keywords internal
#' @noRd
#' 
get_nlf_custom_arg <- function(str, 
                               search, 
                               allowed_nlf_custom_arg = NULL,
                               clean = TRUE) {
  if(str == "fitted") {
    str <- "fittedpower"
  }
  if(str == "mean") {
    str <- "meanpower"
  }
  if(str == "residual") {
    str <- "residualpower"
  }
  search.o <- search
  search <- paste0(search, "=")
  if(!grepl(search, str)) {
    return(str)
  }
  
  if(search.o == "method") {
    method_nlf_custom_arg_msg <- c("basic (ba)",
                                   "varpower (vp)",
                                   "varconstpower (cp)",
                                   "varexp (ve)",
                                   # "fitted (fi)",
                                   "fittedz (fz)",
                                   "fittedpower (fp)",
                                   "fittedexp (fe)",
                                   # "mean (mi)",
                                   "meanpower (mp)",
                                   "meanexp (me)",
                                   # "residual (ri)",
                                   "residualpower (rp)",
                                   "residualexp (re)",
                                   "ls")
  }

  look_for_sigma_method_paran <- 
    replace_string_part(x = str,
                        start = search,
                        end = ")",
                        replace = "",
                        extract = T)
  look_for_sigma_method_comma <- 
    replace_string_part(x = str,
                        start = search,
                        end = ",",
                        replace = "",
                        extract = T)
  look_for_sigma_method_paran <- substr(look_for_sigma_method_paran, 1, 
                                        nchar(look_for_sigma_method_paran) - 1)
  look_for_sigma_method_comma <- substr(look_for_sigma_method_comma, 1, 
                                        nchar(look_for_sigma_method_comma) - 1)
  look_for_sigma_method_paran <- gsub(search, "", 
                                      look_for_sigma_method_paran, fixed = T)
  look_for_sigma_method_comma <- gsub(search, "", 
                                      look_for_sigma_method_comma, fixed = T)
  look_for_sigma_method_paran.o <- look_for_sigma_method_paran
  if(grepl("'", look_for_sigma_method_paran, fixed = T)) {
    look_for_sigma_method_paran <- gsub("'", "", 
                                        look_for_sigma_method_paran, fixed=T)
  }
  if(grepl(",", look_for_sigma_method_paran)) {
    look_for_sigma_method_paran <- 
      strsplit(look_for_sigma_method_paran, ",", fixed = TRUE)[[1]][1]
  }
  if(is_only_letters(look_for_sigma_method_paran)) {
    out <- look_for_sigma_method_paran
  } else if(is_only_letters(look_for_sigma_method_comma)) {
    out <- look_for_sigma_method_comma
  } else {
    out <- NULL
  }
  if(!is.null(out)) {
    if(!is.null(allowed_nlf_custom_arg)) {
      if(!out %in% allowed_nlf_custom_arg) {
        if(search.o == "method") {
          stop2c(paste0("The custom arg '", search.o, "' used in nlf() must",
                      " be one of the following (short form in parenthese):", 
                      "\n  ",
                      collapse_comma(method_nlf_custom_arg_msg),
                      "\n  ",
                      ". Note that 'ls' is location-scale model which does't have ",
                      "\n  ",
                      "any alternative (full) name.",
                      "\n  "))
        } else if(search.o == "prior") {
          stop2c(paste0("The custom arg '", search.o, "' used in nlf() ",
                      "must be one of the following:", 
                      "\n  ",
                      collapse_comma(allowed_nlf_custom_arg),
                      "\n  "))
        } 
      } 
    } 
  } 
  out.org <- out
  if(out == "no")  out <- "none"
  if(out == "ba")  out <- "basic"
  if(out == "fz")  out <- "fittedz"
  if(out == "fi")  out <- "fittedpower"
  if(out == "fp")  out <- "fittedpower"
  if(out == "fe")  out <- "fittedexp"
  if(out == "ve")  out <- "varexp"
  if(out == "vp")  out <- "varpower"
  if(out == "cp")  out <- "varconstpower"
  if(out == "rp")  out <- "residualpower"
  if(out == "re")  out <- "residualexp"
  if(out == "mp")  out <- "meanpower"
  if(out == "me")  out <- "meanexp"
  if(out == "ls")  out <- "ls"
  if(!clean) {
    return(out)
  }
  removeit <- paste0(",", search, "", collapse_comma(out.org))
  if(!is.null(out)) {
    str <- gsub(removeit, "", str, fixed = T)
    out <- c(str, out)
  }
  return(out)
}


#' Function to extract function names and code from a string
#' 
#' @details
#' This will get function names and code specified in nlf() via
#' sigma_formula_manualsi such 
#' code{"lf(z~1+splines2:::bsp(age)+(1++splines2::nsp(age)|55|gr(id,by=NULL)))"}
#' 
#' @param str A string
#' @param replace_ns Replace namespace \code{'::'} and \code{':::'} with 
#' \code{'_'} in the \code{str} (if \code{replace_ns = TRUE})
#' 
#' @return A list
#'
#' @keywords internal
#' @noRd
#' 
get_function_names_code_from_string <- function(str, replace_ns = TRUE) {
  token <- NULL;
  text <- NULL;
  parent <- NULL;
  . <- NULL;
  str <- paste0(gsub_space(str), collapse = ",")
  set_getParseData <- utils::getParseData(parse(text=str, keep.source=TRUE)) 
  packages_included <- set_getParseData %>% 
    dplyr::filter(token=="SYMBOL_PACKAGE") %>% dplyr::pull(text)
  packages_included <- unique(packages_included)
  insight::check_if_installed(packages_included, prompt = FALSE)
  functions_namespace_included_id <- set_getParseData %>% 
    dplyr::filter(token=="NS_GET" | token=="NS_GET_INT") %>% 
    dplyr::pull(parent)
  functions_namespace_included_c <- c()
  functions_namespace_attr_c <- c()
  functions_namespace_included_c_without_ns <- c()
  functions_namespace_str_c <- c()
  if(length(functions_namespace_included_id) > 0) {
    itcx <- 0
    for (i in 1:length(functions_namespace_included_id)) {
      itcx <- itcx + 1
      fid_i <- functions_namespace_included_id[i]
      get_ns <- set_getParseData %>% dplyr::filter(parent == fid_i) %>% 
        dplyr::pull(text) %>% paste0(., collapse = "")
      functions_ns <- environmentName(environment(ept(get_ns)))
      set_getParseData <- set_getParseData %>% dplyr::filter(!parent == fid_i)
      if(functions_ns %in%  allowed_namespace_for_sigma_d1()) {
        get_without_ns <- gsub("::", "_", get_ns, fixed = T)
        get_without_ns <- gsub("_:", "_", get_without_ns, fixed = T)
        functions_namespace_included_c <- c(functions_namespace_included_c, 
                                            get_ns)
        functions_namespace_included_c_without_ns <- 
          c(functions_namespace_included_c_without_ns, get_without_ns)
        functions_namespace_str <- paste(deparse(ept(get_ns)), collapse = "\n")
        functions_namespace_str <- paste0(get_without_ns, "<-", 
                                          functions_namespace_str)
        functions_namespace_attr_c <- c(functions_namespace_attr_c, functions_ns)
        functions_namespace_str <- trimws(functions_namespace_str)
        if(functions_ns == 'splines2') {
          gsub_it <- ".engine"
          gsub_by <- paste0(functions_ns, ":::", gsub_it)
          functions_namespace_str <- gsub(gsub_it, gsub_by, 
                                          functions_namespace_str, fixed = T)
        }
        functions_namespace_str_c <- c(functions_namespace_str_c, 
                                       functions_namespace_str)
        if(replace_ns) {
          str <- gsub(get_ns, get_without_ns, str, fixed = T)
        }
    } 
    if(!functions_ns %in%  allowed_namespace_for_sigma_d1()) {
      get_without_ns <- gsub("::", "::", get_ns, fixed = T)
      get_without_ns <- gsub(":::", ":::", get_without_ns, fixed = T)
      functions_namespace_included_c <- c(functions_namespace_included_c, 
                                          get_ns)
      functions_namespace_included_c_without_ns <- 
        c(functions_namespace_included_c_without_ns, get_without_ns)
      functions_namespace_str <- paste(deparse(ept(get_ns)), collapse = "\n")
      functions_namespace_str <- paste0(get_without_ns, "<-", 
                                        functions_namespace_str)
      functions_namespace_attr_c <- c(functions_namespace_attr_c, functions_ns)
      functions_namespace_str <- trimws(functions_namespace_str)
      if(functions_ns == 'splines2') {
        gsub_it <- ".engine"
        gsub_by <- paste0(functions_ns, ":::", gsub_it)
        functions_namespace_str <- gsub(gsub_it, gsub_by, 
                                        functions_namespace_str, fixed = T)
      }
      functions_namespace_str_c <- c(functions_namespace_str_c, 
                                     functions_namespace_str)
      if(replace_ns) {
        str <- gsub(get_ns, get_without_ns, str, fixed = T)
      }
    } 
   } 
  } 
  functions_global_included_id <- set_getParseData %>% 
    dplyr::filter(token=="SYMBOL_FUNCTION_CALL") %>% dplyr::pull(parent)
  functions_global_included_c <- c()
  functions_global_attr_c <- c()
  functions_global_included_c_without_ns <- c()
  functions_global_str_c <- c()
  if(length(functions_global_included_id) > 0) {
    itcx <- 0
    for (i in 1:length(functions_global_included_id)) {
      itcx <- itcx + 1
      fid_i <- functions_global_included_id[i]
      get_ns <- set_getParseData %>% dplyr::filter(parent==fid_i) %>% 
        dplyr::pull(text) %>% paste0(., collapse = "")
      functions_ns <- environmentName(environment(ept(get_ns)))
      get_without_ns <- get_ns
      if(!functions_ns %in%  exclude_global_for_sigma_d1()) {
        functions_global_included_c <- c(functions_global_included_c, get_ns)
        functions_global_included_c_without_ns <- 
          c(functions_global_included_c_without_ns, get_without_ns)
        functions_global_str <- paste(deparse(ept(get_ns)), collapse = "\n")
        functions_global_str <- paste0(get_without_ns, "<-", 
                                       functions_global_str)
        functions_global_attr_c <- c(functions_global_attr_c, functions_ns)
        functions_global_str <- trimws(functions_global_str)
        functions_global_str_c <- 
          c(functions_global_str_c, functions_global_str)
      } 
    }
  } 
  ns_name          <- functions_namespace_included_c %>% unique()
  ns_name_used     <- functions_namespace_included_c_without_ns %>% unique()
  ns_attr_used     <- functions_namespace_attr_c %>% unique()
  global_name      <- functions_global_included_c %>% unique()
  global_name_used <- functions_global_included_c %>% unique()
  global_attr_used <- functions_global_attr_c %>% unique()
  ns_code_used     <- functions_namespace_str_c %>% unique()
  global_code_used <- functions_global_str_c %>% unique()
  name <- c(ns_name_used, global_name_used)
  code <- c(ns_code_used, global_code_used)
  attr <- c(ns_attr_used, global_attr_used)
  out <- list(str  = str, 
              name = name,
              code = code,
              attr = attr)
  return(out)
} 



#' Function to switch condition and by arguments for marginaleffect based draws
#' 
#' @details
#' Used in get_predictions and get_comparisons
#' 
#' @param arg  A list
#' @param condition A string
#' @param by A string
#' @param rm_by A logical
#' @param rm_condition A logical
#' @param verbose logical
#' 
#' @return A list
#'
#' @keywords internal
#' @noRd
#' 
condition_by_switch_fun <- function(arg, 
                                    condition, 
                                    by, 
                                    rm_by = TRUE,
                                    rm_condition = FALSE,
                                    verbose = FALSE) {
  if(rm_by & rm_condition) {
    stop2c("both 'rm_by' and 'rm_condition' can not be TRUE")
  }
  condition.org <- arg[[condition]]
  by.org        <- arg[[by]]
  if(!is.null(condition.org)) {
    arg[[condition]] <- condition.org
    if(rm_by) {
      arg[[by]]        <- NULL
      if(verbose) {
        message2c("The 'by' argument set as NULL")
      }
    }
  } else if(is.null(condition.org)) {
    if(is.null(by.org)) {
      stop2c("pleasae provide 'condition' argument, both 'by' and 'condition' are NULL")
    } else if(!is.null(by.org)) {
      if(!is.logical(by.org)) {
        arg[[condition]] <- by.org
        if(rm_by) {
          arg[[by]]        <- NULL
          if(verbose) {
            message2c("The 'by' argument set as NULL")
          }
        }
        if(verbose) {
          message2c("The 'condition' argument was NULL, now set same as by argument")
        }
      } 
    } 
  } 
  if(rm_condition) {
    arg[[condition]]        <- NULL
    if(verbose) {
      message2c("The 'condition' argument set as NULL")
    }
  }
  return(arg)
} 






#' A function that asks the user to continue or exit.
#' 
#' @details
#' A wrapper around \code{'utils::menu()'}
#'
#' @return NULL
#'
#' @keywords internal
#' @noRd
#' 
user_prompt <- function() {
  choice <- utils::menu(
    choices = c("Press 1 to continue", "Press 0 or esc to exit")
    , title = "Do you want to proceed?")
  return(choice == 1)
}



#' A function that check if model object already exists in results
#' 
#' @details
#' A helper
#'
#' @return NULL
#'
#' @keywords internal
#' @noRd
#' 
check_model_file_exists <- function(model_str, 
                                    set_getwd, 
                                    results_folder, 
                                    force_refit = FALSE, 
                                    assign_model = TRUE,
                                    invisible_TF = TRUE,
                                    envir = NULL) {
  if(is.null(envir)) envir <- parent.frame()
  model_str_rds <- paste0(model_str, ".", 'RDS')
  if(grepl(results_folder, set_getwd)) {
    set_file.path <- file.path(set_getwd, model_str_rds)
  } else {
    set_file.path <- file.path(set_getwd, results_folder, model_str_rds)
  }
  if(force_refit) {
    if(invisible_TF) return(invisible(FALSE)) else return(FALSE)
  }
  if(!force_refit) {
    file_exists_TF <- file.exists(set_file.path)
    if(file_exists_TF) {
      if(assign_model) assign(model_str, readRDS(set_file.path), envir = envir)
      if(invisible_TF)  return(invisible(TRUE)) else return(TRUE)
    } else if(!file_exists_TF) {
      if(invisible_TF)  return(invisible(FALSE)) else return(FALSE)
    }
  } 
} 




#' Build call arguments
#' @noRd
build_args_call <- function(add_args, fun, verbose = FALSE) {
  defaults_it <- base::as.list(base::formals(fun))
  for (i in names(defaults_it)) {
    if(is.null(add_args[[i]])) add_args[[i]] <- defaults_it[[i]]
  }
  defaults <- base::as.list(base::formals(fun))
  build_args <- utils::modifyList(defaults, add_args)
  return(build_args)
}



