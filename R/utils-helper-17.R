

#' @title Stan Code for Bayesian models
#'
#' @description \code{stancode_custom} is a generic function that can be used to
#'   generate custom Stan code sections. Its original use is within the
#'   \pkg{brms} package, but can be used in \pkg{bsitar} when generating 
#'   custom code for \pkg{data} and \pkg{genqunt} blocks.
#'
#' @param object An object whose class will determine which method to apply.
#'   Usually, it will be some kind of symbolic description of the model
#'   form which Stan code should be generated.
#' @param formula Synonym of \code{object} for use in \code{make_stancode_custom}.
#' @param ... Further arguments passed to the specific method.
#'
#' @return Usually, a character string containing the generated Stan code.
#'   For pretty printing, we recommend the returned object to be of class
#'   \code{c("character", "brmsmodel")}.
#'
#' @keywords internal
#' @noRd
#' 
stancode_custom <- function(object, ...) {
  UseMethod("stancode_custom")
}

#' @keywords internal
#' @noRd
#' 
make_stancode_custom <- function(formula, ...) {
  # became an alias of 'stancode_custom' in 2.20.14
  stancode_custom(formula, ...)
}

#' Stan Code for \pkg{brms} Models
#'
#' Generate Stan code for \pkg{brms} models
#'
#' @inheritParams brms::brm
#' @param object An object of class \code{brmsformula}.
#' @param ... Other arguments for internal usage only.
#'
#' @return A character string containing the fully commented \pkg{Stan} code
#'   to fit a \pkg{brms} model. 
#'
#' @examples
#' stancode_custom(rating ~ treat + period + carry + (1|subject),
#'          data = inhaler, family = "cumulative")
#'
#' stancode_custom(count ~ zAge + zBase * Trt + (1|patient),
#'          data = epilepsy, family = "poisson")
#'
#' @keywords internal
#' @noRd
stancode_custom.default <- function(object, data, family = gaussian(),
                             prior = NULL, autocor = NULL, data2 = NULL,
                             cov_ranef = NULL, sparse = NULL,
                             sample_prior = "no", stanvars = NULL,
                             stan_funs = NULL, knots = NULL,
                             drop_unused_levels = TRUE,
                             threads = getOption("brms.threads", NULL),
                             normalize = getOption("brms.normalize", TRUE),
                             save_model = NULL, ...) {
  
  # get funs
  getFromNamespace_brsm_funs <-
  c('.validate_prior', 'as_one_character', 'as_one_logical', 
    'backend_choices', 'brmsframe', 'collapse_lists', 
    'collapse_stanvars', 'collapse_stanvars_pll_args', 
    'escape_all', 'eval2', 'expand_include_statements', 
    'get_data2_autocor', 'get_data2_cov_ranef', 'get_matches', 
    'get_sample_prior', 'glue', 'parse_model', 'require_old_stan_syntax', 
    'stan_Xme', 'stan_clean_pll_args', 'stan_predictor', 'stan_re', 
    'stan_rngprior', 'stan_unchecked_prior', 'str_if', 'usc', 
    'use_threading', 'validate_data', 'validate_data2', 'validate_formula', 
    'validate_stanvars', 'validate_threads', 'wsp_per_line'
    )

  .validate_prior <- as_one_character <- as_one_logical <- backend_choices <- 
    brmsframe <- collapse_lists <- collapse_stanvars <- 
    collapse_stanvars_pll_args <- escape_all <- eval2 <- 
    expand_include_statements <- get_data2_autocor <- 
    get_data2_cov_ranef <- get_matches <- get_sample_prior <- 
    glue <- parse_model <- require_old_stan_syntax <- stan_Xme <- 
    stan_clean_pll_args <- stan_predictor <- stan_re <- stan_rngprior <- 
    stan_unchecked_prior <- str_if <- usc <- use_threading <- 
    validate_data <- validate_data2 <- validate_formula <- 
    validate_stanvars <- validate_threads <- wsp_per_line <- NULL;
  
  for (i in getFromNamespace_brsm_funs) {
    assign(i, utils::getFromNamespace(i, "brms"))
  }
  
  `c<-` <- NULL;
  `c<-` <- utils::getFromNamespace("c<-", "brms")
  `str_add<-` <- NULL;
  `str_add<-` <- utils::getFromNamespace("str_add<-", "brms")
  
  object <- validate_formula(
    object, data = data, family = family,
    autocor = autocor, sparse = sparse,
    cov_ranef = cov_ranef
  )
  bterms <- brmsterms(object)
  data2 <- validate_data2(
    data2, bterms = bterms,
    get_data2_autocor(object),
    get_data2_cov_ranef(object)
  )
  data <- validate_data(
    data, bterms = bterms,
    data2 = data2, knots = knots,
    drop_unused_levels = drop_unused_levels
  )
  bframe <- brmsframe(bterms, data)
  prior <- .validate_prior(
    prior, bframe = bframe,
    sample_prior = sample_prior
  )
  stanvars <- validate_stanvars(stanvars, stan_funs = stan_funs)
  threads <- validate_threads(threads)

  
 .stancode_custom(
   bframe, prior = prior, stanvars = stanvars, threads = threads,
   normalize = normalize, save_model = save_model, ...
 )
}

# internal work function of 'stancode_custom.default'
# @param parse parse the Stan model for automatic syntax checking
# @param backend name of the backend used for parsing
# @param silent silence parsing messages
.stancode_custom <- function(bterms, prior, stanvars, threads = threading(),
                      normalize = getOption("brms.normalize", TRUE),
                      parse = getOption("brms.parse_stancode_custom", FALSE),
                      backend = getOption("brms.backend", "rstan"),
                      silent = TRUE, save_model = NULL, ...) {
  
  dotslist <- list(...)
  
  if(!is.null(dotslist[['return_gq']])) {
    return_gq <- dotslist[['return_gq']]
  } else {
    return_gq <- FALSE
  }
  
  if(!is.null(dotslist[['return_data']])) {
    return_data <- dotslist[['return_data']]
  } else {
    return_data <- FALSE
  }
  
  # get funs
  getFromNamespace_brsm_funs <-
    c('.validate_prior', 'as_one_character', 'as_one_logical', 
      'backend_choices', 'brmsframe', 'collapse_lists', 
      'collapse_stanvars', 'collapse_stanvars_pll_args', 
      'escape_all', 'eval2', 'expand_include_statements', 
      'get_data2_autocor', 'get_data2_cov_ranef', 'get_matches', 
      'get_sample_prior', 'glue', 'parse_model', 'require_old_stan_syntax', 
      'stan_Xme', 'stan_clean_pll_args', 'stan_predictor', 'stan_re', 
      'stan_rngprior', 'stan_unchecked_prior', 'str_if', 'usc', 
      'use_threading', 'validate_data', 'validate_data2', 'validate_formula', 
      'validate_stanvars', 'validate_threads', 'wsp_per_line'
    )
  
  .validate_prior <- as_one_character <- as_one_logical <- backend_choices <- 
    brmsframe <- collapse_lists <- collapse_stanvars <- 
    collapse_stanvars_pll_args <- escape_all <- eval2 <- 
    expand_include_statements <- get_data2_autocor <- 
    get_data2_cov_ranef <- get_matches <- get_sample_prior <- 
    glue <- parse_model <- require_old_stan_syntax <- stan_Xme <- 
    stan_clean_pll_args <- stan_predictor <- stan_re <- stan_rngprior <- 
    stan_unchecked_prior <- str_if <- usc <- use_threading <- 
    validate_data <- validate_data2 <- validate_formula <- 
    validate_stanvars <- validate_threads <- wsp_per_line <- NULL;
  
  for (i in getFromNamespace_brsm_funs) {
    assign(i, utils::getFromNamespace(i, "brms"))
  }
  
  `c<-` <- NULL;
  `c<-` <- utils::getFromNamespace("c<-", "brms")
  `str_add<-` <- NULL;
  `str_add<-` <- utils::getFromNamespace("str_add<-", "brms")
  
  
  
  

  normalize <- as_one_logical(normalize)
  parse <- as_one_logical(parse)
  backend <- match.arg(backend, backend_choices())
  silent <- as_one_logical(silent)
  scode_predictor <- stan_predictor(
    bterms, prior = prior, normalize = normalize,
    stanvars = stanvars, threads = threads
  )
  scode_re <- stan_re(
    bterms, prior = prior, threads = threads, normalize = normalize
  )
  scode_Xme <- stan_Xme(
    bterms, prior = prior, threads = threads, normalize = normalize
  )

  # extend Stan's likelihood part
  if (use_threading(threads)) {
    # threading is activated
    for (i in seq_along(scode_predictor)) {
      resp <- usc(names(scode_predictor)[i])
      pll_args <- stan_clean_pll_args(
        scode_predictor[[i]][["pll_args"]],
        scode_re[["pll_args"]],
        scode_Xme[["pll_args"]],
        collapse_stanvars_pll_args(stanvars)
      )
      partial_log_lik <- paste0(
        scode_predictor[[i]][["pll_def"]],
        scode_predictor[[i]][["model_def"]],
        collapse_stanvars(stanvars, "likelihood", "start"),
        scode_predictor[[i]][["model_comp_basic"]],
        scode_predictor[[i]][["model_comp_eta_basic"]],
        scode_predictor[[i]][["model_comp_eta_loop"]],
        scode_predictor[[i]][["model_comp_dpar_link"]],
        scode_predictor[[i]][["model_comp_dpar_trans"]],
        scode_predictor[[i]][["model_comp_mix"]],
        scode_predictor[[i]][["model_comp_arma"]],
        scode_predictor[[i]][["model_comp_catjoin"]],
        scode_predictor[[i]][["model_comp_mvjoin"]],
        scode_predictor[[i]][["model_log_lik"]],
        collapse_stanvars(stanvars, "likelihood", "end")
      )
      partial_log_lik <- gsub(" target \\+=", " ptarget +=", partial_log_lik)
      partial_log_lik <- paste0(
        "// compute partial sums of the log-likelihood\n",
        "real partial_log_lik", resp, "_lpmf(array[] int seq", resp,
        ", int start, int end", pll_args$typed, ") {\n",
        "  real ptarget = 0;\n",
        "  int N = end - start + 1;\n",
        partial_log_lik,
        "  return ptarget;\n",
        "}\n"
      )
      partial_log_lik <- wsp_per_line(partial_log_lik, 2)
      scode_predictor[[i]][["partial_log_lik"]] <- partial_log_lik
      static <- str_if(threads$static, "_static")
      scode_predictor[[i]][["model_lik"]] <- paste0(
        "  target += reduce_sum", static, "(partial_log_lik", resp, "_lpmf",
        ", seq", resp, ", grainsize", pll_args$plain, ");\n"
      )
      str_add(scode_predictor[[i]][["tdata_def"]]) <- glue(
        "  array[N{resp}] int seq{resp} = sequence(1, N{resp});\n"
      )
    }
    scode_predictor <- collapse_lists(ls = scode_predictor)
    scode_predictor[["model_lik"]] <- paste0(
      scode_predictor[["model_no_pll_def"]],
      scode_predictor[["model_no_pll_comp_basic"]],
      scode_predictor[["model_no_pll_comp_mvjoin"]],
      scode_predictor[["model_lik"]]
    )
    str_add(scode_predictor[["fun"]]) <-
      "  #include 'fun_sequence.stan'\n"
    str_add(scode_predictor[["data"]]) <-
      "  int grainsize;  // grainsize for threading\n"
  } else {
    # threading is not activated
    scode_predictor <- collapse_lists(ls = scode_predictor)
    scode_predictor[["model_lik"]] <- paste0(
      scode_predictor[["model_no_pll_def"]],
      scode_predictor[["model_def"]],
      collapse_stanvars(stanvars, "likelihood", "start"),
      scode_predictor[["model_no_pll_comp_basic"]],
      scode_predictor[["model_comp_basic"]],
      scode_predictor[["model_comp_eta_basic"]],
      scode_predictor[["model_comp_eta_loop"]],
      scode_predictor[["model_comp_dpar_link"]],
      scode_predictor[["model_comp_dpar_trans"]],
      scode_predictor[["model_comp_mix"]],
      scode_predictor[["model_comp_arma"]],
      scode_predictor[["model_comp_catjoin"]],
      scode_predictor[["model_no_pll_comp_mvjoin"]],
      scode_predictor[["model_comp_mvjoin"]],
      scode_predictor[["model_log_lik"]],
      collapse_stanvars(stanvars, "likelihood", "end")
    )
  }
  
  
  ######################################################################
  # return_gq
  ######################################################################
  
  scode_predictor_gq <- stan_predictor(
    bterms, prior = prior, normalize = normalize,
    stanvars = stanvars, threads = threads
  )
  scode_re_gq <- stan_re(
    bterms, prior = prior, threads = threads, normalize = normalize
  )
  scode_Xme_gq <- stan_Xme(
    bterms, prior = prior, threads = threads, normalize = normalize
  )
  
  scode_predictor_gq <- collapse_lists(ls = scode_predictor_gq)
  
  for (i in names(scode_predictor_gq)) {
    scode_predictor_gq[[i]] <- 
      normalize_stancode_custom(scode_predictor_gq[[i]])
  }
  
  scode_predictor_gq[["model_lik"]] <- paste0(
    scode_predictor_gq[["model_no_pll_def"]],
    scode_predictor_gq[["model_def"]],
    scode_predictor_gq[["model_def_be"]],
    scode_predictor_gq[["model_def_re"]],
    collapse_stanvars(stanvars, "likelihood", "start"),
    scode_predictor_gq[["model_no_pll_comp_basic"]],
    scode_predictor_gq[["model_comp_basic"]],
    scode_predictor_gq[["model_comp_eta_basic"]],
    scode_predictor_gq[["model_comp_eta_basic_be"]],
    scode_predictor_gq[["model_comp_eta_basic_re"]],
    scode_predictor_gq[["model_comp_eta_loop"]],
    scode_predictor_gq[["model_comp_eta_loop_re"]],
    scode_predictor_gq[["model_comp_dpar_link"]],
    scode_predictor_gq[["model_comp_dpar_link_be"]],
    scode_predictor_gq[["model_comp_dpar_link_re"]],
    scode_predictor_gq[["model_comp_dpar_trans"]],
    scode_predictor_gq[["model_comp_mix"]],
    scode_predictor_gq[["model_comp_arma"]],
    scode_predictor_gq[["model_comp_catjoin"]],
    scode_predictor_gq[["model_no_pll_comp_mvjoin"]],
    scode_predictor_gq[["model_comp_mvjoin"]],
    # scode_predictor_gq[["model_log_lik"]],
    collapse_stanvars(stanvars, "likelihood", "end")
  )
  
  
  for (i in names(scode_predictor_gq)) {
    scode_predictor_gq[[i]] <- gsub(";", ";\n", 
                                    scode_predictor_gq[[i]], fixed = T)
    scode_predictor_gq[[i]] <- gsub("{", "{\n", 
                                    scode_predictor_gq[[i]], fixed = T)
    scode_predictor_gq[[i]] <- gsub("}", "}\n", 
                                    scode_predictor_gq[[i]], fixed = T)
  }
  
  
  if(return_gq) {
    return(scode_predictor_gq[["model_lik"]])
  }
  
  ######################################################################
  # end return_gq
  ######################################################################
  
  
  scode_predictor[["model_lik"]] <-
    wsp_per_line(scode_predictor[["model_lik"]], 2)
  
 

  # get all priors added to 'lprior'
  scode_tpar_prior <- paste0(
    scode_predictor[["tpar_prior"]],
    scode_re[["tpar_prior"]],
    scode_Xme[["tpar_prior"]]
  )

  # generate functions block
  scode_functions <- paste0(
    "// generated with brms ", utils::packageVersion("brms"), "\n",
    "functions {\n",
      scode_predictor[["fun"]],
      scode_re[["fun"]],
      collapse_stanvars(stanvars, "functions"),
      scode_predictor[["partial_log_lik"]],
    "}\n"
  )

  # generate data block
  scode_data <- paste0(
    "data {\n",
    "  int<lower=1> N;  // total number of observations\n",
    scode_predictor[["data"]],
    scode_re[["data"]],
    scode_Xme[["data"]],
    "  int prior_only;  // should the likelihood be ignored?\n",
    collapse_stanvars(stanvars, "data"),
    "}\n"
  )
  
  ######################################################################
  # return_data
  ######################################################################
  scode_predictor_data <- scode_re_data <- scode_Xme_data <- list()
  scode_predictor_data[["data"]] <- scode_predictor[["data"]]
  scode_re_data[["data"]]        <- scode_re[["data"]]
  scode_Xme_data[["data"]]       <- scode_Xme[["data"]]
  
  # scode_predictor_datax <<- scode_predictor_data
  # scode_re_datax <<- scode_re_data
  # scode_Xme_datax <<- scode_Xme_data
  
  scode_data_data <- paste0(
    # "data {\n",
     "  int<lower=1> N;\n",
    scode_predictor_data[["data"]],
    scode_re_data[["data"]],
    scode_Xme_data[["data"]] # ,
     # "  int prior_only; \n" ,
     # collapse_stanvars(stanvars, "data") ,
    # "}\n"
  )
  
  # for (i in names(scode_predictor_gq)) {
  #   scode_predictor_gq[[i]] <- normalize_stancode_custom(scode_predictor_gq[[i]])
  # }
  
#  scode_data_data <- normalize_stancode_custom(scode_data_data)
  
  # scode_data_data <- gsub(";", ";\n", scode_data_data, fixed = T)
  # scode_data_data <- gsub("{", "{\n", scode_data_data, fixed = T)
  # scode_data_data <- gsub("}", "}\n", scode_data_data, fixed = T)
  
  # for (i in c(";", "{", "}")) {
  #   scode_data_data <- x_gsubit_gsubby(scode_data_data,
  #                                      gsubit = i, gsubby = "\n",
  #                                      pasteit = TRUE, fixed = TRUE)
  # }
  
  
  if(return_data) return(scode_data_data)
  
  ######################################################################
  # end return_data
  ######################################################################

  
  # generate transformed parameters block
  scode_transformed_data <- paste0(
    "transformed data {\n",
       scode_predictor[["tdata_def"]],
       collapse_stanvars(stanvars, "tdata", "start"),
       scode_predictor[["tdata_comp"]],
       collapse_stanvars(stanvars, "tdata", "end"),
    "}\n"
  )

  # generate parameters block
  scode_parameters <- paste0(
    scode_predictor[["par"]],
    scode_re[["par"]],
    scode_Xme[["par"]]
  )
  # prepare additional sampling from priors
  scode_rngprior <- stan_rngprior(
    tpar_prior = scode_tpar_prior,
    par_declars = scode_parameters,
    gen_quantities = scode_predictor[["gen_def"]],
    special_prior = attr(prior, "special"),
    sample_prior = get_sample_prior(prior)
  )
  scode_parameters <- paste0(
    "parameters {\n",
      scode_parameters,
      scode_rngprior[["par"]],
      collapse_stanvars(stanvars, "parameters"),
    "}\n"
  )

  # generate transformed parameters block
  scode_lprior_def <- "  real lprior = 0;  // prior contributions to the log posterior\n"
  scode_transformed_parameters <- paste0(
    "transformed parameters {\n",
      scode_predictor[["tpar_def"]],
      scode_re[["tpar_def"]],
      scode_Xme[["tpar_def"]],
      str_if(normalize, scode_lprior_def),
      collapse_stanvars(stanvars, "tparameters", "start"),
      scode_predictor[["tpar_prior_const"]],
      scode_re[["tpar_prior_const"]],
      scode_Xme[["tpar_prior_const"]],
      scode_predictor[["tpar_comp"]],
      scode_predictor[["tpar_special_prior"]],
      scode_re[["tpar_comp"]],
      scode_Xme[["tpar_comp"]],
      # lprior cannot contain _lupdf functions in transformed parameters
      # as discussed on github.com/stan-dev/stan/issues/3094
      str_if(normalize, scode_tpar_prior),
      collapse_stanvars(stanvars, "tparameters", "end"),
    "}\n"
  )

  # combine likelihood with prior part
  not_const <- str_if(!normalize, " not")
  scode_model <- paste0(
    "model {\n",
      str_if(!normalize, scode_lprior_def),
      collapse_stanvars(stanvars, "model", "start"),
      "  // likelihood", not_const, " including constants\n",
      "  if (!prior_only) {\n",
      scode_predictor[["model_lik"]],
      "  }\n",
      "  // priors", not_const, " including constants\n",
      str_if(!normalize, scode_tpar_prior),
      "  target += lprior;\n",
      scode_predictor[["model_prior"]],
      scode_re[["model_prior"]],
      scode_Xme[["model_prior"]],
      stan_unchecked_prior(prior),
      collapse_stanvars(stanvars, "model", "end"),
    "}\n"
  )
  # generate generated quantities block
  scode_generated_quantities <- paste0(
    "generated quantities {\n",
      scode_predictor[["gen_def"]],
      scode_re[["gen_def"]],
      scode_Xme[["gen_def"]],
      scode_rngprior[["gen_def"]],
      collapse_stanvars(stanvars, "genquant", "start"),
      scode_predictor[["gen_comp"]],
      scode_re[["gen_comp"]],
      scode_rngprior[["gen_comp"]],
      scode_Xme[["gen_comp"]],
      collapse_stanvars(stanvars, "genquant", "end"),
    "}\n"
  )
  # combine all elements into a complete Stan model
  scode <- paste0(
    scode_functions,
    scode_data,
    scode_transformed_data,
    scode_parameters,
    scode_transformed_parameters,
    scode_model,
    scode_generated_quantities
  )

  scode <- expand_include_statements(scode)
  if (parse) {
    scode <- parse_model(scode, backend, silent = silent)
  }
  # if (backend == "cmdstanr") {
  #   if (requireNamespace("cmdstanr", quietly = TRUE) &&
  #       cmdstanr::cmdstan_version() >= "2.29.0") {
  #     tmp_file <- cmdstanr::write_stan_file(scode)
  #     scode <- .canonicalize_stan_model(tmp_file, overwrite_file = FALSE)
  #   }
  # }
  if (is.character(save_model)) {
    cat(scode, file = save_model)
  }
  class(scode) <- c("character", "brmsmodel")
  scode
}



# print.brmsmodel <- function(x, ...) { 
#   cat(x)
#   invisible(x)
# }



#' Extract Stan code from \code{bgmfit} objects
#'
#' Extract Stan code from a fitted \pkg{brms} model.
#'
#' @param object An object of class \code{bgmfit}.
#' @param version Logical; indicates if the first line containing the \pkg{brms}
#'   version number should be included. Defaults to \code{TRUE}.
#' @param regenerate Logical; indicates if the Stan code should be regenerated
#'   with the current \pkg{brms} version. By default, \code{regenerate} will be
#'   \code{FALSE} unless required to be \code{TRUE} by other arguments.
#' @param threads Controls whether the Stan code should be threaded. See
#'   \code{\link{threading}} for details.
#' @param backend Controls the Stan backend. See \code{brm} for details.
#' @param ... Further arguments passed to  \code{stancode_custom.default}.
#'
#' @return Stan code for further processing.
#'
#' @keywords internal
#' @noRd
stancode_custom.bgmfit <- function(object, version = TRUE, regenerate = NULL,
                             threads = NULL, backend = NULL, ...) {
  
  # get funs
   as_one_logical <- validate_threads <- use_threading <- 
   backend_choices <- get_sample_prior <- 
     require_old_stan_syntax <- NULL;
  
  as_one_logical <- utils::getFromNamespace("as_one_logical", "brms")
  validate_threads <- utils::getFromNamespace("validate_threads", "brms")
  use_threading <- utils::getFromNamespace("use_threading", "brms")
  backend_choices <- utils::getFromNamespace("backend_choices", "brms")
  get_sample_prior <- utils::getFromNamespace("get_sample_prior", "brms")
  require_old_stan_syntax <- utils::getFromNamespace("require_old_stan_syntax", "brms")
  
  

  if (is.null(regenerate)) {
    # determine whether regenerating the Stan code is required
    regenerate <- FALSE
    cl <- match.call()
    if ("threads" %in% names(cl)) {
      threads <- validate_threads(threads)
      if (use_threading(threads) && !use_threading(object$threads) ||
          !use_threading(threads) && use_threading(object$threads)) {
        # threading changed; regenerated Stan code
        regenerate <- TRUE
      }
      object$threads <- threads
    }
    if ("backend" %in% names(cl)) {
      backend <- match.arg(backend, backend_choices())
      # older Stan versions do not support array syntax
      if (require_old_stan_syntax(object, backend, "2.29.0")) {
        regenerate <- TRUE
      }
      object$backend <- backend
    }
  }
  regenerate <- as_one_logical(regenerate)
  if (regenerate) {
    object <- restructure(object)
    out <- make_stancode_custom(
      formula = object$formula,
      data = object$data,
      prior = object$prior,
      data2 = object$data2,
      stanvars = object$stanvars,
      sample_prior = get_sample_prior(object$prior),
      threads = object$threads,
      backend = object$backend,
      ...
    )
  } else {
    # extract Stan code unaltered
    out <- object$model
  }
  if (!version) {
    out <- sub("^[^\n]+[[:digit:]]\\.[^\n]+\n", "", out)
  }
  out
}

# expand '#include' and '#includeR' statements
# For '#include' this could also be done automatically by Stan at compilation time
# but would result in Stan code that is not self-contained until compilation
# @param model Stan code that may contain '#include' and '#includeR' statements
# @return Stan code with '#include' and '#includeR' statements expanded
expand_include_statements <- function(model) {
  get_matches <- NULL;
  escape_all <- NULL;
  eval2 <- NULL;
  get_matches <- utils::getFromNamespace("get_matches", "brms")
  escape_all <- utils::getFromNamespace("escape_all", "brms")
  eval2 <- utils::getFromNamespace("eval2", "brms")
  
  # '#include' statements will be replaced by the content of a file
  path <- system.file("chunks", package = "brms")
  includes <- unique(get_matches("#include '[^']+'", model))
  files <- gsub("(#include )|(')", "", includes)
  for (i in seq_along(includes)) {
    code <- readLines(paste0(path, "/", files[i]))
    code <- paste0(code, collapse = "\n")
    pattern <- paste0(" *", escape_all(includes[i]))
    model <- sub(pattern, code, model)
    # remove all duplicated include statements
    model <- gsub(pattern, "", model)
  }
  # '#includeR' statements will be replaced by the call to an R function
  includes <- unique(get_matches("#includeR `[^`]+`", model))
  calls <- gsub("(#includeR )|(`)", "", includes)
  for (i in seq_along(includes)) {
    code <- eval2(calls[i])
    pattern <- paste0(" *", escape_all(includes[i]))
    model <- sub(pattern, code, model)
    # remove all duplicated include statements
    model <- gsub(pattern, "", model)
  }
  model
}

# check if Stan code includes normalization constants
is_normalized <- function(stancode_custom) {
  !grepl("_lup(d|m)f\\(", stancode_custom)
}

# Normalizes Stan code to avoid triggering refit after whitespace and
# comment changes in the generated code.
# In some distant future, StanC3 may provide its own normalizing functions,
# until then this is a set of regex hacks.
# @param x a string containing the Stan code
normalize_stancode_custom <- function(x) {
  as_one_character <- utils::getFromNamespace("as_one_character", "brms")
  x <- as_one_character(x)
  # Remove single-line comments
  x <- gsub("//[^\n\r]*[\n\r]", " ", x)
  x <- gsub("//[^\n\r]*$", " ", x)
  # Remove multi-line comments
  x <- gsub("/\\*([^*]*(\\*[^/])?)*\\*/", " ", x)
  # Standardize whitespace (including newlines)
  x <- gsub("[[:space:]]+"," ", x)
  trimws(x)
}

