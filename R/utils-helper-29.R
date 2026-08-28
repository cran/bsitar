
#' Get Summary Table(s) from one or more Bayesian SITAR models
#'
#' @description
#' 
#' Create customizable tables to summarize Bayesian SITAR models side-by-side.
#' This function supports producing tables in HTML, LaTeX, Word, Markdown,
#' Typst, PDF, PowerPoint, Excel, RTF, JPG, or PNG. The appearance of the tables
#' can be customized extensively by specifying the `output` argument, and by
#' using functions from one of the supported table customization packages:
#' `tinytable`, `kableExtra`, `gt`, `flextable`, `huxtable`, `DT`. For more
#' information, see [modelsummary::modelsummary()]
#' 
#' @param add_component An indicator to set whether \code{nlpar} such as
#'   \code{nlpar}, \code{a}, \code{b}, \code{c} etc must be returned as column.
#'   Default \code{FALSE}
#'   
#' @param sigma_exp An indicator to set whether exponentiation should be applied
#'   to the distributional parameter \code{sigma} estimates. Default \code{NULL}
#'   
#' @param sort_gr_by An indicator to set whether \code{SD} and \code{Cor}
#'   estimates should be sorted by the group by argument. This is only useful
#'   when user requested separate variance covariance for a factor variable by
#'   specifying \code{gr(..., by =)}. Default \code{TRUE}
#'
#' @param gsub_fixed A paired list to set the substitution of character from the
#'   fixed effects names. This is typically to replace co variate by
#'   multiplication sign \code{x}. Default \code{" x "}
#'   
#' @param gsub_random A paired list to set the substitution of character from
#'   the random effects names. This is typically to replace co variate by
#'   multiplication sign \code{x}. Default \code{list(it = c("::", ":"),
#'   by = c(" = ", " x "))}
#'   
#' @param add_diagnostic A logical or a vector of character string indicating to
#'   include diagnostic in the table such as \code{ESS, Rhat, MCSE} (default
#'   \code{TRUE})
#'   
#' @param add_pd A logical indicating whether to include \code{pd} in the table.
#'   Default \code{TRUE}
#'   
#' @param add_section A logical (default \code{TRUE}) to include a column that
#'   includes information on rows (\code{Fixed}, \code{Random}, and
#'   \code{Residual})
#'   
#' @param sort_corr A logical (default \code{TRUE}) that sort \code{SD} and
#'   \code{Cor} estimates within the \code{Random} effects.
#'   
#' @param custom_rename A logical (default \code{TRUE}). Ignored
#' 
#' @param custom_rename_map A paired list to set the substitution of character
#'   from the fixed and random effects names. This is typically to replace
#'   \code{_b}, \code{_sd} etc.
#'   
#' @inheritParams modelsummary::modelsummary
#' @inheritParams parameters::model_parameters
#'
#' @return An object of the supported classes as descibed in description
#' 
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid time-consuming model estimation, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check if the model fit object 'berkeley_exfit' exists
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Population average distance curve
#' summary_table(model)
#' 
#' }
#' 
#' @keywords internal
#' @noRd
#' 
summary_table <- 
  function(
    # modelsummary::modelsummary
    models,
    output = getOption("modelsummary_output", default = "default"),
    fmt = getOption("modelsummary_fmt", default = 3),
    estimate = getOption("modelsummary_estimate", default = "estimate"),
    statistic = getOption("modelsummary_statistic", default = "std.error"),
    vcov = getOption("modelsummary_vcov", default = NULL),
    conf_level = getOption("modelsummary_conf_level", default = 0.95),
    exponentiate = getOption("modelsummary_exponentiate", default = FALSE),
    stars = getOption("modelsummary_stars", default = FALSE),
    shape = getOption("modelsummary_shape", default = term + statistic ~ model),
    coef_map = getOption("modelsummary_coef_map", default = NULL),
    coef_omit = getOption("modelsummary_coef_omit", default = NULL),
    coef_rename = getOption("modelsummary_coef_rename", default = FALSE),
    gof_map = getOption("modelsummary_gof_map", default = NULL),
    gof_omit = getOption("modelsummary_gof_omit", default = NULL),
    gof_function = getOption("modelsummary_gof_function", default = NULL),
    group_map = getOption("modelsummary_group_map", default = NULL),
    add_columns = getOption("modelsummary_add_columns", default = NULL),
    add_rows = getOption("modelsummary_add_rows", default = NULL),
    align = getOption("modelsummary_align", default = NULL),
    notes = getOption("modelsummary_notes", default = NULL),
    title = getOption("modelsummary_title", default = NULL),
    escape = getOption("modelsummary_escape", default = TRUE),
   
    centrality = "median",
    dispersion = FALSE,
    ci = 0.95,
    ci_method = "eti",
    test = "pd",
    rope_range = "default",
    rope_ci = 0.95,
    bf_prior = NULL,
    diagnostic = c("ESS", "Rhat"),
    priors = FALSE,
    effects = "fixed",
    component = "all",
    standardize = NULL,
    group_level = FALSE,
    keep = NULL,
    drop = NULL,
    verbose = TRUE,
    add_component = FALSE,
    sigma_exp = NULL,
    sort_gr_by = TRUE,
    gsub_fixed = " x ",
    gsub_random = list(it = c("::", ":"), 
                       by = c(" = ", " x ")),
    add_diagnostic = TRUE,
    add_pd = TRUE,
    add_section = TRUE,
    sort_corr = TRUE,
    custom_rename = TRUE,
    custom_rename_map = NULL,
    ...) {
    
    insight::check_if_installed('performance',  minimum_version = '0.17.1')
    insight::check_if_installed('modelsummary', minimum_version = '2.6.0')
    insight::check_if_installed('parameters',   minimum_version = '0.29.2')
    
    
    custom_model_parameters.brmsfit <- function (model, 
                                                 centrality = "median",
                                                 dispersion = FALSE, 
                                                 ci = 0.95, 
                                                 ci_method = "eti", 
                                                 test = "pd", 
                                                 rope_range = "default", 
                                                 rope_ci = 0.95, 
                                                 bf_prior = NULL, 
                                                 diagnostic = c("ESS", "Rhat"), 
                                                 priors = FALSE, 
                                                 effects = "fixed", 
                                                 component = "all", 
                                                 exponentiate = FALSE, 
                                                 standardize = NULL, 
                                                 group_level = FALSE, 
                                                 keep = NULL, 
                                                 drop = NULL, 
                                                 verbose = TRUE, 
                                                 ...) {
      extra_arg <- list(...)
      sigma_exp          <- extra_arg[['sigma_exp']]
      sort_gr_by         <- extra_arg[['sort_gr_by']]
      # pull_bsigma        <- extra_arg[['pull_bsigma']]
      add_diagnostic     <- extra_arg[['add_diagnostic']]
      add_pd             <- extra_arg[['add_pd']]
      add_section        <- extra_arg[['add_section']]
      sort_corr          <- extra_arg[['sort_corr']]
      custom_rename      <- extra_arg[['custom_rename']]
      custom_rename_map  <- extra_arg[['custom_rename_map']]
      gsub_random        <- extra_arg[['gsub_random']]
      gsub_fixed         <- extra_arg[['gsub_fixed']]
      envir              <- extra_arg[['envir']]
      
      if(is.null(sigma_exp)) {
        model_link_sigma <- model$family$link_sigma
        if(model_link_sigma == "log") sigma_exp <- TRUE
      }
      
      sort_random_rows <- function(data, 
                                   select_col = 'Parameter',
                                   by_order = NULL,
                                   verbose = FALSE) {
        
        data$row_id <- seq_len(nrow(data))
        data$is_cor <- grepl("^cor_", data[[select_col]])
        data$sex_group <- vapply(
          data[[select_col]],
          function(x) {
            hit <- by_order[by_order %in% 
                              regmatches(x, 
                                         gregexpr(paste(by_order, 
                                                        collapse = "|"), 
                                                  x))[[1]]]
            if (length(hit) == 0) Inf else match(hit[1], by_order)
          },
          numeric(1)
        )
        out <- data[order(data$is_cor, data$sex_group, data$row_id), ]
        out[c("row_id", "is_cor", "sex_group")] <- NULL
        return(out)
      }

      .add_pretty_names2 <- function (params, model) {
        .safe <- NULL;
        getfrom_ <- c('.safe')
        for (i in getfrom_) {
          assign(i, utils::getFromNamespace(i, 'parameters'))
        }
        attr(params, "model_class") <- class(model)
        cp <- .get_cleaned_parameters(params, model)
        clean_params <- cp[cp$Parameter %in% params$Parameter, ]
        named_clean_params <- stats::setNames(
          clean_params$Cleaned_Parameter[match(params$Parameter, 
                                               clean_params$Parameter)], 
          params$Parameter)
        
        if (!is.null(clean_params$Group) && any(nzchar(clean_params$Group, 
                                                       keepNA = TRUE))) {
        }
        attr(params, "cleaned_parameters") <- named_clean_params
        attr(params, "pretty_names") <- named_clean_params
        return(params)
      } 
      
      transform_names <- function(x, 
                                  prefix = c("id", "study"), 
                                  corsub = 'vs') {
        prefix_pattern <- paste(prefix, collapse = "|")
        search_pattern <- paste0(prefix_pattern, "__")
        matches <- grepl(search_pattern, x)
        x[matches] <- sub(
          paste0("(\\w+\\s+)?\\b(", prefix_pattern, ")__(.*)"),
          "\\1\\3 (\\2)",
          x[matches])
        if(!is.null(corsub)) {
          x <- gsub(
            paste0("Cor (.*?)__(.*?) \\((", prefix_pattern, ")\\)"),
            paste0("Cor \\1", corsub, "\\2 (\\3)"),
            x)
        } 
        return(x)
      }
      
      # by_order
      if(!is_emptyx(model$ranef$by)) {
        gr_by       <- model$ranef$by
        gr_bylevels <- model$ranef$bylevels
        sort_gr_bylevels <- unlist(gr_bylevels) %>% unique()
        by_order    <- paste0(gr_by, unlist(gr_bylevels))
        by_order    <- unique(by_order)
        if(is.null(sort_gr_by)) sort_gr_by  <- sort_gr_by
      } else {
        gr_by    <- NULL
        by_order    <- NULL
        sort_gr_bylevels <- NULL
        if(is.null(sort_gr_by)) sort_gr_by  <- FALSE
      }

      fixed_covs <- model$model_info$covs

      replace_terms_base <- function(df, 
                                     gsub_it, 
                                     gsub_by, 
                                     select_col = 'Parameter', 
                                     group = NULL,
                                     corsub = NULL, 
                                     fixed = FALSE, 
                                     call_custom = FALSE,
                                     add_suffix = F) {
        if(!call_custom) {
          for (i in seq_along(gsub_it)) {
            name_it <- df[[select_col]]
            gsub_iti <- gsub_it[i]
            gsub_byi <- gsub_by[i]
            df[[select_col]] <- gsub(gsub_iti,
                                     gsub_byi,
                                     df[[select_col]],
                                     fixed = fixed)
          } 
        } 
        if(call_custom & !add_suffix) {
          for (i in seq_along(gsub_it)) {
            name_it <- df[[select_col]]
            gsub_iti <- gsub_it[i]
            gsub_byi <- gsub_by[i]
            name_by <- replace_string_part(df[[select_col]],
                                                    start = gsub_iti, 
                                                    end = "", 
                                                    replace = gsub_byi, 
                                                    extract = F,
                                                    exclude_start = F,
                                                    exclude_end = F)
            df[[select_col]] <- name_by
          } 
        } 
        if(!is.null(group)) {
          df[[select_col]] <- transform_names(x = df[[select_col]], 
                                              prefix = group,
                                              corsub = corsub)
        }
        return(df)
      }
      
      extract_sex <- function(x, by) {
        out <- rep(NA_character_, length(x))
        for (byi in by) {
          out[grepl(byi, x)] <- byi
          out[grepl(byi, x)] <- byi
        }
        return(out)
      }

      .extract_parameters_bayesian <- .add_pretty_names <- NULL;
      .add_model_parameters_attributes <- .get_cleaned_parameters <- NULL;
      .exponentiate_parameters <- NULL;
      .model_parameters_brms_meta <- .group_level_total <- NULL; 
      
      getfrom_ <- c('.extract_parameters_bayesian', 
                    '.add_pretty_names',
                    ".add_model_parameters_attributes",
                    ".get_cleaned_parameters",
                    ".model_parameters_brms_meta",
                    ".group_level_total",
                    '.exponentiate_parameters')
      for (i in getfrom_) {
        assign(i, utils::getFromNamespace(i, 'parameters'))
      }
      
      modelinfo <- insight::model_info(model, verbose = FALSE)
      if (!insight::is_multivariate(model) && isTRUE(modelinfo$is_meta)) {
        params <- .model_parameters_brms_meta(model, 
                                              centrality = centrality, 
                                              dispersion = dispersion, 
                                              ci = ci, 
                                              ci_method = ci_method, 
                                              test = test, 
                                              rope_range = rope_range, 
                                              rope_ci = rope_ci, 
                                              diagnostic = diagnostic, 
                                              priors = priors, 
                                              exponentiate = exponentiate, 
                                              standardize = standardize, 
                                              keep_parameters = keep, 
                                              drop_parameters = drop, 
                                              ...)
      }
      else if (effects %in% c("total", "random_total")) {
        params <- .group_level_total(model, centrality, dispersion, 
                                     ci, ci_method, test, 
                                     rope_range, rope_ci, ...)
        params$Effects <- "total"
        class(params) <- c("parameters_coef", "see_parameters_coef", 
                           class(params))
      }
      else {
        if (effects == "random" && group_level) {
          effects <- "grouplevel"
        }
        params <- .extract_parameters_bayesian(model, 
                                               centrality = centrality, 
                                               dispersion = dispersion, 
                                               ci = ci, 
                                               ci_method = ci_method, 
                                               test = test, 
                                               rope_range = rope_range, 
                                               rope_ci = rope_ci, 
                                               bf_prior = bf_prior,
                                               diagnostic = diagnostic, 
                                               priors = priors, 
                                               effects = effects, 
                                               component = component, 
                                               standardize = standardize, 
                                               keep_parameters = keep, 
                                               drop_parameters = drop, 
                                               verbose = verbose, 
                                               ...)
        
        Component <- NULL;
        Parameter <- NULL;
        tempxx_ <- NULL;
        params <- .add_pretty_names2(params, model)
        params <- .exponentiate_parameters(params, model, exponentiate)
        params <- .add_model_parameters_attributes(params, model, 
                                                   ci, exponentiate, 
                                                   ci_method = ci_method, 
                                                   group_level = group_level, 
                                                   modelinfo = modelinfo, 
                                                   verbose = verbose, ...)
        attr(params, "parameter_info") <- .get_cleaned_parameters(params, 
                                                                  model)
        attr(params, "object_name") <- 
          insight::safe_deparse_symbol(substitute(model))
        attr(params, "dpars") <- 
          insight::find_auxiliary(model, verbose = FALSE)
        class(params) <- unique(c("parameters_model", 
                                  "see_parameters_model", 
                                  class(params)))
      }

      if(add_section) {
        # Residual
        Residual <- params %>% dplyr::filter(grepl("^Sigma", Component) | 
                                               grepl("^sigma", Component)) %>%
          dplyr::mutate(section = "Residual")
        
        if(sigma_exp) {
          exp_these <- c("Mean", "Median", "CI", "CI_low", "CI_high")
          for (exp_thesei in exp_these) {
            if(!is.null(Residual[[exp_thesei]])) {
              Residual[[exp_thesei]] <- exp(Residual[[exp_thesei]])
            }
          }
        } 
        
        # Non_Residual
        Non_Residual <- params %>% dplyr::filter(!grepl("^Sigma", Component) &
                                                   !grepl("^sigma", Component)) %>%
          dplyr::mutate(section = "Non_Residual")
        
        
        # ########## Fixed ##########
        Fixed <- Non_Residual %>% dplyr::filter(grepl("^b_", Parameter) | 
                                                  grepl("^b_", Parameter)) %>%
          dplyr::mutate(section = "Fixed")
        
        # Random
        Random <- Non_Residual %>% dplyr::filter(grepl("^sd_", Parameter) | 
                                                   grepl("^cor_", Parameter)) %>%
          dplyr::mutate(section = "Random")
        Random <- sort_random_rows(Random, by_order = by_order)
        
        if(sort_gr_by) {
          Random <- Random %>% dplyr::mutate(
            tempxx_ = extract_sex(Parameter, by_order),
            tempxx_ = factor(tempxx_, levels = by_order)) %>%
            dplyr::arrange(tempxx_) %>% dplyr::select(-tempxx_)
        }
        
        gsub_random_idx <- Random$Parameter != "gof"
        if(!is.null(gsub_random)) {
          gsub_random_it <- gsub_random[[1]]
          gsub_random_by <- gsub_random[[2]]
          for (ixxx in 1:length(gsub_random_it)) {
            Random$Parameter <- ifelse(gsub_random_idx, 
                                       gsub(gsub_random_it[ixxx],
                                            gsub_random_by[ixxx],
                                            Random$Parameter),
                                       Random$Parameter)
          }
        }
        
        # all_section
        all_section <- list()
        if(!  is_emptyx(Fixed)) all_section[['Fixed']] <- Fixed
        if(!  is_emptyx(Random)) all_section[['Random']] <- Random
        if(!  is_emptyx(Residual)) all_section[['Residual']] <-Residual
        params <- do.call(rbind, all_section)
      }
      
      # b_ sd_ cor_
      custom_rename_map[[1]] <- c("^b_", "^sd_", "^cor_")
      custom_rename_map[[2]] <- c("", "SD ", "Cor ")
      params <- replace_terms_base(df = params,
                                   gsub_it = custom_rename_map[[1]],
                                   gsub_by = custom_rename_map[[2]],
                                   select_col = 'Parameter',
                                   group = NULL,
                                   corsub = NULL, 
                                   fixed = FALSE,
                                   call_custom = FALSE,
                                   add_suffix = F)
      
      #a_ b_ c_ d_
      custom_rename_map[[1]] <- c("a_",    "b_",      "c_",         " d_",
                                  "_Intercept", "sigmsize", "sigmSize")
      custom_rename_map[[2]] <- c("Size_", "Timing_", "Intensity_", " Rate_",
                                  "",           'Sigma', "Sigma")
      re_group <- model$ranef$group %>% unique()
      params <- replace_terms_base(df = params,
                                   gsub_it = custom_rename_map[[1]],
                                   gsub_by = custom_rename_map[[2]],
                                   select_col = 'Parameter',
                                   group = re_group,
                                   corsub = ' - ', 
                                   fixed = FALSE,
                                   call_custom = T,
                                   add_suffix = F)
      
      # fixed_covs
      custom_rename_map[[1]] <- paste0("_", fixed_covs)
      custom_rename_map[[2]] <- c(" ")
      if(!is.null(gsub_fixed)) custom_rename_map[[2]] <- gsub_fixed
      
      params <- replace_terms_base(df = params,
                                   gsub_it = custom_rename_map[[1]],
                                   gsub_by = custom_rename_map[[2]],
                                   select_col = 'Parameter',
                                   group = NULL,
                                   corsub = NULL,
                                   fixed = FALSE,
                                   call_custom = FALSE,
                                   add_suffix = F)
      
      
      # residual:
      custom_rename_map[[1]] <- paste0("", ":")
      custom_rename_map[[2]] <- c(" x ")
      
      params <- replace_terms_base(df = params,
                                   gsub_it = custom_rename_map[[1]],
                                   gsub_by = custom_rename_map[[2]],
                                   select_col = 'Parameter',
                                   group = NULL,
                                   corsub = NULL,
                                   fixed = FALSE,
                                   call_custom = FALSE,
                                   add_suffix = F)
      custom_rename_map[[1]] <- paste0("", gr_by)
      custom_rename_map[[2]] <- c("")
      params <- replace_terms_base(df = params,
                                   gsub_it = custom_rename_map[[1]],
                                   gsub_by = custom_rename_map[[2]],
                                   select_col = 'Parameter',
                                   group = NULL,
                                   corsub = NULL,
                                   fixed = FALSE,
                                   call_custom = FALSE,
                                   add_suffix = F)
      
      # add_diagnostic
      if(add_section) custom_add <- c('section') else custom_add <- c()
      if(add_diagnostic) {
        if(!is.null(diagnostic)) {
          if(length(diagnostic) == 1) {
            if(diagnostic == 'all') {
              diagnostic <- c("ESS", "Rhat", "MCSE")
            }
          }
        } 
        custom_add <- c(custom_add, diagnostic)
      } 
      
      if(add_pd) {
        custom_add <- c(custom_add, 'pd')
      }
      attr(params, "custom_add") <- custom_add
      return(params)
    } 
    
    custom_get_gof <- function (model, 
                                gof_function = NULL, 
                                vcov_type = NULL, ...) {
      
      get_gof_broom <- NULL;
      get_gof_parameters <- NULL;
      glance_custom_internal <- NULL;
      format_gof <- NULL;
      bind_cols <- NULL;
      modelsummary_rbind <- NULL;
      settings_init <- NULL;
      sanitize_output <- NULL;

      getfrom_ <- c('get_gof_broom', 
                    "get_gof_parameters",
                    "glance_custom_internal",
                    "bind_cols",
                    "modelsummary_rbind",
                    "settings_init",
                    "sanitize_output",
                    "format_gof")
      for (i in getfrom_) {
        assign(i, utils::getFromNamespace(i, 'modelsummary'))
      }
      
      dots <- list(...)
      if (isTRUE(is.na(dots$gof_map))) 
        return(NULL)
      checkmate::assert_function(gof_function, null.ok = TRUE)
      get_priority <- getOption("modelsummary_get", default = "easystats")
      checkmate::assert_choice(get_priority, choices = c("broom", 
                                                         "easystats", 
                                                         "parameters", 
                                                         "performance", 
                                                         "all"))
      if (get_priority %in% c("all", "broom")) {
        funs <- list(broom = get_gof_broom, parameters = get_gof_parameters)
      }
      else {
        funs <- list(parameters = get_gof_parameters, broom = get_gof_broom)
      }
      warning_msg <- NULL
      gof <- NULL
      for (f in names(funs)) {
        if (get_priority == "all") {
          tmp <- funs[[f]](model, ...)
          if (!is.null(tmp)) {
            attr(tmp, "backend") <- f
          }
          if (inherits(tmp, "data.frame") && inherits(gof, 
                                                      "data.frame")) {
            idx <- !tolower(colnames(tmp)) %in% tolower(colnames(gof))
            tmp <- tmp[, idx, drop = FALSE]
            if (ncol(tmp) > 0) {
              gof <- bind_cols(gof, tmp)
            }
          }
          else if (inherits(tmp, "data.frame")) {
            gof <- tmp
          }
          else {
            warning_msg <- c(warning_msg, tmp)
          }
        }
        else {
          if (!inherits(gof, "data.frame")) {
            gof <- funs[[f]](model, ...)
            if (!is.null(gof)) {
              attr(gof, "backend") <- f
            }
          }
        }
      }
      if (is.character(vcov_type) && !vcov_type %in% c("matrix", 
                                                       "vector", "function")) {
        gof$vcov.type <- vcov_type
      }
      
      assign('glance_custom', utils::getFromNamespace('get_gof_parameters', 
                                                      'modelsummary'))
      
      set_metrics <- NULL
      if (length(gof_map) == 1) {
        if(gof_map == "all") {
          set_metrics <- "all"
        } else if(gof_map == "common") {
          set_metrics <- c("LOOIC", "WAIC", "R2", "RMSE")
        }
      } else if (length(gof_map) > 1) {
        gof_map_c <- c()
        for (gof_mapi in gof_map) {
          if(gof_mapi == "elpd") {
            set_metrics <- "ELPD"
          } else if(gof_mapi == "loo") {
            set_metrics <- "LOOIC"
          } else if(gof_mapi == "looic") {
            set_metrics <- "LOOIC"
          }  else if(gof_mapi == "r2") {
            set_metrics <- "R2"
          }  else if(gof_mapi == "r2_adj") {
            set_metrics <- "R2_adj"
          } else if(gof_mapi == "R2_adj") {
            set_metrics <- "R2_adj"
          }  else if(gof_mapi == "rmse") {
            set_metrics <- "RMSE"
          }  else if(gof_mapi == "waic") {
            set_metrics <- "WAIC"
          } else if(gof_mapi == "score") {
            set_metrics <- "SCORE"
          } else if(gof_mapi == "sigma") {
            set_metrics <- "SIGMA"
          } else if(gof_mapi == "logloss") {
            set_metrics <- "LOGLOSS"
          }
          gof_map_c <- c(gof_map_c, set_metrics)
        } 
        set_metrics <- gof_map_c
      } 

      if(!is.null(set_metrics)) {
        set_metrics   <- unique(set_metrics)
        gof_custom_df <- glance_custom(model, metrics = set_metrics)
        return(gof_custom_df)
      } else {
        return(NULL)
      }

      if (!is.null(gof_custom_df) && is.data.frame(gof)) {
        for (n in colnames(gof_custom_df)) {
          if (is.null(vcov_type) || n != "vcov.type") {
            gof[[n]] <- gof_custom_df[[n]]
          }
        }
      }
      if (is.function(gof_function)) {
        if (!"model" %in% names(formals(gof_function))) {
          msg <- "`gof_function` must accept an argument named 'model'."
          insight::format_error(msg)
        }
        tmp <- try(gof_function(model = model))
        if (!isTRUE(checkmate::check_data_frame(tmp, nrows = 1, 
                                                col.names = "unique"))) {
          msg <- "`gof_function` must be a function which accepts a 
          model and returns a 1-row data frame with unique column names."
          insight::format_error(msg)
        }
        else {
          for (n in names(tmp)) {
            gof[[n]] <- tmp[[n]]
          }
        }
      }
      for (i in rev(seq_along(gof))) {
        if (isTRUE(is.na(gof[[i]]))) {
          gof[[i]] <- NULL
        }
      }
      for (col in colnames(gof)) {
        if (inherits(gof[[col]], "logLik")) {
          gof[[col]] <- as.numeric(gof[[col]])
        }
      }
      if (inherits(gof, "data.frame")) {
        return(gof)
      }
      warning(sprintf("`modelsummary could not extract 
            goodness-of-fit statistics from a model\nof class \"%s\". 
            The package tried a sequence of 2 helper 
            functions:\n\
            nperformance::model_performance(model)\
            nbroom::glance(model)\n\nOne of these 
            functions must return a one-row `data.frame`. 
            The `modelsummary` website explains how to summarize 
            unsupported models or add support for new models 
            yourself:\n\nhttps://modelsummary.com/vignettes/modelsummary.html", 
                      class(model)[1]), call. = FALSE)
    }
    
    custom_modelsummary <- function (
    models, 
    output = getOption("modelsummary_output", default = "default"), 
    fmt = getOption("modelsummary_fmt", default = 3), 
    estimate = getOption("modelsummary_estimate", default = "estimate"), 
    statistic = getOption("modelsummary_statistic", default = "std.error"), 
    vcov = getOption("modelsummary_vcov", default = NULL), 
    conf_level = getOption("modelsummary_conf_level", default = 0.95),
    exponentiate = getOption("modelsummary_exponentiate", default = FALSE), 
    stars = getOption("modelsummary_stars", default = FALSE), 
    shape = getOption("modelsummary_shape", default = term + statistic ~ model), 
    coef_map = getOption("modelsummary_coef_map", default = NULL), 
    coef_omit = getOption("modelsummary_coef_omit", default = NULL), 
    coef_rename = getOption("modelsummary_coef_rename", default = FALSE), 
    gof_map = getOption("modelsummary_gof_map", default = NULL), 
    gof_omit = getOption("modelsummary_gof_omit", default = NULL), 
    gof_function = getOption("modelsummary_gof_function", default = NULL), 
    group_map = getOption("modelsummary_group_map", default = NULL), 
    add_columns = getOption("modelsummary_add_columns", default = NULL),
    add_rows = getOption("modelsummary_add_rows", default = NULL),
    align = getOption("modelsummary_align", default = NULL), 
    notes = getOption("modelsummary_notes", default = NULL), 
    title = getOption("modelsummary_title", default = NULL), 
    escape = getOption("modelsummary_escape", default = TRUE), 
    ...) {
      
      glance_custom <- NULL;
      assign('glance_custom', utils::getFromNamespace('get_gof_parameters', 
                                                      'modelsummary'))
      
      settings_equal <- settings_init <- sanitize_output <- NULL;
      get_span_cbind <- sanitize_estimate <- NULL;
      sanitize_exponentiate <- sanitize_exponentiate <- NULL;
      sanitize_statistic <- sanitize_gof_map <- sanitize_fmt <- NULL;
      sanity_group_map <- map_estimates <- NULL;
      sanitize_conf_level <- sanitize_shape <- NULL;
      sanity_coef <- sanity_stars <- sanity_align <- sanitize_escape <- NULL;
      sanity_ellipsis <- sanitize_models <- sanitize_vcov <- NULL;
      pad <- get_list_of_modelsummary_lists <- format_estimates <- NULL;
      shape_estimates <- map_gof <- bind_est_gof <- settings_get <- NULL;
      factory <- set_span_cbind <- redundant_labels <- format_gof <- NULL;
      settings_rm <- NULL;
      
      modelsummary_rbind <- NULL;
      replace_dict <- NULL;
      make_stars_note <- NULL;
      
      getfrom_ <- c('settings_equal', 
                    "settings_init",
                    "sanitize_output",
                    "get_span_cbind",
                    "sanitize_estimate",
                    "sanitize_exponentiate",
                    "sanitize_statistic", 
                    "sanitize_gof_map", 
                    "sanitize_fmt", 
                    "sanity_group_map", 
                    "sanitize_conf_level", 
                    "sanity_coef", 
                    "sanity_stars", 
                    "sanity_align",
                    "sanitize_escape",
                    "sanity_ellipsis",
                    "sanitize_models",
                    "sanitize_vcov",
                    "sanitize_shape",
                    "pad",
                    "get_list_of_modelsummary_lists",
                    "format_estimates",
                    "map_estimates",
                    "shape_estimates",
                    "map_gof",
                    "bind_est_gof",
                    "settings_get",
                    "factory",
                    "set_span_cbind",
                    "redundant_labels",
                    "format_gof",
                    "modelsummary_rbind",
                    "replace_dict",
                    "make_stars_note",
                    "settings_rm")
      for (i in getfrom_) {
        assign(i, utils::getFromNamespace(i, 'modelsummary'))
      }
      checkmate::assert(checkmate::check_formula(shape), 
                        checkmate::check_choice(shape, 
                                                c("cbind", 
                                                  "rbind",
                                                  "rcollapse")), 
                        checkmate::check_null(shape))
      
      if (isTRUE(checkmate::check_choice(shape, c("rbind", "rcollapse")))) {
        out <- modelsummary_rbind(models,
                                  output = output, 
                                  fmt = fmt, 
                                  estimate = estimate, 
                                  statistic = statistic, 
                                  vcov = vcov, 
                                  conf_level = conf_level, 
                                  exponentiate = exponentiate, 
                                  stars = stars, 
                                  coef_map = coef_map, 
                                  coef_omit = coef_omit, 
                                  coef_rename = coef_rename, 
                                  gof_map = gof_map, 
                                  gof_omit = gof_omit, 
                                  gof_function = gof_function, 
                                  add_columns = add_columns, 
                                  add_rows = add_rows, 
                                  align = align, 
                                  shape = shape, 
                                  group_map = NULL, 
                                  notes = notes, 
                                  title = title, 
                                  escape = escape, 
                                  ...)
        return(out)
      }
      
      extra_arg <- list(...)
      add_component          <- extra_arg[['add_component']]
      dots <- list(...)
      if (!settings_equal("function_called", "modelsummary_rbind")) {
        settings_init(settings = list(function_called = "modelsummary"))
      }
      checkmate::assert_string(gof_omit, null.ok = TRUE)
      tmp <- sanitize_output(output)
      output_format <- tmp$output_format
      output_factory <- tmp$output_factory
      output_file <- tmp$output_file
      tmp <- get_span_cbind(models, shape)
      shape <- tmp$shape
      models <- tmp$models
      span_cbind <- tmp$span_cbind
      sanitize_escape(escape)
      sanity_ellipsis(vcov, ...)
      models <- sanitize_models(models, ...)
      vcov <- sanitize_vcov(vcov, models, ...)
      number_of_models <- max(length(models), length(vcov))
      estimate <- sanitize_estimate(estimate, number_of_models)
      exponentiate <- sanitize_exponentiate(exponentiate, number_of_models)
      shape <- sanitize_shape(shape)
      statistic <- sanitize_statistic(statistic, shape, conf_level)
      gof_map <- sanitize_gof_map(gof_map)
      fmt <- sanitize_fmt(fmt, calling_function = "modelsummary")
      sanity_group_map(group_map)
      conf_level <- sanitize_conf_level(conf_level, estimate, statistic)
      sanity_coef(coef_map, coef_rename, coef_omit)
      sanity_stars(stars)
      sanity_align(align, estimate = estimate, 
                   statistic = statistic, stars = stars)
      checkmate::assert_function(gof_function, null.ok = TRUE)
      if (!any(grepl("conf", c(estimate, statistic)))) {
        conf_level <- NULL
      }
      
      modelsummary_model_labels <- getOption("modelsummary_model_labels", 
                                             default = "(arabic)")
      if (is.null(names(models))) {
        checkmate::assert_choice(modelsummary_model_labels, 
                                 choices = c("model", 
                                             "arabic", 
                                             "letters",
                                             "roman", 
                                             "(arabic)", "(letters)", 
                                             "(roman)"))
        
        if (modelsummary_model_labels == "model") {
          model_names <- paste("Model", 1:number_of_models)
        }
        else if (grepl("arabic", modelsummary_model_labels)) {
          model_names <- as.character(1:number_of_models)
        }
        else if (grepl("letters", modelsummary_model_labels)) {
          model_names <- LETTERS[1:number_of_models]
        }
        else if (grepl("roman", modelsummary_model_labels)) {
          model_names <- as.character(utils::as.roman(1:number_of_models))
        }
        if (grepl("\\(", modelsummary_model_labels)) {
          model_names <- sprintf("(%s)", model_names)
        }
      }
      else {
        model_names <- names(models)
      }
      model_names <- pad(model_names, output_format = output_format)
      if (!settings_equal("function_called", "modelsummary_rbind") && 
          all(grepl("^\\(\\d+\\)$", model_names)) && identical(output_format, 
                                                               "kableExtra")) {
        model_names <- paste0("&nbsp;", model_names)
      }
      
      msl <- get_list_of_modelsummary_lists(models = models, 
                                            conf_level = conf_level, 
                                            vcov = vcov, 
                                            gof_map = gof_map, # 'WAIC',
                                            gof_function = glance_custom,
                                              # modelsummary:::get_gof_parameters, 
                                            shape = shape,
                                            coef_rename = coef_rename,
                                            output_format = output_format, 
                                            ...)
      names(msl) <- model_names
      attr(msl, "backend") <- list()
      for (mod_name in names(msl)) {
        attr(msl, "backend")[[mod_name]] <- attr(msl[[mod_name]], 
                                                 "backend")
      }
      if (identical(output_format, "modelsummary_list")) {
        if (length(msl) == 1) {
          return(msl[[1]])
        }
        else {
          return(msl)
        }
      }

      est <- list()
      for (i in seq_along(msl)) {
        tmp <- format_estimates(est = msl[[i]]$tidy, fmt = fmt, 
                                estimate = estimate[[i]], 
                                estimate_label = names(estimate)[1], 
                                statistic = statistic, 
                                vcov = vcov[[i]], 
                                conf_level = conf_level, 
                                stars = stars, 
                                shape = shape, 
                                group_name = shape$group_name, 
                                exponentiate = exponentiate[[i]], ...)
        
        if(add_component) {
          tmp$component <- msl[[i]]$tidy$component
          tmp <- tmp %>% dplyr::relocate(component, .after = 'term')
        }

        tmp_custom_est <- list()
        for_custom_est <- msl[[i]]$tidy
        if(!is.null(attr(for_custom_est, "custom_add"))) {
          for (for_custom_esti in attr(for_custom_est, "custom_add")) {
            for_custom_esti_lc <- tolower(for_custom_esti)
            if(for_custom_esti_lc == "ess") {
              tmp_custom_est[[for_custom_esti]] <- 
                formatC(msl[[i]]$tidy[[for_custom_esti_lc]],
                        format = "f", digits = 0)
            } else if(for_custom_esti_lc == "rhat") {
              tmp_custom_est[[for_custom_esti]] <- 
                formatC(msl[[i]]$tidy[[for_custom_esti_lc]],
                        format = "f", digits = 2)
            } else {
              tmp_custom_est[[for_custom_esti]] <- 
                fmt(msl[[i]]$tidy[[for_custom_esti_lc]])
            }
          } 
        } 
        
        tmp <- c(tmp, tmp_custom_est)
        tmp <- as.data.frame(tmp)
        tmp <- tmp %>% dplyr::relocate(dplyr::any_of('section'))
        
        tmp <- map_estimates(tmp, 
                             coef_rename = coef_rename, 
                             coef_map = coef_map, 
                             coef_omit = coef_omit, 
                             group_map = group_map, 
                             shape = shape)
        colnames(tmp)[match("modelsummary_value", 
                            colnames(tmp))] <- model_names[i]
        est[[model_names[i]]] <- tmp
      }
      
      
      term_order <- unique(unlist(lapply(est, function(x) x$term)))
      statistic_order <- unique(unlist(lapply(est, function(x) x$statistic)))
      bycols <- c(list(c(shape$group_name, "group", "term", "statistic")), 
                  lapply(est, colnames))
      bycols <- Reduce(intersect, bycols)
      f <- function(x, y) {
        merge(x, y, all = TRUE, sort = FALSE, by = bycols)
      }
      est <- Reduce(f, est)
      if (is.null(shape$group_name)) {
        idx <- paste(est$term, est$statistic)
        if (anyDuplicated(idx) > 0) {
          candidate_groups <- sapply(msl, function(x) colnames(x[["tidy"]]))
          candidate_groups <- unlist(candidate_groups)
          candidate_groups <- setdiff(candidate_groups, c("term", 
                                                          "type", 
                                                          "estimate", 
                                                          "std.error",
                                                          "conf.level", 
                                                          "conf.low", 
                                                          "conf.high",
                                                          "statistic", 
                                                          "df.error", 
                                                          "p.value"))
          msg <- c("There are duplicate term names in the table.", 
          "The `shape` argument of the `modelsummary` function can be used 
          to print related terms together. The `group_map` argument can be used 
          to reorder, subset, and rename group identifiers. See `?modelsummary` 
          for details.", 
          "You can find the group identifier to use in the `
          shape` argument by 
          calling `get_estimates()` on one of your models. Candidates include:", 
          paste(candidate_groups, collapse = ", "))
          
          insight::format_warning(msg)
        }
      }
      
      est <- shape_estimates(est, shape, conf_level = conf_level, 
                             statistic = statistic, estimate = estimate)
      
      est$part <- "estimates"
      est <- est[, unique(c("part", names(est)))]
      est[is.na(est)] <- ""
      if ("term" %in% colnames(est)) {
        if (!is.null(coef_map)) {
          term_order <- coef_map
        }
        est$term <- factor(est$term, unique(term_order))
        if ("group" %in% colnames(est)) {
          if (!is.null(group_map)) {
            est$group <- factor(est$group, group_map)
          }
          else {
            est$group <- factor(est$group, unique(est$group))
          }
        }
      }
      else if ("model" %in% colnames(est)) {
        est$model <- factor(est$model, model_names)
      }
      if ("statistic" %in% colnames(est)) {
        est$statistic <- factor(est$statistic, statistic_order)
      }
      est <- est[do.call(order, as.list(est)), ]
      if (is.numeric(coef_omit)) {
        coef_omit <- unique(round(coef_omit))
        if (length(unique(sign(coef_omit))) != 1) {
          insight::format_error("All elements of `coef_omit` must 
                                have the same sign.")
        }
        if (!"term" %in% shape$lhs) {
          msg <- "`term` must be on the left-hand side of 
          the `shape` formula 
          when `coef_omit` is a numeric vector."
          insight::format_error(msg)
        }
        term_idx <- paste(est$group, est$term)
        if (max(abs(coef_omit)) > length(unique(term_idx))) {
          msg <- sprintf("There are %s unique terms, but `coef_omit` tried to
               omit more than that.", 
                         length(term_idx))
          insight::format_error(msg)
        }
        idx <- !term_idx %in% unique(term_idx)[abs(coef_omit)]
        if (any(coef_omit > 0)) {
          est <- est[idx, , drop = FALSE]
        }
        else {
          est <- est[!idx, , drop = FALSE]
        }
      }
      if (isTRUE(checkmate::check_character(coef_rename, names = "unnamed"))) {
        nterms <- length(unique(est$term))
        if (length(coef_rename) != nterms) {
          msg <- "`coef_rename` must be a named character vector or 
          an unnamed vector of length %s"
          insight::format_error(sprintf(msg, nterms))
        }
        dict <- stats::setNames(coef_rename, as.character(unique(est$term)))
        tmp <- replace_dict(as.character(est$term), dict)
        est$term <- factor(tmp, unique(tmp))
      }
      if (is.null(shape$group_name)) {
        est[["group"]] <- NULL
      }
      cols <- intersect(colnames(est), c("term", shape$group_name, 
                                         "model", "statistic"))
      for (col in cols) {
        est[[col]] <- as.character(est[[col]])
      }
      if (all(est[["term"]] == "cross")) {
        est[["term"]] <- NULL
      }
      gof <- list()
      for (i in seq_along(msl)) {
        if (is.data.frame(msl[[i]]$glance)) {
          gof[[i]] <- format_gof(msl[[i]]$glance, fmt = fmt, 
                                 gof_map = gof_map, ...)
          colnames(gof[[i]])[2] <- model_names[i]
        }
        else {
          gof[[i]] <- NULL
        }
      }
      
      f <- function(x, y) merge(x, y, all = TRUE, sort = FALSE, 
                                by = "term")
      gof <- Reduce(f, gof)
      gof <- map_gof(gof, gof_omit, gof_map)
      tab <- bind_est_gof(est, gof)
      tab[is.na(tab)] <- ""
      if (is.null(coef_map) && isFALSE(coef_rename) && "term" %in% 
          colnames(tab) && !identical(output_format, "rtf")) {
        idx <- tab$part != "gof"
        # tab$term <- ifelse(idx, gsub("::", " = ", tab$term), 
        #                    tab$term)
        # tab$term <- ifelse(idx, gsub(":", " × ", tab$term), 
        #                    tab$term)
      }
      hrule <- match("gof", tab$part)
      if (!is.na(hrule) && !is.null(add_rows) && !is.null(attr(add_rows, 
                                                               "position"))) {
        pos <- attr(add_rows, "position")
        if (identical(pos, "gof_start")) {
          pos <- match("gof", tab$part)
          pos <- pos:(pos + nrow(add_rows) - 1)
          attr(add_rows, "position") <- pos
        }
        else if (identical(pos, "gof_end")) {
          attr(add_rows, "position") <- NULL
        }
        else if (identical(pos, "coef_start")) {
          pos <- seq_len(nrow(add_rows))
          attr(add_rows, "position") <- pos
          hrule <- hrule + sum(length(pos < hrule - 1))
        }
        else if (identical(pos, "coef_end")) {
          pos <- match("gof", tab$part)
          pos <- pos:(pos + nrow(add_rows) - 1)
          attr(add_rows, "position") <- pos
          hrule <- hrule + sum(pos < hrule - 1)
        }
        else {
          hrule <- hrule + sum(pos < hrule - 1)
        }
      }
      if (is.na(hrule)) {
        hrule <- NULL
      }
      
      stars_note <- settings_get("stars_note")
      if (isTRUE(stars_note) && !isFALSE(stars) && !any(grepl("\\{stars\\}", 
                                                              c(estimate, 
                                                                statistic)))) {
        stars_note <- make_stars_note(stars, output_format = output_format, 
                                      output_factory = output_factory)
        if (is.null(notes)) {
          notes <- stars_note
        }
        else {
          notes <- c(stars_note, notes)
        }
      }
      if (settings_equal("function_called", "modelsummary_rbind")) {
        tab <- redundant_labels(tab, "term")
      }
      
      term_label <- getOption("modelsummary_model_labels_term", 
                              default = NULL)
      model_label <- getOption("modelsummary_model_labels_model", 
                               default = NULL)
      group_label <- getOption("modelsummary_model_labels_group", 
                               default = NULL)
      if (!identical(output_format, "dataframe") && 
          !settings_equal("function_called", 
                          "modelsummary_rbind")) {
        dups <- c("term", "model", shape$group_name)
        for (d in dups) {
          tab <- redundant_labels(tab, d)
        }
        tab$statistic <- tab$part <- NULL
        if ("term" %in% colnames(tab)) {
          if (is.null(term_label)) {
            colnames(tab)[colnames(tab) == "term"] <- "       "
          }
          else {
            colnames(tab)[colnames(tab) == "term"] <- term_label
          }
        }
        if ("model" %in% colnames(tab)) {
          if (is.null(model_label)) {
            colnames(tab)[colnames(tab) == "model"] <- "         "
          }
          else {
            colnames(tab)[colnames(tab) == "model"] <- model_label
          }
        }
      }
      if (identical(output_format, "dataframe") || 
          settings_equal("function_called", 
                         "modelsummary_rbind")) {
        if (!is.null(term_label) && "term" %in% colnames(tab)) {
          colnames(tab)[colnames(tab) == "term"] <- term_label
        }
        if (!is.null(model_label) && "model" %in% colnames(tab)) {
          colnames(tab)[colnames(tab) == "model"] <- model_label
        }
        if (!is.null(group_label) && "group" %in% colnames(tab)) {
          colnames(tab)[colnames(tab) == "group"] <- group_label
        }
      }
      if (length(unique(tab$group)) == 1) {
        tab$group <- NULL
      }
      tmp <- setdiff(shape$lhs, c("model", "term"))
      if (length(tmp) == 0) {
        tab$group <- NULL
      }
      else if (!identical(output_format, "dataframe")) {
        if (is.null(group_label)) {
          colnames(tab)[colnames(tab) == "group"] <- "        "
        }
        else {
          colnames(tab)[colnames(tab) == "group"] <- group_label
        }
      }
      if (is.null(align)) {
        stub_labels <- c(" ", shape$group_name)
        if (!is.null(term_label)) {
          stub_labels <- c(stub_labels, term_label)
        }
        if (!is.null(model_label)) {
          stub_labels <- c(stub_labels, model_label)
        }
        if (!is.null(group_label)) {
          stub_labels <- c(stub_labels, group_label)
        }
        n_stub <- sum(grepl("^ *$", colnames(tab))) + sum(colnames(tab) %in% 
                                                            stub_labels)
        align <- paste0(strrep("l", n_stub), strrep("c", ncol(tab) - 
                                                      n_stub))
        if (isTRUE(checkmate::check_data_frame(add_columns))) {
          align <- paste0(align, strrep("c", ncol(add_columns)))
        }
      }
      flag <- any(sapply(models, inherits, c("marginaleffects", "comparisons", 
                                             "marginalmeans")))
      if (isTRUE(flag)) {
        colnames(tab) <- gsub("^value$", " ", colnames(tab))
        colnames(tab) <- gsub("^contrast_", "", colnames(tab))
        colnames(tab) <- gsub("^contrast$", " ", colnames(tab))
      }
      for (i in seq_along(tab)) {
        tab[[i]] <- gsub("\\(\\s*\\)", "", tab[[i]])
        tab[[i]] <- gsub("\\(\\\\num\\{NA\\}\\)", "", tab[[i]])
        tab[[i]] <- gsub("\\[,\\s*\\]", "", tab[[i]])
        tab[[i]] <- gsub("\\[\\\\num\\{NA\\}, \\\\num\\{NA\\}\\]", 
                         "", tab[[i]])
      }
      idx <- apply(tab, 1, function(x) any(x != ""))
      tab <- tab[idx, ]
      out <- factory(tab, align = align, fmt = fmt, hrule = hrule, 
                     gof_idx = hrule, notes = notes, output = output, 
                     title = title, 
                     add_rows = add_rows, add_columns = add_columns, 
                     escape = escape, 
                     output_factory = output_factory, 
                     output_format = output_format, 
                     output_file = output_file, ...)
      attr(out, "backend") <- attr(msl, "backend")
      out <- set_span_cbind(out, span_cbind)
      if (settings_equal("function_called", "modelsummary_rbind")) {
        return(out)
      }
      else if (!is.null(output_file) || isTRUE(output == "jupyter") || 
               (isTRUE(output == "default") && settings_equal("output_default", 
                                                              "jupyter"))) {
        settings_rm()
        return(invisible(out))
      }
      else {
        settings_rm()
        return(out)
      }
    } 
    
    allow_unlock_replace_bind <- FALSE
    if(!inherits(models[[1]], "bgmfit")) {
      if(is.null(models$test_mode)) {
        allow_unlock_replace_bind <- TRUE
      } else if(!models$test_mode) {
        allow_unlock_replace_bind <- TRUE
      }
    } else { 
      test_mode_c <- c()
      for (modelsi in 1:length(models)) {
        test_mode_c <- c(test_mode_c, models[[modelsi]]$test_mode)
      }
      if(is.null(test_mode_c)) test_mode_c <- FALSE
      if(all(!test_mode_c)) allow_unlock_replace_bind <- TRUE
    }

    if(allow_unlock_replace_bind) {
      unlock_replace_bind(package = "parameters", what = "parameters",
                          replacement = custom_model_parameters.brmsfit,
                          ept_str = T)

      unlock_replace_bind(package = "modelsummary", what = "modelsummary",
                          replacement = custom_modelsummary,
                          ept_str = T)

      unlock_replace_bind(package = "modelsummary", what = "get_gof",
                          replacement = custom_get_gof,
                          ept_str = T)
    } 
    
    modelsummary::modelsummary(models, 
                               effects = effects, 
                               output = output,
                               fmt = fmt,
                               estimate = estimate,
                               statistic = statistic,
                               vcov = vcov,
                               conf_level = conf_level,
                               exponentiate = exponentiate,
                               stars = stars,
                               shape = shape,
                               coef_map = coef_map,
                               coef_omit = coef_omit,
                               coef_rename = coef_rename,
                               gof_map = gof_map,
                               gof_omit = gof_omit,
                               # gof_function = custom_gof_function,
                               gof_function = gof_function,
                               group_map = group_map,
                               # add_columns = paramsx$ESS %>% data.frame(),
                               add_columns = add_columns,
                               add_rows = add_rows,
                               align = align,
                               notes = notes,
                               title = title,
                               escape = escape,
                               add_component = add_component,
                               # parameters::model_parameters
                               # model,
                               centrality = centrality,
                               dispersion = dispersion,
                               ci = ci,
                               ci_method = ci_method,
                               test = test,
                               rope_range = rope_range,
                               rope_ci = rope_ci,
                               bf_prior = bf_prior,
                               diagnostic = diagnostic,
                               priors = priors,
                               component = component,
                               standardize = standardize,
                               group_level = group_level,
                               keep = keep,
                               drop = drop,
                               verbose = verbose,
                               # extra_arg
                               sigma_exp = sigma_exp,
                               sort_gr_by = sort_gr_by,
                               gsub_fixed = gsub_fixed,
                               gsub_random = gsub_random,
                               # pull_bsigma = pull_bsigma,
                               add_diagnostic = add_diagnostic,
                               add_pd = add_pd,
                               add_section = add_section,
                               sort_corr = sort_corr,
                               custom_rename = custom_rename,
                               custom_rename_map = custom_rename_map,
                               envir = environment(),
                               ...)
  } # End 


