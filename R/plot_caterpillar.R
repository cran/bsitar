


#' Caterpillar plot for bsitar random effects
#'
#' Create a caterpillar plot of group-level effects from a `bsitar` model using
#' posterior summaries returned by [brms::ranef()]. When multiple plots are
#' requested, plots are arranged and combined into a single display with
#' \pkg{patchwork}.
#'
#' @param model A \code{bsitar} object.
#'
#' @param group A character string specifying the grouping factor to plot. This
#'   is typically the individual identifier, such as \code{id}. If \code{group =
#'   NULL} (the default), it is determined automatically from the model object.
#'   In most cases, \code{idvar} and \code{levels_id} are used to define the
#'   grouping structure, so \code{group} will rarely need to be specified
#'   manually.
#'
#' @param groupby A character string specifying the factor within which
#'   \code{group} is nested. This corresponds to the \code{by} argument in
#'   \code{gr(..., by = xxx)}, where \code{xxx} is the grouping variable. If
#'   \code{groupby = NULL} (the default), it is determined automatically from
#'   the model object.
#'
#' @param effect A character string specifying the random-effect coefficient to
#'   plot. If \code{NULL}, the first available effect is used.
#'
#' @param transform_parameter Optional character vector of parameter names to
#'   transform. Each supplied value must match an available random-effect
#'   parameter.
#'
#' @param transform_fun Optional function or list of functions used to transform
#'   the selected parameters. Accepted forms are:
#'   \itemize{
#'     \item A single function, applied to all targeted combinations.
#'     \item A named list of functions with one function per
#'     \code{transform_parameter}.
#'     \item A named list of functions whose names match the internal
#'     transformation keys used by the method.
#'   }
#'   Each function must accept a numeric vector and return a numeric vector of
#'   the same length.
#'   
#' @param pdata An optional parameter data frame in wide format where parameter
#'   name and key variable such as \code{Estimate} are named as
#'   \code{effect.Estimate}. Currently ignored. Included here for future
#'   expansion.
#'
#' @param probs Numeric vector of length 2 giving the lower and upper credible
#'   interval probabilities passed to [brms::ranef()]. Default is
#'   \code{c(0.025, 0.975)}.
#'
#' @param sort Logical; if \code{TRUE}, levels are ordered by the posterior mean
#'   estimate. If \code{FALSE}, the original order is retained.
#'
#' @param combine Logical; if \code{TRUE} and more than one plot is generated,
#'   plots are combined with [patchwork::wrap_plots()]. If \code{FALSE}, the
#'   function returns the individual plot objects without combining them.
#'
#' @param point_color Color used for the posterior mean points. When \code{NULL}
#'   (default), the color is set internally to \code{"#1f4e79"} if
#'   \code{groupby} is \code{NULL}; otherwise, distinct colors are used for each
#'   level of the \code{groupby} variable.
#'
#' @param point_size Size used for the posterior mean points. When \code{NULL}
#'   (default), the size is set internally to \code{1.5}.
#'
#' @param interval_color Color used for the posterior interval lines. When
#'   \code{NULL} (default), the color is set internally to \code{"#1f4e79"} if
#'   \code{groupby} is \code{NULL}; otherwise, distinct colors are used for each
#'   level of the \code{groupby} variable.
#'
#' @param patch_plot.margin Optional margin specification passed to
#'   \pkg{patchwork}.
#'
#' @param refline Numeric value for the vertical reference line. The default is
#'   \code{NULL}, which is internally treated as \code{0}. Note that the same
#'   transformation is applied to \code{refline} as specified by
#'   \code{transform_parameter}. For example, if the \code{parameter} \code{'c'}
#'   is transformed using \code{'exp'}, then \code{refline = 0} is transformed
#'   to \code{refline = exp(0)}, that is, \code{refline = 1}.
#'
#' @param dodge Width used to control the horizontal thickness of the interval
#'   bars.
#'
#' @param label.title Optional character string used as the plot title. Default
#'   is \code{NULL}.
#'
#' @param ticks.x Optional logical value. If \code{FALSE}, x-axis tick labels
#'   are suppressed. Default is \code{NULL}, which is treated internally as
#'   \code{FALSE}.
#'
#' @param ticks.y Optional logical value. If \code{FALSE}, y-axis tick labels
#'   are suppressed. Default is \code{NULL}, which is treated internally as
#'   \code{FALSE}.
#'   
#' @param facet_ncol Optional integer set the number of columns for the
#'   [ggplot2::facet_wrap()]. Default is \code{NULL}, which allows internally
#'   setting the \code{facet_ncol} based on the number of plots.
#'   
#' @param facet_nrow Optional integer set the number of rows for the
#'   [ggplot2::facet_wrap()]. Default is \code{NULL}, which allows internally
#'   setting the \code{facet_ncol} based on the number of plots.
#'   
#' @param wrap_plots_ncol Optional integer set the number of columns for the
#'   [patchwork::wrap_plots()]. Default is \code{NULL}, which allows internally
#'   setting the \code{wrap_plots_ncol} based on the number of plots.
#'
#' @param wrap_plots_nrow Optional integer set the number of rows for the
#'   [patchwork::wrap_plots()]. Default is \code{NULL}, which allows internally
#'   setting the \code{wrap_plots_nrow} based on the number of plots.
#'
#' @inheritParams plot_curves.bgmfit
#' @inheritParams ggplot2::facet_wrap
#'
#' @return A \code{ggplot} object showing a caterpillar plot of the selected
#'   random effect. 
#'
#' @details
#' The function extracts group-level summaries using [brms::ranef()] with
#' \code{summary = TRUE}, then plots posterior means and credible intervals for
#' the selected grouping factor and effect. If transformations are requested,
#' they are applied before summarizing the draws so that the plotted summaries
#' are on the transformed scale.
#'
#' @rdname plot_caterpillar
#' @export
#' 
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model is fit  
#' # to the 'berkeley_exdata' and saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether the model fit object 'berkeley_exfit' exists
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Individual-specific size (parameter a) 
#' plot_caterpillar(model, group = "id", effect = "a_Intercept")
#' 
#' # Individual-specific timing (parameter b) 
#' 
#' plot_caterpillar(model, group = "id", effect = "b_Intercept")
#' 
#' # Individual-specific intensity (parameter c) 
#' 
#' plot_caterpillar(model, group = "id", effect = "c_Intercept")
#' 
#' 
#' }
#' 
plot_caterpillar.bgmfit <- function(
    model = NULL,
    group = NULL,
    groupby = NULL,
    effect = NULL,
    transform_parameter = NULL,
    transform_fun = NULL,
    pdata = NULL, 
    probs = c(0.025, 0.975),
    conf = NULL,
    robust = FALSE,
    sort = TRUE,
    resp = NULL,
    dpar = NULL,
    levels_id = NULL,
    idvar = NULL,
    layout = NULL,
    linecolor = NULL,
    linecolor1 = NULL,
    linecolor2 = NULL,
    label.x = NULL,
    label.y = NULL,
    ticks.x = NULL,
    ticks.y = NULL,
    label.title = NULL,
    label.subtitle = NULL,
    legendpos = 'bottom',
    linewidth.main = NULL,
    linetype.groupby = NA,
    color.groupby = NULL,
    fill.groupby = NULL,
    point_size = NULL,
    point_color = NULL,
    interval_color = NULL,
    combine = TRUE,
    ncol = NULL,
    wrap_plots_ncol = NULL,
    wrap_plots_nrow = NULL,
    facet_ncol = NULL,
    facet_nrow = NULL,
    scales = "fixed",
    space  = "fixed",
    patch_plot.margin = ggplot2::margin(10, 0, 10, 0),
    each_object = FALSE,
    refline = NULL,
    dodge = 0.4,
    envir = NULL,
    verbose = FALSE,
    ...) {
  
  
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- envir
  }

  gsubfunction <- function(x) {
    xx <- x
    if(grepl("^sigma_", xx)) {
      xx <- gsub("_Intercept", "_Intercept", xx, fixed = T)
    } else if(grepl("^b_", xx)) {
      xx <- gsub("_Intercept", "_Intercept", xx, fixed = T)
    } else if(grepl("^sd_", xx)) {
      xx <- gsub("_Intercept", "", xx, fixed = T)
    } else if(grepl("^cor_", xx)) {
      xx <- gsub("_Intercept", "", xx, fixed = T)
    }
    xx <- gsub("classid", "", xx, fixed = T)
    xx <- gsub("classClassI", "ClassI", xx, fixed = T)
    xx <- gsub("classClassII", "ClassII", xx, fixed = T)
    xx <- gsub("__", "_", xx, fixed = T)
    xx <- gsub("_", " ", xx, fixed = T)
    xx <- gsub(":", ": ", xx, fixed = T)
    xx <- gsub("ClassII", "Class II", xx, fixed = T)
    xx <- gsub("ClassI", "Class I", xx, fixed = T)
    xx
  }
  
  my_labeller <- ggplot2::as_labeller(
    function(x) {
      gsubfunction(x)
    },
    default = ggplot2::label_parsed
  )
  
  summarize_randomeffects <- function(arr, 
                                      idvar, 
                                      transform_parameter, 
                                      transform_fun,
                                      probs = c(0.025, 0.975), 
                                      robust = FALSE,
                                      verbose = FALSE) {
    
    arr <- arr[[idvar]]
    dn <- dimnames(arr)
    nlev <- dim(arr)[2]
    ncoef <- dim(arr)[3]
    transform_map <- make_transform_map(transform_parameter, transform_fun)
    x0 <- arr[, , 1]
    prefix0 <- sub("_.*$", "", dn[[3]][1])
    if (prefix0 %in% names(transform_map)) x0 <- transform_map[[prefix0]](x0)
    s0 <- brms::posterior_summary(x0, probs = probs, robust = robust)
    stat_names <- colnames(s0)
    out <- array(
      NA_real_,
      dim = c(nlev, length(stat_names), ncoef),
      dimnames = list(dn[[2]], stat_names, dn[[3]])
    )
    for (j in seq_len(ncoef)) {
      coef_name <- dn[[3]][j]
      prefix <- sub("_.*$", "", coef_name)
      x <- arr[, , j]
      if (prefix %in% names(transform_map)) {
        x <- transform_map[[prefix]](x)
      }
      s <- brms::posterior_summary(x, probs = probs, robust = robust)
      if (all(dim(s) == c(nlev, length(stat_names)))) {
        out[, , j] <- s
      } else if (all(dim(s) == c(length(stat_names), nlev))) {
        out[, , j] <- t(s)
      } else {
        stop2c(sprintf(
          "Coefficient %s returned %s x %s, expected %s x %s (or transpose)",
          coef_name, dim(s)[1], dim(s)[2], nlev, length(stat_names)
        ))
      }
    }
    outlist <- list()
    outlist[[idvar]] <- out
    return(outlist)
  }
  
  make_transform_map <- function(transform_parameter, transform_fun) {
    id_fun <- function(x) x
    if (is.null(transform_fun)) {
      transform_fun <- stats::setNames(
        rep(list(id_fun), length(transform_parameter)),
        transform_parameter
      )
    } else if (is.function(transform_fun)) {
      transform_fun <- stats::setNames(
        rep(list(transform_fun), length(transform_parameter)),
        transform_parameter
      )
    } else if (!is.list(transform_fun)) {
      stop2c("transform_fun must be NULL, a function, or 
             a named list of functions.")
    }
    stats::setNames(
      lapply(transform_parameter, function(nm) {
        fn <- transform_fun[[nm]]
        if (is.null(fn)) id_fun else get_fun_form(fn)
      }),
      transform_parameter
    )
  }
  
  check_and_set_effect_loop <- function(effect_loop, resp) {
    if(is.null(resp)) return(effect_loop)
    check_resp_prefix <- paste0(resp, "_")
    effect_loopi_c <- c()
    for (effect_loopi in effect_loop) {
      if(!grepl(check_resp_prefix, effect_loopi)) {
        effect_loopiset <- paste0(check_resp_prefix, effect_loopi)
      } else {
        effect_loopiset <- effect_loopi
      }
      effect_loopi_c <- c(effect_loopi_c, effect_loopiset)
    }
    effect_loop <- effect_loopi_c
    return(effect_loop)
  }
  
  
  if(is.null(pdata)) {
    if(is.null(refline)) refline <- 0
    if(is.null(model)) stop2c("Argument 'model' is required")
    if(is.null(dpar)) {
      dpar <- "mu"
    }
    
    model <- getmodel_info(model = model, 
                           dpar = dpar, 
                           resp = resp, 
                           deriv = NULL, 
                           verbose = verbose)
    
    estimate <- NULL;
    idrow <- NULL;
    lower <- NULL;
    upper <- NULL;
    . <- NULL;
    
    xcall <- match.call()
    match.call.list.in <- as.list(match.call())[-1]
    setxcall_ <- match.call()
    post_processing_checks_args <- list()
    post_processing_checks_args[['model']]    <- model
    post_processing_checks_args[['xcall']]    <- setxcall_
    post_processing_checks_args[['resp']]     <- resp
    post_processing_checks_args[['envir']]    <- envir
    post_processing_checks_args[['check_d0']] <- FALSE
    post_processing_checks_args[['check_d1']] <- TRUE
    post_processing_checks_args[['check_d2']] <- FALSE
    o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
    rlang_trace_back <- rlang::trace_back()
    check_trace_back.bgmfit <- grepl(".bgmfit", rlang_trace_back[[1]])
    if(all(!check_trace_back.bgmfit)) {
      # 
    } else {
      rlang_trace_back.bgmfit_i <- min(which(check_trace_back.bgmfit == TRUE))
      rlang_trace_back.bgmfit <- 
        rlang_trace_back[[1]][[rlang_trace_back.bgmfit_i]]
      rlang_call_name <- rlang::call_name(rlang_trace_back.bgmfit)
      xcall <- rlang_call_name
    }
    
    check_if_package_installed(model, xcall = xcall)
    model$xcall            <- xcall
    arguments              <- get_args_(match.call.list.in, xcall)
    arguments$model        <- model
    arguments$dpar         <- dpar
    
    if(is.null(envir)) {
      arguments$envir <- envir <- parent.frame()
    }
    
    if(!is.null(conf)) {
      probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    }
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', 'Est.Error', probtitles)
    
    
    future <- FALSE
    future_session <- 'multisession'
    
    cores <- 1
    get.cores_ <- get.cores(arguments$cores)
    arguments$cores <- cores <-  get.cores_[['max.cores']] 
    .cores_ps <- get.cores_[['.cores_ps']]
    
    if (future) {
      if(is.null(cores)) stop2c("Please set the number of cores for 'future' by  
                            using the the 'cores' argument, e.g. cores = 4")
      if (future_session == 'multisession') {
        future::plan('multisession', workers = cores)
      } else if (future_session == 'multicore') {
        future::plan('multicore', workers = cores)
      }
    }
    
    newdata <- NULL
    if (is.null(newdata)) {
      newdata <- model$model_info$bgmfit.data
    } else {
      newdata <- newdata
    } 
    
    if (is.null(resp)) {
      resp_rev_ <- resp
    } else if (!is.null(resp)) {
      resp_rev_ <- paste0("_", resp)
    }
    
    groupvar_     <- paste0('groupvar', resp_rev_)
    yvar_         <- paste0('yvar', resp_rev_)
    yvar          <- model$model_info[[yvar_]]
    hierarchical_ <- paste0('hierarchical', resp_rev_)
    if(is.null(levels_id) & is.null(idvar)) {
      idvar <- model$model_info[[groupvar_]]
      if (!is.null(model$model_info[[hierarchical_]])) {
        idvar <- model$model_info[[hierarchical_]]
      }
      model$model_info[[groupvar_]] <- idvar # idvar[1]
    } else if (!is.null(levels_id)) {
      idvar <- levels_id
    } else if (!is.null(idvar)) {
      idvar <- idvar
    }
    if(is.null(idvar)) {
      if(is.null(idvar)) {
        if(!is.null(model$model_info[['idvars']])) {
          idvar <- model$model_info[['idvars']]
        }
      }
    }

    if(!is_emptyx(model$ranef$by)) {
      gr_by       <- model$ranef$by
      gr_bylevels <- model$ranef$bylevels
      sort_gr_bylevels <- unlist(gr_bylevels) %>% unique()
      by_order    <- paste0(gr_by, unlist(gr_bylevels))
      by_order    <- unique(by_order)
      gr_by <- unique(gr_by) 
    } else {
      gr_by    <- NULL
      by_order    <- NULL
      sort_gr_bylevels <- NULL
    }
    
    if(!is.null(group)) idvar <- group
   
    if(is.null(groupby)) groupby <- gr_by
    
    
    if(is.null(transform_fun)) {
      transform_fun <- function(x)x
    }
    
    if(is.null(transform_parameter)) {
      re_list <- brms::ranef(model, summary = TRUE, probs = probs, 
                             robust = robust)
    } else {
      re_list <- brms::ranef(model, summary = FALSE)
      re_list <- 
        summarize_randomeffects(arr = re_list, idvar = idvar,
                                probs = probs, robust = robust,
                                transform_parameter = transform_parameter,
                                transform_fun = transform_fun,
                                verbose = FALSE)
    }
    
    if (!idvar %in% names(re_list)) {
      stop2c("`Individual identifier` not found. Available groups: ",
           paste(names(re_list), collapse = ", "))
    }
    
    re_arr <- re_list[[idvar]]
    
    available_effects <- dimnames(re_arr)[[3]]
    
    if (is.null(effect)) {
      default_effect <- "a_Intercept"
      default_effect <- check_and_set_effect_loop(default_effect, resp = resp)
      if (default_effect %in% available_effects) effect_loop <- default_effect
    } else if(length(effect) == 1) {
      if(effect == 'all') {
        effect_loop <- available_effects
      } else {
        effect_loop <- effect
      }
    } else {
      effect_loop <- effect
    }
  } # if(is.null(pdata)) {

  
  if(!is.null(pdata)) {
    probtitles <- names(pdata)[grepl("^Q.*[0-9]$", names(pdata))]
    if(is_emptyx(probtitles)) {
      stop2c("The 'pdata' must have probability tiles that start with Q such 
              'Q2.5' and 'Q97.5'")
    }
    
    if(is.null(effect)) 
      stop2c("For 'pdata', argument 'effect' must be specified")
    if(is.null(group)) 
      stop2c("For 'pdata', argument 'group' must be specified")
    if("Parameter" %in% names(pdata)) {
      namevar <- "Parameter"
    } else if("parameter" %in% names(pdata)) {
      namevar <- "parameter"
    } else {
      stop2c("'pdata' must have column 'Parameter'")
    }
    
    available_effects <- unique(pdata[[namevar]])
    
    if (is.null(effect)) {
      default_effect <- "APGV"
      default_effect <- check_and_set_effect_loop(default_effect, resp = resp)
      if (default_effect %in% available_effects) effect_loop <- default_effect
    } else if(length(effect) == 1) {
      if(effect == 'all') {
        effect_loop <- available_effects
      } else {
        effect_loop <- effect
      }
    } else {
      effect_loop <- effect
    }
    
    idvar   <- group # "id"
    valsvar  <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    
    pdata <- pdata %>%
      tidyr::pivot_wider(
        id_cols = dplyr::all_of(idvar),
        names_from = dplyr::all_of(namevar),
        values_from = dplyr::all_of(valsvar),
        # names_glue = paste0("{", namevar, "}.{.value}")
        names_glue = paste0("{.value}.{", namevar, "}")
      )
    pdata <- pdata %>% dplyr::rename('idrow' = dplyr::all_of(idvar))
    
    transform_parameter <- NULL
    transform_fun <- NULL
    if(is.null(transform_fun)) {
      transform_fun <- function(x)x
    }
  }
  

  if (is.null(ticks.x)) {
    ticks.x     <- FALSE
  }
  if (is.null(ticks.y)) {
    ticks.y     <- TRUE
  }
  
  if (is.null(legendpos)) {
    legendpos <- "none"
  } else if (!is.null(legendpos)) {
    legendpos <- legendpos
  }
  
  if (is.null(linewidth.main)) {
    linewidth.main <- 0.4
  }
  
  if(is.null(point_size)) {
    point_size <- 1.5
  }
  
  groupby_str_d <- NULL
  groupby_str_d <- unique(c(groupby_str_d, groupby))
  
  band.alpha <- NULL
  if (is.null(band.alpha)) {
    band.alpha <- 1
  }
  band.legends <- TRUE
  
  groupby_str_d_o <- groupby_str_d
  groupby_str_d <- c('idrow', groupby_str_d)
  
  if(is.null(color.groupby)) {
    color.groupby <- groupby 
  }
  
  if(is.null(fill.groupby)) {
    fill.groupby <- color.groupby
  } else if(!is.null(fill.groupby)) {
    if(isFALSE(fill.groupby)) fill.groupby <- NA
  }
  
  if(is.null(point_color)) {
    if(is.null(groupby)) {
      point_color <- "#1f4e79" 
    } else {
      point_color <- NULL # color.groupby
    }
  }
  
  if(is.null(interval_color)) {
    if(is.null(groupby)) {
      interval_color <- "#1f4e79" 
    } else {
      interval_color <- NULL #color.groupby
    }
  }

  if(is.null(layout)) {
    if(is.null(groupby)) {
      if (length(effect_loop) > 2) layout <- 'facet' else layout <- 'single'
    } else if(!is.null(groupby)) {
      layout <- 'single'
    }
  }
  

  plot_list <- list()
  for (effect in effect_loop) {
    if (!effect %in% available_effects) {
      stop2c("`effect` not found. Available effects: ", paste(available_effects, 
                                                            collapse = ", "))
    }

    if(is.null(pdata)) {
      dframe <- as.data.frame(re_arr[, , effect, drop = FALSE])
      dframe$idrow <- rownames(re_arr)
      if(!is.null(attr(dframe$idrow, "by"))) {
        dframe <- dframe %>% 
          dplyr::mutate(!!as.name(gr_by) := attr(dframe$idrow, "by"))
      }
      Estimate.effect <- paste0(set_names_[1], '.', effect)
    } else if(!is.null(pdata)) {
      dframe <- pdata %>% dplyr::select(dplyr::ends_with(paste0(".", effect)))
      dframe[['idrow']] <- pdata[['idrow']]
      set_names_ <- valsvar
      gr_by <- groupby
      Estimate.effect <- paste0(set_names_[1], '.', effect)
      refline <- mean(dframe[[Estimate.effect]])
    }

    if(is.null(gr_by)) {
      if (sort) {
        dframe <- dframe %>% 
          dplyr::mutate(idrow = 
                          forcats::fct_reorder(.f = .data[['idrow']], 
                                               .x = .data[[Estimate.effect]], 
                                               .desc = FALSE))
      } else {
        dframe <- dframe %>% 
          dplyr::mutate(idrow = 
                          forcats::fct_reorder(.f = .data[['idrow']], 
                                               .x = .data[[Estimate.effect]], 
                                               .desc = TRUE))
      }
    } else if(!is.null(gr_by)) {
      if (sort) {
        dframe <- dframe %>% 
          dplyr::group_by(.data[[gr_by]]) %>% 
          dplyr::mutate(idrow = 
                          forcats::fct_reorder(.f = .data[['idrow']], 
                                               .x = .data[[Estimate.effect]], 
                                               .desc = FALSE))
      } else {
        dframe <- dframe %>% 
          dplyr::group_by(.data[[gr_by]]) %>% 
          dplyr::mutate(idrow = 
                          forcats::fct_reorder(.f = .data[['idrow']], 
                                               .x = .data[[Estimate.effect]], 
                                               .desc = TRUE))
      }
    }
    
    dframe <-
      dframe %>% dplyr::ungroup() %>% 
      dplyr::mutate(
        groupby = interaction(dplyr::across(dplyr::all_of(groupby_str_d)))
      )
    
    if (is.null(label.x)) {
      label.x     <- paste0("Estimate", "")
    }
    if (is.null(label.y)) {
      label.y     <- paste0(idvar, "")
    }

    label.title     <- effect
   
    if(layout == 'single') {
      if(!is.null(gr_by)) {
        dframe <- dframe %>%
          dplyr::mutate(!!as.name(gr_by) := forcats::fct_relevel(
            .data[[gr_by]])
            ) %>%
          dplyr::arrange(.data[[gr_by]], .data[[Estimate.effect]]) %>%
          dplyr::mutate(idrow = factor(idrow, levels = unique(idrow)))
      }
    }
    
    transform_map <- make_transform_map(transform_parameter, transform_fun)
    
    prefix0 <- sub("_.*$", "", effect)

    if (prefix0 %in% names(transform_map)) {
      refline_mark <- transform_map[[prefix0]](refline)
    } else {
      refline_mark <- refline
    }
   
    
    if(is.null(interval_color)) {
      plot.o.d <- dframe %>%
        ggplot2::ggplot(., ggplot2::aes(x = .data[[Estimate.effect]], 
                                        y = .data[['idrow']])) +
        ggplot2::geom_vline(xintercept = refline_mark, linetype = 2, 
                            color = "grey50") +
        ggplot2::geom_errorbarh(
          ggplot2::aes(
            xmin = .data[[paste0(probtitles[1], '.', effect)]],
            xmax = .data[[paste0(probtitles[2], '.', effect)]],
            group = groupby,
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, gr_by)
          ),
          width = dodge / 2,
          # color = interval_color,
          linewidth = linewidth.main,
          alpha = band.alpha, show.legend = band.legends
        ) +
        ggplot2::geom_point(
          ggplot2::aes(color = line_color_key(.data, gr_by)),
          size = point_size
        ) 
    } else if(!is.null(interval_color)) {
      plot.o.d <- dframe %>%
        ggplot2::ggplot(., ggplot2::aes(x = .data[[Estimate.effect]], 
                                        y = .data[['idrow']])) +
        ggplot2::geom_vline(xintercept = refline_mark, linetype = 2, 
                            color = "grey50") +
        ggplot2::geom_errorbarh(
          ggplot2::aes(
            xmin = .data[[paste0(probtitles[1], '.', effect)]],
            xmax = .data[[paste0(probtitles[2], '.', effect)]],
            group = groupby,
            linetype = line_color_key(.data, linetype.groupby),
            color = line_color_key(.data, gr_by)
          ),
          width = dodge / 2,
          color = interval_color,
          linewidth = linewidth.main,
          alpha = band.alpha, show.legend = band.legends
        ) +
        ggplot2::geom_point(
          color = point_color,
          size = point_size
        ) 
    }
    
    plot.o.d <- plot.o.d +
      ggplot2::labs(
        x = label.x,
        y = label.y,
        title = label.title,
        subtitle = label.subtitle
      ) 
    
    plot.o.d <- plot.o.d +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank()
      )
    
    plot.o.d <- plot.o.d +
      jtools::theme_apa(legend.pos = legendpos) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5))
    
    if(!band.legends) plot.o.d <- plot.o.d + ggplot2::guides(fill = "none")
    
    if(!ticks.x) {
      plot.o.d <- plot.o.d +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
      
    }
    
    if(layout == 'facet') {
      if(!is.null(gr_by)) {
        plot.o.d <- plot.o.d + ggplot2::facet_wrap(~ .data[[gr_by]], 
                                                   ncol = facet_ncol,
                                                   nrow = facet_nrow,
                                                   space = space,
                                                   scales = scales)
      }
    }
    
    plot_list[[effect]] <- plot.o.d
  } # for (effect in effect_loop) {
  
  out <- plot_list
  plots <- names(out)
  
  if (!combine || length(out) <= 1) {
    if(length(out) == 1) return(out[[1]]) else return(out) 
  }
  
  is_patchable <- vapply(
    out,
    function(x) inherits(x, c("gg", "ggplot")),
    logical(1)
  )

  patchable_plots <- out[is_patchable]
  
  if (length(patchable_plots) == 0) {
    return(out)
  }
  
  if (is.null(wrap_plots_ncol)) {
    wrap_plots_ncol <- if (length(patchable_plots) <= 2) 1 else 2
  }

  patchable_plots <- patchable_plots[plots]

  if(is.null(patch_plot.margin)) {
    combined <- patchwork::wrap_plots(patchable_plots, 
                                      ncol = wrap_plots_ncol,
                                      nrow = wrap_plots_nrow)
  } else if(!is.null(patch_plot.margin)) {
    combined <- patchwork::wrap_plots(patchable_plots, 
                                      ncol = wrap_plots_ncol,
                                      nrow = wrap_plots_nrow) &
      ggplot2::theme(plot.margin = patch_plot.margin) 
  }
  
  out_all <- list(combined = combined, plots = out)
  
  if(each_object) {
    return(out_all)
  } else {
    if(is.null(combined)) return(out) else return(combined)
  }
  
 return(invisible(NULL))
}






#' @rdname plot_caterpillar
#' @export
plot_caterpillar <- function(model, ...) {
  UseMethod("plot_caterpillar")
}


