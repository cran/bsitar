

#' QQ plot or Pearson-residual data for a `brms` model with distributional
#' `sigma`
#'
#' Computes Pearson-type standardized residuals for a `brms` model in which
#' `sigma` may depend on predictors and, optionally, group-level terms. The
#' function can return the residual data itself or construct QQ plots based on
#' either a single posterior draw or summaries across multiple posterior draws.
#' 
#' @details
#' Internally, the function:
#' \itemize{
#'   \item Extracts posterior residual draws using
#'   `tidybayes::add_residual_draws()`.
#'   \item Reconstructs observation-specific `sigma` values using
#'     `brms::posterior_linpred(..., dpar = "sigma")` on the link scale and then
#'     transforms them to the response scale.
#'   \item Computes Pearson-style residuals as residual divided by the
#'     observation-specific `sigma`.
#'   \item Optionally summarizes posterior Pearson residuals across selected draws
#'     before creating a QQ plot.
#' }
#'
#' The `summary` argument controls how posterior Pearson residuals are reduced to
#' one value per observation when averaged QQ plots are requested:
#' \itemize{
#'   \item `"mean"`: uses the posterior mean Pearson residual.
#'   \item `"median"`: uses the posterior median Pearson residual.
#'   \item `"robust"`: uses the posterior median Pearson residual and also
#'     computes `MAD` as an auxiliary robustness diagnostic.
#' }
#'
#' The `qq_type` argument controls the returned object:
#' \itemize{
#'   \item `"data"`: returns the full residual object, including posterior
#'     Pearson residuals and extracted `sigma` parameters.
#'   \item `"qq"`: returns a `ggplot2` QQ plot based on a single posterior draw
#'     selected by `draw_ids_select`.
#'   \item `"qq_avg"`: returns a `ggplot2` QQ plot based on summarized Pearson
#'     residuals across posterior draws.
#' }
#'
#' @param model A fitted `brms` model object.
#' @param resp A character string
#' @param data A data frame used to compute residuals. If `NULL`, `model$data`
#'   is used.
#' @param ndraws An integer giving the number of posterior draws to use. If
#'   `NULL`, `brms::ndraws(model)` is used.
#' @param draw_ids An optional integer vector of posterior draw identifiers to
#'   include when summarizing across draws. If `NULL`, all available draws are
#'   used.
#' @param draw_ids_select An integer giving the posterior draw identifier used
#'   for single-draw QQ plots when `qq_type = "qq"`. If `NULL`, defaults to `1`.
#' @param resid_draws Pre computed residual draws
#' @param eta_sigma Pre computed sigma draws
#' @param summary A character string specifying the summary used when
#'   `qq_type = "qq_avg"`. Must be one of `"mean"`, `"median"`, or `"robust"`.
#' @param qq_type A character string specifying the output type. Must be one of
#'   `"data"`, `"qq"`, or `"qq_avg"`.
#' @param verbose A logical indicating whether extra progress information should
#'   be printed. Currently retained for future use.
#'   
#' @inheritParams brms::posterior_linpred
#' @inheritParams tidybayes::add_residual_draws
#'
#' @return
#' Returns one of the following:
#' \itemize{
#'   \item If `qq_type = "data"`, a list with components `data`,
#'     `sigma_parameters`, and `sigma_draws`.
#'   \item If `qq_type = "qq"`, a `ggplot2` object containing a QQ plot of
#'     Pearson residuals from a single posterior draw.
#'   \item If `qq_type = "qq_avg"`, a `ggplot2` object containing a QQ plot of
#'     summarized Pearson residuals across posterior draws.
#' }
#'
#' @keywords internal
#' @noRd
qq_plot_pearson <- function(model,
                            data = NULL,
                            ndraws = NULL, 
                            draw_ids = NULL,
                            draw_ids_select = NULL,
                            resid_draws = NULL,
                            eta_sigma = NULL,
                            summary = c("mean", "median", "robust"),
                            qq_type = 'qq',
                            resp = NULL,
                            transform = FALSE,
                            re_formula = NULL,
                            nlpar = NULL,
                            incl_thres = NULL,
                            sort = FALSE,
                            value = ".epred",
                            seed = 123,
                            category = ".category",
                            verbose = FALSE) {
  
  if(is.null(draw_ids_select)) draw_ids_select <- 1
  if(is.null(ndraws))          ndraws <- brms::ndraws(model)
  if(is.null(draw_ids))        draw_ids <- NULL
  if(is.null(data))            data <- model$data
  summary <- match.arg(summary)
  pearson_sigma_residuals <- function(model, 
                                      data, 
                                      ndraws, 
                                      draw_ids,
                                      draw_ids_select,
                                      resp,
                                      transform,
                                      re_formula,
                                      nlpar,
                                      incl_thres,
                                      sort,
                                      value,
                                      seed,
                                      category) {
    .draw <- NULL;
    .sigma <- NULL;
    draws_df <- posterior::as_draws_df(model)
    sigma_fixed_cols <- grep("^b_sigma", colnames(draws_df), value = TRUE)
    sigma_re_cols <- grep("(^r_.*sigma)|(^sd_.*sigma)|(^cor_.*sigma)", 
                          colnames(draws_df), value = TRUE)
    sigma_cols <- c(sigma_fixed_cols, sigma_re_cols)
    sigma_draws <- draws_df  %>% 
      dplyr::select(dplyr::any_of(c(".chain", ".iteration", ".draw", 
                                    sigma_cols)))
    
    if(is.null(resid_draws)) {
      resid_draws <- tidybayes::add_residual_draws(
        newdata = data,
        object = model,
        ndraws = ndraws,
        draw_ids = draw_ids,
        value = ".residual",
        seed = seed,
        category = category)
    }
    
    if(is.null(eta_sigma)) {
      eta_sigma <- brms::posterior_linpred(
        object = model,
        newdata = data,
        dpar = "sigma",
        ndraws = ndraws,
        draw_ids = draw_ids,
        resp = resp,
        transform = transform,
        re_formula = re_formula,
        nlpar = nlpar,
        incl_thres = incl_thres,
        sort = sort)
    }
    
    sigma_link <- model[['family']] [['link_sigma']]
    sigma_mat <- switch(
      sigma_link,
      log = exp(eta_sigma),
      identity = eta_sigma,
      softplus = log1p(exp(eta_sigma)),
      exp(eta_sigma)
    )
    
    sigma_long <- as.data.frame(sigma_mat) %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(.draw = dplyr::row_number()) %>% 
      tidyr::pivot_longer(
        cols = - .draw,
        names_to = ".row",
        values_to = ".sigma"
      ) %>% 
      dplyr::mutate(.row = as.numeric(sub(".*?([+-]?\\d+\\.?\\d*).*", "\\1",
                                          .row)))

    resid_col <- intersect(c(".residual", "residual", 
                             ".epred_residual", ".value"), 
                           colnames(resid_draws))[1]
    
    out <- resid_draws %>% 
      dplyr::left_join(sigma_long, by = c(".draw", ".row")) %>% 
      dplyr::mutate(
        .pearson = .data[[resid_col]] / .sigma
      )
    
    list(
      data = out,
      sigma_parameters = sigma_cols,
      sigma_draws = sigma_draws
    )
  }

  qq_pearson_sigma <- function(model, 
                               data, 
                               ndraws, 
                               draw_ids,
                               draw_ids_select,
                               resp,
                               transform,
                               re_formula,
                               nlpar,
                               incl_thres,
                               sort) {
    .draw <- NULL;
    .pearson <- NULL;
    res <- pearson_sigma_residuals(model = model, 
                                   data = data, 
                                   ndraws = ndraws,
                                   draw_ids = draw_ids,
                                   resp = resp,
                                   transform = transform,
                                   re_formula = re_formula,
                                   nlpar = nlpar,
                                   incl_thres = incl_thres,
                                   sort = sort,
                                   value = value,
                                   seed = seed,
                                   category = category)
    
    plot_data <- res$data %>%  dplyr::filter(.draw == draw_ids_select)
    
    ggplot2::ggplot(plot_data, ggplot2::aes(sample = .pearson)) +
      ggplot2::stat_qq(distribution = stats::qnorm) +
      ggplot2::stat_qq_line(distribution = stats::qnorm) +
      ggplot2::labs(
        x = "Theoretical Quantiles",
        y = "Pearson Residuals",
        title = "Q-Q Plot of Residuals",
      ) +
      ggplot2::theme_bw()
  }
  
  qq_pearson_avg_sigma <- function(model,
                                   data,
                                   ndraws,
                                   draw_ids,
                                   draw_ids_select,
                                   summary,
                                   resp,
                                   transform,
                                   re_formula,
                                   nlpar,
                                   incl_thres,
                                   sort) {
    .draw <- NULL;
    .pearson <- NULL;
    pearson_mean <- NULL;
    pearson_median <- NULL;
    pearson_plot <- NULL;
    
    res <- pearson_sigma_residuals(
      model = model,
      data = data,
      ndraws = ndraws,
      draw_ids = draw_ids,
      resp = resp,
      transform = transform,
      re_formula = re_formula,
      nlpar = nlpar,
      incl_thres = incl_thres,
      sort = sort,
      value = value,
      seed = seed,
      category = category
    )
    
    all_draws <- res$data %>% 
      dplyr::distinct(.draw) %>% 
      dplyr::pull(.draw)
    
    if (is.null(draw_ids)) {
      draw_ids <- all_draws
    } else {
      draw_ids <- intersect(draw_ids, all_draws)
    }
    
    if (length(draw_ids) == 0) {
      stop("No valid draw_ids found in res$data")
    }
    
    summary_data <- res$data %>% 
      dplyr::filter(.draw %in% draw_ids) %>% 
      dplyr::group_by(.row) %>% 
      dplyr::summarise(
        pearson_mean = mean(.pearson, na.rm = TRUE),
        pearson_median = stats::median(.pearson, na.rm = TRUE),
        pearson_sd = stats::sd(.pearson, na.rm = TRUE),
        pearson_mad = stats::mad(.pearson,
                                 center = stats::median(.pearson, na.rm = TRUE),
                                 na.rm = TRUE),
        n_draws = dplyr::n(),
        .groups = "drop"
      ) %>% 
      dplyr::mutate(
        pearson_plot = if (summary == "mean") {
          pearson_mean
        } else if (summary == "median") {
          pearson_median
        } else if (summary == "robust") {
          pearson_median
        } 
      )
    
    p <- ggplot2::ggplot(summary_data, ggplot2::aes(sample = pearson_plot)) +
      ggplot2::stat_qq(distribution = stats::qnorm) +
      ggplot2::stat_qq_line(distribution = stats::qnorm) +
      ggplot2::labs(
        title = "Q-Q Plot of Residuals",
        x = "Theoretical Quantiles",
        y = "Pearson Residuals"
      ) +
      ggplot2::theme_bw()
    
    list(
      plot = p,
      data = summary_data,
      draw_ids = draw_ids,
      summary = summary
    )
  }

  if(qq_type == 'data') {
    out <- pearson_sigma_residuals(
      model = model,
      data = data,
      ndraws = ndraws,
      draw_ids = draw_ids,
      draw_ids_select = draw_ids_select,
      resp = resp,
      transform = transform,
      re_formula = re_formula,
      nlpar = nlpar,
      incl_thres = incl_thres,
      sort = sort,
      value = value,
      seed = seed,
      category = category)
  } else if(qq_type == 'qq') {
    out <- qq_pearson_sigma(
      model = model,
      data = data,
      ndraws = ndraws,
      draw_ids = draw_ids,
      draw_ids_select = draw_ids_select,
      resp = resp,
      transform = transform,
      re_formula = re_formula,
      nlpar = nlpar,
      incl_thres = incl_thres,
      sort = sort)
  } else if(qq_type == 'qq_avg') {
    out_mean <- qq_pearson_avg_sigma(
      model = model,
      data = data,
      ndraws = ndraws,
      draw_ids = draw_ids,
      draw_ids_select = draw_ids_select,
      summary = summary,
      resp = resp,
      transform = transform,
      re_formula = re_formula,
      nlpar = nlpar,
      incl_thres = incl_thres,
      sort = sort)
    out <- out_mean$plot
  }
  return(out)
}









qq_brms_residuals <- function(
    model,
    residual_type = c("raw", "standardized", "pearson"),
    summary = c("mean", "median"),
    robust = FALSE,
    re_formula = NULL,
    ndraws = NULL,
    draw_ids = NULL,
    main_title = NULL,
    print = FALSE,
    ...) {
  theoretical <- NULL;
  residual_type <- match.arg(residual_type)
  summary <- match.arg(summary)
  .center_draws <- function(x, summary = "mean") {
    if (length(dim(x)) == 1L) {
      return(x)
    }
    if (summary == "mean") {
      base::colMeans(x)
    } else {
      apply(x, 2L, stats::median)
    }
  }
  .get_sigma_i <- function(fit, summary = "mean", robust = FALSE,
                           re_formula = NULL, ndraws = NULL, 
                           draw_ids = NULL, ...) {
    sigma_draws <- brms::posterior_epred(
      fit,
      dpar = "sigma",
      re_formula = re_formula,
      ndraws = ndraws,
      draw_ids = draw_ids,
      ...
    )
    sigma_i <- .center_draws(sigma_draws, summary = summary)
    if (robust) {
      sigma_i <- base::pmax(sigma_i, .Machine$double.eps)
    }
    sigma_i
  }
  # Get ordinary residuals (Y - Yrep)
  raw_res_draws <- stats::residuals(
    model,
    type = "ordinary",
    summary = FALSE,
    re_formula = re_formula,
    ndraws = ndraws,
    draw_ids = draw_ids,
    ...
  )
  raw_res <- .center_draws(raw_res_draws, summary = summary)
  if (residual_type == "raw") {
    res <- raw_res
    y_labs <- "Raw Residuals"
    
  } else if (residual_type == "pearson") {
    # Manually compute Pearson residuals: (Y - Yrep) / SD(Yrep)
    sigma_i <- .get_sigma_i(
      model,
      summary = summary,
      robust = robust,
      re_formula = re_formula,
      ndraws = ndraws,
      draw_ids = draw_ids,
      ...
    )
    if (length(sigma_i) != length(raw_res)) {
      stop("Length of sigma_i does not match number of observations.")
    }
    res <- raw_res / sigma_i
    y_labs <- "Pearson Residuals"
  } else if (residual_type == "standardized") {
    # Standardized residuals: raw residuals / sigma 
    # (observation-specific for distributional models)
    sigma_i <- .get_sigma_i(
      model,
      summary = summary,
      robust = robust,
      re_formula = re_formula,
      ndraws = ndraws,
      draw_ids = draw_ids,
      ...
    )
    if (length(sigma_i) != length(raw_res)) {
      stop("Length of sigma_i does not match number of observations.")
    }
    res <- raw_res / sigma_i
    y_labs <- "Standardized Residuals"
  }
  res <- res[base::is.finite(res)]
  qq_data <- base::data.frame(sample = base::sort(res))
  n <- base::nrow(qq_data)
  qq_data$theoretical <- stats::qnorm(stats::ppoints(n))
  q_sample <- stats::quantile(qq_data$sample, probs = c(0.25, 0.75), na.rm = TRUE)
  q_theory <- stats::qnorm(c(0.25, 0.75))
  slope <- base::diff(q_sample) / base::diff(q_theory)
  intercept <- q_sample[1L] - slope * q_theory[1L]
  p <- ggplot2::ggplot(qq_data, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_point(color = "steelblue", alpha = 0.7, size = 1.8) +
    ggplot2::geom_abline(
      intercept = intercept,
      slope = slope,
      linetype = "dashed",
      color = "red",
      linewidth = 0.7
    ) +
    ggplot2::labs(
      x = "Theoretical Quantiles",
      y = y_labs,
      title = if (!base::is.null(main_title)) main_title else base::paste("QQ Plot", "")
    ) +
    ggplot2::theme_minimal()
  if(print) base::print(p)
  return(p)
}






acf_residuals_custom <- function(
    model,
    residual_type = c("raw", "standardized"),
    summary = c("mean", "median"),
    re_formula = NULL,
    ndraws = NULL,
    draw_ids = NULL,
    idvar = NULL,
    xvar = NULL,
    sort_residuals = TRUE,
    lag_max = NULL,
    plot_type = c("ggplot", "base", "none"),
    ci_level = NULL,
    acf_cutoffs = NULL,
    acf_cutoffs_pm = FALSE,
    print = FALSE,
    main_title = "ACF of Residuals") {
  acfvar <- NULL;
  lagvar <- NULL;
  residual_type <- match.arg(residual_type)
  summary <- match.arg(summary)
  plot_type <- match.arg(plot_type)
  if (is.null(ci_level) && is.null(acf_cutoffs)) {
    acf_cutoffs <- 0.2
  }
  if (!is.null(ci_level)) {
    if (!is.numeric(ci_level) || length(ci_level) != 1L ||
        ci_level <= 0 || ci_level >= 1) {
      stop("ci_level must be a single number strictly between 0 and 1.")
    }
  }
  if (!is.null(acf_cutoffs)) {
    if (!is.numeric(acf_cutoffs)) {
      stop("acf_cutoffs must be numeric.")
    }
    if (length(acf_cutoffs) == 1L) {
      if(acf_cutoffs_pm) acf_cutoffs <- c(-abs(acf_cutoffs), abs(acf_cutoffs))
    }
  }
  dat <- model$data
  if (!is.null(idvar)) {
    if (!is.character(idvar) || length(idvar) != 1L || !idvar %in% names(dat)) {
      stop("idvar must be a valid column name in model$data.")
    }
  }
  if (!is.null(xvar)) {
    if (!is.character(xvar) || length(xvar) != 1L || !xvar %in% names(dat)) {
      stop("xvar must be a valid column name in model$data.")
    }
  }
  .center_draws <- function(x, summary = "mean") {
    if (length(dim(x)) == 1L) return(x)
    if (summary == "mean") colMeans(x) else apply(x, 2L, median)
  }
  .get_sigma_i <- function(fit, summary = "mean",
                           re_formula = NULL, ndraws = NULL, draw_ids = NULL) {
    sigma_draws <- brms::posterior_epred(
      fit,
      dpar = "sigma",
      re_formula = re_formula,
      ndraws = ndraws,
      draw_ids = draw_ids,
      sort = FALSE
    )
    .center_draws(sigma_draws, summary = summary)
  }
  raw_res_draws <- stats::residuals(
    model,
    type = "ordinary",
    summary = FALSE,
    re_formula = re_formula,
    ndraws = ndraws,
    draw_ids = draw_ids,
    sort = FALSE
  )
  raw_res <- .center_draws(raw_res_draws, summary = summary)
  if (residual_type == "raw") {
    res <- raw_res
  } else {
    sigma_i <- .get_sigma_i(
      model,
      summary = summary,
      re_formula = re_formula,
      ndraws = ndraws,
      draw_ids = draw_ids
    )
    if (length(sigma_i) != length(raw_res)) {
      stop("Length of sigma_i does not match number of observations.")
    }
    res <- raw_res / pmax(sigma_i, .Machine$double.eps)
  }
  keep <- is.finite(res)
  res <- res[keep]
  dat <- dat[keep, , drop = FALSE]
  if (length(res) < 2L) 
    stop("Need at least 2 finite residuals to compute autocorrelation.")
  if (is.null(lag_max)) {
    if (!is.null(idvar)) {
      lag_max <- max(table(dat[[idvar]])) - 1L
    } else {
      lag_max <- length(res) - 1L
    }
  }
  lag_max <- max(0L, min(as.integer(lag_max), length(res) - 1L))
  lag_max <- pmin(brms::ndraws(model), lag_max)
  
  if (sort_residuals) {
    if (!is.null(idvar) && !is.null(xvar)) {
      # ord <- order(dat[[idvar]], dat[[xvar]])
      ord <- order(dat[[xvar]])
      res <- res[ord]
      dat <- dat[ord, , drop = FALSE]
    } else if (!is.null(xvar)) {
      ord <- order(dat[[xvar]])
      res <- res[ord]
      dat <- dat[ord, , drop = FALSE]
    } else if (!is.null(idvar)) {
      ord <- order(dat[[idvar]])
      res <- res[ord]
      dat <- dat[ord, , drop = FALSE]
    }
  }
  acf_obj <- stats::acf(
    res,
    lag.max = lag_max,
    plot = FALSE,
    na.action = stats::na.pass
  )
  acf_df <- data.frame(
    lagvar = as.numeric(acf_obj$lag),
    acfvar = as.numeric(acf_obj$acf)
  )
  acf_df <- acf_df[acf_df$lagvar >= 0, , drop = FALSE]
  ci <- NULL
  if (!is.null(ci_level)) {
    z_crit <- stats::qnorm((1 + ci_level) / 2)
    ci <- z_crit / sqrt(length(res))
  }
  if (plot_type == "base") {
    stats::acf(
      res,
      lag.max = lag_max,
      plot = TRUE,
      ci = if (is.null(ci_level)) 0 else ci_level,
      main = main_title,
      na.action = stats::na.pass
    )
    if (!is.null(acf_cutoffs)) {
      for (cc in acf_cutoffs) graphics::abline(h = cc, col = "red", lty = 2)
    }
    
    if (!is.null(ci)) graphics::abline(h = c(-ci, ci), col = "blue", lty = 3)
  }
  if (plot_type == "ggplot") {
    p <- ggplot2::ggplot(acf_df, ggplot2::aes(x = lagvar, y = acfvar)) +
      ggplot2::geom_col(fill = "steelblue", width = 0.8) +
      ggplot2::geom_hline(yintercept = 0, color = "black") +
      ggplot2::labs(x = "Lag", y = "ACF", title = main_title) +
      ggplot2::theme_minimal()
    if (!is.null(acf_cutoffs)) {
      for (cc in acf_cutoffs) {
        p <- p + ggplot2::geom_hline(yintercept = cc, 
                                     linetype = "dashed", color = "red")
      }
    }
    if (!is.null(ci)) {
      p <- p + ggplot2::geom_hline(yintercept = c(-ci, ci), 
                                   linetype = "dotted", color = "blue")
    }
    if(print) print(p)
  }
  return(p)
}









