
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}


if(set_skip_run_ci) {
  skip_run_ci()
}


###############################################################################
# Test bsitar with rcs settings
###############################################################################

test_that("bsitar works fully with rcs settings", {
  skip_on_cran()
  
  testthat::skip_if_not(!test_univariate_fit_cov)
  
  test_fit <- univariate_fit_cov
  
  ##############################################################################
  # expose_model_functions()
  ##############################################################################
  
  outx <- expose_model_functions(test_fit, expose = FALSE, returnobj = FALSE)
  expect_true(is.null(outx))
  
  outx <- expose_model_functions(test_fit, expose = FALSE, returnobj = TRUE)
  expect_true(inherits(outx, 'bgmfit'))
  
  test_fit <- expose_model_functions(test_fit, expose = TRUE, returnobj = TRUE)
  expect_true(inherits(test_fit, 'bgmfit'))
  
  
  ##############################################################################
  # growthparameters()
  ##############################################################################
  
  test_gparms_gp <- growthparameters(test_fit, re_formula = NA,
                                     draw_ids = draw_ids)
  check_val <- 9.12
  compa_val <- mean(test_gparms_gp[["Estimate"]])
  expect_equal(compa_val, check_val,tolerance = 0.01)
  
  
  ##############################################################################
  # get_growthparameters() vs growthparameters()
  ##############################################################################
  
  test_gparms_re <- get_growthparameters(test_fit, by = 'id', 
                                         re_formula = NULL, 
                                         draw_ids = draw_ids)
  test_gparms_gp_re <- growthparameters(test_fit, 
                                        re_formula = NULL, 
                                        draw_ids = draw_ids)
  
  expect_equal(mean(test_gparms_re[["Estimate"]], na.rm = T), 
               mean(test_gparms_gp_re[["Estimate"]], na.rm = T),
               tolerance = 0.01)
  
  check_val <- 9.14
  compa_val <- mean(test_gparms_re[["Estimate"]], na.rm = T)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  
  ##############################################################################
  # modelbased_growthparameters()
  ##############################################################################
  
  mgc_R <- modelbased_growthparameters(test_fit, 
                                       method = 'custom', 
                                       call_function = "R",
                                       draw_ids = draw_ids, 
                                       verbose = FALSE)
  
  mgc_Stan <- modelbased_growthparameters(test_fit, 
                                          method = 'custom', 
                                          call_function = "Stan",
                                          draw_ids = draw_ids,
                                          verbose = FALSE)
  
  expect_equal(mean(mgc_R$estimate), mean(mgc_Stan$estimate), tolerance = 0.01)
  
  check_val <- 11.49426
  compa_val <- mean(mgc_R$estimate)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  ##############################################################################
  # modelbased_growthparameters_call()
  ##############################################################################
  
  mgc <- modelbased_growthparameters_call(test_fit, 
                                          method = 'custom', 
                                          draw_ids = draw_ids)
  check_val <- 148.4
  compa_val <- mean(mgc$distance$Estimate)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  check_val <- 6.14
  compa_val <- mean(mgc$velocity$Estimate)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  ##############################################################################
  # optimize_model
  ##############################################################################
  
  # suppressWarnings(suppressMessages({
    omc <- optimize_model(test_fit,
                          newdata = NULL,
                          optimize_df = 3,
                          optimize_x = list(NULL, log, sqrt),
                          optimize_y = list(NULL, log, sqrt),
                          transform_prior_class = c('beta', 'sd', 
                                                    'rsd', 'sigma', 'dpar'),
                          transform_beta_coef = c('b', 'c'),
                          transform_sd_coef = c('b', 'c'),
                          exclude_default = TRUE,
                          add_fit_criteria = c("loo", "waic", "bayes_R2"),
                          byresp = FALSE,
                          model_name = NULL,
                          overwrite = FALSE,
                          file = NULL,
                          force_save = FALSE,
                          save_each = FALSE,
                          digits = 2,
                          cores = 1,
                          verbose = FALSE,
                          expose_function = FALSE,
                          usesavedfuns = FALSE,
                          clearenvfuns = NULL,
                          envir = NULL)
  # }))
  
  expect_true(is.data.frame(omc$optimize_summary))
  
  omcdf <- omc$optimize_summary %>% tidyr::drop_na('Estimate')
  omcdf_looic <- omcdf %>% dplyr::filter(Parameter == 'looic') %>% 
    dplyr::pull(Estimate) %>% mean()
  omcdf_waic <- omcdf %>% dplyr::filter(Parameter == 'waic') %>% 
    dplyr::pull(Estimate) %>% mean()
  omcdf_bayes_R2 <- omcdf %>% dplyr::filter(Parameter == 'bayes_R2') %>% 
    dplyr::pull(Estimate) %>% mean()
  
  expect_equal(omcdf_looic,   -9293.635, tolerance = 0.01)
  expect_equal(omcdf_waic,    -9293.145, tolerance = 0.01)
  expect_equal(omcdf_bayes_R2,  0.9625, tolerance = 0.01)
  
  
  
  
  ##############################################################################
  # loo_validation()
  ##############################################################################
  
  outx <- loo_validation(omc$models[[1]], omc$models[[2]], draw_ids = draw_ids)
  
  expect_true(is.list(outx))
  
  
  
  ##############################################################################
  # compare_models()
  ##############################################################################
  
  outx <- compare_models(omc$models)
  
  expect_true(is.data.frame(outx))
  
  
  
  ##############################################################################
  # hypothesis_test
  ##############################################################################
  
  set_hypothesis <- c("a_Intercept = 0", "b_Intercept > 1")
  hyp_ex1 <- brms::hypothesis(test_fit, hypothesis = set_hypothesis,
                              draw_ids = draw_ids, seed = 123)
  
  hyp_ex11 <- hypothesis_test(test_fit,
                              draw_ids = draw_ids,
                              hypothesis_str = set_hypothesis,
                              seed = 123)
  
  expect_true(all.equal(hyp_ex11, hyp_ex1))
  
  
  set_parameters <- c("b_a_Intercept", "b_b_Intercept")
  set_range <- list(b_a_Intercept = c(100, 150), b_b_Intercept = c(-2, 2))
  
  
  suppressMessages({
    hyp_ex1 <- bayestestR::equivalence_test(test_fit,
                                            parameters = set_parameters,
                                            range = set_range,
                                            draw_ids = draw_ids)
  })
  
  expect_equal(TRUE, inherits(hyp_ex1, "equivalence_test"))
  
  check_val <- 0.237
  compa_val <- mean(hyp_ex1$ROPE_Percentage)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  hyp_ex11 <- hypothesis_test(test_fit,
                              draw_ids = draw_ids,
                              parameters = set_parameters,
                              range = set_range,
                              plot = FALSE)
  
  check_val <- 65.5
  compa_val <- mean(hyp_ex11$HDI_low)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  set_parameters <- c("b_a_Intercept", "b_b_Intercept")
  set_null <- list(b_a_Intercept = 0, b_b_Intercept = 0)
  hyp_ex1 <- bayestestR::p_direction(test_fit,
                                     parameters = set_parameters,
                                     null = set_null,
                                     draw_ids = draw_ids)
  
  check_val <- 0.996
  compa_val <- mean(hyp_ex1$pd)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  ##############################################################################
  # fitted_draws()
  ##############################################################################
  
  
  outx <- fitted_draws(test_fit, draw_ids = draw_ids, deriv = 0, 
                       model_deriv = TRUE)
  check_val <- 157.665
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  outx <- fitted_draws(test_fit,  draw_ids = draw_ids, deriv = 0, 
                       model_deriv = FALSE)
  check_val <- 157.665
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  outx <- fitted_draws(test_fit, draw_ids = draw_ids, deriv = 1,
                       model_deriv = TRUE)
  check_val <- 4.23
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  outx <- fitted_draws(test_fit,  draw_ids = draw_ids, deriv = 1, 
                       model_deriv = FALSE)
  check_val <- 4.23
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  outx <- fitted_draws(test_fit,  draw_ids = draw_ids, deriv = 1, 
                       model_deriv = FALSE)
  check_val <- 4.23
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  
  
  
  outx <- fitted_draws(test_fit, newdata_fixed  = FALSE ,
                       draw_ids = draw_ids, deriv = 1, model_deriv = FALSE)
  check_val <- 0.340
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  outx <- fitted_draws(test_fit, newdata_fixed  = 0 ,
                       draw_ids = draw_ids, deriv = 1, model_deriv = FALSE)
  check_val <- 0.340
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  outx <- fitted_draws(test_fit, newdata_fixed  = TRUE ,
                       draw_ids = draw_ids, deriv = 1, model_deriv = FALSE)
  check_val <- 0.340
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  outx <- fitted_draws(test_fit, newdata_fixed  = 1 ,
                       draw_ids = draw_ids, deriv = 1, model_deriv = FALSE)
  check_val <- 0.340
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  
  outx <- fitted_draws(test_fit, newdata_fixed  = 2 ,
                       draw_ids = draw_ids, deriv = 1, model_deriv = FALSE)
  check_val <- 4.23
  compa_val <- mean(outx[,1])
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  
  
  
  
  ##############################################################################
  # predict_draws()
  ##############################################################################
  
  outx_predict_draws <- predict_draws(test_fit, draw_ids = draw_ids, deriv = 0,
                                      model_deriv = TRUE)
  check_val_predict_draws <- 157.665
  compa_val_predict_draws <- mean(outx_predict_draws[,1])
  expect_equal(compa_val_predict_draws, check_val_predict_draws, tolerance = 0.01)
  
  
  outx_predict_draws <- predict_draws(test_fit,  draw_ids = draw_ids, deriv = 0,
                                      model_deriv = FALSE)
  check_val_predict_draws <- 157.665
  compa_val_predict_draws <- mean(outx_predict_draws[,1])
  expect_equal(compa_val_predict_draws, check_val_predict_draws, tolerance = 0.01)
  
  
  outx_predict_draws <- predict_draws(test_fit, draw_ids = draw_ids, deriv = 1,
                                      model_deriv = TRUE)
  check_val_predict_draws <- 4.23
  compa_val_predict_draws <- mean(outx_predict_draws[,1])
  expect_equal(compa_val_predict_draws, check_val_predict_draws, tolerance = 0.01)
  
  
  outx_predict_draws <- predict_draws(test_fit,  draw_ids = draw_ids, deriv = 1,
                                      model_deriv = FALSE)
  check_val_predict_draws <- 4.22
  compa_val_predict_draws <- mean(outx_predict_draws[,1])
  expect_equal(compa_val_predict_draws, check_val_predict_draws, tolerance = 0.01)
  
  outx_predict_draws <- predict_draws(test_fit,  draw_ids = draw_ids, deriv = 1,
                                      model_deriv = FALSE)
  check_val_predict_draws <- 4.24
  compa_val_predict_draws <- mean(outx_predict_draws[,1])
  expect_equal(compa_val_predict_draws, check_val_predict_draws, tolerance = 0.01)
  
  
  ##############################################################################
  # plot_ppc()
  ##############################################################################
  
  outx <- plot_ppc(test_fit, draw_ids = draw_ids)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  
  ##############################################################################
  # plot_conditional_effects()
  ##############################################################################
  
  outx <- plot_conditional_effects(test_fit, re_formula = NA,
                                   draw_ids = draw_ids, deriv = 0, 
                                   # effects = c('sex', 'age'),
                                   model_deriv = F)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  check_val <- 160.504
  compa_val <- mean(outx@data$estimate__)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  outx <- plot_conditional_effects(test_fit, re_formula = NA,
                                   draw_ids = draw_ids, deriv = 0, 
                                   model_deriv = T)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  check_val <- 160.504
  compa_val <- mean(outx@data$estimate__)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  
  
  outx <- plot_conditional_effects(test_fit, re_formula = NA,
                                   draw_ids = draw_ids, deriv = 1, 
                                   model_deriv = F)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  check_val <- 4.24
  compa_val <- mean(outx@data$estimate__)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  outx <- plot_conditional_effects(test_fit, re_formula = NA,
                                   draw_ids = draw_ids, deriv = 1, 
                                   model_deriv = T)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  check_val <- 4.24
  compa_val <- mean(outx@data$estimate__)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  #
  outx <- plot_conditional_effects(test_fit, re_formula = NULL,
                                   draw_ids = draw_ids, deriv = 0, 
                                   # effects = c('sex', 'age'),
                                   model_deriv = F)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  check_val <- 160.504
  compa_val <- mean(outx@data$estimate__)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  outx <- plot_conditional_effects(test_fit, re_formula = NULL,
                                   draw_ids = draw_ids, deriv = 0, 
                                   model_deriv = T)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  check_val <- 160.504
  compa_val <- mean(outx@data$estimate__)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  
  
  outx <- plot_conditional_effects(test_fit, re_formula = NULL,
                                   draw_ids = draw_ids, deriv = 1, 
                                   model_deriv = F)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  check_val <- 4.24
  compa_val <- mean(outx@data$estimate__)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  outx <- plot_conditional_effects(test_fit, re_formula = NULL,
                                   draw_ids = draw_ids, deriv = 1, 
                                   model_deriv = T)
  
  expect_true(inherits(outx, 'ggplot'))
  
  
  check_val <- 4.24
  compa_val <- mean(outx@data$estimate__)
  expect_equal(compa_val, check_val, tolerance = 0.01)
  
  
  
  ##############################################################################
  # update_model()
  ##############################################################################
  
  # devtools::load_all()
  
  
  update_model_args <- list()
  update_model_args[['model']] <- test_fit
  update_model_args[['get_stancode']] <- TRUE
  update_model_args[['get_standata']] <- FALSE
  update_model_args[['iter']] <- 10
  
  update_model_args[['refresh']] <- 0
  update_model_args[['silent']]  <- 2
  
  
  axx <- do.call(update_model, update_model_args)
  expect_true(inherits(axx, 'character'))
  
  update_model_args[['get_stancode']] <- FALSE
  update_model_args[['get_standata']] <- TRUE
  
  axx <- do.call(update_model, update_model_args)
  expect_true(inherits(axx, 'list'))
  
  
  update_model_args[['get_stancode']] <- FALSE
  update_model_args[['get_standata']] <- FALSE
  update_model_args[['recompile']] <- TRUE
  
  axx <- do.call(update_model, update_model_args)
  expect_true(inherits(axx, 'bgmfit'))
  
  
  update_formual_args <- list(a_formula = ~ 1,
                              b_formula = ~ 1,
                              c_formula = ~ 1,
                              d_formula = ~ 1,
                              s_formula = ~ 1,
                              a_formula_gr = ~ 1,
                              b_formula_gr = ~ 1,
                              c_formula_gr = ~ 1,
                              d_formula_gr = ~ 1,
                              a_formula_gr_str =  ~ (1|i|id),
                              b_formula_gr_str =  ~ (1|i|id),
                              c_formula_gr_str =  ~ (1|i|id),
                              d_formula_gr_str =  ~ (1|i|id),
                              d_adjusted = FALSE,
                              sigma_formula = ~ 1,
                              sigma_formula_gr = NULL,
                              sigma_formula_gr_str = NULL,
                              sigma_formula_manual = "list(nlf(sigma ~ z) + 
                              lf(z ~ 1 + (1 | gr(id))))",
                              dpar_formula = NULL,
                              autocor_formula =  ~unstr(time, id))
  
  
  
  update_model_args[['get_stancode']] <- FALSE
  update_model_args[['get_standata']] <- FALSE
  update_model_args[['recompile']] <- FALSE
  
  update_model_args[['new_threads']] <- NA
  
  update_model_args[['sigma_formula_manual']] <- 
    gsub_space(update_model_args[['sigma_formula_manual']] )
  
  for (check_formulasi in names(update_formual_args)) {
    update_model_args[[check_formulasi]] <- update_formual_args[[check_formulasi]]
    # axx <- do.call(update_model, update_model_args)
    # expect_true(inherits(axx, 'bgmfit'))
  }
  
  axx <- do.call(update_model, update_model_args)
  expect_true(inherits(axx, 'bgmfit'))
  
  
  ##############################################################################
  # 
  ##############################################################################

  
})





test_that("bsitar works fully with rcs settings - summary_table", {
  skip_on_cran()
  
  testthat::skip_if_not(!test_univariate_fit_cov)
  
  test_fit <- univariate_fit_cov
  
  ##############################################################################
  # summary_table
  ##############################################################################
  model_list <- list()
  model_list[['a']] <- test_fit
  model_list[['b']] <- test_fit
  
  summary_table_mf <- summary_table(model_list,
                                    effects = "all",
                                    gof_map = c('looic', 'waic', 'nobs'),
                                    # gof_map = 'all',
                                    estimate = "{estimate} [{conf.low}, {conf.high}]",
                                    # add_columns = 'ESS',
                                    output = 'dataframe',
                                    diagnostic = 'all',
                                    add_diagnostic = F,
                                    add_pd = F,
                                    add_section = T,
                                    add_component = F,
                                    sort_gr_by = T,
                                    gsub_fixed = " x ",
                                    gsub_random = list(it = c("::", ":"),
                                                       by = c(" = ", " x ")),
                                    statistic = NULL,
                                    sigma_exp = F,
                                    fmt = 2,
                                    centrality = 'mean')
  
  
  rm_these <- c('section.y', 'part', 'statistic')
  
  summary_table_mf <- summary_table_mf %>%
    dplyr::select(-dplyr::any_of(rm_these)) %>%
    dplyr::relocate(dplyr::starts_with('section')) %>%
    # dplyr::select(dplyr::any_of(dplyr::starts_with('section.'))) %>%
    dplyr::rename_with(
      .fn = ~ sub("\\.x$", "", .x),  # removes ".x" at the end
      .cols = dplyr::any_of(c(dplyr::starts_with('section.')))
    )
  
  # writexl::write_xlsx(summary_table_mf, "summary_table_mf.xlsx")
  # flextable::as_flextable(summary_table_mf, max_row = nrow(summary_table_mf))
  
  expect_true(is.data.frame(summary_table_mf))
  
  
  
})





test_that("bsitar works fully with rcs settings - prior_summary_table", {
  skip_on_cran()
  
  testthat::skip_if_not(!test_univariate_fit_cov)
  
  test_fit <- univariate_fit_cov
  
  ##############################################################################
  # prior_summary_table
  ##############################################################################
  
  prior_summary_table <- prior_summary_table(test_fit,
                                          set_width = c(0.95, 0.9999),
                                          set_digits = 1,
                                          empty = "-",
                                          print = FALSE,
                                          return_table = TRUE,
                                          return_file = NULL,
                                          flex_table = TRUE,
                                          path = NULL,
                                          title = NULL,
                                          align = "center",
                                          sheet_name = "table",
                                          draw_samples = 100,
                                          add_range = FALSE,
                                          transform_class = NULL,
                                          transform_parameter = NULL,
                                          transform_fun = NULL,
                                          range_method_arg = NULL,
                                          seed = 123,
                                          verbose = FALSE)
  
  
  
  expect_true(inherits(prior_summary_table, "flextable"))
  
  
  
})
