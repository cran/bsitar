
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}


###############################################################################
# Test plot_diagnostics with rcs settings
###############################################################################

test_that("test plot_diagnostics", {
  skip_on_cran()
 
  suppressWarnings(suppressMessages({
    test_fit <- bsitar(x = age, y = height, id = id,
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
                       get_stancode = FALSE,
                       get_standata = FALSE,
                       chains = 1, cores = 1, iter = 10,
                       backend = "rstan",
                       sample_prior = "no",
                       threads = threading(NULL),
                       # init = '0',
                       init = NULL, # Don't use default random with init_r = 0.5
                       vcov_init_0 = TRUE,
                       refresh = 0, silent = 2,
                       seed = 123)
  }))
  
  # "pairs", - need ggpubr 
  # set wrap_title = FALSE, otherwise need ggtext
  
  plots_c <- c("rvf", "rvp", "qq", "qqn", "qqp",  "acf", "acfp",
               "acfr", "trace", "dens_overlay", "rhat", "rhat_hist", 
               "neff", "ppc_overlay", "ppc_hist", "ppc_scatter", "ppc_stat")
  
  plots_list <- list()
  plots_est_list <- list()
  for (plots_ci in plots_c) {
    plots_list[[plots_ci]] <- 
      plot_diagnostics(test_fit, plots = plots_ci, 
                       variable = c("b_a_Intercept", "b_b_Intercept"),
                       wrap_title = FALSE, combine = FALSE, print = FALSE)
  }
  
  
  for (plots_ci in plots_c) {
    if(length(plots_list[[plots_ci]]) == 0) {
      # 
    } else {
      TFOB <- inherits(plots_list[[plots_ci]], 'ggplot') 
      expect_true(TFOB)
    }
  }
  
  
})



