
.onAttach <- function(libname, pkgname) {
  
  invisible(suppressPackageStartupMessages(
    sapply( c("brms", "future"), requireNamespace, quietly = TRUE)
    ))
  
  if("rethinking" %in% (.packages())){
    packageStartupMessage("Package 'rethinking' detached and unloaded as it",
            " \ncreates conflict with the current rstan version ", 
            utils::packageVersion('rstan'))
    detach("package:rethinking", unload=TRUE) 
  }
  
}


# utils-helper-1-misc-1
# utils-helper-2-misc-2
# utils-helper-3-prepare_data2
# utils-helper-4-prepare_formula
# utils-helper-5-prepare_formula_sigma
# utils-helper-6-export 'nsp' 'nsk' 'rcs' etc string for stan
# utils-helper-7-prepare_function_nsp_rcs
# utils-helper-8-prepare_function_sigma - not used now. instead use for bsp msp isp R_str
# utils-helper-9-set_priors_initials
# utils-helper-10-prepare_priors
# utils-helper-11-prepare_initials
# utils-helper-12-R_support_nsp
# utils-helper-13-get_hypothesis_x
# utils-helper-14-get.newdata
# utils-helper-15-get_idata
# utils-helper-16-plot_curves
# utils-helper-17-stancode_custom
# utils-helper-18-standata_custom
# utils-helper-19-prepare_transformations
# utils-helper-20-Pipe
# utils-helper-21-modelbased_growthparameters_nonS3 - not using
# utils-helper-22-modelbased_growthparameters_call - not using
# utils-helper-23-modelbased_growthparameters_helpers
# utils-helper-24-modelbased_growthparameters - via Stan
# utils-helper-25-fast_nsk_rcsfun_str_get
# utils-helper-26-generate_age_list
# utils-helper-27-knot_search



