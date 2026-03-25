

#' Title
#'
#' @param model 
#' @param newdata 
#' @param ndraws 
#' @param draw_ids 
#' @param peak 
#' @param takeoff 
#' @param trough 
#' @param acgv 
#' @param verbose 
#' @param newdata_fixed 
#' @param ... 
#'
#' @returns A list
#' @keywords internal
#' @noRd
#'
modelbased_growthparameters_nonS3 <- function(model, 
                                        newdata = NULL,
                                        ndraws = NULL,
                                        draw_ids = NULL,
                                        peak = TRUE,
                                        takeoff = FALSE,
                                        trough = FALSE,
                                        acgv = FALSE,
                                        verbose = NULL, 
                                        newdata_fixed = FALSE,
                                        ...) {
  
  if(is.null(newdata)) {
    newdata <- model$model_info$bgmfit.data
  }
  
  if(!is.null(draw_ids)) {
    eval_draw_ids <- eval(draw_ids)
  }
  
  # growthparameters_args <- list()
  # growthparameters_args[['model']] <- model
  # growthparameters_args[['newdata']] <- newdata
  # growthparameters_args[['ndraws']] <- ndraws
  # growthparameters_args[['draw_ids']] <- draw_ids
  # growthparameters_args[['peak']] <- peak
  # growthparameters_args[['takeoff']] <- takeoff
  # growthparameters_args[['trough']] <- trough
  # growthparameters_args[['acgv']] <- acgv
  # growthparameters_args[['get_dv']] <- TRUE
  # growthparameters_args[['summary']] <- FALSE
  
  
  # azz2 <- CustomDoCall(growthparameters, growthparameters_args)

  
  fitted_draws_args <- list()
  fitted_draws_args[['model']] <- model
  fitted_draws_args[['newdata']] <- newdata
  fitted_draws_args[['ndraws']] <- ndraws
  fitted_draws_args[['draw_ids']] <- draw_ids
  fitted_draws_args[['re_formula']] <- NULL
  fitted_draws_args[['summary']] <- FALSE
  
  fitted_draws_args_d0 <- fitted_draws_args
  fitted_draws_args_d1 <- fitted_draws_args
  fitted_draws_args_d0[['deriv']] <- 0
  fitted_draws_args_d1[['deriv']] <- 1
  
  
  get_growthparameters_args <- list()
  get_growthparameters_args[['model']] <- model
  get_growthparameters_args[['newdata']] <- newdata
  get_growthparameters_args[['ndraws']] <- ndraws
  get_growthparameters_args[['draw_ids']] <- draw_ids
  get_growthparameters_args[['pdrawsp']] <- 'return'
  get_growthparameters_args[['parameter']] <- 'apgv'
  get_growthparameters_args[['itransform']] <- ''
  
  get_predictions_args <- get_growthparameters_args
  get_predictions_args[['re_formula']] <- NULL
  get_predictions_args[['peak']] <- NULL
  get_predictions_args[['takeoff']] <- NULL
  get_predictions_args[['trough']] <- NULL
  get_predictions_args[['acgv']] <- NULL
  get_predictions_args[['newdata_fixed']] <- TRUE
  
  get_predictions_args_d0 <- get_predictions_args
  get_predictions_args_d1 <- get_predictions_args
  get_predictions_args_d0[['deriv']] <- 0
  get_predictions_args_d1[['deriv']] <- 1
 
  
  ApvX0 <- CustomDoCall(get_growthparameters, 
                        get_growthparameters_args)
  

  xyadj_curves_args <- list()
  xyadj_curves_args[['model']]  <- model
  # xyadj_curves_args[['newdata']]  <- newdata
  xyadj_curves_args[['ndraws']] <- ndraws
  xyadj_curves_args[['tomean']] <- FALSE
  xyadj_curves_args[['get_dv']] <- TRUE
  
  
  ApvX0[['drawid']] <- eval_draw_ids
  
  ApvX0_for_list_c <- tibble::as_tibble(ApvX0) %>% 
    fastplyr::f_select(-'parameter') %>% 
    fastplyr::f_rename('age' = 'estimate') %>% 
    dplyr::mutate(tempjoinby = 'tempjoinby') %>% 
    dplyr::relocate('age') 
  
  # %>% 
  #   dplyr::mutate(loopcounter = dplyr::row_number()) 
    
  newdata_for_list_c <-tibble::as_tibble(newdata) %>%
    fastplyr::f_select(-'age') %>% 
    dplyr::mutate(tempjoinby = 'tempjoinby')
  
  
  

  
  newdata_list_c <- list()
  for (i in 1:nrow(ApvX0)) {
    newdata_list_c[[i]] <- fastplyr::f_left_join(newdata_for_list_c, 
                                                 ApvX0_for_list_c[i,], 
                                                 by = "tempjoinby") %>% 
      fastplyr::f_select(-'tempjoinby')
  }
  
  # my_array  <- array(NA, dim = c(2, nrow(newdata), 3))
  
  my_matrix <- matrix(NA, nrow = 3, ncol = nrow(newdata))
  analyze_data <- function(data, ...) {
    xyadj_curves_args[['newdata']] <- data
    xyadj_curves_args[['draw_ids']] <- data[['drawid']] [1]
    xyadj_dv <- CustomDoCall(xyadj_curves, xyadj_curves_args)
    # i <-  data[['loopcounter']] [1]
    fitted_draws_args_d0[['draw_ids']] <- data[['drawid']] [1]
    fitted_draws_args_d1[['draw_ids']] <- data[['drawid']] [1]
    fitted_draws_args_d0[['newdata']] [['age']]  <- xyadj_dv %>% t()
    fitted_draws_args_d1[['newdata']] [['age']]  <- xyadj_dv %>% t()
    
    get_predictions_args_d0[['draw_ids']] <- data[['drawid']] [1]
    get_predictions_args_d1[['draw_ids']] <- data[['drawid']] [1]
    get_predictions_args_d0[['newdata']] [['age']]  <- xyadj_dv %>% t()
    get_predictions_args_d1[['newdata']] [['age']]  <- xyadj_dv %>% t()

    yyadj_dv <- CustomDoCall(get_predictions, get_predictions_args_d0)
    vyadj_dv <- CustomDoCall(get_predictions, get_predictions_args_d1)
    
    yyadj_dv <- yyadj_dv[["estimate"]]
    vyadj_dv <- vyadj_dv[["estimate"]]
    
    # yyadj_dv <- CustomDoCall(fitted_draws, fitted_draws_args_d0)
    # vyadj_dv <- CustomDoCall(fitted_draws, fitted_draws_args_d1)
    # my_array <- array(c(xyadj_dv, yyadj_dv, vyadj_dv),
    #                   c(1, ncol(xyadj_dv), 3))
    # my_array
    
    my_matrix[1, ] <- xyadj_dv
    my_matrix[2, ] <- yyadj_dv
    my_matrix[3, ] <- vyadj_dv
    my_matrix
  }
  
  # results_list <- sapply(newdata_list_c, analyze_data, simplify = F)
  
  
  # results_list <- lapply(newdata_list_c, analyze_data)
  
  results_list <- future.apply::future_lapply(newdata_list_c, analyze_data)
  xyvyadj_rows <- CustomDoCall(cbind, results_list)
  xyadj_rows <- matrix(xyvyadj_rows[1,], nrow = length(results_list), byrow = T)
  yyadj_rows <- matrix(xyvyadj_rows[2,], nrow = length(results_list), byrow = T)
  vyadj_rows <- matrix(xyvyadj_rows[3,], nrow = length(results_list), byrow = T)
  
  xyadj_rows <- apply(xyadj_rows, 2, model$model_info$ixfuntransform2)
  xyadj_summary <- brms::posterior_summary(xyadj_rows)
  yyadj_summary <- brms::posterior_summary(yyadj_rows)
  vyadj_summary <- brms::posterior_summary(vyadj_rows)

  out <- list()
  out[['xyadj']] <- xyadj_summary
  out[['yyadj']] <- yyadj_summary
  out[['vyadj']] <- vyadj_summary
  
  return(out)
}

