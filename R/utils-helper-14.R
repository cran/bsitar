

#' An internal function to create data frame for post-processing
#' 
#' @param idata_method A character string to indicate the interpolation method.
#'   Options are \emph{method 1} (specified as  \code{'m1'}, default) and
#'   \emph{method 2} (specified as \code{'m2'}). The \code{'m1'} calls an
#'   internal function \code{idatafunction} whereas \code{'m2'} calls the
#'   \code{get_idata} function for data interpolation. The \emph{method 1}
#'   (\code{'m1'}) is adapted from the the \pkg{iapvbs} and is documented here
#'   <https://rdrr.io/github/Zhiqiangcao/iapvbs/src/R/exdata.R>. The
#'   \emph{method 2} (\code{'m2'}) is adapted from the the \pkg{JMbayes} and is
#'   documented here
#'   <https://github.com/drizopoulos/JMbayes/blob/master/R/dynPred_lme.R>.
#'   If \code{NULL} (default), method \code{'m1'} is automatically set. 
#' 
#' @return A data frame object. 
#' 
#' @inheritParams  growthparameters.bgmfit 
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#'
get.newdata <- function(model,
                        newdata = NULL,
                        xvar = NULL,
                        idvar = NULL,
                        resp = NULL,
                        dpar = NULL,
                        numeric_cov_at = NULL,
                        aux_variables = NULL,
                        levels_id = NULL,
                        ipts = NULL,
                        xrange = NULL,
                        idata_method = NULL,
                        dummy_to_factor = NULL,
                        newdata_fixed = NULL,
                        verbose = FALSE) {
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
 
  if(is.null(dpar)) {
    dpar <- "mu"
  }
  
  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  

  # Initiate non formalArgs()
  `:=` <- NULL
  . <- NULL;
  

  validate_response(model, resp)
  
  list_c <- list()
  
  
  # For sigma
  xvar_      <- paste0('xvar', resp_rev_)
  sigmaxvar_ <- paste0('sigma', xvar_)
  cov_       <- paste0('cov', resp_rev_)
  sigmacov_  <- paste0('sigma', cov_)
  uvarby     <- model$model_info$univariate_by$by
  if(is.null(uvarby)) uvarby <- NA 
 
  
  if(dpar == "mu") {
    if(is.null(xvar)) {
      xvar   <- model$model_info[[xvar_]]
    }
    cov    <- model$model_info[[cov_]]
  } else if(dpar == "sigma") {
    if(!is.na(model$model_info[[sigmaxvar_]])) {
      xvar   <- model$model_info[[sigmaxvar_]]
    } else if(is.na(model$model_info[[sigmaxvar_]]) & 
              !is.null(model$model_info[[xvar_]])) {
      xvar   <- model$model_info[[xvar_]]
    }
    cov    <- model$model_info[[sigmacov_]]
  } # if(dpar == "mu") { else if(dpar == "sigma") {
  
  groupvar_     <- paste0('groupvar', resp_rev_)
  yvar_         <- paste0('yvar', resp_rev_)
  yvar          <- model$model_info[[yvar_]]
  hierarchical_ <- paste0('hierarchical', resp_rev_)
  
  
  if(is.null(levels_id) & is.null(idvar)) {
    idvar <- model$model_info[[groupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      idvar <- model$model_info[[hierarchical_]]
    }
    # 29.08.2025 - re assign idvar to groupvar_ if hierarchical_
    model$model_info[[groupvar_]] <- idvar # idvar[1]
  } else if (!is.null(levels_id)) {
    idvar <- levels_id
  } else if (!is.null(idvar)) {
    idvar <- idvar
  }
  
  cov_       <- paste0('cov', resp_rev_)
  sigmacov_  <- paste0('sigma', cov_)
  
  # When no random effects and hierarchical, IDvar <- NULL problem 02 03 2024
  if(is.null(idvar)) {
    if(is.null(idvar)) {
      if(!is.null(model$model_info[['idvars']])) {
        idvar <- model$model_info[['idvars']]
      }
    }
  }
  
  
  if(dpar == "sigma") {
    sigma_model <- get_sigmamodel_info(model = model,
                                       newdata = NULL, # what = 'model' ignored 
                                       dpar = dpar, 
                                       resp = resp, 
                                       what = 'model',
                                       cov = NULL, 
                                       all = FALSE, 
                                       verbose = verbose)
    }
  
  
  
  if (is.null(newdata)) {
    if(idata_method == 'm1') newdata <- model$model_info$bgmfit.data
    if(idata_method == 'm2') newdata <- model$model_info$bgmfit.data
  } else {
    newdata <- newdata
  }
  
  
  if(data.table::is.data.table(newdata)) {
    newdata <- newdata %>% data.frame()
    newdata_was_data_table <- TRUE
  } else {
    newdata <- newdata
    newdata_was_data_table <- FALSE
  }
  
  
  
  
  # if(dpar == "mu") {
  #   if (is.null(newdata)) {
  #     if(idata_method == 'm1') newdata <- model$model_info$bgmfit.data
  #     if(idata_method == 'm2') newdata <- model$model_info$bgmfit.data
  #   } else {
  #     newdata <- newdata
  #   }
  # } else if(dpar == "sigma") {
  #   if(is.null(newdata)) {
  #     if(!is.null(sigma_model)) {
  #       if(sigma_model != "ls") {
  #         # newdata <- model$model_info$bgmfit.data
  #         newdata <- newdata
  #         if(verbose) {
  #           message("For dpar = 'sigma', the data used for model fitting is used")
  #         }
  #       }
  #     }
  #   }
  # }
  
  newdata.in <- newdata
  
  
  
  
  
  
  if(!is.null(dummy_to_factor)) {
    if(!is.list(dummy_to_factor)) {
      stop("dummy_to_factor must be a named list as follows:",
           "\n ",
           "dummy_to_factor = list(factor.dummy = , factor.name = )",
           "\n ",
           "where factor.dummy is a vector of character strings that will" ,
           "\n ",
           "be converted to a factor variable, and ",
           "\n ",
           "factor.name is a single character string that is used to name the",
           "\n ",
           "newly created factor variable"
      )
    }
    factor.dummy <- dummy_to_factor[['factor.dummy']]
    factor.level <- dummy_to_factor[['factor.level']]
    if(is.null(dummy_to_factor[['factor.name']])) {
      factor.name <- 'factor.dummy'
    } else {
      factor.name <- dummy_to_factor[['factor.name']]
    }
    newdata <- convert_dummy_to_factor(df = newdata, 
                                       factor.dummy = factor.dummy,
                                       factor.level = factor.level,
                                       factor.name = NULL)
    
    colnames(newdata) <- gsub('factor.var', factor.name, names(newdata))
    newdata[[factor.name]] <- as.factor(newdata[[factor.name]])
  }
  
  
  # This is when no random effects and this groupvar is NULL
  # Therefore, an artificial group var created
  # see also changes made to the get_idata function lines 17
  
  # if (is.null(model$model_info$groupvar)) {
  #   name_hypothetical_id <- paste0("id", resp_rev_)
  #   model$model_info$groupvar <- name_hypothetical_id
  #   newdata[[name_hypothetical_id]] <- as.factor("tempid")
  # } else if (!is.null(model$model_info$groupvar)) {
  #   if(length(newdata[[model$model_info$groupvar]]) == 0) {
  #     # name_hypothetical_id <- paste0("hy_id", resp_rev_)
  #     if(length(IDvar) > 1) {
  #       name_hypothetical_id <- IDvar[1] 
  #     } else {
  #       name_hypothetical_id <- IDvar
  #     }
  #     model$model_info$groupvar <- name_hypothetical_id
  #     newdata[[name_hypothetical_id]] <- as.factor("tempid")
  #   }
  # }
  
  newdata <- check_newdata_args(model, newdata, idvar, resp, verbose = verbose)
  
  if(dpar == "sigma") {
    if(is.null(newdata_fixed)) {
      if(!is.null(sigma_model)) {
        if(sigma_model != "ls") {
          # newdata_fixed <- 1
          # idata_method <- 'm2'
          if(verbose) {
            message("For dpar = 'sigma', the data used for model fitting is used")
          }
        }
      }
    }
  }
  
  

  # prepare_data2 with model = model will get all the necessary info
  
  if(dpar == "mu") {
    itransform_set <- 'x'
    transform_set  <- 'x'
  } else if(dpar == "sigma") {
    itransform_set <- 'sigma'
    transform_set  <- 'sigma'
  }
  
  
  
  add_just_list_c <- FALSE
  if(is.null(newdata_fixed)) {
    newdata <- prepare_data2(data = newdata, model = model)
    newdata <- prepare_transformations(data = newdata, 
                                       # itransform = itransform_set, 
                                       # transform = transform_set,
                                       model = model)
  } else if(!is.null(newdata_fixed)) {
    if(newdata_fixed == 0) {
      # do nothing assuming that user has set up for 'uvarby' etc
      # just apply 'dummy_to_factor'
      # only list_c elements will be added
      newdata         <- newdata 
      add_just_list_c <- TRUE
    } else if(newdata_fixed == 1) {
      newdata <- prepare_data2(data = newdata, model = model)
    } else if(newdata_fixed == 2) {
      newdata <- prepare_transformations(data = newdata, model = model)
    } else if(newdata_fixed == 3) {
      return(newdata) # i.e., not even applied 'dummy_to_factor' and return
    } else {
      stop("'newdata_fixed' should be either NULL or an integer, 1, 2, or 3")
    }
  }
  
  
  newdata <- newdata[,!duplicated(colnames(newdata))]
  
  
  
  ##################################################
  # prepare_data2 changes 
  
  if (!is.na(uvarby)) {
    subindicatorsi <- 
      model$model_info$subindicators[grep(resp,
                                          model$model_info$yvars)]
    aux_variables <- subindicatorsi
    list_c[['subindicatorsi']] <- subindicatorsi
    list_c[['uvarby']]         <- uvarby
  } 
  
  cov_vars       <-  model$model_info[[cov_]]
  sigmacov_vars  <-  model$model_info[[sigmacov_]]
  
  # Now instead of NULL, bsitar loop ii return NA when no covar
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
  

  if(is.null(cov_vars) & is.null(sigmacov_vars)) {
    if(!is.null(dummy_to_factor)) {
      warning("There are no covariate(s) but have specified dummy_to_factor",
              "\n ", 
              "Please check if this is an error")
    }
  }
  
  
  if(!is.null(dummy_to_factor)) {
    cov_vars <- c(cov_vars, factor.name)
  }
  
  
  # check if cov is charcater but not factor 
  checks_for_chr_fact <- c(cov_vars, sigmacov_vars)
  checks_for_chr_fact <- unique(checks_for_chr_fact)
  for (cov_varsi in checks_for_chr_fact) {
    if(is.character(newdata[[cov_varsi]])) {
      if(!is.factor(newdata[[cov_varsi]])) {
        if(verbose) {
          message("\nVariable '", cov_varsi, "' used as a covariate in the model ",
                  "\n ",
                  " is a character but not factor. Converting it to factor.")
        }
        newdata[[cov_varsi]] <- as.factor(newdata[[cov_varsi]])
        factor_vars <- names(newdata[sapply(newdata, is.factor)])
      }
    }
  }
  
  factor_vars  <- names(newdata[sapply(newdata, is.factor)])
  numeric_vars <- names(newdata[sapply(newdata, is.numeric)])
  

  cov_factor_vars  <- intersect(cov_vars, factor_vars)
  cov_numeric_vars <- intersect(cov_vars, numeric_vars)
  groupby_fstr     <- c(cov_factor_vars)
  groupby_fistr    <- c(idvar, cov_factor_vars)
  
  sigmacov_factor_vars  <- intersect(sigmacov_vars, factor_vars)
  sigmacov_numeric_vars <- intersect(sigmacov_vars, numeric_vars)
  
  sigmagroupby_fstr     <- c(sigmacov_factor_vars)
  sigmagroupby_fistr    <- c(idvar, sigmacov_factor_vars)

  
  if (identical(cov_factor_vars, character(0)))
    cov_factor_vars <- NULL
  if (identical(cov_numeric_vars, character(0)))
    cov_numeric_vars <- NULL
  
  if (identical(sigmacov_factor_vars, character(0)))
    sigmacov_factor_vars <- NULL
  if (identical(sigmacov_numeric_vars, character(0)))
    sigmacov_numeric_vars <- NULL

  cov_factor_vars  <- c(cov_factor_vars, sigmacov_factor_vars)
  cov_numeric_vars <- c(cov_numeric_vars, sigmacov_numeric_vars)
  
  if (!is.na(uvarby)) {
    if(idata_method == 'm1') groupby_fstr  <- c(uvarby, groupby_fstr)
    if(idata_method == 'm1') groupby_fistr <- c(uvarby, groupby_fistr)
  }
  
  
  if(dpar == "mu") {
    cov_vars         <- cov_vars
    cov_factor_vars  <- cov_factor_vars
    cov_numeric_vars <- cov_numeric_vars
    groupby_fstr     <- groupby_fstr
    groupby_fistr    <- groupby_fistr
  } else if(dpar == "sigma") {
    cov_vars         <- sigmacov_vars
    cov_factor_vars  <- sigmacov_factor_vars
    cov_numeric_vars <- sigmacov_numeric_vars
    groupby_fstr     <- sigmagroupby_fstr
    groupby_fistr    <- sigmagroupby_fistr
  }

  if(add_just_list_c) {
    list_c[['xvar']] <- xvar
    list_c[['yvar']] <- yvar
    list_c[['idvar']] <- idvar
    list_c[['cov_vars']] <- cov_vars
    list_c[['cov_factor_vars']] <- cov_factor_vars
    list_c[['cov_numeric_vars']] <- cov_numeric_vars
    list_c[['groupby_fstr']] <- groupby_fstr
    list_c[['groupby_fistr']] <- groupby_fistr
    if(is.null(list_c[['uvarby']])) {
      list_c[['uvarby']] <- NA
    }
    attr(newdata, 'list_c') <- list_c
    return(newdata)
  }
  
  
  set_numeric_cov_at <- function(x, numeric_cov_at) {
    name_ <- deparse(substitute(x))
    if (is.null((numeric_cov_at[[name_]]))) {
      . <- mean(x, na.rm = T)
    } else if (!is.null((numeric_cov_at[[name_]]))) {
      if (numeric_cov_at[[name_]] == 'mean') {
        . <- mean(x,  na.rm = T)
      } else if (numeric_cov_at[[name_]] == 'median') {
        . <- median(x, na.rm = T)
      } else if (numeric_cov_at[[name_]] == 'min') {
        . <- min(x, na.rm = T)
      } else if (numeric_cov_at[[name_]] == 'max') {
        . <- max(x, na.rm = T)
      } else {
        . <- numeric_cov_at[[name_]]
      }
    }
    round(., 3)
  }
  
  
  
  get.data.grid <- function(data,
                            xvar,
                            yvar,
                            idvar,
                            cov_numeric_vars,
                            numeric_cov_at,
                            aux_variables,
                            uvarby) {
    if (!is.null(idvar))
      relocate_vars <- c(xvar, idvar)
    if (is.null(idvar))
      relocate_vars <- c(xvar)
    if (!is.na(uvarby))
      if(idata_method == 'm1') {
        relocate_vars <- c(relocate_vars, uvarby)
      } else if(idata_method == 'm2') {
        relocate_vars <- c(relocate_vars)
      }
    if (!is.null(cov_numeric_vars)) {
      cov_numeric_vars__ <- cov_numeric_vars
      if (identical(cov_numeric_vars__, character(0)))
        cov_numeric_vars__ <- NULL
      if (!is.null(cov_numeric_vars__)) {
        for (cov_numeric_vars__i in cov_numeric_vars__) {
          if (!is.null(numeric_cov_at)) {
            if (!cov_numeric_vars__i %in% names(numeric_cov_at)) {
              stop(
                "You have used the argument 'numeric_cov_at' to specify the",
                "\n ",
                " value of continous covariate '",
                cov_numeric_vars__i,
                "'.",
                "\n ",
                " However, the name of this covariate is missing from the list",
                "\n ",
                " Please use argument 'numeric_cov_at' correctly as follows:",
                "\n ",
                " numeric_cov_at = list(",
                cov_numeric_vars__i,
                " = xx)"
              )
            }
          }
        }
        data <-
          data %>% dplyr::mutate_at(cov_numeric_vars__,
                                    set_numeric_cov_at,
                                    numeric_cov_at)
        
        # This is good but get.data.grid is called twice - why?
        # hence this cat("\n"... is printed twice
        # cat("Continous covariate(s) set at:\n")
        # for (cov_numeric_vars__i in cov_numeric_vars__) {
        #   cat("\n", cov_numeric_vars__i, "at",
        #       unique(data[[cov_numeric_vars__i]]))
        # }
      }
    }
    
    
    if (!is.null(yvar)) {
      if (yvar %in% colnames(data)) {
        relocate_vars <- c(yvar, relocate_vars)
      }
    }
    
    data %>% 
      dplyr::relocate(dplyr::all_of(relocate_vars)) %>% 
      data.frame()
  }
  
  ########
  

  i_data <-
    function(model,
             newdata,
             xvar = NULL,
             idvar = NULL,
             resp = NULL,
             cov_factor_vars = NULL,
             cov_numeric_vars = NULL,
             aux_variables = NULL,
             levels_id = NULL,
             ipts = NULL,
             xrange = NULL) {
      if (is.null(resp)) {
        resp_rev_ <- resp
      } else if (!is.null(resp)) {
        resp_rev_ <- paste0("_", resp)
      }
      xvar_ <- paste0('xvar', resp_rev_)
      yvar_ <- paste0('yvar', resp_rev_)
      groupvar_ <- paste0('groupvar', resp_rev_)
      if(is.null(xvar)) xvar <- model$model_info[[xvar_]]
      yvar <- model$model_info[[yvar_]]
     
      hierarchical_ <- paste0('hierarchical', resp_rev_)
      
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
      
      if (!is.na(uvarby))
        cov_factor_vars <- c(uvarby, cov_factor_vars)
      
      
      
      # this idatafunction i.e., 'm1'
      if (idata_method == 'm1') {
        idatafunction <- function(.x,
                                  xvar,
                                  idvar,
                                  nmy,
                                  xrange,
                                  set_xrange,
                                  aux_var = NULL) {
          index__x <- NA
          exdata <-
            function(x,
                     id,
                     idmat,
                     nmy,
                     xrange,
                     set_xrange,
                     aux_var) {
              n <- round(nmy * diff(range(x)))
              npt <- n / diff(range(x))
              
              extage <- apply(idmat, 1, function(x1) {
                index__x <-
                  id == x1
                
                if (is.null(xrange)) {
                  id.x <- x[index__x]
                }
                if (!is.null(xrange)) {
                  if (length(xrange) == 1) {
                    if (xrange == 1 & is.null(set_xrange))
                      id.x <- x
                    if (xrange == 2 &
                        !is.null(set_xrange))
                      id.x <- set_xrange
                  }
                  if (length(xrange) == 2) {
                    id.x <- set_xrange
                  }
                }
                
                nt <- floor(npt * diff(range(id.x))) + 1
                newx <- seq(min(id.x), max(id.x), length = nt)
                
                newid <-
                  rep(x1, nt)
                extx <- data.frame(x = newx, id = newid)
                colnames(extx) <- c("x", "id")
                extx
              })
              df <-
                extage[[1]][FALSE,]
              for (dft in extage)
                df <- rbind(df, dft)
              df
            }
          
          inidnull <- FALSE
          if(is.null(.x[[idvar]])) {
            inidnull <- TRUE
            .x[[idvar]] <- unique(levels(newdata[[idvar]]))[1]
          }
          
          out <- exdata(
            x = .x[[xvar]],
            id = .x[[idvar]],
            idmat = matrix(unique(.x[[idvar]], ncol = 1)),
            nmy = nmy,
            xrange = xrange,
            set_xrange = set_xrange,
            aux_var = aux_var
          )
          
          out <- out %>% dplyr::rename(!!idvar := 'id') %>% data.frame()
          
          if(inidnull) out <- out %>% dplyr::select(-dplyr::all_of(idvar))
          
          idxx <- NULL
          if (!is.null(aux_var)) {
            aux_varx <- c(aux_var, idvar)
            newx. <- .x %>% dplyr::select(dplyr::all_of(aux_varx)) %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(idvar))) %>%
              dplyr::mutate(idxx = dplyr::row_number()) %>%
              dplyr::ungroup()
            outx. <- out %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(idvar))) %>%
              dplyr::mutate(idxx = dplyr::row_number()) %>%
              dplyr::ungroup()
            out <- outx. %>% dplyr::left_join(., newx.,
                                              by = c(idvar, 'idxx')) %>%
              dplyr::select(-idxx) %>% data.frame()
          }
          
          out 
        } # end idatafunction -> m1
        
       
        
        if (!is.null(xrange)) {
          if (length(xrange) < 1 | length(xrange) > 2) {
            stop(
              "Argument xrange should be either NULL, numeric value 1 or 2",
              "\n ",
              "or else a paired values indicating the range e.g., c(6, 20)"
            )
          }
        }
        
        
        if (!is.null(xrange)) {
          if (length(xrange) == 1) {
            if (xrange == 1)
              set_xrange <- NULL
            if (xrange == 2)
              set_xrange <- range(newdata[[xvar]])
          }
          if (length(xrange) == 2) {
            set_xrange <- xrange
          }
        }
        
        if (is.null(xrange))
          set_xrange <- NULL
       
        
        

        if (is.null(model$model_info[[hierarchical_]])) {
          if (!is.null(ipts) & is.null(cov_factor_vars)) {
            # 20.03.2025 -> need to keep all variables 
            all_names_x   <- colnames(newdata)
            core_names_x  <- c(idvar, xvar)
            aux_variables <- c(aux_variables, setdiff(all_names_x, core_names_x))
            newdata %>% dplyr::arrange(idvar, xvar) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  idvar = idvar,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              dplyr::mutate(!!idvar := as.factor(eval(parse(text = idvar)))) %>%
              dplyr::relocate(dplyr::all_of(idvar), dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          } else if (!is.null(ipts) & !is.null(cov_factor_vars)) {
            # 20.03.2025 -> need to keep all variables 
            all_names_x   <- colnames(newdata)
            core_names_x  <- c(idvar, xvar, cov_factor_vars)
            aux_variables <- c(aux_variables, setdiff(all_names_x, core_names_x))
            newdata %>% dplyr::arrange(idvar, xvar) %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(cov_factor_vars))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  idvar = idvar,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              dplyr::mutate(!!idvar := as.factor(eval(parse(text = idvar)))) %>%
              dplyr::relocate(dplyr::all_of(idvar), dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          }
        } # if(is.null(model$model_info[[hierarchical_]]))
        
       
        
        multiNewVar <- function(df, df2, varname) {
          df %>% dplyr::mutate(.,!!varname := df2[[varname]])
        }
        
        if (!is.null(model$model_info[[hierarchical_]])) {
          if (!is.null(ipts) & is.null(cov_factor_vars)) {
            IDvar_ <- idvar[1]
            higher_ <- idvar[2:length(idvar)]
            arrange_by <- c(IDvar_, xvar)
            cov_factor_vars_by <- c(higher_, cov_factor_vars)
            newdata_o <- newdata
            # 20.03.2025 -> need to keep all variables 
            all_names_x   <- colnames(newdata)
            core_names_x  <- c(IDvar_, xvar)
            aux_variables <- c(aux_variables, setdiff(all_names_x, core_names_x))
            newdata <-
              newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
              dplyr::group_by(
                dplyr::across(dplyr::all_of(cov_factor_vars_by))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  idvar = IDvar_,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                ),
                .keep = F
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              data.frame()
            
            for (i in idvar) {
              newdata <- newdata %>% multiNewVar(df = .,
                                                 df2 = newdata,
                                                 varname = i)
            }
            
            newdata %>% dplyr::relocate(dplyr::all_of(idvar), 
                                        dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          }
          
          if (!is.null(ipts) & !is.null(cov_factor_vars)) {
            IDvar_ <- idvar[1]
            higher_ <- idvar[2:length(idvar)]
            arrange_by <- c(IDvar_, xvar)
            # cov_factor_vars_by <- c(higher_, cov_factor_vars)
            if(length(idvar) > 1) {
              cov_factor_vars_by <- c(higher_, cov_factor_vars)
            } else {
              cov_factor_vars_by <- c(cov_factor_vars)
            }
            all_names_x   <- colnames(newdata)
            core_names_x  <- c(IDvar_, xvar, cov_factor_vars_by)
            aux_variables <- c(aux_variables, setdiff(all_names_x, core_names_x))
            newdata <-
              newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
              dplyr::group_by(
                dplyr::across(dplyr::all_of(cov_factor_vars_by))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  idvar = IDvar_,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>% data.frame()
            for (i in idvar) {
              newdata <- newdata %>% multiNewVar(df = .,
                                                 df2 = newdata,
                                                 varname = i)
            }
            newdata %>% dplyr::relocate(dplyr::all_of(idvar), 
                                        dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          }
        } # if(!is.null(model$model_info[[hierarchical_]])) {
        
      } # end of if(idata_method == 'm1') {
      
      
      # this is get_idata i.e., 'm2'
      if (idata_method == 'm2') {
        if (!is.null(ipts)) {
          # for 3 or more level data, idvar shoud be first of vector
          if (is.null(model$model_info[[hierarchical_]])) {
            IDvar_for_idata <- idvar
          } else if (!is.null(model$model_info[[hierarchical_]])) {
            IDvar_for_idata <- idvar[1]
          }
          newdata <-
            get_idata(
              model = NULL,
              newdata = newdata,
              idvar = IDvar_for_idata,
              xvar = xvar,
              times = NULL,
              length.out = ipts,
              xrange = xrange,
              keeplevels = FALSE, 
              asdf = FALSE,
              newdata_fixed = 0
            )
        }
        
      } # end if(idata_method == 'm2') {
      
      
      
      if (is.null(ipts)) {
        newdata <- newdata
      }
     
      ##################################################
      newdata.get.data.grid <- get.data.grid(
        data = newdata,
        xvar = xvar,
        yvar = yvar,
        idvar = idvar,
        cov_numeric_vars = cov_numeric_vars,
        numeric_cov_at = numeric_cov_at,
        aux_variables = aux_variables,
        uvarby = uvarby
      )
      
      j_b_names <- intersect(names(newdata), names(newdata.get.data.grid))
      j_b_names__ <- c(j_b_names, cov_numeric_vars)
      j_b_names__ <- unique(j_b_names__)
      
      
      if(idata_method == 'm1') {
        newdata <-
          newdata %>% 
          dplyr::left_join(., 
                           newdata.get.data.grid %>%
                             dplyr::select(dplyr::all_of(j_b_names__)),
                           by = j_b_names)
        
      } else if(idata_method == 'm2') {
        if(length(idvar) > 1) {
          name_hypothetical_id <- idvar[1]
        } else {
          name_hypothetical_id <- idvar
        }
        if(length(unique(newdata[[name_hypothetical_id]])) > 1) {
          newdata <-
            newdata %>% 
            dplyr::left_join(., 
                             newdata.get.data.grid %>%
                               dplyr::select(dplyr::all_of(j_b_names__)),
                             by = j_b_names)
        } else if(length(unique(newdata[[name_hypothetical_id]])) == 1) {
          newdata <- newdata
        }
      } # else if(idata_method == 'm2') {
      
      return(newdata)
    }
  
  
  newdata <- i_data(
    model,
    newdata,
    xvar = xvar,
    idvar = idvar,
    resp = resp,
    cov_factor_vars = cov_factor_vars,
    cov_numeric_vars = cov_numeric_vars,
    aux_variables = aux_variables,
    levels_id = levels_id,
    ipts = ipts,
    xrange = xrange
  )
  
  if(dpar == "mu") {
    cov_vars         <- cov_vars
    cov_factor_vars  <- cov_factor_vars
    cov_numeric_vars <- cov_numeric_vars
    groupby_fstr     <- groupby_fstr
    groupby_fistr    <- groupby_fistr
  } else if(dpar == "sigma") {
    cov_vars         <- sigmacov_vars
    cov_factor_vars  <- sigmacov_factor_vars
    cov_numeric_vars <- sigmacov_numeric_vars
    groupby_fstr     <- sigmagroupby_fstr
    groupby_fistr    <- sigmagroupby_fistr
  }
  
  list_c[['xvar']] <- xvar
  list_c[['yvar']] <- yvar
  list_c[['idvar']] <- idvar
  list_c[['cov_vars']] <- cov_vars
  list_c[['cov_factor_vars']] <- cov_factor_vars
  list_c[['cov_numeric_vars']] <- cov_numeric_vars
  list_c[['groupby_fstr']] <- groupby_fstr
  list_c[['groupby_fistr']] <- groupby_fistr
  if(is.null(list_c[['uvarby']])) {
    list_c[['uvarby']] <- NA
  }
  attr(newdata, 'list_c') <- list_c
  
  if(newdata_was_data_table) {
    newdata <- data.table::as.data.table(newdata)
  }
  
  return(newdata)
} # End get.newdata




