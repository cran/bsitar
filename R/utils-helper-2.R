

#' An internal function to create data frame for post-processing
#'
#' @inheritParams  growthparameters.bgmfit 
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
#' @keywords internal
#' 
#' @return A data frame object. 
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#'
get.newdata <- function(model,
                        newdata = NULL,
                        resp = NULL,
                        numeric_cov_at = NULL,
                        aux_variables = NULL,
                        levels_id = NULL,
                        ipts = NULL,
                        xrange = NULL,
                        idata_method = NULL,
                        dummy_to_factor = NULL,
                        verbose = FALSE,
                        ...) {
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  

  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  
  
  # Initiate non formalArgs()
  `:=` <- NULL
  . <- NULL;
  
  
  validate_response(model, resp)
  
  list_c <- list()
  xvar_ <- paste0('xvar', resp_rev_)
  yvar_ <- paste0('yvar', resp_rev_)
  groupvar_ <- paste0('groupvar', resp_rev_)
  xvar <- model$model_info[[xvar_]]
  yvar <- model$model_info[[yvar_]]
  hierarchical_ <- paste0('hierarchical', resp_rev_)
  
  if (is.null(levels_id)) {
    IDvar <- model$model_info[[groupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      IDvar <- model$model_info[[hierarchical_]]
    } else if (is.null(model$model_info[[hierarchical_]])) {
      # if(!is.null(model$model_info[['ids']])) {
      #   IDvar <- model$model_info[['ids']]
      # }
    }
  } else if (!is.null(levels_id)) {
    IDvar <- levels_id
  }
  
  xfun_ <- paste0('xfun', resp_rev_)
  yfun_ <- paste0('yfun', resp_rev_)
  xfun <- model$model_info[[xfun_]]
  yfun <- model$model_info[[yfun_]]
  
  
  cov_ <- paste0('cov', resp_rev_)
  cov_sigma_ <- paste0('cov_sigma', resp_rev_)
  uvarby <- model$model_info$univariate_by
  
  
  # When no random effects and hierarchical, IDvar <- NULL problem 02 03 2024
  #if(idata_method == 'm2') {
    if(is.null(levels_id)) {
      if(is.null(IDvar)) {
        if(!is.null(model$model_info[['ids']])) {
          IDvar <- model$model_info[['ids']]
        }
      }
    }
  #}
  
 
  
  if (is.null(newdata)) {
    if(idata_method == 'm1') newdata <- model$model_info$bgmfit.data
    if(idata_method == 'm2') newdata <- model$data
  } else {
    newdata <- newdata
  }
  
  
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
  
  newdata <- check_newdata_args(model, newdata, IDvar, resp)

  
  if (!is.na(model$model_info$univariate_by)) {
    if (is.symbol(model$model_info$call.bgmfit$y)) {
      setorgy <- deparse(model$model_info$call.bgmfit$y)
    } else if (is.list(model$model_info$call.bgmfit$y)) {
      setorgy <- unname(unlist(model$model_info$call.bgmfit$y))
      if (is.symbol(setorgy))
        setorgy <- deparse(setorgy)
    } else {
      setorgy <- model$model_info$call.bgmfit$y
    }
  }
  
  if (is.na(model$model_info$univariate_by)) {
    setorgy <- model$model_info$ys
  }
 
  
  if (idata_method == 'm1') {
    newdata <- prepare_data(
      data = newdata,
      x = model$model_info$xs,
      y = setorgy,
      id = model$model_info$ids,
      uvarby = model$model_info$univariate_by,
      mvar = model$model_info$multivariate,
      xfuns = model$model_info$xfuns,
      yfuns = model$model_info$yfuns,
      outliers = model$model_info$outliers)
  }
  
  
 
  newdata <- newdata[,!duplicated(colnames(newdata))]
  
  if (!is.na(model$model_info$univariate_by)) {
    if(idata_method == 'm1') {
      sortbylayer <- NA
      newdata <- newdata %>%
        dplyr::mutate(sortbylayer =
                        forcats::fct_relevel(!!as.name(uvarby),
                                             (levels(
                                               !!as.name(uvarby)
                                             )))) %>%
        dplyr::arrange(sortbylayer) %>%
        dplyr::mutate(!!as.name(IDvar) := factor(!!as.name(IDvar),
                                                 levels =
                                                   unique(!!as.name(IDvar)))) %>%
        dplyr::select(-sortbylayer)
    } # if(idata_method == 'm1') {
    
    subindicatorsi <- model$model_info$subindicators[grep(resp,
                                                          model$model_info$ys)]
    list_c[['subindicatorsi']] <- subindicatorsi
    list_c[['uvarby']] <- uvarby
  }
  
  
  covars_extrcation <- function(str) {
    str <- gsub("[[:space:]]", "", str)
    for (ci in c("*", "+", ":")) {
      str <- gsub(ci, ' ', str, fixed = T)
    }
    str <- strsplit(str, " ")[[1]]
    str
  }
  
  
  cov_vars <-  model$model_info[[cov_]]
  cov_sigma_vars <-  model$model_info[[cov_sigma_]]

  
  if (!is.null(cov_vars)) {
    cov_vars <- covars_extrcation(cov_vars)
  }
  if (!is.null(cov_sigma_vars)) {
    cov_sigma_vars <- covars_extrcation(cov_sigma_vars)
  }
  
  
  if(is.null(cov_vars) & is.null(cov_sigma_vars)) {
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
  checks_for_chr_fact <- c(cov_vars, cov_sigma_vars)
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
    
  factor_vars <- names(newdata[sapply(newdata, is.factor)])
  numeric_vars <- names(newdata[sapply(newdata, is.numeric)])
  
  # if (!is.null(cov_sigma_vars))
  #   cov_sigma_vars <- covars_extrcation(cov_sigma_vars)
  
  
  
  cov_factor_vars <- intersect(cov_vars, factor_vars)
  cov_numeric_vars <- intersect(cov_vars, numeric_vars)
  groupby_fstr <- c(cov_factor_vars)
  groupby_fistr <- c(IDvar, cov_factor_vars)
  
  cov_sigma_factor_vars <- intersect(cov_sigma_vars, factor_vars)
  cov_sigma_numeric_vars <- intersect(cov_sigma_vars, numeric_vars)
  
  if (identical(cov_factor_vars, character(0)))
    cov_factor_vars <- NULL
  if (identical(cov_numeric_vars, character(0)))
    cov_numeric_vars <- NULL
  
  if (identical(cov_sigma_factor_vars, character(0)))
    cov_sigma_factor_vars <- NULL
  if (identical(cov_sigma_numeric_vars, character(0)))
    cov_sigma_numeric_vars <- NULL
  
  # Merge here a b c covariate with sigma co variate
  # IMP: Note that groupby_fstr and groupby_fistr are stil  a b c covariate
  # This way, plot_curves and gparameters will not produce sigam cov specific
  # curves and g parameters
  
  cov_factor_vars <- c(cov_factor_vars, cov_sigma_factor_vars)
  cov_numeric_vars <- c(cov_numeric_vars, cov_sigma_numeric_vars)
  
  if (!is.na(model$model_info$univariate_by)) {
    if(idata_method == 'm1') groupby_fstr <- c(uvarby, groupby_fstr)
    if(idata_method == 'm1') groupby_fistr <- c(uvarby, groupby_fistr)
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
                            IDvar,
                            cov_numeric_vars,
                            numeric_cov_at,
                            aux_variables,
                            uvarby) {
    if (!is.null(IDvar))
      relocate_vars <- c(xvar, IDvar)
    if (is.null(IDvar))
      relocate_vars <- c(xvar)
    if (!is.na(uvarby))
      if(idata_method == 'm1') relocate_vars <- c(relocate_vars, uvarby)
      if(idata_method == 'm2') relocate_vars <- c(relocate_vars)
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
    data %>% dplyr::relocate(dplyr::all_of(relocate_vars)) %>% data.frame()
    
  }
  
  ########
  
  i_data <-
    function(model,
             newdata,
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
      xvar <- model$model_info[[xvar_]]
      yvar <- model$model_info[[yvar_]]
      
      hierarchical_ <- paste0('hierarchical', resp_rev_)
      if (is.null(levels_id)) {
        IDvar <- model$model_info[[groupvar_]]
        if (!is.null(model$model_info[[hierarchical_]])) {
          IDvar <- model$model_info[[hierarchical_]]
        }
      } else if (!is.null(levels_id)) {
        IDvar <- levels_id
      }
      
      uvarby <- model$model_info$univariate_by
      if (!is.na(uvarby))
        cov_factor_vars <- c(uvarby, cov_factor_vars)
      
      
      
      # this idatafunction i.e., 'm1'
      if (idata_method == 'm1') {
        idatafunction <- function(.x,
                                  xvar,
                                  IDvar,
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
          if(is.null(.x[[IDvar]])) {
            inidnull <- TRUE
            .x[[IDvar]] <- unique(levels(newdata[[IDvar]]))[1]
          }
          
          out <- exdata(
            x = .x[[xvar]],
            id = .x[[IDvar]],
            idmat = matrix(unique(.x[[IDvar]], ncol = 1)),
            nmy = nmy,
            xrange = xrange,
            set_xrange = set_xrange,
            aux_var = aux_var
          )
          
          out <- out %>% dplyr::rename(!!IDvar := 'id') %>% data.frame()
          
          if(inidnull) out <- out %>% dplyr::select(-dplyr::all_of(IDvar))
          
          idxx <- NULL
          if (!is.null(aux_var)) {
            aux_varx <- c(aux_var, IDvar)
            newx. <- .x %>% dplyr::select(dplyr::all_of(aux_varx)) %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(IDvar))) %>%
              dplyr::mutate(idxx = dplyr::row_number()) %>%
              dplyr::ungroup()
            outx. <- out %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(IDvar))) %>%
              dplyr::mutate(idxx = dplyr::row_number()) %>%
              dplyr::ungroup()
            out <- outx. %>% dplyr::left_join(., newx.,
                                              by = c(IDvar, 'idxx')) %>%
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
            newdata %>% dplyr::arrange(IDvar, xvar) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  IDvar = IDvar,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              dplyr::mutate(!!IDvar := as.factor(eval(parse(text = IDvar)))) %>%
              dplyr::relocate(dplyr::all_of(IDvar), dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          } else if (!is.null(ipts) & !is.null(cov_factor_vars)) {
            newdata %>% dplyr::arrange(IDvar, xvar) %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(cov_factor_vars))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  IDvar = IDvar,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              dplyr::mutate(!!IDvar := as.factor(eval(parse(text = IDvar)))) %>%
              dplyr::relocate(dplyr::all_of(IDvar), dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          }
        } # if(is.null(model$model_info[[hierarchical_]]))
        
        
        multiNewVar <- function(df, df2, varname) {
          df %>% dplyr::mutate(.,!!varname := df2[[varname]])
        }
        
        if (!is.null(model$model_info[[hierarchical_]])) {
          if (!is.null(ipts) & is.null(cov_factor_vars)) {
            IDvar_ <- IDvar[1]
            higher_ <- IDvar[2:length(IDvar)]
            arrange_by <- c(IDvar_, xvar)
            cov_factor_vars_by <- c(higher_, cov_factor_vars)
            newdata_o <- newdata
            newdata <-
              newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
              dplyr::group_by(
                dplyr::across(dplyr::all_of(cov_factor_vars_by))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  IDvar = IDvar_,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                ),
                .keep = F
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              data.frame()
            
            for (i in IDvar) {
              newdata <- newdata %>% multiNewVar(df = .,
                                                 df2 = newdata,
                                                 varname = i)
            }
            
            newdata %>% dplyr::relocate(dplyr::all_of(IDvar), 
                                        dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          }
          
          if (!is.null(ipts) & !is.null(cov_factor_vars)) {
            IDvar_ <- IDvar[1]
            higher_ <- IDvar[2:length(IDvar)]
            arrange_by <- c(IDvar_, xvar)
            # cov_factor_vars_by <- c(higher_, cov_factor_vars)
            if(length(IDvar) > 1) {
              cov_factor_vars_by <- c(higher_, cov_factor_vars)
            } else {
              cov_factor_vars_by <- c(cov_factor_vars)
            }
            newdata <-
              newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
              dplyr::group_by(
                dplyr::across(dplyr::all_of(cov_factor_vars_by))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  IDvar = IDvar_,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>% data.frame()
            for (i in IDvar) {
              newdata <- newdata %>% multiNewVar(df = .,
                                                 df2 = newdata,
                                                 varname = i)
            }
            newdata %>% dplyr::relocate(dplyr::all_of(IDvar), 
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
            IDvar_for_idata <- IDvar
          } else if (!is.null(model$model_info[[hierarchical_]])) {
            IDvar_for_idata <- IDvar[1]
          }
          newdata <-
            get_idata(
              newdata = newdata,
              idVar = IDvar_for_idata,
              timeVar = xvar,
              times = NULL,
              length.out = ipts,
              xrange = xrange
            )
        }
        
      } # end if(idata_method == 'm2') {
      
      

      if (is.null(ipts)) {
        newdata <- newdata
      }
        
      
      
      if (!is.null(ipts)) {
        # outliers must be NULL
        # Because these has already been taken care of by get.newdata
        if (!is.na(model$model_info$univariate_by)) {
          if (is.symbol(model$model_info$call.bgmfit$y)) {
            setorgy <- deparse(model$model_info$call.bgmfit$y)
          } else if (is.list(model$model_info$call.bgmfit$y)) {
            setorgy <- unname(unlist(model$model_info$call.bgmfit$y))
            if (is.symbol(setorgy))
              setorgy <- deparse(setorgy)
          } else {
            setorgy <- model$model_info$call.bgmfit$y
          }
        }
        
        if (is.na(model$model_info$univariate_by)) {
          setorgy <- model$model_info$ys
        }
        
        if (idata_method == 'm1') {
          newdata <- prepare_data(
            data = newdata,
            x = model$model_info$xs,
            y = setorgy,
            id = model$model_info$ids,
            uvarby = model$model_info$univariate_by,
            mvar = model$model_info$multivariate,
            xfuns = model$model_info$xfuns,
            yfuns = model$model_info$yfuns,
            outliers = NULL) # model$model_info$outliers
        } # if (idata_method == 'm1') {
        
      }
      
      
      if (!is.na(model$model_info$univariate_by)) {
        if(idata_method == 'm1') {
          sortbylayer <- NA
          unique_names <- unique(names(newdata))
          newdata <- newdata %>% dplyr::select(dplyr::all_of(unique_names))
          newdata <- newdata %>%
            dplyr::mutate(sortbylayer =
                            forcats::fct_relevel(!!as.name(uvarby),
                                                 (levels(
                                                   !!as.name(uvarby)
                                                 )))) %>%
            dplyr::arrange(sortbylayer) %>%
            dplyr::mutate(!!as.name(IDvar) :=
                            factor(!!as.name(IDvar),
                                   levels =
                                     unique(!!as.name(IDvar)))) %>%
            dplyr::select(-sortbylayer)
        } # if(idata_method == 'm1') {
      }
      
      newdata.oo <- get.data.grid(
        data = newdata,
        xvar = xvar,
        yvar = yvar,
        IDvar = IDvar,
        cov_numeric_vars = cov_numeric_vars,
        numeric_cov_at = numeric_cov_at,
        aux_variables = aux_variables,
        uvarby = uvarby
      )
      
      
      #         newdata <- check_newdata_args(model, newdata, IDvar, resp)
      
      j_b_names <- intersect(names(newdata), names(newdata.oo))
      j_b_names__ <- c(j_b_names, cov_numeric_vars)
      j_b_names__ <- unique(j_b_names__)
      
      if(idata_method == 'm1') {
        newdata <-
          newdata %>% 
          dplyr::left_join(., 
                           newdata.oo %>%
                             dplyr::select(dplyr::all_of(j_b_names__)),
                           by = j_b_names)
        
      } else if(idata_method == 'm2') {
        if(length(IDvar) > 1) {
          name_hypothetical_id <- IDvar[1]
        } else {
          name_hypothetical_id <- IDvar
        }
        if(length(unique(newdata[[name_hypothetical_id]])) > 1) {
          newdata <-
            newdata %>% 
            dplyr::left_join(., 
                             newdata.oo %>%
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
    resp = resp,
    cov_factor_vars = cov_factor_vars,
    cov_numeric_vars = cov_numeric_vars,
    aux_variables = aux_variables,
    levels_id = levels_id,
    ipts = ipts,
    xrange = xrange
  )
  
 
  list_c[['xvar']] <- xvar
  list_c[['yvar']] <- yvar
  list_c[['IDvar']] <- IDvar
  list_c[['cov_vars']] <- cov_vars
  list_c[['cov_factor_vars']] <- cov_factor_vars
  list_c[['cov_numeric_vars']] <- cov_numeric_vars
  list_c[['groupby_fstr']] <- groupby_fstr
  list_c[['groupby_fistr']] <- groupby_fistr
  
  attr(newdata, 'list_c') <- list_c
  
  return(newdata)
} # get.newdata







#' An internal function to imterpolate data for plotting smooth curves
#' 
#' @param model An object of class \code{bgmfit}. This is optional (default
#'   \code{NULL}) i.e., it is not neccessary to specify the model object. When
#'   \code{model} is specified, then values for \code{newdata}, \code{idVar},
#'   and \code{timeVar} are automatically taken from the \code{model}.
#'
#' @param newdata A data frame. If \code{NULL} (default), data analysed in the
#'   original model fit is used.
#'
#' @param idVar A character string to specify the group identifier. If
#'   \code{NULL} (default), \code{id} from the model fit is used.
#'
#' @param timeVar  A character string to specify the time variable. If
#'   \code{NULL} (default), \code{x} from the model fit is used.
#'
#' @param times  A numeric vector to specify the time range. Currently ignored.
#'
#' @param length.out A numeric value to specify the length of interpolation
#'   points. Default 10.
#'
#' @param xrange An integer to set the predictor range (i.e., age) when
#'   executing the interpolation via \code{ipts}. The default \code{NULL} sets
#'   the individual specific predictor range whereas code \code{xrange = 1} sets
#'   same range for all individuals within the higher order grouping variable
#'   (e.g., study). Code \code{xrange  = 2} sets the identical range
#'   dplyr::across the entire sample. Lastly, a paired numeric values can be
#'   supplied e.g., \code{xrange = c(6, 20)} will set the range between 6 and
#'   20.
#' @param keeplevels A logical in case factor variables other than \code{idVar}
#'   present
#' @return A data frame.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#'
get_idata <-
  function(model = NULL,
           newdata = NULL,
           idVar = NULL,
           timeVar = NULL,
           times = NULL,
           length.out = 10,
           xrange = 1, 
           keeplevels = FALSE, 
           asdf = FALSE) {
    
    if (is.null(newdata)) {
      newdata <- model$data
    } else {
      newdata <- newdata
    }
    
    if(data.table::is.data.table(newdata)) {
      setasdt <- TRUE 
      newdata <- as.data.frame(newdata)
    } else {
      setasdt <- FALSE
    }
    
    if(keeplevels) {
      is.fact <- names(newdata[, sapply(newdata, is.factor)])
      cnames  <- colnames(newdata)
    }
    
    
    if(is.null(model)) {
      if (is.null(idVar)) stop("Specify model or idVar, both can not be NULL")
      if (is.null(timeVar)) stop("Specify model or timeVar, both can't be NULL")
    }
    
    if(!is.null(model)) {
      if (!is.null(idVar)) stop("Specify either model or idVar, not both")
      if (!is.null(timeVar)) stop("Specify either model or timeVar, not both")
    }
    
    if(!is.null(model)) {
      if(length(model$model_info$ids) > 1) {
        stop("Please specify newdata, idVar, and timeVar manullay because 
            currently value for these can not be infered from model with three 
            or more levels of hierarchy")
      }
    }
    
    `.` <- NULL;
    
    if (is.null(idVar)) {
      idVar <- model$model_info$ids
    } else {
      idVar <- idVar
    }
    
    if (is.null(timeVar)) {
      timeVar <- model$model_info$xvar
    } else {
      timeVar <- timeVar
    }
    
    all_times <- TRUE
    if (is.null(xrange))
      xrange <- 1
    else
      xrange <- xrange
    times_orig <- newdata[[timeVar]]
    times_orig <- times_orig[!is.na(times_orig)]
    
    if (is.null(times) || !is.numeric(times)) {
      times <-
        seq(min(times_orig), max(times_orig), length.out = length.out)
    }
    
    # This is when no random effects and groupvar is NULL
    # Therefore, an artificial group var created
    # Check utils-helper function lines 60
    
    if (nlevels(newdata[[idVar]]) == 1) {
      newdata <- newdata %>%
        dplyr::distinct(newdata[[timeVar]], .keep_all = T) %>%
        dplyr::arrange(!!as.name(timeVar))
    }
    
    id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
    
    if(length( unique(newdata[[idVar]])) == 1) {
      if(length.out == 1) stop("The argument 'ipts' should be > 1")
    }
    
    last_time <- tapply(newdata[[timeVar]], id, max)
    first_time <- tapply(newdata[[timeVar]], id, min)
    
    newdata_nomiss <- newdata[complete.cases(newdata),]
    id_nomiss <-
      match(newdata_nomiss[[idVar]], unique(newdata_nomiss[[idVar]]))
    n <- length(unique(id_nomiss))
    
    if (xrange == 1) {
      times_to_pred <- list()
      for (i in 1:length(unique(newdata[[idVar]]))) {
        numx <- as.character(i)
        times_to_pred[[numx]] <-
          seq(first_time[i], last_time[i], length.out = length.out)
      }
    }
    
    if (xrange == 2) {
      times_to_pred <- lapply(last_time, function (t)
        if (all_times)
          times
        else
          times[times > t])
    }
    
    id_pred <- rep(seq_len(n), sapply(times_to_pred, length))
    
    right_rows <- function (data, times, ids, Q_points) {
      fids <- factor(ids, levels = unique(ids))
      if (!is.list(Q_points))
        Q_points <- split(Q_points, row(Q_points))
      ind <- mapply(findInterval, Q_points, split(times, fids))
      ind[ind < 1] <- 1
      rownams_id <- split(row.names(data), fids)
      ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
      data[c(ind),]
    }
    
    newdata_pred <-
      right_rows(newdata, newdata[[timeVar]], id, times_to_pred)
    
    
    if(keeplevels) {
      if(length(setdiff(is.fact, idVar)) > 0) {
        newdata_pred <- newdata_pred %>% droplevels
        newdata_pred <- newdata_pred %>% 
          dplyr::select(-dplyr::all_of(setdiff(is.fact, idVar)))
        
        newdata_is.factx <- newdata %>% 
          dplyr::select(dplyr::all_of(is.fact))
        
        newdata_pred <- newdata_pred %>% 
          dplyr::left_join(., newdata_is.factx, by = idVar,
                    relationship = "many-to-many")
      }
    }
    
   
    
    newdata_pred[[timeVar]] <- unlist(times_to_pred)
    
    if(keeplevels) {
      newdata_pred <- newdata_pred %>% dplyr::select(dplyr::all_of(cnames))
    }
    if(setasdt) newdata_pred <- data.table::as.data.table(newdata_pred)
    if(asdf) out <- as.data.frame(newdata_pred) else out <- newdata_pred 
    # newdata_pred
    out
  }



#' An internal function to edit stancode for tripple logistic model
#' 
#' @param stancode A string character of stan code
#' 
#' @param set_positive_ordered A logical (default \code{TRUE}) to indicate
#'   whether to set \code{transformed parameters} to \code{positive_ordered} or
#'   \code{ordered}. If \code{TRUE}, \code{transformed parameters} are set to
#'   \code{positive_ordered} and to \code{ordered} if \code{FALSE}. if 
#'   \code{NULL}, set to \code{vector}.  
#' 
#' @param constraint A logical (default \code{TRUE}) to indicate whether to 
#' put constraints on the positive ordered vector space.
#' 
#' @param normalize A logical (default \code{TRUE}) to indicate whether to 
#' include the normalizing constant in the prior target density.
#' 
#' @keywords internal
#' @return A character string.
#' @noRd
#'
edit_scode_for_logistic3 <- function(stancode, 
                                     set_positive_ordered = TRUE,
                                     constraint = TRUE,
                                     normalize = TRUE) {
  
  setorder_d <- c(2, 3, 1)
  setorder_v <- c(3, 1, 2)
  setorder_t <- c(1, 2, 3)
  
  # Seems both set_positive_ordered and constraint should be TRUE
  
  true_name_p      <- 'parameters'
  true_name_tp     <- 'transformed parameters'
  true_name_td     <- 'transformed data'
  true_name_model  <- 'model'
  
  tempt_name_p  <- 'ppppppppp'
  tempt_name_tp <- 'tptptptpt'
  tempt_name_td <- 'tdtdtdtdt'
  
  clines_tp <- get_par_names_from_stancode(stancode,
                                           section = true_name_tp,
                                           semicolan = TRUE,
                                           full = TRUE)
  
  clines_p <- get_par_names_from_stancode(stancode,
                                          section = true_name_p,
                                          semicolan = TRUE,
                                          full = TRUE)
  
  clines_m <- get_par_names_from_stancode(stancode,
                                          section = true_name_model,
                                          semicolan = TRUE,
                                          full = TRUE)
  
  
  editedcode    <- stancode 
  editedcode    <- gsub(true_name_tp, tempt_name_tp, editedcode, fixed = T)
  editedcode    <- gsub(true_name_p,  tempt_name_p,  editedcode, fixed = T)
  editedcode    <- gsub(true_name_td,  tempt_name_td,  editedcode, fixed = T)
  
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
  
  
  clines_p2 <- c()
  for (il in clines_p) {
    il <- gsub(pattern = "//", replacement = "//", x = il, fixed = T)
    il <- gsub(pattern = "//[^\\\n]*", replacement = "", x = il)
    if(!grepl('^lprior', gsub_space(il)) & 
       !grepl('^reallprior', gsub_space(il)) & 
       !grepl('^-', gsub_space(il))) {
      if(!is_emptyx(il)) {
        clines_p2 <- c(clines_p2, il) 
      }
    }
  }
  
  clines_p <- clines_p2
  
  
  
  
  b_what_by_pair <- matrix(NA, 9, 3)
  K_name <- "K"
  move_to_tp <- add_move_to_tp <- c()
  for (clines_tpi in clines_p) {
    parameter_name_temp <- sub('.+](.+)', '\\1', clines_tpi)
    parameter_name_temp <- gsub(";", "", parameter_name_temp)
    parameter_name_temp <- gsub_space(parameter_name_temp)
    for (igr in 1:9) {
      if(grepl(paste0("^b", "_"), parameter_name_temp) & 
         grepl(paste0("_", letters[igr]), parameter_name_temp)) {
        move_to_tp <- c(move_to_tp, clines_tpi)
        if(grepl("<lower=", clines_tpi, fixed = T) |
           grepl("<upper=", clines_tpi, fixed = T)) {
          check_bounds_y_n <- 'y'
        } else {
          check_bounds_y_n <- 'n'
        }
        if(letters[igr] == 'a')  poparm <- paste0("raw_d[", setorder_d[1], "]")
        if(letters[igr] == 'd')  poparm <- paste0("raw_d[", setorder_d[2], "]")
        if(letters[igr] == 'g')  poparm <- paste0("raw_d[", setorder_d[3], "]")
        
        if(letters[igr] == 'b')  poparm <- paste0("raw_v[", setorder_v[1], "]")
        if(letters[igr] == 'e')  poparm <- paste0("raw_v[", setorder_v[2], "]")
        if(letters[igr] == 'h')  poparm <- paste0("raw_v[", setorder_v[3], "]")
        
        if(letters[igr] == 'c')  poparm <- paste0("raw_t[", setorder_t[1], "]")
        if(letters[igr] == 'f')  poparm <- paste0("raw_t[", setorder_t[2], "]")
        if(letters[igr] == 'i')  poparm <- paste0("raw_t[", setorder_t[3], "]")
        
        add_move_to_tp_k_dim <- paste0(K_name, "_", letters[igr])
        add_move_to_tp_ <- paste0(" = ", poparm, " + ", "rep_vector(0.0, ",
                                  add_move_to_tp_k_dim, ")")
        add_move_to_tp_2 <- gsub(";", paste0(add_move_to_tp_, ";"),  
                                 clines_tpi, fixed = T)
        
        add_move_to_tp <- c(add_move_to_tp, add_move_to_tp_2)
        parameter_name <- sub('.+](.+)', '\\1', clines_tpi)
        parameter_name <- gsub(";", "", parameter_name)
        parameter_name <- gsub_space(parameter_name)
        b_what_by_pair[igr , 1] <- parameter_name
        b_what_by_pair[igr , 2] <- poparm
        b_what_by_pair[igr , 3] <- check_bounds_y_n
      } 
    }
  } 
  

  
  b_what_it_c <- b_what_by_c <- c()
  for (igr in 1:nrow(b_what_by_pair)) {
    set_check_bounds_y_n_cnt <- 0
    for (igri in 1:2) {
      set_check_bounds_y_n_cnt <- set_check_bounds_y_n_cnt + 1
      b_what_it <- b_what_by_pair[igr , 1]
      b_what_by <- b_what_by_pair[igr , 2]
      check_bounds_y_n <- b_what_by_pair[igr , 3]
      if(check_bounds_y_n == "n") {
        b_what_by <- gsub("[", paste0("[", igri, ","), b_what_by, fixed = T)
        b_what_it <- paste0(b_what_it, "[", igri,"]")
        b_what_it_c <- c(b_what_it_c, b_what_it)
        b_what_by_c <- c(b_what_by_c, b_what_by)
      } else if(check_bounds_y_n == "y") {
        if(set_check_bounds_y_n_cnt == 1) {
          b_what_by <- gsub("[", paste0("[", ","), b_what_by, fixed = T)
          b_what_it <- b_what_it
          b_what_it_c <- c(b_what_it_c, b_what_it)
          b_what_by_c <- c(b_what_by_c, b_what_by)
        }
      } # else if(check_bounds_y_n
    }
  }
  
  
  
  set_b_a <- paste0("b_a = to_vector(raw_dx[", ",", setorder_d[1], "]);")
  set_b_d <- paste0("b_d = to_vector(raw_dx[", ",", setorder_d[2], "]);")
  set_b_g <- paste0("b_g = to_vector(raw_dx[", ",", setorder_d[3], "]);")
  
  set_b_b <- paste0("b_b = to_vector(raw_vx[", ",", setorder_v[1], "]);")
  set_b_e <- paste0("b_e = to_vector(raw_vx[", ",", setorder_v[2], "]);")
  set_b_h <- paste0("b_h = to_vector(raw_vx[", ",", setorder_v[3], "]);")
  
  set_b_c <- paste0("b_c = to_vector(raw_tx[", ",", setorder_t[1], "]);")
  set_b_f <- paste0("b_f = to_vector(raw_tx[", ",", setorder_t[2], "]);")
  set_b_i <- paste0("b_i = to_vector(raw_tx[", ",", setorder_t[3], "]);")
  
  set_b_elements <- paste(set_b_a, set_b_b, set_b_c,
                          set_b_d, set_b_e, set_b_f,
                          set_b_g, set_b_h, set_b_i,
                          sep = "\n   ")
  
  
  if(!is.null(set_positive_ordered)) {
    if(set_positive_ordered) {
      move_to_tp_add_ordered_positive_ordered <- 
        "array[Kedit] positive_ordered[Cedit] raw_dx;
   array[Kedit] positive_ordered[Cedit] raw_vx;
   array[Kedit] positive_ordered[Cedit] raw_tx;"
    }
    
    if(!set_positive_ordered) {
      move_to_tp_add_ordered_positive_ordered <- 
        "array[Kedit] ordered[Cedit] raw_dx;
   array[Kedit] ordered[Cedit] raw_vx;
   array[Kedit] ordered[Cedit] raw_tx;"
    }
  } else if(is.null(set_positive_ordered)) {
    move_to_tp_add_ordered_positive_ordered <- 
      "array[Kedit] vector[Cedit] raw_dx;
   array[Kedit] vector[Cedit] raw_vx;
   array[Kedit] vector[Cedit] raw_tx;"
  }
  
  
  
  if(constraint) {
    move_to_tp_add_constraint <- 
      " for(k in 1:Kedit) {
     raw_dx[k, ] = ordered_lb_ub_lp(raw_d[k,], min_d[k], max_d[k]);
     raw_vx[k, ] = ordered_lb_ub_lp(raw_v[k,], min_v[k], max_v[k]);
     raw_tx[k, ] = ordered_lb_ub_lp(raw_t[k,], min_t[k], max_t[k]);
    } "
  } else if(!constraint) {
    move_to_tp_add_constraint <- 
      " for(k in 1:Kedit) {
     raw_dx[k, ] = raw_d[k,];
     raw_vx[k, ] = raw_v[k,];
     raw_tx[k, ] = raw_t[k,];
    } "
  } 
  
  
  move_to_tp_add <- paste(move_to_tp_add_ordered_positive_ordered,
                          move_to_tp_add_constraint,
                          sep = "\n  ")
  
  
  move_to_tp_add <- paste(move_to_tp_add, 
                          set_b_elements, 
                          sep = "\n   ")
  
  
  move_to_tp2 <- c()
  for (move_to_tpi in move_to_tp) {
    move_to_tp2 <- c(move_to_tp2, paste0("   ", move_to_tpi))
  }
  
  tpcode <- paste0(paste(move_to_tp2, collapse = "\n"), 
                   "\n    ", 
                   move_to_tp_add)
  
  
  pcode <- 
    "array[Kedit] vector[Cedit] raw_d;
  array[Kedit] vector[Cedit] raw_v;
  array[Kedit] vector[Cedit] raw_t;
  "
  
  
  # parameter names - remove from the parameters block
  for (il in move_to_tp) {
    editedcode2 <- gsub(pattern = "//", replacement = "//", 
                        x = editedcode2, fixed = T)
    editedcode2 <- gsub(pattern = "//[^\\\n]*", replacement = "", 
                        x = editedcode2)
    editedcode2 <- gsub(paste0(il, ""), "", editedcode2, fixed = T)
    
  }
  
  
  
  # Remove empty lines
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
      zz_c <- c(zz_c, zz_in)
    }
  }
  editedcode2 <- paste(zz_c, collapse = '\n')
  
  
  p_block_syb_by <- paste0("", tempt_name_tp, " {")
  p_block_syb_it <- paste0(p_block_syb_by, "\n", tpcode)
  editedcode2 <- gsub(paste0("", p_block_syb_by), p_block_syb_it, 
                      editedcode2, fixed=T, perl=F)
  
  
  editedcode2 <- gsub(tempt_name_tp, true_name_tp, editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_p,  true_name_p,  editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_td,  true_name_td,  editedcode2, fixed = T)
  
  
  # Replace parameters in priors
  for (igr in 1:length(b_what_it_c)) {
    parameter_name <- b_what_it_c[igr]
    poparm         <- b_what_by_c[igr]
    editedcode2 <- gsub(paste0("(", parameter_name),
                        paste0("(", poparm),
                        editedcode2, fixed = T)
  }
  
  # Remove empty lines
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
      zz_c <- c(zz_c, zz_in)
    }
  }
  
  editedcode2 <- paste(zz_c, collapse = '\n')
  
  # https://github.com/stan-dev/math/issues/2959
  fcode <- 
    "vector ordered_lb_ub_lp (vector y, real lb, real ub) {
    int N = rows(y);
    vector[N] x;
    x[1] = lb + (ub - lb) * inv_logit(y[1]);
    target += log(ub - lb) + log_inv_logit(y[1]) + log1m_inv_logit(y[1]);
    for (i in 2:N) {
      x[i] = x[i - 1] + (ub - x[i - 1]) * inv_logit(y[i] - log(N+1-i));
      target += log(ub - x[i - 1]) + log_inv_logit(y[i] - log(N+1-i)) + 
        log1m_inv_logit(y[i] - log(N+1-i));
    }
    return x;
  }"
  
  
  
  out <- list(editedcode = editedcode2, 
              tpcode = tpcode, 
              pcode = pcode, 
              fcode = fcode)
  return(out)
}






#' Identify (and remove) outliers with abnormal growth velocity
#'
#' @description The function identifies and remove the putative outliers for y
#'   in data, Codes range from 0 (normal) to 8, where 4 and 6 are conventional
#'   outliers
#'
#' @inherit sitar::velout params description
#' 
#' @inherit sitar::zapvelout params
#' 
#' @param remove A logical (default \code{FALSE}) to indicate whether identified
#'   remove to be removed. When \code{FALSE}, outliers are set as \code{NA}. If
#'   \code{TRUE}, outliers set as \code{NA} are removed.
#'   
#' @param verbose A logical (default \code{FALSE}) to show frequency of
#'   observations with different codes (see [sitar::velout()]), and the number
#'   of observations (outliers) removed when \code{remove=TRUE}.
#'
#' @return A data frame with outliers removed when \code{remove=TRUE}
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
#' @examples
#' 
#' \dontrun{
#' data <- berkeley_mdata
#' 
#' outliers(x = "age", y = "height", id = "id", data = data, remove=TRUE)
#' }
#' 
outliers <-
  function (x,
            y,
            id,
            data,
            icode = c(4:6),
            lag = 1,
            velpower = 0.5,
            limit = 5,
            linearise = FALSE,
            remove = FALSE,
            verbose = TRUE) {
    
    mcall <- match.call()
    if(is.symbol(x)) xx_ <- deparse(substitute(x)) else xx_ <- x
    if(is.symbol(y)) yy_ <- deparse(substitute(y)) else yy_ <- y
    if(is.symbol(id)) idid_ <- deparse(substitute(id)) else idid_ <- id
    
    data <- data %>% dplyr::mutate(order = dplyr::row_number())
    data <- data %>% dplyr::arrange(idid_, xx_)
    
    dc <- data %>% dplyr::select(!!as.symbol(xx_), 
                                 !!as.symbol(yy_),
                                 !!as.symbol(idid_))
    
    colnames_data_ex <- colnames(data)
    colnames_dc_ex <- colnames(dc)
    data_ex <- data %>% dplyr::select(-colnames(dc))
    
    nrow <- nrow(data)
    dc <- na.omit(cbind(dc, count = 1:nrow))
    dc <- dc[order(dc[, 3], dc[, 1]),]
    if (linearise) {
      spline.lm <- loess(dc[, 2] ~ dc[, 1])
      dc[, 2] <- residuals(spline.lm)
    }
    dt1 <- diff(dc[, 1], lag = lag)
    vel1 <- diff(dc[, 2], lag = lag) / dt1 ^ velpower
    dt2 <- diff(dc[, 1], lag = lag * 2)
    vel2 <- diff(dc[, 2], lag = lag * 2) / dt2 ^ velpower
    dt1 <- dt1 == 0
    idlev <- as.numeric(dc[, 3])
    dt1[diff(idlev, lag = lag) != 0] <- FALSE
    vel1[diff(idlev, lag = lag) != 0] <- NA
    vel2[diff(idlev, lag = lag * 2) != 0] <- NA
    vel1 <- trunc(vel1 / mad(vel1, na.rm = TRUE))
    vel2 <- trunc(vel2 / mad(vel2, na.rm = TRUE))
    vel3 <- c(rep(NA, lag), vel2, rep(NA, lag))
    vel2 <- c(vel1, rep(NA, lag))
    vel1 <- c(rep(NA, lag), vel1)
    code <- (as.numeric(abs(vel1) >= limit) + as.numeric(abs(vel2) >=
                                                           limit)) * 2 + 
      as.numeric(abs(vel3) >= limit)
    
    dt2 <- c(dt1, rep(FALSE, lag))
    dt1 <- c(rep(FALSE, lag), dt1)
    code[dt2 | dt1] <- 8
    code[dt2 & !dt1] <- 7
    t <- is.na(vel3) & !(dt1 | dt2)
    code[t] <- (as.numeric(!is.na(vel1[t]) & abs(vel1[t]) >=
                             limit) + as.numeric(!is.na(vel2[t]) &
                                                   abs(vel2[t]) >=
                                                   limit)) * 6
    dc <- cbind(dc[, c(3, 1, 2, 4)], code, vel1, vel2, vel3)
    mat <- as.data.frame(matrix(
      nrow = nrow,
      ncol = dim(dc)[2],
      dimnames = list(row.names(data), dimnames(dc)[[2]])
    ))
    attr(mat, "data") <- deparse(mcall$data)
    mat[dc$count,] <- dc
    if (is.factor(dc[, 1])) {
      mat[, 1] <- as.factor(mat[, 1])
      levels(mat[, 1]) <- levels(dc[, 1])
    }
    mat$count <- NULL
    mat$code <- factor(mat$code)
    if (verbose) {
      cat("code frequencies\n")
      print(summary(mat$code))
    }
    
    mat <- cbind(mat, data_ex)
    mat <- mat %>% dplyr::relocate(dplyr::all_of(colnames_data_ex))
    
    if (remove) {
      zap <- mat$code %in% icode
      if (verbose) {
        cat(sum(zap), yy_, "values set missing\n")
      }
      mat[mat$code %in% icode, yy_] <- NA
      mat <- mat %>% dplyr::select(dplyr::all_of(colnames_data_ex))
      mat <- mat %>% tidyr::drop_na() %>% droplevels()
    }
    mat <- mat %>% dplyr::arrange(order)
    mat <- mat %>% dplyr::select(-order)
    mat
  }



#' An internal function to set default initials for model specific parameters
#'
#' @param cargs A character string specifying the model fitted. 
#' 
#' @param fargs A character string specifying the initials. 
#' 
#' @param dargs A character string specifying the parameter class. Options
#' are \code{'b'}, \code{'sd'} and \code{'cor'}. Default \code{NULL} indicates 
#' that class name in infered automatically. 
#' 
#' @param verbose A logical (default \code{FALSE}) to indicate whetehr to print
#'  \code{warnings} and \code{messages} during the function evaluation.
#'
#' @return A list.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
evaluate_call_args <- function(cargs = NULL,
                            fargs = NULL,
                            dargs = NULL,
                            verbose = FALSE) {
  for (fargsi in names(dargs)) {
    if(is.null(cargs[[fargsi]])) cargs[[fargsi]] <- fargs[[fargsi]]
  }
  for (fargsi in names(fargs)) {
    if(is.null(cargs[[fargsi]])) cargs[[fargsi]] <- fargs[[fargsi]]
  }
 
  cargs
}



#' An internal function to perform checks when calling post-processing functions
#'
#' @description The \code{post_processing_args_sanitize} perform essential 
#'   checks for the arguments passed on to the \code{brms} post-processing
#'   functions.
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
#' @param dots A list passing the \code{...} arguments. Default \code{NULL}.
#' 
#' @param misc A vector of character strings specifying the miscellanous
#'   arguments to be checked. Default \code{NULL}.
#'
#' @param envir Indicator to set the environment of function evaluation. The
#'   default is \code{parent.frame}.
#'   
#' @param verbose A logical (default \code{FALSE}) to indicate whetehr to print
#'  \code{warnings} and \code{messages} during the function evaluation.
#'
#' @return A list with the filtered necessary arguments required for the 
#'   successgull execition of the post-processing function
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
post_processing_args_sanitize <- function(model, 
                                          xcall, 
                                          resp = NULL, 
                                          deriv = NULL,
                                          dots = NULL,
                                          misc = NULL,
                                          envir = NULL, 
                                          verbose = FALSE) {
  
  if(is.null(envir)) envir <- parent.frame()
  if(is.null(deriv)) deriv <- 0
  
  if(!'bgmfit' %in% class(model)) {
    stop("The class of model object should be 'bgmfit' ")
  }
  
  allargs <- c(as.list(xcall), dots)
  
  excall_ <- c("plot_ppc", "loo_validation")
  excall_ <- c(excall_, paste0(excall_, ".", "bgmfit"))
  
  checkwhat_all <- c()
  if (strsplit(deparse((xcall[1])), "\\.")[[1]][1] %in% excall_) {
    # Check deriv
    checkwhat <- ""
    checkwhat <- 'deriv'
    if(!is.null(allargs[[checkwhat]])) {
      checkwhat_all <- c(checkwhat_all, checkwhat)
      if(verbose) {
        message(
          "\nargument 'deriv' is not allowed for the ",
          " post-processing function",  " '",
          strsplit(deparse((xcall[1])), "\\.")[[1]][1], "'",
          "\n ",
          "Therefore, it is set to i.e., deriv = NULL"
        )
      }
    } # if(!is.null(allargs$idata_method)) {
    
    # Check idata_method
    checkwhat <- ""
    checkwhat <- 'idata_method'
    if(!is.null(allargs[[checkwhat]])) {
      checkwhat_all <- c(checkwhat_all, checkwhat)
      if(verbose) {
        message(
          "\nargument 'idata_method' is not allowed for the ",
          " post-processing function",  " '",
          strsplit(deparse((xcall[1])), "\\.")[[1]][1], "'",
          "\n ",
          "Therefore, it is set to i.e., idata_method = NULL"
        )
      }
    } # if(!is.null(allargs$idata_method)) {
    
    # Check re_formula
    checkwhat <- ""
    checkwhat <- 're_formula'
    if(!is.null(allargs[[checkwhat]])) {
      checkwhat_all <- c(checkwhat_all, checkwhat)
      if(verbose) {
        message(
          "\nargument 'idata_method' is not allowed for the ",
          " post-processing function",  " '",
          strsplit(deparse((xcall[1])), "\\.")[[1]][1], "'",
          "\n ",
          "Therefore, it is set to i.e., idata_method = NULL"
        )
      }
    } # if(!is.null(allargs$idata_method)) {
    
    
  } # if (strsplit(deparse((xcall[1])), "\\.")[[1]][1] %in% excall_) {
  
  
  if(length(checkwhat_all) == 0) checkwhat_all <- NULL
  checkwhat_all <- c(checkwhat_all, misc)
  
  if(!is.null(checkwhat_all)) {
    allargs[which(names(allargs)%in%checkwhat_all)] <- NULL
  }
  
  allargs <- allargs[-1]
  names(allargs) <- gsub("model", "object", names(allargs))
 
  return(allargs)
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
                                   verbose = FALSE) {
  
  if(is.null(envir)) envir <- parent.frame()
  if(is.null(deriv)) deriv <- 0
  
  if(!'bgmfit' %in% class(model)) {
    stop("The class of model object should be 'bgmfit' ")
  }
  
  excall_ <- c("plot_ppc", "loo_validation")
  
  check_it <- strsplit(deparse((xcall[1])), "\\.")[[1]][1] 
  check_it <- gsub("\"",  "", check_it)
  
  if (check_it %in% excall_) {
    if(is.null(as.list(xcall)[['deriv']])) deriv <- ''
    if (!is.null(as.list(xcall)[['deriv']])) {
      deriv <- ''
      if(verbose) {
        message(
          "\nargument 'deriv' is not allowed for the ",
          " post-processing function",  " '",
          check_it, "'",
          "\n ",
          "Therefore, it is set to missing i.e., deriv = ''"
        )
      }
    } # if(!is.null(chcallls$idata_method)) {
  }
  
  if (model$model_info$nys == 1 & !is.null(resp)) {
    stop(
      "You have fit a univariate model",
      " but set resp option as: ",
      resp,
      ".",
      "\n ",
      " For univariate model, the resp option should be NULL",
      "\n ",
      " (i.e., resp = NULL)"
    )
  }
  if (model$model_info$nys > 1 & is.null(resp)) {
    if (!is.na(model$model_info$univariate_by)) {
      stop(
        "You have fit a univariate_by model for ",
        model$model_info$univariate_by,
        "\n ",
        " but did not correctly specified the 'resp' option",
        " (which is NULL at present).",
        "\n ",
        " The response options are: ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
    if (model$model_info$multivariate) {
      stop(
        "You have fit a multivariate model ",
        "\n ",
        " but dit not set the the resp options correctly",
        " (which is NULL at present).",
        "\n ",
        " The response options are: ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
  }
  if (is.null(resp)) {
    resp_ <- resp
  } else if (!is.null(resp)) {
    resp_ <- paste0(resp, "_")
  }
  
  
  # assign expose default funs 
  if(model$model_info[['expose_method']] == 'R') {
    assign(paste0(resp_, 
                  model$model_info[['namesexefuns']], 
                  '0'), 
           model$model_info$exefuns[[paste0(resp_, 
                                            model$model_info[['namesexefuns']], 
                                            '0')]], envir = envir)
    
    assign(paste0(resp_, 'getX'), 
           model$model_info$exefuns[[paste0(resp_, 'getX')]], 
           envir = envir)
    
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
  
  return(out)
}




#' An internal function to set up exposed functions and their environment
#'
#' @param o A logical (default \code{FALSE}) to indicate whether to
#' return the object as a character string.
#' @inherit growthparameters.bgmfit params
#' @param ... Additional arguments. Currently ignored.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'

setupfuns <- function(model,
                      resp = NULL,
                      o = NULL,
                      oall = NULL,
                      usesavedfuns = NULL,
                      deriv = NULL,
                      envir = NULL,
                      deriv_model = NULL,
                      verbose = FALSE,
                      ...) {
  
  if(is.null(envir)) {
    envir <- parent.frame()
  }
  
  if (is.null(resp)) {
    resp_ <- resp
  } else if (!is.null(resp)) {
    resp_ <- paste0(resp, "_")
  }
  
  if(is.null(model$xcall)) {
    xcall <- strsplit( deparse(sys.calls()[[sys.nframe()-1]]) , "\\(")[[1]][1]
  } else {
    xcall <- model$xcall
  }
  
  
  excall_ <- c("plot_ppc", "loo_validation")
  
  check_it <- strsplit(deparse((xcall[1])), "\\.")[[1]][1] 
  check_it <- gsub("\"",  "", check_it)
  
  if (check_it %in% excall_) {
    if(is.null(as.list(xcall)[['deriv']])) deriv <- ''
    if (!is.null(as.list(xcall)[['deriv']])) {
      deriv <- ''
      if(verbose) {
        message(
          "\nargument 'deriv' is not allowed for the ",
          " post-processing function",  " '",
          check_it, "'",
          "\n ",
          "Therefore, it is set to missing i.e., deriv = ''"
        )
      }
    } # if(!is.null(chcallls$idata_method)) {
  }
  
  if(!usesavedfuns) {
    if(is.null(check_if_functions_exists(model, o, xcall))) {
      return(invisible(NULL))
    }
  }
  
  if(usesavedfuns) {
    if(is.null(check_if_functions_exists(model, o, xcall,
                                         verbose = F))) {
      envir <- envir
    } else {
      #  envir <- getEnv(o[[1]], geteval = TRUE)
    }
    # envir <- getEnv(o[[1]], geteval = TRUE)
    
    oall <- model$model_info[['exefuns']]
    oalli_c <- names(oall)
    for (oalli in oalli_c) {
      assign(oalli, oall[[oalli]], envir = envir)
    }
  }
  
  if(!is.null(deriv)) {
    if(deriv == 0) {
      assignfun <- paste0(model$model_info[['namesexefuns']], deriv)
      assignfun <- paste0(resp_, assignfun)
      assign(o[[1]], model$model_info[['exefuns']][[assignfun]], envir = envir)
    } else if(deriv > 0) {
      if(deriv_model) {
        assignfun <- paste0(model$model_info[['namesexefuns']], deriv)
        assignfun <- paste0(resp_, assignfun)
      } else if(!deriv_model) {
        assignfun <- paste0(model$model_info[['namesexefuns']], '0')
        assignfun <- paste0(resp_, assignfun)
      }
      assign(o[[1]], model$model_info[['exefuns']][[assignfun]], envir = envir)
    }
  }
  
  if(is.null(deriv)) {
    assignfun <- paste0(model$model_info[['namesexefuns']], "")
    assignfun <- paste0(resp_, assignfun)
    assign(o[[1]], model$model_info[['exefuns']][[assignfun]], envir = envir)
  }
  return(envir)
}







#' An internal function to set default prior for model specific parameters
#'
#' @param select_model A character string specifying the model fitted. 
#' 
#' @param prior A character string specifying the prior. 
#' 
#' @param class A character string specifying the parameter class. Options
#' are \code{'b'}, \code{'sd'} and \code{'cor'}. Default \code{NULL} indicates 
#' that class name in infered automatically. 
#' 
#' @param parameter A character string specifying the parameter name. Options
#'  are \code{'a'}, \code{'b'}, \code{'c'}, \code{'d'}, \code{'e'}, \code{'f'},
#'  \code{'g'}, \code{'h'}, and \code{'i'}. Default \code{NULL} indicates that
#'  parameter name in infered automatically. 
#'
#' @return A character string.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
set_default_priors <- function(select_model,
                               prior,
                               class = NULL,
                               parameter = NULL) {
  if (is.null(parameter)) {
    parameter <- strsplit(deparse(substitute(prior)), "_")[[1]][1]
  }
  if (is.null(class)) {
    get_suffix <- strsplit(deparse(substitute(prior)), "_")[[1]]
    get_suffix <- get_suffix[length(get_suffix)]
    if (grepl("^beta", get_suffix))
      class <- 'b'
    if (grepl("^sd", get_suffix))
      class <- 'sd'
  }
  
  
  ##############################################################
  # class b
  ##############################################################
  
  # parameter a class b
  if (parameter == 'a' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(ymin, ysdxmin, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter b class b
  if (parameter == 'b' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(ymaxs, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0.1, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(ymaxs, ysd, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.5, 0.5, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter c class b
  if (parameter == 'c' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0.1, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(5, 3, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0.1, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0.1, 0.1, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter d class b
  if (parameter == 'd' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(ymeanxmid, ysdxmid, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter e class b
  if (parameter == 'e' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(7, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0.2, 0.1, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter f class b
  if (parameter == 'f' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(13, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(5, 3, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter g class b
  if (parameter == 'g' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- 
            "normal(ymeanxmidxmaxdiff, ysdxmidxmaxdiff, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter h class b
  if (parameter == 'h' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.5, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter i class b
  if (parameter == 'i' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(14, 2, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  
  
  
  ##############################################################
  # class sd
  ##############################################################
  
  # parameter a class sd
  if (parameter == 'a' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, ysdxmin, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter b class sd
  if (parameter == 'b' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, ysd, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter c class sd
  if (parameter == 'c' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 3, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter d class sd
  if (parameter == 'd' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, ysdxmid, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter e class sd
  if (parameter == 'e' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 0.15, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter f class sd
  if (parameter == 'f' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter g class sd
  if (parameter == 'g' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, ysdxmidxmaxdiff, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter h class sd
  if (parameter == 'h' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter i class sd
  if (parameter == 'i' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  return(prior_out)
}







#' An internal function to set default initials for model specific parameters
#'
#' @param select_model A character string specifying the model fitted. 
#' 
#' @param init A character string specifying the initials. 
#' 
#' @param class A character string specifying the parameter class. Options
#' are \code{'b'}, \code{'sd'} and \code{'cor'}. Default \code{NULL} indicates 
#' that class name in infered automatically. 
#' 
#' @param parameter A character string specifying the parameter name. Options
#'  are \code{'a'}, \code{'b'}, \code{'c'}, \code{'d'}, \code{'e'}, \code{'f'},
#'  \code{'g'}, \code{'h'}, and \code{'i'}. Default \code{NULL} indicates that
#'  parameter name in infered automatically. 
#'
#' @return A character string.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
set_default_inits <- function(select_model,
                              init,
                              class = NULL,
                              parameter = NULL) {
  if (is.null(parameter)) {
    parameter <- strsplit(deparse(substitute(init)), "_")[[1]][1]
  }
  if (is.null(class)) {
    get_suffix <- strsplit(deparse(substitute(init)), "_")[[1]]
    get_suffix <- get_suffix[length(get_suffix)]
    if (grepl("^beta", get_suffix))
      class <- 'b'
    if (grepl("^sd", get_suffix))
      class <- 'sd'
  }
  
  
  ##############################################################
  # class b
  ##############################################################
  
  # parameter a class b
  if (parameter == 'a' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "ymean"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "ymean"
      } else if (grepl('^pb', select_model)) {
        init_out <- "ymax"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "ymax"
        }
        if (select_model == 'logistic2') {
          init_out <- "ymax"
        }
        if (select_model == 'logistic3') {
          init_out <- "ymin"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter b class b
  if (parameter == 'b' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- "ymaxs"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 0.1
        }
        if (select_model == 'logistic2') {
          init_out <- "ymaxs"
        }
        if (select_model == 'logistic3') {
          init_out <- 1.5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter c class b
  if (parameter == 'c' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- 0.1
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 5.0
        }
        if (select_model == 'logistic2') {
          init_out <- 0.1
        }
        if (select_model == 'logistic3') {
          init_out <- 0.1
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter d class b
  if (parameter == 'd' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 1.0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- 1.2
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 1.0
        }
        if (select_model == 'logistic2') {
          init_out <- 1.2
        }
        if (select_model == 'logistic3') {
          init_out <- "ymeanxmid"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter e class b
  if (parameter == 'e' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        if(select_model == 'pb1') init_out <- 13
        if(select_model == 'pb2') init_out <- 13
        if(select_model == 'pb3') init_out <- 13
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- 7
        }
        if (select_model == 'logistic3') {
          init_out <- 0.15
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter f class b
  if (parameter == 'f' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        if(select_model == 'pb1') init_out <- NULL
        if(select_model == 'pb2') init_out <- 2
        if(select_model == 'pb3') init_out <- 1
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- 13
        }
        if (select_model == 'logistic3') {
          init_out <- 5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter g class b
  if (parameter == 'g' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 
            "ymeanxmidxmaxdiff"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter h class b
  if (parameter == 'h' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 1.5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter i class b
  if (parameter == 'i' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 13
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  
  
  
  ##############################################################
  # class sd
  ##############################################################
  
  # parameter a class sd
  if (parameter == 'a' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmin, autoscale = 2.5)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter b class sd
  if (parameter == 'b' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, ysd, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter c class sd
  if (parameter == 'c' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 3, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter d class sd
  if (parameter == 'd' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmid, autoscale = 2.5)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter e class sd
  if (parameter == 'e' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 0.15, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter f class sd
  if (parameter == 'f' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter g class sd
  if (parameter == 'g' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmidxmaxdiff, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter h class sd
  if (parameter == 'h' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter i class sd
  if (parameter == 'i' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  return(init_out)
}






#' An internal function to set variance covariance initials values to zero
#'
#' @param xscode A character string to specify the stancode.
#'
#' @param xsdata A character string to specify the standata.
#'
#' @param full A logical (default \code{TRUE}) to indicate whether full names
#'   should be extracted.
#'
#' @param what A character string to specify the variance covariance parameter.
#'   Default \code{L} indicating the correlation parameters. Other options are
#'   \code{sd} and \code{z}.
#'   
#' @param parameterization A character string to specify the 
#' the parameterization, CP or NCP.
#'
#' @param sd_value A numeric value to set intials for standard deviation
#'   parameters. Default \code{1} which is translated to zero initial i.e.,
#'   \code{log(1) = 0}.
#'
#' @param z_value A numeric value (default \code{0}) to set initials for
#'   \code{z} parameter which is part of the non centered parameterisation
#'   implemented in the [brms::brm()].
#'
#' @param L_value A numeric value (default \code{0}) to set initials for
#'   correlation parameter, \code{L}.
#'   
#' @return A prior object.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
set_init_gr_effects <- function(xscode,
                                xsdata,
                                full = TRUE,
                                what = 'L',
                                parameterization = 'ncp',
                                sd_value = 1,
                                z_value = 0,
                                r_value = 0,
                                L_value = 0) {
  xscode <-
    get_par_names_from_stancode(xscode, full = full, what = what)
  
  sdi_c <- c()
  parm_c <- c()
  sd_init_list_c_ <- list()
  getindexl <- 0
  for (sdi in xscode) {
    fullxscode_i <- sdi
    getparm <- tail(strsplit(sdi, split = " ")[[1]], 1)
    parm_c <- getparm # c(parm_c, getparm)
    sdi <- gsub("[", "(", sdi, fixed = T)
    sdi <- gsub("]", ")", sdi, fixed = T)
    str_d <-
      regmatches(sdi, gregexpr("(?=\\().*?(?<=\\))", sdi, perl = T))[[1]]
    # this to remove dim after par name
    parm_cm <- parm_c
    if (grepl("[", parm_cm, fixed = T)) {
      parm_cm <- gsub("[", "(", parm_cm, fixed = T)
      parm_cm <- gsub("]", ")", parm_cm, fixed = T)
      parm_cm_ <-
        regmatches(parm_cm,
                   gregexpr("(?=\\().*?(?<=\\))", parm_cm, perl = T))[[1]]
      parm_cm2 <- gsub(parm_cm_, "", parm_cm, fixed = T)
      parm_c <- parm_cm2
    } else if (!grepl("[", parm_cm, fixed = T)) {
      parm_c <- parm_c
    }
    parm_c <- parm_c
    # sdi_c <- c(sdi_c, str_d)
    sdi <- str_d
    sdi <- gsub("(", "", sdi, fixed = T)
    sdi <- gsub(")", "", sdi, fixed = T)
    str_d <- sdi
    #
    sdi <- gsub("(", "", sdi, fixed = T)
    sdi <- gsub(")", "", sdi, fixed = T)
    str_d <- gsub("[[:space:]]", "", sdi)
    str_d_ <- strsplit(str_d, ",") %>% unlist()
    
    # for student_nu distribution parameter
    if (grepl("sd_nu", parm_c, fixed = T)) {
      set_value <- 3
    }
    
    if (grepl("sd", parm_c, fixed = T)) {
      set_value <- sd_value
    }
    if (grepl("z", parm_c, fixed = T)) {
      set_value <- z_value
    }
    if (grepl("L", parm_c, fixed = T)) {
      set_value <- L_value
    }
    
    # to exclude student_nu distribution parameter
    if (grepl("sd", parm_c, fixed = T)) {
      # to exclude student_nu distribution parameter
      if (!grepl("sd_nu", parm_c, fixed = T)) {
        if (length(str_d_) == 1) {
          dim1 <- xsdata[[str_d_[1]]]
          out <- rep(set_value, dim1)
        } else if (length(str_d_) == 2) {
          dim1 <- xsdata[[str_d_[1]]]
          dim2 <- xsdata[[str_d_[2]]]
          out <- matrix(set_value, dim1, dim2)
        }
        if (is.vector(out)) {
          out <- array(out, dim = length(out))
        }
      } # if(!grepl("sd_nu", parm_c, fixed = T)) {
    } # if(grepl("sd", parm_c, fixed = T)) {
    
    if (grepl("z", parm_c, fixed = T)) {
      if (length(str_d_) == 1) {
        dim1 <- xsdata[[str_d_[1]]]
        out <- rep(set_value, dim1)
      } else if (length(str_d_) == 2) {
        dim1 <- xsdata[[str_d_[1]]]
        dim2 <- xsdata[[str_d_[2]]]
        out <- matrix(set_value, dim1, dim2)
      }
      
      if (is.vector(out)) {
        out <- array(out, dim = length(out))
      }
      
      if (ncol(out) == 1) {
        out <- t(out)
      }
    } # if(grepl("z", parm_c, fixed = T)) {
    
    
    if(parameterization == 'cp') {
      if (grepl("r", parm_c, fixed = T)) {
        if (length(str_d_) == 1) {
          dim1 <- xsdata[[str_d_[1]]]
          out <- rep(set_value, dim1)
        } else if (length(str_d_) == 2) {
          dim1 <- xsdata[[str_d_[1]]]
          dim2 <- xsdata[[str_d_[2]]]
          out <- matrix(set_value, dim1, dim2)
        }
        
        if (is.vector(out)) {
          out <- array(out, dim = length(out))
        }
        
        if (ncol(out) == 1) {
          out <- t(out)
        }
      } # if(grepl("r", parm_c, fixed = T)) {
    } # if(parameterization == 'cp') {
    
   
    
    if (grepl("L", parm_c, fixed = T)) {
      if (length(str_d_) == 1) {
        dim1 <- xsdata[[str_d_[1]]]
        out <- matrix(set_value, dim1, dim1)
        diag(out) <- 1
      }
      if (length(str_d_) == 2) {
        dim1 <- xsdata[[str_d_[1]]] 
        dim2 <- xsdata[[str_d_[2]]] 
        # outxx <- array(0, dim = c(dim2, dim1, dim1))
        outxx <- array(0, dim = c(dim1, dim2, dim2)) # 13 12 23
        for (i in 1:dim(outxx)[2]) {
          outxx[, i, i] <- 1
        }
        out <- outxx
      }
    } # if(grepl("L", parm_c, fixed = T)) {
    sd_init_list_c_[[parm_c]] <- out
  }
  sd_init_list_c_
}




#' An internal function to get initials from pathfinder algorithm
#'
#' @param pthf A pathfinder draw object
#' @param ndraws An integer specifying the number of paths
#' @param init_structure A list of initial to set the appropriate dimensions of 
#'  inits returned from the pathfinder
#' @param variables A character vector specifying the variable names
#' @param model An object of class \code{bgmfit}
#' @param compile_init_model_methods A logical to indicate whether to compile
#' \code{init_model_methods()}
#' @param verbose A logical
#' @return A named list.
#' @keywords internal
#' @noRd
#'
get_pathfinder_init <- function(pthf = NULL, 
                                ndraws = 1, 
                                init_structure = NULL,
                                variables = NULL,
                                model = NULL,
                                compile_init_model_methods = NULL,
                                verbose = FALSE) {
  
  if(is.null(pthf) & is.null(model)) 
    stop('Specify at least pthf or model')
  if(!is.null(pthf) & !is.null(model)) 
    stop('Specify either least pthf or model')
  
  if(!is.null(model)) {
    mod <- attr(model$fit, "CmdStanModel")
    dat <- brms::standata(model)
    ini <- model$stan_args$init
    threads <- model$threads$threads
    if(!is.null(threads)) {
      pth1 <- mod$pathfinder(data = dat, init = ini,  num_threads = threads)
    } else if(!is.null(threads)) {
      pth1 <- mod$pathfinder(data = dat, init = ini)
    }
    if(is.null(compile_init_model_methods)) compile_init_model_methods <- TRUE
  } else {
    if(is.null(compile_init_model_methods)) compile_init_model_methods <- FALSE
    pthf <- pthf
  }
  
  if(compile_init_model_methods) {
    pth1$init_model_methods()
  }
  
  lp__ <- NULL;
  lp_approx__ <- NULL;
  lw <- NULL;
  .draw <- NULL;
  lw.x <- NULL;
  lp__ <- NULL;
  
  
  as_inits <- function(draws, variable=NULL, ndraws=ndraws) {
    ndraws <- min(posterior::ndraws(draws),ndraws)
    if (is.null(draws)) {variable = variables(draws)}
    draws <- draws  %>%  posterior::as_draws_matrix()
    inits <- lapply(1:ndraws,
                    function(drawid) {
                      sapply(variable,
                             function(var) {
                               as.numeric(posterior::subset_draws(draws, 
                                                                  variable=var, 
                                                                  draw=drawid))
                             })
                    })
    if (ndraws==1) { inits[[1]] } else { inits }
  }
  
  if (is.null(variables)) {
    # set variable names to be list of parameter names
    variables <- names(pthf$variable_skeleton(transformed_parameters = FALSE,
                                              generated_quantities = FALSE))
  }
  draws <- pthf$draws(format="df")
  draws <- draws %>% 
    posterior::mutate_variables(lw = lp__ - lp_approx__)
  ndist <- dplyr::n_distinct(posterior::extract_variable(draws,"lw"))
  if (ndist < ndraws) {
    stop(paste0("Not enough distinct draws (", ndist, ") to create inits."))
  }
  if (ndist < 0.95*posterior::ndraws(draws)) {
    # Resampling has been done in Stan, compute weights for distinct draws
    #these are now non Pareto smoothed as we have lost the original information
    draws <- draws %>% 
      dplyr::group_by(lw) %>% 
      dplyr::summarise(.draw=min(.draw)) %>% 
      dplyr::left_join(draws, by = ".draw") %>% 
      posterior::as_draws_df() %>% 
      posterior::mutate_variables(lw = lw.x,
                                  w = exp(lw-max(lw)))
  } else {
    # Resampling was not done in Stan, compute Pareto smoothed weights
    draws <- draws %>% 
      posterior::mutate_variables(w=posterior::pareto_smooth(exp(lw-max(lw)), 
                                                             tail="right"))
  }
  
  out <- 
    draws %>% 
    posterior::weight_draws(weights=posterior::extract_variable(draws,"w"), 
                            log=FALSE) %>% 
    posterior::resample_draws(ndraws=ndraws, method = "simple_no_replace") %>% 
    as_inits(variable=variables, ndraws=ndraws)
  
  
  if(!is.null(init_structure)) {
    path_inits <- out # inits_pathfinder
    init_str_x <- init_structure # fit_m$stan_args$init[[1]]
    for (stri in names(init_str_x)) {
      if(is.array( init_str_x[[stri]] )) {
        if(!is.null(path_inits[[stri]])) {
          path_inits[[stri]] <- array(path_inits[[stri]], dim = dim(init_str_x[[stri]]) )
        }
      } else if(is.vector( init_str_x[[stri]] )) {
      } else if(is.numeric( init_str_x[[stri]] )) {
      }
      out <- path_inits
    }
  } # if(!is.null(init_structure)) {
  
  return(out)
}




# check if a str obj is actually numeric
# @description check if a str obj is actually numeric
# https://stackoverflow.com/questions/13638377/test-for-numeric-elements-in-a-character-string
# is.numeric.like -> changed to check_is_numeric_like -> is called in bsitar.R only
# Thid because of notes in rmd check which says 
# Mismatches for apparent methods not registered
# This perhaps because is.numeric.like sound like is.numeric
#' @param x a str vector, or a factor of str vector, or numeric vector. x will
#'   be coerced and trimws.
#' @param na.strings case sensitive strings that will be treated to NA.
#' @param naAsTrue whether NA (including actual NA and na.strings) will be
#'   treated as numeric like
#' @return a logical vector (vectorized).
#' @note Using regular expression
#' \cr TRUE for any actual numeric c(3,4,5,9.9) or c("-3","+4.4",
#' "-42","4L","9L",   "1.36e4","1.36E4",    NA, "NA", "","NaN", NaN):
#' \cr positive or negative numbers with no more than one decimal c("-3","+4.4")
#' OR
#' \cr positive or negative integers (e.g., c("-42","4L","39L")) OR
#' \cr positive or negative numbers in scientific notation c("1.36e4","1.36E4")
#' \cr NA, or na.strings
#' @keywords internal
#' @noRd
#'
check_is_numeric_like <- function(x, 
                            naAsTrue = TRUE, 
                            na.strings = 
                              c('','.','NA','na','N/A','n/a','NaN','nan')
                            ){
  x = trimws(x,'both')
  x[x %in% na.strings] = NA
  # https://stackoverflow.com/a/21154566/2292993
  result = grepl("^[\\-\\+]?[0-9]+[\\.]?[0-9]*$|^[\\-\\+]?[0-9]+[L]?$|^[\\-\\+]?[0-9]+[\\.]?[0-9]*[eE][0-9]+$",x,perl=TRUE)
  if (naAsTrue) result = result | is.na(x)
  return((result))
}




#' Title
#'
#' @param data data frame
#' @param id id
#' @param outcome outcome
#' @param time time variable such as age
#' @param timeval the value of time beyound which changes made
#' @param nset number of last time point to be flattened
#' @param inc increment for last n time point.
#'
#' @keywords internal
#' @noRd
#' 
flattten_last_time <-
  function (data,
            id,
            outcome,
            time,
            timeval = NULL,
            nset = 1,
            inc = NULL) {
    occtemp <- NULL;
    temdata <-
      data %>% dplyr::group_by_at(id) %>% 
      dplyr::mutate(`:=`("occtemp",
                         dplyr::row_number())) %>% 
      dplyr::mutate(`:=`("nocctemp",
                         max(.data[["occtemp"]])))
    setseq <- seq(1, nset, 1) - 1
    inc <- rev(inc)
    if (is.null(inc)) {
      inc <- 0
      inc <- rep(inc, length(nset))
    }
    else if (length(inc) == 1) {
      inc <- rep(inc, nset)
    }
    else if (length(inc) != nset) {
      stop("lenhth of 'inc' must be either 1 or same as the 'nset'")
    }
    if (is.null(timeval))
      settime <- min(.data[[time]])
    else
      settime <- timeval
    if (nset == 1) {
      j = 0
      for (i in setseq) {
        j <- j + 1
        addinc <- inc[j]
        temdata <- temdata %>% dplyr::group_by_at(id) %>%
          dplyr::mutate(`:=`(
            !!base::as.symbol(outcome),
            dplyr::if_else(
              .data[[time]] > settime & occtemp ==
                max(.data[["occtemp"]]) - i,
              (.data[[outcome]] +
                 addinc),
              .data[[outcome]]
            )
          ))
      }
    }
    else {
      j = 0
      for (i in setseq) {
        j <- j + 1
        addinc <- inc[j]
        temdata <- temdata %>% dplyr::group_by_at(id) %>%
          dplyr::mutate(`:=`(
            !!base::as.symbol(outcome),
            dplyr::if_else(
              .data[[time]] > settime & occtemp ==
                max(.data[["occtemp"]]) - i,
              cummax(.data[[outcome]] +
                       addinc),
              .data[[outcome]]
            )
          ))
      }
    }
    temdata2 <- temdata %>% dplyr::select(-c("occtemp",
                                             "nocctemp"))
    return(temdata2)
  }




#' Save list of ggplot2 objects to single pdf
#'
#' @param list A list of ggplot2 objects.
#' @param filename A character string to name the pdf filr.
#'
#' @return Invisible NULL.
#' @keywords internal
#' @noRd
#'
#' @examples
#' #plot histogram of each numeric variable in iris
#' list_iris = map(names(iris[-5]), ~ggplot(iris, aes_string(.)) + geom_histogram())
#' #save to a single pdf
#' GG_save_pdf(list_iris, "test.pdf")
GG_save_pdf = function(list, filename, 
                       width = 10, height = 7,
                       onefile = TRUE, compress = TRUE) {
  #start pdf
  filename <- paste0(filename, ".", "pdf")
  grDevices::pdf(filename, width = width, height = height,
      onefile = onefile, compress = compress)
  #loop
  for (p in list) {
    print(p)
  }
  #end pdf
  grDevices::dev.off()
  invisible(NULL)
}




#' An internal to check if 'bsitar' argument such as xfun is set or NULL
#'
#' @param x A symbol or a character string 
#'
#' @return A logical TRUE/FALSE
#' @keywords internal
#' @noRd
#'
check_if_arg_set <- function(x) {
  if(is.null(x)) {
    set_x <- FALSE
  } else if(is.character(x)) {
    if(x == "NULL") {
      set_x <- FALSE
    } else {
      set_x <- TRUE
    }
  } else if(is.list(x)) {
    if(length(x) == 1) {
      if(is.null(x[[1]])) {
        set_x <- FALSE
      } else {
        set_x <- TRUE
      }
    } else if(length(x) > 1) {
      if(is.null(x[[1]][1])) {
        set_x <- FALSE
      } else {
        set_x <- TRUE
      }
    }
  }
  return(set_x)
}



#' Check if 'bsitar' argument such as xfun is set or NULL
#'
#' @param x A character string 
#' @param splitat A character at which string to be split default (\code{NULL}) 
#'
#' @return A logical TRUE/FALSE
#' @keywords internal
#' @noRd
#'
remove_between_first_last_parnth <- function(x, splitat = NULL) {
  a <- sub("^.*?\\(", "", x)
  b <- sub(")\\s*$", "", a)
  if(!is.null(splitat)) {
    b <- strsplit(b, splitat)[[1]]
  }
  b
}


