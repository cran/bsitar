

#' An internal function to prepare transformations for \code{bsitar} model
#'
#' @param data An input data frame.
#'
#' @param xvar The predictor (typically age) variables.
#'
#' @param yvar The outcome variables. Must be a single name except when fitting
#'   a multivariate model.
#'
#' @param sigmaxvar The predictor (typically age) variables for \code{sigma}
#'
#' @param xfun Optional name(s) of the transformation function(s) applied to the
#'   predictor variable (typically age). Default \code{NULL}.
#'
#' @param yfun Optional name(s) of the transformation function(s) applied to the
#'   outcome variable. Default \code{NULL}.
#'
#' @param sigmaxfun Optional name(s) of the transformation function(s) applied
#'   to the predictor variable (typically age)  for \code{sigma}. Default
#'   \code{NULL}.
#'
#' @param ixfun Optional name(s) of the inverse transformation function(s)
#'   applied to the predictor variable (typically age). Default \code{NULL}.
#'
#' @param iyfun Optional name(s) of the inverse transformation function(s)
#'   applied to the outcome variable. Default \code{NULL}.
#'
#' @param sigmaixfun Optional name(s) of the inverse transformation function(s)
#'   applied to the predictor variable (typically age)  for \code{sigma}.
#'   Default \code{NULL}.
#' 
#' @param xoffset A real number.
#' 
#' @param sigmaxoffset A real number.
#'    
#' @param transform A character vector to specify variables to be transformed.
#'   currently ignored.
#' 
#' @param itransform A character vector to specify variables to be inverse
#'   transformed. The options are \code{c('x', 'y', 'sigma')} or an empty string
#'   i.e., \code{itransform = ""} (default). The \code{itransform} option is
#'   particularly useful when using just \code{data} and \code{model} arguments
#'   to prepare transformation during post-processing. Since all other
#'   information is retrieved from the \code{model_info}, setting
#'   \code{itransform = ""} will automatically set inverse transformations to
#'   \code{FALSE} for \code{xvar}, \code{yvar}, and \code{sigmaxvar}.
#'   
#' @param restore_decimal Optional logical (default \code{TRUE}) to indicate
#'   whether to restore the decimal places for \code{xvar}, \code{yvar}, and
#'   \code{sigmaxvar} variables after transformation. This does not make any
#'   substantial changes but only to exact match the recovered variable values.
#'   
#' @param envir A logical (default \code{TRUE})
#'
#' @return A data frame with necessary information added a attributes.
#'
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
prepare_transformations <- function(data = NULL,
                                    xvar = NULL, 
                                    yvar = NULL,
                                    sigmaxvar = NULL,
                                    xfun = NULL,
                                    yfun = NULL,
                                    sigmaxfun = NULL,
                                    ixfun = FALSE,
                                    iyfun = FALSE,
                                    sigmaixfun = FALSE,
                                    xoffset = NULL,
                                    sigmaxoffset = NULL,
                                    transform = "",
                                    itransform = "",
                                    model = NULL,
                                    envir = NULL,
                                    verbose = FALSE,
                                    restore_decimal = TRUE) {
  
  
  #######################################################################
  #######################################################################
  # # Functions dictionary:
  # # functions for xvar 
  # 1: xfuns - names of function as string e.g., 'log' 'sqrt'
  # 2: xfun_y - outcome specific names of function as string
  # 3: xfuntransforms - function e.g., function(x)x
  # 4: xfuntransform_y - outcome specific function
  # 4: xfuntransform2s - function with offset included e.g., function(x)x-12
  # 6: xfuntransform2_y - outcome specific function with offset included
  # # functions for yvar 
  # # same as above with prefix y instead of x
  # # functions for sigmaxvar 
  # # same as above with prefix sigmax instead of sigmax
  # # inverse functions: just prefix in before x, y or sigma
  #######################################################################
  #######################################################################
  
  
  if(is.null(envir)) {
    enverr. <- parent.frame()
  } else {
    enverr. <- envir
  }
  
  if(is.null(data) & is.null(model)) {
    stop("specify at least one of the data or model")
  } else if(!is.null(model)) {
    if(is.null(data)) {
      # data <- model$data
      # if(verbose) message("'data' is extracted from the 'model'")
      stop("data must be specified even when model is not NULL")
    } 
    if(is.null(xvar)) {
      xvar <- model$model_info$xvars
      if(verbose) message("'xvar' is extracted from the 'model'")
    } 
    if(is.null(yvar)) {
      yvar <- model$model_info$yvars
      if(verbose) message("'yvar' is extracted from the 'model'")
    } 
    if(is.null(sigmaxvar)) {
      sigmaxvar <- model$model_info$sigmaxvars
      if(verbose) message("'sigmaxvar' is extracted from the 'model'")
    } 
    if(is.null(xfun)) {
      xfun <- model$model_info$xfuntransforms
      if(verbose) message("'xfun' is extracted from the 'model'")
    } 
    if(is.null(yfun)) {
      yfun <- model$model_info$yfuntransforms
      if(verbose) message("'yfun' is extracted from the 'model'")
    } 
    if(is.null(sigmaxfun)) {
      sigmaxfun <- model$model_info$sigmaxfuntransforms
      if(verbose) message("'sigmaxfun' is extracted from the 'model'")
    } 
    if(is.null(xoffset)) {
      xoffset <- model$model_info$xoffsets
      if(verbose) message("'xoffset' is extracted from the 'model'")
    } 
    if(is.null(sigmaxoffset)) {
      sigmaxoffset <- model$model_info$sigmaxoffsets
      if(verbose) message("'sigmaxoffset' is extracted from the 'model'")
    } 
  } # if... else if(!is.null(model)) {
  
  
  if(!is.null(model)) {
    if('x' %in% itransform) {
      ixfun      <- TRUE
    } else if(!'x' %in% itransform) {
      ixfun      <- FALSE
    }
    
    if('y' %in% itransform) {
      iyfun      <- TRUE
    } else if(!'y' %in% itransform) {
      iyfun      <- FALSE
    }
     
    if('sigma' %in% itransform) {
      sigmaixfun <- TRUE
    } else if(!'sigma' %in% itransform) {
      sigmaixfun <- FALSE
    }
  } # ... if(!is.null(model)) {
  
  
  # Over ride when variable itself is NULL
  if(is.null(xvar))      ixfun      <- FALSE
  if(is.null(yvar))      iyfun      <- FALSE
  # if(is.null(sigmaxvar)) sigmaixfun <- FALSE
  
  
  for (j in 1:length(sigmaxvar)) {
    if(is.null(sigmaxvar[j])) {
      sigmaxvar  <- NULL
      sigmaxfun  <- NULL
      sigmaixfun <- FALSE
    } else if(is.na(sigmaxvar[j]) |  sigmaxvar[j] == "NA") {
      sigmaxvar  <- NULL
      sigmaxfun  <- NULL
      sigmaixfun <- FALSE
    } else {
      sigmaxvar[j] <- sigmaxvar[j]
    }
  }
  
  
  #######################################################################
  #######################################################################
  
  if(restore_decimal) {
    # xvar 
    if(!is.null(xvar)) {
      ndecimal_xvar <- c()
      for (ii in 1:length(xvar)) {
        i <- xvar[ii]
        ndecimal_xvar <- c(ndecimal_xvar, 
                           max(get_decimal_places(data, i))
        )
      }
    } else if(is.null(xvar)) {
      ndecimal_xvar <- rep(NA)
    }
    # yvar 
    if(!is.null(yvar)) {
      ndecimal_yvar <- c()
      for (ii in 1:length(yvar)) {
        i <- yvar[ii]
        ndecimal_yvar <- c(ndecimal_yvar, 
                           max(get_decimal_places(data, i))
        )
      }
    } else if(is.null(yvar)) {
      ndecimal_yvar <- rep(NA)
    }
    # sigmaxvar 
    if(!is.null(sigmaxvar)) {
      ndecimal_sigmaxvar <- c()
      for (ii in 1:length(sigmaxvar)) {
        i <- sigmaxvar[ii]
        ndecimal_sigmaxvar <- c(ndecimal_sigmaxvar, 
                                max(get_decimal_places(data, i))
        )
      }
    } else if(is.null(sigmaxvar)) {
      ndecimal_sigmaxvar <- rep(NA)
    }
  } # if(restore_decimal) {
  
  
  #######################################################################
  #######################################################################
  # check and allow function as 'xfun' and also match 'xvar'
  if(!is.null(xfun)) {
    if(!is.list(xfun)) {
      if(is.function(xfun)) {
        tempfuns <- list()
        for (i in 1:length(xvar)) {
          tempfuns[[i]] <- xfun
        }
        xfun <- tempfuns
      }
    } else if(is.list(xfun)) {
      if(length(xfun) == 1) {
        if(!is.function(xfun[[1]])) stop("xfun must be a function")
        tempfuns <- list()
        for (i in 1:length(xvar)) {
          tempfuns[[i]] <- xfun[[1]]
        }
        xfun <- tempfuns
      } else if(length(xfun) != length(xvar)) {
        stop("length of 'xfun' should be either one or same as 'xvar'")
      } else {
        for (i in 1:length(xfun)) {
          if(!is.function(xfun[[i]])) stop("xfun must be a function")
        }
      }
    }
  }
  
  # check and allow function as 'yfun' and also single to match 'yvar'
  if(!is.null(yfun)) {
    if(!is.list(yfun)) {
      if(is.function(yfun)) {
        tempfuns <- list()
        for (i in 1:length(yvar)) {
          tempfuns[[i]] <- yfun
        }
        yfun <- tempfuns
      }
    } else if(is.list(yfun)) {
      if(length(yfun) == 1) {
        if(!is.function(yfun[[1]])) stop("yfun must be a function")
        tempfuns <- list()
        for (i in 1:length(yvar)) {
          tempfuns[[i]] <- yfun[[1]]
        }
        yfun <- tempfuns
      } else if(length(yfun) != length(yvar)) {
        stop("length of 'yfun' should be either one or same as 'yvar'")
      } else {
        for (i in 1:length(yfun)) {
          if(!is.function(yfun[[i]])) stop("yfun must be a function")
        }
      }
    }
  }
  
  # check and allow function as 'sigmaxfun' and also single to match 'sigmaxvar'
  if(!is.null(sigmaxfun)) {
    if(!is.list(sigmaxfun)) {
      if(is.function(sigmaxfun)) {
        tempfuns <- list()
        for (i in 1:length(sigmaxvar)) {
          tempfuns[[i]] <- sigmaxfun
        }
        sigmaxfun <- tempfuns
      }
    } else if(is.list(sigmaxfun)) {
      if(length(sigmaxfun) == 1) {
        if(!is.function(sigmaxfun[[1]])) stop("sigmaxfun must be a function")
        tempfuns <- list()
        for (i in 1:length(sigmaxvar)) {
          tempfuns[[i]] <- sigmaxfun[[1]]
        }
        sigmaxfun <- tempfuns
      } else if(length(sigmaxfun) != length(sigmaxvar)) {
        stop("length of 'sigmaxfun' should be either one or same as 'sigmaxvar'")
      } else {
        for (i in 1:length(sigmaxfun)) {
          if(!is.function(sigmaxfun[[i]])) stop("sigmaxfun must be a function")
        }
      }
    }
  }
  
  #######################################################################
  #######################################################################
  
  # if xoffset = NULL, even though it will be given 0, but name not assigned
  if(is.null(xoffset))      xoffset      <- 0
  if(is.null(sigmaxoffset)) sigmaxoffset <- 0
  
  # This is_emptyx will handle numeric(0) when xoffset was set as NULL
  if(is_emptyx(xoffset))      xoffset      <- 0
  if(is_emptyx(sigmaxoffset)) sigmaxoffset <- 0
  
  
  # Assign xoffset names same as xvar 
  if(!is.null(xoffset) & !is.null(xvar)) {
    templist <- itemplist <- list()
    for (i in 1:length(xvar)) {
      templist[[xvar[[i]]]] <- xoffset[[i]]
    }
    xoffset <- templist
    if(is_emptyx(xoffset)) xoffset <- 0
  } else {
    xoffset <- 0
  }
  
  
  # Assign sigmaxoffset names same as sigmaxvar 
  if(!is.null(sigmaxoffset) & !is.null(sigmaxvar)) {
    templist <- itemplist <- list()
    for (i in 1:length(sigmaxvar)) {
      templist[[sigmaxvar[[i]]]] <- sigmaxoffset[[i]]
    }
    sigmaxoffset <- templist
    if(is_emptyx(sigmaxoffset) | is.null(sigmaxoffset)) sigmaxoffset <- 0
  } else {
    sigmaxoffset <- 0
  }
  
 
  #######################################################################
  #######################################################################
  
  # check and tranform TRUE/FALSE ixfun iyfun sigmaixfun to NULL or funs
  
  # Although below we keep option of setting ifuns, for now we allow only T/F
  if(missing(ixfun)) {
    ixfun <- FALSE
    if(verbose) message("missing 'ixfun' set as FALSE")
  }
  if(missing(iyfun)) {
    iyfun <- FALSE
    if(verbose) message("missing 'iyfun' set as FALSE")
  }
  if(missing(sigmaixfun)) {
    iyfun <- FALSE
    if(verbose) message("missing 'sigmaixfun' set as FALSE")
  }
  
  if(!is.logical(ixfun)) {
    stop("'ixfun' must be a logical, TRUE or FALSE")
  }
  if(!is.logical(iyfun)) {
    stop("'iyfun' must be a logical, TRUE or FALSE")
  }
  if(!is.logical(sigmaixfun)) {
    stop("'sigmaixfun' must be a logical, TRUE or FALSE")
  }
  
  
  #######################################################################
  #######################################################################

  if(!is.null(xvar)) {
    templist <- itemplist <- list()
    for (i in 1:length(xvar)) {
      templist[[xvar[[i]]]] <- xfun[[i]]
      # itemplist[[xvar[[i]]]] <- ixfun[[i]]
    }
    xfun <- templist
    if(is_emptyx(xfun)) xfun <- NULL
  }

  if(!is.null(yvar)) {
    templist <- itemplist <- list()
    for (i in 1:length(yvar)) {
      templist[[yvar[[i]]]] <- yfun[[i]]
    }
    yfun <- templist
    if(is_emptyx(yfun)) yfun <- NULL
  }

  if(!is.null(sigmaxvar)) {
    templist <- itemplist <- list()
    for (i in 1:length(sigmaxvar)) {
      templist[[sigmaxvar[[i]]]] <- sigmaxfun[[i]]
    }
    sigmaxfun <- templist
    if(is_emptyx(sigmaxfun)) sigmaxfun <- NULL
  }
  
  
  #######################################################################
  #######################################################################
  # Re construct 'xfun' functions by including xoffset
  if(!is.null(xfun)) {
    for (i in names(xfun)) {
      evalxoffset             <- xoffset[[i]]
      bodyoffun               <- deparse(body(xfun[[i]]))
      addtobodyoffun          <- paste0("-", evalxoffset)
      bodyoffun2              <- paste0(bodyoffun, addtobodyoffun)
      body(xfun[[i]])        <- str2lang(bodyoffun2)
    }
  }
  
  
  # Re construct 'sigmaxfun' functions by including sigmaxoffset
  if(!is.null(sigmaxfun)) {
    for (i in names(sigmaxfun)) {
      evalxoffset <- sigmaxoffset[[i]]
      bodyoffun               <- deparse(body(sigmaxfun[[i]]))
      addtobodyoffun          <- paste0("-", evalxoffset)
      bodyoffun2              <- paste0(bodyoffun, addtobodyoffun)
      body(sigmaxfun[[i]])   <- str2lang(bodyoffun2)
    }
  }
  
  
  #######################################################################
  #######################################################################
  
  # Since we are allowing only T/F, when ifuns = TRUE, set fun to NULL
  # Generate inverse functions if NULL 
  if(!is.logical(ixfun)) { # i.e., ixfun not T/F
    if(is.null(ixfun)) {
      if(!is.null(xfun)) {
        ixfun <- list()
        for (i in names(xfun)) {
          ixfun[[i]] <- inverse_transform(base::body(xfun[[i]]))
        }
      } else {
        # stop("Please specify 'ixfun'")
      }
    }
  } else if(is.logical(ixfun)) { # i.e., ixfun T/F
    if(ixfun) {
      if(!is.null(xfun)) {
        ixfun <- list()
        for (i in names(xfun)) {
          ixfun[[i]] <- inverse_transform(base::body(xfun[[i]]))
        }
      }
      xfun <- NULL
    } else if(!ixfun) {
      ixfun <- NULL
    }
  }
  
  if(!is.logical(iyfun)) {
    if(is.null(iyfun)) {
      if(!is.null(yfun)) {
        iyfun <- list()
        for (i in names(yfun)) {
          iyfun[[i]] <- inverse_transform(base::body(yfun[[i]]))
        }
      } else {
        # stop("Please specify 'iyfun'")
      }
    }
  } else if(is.logical(iyfun)) { # i.e., ixfun T/F
    if(iyfun) {
      if(!is.null(yfun)) {
        iyfun <- list()
        for (i in names(yfun)) {
          iyfun[[i]] <- inverse_transform(base::body(yfun[[i]]))
        }
      }
      yfun <- NULL
    } else if(!iyfun) {
      iyfun <- NULL
    }
  }
  
  
  if(!is.logical(sigmaixfun)) {
    if(is.null(sigmaixfun)) {
      if(!is.null(sigmaxfun)) {
        sigmaixfun <- list()
        for (i in names(sigmaxfun)) {
          sigmaixfun[[i]] <- inverse_transform(base::body(sigmaxfun[[i]]))
        }
      } else {
        # stop("Please specify 'sigmaixfun'")
      }
    }
  } else if(is.logical(sigmaixfun)) { # i.e., ixfun T/F
    if(sigmaixfun) {
      if(!is.null(sigmaxfun)) {
        sigmaixfun <- list()
        for (i in names(sigmaxfun)) {
          sigmaixfun[[i]] <- inverse_transform(base::body(sigmaxfun[[i]]))
        }
      }
      sigmaxfun <- NULL
    } else if(!sigmaixfun) {
      sigmaixfun <- NULL
    }
  }
  
 
  
  #######################################################################
  #######################################################################
  
  transform <- itransform <- ""
  if(!is.null(xfun)) {
    if(is.null(ixfun)) transform <- c(transform, 'x')
  }
  if(!is.null(yfun)) {
    if(is.null(iyfun)) transform <- c(transform, 'y')
  }
  if(!is.null(sigmaxfun)) {
    if(is.null(sigmaixfun)) transform <- c(transform, 'sigma')
  }
  
  if(!is.null(ixfun)) {
    itransform <- c(itransform, 'x')
  }
  if(!is.null(iyfun)) {
    itransform <- c(itransform, 'y')
  }
  if(!is.null(sigmaixfun)) {
    itransform <- c(itransform, 'sigma')
  }
  
  transform  <- transform[  transform!=""]
  itransform <- itransform[itransform!=""]
 
  #######################################################################
  #######################################################################
  ## Transform
  if('x' %in% transform & !is.null(xfun)) {
    if(is.null(xfun)) {
      stop("you requested transformation of 'x' but xfun is 'NULL'")
    } else if(is.null(xvar)) {
      stop("you requested transformation of 'x' but xvar is 'NULL'")
    } else {
      if(length(xvar) != length(xfun)) 
        stop("Number of 'xvar' and 'xfun' differs")
      for (i in names(xfun)) { 
        evalfuns <- xfun[[i]]
        # check and do nothing if matrix, typically posterior draws
        if(tibble::is_tibble(data) | is.data.frame(data)) {
          if(!is.null(data[[i]])) data[[i]] <- evalfuns(data[[i]])
        }
        evalfuns <- NULL
      }
    }
  }
 

  if('y' %in% transform & !is.null(yfun)) {
    if(is.null(yfun)) {
      stop("you requested transformation of 'y' but yfun is 'NULL'")
    } else if(is.null(yvar)) {
      stop("you requested transformation of 'y' but yvar is 'NULL'")
    } else {
      if(length(yvar) != length(yfun)) 
        stop("Number of 'yvar' and 'yfun' differs")
      for (i in names(yfun)) {
        evalfuns <- yfun[[i]]
        if(tibble::is_tibble(data) | is.data.frame(data)) {
          if(!is.null(data[[i]])) data[[i]] <- evalfuns(data[[i]])
        }
        evalfuns <- NULL
      }
    }
  }
  
  
  if('sigma' %in% transform & !is.null(sigmaxfun)) {
    if(is.null(sigmaxfun)) {
      stop("you requested transformation of 'sigma' but sigmaxfun is 'NULL'")
    } else if(is.null(sigmaxvar)) {
      stop("you requested transformation of 'sigma' but sigmaxvar is 'NULL'")
    } else {
      if(length(sigmaxvar) != length(sigmaxvar)) 
        stop("Number of 'sigmaxvar' and 'sigmaxvar' differs")
      for (i in names(sigmaxfun)) {
        evalfuns <- sigmaxfun[[i]]
        if(tibble::is_tibble(data) | is.data.frame(data)) {
          if(!is.null(data[[i]])) data[[i]] <- evalfuns(data[[i]])
        }
        evalfuns <- NULL
      }
    }
  }
  
  #######################################################################
  #######################################################################
  ## itransform -> inverse
  if('x' %in% itransform & !is.null(ixfun)) {
    if(is.null(ixfun)) {
      stop("you requested transformation of 'x' but ixfun is 'NULL'")
    } else if(is.null(xvar)) {
      stop("you requested transformation of 'x' but xvar is 'NULL'")
    } else {
      if(length(xvar) != length(ixfun)) 
        stop("Number of 'xvar' and 'ixfun' differs")
      for (i in names(ixfun)) { # don't use xvar, it might reult in duplicate exe
        evalfuns <- ixfun[[i]]
        if(tibble::is_tibble(data) | is.data.frame(data)) {
          if(!is.null(data[[i]])) data[[i]] <- evalfuns(data[[i]])
        }
        evalfuns <- NULL
      }
    }
  }
 
  if('y' %in% itransform & !is.null(iyfun)) {
    if(is.null(iyfun)) {
      stop("you requested transformation of 'y' but iyfun is 'NULL'")
    } else if(is.null(yvar)) {
      stop("you requested transformation of 'y' but yvar is 'NULL'")
    } else {
      if(length(yvar) != length(iyfun)) 
        stop("Number of 'yvar' and 'iyfun' differs")
      for (i in names(iyfun)) {
        evalfuns <- iyfun[[i]]
        if(tibble::is_tibble(data) | is.data.frame(data)) {
          if(!is.null(data[[i]])) data[[i]] <- evalfuns(data[[i]])
        }
        evalfuns <- NULL
      }
    }
  }
  
  if('sigma' %in% itransform & !is.null(sigmaixfun)) {
    if(is.null(sigmaixfun)) {
      stop("you requested transformation of 'sigma' but sigmaixfun is 'NULL'")
    } else if(is.null(sigmaxvar)) {
      stop("you requested transformation of 'sigma' but sigmaxvar is 'NULL'")
    } else {
      if(length(sigmaxvar) != length(sigmaixfun)) 
        stop("Number of 'sigmaxvar' and 'sigmaixfun' differs")
      for (i in names(sigmaixfun)) {
        evalfuns <- sigmaixfun[[i]]
        if(tibble::is_tibble(data) | is.data.frame(data)) {
          if(!is.null(data[[i]])) data[[i]] <- evalfuns(data[[i]])
        }
        evalfuns <- NULL
      }
    }
  }
  
  
  #######################################################################
  #######################################################################
  
  if(restore_decimal) {
    # xvar
    if(!is.null(xvar)) {
      for (ii in 1:length(xvar)) {
        i <- xvar[ii]
        if(!is.na(ndecimal_xvar[i])) {
          data[[i]] <- round(data[[i]], ndecimal_xvar[i]) 
        }
      }
    }
    # yvar
    if(!is.null(yvar)) {
      for (ii in 1:length(yvar)) {
        i <- yvar[ii]
        if(!is.na(ndecimal_yvar[i])) {
          data[[i]] <- round(data[[i]], ndecimal_yvar[i]) 
        }
      }
    }
    # sigmaxvar
    if(!is.null(sigmaxvar)) {
      for (ii in 1:length(sigmaxvar)) {
        i <- sigmaxvar[ii]
        if(!is.na(ndecimal_sigmaxvar[i])) {
          data[[i]] <- round(data[[i]], ndecimal_sigmaxvar[i]) 
        }
      }
    }
  } # if(restore_decimal) {
  
  return(data)
}




