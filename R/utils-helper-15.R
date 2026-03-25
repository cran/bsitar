


#' An internal function to interpolate data for plotting smooth curves
#' 
#' @param model An object of class \code{bgmfit}. This is optional (default
#'   \code{NULL}) i.e., it is not neccessary to specify the model object. When
#'   \code{model} is specified, then values for \code{newdata}, \code{idvar},
#'   and \code{xvar} are automatically taken from the \code{model}.
#'
#' @param newdata A data frame. If \code{NULL} (default), data analysed in the
#'   original model fit is used.
#'
#' @param idvar A character string to specify the group identifier. If
#'   \code{NULL} (default), \code{id} from the model fit is used.
#'
#' @param xvar  A character string to specify the time variable. If
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
#'   
#' @param keeplevels A logical in case factor variables other than \code{idvar}
#'   present
#'   
#' @param asdf A logical to indicate whether to return the value as data frame.
#' 
#' @param newdata_fixed
#' 
#' @param verbose
#'   
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
           idvar = NULL,
           xvar = NULL,
           times = NULL,
           length.out = 10,
           xrange = 1, 
           keeplevels = FALSE, 
           asdf = FALSE,
           newdata_fixed = NULL,
           verbose = FALSE) {
    
    `.` <- NULL;
    
    if (is.null(newdata)) {
      newdata <- model$model_info$bgmfit.data
    } else {
      newdata <- newdata
    }
    
    
    if(is.null(model)) {
      if(is.null(newdata_fixed)) {
        stop2c("The 'model' must be specified when newdata_fixed = 'NULL'")
      } else {
        if(newdata_fixed == 1) {
          stop2c("The 'model' must be specified when newdata_fixed = 1")
        }
      }
    }
    
    add_just_list_c <- FALSE
    if(is.null(newdata_fixed)) {
      newdata <- prepare_data2(data = newdata, model = model)
      newdata <- prepare_transformations(data = newdata, model = model)
    } else if(!is.null(newdata_fixed)) {
      if(newdata_fixed == 0) {
        # do nothing assuming that user has set up for 'uvarby' etc
        # not even applied 'dummy_to_factor'. only list_c elements will be added
        newdata <- newdata 
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
      if (is.null(idvar)) stop("Specify model or idvar")
      if (is.null(xvar)) stop("Specify model or xvar")
    }
    
    if (is.null(idvar)) {
      idvar <- model$model_info$idvars
    } else {
      idvar <- idvar
    }
    
    if (is.null(xvar)) {
      xvar <- model$model_info$xvar
    } else {
      xvar <- xvar
    }
    
    all_times <- TRUE
    if (is.null(xrange))
      xrange <- 1
    else
      xrange <- xrange
    times_orig <- newdata[[xvar]]
    times_orig <- times_orig[!is.na(times_orig)]
    
    
    if (is.null(times) || !is.numeric(times)) {
      times <-
        seq(min(times_orig), max(times_orig), length.out = length.out)
    }
    
    # This is when no random effects are specified or groupvar is NULL
    # A dummy group var is created. Check utils-helper function lines 60
    if (nlevels(newdata[[idvar]]) == 1) {
      newdata <- newdata %>%
        dplyr::distinct(newdata[[xvar]], .keep_all = T) %>%
        dplyr::arrange(!!as.name(xvar))
    }
    
    id <- match(newdata[[idvar]], unique(newdata[[idvar]]))
    
    if(length( unique(newdata[[idvar]])) == 1) {
      if(length.out == 1) stop("The argument 'ipts' should be > 1")
    }
    
    last_time  <- tapply(newdata[[xvar]], id, max)
    first_time <- tapply(newdata[[xvar]], id, min)
    
    newdata_nomiss <- newdata[complete.cases(newdata),]
    id_nomiss <-
      match(newdata_nomiss[[idvar]], unique(newdata_nomiss[[idvar]]))
    n <- length(unique(id_nomiss))
    
    if (xrange == 1) {
      times_to_pred <- list()
      for (i in 1:length(unique(newdata[[idvar]]))) {
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
    
    # Why this id_pred here??
    # id_pred <- rep(seq_len(n), sapply(times_to_pred, length))
    
    right_rows <- function (data, times, ids, Q_points) {
      fids <- factor(ids, levels = unique(ids))
      if (!is.list(Q_points))
        Q_points <- split(Q_points, row(Q_points))
      if(base::is.unsorted(times)) {
        times <- sort(times) # Faced issue newdata[[xvar]] not sorted 
      }
      ind <- mapply(findInterval, Q_points, split(times, fids))
      ind[ind < 1] <- 1
      rownams_id <- split(row.names(data), fids)
      ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
      data[c(ind),]
    }
    
    if(is.null(length.out)) {
      newdata_pred <- newdata
    } else {
      if(length.out == 1) {
        stop("The 'ipts' should be either 'NULL' or > 1")
      } else {
        newdata_pred <- right_rows(newdata, newdata[[xvar]], id, 
                                   times_to_pred)
      }
    }
    
    # newdata_pred <- right_rows(newdata, newdata[[xvar]], id, times_to_pred)
    if(keeplevels) {
      if(length(setdiff(is.fact, idvar)) > 0) {
        newdata_pred <- newdata_pred %>% droplevels
        newdata_pred <- newdata_pred %>% 
          dplyr::select(-dplyr::all_of(setdiff(is.fact, idvar)))
        
        newdata_is.factx <- newdata %>% 
          dplyr::select(dplyr::all_of(is.fact))
        
        newdata_pred <- newdata_pred %>% 
          dplyr::left_join(., newdata_is.factx, by = idvar,
                    relationship = "many-to-many")
      }
    }
    
    newdata_pred[[xvar]] <- unlist(times_to_pred)
    
    if(keeplevels) {
      newdata_pred <- newdata_pred %>% dplyr::select(dplyr::all_of(cnames))
    }
    
    if(setasdt) {
      newdata_pred <- data.table::as.data.table(newdata_pred)
    }
    
    if(asdf) {
      out <- as.data.frame(newdata_pred)  
    } else {
      out <- newdata_pred 
    }
    
   return(out)
  }



