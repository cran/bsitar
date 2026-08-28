


#' An internal function to interpolate data for plotting smooth growth curves
#' 
#' @noRd
#'
get_idata <- function(model = NULL,
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
    right_rows <- function (data, times, ids, Q_points) {
      fids <- factor(ids, levels = unique(ids))
      if (!is.list(Q_points))
        Q_points <- split(Q_points, row(Q_points))
      if(base::is.unsorted(times)) {
        times <- sort(times) 
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



