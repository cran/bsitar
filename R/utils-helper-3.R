

#' An internal function to set up data for \code{bsitar} model
#'
#' @param data An input data frame.
#'
#' @param xvar The predictor (typically age) variables.
#'
#' @param yvar The outcome variables. Must be a single name except when fitting a
#'   multivariate model.
#'
#' @param idvar The group identifier.
#'
#' @param nys An integer indicating the number of outcomes.
#'
#' @param univariate_by An optional (default \code{NA}) to specify the indicator
#'   variable for fitting univariate-by-subgroup model. See \code{univariate_by}
#'   argument in the [bsitar::bsitar()] function. If not \code{NA}, then it
#'   should be a valid factor variable present in the \code{data}.
#'
#' @param multivariate A logical (default \code{FALSE}) to specify the the
#'   multivariate model. See \code{multivariate} argument in the
#'   [bsitar::bsitar()] function.
#'
#' @param outliers An optional (default \code{NULL}) to remove velocity
#'   outliers. The argument should be a named list to pass options to the
#'   [bsitar::outliers()] function. See [bsitar::outliers()] for details.
#'
#' @param subset A logical (default \code{TRUE}) to indicate whether to create
#'   data for each level of the \code{univariate_by} variable, or only for a
#'   subset of levels. The \code{subset = TRUE} is typically used during model
#'   fit and \code{subset = FALSE} during post processing of each sub model. The
#'   argument \code{subset} is ignored when \code{univariate_by} is \code{NA} or
#'   \code{NULL}.
#'
#' @param sigmaxvar An optional (default \code{NULL}) predictor (typically age)
#'   variable for \code{sigma}
#'   
#' @param returnys A logical (default \code{FALSE})
#' 
#' @param model A placeholder Ignored (default \code{NULL})
#'
#' @param resp A placeholder Ignored (default \code{NULL})
#'
#' @param verbose A logical (default \code{FALSE})
#'
#' @param displayit A character string logical (default \code{"})
#'
#' @param setcolb A character string logical (default \code{"})
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
prepare_data2 <- function(data = NULL,
                          xvar = NULL,
                          yvar = NULL,
                          idvar = NULL,
                          nys = NULL,
                          univariate_by = NULL,
                          multivariate = NULL,
                          outliers = NULL,
                          subset = NULL,
                          sigmaxvar = NULL,
                          returnys = FALSE,
                          verbose = FALSE,
                          displayit = "",
                          setcolb = "",
                          model = NULL,
                          resp = NULL,
                          envir = NULL) {

  . <- NULL;
  if(is.null(subset)) {
    subset <- FALSE
  }
  if(is.null(envir)) {
    enverr. <- parent.frame()
  } else {
    enverr. <- envir
  }
  
  if(is.null(data) & is.null(model)) {
    stop("specify at least one of the data or model")
  }  
  
  if(!is.null(model)) {
    if(is.null(data)) {
       data <- model$model_info$bgmfit.data
      if(verbose) message("'data' is extracted from the 'model'")
    } 
    if(is.null(xvar)) {
      xvar <- extract_names_from_call(model = model, arg = "x")
      if(verbose) message("'xvar' is extracted from the 'model'")
    } 
    if(is.null(yvar)) {
      yvar <- extract_names_from_call(model = model, arg = "y")
      if(verbose) message("'yvar' is extracted from the 'model'")
    } 
    if(is.null(idvar)) {
      idvar <- extract_names_from_call(model = model, arg = "id")
      if(verbose) message("'idvar' is extracted from the 'model'")
    } 
    if(is.null(sigmaxvar)) {
      sigmaxvar <- model$model_info$sigmaxvars
      if(verbose) message("'sigmaxvar' is extracted from the 'model'")
    } 
    if(is.null(univariate_by)) {
      univariate_by <- model$model_info$univariate_by
      if(verbose) message("'univariate_by' is extracted from the 'model'")
    } 
    if(is.null(multivariate)) {
      multivariate <- model$model_info$multivariate
      if(verbose) message("'multivariate' is extracted from the 'model'")
    } 
    if(is.null(nys)) {
      nys <- model$model_info$nys
      if(verbose) message("'nys' is extracted from the 'model'")
    } 
  } # if... else if(!is.null(model)) {
  
  data   <- data %>% droplevels()
  uvarby <- univariate_by$by
  mvar   <- multivariate$mvar
  
  if(is.null(uvarby)) {
    uvarby <- NA
  }

  if (!(is.na(uvarby) | uvarby == "NA")) {
    check_if_any_varibale_all_NA(data, factor_var = uvarby)
    ys <- c()
    for (l in levels(data[[uvarby]])) {
      # why subset here?
      if(!subset) data[[l]] <- data[[yvar[1]]]
      if( subset) data[[l]] <- data[[l]]
      if(!subset)        ys <- c(ys, l[1])
      if( subset)        ys <- c(ys, l)
    }
  } else { # end need to set ys for univariate by levels
    check_if_any_varibale_all_NA(data, factor_var = NULL)
    ys      <- yvar
  }

  xs      <- xvar
  ids     <- idvar
  
  if (nys > 1) {
    unique_xs  <- unique(xs)
    unique_ys  <- unique(ys)
    unique_ids <- unique(ids)
    if(length(unique_xs) == 1) {
      xs <- c()
      for (j in unique_ys) {
        tempname <- paste0(unique_xs, "_", j)
        if(is.null(data[[tempname]])) {
          data[[tempname]] <- data[[unique_xs]]
        }
        xs <- c(xs, tempname)
        rm('tempname')
      }
    } else if(length(unique_xs) != nys) { 
      stop2c("The number of 'xvar' variables shoud be either 1 or 
             same as the 'yvar'")
    } # if(length(unique_xs) == 1) { else if(length(unique_xs) == 1) {
  } # if (nys > 1) {
  

  
  if(!is.null(model)) {
    if (is.null(resp)) {
      resp_    <- resp
      revresp_ <- ""
    } else if (!is.null(resp)) {
      resp_    <- paste0(resp, "_")
      revresp_ <- paste0("_", resp)
    }
    
    sigma_model_      <- paste0('sigmamodel', revresp_)
    sigma_model_name_ <- paste0('sigmabasicfunname', revresp_)
    sigma_model_attr_ <- paste0('sigmabasicfunattr', revresp_)
    
    sigma_model       <- model$model_info[[sigma_model_]]
    sigma_model_name  <- model$model_info[[sigma_model_name_]]
    sigma_model_attr  <- model$model_info[[sigma_model_attr_]]
  }
  

  sigmaxs_c <- c()
  for (j in 1:nys) {
    if(is.na(sigmaxvar[j]) |  sigmaxvar[j] == "NA") {
      addsigmaxvar <- NA
    } else if(isTRUE(sigmaxvar[j]) |  sigmaxvar[j] == "TRUE") {
      addsigmaxvar <- paste0('sigma', xs[j])
      data[[addsigmaxvar]] <- data[[xs[j]]]
    } else {
      addsigmaxvar <- tempbane <- sigmaxvar[j]
      if(!addsigmaxvar %in% names(data)) {
        stop2c("The variable ", collapse_comma(addsigmaxvar), " used as an ",
             "\n ",
             " argument for 'sigmax' is not available in the dataframe")
      }
      if(nys > 1) {
        addsigmaxvar <- paste0(addsigmaxvar, "_", ys[j])
      } else {
        addsigmaxvar <- addsigmaxvar
      }
      # addsigmaxvar <- paste0(addsigmaxvar, "_", ys[j])
      data[[addsigmaxvar]] <- data[[tempbane]]
    }
    sigmaxs_c <- c(sigmaxs_c, addsigmaxvar)
  } # for (j in 1:nys) {
  
  sigmaxs <- sigmaxs_c
  
  xvar      <- xs
  yvar      <- ys
  idvar     <- ids
  org.data  <- data
  
  
  if (!is.null(outliers)) {
    remove_ <- outliers$remove
    icode_ <- outliers$icode
    icode_ <- deparse(substitute(icode_))
    limit_ <- outliers$limit
    velpower_ <- outliers$velpower
    lag_ <- outliers$lag
    linearise_ <- outliers$linearise
    verbose_ <- outliers$verbose

    for (yi in 1:length(yvar)) {
      if (!yvar[yi] %in% colnames(data)) {
        stop(
          "When model is fit with argument outliers (outliers not NULL), ",
          "\n",
          "  then outcome variable should be part of the newdata specified.",
          "\n",
          "  please check the missing outcome varibale: ",
          yvar[yi]
        )
      }
      if (!xvar[yi] %in% colnames(data)) {
        stop(
          "When model is fit with argument outliers (i.e., outliers not NULL),",
          " \n ",
          "  then predictor variable should be part of the newdata specified.",
          "\n",
          "  please check the missing predictor varibale: ",
          xvar[yi]
        )
      }
      if (!idvar[yi] %in% colnames(data)) {
        stop(
          "When model is fit with argument outliers
          (i.e., outliers not NULL), ",
          "\n",
          "  then group identifier variable should be
          part of the newdata specified.",
          "\n",
          "  please check the missing group identifier varibale: ",
          idvar[yi]
        )
      }
      data <-
        outliers(
          x = xvar[yi],
          y =  yvar[yi],
          id = idvar[yi],
          data = data,
          icode = icode_,
          lag = lag_,
          velpower = velpower_,
          limit = limit_,
          linearise = linearise_,
          remove = remove_,
          verbose = verbose_
        )
    }
  } # if(!is.null(outliers)) {

  # Internal argument 'uvarby_method2' to set data for multivariate framework 
  # for uvarby. Did not work. The 'uvarby_method1' is the traditional and 
  # correct approach.
  uvarby_method <- 'uvarby_method1'
  
  if (!(is.na(uvarby) | uvarby == "NA")) {
    if (!uvarby %in% colnames(data)) {
      stop(paste(
        "\nvariable",
        uvarby,
        "used for setting univariate submodels is missing"
      ))
    }
    if (!is.factor(data[[uvarby]])) {
      stop("subset by variable '",
           uvarby,
           "' should be a factor variable")
    }
    if(uvarby_method == 'uvarby_method1') {
      for (l in levels(data[[uvarby]])) {
        if(!subset) data[[l]] <- data[[yvar[1]]]
        if(subset) data[[l]]  <- data[[l]]
      }
      unibyimat <- model.matrix(~ 0 + eval(parse(text = uvarby)), data)
      subindicators <- paste0(uvarby, levels(data[[uvarby]]))
      colnames(unibyimat) <- subindicators
      unibyimat <- unibyimat %>% data.frame()
      unibyimat <- sapply(unibyimat, as.integer ) %>% data.frame()
      yvar <- levels(data[[uvarby]])
      data <- as.data.frame(cbind(data, unibyimat))
    }
    if(uvarby_method == 'uvarby_method2') {
      id_colsx <- setdiff(colnames(data), c(yvar, uvarby))
      uvarbyx  <- levels(data[[uvarby]])
      data <-
        data %>% data.frame() %>% 
        tidyr::pivot_wider(., id_cols= dplyr::all_of(id_colsx), 
                           names_from = uvarby, values_from = yvar) %>% 
        dplyr::mutate(dplyr::across(dplyr::all_of(uvarbyx), 
                             ~ dplyr::if_else(is.na(.x), FALSE, TRUE), 
                             .names = "{.col}s"))
      data[[uvarby]] <- org.data[[uvarby]]
      subindicators <- paste0(uvarbyx, 's')
      yvar <- uvarbyx
    } # if(uvarby_method == 'uvarby_method2') {
    

    if (!is.na(uvarby) & univariate_by$verbose) {
      resvcts_ <- levels(data[[uvarby]])
      resvcts <- paste0(resvcts_, collapse = " ")
      setmsgtxt <- paste0(
        "\n For univariate-by-subgroup model fitting for variable '",
        uvarby,
        "'",
        " (specified via 'univariate_by' argument)",
        "\n ",
        resvcts,
        " response vectors created based on the factor levels",
        "\n\n ",
        "Please check corresponding arguments list.",
        " E.g, df = list(4, 5) denotes that\n df = 4 is for ",
        resvcts_[1],
        ", and  df = 5 is for ",
        resvcts_[2],
        " (and similalry knots, priors, initials etc)",
        "\n\n ",
        "If it does't correspond correctly, then either reverse the list ",
        "arguments\n such as df = list(5, 4),",
        " or else reverse sort the order of factor levels"
      )
      if (displayit == 'msg') {
        message(setmsgtxt)
      } else if (displayit == 'col') {
        col <- setcolb
        cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
      }
    }
    
    
    if (!(is.na(uvarby) | uvarby == "NA")) {
      sortbylayer <- NA
      data <- data %>%
        dplyr::mutate(sortbylayer =
                        forcats::fct_relevel(!!as.name(uvarby),
                                             (levels(
                                               !!as.name(uvarby)
                                             )))) %>%
        dplyr::arrange(sortbylayer) %>%
        dplyr::select(-sortbylayer)
    }
    
    
    if(returnys) return(yvar)
    
    attr(data, "multivariate")  <- FALSE
    attr(data, "uvarby")        <- uvarby
    attr(data, "subindicators") <- subindicators
    # data_out <- data
  } else if (mvar) {
    data <- org.data
    attr(data, "multivariate")  <- TRUE
    attr(data, "uvarby")        <- NULL
    attr(data, "subindicators") <- NULL
  } else {
    data <- org.data
    attr(data, "multivariate")  <- FALSE
    attr(data, "uvarby")        <- NULL
    attr(data, "subindicators") <- NULL
  }
  attr(data, "xs") <- xvar
  attr(data, "ys") <- yvar
  attr(data, "ids") <- idvar
  attr(data, "sigmaxs") <- sigmaxs # sigmaxvar
  return(data)
}




