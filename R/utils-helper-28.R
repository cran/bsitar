



#' Select rows by parameter values
#'
#' Filters a data.frame or data.table to rows where a specified column matches
#' given values. Supports optional case-insensitive matching.
#'
#' @param dt A \code{data.frame} or \code{data.table} to filter.
#' @param values Character vector of values to match (e.g., \code{c("apgv", "pgv")}).
#' @param ignore_case \code{TRUE} (default) for case-insensitive matching using
#'   \code{toupper()}; \code{FALSE} for exact match.
#' @param col Character name of column to filter on (default: \code{"parameter"}).
#'
#' @return A \code{data.frame} or \code{data.table} (same class as input) with
#'   matching rows.
#'
#' @examples
#' DT <- data.frame(
#'   parameter = c("APGV", "apgv", "PGV", "pgv", "other"),
#'   sex = rep(c("Male", "Female"), 3)[1:5],
#'   estimate = 1:5
#' )
#'
#' # Case-insensitive (default)
#' get_selected_rows(DT)
#'
#' # Case-sensitive
#' get_selected_rows(DT, ignore_case = FALSE)
#'
#' @keywords internal
#' @noRd
#' 
get_selected_rows <- function(dt, 
                              values = c("apgv", "pgv"), 
                              ignore_case = TRUE, 
                              col = "parameter") {
  if (data.table::is.data.table(dt)) {
    param_col <- dt[[col]]
    if (ignore_case) {
      param_col <- toupper(param_col)
      values <- toupper(values)
    }
    return(dt[param_col %in% values, ])
  } else {
    param_col <- dt[[col]]
    if (ignore_case) {
      param_col <- toupper(param_col)
      values <- toupper(values)
    }
    return(dt[param_col %in% values, ])
  }
}


#' Remove groups with missing values
#'
#' @description
#' 
#' Removes all rows for groups where any of the selected columns contain missing
#' values, and optionally issues a warning listing the affected groups.
#' 
#' @param DT A data.table or data.frame. Data to be filtered.
#' @param variable Character vector of column names to check for missing values.
#'   If `NULL`, all columns in `DT` are checked.
#' @param group Character scalar giving the grouping variable name
#'   (e.g. `"drawid"`). Groups are removed if any of the `variable`
#'   columns are missing within that group.
#' @param verbose Logical; if `TRUE`, a warning is emitted listing the
#'   group values removed and the columns that were checked.
#'
#' @return An object of the same class as `DT` (data.table or data.frame)
#'   containing only groups for which the selected columns have no
#'   missing values. 
#'
#' @examples
#' clean_draws(DT, variable = "draw", group = "drawid")
#'
#' @keywords internal
#' @noRd
#' 
clean_draws <- function(DT, 
                        variable = NULL, 
                        group = "drawid", 
                        verbose = FALSE) {
  
  if(!is.null(DT)) {
    if(nrow(DT) == 0) {
      return(DT)
    }
  } else if(is.null(DT)) {
    return(DT)
  }
  
  has_na <- NULL;
  . <- NULL;
  is_dt <- data.table::is.data.table(DT)
  
  if (!is_dt) {
    DT <- data.table::as.data.table(DT)
  }
  
  if (is.null(variable)) {
    variable <- base::names(DT)
  }
  
  bad_ids <- DT[
    ,
    .(has_na = base::any(base::sapply(.SD, base::anyNA))),
    by = group,
    .SDcols = variable
  ][has_na == TRUE, base::get(group)]
  
  out <- DT[!(base::get(group) %in% bad_ids)]
  
  if (verbose && base::length(bad_ids)) {
    base::warning(
      "Removed ", group, " value(s) due to NA in: ",
      base::paste(variable, collapse = ", "),
      " | ", group, "s. removed group ids are: ",
      base::paste(bad_ids, collapse = ", ")
    )
  }
  
  # convert back to data.frame if original was not data.table
  if (!is_dt) {
    out <- base::as.data.frame(out)
  }
  
  return(out)
}




#' The main function get_comparison_hypothesis for hypothesis test
#'
#' @param data A fitted model object of class \code{bgmfit}, or a posterior data
#'   frame. See Details for engine-specific support.
#' @param full.args See [hypothesis_test()] and [get_growthparameters()]
#'   for details.
#' @param parameter See [hypothesis_test()] and [get_growthparameters()]
#'   for details.
#' @param by See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param comparison_args See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param hypothesis_args See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param comparison_by See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param hypothesis_by See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param evaluate_comparison See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param evaluate_hypothesis See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param rope_test See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param pd_test See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param comparison_range_null See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param hypothesis_range_null See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param comparison_range See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param hypothesis_range See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param comparison_null See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param hypothesis_null See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param method See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param as_p See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param remove_na See [hypothesis_test()] and [get_growthparameters()]
#'   for details.
#' @param rvar_col See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param conf_level See [hypothesis_test()] and [get_growthparameters()]
#'   for details.
#' @param probs See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param pd_as_percent See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param rope_as_percent See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param nthreads See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param estimate_center See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param estimate_interval See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param na.rm See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param full_frame See [hypothesis_test()] and [get_growthparameters()]
#'   for details.
#' @param return_draws See [hypothesis_test()] and [get_growthparameters()]
#'   for details.
#' @param get_range_null_form See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param get_range_null_value See [hypothesis_test()] and
#'   [get_growthparameters()] for details.
#' @param digits See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param format See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#' @param verbose See [hypothesis_test()] and [get_growthparameters()] for
#'   details.
#'
#' @returns A data frame
#' 
#' @keywords internal
#' @noRd
#' 
get_comparison_hypothesis <- function(data,
                                      full.args = NULL,
                                      parameter = c('apgv', 'pgv'),
                                      by = NULL,
                                      comparison_args = NULL,
                                      hypothesis_args = NULL,
                                      comparison_by = NULL,
                                      hypothesis_by = NULL,
                                      evaluate_comparison = TRUE,
                                      evaluate_hypothesis = TRUE,
                                      rope_test = TRUE,
                                      pd_test = TRUE,
                                      comparison_range_null = NULL,
                                      hypothesis_range_null = NULL,
                                      comparison_range = NULL,
                                      hypothesis_range = NULL,
                                      comparison_null = NULL,
                                      hypothesis_null = NULL,
                                      method = "direct",
                                      as_p = FALSE,
                                      remove_na = TRUE,
                                      rvar_col = NULL,
                                      conf_level = NULL,
                                      probs = c(0.025, 0.975),
                                      pd_as_percent = TRUE,
                                      rope_as_percent = TRUE,
                                      nthreads = NULL,
                                      estimate_center = NULL,
                                      estimate_interval = NULL,
                                      na.rm = TRUE,
                                      full_frame = FALSE,
                                      return_draws = FALSE,
                                      get_range_null_form = FALSE,
                                      get_range_null_value = FALSE,
                                      digits = 2,
                                      parms_sat_elements = NULL,
                                      format = FALSE,
                                      verbose = FALSE) {
  
  string_sat <- 'sat'
  if(!is_emptyx(parms_sat_elements)) {
    get_parms_size     <- parms_sat_elements[['get_parms_size']]    
    parameter_sat      <- parms_sat_elements[['parameter_sat']]
    string_sat         <- parms_sat_elements[['string_sat']] 
    numeric_sat        <- parms_sat_elements[['numeric_sat']]       
    string_numeric_sat <- parms_sat_elements[['string_numeric_sat']] 
  }
 
  
  data <- clean_draws(data,
                      variable = "draw", 
                      group = "drawid", 
                      verbose = verbose)

  if(is.null(full.args)) {
    if(evaluate_comparison) {
      if(is.null(comparison_args)) {
        stop2c("when evaluate_comparison = TRUE, 
               full.args or comparison_args are required")
      }
    }
    if(evaluate_hypothesis) {
      if(is.null(hypothesis_args)) {
        stop2c("when evaluate_comparison = TRUE, 
               full.args or hypothesis_args are required")
      }
    }
  } else if(!is.null(full.args)) {
    if(is.null(comparison_args)) {
      if(is.null(full.args[['comparison_by']])) {
        full.args[['comparison_by']] <- by
      }
      comparison_args <- list()
      if(evaluate_comparison) {
        comparison_args_names <- c('comparison_by')
        for (i in comparison_args_names) {
          comparison_args[[i]] <- full.args[[i]]
        }
        if(is_emptyx(comparison_args)) comparison_args <- NULL
      }
    } else if(!is.null(comparison_args)) {
      if(is.null(comparison_args[['by']])) {
        if(!is.null(comparison_by)) {
          comparison_args[['by']] <- comparison_by
        } else if(!is.null(by)) {
          comparison_args[['by']] <- by
        }
      }
    } # if(is.null(comparison_args)) { else if(!is.null(comparison_args)) {
    if(is.null(hypothesis_args)) {
      if(is.null(full.args[['hypothesis_by']])) {
        full.args[['hypothesis_by']] <- by
      }
      hypothesis_args <- list()
      if(evaluate_hypothesis) {
        hypothesis_args_names <- c('hypothesis_by', 'hypothesis')
        hypothesis_args <- list()
        for (i in hypothesis_args_names) {
          hypothesis_args[[i]] <- full.args[[i]]
        }
        if(is_emptyx(hypothesis_args)) hypothesis_args <- NULL
      }
    } else if(is.null(hypothesis_args)) {
      if(!is.null(hypothesis_by)) {
        hypothesis_args[['by']] <- hypothesis_by
      } else if(!is.null(by)) {
        hypothesis_args[['by']] <- by
      }
    } # if(is.null(hypothesis_args)) { else if(is.null(hypothesis_args)) {
  } # if(is.null(full.args)) { else if(!is.null(full.args)) {
  

  if(!is.null(comparison_args)) {
    if(!is.null(comparison_args[['comparison_by']])) {
      if(is.null(comparison_args[['by']])) {
        comparison_args[['by']] <- comparison_args[['comparison_by']]
      }
      comparison_args[['comparison_by']] <- NULL
    }
  }
  if(!is.null(hypothesis_args)) {
    if(!is.null(hypothesis_args[['hypothesis_by']])) {
      if(is.null(hypothesis_args[['by']])) {
        hypothesis_args[['by']] <- hypothesis_args[['hypothesis_by']]
      }
      hypothesis_args[['hypothesis_by']] <- NULL
    } 
  }
  
  
  collapse::set_collapse(nthreads = parallel::detectCores() - 1)
  
  data.table::setDTthreads(threads = parallel::detectCores() - 1) 
  
  ##########################################
  
  if(is.null(probs) & is.null(conf_level)) {
    stop2c("Please specify either 'probs' or 'conf_level'")
  } else if(!is.null(conf_level)) {
    conf <- ci <- conf_level
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  } else if(!is.null(probs)) {
    conf <- ci <- conf_level <- probs[2] - probs[1] 
    probs <- probs
  }
  
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', probtitles)
  
  
  ##########################################
  
  if(!is.null(estimate_center)) {
    ec_ <- getOption("marginaleffects_posterior_center")
    options("marginaleffects_posterior_center" = estimate_center)
    on.exit(options("marginaleffects_posterior_center" = ec_), add = TRUE)
  }
  if(!is.null(estimate_interval)) {
    ei_ <- getOption("marginaleffects_posterior_interval")
    options("marginaleffects_posterior_interval" = estimate_interval)
    on.exit(options("marginaleffects_posterior_interval" = ei_), add = TRUE)
  }
  ec_agg <- getOption("marginaleffects_posterior_center")
  ei_agg <- getOption("marginaleffects_posterior_interval")
  if(is.null(ec_agg)) ec_agg <- "mean"
  if(is.null(ei_agg)) ei_agg <- "eti"
  
  ##########################################
  
  if(get_range_null_form & get_range_null_value) {
    stop2c(
      "Specify either 'get_range_null_form' or 'get_range_null_value', not both. 
     While 'get_range_null_form' returns the expected structure for the ROPE 
     range and the null for p_direction, the 'get_range_null_value' return the  
     actuall values that will be used for comparison and hypothesis testing")
  }
  
  ##########################################
  
  comparison_equivalence_test_arg <- list()
  comparison_equivalence_test_arg[['ci']] <- ci
  comparison_equivalence_test_arg[['rvar_col']] <- rvar_col
  comparison_equivalence_test_arg[['verbose']] <- verbose
  
  comparison_p_direction_arg <- list()
  comparison_p_direction_arg[['method']] <- method
  comparison_p_direction_arg[['as_p']] <- FALSE
  comparison_p_direction_arg[['remove_na']] <- TRUE
  comparison_p_direction_arg[['rvar_col']] <- NULL
  
  ##########################################
  
  hypothesis_equivalence_test_arg <- comparison_equivalence_test_arg
  hypothesis_p_direction_arg <- comparison_p_direction_arg
  
  ##########################################

  # If somehow hypothesis_args are set by hypothesis itself is NULL
  if(!is.null(hypothesis_args)) {
    if(is.null(hypothesis_args[['hypothesis']])) {
      evaluate_hypothesis  <- FALSE
    }
  }
  
  get_depth <- function(x) {
    if (is.list(x) && length(x) > 0) {
      1 + max(sapply(x, get_depth))
    } else {
      0
    }
  }
  
  list_range_null_to_df <- function(range_null, 
                                    parameter,
                                    data = NULL, 
                                    by = NULL, 
                                    evaluate_comparison = NULL, 
                                    evaluate_hypothesis = NULL, 
                                    comparison_args = NULL,
                                    hypothesis_args = NULL,
                                    what = NULL) {

    
    if(is.list(range_null)) {
      range_null <- range_null
    } else if(is.numeric(range_null)) {
      range_null_temp <- list()
      for (i in parameter) {
        range_null_temp[[i]] <- range_null
      }
      range_null <- range_null_temp
    } else if(is.null(range_null)) {
      range_null <- "default"
    }
    
  
    if(get_depth(range_null) > 1) {
      range_null <- range_null
    } else {
      range_null <- list(range_null)
    }
   
    create_range_lists_pair_args <- list()
    create_range_lists_pair_args[['data']] <- data
    create_range_lists_pair_args[['by']] <- by
    create_range_lists_pair_args[['set_grid']] <- T
    create_range_lists_pair_args[['full_frame']] <- FALSE
    create_range_lists_pair_args[['evaluate_comparison']] <- evaluate_comparison
    create_range_lists_pair_args[['evaluate_hypothesis']] <- evaluate_hypothesis
    create_range_lists_pair_args[['comparison_args']] <- comparison_args
    create_range_lists_pair_args[['hypothesis_args']] <- hypothesis_args
    create_range_lists_pair_args[['parms_sat_elements']] <- parms_sat_elements
    if(!is.null(what)) {
      if(!is.null(range_null)) {
        create_range_lists_pair_args[['what']] <- what
        create_range_lists_pair_args <- c(range_null, 
                                          create_range_lists_pair_args)
        out <- do.call(get_test_range_null, create_range_lists_pair_args)
        return(out)
      } else if(is.null(range_null)) {
        return(NULL)
      }
    } else if(is.null(what)) {
      create_range_lists_pair_args[['what']] <- 'null'
      create_range_lists_pair_args <- c(range_null, 
                                        create_range_lists_pair_args)
      
      test_null <- do.call(get_test_range_null, 
                                      create_range_lists_pair_args) 
      create_range_lists_pair_args[['what']] <- 'range'
      create_range_lists_pair_args <- c(range_null, 
                                        create_range_lists_pair_args)
      test_range <- do.call(get_test_range_null, 
                           create_range_lists_pair_args) 
      test_null_range <- join_df_or_lists(test_null, 
                                          test_range, 
                                          join_on = NULL, 
                                          remove_duplicate = "both")
      return(test_null_range)
    } # if(!is.null(what)) { else if(is.null(what)) {
    
  } # End of list_range_null_to_df
  
  list_comparison_range <- list_comparison_null <- NULL
  list_hypothesis_range <- list_hypothesis_null <- NULL
  
  if(NullFALSE(evaluate_comparison)) {
    if(!is.null(comparison_args) & !is.null(comparison_args$range)) {
      list_comparison_range <- comparison_args$range
      list_comparison_null  <- comparison_args$null
    } else if(!is.null(full.args$equivalence_test)) {
      if(!is.null(full.args$equivalence_test$range)) {
        list_comparison_range <- full.args$equivalence_test$range
      } else {
        list_comparison_range <- NULL
      }
      if(!is.null(full.args$p_direction$null)) {
        list_comparison_null <- full.args$p_direction$null
      } else {
        list_comparison_null <- NULL
      }
    } # if(!is.null(comparison_args) & !is.null(comparison_args$range)) { else .
  }# if(NullFALSE(evaluate_comparison)) {
  
  if(NullFALSE(evaluate_hypothesis)) {
    if(!is.null(hypothesis_args) & !is.null(hypothesis_args$range)) {
      list_hypothesis_range <- hypothesis_args$range
      list_hypothesis_null  <- hypothesis_args$null
    } else if(!is.null(full.args$equivalence_test)) {
      if(!is.null(full.args$equivalence_test$range)) {
        list_hypothesis_range <- full.args$equivalence_test$range
      } else {
        list_hypothesis_range <- NULL
      }
      if(!is.null(full.args$p_direction$null)) {
        list_hypothesis_null <- full.args$p_direction$null
      } else {
        list_hypothesis_null <- NULL
      }
      # list_hypothesis_range <- full.args$equivalence_test$range
      # list_hypothesis_null <- full.args$p_direction$null
    } # if(!is.null(hypothesis_args) & !is.null(hypothesis_args$range)) { else .
  } # if(NullFALSE(evaluate_hypothesis)) {
  
  if(is.null(comparison_range)) {
    if(!is.null(list_comparison_range)) {
      comparison_range <- list_range_null_to_df(list_comparison_range,
                                                parameter = parameter,
                                                data = as.data.frame(data),
                                                by = comparison_args[['by']],
                                                evaluate_comparison = 
                                                  evaluate_comparison,
                                                evaluate_hypothesis = FALSE,
                                                comparison_args = 
                                                  comparison_args,
                                                hypothesis_args = NULL,
                                                what = 'range')
    } else {
      comparison_range <- NULL
    }
  }
  if(is.null(comparison_null)) {
    if(!is.null(list_comparison_null)) {
      comparison_null <- list_range_null_to_df(list_comparison_null,
                                               parameter = parameter,
                                               data = as.data.frame(data),
                                               by = comparison_args[['by']],
                                               evaluate_comparison = 
                                                 evaluate_comparison,
                                               evaluate_hypothesis = FALSE,
                                               comparison_args = 
                                                 comparison_args,
                                               hypothesis_args = NULL,
                                               what = 'null')
    } else {
      comparison_null <- NULL
    }
  }
  
  
  
  if(is.null(hypothesis_range)) {
    if(!is.null(list_hypothesis_range)) {
      hypothesis_range <- list_range_null_to_df(list_hypothesis_range,
                                                parameter = parameter,
                                                data = as.data.frame(data),
                                                by = hypothesis_args[['by']],
                                                evaluate_comparison = FALSE,
                                                evaluate_hypothesis = 
                                                  evaluate_hypothesis,
                                                comparison_args = NULL,
                                                hypothesis_args = 
                                                  hypothesis_args,
                                                what = 'range')
    } else {
      hypothesis_range <- NULL
    }
  }
  if(is.null(hypothesis_null)) {
    if(!is.null(list_hypothesis_null)) {
      hypothesis_null <- list_range_null_to_df(list_hypothesis_null,
                                               parameter = parameter,
                                               data = as.data.frame(data),
                                               by = hypothesis_args[['by']],
                                               evaluate_comparison = FALSE,
                                               evaluate_hypothesis = 
                                                 evaluate_hypothesis,
                                               comparison_args = NULL,
                                               hypothesis_args = 
                                                 hypothesis_args,
                                               what = 'null')
    } else {
      hypothesis_null <- NULL
    }
  }
  

  comparison_range_null <- join_df_or_lists(comparison_null, 
                                            comparison_range, 
                                            join_on = NULL, 
                                            remove_duplicate = "both")
  
  
  hypothesis_range_null <- join_df_or_lists(hypothesis_null, 
                                            hypothesis_range, 
                                            join_on = NULL, 
                                            remove_duplicate = "both")
 

  
  full.args_equivalence_test_range <- full.args$equivalence_test$range

  # get_test_range_null_wrapper_docall_ars_as_list
  
  ####################################################
  # Initiate  create_range_lists_pair_args
  ####################################################
  create_range_lists_pair_args <- list()
  create_range_lists_pair_args[['data']]          <- as.data.frame(data)
  create_range_lists_pair_args[['parameter']]     <- parameter
  create_range_lists_pair_args[['full_frame']]    <- full_frame
  create_range_lists_pair_args[['get_range_null_form']] <- get_range_null_form
  
  
  
  ###########################################################
  # Check - Set comparison and hypothesis -> null range
  ###########################################################
  comparison_test_null <- comparison_test_range <- NULL
  hypothesis_test_null <- hypothesis_test_range <- NULL
  
  if(evaluate_comparison) {
    create_range_lists_pair_args[['by']]  <- comparison_args[['by']]
    create_range_lists_pair_args[['comparison_args']]     <- comparison_args
    create_range_lists_pair_args[['evaluate_comparison']] <- evaluate_comparison
    create_range_lists_pair_args[['hypothesis_args']]     <- NULL
    create_range_lists_pair_args[['evaluate_hypothesis']] <- FALSE
    create_range_lists_pair_args[['parms_sat_elements']] <- parms_sat_elements
    if(pd_test) {
      create_range_lists_pair_args[['what']] <- 'null'
      comparison_test_null <- do.call(get_test_range_null, 
                                      create_range_lists_pair_args)
    }
    if(NullFALSE(rope_test)) {
      create_range_lists_pair_args[['what']] <- 'range'
      comparison_test_range <- do.call(get_test_range_null, 
                                       create_range_lists_pair_args)
    }
  } # if(evaluate_comparison) {
  
  


  if(evaluate_hypothesis) {
    create_range_lists_pair_args[['by']] <- hypothesis_args[['by']]
    create_range_lists_pair_args[['comparison_args']]     <- NULL
    create_range_lists_pair_args[['evaluate_comparison']] <- FALSE
    create_range_lists_pair_args[['hypothesis_args']]     <- hypothesis_args
    create_range_lists_pair_args[['evaluate_hypothesis']] <- evaluate_hypothesis
    if(pd_test) {
      create_range_lists_pair_args[['what']] <- 'null'
      hypothesis_test_null <- do.call(get_test_range_null, 
                                      create_range_lists_pair_args)
    }

    if(NullFALSE(rope_test)) {
      create_range_lists_pair_args[['what']] <- 'range'
      hypothesis_test_range <- do.call(get_test_range_null, 
                                       create_range_lists_pair_args)
    }
  } # if(evaluate_hypothesis) {
  
  

  
  comparison_test_null_range <- hypothesis_test_null_range <- NULL
  comparison_hypothesis_test_null <- list()
  if(!is.null(comparison_test_null) | !is.null(comparison_test_range)) {
    comparison_test_null_range <- join_df_or_lists(comparison_test_null, 
                                                   comparison_test_range, 
                                                   join_on = NULL, 
                                                   remove_duplicate = "both")
    comparison_hypothesis_test_null[['comparison']]<-comparison_test_null_range
  }
  if(!is.null(hypothesis_test_null) | !is.null(hypothesis_test_range)) {
    hypothesis_test_null_range <- join_df_or_lists(hypothesis_test_null, 
                                                   hypothesis_test_range, 
                                                   join_on = NULL, 
                                                   remove_duplicate = "both")
    comparison_hypothesis_test_null[['hypothesis']]<-hypothesis_test_null_range
  }
  
  if(get_range_null_form) {
    return(comparison_hypothesis_test_null)
  } else if(get_range_null_value) {
    return(comparison_hypothesis_test_null)
  }
  
  
  #############################################################################
  
  check_names_exits <- function(data, names) {
    if(is.null(data)) return(invisible(NULL)) # This when no hypothesis 
    `%chin%` <- data.table::`%chin%`
    if(data.table::is.data.table(data)) {
      if(!all(names %chin% names(data))) {
        stop2c("The following variables are missing the data ",
               collapse_comma(setdiff(names , names(data))))
      }
    } else if(is.data.frame(data)) {
      if(!all(names %in% names(data))) {
        stop2c("The following variables are missing the data ",
               collapse_comma(setdiff(names , names(data))))
      }
    } else {
      stop("data should be a data frame")
    }
    return(invisible(NULL))
  }
  
  
  # Function that triggers the error
  check_range_null_structure_rows <- function(out_range_null, set_range_null) {
    
    if(is.null(set_range_null)) return(invisible(NULL))
    if(nrow(out_range_null) != nrow(set_range_null)) {
      set_range_null <- set_range_null[, 'null' := NA]
      set_range_null <- set_range_null[, 'range' := NA]
      set_range_null <- set_range_null %>% data.frame()
      formatted_df <- paste(utils::capture.output(print(set_range_null)), 
                            collapse = "\n")
      stop("The number of rows for specified range and null should be ", 
           nrow(set_range_null), " but is currently ",nrow(out_range_null), ".",
           " Below is the expected structure.",
           " Please use similar structure to set values for range and null.",
           "\n\n",
           formatted_df, call. = FALSE)
    }
  } # check_range_null_structure_rows
  
  
  set_test_null_range_fun <- function(range_null, 
                                      range, 
                                      null, 
                                      set_range_null,
                                      rope_test,
                                      pd_test) {
    check_names_null <- c('null')
    check_names_range <- c('range')
    if(NullFALSE(rope_test) & NullFALSE(pd_test)) {
      check_names_range_null <- c('range', 'null')
    } else if(NullFALSE(rope_test)) {
      check_names_range_null <- c('range')
    } else if(NullFALSE(pd_test)) {
      check_names_range_null <- c('null')
    }
    
    # if(is_emptyx(range_null)) range_null <- NULL
    # 
    # if(is.null(range_null) & is.null(range) & is.null(null)) {
    #   return(NULL)
    # }
    
    if(!is.null(range_null)) {
      check_names_exits(range_null, check_names_range_null)
      out_range_null <- range_null
    } else if(!is.null(range) & is.null(null)) { 
      check_names_exits(range, check_names_range)
      out_range_null <- set_range_null
      out_range_null[['range']] <- range[['range']]
    } else if(is.null(range) & !is.null(null)) { 
      check_names_exits(null, check_names_null)
      out_range_null <- set_range_null
      out_range_null[['null']] <- null[['null']]
    } else if(is.null(range) & is.null(null)) { 
      check_names_exits(set_range_null, check_names_range_null)
      out_range_null <- set_range_null
    }
    check_range_null_structure_rows(out_range_null, set_range_null)
    return(out_range_null)
  } # set_test_null_range_fun
  
  

  
  # Update comparison_test_null_range with user specified range_null/range/null
  if(!is_emptyx(comparison_range_null)) {
    comparison_test_null_range <- 
      set_test_null_range_fun(range_null = comparison_range_null, 
                              range = comparison_range, 
                              null = comparison_null, 
                              set_range_null = comparison_test_null_range,
                              rope_test = rope_test,
                              pd_test = pd_test)
  } # if(!is_emptyx(comparison_range_null)) {
  
  
  if(!is_emptyx(hypothesis_test_null_range)) {
    hypothesis_test_null_range <- 
      set_test_null_range_fun(range_null = hypothesis_range_null, 
                              range = hypothesis_range, 
                              null = hypothesis_null, 
                              set_range_null = hypothesis_test_null_range,
                              rope_test = rope_test,
                              pd_test = pd_test)
  } # if(!is_emptyx(hypothesis_test_null_range)) {
  
  
  
  
  
  ####################################################
  # Check - Set comparison -> range
  ####################################################
  
  comparison_args[['equivalence_test']] <- comparison_equivalence_test_arg
  comparison_args[['p_direction']]      <- comparison_p_direction_arg
  comparison_args[['range_null']]       <- comparison_test_null_range
  
  hypothesis_args[['equivalence_test']] <- hypothesis_equivalence_test_arg
  hypothesis_args[['p_direction']]      <- hypothesis_p_direction_arg
  hypothesis_args[['range_null']]       <- hypothesis_test_null_range
  
  
  call_equivalence_test_p_direction_args <- list()
  call_equivalence_test_p_direction_args[['data']] <- data
  call_equivalence_test_p_direction_args[['by']] <- by
  call_equivalence_test_p_direction_args[['evaluate_comparison']] <- 
    evaluate_comparison
  call_equivalence_test_p_direction_args[['evaluate_hypothesis']] <- 
    evaluate_hypothesis
  call_equivalence_test_p_direction_args[['hypothesis_args']] <- 
    hypothesis_args
  call_equivalence_test_p_direction_args[['comparison_args']] <- 
    comparison_args
  
  call_equivalence_test_p_direction_args[['rope_test']] <- 
    rope_test
  call_equivalence_test_p_direction_args[['pd_test']] <- 
    pd_test
  
  
  call_equivalence_test_p_direction_args[['conf_level']] <- conf_level
  call_equivalence_test_p_direction_args[['probs']] <- probs
  call_equivalence_test_p_direction_args[['pd_as_percent']] <- pd_as_percent
  call_equivalence_test_p_direction_args[['rope_as_percent']] <- rope_as_percent
  call_equivalence_test_p_direction_args[['nthreads']] <- nthreads
  call_equivalence_test_p_direction_args[['ec_agg']] <- ec_agg
  call_equivalence_test_p_direction_args[['ei_agg']] <- ei_agg
  call_equivalence_test_p_direction_args[['na.rm']] <- na.rm
  call_equivalence_test_p_direction_args[['digits']] <- digits
  call_equivalence_test_p_direction_args[['return_draws']] <- return_draws
  call_equivalence_test_p_direction_args[['get_range_null_form']] <- 
    get_range_null_form
  call_equivalence_test_p_direction_args[['digits']] <- digits
  call_equivalence_test_p_direction_args[['verbose']] <- verbose
  

  comparison_hypothesis_results<-do.call(call_equivalence_test_p_direction, 
                                         call_equivalence_test_p_direction_args)
  
  
  if(is_emptyx(comparison_hypothesis_results$comparison)) {
    comparison_hypothesis_results$comparison <- NULL
  }
  if(is_emptyx(comparison_hypothesis_results$hypothesis)) {
    comparison_hypothesis_results$hypothesis <- NULL
  } 
  
  
  if(!is.null(string_sat)) {
    if(!is.null(comparison_hypothesis_results$comparison)) {
      comparison_hypothesis_results$comparison <- 
        rename_vector_in_column_dt(comparison_hypothesis_results$comparison, 
                                           column = 'parameter', 
                                           it = string_sat,
                                           by = string_numeric_sat)
    }
    if(!is.null(comparison_hypothesis_results$hypothesis)) {
      comparison_hypothesis_results$hypothesis <- 
        rename_vector_in_column_dt(comparison_hypothesis_results$hypothesis, 
                                   column = 'parameter', 
                                   it = string_sat,
                                   by = string_numeric_sat)
    }
  } # if(!is.null(string_sat)) {
  
  
  
  # Remove paranthesis () from the hypothesis terms
  if(!is.null(comparison_hypothesis_results$hypothesis)) {
    comparison_hypothesis_results$hypothesis <- 
      comparison_hypothesis_results$hypothesis[, 
                                               hypothesis := 
                                                 gsub("\\(|\\)", "", 
                                                      hypothesis)][,
                                                          hypothesis := 
                                                          trimws(gsub("\\(|\\)", 
                                                          "", hypothesis))]
  }
  
  
  if(format) {
    merge_ranges_eqpd_args <- list()
    merge_ranges_eqpd_args[['x']] <- comparison_hypothesis_results
    merge_ranges_eqpd_args[['string']] <- T
    merge_ranges_eqpd_args[['sep']] <- "; "
    merge_ranges_eqpd_args[['bracket']] <- 'square'
    merge_ranges_eqpd_args[['verbose']] <- FALSE
    comparison_hypothesis_results <- do.call(merge_ranges_eqpd, 
                                             merge_ranges_eqpd_args)
  } # if(format) {
  
  ##########################################
  
  if(!is.null(comparison_hypothesis_results)) {
    if(length(comparison_hypothesis_results) == 1) {
      comparison_hypothesis_results <- comparison_hypothesis_results[[1]]
    }
  }
  
  return(comparison_hypothesis_results)
}


#' get_hypothesis_group_fun
#'
#' @param hypothesis marginaleffect hypothesis argument
#' @keywords internal
#' @noRd
#' 
get_hypothesis_group_fun <- function(hypothesis) {
  sanitize_hypothesis_formula <- get_labels <- NULL;
  getfrom_ <- c('sanitize_hypothesis_formula', 'get_labels')
  for (i in getfrom_) {
    assign(i, utils::getFromNamespace(i, 'marginaleffects'))
  }
  hypothesis_group <- NULL
  if(!is.null(hypothesis)) {
    if(rlang::is_formula(hypothesis)) {
      form <- sanitize_hypothesis_formula(hypothesis)
      hypothesis_group <- form$group
      if(is_emptyx(hypothesis_group)) hypothesis_group <- NULL
    }
  }
  return(hypothesis_group)
}



#' Create Hierarchical Parameter Grids with Path-Aware Inheritance
#'
#' @description
#' Generates fully resolved parameter grids for hierarchical structures (e.g.,
#' Parameter -> Sex -> Study -> Cohort). The function produces two corresponding
#' data frames: one containing the **names** (labels) of the hierarchy levels,
#' and another containing the **resolved numeric values** assigned to each node.
#' Generates fully resolved parameter grids where every value is a **pair**
#' (vector of length 2).
#' - **Single Value**: Automatically converted to a symmetric pair `c(-x, x)`.
#' - **Paired Value**: Preserved exactly as provided `c(min, max)`.
#'
#' It supports advanced value assignment rules, including:
#' \itemize{
#'   \item **Scalar Assignment**: Assigns a single value to all branches of a
#'   node (e.g., `S2 = 20`).
#'   \item **Named List Assignment**: Assigns specific values to specific root
#'   parameters within a branch (e.g., `Male = list(apgv = c(1, 2))`).
#'   Unspecified parameters inherit.
#'   \item **Root-Matching Assignment**: Assigns values matching the number of
#'   root parameters (e.g., `Male = c(0.2, 8)` maps to `apgv` and `pgv`).
#'   \item **Path-Aware Assignment**: Assigns values based on the full path
#'   leading to the current node. For example, if `S1` follows `Male` and
#'   `Female`, providing a vector of length 4 (`c(0.2, 8, 3, 4)`) correctly maps
#'   values to `Male-apgv`, `Male-pgv`, `Female-apgv`, and `Female-pgv`
#'   respectively.
#'   \item **Inheritance**: If a level does not explicitly define a value, it
#'   inherits the value from its parent level. Deepest explicit definitions
#'   override higher ones (e.g., a `Cohort` value overrides a `Study` value).
#' }
#'
#' @param parameter A named list defining the root parameters and their default
#'   values. Example: `list(apgv = 1, pgv = 0)`. These values serve as the
#'   top-level defaults.
#'
#' @param parameter_value (Optional) A numeric vector or list providing values
#'   for the `parameter` argument. Rarely used if `parameter` is already a named
#'   list of values.
#'
#' @param data (Optional) A data frame. If a level is provided as an empty list
#'   (e.g., `study = list()`), the function attempts to infer the level names
#'   (e.g., "S1", "S2") from the corresponding column in this data frame.
#'
#' @param ... Additional stratifying variables representing deeper levels of the
#'   hierarchy (e.g., `sex`, `study`, `cohort`). Each argument should be a named
#'   list where names correspond to the level's items (e.g., "Male", "Female")
#'   and values correspond to the assigned parameters.
#'   \itemize{
#'     \item **Scalar**: `list(C1 = 100)` assigns 100 to C1 regardless of path.
#'     \item **Root-Match Vector**: `list(Male = c(0.2, 8))` assigns 0.2 to the
#'     first root param (`apgv`) and 8 to the second (`pgv`).
#'     \item **Full-Path Vector**: `list(S1 = c(0.2, 8, 3, 4))` assigns values
#'     sequentially based on the expanded grid_str of previous levels (Parameter x
#'     Sex).
#'   }
#'   
#' @param get_grid (Optional) A logical (default \code{FALSE}) to get grid_str 
#' returned from the \code{generate_grid_list} function.
#' 
#' @param set_grid (Optional) A logical (default \code{FALSE}) to set parameter
#'   and parameter values using the grid_str list returned from the
#'   \code{generate_grid_list} function.
#' 
#' @param by (Optional) A character string or a vector that is used in setting
#' up the grid_str via the \code{generate_grid_list} function.
#' 
#' @param evaluate_comparison (Optional)
#' 
#' @param evaluate_hypothesis (Optional)
#'  
#' @param return (Optional) A character string indicating whether to return the
#' grid_str as data frame or a list.
#'   
#' @param parameter_grid_value A named list or logical to set the parameter
#'   values for the \code{'apgv'}, \code{'pgv'}, \code{'spgv'}, \code{'atgv'}
#'   \code{'tgv'}, \code{'stgv'}, \code{'acgv'}, \code{'cgv'}, \code{'acgv'},
#'   when logical \code{TRUE}, then each of the nine parameters are assigned
#'   default values which are: \code{1} for age parameters, \code{o.5} for
#'   velocity parameters, and \code{5.0} for size parameters. when
#'   \code{parameter_grid_value} is a list, it must be a named. whatever names
#'   are missing, those parameters are assigned the default values.
#'   
#' @param full_frame A logical to indicate whether to return full frame of or
#'   only the null and range columns. When both data and by are null i.e., when
#'   list structure is used to set up the null and range columns, then
#'   full_frame should be TRUE, otherwise FALSE. Default \code{full_frame =
#'   NULL} which allows automatic setting of \code{full_frame = TRUE} or
#'   \code{full_frame = FALSE}.
#'   
#' @return A named list containing two data frames:
#' \item{final_grid}{A data frame where columns correspond to the hierarchy
#' levels. Cells contain the **resolved numeric values** for that specific path.
#' This table acts as a "trace" of the values, showing inheritance and
#' overrides.}
#' \item{final_names}{A corresponding data frame where cells contain the
#' **names/labels** (e.g., "apgv", "Male", "S1") for the path.}
#' 
#' @examples
#' # Example: Defining a hierarchy with Parameter -> Sex -> Study -> Cohort
#' # S1 uses path-aware assignment (different values for Male/Female x apgv/pgv)
#' # Cohort overrides everything with fixed values 100 and 200
#'
#' # --- Example 1: Basic Usage with Scalars and Pairs ---
#' # - Scalars (e.g., 0, 1) are converted to symmetric pairs c(0, 0), c(-1, 1).
#' # - Explicit pairs (e.g., c(1, 2)) are preserved.
#' # - Unspecified items inherit from parent.
#' 
#' res <- get_test_range_null(
#'   parameter = list(apgv = c(1, 2), pgv = 0),
#'   sex       = list(Male = NULL, Female = 1), 
#'   study     = list(S1 = 10, S2 = 20),
#'   cohort    = list(C1 = 100, C2 = 200)
#' )
#' print(res)
#' 
#' 
#' # --- Example 2: Named List Assignment (Granular Control) ---
#' # Here, we explicitly set 'apgv' for 'Male' using a named list.
#' # 'pgv' for 'Male' is NOT specified in the list, so it inherits 
#' # from the 'parameter' level (pgv = 0).
#' 
#' # Note that for sex = list(Male = list(apgv = c(4, 9)), Female = 1), the
#' # c(4, 9) is treated as pair for apgv whereas the pgv is assigned the default
#' # value. However, as shown in the next example below, the
#' # list(Male = c(0.2, 8), Female = c(2, 1)), the c(0.2, 8) is considered as  
#' # values specified for apgv (0.2) and pgv(8) respectively and hence it will  
#' # be translated into apgv = c(-0.2, 0.2), and pgv = c(-8, 8)
#' 
#' res_named <- get_test_range_null(
#'   parameter = list(apgv = c(1, 2), pgv = 0),
#'   sex       = list(Male = list(apgv = c(4, 9)), Female = 1),
#'   study     = list(S1 = 10, S2 = 20),
#'   cohort    = list(C1 = 100, C2 = 200)
#' )
#' print(res_named)
#' 
#' 
#' # --- Example 3: Path-Aware Vector Assignment ---
#' # S1 is assigned a vector of length 4: c(0.2, 8, 3, 4).
#' # This maps sequentially to the expansion of previous levels (Parameter x Sex).
#' # Order: (apgv-Male, pgv-Male, apgv-Female, pgv-Female)
#' 
#' res_path <- get_test_range_null(
#'   parameter = list(apgv = 1, pgv = 0),
#'   sex       = list(Male = c(0.2, 8), Female = c(2, 1)),
#'   study     = list(S1 = c(0.2, 8, 3, 4), S2 = 20),
#'   cohort    = list(C1 = 100, C2 = 200)
#' )
#' print(res_path)
#' 
#' 
#' # --- Example 4: Handling NULLs for Inheritance ---
#' # When a level is explicitly set to NULL (or omitted in list), 
#' # it inherits the value from the level above it.
#' 
#' res_null <- get_test_range_null(
#'   parameter = list(apgv = 1, pgv = 0),
#'   sex       = list(Male = c(0.2, 8), Female = c(2, 1)),
#'   study     = list(S1 = 10, S2 = NULL), # S2 inherits from Sex
#'   cohort    = list(C1 = NULL, C2 = NULL) # Cohorts inherit from Study
#' )
#' print(res_null)
#' 
#' 
#' # --- Example 5: Using Grid Getter and Setter ---
#' 
#' # 5a. Generate the underlying grid_str structure without values
#' df_mock <- data.frame(parameter="apgv", sex="Male", study="S1", cohort="C1")
#' by <- c("parameter", "sex", "study", "cohort")
#' 
#' the_grid <- get_test_range_null(
#'   get_grid = TRUE,
#'   by = by,
#'   data = df_mock
#' )
#' print(the_grid)
#' 
#' # Use the generated grid_str structure to initialize a new call (set_grid)
#' # This creates the structure defined in 'the_grid' (as list) 
#' # and allows populating it via '...'.
#' 
#' the_sgrid <- get_test_range_null(
#'   set_grid = TRUE,
#'   by = by,
#'   data = df_mock,
#'   # Now we can provide values for the levels defined in the grid_str
#'   parameter = list(apgv = 1),
#'   sex = list(Male = 2)
#' )
#' print(the_sgrid)
#'
#' @keywords internal
#' @noRd
#' 
get_test_range_null <- function(parameter = NULL,
                                parameter_value = NULL, 
                                data = NULL, 
                                ...,
                                set_grid = FALSE,
                                get_grid = FALSE,
                                get_range_null_form = FALSE,
                                get_range_null_value = FALSE,
                                comparison_args = NULL,
                                hypothesis_args = NULL,
                                evaluate_comparison = NULL,
                                evaluate_hypothesis = NULL,
                                what = 'range',
                                by = NULL,
                                return = 'df',
                                parameter_grid_value = NULL, 
                                full_frame = NULL,
                                parms_sat_elements = NULL, 
                                verbose = FALSE) {
  
  # Default, will get updated by parms_sat_elements[['string_sat']] 
  string_sat <- 'sat'
  if(!is_emptyx(parms_sat_elements)) {
    get_parms_size     <- parms_sat_elements[['get_parms_size']]    
    parameter_sat      <- parms_sat_elements[['parameter_sat']]
    string_sat         <- parms_sat_elements[['string_sat']] 
    numeric_sat        <- parms_sat_elements[['numeric_sat']]       
    string_numeric_sat <- parms_sat_elements[['string_numeric_sat']] 
  } 
  
  
  if(!is.null(by)) {
    if(is.logical(by)) by <- NULL
  }
  
  
  if(is.null(what)) what <- 'range'
  checkmate::assert_choice(what, choices = c('range', 'null'), null.ok = FALSE)
  range <-  null <- FALSE  
  if(what == 'range') range <- TRUE 
  if(what == 'null')  null  <- TRUE
  if(!range & !null) {
    if(verbose) {
      message2c("'null' set as TRUE")
    }
  }
  if(range & null) {
    stop2c("Please set either 'range' or 'null' as TRUE, not both")
  }
  if(range) setpair <- TRUE
  if(null)  setpair <- FALSE

  # when both data and by are NULL, it means all list elements used 
  if(is.null(data) & is.null(by)) {
    call_standalone_set_grid_fun <- TRUE
    set_grid <- TRUE
    call_evaluate_hypothesis_fun <- call_evaluate_comparison_fun <- FALSE
    if(is.null(full_frame)) full_frame <- TRUE
  } # if(is.null(data) & is.null(data)) {
  
  # when either data or by specified
  if(!is.null(data) | !is.null(by)) {
    if(get_range_null_form) {
      get_grid <- TRUE
    }
    
    # call_standalone_set_grid_fun <- FALSE
    if(is.null(by)) {
      call_standalone_set_grid_fun <- TRUE
    } else if(!is.null(by)) {
      call_standalone_set_grid_fun <- FALSE
    }
    
    if(is.null(full_frame)) full_frame <- FALSE
    call_evaluate_hypothesis_fun <- call_evaluate_comparison_fun <- FALSE
    if(is_emptyx(comparison_args)) comparison_args <- NULL
    if(is_emptyx(hypothesis_args)) hypothesis_args <- NULL
    if(is.null(comparison_args) & is.null(hypothesis_args)) {
      if(get_range_null_form) {
        stop2c("When 'get_range_null_form' = TRUE, then either comparison_args 
           or hypothesis_args must be specified")
      }
    }
    if(!is.null(comparison_args)) {
      if(!is.list(comparison_args)) stop2c("'comparison_args' must be a list")
      call_evaluate_comparison_fun <- TRUE
    }
    if(!is.null(hypothesis_args)) {
      if(!is.list(hypothesis_args)) stop2c("'hypothesis_args' must be a list")
      call_evaluate_hypothesis_fun  <- TRUE
      # If somehow hypothesis_args are set by hypothesis itself is NULL
      if(is.null(hypothesis_args[['hypothesis']])) {
        call_evaluate_hypothesis_fun  <- FALSE
      }
    } # if(!is.null(hypothesis_args)) {
    
   
    if(NullFALSE(evaluate_comparison)) {
      call_evaluate_comparison_fun  <- TRUE
    }
    if(NullFALSE(evaluate_hypothesis)) {
      call_evaluate_hypothesis_fun  <- TRUE
    }
    
    
    if(get_grid | set_grid) {
      if(is.null(data)) stop2c("'data' required when get_grid = TRUE")
    } else if(!get_grid & !set_grid) {
      if(is.null(parameter)) stop2c("'parameter' 
                                    required when get_grid = FALSE")
    }
    
    
    if(!is.null(by)) {
      if(!is.logical(by)) {
        if(!is.character(by)) {
          stop2c("'by' must be a character vector")
        }
      }
    }
  } # if(!is.null(data) & !is.null(data)) {
  
  
  
  # if(!is.null(hypothesis_args)) {
  #   if(!is.list(hypothesis_args)) stop2c("'hypothesis_args' must be a list")
  #   call_evaluate_hypothesis_fun  <- TRUE
  # }
  
  if(is.null(parameter)) {
    if(is.null(parameter_value)) {
      if(!is.null(parameter_grid_value)) {
        if(!is.list(parameter_grid_value)) {
          stop("'parameter_grid_value' must be a named list")
        }
        parameter <- names(parameter_grid_value)
        # parameter_value <- parameter_grid_value
        #parameter_grid_value <- NULL
      }
    }
  } else if(!is.null(parameter)) {
    if(is.list(parameter)) {
      parameter_grid_value <- parameter
      parameter <- names(parameter)
    } # if(is.list(parameter)) {
  } # if(is.null(parameter)) { else if(!is.null(parameter)) {
  
  
  
  # if parameter <- c('apgv', 'pgv') and no parameter_value/parameter_grid_value
  only_parameter_no_parameter_or_grid_value <- FALSE
  if(!is.list(parameter)) {
    if(is.vector(parameter)) {
      checkmate::assert_character(parameter)
      if(is.null(parameter_value) & is.null(parameter_grid_value)) {
        only_parameter_no_parameter_or_grid_value <- TRUE
        if(!get_grid & !set_grid) set_grid <- TRUE
      }
    } # if(is.vector(parameter)) {
  } # if(!is.list(parameter)) {
  
  
  
  get_depth <- function(x) {
    if (is.list(x) && length(x) > 0) {
      1 + max(sapply(x, get_depth))
    } else {
      0
    }
  }
  
  copy_list_info <- function(dots, dots_in) {
     if(is_emptyx(dots)) return(dots_in) # for call_standalone_set_grid_fun
    dots[names(dots_in)] <- mapply(
      utils::modifyList,                  # or purrr::list_modify
      dots[names(dots_in)],
      dots_in,
      SIMPLIFY = FALSE
    )
    return(dots)
  }
  
  
  get_eqpd_form_rename <- function(x, setpair) {
    if(!is.null(x)) {
      x <- data.table::setDF(x)
      if(setpair)  data.table::setnames(x, 'estimate', 'range')
      if(!setpair) data.table::setnames(x, 'estimate', 'null')
    }
    return(x)
  }
 
  
  sanitize_dots_in <- function(dots_in) {
    rm_dots_in <- c('evaluate_comparison', 'evaluate_hypothesis')
    for (i in rm_dots_in) {
      dots_in[[i]] <- NULL
    }
    return(dots_in)
  }
 
  hypothesis_group_by_msg <- 
    paste0("'hypothesis_group' should be a subset of 'by' variables. 
           Furthermore, make sure that 'hypothesis_group' is nested. E.g., 
           for by=c('sex', 'id'), the 'hypothesis_group' should be 'sex' and  
             not 'id' i.e., ... | sex")
  
  
  # Even when set_grid, grid_str_hypothesis is needed
  if(call_evaluate_hypothesis_fun) {
    hypothesis_group<-get_hypothesis_group_fun(hypothesis_args[['hypothesis']])
    if(!is.null(hypothesis_group)) {
      hypothesis_group_by_str <- paste0("hypothesis_group =  ", 
                                        collapse_comma(hypothesis_group),
                                        " and by = ", 
                                        collapse_comma(by))
      hypothesis_group_by_msg <- paste0(hypothesis_group_by_msg, ". Currently ",
                                        hypothesis_group_by_str)
      if(identical(hypothesis_group, by)) {
        stop2c(hypothesis_group_by_msg)
      } else if(any(!hypothesis_group %in% by)) {
        stop2c(hypothesis_group_by_msg)
      } else if(length(hypothesis_group) >= length(by)) {
        stop2c(hypothesis_group_by_msg)
      }
    } # if(!is.null(hypothesis_group)) {
    grid_by <- c('parameter', by)
    grid_str <- generate_grid_list(df = data, by = grid_by, return = return)
    grid_str_hypothesis <-  evaluate_hypothesis_fun(data = grid_str, 
                                                    hypothesis_args = 
                                                      hypothesis_args,
                                                    set_grid = TRUE,
                                                    get_range_null_form = TRUE)

    converfactor <- c('hypothesis')
    grid_str_hypothesis <- grid_str_hypothesis[, 
                                               (converfactor) := 
                                                 lapply(.SD, as.factor), 
                                               .SDcols = converfactor]
    
    
  } # if(call_evaluate_hypothesis_fun) {
  
  
  
  
  if(get_grid) {
    grid_by <- c('parameter', by)
    grid_str <- generate_grid_list(df = data, by = grid_by, return = return)
    if(get_range_null_form) {
      grid_str_list <- list()
      if(call_evaluate_comparison_fun) {
        grid_str_comparison <-  evaluate_comparison_fun(data = grid_str, 
                                             comparison_args = comparison_args,
                                             set_grid = FALSE,
                                             get_range_null_form = 
                                               get_range_null_form)
        grid_str_comparison <- get_eqpd_form_rename(grid_str_comparison, 
                                                    setpair)
        grid_str_list[['comparison']] <- grid_str_comparison
      }
      if(call_evaluate_hypothesis_fun) {
        grid_str_hypothesis <-  evaluate_hypothesis_fun(data = grid_str, 
                                             hypothesis_args = hypothesis_args,
                                             set_grid = FALSE,
                                             get_range_null_form = 
                                               get_range_null_form)
        grid_str_hypothesis <- get_eqpd_form_rename(grid_str_hypothesis, 
                                                    setpair)
        grid_str_list[['hypothesis']] <- grid_str_hypothesis
      }
      return(grid_str_list)
    } # if(get_range_null_form) {
    if(!data.table::is.data.table(grid_str)) {
      if(setpair)  grid_str <- grid_str %>% dplyr::mutate('range' = NA)
      if(!setpair) grid_str <- grid_str %>% dplyr::mutate('null' = NA)
    }
    if(data.table::is.data.table(grid_str)) {
      if(setpair)  grid_str <- grid_str[, 'range' := NA]
      if(!setpair) grid_str <- grid_str[, 'null' := NA]
    }
   return(grid_str)
  } else if(set_grid) {
    if(call_evaluate_comparison_fun) {
      grid_by <- c('parameter', by)
      grid_str <- generate_grid_list(df = data, by = grid_by, return = 'list')
    }
    if(call_evaluate_hypothesis_fun) {
      grid_by <- c('parameter', 'hypothesis', hypothesis_group)
      grid_str <- generate_grid_list(df = grid_str_hypothesis, by = grid_by, 
                                     return = 'list')
    }
    if(call_standalone_set_grid_fun) {
      grid_by <- c('parameter', by)
      grid_str <- generate_grid_list(df = data, by = grid_by, 
                                     parameter = parameter, return = 'list')
    }
    parameter <- grid_str$parameter
    nonparameter <- grid_str
    nonparameter$parameter <- NULL
    if(!is.null(parameter_value) & is.null(parameter_grid_value)) {
      parameter_grid_value <- parameter_value
    }
    dots_in <- list(...)
    dots_in <- sanitize_dots_in(dots_in)
    dots <- nonparameter
    dots <- copy_list_info(dots, dots_in)
  } else {
    dots <- list(...)
  }
  
  
  
  
  
  # ... [Parameter Grid Value Defaults - Same as before] ...
  defaut_parameter_grid_value <- list()
  if(range) {
    defaut_parameter_grid_value[['apgv']] <- 1.0
    defaut_parameter_grid_value[['pgv']]  <- 0.5
    defaut_parameter_grid_value[['spgv']] <- 5.0
    defaut_parameter_grid_value[['atgv']] <- 1.0
    defaut_parameter_grid_value[['tgv']]  <- 0.5
    defaut_parameter_grid_value[['stgv']] <- 5.0
    defaut_parameter_grid_value[['acgv']] <- 1.0
    defaut_parameter_grid_value[['cgv']]  <- 0.5
    defaut_parameter_grid_value[['scgv']] <- 5.0
    defaut_parameter_grid_value[[string_sat]]  <- 5.0
  }
  if(null) {
    defaut_parameter_grid_value[['apgv']] <- 0.0
    defaut_parameter_grid_value[['pgv']]  <- 0.0
    defaut_parameter_grid_value[['spgv']] <- 0.0
    defaut_parameter_grid_value[['atgv']] <- 0.0
    defaut_parameter_grid_value[['tgv']]  <- 0.0
    defaut_parameter_grid_value[['stgv']] <- 0.0
    defaut_parameter_grid_value[['acgv']] <- 0.0
    defaut_parameter_grid_value[['cgv']]  <- 0.0
    defaut_parameter_grid_value[['scgv']] <- 0.0
    defaut_parameter_grid_value[[string_sat]]  <- 0.0
  }
  
  
  allowed_parameter_grid_value_names <- c('apgv', 'pgv', 'spgv', 
                                          'atgv', 'tgv', 'stgv', 
                                          'acgv', 'cgv', 'scgv')
  
  
  if(only_parameter_no_parameter_or_grid_value & !set_grid) {
    parameter_vector <- parameter
    parameter <- as.list(parameter)
    names(parameter) <- parameter_vector
    for (i in names(parameter)) {
      parameter[[i]] <- NULL
    }
    parameter_grid_value <- TRUE
  } else if(only_parameter_no_parameter_or_grid_value & set_grid) {
    parameter_grid_value <- TRUE
  }
  
  
  
  if(!is.null(parameter_grid_value)) {
    if(is.logical(parameter_grid_value)) {
      if(parameter_grid_value) {
        for (i in names(parameter)) {
          parameter[[i]] <- defaut_parameter_grid_value[[i]]
        } 
      } 
    } else { 
      if(!is.list(parameter_grid_value)) {
        stop2c("parameter_grid_value must be a list")
      }
      if(length(parameter_grid_value) != length(names(parameter_grid_value))) {
        stop2c("parameter_grid_value must be a named list")
      }
      for (i in names(parameter_grid_value)) {
        if(!i %in% allowed_parameter_grid_value_names) {
          stop2c(i, " is not a valid parameter name.
                 Allowed parameter names are: ", 
                 collapse_comma(allowed_parameter_grid_value_names))
        }
      }
      for (i in names(parameter)) {
        if(!is.null(parameter_grid_value[[i]])) {
          parameter[[i]] <- parameter_grid_value[[i]]
        } else {
          parameter[[i]] <- defaut_parameter_grid_value[[i]]
        }
      }
    } 
  }
  
  
  # Drop first layer (one level flatten)
  if(get_depth(parameter) > 1) {
    parameter <- unlist(parameter, recursive = FALSE, use.names = TRUE)
    names(parameter) <- sub(".*\\.", "", names(parameter))
  }
  
  if(is_emptyx(parameter)) {
    stop2c("'parameter' values not found")
  }
  
  
  
  # --- Helper: Transform Value to Pair ---
  transform_to_pair <- function(val, item_name, setpair) {
    if(!setpair) {
      if(!is.null(val)) {
        if (length(val) > 1) {
          stop2c("For 'null', the length should be one.", 
                 " Check the below values ", collapse_comma(val),
                 " specified for ", collapse_comma(item_name))
        }
      }
      return(val)
    } # if(setpair) {
    
    if(!is.null(val)) {
      if (length(val) > 2) {
        stop2c("For 'range', the length should be either one or two.", 
               " Check the below values ", collapse_comma(val),
               " specified for ", collapse_comma(item_name))
      }
    }
    
    if (length(val) == 1) {
      return(c(-val, val))
    } else if (length(val) == 2) {
      return(val)
    } else {
      return(val) 
    }
  }
  
  # --- 1. PARSE LEVELS ---
  parse_level <- function(main_input, val_input, label) {
    names_vec <- NULL; vals_list <- list(); is_inherited <- FALSE
    if((is.list(main_input) && length(main_input)==0) || is.null(main_input)) {
      if (is.null(data)) {
        stop2c(sprintf("Argument '%s' is empty but no 'data' frame provided.", 
                       label))
      }
      if (!label %in% names(data)) {
        stop2c(sprintf("Column '%s' not found in 'data'.", label))
      }
      names_vec <- as.character(unique(data[[label]]))
      names_vec <- names_vec[!is.na(names_vec)]
      is_inherited <- TRUE
    } else if (is.list(main_input) && !is.null(names(main_input))) {
      names_vec <- names(main_input); vals_list <- main_input 
    } else {
      names_vec <- main_input
      if (is.null(val_input)) is_inherited <- TRUE 
      else vals_list <- as.list(rep(val_input, length.out = length(names_vec)))
    }
    return(list(names = names_vec, values = vals_list, label = label, 
                inherited = is_inherited))
  }
  
  processed_levels <- list()
  processed_levels[[1]] <- parse_level(parameter, parameter_value, "parameter")
  
  dot_names <- names(dots)
  primary_keys <- dot_names[!grepl("_value$", dot_names)]
  for (key in primary_keys) {
    main_dat <- dots[[key]]
    val_dat  <- dots[[paste0(key, "_value")]]
    processed_levels[[length(processed_levels) + 1]] <- parse_level(main_dat, 
                                                                    val_dat, 
                                                                    key)
  }
  
  # --- 2. CREATE MASTER GRID (Names) ---
  level_names_list <- lapply(processed_levels, function(x) x$names)
  names(level_names_list) <- sapply(processed_levels, function(x) x$label)
  
  full_grid_names <- do.call(expand.grid, c(level_names_list, 
                                            list(stringsAsFactors = FALSE)))
  
  # Create working grid_str for values (Initialize as LIST columns)
  full_grid_vals <- full_grid_names
  for (i in seq_along(processed_levels)) {
    full_grid_vals[[paste0("val_L", i)]] <- vector("list", nrow(full_grid_vals))
  }
  
  # --- 3. POPULATE VALUES PER LEVEL ---
  root_params <- processed_levels[[1]]$names
  
  for (depth in seq_along(processed_levels)) {
    lvl <- processed_levels[[depth]]
    col_name <- paste0("val_L", depth)
    
    if (!lvl$inherited) {
      if (depth == 1) context_grid <- data.frame(dummy = 1)
      else {
        prev_levels <- level_names_list[1:(depth-1)]
        context_grid <- do.call(expand.grid, c(prev_levels, 
                                               list(stringsAsFactors = FALSE)))
      }
      expected_vec_len <- nrow(context_grid)
      
      for (item_name in lvl$names) {
        raw_val <- lvl$values[[item_name]]
        row_indices <- which(full_grid_vals[[lvl$label]] == item_name)
        
        vals_to_assign <- list()
        
        # NEW LOGIC: Named List Handling (e.g., Male = list(apgv = c(1,2)))
        if (is.list(raw_val) && !is.null(names(raw_val))) {
          # The user provided values for specific parameters for this item
          # E.g. raw_val = list(apgv = c(1, 2))
          
          # Initialize all rows for this item as NULL (they inherit by default)
          # (Already NULL by default initialization)
          
          for (param_key in names(raw_val)) {
            # Check if this param_key is valid
            if (!param_key %in% root_params) {
              stop2c("The ", collapse_comma(param_key), 
                     "is not a valid parameter name. The options are: ",
                     collapse_comma(root_params))
            }
            
            val_content <- raw_val[[param_key]]
            pair <- transform_to_pair(val_content, item_name, setpair)
            param_label <- processed_levels[[1]]$label
            
            specific_indices <- 
              row_indices[full_grid_vals[[param_label]][row_indices]==param_key]
            
            if (length(specific_indices) > 0) {
              full_grid_vals[[col_name]][specific_indices] <- 
                rep(list(pair), length(specific_indices))
            }
          }
          # We do NOT assign to `vals_to_assign` generic logic below. 
          # We are done for this item.
          next 
        }
        
        # --- STANDARD LOGIC (Scalar, Vector, Pair) ---
        
        # 1. Scalar (Length 1) -> Convert to Pair c(-x, x)
        if (length(raw_val) == 1) {
          pair <- transform_to_pair(raw_val, item_name, setpair)
          vals_to_assign <- rep(list(pair), length(row_indices))
          
          # 2. Explicit Pair (Length 2) at Root Level (Depth 1) -> Keep as is
        } else if (depth == 1 && length(raw_val) == 2) {
          pair <- raw_val 
          vals_to_assign <- rep(list(pair), length(row_indices))
          
          # 3. Explicit Pair (Length 2) matching Root Param Count -> 
          # Treat as Vector of Scalars
        } else if (length(raw_val) == length(root_params)) {
          n_repeats <- length(row_indices) / length(raw_val)
          base_pairs <- lapply(raw_val, transform_to_pair, item_name, setpair)
          vals_to_assign <- rep(base_pairs, times = n_repeats)
          
          # 4. Explicit Pair (Length 2) but NOT matching Root Param Count -> 
          # Treat as single pair
        } else if (length(raw_val) == 2) {
          pair <- raw_val
          vals_to_assign <- rep(list(pair), length(row_indices))
          
          # 5. Full Vector (Length N) -> Transform elements individually
        } else if (length(raw_val) == expected_vec_len) {
          n_repeats <- length(row_indices) / expected_vec_len
          base_pairs <- lapply(raw_val, transform_to_pair, item_name, setpair)
          vals_to_assign <- rep(base_pairs, times = n_repeats)
        }
        
        if(length(vals_to_assign) > 0) {
          full_grid_vals[[col_name]][row_indices] <- vals_to_assign
        }
      }
    }
  }
  
  # --- 4. COALESCE VALUES ---
  final_grid_vals <- data.frame(row.names = 1:nrow(full_grid_vals))
  
  for (depth in seq_along(processed_levels)) {
    label <- processed_levels[[depth]]$label
    current_col <- full_grid_vals[[paste0("val_L", depth)]]
    
    if (depth > 1) {
      prev_final_col <- final_grid_vals[[processed_levels[[depth-1]]$label]]
      combined_col <- vector("list", length(current_col))
      for(k in seq_along(current_col)) {
        if (!is.null(current_col[[k]])) combined_col[[k]] <- current_col[[k]]
        else combined_col[[k]] <- prev_final_col[[k]]
      }
      current_col <- combined_col
    }
    final_grid_vals[[label]] <- I(current_col)
  }
  
  values_for_parms       <- names(final_grid_vals)
  names(final_grid_vals) <- paste0(values_for_parms, "_", 'value')
  
  combination <- tidyr::unite(full_grid_names,
                              !!as.symbol(names(full_grid_names)),
                              remove = TRUE)
  names(combination) <- 'combination'
  out <- cbind(combination, full_grid_names, final_grid_vals)
  
  if(full_frame) {
    attr(out, 'combination') <- values_for_parms
    return(out)
  } else if(!full_frame) {
    grid_by <- c(grid_by, dot_names)
    grid_by <- unique(grid_by)
    if(data.table::is.data.table(out)) {
      out <- out
    } else if(!data.table::is.data.table(out)) {
      out <- data.table::setDF(out)
    }
    last_col <- ncol(out)
    out <- out[, c(grid_by, names(out)[last_col])]
    if(setpair)  out <- data.table::setnames(out, ncol(out), "range")
    if(!setpair) out <- data.table::setnames(out, ncol(out), "null")
    out <- data.table::setDT(out)[, (grid_by) := lapply(.SD, as.factor), 
                                  .SDcols = grid_by]
    attr(out, 'combination') <- values_for_parms
  } # if(full_frame) { else if(!full_frame) {
  
  
  return(out)
}


##################################################################

#' get_test_range_null_wrapper_docall_ars_as_list
#' @param x vectors, lists, rows, etc
#' @keywords internal
#' @noRd
#' 
get_test_range_null_wrapper_docall_ars_as_list <- function(x) {
  # Works for vectors, lists, rows, etc.
  do.call(get_test_range_null, as.list(x)) 
}


##################################################################
#' NullFALSE
#' @param x vectors, lists, rows, etc
#' @keywords internal
#' @noRd
#' 
NullFALSE <- function(x) {
  
  if(is.null(x)) {
    return(FALSE)
  } else if(is.logical(x)) {
    return(x)
  } else if(is.character(x)) {
    if(x == "F" | x == "FALSE") return(FALSE)
    if(x == "T" | x == "TRUE") return(TRUE)
  } else if(is.list(x)) {
    # list 
    return(TRUE)
    # stop("'x' muste be either NULL or logical")
  }
}

##################################################################

#' Generalized Join for Lists of Data Frames or Single Data Frames
#'
#' This function merges multiple lists of data.frames (or single data.frames)
#' sequentially using a `Map` -> `Reduce` strategy. It supports:
#' \itemize{
#'   \item Merging corresponding elements from multiple lists (e.g.,
#'   \code{list1[[1]] + list2[[1]]}).
#'   \item Merging a single data.frame against a list (recycling the single
#'   data.frame).
#'   \item Merging multiple single data.frames into one result.
#'   \item Automatic or explicit specification of merge columns.
#'   \item Flexible removal of duplicate columns created during joins (e.g.,
#'   \code{i.variable}, \code{variable.1}).
#' }
#'
#' @param ... Input arguments. Can be any combination of:
#'   \itemize{
#'     \item \code{list} of \code{data.frame}s or \code{data.table}s.
#'     \item Single \code{data.frame} or \code{data.table}.
#'   }
#'   If a mix of lists and single data.frames is provided, the single
#'   data.frames are recycled to match the length of the lists.
#'
#' @param join_on Character vector. The column name(s) to join by. If
#'   \code{NULL} (default), the function automatically uses the intersection of
#'   column names between the two tables currently being joined
#'   (\code{intersect(names(x), names(y))}).
#' 
#' @param remove_duplicate Character string specifying how to handle duplicate
#'   columns generated during the join (e.g., \code{i.name} from data.table
#'   joins or \code{name.1} suffixes). Options are:
#'   \itemize{
#'     \item \code{"both"} (Default): Removes duplicates matching both "i."
#'     prefix and ".N" numeric suffix.
#'     \item \code{"string"}: Removes only columns starting with "i." (standard
#'     data.table behavior).
#'     \item \code{"numeric"}: Removes only columns ending in ".1", ".2", etc.
#'     \item \code{"none"}: Retains all duplicate columns.
#'   }
#'
#' @return A \code{list} of \code{data.table}s.
#'   \itemize{
#'     \item If inputs are lists of length N, returns a list of length N
#'     containing the merged results.
#'     \item If inputs are single data.frames, returns a list of length 1
#'     containing the single merged result.
#'   }
#'
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' library(magrittr) # for the pipe %>%
#'
#' # --- Setup Test Data ---
#' # List 1: Range data
#' test_rangex <- list(
#'   data.table(id = 1:3, range = c(10, 20, 30), group = "A"),
#'   data.table(id = 4:6, range = c(40, 50, 60), group = "B")
#' )
#'
#' # List 2: Null data (to be joined)
#' test_nullx <- list(
#'   data.table(id = 1:3, null_val = c(0, 0, 0), extra = "X"),
#'   data.table(id = 4:6, null_val = c(1, 1, 1), extra = "Y")
#' )
#'
#' # --- Example 1: Joining Lists ---
#' # Merges test_rangex[[1]] with test_nullx[[1]], and [[2]] with [[2]]
#' list_result <- join_all_lists(
#'   test_rangex, 
#'   test_nullx, 
#'   join_on = "id"
#' )
#' print(list_result)
#'
#' # --- Example 2: Joining Single Data Frames ---
#' # Merges a single pair of data frames.
#' # Note: We convert to base data.frame to demonstrate the function handles 
#' # them automatically.
#' df1 <- test_rangex[[1]] %>% as.data.frame()
#' df2 <- test_nullx[[1]] %>% as.data.frame()
#'
#' single_result <- join_all_lists(
#'   df1, 
#'   df2, 
#'   join_on = "id"
#' )
#' print(single_result)
#'
#' # --- Example 3: Recycling (List + Single DF) ---
#' # Merges 'df2' into EVERY element of 'test_rangex'
#' recycle_result <- join_all_lists(
#'   test_rangex, 
#'   df2, 
#'   join_on = "id"
#' )
#' }
#' 
#' @keywords internal
#' @noRd
#' 
join_df_or_lists <- function(..., join_on = NULL, remove_duplicate = "both") {
  
  # 1. Capture and Standardize Inputs
  input_args <- list(...)
  
  input_args <- input_args[!sapply(input_args, is.null)]
  
  standardized_lists <- lapply(input_args, function(arg) {
    if (is.data.frame(arg)) {
      # Wrap single DF in list so it is treated as a recyclable unit by Map
      return(list(arg))
    } else if (is.list(arg) && !is.data.frame(arg)) {
      return(arg)
    } else {
      stop("All arguments must be data.frames or lists of data.frames")
      # return(arg)
    }
  })
  
  # 2. Regex Definitions for duplicates
  pat_i <- "^i\\."
  pat_num <- "\\.\\d+$"
  
  # 3. Main Logic: Map + Reduce
  result_list <- do.call(Map, c(list(f = function(...) {
    
    tables <- list(...)
    
    final_dt <- Reduce(function(x, y) {
      data.table::setDT(x)
      data.table::setDT(y)
      
      # Determine merge columns
      merge_cols <- if(!is.null(join_on)) join_on else base::intersect(names(x), 
                                                                       names(y))
      
      # Perform Join (data.table update join syntax)
      dt_joined <- x[y, on = merge_cols]
      
      # --- Duplicate Removal Logic ---
      if (!is.null(remove_duplicate) && remove_duplicate != "none") {
        
        cols_to_remove <- c()
        
        # Case 1: Remove "i." prefixed columns
        if (remove_duplicate %in% c("string", "both")) {
          cols_to_remove <- c(cols_to_remove, grep(pat_i, names(dt_joined), 
                                                   value = TRUE))
        }
        
        # Case 2: Remove numeric suffixed columns (.1, .2)
        if (remove_duplicate %in% c("numeric", "both")) {
          cols_to_remove <- c(cols_to_remove, grep(pat_num, names(dt_joined), 
                                                   value = TRUE))
        }
        
        # Execute removal efficiently by reference
        if (length(cols_to_remove) > 0) {
          cols_to_remove <- unique(cols_to_remove)
          for (col in cols_to_remove) {
            data.table::set(dt_joined, j = col, value = NULL)
          }
        }
      }
      
      return(dt_joined)
    }, tables)
    
    return(final_dt)
  }), standardized_lists))
  
  if(length(result_list) == 1) result_list <- result_list[[1]]
  
  return(result_list)
}




##################################################################




#' generate_grid_list
#'
#' @param df A data frame
#' @param by A character vector
#' @param return lA data frame \code{'df'} or a list \code{'list'}
#' @param verbose A logical
#' 
#' @return A function
#' @keywords internal
#' @noRd
#'
generate_grid_list <- function(df, 
                               by, 
                               parameter = NULL,
                               return = 'df', 
                               verbose = FALSE) {
  

  if(missing(by)) {
    if(is.null(parameter)) {
      stop2c("Argument 'by' is required when constructing the grid_str")
    }
  } else if(is.null(by)) {
    if(is.null(parameter)) {
      stop2c("Argument 'by' is required when constructing the grid_str")
    }
  }
  
  
  if(is.null(df)) {
    parameter_null   <- FALSE
    parameter_as_var <- FALSE
  }
  
  if(!is.null(df)) {
    parameter_null   <- TRUE
    parameter_as_var <- FALSE
    if(!is.null(parameter)) {
      checkmate::assert_character(parameter)
      parameter_null <- FALSE
      if('parameter' %in% names(df)) {
        parameter_as_var <- TRUE
        by <- c(parameter, by)
      }
    }
    
    # df <- data.table::setDT(df)
    
    # 1. Keep only columns that exist in df - e.g., remove 'hypothesis'
    by <- intersect(by, names(df))
    if (length(by) == 0) {
      stop("None of the specified 'by' variables exist in the dataset.")
    }
    
    if(data.table::is.data.table(df)) {
      # expand.grid_args <- do.call(CJ, c(df[, ..by], unique = TRUE))
      # expand.grid_args <- lapply(df[, ..by], collapse::funique)
      expand.grid_args <- lapply(df[, mget(by)], collapse::funique)
    } else {
      expand.grid_args <- c(lapply(df[by], unique), 
                            list(stringsAsFactors = FALSE)  )
    }
  } # if(!is.null(df)) {
  
 
  
  

  if(!parameter_null) {
    if(!parameter_as_var) {
      expand.grid_args <- list()
      expand.grid_args[['parameter']] <- parameter
      # Move 'parameter' to the first position
      expand.grid_args <- expand.grid_args[c("parameter", 
                                             setdiff(names(expand.grid_args), 
                                                     "parameter"))]
    }
  }
  
  
  
  grid_data <- do.call(expand.grid, expand.grid_args)
  na_list <- lapply(grid_data, function(col) {
    lvl_names <- unique(col)
    setNames(as.list(rep(NA, length(lvl_names))), lvl_names)
  })
  null_list <- lapply(grid_data, function(col) {
    lvl_names <- unique(col)
    setNames(vector("list", length(lvl_names)), lvl_names)
  })
  if(return == 'df') return(grid_data)
  if(return == 'list') return(null_list)
  if(return == 'both') {
    return(list(grid_data = grid_data, null_list = null_list))
  }
}




#' comparison / hypothesis grid_str structure
#'
#' @param df A data frame
#' @param by A character vector
#' @param hypothesis see marginaleffects
#' @param verbose A logical
#' 
#' @return A function
#' @keywords internal
#' @noRd
#'
get_comparison_hypothesis_grid <- function(df, 
                                by = NULL, 
                                parameter = NULL,
                                # hypothesis = difference ~ pairwise,
                                hypothesis = NULL, 
                                verbose = FALSE) {
  
  estimate <- NULL;
  
  if(missing(by)) {
    if(is.null(parameter)) {
      stop2c("Argument 'by' is required when constructing the grid_str")
    }
  } else if(is.null(by)) {
    if(is.null(parameter)) {
      stop2c("Argument 'by' is required when constructing the grid_str")
    }
  }
  
  gdout <- generate_grid_list(df, 
                              parameter = parameter, 
                              by = by,
                              return = 'df',
                              verbose = verbose)
  
  if(is.null(hypothesis)) {
    out <- gdout
  } else if(!is.null(hypothesis)) {
    if(is.null(parameter)) h_out_by <- NULL else h_out_by <- "parameter"
    out <- data.table::setDT(gdout)
    marginaleffects_hypothesis_fun <- utils::getFromNamespace('get_hypothesis', 
                                                              'marginaleffects')
    out <- out[, estimate := 1][, 
                                marginaleffects_hypothesis_fun(.SD, 
                                                                 by = by, 
                                                                 hypothesis = 
                                                                   hypothesis), 
                                by = h_out_by][, 
                                               estimate := NA] %>% 
      data.table::setnames("estimate", "range") %>% 
      data.table::setDF()
  } # if(is.null(hypothesis)) { else if(!is.null(hypothesis)) {
  
  return(out)
} # get_comparison_str


# get_comparison_hypothesis_grid(fit$data, parameter = 'apgv')

# get_comparison_hypothesis_grid(fit$data, parameter= c('apgv', 'pgv'),by='sex')





#' custom_marginaleffects_equivalence for for \code{insight get_data}
#'
#' @param x A brms objects
#' @param effects levels of group levels - see [insight::get_data()]
#' @param component levels of group levels - see [insight::get_data()]
#' @param source levels of group levels - see [insight::get_data()]
#' @param verbose levels of group levels - see [insight::get_data()]
#' @param ... Other arguments
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
custom_marginaleffects_equivalence <- function (x, 
                                                equivalence  = NULL, 
                                                df = Inf, 
                                                draws  = NULL, 
                                                verbose = FALSE, 
                                                ...) {
  # ..filterby <- NULL;
  target_n <- NULL;
  ROPE_Percentage <- NULL;
  pd <- NULL;
  . <- NULL;
  
  
  dots <- list(...)
  equivalence_test_args       <- dots[['equivalence_test']]
  p_direction_args            <- dots[['p_direction']]
  
  if(is.null(equivalence_test_args) & is.null(p_direction_args)) {
    return(x)
  }
  
  

  equivalence_test_args_by    <-  equivalence_test_args[['by']]
  p_direction_args_by         <-  p_direction_args[['by']]
  
  equivalence_test_args_range <- equivalence_test_args[['range']]
  p_direction_args_null       <-  p_direction_args[['null']]
  
  format_eq <- equivalence_test_args[['format']]
  format_pd <- p_direction_args[['format']]
  format <- c(format_eq, format_pd)
  if(all(is.null(format))) {
    format <- TRUE
  } else if(any(format)) {
    format <- TRUE
  } else {
    format <- FALSE
  }
  
  digits_eq <- equivalence_test_args[['digits']]
  digits_pd <- p_direction_args[['digits']]
  digits_eqpd <- c(digits_eq, digits_pd)
  if(all(is.null(digits_eqpd))) {
    digits <- NULL
  } else if(any(!is.null(digits_eqpd))) {
    digits <- digits_eqpd[1]
  } else {
    # digits <- 2
  } 
  
  
  as_percent_eq <- equivalence_test_args[['as_percent']]
  as_percent_pd <- p_direction_args[['as_percent']]
  if(is.null(as_percent_eq)) as_percent_eq <- TRUE
  if(is.null(as_percent_pd)) as_percent_pd <- TRUE
  if(is.null(equivalence_test_args)) as_percent_eq <- FALSE
  if(is.null(p_direction_args)) as_percent_pd <- FALSE
  
  
  equivalence_test_args[['by']] <- NULL
  p_direction_args[['by']]      <- NULL
  equivalence_test_args[['by']]                 <- NULL
  p_direction_args[['by']]                      <- NULL
  equivalence_test_args[['format']]             <- NULL
  p_direction_args[['format']]                  <- NULL
  equivalence_test_args[['digits']]             <- NULL
  p_direction_args[['digits']]                  <- NULL
  equivalence_test_args[['as_percent']]         <- NULL
  p_direction_args[['as_percent']]              <- NULL
  
  
  
  
 
  
  if(!is.null(equivalence_test_args_range)) { # equivalence_test_args
    hypothesis  <- equivalence_test_args[['hypothesis']]
    hypothesis_group <- 
      get_hypothesis_group_fun(equivalence_test_args[['hypothesis']])
    df_filtered <- data.table::as.data.table(equivalence_test_args_range)
    if(is.null(x$hypothesis)) {
      if(!is.null(equivalence_test_args_by)) {
        filterby <- equivalence_test_args_by
        df_filtered[, (filterby) := lapply(.SD, as.factor), .SDcols = filterby]
        # my_levels <- lapply(df_filtered[, ..filterby], levels)
        my_levels <- lapply(df_filtered[, mget(filterby)], levels)
        df_filtered <- df_filtered[my_levels, on = filterby, nomatch = NULL]
        df_filtered <- df_filtered[, .SD[[ncol(.SD)]]]
      } 
    } else if(!is.null(x$hypothesis)) {
      if(!is.null(hypothesis_group)) {
        filterby <- hypothesis_group
        counts <- x[, .(target_n = .N), by = filterby]
        df_filtered <- df_filtered[counts, on = filterby,
                                   head(.SD, target_n), by = .EACHI]
      } else if(is.null(hypothesis_group)) {
        counts <- nrow(x)
        df_filtered <- head(df_filtered, counts)
      }
      
      df_filtered <- data.table::setDT(df_filtered)[, .SD[[ncol(.SD)]]]
      # For hypothesis, it's too complicated to assign a uniue range
      # work out is to just select the first range that same will be replicated 
      counts <- 1
      df_filtered <- head(df_filtered, counts)
    } # if(is.null(x$hypothesis)) { else if(!is.null(x$hypothesis)) {
    
    equivalence_test_args_range <- df_filtered
    equivalence_test_args$range <- equivalence_test_args_range
  } # if(!is.null(equivalence_test_args)) {
  
  
 
  eqp <- call_bayestest_eq(x = draws, 
                           equivalence_test = equivalence_test_args, 
                           p_direction = p_direction_args, 
                           via = 'marginals',
                           verbose = verbose) 
  
  if(!is_emptyx(eqp)) {
    eqp <- do.call(cbind, eqp)
    out <- cbind(x, eqp)
  } else {
    out <- x
  }
 
  if(as_percent_eq) {
    out <- out[, ROPE_Percentage := ROPE_Percentage * 100]
  }
  if(as_percent_pd) {
    out <- out[, pd := pd * 100]
  }
  if(!is.null(digits)) {
    out <- out[, lapply(.SD, function(z) 
      if (is.numeric(z)) round(z, digits) else z)]
  }
  if(format) {
    merge_ranges_eqpd_args <- list()
    merge_ranges_eqpd_args[['x']]       <- out
    merge_ranges_eqpd_args[['string']]  <- TRUE
    merge_ranges_eqpd_args[['sep']]     <- "; "
    merge_ranges_eqpd_args[['bracket']] <- 'square'
    merge_ranges_eqpd_args[['verbose']] <- FALSE
    out <- do.call(merge_ranges_eqpd, merge_ranges_eqpd_args)
  }
  return(out)
} # custom_marginaleffects_equivalence





#' call_bayestest_eq
#'
#' @param x A brms objects, data frame or a matrix
#' @param equivalence_test A list
#' @param via A character string
#' @param verbose Logical
#' @param ... Other arguments
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
call_bayestest_eq <- function(x,
                              equivalence_test = NULL, 
                              p_direction = NULL, 
                              via = 'marginals',
                              verbose = FALSE) {
  
  pd <- NULL;
  hby <- NULL;
  draw <- NULL;
  digits <- NULL;
  
  ##############################################################################
  # set argumenst
  ##############################################################################
  equivalence_test_args <- equivalence_test
  p_direction_args      <- p_direction
  
  equivalence_test_args_reformat_args_names   <- c('digits', 'reformat')
  equivalence_test_args_reformat_args         <- list()
  
  p_direction_args_reformat_args_names   <- c('digits', 'reformat', 'percent')
  p_direction_args_reformat_args         <- list()
  
  if(!is.null(equivalence_test_args)) {
    if(is.null(equivalence_test_args[['range']])) {
      equivalence_test_args[['range']] <- "default"
    }
    if(is.null(equivalence_test_args[['ci']])) {
      equivalence_test_args[['ci']] <- 0.95
    }
    if(is.null(equivalence_test_args[['verbose']])) {
      equivalence_test_args[['verbose']] <- FALSE
    }
    if(is.null(equivalence_test_args[['rvar_col']])) {
      equivalence_test_args[['rvar_col']] <- NULL
    }
    equivalence_test_args[['verbose']] <- verbose
  } # if(!is.null(equivalence_test_args)) {
  
  if(!is.null(p_direction_args)) {
    if(is.null(p_direction_args[['method']])) {
      p_direction_args[['method']] <- "direct"
    }
    if(is.null(p_direction_args[['null']])) {
      p_direction_args[['null']] <- 0
    }
    if(is.null(p_direction_args[['as_p']])) {
      p_direction_args[['as_p']] <- FALSE
    }
    if(is.null(p_direction_args[['remove_na']])) {
      p_direction_args[['remove_na']] <- TRUE
    }
    if(is.null(p_direction_args[['rvar_col']])) {
      p_direction_args[['rvar_col']] <- NULL
    }
  } # if(!is.null(p_direction_args)) {
  
  
  if(!is.null(equivalence_test_args)) {
    for (i in equivalence_test_args_reformat_args_names) {
      if(!is.null(equivalence_test_args[[i]])) {
        equivalence_test_args_reformat_args[[i]] <- equivalence_test_args[[i]] 
        equivalence_test_args[[i]]               <- NULL
      }
    } # for (i in equivalence_test_args_reformat_args_names) {
    if(is.null(equivalence_test_args_reformat_args[['digits']])) {
      equivalence_test_args_reformat_args[['digits']] <- 2
    }
    if(is.null(equivalence_test_args_reformat_args[['reformat']])) {
      equivalence_test_args_reformat_args[['reformat']] <- FALSE
    }
  } # if(!is.null(equivalence_test_args)) {
  
  
  
  
  if(!is.null(p_direction_args)) {
    for (i in p_direction_args_reformat_args_names) {
      if(!is.null(p_direction_args[[i]])) {
        p_direction_args_reformat_args[[i]] <- p_direction_args[[i]] 
        p_direction_args[[i]]               <- NULL
      }
    } # for (i in p_direction_args_reformat_args_names) {
    if(is.null(p_direction_args_reformat_args[['percent']])) {
      p_direction_args_reformat_args[['percent']] <- FALSE
    }
    if(is.null(p_direction_args_reformat_args[['digits']])) {
      p_direction_args_reformat_args[['digits']] <- 2
    }
    if(is.null(p_direction_args_reformat_args[['reformat']])) {
      p_direction_args_reformat_args[['reformat']] <- FALSE
    }
  } # if(!is.null(p_direction_args)) {
  
  
  
  
  
  
  
  ##############################################################################
  # set equivalence_test_args_run_marginals
  ##############################################################################
  equivalence_test_args_run_marginals <- function(x,
                                                  set_args,
                                                  reformat_args) {
    dimx1           <- dim(x)[1]
    verbose         <- set_args[['verbose']]
    set_args[['x']] <- x
    if(length(set_args$range) == dimx1) {
      set_args$range <- set_args$range
    } else if(length(set_args$range) == 1) {
      set_args$range <- rep(set_args$range, dimx1)
    } else {
      stop("range should be one or same as the length of comparisons")
    }
    
    exe_eq <- function(x) {
      current_args       <- set_args
      current_args$range <- set_args$range[[x]]
      current_args$x     <- set_args$x[x,] 
      do.call(bayestestR::equivalence_test, current_args)
    }
    
    # lapply(1:dimx1, exe_eq) %>% do.call(rbind, .)
    do.call(rbind, lapply(1:dimx1, exe_eq))
  } # equivalence_test_args_run_marginals
  
  
  
  ##############################################################################
  # set p_direction_args_args_run_marginals
  ##############################################################################
  p_direction_args_run_marginals <- function(x,
                                             set_args,
                                             reformat_args) {
    dimx1           <- dim(x)[1]
    verbose         <- set_args[['verbose']]
    set_args[['x']] <- x
    if(length(set_args$null) == dimx1) {
      set_args$null <- set_args$null
    } else if(length(set_args$null) == 1) {
      set_args$null <- rep(set_args$null, dimx1)
    } else {
      stop("null should be one or same as lenght of comparisons")
    }
    
    exe_eq <- function(x) {
      current_args       <- set_args
      current_args$null  <- set_args$null[[x]]
      current_args$x     <- set_args$x[x,] 
      do.call(bayestestR::p_direction, current_args)
    }
    pd <- NULL;
    out <- do.call(rbind, lapply(1:dimx1, exe_eq)) 
    # out %>% dplyr::mutate(pd = pd * 100) %>% dplyr::select(pd)
    out %>% dplyr::select(pd)
  } # p_direction_args_run_marginals
  
  
  
  
  ##############################################################################
  # via != 'marginals'
  ##############################################################################
  
  
  if(via != 'marginals') {
    allowed_eqtestby <- equivalence_test_args$by 
    if(is.null(equivalence_test_args$by)) {
      equivalence_test_args$by <- 'drawid'
    } else {
      n_eqtestby <- equivalence_test_args$by
      for (i in n_eqtestby) {
        if(! i %in% allowed_eqtestby) {
          stop("Allowed by options are: ", collapse_comma(allowed_eqtestby), 
               " , not ", i)
        }
      }
    }
    
    if(length(equivalence_test_args$by) == 1) {
      eqtestby <- equivalence_test_args$by
    } else {
      eqtestby <- 'eqby'
      x <- x %>% tidyr::unite(eqtestby, equivalence_test_args$by, 
                              sep = "_", remove = F)
    }
    
    if(is.null(equivalence_test_args$range)) {
      equivalence_test_args$range <- "default"
    } else if(!is.null(equivalence_test_args$range)) {
      if(is.character(equivalence_test_args$range)) {
        if(equivalence_test_args$range == "default") {
          equivalence_test_args$range <- "default"
        } else {
          stop("Range must be either character 'default', or a list")
        }
      } else {
        for (i in 1:length(equivalence_test_args$range)) {
          if(length(equivalence_test_args$range[[i]]) == 1) {
            equivalence_test_args$range[[i]] <- 
              c(-equivalence_test_args$range[[i]], 
                equivalence_test_args$range[[i]])
          } else if(length(equivalence_test_args$range[[i]]) == 2) {
            equivalence_test_args$range[[i]] <- equivalence_test_args$range[[i]]
          } else {
            stop2c("The length of range should be either 1, or 2, not ", 
                   length(equivalence_test_args$range[[i]]), ".",
                   " Check the following range: ",
                   deparse(equivalence_test_args$range[[i]]))
          }
        }
      } # if(equivalence_test_args$range == 'default') {
    } # if(is.null(equivalence_test_args$range)) {
    
    if(length(equivalence_test_args$range) == 
       length(unique(x[[eqtestby]]))) {
      equivalence_test_args$range <- equivalence_test_args$range
    } else {
      equivalence_test_args$range <- rep(equivalence_test_args$range, 
                                         length(unique(x[[eqtestby]])) )
    }
    
    if(length(equivalence_test_args$range) == 
       length(names(equivalence_test_args$range))) {
      equivalence_test_args$range <- equivalence_test_args$range
    } else {
      names(equivalence_test_args$range) <- unique(x[[eqtestby]])
    }
    
    
    hbyeqtestby <- c(hby, eqtestby)
    result <- setDT(x)[, {
      current_args <- equivalence_test_args
      # This below will work on unnamed equivalence_test_args$range
      # current_args$range <- equivalence_test_args$range[[.GRP]]
      # This below will work on named equivalence_test_args$range
      current_args$range <- equivalence_test_args$range[[ .BY[[1]] ]] 
      current_args$x <- draw
      do.call(bayestestR::equivalence_test, current_args)
    }, by = hbyeqtestby][
      # Chain a second call to round numeric columns
      , lapply(.SD, function(x) if(is.numeric(x)) round(x, digits) else x)
    ]
    
    result <- result %>% data.frame()
  } # if(via != 'marginals') {
  
  
  
  
  ##############################################################################
  # via == 'marginals'
  ##############################################################################
  
  if(via == 'marginals') {
    
    # set equivalence_test_args$range
    if(!is.null(equivalence_test_args)) {
      if(is.null(equivalence_test_args$range)) {
        equivalence_test_args$range <- "default"
      } else if(!is.null(equivalence_test_args$range)) {
        if(is.character(equivalence_test_args$range)) {
          if(equivalence_test_args$range == "default") {
            equivalence_test_args$range <- "default"
          } else {
            stop("Range must be either character 'default', or a list")
          }
        } else {
          for (i in 1:length(equivalence_test_args$range)) {
            if(length(equivalence_test_args$range[[i]]) == 1) {
              equivalence_test_args$range[[i]] <- 
                c(-equivalence_test_args$range[[i]], 
                  equivalence_test_args$range[[i]])
            } else if(length(equivalence_test_args$range[[i]]) == 2) {
              equivalence_test_args$range[[i]] <- 
                equivalence_test_args$range[[i]]
            } else {
              stop2c("The length of range should be either 1, or 2, not ", 
                     length(equivalence_test_args$range[[i]]), ".",
                     " Check the following range: ", 
                     deparse(equivalence_test_args$range[[i]]))
            }
          }
        } # if(equivalence_test_args$range == 'default') {
      } # if(is.null(equivalence_test_args$range)) {
    } # if(!is.null(equivalence_test_args)) {
    # End # set equivalence_test_args$range
    
  
    result <- list()
    if(!is.null(equivalence_test_args) & !is.null(p_direction_args)) {
      result[[1]] <- 
        equivalence_test_args_run_marginals(
          x = x,
          set_args = equivalence_test_args,
          reformat_args = equivalence_test_args_reformat_args)
      result[[2]] <- p_direction_args_run_marginals(
        x = x, 
        set_args = p_direction_args,
        reformat_args = p_direction_args_reformat_args)
    } else if(!is.null(equivalence_test_args)) {
      result[[1]] <- 
        equivalence_test_args_run_marginals(
          x = x,
          set_args = equivalence_test_args,
          reformat_args = equivalence_test_args_reformat_args)
    } else if(!is.null(p_direction_args)) {
      result[[1]] <- p_direction_args_run_marginals(
        x = x, 
        set_args = p_direction_args,
        reformat_args = p_direction_args_reformat_args)
    } else {
      result <- NULL
    }
  } # if(via == 'marginals') {
  
  
  return(result)
}





#' evaluate_comparison_fun
#'
#' @param data A brms objects, data frame or a matrix
#' @param hypothesis_args A list
#' @param return_draws A character string
#' @param nthreads A character string
#' @param probs A character string
#' @param ec_agg A character string
#' @param ei_agg A character string
#' @param na.rm A character string
#' @param verbose Logical
#' @param ... Other arguments
#'
#' @return A list
#' @keywords internal
#' @noRd
#'
evaluate_comparison_fun <- function(data,
                                    comparison_args,
                                    comparison_test = FALSE,
                                    conf_level = NULL,
                                    probs = c(0.025, 0.975),
                                    nthreads = 1,
                                    ec_agg = 'mean',
                                    ei_agg = "eti",
                                    na.rm = TRUE,
                                    digits = NULL,
                                    return_draws = FALSE,
                                    get_range_null_form = FALSE,
                                    set_grid = FALSE,
                                    verbose = FALSE) {
  
  
  . <- NULL
  estimate <- NULL
  ##########################################
  
  if(is.null(probs) & is.null(conf_level)) {
    stop2c("Please specify either 'probs' or 'conf_level'")
  } else if(!is.null(conf_level)) {
    conf <- ci <- conf_level
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  } else if(!is.null(probs)) {
    conf <- ci <- conf_level <- probs[2] - probs[1] 
    probs <- probs
  }
  
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', probtitles)
  
  ##########################################
  
  data <- data.table::setDT(data)
  `%chin%` <- data.table::`%chin%`
  if(get_range_null_form) {
    required_cols <- c("drawid", "draw")
    for (col in required_cols) {
      if (!col %in% names(data)) {
        if(col == "drawid") value <- as.factor(1)
        if(col == "draw") value <- 0
        set(data, j = col, value = value) 
      }
    }
  } # if(get_range_null_form) {
  
  
  if(!is.null(comparison_args[['by']])) {
    if(is.logical(comparison_args[['by']])) {
      comparison_args[['by']] <- NULL
    }
  }
  
  parameter_by             <- c('parameter', comparison_args[['by']])
  drawid_parameter         <- c('drawid', 'parameter')
  drawid_parameter_draw_by <- c(drawid_parameter, 'draw', 
                                comparison_args[['by']])
  
  if(!all(drawid_parameter_draw_by %chin% names(data))) {
    stop2c("The following variables are missing the data ",
           collapse_comma(setdiff(drawid_parameter_draw_by , names(data)))
    )
  }
  
  # ..drawid_parameter_draw_by ..setdrawidh_draw_estimate
  # data[, .SD[, drawid_parameter, with = FALSE]]  # if inside groups/by
  # # or
  # data[, mget(drawid_parameter)]  # returns list; use setDT() if needed
  
  if(get_range_null_form) {
    if(!set_grid) data <- unique(data[, mget(drawid_parameter_draw_by)])
  }

  comparison_draws <- data[, mget(drawid_parameter_draw_by)] %>% 
    data.table::setnames(., old = "draw", new = "estimate")
  
  if(get_range_null_form) {
    setdrawidh_draw_estimate <- setdiff(c(drawid_parameter_draw_by, "estimate"), 
                                        c('draw', 'drawid'))
    return(comparison_draws[, mget(setdrawidh_draw_estimate)][,
                                                          "estimate" := NA])
  }
  
  if(return_draws) {
    return(comparison_draws)
  }
  
  # New
  comparison_results <- comparison_draws[,
                          as.list(get_pe_ci_collapse(estimate, 
                                                     ec_agg= ec_agg,
                                                     ei_agg = ei_agg,
                                                     na.rm =na.rm,                 
                                                     nthreads = nthreads, 
                                                     probs = probs)),
                          by = parameter_by]
  
  if(is_emptyx(comparison_results)) return(NA)
  
  comparison_results <- comparison_results[,
                                           setnames(.SD, c("V1", "V2", "V3"),
                                                    c("estimate", "conf.low", 
                                                      "conf.high"))]
  
  # comparison_results <-
  #   comparison_draws[,
  #                   as.list(get_pe_ci_collapse(estimate,
  #                                              ec_agg= ec_agg,
  #                                              ei_agg = ei_agg,
  #                                              na.rm =na.rm,
  #                                              nthreads = nthreads,
  #                                              probs = probs)),
  #                   by = parameter_by
  #   ][,
  #     setnames(.SD, c("V1", "V2", "V3"),
  #              c("estimate", "conf.low", "conf.high"))
  #   ]
  
  if(!is.null(digits)) {
    comparison_results <- comparison_results[,
                                             lapply(.SD, 
                                                    function(z) if 
                                                    (is.numeric(z)) 
                                                      round(z, digits) 
                                                    else z)]
  }
  
  if(comparison_test) {
    attr(comparison_results, 'comparison_draws') <- comparison_draws
    attr(comparison_results, 'parameter_by') <- parameter_by
  }
  
  return(comparison_results)
}




#' evaluate_hypothesis_fun
#'
#' @param data A brms objects, data frame or a matrix
#' @param hypothesis_args A list
#' @param return_draws A character string
#' @param nthreads A character string
#' @param probs A character string
#' @param ec_agg A character string
#' @param ei_agg A character string
#' @param na.rm A character string
#' @param verbose Logical
#' @param ... Other arguments
#'
#' @return A list
#' @keywords internal
#' @noRd
#'
evaluate_hypothesis_fun <- function(data,
                                    hypothesis_args,
                                    hypothesis_test = FALSE,
                                    conf_level = NULL,
                                    probs = c(0.025, 0.975),
                                    nthreads = 1,
                                    ec_agg = 'mean',
                                    ei_agg = "eti",
                                    na.rm = TRUE,
                                    digits = NULL,
                                    return_draws = FALSE,
                                    get_range_null_form = FALSE,
                                    set_grid = FALSE,
                                    verbose = FALSE) {
  
  
  . <- NULL
  estimate <- NULL
  
  if(is.null(hypothesis_args[['hypothesis']])) {
    return(NULL)
  }
  
  
  
  ##########################################
  
  if(is.null(probs) & is.null(conf_level)) {
    stop2c("Please specify either 'probs' or 'conf_level'")
  } else if(!is.null(conf_level)) {
    conf <- ci <- conf_level
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  } else if(!is.null(probs)) {
    conf <- ci <- conf_level <- probs[2] - probs[1] 
    probs <- probs
  }
  
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', probtitles)
  
  ##########################################

  data <- data.table::setDT(data)
  `%chin%` <- data.table::`%chin%`
  if(get_range_null_form) {
    required_cols <- c("drawid", "draw")
    for (col in required_cols) {
      if (!col %in% names(data)) {
        if(col == "drawid") value <- as.factor(1)
        if(col == "draw") value <- 0
        set(data, j = col, value = value) 
      }
    }
  } # if(get_range_null_form) {
  
  
  hypothesis_group <- get_hypothesis_group_fun(hypothesis_args[['hypothesis']])

  hypothesis_by <- hypothesis_args[['by']]
  
  hypothesis_groupby       <- c('hypothesis', hypothesis_group)
  
  # When getting get_grid_form = TRUE, 'hypothesis' is not available 
  hypothesis_groupby   <- intersect(hypothesis_groupby, names(data))
  
  drawid_parameter         <- c('drawid', 'parameter')
  drawid_parameter_draw_by <- c(drawid_parameter, 'draw')
    
  drawid_parameter_draw_by <- c(drawid_parameter_draw_by, hypothesis_by)
  
  
  if(!all(drawid_parameter_draw_by %chin% names(data))) {
    stop2c("The following variables are missing the data ",
           collapse_comma(setdiff(drawid_parameter_draw_by , names(data)))
           )
  }
  
  if(get_range_null_form) {
    if(!set_grid) data <- unique(data[, mget(drawid_parameter_draw_by)])
  }
  
  if(!is.null(hypothesis_by)) {
    parameter_hypothesis <- c('parameter', hypothesis_group, 'hypothesis')
    marginaleffects_hypothesis_fun <- utils::getFromNamespace('get_hypothesis', 
                                                              'marginaleffects')
  } else if(is.null(hypothesis_by)) {
    parameter_hypothesis <- c('parameter', hypothesis_group)
    marginaleffects_hypothesis_fun <- utils::getFromNamespace('get_hypotheses', 
                                                              'marginaleffects')
  }

    hypothesis_data <- data[, mget(drawid_parameter_draw_by)] 
    
    hypothesis_data <- hypothesis_data %>%
      data.table::setnames(., old = "draw", new = "estimate")
  

  hypothesis_args[['equivalence_test']] <- NULL
  hypothesis_args[['p_direction']] <- NULL
  hypothesis_args[['range_null']] <- NULL
  
  hypothes_draws <- hypothesis_data[, do.call(marginaleffects_hypothesis_fun, 
                                              c(list(.SD), 
                                                hypothesis_args)), 
                                    by = drawid_parameter]
 
  if(is.null(hypothesis_by)) {
    hypothes_draws <- hypothes_draws[, 
                                     data.table::setnames(.SD, 
                                                          old = "V1", 
                                                          new = "estimate")]
  }
  
  if(get_range_null_form) {
    setdrawidh_draw_estimate <-setdiff(c(drawid_parameter_draw_by, "hypothesis",
                                          "estimate"), 
                                        c('draw', 'drawid', 
                                          hypothesis_by))
    setdrawidh_draw_estimate <- c(setdrawidh_draw_estimate,hypothesis_groupby)
    return(hypothes_draws[, mget(setdrawidh_draw_estimate)][,
                                                        "estimate" := NA])
  }
  
  if(return_draws) return(hypothes_draws)
  
  
  # New
  hypothesis_results <- 
    hypothes_draws[,
                   as.list(get_pe_ci_collapse(estimate, 
                                              ec_agg= ec_agg,
                                              ei_agg = ei_agg,
                                              na.rm =na.rm,                 
                                              nthreads = nthreads, 
                                              probs = probs)),
                   by = parameter_hypothesis]
  
  if(is_emptyx(hypothesis_results)) return(NA)
  
  hypothesis_results <- hypothesis_results[,
                                           setnames(.SD, c("V1", "V2", "V3"),
                                                    c("estimate", "conf.low", 
                                                      "conf.high"))]
  
  # hypothesis_results <- 
  #   hypothes_draws[,
  #                  as.list(get_pe_ci_collapse(estimate, 
  #                                             ec_agg= ec_agg,
  #                                             ei_agg = ei_agg,
  #                                             na.rm =na.rm,                 
  #                                             nthreads = nthreads, 
  #                                             probs = probs)),
  #                  by = parameter_hypothesis
  #   ][,
  #     setnames(.SD, c("V1", "V2", "V3"),
  #              c("estimate", "conf.low", "conf.high"))
  #   ]
  
  
  if(!is.null(digits)) {
    hypothesis_results <- hypothesis_results[,
                                             lapply(.SD, 
                                                    function(z) if 
                                                    (is.numeric(z)) 
                                                      round(z, digits) 
                                                    else z)]
  }
  
  if(hypothesis_test) {
    attr(hypothesis_results, 'hypothes_draws') <- hypothes_draws
    attr(hypothesis_results, 'parameter_hypothesis') <- parameter_hypothesis
  }
  
  return(hypothesis_results)
}





#' call_equivalence_test_p_direction
#'
#' @param data A brms objects, data frame or a matrix
#' @param hypothesis_args A list
#' @param return_draws A character string
#' @param nthreads A character string
#' @param probs A character string
#' @param ec_agg A character string
#' @param ei_agg A character string
#' @param na.rm A character string
#' @param verbose Logical
#' @param ... Other arguments
#'
#' @return A list
#' @keywords internal
#' @noRd
#' 
call_equivalence_test_p_direction <- function(data,
                                               comparison_args = NULL,
                                               hypothesis_args = NULL,
                                               evaluate_comparison = NULL,
                                               evaluate_hypothesis = NULL,
                                              rope_test = NULL,
                                              pd_test = NULL,
                                               by = NULL,
                                              conf_level = NULL,
                                              probs = c(0.025, 0.975),
                                              pd_as_percent = FALSE,
                                              rope_as_percent = FALSE,
                                              nthreads = 1,
                                               ec_agg = 'mean',
                                               ei_agg = "eti",
                                               na.rm = TRUE,
                                              return_draws = FALSE,
                                              get_range_null_form = FALSE,
                                               digits = 2,
                                               verbose = FALSE) {
  
  
  estimate <- NULL;
  ROPE_Percentage <- NULL;
  null <- NULL;
  . <- NULL;
  pd <- NULL;
  ##########################################
  
  if(is.null(probs) & is.null(conf_level)) {
    stop2c("Please specify either 'probs' or 'conf_level'")
  } else if(!is.null(conf_level)) {
    conf <- ci <- conf_level
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  } else if(!is.null(probs)) {
    conf <- ci <- conf_level <- probs[2] - probs[1] 
    probs <- probs
  }
  
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', probtitles)
  
  
  
  ################################################################
  

  # library(data.table)
  # keep_side = "x"  -> keep columns from x, drop RHS duplicates
  # keep_side = "i"  -> keep columns from i, drop LHS duplicates
  join_nodup <- function(x, i, by, keep_side = c("x", "i"), ...) {
    ..keep_cols <- NULL;
    keep_side <- match.arg(keep_side)
    out <- x[i, on = by, ...]  # allow passing nomatch, mult, etc.
    # .1 columns created by join
    dup1  <- grep("\\.1$", names(out), value = TRUE)
    if (length(dup1) == 0L) {
      # no duplicate columns -> nothing to fix
      return(out)
    }
    base  <- sub("\\.1$", "", dup1)
    if (keep_side == "x") {
      # keep x’s (unsuffixed), drop ".1"
      out[, (dup1) := NULL]
    } else {
      # keep i’s: drop base, keep ".1" and then strip suffix
      drop_base <- base
      keep_cols <- c(setdiff(names(out), drop_base), dup1)
      out <- out[, ..keep_cols]
      setnames(out, sub("\\.1$", "", names(out)))
    }
    out[]
  }
  
  reorder_num_last <- function(dt) {
    stopifnot(is.data.table(dt))
    is_num <- sapply(dt, is.numeric)
    nn  <- names(dt)[!is_num]   # non-numeric
    num <- names(dt)[ is_num]   # numeric
    setcolorder(dt, c(nn, num))
    invisible(dt)
  }
  
  ################################################################
  
  if(!data.table::is.data.table(data)) {
    data <- data.table::setDT(data)
  } 
  
 
  
  if(is.null(comparison_args[['by']])) {
    comparison_by <- by
  } else {
    comparison_by <- comparison_args[['by']]
  }
  
  
  if(is.null(hypothesis_args[['by']])) {
    hypothesis_by <- by
  } else {
    hypothesis_by <- hypothesis_args[['by']]
  }
  
  
  
  if(is.null(evaluate_comparison)) {
    if(is_emptyx(comparison_args)) {
      evaluate_comparison <- FALSE
    } else {
      evaluate_comparison <- TRUE
    }
  }
  
  
  if(is.null(evaluate_hypothesis)) {
    if(is_emptyx(hypothesis_args)) {
      evaluate_hypothesis <- FALSE
    } else {
      evaluate_hypothesis <- TRUE
      if(is.null(hypothesis_args[['hypothesis']])) {
        evaluate_hypothesis <- FALSE
      }
    }
  }
  
  
  
  evaluate_comparison_equivalence_test <- FALSE
  evaluate_comparison_p_direction <- FALSE
  if(evaluate_comparison) {
    comparison_equivalence_test_arg <- comparison_args[['equivalence_test']]
    comparison_p_direction_arg      <- comparison_args[['p_direction']]
    comparison_range_null           <- comparison_args[['range_null']]
    if(is_emptyx(comparison_equivalence_test_arg)) {
      evaluate_comparison_equivalence_test <- FALSE
    } else {
      evaluate_comparison_equivalence_test <- TRUE
    }
    if(is_emptyx(comparison_p_direction_arg)) {
      evaluate_comparison_p_direction <- FALSE
    } else {
      evaluate_comparison_p_direction <- TRUE
    }
  } # if(evaluate_comparison) {
  
  
  evaluate_hypothesis_equivalence_test <- FALSE
  evaluate_hypothesis_p_direction <- FALSE
  if(evaluate_hypothesis) {
    hypothesis_equivalence_test_arg <- hypothesis_args[['equivalence_test']]
    hypothesis_p_direction_arg      <- hypothesis_args[['p_direction']]
    hypothesis_range_null           <- hypothesis_args[['range_null']]
    if(is_emptyx(hypothesis_equivalence_test_arg)) {
      evaluate_hypothesis_equivalence_test <- FALSE
    } else {
      evaluate_hypothesis_equivalence_test <- TRUE
    }
    if(is_emptyx(hypothesis_p_direction_arg)) {
      evaluate_hypothesis_p_direction <- FALSE
    } else {
      evaluate_hypothesis_p_direction <- TRUE
    }
  } # if(evaluate_hypothesis) {
  
  
  
  
  
  
  if(!evaluate_comparison & !evaluate_hypothesis) {
    return(NULL)
  }
  
  
  if(!evaluate_comparison) {
    comparison_results <- NULL
    evaluate_comparison_equivalence_test <- FALSE
    evaluate_comparison_p_direction <- FALSE
  }
  
  
  if(!evaluate_hypothesis) {
    hypothesis_results <- NULL
    evaluate_hypothes_equivalence_test <- FALSE
    evaluate_hypothes_p_direction <- FALSE
  }
  
  
  if(!evaluate_comparison_equivalence_test & 
     !evaluate_comparison_p_direction) {
    comparison_eqpd_results <- NULL
  }
  
  if(!evaluate_hypothesis_equivalence_test & 
     !evaluate_hypothesis_p_direction) {
    hypothesis_eqpd_results <- NULL
  }
  
  
  # override above setting when rope_test and pd_test
  # evaluate_hypothesis_equivalence_test & evaluate_hypothesis_p_direction
  if(!is.null(rope_test)) {
    evaluate_hypothesis_equivalence_test <- NullFALSE(rope_test)
    evaluate_comparison_equivalence_test <- NullFALSE(rope_test)
  }
  if(!is.null(pd_test)) {
    evaluate_hypothesis_p_direction <- pd_test
    evaluate_comparison_p_direction <- pd_test
  }
  
  
  evaluate_comparison_fun_arg <- list()
  evaluate_comparison_fun_arg[['data']] <- data
  evaluate_comparison_fun_arg[['conf_level']] <- conf_level
  evaluate_comparison_fun_arg[['probs']] <- probs
  evaluate_comparison_fun_arg[['nthreads']] <- nthreads
  evaluate_comparison_fun_arg[['ec_agg']] <- ec_agg
  evaluate_comparison_fun_arg[['ei_agg']] <- ei_agg
  evaluate_comparison_fun_arg[['digits']] <- digits
  evaluate_comparison_fun_arg[['na.rm']] <- na.rm
  evaluate_comparison_fun_arg[['return_draws']] <- return_draws
  evaluate_comparison_fun_arg[['get_range_null_form']] <- get_range_null_form
  evaluate_comparison_fun_arg[['verbose']] <- verbose
  
  evaluate_hypothesis_fun_arg <- evaluate_comparison_fun_arg
  
  evaluate_comparison_fun_arg[['comparison_args']] <- comparison_args
  evaluate_hypothesis_fun_arg[['hypothesis_args']] <- hypothesis_args
  
  run_comparison_test <- FALSE
  run_hypothesis_test <- FALSE
  
  if(evaluate_comparison) {
    if(evaluate_comparison_equivalence_test | evaluate_comparison_p_direction) {
      run_comparison_test <- TRUE
    }
  }
  
  if(evaluate_hypothesis) {
    if(evaluate_hypothesis_equivalence_test | evaluate_hypothesis_p_direction) {
      run_hypothesis_test <- TRUE
    }
  }
  
  
  evaluate_comparison_fun_arg[['comparison_test']] <- run_comparison_test
  evaluate_hypothesis_fun_arg[['hypothesis_test']] <- run_hypothesis_test
  
  
  #########################################################################
  # Evaluate comparison
  #########################################################################
  if(evaluate_comparison) {
    comparison_results <- do.call(evaluate_comparison_fun,
                               evaluate_comparison_fun_arg)
    comparison_draws <- attr(comparison_results, 'comparison_draws')
    parameter_by     <- attr(comparison_results, 'parameter_by')
  } # if(evaluate_comparison) {

  
  #########################################################################
  # Evaluate hypothesis
  #########################################################################
  if(evaluate_hypothesis) {
    evaluate_hypothesis_fun_arg[['comparison_args']] <- NULL
    hypothesis_results <- do.call(evaluate_hypothesis_fun,
                               evaluate_hypothesis_fun_arg)
    hypothes_draws       <- attr(hypothesis_results, 'hypothes_draws')
    parameter_hypothesis <- attr(hypothesis_results, 'parameter_hypothesis')
  } # if(evaluate_hypothesis) {
  
  if(is.null(hypothesis_results)) {
    run_hypothesis_test <- FALSE
  }
  
 
  
  comparison_eqpd_results <- NULL
  hypothesis_eqpd_results <- NULL
 
  #########################################################################
  # Evaluate equivalence_test and p_direction for comparison
  #########################################################################
  if(run_comparison_test) {
    eqpdby <- parameter_by # c( "parameter", "sex")
    if(is.null(comparison_range_null)) {
      set_data_dt <- comparison_draws
    } else if(!is.null(comparison_range_null)) {
      set_data_dt <- join_nodup(comparison_draws, comparison_range_null, 
                                by = eqpdby, keep_side = "x")
      set_data_dt <- reorder_num_last(set_data_dt)
    } # if(is.null(comparison_range_null)) { else if
    
    set_data_dt <- stats::na.omit(set_data_dt, cols = 'estimate')
    
    comparison_eqpd_results <- set_data_dt[
      ,
      {
        if(evaluate_comparison_equivalence_test) {
          current_equivalence_test_args <- comparison_equivalence_test_arg
          grp_range <- range[1L]
          if (is.character(grp_range)) {
            grp_range <- as.numeric(strsplit(grp_range, ",")[[1L]])
          } else {
            grp_range <- grp_range
          }
          
          # Check and set up the grp_range
          if(any(is.na(grp_range))) {
            grp_range <- "default"
          } else if(is.list(grp_range)) {
            # grp_range <- "default"
          } else if(is.vector(grp_range)) {
            if(length(grp_range) == 1) {
              grp_range <- c(-grp_range, grp_range)
            } else if(length(grp_range) > 2) {
              stop2c("The range should be 'default' or a vector of 2 
                     numeric values (e.g., c(-0.1,0.1)).")
            }
          }
          
          
          if(is.list(grp_range)) grp_range <- unlist(grp_range)
          current_equivalence_test_args$range <- grp_range
          current_equivalence_test_args$x <- estimate
          eqres <- do.call(bayestestR::equivalence_test, 
                           current_equivalence_test_args)
          eq_dt <- as.data.table(eqres)
          if(rope_as_percent) {
            eq_dt <- eq_dt[, ROPE_Percentage := ROPE_Percentage * 100]
          }
        }
        
        if(evaluate_comparison_p_direction) {
          current_p_direction_args <- comparison_p_direction_arg
          grp_null <- null[1L]
          if (is.character(grp_null)) {
            grp_null <- as.numeric(strsplit(grp_null, ",")[[1L]])
          } else {
            grp_null <- grp_null
          }
          if(is.list(grp_null)) grp_null <- unlist(grp_null)
          
          
          # Check and set up the grp_null
          if(any(is.na(grp_null))) {
            grp_null <- 0
          } else if(is.list(grp_null)) {
            # grp_null <- "default"
          } else if(is.vector(grp_null)) {
            if(length(grp_null) == 1) {
              grp_null <- grp_null
            } else if(length(grp_null) > 1) {
              stop2c("The null should be a vector of a 
                     numeric values (e.g., c(0)).")
            }
          }
          
          
          current_p_direction_args$null <- grp_null
          current_p_direction_args$x <- estimate
          pdres <-do.call(bayestestR::p_direction, current_p_direction_args)
          pd_dt <- as.data.table(pdres)
          if(pd_as_percent) pd_dt <- pd_dt[, .(pd = pd * 100)]
          pd_dt <- pd_dt[, 'pd_null' := grp_null]
        }
        
        if(evaluate_comparison_equivalence_test & 
           evaluate_comparison_p_direction) {
          out <- cbind(eq_dt, pd_dt)
        } else if(evaluate_comparison_equivalence_test & 
                  !evaluate_comparison_p_direction) {
          out <- eq_dt
        } else if(!evaluate_comparison_equivalence_test & 
                  evaluate_comparison_p_direction) {
          out <- pd_dt
        } else if(!evaluate_comparison_equivalence_test & 
                  !evaluate_comparison_p_direction) {
          out <- NULL
        }
        out
      },
      by = eqpdby
    ][
      ,
      lapply(.SD, function(z) if (is.numeric(z)) round(z, digits) else z)
    ][]

  } # ifrun_comparison_test) {
  
  
  #########################################################################
  # Evaluate equivalence_test and p_direction for hypothesis
  #########################################################################
  if(run_hypothesis_test) {
    eqpdby <- parameter_hypothesis
    if(is.null(hypothesis_range_null)) {
      set_data_dt <- hypothes_draws
    } else if(!is.null(hypothesis_range_null)) {
      set_data_dt <- join_nodup(hypothes_draws, hypothesis_range_null, 
                                by = eqpdby, keep_side = "x")
      set_data_dt <- reorder_num_last(set_data_dt)
    } # if(is.null(hypothesis_range_null)) { else if
    
    set_data_dt <- stats::na.omit(set_data_dt, cols = 'estimate')
    
    hypothesis_eqpd_results <- set_data_dt[
      ,
      {
        if(evaluate_hypothesis_equivalence_test) {
          current_equivalence_test_args <- hypothesis_equivalence_test_arg
          grp_range <- range[1L]
          if (is.character(grp_range)) {
            grp_range <- as.numeric(strsplit(grp_range, ",")[[1L]])
          } else {
            grp_range <- grp_range
          }
          if(is.list(grp_range)) grp_range <- unlist(grp_range)
          
          # Check and set up the grp_range
          if(any(is.na(grp_range))) {
            grp_range <- "default"
          } else if(is.list(grp_range)) {
            # grp_range <- "default"
          } else if(is.vector(grp_range)) {
            if(length(grp_range) == 1) {
              grp_range <- c(-grp_range, grp_range)
            } else if(length(grp_range) > 2) {
              stop2c("The range should be 'default' or a vector of 2 
                     numeric values (e.g., c(-0.1,0.1)).")
            }
          }
          
          
          current_equivalence_test_args$range <- grp_range
          current_equivalence_test_args$x <- estimate
          eqres <- do.call(bayestestR::equivalence_test, 
                           current_equivalence_test_args)
          eq_dt <- as.data.table(eqres)
          if(rope_as_percent) {
            eq_dt <- eq_dt[, ROPE_Percentage := ROPE_Percentage * 100]
          }
        }
        
        if(evaluate_hypothesis_p_direction) {
          current_p_direction_args <- hypothesis_p_direction_arg
          grp_null <- null[1L]
          if (is.character(grp_null)) {
            grp_null <- as.numeric(strsplit(grp_null, ",")[[1L]])
          } else {
            grp_null <- grp_null
          }
          if(is.list(grp_null)) grp_null <- unlist(grp_null)
          
          
          # Check and set up the grp_null
          if(any(is.na(grp_null))) {
            grp_null <- 0
          } else if(is.list(grp_null)) {
            # grp_null <- "default"
          } else if(is.vector(grp_null)) {
            if(length(grp_null) == 1) {
              grp_null <- grp_null
            } else if(length(grp_null) > 1) {
              stop2c("The null should be a vector of a 
                     numeric values (e.g., c(0)).")
            }
          }
          
          
          current_p_direction_args$null <- grp_null
          current_p_direction_args$x <- estimate
          pdres <-do.call(bayestestR::p_direction, current_p_direction_args)
          pd_dt <- as.data.table(pdres)
          if(pd_as_percent) pd_dt <- pd_dt[, .(pd = pd * 100)]
          pd_dt <- pd_dt[, 'pd_null' := grp_null]
        }
        
        if(evaluate_hypothesis_equivalence_test & 
           evaluate_hypothesis_p_direction) {
          out <- cbind(eq_dt, pd_dt)
        } else if(evaluate_hypothesis_equivalence_test & 
                  !evaluate_hypothesis_p_direction) {
          out <- eq_dt
        } else if(!evaluate_hypothesis_equivalence_test & 
                  evaluate_hypothesis_p_direction) {
          out <- pd_dt
        } else if(!evaluate_hypothesis_equivalence_test & 
                  !evaluate_hypothesis_p_direction) {
          out <- NULL
        }
        out
      },
      by = eqpdby
    ][
      ,
      lapply(.SD, function(z) if (is.numeric(z)) round(z, digits) else z)
    ][]

  } # if(run_hypothesis_test) {
  
  
  comparison_hypothesis_results <- list()
  
  comparison_results <- join_df_or_lists(comparison_results, 
                                         comparison_eqpd_results, 
                                         join_on = NULL, 
                                         remove_duplicate = "both") 
  
  comparison_hypothesis_results[['comparison']] <- comparison_results 
  
  
  hypothesis_results <- join_df_or_lists(hypothesis_results, 
                                         hypothesis_eqpd_results, 
                                                 join_on = NULL, 
                                                 remove_duplicate = "both")  
  
  comparison_hypothesis_results[['hypothesis']] <- hypothesis_results 
  
  return(comparison_hypothesis_results)
}





#' merge_ranges_eqpd
#'
#' @param X A 
#' @param string A list
#' @param sep A character string
#' @param bracket A character string
#' @param verbose Logical
#'
#' @return A list
#' @keywords internal
#' @noRd
#' 
merge_ranges_eqpd <- function(x, 
                              string = FALSE, 
                              sep = " - ", 
                              bracket = 'square',
                              verbose = FALSE) {
  is_list <- is.list(x) && !is.data.frame(x) && !"data.table" %in% class(x)
  if(bracket == 'square') setbracket <- c("[", "]")
  if(bracket == 'round') setbracket <- c("(", ")")
  inner_merge <- function(dt, string, sep) {
    . <- NULL;
    ROPE_range <- NULL;
    ROPE_low <- NULL;
    ROPE_high <- NULL;
    HDI_range <- NULL;
    HDI_low <- NULL;
    HDI_high <- NULL;
    if (!"data.table" %in% class(dt)) data.table::setDT(dt)
    rope_pos <- match("ROPE_low", names(dt))
    hdi_pos <- match("HDI_low", names(dt))
    if (!is.na(rope_pos)) {
      if (string) {
        dt[, ROPE_range := paste0(setbracket[1], ROPE_low, sep, ROPE_high,
                                  setbracket[2])]
      } else {
        # Correct: cbind creates matrix → list-column of 2-element vectors
        dt[, ROPE_range := list(.(ROPE_low, ROPE_high))]
      }
      dt[, c("ROPE_low", "ROPE_high") := NULL]
      # Build order: left+ROPE_range+right (exclude ROPE_range from right)
      nms <- names(dt)
      left <- nms[1:(rope_pos-1)]
      right <- setdiff(nms[rope_pos:length(nms)], "ROPE_range")
      data.table::setcolorder(dt, c(left, "ROPE_range", right))
    }
    if (!is.na(hdi_pos)) {
      # Recompute hdi_pos after ROPE (if processed)
      hdi_pos <- match("HDI_low", names(dt))
      if (is.na(hdi_pos)) return(dt)  # Already processed
      if (string) {
        dt[, HDI_range := paste0(setbracket[1], HDI_low, sep, HDI_high,
                                 setbracket[2])]
      } else {
        dt[, HDI_range := list(.(HDI_low, HDI_high))]
      }
      dt[, c("HDI_low", "HDI_high") := NULL]
      # Build order: left+HDI_range+right (exclude HDI_range from right)
      nms <- names(dt)
      left <- nms[1:(hdi_pos-1)]
      right <- setdiff(nms[hdi_pos:length(nms)], "HDI_range")
      data.table::setcolorder(dt, c(left, "HDI_range", right))
    }
    return(dt)
  }
  if (is_list) {
    out <- lapply(x, inner_merge, string, sep)
  } else {
    out <- inner_merge(x, string = string, sep = sep)
  }
  return(out)
}





#' set_up_equivalence_test_p_direction_args
#'
#' @param inbound_arguments A list
#' @param checking_inline A logical, TRUE for marginal_* 
#' and FALSE for hypothesis test
#' @param verbose levels of group levels - see [insight::get_data()]
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
set_up_equivalence_test_p_direction_args <- function(inbound_arguments,
                                                     checking_inline,
                                                     xcall,
                                                     verbose = FALSE) {
  
  # When set_up_equivalence_test_p_direction_args called from 
  # marginal_* that too is from hypothesis_test()
  if(grepl("^hypothesis_test", xcall) |
     grepl("^hypothesis_test", xcall)) {
    if(!is.list(inbound_arguments$equivalence_test)) {
      inbound_arguments$equivalence_test <- NULL
    }
    if(!is.list(inbound_arguments$equivalence_test)) {
      inbound_arguments$p_direction <- NULL
    }
    format_eq <- format_pd <- NULL
    get_eq_form <- get_pd_form <- NULL
    get_eq_value <- get_pd_value <- NULL
    rope_test <- pd_test <- FALSE
  }
 
  
  
  
  if(!is.null(inbound_arguments$equivalence_test)) {
    format_eq <- inbound_arguments$equivalence_test[['format']]
    get_eq_form  <- inbound_arguments$equivalence_test[['get_form']]
    get_eq_value <- inbound_arguments$equivalence_test[['get_value']]
    rope_test <- TRUE
  } else {
    format_eq <- NULL
    get_eq_form <- NULL
    get_eq_value <- NULL
    rope_test <- FALSE
  }
  if(!is.null(inbound_arguments$p_direction)) {
    format_pd <- inbound_arguments$p_direction[['format']]
    get_pd_form <- inbound_arguments$p_direction[['get_form']]
    get_pd_value <- inbound_arguments$p_direction[['get_value']]
    pd_test  <- TRUE
  } else {
    format_pd <- NULL
    get_pd_form <- NULL
    get_pd_value <- NULL
    pd_test  <- FALSE
  }
  
  format <- c(format_eq, format_pd)
  if(all(is.null(format))) {
    format <- TRUE
  } else if(any(format)) {
    format <- TRUE
  } else {
    format <- FALSE
  }
  
  get_range_null_form  <- c(get_eq_form, get_pd_form)
  if(all(is.null(get_range_null_form))) {
    get_range_null_form  <- FALSE
  } else if(any(get_range_null_form)) {
    get_range_null_form  <- TRUE
  } else {
    get_range_null_form  <- FALSE
  }
  
  get_range_null_value  <- c(get_eq_value, get_pd_value)
  if(all(is.null(get_range_null_value))) {
    get_range_null_value  <- FALSE
  } else if(any(get_range_null_value)) {
    get_range_null_value  <- TRUE
  } else {
    get_range_null_value  <- FALSE
  }
  
  check_equivalence_test_full.args <- FALSE
  if(!is.null(inbound_arguments$equivalence_test)) {
    if(is.null(inbound_arguments$equivalence_test$inline)) {
      inbound_arguments$equivalence_test$inline <- FALSE
      if(checking_inline) inbound_arguments$equivalence_test <- NULL
      check_equivalence_test_full.args <- TRUE
    } else if(!is.null(inbound_arguments$equivalence_test$inline)) {
      if(!inbound_arguments$equivalence_test$inline) {
        if(checking_inline) inbound_arguments$equivalence_test <- NULL
        check_equivalence_test_full.args <- TRUE
      } else if(inbound_arguments$equivalence_test$inline) {
        check_equivalence_test_full.args <- FALSE
      }
    }
  } # if(!is.null(inbound_arguments$equivalence_test)) {
  
  check_p_direction_full.args <- FALSE
  if(!is.null(inbound_arguments$p_direction)) {
    if(is.null(inbound_arguments$p_direction$inline)) {
      inbound_arguments$p_direction$inline <- FALSE
      if(checking_inline) inbound_arguments$p_direction <- NULL
      check_p_direction_full.args <- TRUE
    } else if(!is.null(inbound_arguments$p_direction$inline)) {
      if(!inbound_arguments$p_direction$inline) {
        if(checking_inline) inbound_arguments$p_direction <- NULL
        check_p_direction_full.args <- TRUE
      } else if(inbound_arguments$p_direction$inline) {
        check_p_direction_full.args <- FALSE
      }
    }
  } # if(!is.null(inbound_arguments$p_direction)) {
  
  
  if(!checking_inline) {
    check_equivalence_test_full.args <- TRUE
    check_p_direction_full.args <- TRUE
  }
  
  
  if(!check_equivalence_test_full.args & !check_p_direction_full.args) {
    unlock_replace_bind(package = "marginaleffects", what = "equivalence",
                        replacement = custom_marginaleffects_equivalence, 
                        ept_str = T)
  }
  
  if(check_equivalence_test_full.args | check_p_direction_full.args) {
    if(checking_inline) inbound_arguments[['hypothesis']] <- NULL
  }
  
  
  
  out <- list(format = format, get_range_null_form = get_range_null_form,
              get_range_null_value = get_range_null_value,
              check_equivalence_test_full.args = 
                check_equivalence_test_full.args,
              check_p_direction_full.args = check_p_direction_full.args,
              inbound_arguments = inbound_arguments,
              rope_test = rope_test,
              pd_test = pd_test)
  
  
  return(out)
}


#' suppresswar_equivalence_test_p_direction_args
#'
#' @param inbound_arguments A list
#' @param parm A string
#' @param verbose levels of group levels - see [insight::get_data()]
#' 
#' @return A list
#' @keywords internal
#' @noRd
#'
suppresswar_equivalence_test_p_direction_args <- function(inbound_arguments,
                                                          parm,
                                                          verbose = FALSE) {
  
  parameter <- NULL;
  
  suppresswar <- FALSE
  
  temprangecheck <- inbound_arguments$equivalence_test$range 
  
  if(!is.null(temprangecheck)) {
    if(!is.data.frame(temprangecheck) | 
       !data.table::is.data.table(temprangecheck)) {
      if(is.character(temprangecheck)) {
        if(temprangecheck == "default") {
          temprangecheck <- NULL
        }
      }
    }
  } # if(!is.null(temprangecheck)) {
  
  inbound_arguments$equivalence_test$range <- temprangecheck
  
  
  if(!is.null(inbound_arguments$equivalence_test)) {
    suppresswar <- TRUE
    inbound_arguments$equivalence_test$by <- inbound_arguments[['by']]
    if(!is.null(inbound_arguments$equivalence_test$range)) {
      inbound_arguments$equivalence_test$range <-
        inbound_arguments$equivalence_test$range %>% 
        dplyr::filter(parameter == parm) 
    } # if(!is.null(inbound_arguments$equivalence_test$range)) {
  }
  
  if(!is.null(inbound_arguments$p_direction)) {
    suppresswar <- TRUE
    inbound_arguments$p_direction$by <- inbound_arguments[['by']]
  } 
  
  
  out <- list(suppresswar = suppresswar, inbound_arguments = inbound_arguments)
  return(out)
}



#' DT_to_data_frames
#'
#' @param x A list of list of data.tables or a data.table
#' @return A list of of data.frames or a data.frame
#' @keywords internal
#' @noRd
#'
DT_to_data_frames <- function(x) {
  if (is.list(x) && all(sapply(x, inherits, "data.table"))) {
    return(lapply(x, as.data.frame))
  } else if (inherits(x, "data.table")) {
    return(as.data.frame(x))
  } else if (inherits(x, "data.frame")) { # new layes
    return(x)
  }
  stop("Input must be a data.table or list of data.tables")
}


#' marginalstyle_reformat
#'
#' @param out A data.frame or a data.table
#' @param set_names_ Additional arguments
#' @return A data.frame 
#' @keywords internal
#' @noRd
#'
marginalstyle_reformat <- function(out, set_names_) {
  . <- NULL;
  out <- out %>% 
    dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
    dplyr::rename(!!as.symbol(set_names_[1]) := 
                    dplyr::all_of('estimate')) %>% 
    dplyr::rename(!!as.symbol(set_names_[2]) := 
                    dplyr::all_of('conf.low')) %>% 
    dplyr::rename(!!as.symbol(set_names_[3]) := 
                    dplyr::all_of('conf.high')) %>% 
    dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
    data.frame()
  return(out)
}




#' Populate NA elements in list with template
#' \code{data.frame}/\code{data.table}
#'
#' Finds the first non-NA \code{data.frame} or \code{data.table} in the list as
#' a template. Populates scalar \code{NA} elements by copying the template and
#' setting specified columns to \code{NA}. Preserves input class:
#' \code{data.frame} inputs yield \code{data.frame} outputs; \code{data.table}
#' inputs yield \code{data.table} outputs.
#'
#' @param lst A named list containing \code{data.frame}s, \code{data.table}s,
#'   and/or scalar \code{NA}s.
#' @param col_names Character vector of column names to set to \code{NA}.
#'   Defaults to \code{"draw"}.
#'
#' @return Modified list matching input class.

#' @examples
#' lst <- list(pgv = data.frame(a = 1:2), atgv = NA)
#' populate_na_elements(lst, "a")
#' 
#' @keywords internal
#' @noRd
#' 
populate_na_elements <- function(lst, col_names = "draw") {
  stopifnot(is.list(lst))
  
  for (i in 1:length(lst)) {
    if(is_emptyx(lst[[i]])) lst[[i]] <- NA
  }
  
  template <- NULL
  class_template <- NULL
  
  # Find first data.frame or data.table
  for (el in lst) {
    if (is.data.frame(el) || methods::is(el, "data.table")) {
      template <- el
      class_template <- class(el)
      break
    }
  }
  if (is.null(template)) {
    # stop("No data.frame or data.table found in list")
    return(template)
  }
  
  # Populate scalar NA elements only
  for (i in seq_along(lst)) {
    el <- lst[[i]]
    if (length(el) == 1 && is.na(el)) {
      new_obj <- template
      for (col in col_names) {
        if (col %in% names(new_obj)) {
          new_obj[[col]] <- NA
        }
      }
      class(new_obj) <- class_template
      lst[[i]] <- new_obj
    }
  }
  return(lst)
}



#' Replace NA with integer
#' 
#' @description
#' Useful in \code{marginaleffects method = 'pkg'} \code{get_growthparameters}
#' as NA are not allowedNot used anywhere till \code{10/03/2026}
#' 
#'
#' @param x numeric value
#' @param val integer
#' @param verbose logical
#' @return A data.frame 
#' @keywords internal
#' @noRd
#'
NAtoint <- function(x, val = 999999L, verbose = FALSE) {
  if(is.null(val)) val <- 999999L
  if(!is.integer(val)) stop2c("'val' must be an integer")
  if(is.na(x)) {
    if(verbose) {
      msg <- sprintf("All '%s' values replaced with integer '%s'", NA, val)
      message2c(msg)
    }
    return(as.integer(val))
  } else if (x == val) {
    if(verbose) {
      msg <- sprintf("All '%s' values replaced with '%s'", val, NA)
      message2c(msg)
    }
    return(NA)
  } else {
    return(x)
  }
} # NAtoint



#' replace_int_with_na Replace integer with NA
#' 
#' @description
#' Useful in \code{marginaleffects method = 'pkg'} \code{get_growthparameters}
#' as NA are not allowedNot used anywhere till \code{10/03/2026}
#'
#' @param x numeric value
#' @param val integer
#' @param int_type integer
#' @param verbose logical
#' @return A data.frame 
#' @keywords internal
#' @noRd
#'
replace_int_with_na <- function(x, val = 999999L, int_type = NA_integer_,
                                verbose = FALSE) {
  if(is.null(val)) val <- 999999L
  if(!is.integer(val)) stop2c("'val' must be an integer")
  
  if(is.null(int_type)) int_type <- NA_integer_
  
  if (inherits(x, "data.table")) {
    for (j in seq_along(x)) {
      col <- x[[j]]
      if (is.integer(col)) {
        idx <- which(col == val)
        if (length(idx)) {
          data.table::set(x, i = idx, j = j, value = int_type)
        }
      }
    }
    return(x)
  }
  
  if (inherits(x, "data.frame")) {
    x[] <- lapply(x, function(col) {
      if (is.integer(col)) {
        col[col == val] <- int_type
      }
      col
    })
    return(x)
  }
  
  if(inherits(x, "numeric")) {
    return(int_type)
  }
  
  stop("x must be a data.frame, tibble, or data.table")
}



