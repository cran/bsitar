

#' Flag parameters with possible prior-data conflict
#'
#' Screen [prior_sensitivity()] results for parameters showing notable
#' prior-data conflict, likelihood conflict, or both. Parameters flagged for
#' both forms of sensitivity may indicate possible prior-data conflict in the
#' \pkg{priorsense} framework, while parameters flagged primarily for prior
#' sensitivity may reflect weak likelihood informativeness. See
#' [prior_sensitivity()] for details.
#' 
#' @param x An object returned by [prior_sensitivity()].
#'
#' @param threshold_prior Numeric threshold used to flag prior sensitivity.
#'
#' @param threshold_lik Numeric threshold used to flag likelihood sensitivity.
#'
#' @param empty A string that is used to replace \code{NULL} values.
#'
#' @param print A logical to print the object.
#'
#' @param return_table A logical indicating whether to return the table
#'   \code{return_table = TRUE} or the list \code{return_table = FALSE}.
#'
#' @param tibble_table A logical indicating whether to return the table as a
#'   data frame (\code{tibble_table = FALSE}) or as a \code{tibble}
#'   (\code{tibble_table = TRUE}). Default \code{FALSE}. Ignored if
#'   \code{return_table = FALSE}.
#'
#' @param print_table A logical indicating whether to print the table using
#'   [knitr::kable()] by setting \code{print_table = TRUE}, or return it as an
#'   object by setting \code{print_table = FALSE}. If \code{print_table = TRUE},
#'   the table is printed but not returned; the function returns \code{NULL}
#'   invisibly. The default is \code{FALSE}. This argument is ignored when
#'   \code{return_table = FALSE}.
#'
#' @param flex_table A logical indicating whether to return the data frame
#'   (\code{flex_table = FALSE}) or a \code{flextable} object (\code{flex_table
#'   = TRUE}).
#'
#' @param return_file An optional character string specifying the output format.
#'   Supported values are \code{"word"}, \code{"docx"}, \code{"html"},
#'   \code{"png"}, \code{"pdf"}, \code{"svg"}, and \code{"xlsx"}. If
#'   \code{NULL}, the function returns either a data frame or a \code{flextable}
#'   object, depending on the value of \code{flex_table}. This argument is
#'   ignored when \code{path} is specified.
#'
#' @param path An optional character string specifying the output file path. If
#'   \code{NULL}, a default file name is generated based on the value of
#'   \code{return_file}; for example, \code{"table_output.docx"} or
#'   \code{"table_output.xlsx"}. If no file extension is specified and
#'   \code{return_file = NULL}, the default extension is \code{".xlsx"}.
#'   Otherwise, the file extension is set to match \code{return_file}.
#'
#' @param title Optional title used in exported output where supported. For Word
#'   and HTML output, named \code{flextable} objects can be used as document
#'   titles or section titles.
#'
#' @param align Alignment used for Word export. Must be one of \code{"left"},
#'   \code{"center"}, or \code{"right"}.
#'
#' @param sheet_name Character string giving the worksheet name for Excel
#'   output. Default is \code{"table"}.
#'
#' @param ... Ignored.
#'
#' @details
#' This helper is designed as a lightweight screening step after
#' [prior_sensitivity()]. It does not replace graphical inspection. The
#' recommended workflow is:
#' \enumerate{
#'   \item Run \code{prior_sensitivity()} to compute sensitivity.
#'   \item Use \code{prior_conflict()} to identify parameters of
#'         potential concern.
#'   \item Re-run \code{prior_sensitivity()} with \code{plot = "dens"},
#'         \code{plot = "ecdf"}, or \code{plot = "quantities"} for the
#'         flagged variables.
#' }
#' 
#' Else, user can call [prior_conflict()] directly from within the
#' [prior_sensitivity()] by setting the argument \code{return_conflicts} as
#' \code{TRUE}. See [prior_sensitivity()] for details and examples.
#' 
#' @return A list when \code{return_table = FALSE} or a \code{flextable} if
#'   \code{return_table = TRUE}
#' \describe{
#'   \item{prior_flagged}{Character vector of variables flagged for
#'     prior sensitivity.}
#'   \item{likelihood_flagged}{Character vector of all variables flagged for
#'     likelihood sensitivity.}
#'   \item{prior_likelihood_flagged}{Character vector of variables flagged for 
#'     both prior and likelihood sensitivity}
#' }
#'
#' @seealso
#' [prior_sensitivity()]
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid model estimation which can take time, the Bayesian SITAR model fit
#' # to the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' model <- getNsObject(berkeley_exfit)
#' 
#' # Note return_table = FALSE which ensures correct output is returned from
#' # the prior_sensitivity()
#' ps <- prior_sensitivity(
#'   model = model,
#'   variable = c("b_a_Intercept", "sigma"),
#'   return_table = FALSE
#' )
#' 
#' conflict <- prior_conflict(ps, return_table = TRUE)
#' print(conflict)
#' 
#' # Note that you can call prior_conflict() from within the prior_sensitivity()
#' # with return_conflict = TRUE that resturns the same results as above, 
#' conflict2 <- prior_sensitivity(
#'   model = model,
#'   variable = c("b_a_Intercept", "sigma"),
#'   return_conflict = TRUE, return_table = TRUE
#' )
#' stopifnot(identical(conflict, conflict2))
#' }
#'
#' @rdname prior_conflict
#' @export
#' 
prior_conflict.bgmfit <- function(x,
                                  threshold_prior = 0.05,
                                  threshold_lik = 0.05,
                                  empty = "-",
                                  print = FALSE,
                                  return_table = TRUE,
                                  tibble_table = FALSE,
                                  print_table = FALSE,
                                  return_file = NULL,
                                  flex_table = FALSE,
                                  path = NULL,
                                  title = NULL,
                                  align = "center",
                                  sheet_name = "table",
                                  ...) {
  
  if(!inherits(x, "prior_sensitivity")) {
    stop2c("Please run `prior_sensitivity()` with 
           'return_table = FALSE' before calling the 
           `prior_conflict()`")
  }
  
  sens <- x$sensitivity
  
  prior_flagged <- sens %>% 
    dplyr::filter(
      abs(.data[["prior"]]) > threshold_prior
    ) %>% 
    dplyr::pull(.data$variable) %>% 
    unique()
  
  likelihood_flagged <- sens %>% 
    dplyr::filter(
      abs(.data[["likelihood"]]) < threshold_lik
    ) %>% 
    dplyr::pull(.data$variable) %>% 
    unique()
  
  prior_only_params <- setdiff(prior_flagged, likelihood_flagged)
  prior_likelihood_flagged <-  intersect(prior_flagged, likelihood_flagged)
  
  if(is_emptyx(prior_flagged)) prior_flagged <- NULL
  if(is_emptyx(likelihood_flagged)) likelihood_flagged <- NULL
  if(is_emptyx(prior_likelihood_flagged)) prior_likelihood_flagged <- NULL
  
  out_list <- list(
    prior_flagged = prior_flagged,
    likelihood_flagged = likelihood_flagged,
    prior_likelihood_flagged = prior_likelihood_flagged)
  
  if(!return_table) {
    return(out_list)
  }
  
  out_flex <- list_to_flextable(out_list, align = align, empty = empty)
  
  if (!is.null(title)) {
    out_flex <- flextable::set_caption(out_flex, caption = title)
  }
  
  if(print) print(out_flex$body$dataset)
  
  out <- export_flextable(ft = out_flex,
                          return_file = return_file,
                          path = path,
                          title = title,
                          align = align,
                          sheet_name = sheet_name)
  
  if(!is.null(return_file)) {
    return_table <- FALSE
  } 
  
  if(return_table) {
    if(!flex_table) {
      out <- out$body$dataset
      if(tibble_table) out <- out %>% tibble::as_tibble()
      if(print_table) {
        print(knitr::kable(out))
        return(invisible(NULL))
      } else {
         return(out)
      }
    } else if( flex_table) {
      if (!is.null(title)) {
        out <- flextable::set_caption(out, caption = title)
      }
      return(out)
    }
  }
  
  if(!return_table) {
    export_flextable(ft = out_flex,
                     return_file = return_file,
                     path = path,
                     title = title,
                     align = align,
                     sheet_name = sheet_name)
  }

}





#' @rdname prior_conflict
#' @export
prior_conflict <- function(x, ...) {
  UseMethod("prior_conflict")
}


