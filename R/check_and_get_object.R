

#' Check and get namespace object if exists
#'
#' @param object An object to be retrieved. Note that \code{object} must be a
#'   symbol and not a character string.
#'   
#' @param namespace A character string specifying the namespace to be checked.
#' 
#' @param envir An environment to be used (default \code{NULL}).
#' 
#' @return An object of same class as input \code{object}.
#'
#' @export
#'
#' @inherit berkeley author
#'
#' @examples
#'
#' \donttest{
#' # Check whether model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' }
#'
getNsObject <- function(object,
                        namespace = NULL,
                        envir = NULL) {
  if (is.null(namespace)) {
    namespace <- "bsitar"
  }
  if (grepl("\"", deparse(substitute(object)), fixed = T)) {
    stop("object must be a symbol and not a string")
  }
  object_str <- deparse(substitute(object))
  if (is.null(envir)) {
    if (exists(object_str)) {
      out <- object
    } else {
      out <- utils::getFromNamespace(object_str, namespace)
    }
  }
  
  if (!is.null(envir)) {
    if (exists(object_str, envir = envir)) {
      out <- assign(object_str, object, envir = envir)
    } else {
      out <- utils::getFromNamespace(object_str, namespace, envir = envir)
    }
  }
  
  return(out)
}
