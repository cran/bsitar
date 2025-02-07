

#' @title Check and Get Namespace Object If Exists
#'
#' @description
#' This function checks if an object exists within a specified namespace and
#' returns the object if it exists. The object must be provided as a symbol
#' (not a character string). The function is designed to facilitate the
#' retrieval of model or other objects from a specified environment or
#' namespace. This function is mainly for internal purposes.
#'
#' @param object A symbol representing the object to be retrieved. The input
#'   must be a symbol (i.e., not a character string) corresponding to an
#'   existing object within the specified namespace.
#'
#' @param namespace A character string specifying the namespace to check for the
#'   object. If the object exists within the given namespace, it will be
#'   returned.
#'
#' @param envir An environment in which to search for the object. If set to
#'   \code{NULL} (default), the function uses the global environment.
#'
#' @return The object of the same class as the input \code{object}, if it
#'   exists. If the object doesn't exist in the specified namespace or
#'   environment, an error is raised.
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
