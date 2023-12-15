

#' Berkeley Child Guidance Study Data for males
#'
#' @description A subset of the [bsitar::berkeley] data that contains
#'   longitudinal growth data for 20 randomly selected males (6 to 20 
#'   years of age).
#'
#' @details A detailed description of the full data including the frequency of
#'   measurements per year is provided in the [bsitar::berkeley] data.
#'   
#' @name berkeley_mdata
#' @docType data
#' @format A data frame with 902 observations on the following 3 variables:
#' \describe{
#' \item{id}{factor variable}
#' \item{age}{years, numeric vector}
#' \item{height}{cm, numeric vector}
#' }
#' @references
#'  \insertAllCited{}
#'  
#' @keywords datasets
#' 
#' @return A data frame with 3 columns.
#' @inherit berkeley author
"berkeley_mdata"