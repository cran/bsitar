

#' @title Berkeley Child Guidance Study Data
#'
#' @description Longitudinal growth records for 136 children.
#'
#' @details The data was originally included as an appendix in the book
#'   *Physical Growth of California Boys and Girls from Birth to Eighteen Years*
#'   authored by \insertCite{Tuddenham1954;textual}{bsitar}. The dataset was
#'   later used as an example in the \pkg{sitar} \insertCite{R-sitar}{bsitar}
#'   package after correcting transcription errors.
#' 
#'   A more detailed description, including the frequency of measurements per
#'   year, is provided in the \pkg{sitar} package \insertCite{R-sitar}{bsitar}.
#'   Briefly, the data consists of repeated growth measurements made on 66 boys
#'   and 70 girls (ages 0 to 21). The children were born in 1928-29 in Berkeley,
#'   California, and were of northern European ancestry. Measurements were taken
#'   at the following ages:
#'   - 0 years (at birth), 
#'   - 0.085 years, 
#'   - 0.25 to 2 years (every 3 months), 
#'   - 2 to 8 years (annually), 
#'   - 8 to 21 years (every 6 months).
#'
#'   The data includes measurements for height, weight (undressed), stem length,
#'   biacromial diameter, bi-iliac diameter, leg circumference, and
#'   \code{dynamometric} strength.
#' 
#' @name berkeley
#' @docType data
#' @format A data frame with 4884 observations on the following 10 variables:
#' \describe{
#' \item{id}{Factor variable with levels 201-278 for males and 301-385 for females.}
#' \item{age}{Age in years (numeric vector).}
#' \item{height}{Height in cm (numeric vector).}
#' \item{weight}{Weight in kg (numeric vector).}
#' \item{stem.length}{Stem length in cm (numeric vector).}
#' \item{bi.acromial}{Biacromial diameter in cm (numeric vector).}
#' \item{bi.iliac}{Bi-iliac diameter in cm (numeric vector).}
#' \item{leg.circ}{Leg circumference in cm (numeric vector).}
#' \item{strength}{\code{Dynamometric} strength in pounds (numeric vector).}
#' \item{sex}{Factor variable with level 1 for male and level 2 for female.}
#' }
#' @references
#'  \insertAllCited{}
#'  
#' @keywords datasets
#' 
#' @return A data frame with 10 columns.
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
"berkeley"
