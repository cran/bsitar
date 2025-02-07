

#' title The Berkeley Child Guidance Study
#'
#' @description The Berkeley Child Guidance Study dataset contains longitudinal
#' anthropometry data for 136 children from birth to 21 years.
#'
#' The data were collected from 66 boys and 70 girls from Berkeley, California,
#' born in 1928-29, of north European ancestry, and followed from birth to 21
#' years. Measurements were taken at the following intervals:
#' - 0, 0.085, 0.25 to 2 (3-monthly),
#' - 2 to 8 years (annually),
#' - 8 to 21 years (6-monthly).
#'
#' The dataset includes measurements for several anthropometric variables,
#' including:
#' - height,
#' - weight (undressed),
#' - stem length,
#' - biacromial diameter,
#' - bi-iliac diameter,
#' - leg circumference,
#' - dynamometric strength.
#'
#' The data were originally provided as an appendix to the book by Tuddenham and
#' Snyder (1954), with a few transcription errors corrected here. Additionally,
#' 19 errors in height and weight reported in issue #7 of the \code{sitar}
#' package have been corrected.
#'
#' The \code{growth} dataset in the \code{fda} package uses heights from the
#' same study.
#'
#' @details This dataset provides valuable information for studying growth
#' patterns in children, and it can be used for various analyses of
#' anthropometric data, including modeling growth curves and studying the
#' effects of early life factors on long-term development.
#' 
#' @name berkeley
#' @docType data
#' @format A data frame with 4884 observations on the following 10 variables:
#' \describe{
#' \item{id}{factor with levels 201-278 male and 301-385 female}
#' \item{age}{years, numeric vector}
#' \item{height}{cm, numeric vector}
#' \item{weight}{kg, numeric vector}
#' \item{stem.length}{cm, numeric vector}
#' \item{bi.acromial}{cm, numeric vector}
#' \item{bi.iliac}{cm, numeric vector}
#' \item{leg.circ}{cm, numeric vector}
#' \item{strength}{lb, numeric vector}
#' \item{sex}{factor with level 1 male and level 2 female}
#' }
#' 
#' @references Tuddenham RD, Snyder MM. Physical growth of California boys and
#' girls from birth to eighteen years. University of California Publications in
#' Child Development 1954;1:183-364.
#'
#' @keywords datasets
#' @export
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' "berkeley"

