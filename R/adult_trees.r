#' Adult trees data set
#'
#' Adult trees data set
#'
#'
#' A pattern of large trees (height > 25 m) originating from an uneven aged multi-species
#' broadleaf nonmanaged forest in Kaluzhskie Zaseki, Russia.
#'
#' The pattern is a sample part of data collected over 10 ha plot as a part of a research
#' program headed by project leader Prof. O.V. Smirnova.
#'
#' @format A \code{data.frame} containing the locations (x- and y-coordinates) of 67 trees
#' in an area of 75 m x 75 m.
#'
#' @usage data("adult_trees")
#' @references
#' Grabarnik, P. and Chiu, S. N. (2002) Goodness-of-fit test for complete spatial randomness against
#' mixtures of regular and clustered spatial point processes. Biometrika, 89, 411–421.
#'
#' van Lieshout, M.-C. (2010) Spatial point process theory. In Handbook of Spatial Statistics (eds. A. E.
#' Gelfand, P. J. Diggle, M. Fuentes and P. Guttorp), Handbooks of Modern Statistical Methods. Boca
#' Raton: CRC Press.
#'
#' @keywords datasets
#' @keywords spatial
#' @name adult_trees
#' @docType data
#' @seealso \code{\link{saplings}}
#' @examples
#' if(require("spatstat.geom", quietly=TRUE)) {
#'   data("adult_trees")
#'   adult_trees <- as.ppp(adult_trees, W = square(75))
#'   plot(adult_trees)
#' }
#'
NULL
