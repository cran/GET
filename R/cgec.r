#' Centred government expenditure centralization ratios
#'
#' Centred government expenditure centralization (GEC) ratios
#'
#'
#' The data includes the government expenditure centralization (GEC) ratio in percent that has been
#' centred with respect to country average in order to remove the differences in absolute values of
#' GEC.
#' The GEC ratio is the ratio of central government expenditure to the total general government expenditure.
#' Data were collected from the Eurostat (2018) database.
#' Only those European countries were included, where the data were available from 1995 to 2016 without
#' interruption. Finally, 29 countries were classified into three groups in the following way:
#' \itemize{
#' \item Group 1: Countries joining EC between 1958 and 1986 (Belgium, Denmark, France, Germany
#' (until 1990 former territory of the FRG), Greece, Ireland, Italy, Luxembourg, Netherlands, Portugal,
#' Spain, United Kingdom. These countries have long history of European integration, representing the
#' core of integration process.
#' \item Group 2: Countries joining the EU in 1995 (Austria, Sweden, Finland) and 2004 (Malta, Cyprus),
#' except CEEC (separate group), plus highly economically integrated non-EU countries, EFTA members
#' (Norway, Switzerland). Countries in this group have been, or in some case even still are standing
#' apart from the integration mainstream. Their level of economic integration is however very high.
#' \item Group 3: Central and Eastern European Countries (CEEC), having similar features in political
#' end economic history. The process of economic and political integration have been initiated by
#' political changes in 1990s. CEEC joined the EU in 2004 and 2007 (Bulgaria, Czech Republic, Estonia,
#' Hungary, Latvia, Lithuania, Poland, Romania, Slovakia, Slovenia, data for Croatia joining in 2013 are
#' incomplete, therefore not included).
#' }
#' This grouping is used in examples.
#'
#' @format A list of two components. The first one is the \code{curve_set} object containing the observed
#' values of centred GEC observed in year 1995-2016 for the above countries.
#' The second component \code{group} gives the grouping.
#'
#' @usage data("cgec")
#' @references
#' Eurostat (2018). "Government revenue, expenditure and main aggregates (gov10amain)”. Retrieved from https://ec.europa.eu/eurostat/data/database(26/10/2018).
#'
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2020) A one-way ANOVA test for functional data with graphical interpretation. Kybernetika 56 (3), 432-458. doi: 10.14736/kyb-2020-3-0432
#' @keywords datasets
#' @name cgec
#' @docType data
#' @seealso \code{\link{graph.fanova}}
#' @examples
#' data("cgec")
#' # Plot data in groups
#' for(i in 1:3)
#'   assign(paste0("p", i), plot(subset(cgec$cgec, cgec$group == i)) +
#'     ggplot2::labs(title=paste("Group ", i, sep=""), y="Centred GEC"))
#' p3
#' if(require("patchwork", quietly=TRUE))
#'   p1 + p2 + p3
NULL
