#' GDP
#'
#' Gross domestic product (GDP)
#'
#'
#' The data here was constructed based on the following data:
#' The GDP data are publicly available at
#' \url{https://data.worldbank.org/indicator/NY.GDP.PCAP.CD}.
#' The excel file that we downloaded was called
#' \code{API_NY.GDP.PCAP.CD_DS2_en_excel_v2_3358980.xls}.
#' The inflation rates are publicly available at
#' \url{https://data.worldbank.org/indicator/NY.GDP.DEFL.KD.ZG}.
#' The excel file that we downloaded was called
#' \code{API_NY.GDP.DEFL.KD.ZG_DS2_en_excel_v2_3469555.xls},
#' from there we took only the inflation rates for United States.
#' Both are distributed under the CC-BY 4.0 license (see
#' https://datacatalog.worldbank.org/public-licenses#cc-by).
#'
#' Then we discounted the GDP of every country in the study to the 1960 USD,
#' and we extrapolated the missing values of the GDP of a country using the
#' closest known ratio of the GDP of the country and the median GDP in that year.
#' Further, the missing values of GDP were interpolated using linear
#' interpolation of the two closest ratios. Appendix of the
#' \code{vignette(FDRenvelopes)} includes the code to prepare the
#' \code{curve_set} object.
#'
#'
#' @format A \code{curve_set} object containing the GDP values for different countries.
#'
#' @usage data("GDP")
#' @keywords datasets
#' @name GDP
#' @docType data
NULL
