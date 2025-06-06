% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/appl_hotspots.r
\name{hotspots.MatClustlpp}
\alias{hotspots.MatClustlpp}
\title{Hot spots detection for Matern point process
- parallel version}
\usage{
hotspots.MatClustlpp(
  PP,
  formula,
  R,
  alpha,
  data,
  sigma = 250,
  nsim = 10000,
  ncores = 1L,
  ...
)
}
\arguments{
\item{PP}{Input, a point pattern object (ppp) of spatstat.}

\item{formula}{An R formula to estimate the first order model.
This formula can contain objects of full size. \code{PP} should be on the
left side of the formula.}

\item{R}{Cluster radius parameter of the Matern cluster process.}

\item{alpha}{Parameter related to mean number of points per cluster.}

\item{data}{Data from where the formula takes objects.
Must be acceptable by the function lppm of spatstat.linnet.}

\item{sigma}{To be passed to density.lpp.}

\item{nsim}{Number of simulations to be performed.}

\item{ncores}{Number of cores used for computations. Default to 1.
If NULL, all available cores are used.}

\item{...}{Additional parameters to be passed to \code{\link{fdr_envelope}}.}
}
\description{
See the hotspots vignette available by starting R, typing
\code{library("GET")} and \code{vignette("GET")}.
}
\references{
Mrkvička et al. (2023). Hotspot detection on a linear network in the presence of covariates: A case study on road crash data. DOI: 10.2139/ssrn.4627591
}
