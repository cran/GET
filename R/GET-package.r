#' Global Envelopes
#'
#' The \pkg{GET} package provides implementation of global envelopes
#' for a set of general d-dimensional vectors T in various applications.
#' A 100(1-alpha)% global envelope is a band bounded by two vectors such that
#' the probability that T falls outside this envelope in any of the d points is
#' equal to alpha. Global means that the probability is controlled simultaneously
#' for all the d elements of the vectors.
#' The global envelopes can be used for central regions of functional or
#' multivariate data (e.g. outlier detection, functional boxplot),
#' for graphical Monte Carlo and permutation tests where the test statistic
#' is a multivariate vector or function (e.g. goodness-of-fit testing for point
#' patterns and random sets, functional ANOVA, functional GLM, n-sample test of
#' correspondence of distribution functions), and for global confidence and
#' prediction bands (e.g. confidence band in polynomial regression,
#' Bayesian posterior prediction).
#'
#'
#' The \pkg{GET} package provides central regions (i.e. global envelopes) and
#' global envelope tests with intrinsic graphical interpretation.
#' The central regions can be constructed from (functional) data.
#' The tests are Monte Carlo or permutation tests, which demand simulations
#' from the tested null model. The methods are applicable for any multivariate
#' vector data and functional data (after discretization).
#'
#' To get an overview of the package, start R and type \code{library("GET")} and \code{vignette("GET")}.
#'
#' To get examples of point pattern analysis, start R and type \code{library("GET")} and \code{vignette("pointpatterns")}.
#'
#' To get examples of Mrkvička and Myllymäki (2022), start R and type \code{library("GET")} and \code{vignette("FDRenvelopes")}.
#'
#' @section Key functions in \pkg{GET}:
#' \itemize{
#' \item \emph{Central regions} or \emph{global envelopes} or \emph{confidence bands}:
#' \code{\link{central_region}}.
#' E.g. 50\% central region of growth curves of girls \code{\link[fda]{growth}}.
#' \itemize{
#'            \item First create a curve_set of the growth curves, e.g.
#'
#'                  \code{
#'                    cset <- curve_set(r = as.numeric(row.names(growth$hgtf)),
#'                                      obs = growth$hgtf)
#'                  }
#'            \item Then calculate 50\% central region (see \code{\link{central_region}} for further arguments)
#'
#'                  \code{
#'                    cr <- central_region(cset, coverage = 0.5)
#'                  }
#'            \item Plot the result (see \code{\link{plot.global_envelope}} for plotting options)
#'
#'                  \code{
#'                    plot(cr)
#'                  }
#' }
#' It is also possible to do combined central regions for several sets of curves provided in a list
#' for the function, see examples in \code{\link{central_region}}.
#'
#' \item \emph{Global envelope tests}: \code{\link{global_envelope_test}} is the main function.
#' E.g. A test of complete spatial randomness (CSR) for a point pattern \code{X}:
#'
#' \code{X <- spruces # an example pattern from spatstat}
#'
#' \itemize{
#'            \item Use the function \code{envelope} of \pkg{spatstat} to create nsim simulations
#'                  under CSR and to calculate the functions you want (below K-functions by Kest).
#'                  Important: use the option 'savefuns=TRUE' and
#'                  specify the number of simulations \code{nsim}.
#'
#'                  \code{
#'                    env <- envelope(X, nsim=999, savefuns = TRUE, fun = Kest, simulate = expression(runifpoint(ex = X)))
#'                  }
#'            \item Perform the test (see \code{\link{global_envelope_test}} for further arguments)
#'
#'                  \code{
#'                    res <- global_envelope_test(env)
#'                  }
#'            \item Plot the result (see \code{\link{plot.global_envelope}} for plotting options)
#'
#'                  \code{
#'                    plot(res)
#'                  }
#' }
#' It is also possible to do combined global envelope tests for several sets of curves provided in a list
#' for the function, see examples in \code{\link{global_envelope_test}}.
#' To obtain false discovery rate envelopes of Mrkvička and Myllymäki (2023) use the argument \code{typeone = "fdr"}.
#' }
#'
#' \itemize{
#'  \item \emph{Functional ordering}: \code{\link{central_region}} and \code{\link{global_envelope_test}}
#'  are based on different measures for ordering the functions (or vectors) from
#'  the most extreme to the least extreme ones. The core functionality of calculating the measures
#'  is in the function \code{\link{forder}}, which can be used to obtain different measures for sets of
#'  curves. Usually there is no need to call \code{\link{forder}} directly.
#' \item \emph{Functional boxplots}: \code{\link{fBoxplot}}
#' \item \emph{Adjusted} global envelope tests for composite null hypotheses
#' \itemize{
#'   \item \code{\link{GET.composite}}, see a detailed example in \code{\link{saplings}}
#' }
#' \item \emph{One-way functional ANOVA}:
#'  \itemize{
#'   \item \emph{Graphical} functional ANOVA tests: \code{\link{graph.fanova}}
#'   \item Global rank envelope based on F-values: \code{\link{frank.fanova}}
#'  }
#' \item \emph{Functional general linear model (GLM)}:
#'  \itemize{
#'   \item \emph{Graphical} functional GLM: \code{\link{graph.flm}}
#'   \item Global rank envelope based on F-values: \code{\link{frank.flm}}
#'   \item For large data (not fitting comfortably in memory): \code{\link{partial_forder}}
#'  }
#' \item \emph{Functional clustering}: \code{\link{fclustering}}
#' \item \emph{Global quantile regression}: \code{\link{global_rq}}
#' \item Functions for performing global envelopes for other specific purposes:
#'  \itemize{
#'   \item Graphical n sample test of correspondence of distribution functions: \code{\link{GET.distrequal}}
#'   \item Permutation-based tests of independence to samples from any bivariate distribution: \code{\link{GET.distrindep}}
#'   \item Testing global and local dependence of point patterns on covariates: \code{\link{GET.spatialF}}
#'   \item Testing local correlations: \code{\link{GET.localcor}}
#'   \item Variogram and residual variogram with global envelopes: \code{\link{GET.variogram}}
#'  }
#'
#' \item Deviation tests (for simple hypothesis): \code{\link{deviation_test}} (no graphical
#' interpretation)
#' \item Most functions accept the curves provided in a \code{curve_set} object.
#' Use \code{\link{curve_set}} to create a \code{curve_set} object from the
#' functions. Other formats to provide the curves to the above functions are
#' also accepted, see the information on the help pages.
#' }
#' See the help files of the functions for examples.
#'
#' @section Workflow for (single hypothesis) tests based on single functions:
#'
#' To perform a test you always first need to obtain the test function \eqn{T(r)}{T(r)}
#' for your data (\eqn{T_1(r)}{T_1(r)}) and for each simulation
#' (\eqn{T_2(r), \dots, T_{s+1}(r)}{T_2(r), ..., T_{s+1}(r)}) in one way or another.
#' Given the set of the functions \eqn{T_i(r), i=1, \dots, s+1}{T_i(r), i=1,...,s+1},
#' you can perform a test by \code{\link{global_envelope_test}}.
#'
#' 1) The workflow when using your own programs for simulations:
#'
#' \itemize{
#' \item (Fit the model and) Create \eqn{s}{s} simulations from the (fitted) null model.
#' \item Calculate the functions \eqn{T_1(r), T_2(r), \dots, T_{s+1}(r)}{T_1(r), T_2(r), ..., T_{s+1}(r)}.
#' \item Use \code{\link{curve_set}} to create a \code{curve_set} object
#'       from the functions \eqn{T_i(r), i=1, \dots, s+1}{T_i(r), i=1,...,s+1}.
#' \item Perform the test
#'
#'       \code{res <- global_envelope_test(curve_set)}
#'
#' where \code{curve_set} is the 'curve_set'-object you created, and plot the result
#'
#'       \code{plot(res)}
#' }
#'
#' 2) The workflow utilizing \pkg{spatstat}: start R, type \code{library("GET")} and \code{vignette("pointpatterns")},
#' which explains the workflow and gives many examples of point pattern analysis
#'
#' @section Functions for modifying sets of functions:
#' It is possible to modify the curve set \eqn{T_1(r), T_2(r), \dots, T_{s+1}(r)}{T_1(r), T_2(r), ..., T_{s+1}(r)} for the test.
#'
#' \itemize{
#' \item You can choose the interval of distances \eqn{[r_{\min}, r_{\max}]}{[r_min, r_max]} by \code{\link{crop_curves}}.
#' \item For better visualisation, you can take \eqn{T(r)-T_0(r)}{T(r)-T_0(r)} by \code{\link{residual}}.
#' Here \eqn{T_0(r)}{T_0(r)} is the expectation of \eqn{T(r)}{T(r)} under the null hypothesis.
#' }
#'
#' @section Example data (see references on the help pages of each data set):
#' \itemize{
#'  \item \code{\link{abide_9002_23}}: see help page
#'  \item \code{\link{adult_trees}}: a point pattern of adult rees
#'  \item \code{\link{cgec}}: centred government expenditure centralization (GEC) ratios (see \code{\link{graph.fanova}})
#'  \item \code{\link{fallen_trees}}: a point pattern of fallen trees
#'  \item \code{\link{GDPtax}}: GDP per capita with country groups and other covariates
#'  \item \code{\link{imageset3}}: a simulated set of images
#'  \item \code{\link{rimov}}: water temperature curves in 365 days of the 36 years
#'  \item \code{\link{saplings}}: a point pattern of saplings (see \code{\link{GET.composite}})
#' }
#' The data sets are used to show examples of the functions of the library.
#'
#' @section Number of functions:
#'
#' If the number of functions is low, the choice of the measure (or type or depth) playes a role,
#' as explained in \code{vignette("GET")} (Section 2.4).
#'
#' Note that the recommended minimum number of simulations for the rank
#' envelope test (Myllymäki et al., 2017) based on a single function in spatial statistics is nsim=2499.
#' When the number of argument values is large, also larger number simulations is needed in order to
#' have a narrow p-interval.
#' The "erl", "cont", "area", "qdir" and "st" global envelope tests and deviation tests can be
#' used with a lower number of simulations, although the Monte Carlo error is obviously larger
#' with a lower number of simulations.
#' For increasing the number of simulations, all the global rank envelopes approach the same curves.
#'
#' Mrkvička et al. (2017) discussed the number of simulations for tests based on many functions.
#'
#' @section Documentation:
#' Myllymäki and Mrkvička (2024) provides description of the package.
#' The material can also be found in the corresponding vignette, which is available by
#' starting R and typing \code{library("GET")} and \code{vignette("GET")}.
#'
#' In the special case of spatial processes (spatial point processes, random sets),
#' the functions are typically estimators of summary functions. The package supports
#' the use of the R package \pkg{spatstat} for generating simulations and calculating
#' estimators of the chosen summary function, but alternatively these can be done by
#' any other way, thus allowing for any user-specified models/functions.
#' To see examples of global envelopes for analysing point pattern data,
#' start R, type \code{library("GET")} and \code{vignette("pointpatterns")}.
#'
#' Mrkvička and Myllymäki (2023) developed false discovery rate (FDR) envelopes.
#' Examples can be found by in associated vignette: start R, and type
#' \code{library("GET")} and \code{vignette("pointpatterns")}.
#'
#' Mrkvička et al. (2023a) proposed global quantile regression. An example of
#' global quantile regression is given in the vignette \code{vignette("QuantileRegression")}.
#'
#' The vignette \code{vignette("HotSpots")} illustrates the methodology proposed by
#' Mrkvička et al. (2023b) for detecting hotspots on a linear network.
#'
#' Type citation("GET") to get a full list of references.
#'
#' @section Acknowledgements:
#'
#' Mikko Kuronen has made substantial contributions of code.
#' Additional contributions and suggestions from Jiří Dvořák, Pavel Grabarnik,
#' Ute Hahn, Michael Rost and Henri Seijo.
#'
#' @author
#' Mari Myllymäki (mari.myllymaki@@luke.fi, mari.j.myllymaki@@gmail.com) and
#' Tomáš Mrkvička (mrkvicka.toma@@gmail.com)
#'
#' @references
#' Dai, W., Athanasiadis, S., Mrkvička, T. (2021) A new functional clustering method with combined dissimilarity sources and graphical interpretation. Intech open, London, UK. doi: 10.5772/intechopen.100124
#'
#' Dvořák, J. and Mrkvička, T. (2022). Graphical tests of independence for general distributions. Computational Statistics 37, 671--699.
#'
#' Konstantinou, K., Mrkvička, T. and Myllymäki, M. (2024) The power of visualizing distributional differences: formal graphical n-sample tests. Computational Statistics. doi: 10.1007/s00180-024-01569-z
#'
#' Mrkvička, T., Konstantinou, K., Kuronen, M. and Myllymäki, M. (2023a) Global quantile regression. arXiv:2309.04746 [stat.ME]. https://doi.org/10.48550/arXiv.2309.04746
#'
#' Mrkvička, T., Kraft, S., Blažek, V. and Myllymäki, M. (2023b) Hotspots detection on a linear network with presence of covariates: a case study on road crash data. Available at SSRN: http://dx.doi.org/10.2139/ssrn.4598454
#'
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017) Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27(5), 1239-1255. doi: 10.1007/s11222-016-9683-9
#'
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2020) A one-way ANOVA test for functional data with graphical interpretation. Kybernetika 56(3), 432-458. doi: 10.14736/kyb-2020-3-0432
#'
#' Mrkvička, T., Myllymäki, M., Kuronen, M. and Narisetty, N. N. (2022) New methods for multiple testing in permutation inference for the general linear model. Statistics in Medicine 41(2), 276-297. doi: 10.1002/sim.9236
#'
#' Mrkvička, T., Myllymäki, M. (2023) False discovery rate envelopes. Statistics and Computing 33, 109. https://doi.org/10.1007/s11222-023-10275-7
#'
#' Mrkvička, T., Roskovec, T. and Rost, M. (2021) A nonparametric graphical tests of significance in functional GLM. Methodology and Computing in Applied Probability 23, 593-612. doi: 10.1007/s11009-019-09756-y
#'
#' Mrkvička, T., Soubeyrand, S., Myllymäki, M., Grabarnik, P., and Hahn, U. (2016) Monte Carlo testing in spatial statistics, with applications to spatial residuals. Spatial Statistics 18, Part A, 40-53. doi: 10.1016/j.spasta.2016.04.005
#'
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015) Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11, 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017) Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 79, 381-404. doi: 10.1111/rssb.12172
#'
#' Myllymäki, M. and Mrkvička, T. (2024). GET: Global envelopes in R. Journal of Statistical Software 111(3), 1-40. doi: 10.18637/jss.v111.i03
#'
#' Myllymäki, M., Kuronen, M. and Mrkvička, T. (2020). Testing global and local dependence of point patterns on covariates in parametric models. Spatial Statistics 42, 100436. doi: 10.1016/j.spasta.2020.100436
#'
#' @name GET-package
#' @aliases GET GET-package
NULL
