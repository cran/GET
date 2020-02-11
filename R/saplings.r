#' Saplings data set
#'
#' Saplings data set
#'
#'
#' A pattern of small trees (height <= 15 m) originating from an uneven aged multi-species
#' broadleaf nonmanaged forest in Kaluzhskie Zaseki, Russia.
#'
#' The pattern is a sample part of data collected over 10 ha plot as a part of a research
#' program headed by project leader Prof. O.V. Smirnova.
#'
#' @format An object of class \code{\link[spatstat]{ppp.object}} representing the point
#' pattern of tree locations.
#'
#' @usage data(saplings)
#' @references
#' Grabarnik, P. and Chiu, S. N. (2002) Goodness-of-fit test for complete spatial randomness against
#' mixtures of regular and clustered spatial point processes. \emph{Biometrika}, \bold{89}, 411–421.
#'
#' van Lieshout, M.-C. (2010) Spatial point process theory. In Handbook of Spatial Statistics (eds. A. E.
#' Gelfand, P. J. Diggle, M. Fuentes and P. Guttorp), Handbooks of Modern Statistical Methods. Boca
#' Raton: CRC Press.
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381-404. doi: 10.1111/rssb.12172
#'
#' @keywords datasets
#' @keywords spatial
#' @name saplings
#' @docType data
#' @examples
#' # This is an example analysis of the saplings data set
#' #=====================================================
#' # Example of Myllymaki et al. (2017, Supplement S4).
#' if(require("spatstat", quietly=TRUE)) {
#'   data(saplings)
#'
#'   # First choose the r-distances for L (r) and J (rJ) functions, respectively.
#'   nr <- 500
#'   rmin <- 0.3; rminJ <- 0.3
#'   rmax <- 10; rmaxJ <- 6
#'   rstep <- (rmax-rmin)/nr; rstepJ <- (rmaxJ-rminJ)/nr
#'   r <- seq(0, rmax, by=rstep)
#'   rJ <- seq(0, rmaxJ, by=rstepJ)
#'
#'   #-- CSR test --# (a simple hypothesis)
#'   #--------------#
#'   # First, a CSR test using the L(r)-r function:
#'   # Note: CSR is simulated by fixing the number of points and generating nsim simulations
#'   # from the binomial process, i.e. we deal with a simple hypothesis.
#'   \donttest{nsim <- 999 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations for testing}
#'   env <- envelope(saplings, nsim=nsim,
#'    simulate=expression(runifpoint(saplings$n, win=saplings$window)), # Simulate CSR
#'    fun="Lest", correction="translate", # T(r) = estimator of L with translational edge correction
#'    transform = expression(.-r),        # Take the L(r)-r function instead of L(r)
#'    r=r,                                # Specify the distance vector
#'    savefuns=TRUE)                      # Save the estimated functions
#'   # Crop the curves to the interval of distances [rmin, rmax]
#'   # (at the same time create a curve_set from 'env')
#'   curve_set <- crop_curves(env, r_min = rmin, r_max = rmax)
#'   # Perform a global envelope test
#'   res <- global_envelope_test(curve_set, type="erl") # type="rank" and larger nsim was used in S4.
#'   # Plot the result.
#'   plot(res, ylab=expression(italic(hat(L)(r)-r)))
#'
#'   # -> The CSR hypothesis is clearly rejected and the rank envelope indicates clear
#'   # clustering of saplings. Next we explore the Matern cluster process as a null model.
#' }
#' \donttest{
#' if(require("spatstat", quietly=TRUE)) {
#'   #-- Testing the Matern cluster process --# (a composite hypothesis)
#'   #----------------------------------------#
#'   # Fit the Matern cluster process to the pattern (using minimum contrast estimation with the pair
#'   # correction function)
#'   fitted_model <- kppm(saplings~1, clusters = "MatClust", statistic="pcf")
#'   summary(fitted_model)
#'
#'   nsim <- 19 # 19 just for experimenting with the code!!
#'   #nsim <- 499 # 499 is ok for type = 'qdir' (takes > 1 h)
#'
#'   # Make the adjusted directional quantile global envelope test using the L(r)-r function
#'   # (For the rank envelope test, choose type = "rank" instead and increase nsim.)
#'   adjenvL <- GET.composite(X = fitted_model,
#'                      fun="Lest", correction="translate",
#'                      transform = expression(.-r), r=r,
#'                      type = "qdir", nsim = nsim, nsimsub = nsim,
#'                      r_min=rmin, r_max=rmax)
#'   # Plot the test result
#'   plot(adjenvL, ylab=expression(italic(L(r)-r)))
#'
#'   # From the test with the L(r)-r function, it appears that the Matern cluster model would be
#'   # a reasonable model for the saplings pattern.
#'   # To further explore the goodness-of-fit of the Matern cluster process, test the
#'   # model with the J function:
#'   # This takes quite some time if nsim is reasonably large.
#'   adjenvJ <- GET.composite(X = fitted_model,
#'                      fun="Jest", correction="none", r=rJ,
#'                      type = "qdir", nsim = nsim, nsimsub = nsim,
#'                      r_min=rminJ, r_max=rmaxJ)
#'   # Plot the test result
#'   plot(adjenvJ, ylab=expression(italic(J(r))))
#'   # -> the Matern cluster process not adequate for the saplings data
#'
#'   # Test with the two test functions jointly
#'   adjenvLJ <- GET.composite(X = fitted_model,
#'                      testfuns = list(L = list(fun="Lest", correction="translate",
#'                                          transform = expression(.-r), r=r),
#'                                      J = list(fun="Jest", correction="none", r=rJ)),
#'                      type = "erl", nsim = nsim, nsimsub = nsim,
#'                      r_min=c(rmin, rminJ), r_max=c(rmax, rmaxJ),
#'                      save.cons.envelope=TRUE)
#'   plot(adjenvLJ)
#' }}
NULL
