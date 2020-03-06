% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm.R
\name{graph.flm}
\alias{graph.flm}
\title{Graphical functional GLM}
\usage{
graph.flm(
  nsim,
  formula.full,
  formula.reduced,
  curve_sets,
  factors = NULL,
  contrasts = FALSE,
  savefuns = FALSE,
  ...,
  GET.args = NULL,
  mc.cores = 1L,
  mc.args = NULL,
  cl = NULL,
  fast = TRUE
)
}
\arguments{
\item{nsim}{The number of random permutations.}

\item{formula.full}{The formula specifying the general linear model,
see \code{formula} in \code{\link[stats]{lm}}.}

\item{formula.reduced}{The formula of the reduced model with nuisance factors only. This model
should be nested within the full model.}

\item{curve_sets}{A named list of sets of curves giving the dependent variable (Y), and
possibly additionally all the factors. The dimensions of the elements should
match with each other, i.e. the factor values should be given for each argument value
and each function. If factors are given in the argument \code{factors}, then can also be just
the curve set representing Y. Also \code{\link[fda.usc]{fdata}} objects allowed.}

\item{factors}{A data frame of factors. An alternative way to specify factors when they
are constant for all argument values. The number of rows of the data frame should be equal
to the number of curves. Each column should specify the values of a factor.}

\item{contrasts}{Logical. FALSE and TRUE specify the two test functions as described in
description part of this help file.}

\item{savefuns}{Logical. If TRUE, then the functions from permutations are saved to the attribute
simfuns.}

\item{...}{Additional arguments to be passed to \code{\link[stats]{lm}}. See details.}

\item{GET.args}{A named list of additional arguments to be passed to \code{\link{global_envelope_test}}.}

\item{mc.cores}{The number of cores to use, i.e. at most how many child processes will be run simultaneously.
Must be at least one, and parallelization requires at least two cores. On a Windows computer mc.cores must be 1
(no parallelization). For details, see \code{\link{mclapply}}, for which the argument is passed.
Parallelization can be used in generating simulations and in calculating the second stage tests.}

\item{mc.args}{A named list of additional arguments to be passed to \code{\link{mclapply}}.
Only relevant if \code{mc.cores} is more than 1.}

\item{cl}{Allows parallelization through the use of \code{\link{parLapply}} (works also
in Windows), see the argument \code{cl} there, and examples.}

\item{fast}{Logical. See details.}
}
\value{
A \code{global_envelope} or \code{combined_global_envelope} object,
which can be printed and plotted directly.
}
\description{
Non-parametric graphical tests of significance in functional general linear model (GLM)
}
\details{
The function \code{graph.flm} performs the graphical functional GLM of Mrkvička et al. (2019).
This is a nonparametric graphical test of significance of a covariate in functional GLM.
The test is able to find not only if the factor of interest is significant, but also which
functional domain is responsible for the potential rejection.
In the case of functional multi-way main effect ANOVA or functional main effect ANCOVA models,
the test is able to find which groups differ (and where they differ).
In the case of functional factorial ANOVA or functional factorial ANCOVA models,
the test is able to find which combination of levels (which interactions) differ (and where they differ).
The described tests are global envelope tests applied in the context of GLMs.
The Freedman-Lane algorithm (Freedman and Lane, 1983) is applied to permute the functions
(to obtain the simulations under the null hypothesis of "no effects");
consequently, the test approximately achieves the desired significance level.

The specification of the full and reduced formulas is important. The reduced model should be
nested within the reduced model. The full model should include in addition to the reduced
model the interesting factors whose effects are under investigation.
The implementation to find the coefficients of the interesting factors is based on dummy.coef and
the restrictions there apply.

There are different versions of the implementation depending on the application.
Given that the argument \code{fast} is TRUE, then
\itemize{
\item If all the covariates are constant across the functions, i.e. they can be provided in the
argument \code{factors}, then a linear model is fitted separately by least-squares estimation to
the data at each argument value of the functions fitting a multiple linear model by \code{\link[stats]{lm}}.
The possible extra arguments passed in \code{...} to \code{\link[stats]{lm}} must be of the form that
\code{\link[stats]{lm}} accepts for fitting a multiple linear model. In the basic case, no extra arguments are
needed.
\item If some of the covariates vary across the space and there are user specified extra arguments given in
\code{...}, then the implementation fits a linear model at each argument value of the functions using
\code{\link[stats]{lm}}, which can be rather slow. The arguments \code{...} are passed to \code{\link[stats]{lm}}
for fitting each linear model.
}
By setting \code{fast = FALSE}, it is possible to use the slow version for any case. Usually this is not desired.
}
\examples{
data(rimov)
\donttest{
res <- graph.flm(nsim=19, # Increase the number of simulations for serious analysis!
                 formula.full = Y~Year,
                 formula.reduced = Y~1,
                 curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014))
}
\dontshow{
res <- graph.flm(nsim = 3,
                 formula.full = Y~Year,
                 formula.reduced = Y~1,
                 curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014),
                 GET.args = list(alpha=0.25))
}
plot(res)

# Test if there is a change in the slope in 1994,
# i.e. the full model is T = a + b*year + c*year:group,
res <- graph.flm(nsim = 19, # Increase the number of simulations for serious analysis!
                 formula.full = Y ~ Year + Year:Group,
                 formula.reduced = Y ~ Year,
                 curve_sets = list(Y=rimov),
                 factors = data.frame(Year = 1979:2014,
                                     Group = factor(c(rep(1,times=24), rep(2,times=12)),
                                                    levels=1:2)),
                 contrasts = FALSE)
plot(res)

\donttest{
data(GDPtax)
factors.df <- data.frame(Group = GDPtax$Group, Tax = GDPtax$Profittax)
res.tax_within_group <- graph.flm(nsim = 999,
                                  formula.full = Y~Group+Tax+Group:Tax,
                                  formula.reduced = Y~Group+Tax,
                                  curve_sets = list(Y=GDPtax$GDP),
                                  factors = factors.df)
plot(res.tax_within_group)
}
}
\references{
Mrkvička, T., Roskovec, T. and Rost, M. (2019) A nonparametric graphical tests of significance in functional GLM. Methodology and Computing in Applied Probability. doi: 10.1007/s11009-019-09756-y

Freedman, D., & Lane, D. (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298. doi:10.2307/1391660
}