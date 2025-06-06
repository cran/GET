% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom_central_region.R
\name{geom_central_region}
\alias{geom_central_region}
\alias{stat_central_region}
\title{Central region plot}
\usage{
geom_central_region(
  mapping = NULL,
  data = NULL,
  stat = "CentralRegion",
  position = "identity",
  ...,
  coverage = 0.5,
  type = "erl",
  filled = TRUE,
  drawcenterline = TRUE,
  colours = grey.colors(length(coverage), start = 0.9, end = 0.5),
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)

stat_central_region(
  mapping = NULL,
  data = NULL,
  position = "identity",
  ...,
  coverage = 0.5,
  type = "erl",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{\link[ggplot2:aes]{aes()}}. If specified and
\code{inherit.aes = TRUE} (the default), it is combined with the default mapping
at the top level of the plot. You must supply \code{mapping} if there is no plot
mapping.}

\item{data}{The data to be displayed in this layer. There are three
options:

If \code{NULL}, the default, the data is inherited from the plot
data as specified in the call to \code{\link[ggplot2:ggplot]{ggplot()}}.

A \code{data.frame}, or other object, will override the plot
data. All objects will be fortified to produce a data frame. See
\code{\link[ggplot2:fortify]{fortify()}} for which variables will be created.

A \code{function} will be called with a single argument,
the plot data. The return value must be a \code{data.frame}, and
will be used as the layer data. A \code{function} can be created
from a \code{formula} (e.g. \code{~ head(.x, 10)}).}

\item{stat}{The statistical transformation to use on the data for this layer.
When using a \verb{geom_*()} function to construct a layer, the \code{stat}
argument can be used the override the default coupling between geoms and
stats. The \code{stat} argument accepts the following:
\itemize{
\item A \code{Stat} ggproto subclass, for example \code{StatCount}.
\item A string naming the stat. To give the stat as a string, strip the
function name of the \code{stat_} prefix. For example, to use \code{stat_count()},
give the stat as \code{"count"}.
\item For more information and other ways to specify the stat, see the
\link[ggplot2:layer_stats]{layer stat} documentation.
}}

\item{position}{A position adjustment to use on the data for this layer. This
can be used in various ways, including to prevent overplotting and
improving the display. The \code{position} argument accepts the following:
\itemize{
\item The result of calling a position function, such as \code{position_jitter()}.
This method allows for passing extra arguments to the position.
\item A string naming the position adjustment. To give the position as a
string, strip the function name of the \code{position_} prefix. For example,
to use \code{position_jitter()}, give the position as \code{"jitter"}.
\item For more information and other ways to specify the position, see the
\link[ggplot2:layer_positions]{layer position} documentation.
}}

\item{...}{Other arguments passed on to \code{\link[ggplot2:layer]{layer()}}'s \code{params} argument. These
arguments broadly fall into one of 4 categories below. Notably, further
arguments to the \code{position} argument, or aesthetics that are required
can \emph{not} be passed through \code{...}. Unknown arguments that are not part
of the 4 categories below are ignored.
\itemize{
\item Static aesthetics that are not mapped to a scale, but are at a fixed
value and apply to the layer as a whole. For example, \code{colour = "red"}
or \code{linewidth = 3}. The geom's documentation has an \strong{Aesthetics}
section that lists the available options. The 'required' aesthetics
cannot be passed on to the \code{params}. Please note that while passing
unmapped aesthetics as vectors is technically possible, the order and
required length is not guaranteed to be parallel to the input data.
\item When constructing a layer using
a \verb{stat_*()} function, the \code{...} argument can be used to pass on
parameters to the \code{geom} part of the layer. An example of this is
\code{stat_density(geom = "area", outline.type = "both")}. The geom's
documentation lists which parameters it can accept.
\item Inversely, when constructing a layer using a
\verb{geom_*()} function, the \code{...} argument can be used to pass on parameters
to the \code{stat} part of the layer. An example of this is
\code{geom_area(stat = "density", adjust = 0.5)}. The stat's documentation
lists which parameters it can accept.
\item The \code{key_glyph} argument of \code{\link[ggplot2:layer]{layer()}} may also be passed on through
\code{...}. This can be one of the functions described as
\link[ggplot2:draw_key]{key glyphs}, to change the display of the layer in the legend.
}}

\item{coverage}{A number between 0 and 1. The 100*coverage\% central region will be calculated.
A vector of values can also be provided, leading to the corresponding number of central regions.}

\item{type}{The options and details for \code{type} are given in \code{\link{central_region}}.}

\item{filled}{Boolean. Should the ribbon be filled?}

\item{drawcenterline}{Boolean. Should the center line be drawn?}

\item{colours}{Colours for different coverage levels}

\item{na.rm}{If \code{FALSE}, the default, missing values are removed with
a warning. If \code{TRUE}, missing values are silently removed.}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{\link[ggplot2:borders]{borders()}}.}
}
\description{
\code{geom_central_region} and \code{stat_central_region} can be used to compute
and plot \code{central_region} from data arranged in a \code{data.frame}.
}
\details{
Plots of central regions (global envelopes) with the specified \code{coverage}
and \code{type} (see \code{\link{central_region}}).
When splitting the set of functions to groups by aesthetics or facets, see
examples, the central regions are constructed separately for each group,
each having the specified \code{coverage}.

If Nfunc*(1-coverage) < 1, where Nfunc is the number of functions/curves,
the curves are plotted instead of any region.
}
\section{Aesthetics}{

\code{geom_central_region} requires \code{x}, \code{y} and \code{curveid}.
Additionally \code{geom_central_region} uses the same aesthetics as
\code{\link[ggplot2]{geom_ribbon}} if \code{filled==TRUE} and
\code{\link[ggplot2]{geom_line}} otherwise.
For multiple coverage values additional aesthetics are not currently supported.
}

\section{Computed variables}{

\code{stat_central_region} computes
\code{after_stat(ymax)} and \code{after_stat(ymin)} for the high and low value of the central region.

For multiple coverages the variables use the same names as \code{\link{central_region}},
i.e. \code{hi.95} and \code{lo.95} for the region with 95\% coverage.
}

\examples{
require("ggplot2")
## Generate some data
#------------------------------------------------------
# Simulate regression data according to the cubic model
# f(x) = 0.8x - 1.8x^2 + 1.05x^3 for x in [0,1]
par <- c(0,0.8,-1.8,1.05) # Parameters of the true polynomial model
res <- 100 # Resolution
x <- seq(0, 1, by=1/res); x2=x^2; x3=x^3;

f <- par[1] + par[2]*x + par[3]*x^2 + par[4]*x^3 # The true function
d <- f + rnorm(length(x), 0, 0.04) # Data

# Estimate polynomial regression model
reg <- lm(d ~ x + x2 + x3)
ftheta <- reg$fitted.values
resid0 <- reg$residuals

# Bootstrap regression
B <- 200 # Number of bootstrap samples
df <- NULL
for(i in 1:B) {
  u <- sample(resid0, size=length(resid0), replace=TRUE)
  reg1 <- lm((ftheta+u) ~ x + x2 + x3)
  df <- rbind(df, data.frame(y=reg1$fitted.values, x=x, i=i,
    g=ifelse(i<14, "A", "B"), g2=ifelse(i<100, "A", "B")))
}

ggplot(df) + geom_line(aes(x, y, group=i))
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=0.50)
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=0.50, filled=FALSE)
# Central regions for two groups as specified by 'g2'
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i, col=g2), coverage=0.90, filled=FALSE)
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=0.90) + facet_wrap(vars(g2))
\dontshow{
# If nr. of functions < 20, then the functions are drawn; otherwise the 100*coverage\% central region
ggplot(df[df$i < 10,]) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=0.90)
# Central regions for two groups split by 'g'; <20 functions in the first group
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i, col=g, fill=g), coverage=0.90)
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i, col=g), coverage=0.90, filled=FALSE)
}

# Central regions with multiple coverage levels
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=c(0.2,0.4,0.6)) +
  theme_minimal()
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=seq(0.1, 0.9, length=20),
  colours=rainbow(20))

\donttest{
# Colors for multiregions are not supported
ggplot(df) + geom_central_region(aes(x=x, y=y+0.1*(g2=="B"),
  curveid=i, col=as.factor(g2)), coverage=c(0.05, 0.2,0.4,0.6))
}
ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i),
  coverage=c(0.05, 0.2,0.4,0.6)) + facet_wrap(vars(g2))

# Using stat_central_region with geom_linerange and geom_rect
ggplot(df) +
  geom_linerange(aes(curveid=i, x=x, y=y, ymax=after_stat(ymax), ymin=after_stat(ymin),
               group=g2, col=factor(g2)),
               stat="central_region", coverage = 0.90, position=position_dodge(0.01))
ggplot(within(df, {x = x+0.004*(g2=="B")})) +
  geom_rect(aes(curveid=i, x=x, y=y, xmax=after_stat(x), xmin=after_stat(x+0.004),
              ymax=after_stat(ymax), ymin=after_stat(ymin), group=g2, fill=factor(g2)),
              stat="central_region", coverage = 0.90)
\donttest{
# Non-finite values are not supported
ggplot(within(df, {y = ifelse(runif(length(y)) < 0.001, Inf, y)})) +
  geom_central_region(aes(x=x, y=y, curveid=i))
}
}
\seealso{
\code{\link{central_region}} for the basic computation and,
\code{\link[ggplot2]{geom_ribbon}} for the default base geom.
}
