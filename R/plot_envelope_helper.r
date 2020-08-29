# A helper function for avoiding repetition of code (in envelope tests).
pick_attributes <- function(curve_set, alternative, type) {
    # saving for attributes / plotting purposes
    lo.name <- "lower critical boundary for %s"
    hi.name <- "upper critical boundary for %s"
    switch(alternative,
            two.sided = {},
            less = {
                hi.name <- "infinite upper boundary"
            },
            greater = {
                lo.name <- "infinite lower boundary"
            })
    if(inherits(curve_set, 'envelope')) {
        names <- c("argu", "fname", "labl", "ylab", "yexp", "einfo")
        for(i in 1:length(names)) assign(names[i], attr(curve_set, names[i]))
        desc <- attr(curve_set, "desc")
        desc[4] <- lo.name
        desc[5] <- hi.name
        einfo$global <- TRUE
        einfo$alternative <- alternative
        einfo$scale <- NULL
        einfo$clamp <- NULL
        einfo$nrank <- NULL
        einfo$VARIANCE <- (type == "st")
        einfo$nSD <- NULL
    }
    else {
        fname <- "T"
        argu <- "r"
        if(curve_set_is1obs(curve_set)) {
          labl <- c("r", "T[obs](r)", "T[0](r)", "T[lo](r)", "T[hi](r)")
          desc <- c("distance argument r",
                    "observed value of %s for data pattern",
                    "central curve under the null hypothesis",
                    lo.name, hi.name)
        }
        else {
          labl <- c("r", "T[0](r)", "T[lo](r)", "T[hi](r)")
          desc <- c("distance argument r",
                    "central curve under the null hypothesis",
                    lo.name, hi.name)
        }
        ylab <- "T(r)"
        yexp <- quote(T(r))
        # tack on envelope information
        einfo <- list(alternative = alternative)
    }
    list(argu=argu, fname=fname, labl=labl, desc=desc, ylab=ylab, yexp=yexp,
         xlab=argu, xexp=quote(r), einfo=einfo)
}

# A helper function to check whether the xaxis needs to be reticked with new values due to
# combined tests. Called also from plot.global_envelope for checking if fv plot style is available.
retick_xaxis <- function(x) {
  if(!inherits(x, "list")) x <- list(x)
  if(any(sapply(x, FUN=function(x) { !(inherits(x, c("global_envelope", "fboxplot", "curve_set"))) })))
    stop("x should consist of global_envelope objects.")
  r_values_ls <- lapply(x, FUN=function(x) x$r)
  r_values <- do.call(c, r_values_ls, quote=FALSE)
  nr <- length(r_values)
  list(retick_xaxis = !(length(x) == 1 & all(r_values[-1] - r_values[-nr] > 0)),
       r_values_ls = r_values_ls, r_values = r_values)
}

# Define breaking r values and labels on x-axis for plotting several
# global envelopes or curve sets jointly.
# @param x A global_envelope object or a list of global_envelope objects.
# Also curve_set objects allowed with a restricted use.
# @param nticks Number of ticks per a sub test.
combined_global_envelope_rhelper <- function(x, nticks = 5) {
  if(!inherits(x, "list")) x <- list(x)
  retick <- retick_xaxis(x)
  r_values_ls <- retick$r_values_ls
  r_values <- retick$r_values
  nr <- length(r_values)
  if(!retick$retick_xaxis) {
      new_r_values <- NULL
      r_break_values <- NULL
      loc_break_values <- NULL
      r_values_newstart_id <- NULL
  }
  else {
    if(length(x) > 1 & any(unlist(lapply(r_values_ls, function(x) { !all(x[-1] - x[-length(x)] > 0) })))) {
      warning(paste("Something strange. The r values are not increasing in a", class(x[[1]])[1], "object.\n", sep=""))
    }
    new_r_values <- 1:nr # to be used in plotting
    # Define where the functions start when they are put all together
    if(length(x) == 1) { # Find from the r-values
      r_values_newstart_id <- which(!(r_values[1:(nr-1)] < r_values[2:nr])) + 1
    }
    else { # Define directly from the r_values_ls
      r_values_newstart_id <- NULL
      r_values_newstart_id[1] <- length(r_values_ls[[1]]) + 1
      if(length(r_values_ls) > 2) {
        for(i in 2:(length(r_values_ls)-1))
          r_values_newstart_id <- c(r_values_newstart_id, r_values_newstart_id[i-1] + length(r_values_ls[[1]]))
      }
    }
    # r-values for labeling ticks
    r_starts <- r_values[c(1, r_values_newstart_id)]
    r_ends <- r_values[c(r_values_newstart_id - 1, nr)]
    r_break_values <- NULL
    # indeces for ticks in the running numbering from 1 to nr
    loc_starts <- (1:nr)[c(1, r_values_newstart_id)]
    loc_ends <- (1:nr)[c(r_values_newstart_id - 1, nr)]
    loc_break_values <- NULL
    nslots <- length(r_starts) # number of combined tests/slots
    for(i in 1:(nslots-1)) {
      r_break_values <- c(r_break_values, seq(r_starts[i], r_ends[i], length=nticks)[1:(nticks-1)])
      loc_break_values <- c(loc_break_values, seq(loc_starts[i], loc_ends[i], length=nticks)[1:(nticks-1)])
    }
    r_break_values <- c(r_break_values, seq(r_starts[nslots], r_ends[nslots], length=nticks))
    loc_break_values <- c(loc_break_values, seq(loc_starts[nslots], loc_ends[nslots], length=nticks))
  }
  if(inherits(x[[1]], c("global_envelope", "fboxplot"))) {
    if(!is.null(x[[1]]$obs))
      x_vec <- data.frame(r = r_values,
                          obs = do.call(c, lapply(x, FUN = function(x) x$obs), quote=FALSE),
                          central = do.call(c, lapply(x, FUN = function(x) x$central), quote=FALSE),
                          lo = do.call(c, lapply(x, FUN = function(x) x$lo), quote=FALSE),
                          hi = do.call(c, lapply(x, FUN = function(x) x$hi), quote=FALSE))
    else
      x_vec <- data.frame(r = r_values,
                          central = do.call(c, lapply(x, FUN = function(x) x$central), quote=FALSE),
                          lo = do.call(c, lapply(x, FUN = function(x) x$lo), quote=FALSE),
                          hi = do.call(c, lapply(x, FUN = function(x) x$hi), quote=FALSE))
  }
  else x_vec <- NULL

  list(x_vec = x_vec,
       retick_xaxis = retick$retick_xaxis,
       new_r_values = new_r_values,
       r_break_values = r_break_values, loc_break_values = loc_break_values,
       r_values_newstart_id = r_values_newstart_id)
}

# An internal GET function for setting the default main for a global envelope plot.
# @param x An 'global_envelope' object.
env_main_default <- function(x, digits=3, alternative=get_alternative(x)) {
  if(!is.null(attr(x, "p_interval"))) {
    if(alternative == "two.sided")
      main <- paste(attr(x, "method"), ": p-interval = (",
                    round(attr(x, "p_interval")[1], digits=digits),", ",
                    round(attr(x, "p_interval")[2], digits=digits), ")", sep="")
    else
      main <- paste(attr(x, "method"), ": p-interval = (",
                    round(attr(x, "p_interval")[1], digits=digits),", ",
                    round(attr(x, "p_interval")[2], digits=digits), ") \n",
                    "Alternative = \"", alternative, "\"\n", sep="")
  }
  else {
    if(!is.null(attr(x, "p"))) {
      p <- round(attr(x, "p"), digits=digits)
      if(p > 0) main <- paste(attr(x, "method"), ": p = ", p, sep="")
      else main <- paste(attr(x, "method"), ": p < ", 10^(-digits), sep="")
      if(alternative != "two.sided")
        main <- paste(main, "\n",
                      "Alternative = \"", alternative, "\"\n", sep="")
    }
    else {
      if(inherits(x, c("fboxplot", "combined_fboxplot")))
        main <- paste(attr(x, "method"), " based on ", 100*(1-attr(x, "alpha")), "% central region (", attr(x, "type"), ")", sep="")
      else if(inherits(x, c("global_envelope")))
         main <- paste(100*(1-attr(x, "alpha")), "% central region (", attr(x, "type"), ")", sep="")
      else main <- NULL
    }
  }
  main
}

# An internal GET function for setting the default ylim for a global envelope plot.
# @param x An 'global_envelope' object or a list of them.
# @param use_ggplot2 TRUE/FALSE, If TRUE, then default ylim are for \code{\link{env_ggplot}}.
# Otherwise the ylim are for \code{\link{env_basic_plot}}.
env_ylim_default <- function(x, use_ggplot2) {
  if(!inherits(x, "list")) x <- list(x)
  if(!use_ggplot2)
    ylim <- lapply(x, {
      function(y) {
        switch(get_alternative(y),
               two.sided = {
                 ylim <- c(min(y[['obs']],y[['lo']],y[['hi']],y[['central']]),
                           max(y[['obs']],y[['lo']],y[['hi']],y[['central']]))
               },
               less = {
                 ylim <- c(min(y[['obs']],y[['lo']],y[['central']]),
                           max(y[['obs']],y[['lo']],y[['central']]))
               },
               greater = {
                 ylim <- c(min(y[['obs']],y[['hi']],y[['central']]),
                           max(y[['obs']],y[['hi']],y[['central']]))
               })
      }
    })
  else ylim <- NULL
  if(length(ylim)==1) ylim <- unlist(ylim)
  ylim
}

# An internal GET function for making a dotplot style "global envelope plot".
#
# An internal GET function for making a dotplot style "global envelope plot".
#
# @param x An 'global_envelope' object.
# @param main See \code{\link{plot.default}}.
# @param ylim See \code{\link{plot.default}}.
# @param xlab See \code{\link{plot.default}}.
# @param ylab See \code{\link{plot.default}}.
# @param color_outside Logical. Whether to color the places where the data function goes
# outside the envelope. Currently red color is used.
# @param labels Labels for the tests at x-axis.
# @param add Whether to add the plot to an existing plot (TRUE) or to draw a new plot (FALSE).
# @param arrows.col Color for the doplot arrows. If not given, 1 (black) is used.
# @param curve_sets If provided, then curves going outside the envelope are plotted.
# @param ... Additional parameters to be passed to the function \code{\link{plot}}.
#' @importFrom graphics plot
#' @importFrom graphics arrows
#' @importFrom graphics points
#' @importFrom graphics axis
env_dotplot <- function(x, main, ylim, xlab, ylab, color_outside = TRUE,
                        labels = NULL, add = FALSE, arrows.col, curve_sets = NULL, ...) {
    nr <- length(x[['r']])
    if(is.null(labels)) labels <- paste(round(x[['r']], digits=2))
    if(missing(arrows.col)) arrows.col <- 1
    if(nr > 10) warning("Dotplot style meant for low dimensional test vectors.")
    if(!add) graphics::plot(1:nr, x[['central']], main=main, ylim=ylim, xlab=xlab, ylab=ylab, cex=0.5, pch=16, xaxt="n", ...)
    else graphics::points(1:nr, x[['central']], main=main, ylim=ylim, xlab=xlab, ylab=ylab, cex=0.5, pch=16, xaxt="n", ...)
    if(get_alternative(x)!="greater")
        graphics::arrows(1:nr, x[['lo']], 1:nr, x[['central']], code = 1, angle = 75, length = .1, col=arrows.col)
    else
        graphics::arrows(1:nr, x[['lo']], 1:nr, x[['central']], code = 1, angle = 75, length = .1, col='grey80')
    if(get_alternative(x)!="less")
        graphics::arrows(1:nr, x[['hi']], 1:nr, x[['central']], code = 1, angle = 75, length = .1, col=arrows.col)
    else
        graphics::arrows(1:nr, x[['hi']], 1:nr, x[['central']], code = 1, angle = 75, length = .1, col='grey80')
    graphics::axis(1, 1:nr, labels=labels)
    if(!is.null(x[['obs']])) {
      graphics::points(1:nr, x[['obs']], pch='x')
      if(color_outside) {
        outside <- x[['obs']] < x[['lo']] | x[['obs']] > x[['hi']]
        graphics::points((1:nr)[outside], x[['obs']][outside], pch='x', col="red")
      }
    }
    if(!is.null(curve_sets)) {
      funcs <- curve_set_funcs(curve_sets)
      for(i in 1:ncol(funcs)) {
        if(any(funcs[,i] < x[['lo']] | funcs[,i] > x[['hi']])) {
          graphics::points(1:nr, funcs[,i], pch='x', col='grey70', type="b")
        }
      }
    }
}

# Turn an envelope into a 'dotplot' with ggplot2.
#' @importFrom ggplot2 arrow ggplot geom_segment aes .data geom_point scale_color_identity scale_x_discrete
env_dotplot_ggplot <- function(x, labels=NULL) {
  if(is.null(labels) && !is.null(x[['r']])) labels <- paste(round(x[['r']], digits=2))
  df <- as.data.frame(x)
  arrow <- arrow(angle=75)
  g <- ggplot(df) + geom_segment(aes(x=.data$r, y=.data$central, xend=.data$r, yend=.data$hi), arrow=arrow) +
    geom_segment(aes(x=factor(.data$r), y=.data$central, xend=.data$r, yend=.data$lo), arrow=arrow)
  if(!is.null(x[['obs']])) g <- g + geom_point(aes(x=factor(.data$r), y=.data$obs, col=ifelse(.data$obs > .data$hi | .data$obs < .data$lo, "red", "black")), shape="x", size=5)
  g <- g + geom_point(aes(x=factor(.data$r), y=.data$central)) +
    scale_color_identity() +
    scale_x_discrete(breaks=paste(x[['r']]), labels=labels)
}

# An internal GET function for making a basic "global envelope plot".
#
# An internal GET function for making a basic "global envelope plot".
#
# @inheritParams env_dotplot
# @param x An 'global_envelope' object or a list of them.
# @param max_ncols_of_plots The maximum number of columns for figures. Default 2.
# @param env.col The color for the envelope lines. Default 1 (black).
# @param nticks The number of ticks on the xaxis, if the xaxis is re-ticked for combined tests.
# @param obs.type The line type of the observed function (if exists).
# @param ... Additional parameters to be passed to the function \code{\link{plot}}.
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics lines
#' @importFrom graphics axis
#' @importFrom graphics abline
env_basic_plot <- function(x, main, ylim, xlab, ylab, color_outside=TRUE,
                           max_ncols_of_plots = 2, add = FALSE, env.col = 1,
                           nticks = 5, curve_sets = NULL, obs.type = "l", ...) {
    if(!inherits(x, "list")) x <- list(x)
    Nfunc <- length(x)
    # Handle combined tests; correct labels on x-axis
    # a) if x is a list of global_envelope objects
    # b) if x[['r']] contains repeated values (when length(x) == 1)
    rdata <- combined_global_envelope_rhelper(x, nticks=nticks)
    alt <- get_alternative(x[[1]])
    x <- rdata$x_vec
    # Plot
    if(Nfunc == 1 & is.null(rdata$r_values_newstart_id)) {
        if(!rdata$retick_xaxis) {
            if(!add) plot(x[['r']], x[['central']], main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                                    type="l", lty=3, lwd=2, ...)
            else lines(x[['r']], x[['central']], lty=3, lwd=2, ...)
        }
        else {
          x[['r']] <- 1:length(x[['r']])
          if(!add) plot(x[['r']], x[['central']], main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                                  type="l", lty=3, lwd=2, xaxt="n", ...)
          else lines(x[['r']], x[['central']], lty=3, lwd=2, xaxt="n", ...)
          if(!is.null(curve_sets) && class(curve_sets)[1] == "list") curve_sets <- combine_curve_sets(curve_sets, equalr=FALSE)
        }
        if(alt != "greater") lines(x[['r']], x[['lo']], lty=2, col=env.col) else lines(x[['r']], x[['lo']], lty=2, col='grey80')
        if(alt != "less") lines(x[['r']], x[['hi']], lty=2, col=env.col) else lines(x[['r']], x[['hi']], lty=2, col='grey80')
        if(!is.null(x[['obs']])) {
          lines(x[['r']], x[['obs']], lty=1, type=obs.type)
          if(color_outside) {
              outside <- x[['obs']] < x[['lo']] | x[['obs']] > x[['hi']]
              graphics::points(x[['r']][outside], x[['obs']][outside], col="red")
          }
        }
        # curves/outliers
        if(!is.null(curve_sets)) {
          funcs <- curve_set_funcs(curve_sets)
          for(i in 1:ncol(funcs)) {
            if(any(funcs[,i] < x[['lo']] | funcs[,i] > x[['hi']]))
              lines(x[['r']], funcs[,i], col='grey70')
          }
        }
        if(rdata$retick_xaxis) {
            axis(1, rdata$loc_break_values, labels=paste(round(rdata$r_break_values, digits=2)))
            abline(v = rdata$new_r_values[rdata$r_values_newstart_id], lty=3)
        }
    }
    else {
        if(Nfunc == 1) warning("The r-values are non-increasing in the given object. Splitting to several plots.")
        n_of_plots <- as.integer(1 + length(rdata$r_values_newstart_id))
        ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
        nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)
        opar <- par(mfrow=c(nrows_of_plots, ncols_of_plots))
        on.exit(par(opar))
        tmp_indeces <- c(1, rdata$r_values_newstart_id, length(rdata$new_r_values)+1)
        if(is.null(ylim) | !(inherits(ylim, "list") && length(ylim) == n_of_plots)) {
          if(is.vector(ylim) & length(ylim)==2) ylim <- rep(list(ylim), times=n_of_plots)
          else {
            ylim <- list()
            for(i in 1:n_of_plots)
              ylim[[i]] <- c(min(x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                                 x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                                 x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)]),
                             max(x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                                 x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                                 x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)]))
          }
        }
        if(length(main) != n_of_plots) { main <- NULL; message("Note: \"main\" Ignored.") }
        if(!inherits(xlab, "list") && length(xlab) == 1) xlab <- rep(list(xlab), times=n_of_plots)
        if(!inherits(ylab, "list") && length(ylab) == 1) ylab <- rep(list(ylab), times=n_of_plots)
        if(!is.null(curve_sets)) {
          curve_sets <- combine_curve_sets(curve_sets, equalr=FALSE)
        }
        for(i in 1:n_of_plots) {
            if(!add)
              plot(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             x[['central']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             main=main[i], xlab=xlab[[i]], ylab=ylab[[i]],
                             type="l", lty=3, lwd=2, ylim=ylim[[i]], ...)
            else
              lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             x[['central']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             lty=3, lwd=2, ...)
            if(alt != "greater")
              lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2)
            else
              lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2, col='grey80')
            if(alt != "less")
              lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2)
            else
              lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2, col='grey80')
            if(!is.null(x[['obs']])) {
              lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=1, type=obs.type)
              if(color_outside) {
                  outside <- x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)] < x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)] | x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)] > x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)]
                  points(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)][outside], x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)][outside], col="red")
              }
            }
            # curves/outliers
            if(!is.null(curve_sets)) {
              funcs <- curve_set_funcs(curve_sets)
              for(j in 1:ncol(funcs)) {
                if(any(funcs[,j] < x[['lo']] | funcs[,j] > x[['hi']]))
                  lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                        funcs[tmp_indeces[i]:(tmp_indeces[i+1]-1),j], col='grey70')
              }
            }
        }
    }
}


# An internal GET function for making a ggplot2 style "global envelope plot".
#
# An internal GET function for making a ggplot2 style "global envelope plot".
#
# @param x An 'global_envelope' object or a list of them.
# @param base_size Base font size, to be passed to theme style when \code{plot_type="ggplot2"}.
# @param main See \code{\link{plot.default}}.
# @param ylim See \code{\link{plot.default}}.
# @param xlab See \code{\link{plot.default}}.
# @param ylab See \code{\link{plot.default}}.
# @param labels Labels for the separate plots.
# @param max_ncols_of_plots The maximum number of columns for figures. Default 2.
# @param labels Labels for components of the combined tests.
# @param nticks The number of ticks on the xaxis, if the xaxis is re-ticked for combined tests.
# @param curve_sets If provided, then curves going outside the envelope are plotted.
# @param x2 Another 'global_envelope' object, which is plotted within x, i.e. x2 is assumed to be narrower
# of the two envelopes.
# @param legend Logical. If FALSE, then legend is removed.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 aes_
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_linetype_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom utils tail
env_ggplot <- function(x, base_size, main, ylim, xlab, ylab,
                       max_ncols_of_plots = 2,
                       labels = NULL, nticks = 5, curve_sets = NULL, x2 = NULL,
                       legend = TRUE, color_outside=TRUE, sign.col="red") {
    if(!inherits(x, "list")) x <- list(x)
    Nfunc <- length(x)
    if(!is.null(x2)) {
      if(!inherits(x2, "list")) x2 <- list(x2)
      if(length(x) != length(x2)) {
        warning("Unsuitable x2. Setting it to NULL.")
        x2 <- NULL
      }
      else {
        for(i in 1:length(x)) {
          if(!all(x[[i]][['r']] == x2[[i]][['r']])) stop("The two envelopes are for different r-values.")
          if(!all(x[[i]][['central']] == x2[[i]][['central']])) warning("The two envelopes have different central functions!")
        }
      }
      rdata <- combined_global_envelope_rhelper(x2, nticks=nticks)
      x2 <- rdata$x_vec
    }
    # Handle combined tests; correct labels on x-axis
    # a) if x is a list of global_envelope objects
    # b) if x[['r']] contains repeated values (when length(x) == 1)
    rdata <- combined_global_envelope_rhelper(x, nticks=nticks)
    alt <- get_alternative(x[[1]])
    x <- rdata$x_vec

    linetype.values <- c('dashed', 'solid')
    size.values <- c(0.2, 0.2)

    outliers <- outliers_id <- NULL
    if(!is.null(curve_sets)) {
      if(inherits(curve_sets, "list")) curve_sets <- combine_curve_sets(curve_sets, equalr=FALSE)
      funcs <- curve_set_funcs(curve_sets)
      for(j in 1:ncol(funcs)) {
        if(any(funcs[,j] < x[['lo']] | funcs[,j] > x[['hi']])) {
          outliers <- c(outliers, funcs[,j])
          outliers_id <- c(outliers_id, j)
        }
      }
    }

    if(Nfunc == 1 & is.null(rdata$r_values_newstart_id)) {
      if(rdata$retick_xaxis) x[['r']] <- 1:length(x[['r']])
      if(is.null(x[['obs']])) {
        df <- data.frame(r = x[['r']],
                         curves = x[['central']],
                         type = factor("Central function", levels = "Central function"),
                         lower = x[['lo']],
                         upper = x[['hi']],
                         main = main)
        if(!is.null(x2)) {
          df$lower2 <- x2[['lo']]
          df$upper2 <- x2[['hi']]
        }
      }
      else {
        df <- data.frame(r = rep(x[['r']], times=2),
                         curves = c(x[['obs']], x[['central']]),
                         type = factor(rep(c("Data function", "Central function"),
                                           each=length(x[['r']])),
                                       levels=c("Central function", "Data function")),
                         lower = rep(x[['lo']], times=2),
                         upper = rep(x[['hi']], times=2),
                         main = main)
        if(!is.null(x2)) {
          df$lower2 <- rep(x2[['lo']], times=2)
          df$upper2 <- rep(x2[['hi']], times=2)
        }
      }
      if(is.null(x2)) {
        p <- ( ggplot2::ggplot()
               + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower, ymax = ~upper),
                                      fill = 'grey59', alpha = 1)
        )
      }
      else {
        p <- ( ggplot2::ggplot()
               + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower, ymax = ~upper),
                                      fill = 'grey80', alpha = 1)
               + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower2, ymax = ~upper2),
                                      fill = 'grey59', alpha = 1)
        )
      }
      p <- ( p + ggplot2::geom_line(data = df, ggplot2::aes_(x = ~r, y = ~curves, group = ~type,
                                                             linetype = ~type, size = ~type))
             + ggplot2::facet_grid('~ main', scales = 'free')
             + ggplot2::scale_y_continuous(name = ylab, limits = ylim)
             + ggplot2::scale_linetype_manual(values = linetype.values, name = '')
             + ggplot2::scale_size_manual(values = size.values, name = '')
             + ThemePlain(base_size=base_size)
      )
      if(is.null(x[['obs']])) p <- p + ggplot2::guides(linetype = "none", size = "none")
      if(!is.null(outliers)) {
        outliers.df <- data.frame(r = rep(x[['r']], times=length(outliers_id)),
                                  curves = outliers,
                                  id = rep(outliers_id, each=length(x[['r']])))
        last_points <- outliers.df[outliers.df$r==max(outliers.df$r),]
        extrasp <- (max(outliers.df$r)-min(outliers.df$r))/50
        p <- p + ggplot2::geom_line(data = outliers.df, ggplot2::aes_(x = ~r, y = ~curves, group = ~id)) +
          ggplot2::geom_text(data = last_points, ggplot2::aes_(x=~r+extrasp, y=~curves, label = ~id))
      }
      if(rdata$retick_xaxis) {
        p <- p + ggplot2::scale_x_continuous(name = xlab,
                                             breaks = rdata$loc_break_values,
                                             labels = paste(round(rdata$r_break_values, digits=2)),
                                             limits = range(rdata$new_r_values))
        p <- p + ggplot2::geom_vline(xintercept = rdata$new_r_values[rdata$r_values_newstart_id],
                                     linetype = "dotted")
      }
      else p <- p + ggplot2::scale_x_continuous(name = xlab)
    }
    else {
      if(Nfunc == 1) warning("The r-values are non-increasing in the given object. Splitting to several plots.")
      n_of_plots <- as.integer(1 + length(rdata$r_values_newstart_id))
      ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
      nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)
      if(is.null(labels)) labels <- paste(1:n_of_plots)
      if(length(labels)!=n_of_plots) {
        if(length(labels)==1) {
          labels <- paste(labels, " - ", 1:n_of_plots, sep="")
          warning(paste("Consider giving labels as a vector of length ", n_of_plots,
                        " containing the label for each test function/vector used.", sep=""))
        }
        else {
          warning("The length of the vector labels is unreasonable. Setting labels to empty.")
          labels <- rep("", times=n_of_plots)
        }
      }

      tmp_indeces <- c(1, rdata$r_values_newstart_id, length(rdata$new_r_values)+1)
      func_labels <- NULL
      for(i in 1:(length(tmp_indeces)-1)) {
        func_labels <- c(func_labels, rep(labels[i], times=tmp_indeces[i+1]-tmp_indeces[i]))
      }

      if(is.null(x[['obs']])) {
        df <- data.frame(r = x[['r']],
                         curves =x[['central']],
                         type = factor("Central function", levels="Central function"),
                         lower = x[['lo']],
                         upper = x[['hi']],
                         main = main,
                         test_function = factor(func_labels, levels=labels))
        if(!is.null(x2)) {
          df$lower2 <- x2[['lo']]
          df$upper2 <- x2[['hi']]
        }
      }
      else {
        df <- data.frame(r = rep(x[['r']], times=2),
                         curves = c(x[['obs']], x[['central']]),
                         type = factor(rep(c("Data function", "Central function"),
                                           each=length(x[['r']])),
                                       levels=c("Central function", "Data function")),
                         lower = rep(x[['lo']], times=2),
                         upper = rep(x[['hi']], times=2),
                         main = main,
                         test_function = factor(func_labels, levels=labels))
        if(!is.null(x2)) {
          df$lower2 <- rep(x2[['lo']], times=2)
          df$upper2 <- rep(x2[['hi']], times=2)
        }
      }
      if(is.null(x2)) {
        p <- ( ggplot2::ggplot()
               + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower, ymax = ~upper),
                                      fill = 'grey59', alpha = 1)
        )
      }
      else {
        p <- ( ggplot2::ggplot()
               + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower, ymax = ~upper),
                                      fill = 'grey80', alpha = 1)
               + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower2, ymax = ~upper2),
                                      fill = 'grey59', alpha = 1)
        )
      }
      p <- (p + ggplot2::geom_line(data = df, ggplot2::aes_(x = ~r, y = ~curves, group = ~type,
                                                            linetype = ~type, size = ~type))
            + ggplot2::facet_wrap(~ test_function, scales="free",
                                  nrow=nrows_of_plots, ncol=ncols_of_plots)
            + ggplot2::scale_y_continuous(name = ylab, limits = ylim)
            + ggplot2::scale_linetype_manual(values = linetype.values, name = '')
            + ggplot2::scale_size_manual(values = size.values, name = '')
            + ThemePlain(base_size=base_size)
            + ggplot2::labs(title=main)
      )
      if(is.null(x[['obs']])) p <- p + ggplot2::guides(linetype = "none", size = "none")
      if(!is.null(outliers)) {
        outliers.df <- data.frame(r = rep(x[['r']], times=length(outliers_id)),
                                  curves = outliers,
                                  id = rep(outliers_id, each=length(x[['r']])),
                                  test_function = factor(func_labels, levels=labels))
        last_points <- outliers.df[outliers.df$r==max(outliers.df$r),]
        extrasp <- (max(outliers.df$r)-min(outliers.df$r))/50
        p <- p + ggplot2::geom_line(data = outliers.df, ggplot2::aes_(x = ~r, y = ~curves, group = ~id)) +
          ggplot2::geom_text(data = last_points, ggplot2::aes_(x=~r+extrasp, y=~curves, label = ~id))
      }
      p <- p + ggplot2::scale_x_continuous(name = xlab)
    }
    if(!legend) p <- p + theme(legend.position = "none")
    if(!is.null(x[['obs']])) {
      if(color_outside) {
        df.outside <- df[df$type == "Data function",]
        df.outside <- df.outside[df.outside$curves < df.outside$lower | df.outside$curves > df.outside$upper,]
        p <- p + geom_point(data=df.outside, ggplot2::aes_(x = ~r, y = ~curves), color=sign.col, size=1)
      }
    }
    # Return
    p
}
