# A helper function for extracting information from the curve set
# (envelope, fdata, curve_set) object given for global envelope calculations
pick_attributes <- function(curve_set, alternative = "two.sided") {
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
    res <- list()
    if(inherits(curve_set, 'envelope')) {
        names <- c("fname", "argu", "labl", "ylab", "yexp", "desc")
        for(i in seq_along(names)) res[[names[i]]] <- attr(curve_set, names[i])
        res[['desc']][4] <- lo.name
        res[['desc']][5] <- hi.name
    } else if(inherits(curve_set, 'fdata')) {
      if(!is.null(curve_set$names[['xlab']]))
        res[['xlab']] <- curve_set$names[['xlab']]
      else
        res[['xlab']] <- expression(italic(r))
      if(!is.null(curve_set$names[['ylab']]))
        res[['ylab']] <- curve_set$names[['ylab']]
      else
        res[['ylab']] <- expression(italic(T(r)))
    } else {
        res[['xlab']] <- expression(italic(r))
        res[['ylab']] <- expression(italic(T(r)))
    }
    res[['alternative']] <- alternative
    res
}

# Reset attributes xlab, argu, xecp, ylab, yexp used for plotting purposes
# (Used by special functions utilizing global envelope tests.)
envelope_set_labs <- function(x, xlab, ylab) {
  if(!missing(xlab)) attr(x, "xlab") <- xlab
  if(!missing(ylab)) attr(x, "ylab") <- ylab
  x
}

# An internal GET function for setting the default main for a global envelope plot.
# @param x An 'global_envelope' object.
env_main_default <- function(x, digits=3, alternative=get_alternative(x)) {
  if(!is.null(attr(x, "p_interval"))) {
    if(alternative == "two.sided")
      main <- paste0(attr(x, "method"), ": p-interval = (",
                     round(attr(x, "p_interval")[1], digits=digits),", ",
                     round(attr(x, "p_interval")[2], digits=digits), ")")
    else
      main <- paste0(attr(x, "method"), ": p-interval = (",
                     round(attr(x, "p_interval")[1], digits=digits),", ",
                     round(attr(x, "p_interval")[2], digits=digits), ") \n",
                     "Alternative = \"", alternative, "\"\n")
  }
  else {
    if(!is.null(attr(x, "p", exact=TRUE))) {
      p <- round(attr(x, "p", exact=TRUE), digits=digits)
      if(p > 0) main <- paste(attr(x, "method"), ": p = ", p, sep="")
      else main <- paste0(attr(x, "method"), ": p < ", 10^(-digits))
      if(alternative != "two.sided")
        main <- paste0(main, "\n",
                       "Alternative = \"", alternative, "\"\n")
    }
    else {
      if(inherits(x, c("fdr_envelope", "combined_fdr_envelope"))) {
        main <- paste0(100*(1-attr(x, "alpha")), "% FDR envelope")
      }
      else {
        if(!is.null(attr(x, "alpha"))) {
          if(inherits(x, c("fboxplot", "combined_fboxplot")))
            main <- paste0(attr(x, "method"), " based on ", 100*(1-attr(x, "alpha")), "% central region (", attr(x, "type"), ")")
          else if(inherits(x, c("global_envelope", "combined_global_envelope"))) {
            main <- NULL
            for(i in seq_along(attr(x, "alpha")))
              main <- paste0(main, 100*(1-attr(x, "alpha")[i]), "% ")
            main <- paste0(main, "central region (", attr(x, "type"), ")")
          }
          else
            main <- NULL
        }
        else
          main <- paste0(attr(x, "method"), " (", attr(x, "type"), ")")
      }
    }
  }
  main
}

# Labels for n plots or for a dotplot style envelope plot
default_labels <- function(x, labels) {
  if(inherits(x, "list")) { # Case: n plots
    n <- length(x)
    # Define labels
    if(missing(labels)) {
      if(!is.null(attr(x, "labels")))
        labels <- attr(x, "labels")
      else {
        if(!is.null(names(x)))
          labels <- names(x)
        else {
          labels <- sapply(x, function(y) attr(y, "ylab"), simplify=TRUE)
          if(all(sapply(labels, FUN=identical, y=labels[[1]])))
            labels <- paste0(1:n)
        }
      }
    }
    # Check and edit length
    if(length(labels)!=n) {
      if(length(labels)==1) {
        labels <- paste(labels, " - ", 1:n, sep="")
        warning(paste("Consider giving labels as a vector of length ", n,
                      " containing the label for each test function/vector used.", sep=""))
      }
      else {
        warning("The length of the vector labels is unreasonable. Setting labels to running numbers.")
        labels <- paste0(1:n)
      }
    }
  }
  else { # Case: dotplot
    if(missing(labels)) {
      if(!is.null(attr(x, "labels")))
        labels <- attr(x, "labels")
      else
        labels <- NULL
    }
  }
  labels
}

plotdefaultlabs <- function(x) {
  choice <- function(attrname) {
    if(is.expression(attr(x, attrname)))
      substitute(i, list(i=attr(x, attrname)))
    else
      substitute(italic(i), list(i=attr(x, attrname)))
  }
  if(!is.null(attr(x, "xlab"))) {
    xlab <- choice("xlab")
  } else if(!is.null(attr(x, "argu"))) {
    xlab <- choice("argu")
  } else {
    xlab <- expression(italic(r))
  }
  if(!is.null(attr(x, "yexp"))) {
    ylab <- choice("yexp")
  } else if(!is.null(attr(x, "ylab"))) {
    ylab <- choice("ylab")
  } else {
    ylab <- expression(italic(T(r)))
  }
  list(xlab=xlab, ylab=ylab)
}

# An inner function for a 'dotplot' style envelope plot with ggplot2.
#' @importFrom ggplot2 ggplot geom_segment aes .data geom_point scale_color_identity scale_x_discrete
#' @importFrom grid arrow
env_dotplot_ggplot <- function(x, labels=NULL, sign.col="red") {
  if(is.null(labels) && !is.null(x[['r']])) labels <- paste(round(x[['r']], digits=2))
  df <- as.data.frame(x)
  df$r <- factor(df$r)
  arrow <- arrow(angle=75)
  if(length(attr(x, "alpha")) > 1) message("Note: dotplot shows only the largest envelope.")
  loname <- env_loname(attr(x, "alpha"), largest=TRUE)
  hiname <- env_hiname(attr(x, "alpha"), largest=TRUE)
  g <- ggplot(df) + geom_segment(aes(x=.data$r, y=.data$central, xend=.data$r, yend=.data[[hiname]]), arrow=arrow) +
    geom_segment(aes(x=.data$r, y=.data$central, xend=.data$r, yend=.data[[loname]]), arrow=arrow)
  if(!is.null(x[['obs']])) {
    if(is.null(sign.col)) sign.col <- "black"
    g <- g + geom_point(aes(x=.data$r, y=.data$obs, col=ifelse(.data$obs > .data[[hiname]] | .data$obs < .data[[loname]], sign.col, "black")), shape="x", size=5)
  }
  g <- g + geom_point(aes(x=.data$r, y=.data$central)) +
    scale_color_identity() +
    scale_x_discrete(breaks=paste(x[['r']]), labels=labels)
  g
}

# Global envelope plots
#----------------------

# Construct a data.frame for the envelope (gg)plot
env_df_construction <- function(x, main) {
  n <- names(x)
  df <- data.frame(r = x[['r']],
                   curves = x[['central']],
                   type = factor("Central function", levels = "Central function"))
  if('obs' %in% n) {
    df_obs <- data.frame(r = x[['r']],
                         curves = x[['obs']],
                         type = factor("Data function", levels = "Data function"))
    df <- rbind(df, df_obs)
  }
  for(w in setdiff(n, c("r", "obs", "central"))) {
    df[, w] <- x[[w]]
  }
  df$plotmain <- main # Needed for combined plots
  df
}

linetype_values <- function() { c('dashed', 'solid') }

# Basic elements of the envelope (gg)plot
# level = The significance level(s) of the data that are found in df.
# If 0 (or other single number), then "lo" and "hi" should be found in df.
# Otherwise several lo.xx and hi.xx, where xx represent different levels.
#' @importFrom ggplot2 geom_ribbon aes_string geom_line
#' @importFrom ggplot2 labs scale_linetype_manual
basic_stuff_for_env_ggplot <- function(df, xlab, ylab, main, level=0) {
  pE <- vector("list", length(level))
  lonames <- env_loname(level)
  hinames <- env_hiname(level)
  if(length(level) ==  1) cols <- 'grey59'
  else cols <- paste0('grey', floor(seq(80, 59, length=length(level))))
  for(i in seq_along(level)) {
    pE[[i]] <- geom_ribbon(data = df, aes_string(x = "r", ymin = lonames[i], ymax = hinames[i]),
                     fill = cols[i], alpha = 1)
  }
  c(pE, list(
       geom_line(data = df, aes_(x = ~r, y = ~curves, group = ~type,
                                 linetype = ~type)), # , size = 0.2
       labs(title = main, x = xlab, y = ylab),
       scale_linetype_manual(values = linetype_values(), name = '')))
}
basic_stuff_for_fclustplot <- function(df, xlab, ylab, main, fillcolor = 'grey59', alpha = 0.5, size=0.3) {
  list(geom_ribbon(data = df, aes_(x = ~r, ymin = ~lo, ymax = ~hi),
                   fill = fillcolor, alpha = alpha),
       geom_line(data = df, aes_(x = ~r, y = ~curves, group = ~type,
                                 linetype = ~type), size = size),
       labs(title = main, x = xlab, y = ylab),
       scale_linetype_manual(values = 'solid', name = ''))
}

# An internal function for making a ggplot2 style "global envelope plot"
# @param x An 'global_envelope' object or a list of them.
# @param main See \code{\link{plot.default}}.
# @param ylim See \code{\link{plot.default}}.
# @param xlab See \code{\link{plot.default}}.
# @param ylab See \code{\link{plot.default}}.
# @param color_outside Logical, whether to use sign.col.
# @param sign.col Color for the observed curve outside the envelope.
#' @importFrom ggplot2 ggplot theme guides geom_point aes_
env_ggplot <- function(x, main, xlab, ylab, sign.col="red") {
  if(!inherits(x, "global_envelope")) stop("Internal error.")
  df <- env_df_construction(x, NULL)
  p <- ( ggplot()
         + basic_stuff_for_env_ggplot(df, xlab, ylab, main, attr(x, "alpha"))
         + set_envelope_legend_position() )
  if(length(levels(df$type)) < 2) p <- p + guides(linetype = "none")
  if("Data function" %in% levels(df$type)) {
    if(!is.null(sign.col)) {
      df.outside <- df[df$type == "Data function",
                       c("r", "curves",
                         env_loname(attr(x, "alpha"), largest=TRUE),
                         env_hiname(attr(x, "alpha"), largest=TRUE))]
      names(df.outside)[3:4] <- c("lo", "hi")
      df.outside <- df.outside[df.outside$curves < df.outside$lo | df.outside$curves > df.outside$hi,]
      p <- p + geom_point(data=df.outside, ggplot2::aes_(x = ~r, y = ~curves), color=sign.col, size=1)
    }
  }
  p
}

# An internal function for making a ggplot2 style "functional boxplot"
#' @importFrom viridisLite viridis
#' @importFrom ggplot2 ggplot geom_ribbon aes_ guides geom_line scale_color_identity
fboxplot_ggplot <- function(x, main, xlab, ylab, plot_outliers = TRUE) {
    if(!inherits(x, "fboxplot")) stop("x should have class fboxplot. Possibly internal error.")

    # Basic df
    df <- env_df_construction(x, NULL)
    p <- ( ggplot2::ggplot()
           + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~whisker.lo, ymax = ~whisker.hi),
                                  fill = 'grey80', alpha = 1)
           + basic_stuff_for_env_ggplot(df, xlab, ylab, main)
           + guides(linetype = "none") )

    if(plot_outliers & !is.null(attr(x, "outliers"))) {
      out <- attr(x, "outliers")
      col <- viridis(ncol(out))
      out.df <- data.frame(r = rep(x[['r']], times=ncol(out)),
                           curves = c(out),
                           id = rep(colnames(out), each=length(x[['r']])),
                           col = rep(col, each=length(x[['r']])))
      p <- ( p + geom_line(data = out.df, ggplot2::aes_(x = ~r, y = ~curves, group = ~id, col=~col))
               + scale_color_identity("", labels = colnames(out), guide = "legend") )
    }
    p
}


# Combined envelope plots
#------------------------

combined_df_construction <- function(x, labels) {
  n <- names(x[[1]])

  dfs <- vector("list", length(x))
  for(i in seq_along(x)) {
    dfs[[i]] <- env_df_construction(x[[i]], labels[i])
  }
  df <- do.call(rbind, dfs)
  df$plotmain <- factor(df$plotmain, levels = labels)
  df
}

# An internal function for making a ggplot2 style "combined global envelope plot"
# @param labels Labels for components of the combined tests.
# @param max_ncols_of_plots The maximum number of columns for figures. Default 2.
#' @importFrom ggplot2 ggplot facet_wrap guides geom_point aes_ theme
env_combined_ggplot <- function(x, main, xlab, ylab, labels, scales = "free",
                       max_ncols_of_plots = 2, sign.col="red") {
  if(!inherits(x, "list")) stop("Internal error. x is not a list.")
  Nfunc <- length(x)

  n_of_plots <- as.integer(Nfunc)
  ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
  nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)

  df <- combined_df_construction(x, labels=labels)
  p <- ( ggplot()
         + basic_stuff_for_env_ggplot(df, xlab, ylab, main, attr(x, "alpha"))
         + facet_wrap(~ plotmain, scales=scales,
                      nrow=nrows_of_plots, ncol=ncols_of_plots)
         + set_envelope_legend_position() )
  if(length(levels(df$type)) < 2)
    p <- p + guides(linetype = "none")
  if("Data function" %in% levels(df$type)) {
    if(!is.null(sign.col)) {
      df.outside <- df[df$type == "Data function",
                       c("r", "curves",
                         env_loname(attr(x, "alpha"), largest=TRUE),
                         env_hiname(attr(x, "alpha"), largest=TRUE),
                         "plotmain")]
      names(df.outside)[3:4] <- c("lo", "hi")
      df.outside <- df.outside[df.outside$curves < df.outside$lo | df.outside$curves > df.outside$hi,]
      p <- p + geom_point(data=df.outside, ggplot2::aes_(x=~r, y=~curves), color=sign.col, size=1)
    }
  }
  p
}
# An internal function for making a dotplot style "combined global envelope plot"
#' @importFrom ggplot2 ggplot geom_segment scale_color_identity facet_wrap geom_point scale_x_discrete
#' @importFrom grid arrow
env_combined_dotplot <- function(x, main, xlab, ylab, labels, scales = "free",
                                max_ncols_of_plots = 2, sign.col="red", ...) {
  if(!inherits(x, "list")) stop("Internal error. x is not a list.")

  n_of_plots <- as.integer(length(x))
  ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
  nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)

  df <- NULL
  plot_labels <- default_labels(x)
  for(i in 1:length(x))
    df <- rbind(df, data.frame(x[[i]], plotmain=plot_labels[i]))

  arrow <- arrow(angle=75, ...)
  if(length(attr(x, "alpha")) > 1) message("Note: dotplot shows only the largest envelope.")
  loname <- env_loname(attr(x, "alpha"), largest=TRUE)
  hiname <- env_hiname(attr(x, "alpha"), largest=TRUE)
  df$r <- factor(df$r)

  g <- ( ggplot(df)
         + geom_segment(aes(x=.data$r, y=.data$central, xend=.data$r, yend=.data[[hiname]]), arrow=arrow)
         + geom_segment(aes(x=.data$r, y=.data$central, xend=.data$r, yend=.data[[loname]]), arrow=arrow)
         + set_envelope_legend_position()
         + scale_color_identity()
         + facet_wrap(~ plotmain, scales=scales,
                      nrow=nrows_of_plots, ncol=ncols_of_plots) )
  if(!is.null(x[[1]][['obs']])) {
    if(is.null(sign.col)) sign.col <- "black"
    g <- g + geom_point(aes(x=.data$r, y=.data$obs, col=ifelse(.data$obs > .data[[hiname]] | .data$obs < .data[[loname]], sign.col, "black")), shape="x", size=5)
  }
  g <- g + geom_point(aes(x=.data$r, y=.data$central))
  if(!missing(labels) && all(sapply(x, FUN = function(y) { length(y$r) == length(x[[1]]$r) }, simplify=TRUE)))
    g <- g + scale_x_discrete(breaks=paste(x[[1]][['r']]), labels=labels)
  g
}

# An internal function for making a ggplot2 style "combined functional boxplot"
#' @importFrom viridisLite viridis
#' @importFrom ggplot2 ggplot geom_ribbon aes_ facet_wrap guides geom_line scale_color_identity
fboxplot_combined_ggplot <- function(x, main, xlab, ylab, labels, scales = "free",
                                max_ncols_of_plots = 2, plot_outliers = TRUE) {
  if(!inherits(x, "list")) stop("Internal error. x is not a list.")
  Nfunc <- length(x)

  n_of_plots <- as.integer(Nfunc)
  ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
  nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)

  df <- combined_df_construction(x, labels=labels)
  p <- ( ggplot2::ggplot()
         + ggplot2::geom_ribbon(data = df,
                                ggplot2::aes_(x = ~r, ymin = ~whisker.lo, ymax = ~whisker.hi),
                                fill = 'grey80', alpha = 1)
         + basic_stuff_for_env_ggplot(df, xlab, ylab, main)
         + facet_wrap(~ plotmain, scales=scales,
                      nrow=nrows_of_plots, ncol=ncols_of_plots)
         + guides(linetype = "none") )
  if(plot_outliers & !is.null(attr(x, "outliers"))) {
    out <- attr(x, "outliers")
    col_values <- viridis(ncol(out[[1]]))
    names(col_values) <- colnames(out[[1]])
    out.df <- do.call(rbind, lapply(seq_along(out), FUN = function(i) {
      data.frame(r = rep(x[[i]][['r']], times=ncol(out[[i]])),
                 curves = c(out[[i]]),
                 id = factor(rep(colnames(out[[i]]), each=length(x[[i]][['r']])), levels = colnames(out[[i]])),
                 col = rep(col_values, each=length(x[[i]][['r']])),
                 plotmain = factor(rep(labels[i], each=length(x[[i]][['r']])), levels=labels))
    }))
    p <- ( p + ggplot2::geom_line(data = out.df, ggplot2::aes_(x = ~r, y = ~curves, group = ~id, col = ~col))
             + scale_color_identity("Outliers", labels = colnames(out[[1]]), guide = "legend") )
  }
  p
}
