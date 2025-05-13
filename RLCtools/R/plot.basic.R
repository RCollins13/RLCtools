#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Basic plotting functions to generate entire plots
# See plot.helpers.R for smaller subroutines used to generate certain plot elements


#' Two-dimensional scatterplot
#'
#' Generate a scatterplot of two numeric values with options for coloring
#'
#' @param X Numeric vector of values to be plotted on the X-axis (see `Details`)
#' @param Y Numeric vector of values to be plotted on the Y-axis (see `Details`)
#' @param colors Vector of colors. Must be supplied in the same order as `X` and `Y` \[default: "gray70"\]
#' @param title Main title for plot \[default: NULL\]
#' @param x.label.line Line for X-axis labels (`label.line` parameter for [RLCtools::clean.axis])
#' @param x.title Title for X-axis
#' @param x.title.line Line for X-axis title (`title.line` parameter for [RLCtools::clean.axis])
#' @param xlims Limits for X-axis
#' @param y.label.line Line for X-axis labels (`label.line` parameter for [RLCtools::clean.axis])
#' @param y.title Title for Y-axis
#' @param y.title.line Line for axis titles (`title.line` parameter for [RLCtools::clean.axis])
#' @param ylims Limits for Y-axis
#' @param legend.vals Named vector mapping category names to colors \[default: NULL\]
#' @param legend.labels Optional vector to overwrite names of `legend.vals`
#' @param cex Character expansion factor for individual points \[default: 0.3\]
#' @param add Should this scatterplot be added on top of an existing graphics
#' device? \[default: FALSE\]
#' @param parmar Numeric vector of values to pass to `par(mar)`
#'
#' @details If `legend.vals` is not provided (or is `NULL`), no legend will be added.
#'
#' @export scatterplot
#' @export
scatterplot <- function(X, Y, colors=NULL, title=NULL,
                        x.label.line=NULL, x.title=NULL, x.title.line=0.5, xlims=NULL,
                        y.label.line=NULL, y.title=NULL, y.title.line=0.5, ylims=NULL,
                        legend.vals=NULL, legend.labels=NULL,
                        cex=0.3, add=FALSE, parmar=c(2.5, 2.5, 1, 1)){
  x <- as.numeric(X)
  y <- as.numeric(Y)
  if(is.null(xlims)){
    xlims <- range(x, na.rm=T)
  }
  if(is.null(ylims)){
    ylims <- range(y, na.rm=T)
  }
  if(is.null(colors)){
    colors <- rep("gray70", length(x))
  }

  # Prepare plot area
  if(!add){
    prep.plot.area(xlims, ylims, parmar=parmar, xaxs="r", yaxs="r")
  }
  mtext(3, text=title)

  # Add X-axis
  if(is.null(x.title)){
    x.title <- "X value"
  }
  clean.axis(1, title=x.title, infinite=T, label.line=x.label.line, title.line=x.title.line)

  # Add Y-axis
  if(is.null(y.title)){
    y.title <- "Y value"
  }
  clean.axis(2, title=y.title, infinite=T, label.line=y.label.line, title.line=y.title.line)

  # Add points
  points(x, y, pch=19, cex=cex, col=colors, xpd=T)

  # Add legend, if optioned
  if(!is.null(legend.vals)){
    if(!is.null(legend.labels)){
      names(legend.vals) <- legend.labels
    }
    quad.counts <- table(x > mean(par("usr")[1:2]), y < mean(par("usr")[3:4]))
    least.dense.quad <- head(which(quad.counts == min(quad.counts)), 1)
    legend.pos <- if(least.dense.quad == 1){
      "topleft"
    }else if(least.dense.quad == 2){
      "topright"
    }else if(least.dense.quad == 3){
      "bottomleft"
    }else if(least.dense.quad == 4){
      "bottomright"
    }
    legend(legend.pos, legend=names(legend.vals), pch=21, pt.bg=legend.vals,
           pt.cex=1.5, bty="n", border=NA, cex=5/6)
  }
}


#' Ridgeplot
#'
#' Generate a ridgeplot using base R syntax
#'
#' @param data List of numeric vectors or [density()] objects to plot.
#' See `Details`.
#' @param bw.adj Numeric vector of `adjust` values passed to [density()]
#' if `data` is provided as numeric vectors rather than pre-computed
#' [density()] objects. See `Details`.
#' @param breaks Numeric vector of histogram bin breaks. If specified,
#' this will override the value of `bw.adj` and cause the plot to be rendered
#' as a histogram. See `Details`.
#' @param names Optional list of names for Y axis \[default: take names from data\]
#' @param hill.overlap Relative fraction of overlap bewtween adjacent
#' hills \[default: 0.35\]
#' @param x.axis.side `side` value for x-axis; `NA` will disable X-axis
#' plotting. \[default: 1\]
#' @param x.title Title for X axis \[default: "Values"\]
#' @param x.title.line Value of `title.line` passed to [RLCtools::clean.axis()]
#' @param x.label.line Value of `label.line` passed to [RLCtools::clean.axis()]
#' @param max.x.ticks Value of `max.ticks` passed to [RLCtools::clean.axis()]
#' \[default: 6\]
#' @param x.tick.len Value of `tck` passed to [RLCtools::clean.axis()]
#' \[default: -0.025\]
#' @param xlims Custom X axis limits
#' @param y.axis Should groups be labeled on the y-axis? \[default: TRUE\]
#' @param ylims Custom Y axis limits
#' @param yaxs Value of `yaxs` passed to [plot()] \[default: "r"\]
#' @param fill Vector of polygon fill colors \[default: "grey70"\]
#' @param border Vector of hill border colors \[default: "grey35"\]
#' @param border.lwd Line width for hill borders \[default: 2\]
#' @param hill.bottom Relative value (0-1) indicating the vertical starting
#' position of each hill should be \[default: 0\]
#' @param fancy.hills Should hills be shaded & marked like a boxplot? Can only
#' be computed if `data` is provided as raw values, not pre-computed densities.
#' \[default: TRUE\]
#' @param fancy.light.fill Vector of colors to use for the sections of the
#' distributions beyond the IQR. Only used if `fancy.hills` is `TRUE`
#' \[default: "grey85"\]
#' @param fancy.quartile.color Vector of colors to be used for quartile
#' indicators if `fancy.hills` is `TRUE` \[default: "white"\]
#' @param fancy.quartile.lwd Value of `lwd` to use for fancy median and
#' IQR lines \[default: 2\]
#' @param fancy.quartile.lend Value of `lend` to use for fancy median line
#' \[default: "square"\]
#' @param parmar Vector of values passed to par(mar)
#'
#' @details `data` can be provided either as a list of numeric vectors (one per
#' hill to be plotted) or as a list of pre-computed [density()] objects.
#'
#' The type of `data` will be automatically checked prior to plotting and will
#' be converted to [density()] if needed.
#'
#' Optionally, `bw.adj` can be supplied alongside a numeric vector-style `data`
#' to customize the bandwidth for each hill. This vector will be recycled in the
#' usual way following R conventions.
#'
#' Lastly, if `breaks` are specified and if `data` is provided as a list of
#' numeric vectors, no [density()] will be computed and all hills will be
#' rendered as histograms instead.
#'
#' @seealso [density()], [hist()]
#'
#' @export ridgeplot
#' @export
ridgeplot <- function(data, bw.adj=NULL, breaks=NULL, names=NULL, hill.overlap=0.35,
                      x.axis.side=1, x.title="Values", x.title.line=0.3,
                      x.label.line=-0.65, max.x.ticks=6, x.tick.len=-0.025,
                      xlims=NULL, y.axis=TRUE, ylims=NULL, yaxs="r", fill=NULL,
                      border=NULL, border.lwd=2, hill.bottom=0, fancy.hills=TRUE,
                      fancy.light.fill=NULL, fancy.quartile.color=NULL,
                      fancy.quartile.lwd=2, fancy.quartile.lend="square",
                      parmar=c(2.5, 3, 0.25, 0.25)){
  # Get names before manipulating data
  if(is.null(names)){
    names <- names(data)
    if(is.null(names)){
      names <- 1:length(data)
    }
  }

  # Collect distribution statistics for fancy hills, if optioned & possible
  if(fancy.hills & is.numeric(data[[1]])){
    meds <- sapply(data, median, na.rm=T)
    q1s <- sapply(data, quantile, prob=0.25, na.rm=T)
    q3s <- sapply(data, quantile, prob=0.75, na.rm=T)
  }else{
    fancy.hills <- FALSE
  }

  # Coerce data to KDE or hist, if needed
  if(is.numeric(data[[1]])){
    if(!is.null(breaks)){
      as.hist <- TRUE
      data <- lapply(1:length(data), function(i){
        if(length(data[[i]]) > 1){
          data[[i]][which(data[[i]] < min(breaks))] <- min(breaks)
          data[[i]][which(data[[i]] > max(breaks))] <- max(breaks)
          hist(data[[i]], breaks=breaks, plot=F)$density
        }else{
          NULL
        }
      })
      names(data) <- names
    }else{
      as.hist <- FALSE
      if(is.null(bw.adj)){
        bw.adj <- rep(1, length(data))
      }
      if(length(bw.adj) < length(data)){
        bw.adj <- rep(bw.adj, length(data) / length(bw.adj))
      }
      data <- lapply(1:length(data), function(i){
        if(length(data[[i]]) > 1){
          density(data[[i]], adjust=bw.adj[i])
        }else{
          NULL
        }
      })
      names(data) <- names
    }
  }

  # Scale Y values of data to [hill.bottom, hill.overlap]
  for(i in 1:length(data)){
    if(!is.null(data[[i]])){
      y <- if(as.hist){data[[i]]}else{data[[i]]$y}
      y <- y - min(y)
      if(as.hist){
        data[[i]] <- ((1 + hill.overlap - hill.bottom) * (y / max(y))) + hill.bottom
      }else{
        data[[i]]$y <- ((1 + hill.overlap - hill.bottom) * (y / max(y))) + hill.bottom
      }
    }
  }

  # Get plot dimensions
  if(is.null(xlims)){
    if(as.hist){
      xlims <- range(breaks)
    }else{
      xlims <- c(min(sapply(data, function(d){min(d$x, na.rm=T)})),
                 max(sapply(data, function(d){max(d$x, na.rm=T)})))
    }
  }
  if(is.null(ylims)){
    ylims <- c(0, length(data) + max(c(0, hill.overlap)))
  }

  # Get ridge colors
  if(is.null(fill)){
    fill <- rep("grey70", length(data))
  }
  if(length(fill) < length(data)){
    fill <- rep(fill, ceiling(length(data) / length(fill)))
  }
  if(is.null(border)){
    border <- rep("grey35", length(data))
  }
  if(length(border) < length(data)){
    border <- rep(border, length(data) / length(border))
  }
  if(fancy.hills){
    if(is.null(fancy.light.fill)){
      fancy.light.fill <- rep("grey85", length(data))
    }
    if(length(fancy.light.fill) < length(data)){
      fancy.light.fill <- rep(fancy.light.fill,
                              length(data) / length(fancy.light.fill))
    }
    if(is.null(fancy.quartile.color)){
      fancy.quartile.color <- rep("white", length(data))
    }
    if(length(fancy.quartile.color) < length(data)){
      fancy.quartile.color <- rep(fancy.quartile.color,
                                  length(data) / length(fancy.quartile.color))
    }
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="i", yaxs=yaxs)
  if(!is.na(x.axis.side)){
    clean.axis(x.axis.side, title=x.title, infinite=TRUE, max.ticks=max.x.ticks,
               label.line=x.label.line, title.line=x.title.line, tck=x.tick.len)
  }
  if(y.axis){
    axis(2, at=(1:length(data)) - 0.5, tick=F, las=2, line=-0.8, labels=names)
  }

  # Add hills
  sapply(length(data):1, function(i){
    abline(h=i-1+hill.bottom, col="gray85")
    if(is.null(data[[i]])){
      text(x=mean(par("usr")[1:2]), y=mean(c(i-1+hill.bottom, i)),
           cex=5/6, labels="No data", font=3, col="gray80")
    }else{
      if(as.hist){
        steps <- step.function(breaks[-length(breaks)], data[[i]], offset=1)
        x <- c(steps$x, rev(steps$x))
        y <- c(steps$y, rep(hill.bottom, times=length(steps$y)))+i-1
      }else{
        x <- c(data[[i]]$x, rev(data[[i]]$x))
        y <- c(data[[i]]$y, rep(hill.bottom, times=length(data[[i]]$y)))+i-1
      }
      if(fancy.hills){
        if(as.hist){
          left.break.idx <- which(breaks < q1s[i])
          left.step <- step.function(as.numeric(c(breaks[left.break.idx], q1s[i])),
                                     data[[i]][c(left.break.idx,
                                                 max(left.break.idx))],
                                     offset=1)
          left.x <- c(left.step$x, rev(left.step$x))
          left.y <- c(left.step$y, rep(hill.bottom,
                                       times=length(left.step$y)))+i-1
          mid.break.idx <- which(breaks >= q1s[i] & breaks < q3s[i])
          mid.step <- step.function(as.numeric(c(q1s[i], breaks[mid.break.idx])),
                                    data[[i]][c(min(mid.break.idx)-1,
                                                mid.break.idx)],
                                    offset=1)
          mid.step$x[length(mid.step$x)] <- q3s[i]
          mid.x <- c(mid.step$x, rev(mid.step$x))
          mid.y <- c(mid.step$y, rep(hill.bottom,
                                     times=length(mid.step$y)))+i-1
          right.break.idx <- which(breaks >= q3s[i])
          right.break.idx <- right.break.idx[-length(right.break.idx)]
          right.step <- step.function(as.numeric(c(q3s[i],
                                                   breaks[right.break.idx])),
                                      data[[i]][c(min(right.break.idx)-1,
                                                  right.break.idx)], offset=1)
          right.x <- c(right.step$x, rev(right.step$x))
          right.y <- c(right.step$y, rep(hill.bottom,
                                         times=length(right.step$y)))+i-1
        }else{
          left.idx <- which(x < q1s[i])
          left.x <- x[left.idx]
          left.y <- y[left.idx]
          mid.idx <- which(x >= q1s[i] & x <= q3s[i])
          mid.x <- x[mid.idx]
          mid.y <- y[mid.idx]
          right.idx <- which(x > q3s[i])
          right.x <- x[right.idx]
          right.y <- y[right.idx]
        }
        polygon(x, y, border=NA, col="white", xpd=T)
        polygon(left.x, left.y, border=fancy.light.fill[i],
                col=fancy.light.fill[i], xpd=T, lwd=0.5)
        polygon(right.x, right.y, border=fancy.light.fill[i],
                col=fancy.light.fill[i], xpd=T, lwd=0.5)
        polygon(mid.x, mid.y, border=fill[i], col=fill[i], xpd=T, lwd=0.5)
        segments(x0=meds[i], x1=meds[i], y0=i-1+hill.bottom,
                 y1=if(as.hist){
                   data[[i]][max(which(meds[i] > breaks))]+i-1
                 }else{
                   y[which.min(abs(meds[i] - x))]
                 },
                 lwd=fancy.quartile.lwd, col=fancy.quartile.color[i], xpd=T,
                 lend=fancy.quartile.lend)
        segments(x0=c(q1s[i], q3s[i]), x1=c(q1s[i], q3s[i]),
                 y0=i-1+hill.bottom,
                 y1=if(as.hist){
                   data[[i]][c(max(which(q1s[i] > breaks)),
                               max(which(q3s[i] > breaks)))]+i-1
                 }else{
                   y[c(which.min(abs(q1s[i] - x)), which.min(abs(q3s[i] - x)))]
                 },
                 lwd=fancy.quartile.lwd, col=fancy.quartile.color[i], xpd=T,
                 lend=fancy.quartile.lend)
      }else{
        polygon(x, y, border=fill[i], col=fill[i], lwd=border.lwd, xpd=T)
      }
      points(if(as.hist){steps$x}else{data[[i]]$x},
             if(as.hist){steps$y+i-1}else{data[[i]]$y+i-1},
             type="l", lwd=border.lwd, col=border[i], xpd=T)
    }
  })
}


#' Quantile-quantile plot
#'
#' Generate a Q-Q plot of association stats vs. uniform null
#'
#' @param pvals Numeric vector of untransformed P-values
#' @param cutoff P-value threshold for significance
#' @param do.fdr Should points meeting `fdr.cutoff` be accented?
#' \[default: TRUE\]
#' @param fdr.cutoff Cutoff for highlighting FDR-significant points
#' \[default: 0.01\]
#' @param print.stats Should genomic inflation statistic be printed on
#' plot? \[default: FALSE\]
#' @param title (Optional) Plot title
#' @param title.line Line for `title` \[default: 0\]
#' @param xmax Maximum X-value to plot \[default: include all points\]
#' @param ymax Maximum Y-value to plot \[default: smallest P-value\]
#' @param cap.pvals Should P-values more significant than `ymax` be capped
#' at `ymax`? \[default: FALSE\]
#' @param x.title Text or expression for x axis title
#' \[default: `Expected -log10 P`\]
#' @param y.title Text or expression for x axis title
#' \[default: `Observed -log10 P`\]
#' @param x.label.line Value of `label.line` passed to
#' [`RLCtools::clean.axis()`]  \[default: -0.75]
#' @param x.title.line Value of `title.line` passed to
#' [`RLCtools::clean.axis()`] \[default: 0.35]
#' @param y.label.line Value of `label.line` passed to
#' [`RLCtools::clean.axis()`]  \[default: -0.65]
#' @param y.title.line Value of `title.line` passed to
#' [`RLCtools::clean.axis()`] \[default: 0.15]
#' @param label.cex Scaling factor for label text \[default: 5/6\]
#' @param ax.title.cex Scaling factor for axis title text \[default: 1\]
#' @param title.cex Scaling factor for label text \[default: 1\]
#' @param axis.tck Value of `tck` passed to [`RLCtools::clean.axis()`]
#' \[default: -0.025].
#' @param pt.color Color for all points; see `Details` \[default: "grey35"\]
#' @param pt.cex Scaling factor for all points \[default: 0.35\]
#' @param fdr.color Color for FDR-significant points \[default: "grey5"\]
#' @param fdr.cex Scaling factor for FDR-significant points \[default: 0.7\]
#' @param plot.ci Should a shaded confidence interval be added to the plot?
#' \[default: TRUE\]
#' @param ci.color Color for shaded confidence interval \[default: "gray90"\]
#' @param oe.line.color Color for observed ~ expected line \[default: "gray50"\]
#' @param plot Should a QQ plot be generated? \[default: TRUE\]
#' @param add Should points be added to an existing graphics device?
#' \[default: generate new plot\]
#' @param parmar Value of `mar` passed to `par()`
#' @param return.xy Should \(x,y\) coordinates of QQ points be returned?
#' \[default: FALSE\]
#'
#' @details
#' `cutoff` should be specified in untransformed P-value units; i.e., not -log10
#'
#' `ymax` should be specified in transformed P-value units; i.e.,
#' -log10\(min desired P\)
#'
#' `pt.color` can be specified as a single color or a vector of colors
#' for pointwise color assignment. If the length of `pt.color` does not match
#' the length of `pvals`, the values in the `pt.color` vector will be recycled
#' according to R conventions
#'
#' @returns Dependent on the value of `return.xy`:
#' - When `TRUE`, returns a data.frame
#' - When `FALSE` \(default\), returns NULL
#'
#' @export plot.qq
#' @export
plot.qq <- function(pvals, cutoff=NULL, do.fdr=TRUE, fdr.cutoff=0.01,
                    print.stats=FALSE, title=NULL, title.line=0, xmax=NULL,
                    ymax=NULL, cap.pvals=FALSE, x.title=NULL, y.title=NULL,
                    x.label.line=-0.75, y.label.line=-0.65, x.title.line=0.35,
                    y.title.line=0.15, label.cex=5/6, ax.title.cex=1,
                    axis.tck=-0.025, pt.color="grey35", pt.cex=0.35,
                    fdr.color="grey5", fdr.cex=0.7, plot.ci=TRUE,
                    ci.color="gray90", oe.line.color="gray50", plot=TRUE,
                    add=FALSE, parmar=c(2.25, 2.5, 0.25, 0.25),
                    return.xy=FALSE){
  # Format P-values and pointwise colors
  if (!is.numeric(pvals)){
    stop("P values must be numeric.")
  }
  if(length(pt.color) <= length(pvals)){
    colors <- rep(pt.color, ceiling(length(pvals) / length(pt.color)))
  }
  colors <- colors[1:length(pvals)]
  keep.idxs <- which(!is.na(pvals) & !is.nan(pvals) & !is.null(pvals) &
                       is.finite(pvals) & pvals <= 1 & pvals >= 0)
  pvals <- pvals[keep.idxs]
  colors <- colors[keep.idxs]
  pw.cex <- rep(pt.cex, length(pvals))
  p.order <- order(pvals)
  pvals <- pvals[p.order]
  colors <- colors[p.order]
  if(do.fdr){
    fdr.idx <- which(p.adjust(pvals, "fdr") < fdr.cutoff)
    if(length(fdr.idx) > 0){
      colors[fdr.idx] <- fdr.color
      pw.cex[fdr.idx] <- fdr.cex
    }
  }

  # Compute expected quantiles
  expected <- ppoints(length(pvals))
  qqconf <- function (p.expected, reflection=F){
    n <- length(p.expected)
    mpts <- matrix(nrow=n * 2, ncol=2)
    for (i in seq(from=1, to=n)) {
      mpts[i, 1] <- -log10((i - 0.5)/n)
      mpts[i, 2] <- -log10(qbeta(0.975, i, n - i))
      mpts[n * 2 + 1 - i, 1] <- -log10((i - 0.5)/n)
      mpts[n * 2 + 1 - i, 2] <- -log10(qbeta(0.025, i, n - i))
    }
    mpts <- as.data.frame(mpts)
    return(mpts)
  }
  conf.int <- qqconf(expected, reflection)
  lambda <- dchisq(median(pvals), df=1)/dchisq(median(expected), df=1)

  # Get other plot parameters
  if(is.null(cutoff)){
    cutoff <- 0.05/length(pvals)
  }
  if(is.null(ymax)){
    maxp <- max(-log10(pvals[which(!(is.infinite(-log10(pvals))))]))
    ymax <- max(c(maxp, -log10(cutoff) + 2))
  }
  expected <- -log10(expected)
  if(is.null(xmax)){
    xmax <- max(expected, na.rm=T)
  }
  pvals <- -log10(pvals)
  log.cutoff <- -log10(cutoff)

  # Floor large P-values, if optioned
  any.pvals.rounded <- FALSE
  if(cap.pvals){
    big.pvals <- which(pvals > ymax)
    if(length(big.pvals) > 0){
      pvals[big.pvals] <- ymax
      any.pvals.rounded <- TRUE
    }
  }

  # Generate plot if optioned
  if(plot){
    # Prep plot area and add null confidence interval
    if(!add){
      prep.plot.area(c(0, 1.1 * xmax), c(0, 1.1 * ymax), parmar=parmar)
      if(plot.ci){
        polygon(x=conf.int[, 1], y=conf.int[, 2], col=ci.color, border=NA)
        abline(0, 1, col=oe.line.color)
        if(!is.null(cutoff)){
          abline(h=log.cutoff, lty=2)
        }
      }
    }

    # Annotate genomic control if optioned
    if (print.stats == T){
      text(par("usr")[1] - (0.025 * (par("usr")[2] - par("usr")[1])),
           0.9 * par("usr")[4], pos=4, font=2, cex=0.8,
           labels=bquote(lambda ~ .(paste("=", sprintf("%.2f", lambda), sep=""))))
    }

    # Add points
    points(x=expected, y=pvals, pch=19, col=colors, cex=pw.cex)

    # Add axes & title
    if(!add){
      if(is.null(x.title)){
        x.title <- expression(Expected ~ ~-log[10] ~ italic(P))
      }
      clean.axis(1, title=x.title,
                 title.line=x.title.line, cex.title=ax.title.cex,
                 label.line=x.label.line, cex.axis=label.cex,
                 infinite.positive=TRUE, tck=axis.tck)
      y.labs <- y.at <- axTicks(2)
      if(any.pvals.rounded){
        y.at <- unique(c(y.at[-length(y.at)], ymax))
        y.labs <- c(y.at[-length(y.at)],
                    paste("\"\" > ", y.at[length(y.at)], sep=""))
      }
      if(is.null(y.title)){
        y.title <- expression(Observed ~ ~-log[10] ~ italic(P))
      }
      clean.axis(2, at=y.at, labels=y.labs, parse.labels=TRUE, title=y.title,
                 title.line=y.title.line, cex.title=ax.title.cex,
                 label.line=y.label.line, cex.axis=label.cex,
                 infinite.positive=TRUE, tck=axis.tck)
      mtext(3, line=title.line, text=title)
    }
  }

  # Return point coordinates, if optioned
  if(return.xy){
    return(data.frame("x"=expected, "y"=pvals))
  }
}


#' Manhattan plot
#'
#' Plot genome-wide summary statistics in the "Manhattan" style
#'
#' @param stats data.frame with at least two columns: "pos" and `p.column`
#' @param p.columns Name of column containing P values
#' @param chrom.colors Named vector of point colors for each chromosome
#' @param gw.sig P-value threshold for genome-wide significance
#' @param pt.cex Value of `cex` passed to [points()] \[default: 0.3\]
#' @param title (Optional) title
#' @param x.title Title for x axis \[default: "Genomic coordinate"\]
#' @param parmar Value of `mar` passed to [par()]
#'
#' @export manhattan
#' @export
manhattan <- function(stats, p.column, chrom.colors, gw.sig, pt.cex=0.3,
                      title=NULL, x.title="Genomic coordinate",
                      parmar=c(2, 2.25, 0.25, 0.3)){
  # Prep plot area
  xlims <- range(stats$pos)
  ylims <- c(0, ceiling(max(c(stats[, p.column], gw.sig), na.rm=T) + 1))
  prep.plot.area(xlims, ylims, parmar=parmar)
  abline(h=gw.sig, lty=5)

  # Add points
  set.seed(2024)
  stats <- stats[sample(1:nrow(stats), size=nrow(stats)), ]
  points(stats[, c("pos", p.column)], pch=19,
         cex=pt.cex, col=chrom.colors[stats$chrom])

  # Add axes
  clean.axis(1, at=c(0, cumsum(contig.lengths)), tck=0.025, infinite=TRUE,
             labels=NA, title=x.title, title.line=0)
  sapply(1:length(contig.lengths), function(x){
    axis(1,
         at=((c(0, cumsum(contig.lengths[-24])) + cumsum(contig.lengths))/2)[x],
         tick=F, cex.axis=4/6,
         labels=gsub("^chr", "", names(contig.lengths)[x]),
         col.axis=chrom.colors[x],
         line=c(-0.8, -1.3)[as.integer((x %% 2) == 1) + 1])
  })
  clean.axis(2, title=bquote(-log[10] ~ italic(P)),
             infinite.positive=T, title.line=0.2)
  mtext(3, text=title, line=0)
}


#' Scale-Aware Barplot Cluster
#'
#' Plot a set of stacked barplots scaled proportional to set size
#'
#' @param values List of character vectors of values to plot
#' @param colors Named vector of colors to map to `values`
#' @param group.names (Optional) group names to assign to each list element
#' in `values`
#' @param sep.wex Relative width scalar for whitespace on the X-axis
#' between groups \[default: 0.05\]
#' @param title (Optional) title to be printed in top-right corner
#' @param legend Should a legend be plotted?
#' @param legend.names (Optional) mapping of `values` to labels for legend
#' @param parmar Margin values passed to par()
#'
#' @details If `values` is supplied as a named list, those names will be used as
#' `group.names` unless `group.names` is explicitly specified. Otherwise,
#' `group.names` will be set to the ordinal value of each group in `values`.
#'
#' @seealso [RLCtools::scaled.swarm], [RLCtools::clean.axis]
#'
#' @export scaled.bars
#' @export
scaled.bars <- function(values, colors, group.names=NULL, sep.wex=0.05,
                        title=NULL, legend=TRUE, legend.names=NULL,
                        parmar=c(1, 2, 1.25, 4.5)){
  # Summarize plotting data
  values <- lapply(values, function(v){v[!is.na(v)]})
  elig.vals <- unique(unlist(values))
  if(!all(elig.vals %in% names(colors))){
    stop("Names of 'colors' must cover all values of 'values'")
  }
  colors <- colors[which(names(colors) %in% elig.vals)]
  legend.names <- legend.names[which(names(legend.names) %in% elig.vals)]
  if(is.null(group.names)){
    group.names <- names(values)
  }
  counts <- lapply(values,
                   function(v){sapply(names(colors),
                                      function(k){length(which(v == k))})})
  counts.df <- do.call("cbind", counts)
  n.elig.vals <- length(colors)
  n.groups <- length(values)
  group.size <- apply(counts.df, 2, sum)
  pct.df <- do.call("cbind", lapply(1:n.groups, function(i){
    as.numeric(counts.df[, i] / group.size[i])
    }))
  rownames(pct.df) <- names(colors)
  colnames(pct.df) <- group.names
  cumpct.df <- apply(pct.df, 2, cumsum)

  # Get plot dimensions
  group.widths <- group.size / sum(group.size)
  group.lefts <- c(0, cumsum(group.widths)[-n.groups]) + (sep.wex * (0:(n.groups-1)))
  group.rights <- group.widths + group.lefts
  group.mids <- (group.lefts + group.rights) / 2
  xlims <- c(-sep.wex, max(group.rights))

  # Prep plot area
  prep.plot.area(xlims, c(1, 0), parmar, xaxs="i", yaxs="i")

  # Add X axis
  x.axis.at <- smart.spacing(group.mids, min.dist=0.15)
  axis(1, at=x.axis.at, line=-1, tick=F, labels=group.names, xpd=T)

  # Add Y axis
  clean.axis(2, at=seq(0, 1, 0.25),
             labels=rev(paste(seq(0, 100, 25), "%", sep="")),
             cex.axis=5/6, infinite=FALSE, label.line=-0.7, title=NULL)

  # Add top scale bar
  scale.k <- floor(log10(sum(group.size)))
  bar.len <- (10^scale.k) / sum(group.size)
  if(bar.len > 0.5){
    scale.k <- scale.k - 1
    bar.len <- (10^scale.k) / sum(group.size)
  }
  if(bar.len < 0.25){
    bar.stretch <- ceiling(0.25 / bar.len)
    bar.len <- bar.len * bar.stretch
  }else{
    bar.stretch <- 1
  }
  segments(x0=mean(xlims) - (0.5 * bar.len),
           x1=mean(xlims) + (0.5 * bar.len),
           y0=-0.03, y1=-0.03, xpd=T, lwd=4, lend="butt")
  text(x=mean(xlims), y=-0.01, pos=3, xpd=T, cex=0.85,
       labels=prettyNum(bar.stretch * (10^scale.k), big.mark=","))

  # Add title
  text(x=xlims[2] - sep.wex, y=-0.05, pos=4, font=2, labels=title, xpd=T)

  # Add legend, if optioned
  if(legend){
    if(is.null(legend.names)){
      legend.names <- names(colors)
    }
    leg.mids <- (c(0, cumpct.df[-n.elig.vals, n.groups]) + cumpct.df[, n.groups])/2
    yaxis.legend(legend.names, xlims[2], leg.mids, sep.wex,
                 min.label.spacing=0.075, lower.limit=0.025, upper.limit=0.975,
                 colors)
  }

  # Add rectangles
  sapply(1:n.groups, function(i){
    rect(xleft=group.lefts[i], xright=group.rights[i],
         ybottom=c(0, cumpct.df[-n.elig.vals, i]), ytop=cumpct.df[, i],
         col=colors, border=NA, bty="n")
    rect(xleft=group.lefts[i], xright=group.rights[i],
         ybottom=0, ytop=1, col=NA, xpd=T)
  })
}


#' Scale-Aware Beeswarm Cluster
#'
#' Plot a set of beeswarm distributions scaled proportional to set size
#'
#' @param values List of numeric vectors of values to plot
#' @param colors Vector of colors for the list elements in `values`
#' @param group.names (Optional) group names to assign to each list element
#' in `values`
#' @param group.widths (Optional) vector of custom relative width assignments
#' for all groups \[default: scale widths proportional to group size\]
#' @param sep.wex Relative width scalar for whitespace on the X-axis
#' between groups \[default: 0.05\]
#' @param pch Value of `pch` passed to `beeswarm`
#' @param pt.cex Value of `cex` passed to `beeswarm`
#' @param title Custom title for top axis
#' @param title.line Line for main title on top axis \[default: 0\]
#' @param title.cex `cex` parameter for main title \[default: 1\]
#' @param add.y.axis Should the Y \(left\) axis be added? \[default: TRUE\]
#' @param y.title Title of Y-axis
#' @param y.title.line Value of `line` for `y.title`
#' @param y.axis.at Custom Y-axis tick positions, if desired
#' @param y.axis.labels Custom Y-axis tick labels, if desired
#' @param y.axis.tck Value of `tck` passed to [RLCtools::clean.axis()]
#' \[default: -0.025\]
#' @param parmar Margin values passed to par()
#'
#' @details If `values` is supplied as a named list, those names will be used as
#' `group.names` unless `group.names` is explicitly specified. Otherwise,
#' `group.names` will be set to the ordinal value of each group in `values`.
#'
#' @seealso [RLCtools::scaled.bars], [RLCtools::clean.axis]
#'
#' @export scaled.swarm
#' @export
scaled.swarm <- function(values, colors, group.names=NULL, group.widths=NULL,
                         sep.wex=0.05, pch=19, pt.cex=0.2, title=NULL,
                         title.line=0, title.cex=1, add.y.axis=TRUE,
                         y.title=NULL, y.title.line=0.5, y.axis.at=NULL,
                         y.axis.labels=NULL, y.axis.tck=-0.025,
                         parmar=c(1, 2.5, 0.25, 0.25)){
  # Ensure beeswarm & vioplot are loaded
  require(beeswarm, quietly=T)
  require(vioplot, quietly=T)

  # Summarize plotting data
  values <- lapply(values, function(v){as.numeric(v[!is.na(v)])})
  if(is.null(group.names)){
    group.names <- names(values)
  }
  drop.groups <- sapply(values, length) == 0
  if(any(drop.groups)){
    drop.group.idx <- which(drop.groups)
    values <- values[-drop.group.idx]
    colors <- colors[-drop.group.idx]
    group.names <- group.names[-drop.group.idx]
  }
  n.groups <- length(values)
  group.size <- sapply(values, length)

  # Get plot dimensions
  if(!is.null(group.widths)){
    if(length(group.widths) != length(values)){
      stop(paste("Length of custom `group.widths` does not match number of",
                 "groups in `values`"))
    }
    # Re-normalize custom widths
    group.widths <- group.widths / sum(group.widths)
  }else{
    group.widths <- group.size / sum(group.size)
  }
  group.lefts <- c(0, cumsum(group.widths)[-n.groups]) + (sep.wex * (0:(n.groups-1)))
  group.rights <- group.widths + group.lefts
  group.mids <- (group.lefts + group.rights) / 2
  xlims <- c(-sep.wex, max(group.rights))
  ylims <- range(unlist(values), na.rm=T)

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="r", yaxs="r")

  # Add X axis
  x.axis.at <- smart.spacing(group.mids, min.dist=0.15)
  axis(1, at=x.axis.at, line=-1, tick=F, labels=group.names, xpd=T)

  # Add Y axis
  if(add.y.axis){
    clean.axis(2, at=y.axis.at, labels=y.axis.labels,
               cex.axis=5/6, infinite=TRUE, label.line=-0.7,
               title.line=y.title.line, title=y.title, tck=y.axis.tck)
  }

  # Add title
  mtext(3, text=title, line=title.line, cex=title.cex)

  # Add boxplots and swarms
  sapply(1:n.groups, function(i){
    if(length(values[[i]] > 0)){
      box.buffer <- min(group.widths[i] / 8, sep.wex / 2)
      metrics <- summary(values[[i]])
      box.x <- (group.lefts + group.mids)[i]/2
      q1q3 <- metrics[c("1st Qu.", "3rd Qu.")]
      segments(x0=box.x, x1=box.x,
               y0=max(metrics["Min."], q1q3[1] - (1.5 * diff(q1q3))),
               y1=min(metrics["Max."], q1q3[2] + (1.5 * diff(q1q3))),
               lwd=3, lend="butt", col=colors[i])
      rect(xleft=group.lefts[i] + box.buffer,
           xright=group.mids[i] - box.buffer,
           ybottom=q1q3[1], ytop=q1q3[2],
           border=NA, bty="n", col=colors[i])
      segments(x0=group.lefts[i],  x1=group.mids[i],
               y0=metrics["Median"], y1=metrics["Median"],
               col="white", lend="butt", lwd=1.5)
      beeswarm(values[[i]], at=group.mids[i], side=1,
               pch=pch, cex=pt.cex, col=colors[i],
               corral="wrap", corralWidth=group.widths[i]/2,
               method="swarm", priority="density", add=T)
    }
  })
}


#' Kaplan-Meier Curves
#'
#' Plot Kaplan-Meier curves for one or more datasets or strata
#'
#' @param surv.models List of one or more [`survival::summary.survfit`] objects
#' @param colors Vector of colors for the list elements in `surv.models`
#' @param group.names (Optional) group names to assign to each list element
#' in `surv.models`
#' @param km.lwd Line width for Kaplan-Meier curves \[default: 3\]
#' @param ci.alpha Transparency value `alpha` for confidence interval
#' shading \[default: 0.1\]
#' @param legend Should a legend be plotted?
#' @param legend.names (Optional) mapping of `values` to labels for legend
#' @param legend.label.spacing Minimum vertical spacing between legend
#' labels \[default: 0.075\]
#' @param legend.label.cex Character expansion value for text labels in
#' legend \[default: 1\]
#' @param title (Optional) Title for plot
#' @param y.title Title for Y-axis \[default: "Survival Probability"\]
#' @param xlims (Optional) two-element vector of start and stop values for
#' X-axis, specified in days (or years if `time.is.days` is `FALSE`)
#' @param x.label.line Line for X-axis labels \[default: -0.75\]
#' @param x.title.line Line for X-axis title \[default: 0\]
#' @param x.tck X-axis tick length \[default: -0.0175\]
#' @param time.is.days Should time be interpreted as days \[defualt: TRUE\]
#' @param parmar Margin values passed to par()
#'
#' @seealso [`survival::Surv`], [`survival::survfit`], [`survival::summary.survfit`]
#'
#' @export km.curve
#' @export
km.curve <- function(surv.models, colors, group.names=NULL, km.lwd=3, ci.alpha=0.1,
                     legend=TRUE, legend.names=NULL, legend.label.spacing=0.075,
                     legend.label.cex=1, title=NULL, y.title="Survival Probability",
                     xlims=NULL, x.label.line=-0.75, x.title.line=0, x.tck=-0.0175,
                     time.is.days=TRUE, parmar=c(2, 3, 0.25, 4)){
  # Ensure survival library is loaded within function scope
  require(survival, quietly=TRUE)

  # Get plotting values
  if(is.null(group.names)){
    group.names <- names(surv.models)
  }
  n.groups <- length(surv.models)
  if(is.null(legend.names)){
    legend.names <- names(surv.models)
  }
  if(is.null(xlims)){
    xlims <- c(0, max(sapply(surv.models, function(ss){max(ss$time, na.rm=T)})))
  }

  # Prep plot area
  prep.plot.area(xlims, c(0, 1.025), parmar)
  if(time.is.days){
    x.ax.step <- max(c(floor(xlims[2] / (365*6)), 1))
    x.ax.years <- seq(0, xlims[2]/365, by=x.ax.step)
  }else{
    x.ax.step <- max(c(floor(xlims[2] / 6), 1))
    x.ax.years <- seq(0, xlims[2], by=x.ax.step)
  }

  # Add confidence intervals
  # Loop over this twice: first to lay white backgrounds, then add colors
  for(layer in c("white", "colors")){
    sapply(1:n.groups, function(i){
      n.times <- length(surv.models[[i]]$time)
      if(n.times > 1){
        x.bottom <- c(0, RLCtools::stretch.vector(surv.models[[i]]$time, 2)[-2*n.times])
        x.top <- rev(x.bottom)
        y.bottom <- c(1, 1, RLCtools::stretch.vector(surv.models[[i]]$lower, 2)[-c(2*n.times-c(0, 1))])
        y.top <- rev(c(1, 1, RLCtools::stretch.vector(surv.models[[i]]$upper, 2)[-c(2*n.times-c(0, 1))]))
        if(layer == "white"){
          ci.col <- "white"
        }else{
          ci.col <- adjustcolor(colors[[i]], alpha=ci.alpha)
        }
        polygon(x=c(x.bottom, x.top), y=c(y.bottom, y.top),
                border=NA, bty="n",
                col=ci.col)
      }
    })
  }

  # Add K-M curves
  sapply(1:n.groups, function(i){
    n.times <- length(surv.models[[i]]$time)
    # If summary.survfit returns no data, this is either because
    # there are no patients in this group or nobody died.
    # If the latter, we can plot as a flat line at Y=1 until rmean.endtime (I think?)
    if(surv.models[[i]]$n > 0){
      if(n.times == 0){
        x <- c(0, surv.models[[i]]$rmean.endtime)
        y <- c(1, 1)
      }else{
        x <- c(0, RLCtools::stretch.vector(surv.models[[i]]$time, 2))
        y <- c(1, 1, RLCtools::stretch.vector(surv.models[[i]]$surv, 2))[1:length(x)]
      }
      points(x, y, type="l", col=colors[[i]], lwd=km.lwd)
    }
  })

  # Add axes
  clean.axis(1, at=if(time.is.days){x.ax.years*365}else{x.ax.years},
             labels=x.ax.years, infinite=TRUE,
             title="Years", label.line=x.label.line,
             title.line=x.title.line, tck=x.tck)
  clean.axis(2, title=y.title, infinite=FALSE, tck=-0.0175)
  mtext(title, side=3, line=0)

  # Add legend
  if(legend){
    final.y <- sapply(surv.models, function(ss){
      # In the case of no events, ss will have zero rows
      # We can default to Y=1 in this case
      if(length(ss$time) == 0){1}else{
        dist.to.rb <- ss$time - xlims[2]
        if(any(dist.to.rb > 0)){
          closest <- which(dist.to.rb == min(dist.to.rb[which(dist.to.rb >= 0)],
                                             na.rm=T) & dist.to.rb >= 0)
          closest <- max(c(1, closest-1))
        }else{
          closest <- which(dist.to.rb == max(dist.to.rb, na.rm=T))
          closest <- max(c(1, closest))
        }
        ss$surv[closest]
      }
    })
    yaxis.legend(legend.names[order(final.y)], x=xlims[2] + (0.05*diff(xlims)),
                 y.positions=final.y[order(final.y)], label.cex=legend.label.cex,
                 min.label.spacing=legend.label.spacing,
                 sep.wex=0.05*diff(xlims), colors=colors[order(final.y)])
  }
}


#' Density with Outlier Rug
#'
#' Plot one-dimensional density with a marginal rug of ticks for outlier observations
#'
#' @param vals Numeric values to plot
#' @param style Style of density plot; either "density" (default) for kernel
#' density estimator or "hist"/"histogram" for histogram
#' @param min.complexity Minimum number of unique values in `vals` before
#' automatically defaulting to `style == "hist"` \[default: 30\]
#' @param bw.adj Bandwidth adjustment for density estimation. See `adjust`
#' in [stats::density()].
#' @param min.bin.width Minimum permissible bin width. Only used if
#' `style == "hist"`.
#' \[default: automatically determine optimal bin width\]
#' @param outlier.lower.bound Threshold below which an observation is treated
#' as an outlier \[default: Q1 - 3*IQR\]
#' @param outlier.upper.bound Threshold below which an observation is treated
#' as an outlier \[default: Q3 + 3*IQR\]
#' @param outlier.tick.hex Relative height expansion of outlier ticks
#' \[default: 0.02\]
#' @param color Color for density area \[default: "gray80"\]
#' @param border Color for density border \[default: "black"\]
#' @param outlier.color Color for outlier ticks \[default: same value is `border`\]
#' @param title (Optional) Title for plot
#' @param title.line Line for `title` \[default: 0.1\]
#' @param xlims Enforce custom X-axis limits \[default: plot entire range of values\]
#' @param x.title Title for X-axis \[default: no title\]
#' @param x.title.line Line for `x.title` \[default: 0.5\]
#' @param x.label.units Units for X-axis labels; passed to [RLCtools::clean.axis()]
#' @param max.x.ticks Maximum number of ticks for X-axis \[default: 5\]
#' @param add.y.axis Should a Y-axis be plotted? \[default: TRUE\]
#' @param y.title Title for Y-axis \[default: no title\]
#' @param y.title.line Line for `y.title` \[default: 0.5\]
#' @param outlier.lwd Line width for outlier ticks \[default: 1/3\]
#' @param parmar Margin values passed to par()
#'
#' @details When `style == "density"`, this function will attempt to approximate
#' the equivalent count values for the Y-axis, although these values will depend
#' on `bw.adj` and can sometimes be misleading
#'
#' @seealso [stats::density()], [RLCtools::clean.axis()]
#'
#' @export density.w.outliers
#' @export
density.w.outliers <- function(vals, style="density", min.complexity=30, bw.adj=1,
                               min.bin.width=NULL, outlier.lower.bound=NULL,
                               outlier.upper.bound=NULL, outlier.tick.hex=0.02,
                               color="gray70", border="black", outlier.color=NULL,
                               title=NULL, title.line=0.1, xlims=NULL,
                               x.title=NULL, x.title.line=0.5, x.label.units=NULL,
                               max.x.ticks=6, add.y.axis=TRUE, y.title=NULL,
                               y.title.line=0.5, outlier.lwd=1/3,
                               parmar=c(2, 2, 0.35, 0.35)){
  # Determine outlier points
  vals <- as.numeric(vals)
  n.unique.vals <- length(unique(vals))
  if(n.unique.vals < 30){
    is.out <- rep(FALSE, length(vals))
  }else{
    is.out <- label.outliers(vals,
                             fixed.min=outlier.lower.bound,
                             fixed.max=outlier.upper.bound)
  }
  n.out <- length(which(is.out))

  # Scale axis limits according to units of values
  xmin <- min(vals, na.rm=T)
  xmax <- max(vals, na.rm=T)
  if(!is.null(x.label.units)){
    if(x.label.units == "percent"){
      xmin <- max(c(0, min(vals, na.rm=T)))
      xmax <- min(c(1, max(vals, na.rm=T)))
    }
  }
  if(is.null(xlims)){
    xlims <- c(xmin, xmax)
  }

  # Compute density over non-outlier points
  histogram <- style %in% c("hist", "histogram") | n.unique.vals < 30
  if(histogram){
    if(is.null(min.bin.width)){
      breaks <- seq(xmin, xmax, length.out=50 / bw.adj)
    }else{
      breaks <- seq(xmin, xmax, by=min.bin.width)
      breaks <- c(breaks, max(breaks, na.rm=T) + min.bin.width)
    }
    v.dens <- hist(vals[which(!is.out)], plot=F, breaks=breaks)
    ymax.label <- ymax <- max(v.dens$counts, na.rm=T)
  }else{
    v.dens <- density(vals[which(!is.out)], adjust=bw.adj)
    # Approximate count estimates with hist
    max.count <- max(hist(vals[which(!is.out)],
                          breaks=seq(min(vals[which(!is.out)], na.rm=T),
                                     max(vals[which(!is.out)], na.rm=T),
                                     length.out=length(v.dens$x)), plot=F)$counts)
    v.dens$y <- max.count * v.dens$y / max(v.dens$y, na.rm=T)
    ymax <- max(v.dens$y, na.rm=T)
  }

  # Prepare plot area
  prep.plot.area(xlims, c(0, ymax), parmar=parmar, xaxs="r")

  # Add ticks for outliers
  if(is.null(outlier.color)){
    outlier.color <- border
  }
  if(n.out > 0){
    out.tick.y1 <- outlier.tick.hex * diff(par("usr")[3:4])
    segments(x0=vals[is.out], x1=vals[is.out], y0=rep(0, n.out),
             y1=rep(out.tick.y1, n.out), xpd=T, lwd=outlier.lwd,
             col=outlier.color)
  }

  # Add density
  if(histogram){
    n.breaks <- length(v.dens$breaks)
    rect(xleft=v.dens$breaks[-n.breaks], xright=v.dens$breaks[-1],
         ybottom=rep(0, length=n.breaks-1), ytop=v.dens$counts,
         col=color, border=border, xpd=T)
  }else{
    polygon(x=c(v.dens$x, rev(v.dens$x)),
            y=c(v.dens$y, rep(0, length(v.dens$x))),
            col=color, border=border, xpd=T)
  }

  # Add axes & title
  clean.axis(1, title=x.title, title.line=x.title.line, max.ticks=max.x.ticks,
             label.units=x.label.units, infinite=TRUE, label.line=-0.9)
  if(add.y.axis){
    clean.axis(2, label.units=y.title, title=y.title,
               title.line=y.title.line, infinite=TRUE)
  }
  mtext(3, text=title, line=title.line)
}


#' Stacked barplot
#'
#' Produce a stacked barplot of one or two categorical variables
#'
#' @param major.values "Major" axis values; used for Y-axis groupings
#' @param minor.values "Minor" axis values; if provided, will be depicted as
#' stacked bars within each major axis group
#' @param colors Color assignments for "minor" axis values; see `Details`
#' @param inner.borders Color assignments for minor segment inner.borders
#' \[default: "white"\]
#' @param outer.borders Color assignment for outer border \[default: "black\]
#' @param outer.border.lwd `lwd` parameter for outer border \[default: 1\]
#' @param as.proportion Should values be plotted on a relative (proportional /
#' percentage) scale? \[default: FALSE, i.e., plot raw counts\]
#' @param add.major.labels Should major (y-axis) groups be labeled in the margin?
#' \[default: TRUE\]
#' @param x.axis.side `side` value for x-axis; `NA` will disable X-axis plotting.
#' \[default: 3 (top axis)\]
#' @param x.title Optional title for X axis
#' @param x.title.line Value of `title.line` passed to [RLCtools::clean.axis()]
#' @param x.label.line Value of `label.line` passed to [RLCtools::clean.axis()]
#' @param x.axis.tck Value of `tck` passed to [RLCtools::clean.axis()]
#' @param y.label.cex Cex parameter for Y-axis "major" group labels
#' @param bar.hex Width of bars relative to size of gap between bars. Setting
#' `bar.hex == 1` will leave no gap between the bars
#' @param add.legend Should a legend of minor value colors be added to the
#' bottom-right of the plot? \[default: add legend\]
#' @param legend.xadj Legend x-position adjustment, in relative user units.
#' Only relevant if `add.legend` is `TRUE`. \[default: -0.075]
#' @param major.legend Should a major color legend be added between Y axis and
#' labels? \[default: FALSE\]
#' @param major.legend.colors Named vector mapping `major.values` to colors for
#' `major.legend` \[default: uniform greyscale\]
#' @param major.legend.xadj X adjustment scalar for major legend
#' \[default: -0.04\]
#' @param minor.labels.on.bars Should minor value labels be printed on bars,
#' where space is permitting? \[default: FALSE\]
#' @param minor.label.letter.width Proportion of total X-axis to apportion per
#' letter for minor labels. Only used if `minor.labels.on.bars` is `TRUE`.
#' \[default: 0.05\]
#' @param minor.label.color Color for minor label text. If `NULL`, color will be
#' dynamically optimized to be most visible against each bar's background color.
#' @param minor.label.cex Character expansion parameter for minor label text
#' \[default: 5/6\]
#' @param annotate.counts Should exact bar counts be annotated at the tip of
#' each bar? \[default: no annotations\]
#' @param end.label.xadj End-label x-position adjustment, in relative user units.
#' Only relevant if `annotate.counts` is `TRUE`. \[default: -0.025]
#' @param end.label.cex Character expansion parameter for end label text
#' \[default: 5/6\]
#' @param orient Should the bar length be increasing to the `right` or
#' `left`? \[default: `right`\]
#' @param custom.major.order Specific order of major values to use for bars
#' @param custom.minor.order Specific order of minor values to use for bars
#' @param sort.minor Should minor labels be sorted by abundance prior to plotting?
#' Only used if `custom.minor.order` is `NULL`. \[default: alphabetically order
#' minor labels]
#' @param parmar Margin values passed to par()
#'
#' @param details
#' By default, `colors` will uniformly sample a grayscale palette and
#' assign one color to each unique value present in `minor.values` (or
#' `major.values`, if `minor.values` is not specified). Alternatively, custom
#' color assignments can be specified in the following ways:
#'
#' 1. As a function to generate a palette, which will be used in place of the
#' grayscale palette described above
#'
#' 2. As a named vector of color codes recognized by R, where the names of the
#' vector elements map onto all unique values of `minor.values`
#'
#' @examples
#' set.seed(2024)
#' x <- sample(1:5, 100, replace=TRUE)
#' y <- sample(letters[1:4], 100, replace=TRUE)
#' stacked.barplot(x, y)
#'
#' @export stacked.barplot
#' @export
stacked.barplot <- function(major.values, minor.values=NULL, colors=NULL,
                            inner.borders=NULL, outer.borders=NULL,
                            outer.border.lwd=1, as.proportion=FALSE,
                            add.major.labels=TRUE, x.axis.side=3, x.title=NULL,
                            x.title.line=0.3, x.label.line=-0.65,
                            x.axis.tck=-0.025, y.label.cex=5/6, bar.hex=0.8,
                            add.legend=TRUE, legend.xadj=-0.075,
                            major.legend=FALSE, major.legend.colors=NULL,
                            major.legend.xadj=-0.04, minor.labels.on.bars=FALSE,
                            minor.label.letter.width=0.05,
                            minor.label.color=NULL, minor.label.cex=5/6,
                            annotate.counts=FALSE,end.label.xadj=-0.025,
                            end.label.cex=5/6, orient="right",
                            custom.major.order=NULL, custom.minor.order=NULL,
                            sort.minor=FALSE, parmar=c(0.5, 3, 2.5, 0.5)){
  # Check if minor values are provided
  no.minor <- is.null(minor.values)
  if(no.minor){
    minor.values <- major.values
  }

  # Handle NAs
  # Minor NAs will be filled
  minor.values[which(is.na(minor.values))] <- "N.S."
  # Major NAs will be dropped outright
  drop.idx <- is.na(major.values)
  if(any(drop.idx)){
    major.values <- major.values[-which(drop.idx)]
    minor.values <- minor.values[-which(drop.idx)]
  }

  # Organize plot data
  major.table <- sort(table(major.values), decreasing=TRUE)
  minor.table <- table(sort(minor.values))
  if(!is.null(custom.minor.order)){
    if(!all(length(union(names(minor.table), custom.minor.order)) %in% c(length(minor.table), length(custom.minor.order)))){
      warning(paste("Not all values of `custom.minor.order` appear in minor",
                    "values (or vice versa)"))
    }
    minor.table <- minor.table[custom.minor.order]
  }else if(sort.minor){
    minor.table <- sort(minor.table, decreasing=TRUE)
  }
  if(!is.null(custom.major.order)){
    if(!all(length(union(names(major.table), custom.major.order)) %in% c(length(major.table), length(custom.major.order)))){
      warning(paste("Not all values of `custom.major.order` appear in major",
                    "values (or vice versa)"))
    }
    major.table <- major.table[custom.major.order]
    if(no.minor){
      minor.table <- minor.table[custom.major.order]
    }
  }
  plot.df <- do.call("rbind", lapply(names(major.table), function(major){
    sapply(names(minor.table), function(minor){
      length(which(major.values == major & minor.values == minor))
    })
  }))
  if(as.proportion){
    plot.df <- t(apply(plot.df, 1, function(vals){vals / sum(vals)}))
  }
  plot.df <- as.data.frame(plot.df)
  rownames(plot.df) <- names(major.table)
  bar.hex <- min(c(bar.hex, 1))

  # Handle color assignment
  if(is.null(colors)){
    colors <- greyscale.palette(length(minor.table))
    names(colors) <- names(minor.table)
  }
  if(is.null(inner.borders)){
    inner.borders <- rep("white", length(minor.table))
    names(inner.borders) <- names(minor.table)
  }
  if(length(inner.borders) < length(minor.table)){
    inner.borders <- rep(inner.borders,
                         length(minor.table) / length(inner.borders))
    names(inner.borders) <- names(minor.table)
  }
  if(is.null(outer.borders)){
    outer.borders <- rep("black", length(major.table))
    names(outer.borders) <- names(major.table)
  }
  if(length(outer.borders) < length(major.table)){
    outer.borders <- rep(outer.borders,
                         length(major.table) / length(outer.borders))
    names(outer.borders) <- names(major.table)
  }

  # Prepare plot area
  xlims <- c(0, if(as.proportion){1}else{max(major.table, na.rm=T)})
  ylims <- c(length(major.table), 0)
  if(orient == "left"){
    xlims <- rev(xlims)
    ylims <- rev(ylims)
  }
  prep.plot.area(xlims, ylims, parmar=parmar)
  if(add.major.labels){
    sapply(1:nrow(plot.df), function(y){
      axis(if(orient == "left"){4}else{2},
           at=y-0.5, tick=F, las=2, cex.axis=y.label.cex,
           labels=rownames(plot.df)[y], line=if(major.legend){-0.4}else{-0.9})
    })
  }
  if(!is.na(x.axis.side)){
    clean.axis(x.axis.side,
               label.units=if(as.proportion){"percent"}else{"count"},
               infinite.positive=TRUE, title=x.title, title.line=x.title.line,
               label.line=x.label.line, tck=x.axis.tck)
  }

  # Add bars
  r.bar <- bar.hex / 2
  sapply(1:length(major.table), function(major.idx){
    rect(xleft=c(0, cumsum(unlist(plot.df[major.idx, ])))[-(ncol(plot.df)+1)],
         xright=cumsum(unlist(plot.df[major.idx, ])),
         ybottom=major.idx - 0.5 - r.bar,
         ytop=major.idx -0.5 + r.bar,
         col=colors[colnames(plot.df)],
         border=inner.borders[colnames(plot.df)], lwd=0.5)
    rect(xleft=0, xright=apply(plot.df, 1, sum, na.rm=T)[major.idx],
         ybottom=major.idx - 0.5 - r.bar,
         ytop=major.idx -0.5 + r.bar,
         col=NA, xpd=T, border=outer.borders[rownames(plot.df)[major.idx]],
         lwd=outer.border.lwd)
  })

  # If optioned, print key codes on physical bars where space permits
  if(minor.labels.on.bars){
    prop.df <- floor((plot.df / max(xlims)) / minor.label.letter.width)
    sapply(1:ncol(prop.df), function(k){
      if(any(prop.df[, k] > 0)){
        longest.idx <- head(which.max(prop.df[, k]), 1)
        longest.x <- mean(c(0, cumsum(unlist(plot.df[longest.idx, ])))[c(k, k+1)])
        label <- colnames(plot.df)[k]
        room <- prop.df[longest.idx, k]
        if(nchar(label) > room){
          if(room == 1){
            label <- substr(label, 1, 1)
          }else{
            label <- paste(substr(label, 1, room), ".", sep="")
          }
        }
        lab.col <- minor.label.color
        if(is.null(lab.col)){
          lab.col <- optimize.label.color(colors[colnames(plot.df)[k]])
        }
        text(x=longest.x, y=longest.idx-0.5, cex=minor.label.cex,
             labels=label, col=lab.col)
      }
    })
  }

  # Add minor legend, if optioned
  if(add.legend){
    legend.colors <- colors[intersect(names(sort(-minor.table)), names(colors))]
    if(orient == "left"){
      legend(x=par("usr")[1] + (legend.xadj * diff(par("usr")[1:2])),
             y=par("usr")[4],
             names(legend.colors), fill=legend.colors, cex=5/6, bty="n", xpd=T)
    }else{
      legend("bottomright", names(legend.colors), fill=legend.colors,
             cex=5/6, bty="n", xpd=T)
    }
  }

  # Add major legend in margin, if optioned
  if(major.legend){
    if(is.null(major.legend.colors)){
      major.legend.colors <- greyscale.palette(length(major.table))
      names(major.legend.colors) <- rownames(plot.df)
    }
    if(orient == "left"){
      maj.leg.xadj <- -major.legend.xadj
    }else{
      maj.leg.xadj <- major.legend.xadj
    }
    points(x=rep(maj.leg.xadj*diff(par("usr")[1:2]), nrow(plot.df)),
           y=(1:nrow(plot.df)) - 0.5, pch=23, xpd=T, col="black",
           bg=major.legend.colors[rownames(plot.df)])
  }

  # Add count labels, if optioned
  if(annotate.counts){
    bar.ends <- apply(plot.df, 1, sum, na.rm=T)
    end.xadj <- if(orient == "left"){-end.label.xadj}else{end.label.xadj}
    text(x=bar.ends + (end.xadj * diff(par("usr")[1:2])),
         y=(1:nrow(plot.df)) - 0.5, labels=prettyNum(bar.ends, big.mark=","),
         cex=end.label.cex, pos=if(orient == "left"){2}else{4}, xpd=T)
  }
}
