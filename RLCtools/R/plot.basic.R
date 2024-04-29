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
#' @param X PC to be plotted on the X-axis (see `Details`)
#' @param Y PC to be plotted on the Y-axis (see `Details`)
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
#' @param parmar Numeric vector of values to pass to `par(mar)`
#'
#' @details If `legend.vals` is not provided (or is `NULL`), no legend will be added.
#'
#' @export scatterplot
#' @export
scatterplot <- function(pcs, X, Y, colors=NULL, title=NULL,
                        x.label.line=NULL, x.title=NULL, x.title.line=0.5, xlims=NULL,
                        y.label.line=NULL, y.title=NULL, y.title.line=0.5, ylims=NULL,
                        legend.vals=NULL, legend.labels=NULL,
                        cex=0.3, parmar=c(2.5, 2.5, 1, 1)){
  x <- as.numeric(X)
  y <- as.numeric(Y)
  if(is.null(xlims)){
    xlims <- range(x)
  }
  if(is.null(ylims)){
    ylims <- range(y)
  }
  if(is.null(colors)){
    colors <- rep("gray70", length(x))
  }

  # Prepare plot area
  prep.plot.area(xlims, ylims, parmar=parmar, xaxs="r", yaxs="r")
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
#' @param data list of [density] objects to plot
#' @param names optional list of names for Y axis \[default: take names from data\]
#' @param hill.overlap relative fraction of overlap bewtween adjacent hills \[default: 0.35\]
#' @param xlims custom X axis limits
#' @param ylims custom Y axis limits
#' @param fill vector of polygon fill colors \[default: "grey70"\]
#' @param border vector of hill border colors \[default: "grey35"\]
#' @param border.lwd line width for hill borders \[default: 2\]
#' @param parmar vector of values passed to par(mar)
#'
#' @export ridgeplot
#' @export
ridgeplot <- function(data, names=NULL, hill.overlap=0.35, xlims=NULL, x.axis=TRUE,
                      fill="grey70", border="grey35", border.lwd=2,
                      parmar=c(2.5, 3, 0.25, 0.25)){
  # Get names before manipulating data
  if(is.null(names)){
    names <- names(data)
    if(is.null(names)){
      names <- 1:length(data)
    }
  }

  # Scale Y values of data to [0, hill.overlap]
  for(i in 1:length(data)){
    y <- data[[i]]$y
    data[[i]]$y <- (1 + hill.overlap) * (y / max(y))
  }

  # Get plot dimensions
  if(is.null(xlims)){
    xlims <- c(min(sapply(data, function(d){min(d$x)})),
               max(sapply(data, function(d){max(d$x)})))
  }
  ylims <- c(0, length(data) + hill.overlap)

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="i", yaxs="r")
  median(unlist(sapply(data, function(df){df$y})), na.rm=T)
  if(x.axis){
    clean.axis(1, title="Values", infinite=TRUE)
  }
  axis(2, at=(1:length(data)) - 0.5, tick=F, las=2, line=-0.8, labels=names)

  # Add hills
  sapply(length(data):1, function(i){
    abline(h=i-1, col="gray85")
    x <- c(data[[i]]$x, rev(data[[i]]$x))
    y <- c(data[[i]]$y, rep(0, times=length(data[[i]]$y)))+i-1
    polygon(x, y, border=fill, col=fill, lwd=border.lwd)
    points(data[[i]]$x, data[[i]]$y+i-1, type="l", lwd=border.lwd, col=border)
  })
}


#' Quantile-quantile plot
#'
#' Generate a Q-Q plot of association stats vs. uniform null
#'
#' @param pvals Numeric vector of untransformed P-values
#' @param cutoff P-value threshold for significance
#' @param do.fdr Should points meeting `fdr.cutoff` be accented? \[default: TRUE\]
#' @param fdr.cutoff Cutoff for highlighting FDR-significant points \[default: 0.01\]
#' @param print.stats Should genomic inflation statistic be printed on plot? \[default: FALSE\]
#' @param title (Optional) Plot title
#' @param title.line Line for `title` \[default: 0\]
#' @param xmax Maximum X-value to plot \[default: include all points\]
#' @param ymax Maximum Y-value to plot \[default: smallest P-value\]
#' @param label.cex Scaling factor for label text \[default: 1\]
#' @param pt.color Color for all points \[default: "grey35"\]
#' @param pt.cex Scaling factor for all points \[default: 0.35\]
#' @param fdr.color Color for FDR-significant points \[default: "grey5"\]
#' @param fdr.cex Scaling factor for FDR-significant points \[default: 0.7\]
#' @param parmar Value of `mar` passed to `par()`
#'
#' @return None
#'
#' @export plot.qq
#' @export
plot.qq <- function(pvals, cutoff=NULL, do.fdr=TRUE, fdr.cutoff=0.01, print.stats=FALSE, 
                    title=NULL, title.line=0, xmax=NULL, ymax=NULL, label.cex=1, 
                    pt.color="grey35", pt.cex=0.35, fdr.color="grey5", fdr.cex=0.7, 
                    parmar=c(2.25, 2.5, 0.25, 0.25)){
  # Format P-values
  if (!is.numeric(pvals)){
    stop("P values must be numeric.")
  }
  keep.idxs <- which(!is.na(pvals) & !is.nan(pvals) & !is.null(pvals) &
                       is.finite(pvals) & pvals <= 1 & pvals >= 0)
  pvals <- pvals[keep.idxs]
  colors <- rep(pt.color, length(pvals))
  pw.cex <- rep(pt.cex, length(pvals))
  pvals <- pvals[order(pvals)]
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

  # Prep plot area and add null confidence interval
  prep.plot.area(c(0, 1.1 * xmax), c(0, 1.1 * ymax), parmar=parmar)
  polygon(x=conf.int[, 1], y=conf.int[, 2], col="gray90", border=NA)
  abline(0, 1, col="gray50")
  if(!is.null(cutoff)){
    abline(h=log.cutoff, lty=2)
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
  clean.axis(1, title=expression(Expected ~ ~-log[10] ~ italic(P)),
             label.line=-0.75, title.line=0.35, infinite=TRUE)
  clean.axis(2, title=expression(Observed ~ ~-log[10] ~ italic(P)),
             title.line=0.15, infinite=TRUE)
  mtext(3, line=title.line, text=title)
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
#' @param parmar Value of `mar` passed to [par()]
#'
#' @export manhattan
#' @export
manhattan <- function(stats, p.column, chrom.colors, gw.sig, pt.cex=0.3,
                      title=NULL, parmar=c(2, 2.25, 0.25, 0.3)){
  # Prep plot area
  xlims <- range(stats$pos)
  ylims <- c(0, ceiling(max(c(stats[, p.column], gw.sig), na.rm=T) + 1))
  prep.plot.area(xlims, ylims, parmar=parmar)
  abline(h=gw.sig, lty=5)

  # Add points
  set.seed(2024)
  stats <- stats[sample(1:nrow(stats), size=nrow(stats)), ]
  points(stats[, c("pos", p.column)], pch=19, cex=pt.cex, col=chrom.colors[stats$chrom])

  # Add axes
  clean.axis(1, at=c(0, cumsum(contig.lengths)), tck=0.025, infinite=TRUE,
             labels=NA, title="Genomic coordinate", title.line=0)
  sapply(1:length(contig.lengths), function(x){
    axis(1, at=((c(0, cumsum(contig.lengths[-24])) + cumsum(contig.lengths))/2)[x],
         tick=F, cex.axis=4/6, labels=gsub("^chr", "", names(contig.lengths)[x]),
         col.axis=chrom.colors[x],
         line=c(-0.8, -1.3)[as.integer((x %% 2) == 1) + 1])
  })
  clean.axis(2, title=bquote(-log[10] ~ italic(P)), infinite.positive=T, title.line=0.2)
  mtext(3, text=title, line=0)
}

