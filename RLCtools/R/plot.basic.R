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
                      fill=NULL, border=NULL, border.lwd=2,
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
    polygon(x, y, border=fill[i], col=fill[i], lwd=border.lwd)
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


#' Scale-Aware Barplot Cluster
#'
#' Plot a set of stacked barplots scaled proportional to set size
#'
#' @param values List of character vectors of values to plot
#' @param colors Named vector of colors to map to `values`
#' @param group.names (Optional) group names to assign to each list element in `values`
#' @param sep.wex Relative width scalar for whitespace on the X-axis
#' between groups \[default: 0.05\]
#' @param title (Optional) title to be printed in top-right corner
#' @param legend Should a legend be plotted?
#' @param legend.names (Optional) mapping of `values` to labels for legend
#' @param parmar Margin values passed to par()
#'
#' @details If `values` is supplied as a named list, those names will be used as
#' `group.names` unless `group.names` is explicitly specified. Otherwise, `group.names`
#' will be set to the ordinal value of each group in `values`.
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
  pct.df <- do.call("cbind", lapply(1:n.groups, function(i){as.numeric(counts.df[, i] / group.size[i])}))
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
  clean.axis(2, at=seq(0, 1, 0.25), labels=rev(paste(seq(0, 100, 25), "%", sep="")),
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
    rect(xleft=group.lefts[i], xright=group.rights[i], ybottom=0, ytop=1, col=NA, xpd=T)
  })
}


#' Scale-Aware Beeswarm Cluster
#'
#' Plot a set of beeswarm distributions scaled proportional to set size
#'
#' @param values List of numeric vectors of values to plot
#' @param colors Vector of colors for the list elements in `values`
#' @param group.names (Optional) group names to assign to each list element in `values`
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
#' @param parmar Margin values passed to par()
#'
#' @details If `values` is supplied as a named list, those names will be used as
#' `group.names` unless `group.names` is explicitly specified. Otherwise, `group.names`
#' will be set to the ordinal value of each group in `values`.
#'
#' @seealso [RLCtools::scaled.bars], [RLCtools::clean.axis]
#'
#' @export scaled.swarm
#' @export
scaled.swarm <- function(values, colors, group.names=NULL, sep.wex=0.05,
                         pch=19, pt.cex=0.2, title=NULL, title.line=0,
                         title.cex=1, add.y.axis=TRUE, y.title=NULL,
                         y.title.line=0.5, y.axis.at=NULL, y.axis.labels=NULL,
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
  group.widths <- group.size / sum(group.size)
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
               title.line=y.title.line, title=y.title)
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
#' @param group.names (Optional) group names to assign to each list element in `surv.models`
#' @param ci.alpha Transparency value `alpha` for confidence interval shading \[default: 0.15\]
#' @param legend Should a legend be plotted?
#' @param legend.names (Optional) mapping of `values` to labels for legend
#' @param legend.label.spacing Minimum vertical spacing between legend labels \[default: 0.075\]
#' @param title (Optional) Title for plot
#' @param y.title Title for Y-axis \[default: "Survival Probability"\]
#' @param xlims (Optional) two-element vector of start and stop values for X-axis, in days
#' @param parmar Margin values passed to par()
#'
#' @seealso [`survival::Surv`], [`survival::survfit`], [`survival::summary.survfit`]
#'
#' @export km.curve
#' @export
km.curve <- function(surv.models, colors, group.names=NULL, ci.alpha=0.15,
                     legend=TRUE, legend.names=NULL, legend.label.spacing=0.075,
                     title=NULL, y.title="Survival Probability",
                     xlims=NULL, parmar=c(2, 3, 0.25, 4)){
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
  x.ax.step <- max(c(floor(xlims[2] / (365*6)), 1))
  x.ax.years <- seq(0, xlims[2]/365, by=x.ax.step)

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
        polygon(x=c(x.bottom, x.top), y=c(y.bottom, y.top), border=NA, bty="n",
                col=if(layer == "white"){"white"}else{adjustcolor(colors[[i]], alpha=ci.alpha)})
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
      points(x, y, type="l", col=colors[[i]], lwd=3)
    }
  })

  # Add axes
  clean.axis(1, at=x.ax.years*365, labels=x.ax.years, infinite=TRUE,
             title="Years", label.line=-0.75, title.line=0, tck=-0.0175)
  clean.axis(2, title=y.title, infinite=FALSE, tck=-0.0175)
  mtext(title, side=3, line=0, font=2)

  # Add legend
  if(legend){
    final.y <- sapply(surv.models, function(ss){
      # In the case of no events, ss will have zero rows
      # We can default to Y=1 in this case
      if(length(ss$time) == 0){1}else{
        dist.to.rb <- ss$time - xlims[2]
        if(any(dist.to.rb > 0)){
          closest <- which(dist.to.rb == min(dist.to.rb[which(dist.to.rb >= 0)], na.rm=T) & dist.to.rb >= 0)
          closest <- max(c(1, closest-1))
        }else{
          closest <- which(dist.to.rb == max(dist.to.rb, na.rm=T))
          closest <- max(c(1, closest))
        }
        ss$surv[closest]
      }
    })
    yaxis.legend(legend.names[order(final.y)], x=xlims[2] + (0.05*diff(xlims)),
                 y.positions=final.y[order(final.y)],
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
#' @param bw.adj Bandwidth adjustment for density estimation. See `adjust` in [stats::density()].
#' @param outlier.lower.bound Threshold below which an observation is treated
#' as an outlier \[default: Q1 - 3*IQR\]
#' @param outlier.upper.bound Threshold below which an observation is treated
#' as an outlier \[default: Q3 + 3*IQR\]
#' @param outlier.tick.hex Relative height expansion of outlier ticks \[default: 0.02\]
#' @param color Color for density area \[default: "gray80"\]
#' @param title (Optional) Title for plot
#' @param x.title Title for X-axis \[default: no title\]
#' @param x.label.units Units for X-axis labels; passed to [RLCtools::clean.axis()]
#' @param y.title Title for Y-axis \[default: no title\]
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
                               outlier.lower.bound=NULL, outlier.upper.bound=NULL,
                               outlier.tick.hex=0.02, color="gray80", title=NULL,
                               x.title=NULL, x.label.units=NULL, y.title=NULL,
                               outlier.lwd=1/3, parmar=c(2, 2, 0.35, 0.35)){
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
  xlims <- c(xmin, xmax)

  # Compute density over non-outlier points
  histogram <- style %in% c("hist", "histogram") | n.unique.vals < 30
  if(histogram){
    v.dens <- hist(vals[which(!is.out)], plot=F,
                   breaks=seq(xmin, xmax, length.out=50 / bw.adj))
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
  clean.axis(1, title=x.title, label.units=x.label.units, infinite=TRUE, label.line=-0.9)
  clean.axis(2, label.units=y.title, title=y.title, infinite=TRUE)
  mtext(3, text=title, line=0.1)

  # Add ticks for outliers
  if(n.out > 0){
    out.tick.y1 <- outlier.tick.hex * diff(par("usr")[3:4])
    segments(x0=vals[is.out], x1=vals[is.out], y0=rep(0, n.out),
             y1=rep(out.tick.y1, n.out), xpd=T, lwd=outlier.lwd)
  }

  # Add density
  if(histogram){
    n.breaks <- length(v.dens$breaks)
    rect(xleft=v.dens$breaks[-n.breaks], xright=v.dens$breaks[-1],
         ybottom=rep(0, length=n.breaks-1), ytop=v.dens$counts,
         col=color, xpd=T)
  }else{
    polygon(x=c(v.dens$x, rev(v.dens$x)),
            y=c(v.dens$y, rep(0, length(v.dens$x))),
            col=color, xpd=T)
  }
}


#' Stacked barplot
#'
#' Produce a stacked barplot of two categorical variables
#'
#' @param major.values "Major" axis values; used for Y-axis groupings
#' @param minor.values "Minor" axis values; will be depicted as stacked bars
#' within each major axis group
#' @param minor.colors Color assignments for "minor" axis values; see `Details`
#' @param x.title Optional title for top X axis
#' @param y.label.cex Cex parameter for Y-axis "major" group labels
#' @param bar.hex Width of bars relative to size of gap between bars. Setting
#' `bar.hex == 1` will leave no gap between the bars
#' @param add.legend Should a legend of minor value colors be added to the
#' bottom-right of the plot? \[default: add legend\]
#' @param parmar Margin values passed to par()
#'
#' @param details
#' By default, `minor.colors` will uniformly sample a grayscale palette and
#' assign one color to each unique value present in `minor.values`. Alternatively,
#' custom color assignments can be specified in the following ways:
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
stacked.barplot <- function(major.values, minor.values, minor.colors=NULL,
                            x.title=NULL, y.label.cex=5/6, bar.hex=0.8,
                            add.legend=TRUE, parmar=c(0.5, 3, 2.5, 0.5)){
  # Organize plot data
  major.table <- sort(table(major.values), decreasing=TRUE)
  minor.table <- table(sort(minor.values))
  plot.df <- do.call("rbind", lapply(names(major.table), function(major){
    sapply(names(minor.table), function(minor){
      length(which(major.values == major & minor.values == minor))
    })
  }))
  rownames(plot.df) <- names(major.table)
  bar.hex <- min(c(bar.hex, 1))

  # Handle color assignment
  if(is.null(minor.colors)){
    minor.colors <- colorRampPalette(c("black", "white"))(length(minor.table) + 2)[-c(1, length(minor.table) + 2)]
    names(minor.colors) <- names(minor.table)
  }

  # Prepare plot area
  prep.plot.area(c(0, max(major.table)), c(length(major.table), 0),
                 parmar=parmar)
  axis(2, at=(1:nrow(plot.df)) - 0.5, tick=F, las=2, cex=y.label.cex,
       labels=rownames(plot.df), line=-0.9)
  clean.axis(3, label.units="count", infinite.positive=TRUE, title=x.title)

  # Add bars
  r.bar <- bar.hex / 2
  sapply(1:length(major.table), function(major.idx){
    rect(xleft=c(0, cumsum(plot.df[major.idx, ]))[-ncol(plot.df)],
         xright=cumsum(plot.df[major.idx, ]),
         ybottom=major.idx - 0.5 - r.bar,
         ytop=major.idx -0.5 + r.bar,
         col=minor.colors[colnames(plot.df)],
         border=NA, bty="n")
    rect(xleft=0, xright=major.table[major.idx],
         ybottom=major.idx - 0.5 - r.bar,
         ytop=major.idx -0.5 + r.bar,
         col=NA)
  })

  # Add legend, if optioned
  if(add.legend){
    legend("bottomright", names(minor.colors), fill=minor.colors, cex=5/6, bty="n")
  }
}

