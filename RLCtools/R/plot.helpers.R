#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Plotting helper functions


#' Prepare Plot Area
#'
#' Prepare a standardized & formatted plot area
#'
#' @param xlims Range of values for X axis
#' @param ylims Range of values for Y axis
#' @param parmar Margin values passed to par()
#' @param xaxs Value of `xaxs` passed to plot()
#' @param yaxs Value of `yaxs` passed to plot()
#'
#' @examples
#' prep.plot.area(xlims=c(0, 5), ylims=(-10, 10), parmar=rep(3, 4));
#'
#' @export prep.plot.area
#' @export
prep.plot.area <- function(xlims, ylims, parmar, xaxs="i", yaxs="i"){
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims, type="n",
       xaxs=xaxs, xlab="", xaxt="n",
       yaxs=yaxs, ylab="", yaxt="n")
}


#' Clean axis
#'
#' Print a clean axis using visually pleasing defaults
#'
#' @param side Value passed to `axis()`. See `?axis` for details.
#' @param at Positions where axis ticks should be plotted \[default: `axTicks(side)`\]
#' @param labels Labels for axis ticks \[default: values of `at`\]
#' @param labels.at Positions for axis labels \[default: values of `at`\]
#' @param label.units Specify custom units for the axis label. Default of `NULL`
#' will display numeric values. Options currently include "percent" for percentages
#' and "count" for counts that will abbreviated with "k" for thousands and "M"
#' for millions. Can be overridden by supplying `labels` directly.
#' @param parse.labels Should `labels` be parsed as R expressions? \[default: FALSE\]
#' @param max.ticks Maximum number of axis ticks. Will be overridden by `at` \[default: 6\]
#' @param title Axis title
#' @param tck Value passed to `axis()`. See `?axis` for details. \[default: -0.025\]
#' @param cex.axis Value passed to `axis()`. See `?axis` for details. \[default: 5/6\]
#' @param axis.line `line` parameter for overall axis \[default: 0\]
#' @param label.line `line` parameter for axis labels \[default: -0.65\]
#' @param title.line `line` parameter for axis title \[default: 0.5\]
#' @param infinite Indicator for the axis to be extended infinitely (without ticks) \[default: FALSE\]
#' @param infinite.positive Indicator for the axis to be extended infinitely
#' in the positive direction (without ticks) \[default: FALSE\]
#' @param infinite.negative Indicator for the axis to be extended infinitely
#' in the negative direction (without ticks) \[default: FALSE\]
#'
#' @returns NULL
#'
#' @export clean.axis
#' @export
clean.axis <- function(side, at=NULL, labels=NULL, labels.at=NULL, label.units=NULL,
                       parse.labels=FALSE, max.ticks=6, title=NULL, tck=-0.025,
                       cex.axis=5/6, line=0, label.line=-0.65, title.line=0.5,
                       infinite=FALSE, infinite.positive=FALSE, infinite.negative=FALSE){
  if(infinite){axis(side, at=c(-10e10, 10e10), tck=0, labels=NA, line=line)}
  if(is.null(at)){
    at <- axTicks(side)
    if(length(at) > max.ticks){
      at <- at[c(TRUE, FALSE)]
    }
  }
  if(infinite.positive){axis(side, at=c(at[1], 10e10), tck=0, labels=NA)}
  if(infinite.negative){axis(side, at=c(-10e10, at[length(at)]), tck=0, labels=NA)}
  if(is.null(labels)){
    labels <- at
    if(!is.null(label.units)){
      if(label.units == "percent"){
        labels <- paste(100 * labels, "%", sep="")
      }
      if(label.units == "count"){
        lab.logs <- floor(log10(labels))
        raw.best <- length(which(lab.logs < 3))
        k.best <- length(which(lab.logs >= 3 & lab.logs < 6))
        M.best <- length(which(lab.logs > 5))
        if(M.best > 0){
          at <- sort(1000000 * unique(round(as.numeric(labels) / 1000000, 1)))
          labels <- prettyNum(at / 1000000, big.mark=",")
          labels <- paste(labels, "M", sep="")
        }else if(k.best >= raw.best){
          at <- sort(1000 * unique(round(as.numeric(labels) / 1000, 1)))
          labels <- prettyNum(at / 1000, big.mark=",")
          labels <- paste(labels, "k", sep="")
        }else{
          labels <- prettyNum(labels, big.mark=",")
        }
      }
    }
  }
  if(paste(labels, collapse="") > 15 | length(labels) > 4){
    cex.axis <- cex.axis - 0.5/6
  }
  if(is.null(labels.at)){labels.at <- at}
  if(side %in% c(1, 3)){
    las <- 1
    title.at <- mean(par("usr")[1:2])
  }else{
    las <- 2
    title.at <- mean(par("usr")[3:4])
  }
  axis(side, at=at, labels=NA, tck=tck, line=line)
  sapply(1:length(labels.at), function(i){
    if(parse.labels){
      label <- parse(text=labels[i])
    }else if(is.numeric(labels[i])){
      label <- prettyNum(labels[i], big.mark=",")
    }else{
      label <- labels[i]
    }
    axis(side, at=labels.at[i], labels=label, tick=F, cex.axis=cex.axis,
         las=las, line=line+label.line)
  })
  if(!is.null(title)){
    axis(side, at=title.at, tick=F, labels=title, line=line+title.line, xpd=T)
  }
}


#' Convert colors to greyscale
#'
#' Convert one or more HEX-encoded color(s) to its greyscale HEX equivalent(s)
#'
#' @param in.colors
#'
#' @return character vector of HEX color
#'
#' @examples
#' hex2grey(c("#FF0000", "#00008B", "#ffffe0"))
#'
#' @export hex2grey
#' @export
hex2grey <- function(in.colors){
  mean.rgbs <- round(apply(col2rgb(in.colors), 2, mean), 0)
  sapply(mean.rgbs, function(k){rgb(k, k, k, maxColorValue=255)})
}


#' Format P-value
#'
#' Format P-value for plotting
#'
#' @param p P-value
#' @param nsmall number of digits after the decimal to retain for scientific
#' notification \[default: 2\]
#' @param max.decimal convert all P-values requiring more digits after the decimal
#' to be converted to scientific notation \[default: 3\]
#' @param equality equality symbol to print after `P` \[default: '='\]
#' @param min.neg.log10.p minimum order of magnitude to process before considering
#' P-value to be arbitrarily/meaninglessly small \[default: 30\]
#'
#' @return formatted P-value as character
#'
#' @export format.pval
#' @export
format.pval <- function(p, nsmall=2, max.decimal=3, equality="=", min.neg.log10.p=30){
  if(-log10(p) > min.neg.log10.p){
    bquote(italic(P) < 10 ^ -.(min.neg.log10.p))
  }else if(ceiling(-log10(p)) > max.decimal){
    parts <- unlist(strsplit(format(p, scientific=T), split="e"))
    base <- gsub(" ", "", format(round(as.numeric(parts[1]), nsmall), digits=1+nsmall), fixed=T)
    exp <- gsub(" ", "", as.numeric(parts[2]), fixed=T)
    if(base %in% c("1", "10")){
      if(base == "10"){
        exp <- as.character(as.numeric(exp) + 1)
      }
      bquote(italic(P) ~ .(equality) ~ 10 ^ .(exp))
    }else{
      bquote(italic(P) ~ .(equality) ~ .(base) ~ "x" ~ 10 ^ .(exp))
    }
  }else{
    bquote(italic(P) ~ .(equality) ~ .(formatC(round(p, max.decimal), digits=max.decimal)))
  }
}


#' Color points by density
#'
#' Generate colors for XY scatterplot based on point density
#'
#' @param x independent variable vector
#' @param y dependent variable vector
#' @param palette 256-color palette to be applied based on density \[default: viridis()\]
#' @param bandwidth `bandwidth` parameter passed to [densCols()]
#'
#' @details Inspired by heatscatter.R from Colby Chiang:
#'  https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R
#'
#' @return dataframe of values to be plotted with density and colors
#'
#' @seealso [viridis()], [densCols()]
#'
#' @export
color.points.by.density <- function(x, y, palette=NULL, bandwidth=1){
  # Based on heatscatter.R from Colby Chiang
  # (https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R)
  plot.df <- data.frame("x"=x, "y"=y)
  plot.df <- plot.df[which(!is.infinite(plot.df$x) & !is.infinite(plot.df$y)
                           & !is.na(plot.df$x) & !is.na(plot.df$y)), ]
  dens <- densCols(plot.df$x, plot.df$y, bandwidth=bandwidth,
                   colramp=colorRampPalette(c("black", "white")))
  plot.df$dens <- col2rgb(dens)[1, ] + 1L
  if(is.null(palette)){
    require(viridis, quietly=TRUE)
    palette <- viridis(256)
  }
  plot.df$col <- palette[plot.df$dens]
  plot.df[order(plot.df$dens), ]
}


#' Smart Element Spacing
#'
#' Enforce minimum distance between elements in a vector
#'
#' @param ideal.values Numeric values to be spaced
#' @param min.dist Minimum distance between values
#' @param lower.limit Values are not allowed to be placed below this limit
#' \[default: no limits\]
#' @param upper.limit Values are not allowed to be placed above this limit
#' \[default: no limits\]
#'
#' @return Numeric vector
#'
#' @export smart.spacing
#' @export
smart.spacing <- function(ideal.values, min.dist, lower.limit=-Inf, upper.limit=Inf){
  # Curate input values
  val.order <- order(ideal.values)
  n.vals <- length(ideal.values)
  if(n.vals * min.dist > abs(upper.limit - lower.limit)){
    warning(paste("Number of points incompatible with 'min.dist' and",
                  "specified limits. Reverting to equidistant spacing."))
    if(is.infinite(lower.limit)){
      lower.limit <- min(ideal.values)
    }
    if(is.infinite(upper.limit)){
      upper.limit <- max(ideal.values)
    }
    return(seq(lower.limit, upper.limit, length.out=n.vals))
  }
  ideal.sorted <- as.numeric(ideal.values[val.order])

  # Enforce lower & upper limits on raw values before spacing
  ideal.sorted[which(ideal.sorted < lower.limit)] <- lower.limit
  ideal.sorted[which(ideal.sorted > upper.limit)] <- upper.limit

  # Iteratively update spacing until best balance is reached
  calc.spacing <- function(ideal.sorted, sig.digits=10){
    round(ideal.sorted[-1] - ideal.sorted[-length(ideal.sorted)], sig.digits)
  }
  spacing <- calc.spacing(ideal.sorted)
  while(any(spacing < min.dist)){
    # Save old spacing info to check for convergence
    old.spacing <- spacing

    # Pick closest two points
    # (break ties by taking pair of points with greatest room to be moved)
    spacing.w.limits <- calc.spacing(c(lower.limit, ideal.sorted, upper.limit))
    room.to.move <- sapply(1:length(spacing), function(i){sum(spacing.w.limits[c(i, i+2)])})
    can.move <- which(room.to.move > 0)
    move.priority <- intersect(order(room.to.move, decreasing=TRUE), can.move)
    smallest <- min(spacing)
    smaller.idx <- head(intersect(move.priority, which(spacing == smallest)), 1)
    larger.idx <- smaller.idx + 1

    # Slide points away from each other equally s/t they are exactly min.dist apart
    # Never allow points to be placed beyond limits or jump another point in order
    smaller.max.move <- spacing.w.limits[smaller.idx]
    larger.max.move <- spacing.w.limits[larger.idx+1]
    pad.each <- (min.dist - smallest) / 2
    move.smaller <- min(pad.each, smaller.max.move)
    move.larger <- min(pad.each, larger.max.move)
    ideal.sorted[smaller.idx] <- ideal.sorted[smaller.idx] - move.smaller
    ideal.sorted[larger.idx] <- ideal.sorted[larger.idx] + move.larger

    # Update spacing from ideal.sorted
    spacing <- calc.spacing(ideal.sorted)

    # Break if spacing has converged to optimal locations
    if(setequal(spacing, old.spacing)){
      break
    }else{
      old.spacing <- spacing
    }
  }

  return(ideal.sorted[val.order])
}


#' Legend-Axis Hybrid
#'
#' Add a legend to the right Y-axis with lines connecting legend labels to
#' specified Y positions
#'
#' @param legend.names Labels to be printed in legend
#' @param x Where the legend will connect to the rest of the plot (in X-axis units)
#' @param y.positions Where should the legend labels be placed (in Y-axis units)
#' @param sep.wex Width expansion term for text relative to `x`
#' @param min.label.spacing Minimum distance between any two labels (in Y-axis units) \[default: 0.1\]
#' @param label.cex Value of `cex` to be used for legend text
#' @param lower.limit No label will be placed below this value on the Y-axis \[default: `par("usr")[3]`\]
#' @param upper.limit No label will be placed above this value on the Y-axis \[default: `par("usr")[4]`\]
#' @param colors Line colors connecting labels to plot body \[default: all black\]
#' @param lwd Width of line connecting labels to plot body \[default: 3\]
#'
#' @return NULL
#'
#' @export yaxis.legend
#' @export
yaxis.legend <- function(legend.names, x, y.positions, sep.wex,
                         min.label.spacing=0.1, label.cex=1,
                         lower.limit=NULL, upper.limit=NULL,
                         colors=NULL, lwd=3){
  if(is.null(colors)){
    colors <- "black"
  }
  if(is.null(lower.limit)){
    lower.limit <- par("usr")[3]
  }
  if(is.null(upper.limit)){
    upper.limit <- par("usr")[4]
  }
  leg.at <- smart.spacing(y.positions, min.dist=min.label.spacing,
                          lower.limit=lower.limit, upper.limit=upper.limit)
  text(x=x + sep.wex, y=leg.at, labels=legend.names, xpd=T, pos=4, cex=label.cex)
  segments(x0=x, x1=x + (1.5*sep.wex), y0=y.positions, y1=leg.at,
           lwd=lwd, col=colors, xpd=T, lend="round")
}


#' HSV-based palette interpolation
#'
#' Generate a color palette based on uniform sampling across HSV space
#'
#' @param hues Hue values in \[0, 1\]
#' @param sats Saturation values in \[0, 1\]
#' @param vals Brightness values in \[0, 1\]
#'
#' @details All of `hues`, `sats`, and `vals` must have the same lengths. If
#' any of these vectors is detected to be shorter than the others, its values
#' will be recycled as per normal R convention.
#'
#' @examples
#' set.seed(1234)
#' colors <- hsv.palette(hues=(0:10)/10, sats=runif(10, 0.5, 0.8), vals=runif(10, 0.4, 0.9))
#' plot(1:10, 1:10, pch=19, cex=3, col=colors)
#'
#' @returns Vector of hex color codes
#'
#' @export hsv.palette
#' @export
hsv.palette <- function(hues, sats, vals){
  # Ensure all vectors are the same length
  comps <- list(hues, sats, vals)
  n.colors <- max(sapply(comps, length))
  comp.df <- do.call("cbind", lapply(comps, function(v){
    l.v <- length(v)
    if(l.v < n.colors){rep(v, times=ceiling(n.colors / l.v))[1:n.colors]}else{v}
  }))
  hsv(h=comp.df[, 1], s=comp.df[, 2], v=comp.df[, 3], alpha=1)
}


#' Categorical rainbow palette
#'
#' Generate a rainbow color palette suitable for N categorical variables
#'
#' @param n Number of colors in palette
#' @param hue.range Two-element numeric vector specifying range of eligible hues \[default: c(0, 1)\]
#' @param saturation.range Two-element numeric vector specifying range of
#' eligible saturations  \[default: c(0, 1)\]
#' @param value.range  Two-element numeric vector specifying range of
#' eligible values  \[default: c(0, 1)\]
#' @param period Scalar for how many complete sine rotations (2*pi) should be sampled
#' for saturation and value \[default: 1\]
#'
#' @details Colors will be uniformly sampled across `hue.range`, with saturation
#' and value being sampled approximately according to `sin(x)` and `-sin(x)`,
#' for x uniformly spaced elements in \[0, `period`*pi\], respectively.
#'
#' @seealso [RLCtools::hsv.palette()]
#'
#' @returns Character vector of hex colors
#'
#' @export categorical.rainbow
#' @export
categorical.rainbow <- function(n, hue.range=c(0, 1), saturation.range=c(0.5, 0.8),
                                value.range=c(0.3, 0.95), period=1){
  # Ensure there's no overlapping colors at the ends of the range
  n <- n+1

  # Determine hues
  hues <- seq(hue.range[1], hue.range[2], length.out=n)

  # Generate opposing sine samples for saturation and value
  s.points <- (sin(pi*(seq(0, 2*period, length.out=n))) + 1)/2
  sats <- saturation.range[1] + (s.points*(saturation.range[2] - saturation.range[1]))
  v.points <- (-sin(pi*(seq(0, 2*period, length.out=n))) + 1)/2
  vals <- value.range[1] + (v.points*(value.range[2] - value.range[1]))

  # Generate palette
  hsv.palette(hues, sats, vals)[1:(n-1)]
}


#' Greyscale palette
#'
#' Generate a sequential greyscale palette that does not end at black & white
#'
#' @param n Number of colors in palette
#' @param oscillate Should colors be interleaved to maximize contrast?
#' \[default: `FALSE`, which provides a sequential greyscale palette\]
#'
#' @returns Character vector of hex colors
#'
#' @export greyscale.palette
#' @export
greyscale.palette <- function(n, oscillate=FALSE){
  pal <- colorRampPalette(c("black", "white"))(n + 2)[-c(1, n + 2)]
  if(oscillate){
    split.idx <- ceiling(n / 2)
    return(pal[interleave(1:split.idx, n:(split.idx+1))])
  }else{
    return(pal)
  }
}


#' Optimize label color
#'
#' Automatically determine whether a label should be white or black based on the
#' brightness of the background color upon which the label will be added
#'
#' @param bg.color Hex code for background color upon which the label will be printed
#' @param cutoff Brightness threshold for switching from white to black \[default: 0.7\]
#'
#' @returns Either "black" or "white"
#'
#' @export optimize.label.color
#' @export
optimize.label.color <- function(bg.color, cutoff=0.7){
  v <- DescTools::ColToHsv(bg.color)[3, 1]
  if(v >= cutoff){"black"}else{"white"}
}


#' Adjust color brightness
#'
#' Adjust brightness of one or more colors via conversion to HSV colorspace
#'
#' @param colors Character vector of one or more colors in hex format
#' @param d Adjustment for brightness
#'
#' @details Final brightness will be bounded within \[0, 1\]
#'
#' @export adjust.brightness
#' @export
adjust.brightness <- function(colors, d){
  require(DescTools)
  c.v <- DescTools::ColToHsv(colors)
  c.v[3, ] <- sapply(c.v[3, ] + d, function(v){max(c(0, min(c(1, v))))})
  apply(c.v, 2, function(cc){hsv(cc[1], cc[2], cc[3])})
}


#' Title case conversion
#'
#' Convert a character vector to title (or, optionally, sentence) case
#'
#' @param x Character vector to be converted
#' @param case Either `title` or `sentence`
#'
#' @returns Character vector `x` in desired case
#'
#' @export title.case
#' @export
title.case <- function(x, case="title"){
  as.character(sapply(x, function(s){
    if(case == "sentence"){
      paste(toupper(substr(s, 1, 1)),
            tolower(substr(s, 2, nchar(s))),
            sep="")
    }else{
      paste(sapply(strsplit(s, split=" "), function(ss){
        paste(toupper(substr(ss, 1, 1)),
              tolower(substr(ss, 2, nchar(ss))),
              sep="")}), collapse=" ")
    }
  }))
}
