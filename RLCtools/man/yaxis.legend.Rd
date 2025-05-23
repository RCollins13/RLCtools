% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.helpers.R
\name{yaxis.legend}
\alias{yaxis.legend}
\title{Legend-Axis Hybrid}
\usage{
yaxis.legend(
  legend.names,
  x,
  y.positions,
  sep.wex,
  parse.labels = FALSE,
  min.label.spacing = 0.1,
  label.cex = 1,
  lower.limit = NULL,
  upper.limit = NULL,
  colors = NULL,
  lwd = 3,
  label.colors = NULL,
  label.font = 1,
  return.label.pos = FALSE,
  direct = FALSE,
  lend = "round"
)
}
\arguments{
\item{legend.names}{Labels to be printed in legend}

\item{x}{Where the legend will connect to the rest of the plot (in X-axis units)}

\item{y.positions}{Where should the legend labels be placed (in Y-axis units)}

\item{sep.wex}{How far should the labels be placed to the right of the plot
margin? (in X-axis units).}

\item{parse.labels}{Should \code{legend.names} be parsed as expressions? [default: FALSE]}

\item{min.label.spacing}{Minimum distance between any two labels (in Y-axis units) [default: 0.1]}

\item{label.cex}{Value of \code{cex} to be used for legend text}

\item{lower.limit}{No label will be placed below this value on the Y-axis
[default: \code{par("usr")[3]}]}

\item{upper.limit}{No label will be placed above this value on the Y-axis
[default: \code{par("usr")[4]}]}

\item{colors}{Line colors connecting labels to plot body [default: all black]}

\item{lwd}{Width of line connecting labels to plot body [default: 3]}

\item{label.colors}{Colors for text labels [default: all black]}

\item{label.font}{Parameter of \code{font} passed to \code{\link[=text]{text()}} [default: 1]}

\item{return.label.pos}{Should label Y positions be returned? [default: FALSE]}

\item{direct}{Draw lines directly from \code{x} and labels [default: snake connectors]}

\item{lend}{Value of \code{lend} passed to \code{\link[=segments]{segments()}} [default: "round"]}
}
\value{
NULL, unless \code{return.label.pos} is \code{TRUE}, in which case the return
will be a numeric vector of Y-position values for the legend labels
}
\description{
Add a legend to the right Y-axis with lines connecting legend labels to
specified Y positions
}
