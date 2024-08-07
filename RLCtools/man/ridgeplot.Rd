% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.basic.R
\name{ridgeplot}
\alias{ridgeplot}
\title{Ridgeplot}
\usage{
ridgeplot(
  data,
  names = NULL,
  hill.overlap = 0.35,
  xlims = NULL,
  x.axis = TRUE,
  fill = NULL,
  border = NULL,
  border.lwd = 2,
  parmar = c(2.5, 3, 0.25, 0.25)
)
}
\arguments{
\item{data}{list of \link{density} objects to plot}

\item{names}{optional list of names for Y axis [default: take names from data]}

\item{hill.overlap}{relative fraction of overlap bewtween adjacent hills [default: 0.35]}

\item{xlims}{custom X axis limits}

\item{fill}{vector of polygon fill colors [default: "grey70"]}

\item{border}{vector of hill border colors [default: "grey35"]}

\item{border.lwd}{line width for hill borders [default: 2]}

\item{parmar}{vector of values passed to par(mar)}

\item{ylims}{custom Y axis limits}
}
\description{
Generate a ridgeplot using base R syntax
}
