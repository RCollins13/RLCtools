% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{clean.numeric.labels}
\alias{clean.numeric.labels}
\title{Clean labels for numeric values}
\usage{
clean.numeric.labels(
  vals,
  suffix.delim = "",
  acceptable.decimals = 0,
  min.label.length = 1,
  return.rounded.vals = FALSE
)
}
\arguments{
\item{vals}{Vector of numeric values to translate to labels}

\item{suffix.delim}{Character delimiter (i.e., separator) between value and
suffix [default: no delimiter, ""]}

\item{acceptable.decimals}{Number of significant decimal digits to permit
before jumping to a greater log suffix [default: 0]}

\item{min.label.length}{Minimum total count of digits to include in each
label; labels shorter than this length will have their decimal places expanded
[default: 1]}

\item{return.rounded.vals}{Should the exact rounded numeric values also
be returned? [default: FALSE]}
}
\description{
Convert large numeric values to simple, human-readable labels with log suffixes
}
\details{
Note that all labels will be converted to the same scale / suffix scheme

If element-wise conversion is desired, it's recommended to apply over \code{vals} instead
}
\examples{

vals <- 10^(0:6)

# Standard invocation:
clean.numeric.labels(vals)

# Element-wise implementation
sapply(vals, clean.numeric.labels)

}
