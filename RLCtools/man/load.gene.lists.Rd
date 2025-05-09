% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{load.gene.lists}
\alias{load.gene.lists}
\title{Load gene list(s)}
\usage{
load.gene.lists(input.tsv)
}
\arguments{
\item{input.tsv}{Path to two-column .tsv of gene list prefix(es) and path to gene list(s)}
}
\value{
list of gene symbol character vectors
}
\description{
Load one or more list(s) of gene symbols
}
\details{
Will deduplicate the list before loading if there are any strictly
identical rows in \code{input.tsv}
}
