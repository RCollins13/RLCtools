#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Mathematical helper functions


#' Arithmetic Mode
#'
#' Compute the arithmetic mode of a vector
#'
#' @param values Vector of values
#' @param break.ties Specify how ties should be broken. See `Details`. \[default: "all"\]
#' @param seed Random seed to break ties, if needed. \[default: 2022\]
#'
#' @details Recognized values for `break.ties` include:
#' * `all` : return all tied modes
#' * `first` : return the first mode value encountered
#' * `last` : return the last mode value encountered
#' * `random` : return a randomly selected mode value
#'
#' @return mode(s) of `values`
#'
#' @export mode
#' @export
mode <- function(values, break.ties="all", seed=2022){
  vt <- sort(table(values), decreasing=TRUE)
  vmax <- max(vt, na.rm=TRUE)
  modes <- names(vt)[which(vt == vmax)]
  if(break.ties == "all"){
    return(modes)
  }else if(break.ties == "first"){
    return(head(modes, 1))
  }else if(break.ties == "last"){
    return(tail(modes, 1))
  }else if(break.ties == "random"){
    set.seed(seed)
    return(sample(modes, 1))
  }
}


#' Define Outlier Values
#'
#' Define outlier observations from a numeric vector using the standard
#' Q1/Q3 + N*IQR approach
#'
#' @param values Vector of numeric values
#' @param n.iqr Multiple of IQRs beyond which a point will be treated
#' as an outlier \[default: 3\]
#' @param fixed.min Optional. If specified, will override the lower bound below
#' which observations will be treated as outliers.
#' @param fixed.min Optional. If specified, will override the upper bound below
#' which observations will be treated as outliers.
#'
#' @return Logical vector of outlier assignments for each element in `values`
#'
#' @export label.outliers
#' @export
label.outliers <- function(values, n.iqr=3, fixed.min=NULL, fixed.max=NULL){
  values <- as.numeric(values)
  v.iqr <- IQR(values, na.rm=T)
  v.q1q3 <- quantile(values, probs=c(0.25, 0.75), na.rm=T)
  v.out.bounds <- v.q1q3 + c(n.iqr * c(-1, 1) * v.iqr)
  if(!is.null(fixed.min)){
    v.out.bounds[1] <- as.numeric(fixed.min)
  }
  if(!is.null(fixed.max)){
    v.out.bounds[2] <- as.numeric(fixed.max)
  }
  (values < v.out.bounds[1] | values > v.out.bounds[2])
}


