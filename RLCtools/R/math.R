#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Mathematical helper functions


#' Chi-square power
#'
#' Calculate the power of a 2x2 Chi-square test given a known sample size and prespecified effect odds ratio
#'
#' @param N Total observations
#' @param P Expected proportion of observations with affirmative outcome
#' (e.g., disease|event) under a null baseline
#' @param OR Target odds ratio
#' @param alpha Desired type I error rate \[default: 0.05\]
#'
#' @export chisq.power.from.or
#' @export
chisq.power.from.or <- function(N, P, OR){
  Eo <- N * (1 - P)
  Ea <- N * P
  R <- OR * (Ea / Eo)
  Oa <- (R * N) / (1 + R)
  Oo <- N - Oa
  NC <- ((Oo - Eo)^2 / Eo) + ((Oa - Ea)^2 / Ea)
  Xc <- qchisq(1 - alpha, 1)
  1 - pchisq(Xc, df=1, ncp=NC)
}


#' Dynamic range
#'
#' Calculate the dynamic range of a numeric vector; that is, the ratio of the
#' maximum to minimum non-NA values in the vector
#'
#' @param values Vector of numeric values
#'
#' @returns numeric ratio of max / min
#'
#' @export dynamic.range
#' @export
dynamic.range <- function(values){
  values <- as.numeric(values)
  max(values, na.rm=T) / min(values, na.rm=T)
}


#' Inverse-variance weighted meta-analysis
#'
#' Conduct a simple inverse-variance weighted meta-analysis of point estimates
#'
#' @param estimates Numeric point estimates (e.g., log-odds ratios)
#' @param vars Variances for estimates
#' @param conf Confidence interval width \[default: 0.95\]
#'
#' @returns Named numeric vector with pooled point estimate, pooled variance,
#' lower CI bound, and upper CI bound
#'
#' @export ivw.meta
#' @export
ivw.meta <- function(estimates, vars, conf=0.95){
  keep.idx <- which(complete.cases(data.frame(estimates, vars)))
  numerator <- sum((estimates/vars)[keep.idx])
  denominator <- sum(1/vars[keep.idx])
  avg <- numerator/denominator
  pooled.var <- 1/denominator
  pooled.se <- sqrt(pooled.var)
  lower.ci <- avg + qnorm((1-conf)/2)*pooled.se
  upper.ci <- avg + qnorm(conf+(1-conf)/2)*pooled.se
  return(c("estimate" = avg,
           "variance" = pooled.var,
           "lower.ci" = lower.ci,
           "upper.ci" = upper.ci))
}


#' Position-interval intersection
#'
#' Logical check of whether one integer is inside an interval
#' defined by two other integers
#'
#' @param pos Query integer
#' @param interval Two-element numeric vector defining the bounds of the target interval
#'
#' @returns Logical
#'
#' @examples
#' is.inside(2, c(1, 3))
#' # TRUE
#'
#' is.inside(5, c(1, 3))
#' # FALSE
#'
#' @export is.inside
#' @export
is.inside <- function(pos, interval){
  pos <- as.numeric(pos)
  pos >= min(interval, na.rm=T) & pos <= max(interval, na.rm=T)
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


#' Arithmetic mode
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


#' Welch's t-test for summary statistics
#'
#' Implementation of Welch's two-sample t-test for unequal variances based on
#' summary statistics rather than observation-level data
#'
#' @param means Two-element numeric vector of sample means
#' @param sds Two-element numeric vector of sample standard deviations
#' @param ns Two-element integer vector of sample sizes
#'
#' @details Based on an implementation on StackExchange:
#' https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
#'
#' @returns data.frame with summary of test, including difference in means,
#' standard error, t statistic, and P-value
#'
#' @export welch.test.from.sumstats
#' @export
welch.test.from.sumstats <- function(means, sds, ns){
  se <- sqrt( (sds[1]^2/ns[1]) + (sds[2]^2/ns[2]) )
  df <- ( (sds[1]^2/ns[1] + sds[2]^2/ns[2])^2 )/( (sds[1]^2/ns[1])^2/(ns[1]-1) + (sds[2]^2/ns[2])^2/(ns[2]-1) )
  t <- (means[1]-means[2])/se
  dat <- c(means[1]-means[2], se, t, 2*pt(-abs(t),df))
  names(dat) <- c("Difference of means", "Std Error", "t", "P.value")
  return(dat)
}
