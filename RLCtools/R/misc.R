#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Miscellaneous helper functions


#' Stretch Vector
#'
#' Expand the length of a vector by repeating (or "stuttering") values
#'
#' @param values Vector of values to be stretched
#' @param k Number of times to duplicate each element in `values`
#'
#' @examples
#' stretch.vector(values=c(1, 2, 3), k=4)
#'
#' @export stretch.vector
#' @export
stretch.vector <- function(values, k){
  as.vector(sapply(values, rep, times=k))
}


#' Interleave two vectors
#'
#' Interleave values of two vectors
#'
#' @param v1 First vector
#' @param v2 Second vector
#'
#' @export interleave
#' @export
interleave <- function(v1, v2){
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1, v2)[order(c(ord1, ord2))]
}


#' Impute missing values
#'
#' Impute missing values into a dataframe
#'
#' @param df Dataframe to impute
#' @param fill.missing Behavior for handling missing values. See `Details`. \[default: "mean"\]
#' @param fill.columns Specify in which columns missing values should be filled.
#'
#' @details Recognized values for `fill.missing` include:
#' * `"mean"` : impute missing values as the mean of all non-missing values
#' * `"median"` : impute missing values as the median of all non-missing values
#' * `"mode"` : impute missing values as the mode of all non-missing values
#'
#' If one of the above values is not specified, the value of `fill.missing` will
#' be taken literally and will be used to fill all missing cells
#'
#' Note that `fill.missing` values of `"mean"` and `"median"` are only applicable to
#' columns with numeric values. For non-numeric columns, the behavior of
#' `fill.missing = "mode"` will always be applied.
#'
#' All columns will be subjected to automated missing value imputation unless a
#' non-`NULL` value is supplied to `fill.columns`, in which case only the column
#' names passed to `fill.columns` will be filled.
#'
#' @return data.frame
#'
#' @export impute.missing.values
#' @export
impute.missing.values <- function(df, fill.missing="mean", fill.columns=NULL){
  categorical.action <- RLCtools::mode
  fixed.fill <- FALSE
  if(fill.missing == "mean"){
    numeric.action <- mean
  }else if(fill.missing == "median"){
    numeric.action <- median
  }else if(fill.missing == "mode"){
    numeric.action <- mode
  }else{
    fixed.fill <- TRUE
  }

  if(is.null(fill.columns)){
    fill.columns <- colnames(df)[-1]
  }
  for(col in fill.columns){
    na.idxs <- which(is.na(df[, col]))
    if(length(na.idxs) > 0 & length(na.idxs) < nrow(df)){
      if(fixed.fill){
        df[na.idxs, col] <- fill.missing
      }else if(is.numeric(df[, col])){
        df[na.idxs, col] <- numeric.action(df[-na.idxs, col])
      }else{
        df[na.idxs, col] <- categorical.action(df[-na.idxs, col], break.ties="first")
      }
    }
  }

  return(df)
}


#' Keyed vector remapping
#'
#' Partial remapping of vectors that can tolerate data missingness
#'
#' @param x Vector of values (keys) to be remapped
#' @param map Named vector mapping keys (vector names) to values
#' @param default.value Optional default value to assign for elements of `x` not
#' present in `names(map)`. By default, values of `x` that fail to remap into
#' `map` will be left unchanged.
#'
#' @details Inspired by Pandas map() function:
#' https://pandas.pydata.org/docs/reference/api/pandas.Series.map.html
#'
#' @export remap
#' @export
remap <- function(x, map, default.value=NULL){
  if(!is.null(default.value)){
    x[which(!(x %in% c("NA", names(map))))] <- default.value
  }
  for(key in names(map)){
    x[which(x == key)] <- map[key]
  }
  if("NA" %in% names(map) & any(is.na(x))){
    x[which(is.na(x))] <- map["NA"]
  }
  return(x)
}

