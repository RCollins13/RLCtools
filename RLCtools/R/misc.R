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

