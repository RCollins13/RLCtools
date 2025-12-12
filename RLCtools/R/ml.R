#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2025-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Machine learning-related helper functions


#' ReLU activation function
#'
#' Applies the Rectified Linear Unit (ReLU) activation function to a numeric input vector
#'
#' @param x Numeric values to be processed
#'
#' @export relu
#' @export
relu <- function(x){
  sapply(as.numeric(x), function(x){max(c(0, x))})
}
