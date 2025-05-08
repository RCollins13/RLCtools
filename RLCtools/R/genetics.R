#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Helper functions pertaining to genetics or genetic data


#' Convert genotype frequencies to Hardy-Weinberg coordinates
#'
#' Convert a pair of heterozgous and homozygous alt genotype frequencies to ternary
#' (x, y) coordinates for Hardy-Weinberg visualization
#'
#' @param f.het Frequency of heterozygous genotypes
#' @param f.hom Frequency of homozygous alternate genotypes
#'
#' @returns Numeric vector of (x, y) coordinates
#'
#' @export calc.hwe.xy
#' @export
calc.hwe.xy <- function(f.het, f.hom){
  if(is.na(f.het) | is.na(f.hom)){
    return(c(NA, NA))
  }
  # Assuming equilateral triangle with l(sides) = 1
  x <- f.hom + (0.5*f.het)
  y <- sin(60 * pi / 180) * f.het
  c(x, y)
}
