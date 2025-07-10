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


#' Simulate genotypes
#'
#' Simulate a vector of biallelic genotypes according to Hardy-Weinberg Equilibrium
#'
#' @param af Minor allele frequency
#' @param n Number of genotypes to simulate
#' @param seed Optional random seed
#' @param return.dosage Return integer allele dosages \[default: FALSE\]
#'
#' @returns Either a character vector of simulated VCF-style genotypes or a
#' numeric vector of allele dosages, depending on `return.dosage`
#'
#' @export simulate.gts
#' @export
simulate.gts <- function(af, n, seed=NULL, return.dosage=FALSE){
  # Convert AF into classical HWE terms
  af <- as.numeric(af)
  p <- af
  q <- 1 - af

  # Genotype probabilities
  p.hom <- p^2
  p.het <- 2*p*q
  p.ref <- q^2

  # Sample genotypes
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(return.dosage){
    sample(0:2, n, replace=TRUE, prob=c(p.ref, p.het, p.hom))
  }else{
    sample(c("0/0", "0/1", "1/1"), n, replace=TRUE, prob=c(p.ref, p.het, p.hom))
  }
}

