#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Statistical and helper functions for association testing


#' Z-score saddlepoint approximation
#'
#' Apply saddlepoint approximation to vector of Z-scores to generate adjusted P-values
#'
#' @param zscores numeric vector of Z-scores
#' @param alternative alternative hypothesis. Options include "two.sided", "less", or "greater". \[default: "two.sided"\]
#'
#' @return data.frame with two columns:
#' * `$zscores` for corrected Z-scores
#' * `$pvalues` for P-values corresponding to corrected Z-scores
#'
#' @export
saddlepoint.adj <- function(zscores, alternative="two.sided"){
  require(EQL)
  zscores.orig <- zscores
  mu.hat <- mean(zscores, na.rm=T)
  sd.hat <- sd(zscores, na.rm=T)
  cumuls <- gaussianCumulants(mu.hat, sd.hat)
  dx <- 0.01
  x <- seq(-40, 40, dx)
  saddle.pdf <- saddlepoint(x, 1, cumuls)$approx
  # Dev note: must infer parameters of saddlepoint-approximated normal for precise extreme P-values with pnorm()
  mu.saddle <- sum(x * saddle.pdf) * dx
  sd.saddle <- sqrt(sum(saddle.pdf * dx * (x - mu.saddle)^2))

  # Compute new Z-scores and P-values
  new.zscores <- (zscores.orig - mu.saddle) / sd.saddle
  if(alternative == "two.sided"){
    new.pvals <- 2*mapply(pnorm, new.zscores, lower.tail=new.zscores<=0)
  }else{
    new.pvals <- pnorm(new.zscores, lower.tail=(alternative=="less"))
  }
  return(data.frame("zscores" = new.zscores, "pvalues" = new.pvals))
}
