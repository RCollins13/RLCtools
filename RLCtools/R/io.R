#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Input/output functions for common file types


#' Load gene list(s)
#'
#' Load one or more list(s) of gene symbols
#'
#' @param input.tsv Path to two-column .tsv of gene list prefix(es) and path to gene list(s)
#'
#' @returns list of gene symbol character vectors
#'
#' @details Will deduplicate the list before loading if there are any strictly
#' identical rows in `input.tsv`
#'
#' @export load.gene.lists
#' @export
load.gene.lists <- function(input.tsv){
  info.df <- read.table(input.tsv, header=F, sep="\t")
  info.df <- info.df[!duplicated(info.df), ]
  set.names <- info.df[, 1]
  sets <- lapply(info.df[, 2], function(path){sort(unique(read.table(path, header=F)[, 1]))})
  names(sets) <- set.names
  return(sets)
}
