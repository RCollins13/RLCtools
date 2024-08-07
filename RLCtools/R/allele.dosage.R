#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2024-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Functions for manipulating allele dosage matrixes


#' Query allele dosage matrix
#'
#' Subsets and optionally compresses an allele dosage matrix
#'
#' @param ad Either a path to an allele dosage .bed.gz matrix or a data.frame
#' from a previous call of [query.ad.matrix]
#' @param query.regions List of tripartite coordinate vectors to query; see `Details`
#' @param query.ids Vector of variant IDs to query
#' @param action Action to apply to query rows; see `Details`
#' @param weights Optional named vector of weights for each variant. Variants
#' not present in this vector will be assigned weight = 1. \[default: no weighting\]
#' @param na.behavior Specify how `NA` entries in `ad.df` should be treated when
#' compressing. See [compress.ad.matrix] for more information.
#' @param na.frac Fraction of `NA` entries allowed before failing a sample.
#' See [compress.ad.matrix] for more information.
#' @param keep.samples Optional. Either character vector of sample IDs to retain
#' or path to single-column flat text file of sample IDs to retain.
#' \[default: keep all samples\]
#'
#' @details The value for `query.regions` is expected to be a list of vectors,
#' where each vector must either be a single chromosome or have exactly
#' three elements in the format (chromosome, start, end).
#'
#' Recognized values for `action` include:
#' * `"verbose"` : return the full query matrix \[default\]
#' * `"any"` : return numeric indicator if any of the query rows are non-zero
#' for each sample
#' * `"all"` : return numeric indicator if all of the query rows are non-zero
#' for each sample
#' * `"count"` : return the count of non-reference variants per sample
#' irrespective of genotype
#' * `"sum"` : return the sum of allele dosages for all query rows per sample
#' * `"max"` : return the maximum value for all query rows per sample
#' * `"mean"` : return the mean of allele dosages for all query rows per sample
#'
#' @return numeric vector or data.frame, depending on `action`
#'
#' @importFrom bedr tabix
#'
#' @seealso [compress.ad.matrix]
#'
#' @export query.ad.matrix
#' @export
query.ad.matrix <- function(ad, query.regions=NULL, query.ids=NULL,
                            action="verbose", weights=NULL,
                            na.behavior="threshold", na.frac=0.05,
                            keep.samples=NULL){

  # Load region(s) of interest or whole matrix if query.regions is NULL
  if(inherits(ad, "character")){
    if(is.null(query.regions)){
      ad <- read.table(ad, header=T, comment.char="", check.names=F,
                       stringsAsFactors=F)
      ad.header <- NULL
    }else{
      require(bedr, quietly=TRUE)
      tabix.query <- sapply(query.regions, function(coords){
        if(length(coords) == 1){
          coords <- c(coords, 0, 400000000)
        }
        coords <- as.character(coords)
        paste(coords[1], ":", coords[2], "-", coords[3], sep="")
      })
      tabix.query <- bedr.merge.region(bedr.sort.region(tabix.query, verbose=FALSE),
                                       verbose=FALSE)
      incon <- gzfile(ad, open="r")
      ad.header <- unlist(strsplit(readLines(incon, 1), split="\t"))
      close(incon)
      ad <- as.data.frame(bedr::tabix(region=tabix.query, file.name=ad),
                          check.names=F)
    }
    if(!is.null(query.ids)){
      ad <- ad[which(ad[, 4] %in% query.ids), ]
    }
    ad <- ad[!duplicated(ad[, 1:4]), ]
    if(!is.null(ad.header)){
      colnames(ad) <- ad.header
    }
    rownames(ad) <- ad[, 4]
    ad <- ad[, -c(1:4)]
  }else{
    # Subset rows based on query IDs
    # This is necessary for use cases when ad is passed as a data.frame
    if(!is.null(query.ids)){
      ad <- ad[query.ids, ]
    }
  }

  # Subset columns based on sample IDs
  if(!is.null(keep.samples)){
    if(length(keep.samples) == 1){
      if(file.exists(keep.samples)){
        keep.samples <- unique(read.table(keep.samples, header=F)[, 1])
      }
    }
    ad <- ad[, keep.samples]
  }

  # Ensure all values are numeric
  ad <- as.data.frame(do.call("cbind", lapply(ad, as.numeric)),
                      row.names=rownames(ad), check.names=F)

  # Lastly, apply various compression strategies prior to returning
  compress.ad.matrix(ad[which(!is.na(colnames(ad)))],
                     action, weights, na.behavior, na.frac)
}


#' Compress allele dosage matrix
#'
#' Compresses information in an allele dosage matrix
#'
#' @param ad.df Data.frame of allele dosages loaded by [query.ad.matrix]
#' @param action Action to apply to query rows; see `Details` in [query.ad.matrix]
#' @param weights Optional named vector of weights for each variant. Variants
#' not present in this vector will be assigned weight = 1. \[default: no weighting\]
#' @param na.behavior Specify how `NA` entries in `ad.df` should be treated when
#' compressing. See `Details`.
#' @param na.frac Fraction of `NA` entries allowed before failing a sample. See `Details`.
#' @param keep.vids Optional vector of variant IDs to keep \[default: keep all rows in `ad.df`\]
#' @param elig.controls Sample IDs eligible to be reported as 0 \[default: all samples\]
#'
#' @return numeric or logical vector, depending on `action`
#'
#' @details The `na.behavior` argument recognizes two options:
#' * `all`: only return `NA` for a sample if none of the rows in `ad.df` are non-NA
#' * `any`: return `NA` for a sample if any of the rows in `ad.df` are NA
#' * `threshold`: return `NA` for a sample if more than `na.frac` genotypes in that
#' sample are NA
#'
#' @seealso [query.ad.matrix]
#'
#' @export compress.ad.matrix
#' @export
compress.ad.matrix <- function(ad.df, action, weights=NULL,
                               na.behavior="threshold", na.frac=0.05,
                               keep.vids=NULL, elig.controls=NULL){
  if(nrow(ad.df) < 2){
    if(nrow(ad.df) == 0 | action == "verbose"){
      return(ad.df)
    }else{
      ad.vals <- as.numeric(ad.df[1, ])
      names(ad.vals) <- colnames(ad.df)
      return(ad.vals)
    }
  }
  if(!is.null(keep.vids)){
    ad.df <- ad.df[intersect(rownames(ad.df), keep.vids), ]
  }
  if(na.behavior == "all"){
    na.fx <- all
  }else if(na.behavior == "any"){
    na.fx <- any
  }else if(na.behavior == "threshold"){
    na.fx <- function(tf){length(which(tf)) / length(tf) > na.frac}
  }else{
    stop(paste("compress.ad.matrix does not recognize", na.behavior,
               "as a valid option for 'na.behavior'."))
  }
  col.na <- apply(ad.df, 2, function(vals){
    na.fx(is.na(as.numeric(vals)))
  })

  # Reorder weights, fill missing weights, and apply weights to all ACs
  ordered.weights <- rep(1, times=nrow(ad.df))
  if(!is.null(weights)){
    for(svid in names(weights)){
      hits <- which(rownames(ad.df) == svid)
      if(length(hits) > 0){
        ordered.weights[hits] <- weights[svid]
      }
    }
  }
  ad.df <- as.data.frame(apply(ad.df, 2, function(vals){vals * ordered.weights}),
                         check.names=F)

  if(action == "verbose"){
    query.res <- ad.df
  }else if(action == "any"){
    query.res <- as.numeric(apply(ad.df, 2, function(vals){
      any(as.logical(as.numeric(vals)), na.rm=T)
    }))
  }else if(action == "all"){
    query.res <- as.numeric(apply(ad.df, 2, function(vals){
      all(as.logical(as.numeric(vals)), na.rm=T)
    }))
  }else if(action == "count"){
    query.res <- apply(ad.df, 2, function(vals){
      sum(as.numeric(vals != 0), na.rm=T)
    })
  }else if(action == "sum"){
    query.res <- apply(ad.df, 2, function(vals){
      sum(as.numeric(vals), na.rm=T)
    })
  }else if(action == "max"){
    query.res <- apply(ad.df, 2, function(vals){
      max(as.numeric(vals), na.rm=T)
    })
  }else if(action == "mean"){
    query.res <- apply(ad.df, 2, function(vals){
      mean(as.numeric(vals), na.rm=T)
    })
  }else{
    stop(paste("action ", action, " not recognized as a viable option for compress.ad.matrix", sep="'"))
  }

  if(length(col.na) > 0 & action != "verbose"){
    query.res[col.na] <- NA
  }
  names(query.res) <- colnames(ad.df)
  if(!is.null(elig.controls)){
    query.res[which(query.res == 0 & !(names(query.res) %in% elig.controls))] <- NA
  }
  return(query.res[which(!is.na(colnames(ad.df)))])
}

