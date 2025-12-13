#!/usr/bin/env R

##################
#    RLCtools    #
##################

# Copyright (c) 2025-Present Ryan L. Collins and the Dana-Farber Cancer Institute
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>


# Machine learning-related helper functions


#' Train an elastic net regression
#'
#' Optimize an elastic net regression with k-fold nested cross-validation
#'
#' @param X Matrix of features
#' @param Y Vector of numeric targets
#' @param k Number of cross-validation folds \[default: 5\]
#' @param seed Random seed to use for cross-validation fold assignment \[default: 2025\]
#' @param fold.indexes Optional list of numeric indexes corresponding to
#' cross-validation folds. If provided, will override `k` and `seed`.
#' \[default: use random CV\]
#' @param inner.k How many inner cross-validations should be performed within
#' each fold? \[default: same as `k`\]
#' @param tune.length Resolution of grid search to optimize elastic net params \[default: 25\]
#'
#' @seealso [caret::trainControl]
#'
#' @return A fitted [caret::train] list object
#'
#' @export train.elastic.net.cv
#' @export
train.elastic.net.cv <- function(X, Y, k=5, seed=2025, fold.indexes=NULL,
                                 inner.k=k, tune.length=25){
  require(glmnet, quietly=TRUE)
  require(caret, quietly=TRUE)

  # Define fold membership if not specified by user
  n.obs <- nrow(X)
  n.labels <- length(Y)
  if(n.obs != n.labels){
    stop("Number of observations does not match number of labels in optimize.elnet.params.cv()")
  }
  if(is.null(fold.indexes)){
    set.seed(seed)
    fold.indexes <- createFolds(1:nrow(X), k=k)
  }

  # Set training parameters
  control <- trainControl(method = "repeatedcv",
                          number = k,
                          repeats = inner.k,
                          search = "random")

  # Train elastic net
  train.df <- cbind(X, Y)
  colnames(train.df)[ncol(train.df)] <- "Y"
  set.seed(seed)
  train(Y ~ .,
        data = train.df,
        method = "glmnet",
        preProcess = c("center", "scale"),
        tuneLength = tune.length,
        trControl = control)
}


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
