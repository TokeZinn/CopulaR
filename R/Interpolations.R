#¤ New Interpolation Setup ¤#
library(purrr)
library(rlang)


# Function NA - similar to NA_real_ - not necessary but convenient.
NA_function_ <- function() NA
class(NA_function_) <- c(class(NA_function_),"NA_function_")
is.na.NA_function_ = function(x){T}
is.na(NA_function_)


# Partition of the real interval [a,b]
Partition <- setClass("Partition",
                      representation(
                        partition = "numeric",
                        index = "ANY",
                        value = "ANY"
                      ),
                      prototype(
                        partition = NA_real_,
                        index = NA_function_,
                        value = NA_function_
                      ))


EquidistantPartition <- setClass("EquidistantPartition",
                   representation(
                     a = "numeric",
                     b = "numeric",
                     n = "numeric"
                   ),
                   prototype(
                     a = NA_real_,
                     b = NA_real_,
                     n = NA_real_
                   ),
                   contains = "Partition")

Mesh <- setClass("Mesh",
                 representation(
                   end_points = "list",
                   mesh = "list",
                   values = "array"
                 ),
                 prototype(
                   end_points = list(),
                   mesh = list(),
                   values = array()
                 ))


setMethod("initialize", signature(.Object = "EquidistantPartition"), function(.Object, ...){
  .Object <- callNextMethod(.Object, ...)

  if((!any(is.na(c(.Object@a, .Object@b, .Object@n)))) & is.na(.Object@partition)){
    .Object@partition <- seq(from = .Object@a, to = .Object@b, length.out = .Object@n)
  }else if((!is.na(.Object@partition)) & all(is.na(c(.Object@a, .Object@b, .Object@n)))){
      .Object@a <- .Object@partition[1]
      .Object@n <- length(.Object@partition)
      .Object@b <- .Object@partition[.Object@n]
  }else{
    stop("Provide either endpoints 'a' and 'b', and length 'n', or a partition.")
  }


  .Object@index <- approxfun(x = .Object@partition, y = seq(from = 1, to = .Object@n), method = "constant")
  .Object@value <- approxfun(x = seq(from = 1, to = .Object@n), y = .Object@partition, method = "constant")


  .Object
})
setGeneric("Index", function(...,P){standardGeneric("Index")})
setGeneric("LongIndex", function(...,P){standardGeneric("LongIndex")})
setMethod("Index", signature(P = "Partition"), function(x,P){P@index(x)})
setMethod("Index", signature(P = "Mesh"), function(x,..., P){

  Points <- list(x,...)
  Functions <- lapply(P@mesh, function(m){m@index})

  if(length(Points) != length(Functions)) stop("Number of vectors is not equal to the number of functions")

  res <- lapply(1:length(Points), function(i){Functions[[i]](Points[[i]])})

  res
})
setMethod("LongIndex", signature(P = "Mesh"), function(x,..., P){

  Points <- list(x,...)
  M <- sapply(P@mesh,function(m){m@n}) %>% lag() %>% replace_na(1) %>% cumprod

  Functions <- lapply(P@mesh, function(m){m@index})

  if(length(Points) != length(Functions)) stop("Number of vectors is not equal to the number of functions")

  res <- lapply(1:length(Points), function(i){Functions[[i]](Points[[i]])-1})
  res <- lapply(1:length(res), function(i){M[i]*res[[i]]})
  res <- reduce(res, `+`) + 1

  res
})


setGeneric("Value", function(n,.Object){standardGeneric("Value")})
setMethod("Value", signature(n = "numeric", .Object = "Partition"), function(n,.Object){.Object@value(n)})
setGeneric("Interval", function(x,.Object){standardGeneric("Interval")})
setMethod("Interval", signature(x = "numeric", .Object = "Partition"), function(x,.Object){
  indices = .Object@index(x)
  return(cbind(indices, indices + 1))})
setGeneric("MakeMesh", function(P, ...){standardGeneric("MakeMesh")})
setMethod("MakeMesh", signature(P = "Partition"), function(P,...){
  M = Mesh()

  Partitions <- list(P, ...)

  M@mesh <- Partitions
  M@end_points <- lapply(Partitions, function(p){c(p@a,p@b)})

  M
  })
setGeneric("FillValues", function(f, M,...){standardGeneric("FillValues")})
setMethod("FillValues", signature(f = "function", M = "Mesh"), function(f,M,...){
  dim = length(M@mesh)

  dimsize = sapply(1:dim, function(i){M@mesh[[i]]@n})

  Partitions = lapply(1:dim, function(i){M@mesh[[i]]@partition})
  Partitions = expand.grid(Partitions, stringsAsFactors = F)

  M@values <- array(pmap_dbl(.l = unname(Partitions),.f = f,...),dim = dimsize)

  M
})
setGeneric("LinearInterpolation", function(x, ..., M){standardGeneric("LinearInterpolation")})
setMethod("LinearInterpolation", signature(x = "numeric", M = "Mesh" ), function(x, ..., M){



})

LinInt = function(x, x_1, x_0, f_1, f_0){
  alpha = (x-x_0)/(x_1-x_0)
  f_0*(1-alpha) + f_1*alpha
}


P = EquidistantPartition(a = 0, b = 1, n = 100)
M = MakeMesh(P,P)

Upper_Index = sapply(M@mesh, function(m){m@n})

x = y = seq(from = 0, to = 1, length.out = 11)

Index(x,y,P = M)







