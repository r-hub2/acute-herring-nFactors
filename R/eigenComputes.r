#' Computes Eigenvalues According to the Data Type
#'
#' The \code{eigenComputes} function computes eigenvalues from the identified data
#' type. It is used internally in many
#' fonctions of the \pkg{nFactors} package in order to apply these to a vector of
#' eigenvalues, a matrix of correlations or covariance or a data frame.
#' @param x numeric: a \code{vector} of eigenvalues, a \code{matrix} of
#' correlations or of covariances or a \code{data.frame} of data
#' @param cor logical: if \code{TRUE} computes eigenvalues from a correlation
#' matrix, else from a covariance matrix
#' @param model character: \code{"components"} or \code{"factors"}
#' @param ... variable: additionnal parameters to give to the \code{cor} or
#'  \code{cov} functions
#' @return numeric: return a vector of eigenvalues
#'
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' \cr \cr David Magis \cr Departement de mathematiques \cr Universite de Liege
#' \cr \email{David.Magis@@ulg.ac.be}
#' @export
#' @importFrom stats cor cov cov2cor
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' # .......................................................
#' # Different data types
#' # Vector of eigenvalues
#' data(dFactors)
#' x1 <- dFactors$Cliff1$eigenvalues
#' eigenComputes(x1)
#'
#' # Data from a data.frame
#' x2 <- data.frame(matrix(20*rnorm(100), ncol=5))
#' eigenComputes(x2, cor=TRUE,  use="everything")
#' eigenComputes(x2, cor=FALSE, use="everything")
#' eigenComputes(x2, cor=TRUE,  use="everything", method="spearman")
#' eigenComputes(x2, cor=TRUE,  use="everything", method="kendall")
#'
# From a covariance matrix
#' x3 <- cov(x2)
#' eigenComputes(x3, cor=TRUE,  use="everything")
#' eigenComputes(x3, cor=FALSE, use="everything")
#'
# From a correlation matrix
#' x4 <- cor(x2)
#' eigenComputes(x4, use="everything")
#' # .......................................................
#'  }
#' }
eigenComputes <-
function(x, cor=TRUE, model="components", ...) {
 dataType <- eigenFrom(x)

 if (model == "components") {
  res <- switch(dataType,
   eigenvalues = as.vector(x),
   correlation = {if (cor == FALSE) eigen(x)$values           else  eigen(stats::cov2cor(x))$values},
   data        = {if (cor == TRUE)  eigen(stats::cor(x, ...))$values else  eigen(stats::cov(x, ...))$values}
   )
  }

 if (model == "factors") {
  res <- switch(dataType,
   eigenvalues = as.vector(x),
   correlation = {if (cor == FALSE) eigen(corFA(x, method="ginv"))$values else   eigen(stats::cov2cor(corFA(x, method="ginv")))$values},
   data        = {if (cor == TRUE)  eigen(corFA(stats::cor(x, ...), method="ginv"))$values else  eigen(corFA(stats::cov(x, ...), method="ginv"))$values}
   )
  }
 return(res)
 }
