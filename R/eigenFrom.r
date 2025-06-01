#' Identify the Data Type to Obtain the Eigenvalues
#'
#' The \code{eigenFrom} function identifies the data type from which to obtain the
#' eigenvalues. The function is used internally in many functions of
#' the \pkg{nFactors} package to be able to apply these to a vector of eigenvalues,
#' a matrix of correlations or covariance or a \code{data.frame}.
#' @param x numeric: a \code{vector} of eigenvalues, a \code{matrix} of correlations or of covariances or a \code{data.frame} of data
#' @return character: return the data type to obtain the eigenvalues: \code{"eigenvalues"}, \code{"correlation"} or \code{"data"}
#'
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' \cr \cr David Magis \cr Departement de mathematiques \cr Universite de Liege
#' \cr \email{David.Magis@@ulg.ac.be}
#' @export
# #' @import methods
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' # .......................................................
#' # Different data types
#' # Examples of adequate data sources
#' # Vector of eigenvalues
#' data(dFactors)
#' x1 <- dFactors$Cliff1$eigenvalues
#' eigenFrom(x1)
#'
#' # Data from a data.frame
#' x2 <- data.frame(matrix(20*rnorm(100), ncol=5))
#' eigenFrom(x2)
#'
#' # From a covariance matrix
#' x3 <- cov(x2)
#' eigenFrom(x3)
#'
#' # From a correlation matrix
#' x4 <- cor(x2)
#' eigenFrom(x4)
#'
#' # Examples of inadequate data sources: not run because of errors generated
#' # x0 <- c(2,1)             # Error: not enough eigenvalues
#' # eigenFrom(x0)
#' # x2 <- matrix(x1, ncol=5) # Error: non a symetric covariance matrix
#' # eigenFrom(x2)
#' # eigenFrom(x3[,(1:2)])    # Error: not enough variables
#' # x6 <- table(x5)          # Error: not a valid data class
#' # eigenFrom(x6)
#' # .......................................................
#'  }
#' }
eigenFrom <-
function(x) {
 #classType <- methods::class1(x)
 classType <- data.class(x)
 res <- switch (classType,
  data.frame  = "data",
  matrix      = "correlation",
  numeric     = "eigenvalues",
  stop("Not a data.frame, a matrix, or a numeric vector")
  )

 switch (res,
  data        = if (dim(x)[2] <= 2) stop("At least 3 variables must be supplied"),
  correlation = if (dim(x)[2] <= 2) stop("At least 3 variables must be supplied"),
  eigenvalues = if (length(x) <= 2) stop("A vector of 3 eigenvalues or more must be supplied")
  )

 if (res == "correlation") if (any(x[lower.tri(x)] != t(x)[lower.tri(t(x))])) {
  stop("A correlation/covariance matrix must be symetric, empirical data must
        come from a data.frame, or eigenvalues must directly come from a vector.
        Verify the documentation about the eigenFrom function.")
  }

 invisible(res)
 }




