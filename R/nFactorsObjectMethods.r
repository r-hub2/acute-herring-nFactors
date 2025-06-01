#' Utility Functions for nFactors Class Objects
#'
#' Utility functions for \code{nFactors} class objects.
#'
#'
# #' @aliases is.nFactors print.nFactors summary.nFactors

#' @rdname nFactorsObjectMethods
#'
#' @param x nFactors: an object of the class nFactors
#' @param ...  variable: additionnal parameters to give to the \code{print}
#' function with \code{print.nFactors} or to the \code{summary} function with
#' \code{summary.nFactors}
#' @return Generic functions for the nFactors class:
#'
#' \item{is.nFactors}{ logical: is the object of the class nFactors? }
#' \item{print.nFactors }{ numeric: vector of the number of components/factors
#' to retain: same as the \code{nFactors} vector from the \code{nFactors}
#' object} \item{summary.nFactors }{ data.frame: details of the results from a
#' nFactors object: same as the \code{details} data.frame from the
#' \code{nFactors} object, but with easier control of the number of decimals
#' with the \code{digits} parameter}
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{nBentler}}, \code{\link{nBartlett}},
#' \code{\link{nCng}}, \code{\link{nMreg}}, \code{\link{nSeScree}}
#' @references
#' Raiche, G., Walls, T. A., Magis, D., Riopel, M. and Blais, J.-G. (2013). Non-graphical solutions
#' for Cattell's scree test. Methodology, 9(1), 23-29.
#'
#' @export
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## SIMPLE EXAMPLE
#'  data(dFactors)
#'  eig      <- dFactors$Raiche$eigenvalues
#'  N        <- dFactors$Raiche$nsubjects
#'
#'  res <- nBartlett(eig,N); res; is.nFactors(res); summary(res, digits=2)
#'  res <- nBentler(eig,N);  res; is.nFactors(res); summary(res, digits=2)
#'  res <- nCng(eig);        res; is.nFactors(res); summary(res, digits=2)
#'  res <- nMreg(eig);       res; is.nFactors(res); summary(res, digits=2)
#'  res <- nSeScree(eig);    res; is.nFactors(res); summary(res, digits=2)
#'
#' ## SIMILAR RESULTS, BUT NOT A nFactors OBJECT
#'  res <- nScree(eig);      res; is.nFactors(res); summary(res, digits=2)
#'
## .................................................................
#'  }
#' }
is.nFactors <-
function(x) {
 if (any(class(x) == "nFactors")) return(TRUE) else return(FALSE)
 }
## .................................................................

## .................................................................
#' @rdname nFactorsObjectMethods
#' @export
print.nFactors <-
function(x, ...) {
 if (!is.nFactors(x)) stop("Not a nFactors object")
 res <- x$nFactors
 print(res, ...)
 }
## .................................................................

## .................................................................
#' @rdname nFactorsObjectMethods
#' @param object nFactors: an object of the class nFactors
#' @export
summary.nFactors <-
function(object, ...) {
 if (!is.nFactors(object)) stop("Not a nFactors object")
 cat("Report For a nFactors Class \n\n")
 NextMethod()
 cat(paste("Details:","\n\n"))
 print(object$detail, ...)
 cat(paste("\n\n Number of factors retained by index","\n\n"))
 print(object$nFactors)
 }
## .................................................................

