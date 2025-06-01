#' Create a Full Correlation/Covariance Matrix from a Matrix With Lower Part Filled and Upper Part With Zeros
#'
#' This function creates a full correlation/covariance matrix from a matrix with
#' lower part filled and upper part with zeros.
#' @param x      numeric: matrix
#' @return numeric: full correlation matrix
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{plotuScree}}, \code{\link{nScree}}, \code{\link{plotnScree}}, \code{\link{plotParallel}}
#' @export
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## ................................................
#'## LOWER CORRELATION MATRIX WITH ZEROS ON UPPER PART
#'## From Gorsuch (table 1.3.1)
#'gorsuch <- c(
#'  1,0,0,0,0,0,0,0,0,0,
#'  .6283, 1,0,0,0,0,0,0,0,0,
#'  .5631, .7353, 1,0,0,0,0,0,0,0,
#'  .8689, .7055, .8444, 1,0,0,0,0,0,0,
#'  .9030, .8626, .6890, .8874, 1,0,0,0,0,0,
#'  .6908, .9028, .9155, .8841, .8816, 1,0,0,0,0,
#' .8633, .7495, .7378, .9164, .9109, .8572, 1,0,0,0,
#'  .7694, .7902, .7872, .8857, .8835, .8884, .7872, 1,0,0,
#'  .8945, .7929, .7656, .9494, .9546, .8942, .9434, .9000, 1,0,
#'  .5615, .6850, .8153, .7004, .6583, .7720, .6201, .6141, .6378, 1)
#'
#'## UPPER CORRELATION MATRIX FILLED WITH UPPER CORRELATION MATRIX
#'gorsuch <- makeCor(gorsuch)
#'gorsuch
#'  }
#' }
makeCor <-
function(x) {
 if (is.matrix(x)) stop("x is not a vector.")
 upper       <- matrix(x,ncol=10, byrow=FALSE)
 diag(upper) <- 0
 lower <- matrix(x,ncol=10, byrow=TRUE)
 res <- lower + upper
 return(res)
 }
