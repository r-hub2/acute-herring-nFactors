#' Utility Functions for nScree Class Objects
#'
#' Utility functions for \code{nScree} class objects. Some of these functions
#' are already implemented in the \code{nFactors} package, but are easier to
#' use with generic functions like these.
#'
#' @rdname nScreeObjectMethods
#' @param object nScree: an object of the class \code{nScree}
#' @param ...  variable: additionnal parameters to give to the \code{print}
#' function with \code{print.nScree}, the \code{plotnScree} with
#' \code{plot.nScree} or to the \code{summary} function with
#' \code{summary.nScree}
#'
#' @return Generic functions for the nScree class:
#' \item{is.nScree}{ logical: is the object of the class \code{nScree}? }
#' \item{plot.nScree }{ graphic: plots a figure according to the
#' \code{plotnScree} function}
#' \item{print.nScree }{ numeric: vector of the
#' number of components/factors to retain: same as the \code{Components} vector
#' from the \code{nScree} object}
#' \item{summary.nScree }{ data.frame: details
#' of the results from a nScree analysis: same as the \code{Analysis}
#' data.frame from the \code{nScree} object, but with easier control of the
#' number of decimals with the \code{digits} parameter}
#'
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
# #' @family nScree
# #' @seealso \code{\link{plotuScree}}, \code{\link{plotnScree}},
# #' \code{\link{parallel}}, \code{\link{plotParallel}},
#' @references
#' Raiche, G., Walls, T. A., Magis, D., Riopel, M. and Blais, J.-G. (2013). Non-graphical solutions
#' for Cattell's scree test. Methodology, 9(1), 23-29.
#' @export
#' @importFrom stats coef
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## INITIALISATION
#'  data(dFactors)                      # Load the nFactors dataset
#'  attach(dFactors)
#'  vect         <- Raiche              # Use the example from Raiche
#'  eigenvalues  <- vect$eigenvalues    # Extract the observed eigenvalues
#'  nsubjects    <- vect$nsubjects      # Extract the number of subjects
#'  variables    <- length(eigenvalues) # Compute the number of variables
#'  rep          <- 100                 # Number of replications for the parallel analysis
#'  cent         <- 0.95                # Centile value of the parallel analysis
#'
#' ## PARALLEL ANALYSIS (qevpea for the centile criterion, mevpea for the mean criterion)
#'  aparallel    <- parallel(var     = variables,
#'                           subject = nsubjects,
#'                           rep     = rep,
#'                           cent    = cent
#'                           )$eigen$qevpea  # The 95 centile
#'
#' ## NOMBER OF FACTORS RETAINED ACCORDING TO DIFFERENT RULES
#'  results      <- nScree(x=eigenvalues, aparallel=aparallel)
#'
#'  is.nScree(results)
#'  results
#'  summary(results)
#'
#' ## PLOT ACCORDING TO THE nScree CLASS
#'  plot(results)
#'  }
#' }


summary.nScree <-
function(object, ...) {
 if (!is.nScree(object)) stop("Not a nScree object")
 cat("Report For a nScree Class \n\n")
 #digits <- 2
 NextMethod()
 cat(paste("Details:",object$Model,"\n\n"))
 object$Analysis[,c(1:5,7)] <- round(object$Analysis[,c(1:5,7)], ...)
 print(object[[2]])
 cat(paste("\n\n Number of factors retained by index","\n\n"))
 print(object[[1]])
 }
## .................................................................

#' @rdname nScreeObjectMethods
#' @param x Results of a previous \code{nScree} analysis
# #' @family nScree
#' @export
# #' @method print print.nScree
# #' @S3method print print.nScree
## .................................................................
print.nScree <-
function(x, ...) {
 res <- x[[1]]
 print(res, ...)
 }
## .................................................................

#' @rdname nScreeObjectMethods
# #' @param x Results of a previous \code{nScree} analysis
# #' @family nScree
#' @export
# #' @method plot plot.nScree
# #' @S3method plot plot.nScree
## .................................................................
plot.nScree <-
function(x, ...) {
 plotnScree(x, ...)
 }
## .................................................................

#' @rdname nScreeObjectMethods
# #' @family nScree
#' @export
# #' @method is is.nScree
# #' @S3method is is.nScree
## .................................................................
is.nScree <-
function(object) {
 if (inherits(object, "nScree")) return(TRUE) else return(FALSE)
 }
## .................................................................
