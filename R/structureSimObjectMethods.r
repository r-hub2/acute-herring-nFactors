#' Utility Functions for nScree Class Objects
#'
#' Utility functions for \code{structureSim} class objects. Note that with the
#' \code{plot.structureSim} a dotted black vertical line shows the median
#' number of factors retained by all the different indices.
#' @rdname structureSimObjectMethods
#'
#' @aliases boxplot.structureSim is.structureSim plot.structureSim
#' print.structureSim summary.structureSim
#' @param eigenSelect numeric: vector of the index of the selected eigenvalues
#' @param index numeric: vector of the index of the selected indices
#' @param main character: main title
#' @param nFactors numeric: if known, number of factors
#' @param object structureSim: an object of the class \code{structureSim}
#' @param vLine character: color of the vertical indicator line of the initial
#' number of factors in the eigen boxplot
#' @param x structureSim: an object of the class \code{structureSim}
#' @param xlab character: x axis label
#' @param ylab character: y axis label
#' @param ...  variable: additionnal parameters to give to the \code{boxplot},
#' \code{plot}, \code{print} and \code{summary functions.}
#' @return Generic functions for the \code{structureSim} class:
#' \item{boxplot.structureSim }{ graphic: plots an eigen boxplot }
#' \item{is.structureSim}{ logical: is the object of the class
#' \code{structureSim}? } \item{plot.structureSim }{ graphic: plots an index
#' acuracy plot} \item{print.structureSim }{ numeric: data.frame of statistics
#' about the number of components/factors to retain according to different
#' indices following a \code{structureSim} simulation}
#' \item{summary.structureSim }{ list: two data.frame, the first with the
#' details of the simulated eigenvalues, the second with the details of the
#' simulated indices}
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{nFactors-package}}
#' @references
#'
#' Raiche, G., Walls, T. A., Magis, D., Riopel, M. and Blais, J.-G. (2013). Non-graphical solutions
#' for Cattell's scree test. Methodology, 9(1), 23-29.
#'
#' @export
#' @importFrom graphics boxplot abline lines
#' @importFrom stats median
#' @keywords multivariate
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#' ## INITIALISATION
#'  library(xtable)
#'  library(nFactors)
#'  nFactors  <- 3
#'  unique    <- 0.2
#'  loadings  <- 0.5
#'  nsubjects <- 180
#'  repsim    <- 10
#'  var       <- 36
#'  pmjc      <- 12
#'  reppar    <- 10
#'  zwick     <- generateStructure(var=var, mjc=nFactors, pmjc=pmjc,
#'                                 loadings=loadings,
#'                                 unique=unique)
#'
#' ## SIMULATIONS
#' mzwick    <-  structureSim(fload=as.matrix(zwick), reppar=reppar,
#'                            repsim=repsim, details=TRUE,
#'                            N=nsubjects, quantile=0.5)
#'
#' ## TEST OF structureSim METHODS
#'  is(mzwick)
#'  summary(mzwick, index=1:5, eigenSelect=1:10, digits=3)
#'  print(mzwick, index=1:10)
#'  plot(x=mzwick, index=c(1:10), cex.axis=0.7, col="red")
#'  graphics::boxplot(x=mzwick, nFactors=3, vLine="blue", col="red")
#'   }
#'  }
#'
## .................................................................
summary.structureSim <- function(object, index=c(1:15), eigenSelect=NULL, ...) {

 if (!is.structureSim(object)) stop("Not a structureSim object")
 if (is.null(eigenSelect)) eigenSelect <- c(1:dim(object$details$eigenvalues)[2])

 cat("Report For a structureSim Class \n\n")
 NextMethod()
 cat(paste("Simulated eigenvalues","\n\n"))
 object$details$eigenvalues <- round(object$details$eigenvalues[,eigenSelect], ...)
 colnames(object$details$eigenvalues) <- paste("E",eigenSelect,sep="")
 print(object$details$eigenvalues)
 cat(paste("\n\n Number of factors retained by each index for each simulation","\n\n"))
 object$details$components <- round(object$details$components[,index], ...)
 print(object$details$components)
 }
 # summary(mzwick, index=1:5, eigenSelect=1:10, digits=2)
 # summary.structureSim(x)
 # summary(x)
## .................................................................

#' @rdname structureSimObjectMethods
#' @export
## .................................................................
print.structureSim <- function(x, index=NULL, ...) {

 if (!is.structureSim(x)) stop("Not a structureSim object")
 if (is.null(index)) index <- c(1:dim(x$nFactors)[2])

 res <- x$nFactors[,index]
 print(res, ...)
 }
 # print(mzwick, index=c(1:13), 2)
 # print.structureSim(x)
 # print(x)
## .................................................................

#' @rdname structureSimObjectMethods
#' @export
## .................................................................
boxplot.structureSim <- function(x, nFactors=NULL, eigenSelect=NULL,
                                 vLine="green", xlab="Factors",
                                 ylab="Eigenvalues", main="Eigen Box Plot", ...) {

 if (!is.structureSim(x)) stop("Not a structureSim object")
 if (is.null(eigenSelect)) eigenSelect <- c(1:dim(x$details$eigenvalues)[2])

 graphics::boxplot(x$details$eigenvalues[,eigenSelect], xlab=xlab, ylab=ylab, main=main, ...)
 graphics::abline(v=nFactors, lty=2, col=vLine)
 }
 # boxplot(mzwick, nFactors=3, eigenSelect=1:5, vLine="blue", col="red")
 # boxplot.structureSim(x)
 # boxplot(x)
## .................................................................

#' @rdname structureSimObjectMethods
#' @export
## .................................................................
plot.structureSim <- function(x, nFactors=NULL, index=NULL, main="Index Acuracy Plot", ...) {

 if (!is.structureSim(x)) stop("Not a structureSim object")
 if (is.null(index)) index <- c(1:dim(x$details$components)[2])

 if (!exists("col")  == TRUE) col  <- "black"
 ylab              <- "Average Number of Factors Retained"
 tx                <- t(x[[2]][,index])
 tx                <- data.frame(Index=rownames(tx),tx)
 colnames(tx)[2]   <- "Mean"
 tx                <- tx[order(tx[,1]),]
 plot(Mean ~ Index, type="n", data=tx, main=main, ...)
 #plot(Mean ~ Index, data=tx, cex.lab=1, cex.axis=0.7, type="n", ylab=ylab)
 graphics::abline(h=nFactors, ...)
 graphics::abline(h=stats::median(tx[2,], na.rm=TRUE), lty=2, col="black")
 for (i in 1:length(tx[,2])) graphics::lines(y=c(0,tx[i,2]), x=c(i,i), lty=2)
 }
 # plot.structureSim(x=mzwick, nFactors=3, index=c(1:10), cex.axis=0.7, col="red")
 # plot.structureSim(x)
 # plot(x)
## .................................................................

#' @rdname structureSimObjectMethods
#' @export
## .................................................................
is.structureSim <- function(object) {
 if (inherits(object, "structureSim")) return(TRUE) else return(FALSE)
 }
 # is.structureSim(mzwick)
 # is.structureSim(x)
## .................................................................

