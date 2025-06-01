#' Statistical Summary of a Data Frame
#'
#' This function produces another summary of a \code{data.frame}. This function
#' was proposed in order to apply some functions globally on a \code{data.frame}:
#' \code{quantile}, \code{median}, \code{min} and \code{max}. The usual \emph{R}
#' version cannot do so.
#'
#' @param x          numeric: matrix or \code{data.frame}
#' @param quantile   numeric: quantile of the distribution
#' @param show       logical: if \code{TRUE} prints the quantile choosen
#' @return numeric: \code{data.frame} of statistics: mean, median, quantile, standard deviation, minimum and maximum
#'
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#'
#' @seealso \code{\link{plotuScree}}, \code{\link{nScree}}, \code{\link{plotnScree}}, \code{\link{plotParallel}}

#' @export
#' @importFrom stats sd median
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## ................................................
#' ## GENERATION OF A MATRIX OF 100 OBSERVATIONS AND 10 VARIABLES
#' x   <- matrix(rnorm(1000),ncol=10)
#'
#' ## STATISTICS
#' res <- moreStats(x, quantile=0.05, show=TRUE)
#' res
#'  }
#' }
moreStats <-
function(x, quantile=0.95, show=FALSE) {
 cent  <- quantile    # The old parameter was labeled cent
 x     <- data.frame(x)
 xMean <- sapply(x, mean) # mean(x)
 xSd   <- sapply(x, stats::sd)   # sd(x)
 xMin  <- xMax <- xMedian <- xQuantile <- numeric(ncol(x))
 for (i in 1:ncol(x)) {
  xMin[i]    <- min(x[,i])
  xMax[i]    <- max(x[,i])
  xMedian[i] <- stats::median(x[,i])
  xQuantile[i]  <- quantile(x[,i],probs=cent,names=FALSE, na.rm=TRUE) # quantile(rnorm(1000),probs=cent)
  }
 names       <- colnames(x)
 results     <- rbind(mean=xMean, median=xMedian, quantile=xQuantile, sd=xSd, min=xMin, max=xMax)
 if (show==TRUE) {
  cat("------------------------ \n")
  cat("Quantile specified:", cent, "\n")
  cat("------------------------ \n")
  }
 return(results)
 }
