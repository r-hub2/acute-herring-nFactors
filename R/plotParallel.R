#' Plot a Parallel Analysis Class Object
#'
#' Plot a scree plot adding information about a parallel analysis.
#'
#' If \code{eig} is \code{FALSE} the plot shows only the parallel analysis
#' without eigenvalues.
#'
#' @param parallel numeric: vector of the results of a previous parallel
#' analysis
#' @param eig depreciated parameter: eigenvalues to analyse (not used if x is
#' used, recommended)
#' @param x numeric: a \code{vector} of eigenvalues, a \code{matrix} of
#' correlations or of covariances or a \code{data.frame} of data
#' @param model character: \code{"components"} or \code{"factors"}
#' @param main character: title of the plot
#' @param xlab character: label of the x axis
#' @param ylab character: label of the y axis
#' @param legend logical: indicator of the presence or not of a legend
#' @param ...  variable: additionnal parameters to give to the \code{cor} or
#' \code{cov} functions
#' @return Nothing returned.
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{plotuScree}}, \code{\link{nScree}},
#' \code{\link{plotnScree}}, \code{\link{parallel}}
#' @references
#' Raiche, G., Walls, T. A., Magis, D., Riopel, M. and Blais, J.-G. (2013). Non-graphical solutions
#' for Cattell's scree test. Methodology, 9(1), 23-29.
#'
#' @export
#' @importFrom graphics plot.default lines
#' @keywords Graphics
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## SIMPLE EXAMPLE OF A PARALLEL ANALYSIS
#' ## OF A CORRELATION MATRIX WITH ITS PLOT
#'  data(dFactors)
#'  eig      <- dFactors$Raiche$eigenvalues
#'  subject  <- dFactors$Raiche$nsubjects
#'  var      <- length(eig)
#'  rep      <- 100
#'  cent     <- 0.95
#'  results  <- parallel(subject,var,rep,cent)
#'
#'  results
#'
#'
#' ## PARALLEL ANALYSIS SCREE PLOT
#'  plotParallel(results, x=eig)
#'  plotParallel(results)
#'  }
#' }
#'
"plotParallel" <-
function(parallel,
         eig    = NA,
         x      = eig,
         model  = "components",
         legend = TRUE,
         ylab   = "Eigenvalues",
         xlab   = "Components",
         main   = "Parallel Analysis",
         ...
                           ) {
  if (any(!is.na(x))) eig <- eigenComputes(x, ...)
  if (!inherits(parallel, "parallel")) stop("Method is only for parallel objects")
  if (model == "factors") xlab <- "Factors"
  var        <- length(parallel$eigen$qevpea)
  if (length(eig) == 1) {
   Component <- var:1
   Location  <- seq(from = 0, to = max(parallel$eigen$qevpea)*3, length.out = var)
   graphics::plot.default(as.numeric(Component),
                as.numeric(Location),
                type = "n",
                main = main,
                xlab = xlab,
                ylab = ylab)
    }

  if (length(eig) > 1) {plotuScree(eig, main = main, xlab = xlab, ylab = ylab) }
  graphics::lines(1:var, parallel$eigen$qevpea , col = "green", type = "p", pch = 2)
  graphics::lines(1:var, parallel$eigen$mevpea,  col = "red")
  if (legend == TRUE) {
   if (length(eig) == 1) {
     leg <-  c("Mean Eigenvalues", "Centiles of the Eigenvalues")
     tco <-  c("red", "green")
     co  <-  c("red", "green")
     pc  <-  c(NA, 2)
   }
   if (length(eig) > 1) {
     leg <-  c("Eigenvalues", "Mean Eigenvalues", "Centiles of the Eigenvalues")
     tco <-  c("black", "red", "green")
     co  <-  c("black", "red", "green")
     pc  <-  c(1, NA, 2)
   }
   legend("topright",
          legend   = leg,
          text.col = tco,
          col      = co,
          pch      = pc
          )
    }
  }

