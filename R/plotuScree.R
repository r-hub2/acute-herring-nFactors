#' Plot of the Usual Cattell's Scree Test
#'
#' \code{uScree} plot a usual scree test of the eigenvalues of a correlation
#' matrix.
#'
#'
#' @param Eigenvalue depreciated parameter: eigenvalues to analyse (not used if
#' x is used, recommended)
#' @param x numeric: a \code{vector} of eigenvalues, a \code{matrix} of
#' correlations or of covariances or a \code{data.frame} of data
#' @param model character: \code{"components"} or \code{"factors"}
#' @param main character: title of the plot (default is \code{Scree Plot})
#' @param xlab character: label of the x axis (default is \code{Component})
#' @param ylab character: label of the y axis (default is \code{Eigenvalue})
#' @param ...  variable: additionnal parameters to give to the
#' \code{eigenComputes} function
#' @return Nothing returned with this function.
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{nScree}}, \code{\link{parallel}}
#' @references Cattell, R. B. (1966). The scree test for the number of factors.
#' \emph{Multivariate Behavioral Research, 1}, 245-276.
#' @export
#' @importFrom graphics lines par text plot.default
#' @importFrom stats cor
#' @keywords Graphics
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## SCREE PLOT
#'  data(dFactors)
#'  attach(dFactors)
#'  eig = Cliff1$eigenvalues
#'  plotuScree(x=eig)
#'  }
#' }
"plotuScree" <-
function(Eigenvalue, x=Eigenvalue, model  = "components",
         ylab   = "Eigenvalues",
         xlab   = "Components",
         main   = "Scree Plot" ,
         ...) {
 Eigenvalue  <- eigenComputes(x, ...)
 if (!inherits(Eigenvalue, "numeric")) stop("use only with \"numeric\" objects")
 if (model == "factors") xlab <- "Factors"
 graphics::par(mfrow = c(1,1))
 nk          <- length(Eigenvalue)
 Component   <- 1:nk
 graphics::plot.default(as.numeric(Component),
              as.numeric(Eigenvalue),
              type = 'b',col = "black", pch = 1,
              ylab = ylab,
              xlab = xlab,
              main = main
              )
 }

