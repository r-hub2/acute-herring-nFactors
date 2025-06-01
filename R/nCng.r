#' Cattell-Nelson-Gorsuch CNG Indices
#'
#' This function computes the \emph{CNG} indices for the eigenvalues of a
#' correlation/covariance matrix (Gorsuch and Nelson, 1981; Nasser, 2002, p.
#' 400; Zoski and Jurs, 1993, p. 6).
#'
#' Note that the \code{nCng} function is only valid when more than six
#' eigenvalues are used and that these are obtained in the context of a
#' principal component analysis. For a factor analysis, some eigenvalues could
#' be negative and the function will stop and give an error message.
#'
#' The slope of all possible sets of three adjacent eigenvalues are compared,
#' so \emph{CNG} indices can be applied only when more than six eigenvalues are
#' used. The eigenvalue at which the greatest difference between two successive
#' slopes occurs is the indicator of the number of components/factors to
#' retain.
#'
#' @param x numeric: a \code{vector} of eigenvalues, a \code{matrix} of
#' correlations or of covariances or a \code{data.frame} of data
#' @param cor logical: if \code{TRUE} computes eigenvalues from a correlation
#' matrix, else from a covariance matrix
#' @param model character: \code{"components"} or \code{"factors"}
#' @param details logical: if \code{TRUE} also returns detains about the
#' computation for each eigenvalue.
#' @param ...  variable: additionnal parameters to give to the
#' \code{eigenComputes} function
#' @return \item{nFactors}{ numeric: number of factors retained by the CNG
#' procedure. } \item{details}{ numeric: matrix of the details for each index.}
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{plotuScree}}, \code{\link{nScree}},
#' \code{\link{plotnScree}}, \code{\link{plotParallel}}
#' @references Gorsuch, R. L. and Nelson, J. (1981). \emph{CNG scree test: an
#' objective procedure for determining the number of factors}. Presented at the
#' annual meeting of the Society for multivariate experimental psychology.
#'
#' Nasser, F. (2002). The performance of regression-based variations of the
#' visual scree for determining the number of common factors. \emph{Educational
#' and Psychological Measurement, 62(3)}, 397-419.
#'
#' Zoski, K. and Jurs, S. (1993). Using multiple regression to determine the
#' number of factors to retain in factor analysis. \emph{Multiple Linear
#' Regression Viewpoints, 20}(1), 5-9.
#' @export
#' @importFrom stats lm
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## SIMPLE EXAMPLE OF A CNG ANALYSIS
#'
#'  data(dFactors)
#'  eig      <- dFactors$Raiche$eigenvalues
#'
#'  results  <- nCng(eig, details=TRUE)
#'  results
#'
#'  plotuScree(eig, main=paste(results$nFactors,
#'                             " factors retained by the CNG procedure",
#'                             sep=""))
#'  }
#' }
nCng <-
function(x, cor=TRUE, model="components", details=TRUE, ...) {
 x       <- eigenComputes(x, cor=cor, model=model, ...)
 detail  <- NULL
 nlength <- 2
 n       <- length(x)
 if (n < 6) stop("The number of variables must be at least 6.")
 i       <- 1
 cng     <- numeric(n-5)
 while ((i+2*nlength+1) <= n) {
  xa     <- c(i:(i+nlength))
  ya     <- x[i:(i+nlength)]
  compa  <- stats::lm(ya ~ xa)$coef[2]
  xb     <- c((i+1+nlength):(i+2*nlength+1))
  yb     <- x[(i+1+nlength):(i+1+2*nlength)]
  compb  <- stats::lm(yb ~ xb)$coef[2]
  cng[i] <-  compb - compa
  i      <- i + 1
  }
 if (details == TRUE) detail  <- data.frame(v=(1:(n-5)),values=x[1:(n-5)], cng)
 cng        <- as.numeric(which(cng==max(cng, na.rm=TRUE))+nlength)
 res        <- list(detail=detail, nFactors=c(cng))
 class(res) <- c("nFactors","list")
 return(res)
 }

