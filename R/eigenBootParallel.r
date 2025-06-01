# PERMETTRE LES CALCULS AVEC DES DONNEES DISCRETES AUSSI

#' Bootstrapping of the Eigenvalues From a Data Frame
#'
#' The \code{eigenBootParallel} function samples observations from a
#' \code{data.frame} to produce correlation or covariance matrices from which
#' eigenvalues are computed. The function returns statistics about these
#' bootstrapped eigenvalues. Their means or their quantile could be used later
#' to replace the eigenvalues inputted to a parallel analysis.  The
#' \code{eigenBootParallel} can also compute random eigenvalues from empirical
#' data by column permutation (Buja and Eyuboglu, 1992).
#'
#'
#' @param x data.frame: data from which a correlation matrix will be obtained
#' @param quantile numeric: eigenvalues quantile to be reported
#' @param nboot numeric: number of bootstrap samples
#' @param option character: \code{"permutation"} or \code{"bootstrap"}
#' @param cor logical: if \code{TRUE} computes eigenvalues from a correlation
#' matrix, else from a covariance matrix (\code{eigenComputes})
#' @param model character: bootstraps from a principal component analysis
#' (\code{"components"}) or from a factor analysis (\code{"factors"})
#' @param ...  variable: additionnal parameters to give to the \code{cor} or
#' \code{cov} functions
#' @return \item{values}{ data.frame: mean, median, quantile, standard
#' deviation, minimum and maximum of bootstrapped eigenvalues }
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{principalComponents}},
#' \code{\link{iterativePrincipalAxis}}, \code{\link{rRecovery}}
#' @references Buja, A. and Eyuboglu, N. (1992). Remarks on parallel analysis.
#' \emph{Multivariate Behavioral Research, 27}(4), 509-540.
#'
#' Zwick, W. R. and Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain.  \emph{Psychological
#' bulletin, 99}, 432-442.
#' @keywords multivariate
#' @export
#' @importFrom stats cov cor
#' @examples
#' \dontrun{
#' if(interactive()){
#' # .......................................................
#' # Example from the iris data
#'  eigenvalues <- eigenComputes(x=iris[,-5])
#'
#' # Permutation parallel analysis distribution
#'  aparallel   <- eigenBootParallel(x=iris[,-5], quantile=0.95)$quantile
#'
#' # Number of components to retain
#'  results     <- nScree(x = eigenvalues, aparallel = aparallel)
#'  results$Components
#'  plotnScree(results)
#' # ......................................................
#'
#' # ......................................................
#' # Bootstrap distributions study of the eigenvalues from iris data
#' # with different correlation methods
#'  eigenBootParallel(x=iris[,-5],quantile=0.05,
#'                    option="bootstrap",method="pearson")
#'  eigenBootParallel(x=iris[,-5],quantile=0.05,
#'                    option="bootstrap",method="spearman")
#'  eigenBootParallel(x=iris[,-5],quantile=0.05,
#'                    option="bootstrap",method="kendall")
#'  }
#' }
#'
"eigenBootParallel" <-
function(x, quantile=0.95, nboot=30, option="permutation", cor=TRUE, model="components", ...)
 {

 if (eigenFrom(x) != "data") stop("Only data from a data.frame must be used as input")

 x          <- data.frame(x)
 res        <- data.frame(matrix(NA, ncol=dim(x)[2], nrow=nboot))
 if (model == "components") { names(res) <- paste("C", 1:dim(x)[2], sep="")
  } else names(res) <- paste("F", 1:dim(x)[2], sep="")

 if (option == "permutation") {
  for (i in 1:nboot) {
   rPerm   <- apply(x,2,sample, replace=TRUE)
   if (cor == TRUE)        corY <- stats::cor(rPerm, ...)
   if (cor == FALSE)       corY <- stats::cov(rPerm, ...)
   if (model == "factors") corY <- corFA(corY, method="ginv")
   res[i,] <- eigen(corY, only.values=TRUE)$values
   }
  }

 if (option == "bootstrap") {
  for (i in 1:nboot) {
   rBoot   <- sample(1:dim(x)[1], dim(x)[1], replace=TRUE)
   if (cor == TRUE)        corY <- stats::cor(x[rBoot,], ...)
   if (cor == FALSE)       corY <- stats::cov(x[rBoot,], ...)
   if (model == "factors") corY <- corFA(corY, method="ginv")
   res[i,] <- eigen(corY, only.values=TRUE)$values
   #if (cor == TRUE)  res[i,] <- eigen(stats::cor(x[rBoot,], ...), only.values=TRUE)$values
   #if (cor == FALSE) res[i,] <- eigen(stats::cov(x[rBoot,], ...), only.values=TRUE)$values
   }
  }

 res        <- data.frame(t(moreStats(res, quantile=quantile)))
 return(res)
 }
