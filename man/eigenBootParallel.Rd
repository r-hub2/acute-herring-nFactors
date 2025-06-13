% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eigenBootParallel.r
\name{eigenBootParallel}
\alias{eigenBootParallel}
\title{Bootstrapping of the Eigenvalues From a Data Frame}
\usage{
eigenBootParallel(
  x,
  quantile = 0.95,
  nboot = 30,
  option = "permutation",
  cor = TRUE,
  model = "components",
  ...
)
}
\arguments{
\item{x}{data.frame: data from which a correlation matrix will be obtained}

\item{quantile}{numeric: eigenvalues quantile to be reported}

\item{nboot}{numeric: number of bootstrap samples}

\item{option}{character: \code{"permutation"} or \code{"bootstrap"}}

\item{cor}{logical: if \code{TRUE} computes eigenvalues from a correlation
matrix, else from a covariance matrix (\code{eigenComputes})}

\item{model}{character: bootstraps from a principal component analysis
(\code{"components"}) or from a factor analysis (\code{"factors"})}

\item{...}{variable: additionnal parameters to give to the \code{cor} or
\code{cov} functions}
}
\value{
\item{values}{ data.frame: mean, median, quantile, standard
deviation, minimum and maximum of bootstrapped eigenvalues }
}
\description{
The \code{eigenBootParallel} function samples observations from a
\code{data.frame} to produce correlation or covariance matrices from which
eigenvalues are computed. The function returns statistics about these
bootstrapped eigenvalues. Their means or their quantile could be used later
to replace the eigenvalues inputted to a parallel analysis.  The
\code{eigenBootParallel} can also compute random eigenvalues from empirical
data by column permutation (Buja and Eyuboglu, 1992).
}
\examples{
\dontrun{
if(interactive()){
# .......................................................
# Example from the iris data
 eigenvalues <- eigenComputes(x=iris[,-5])

# Permutation parallel analysis distribution
 aparallel   <- eigenBootParallel(x=iris[,-5], quantile=0.95)$quantile

# Number of components to retain
 results     <- nScree(x = eigenvalues, aparallel = aparallel)
 results$Components
 plotnScree(results)
# ......................................................

# ......................................................
# Bootstrap distributions study of the eigenvalues from iris data
# with different correlation methods
 eigenBootParallel(x=iris[,-5],quantile=0.05,
                   option="bootstrap",method="pearson")
 eigenBootParallel(x=iris[,-5],quantile=0.05,
                   option="bootstrap",method="spearman")
 eigenBootParallel(x=iris[,-5],quantile=0.05,
                   option="bootstrap",method="kendall")
 }
}

}
\references{
Buja, A. and Eyuboglu, N. (1992). Remarks on parallel analysis.
\emph{Multivariate Behavioral Research, 27}(4), 509-540.

Zwick, W. R. and Velicer, W. F. (1986). Comparison of five rules for
determining the number of components to retain.  \emph{Psychological
bulletin, 99}, 432-442.
}
\seealso{
\code{\link{principalComponents}},
\code{\link{iterativePrincipalAxis}}, \code{\link{rRecovery}}
}
\author{
Gilles Raiche \cr Centre sur les Applications des Modeles de
Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
\email{raiche.gilles@uqam.ca}
}
\keyword{multivariate}
