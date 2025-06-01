#' Principal Component Analysis
#'
#' The \code{principalComponents} function returns a principal component
#' analysis.  Other R functions give the same results, but
#' \code{principalComponents} is customized mainly for the other factor
#' analysis functions available in the \pkg{nfactors} package. In order to
#' retain only a small number of components the \code{componentAxis} function
#' has to be used.
#'
#'
#' @param R numeric: correlation or covariance matrix
#' @return \item{values}{ numeric: variance of each component }
#' \item{varExplained}{ numeric: variance explained by each component }
#' \item{varExplained}{ numeric: cumulative variance explained by each
#' component } \item{loadings}{ numeric: loadings of each variable on each
#' component }
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{componentAxis}}, \code{\link{iterativePrincipalAxis}},
#' \code{\link{rRecovery}}
#' @references Joliffe, I. T. (2002). \emph{Principal components analysis} (2th
#' Edition).  New York, NJ: Springer-Verlag.
#'
#' Kim, J.-O. and Mueller, C. W. (1978). \emph{Introduction to factor analysis.
#' What it is and how to do it}. Beverly Hills, CA: Sage.
#'
#' Kim, J.-O. and Mueller, C. W. (1987). \emph{Factor analysis. Statistical
#' methods and practical issues}. Beverly Hills, CA: Sage.
#' @export
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' # .......................................................
#' # Example from Kim and Mueller (1978, p. 10)
#' # Population: upper diagonal
#' # Simulated sample: lower diagnonal
#'  R <- matrix(c( 1.000, .6008, .4984, .1920, .1959, .3466,
#'                 .5600, 1.000, .4749, .2196, .1912, .2979,
#'                 .4800, .4200, 1.000, .2079, .2010, .2445,
#'                 .2240, .1960, .1680, 1.000, .4334, .3197,
#'                 .1920, .1680, .1440, .4200, 1.000, .4207,
#'                 .1600, .1400, .1200, .3500, .3000, 1.000),
#'                 nrow=6, byrow=TRUE)
#'
#' # Factor analysis: Principal component -
#' # Kim et Mueller (1978, p. 21)
#' # Replace upper diagonal with lower diagonal
#'  RU <- diagReplace(R, upper=TRUE)
#'  principalComponents(RU)
#'
#' # Replace lower diagonal with upper diagonal
#'  RL <- diagReplace(R, upper=FALSE)
#'  principalComponents(RL)
#' # .......................................................
#'  }
#' }
"principalComponents" <-
function(R) {
 nVar   <- dim(R)[2]
 acp             <- eigen(R)
 values          <- acp$values
 vectors         <- acp$vectors # Normed vecteurs to 1
 varExplained    <- round((values/nVar)*100,    2)
 cumVarExplained <- round(cumsum(varExplained), 2)
 loadings        <- vectors  %*% diag(values^0.5)  # F1 * diag(E)
 acp             <- list(values          = values,
                         varExplained    = varExplained,
                         cumVarExplained = cumVarExplained,
                         vectors         = vectors,
                         loadings        = loadings)
 return(acp)
 }

