#' Iterative Principal Axis Analysis
#'
#'  The \code{iterativePrincipalAxis} function returns a principal axis analysis with
#'  iterated communality estimates. Four different choices of initial communality
#'  estimates are given: maximum correlation, multiple correlation (usual and
#'  generalized inverse) or estimates based
#'  on the sum of the squared principal component analysis loadings. Generally, statistical
#'  packages initialize the communalities at the multiple correlation value.
#'  Unfortunately, this strategy cannot always deal with singular correlation or
#'  covariance matrices.
#'  If a generalized inverse, the maximum correlation or the estimated communalities
#'  based on the sum of loadings
#'  are used instead, then a solution can be computed.
#'
#'
#' @param R             numeric: correlation or covariance matrix
#' @param nFactors      numeric: number of factors to retain
#' @param communalities character: initial values for communalities (\code{"component", "maxr", "ginv" or "multiple"})
#' @param iterations    numeric: maximum number of iterations to obtain a solution
#' @param tolerance     numeric: minimal difference in the estimated communalities after a given iteration
#'
#' @return   values       numeric: variance of each component
#' @return   varExplained numeric: variance explained by each component
#' @return   varExplained numeric: cumulative variance explained by each component
#' @return   loadings     numeric: loadings of each variable on each component
#' @return   iterations   numeric: maximum number of iterations to obtain a solution
#' @return   tolerance    numeric: minimal difference in the estimated communalities after a given iteration
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' \cr \cr David Magis \cr Departement de mathematiques \cr Universite de Liege
#' \cr \email{David.Magis@@ulg.ac.be}
#'
#' @references
#' Kim, J.-O. and Mueller, C. W. (1978). \emph{Introduction to factor analysis. What it
#'   is and how to do it}. Beverly Hills, CA: Sage.
#'
#' Kim, J.-O. and Mueller, C. W. (1987). \emph{Factor analysis. Statistical methods and
#'   practical issues}. Beverly Hills, CA: Sage.
#'
#' @export
#' @importFrom MASS ginv
#' @keywords multivariate
#' @seealso \code{\link{componentAxis}}, \code{\link{principalAxis}}, \code{\link{rRecovery}}
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## ................................................
#' # Example from Kim and Mueller (1978, p. 10)
#' # Population: upper diagonal
#' # Simulated sample: lower diagnonal
#' R <- matrix(c( 1.000, .6008, .4984, .1920, .1959, .3466,
#'                .5600, 1.000, .4749, .2196, .1912, .2979,
#'                .4800, .4200, 1.000, .2079, .2010, .2445,
#'                .2240, .1960, .1680, 1.000, .4334, .3197,
#'                .1920, .1680, .1440, .4200, 1.000, .4207,
#'                .1600, .1400, .1200, .3500, .3000, 1.000),
#'             nrow=6, byrow=TRUE)
#'
#' # Factor analysis: Principal axis factoring with iterated communalities
#' # Kim and Mueller (1978, p. 23)
#' # Replace upper diagonal with lower diagonal
#' RU         <- diagReplace(R, upper=TRUE)
#' nFactors   <- 2
#' fComponent <- iterativePrincipalAxis(RU, nFactors=nFactors,
#'                                      communalities="component")
#' fComponent
#' rRecovery(RU,fComponent$loadings, diagCommunalities=FALSE)
#'
#' fMaxr      <- iterativePrincipalAxis(RU, nFactors=nFactors,
#'                                      communalities="maxr")
#' fMaxr
#' rRecovery(RU,fMaxr$loadings, diagCommunalities=FALSE)
#'
#' fMultiple  <- iterativePrincipalAxis(RU, nFactors=nFactors,
#'                                      communalities="multiple")
#' fMultiple
#' rRecovery(RU,fMultiple$loadings, diagCommunalities=FALSE)
#' # .......................................................
#'  }
#' }
#'
iterativePrincipalAxis <-
function(R, nFactors=2, communalities="component", iterations=20, tolerance=0.001) {
 if (communalities == "component")            diag(R)  <- componentAxis(R)$communalities
 if (communalities == "maxr")      { RT <- R; diag(RT) <- 0; diag(R) <- apply(RT, 1, max)}
 if (communalities == "ginv")                 diag(R)  <- sqrt(1-1/diag(MASS::ginv(R)))
 if (communalities == "multiple")  {
  if (all(eigen(R)$values > 0)) diag(R) <- sqrt(1-1/diag(solve(R)))  # Gorsuch (1983, p. 106)
  else return("Not all eigenvalues are grater than 0") # Verication of positive definiteness
  }
  iter <- 1; tol <- 1
  while ((iter < iterations) && (tol > tolerance)) {     # for (i in (1:iterations))
   oldR    <- diag(R)
   diag(R) <- componentAxis(R, nFactors)$communalities
   tol     <- max(abs(diag(R) - oldR))
   iter    <- iter + 1
  }
 if (tol > tolerance) warning("Maximum number of iterations needed before the desired tolerance: cautious solution.")
 iapa <- componentAxis(R, nFactors)
 iapa <- list(values          = iapa$values,
              varExplained    = iapa$varExplained,
              cumVarExplained = iapa$cumVarExplained,
              loadings        = iapa$loadings,
              iterations      = iter,
              tolerance       = tol)
 return(iapa)
 }

