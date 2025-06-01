#' Bentler and Yuan's Computation of the LRT Index and the Linear Trend
#' Coefficients
#'
#' This function computes the Bentler and Yuan's (1996, 1998) \emph{LRT} index
#' for the linear trend in eigenvalues of a covariance matrix.  The related
#' \eqn{\chi^2} and \emph{p}-value are also computed.  This function is
#' generally called from the \code{nBentler} function.  But it could be of use
#' for graphing the linear trend function and to study it's behavior.
#'
#' The implemented Bentler and Yuan's procedure must be used with care because
#' the minimized function is not always stable. In many cases, constraints must
#' applied to obtain a solution.  The actual implementation did, but the user
#' can modify these constraints.
#'
#' The hypothesis tested (Bentler and Yuan, 1996, equation 10) is: \cr \cr
#'
#' (1) \eqn{\qquad \qquad H_k: \lambda_{k+i} = \alpha + \beta x_i, (i = 1,
#' \ldots, q)} \cr
#'
#' The solution of the following simultaneous equations is needed to find
#' \eqn{(\alpha, \beta) \in} \cr
#'
#' (2) \eqn{\qquad \qquad f(x) = \sum_{i=1}^q \frac{ [ \lambda_{k+j} - N \alpha
#' + \beta x_j ] x_j}{(\alpha + \beta x_j)^2} = 0} \cr \cr
#'
#' and \eqn{\qquad \qquad g(x) = \sum_{i=1}^q \frac{ \lambda_{k+j} - N \alpha +
#' \beta x_j x_j}{(\alpha + \beta x_j)^2} = 0} \cr
#'
#' The solution to this system of equations was implemented by minimizing the
#' following equation: \cr
#'
#' (3) \eqn{\qquad \qquad (\alpha, \beta) \in \inf{[h(x)]} = \inf{\log{[f(x)^2
#' + g(x)^2}}]} \cr
#'
#' The likelihood ratio test \eqn{LRT} proposed by Bentler and Yuan (1996,
#' equation 7) follows a \eqn{\chi^2} probability distribution with \eqn{q-2}
#' degrees of freedom and is equal to: \cr
#'
#' (4) \eqn{\qquad \qquad LRT = N(k - p)\left\{ {\ln \left( {{n \over N}}
#' \right) + 1} \right\} - N\sum\limits_{j = k + 1}^p {\ln \left\{ {{{\lambda
#' _j } \over {\alpha + \beta x_j }}} \right\}} + n\sum\limits_{j = k + 1}^p
#' {\left\{ {{{\lambda _j } \over {\alpha + \beta x_j }}} \right\}} } \cr
#'
#' With \eqn{p} beeing the number of eigenvalues, \eqn{k} the number of
#' eigenvalues to test, \eqn{q} the \eqn{p-k} remaining eigenvalues, \eqn{N}
#' the sample size, and \eqn{n = N-1}.  Note that there is an error in the
#' Bentler and Yuan equation, the variables \eqn{N} and \eqn{n} beeing inverted
#' in the preceeding equation 4.
#'
#' A better strategy proposed by Bentler an Yuan (1998) is to use a minimized
#' \eqn{\chi^2} solution. This strategy will be implemented in a future version
#' of the \pkg{nFactors} package.
#'
#' @param x numeric: a \code{vector} of eigenvalues, a \code{matrix} of
#' correlations or of covariances or a \code{data.frame} of data
#' @param N numeric: number of subjects.
#' @param nFactors numeric: number of components to test.
#' @param log logical: if \code{TRUE} the minimization is applied on the log
#' values.
#' @param cor logical: if \code{TRUE} computes eigenvalues from a correlation
#' matrix, else from a covariance matrix
#' @param minPar numeric: minimums for the coefficient of the linear trend.
#' @param maxPar numeric: maximums for the coefficient of the linear trend.
#' @param resParx numeric: restriction on the \eqn{\alpha} coefficient (x) to
#' graph the function to minimize.
#' @param resPary numeric: restriction on the \eqn{\beta} coefficient (y) to
#' graph the function to minimize.
#' @param graphic logical: if \code{TRUE} plots the minimized function
#' \code{"wireframe"}, \code{"contourplot"} or \code{"levelplot"}.
#' @param resolution numeric: resolution of the 3D graph (number of points from
#' \eqn{\alpha} and from \eqn{\beta}).
#' @param typePlot character: plots the minimized function according to a 3D
#' plot: \code{"wireframe"}, \code{"contourplot"} or \code{"levelplot"}.
#' @param ...  variable: additionnal parameters from the \code{"wireframe"},
#' \code{"contourplot"} or \code{"levelplot"} \code{lattice} functions. Also
#' additionnal parameters for the \code{eigenFrom} function.
#'
#' @return \item{nFactors}{ numeric: vector of the number of factors retained
#' by the Bentler and Yuan's procedure. } \item{details}{ numeric: matrix of
#' the details of the computation.}
#'
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' \cr \cr David Magis \cr Departement de mathematiques \cr Universite de Liege
#' \cr \email{David.Magis@@ulg.ac.be}
#' @seealso \code{\link{nBartlett}}, \code{\link{nBentler}}
#' @references
#' Bentler, P. M. and Yuan, K.-H. (1996). Test of linear trend in
#' eigenvalues of a covariance matrix with application to data analysis.
#' \emph{British Journal of Mathematical and Statistical Psychology, 49},
#' 299-312.
#'
#' Bentler, P. M. and Yuan, K.-H. (1998). Test of linear trend in the smallest
#' eigenvalues of the correlation matrix. \emph{Psychometrika, 63}(2), 131-144.
#' @importFrom lattice wireframe contourplot levelplot
#' @importFrom stats nlminb lm pchisq
#' @export
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## ................................................
#' ## SIMPLE EXAMPLE OF THE BENTLER AND YUAN PROCEDURE
#' ## #' @importFrom graphics abline
#'
#' # Bentler (1996, p. 309) Table 2 - Example 2 .............
#' n=649
#' bentler2<-c(5.785, 3.088, 1.505, 0.582, 0.424, 0.386, 0.360, 0.337, 0.303,
#'             0.281, 0.246, 0.238, 0.200, 0.160, 0.130)
#'
#' results  <- nBentler(x=bentler2, N=n,  details=TRUE)
#' results
#'
#' # Two different figures to verify the convergence problem identified with
#' # the 2th component
#' bentlerParameters(x=bentler2, N=n, nFactors= 2, graphic=TRUE,
#'                   typePlot="contourplot",
#'                   resParx=c(0,9), resPary=c(0,9), cor=FALSE)
#'
#' bentlerParameters(x=bentler2, N=n, nFactors= 4, graphic=TRUE, drape=TRUE,
#'                   resParx=c(0,9), resPary=c(0,9),
#'                   scales = list(arrows = FALSE) )
#'
#' plotuScree(x=bentler2, model="components",
#'   main=paste(results$nFactors,
#'   " factors retained by the Bentler and Yuan's procedure (1996, p. 309)",
#'   sep=""))
#' # ........................................................
#'
#' # Bentler (1998, p. 140) Table 3 - Example 1 .............
#' n        <- 145
#' example1 <- c(8.135, 2.096, 1.693, 1.502, 1.025, 0.943, 0.901, 0.816,
#'               0.790,0.707, 0.639, 0.543,0.533, 0.509, 0.478, 0.390,
#'               0.382, 0.340, 0.334, 0.316, 0.297,0.268, 0.190, 0.173)
#'
#' results  <- nBentler(x=example1, N=n,  details=TRUE)
#' results
#'
#' # Two different figures to verify the convergence problem identified with
#' # the 10th component
#' bentlerParameters(x=example1, N=n, nFactors= 10, graphic=TRUE,
#'                   typePlot="contourplot",
#'                   resParx=c(0,0.4), resPary=c(0,0.4))
#'
#' bentlerParameters(x=example1, N=n, nFactors= 10, graphic=TRUE, drape=TRUE,
#'                   resParx=c(0,0.4), resPary=c(0,0.4),
#'                   scales = list(arrows = FALSE) )
#'
#' plotuScree(x=example1, model="components",
#'    main=paste(results$nFactors,
#'    " factors retained by the Bentler and Yuan's procedure (1998, p. 140)",
#'    sep=""))
#' # ........................................................
#'  }
#' }
#'
bentlerParameters <-
function(x, N, nFactors, log=TRUE, cor=TRUE,
         minPar=c(min(lambda) - abs(min(lambda)) +.001, 0.001),
         maxPar=c(max(lambda), stats::lm(lambda ~ I(length(lambda):1))$coef[2]),
         resParx=c(0.01, 2), resPary=c(0.01, 2),
         graphic=TRUE, resolution=30, typePlot="wireframe", ...){
 stopMessage  <- paste("\n These indices are only valid with a principal component solution.\n",
                       " ...................... So, only positive eigenvalues are permitted.\n",
                       sep="")
 lambda       <- eigenComputes(x, cor=cor, ...)
 if (length(which(lambda <0 )) > 0) {cat(stopMessage);stop()}

 k     <- nFactors
 p     <- length(lambda)
 q     <- p-k
 i     <- 1:q
 x     <- q-i
 l     <- lambda[k+i]
 n     <- N - 1

 # Bentler (1996, p. 133) maximization of equations 8 and 9
 f1    <- function(n,l,x,alpha,beta) sum((n*l-(n+1)*(alpha+beta*x))/((alpha+beta*x)^2))
 f2    <- function(n,l,x,alpha,beta) sum((n*l-(n+1)*(alpha+beta*x))*x/((alpha+beta*x)^2))
 f     <- function(alpha,beta) f1(n,l,x,alpha,beta)^2+f2(n,l,x,alpha,beta)^2
 if (log == FALSE)  F <- function(y) f(y[1],y[2])  else  F <- function(y) log(f(y[1],y[2]))

 figure <- NULL
 if (graphic == TRUE) {
  p1        <- seq(resParx[1], resParx[2], length=resolution)
  p2        <- seq(resPary[1], resPary[2], length=resolution)
  data      <- expand.grid(Alpha = p1, Beta = p2)
  data      <- data.frame(data, y=numeric(length(data$Alpha)))
  for( i in 1:length(data$Alpha)) data$y[i] <- F(c(data$Alpha[i],data$Beta[i]))

  if (log == FALSE) zlab <- "y" else zlab <- "log(y)"
  if (typePlot == "wireframe")   figure    <- lattice::wireframe(  y ~ Alpha * Beta, data=data, zlab=zlab, ...)
  if (typePlot == "contourplot") figure    <- lattice::contourplot(y ~ Alpha * Beta, data=data, region=TRUE, ...)
  if (typePlot == "levelplot")   figure    <- lattice::levelplot(  y ~ Alpha * Beta, data=data, region=TRUE, ...)
  }

 res   <- stats::nlminb(objective=F,start=stats::lm(l~x)$coefficients,lower=c(minPar[1],minPar[2]),upper=c(maxPar[1],maxPar[2]))
 para  <- res$par[1]
 parb  <- res$par[2]
 # Bentler (1996, p. 133) equation 7
 # !!! Warning: Bentler and Yuan (1998) were in error for the definition of LRT !!!
 # !!! So N and n must be inversed in the first logarithm                       !!!
 lrt   <- N*(k-p)*(log(n/N)+1)-N*sum(log(lambda[(k+1):p]/(para+parb*x))) + n*sum(lambda[(k+1):p]/(para+parb*x))
 df    <- q-2
 resp  <- list(convergence=res$convergence, figure=figure, coefficients=res$par,
              lrt=lrt, df=df,k=k,p.value=1-stats::pchisq(lrt,df))
 names(resp$coefficients)<-c("alpha","beta")
 return(resp)
 }
