#' Bentler and Yuan's Procedure to Determine the Number of Components/Factors
#'
#' This function computes the Bentler and Yuan's indices for determining the
#' number of components/factors to retain.
#'
#' The implemented Bentler and Yuan's procedure must be used with care because
#' the minimized function is not always stable, as Bentler and Yan (1996, 1998)
#' already noted. In many cases, constraints must applied to obtain a solution,
#' as the actual implementation did, but the user can modify these constraints.
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
#' + \beta x_j ] x_j}{(\alpha + \beta x_j)^2} = 0} \cr \cr and \eqn{\qquad
#' \qquad g(x) = \sum_{i=1}^q \frac{ \lambda_{k+j} - N \alpha + \beta x_j
#' x_j}{(\alpha + \beta x_j)^2} = 0} \cr
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
#' A better strategy proposed by Bentler an Yuan (1998) is to used a minimized
#' \eqn{\chi^2} solution. This strategy will be implemented in a future version
#' of the \pkg{nFactors} package.
#'
#' @param x numeric: a \code{vector} of eigenvalues, a \code{matrix} of
#' correlations or of covariances or a \code{data.frame} of data
#' @param N numeric: number of subjects.
#' @param log logical: if \code{TRUE} does the maximization on the log values.
#' @param alpha numeric: statistical significance level.
#' @param cor logical: if \code{TRUE} computes eigenvalues from a correlation
#' matrix, else from a covariance matrix
#' @param details logical: if \code{TRUE} also returns detains about the
#' computation for each eigenvalue.
#' @param minPar numeric: minimums for the coefficient of the linear trend to
#' maximize.
#' @param maxPar numeric: maximums for the coefficient of the linear trend to
#' maximize.
#' @param ...  variable: additionnal parameters to give to the \code{cor} or
#' \code{cov} functions
#' @return \item{nFactors}{ numeric: vector of the number of factors retained
#' by the Bentler and Yuan's procedure. } \item{details}{ numeric: matrix of
#' the details of the computation.}
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' \cr \cr David Magis \cr Departement de mathematiques \cr Universite de Liege
#' \cr \email{David.Magis@@ulg.ac.be}
#' @seealso \code{\link{nBartlett}}, \code{\link{bentlerParameters}}
#' @references Bentler, P. M. and Yuan, K.-H. (1996). Test of linear trend in
#' eigenvalues of a covariance matrix with application to data analysis.
#' \emph{British Journal of Mathematical and Statistical Psychology, 49},
#' 299-312.
#'
#' Bentler, P. M. and Yuan, K.-H. (1998). Test of linear trend in the smallest
#' eigenvalues of the correlation matrix. \emph{Psychometrika, 63}(2), 131-144.
#' @export
#' @importFrom stats lm
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## ................................................
#' ## SIMPLE EXAMPLE OF THE BENTLER AND YUAN PROCEDURE
#'
#' # Bentler (1996, p. 309) Table 2 - Example 2 .............
#' n=649
#' bentler2<-c(5.785, 3.088, 1.505, 0.582, 0.424, 0.386, 0.360, 0.337, 0.303,
#'             0.281, 0.246, 0.238, 0.200, 0.160, 0.130)
#'
#' results  <- nBentler(x=bentler2, N=n)
#' results
#'
#' plotuScree(x=bentler2, model="components",
#'     main=paste(results$nFactors,
#'     " factors retained by the Bentler and Yuan's procedure (1996, p. 309)",
#'     sep=""))
#' # ........................................................
#'
#' # Bentler (1998, p. 140) Table 3 - Example 1 .............
#' n        <- 145
#' example1 <- c(8.135, 2.096, 1.693, 1.502, 1.025, 0.943, 0.901, 0.816, 0.790,
#'               0.707, 0.639, 0.543,
#'               0.533, 0.509, 0.478, 0.390, 0.382, 0.340, 0.334, 0.316, 0.297,
#'               0.268, 0.190, 0.173)
#'
#' results  <- nBentler(x=example1, N=n)
#' results
#'
#' plotuScree(x=example1, model="components",
#'    main=paste(results$nFactors,
#'    " factors retained by the Bentler and Yuan's procedure (1998, p. 140)",
#'    sep=""))
#' # ........................................................
#'  }
#' }
nBentler <-
function(x, N, log=TRUE, alpha=0.05, cor=TRUE, details=TRUE,
         minPar=c(min(lambda) - abs(min(lambda)) +.001, 0.001),
         maxPar=c(max(lambda), stats::lm(lambda ~ I(length(lambda):1))$coef[2]),
         ...) {
 stopMessage  <- paste("\n These indices are only valid with a principal component solution.\n",
                       " ...................... So, only positive eugenvalues are permitted.\n",
                       sep="")
 lambda       <- eigenComputes(x, cor=cor, ...)
 if (length(which(lambda <0 )) > 0) {cat(stopMessage);stop()}

 n            <- N
 significance <- alpha
 min.k        <- 3
 LRT          <- data.frame(q=numeric(length(lambda)-min.k), k=numeric(length(lambda)-min.k),
                            LRT=numeric(length(lambda)-min.k), a=numeric(length(lambda)-min.k),
                            b=numeric(length(lambda)-min.k),
                            p=numeric(length(lambda)-min.k),
                            convergence=numeric(length(lambda)-min.k))
 bentler.n    <- 0
 for (i in 1:(length(lambda)-min.k)) {
  temp     <- bentlerParameters(x=lambda, N=n, nFactors=i, log=log, cor=cor, minPar=minPar, maxPar=maxPar)
  LRT[i,3] <- temp$lrt
  LRT[i,4] <- ifelse(is.null(temp$coef[1]),     NA, temp$coef[1])
  LRT[i,5] <- ifelse(is.null(temp$coef[2]),     NA, temp$coef[2])
  LRT[i,6] <- ifelse(is.null(temp$p.value),     NA, temp$p.value)
  LRT[i,7] <- ifelse(is.null(temp$convergence), NA, temp$convergence)
  LRT[i,2] <- i
  LRT[i,1] <- length(lambda) - i
  }
 #LRT     <- LRT[order(LRT[,1],decreasing = TRUE),]
 for (i in 1:(length(lambda)-min.k)) {
  if (i == 1)                         bentler.n <- bentler.n + as.numeric(LRT$p[i] <= significance)
  if (i > 1) {if(LRT$p[i-1] <= 0.05)  bentler.n <- bentler.n + as.numeric(LRT$p[i] <= significance)}
  }
 if (bentler.n == 0)  bentler.n <- length(lambda)
 if (details == TRUE) details <- LRT else details <- NULL
 res        <- list(detail=details, nFactors=bentler.n)
 class(res) <- c("nFactors","list")
 return(res)
 }
