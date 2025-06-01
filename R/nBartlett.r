#'
#' Bartlett, Anderson and Lawley Procedures to Determine the Number of Components/Factors
#'
#' This function computes the Bartlett, Anderson and Lawley indices for determining the
#' number of components/factors to retain.
#' @details Note: the latex formulas are available only in the pdf version of this help file.
#'
#' The hypothesis tested is: \cr
#'
#' (1)  \eqn{\qquad \qquad H_k: \lambda_{k+1} = \ldots = \lambda_p} \cr
#'
#' This hypothesis is verified by the application of different version of a
#' \eqn{\chi^2} test with different values for the degrees of freedom.
#' Each of these tests shares the compution of a \eqn{V_k} value: \cr
#'
#' (2) \eqn{\qquad \qquad V_k  =
#'   \prod\limits_{i = k + 1}^p {\left\{ {{{\lambda _i }
#'     \over {{\raise0.7ex\hbox{$1$} \!\mathord{\left/
#'     {\vphantom {1 q}}\right.\kern-\nulldelimiterspace}
#'       \!\lower0.7ex\hbox{$q$}}\sum\limits_{i = k + 1}^p {\lambda _i } }}} \right\}}
#' }
#'
#' \eqn{p} is the number of eigenvalues, \eqn{k} the number of eigenvalues to test,
#' and \eqn{q} the \eqn{p-k} remaining eigenvalues. \eqn{n} is equal to the sample size
#' minus 1 (\eqn{n = N-1}). \cr
#'
#' The Anderson statistic is distributed as a \eqn{\chi^2} with \eqn{(q + 2)(q - 1)/2} degrees
#' of freedom and is equal to: \cr
#'
#' (3) \eqn{\qquad \qquad - n\log (V_k ) \sim \chi _{(q + 2)(q - 1)/2}^2 } \cr
#'
#' An improvement of this statistic from Bartlett (Bentler, and Yuan, 1996, p. 300;
#'                                                 Horn and Engstrom, 1979, equation 8) is distributed as a \eqn{\chi^2}
#' with \eqn{(q)(q - 1)/2} degrees of freedom and is equal to: \cr
#'
#' (4) \eqn{\qquad \qquad - \left[ {n - k - {{2q^2 q + 2} \over {6q}}}
#'                                 \right]\log (V_k ) \sim \chi _{(q + 2)(q - 1)/2}^2 }  \cr
#'
#' Finally, Anderson (1956) and James (1969) proposed another statistic. \cr
#'
#' (5) \eqn{\qquad \qquad - \left[ {n - k - {{2q^2 q + 2} \over {6q}}
#'   + \sum\limits_{i = 1}^k {{{\bar \lambda _q^2 } \over {\left( {\lambda _i
#'     - \bar \lambda _q } \right)^2 }}} } \right]\log (V_k ) \sim \chi _{(q + 2)(q - 1)/2}^2 } \cr
#'
#' Bartlett (1950, 1951) proposed a correction to the degrees of freedom of these \eqn{\chi^2} after the
#' first significant test: \eqn{(q+2)(q - 1)/2}. \cr
#'
#' @param x          numeric: a \code{vector} of eigenvalues, a \code{matrix} of correlations or of covariances or a \code{data.frame} of data (eigenFrom)
#' @param N          numeric: number of subjects
#' @param alpha      numeric: statistical significance level
#' @param cor        logical: if \code{TRUE} computes eigenvalues from a correlation matrix, else from a covariance matrix
#' @param details    logical: if \code{TRUE} also returns detains about the computation for each eigenvalue
#' @param correction logical: if \code{TRUE} uses a correction for the degree of freedom after the first eigenvalue
#' @param ...        variable: additionnal parameters to give to the \code{cor} or \code{cov} functions
#' @return \item{nFactors}{numeric: vector of the number of factors retained by the Bartlett, Anderson and Lawley procedures.}
#' @return \item{details}{numeric: matrix of the details for each index.}
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{plotuScree}}, \code{\link{nScree}}, \code{\link{plotnScree}}, \code{\link{plotParallel}}
#'
#' @references
#'  Anderson, T. W. (1963). Asymptotic theory for principal component analysis. \emph{Annals of Mathematical Statistics, 34}, 122-148.
#'
#'  Bartlett, M. S. (1950). Tests of significance in factor analysis. \emph{British Journal of Psychology, 3}, 77-85.
#'
#'  Bartlett, M. S. (1951). A further note on tests of significance. \emph{British Journal of Psychology, 4}, 1-2.
#'
#'  Bentler, P. M. and Yuan, K.-H. (1996). Test of linear trend in eigenvalues of a covariance matrix with application to data analysis.
#'  \emph{British Journal of Mathematical and Statistical Psychology, 49}, 299-312.
#'
#'  Bentler, P. M. and Yuan, K.-H. (1998). Test of linear trend in the smallest
#'  eigenvalues of the correlation matrix. \emph{Psychometrika, 63}(2), 131-144.
#'
#'  Horn, J. L. and Engstrom, R. (1979). Cattell's scree test in relation to
#'  Bartlett's chi-square test and other observations on the number of factors
#'  problem. \emph{Multivariate Behavioral Reasearch, 14}(3), 283-300.
#'
#'  James, A. T. (1969). Test of equality of the latent roots of the covariance
#'  matrix. \emph{In} P. K. Krishna (Eds): \emph{Multivariate analysis, volume 2}.New-York, NJ: Academic Press.
#'
#'  Lawley, D. N. (1956). Tests of significance for the latent roots of covarianceand correlation matrix. \emph{Biometrika, 43}(1/2), 128-136.
#'
#' @export
#' @importFrom stats pchisq
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## ................................................
#' ## SIMPLE EXAMPLE OF A BARTLETT PROCEDURE
#'
#' data(dFactors)
#' eig      <- dFactors$Raiche$eigenvalues
#'
#' results  <- nBartlett(x=eig, N= 100, alpha=0.05, details=TRUE)
#' results
#'
#' plotuScree(eig, main=paste(results$nFactors[1], ", ",
#'                            results$nFactors[2], " or ",
#'                            results$nFactors[3],
#'                            " factors retained by the LRT procedures",
#'                            sep=""))
#'  }
#' }
nBartlett <-
function(x, N, alpha=0.05, cor=TRUE, details=TRUE, correction=TRUE, ...) {
 stopMessage  <- paste("\n These indices are only valid with a principal component solution.\n",
                       " ...................... So, only positive eugenvalues are permitted.\n",
                       sep="")
 x            <- eigenComputes(x, cor=cor, ...)
 if (length(which(x<0)) > 0) {cat(stopMessage);stop()}

 n            <- length(x)
 detail       <- NULL
 bartlett.n   <- anderson.n   <- lawley.n                   <- 0
 bartlett     <- bartlett.chi <- bartlett.df <- bartlett.p  <- numeric(n)
 anderson.chi <- anderson.df  <- anderson.p                 <- numeric(n)
 lawley.chi   <- lawley.df    <- lawley.p                   <- numeric(n)
 for (k in 0:(n-1)) {
  i <- k+1
  bartlett[i]     <- prod(x[(k+1):n]) /  (sum(x[(k+1):n])/(n-k))^(n-k) # From Horn et Engstrom (1979)
  bartlett.chi[i] <- -(N - 1 - ((2*n+5)/6) - ((2*k)/3)) * log(bartlett[i])
  bartlett.df[i]  <- .5 * (n-k) * (n-k-1)   # Bartlett without correction, from Horn and Engstrom (1979. p. 291, equation 8)
  if (correction==TRUE & bartlett.n > 0) bartlett.df[i]  <- .5 * (n-k+2) * (n-k-1)  # From Bentler and Yuan (1996, p. 300)
  bartlett.p[i]   <- stats::pchisq(bartlett.chi[i] , bartlett.df[i], lower.tail = FALSE)
  # Conditions to stop when non significant test are obtained
  anderson.chi[i] <- -N * log(bartlett[i])  # From Bentler and Yuan (1996, p. 300, equations 3-4)
  anderson.df[i]  <- .5 * (n-k+2) * (n-k-1) # From Bentler and Yuan (1996, p. 300)
  anderson.p[i]   <- stats::pchisq(anderson.chi[i] , anderson.df[i], lower.tail = FALSE)
  lMean           <- mean(x[(k+1):n])
  lawley.chi[i]   <- -(N - 1 - ((2*n+5)/6) - ((2*k)/3) + sum((lMean^2)/((x[k]+lMean)^2))) * log(bartlett[i]) # From Bentler and Yuan (1996, p. 300, equation 6)
  lawley.df[i]    <- .5 * (n-k) * (n-k-1) # From Horn and Engstrom (1979. p. 291, equation 8)
  lawley.p[i]     <- stats::pchisq(lawley.chi[i] , lawley.df[i], lower.tail = FALSE)
# print(c(bartlett[i], bartlett.chi[i], bartlett.df[i], bartlett.p[i]),2)  ############ TEST #############
  if (i == 1) {
   bartlett.n <- bartlett.n + as.numeric(bartlett.p[i] <= alpha)
   anderson.n <- anderson.n + as.numeric(anderson.p[i] <= alpha)
   lawley.n   <- lawley.n   + as.numeric(lawley.p[i]   <= alpha)
      }
  if (i > 1)  {
   if(bartlett.p[i-1] <= 0.05) bartlett.n <- bartlett.n + as.numeric(bartlett.p[i] <= alpha)
   if(anderson.p[i-1] <= 0.05) anderson.n <- anderson.n + as.numeric(anderson.p[i] <= alpha)
   if(lawley.p[i-1]   <= 0.05) lawley.n   <- lawley.n   + as.numeric(lawley.p[i]   <= alpha)
   }
  }
 if (bartlett.n == 0) bartlett.n <- n # If no test if significant, retain all components
 if (anderson.n == 0) anderson.n <- n
 if (lawley.n   == 0) lawlwy.n   <- n
 if (details == TRUE) detail    <- data.frame(v=(1:(n)),values=x[1:(n)],
                                               bartlett, bartlett.chi, bartlett.df, bartlett.p,
                                               anderson.chi, anderson.df, anderson.p,
                                               lawley.chi,   lawley.df,   lawley.p)
 res        <- list(detail=detail,
                    nFactors=c(bartlett=bartlett.n, anderson=anderson.n, lawley=lawley.n))
 class(res) <- c("nFactors","list")
 return(res)
 }

