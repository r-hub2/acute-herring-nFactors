#' Parallel Analysis of a Correlation or Covariance Matrix
#'
#' This function gives the distribution of the eigenvalues of correlation or a
#' covariance matrices of random uncorrelated standardized normal variables.
#' The mean and a selected quantile of this distribution are returned.
#'
#' Note that if the decision is based on a quantile value rather than on the
#' mean, care must be taken with the number of replications (\code{rep}). In
#' fact, the smaller the quantile (\code{cent}), the bigger the number of
#' necessary replications.
#'
#' @param subject numeric: nmber of subjects (default is 100)
#' @param var numeric: number of variables (default is 10)
#' @param rep numeric: number of replications of the correlation matrix
#' (default is 100)
#' @param cent depreciated numeric (use quantile instead): quantile of the
#' distribution on which the decision is made (default is 0.05)
#' @param quantile numeric: quantile of the distribution on which the decision
#' is made (default is 0.05)
#' @param model character: \code{"components"} or \code{"factors"}
#' @param sd numeric: vector of standard deviations of the simulated variables
#' (for a parallel analysis on a covariance matrix)
#' @param ...  variable: other parameters for the \code{"mvrnorm"}, \code{corr}
#' or \code{cov} functions
#' @return \item{eigen}{ Data frame consisting of the mean and the quantile of
#' the eigenvalues distribution } \item{eigen$mevpea}{ Mean of the eigenvalues
#' distribution} \item{eigen$sevpea}{ Standard deviation of the eigenvalues
#' distribution} \item{eigen$qevpea}{ quantile of the eigenvalues distribution}
#' \item{eigen$sqevpea}{ Standard error of the quantile of the eigenvalues
#' distribution} \item{subject}{ Number of subjects} \item{variables}{ Number
#' of variables} \item{centile}{ Selected quantile} Otherwise, returns a
#' summary of the parallel analysis.
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{plotuScree}}, \code{\link{nScree}},
#' \code{\link{plotnScree}}, \code{\link{plotParallel}}
#' @references Drasgow, F. and Lissak, R. (1983) Modified parallel analysis: a
#' procedure for examining the latent dimensionality of dichotomously scored
#' item responses. \emph{Journal of Applied Psychology, 68}(3), 363-373.
#'
#' Hoyle, R. H. and Duvall, J. L. (2004). Determining the number of factors in
#' exploratory and confirmatory factor analysis.  In D. Kaplan (Ed.): \emph{The
#' Sage handbook of quantitative methodology for the social sciences}. Thousand
#' Oaks, CA: Sage.
#'
#' Horn, J. L. (1965). A rationale and test of the number of factors in factor
#' analysis. \emph{Psychometrika, 30}, 179-185.
#' @export
#' @importFrom MASS ginv mvrnorm
#' @importFrom stats cov dnorm qnorm
#' @keywords multivariate
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
#'  quantile <- 0.95
#'  results  <- parallel(subject, var, rep, quantile)
#'
#'  results
#'
#' ## IF THE DECISION IS BASED ON THE CENTILE USE qevpea INSTEAD
#' ## OF mevpea ON THE FIRST LINE OF THE FOLLOWING CALL
#'  plotuScree(x    = eig,
#'             main = "Parallel Analysis"
#'             )
#'
#'  lines(1:var,
#'        results$eigen$qevpea,
#'        type="b",
#'        col="green"
#'        )
#'
#'
#' ## ANOTHER SOLUTION IS SIMPLY TO
#'  plotParallel(results)
#'  }
#' }
"parallel" <-
function(subject=100, var=10, rep=100, cent=0.05, quantile=cent, model="components", sd=diag(1,var), ...)
 {
  r             <- subject
  c             <- var
  y             <- matrix(c(1:r*c), nrow=r, ncol=c)
  ycor          <- matrix(c(1:c*c), nrow=c, ncol=c)
  evpea         <- NULL
  leg.txt       <- "Pearson"

  # Simulation of k samples to obtain k random eigenvalues vectors
  # for Pearson correlation coefficients
  for (k in c(1:rep)) {
   # y              <- rnorm(y, sd=sqrt(mean(diag(sd))))  # Old version without covariance
   # y              <- matrix(y, nrow=r, ncol=c)          # Old version without covariance
   y <- MASS::mvrnorm(n = r, mu=rep(0,var), Sigma=sd, empirical=FALSE)
   corY           <- stats::cov(y, ...) # The previous version was only cor(y)
   if (model == "components") diag(corY) <- diag(sd) # To constraint the diagonal to sd for PCA
   if (model == "factors") corY <- corY - MASS::ginv(diag(diag(MASS::ginv(corY)))) # To constraint the diagonal to communalities for FCA
   evpea          <- rbind(evpea, eigen(corY)[[1]])
   }
  # Temporay function to compute the standard error of a quantile
  SEcentile <- function(sd, n = 100, p = 0.95) {return(sd/sqrt(n) * sqrt(p*(1-p))/stats::dnorm(stats::qnorm(p))) }

  # Summary statistics
  sprob         <- c(cent)
  mevpea        <- sapply(as.data.frame(evpea),  mean)                          # Eigenvalues means
  sevpea        <- sapply(as.data.frame(evpea),  sd  )                          # Eigenvalues Standard deviations
  qevpea        <- moreStats(evpea, quantile=quantile)[3,]                      # Would be more in line with version 2.3
  #quant         <- function(x, sprobs = sprobs) {return(as.vector(quantile(x, probs = sprob))) }
  #qevpea        <- sapply(as.data.frame(evpea),  quant)                         # Eigenvalues centiles
  sqevpea       <- sevpea
  sqevpea       <- sapply(as.data.frame(sqevpea), SEcentile, n = rep, p = cent) # Standard error of the centiles

  # List of results return
  result        <- list(eigen     = data.frame(mevpea, sevpea, qevpea, sqevpea),
                        subject   = r,
                        variables = c,
                        centile   = cent
                        )
  class(result) <- 'parallel'                                                  # For future use
  return(result)
 }

