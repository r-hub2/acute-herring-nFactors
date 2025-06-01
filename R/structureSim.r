#' Population or Simulated Sample Correlation Matrix from a Given Factor
#' Structure Matrix
#'
#' The \code{structureSim} function returns a population and a sample
#' correlation matrices from a predefined congeneric factor structure.
#'
#'
#' @param fload matrix: loadings of the factor structure
#' @param reppar numeric: number of replications for the parallel analysis
#' @param repsim numeric: number of replications of the matrix correlation
#' simulation
#' @param N numeric: number of subjects
#' @param quantile numeric: quantile for the parallel analysis
#' @param model character: \code{"components"} or \code{"factors"}
#' @param adequacy logical: if \code{TRUE} prints the recovered population
#' matrix from the factor structure
#' @param details logical: if \code{TRUE} outputs details of the \code{repsim}
#' simulations
#' @param r2limen numeric: R2 limen value for the R2 Nelson index
#' @param all logical: if \code{TRUE} computes the Bentler and Yuan index (very
#' long computing time to consider)
#' @return \item{values}{ the output depends of the logical value of details.
#' If \code{FALSE}, returns only statistics about the eigenvalues: mean,
#' median, quantile, standard deviation, minimum and maximum. If \code{TRUE},
#' returns also details about the \code{repsim} simulations.  If
#' \code{adequacy} = \code{TRUE} returns the recovered factor structure}
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{principalComponents}},
#' \code{\link{iterativePrincipalAxis}}, \code{\link{rRecovery}}
#' @references
#' Raiche, G., Walls, T. A., Magis, D., Riopel, M. and Blais, J.-G. (2013). Non-graphical solutions
#' for Cattell's scree test. Methodology, 9(1), 23-29.
#'
#' Zwick, W. R. and Velicer, W. F. (1986). Comparison of five rules
#' for determining the number of components to retain. \emph{Psychological
#' Bulletin, 99}, 432-442.
#' @export
#' @importFrom stats median factanal
#' @importFrom graphics boxplot plot abline lines
#' @importFrom psych sim.structure
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' # .......................................................
#' # Example inspired from Zwick and Velicer (1986, table 2, p. 437)
#' ## ...................................................................
#'  nFactors  <- 3
#'  unique    <- 0.2
#'  loadings  <- 0.5
#'  nsubjects <- 180
#'  repsim    <- 30
#'  zwick     <- generateStructure(var=36, mjc=nFactors, pmjc=12,
#'                                 loadings=loadings,
#'                                 unique=unique)
#' ## ...................................................................
#'
#' # Produce statistics about a replication of a parallel analysis on
#' # 30 sampled correlation matrices
#'
#'  mzwick.fa <-  structureSim(fload=as.matrix(zwick), reppar=30,
#'                             repsim=repsim, N=nsubjects, quantile=0.5,
#'                             model="factors")
#'
#'  mzwick    <-  structureSim(fload=as.matrix(zwick), reppar=30,
#'                             repsim=repsim, N=nsubjects, quantile=0.5, all=TRUE)
#'
#' # Very long execution time that could be used only with model="components"
#' # mzwick    <-  structureSim(fload=as.matrix(zwick), reppar=30,
#' #                            repsim=repsim, N=nsubjects, quantile=0.5, all=TRUE)
#'
#'  par(mfrow=c(2,1))
#'  graphics::plot(x=mzwick,    nFactors=nFactors, index=c(1:14), cex.axis=0.7, col="red")
#'  graphics::plot(x=mzwick.fa, nFactors=nFactors, index=c(1:11), cex.axis=0.7, col="red")
#'  par(mfrow=c(1,1))
#'
#'  par(mfrow=c(2,1))
#'  graphics::boxplot(x=mzwick,    nFactors=3, cex.axis=0.8, vLine="blue", col="red")
#'  graphics::boxplot(x=mzwick.fa, nFactors=3, cex.axis=0.8, vLine="blue", col="red",
#'          xlab="Components")
#'  par(mfrow=c(1,1))
#' # ......................................................
#'  }
#' }
structureSim <-
function(fload, reppar=30, repsim=100, N, quantile=0.95, model="components",
         adequacy=FALSE, details=TRUE, r2limen=0.75, all=FALSE) {
 simulation   <- psych::sim.structure(fx=fload, n=N, raw=TRUE)
 if (adequacy == TRUE) print(stats::factanal(covmat=simulation$model, factors=dim(fload)[2])) # Verification of the adequacy of the model
 eigenvalues  <- eigenComputes(simulation$r, cor=TRUE, model=model)
 variables    <- length(eigenvalues) # Compute the number of variables
 aparallel    <- parallel(var=dim(fload)[1],subject=N,rep=reppar,cent=quantile,model=model)$eigen$qevpea  # The percentile
 components   <- matrix(NA, ncol=15,nrow=repsim)
 analysis     <- NA
 values       <- matrix(NA, ncol=length(eigenvalues),nrow=repsim)
 for (i in 1:repsim) {
  simulation             <- psych::sim.structure(fx=fload, n=N, raw=TRUE)
  aparallel              <- parallel(var=dim(fload)[1],subject=N,rep=reppar,cent=quantile,model=model)$eigen$qevpea
  eigenvalues            <- eigenComputes(simulation$r, cor=TRUE, model=model)
  values[i,]             <- eigenvalues
  results                <- nScree(x=eigenvalues,aparallel = aparallel, cor=TRUE, model=model)
  components[i,(1:4)]    <- t(results$Components)
  ### PERMUTATIONS
  if (eigenFrom(data.frame(simulation$observed)) == "data")  {
   permutation <- eigenBootParallel(x=data.frame(simulation$observed), quantile=quantile, model=model)$quantile
   }
  results                <- nScree(x=eigenvalues,aparallel = permutation, cor=TRUE, model=model)
  components[i, 5]       <- results$Components$nparallel
  ### ...
  components[i, 6]       <- nCng(x=eigenvalues, model=model)$nFactors
  components[i, (7:9)]   <- nMreg(x=eigenvalues, model=model)$nFactors
  components[i, (10:11)] <- nSeScree(x=eigenvalues, model=model, r2limen=r2limen)$nFactors

  if (model == "components") {
   components[i, (12:14)] <- nBartlett(x=eigenvalues, N=N, alpha=1-quantile, cor=TRUE, correction=TRUE)$nFactors
   if (all == TRUE) {
    cat(paste("-- repsim = ", i, "/",repsim,"\n", sep=""))
    components[i, (15)] <- nBentler(x=eigenvalues, N=N, alpha=1-quantile, cor=TRUE)$nFactors
    }
   }
  # analysis       <- rbind(analysis, results$Analysis)
  #components[2,] <- t(results$Components);components
  }

 names                <- colnames(results$Components)
 names                <- c("oc", "af", "par", "mean.eig", "per")
 components           <- data.frame(components)
 colnames(components) <- c(names,"cng","b","t.b","p.b","sescree","R2","Bartlett","Anderson","Lawley","Bentler")
 if (details == TRUE) analysis <- list(components=components, eigenvalues=values)
 if (repsim > 1)      components <- moreStats(components, quantile=quantile) else components <- NA

 res <- list(details=analysis, nFactors=components)
 class(res) <- 'structureSim'
 return(res)
 }
## LIGNE 21 MODIFIEE: ETAIT quantile=0.95
## LIGNE 42 MODIFIEE: EAIT c("oc", "af", "par", "per", "mean.eig")
