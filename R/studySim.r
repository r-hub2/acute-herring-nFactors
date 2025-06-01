#' Simulation Study from Given Factor Structure Matrices and Conditions
#'
#' The \code{structureSim} function returns statistical results from
#' simulations from predefined congeneric factor structures. The main ideas
#' come from the methodology applied by Zwick and Velicer (1986).
#'
#'
#' @param var numeric: vector of the number of variables
#' @param nFactors numeric: vector of the number of components/factors
#' @param pmjc numeric: vector of the number of major loadings on each
#' component/factor
#' @param loadings numeric: vector of the major loadings on each
#' component/factor
#' @param unique numeric: vector of the unique loadings on each
#' component/factor
#' @param N numeric: vector of the number of subjects/observations
#' @param repsim numeric: number of replications of the matrix correlation
#' simulation
#' @param reppar numeric: number of replications for the parallel and
#' permutation analysis
#' @param stats numeric: vector of the statistics to return: mean(1),
#' median(2), sd(3), quantile(4), min(5), max(6)
#' @param quantile numeric: quantile for the parallel and permutation analysis
#' @param model character: \code{"components"} or \code{"factors"}
#' @param r2limen numeric: R2 limen value for the R2 Nelson index
#' @param all logical: if \code{TRUE} computes the Bentler and Yuan index (very
#' long computing time to consider)
#' @param dir character: directory where to save output. Default to NA
#' @param trace logical: if \code{TRUE} outputs details of the status of the
#' simulations
#' @return \item{values}{ Returns selected statistics about the number of
#' components/factors to retain: mean, median, quantile, standard deviation,
#' minimum and maximum.}
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{generateStructure}}, \code{\link{structureSim}}
#' @references
#' Raiche, G., Walls, T. A., Magis, D., Riopel, M. and Blais, J.-G. (2013). Non-graphical solutions
#' for Cattell's scree test. Methodology, 9(1), 23-29.
#'
#' Zwick, W. R. and Velicer, W. F. (1986). Comparison of five rules
#' for determining the number of components to retain. \emph{Psychological
#' Bulletin, 99}, 432-442.
#' @export
#' @keywords multivariate
#' @examples
#'
#' \dontrun{
#' # ....................................................................
#' # Example inspired from Zwick and Velicer (1986)
#' # Very long computimg time
#' # ...................................................................
#'
#' # 1. Initialisation
#' # reppar    <- 30
#' # repsim    <- 5
#' # quantile  <- 0.50
#'
#' # 2. Simulations
#' # X         <- studySim(var=36,nFactors=3, pmjc=c(6,12), loadings=c(0.5,0.8),
#' #                       unique=c(0,0.2), quantile=quantile,
#' #                       N=c(72,180), repsim=repsim, reppar=reppar,
#' #                       stats=c(1:6))
#'
#' # 3. Results (first 10 results)
#' # print(X[1:10,1:14],2)
#' # names(X)
#'
#' # 4. Study of the error done in the determination of the number
#' #    of components/factors. A positive value is associated to over
#' #    determination.
#' # results   <- X[X$stats=="mean",]
#' # residuals <- results[,c(11:25)] - X$nfactors
#' # BY        <- c("nsubjects","var","loadings")
#' # round(aggregate(residuals, by=results[BY], mean),0)
#'  }
#'
studySim <- function(var, nFactors, pmjc, loadings, unique, N, repsim, reppar,
                     stats=1, quantile=0.5, model="components", r2limen=0.75,
                     all=FALSE, dir=NA, trace=TRUE) {
 nsubjects <- N
 result    <- NULL
 id        <- 0
 nid       <- length(nFactors) * length(loadings) * length(pmjc) * length(var) * length(unique) * length(nsubjects)
 for (i in 1:length(nFactors))  {
  for (j in 1:length(loadings))  {
   for (l in 1:length(pmjc))      {
    for (n in 1:length(var))       {
     for (k in 1:length(unique))    {
      for (m in 1:length(nsubjects)) {
       id    <- id + 1
       kid  <- paste(id,"/",nid,sep="")
       ident <- c(nFactors=nFactors[i], loadings=loadings[j], unique=unique[k], quantile=quantile,
                  pmjc=pmjc[l], nsubjects=nsubjects[m], var=var[n], reppar=reppar,
                  repsim=repsim, id=kid, model=model)
       if (trace == TRUE) print(ident)
       fStruct <- generateStructure(var=var[n], mjc=nFactors[i], pmjc=pmjc[l], loadings=loadings[j], unique=unique[k])
       fSim    <- structureSim(fload=as.matrix(fStruct), reppar=reppar, repsim=repsim, details=FALSE, all=all,
                               N=nsubjects[m], quantile=quantile, model=model, r2limen=r2limen)[[2]][stats,]
       if (length(stats) == 1) {
        fSim  <- data.frame(var=var[n], nsubjects=nsubjects[m], nfactors=nFactors[i], pmjc=pmjc[l],
                            loadings=loadings[j], unique=unique[k], t(fSim), repsim=repsim, reppar=reppar)
        }
       if (length(stats) > 1) {
        ls    <- length(stats)
        info  <-  data.frame(stats   =rownames(fSim),       id       =rep(id, ls),
                             var     =rep(var[n], ls),      nsubjects=rep(nsubjects[m], ls),
                             nfactors=rep(nFactors[i], ls), pmjc     =rep(pmjc[l], ls),
                             loadings=rep(loadings[j], ls), unique   =rep(unique[k], ls),
                             repsim  =rep(repsim, ls),      reppar   =rep(reppar, ls))
         fSim <- data.frame(info, fSim)
         }
        result           <- rbind(result, fSim)
        rownames(result) <- 1:dim(result)[1]
        fString          <- paste("RES_", paste(ident,"_", sep="", collapse=""), sep="")
        # if (!is.na(dir)) save("fSim", file=paste(dirPack, fString,".Rdata", sep=""))  # Old erroneous code
        if (!is.na(dir)) save("fSim", file=paste(dir, fString,".Rdata", sep=""))
  }}}}}}
 return(result)
 }


