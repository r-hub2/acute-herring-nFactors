#' Scree Plot According to a nScree Object Class
#'
#' Plot a scree plot adding information about a non graphical \code{nScree}
#' analysis.
#'
#'
#' @param nScree Results of a previous \code{nScree} analysis
#' @param legend Logical indicator of the presence or not of a legend
#' @param xlab Label of the x axis (default to \code{"Component"})
#' @param ylab Label of the y axis (default to \code{"Eigenvalue"})
#' @param main Main title (default to \code{"Non Graphical Solutions to the
#' Scree Test"})
#' @return Nothing returned.
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{plotuScree}}, \code{\link{nScree}},
#' \code{\link{plotParallel}}, \code{\link{parallel}}
#' @references
#' Raiche, G., Walls, T. A., Magis, D., Riopel, M. and Blais, J.-G. (2013). Non-graphical solutions
#' for Cattell's scree test. Methodology, 9(1), 23-29.
#' @export
#' @importFrom graphics lines par text plot.default
#' @importFrom stats lm coef
#' @keywords Graphics
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## INITIALISATION
#'  data(dFactors)                      # Load the nFactors dataset
#'  attach(dFactors)
#'  vect         <- Raiche              # Use the second example from Buja and Eyuboglu
#'                                      # (1992, p. 519, nsubjects not specified by them)
#'  eigenvalues  <- vect$eigenvalues    # Extract the observed eigenvalues
#'  nsubjects    <- vect$nsubjects      # Extract the number of subjects
#'  variables    <- length(eigenvalues) # Compute the number of variables
#'  rep          <- 100                 # Number of replications for the parallel analysis
#'  cent         <- 0.95                # Centile value of the parallel analysis
#'
#' ## PARALLEL ANALYSIS (qevpea for the centile criterion, mevpea for the mean criterion)
#'  aparallel    <- parallel(var     = variables,
#'                           subject = nsubjects,
#'                           rep     = rep,
#'                           cent    = cent)$eigen$qevpea  # The 95 centile
#'
#' ## NOMBER OF FACTORS RETAINED ACCORDING TO DIFFERENT RULES
#'  results <- nScree(eig       = eigenvalues,
#'                    aparallel = aparallel
#'                    )
#'
#'  results
#'
#' ## PLOT ACCORDING TO THE nScree CLASS
#'  plotnScree(results)
#'  }
#' }
"plotnScree" <-
function (nScree,
          legend = TRUE,
          ylab   = "Eigenvalues",
          xlab   = "Components",
          main   = "Non Graphical Solutions to Scree Test")
          {
   if (!inherits(nScree, "nScree"))  stop("Method is only for nScree objects")
   #if (!exists("legend", mode="logical") ) legend <- TRUE                                   # To develop
   #if (!exists("ylab"))                    ylab <- "Eigenvalues"                            # To develop
   #if (!exists("xlab"))                    xlab <- "Components"                             # To develop
   #if (!exists("main"))                    main <- "Non Graphical Solutions to Scree Test"  # To develop
   if (nScree$Model == "components") nkaiser = "Eigenvalues (>mean  = " else nkaiser = "Eigenvalues (>0 = "
   if (nScree$Model == "factors")  xlab   = "Factors"
   graphics::par(col   = 1, pch = 1)     # Color and symbol for usual scree
   graphics::par(mfrow = c(1,1))
   eig        <- nScree$Analysis$Eigenvalues
   k          <- 1:length(eig)
   #plotuScree(x=eig, ...)                                                                   # To develop
   plotuScree(x=eig, main=main, xlab=xlab, ylab=ylab)
   nk         <- length(eig)
   noc        <- nScree$Components$noc
   vp.p       <- stats::lm(eig[c(noc+1,nk)] ~ k[c(noc+1,nk)])
   x          <- sum(c(1,1) * stats::coef(vp.p))
   y          <- sum(c(1,nk)* stats::coef(vp.p))
   graphics::par(col = 10)            # Color for optimal coordinates
   graphics::lines(k[c(1,nk)],c(x,y))
   graphics::par(col = 11,pch=2)            # Color and symbol for parallel analysis
   graphics::lines(1:nk, nScree$Analysis$Par.Analysis, type = "b")
   if (legend == TRUE) {
     leg.txt  <- c(paste(nkaiser,nScree$Components$nkaiser,")"),
                 c(paste("Parallel Analysis (n = ",nScree$Components$nparallel,")")),
                 c(paste("Optimal Coordinates (n = ",nScree$Components$noc,")")),
                 c(paste("Acceleration Factor (n = ",nScree$Components$naf,")")) )
     legend("topright",
            legend   = leg.txt,
            pch      = c(1,2,NA,NA),
            text.col = c(1,3,2,4), col = c(1,3,2,4)
            )
     }
   naf        <-   nScree$Components$naf
   graphics::text(x = noc ,    y = eig[noc],     label = " (OC)", cex = .70, adj = c(0,0), col = 2)
   graphics::text(x = naf + 1, y = eig[naf + 1], label = " (AF)", cex = .70, adj = c(0,0), col = 4)
   }

