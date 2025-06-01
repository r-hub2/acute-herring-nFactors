#' Multiple Regression Procedure to Determine the Number of Components/Factors
#'
#' This function computes the \eqn{\beta} indices, like their associated
#' Student \emph{t} and probability (Zoski and Jurs, 1993, 1996, p. 445). These
#' three values can be used as three different indices for determining the
#' number of components/factors to retain.
#'
#' When the associated Student \emph{t} test is applied, the following
#' hypothesis is considered: \cr
#'
#' (1) \eqn{\qquad \qquad H_k: \beta (\lambda_1 \ldots \lambda_k) - \beta
#' (\lambda_{k+1} \ldots \lambda_p), (k = 3, \ldots, p-3) = 0} \cr
#'
#'
#' @param x numeric: a \code{vector} of eigenvalues, a \code{matrix} of
#' correlations or of covariances or a \code{data.frame} of data (eigenFrom)
#' @param cor logical: if \code{TRUE} computes eigenvalues from a correlation
#' matrix, else from a covariance matrix
#' @param model character: \code{"components"} or \code{"factors"}
#' @param details logical: if \code{TRUE} also returns details about the
#' computation for each eigenvalue.
#' @param ...  variable: additionnal parameters to give to the
#' \code{eigenComputes} and \code{cor} or \code{cov} functions
#' @return \item{nFactors}{ numeric: number of components/factors retained by
#' the \emph{MREG} procedures. } \item{details}{ numeric: matrix of the details
#' for each indices.}
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @seealso \code{\link{plotuScree}}, \code{\link{nScree}},
#' \code{\link{plotnScree}}, \code{\link{plotParallel}}
#' @references Zoski, K. and Jurs, S. (1993). Using multiple regression to
#' determine the number of factors to retain in factor analysis. \emph{Multiple
#' Linear Regression Viewpoints, 20}(1), 5-9.
#'
#' Zoski, K. and Jurs, S. (1996). An objective counterpart to the visual scree
#' test for factor analysis: the standard error scree test.  \emph{Educational
#' and Psychological Measurement, 56}(3), 443-451.
#' @export
#' @importFrom stats sd lm pt
#' @keywords multivariate
#' @examples
#' \dontrun{
#' if(interactive()){
#' ## SIMPLE EXAMPLE OF A MREG ANALYSIS
#'
#'  data(dFactors)
#'  eig      <- dFactors$Raiche$eigenvalues
#'
#'  results  <- nMreg(eig)
#'  results
#'
#'  plotuScree(eig, main=paste(results$nFactors[1], ", ",
#'                             results$nFactors[2], " or ",
#'                             results$nFactors[3],
#'                             " factors retained by the MREG procedures",
#'                             sep=""))
#'  }
#' }
nMreg <-
function(x, cor=TRUE, model="components", details=TRUE, ...) {
 x       <- eigenComputes(x, cor=cor, model=model, ...)
 nlength <- 3
 detail  <- NULL
 n       <- length(x)
 if (n < 6) stop("The number of variables must be at least 6.")
 i       <- 1
 mreg    <- tmreg <- tmreg2 <-pmreg <- numeric(n-5)
 while (i <= (length(x)-5)) {
  xa        <- c(1:(i+2))
  ya        <- x[1:(i+2)]
  ma        <- stats::lm(ya ~ xa)
  Syx.a     <- stats::sd(ya)*sqrt((1-summary(ma)$r.squared) * ((length(ya)-1)/(length(ya)-2))) # Howell(2008, p. 253)
  compa     <- ma$coef[2]
  seCompa   <- summary(ma)$coef[2,2]

  xb        <- c((i+1+nlength):length(x))
  yb        <- x[(i+1+nlength):length(x)]
  mb        <- stats::lm(yb ~ xb)
  Syx.b     <- stats::sd(yb)*sqrt((1-summary(mb)$r.squared) * ((length(yb)-1)/(length(yb)-2))) # Howell(2008, p. 253)
  compb     <- mb$coef[2]
  seCompb   <- summary(mb)$coef[2,2]

  mreg[i]   <- compb - compa
  semreg    <- sqrt((Syx.a^2)/((length(xa)-1)*stats::sd(xa)^2) + (Syx.b^2)/((length(xb)-1)*stats::sd(xb)^2))     # Se_dif_b -> Howell(2008, p. 259, 266)
  tmreg[i]  <- (compb - compa)/(semreg)
  tmreg2[i] <- (mreg[i])/sqrt(seCompa^2 + seCompb^2) # Il semble, selon moi, qu'il y aurait une erreur dans la formule de Zoski et Just. Et ce serait la bonne formul, comme celle plu shaut, mais plus rapide de calcul.
  pmreg[i]  <- stats::pt(tmreg[i],(length(xa)-1) + (length(xb)-1) - 4, lower.tail=FALSE, log.p=TRUE)
  i         <- i + 1
  }
 if (details == TRUE) detail  <- data.frame(v=(1:(n-5)),values=x[1:(n-5)], mreg=mreg, tmreg=tmreg, pmreg=pmreg)
 mreg  <- as.numeric(which(mreg ==max( mreg, na.rm=TRUE)) + nlength)
 tmreg <- as.numeric(which(tmreg==max(tmreg, na.rm=TRUE)))
 pmreg <- as.numeric(which(pmreg==min(pmreg, na.rm=TRUE)))
 res   <- list(detail=detail, nFactors=c(b=mreg,t.p=tmreg,p.b=pmreg))
 class(res) <- c("nFactors","list")
 return(res)
 }
