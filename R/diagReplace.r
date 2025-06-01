#' Replacing Upper or Lower Diagonal of a Correlation or Covariance Matrix
#'
#' The \code{diagReplace} function returns a modified correlation or covariance
#' matrix by replacing upper diagonal with lower diagonal, or lower diagonal
#' with upper diagonal.
#'
#' @param R numeric: correlation or covariance matrix
#' @param upper logical: if \code{TRUE} upper diagonal is replaced with lower
#' diagonal. If \code{FALSE}, lower diagonal is replaced with upper diagonal.
#' @return \item{R }{ numeric: correlation or covariance matrix }
#' @author Gilles Raiche \cr Centre sur les Applications des Modeles de
#' Reponses aux Items (CAMRI) \cr Universite du Quebec a Montreal\cr
#' \email{raiche.gilles@@uqam.ca}
#' @keywords manip
#' @export
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
#' # Replace upper diagonal with lower diagonal
#'  RU <- diagReplace(R, upper=TRUE)
#'
#' # Replace lower diagonal with upper diagonal
#'  RL <- diagReplace(R, upper=FALSE)
#' # .......................................................
#'  }
#' }

diagReplace <- function(R, upper=TRUE) {
 RT <- R
 if (upper == TRUE) {
  Rtranspose         <- t(RT)
  # Replacing upper diagonal with lower diagonal
  RT[upper.tri(RT)]  <- Rtranspose[upper.tri(Rtranspose)]
  return(RT)
  }
 if (upper == FALSE) {
  Rtranspose         <- t(RT)
  # Replacing lower diagonal with upper diagonal
  RT[lower.tri(RT)]  <- Rtranspose[lower.tri(Rtranspose)]
  return(RT)
  }
 }


