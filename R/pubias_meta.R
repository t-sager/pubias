#' Title
#'
#' @param data
#' @param studynames
#' @param symmertic
#' @param cutoffs
#' @param GMM
#'
#' @return
#' @export
#'
#' @examples
pubias_meta <- function(data, studynames, symmetric = 1, symmetric_p = 1,cutoffs = 1.96, GMM = FALSE){

  # Preliminaries
  X <- as.matrix(data[,1])
  sigma <- as.matrix(data[,2])
  cluster_ID <- as.matrix(data[,3])
  n <- length(X)
  C <- matrix(1,n,1)
  includeinfigure <<- as.logical(rep(1,n))
  includeinestimation <- as.logical(rep(1,n))
  identificationapproach <- 2

  if (GMM == TRUE) {
    name <- 'GMM'
    gmm_meta(X, sigma, symmetric, cluster_ID, cutoffs, studynames)
  } else {
    name <- 'MLE'
    mle_meta(X, sigma, symmetric, symmetric_p, cluster_ID, cutoffs, studynames, includeinestimation)
  }

}
