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
pubias_replication <- function(data, studynames, symmetric = 1, cutoffs = 1.96, GMM = FALSE){

  # Preliminaries
  Z <- cbind(data[,1]/data[,2], data[,3]/data[,2])
  sigmaZ2 <- as.matrix(data[,4]/data[,2])
  X <- as.matrix(data[,1])
  sigma <- as.matrix(data[,2])
  cluster_ID <- as.matrix(data[,3])
  n <- length(X)
  C <- matrix(1,n,1)
  includeinfigure <<- as.logical(rep(1,n))
  includeinestimation <- as.logical(rep(1,n))
  identificationapproach <- 1

  if (GMM == TRUE) {
    name <- 'GMM'
    gmm_replication(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, studynames)
  } else {
    name <- 'MLE'
    mle_replication(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, studynames)
  }

}
