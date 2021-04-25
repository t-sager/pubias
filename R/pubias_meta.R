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
pubias_meta <- function(data, studynames, symmetric = 1, symmetric_p = 1, cutoffs = 1.96, GMM = FALSE){

  # Preliminaries
  X <<- as.matrix(data[,1])
  sigma <- as.matrix(data[,2])
  cluster_ID <<- as.matrix(data[,3])
  n <- nrow(X)
  C <- matrix(1,n,1)
  includeinfigure <<- as.logical(rep(1,n))
  includeinestimation <<- as.logical(rep(1,n))
  identificationapproach <- 2

  if (GMM == TRUE) {
    name <- 'GMM_Meta'
    result <<- gmm_meta(X, sigma, symmetric, cluster_ID, cutoffs, studynames)
    descriptive_stats(X, sigma, identificationapproach, name, symmetric, cluster_ID)
    corrected_estimates <<- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p,identificationapproach,GMM)
    plot_correction(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,studynames,identificationapproach, corrected_estimates)
  } else {
    name <- 'MLE_Meta'
    result <<- mle_meta(X, sigma, symmetric, symmetric_p, cluster_ID, cutoffs, studynames, C)
    descriptive_stats(X, sigma, identificationapproach, name, symmetric, cluster_ID)
    corrected_estimates <<- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p,identificationapproach,GMM)
    plot_correction(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,studynames,identificationapproach, corrected_estimates)
  }
}
