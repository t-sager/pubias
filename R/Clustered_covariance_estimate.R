#' Title
#'
#' @param g
#' @param cluster_ID
#'
#' @return
#' @export
#'
#' @examples
clustered_covariance_estimate <- function(g,cluster_ID) {
  sorted <- sort(cluster_ID,decreasing=FALSE,index.return=TRUE)
  cluster_ID <- as.matrix(sorted$x)
  I <- as.matrix(sorted$ix)
  g <- as.matrix(g[I,])
  g <- g-matrix(rep(colMeans(g),nrow(g)),ncol=ncol(g),byrow=TRUE)
  gsum <- matrix(cumsum(g),nrow(g),ncol(g))
  index_diff <- as.matrix(ifelse(cluster_ID[2:length(cluster_ID),1]!=cluster_ID[1:(length(cluster_ID)-1),1],1,0))
  index_diff <- rbind(index_diff,1)
  gsum <- subset(gsum,index_diff==1)
  gsum <- rbind(gsum[1,],diff(gsum))
  Sigma <- (1/(nrow(g)-1))*(t(gsum)%*%gsum)
  return(Sigma)
}
