Clustered_covariance_estimate <- function(g,cluster_index) {
  #given a matrix of moment condition values g, compute a clustering-robust
  #estimate of the covariance matrix Sigma
  sorted=sort(cluster_index,decreasing=FALSE,index.return=TRUE)
  cluster_index=as.matrix(sorted$x)
  I=as.matrix(sorted$ix)
  g=as.matrix(g[I,])
  g=g-matrix(rep(colMeans(g),dim(g)[1]),ncol=dim(g)[2],byrow=TRUE)
  gsum=matrix(cumsum(g),dim(g)[1],dim(g)[2])
  index_diff=as.matrix(ifelse(cluster_index[2:length(cluster_index),1]!=cluster_index[1:(length(cluster_index)-1),1],1,0))
  index_diff=rbind(index_diff,1)
  gsum=subset(gsum,index_diff==1)
  gsum=rbind(gsum[1,],diff(gsum))
  Sigma=(1/(dim(g)[1]-1))*(t(gsum)%*%gsum)
  return(Sigma)
}
