MetaRegressionTable<-function(pathname,Psihat,se_robust,name,Zstats) {
  if (Zstats==0) {
    FID <-sink(paste0(pathname,"/FiguresandTables/" ,name, "MetaRegBeta.tex"));
  } else {
    FID <- sink(paste0(pathname,"/FiguresandTables/" ,name, "MetaRegGamma.tex"));
  }
  
  if (Zstats==0) {
    cat(sprintf("\\begin{tabular}{cc} \n"));
    cat(sprintf("$\\gamma_0$ & $\\gamma_1$ \n"));
  } else {
    cat(sprintf("\\begin{tabular}{cc} \n"));
    cat(sprintf("$\\beta_0$ & $\\beta_1$ \n"));

  }
  cat(sprintf( "%.3f ", Psihat[1]));
      for (j in (2:length(Psihat))) {
        cat(sprintf( "& %.3f ", Psihat[j]));
      }
      
  cat(sprintf( "\\\\  \n"));
  cat(sprintf( "(%.3f) ", se_robust[1]));
  for (j in (2:length(Psihat))) {
                cat(sprintf( "& (%.3f) ", se_robust[j]));
  }
              
  cat(sprintf( "\\\\  \n"));
  cat(sprintf( "\\end{tabular}\n"));
  sink();
}

#if (test) {
#  Psihat<-c(1,2)
#  se_robust<-c(1,2)
#  name<-"Test"
#  symmetric<-1
#  notilde<-FALSE
#  Zstats<-0
#  pathname<-"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication"
  
#  MetaRegressionTable(pathname,Psihat,se_robust,name,Zstats)
#} 
