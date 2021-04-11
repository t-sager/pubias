SelectionTable<-function(pathname,Psihat,se_robust,name,notilde,symmetric) {
  FID = sink(paste0(pathname,"/FiguresandTables/", name, "SelectionModel.tex"));
  
  
  if (symmetric==0) {
    len<-length(Psihat)-2;

  cat(sprintf( paste0("\\begin{tabular}{cc |", rep("c", 1, len),"}")),"\n");
  cat(sprintf( "$\\bar \\theta$ & "));
  } else {
    len=length(Psihat)-1;
    cat(sprintf( paste0("\\begin{tabular}{c |", rep("c", 1, len),"}")),"\n");
  }
   

  
  
  if (notilde) {
    cat(sprintf( "$\\tau$ "));
  } else {
    cat(sprintf( "$\\tilde \\tau$ "));
  }
   
 
  
  if (len==1) {
    cat(sprintf( "& $\\beta_p$ \\\\ \n \\hline \n"));
  } else {
    for (j in 1:len) {
      cat(sprintf( "& $\\beta_{p,%i}$ ", j));
    }
        cat(sprintf( "\\\\ \n \\hline \n"));
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
  sink(FID);
  
  
}
 
#if (test) {
#  Psihat<-c(1,2)
 # se_robust<-c(1,2)
 # name<-"Test"
 # symmetric<-1
 # notilde<-FALSE
 # pathname<-"/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication"
  
 # SelectionTable(pathname,Psihat,se_robust,name,notilde,symmetric)
#} 
