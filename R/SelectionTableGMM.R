SelectionTableGMM <- function(pathname,Psihat,se_robust,name,S,Psigrid,dof){
     FID = sink(paste0(pathname,"/FiguresandTables/", name, "SelectionModel.tex"));
     
    len <- length(Psihat)
    cat(sprintf( paste0("\\begin{tabular}{cc |", rep("c", 1, len),"}")),"\n");

#################    

if (len==1){
    cat(sprintf( "& $\\beta_p$ \\\\ \n \\hline \n"));
} else {
    cat(sprintf( "& $\\beta_{p,1}$ "));
    for (j in 2:len){
        cat(sprintf( "& $\\beta_{p,%i}$ ", j));
    }
    cat(sprintf( "\\\\ \n \\hline \n"));
}

################# 
    
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
    
#################    

if (len==1){
    CS_LB <- min(Psigrid(S<qchisq(0.95,dof),1,0))
    CS_UB <- max(Psigrid(S<qchisq(0.95,dof),1,0))

    FID = sink(paste0(pathname,"/FiguresandTables/", name, "SelectionModelCS.tex"));

    cat(sprintf( paste0("\\begin{tabular}{cc |", rep("c", 1, len),"}")),"\n");
    cat(sprintf( "$\\beta_p$ Lower Bound "));
    cat(sprintf( "& $\\beta_p$ Upper Bound "));
    cat(sprintf( "\\\\ \n \\hline \n "));

    if (CS_LB==min(Psigrid)){
        cat(sprintf( "%.3f ", 0))
    } else {
        cat(sprintf( "%.3f ", CS_LB))
    }

    if (CS_UB==max(Psigrid)){
        cat(sprintf( "& $\\infty$ "))
    } else {
        cat(sprintf( "%.3f ", CS_UB))
    }
    cat(sprintf( "\\\\  \n"));
    cat(sprintf( "\\end{tabular}\n"));
    sink(FID);

} else if (len==2){
        color<-brewer.pal(8,"Greys");
        GREY=color[3];
        S_CS = ifelse(S_store<=qchisq(0.95,dof),1,0)
        S_CS=as.matrix(S_CS)
        S_CS=t(S_CS)
        pdf(here("FiguresandTables","GMMresultsCS.pdf"), width = 4, height = 4)
        contour(Psi_grid[,1],Psi_grid[,2],S_CS,xlab="$\\beta_{1}$",ylab="$\\beta_{2}$");
        dev.off()
    }

}
