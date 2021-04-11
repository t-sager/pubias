EstimatingSelection <- function(X,sigma,symmetric,cluster_ID,cutoffs,Studynames, identificationapproach, GMMapproach) {
    GMM_obj <- function(Psi){
        MetastudyGMMObjective(Psi, cutoffs, symmetric, X, sigma, cluster_ID)+max(-min(Psi),0)*10^5
    }
    
    mini<-optim(par=Psihat0,fn=GMM_obj,method="BFGS",control = list(abstol=10^-8,maxit=10^5));
    
    #For more accurate optimization, use the following:
    #Psihat1=mini$par;
    #Objval=mini$value;
    
    #mini<-optim(par=Psihat1,fn=GMM_obj);
    
    Psihat=mini$par;
    Objval=mini$value;
    
    GMM_obj_new<-function(Psi){
        #all_in_one=cbind(Psi,Psihat)
        ComputingMetastudyGMMObjective(Psi,Psihat,cutoffs,symmetric,X,sigma,cluster_ID)
    }
    
    mini<-optim(par=Psihat0_theta,fn=GMM_obj_new,method="BFGS",control = list(abstol=10^-8,maxit=10^5));
    
    #For more accurate optimization, use the following:
    #Psihat1_theta=mini$par;
    #Objval_theta=mini$value;
    
    #mini<-optim(par=Psihat1_theta,fn=GMM_obj_new);
    
    Psihat_theta=mini$par;
    Objval_theta=mini$value;
    
    Psihat_theta[2]=abs(Psihat_theta[2])
    
    Psihat_all = cbind(Psihat_theta,Psihat);
    
    mom=ComputingMetastudyMoments(cbind(Psihat_all,1),cutoffs,symmetric,X,sigma)
    moments=mom$moment_mean;
    Sigma_temp=mom$moment_var;
    rhat=mom$raw_moments;
    
    Sigma_hat=Clustered_covariance_estimate(rhat,cluster_ID);
    
    stepsize=10^-3;
    G=matrix(0,dim(moments)[2],length(Psihat_all));
    for (n1 in c(1:length(Psihat_all))) {
        beta_plus=Psihat_all;
        beta_plus[n1]=beta_plus[n1]+stepsize
        mom=ComputingMetastudyMoments(c(beta_plus,1),cutoffs,symmetric,X,sigma);
        moments_plus=mom$moment_mean;
        beta_minus=Psihat_all;
        beta_minus[n1]=beta_minus[n1]-stepsize
        mom=ComputingMetastudyMoments(cbind(beta_minus,1),cutoffs,symmetric,X,sigma);
        moments_minus=mom$moment_mean;
        G[,n1]=t(moments_plus-moments_minus)/(2*stepsize);
    }
    
    #In any case you get an error on an almost numerically singular matrix,
    #try using ``MASS::ginv'' instead of ``solve''. This command will run a generalized
    #inverse matrix and usually deals with this error.
    
    Varhat_all=solve(t(G)%*%solve(Sigma_hat)%*%G)/length(X);
    se_robust=sqrt(diag(Varhat_all));
    dof=dim(moments)[2];
    Varhat=Varhat_all[3:dim(Varhat_all)[1],3:dim(Varhat_all)[2]];
    
    if (length(cutoffs)==1) {
        Psi_grid=c(seq(0.001,5,0.001),10,10^5)
        S_store=matrix(0,length(Psi_grid),1)
        for (m in c(1:length(Psi_grid))) {
            S_store[m,1]=GMM_obj(Psi_grid[m])
        }
        
        CS_LB=min(Psi_grid[S_store<qchisq(0.95,1)]);
        CS_UB=max(Psi_grid[S_store<qchisq(0.95,1)]);
    } else if (length(cutoffs)==2) {
        gridsteps=200;
        Psi_grid1=seq(0.01,0.25,(0.25-0.01)/(gridsteps-1));
        Psi_grid2=seq(0.05,5,(5-0.05)/(gridsteps-1));
        S_store=matrix(0,length(Psi_grid1),length(Psi_grid2));
        for (m1 in c(1:length(Psi_grid1))) {
            for (m2 in c(1:length(Psi_grid2))) {
                Psi_temp=c(Psi_grid1[m1],Psi_grid2[m2]);
                S_store[m1,m2]=GMM_obj(Psi_temp);
            }
            
        }
        Psi_grid=cbind(t(Psi_grid1),t(Psi_grid2));
    } else {
        Psi_grid=c(seq(0.001,5,0.001),10,10^5)
        S_store=matrix(0,length(Psi_grid),1)
    }
    
    Psi_grid=as.matrix(Psi_grid)
    
    # Ab hier SelectionTableGMM
    sink(here("FiguresandTables","GMMEstimatesAllSelectionModel.tex"))
    columns_names<-c("$\\Theta$","$\\Sigma$")
    if (length(Psihat)==1) {
        columns_names<-c(columns_names,"$\\beta_{p}$")
    } else {
        for (i in c(1:length(Psihat))) {
            columns_names<-c(columns_names,paste0("$\\beta_{p,",i,"}$"))
        }
        
    }
    results<-matrix(rbind(format(Psihat_all,digits=2),format(se_robust,digits=2)),nrow=2,dimnames = list(c(1:2),columns_names))
    results[2,]<-paste0("(",results[2,],")")
    results<-as.matrix(results)
    alignment <- c("c")
    for (i in c(1:(length(Psihat_all)))) {
        alignment<-c(alignment,"c")
    }
    print(xtable(results,align = alignment), sanitize.rownames.function = function(a) {''}, sanitize.colnames.function = function(a) {a})
    sink()
    
    sink(here("FiguresandTables","GMMresultsSelectionModel.tex"))
    if (length(Psihat)>1) {
        columns_names<-c("$\\beta_{p,1}$")
        for (i in c(2:length(Psihat))) {
            columns_names<-c(columns_names,paste0("$\\beta_{p,",i,"}$"))
        }
    } else {
        columns_names<-c("$\\beta_{p}$")
    }
    
    results<-matrix(rbind(format(Psihat,digits=2),format(se_robust[3:length(se_robust)],digits=2)),nrow=2,dimnames = list(c(1:2),columns_names))
    results[2,]<-paste0("(",results[2,],")")
    results<-as.matrix(results)
    alignment <- c("c")
    if (length(Psihat)>=1) {
        for (i in c(1:(length(Psihat)))) {
            alignment<-c(alignment,"c")
        }
    }
    
    print(xtable(results,align = alignment), sanitize.rownames.function = function(a) {''}, sanitize.colnames.function = function(a) {a})
    sink()
    
    dof=length(moments[3:length(moments)])
    len=length(Psihat)
    if (dorobust==1) {
        if (len==1) {
            CS_LB = min(Psi_grid[S_store<qchisq(0.95,dof)])
            CS_UB = max(Psi_grid[S_store<qchisq(0.95,dof)])
            
            if (CS_LB==min(Psi_grid)) {
                CS_LB=0
            }
            if (CS_UB==max(Psi_grid)) {
                CS_UB=c("$\\infty$")
            }
            sink(here("FiguresandTables","GMMresultsSelectionModelCS.tex"))
            columns_CS_names<-c(paste0("$\\beta_p$", " Lower Bound"),paste0("$\\beta_p$"," Upper Bound"))
            CS_bounds<-matrix(c(CS_LB,CS_UB),nrow=1,dimnames = list(c(1),columns_CS_names))
            alignment <- c("c","c","c")
            print(xtable(CS_bounds,align = alignment), sanitize.text.function = function(a) {a}, sanitize.rownames.function = function(a) {''}, sanitize.colnames.function = function(a) {a})
            sink()
        } else if (len==2) {
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
    return(list("Psihat"= Psihat, "Varhat" = Varhat))
}


