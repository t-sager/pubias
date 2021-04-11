#Potential issues: Z^r, sigma
# Figure size
# Stopped at histogram

DescriptiveStats<-function(X, sigma, identificationapproach,name,symmetric,cluster_ID) {
    #% for approach 1:
    #  % X has 2 columns, initial estimate Zstat, and replication estimate
    #  % normalized by sigma1; sigma is SE2/SE1
    # % for approach 2:
    #   % X has 1 column, initial estimate unscaled, sigma=sigma1
    
    # Significance
    critval<-1.96;
    
    if (identificationapproach==1) {
        ll<-floor(min(min(X)));
        uu<-ceiling(max(max(X)));
        significant<-(abs(X[,1])>critval);
        dat<-data.frame(xvar = X[,1],
                        yvar = X[,2],
                        significant = significant)
        
        # Plot scatter and scatterslides (with letters A,B)    
        p<- ggplot(dat, aes(x=xvar, y=yvar)) +geom_vline(xintercept =critval,color="grey")+
            geom_hline(yintercept =critval,color="grey")+
            geom_abline(intercept = 0,slope=1,color="grey")+
            xlab("Z")+
            ylab("ZRep" )+
            xlim(c(max(ll,0),uu))+
            ylim(c(ll,uu))+
            geom_point(shape=21,size = 2,aes(fill = significant))
        filepath<-paste0(pathname,'/FiguresandTables/',name, 'Scatter.pdf')
        pdf(filepath, width=5, height=5)
        print(p)
        dev.off()
        q<-p+ annotate("text", x = 4, y = 1, label = "A")+
            annotate("text", x = 1, y = 4, label = "B")
        filepath<-paste0(pathname,'/FiguresandTables/',name, 'SlidesScatter.pdf')
        pdf(filepath, width=5, height=5)
        print(q)
        dev.off()
        
    } else  if (identificationapproach==2)  {
        significant<-(abs(X/sigma)>critval);
        nooutlier<-sigma<50;
        dat<-data.frame(xvar=X,
                        yvar=sigma,
                        significant = significant &nooutlier )
        rangeX=1.1*max(max(abs(X)), max(abs(sigma[nooutlier]))*critval);
        p<-ggplot(dat, aes(x=xvar,y=yvar)) +
            xlim(c(-rangeX,rangeX))+
            ylim(c(min(sigma)-0.03,rangeX/critval))+
            xlab("X")+
            ylab(expression("sigma" ))+
            geom_abline(intercept = 0,slope=1/critval,color="grey")+ 
            geom_abline(intercept = 0,slope=-1/critval,color="grey")+
            geom_point(shape=21,size = 2,aes(fill = significant))
        
        filepath<-paste0(pathname,'/FiguresandTables/',name, 'Scatter.pdf')
        pdf(filepath, width=5, height=5)
        print(p)
        dev.off()
        
    }
    # Create histogram plot for both identifiactionapproach = 1 and 2
    
    if (identificationapproach==1) {
        Zuse = Z[,1];
        
    } else {
        Zuse = X/sigma;
        if (max(abs(Zuse))>10) {
            Zuse<-Zuse[abs(Zuse)<6];
        } 
    }
    ll=floor(min(Zuse));
    lleven=floor(ll/2)*2;
    uu=ceiling(max(Zuse));
    
    if (symmetric == 0) {
        if (n>=30) {
            uu2<-ceiling((uu-.36)/.32)*.32+.36;
            ll2<-floor((ll+.36)/.32)*.32-.36;
            edges<-c(seq(from=ll2,
                         to=-0.36,
                         by=0.32), 0, seq(from=0.36,
                                          to=uu2,
                                          by=0.32));
        } else {
            uu2<-ceiling((uu-.68)/.64)*.64+.68;
            ll2<-floor((ll+.68)/.64)*.64-.68;
            edges<-c(seq(from=ll2,
                         to=-0.68,
                         by=0.64), 0, seq(from=0.68,
                                          to=uu2,
                                          by=0.64));
        }
    } else if (symmetric==1){
        if (n>=30) {
            uu2=ceiling((uu-.36)/.32)*.32+.36;
            edges = c(0,seq(from=.36,
                            to=1.64,
                            by=0.32),
                      seq(from=1.96,
                          to=uu2,
                          by=0.32)
            )
        }  else {
            uu2=ceiling((uu-.68)/.64)*.64+.68;
            edges= c(0, seq(from=.68,
                            to=uu2,
                            by=0.64))
        }
        
    }
    
    
    
    if (symmetric == 0) {
        
        h<-ggplot(data = as.data.frame(Zuse), aes(Zuse))+
            geom_histogram(aes(y = ..density..),
                           fill = 'blue',
                           breaks=edges)+
            geom_vline(xintercept =-1.96,color='grey')+
            geom_vline(xintercept =1.96, color='grey')+
            xlab('X/sigma')+
            ylab('Density')+
            xlim(c(min(edges),max(edges)))
        
        
    } else {
        h<-ggplot(data = as.data.frame(Zuse), aes(Zuse))+
            geom_histogram(aes(y = ..density..),
                           fill = 'blue',
                           breaks=edges)+
            geom_vline(xintercept =1.96, color='grey')+
            xlab('|X|/sigma')+
            ylab('Density')+
            xlim(c(min(edges),max(edges)))
    }
    # If  histogram is for both id approachaes, just remove if clause below
    
    if ( identificationapproach == 2) {
        # handp<-multiplot(h,p,cols=2)
        handp<-grid.arrange(h, p, ncol = 2)
        filepath<-paste0(pathname,'/FiguresandTables/',name, 'ScatterHist.pdf')
        save_plot(filepath,handp,ncol=2,base_width = 4, base_height=3)
    }
    
    
    
    
    
    #  %Meta-regression estimates
    if (identificationapproach==2 && symmetric==0) {
        source("MetaRegressionTable.R")
        source("Clustered_covariance_estimate.R")
        do_metatable<-function(sigma,Zstats) {
            R = cbind(rep(1,length(X)),sigma);
            betahat<-solve(t(R)%*%R)%*%t(R)%*%X
            ehat<-X-R%*%betahat
            
            Sigma_base<-R*matrix(rep(ehat,dim(R)[2]),ncol=dim(R)[2],byrow=FALSE)
            Sigma<-Clustered_covariance_estimate(Sigma_base,cluster_ID);
            Vhat<-n*solve((t(R)%*%R))%*%Sigma%*%solve((t(R)%*%R))
            se_robust<-sqrt(diag(Vhat))
            
            betahat<-as.vector(betahat)
            MetaRegressionTable(pathname,betahat,se_robust,name,Zstats)
        }
        
        do_metatable(sigma,0)
        
        do_metatable(sigma^(-1),1)
        
    }
    
    
}

