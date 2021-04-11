DescriptiveStatsCombined <-function(X_rep,X_meta,sigma_meta,name,symmetric) {
  
  critval<-1.96
  significant<-(abs(X_rep[,1])>critval);
  ll<-floor(min(min(X_rep)));
  uu<-ceiling(max(max(X_rep)));
 
  dat<-data.frame(xvar = X_rep[,1],
                  yvar = X_rep[,2],
                  significant = significant)
  
  # Plot scatter and scatterslides (with letters A,B)    
  p1<- ggplot(dat, aes(x=xvar, y=yvar)) +geom_vline(xintercept =critval,color="grey")+
    geom_hline(yintercept =critval,color="grey")+
    geom_abline(intercept = 0,slope=1,color="grey")+
    xlab("W")+
    ylab("WRep" )+
    xlim(c(max(ll,0),uu))+
    ylim(c(ll,uu))+
    geom_point(shape=21,size = 2,aes(fill = significant))
  
  significant<-(abs(X_meta/sigma_meta)>critval);
  #%added to exclude psych outlier
  nooutlier<-sigma<50;
  dat<-data.frame(xvar=X_meta,
                  yvar=sigma_meta,
                  significant = significant &nooutlier )
  rangeX=1.1*max(max(abs(X_meta)), max(abs(sigma_meta[nooutlier]))*critval);
  p2<-ggplot(dat, aes(x=xvar,y=yvar)) +
    xlim(c(-rangeX,rangeX))+
    ylim(c(min(sigma)-0.03,rangeX/critval))+
    xlab("|X|")+
    ylab(expression("sigma" ))+
    geom_abline(intercept = 0,slope=1/critval,color="grey")+ 
    geom_abline(intercept = 0,slope=-1/critval,color="grey")+
    geom_point(shape=21,size = 2,aes(fill = significant))
  
  p12<-grid.arrange(p1, p2, ncol = 2)
  filepath<-paste0(pathname,'/FiguresandTables/',name, 'CombinedScatter.pdf')
  save_plot(filepath,p12,ncol=2,base_width = 4, base_height=3)
  dev.off()
  #####
  Zuse<-X_meta/sigma_meta;
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
  
  p123<-grid.arrange(h,p1, p2, ncol = 3)
  filepath<-paste0(pathname,'/FiguresandTables/',name, 'CombinedScatterHist.pdf')
  save_plot(filepath,p123,ncol=3,base_width = 4, base_height=3)
  dev.off()
}