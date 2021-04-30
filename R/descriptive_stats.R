#' Title
#'
#' @param X
#' @param sigma
#' @param identificationapproach
#' @param name
#' @param symmetric
#' @param cluster_ID
#'
#' @return
#' @export
#'
#' @examples
descriptive_stats <- function(X, sigma, identificationapproach, name, symmetric, cluster_ID) {

    # for approach 1:
    #   X has 2 columns, initial estimate Zstat, and replication estimate
    #   normalized by sigma1; sigma is SE2/SE1
    # for approach 2:
    #   X has 1 column, initial estimate unscaled, sigma=sigma1
    n <- nrow(X)
    # Significance
    critval <- 1.96

    # Define blue colors
    color <- brewer.pal(8,"Blues")
    blue1 <- color[8]
    blue2 <- color[5]
    blue3 <- color[3]

    if (identificationapproach==1) {
        ll<-floor(min(min(X)));
        uu<-ceiling(max(max(X)));
        significant<-(abs(X[,1])>critval);
        dat<-data.frame(xvar = X[,1],
                        yvar = X[,2],
                        significant = significant)

        # Plot scatter
        p <- ggplot(dat, aes(x=xvar, y=yvar)) +geom_vline(xintercept =critval,color="grey")+
            geom_hline(yintercept =critval,color="grey")+
            geom_abline(intercept = 0,slope=1,color="grey")+
            xlab("Z")+
            ylab("ZRep" )+
            xlim(c(max(ll,0),uu))+
            ylim(c(ll,uu))+
            geom_point(shape=21,size = 2,aes(fill = significant)) +
            scale_fill_manual(values=c("#FF0000", blue2))+
            theme_minimal()+ theme(axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.grid.major.x = element_blank())

        ggsave(paste0(getwd(),"/Scatter.pdf"), width = 10, height = 6.5)

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
            geom_point(shape=21,size = 2,aes(fill = significant))+
            scale_fill_manual(values=c("#FF0000", blue2))+
            theme_minimal()+ theme(axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.grid.major.x = element_blank())

        ggsave(paste0(getwd(),"/Scatter.pdf"), width = 10, height = 6.5)

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
                           fill = blue1,
                           breaks=edges)+
            geom_vline(xintercept =-1.96,color='grey')+
            geom_vline(xintercept =1.96, color='grey')+
            xlab('X/sigma')+
            ylab('Density')+
            xlim(c(min(edges),max(edges)))+
            theme_minimal()+ theme(axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.grid.major.x = element_blank())


    } else {
        h<-ggplot(data = as.data.frame(Zuse), aes(Zuse))+
            geom_histogram(aes(y = ..density..),
                           fill = blue1,
                           breaks=edges)+
            geom_vline(xintercept =1.96, color='grey')+
            xlab('|X|/sigma')+
            ylab('Density')+
            xlim(c(min(edges),max(edges)))+
            theme_minimal()+ theme(axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.grid.major.x = element_blank())
    }

    handp<-grid.arrange(h, p, ncol = 2)
    ggsave(plot = handp, paste0(getwd(),"/Scatter_Hist.pdf"), width = 10, height = 6.5)

    h <- ggplotly(h)
    p <- ggplotly(p)

    return(list("hist"= h,"scatter"= p))

}

