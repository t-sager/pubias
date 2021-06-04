descriptive_stats <- function(X, sigma, identificationapproach, name, symmetric, cluster_ID, cutoffs) {

    n <- nrow(X)

    # Cutoff
    critval <- cutoffs

    # Define blue colors
    color <- brewer.pal(8,"Blues")
    blue1 <- color[8]
    blue2 <- color[5]
    blue3 <- color[3]

    if (identificationapproach == 1) {
        ll <- floor(min(min(X)))
        uu <- ceiling(max(max(X)))
        significant<-(abs(X[,1])>critval)
        dat <- data.frame(xvar = X[,1],
                        yvar = X[,2],
                        significant = significant)

        # Scatter plot if replication
        p <- ggplot(dat, aes(x=xvar, y=yvar)) +geom_vline(xintercept =critval,color="grey")+
            geom_hline(yintercept =critval,color="grey")+
            geom_abline(intercept = 0,slope=1,color="grey")+
            xlab("Z")+
            ylab("ZRep")+
            xlim(c(max(ll,0),uu))+
            ylim(c(ll,uu))+
            geom_point(shape=21,size = 2,aes(fill = significant)) +
            scale_fill_manual(values=c("#FF0000", blue2))+
            theme_minimal()+ theme(panel.grid.minor = element_blank(),
                                   panel.grid.major.x = element_blank())

        ggsave(paste0(getwd(),"/Scatter.pdf"), width = 10, height = 5)

    } else  if (identificationapproach == 2)  {

        # Scatter plot if meta
        significant <- (abs(X/sigma)>critval)
        nooutlier <- sigma<50
        dat <- data.frame(xvar=X,
                        yvar=sigma,
                        significant = significant &nooutlier )
        rangeX <- 1.1*max(max(abs(X)), max(abs(sigma[nooutlier]))*critval)
        p <- ggplot(dat, aes(x=xvar,y=yvar)) +
            xlim(c(-rangeX,rangeX))+
            ylim(c(min(sigma)-0.03,rangeX/critval))+
            xlab("X")+
            ylab("sigma" )+
            geom_abline(intercept = 0,slope=1/critval,color="grey")+
            geom_abline(intercept = 0,slope=-1/critval,color="grey")+
            geom_point(shape=21,size = 2,aes(fill = significant))+
            scale_fill_manual(values=c("#FF0000", blue2))+
            theme_minimal()+ theme(panel.grid.minor = element_blank(),
                                   panel.grid.major.x = element_blank())

        ggsave(paste0(getwd(),"/Scatter.pdf"), width = 10, height = 5)

    }

    # Create histogram plot for both identifiactionapproach = 1 and 2
    if (identificationapproach==1) {
        Zuse = Z[,1];

    } else {
        Zuse <- X/sigma
        if (max(abs(Zuse))>10) {
            Zuse <- Zuse[abs(Zuse)<6]
        }
    }

    # Axis scaling
    ll <- floor(min(Zuse))
    lleven <- floor(ll/2)*2
    uu <- ceiling(max(Zuse))

    if (symmetric == FALSE) {
        if (n >= 30) {
            uu2 <- ceiling((uu-.36)/.32)*.32+.36
            ll2 <- floor((ll+.36)/.32)*.32-.36
            edges<-c(seq(from = ll2,
                         to = -0.36,
                         by = ifelse(ll2>0,-0.32-0.32)), 0, seq(from = 0.36,
                                          to = uu2,
                                          by = 0.32))
        } else {
            uu2 <- ceiling((uu-.68)/.64)*.64+.68
            ll2 <- floor((ll+.68)/.64)*.64-.68
            edges<-c(seq(from = ll2,
                         to = -0.68,
                         by = ifelse(ll2>0,-0.64, 0.64)), 0, seq(from = 0.68,
                                          to = uu2,
                                          by = 0.64))
        }
    }

    edges <- seq(from = ll,
                    to = uu,
                    by = 0.64)



    if (symmetric == FALSE) {

        h <- ggplot(data = as.data.frame(Zuse), aes(Zuse))+
            geom_histogram(aes(y = ..density..),
                           fill = blue1,
                           binwidth = 0.3)+
            geom_vline(xintercept =-cutoffs,color='grey')+
            geom_vline(xintercept =cutoffs, color='grey')+
            geom_vline(xintercept =0, color='grey', linetype="dotted")+
            xlab('Z')+
            ylab('Density')+
            xlim(c(min(edges),max(edges)))+
            theme_minimal()+ theme(panel.grid.minor = element_blank(),
                                   panel.grid.major.x = element_blank())


    } else {
        h <- ggplot(data = as.data.frame(Zuse), aes(Zuse))+
            geom_histogram(aes(y = ..count..),
                           fill = blue1,
                           binwidth = 0.3)+
            geom_vline(xintercept =-cutoffs, color='grey')+
            geom_vline(xintercept =cutoffs, color='grey')+
            geom_vline(xintercept =0, color='grey', linetype="dotted")+
            xlab('Z')+
            ylab('Density')+
            xlim(c(min(edges),max(edges)))+
            theme_minimal()+ theme(panel.grid.minor = element_blank(),
                                   panel.grid.major.x = element_blank())
    }

    comb <- grid.arrange(h, p, ncol = 2)
    ggsave(plot = comb, paste0(getwd(),"/Scatter_Hist.pdf"), width = 10, height = 5)

    h <- ggplotly(h)
    p <- ggplotly(p)

    return(list("hist"= h,"scatter"= p))

}

