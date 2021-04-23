#' Title
#'
#' @param X
#' @param sigma
#' @param Psihat
#' @param Varhat
#' @param cutoffs
#' @param symmetric
#' @param symmetric_p
#' @param Studynames
#' @param identificationapproach
#'
#' @return
#' @export
#'
#' @examples
plot_correction <- function(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,Studynames,identificationapproach) {

  color <- brewer.pal(8,"Greys")
  GREY <- color[3]
  SLIDE_GREY <- color[5]
  color <- brewer.pal(8,"Purples")
  BLUE <- color[6]

  # Plot original and adjusted confidence sets ------------------------------

  Rl=min(min(Z1_L),min(Z1)-2)-0.5
  Ru=max(max(Z1_U),max(Z1)+2)+0.5
  R=Ru-Rl

  W <- 6
  H <- 3.5
  ytick<-Studynames
  xtick<-seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  pdf(here("FiguresandTables","OriginalAndAdjusted.pdf"), width = 6, height = 3.5)
  par(mgp = c(3,0.5,0), mar = c(1, 13, 2, 1))
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(0,H), col = BLUE,cex = 0.3); #original
  par(new=TRUE)
  lines(c(0,0),c(0,n+1),col=GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96/R*W,1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  lines(c(-1.96/R*W,-1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(0,H), col = BLUE,cex = 0.3); #original
  par(new=TRUE)
  plot(x = Z1_M/R*W, y = (c(1:n)-0.3)/n*H,yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(0,H),cex = 0.3); #adjusted
  par(new=TRUE)
  for (i in c(1:n)) {
    lines(c(Z1[i]-1.96,Z1[i]+1.96)/R*W,c(i,i)/n*H, xlab="", ylab="",col=BLUE,lwd=1.5);
    lines(c(Z1_L[i],Z1_U[i])/R*W,(c(i,i)-0.3)/n*H, xlab="", ylab="",lwd=1.5);
  }

  axis(2, at = seq(1/n,H,H/n), labels = ytick, las = 1, cex.axis=0.75)
  axis(3, at=seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels = xtick, cex.axis=0.6)

  legend(2,1,legend=c("Original Estimates", "Adjusted Estimates"),col=c(BLUE,"black"), pch=c(4,19), cex=0.5)

  dev.off()


# Plot Original and adjusted confidence sets, including bonferroni  -------

  Rl=min(min(Z1_LB),min(Z1)-2)-0.5
  Ru=max(max(Z1_UB),max(Z1)+2)+0.5
  R=Ru-Rl

  W <- 6
  H <- 3.5
  ytick<-Studynames
  xtick<-seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  pdf(here("FiguresandTables","OriginalAndAdjustedBonferroni.pdf"), width = 6, height = 3.5)
  par(mgp = c(3,0.5,0), mar = c(1, 13, 2, 1))
  par(new=TRUE)
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(1,n)*H/n, col = BLUE,cex = 0.3); #original
  lines(c(0,0),c(0,n+1),col=GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96/R*W,1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  lines(c(-1.96/R*W,-1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(1,n)*H/n, col = BLUE,cex = 0.3,add=TRUE); #original
  par(new=TRUE)
  plot(x = Z1_M/R*W, y = (c(1:n)-0.35)/n*H,yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(1,n)*H/n ,cex = 0.3); #adjusted
  par(new=TRUE)
  for (i in c(1:n)) {
    lines(c(Z1[i]-1.96,Z1[i]+1.96)/R*W,c(i,i)/n*H, xlab="", ylab="",col=BLUE,lwd=1.5);
    lines(c(Z1_L[i],Z1_U[i])/R*W,(c(i,i)-0.3)/n*H, xlab="", ylab="",lwd=1.5);

    lines(c(1,1)*Z1_LB[i]/R*W,c(i-0.4,i-0.2)/n*H, xlab="", ylab="");
    lines(c(1,1)*Z1_UB[i]/R*W,(c(i-0.4,i-0.2))/n*H, xlab="", ylab="");
  }

  axis(2, at = seq(1/n,H,H/n), labels = ytick, las = 1,cex.axis=0.75) # TSA CHECK
  axis(3, at=seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels = xtick,cex.axis=0.6) # TSA CHECK

  legend(2,1,legend=c("Original Estimates", "Adjusted Estimates"),col=c(BLUE,"black"), pch=c(4,19), cex=0.5)

  dev.off()

# Plot original and adjusted confidence sets for slides ---------------------------
  Rl <- min(min(Z1_LB),min(Z1)-2)-0.5
  Ru <- max(max(Z1_UB),max(Z1)+2)+0.5
  R <- Ru-Rl

  ytick<-Studynames
  xtick<-seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  pdf(here("FiguresandTables","OriginalAndAdjustedSlides.pdf"), width = 6, height = 3.5)
  par(mgp = c(3,0.5,0), mar = c(1, 13, 2, 1))
  par(new=TRUE)
  plot(x = Z1, y = (1:n),
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl, Ru), ylim = c(1,n+1), col = BLUE,cex = 0.3); #original
  lines(c(0,0),c(0,n+1),col=SLIDE_GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96, 1.96), c(1,n+1), xlab="", ylab="",col=SLIDE_GREY);
  par(new=TRUE)
  lines(c(-1.96,-1.96),c(1,n+1), xlab="", ylab="",col=SLIDE_GREY);
  par(new=TRUE)
  plot(x = Z1, y = (1:n),
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl,Ru), ylim = c(1,n+1), col = BLUE,cex = 0.3); #original
  par(new=TRUE)
  plot(x = Z1_M, y = (1:n)-0.35,yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
       xlim = c(Rl,Ru), ylim = c(1,n+1) ,cex = 0.3); #adjusted
  par(new=TRUE)
  for (i in c(1,n+1)) {
    lines(c(Z1[i]-1.96,Z1[i]+1.96),c(i,i), xlab="", ylab="",col=BLUE,lwd=1.5) # original
    lines(c(Z1_L[i],Z1_U[i]),(c(i,i)-0.2), xlab="", ylab="",lwd=1.5) # adjusted
  }

  axis(2, at = seq(1/n,H,H/n), labels = ytick, las = 1,cex.axis=0.75) # TSA CHECK
  axis(3, at=seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels = xtick,cex.axis=0.6) # TSA CHECK

  dev.off()

# Plot original, adjusted, and replication point estimates ---------------------------
  if (identificationapproach == 1) {
  Z2 <- Z[ ,2]
  Z2_rescaled <- Z2

  H <- 3.5
  W <- 6

  ytick<-Studynames
  xtick<-seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  pdf(here("FiguresandTables","OriginalReplicationAndAdjusted.pdf"), width = 6, height = 3.5)
  par(mgp = c(3,0.5,0), mar = c(1, 13, 2, 1))
  par(new=TRUE)
  plot(x = Z1/R*W, y = seq(1/n,H,H/n),
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = H, col = BLUE,cex = 0.3); #original
  lines(c(0,0),c(0,n+1),col=GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96/R*W,1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  lines(c(-1.96/R*W,-1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  plot(x = Z1/R*W, y = seq(1/n,H,H/n),
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = H, col = BLUE,cex = 0.3, add = TRUE); #original
  par(new=TRUE)
  plot(x = Z1_M/R*W, y = seq(1/n,H,H/n),yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = H, col = GREY ,cex = 0.3, add = TRUE); #adjusted
  par(new=TRUE)
  plot(x = Z2_rescaled/R*W, y = seq(1/n,H,H/n),yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = H ,cex = 0.3, add = TRUE); #adjusted

  axis(2, at = seq(1/n,H,H/n), labels = ytick, las = 1,cex.axis=0.75)
  axis(3, at = seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels = xtick,cex.axis=0.6)

  dev.off()

}

# Plot original, adjusted, and replication point estimates for slides ---------------------------
if (identificationapproach == 1) {
  Z2 <- Z[ ,2]
  Z2_rescaled <- Z2

  ytick<-Studynames
  xtick<-seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  pdf(here("FiguresandTables","OriginalReplicationAndAdjustedSlides.pdf"), width = 6, height = 3.5)
  par(mgp = c(3,0.5,0), mar = c(1, 13, 2, 1))
  par(new=TRUE)
  plot(x = Z1, y = 1:n,
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = H, col = BLUE,cex = 0.3); #original
  lines(c(0,0),c(0,n+1),col=GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96,1.96),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  lines(c(-1.96,-1.96),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  plot(x = Z1, y = 1:n,
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = H, col = BLUE,cex = 0.3, add = TRUE); #original
  par(new=TRUE)
  plot(x = Z1_M, y = 1:n,yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = H, col = GREY ,cex = 0.3, add = TRUE); #adjusted
  par(new=TRUE)
  plot(x = Z2_rescaled, y = 1:n,yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = H ,cex = 0.3, add = TRUE); #adjusted

  axis(2, at = 1:n, labels = ytick, las = 1,cex.axis=0.75)
  axis(3, at = seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels = xtick,cex.axis=0.6)

  dev.off()

}


# Plot corrected estimators and confidence sets - Correction Plot---------------------------

  xmin=min(xgrid)
  xmax=max(xgrid)

  ymin=min(min(Theta_LB_store,xmin-1.96))
  ymax=max(max(Theta_UB_store,xmax+1.96))

  pdf(here("FiguresandTables","GMMresultsCorrection_plot.pdf"), width = 4, height = 4)
  par(mar = c(3, 3, 1, 1))
  plot(x = xgrid, y = Theta_L_store,
       lty = 1, xlab='', ylab='',
       xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       col = BLUE,cex=0.1)
  title(ylab=TeX('Estimation $\\theta$'), xlab=TeX('$X$'), line=1.75, cex.lab=0.8)
  par(new=TRUE)
  for (n in seq(ceiling(xmin)+1,floor(xmax)-1,1)) {
    lines(c(n,n),c(ymin,ymax),lty=3,col=GREY,lwd=1.5)
  }
  par(new=TRUE)
  for (n in seq(ceiling(ymin),floor(ymax),1)) {
    lines(c(xmin,xmax),c(n,n),lty=3,col=GREY,lwd=1.5)
  }
  par(new=TRUE)
  lines(xgrid,xgrid-1.96,col=GREY, xlab="", ylab="",lwd=1.5);
  par(new=TRUE)
  lines(xgrid,xgrid, xlab="", ylab="",col=GREY,lwd=1.5);
  par(new=TRUE)
  lines(xgrid,xgrid+1.96, xlab="", ylab="",col=GREY,lwd=1.5);
  par(new=TRUE)
  lines(x = xgrid, y = Theta_L_store,
        lty = 1,   xlab=TeX('$X$'), ylab=TeX('Estimation $\\theta$'),
        xlim=c(xmin,xmax), ylim=c(ymin,ymax),
        col = BLUE,lwd=1.5);
  par(new=TRUE)
  lines(x = xgrid, y = Theta_M_store,
        lty = 1,   xlab="", ylab="",
        col = BLUE,lwd=3);
  par(new=TRUE)
  lines(x = xgrid, y = Theta_LB_store,
        lty = 1,   xlab="", ylab="",
        col = "black");
  par(new=TRUE)
  lines(x = xgrid, y = Theta_UB_store,
        lty = 1,   xlab="", ylab="",
        col = "black");
  lines(x = xgrid, y = Theta_U_store,
        lty = 1,   xlab="", ylab="",
        col = BLUE,lwd=1.5);

  legend(0,7,legend=c("95% confidence Bounds", "Bonferroni Corrected 95% Bounds", "Median Unbiased Estimator"),col=c(BLUE,"black",BLUE), lty=1, lwd=c(1.5,1,3), cex=0.6)

  dev.off()

}
