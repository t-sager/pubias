plot_correction <- function(X, corrected_estimates,cutoffs,symmetric,studynames,identificationapproach) {

  # Get results from correction
  original <- corrected_estimates$original
  adj_U <- corrected_estimates$adj_U
  adj_L <- corrected_estimates$adj_L
  adj_estimates <- corrected_estimates$adj_estimates
  adj_UB <- corrected_estimates$adj_UB
  adj_LB <- corrected_estimates$adj_LB
  Theta_U_store <- corrected_estimates$Theta_U_store
  Theta_L_store <- corrected_estimates$Theta_L_store
  Theta_M_store <- corrected_estimates$Theta_M_store
  Theta_UB_store <- corrected_estimates$Theta_UB_store
  Theta_LB_store <- corrected_estimates$Theta_LB_store
  xgrid <- corrected_estimates$xgrid

  # n
  n <- nrow(original)

  # Define blue colors
  color <- brewer.pal(8,"Blues")
  blue1 <- color[8]
  blue2 <- color[5]
  blue3 <- color[3]

  # GGPLOT: Plot original and adjusted confidence sets ------------------------------
  # Preliminaries
  Rl<-min(min(adj_L),min(original)-2)-0.5
  Ru<-max(max(adj_U),max(original)+2)+0.5
  R<-Ru-Rl
  W <- 7
  H <- 6.5
  ylabel <- studynames
  xlabel <- seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)
  cols <- c("Adjusted Estimates"=color[5],"Original Estimates"=color[8])

  # CIs
  segment_data_og <- data.frame(
    x = (original-1.96)/R*W,
    xend = (original+1.96)/R*W,
    y = c(1:n)/n*H+0.1,
    yend = c(1:n)/n*H+0.1
  )

  segment_data_adj <- data.frame(
    x = adj_L/R*W,
    xend = adj_U/R*W,
    y = c(1:n)/n*H,
    yend = c(1:n)/n*H
  )

  # Plot
  g <- ggplot()+
    geom_point(aes(x = original/R*W, y = c(1:n)/n*H+0.1, colour = "Original Estimates")) + # original
    geom_point(aes(x = adj_estimates/R*W, y = c(1:n)/n*H, colour = "Adjusted Estimates")) + # adjusted
    geom_segment(data = segment_data_og, aes(x = x, y = y, xend = xend, yend = yend), color = blue1) + # original
    geom_segment(data = segment_data_adj, aes(x = x, y = y, xend = xend, yend = yend), color = blue2) + # adjusted
    geom_line(aes(x = c(cutoffs/R*W,cutoffs/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    geom_line(aes(x = c(-cutoffs/R*W,-cutoffs/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    scale_y_continuous(breaks = c(1:n)/n*H+0.05, labels=ylabel) +
    scale_x_continuous(breaks = seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels=xlabel) +
    scale_colour_manual(values=cols)+
    xlab("Z")+
        theme_minimal()+ theme(axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.grid.major.x = element_blank(),
                               legend.position="top",
                               legend.margin=margin(0,0,0,0),
                               legend.box.margin=margin(5,5,-17,-17),
                               legend.title = element_blank())
  g

  ggsave(paste0(getwd(),"/OriginalAndAdjusted.pdf"), width = 10, height = 5)

  OriginalAndAdjusted <- ggplotly(g)

  # GGPLOT: Original and adjusted confidence sets, including bonferroni ------------------------------
  # Preliminaries
  Rl<-min(min(adj_LB),min(original)-2)-0.5
  Ru<max(max(adj_UB),max(original)+2)+0.5
  R<-Ru-Rl
  W <- 7
  H <- 6.5
  ylabel <- studynames
  xlabel <- seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)
  cols <- c("Adjusted Estimates"=color[5],"Original Estimates"=color[8])

  # CIs
  segment_data_og <- data.frame(
    x = (original-1.96)/R*W,
    xend = (original+1.96)/R*W,
    y = c(1:n)/n*H+0.1,
    yend = c(1:n)/n*H+0.1
  )

  segment_data_adj <- data.frame(
    x = adj_LB/R*W,
    xend = adj_UB/R*W,
    y = c(1:n)/n*H,
    yend = c(1:n)/n*H
  )

  # Plot
  g <- ggplot()+
    geom_point(aes(x = original/R*W, y = c(1:n)/n*H+0.1, colour = "Original Estimates")) + # original
    geom_point(aes(x = adj_estimates/R*W, y = c(1:n)/n*H, colour = "Adjusted Estimates")) + # adjusted
    geom_segment(data = segment_data_og, aes(x = x, y = y, xend = xend, yend = yend), color = blue1) + # original
    geom_segment(data = segment_data_adj, aes(x = x, y = y, xend = xend, yend = yend), color = blue2) + # adjusted
    geom_line(aes(x = c(cutoffs/R*W,cutoffs/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    geom_line(aes(x = c(-cutoffs/R*W,-cutoffs/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    scale_y_continuous(breaks = c(1:n)/n*H+0.05, labels=ylabel) +
    scale_x_continuous(breaks = seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels=xlabel) +
    scale_colour_manual(values=cols)+
    xlab("Z")+
    theme_minimal()+ theme(axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.grid.major.x = element_blank(),
                           legend.position="top",
                           legend.margin=margin(0,0,0,0),
                           legend.box.margin=margin(5,5,-17,-17),
                           legend.title = element_blank())
  g

  ggsave(paste0(getwd(),"/OriginalAndAdjustedBonferroni.pdf"), width = 10, height = 5)

  OriginalAndAdjustedBonferroni <- ggplotly(g)

  # GGPLOT: Original, adjusted & replication confidence sets, including Bonferroni ------------------------------
  if (identificationapproach == 1) {

  # Preliminaries
  Z2 <- X[ ,2]
  Z2_rescaled <- Z2

  Rl<-min(min(adj_LB),min(Z2_rescaled)-2)-0.5
  Ru<-max(max(adj_UB),max(Z2_rescaled)+2)+0.5
  R<-Ru-Rl

  # Uses Bonferroni
  W <- 7
  H <- 6.5
  ylabel <- studynames
  xlabel <- seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)
  cols <- c("Adjusted Estimates"=color[5],"Original Estimates"=color[8], "Replication Estimates"=color[3])

  # CIs
  segment_data_og <- data.frame(
    x = (original-1.96)/R*W,
    xend = (original+1.96)/R*W,
    y = c(1:n)/n*H+0.1,
    yend = c(1:n)/n*H+0.1
  )

  segment_data_adj <- data.frame(
    x = adj_LB/R*W,
    xend = adj_UB/R*W,
    y = c(1:n)/n*H,
    yend = c(1:n)/n*H
  )

  segment_data_rep <- data.frame(
    x = (Z2_rescaled-1.96)/R*W,
    xend = (Z2_rescaled+1.96)/R*W,
    y = c(1:n)/n*H+0.2,
    yend = c(1:n)/n*H+0.2
  )

  # Plot
  g <- ggplot()+
    geom_point(aes(x = original/R*W, y = c(1:n)/n*H+0.1, colour = "Original Estimates")) + # original
    geom_point(aes(x = adj_estimates/R*W, y = c(1:n)/n*H, colour = "Adjusted Estimates")) + # adjusted
    geom_point(aes(x = Z2_rescaled/R*W, y = c(1:n)/n*H+0.2, colour = "Replication Estimates")) + # replication
    geom_segment(data = segment_data_og, aes(x = x, y = y, xend = xend, yend = yend, colour = "Original Estimates")) + # original
    geom_segment(data = segment_data_adj, aes(x = x, y = y, xend = xend, yend = yend, colour = "Adjusted Estimates")) + # adjusted
    geom_segment(data = segment_data_rep, aes(x = x, y = y, xend = xend, yend = yend, colour = "Replication Estimates")) + # replication
    geom_line(aes(x = c(cutoffs/R*W,cutoffs/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    geom_line(aes(x = 0,y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    scale_y_continuous(breaks = c(1:n)/n*H+0.1, labels=ylabel) +
    scale_x_continuous(breaks = seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels=xlabel) +
    scale_colour_manual(values=cols)+
    xlab("Z")+
    theme_minimal()+ theme(axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.grid.major.x = element_blank(),
                           legend.position="top",
                           legend.margin=margin(0,0,0,0),
                           legend.box.margin=margin(5,5,-17,-17),
                           legend.title = element_blank())
  g

  ggsave(paste0(getwd(),"/OriginalReplicationAndAdjusted.pdf"), width = 10, height = 5)

  OriginalReplicationAndAdjusted <- ggplotly(g)

  }

  # GGPLOT: Correction Plot  ------------------------------

  # Preliminaries
  xmin<-min(xgrid)
  xmax<-max(xgrid)
  ymin<-min(min(Theta_LB_store,xmin-cutoffs))
  ymax<-max(max(Theta_UB_store,xmax+cutoffs))
  cols <- c("95% confidence Bounds"=color[5],"Bonferroni Corrected 95% Bounds"=color[8], "Median Unbiased Estimator"=color[3])

  # Plot
  g <- ggplot()+
    geom_line(aes(x = xgrid, y = Theta_L_store,colour = "95% confidence Bounds"), lwd = 1) +
    geom_line(aes(x = xgrid, y = Theta_M_store,colour = "Median Unbiased Estimator"), lwd = 1) +
    geom_line(aes(x = xgrid, y = Theta_LB_store,colour = "Bonferroni Corrected 95% Bounds"), lwd = 1) +
    geom_line(aes(x = xgrid, y = Theta_UB_store,colour = "Bonferroni Corrected 95% Bounds"), lwd = 1) +
    geom_line(aes(x = xgrid, y = Theta_U_store,colour = "95% confidence Bounds"), lwd = 1) +
    geom_line(aes(x = xgrid, y = xgrid - cutoffs), color = "grey") +
    geom_line(aes(x = xgrid, y = xgrid), color = "grey") +
    geom_line(aes(x = xgrid, y = xgrid + cutoffs), color = "grey") +
    geom_vline(xintercept = -cutoffs, color='grey', linetype="dotted") +
    geom_vline(xintercept = cutoffs, color='grey', linetype="dotted") +
    xlab('Z') +
    ylab('Estimation') +
    scale_colour_manual(values=cols) +
    xlim(c(xmin, xmax)) +
    ylim(c(ymin, ymax)) +
    theme_minimal() + theme(
      panel.grid.minor.y = element_blank(),
      legend.position = "top",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(5, 5, -10, -10),
      legend.title = element_blank()
    )
  g

  ggsave(paste0(getwd(),"/CorrectionPlot.pdf"), width = 10, height = 5)

  CorrectionPlot <- ggplotly(g)

  if (identificationapproach == 1) {
  return(list("OriginalAndAdjusted"= OriginalAndAdjusted,"OriginalAndAdjustedBonferroni"= OriginalAndAdjustedBonferroni, "OriginalReplicationAndAdjusted" = OriginalReplicationAndAdjusted, "CorrectionPlot" = CorrectionPlot))
  } else {
    return(list("OriginalAndAdjusted"= OriginalAndAdjusted,"OriginalAndAdjustedBonferroni"= OriginalAndAdjustedBonferroni, "CorrectionPlot" = CorrectionPlot))
}

    }
