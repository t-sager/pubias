#' Title
#'
#' @param X
#' @param sigma
#' @param Psihat
#' @param Varhat
#' @param cutoffs
#' @param symmetric
#' @param symmetric_p
#' @param studynames
#' @param identificationapproach
#'
#' @return
#' @export
#'
#' @examples
plot_correction <- function(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,studynames,identificationapproach, corrected_estimates) {

  Z1 <- corrected_estimates$Z1
  Z1_U <- corrected_estimates$Z1_U
  Z1_L <- corrected_estimates$Z1_L
  Z1_M <- corrected_estimates$Z1_M
  Z1_UB <- corrected_estimates$Z1_UB
  Z1_LB <- corrected_estimates$Z1_LB
  Theta_U_store <- corrected_estimates$Theta_U_store
  Theta_L_store <- corrected_estimates$Theta_L_store
  Theta_M_store <- corrected_estimates$Theta_M_store
  Theta_UB_store <- corrected_estimates$Theta_UB_store
  Theta_LB_store <- corrected_estimates$Theta_LB_store
  xgrid <- corrected_estimates$xgrid

  n <- nrow(Z1)

  # Define blue colors
  color <- brewer.pal(8,"Blues")
  blue1 <- color[8]
  blue2 <- color[5]
  blue3 <- color[3]

  # GGPLOT: Plot original and adjusted confidence sets ------------------------------
  Rl=min(min(Z1_L),min(Z1)-2)-0.5
  Ru=max(max(Z1_U),max(Z1)+2)+0.5
  R=Ru-Rl
  W <- 7
  H <- 6.5
  ylabel <- studynames
  xlabel <- seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  segment_data_og <- data.frame(
    x = (Z1-1.96)/R*W,
    xend = (Z1+1.96)/R*W,
    y = c(1:n)/n*H+0.1,
    yend = c(1:n)/n*H+0.1
  )

  segment_data_adj <- data.frame(
    x = Z1_L/R*W,
    xend = Z1_U/R*W,
    y = c(1:n)/n*H,
    yend = c(1:n)/n*H
  )

cols <- c("Adjusted Estimates"=color[5],"Original Estimates"=color[8])

  g <- ggplot()+
    geom_point(aes(x = Z1/R*W, y = c(1:n)/n*H+0.1, colour = "Original Estimates")) + # original
    geom_point(aes(x = Z1_M/R*W, y = c(1:n)/n*H, colour = "Adjusted Estimates")) + # adjusted
    geom_segment(data = segment_data_og, aes(x = x, y = y, xend = xend, yend = yend), color = blue2) + # original
    geom_segment(data = segment_data_adj, aes(x = x, y = y, xend = xend, yend = yend), color = blue1) + # adjusted
    geom_line(aes(x = c(1.96/R*W,1.96/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    geom_line(aes(x = c(-1.96/R*W,-1.96/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    scale_y_continuous(breaks = c(1:n)/n*H+0.05, labels=ylabel) +
    scale_x_continuous(breaks = seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels=xlabel) +
    scale_colour_manual(values=cols)+
    xlab("")+
        theme_minimal()+ theme(axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.grid.major.x = element_blank(),
                               legend.position="top",
                               legend.margin=margin(0,0,0,0),
                               legend.box.margin=margin(5,5,-17,-17),
                               legend.title = element_blank())
  g

  ggsave("OriginalAndAdjusted.pdf", width = 10, height = 6.5)


  # GGPLOT: Original and adjusted confidence sets, including bonferroni ------------------------------
  Rl=min(min(Z1_LB),min(Z1)-2)-0.5
  Ru=max(max(Z1_UB),max(Z1)+2)+0.5
  R=Ru-Rl
  W <- 7
  H <- 6.5
  ylabel <- studynames
  xlabel <- seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  segment_data_og <- data.frame(
    x = (Z1-1.96)/R*W,
    xend = (Z1+1.96)/R*W,
    y = c(1:n)/n*H+0.1,
    yend = c(1:n)/n*H+0.1
  )

  segment_data_adj <- data.frame(
    x = Z1_LB/R*W,
    xend = Z1_UB/R*W,
    y = c(1:n)/n*H,
    yend = c(1:n)/n*H
  )

  cols <- c("Adjusted Estimates"=color[5],"Original Estimates"=color[8])

  g <- ggplot()+
    geom_point(aes(x = Z1/R*W, y = c(1:n)/n*H+0.1, colour = "Original Estimates")) + # original
    geom_point(aes(x = Z1_M/R*W, y = c(1:n)/n*H, colour = "Adjusted Estimates")) + # adjusted
    geom_segment(data = segment_data_og, aes(x = x, y = y, xend = xend, yend = yend), color = blue2) + # original
    geom_segment(data = segment_data_adj, aes(x = x, y = y, xend = xend, yend = yend), color = blue1) + # adjusted
    geom_line(aes(x = c(1.96/R*W,1.96/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    geom_line(aes(x = c(-1.96/R*W,-1.96/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    scale_y_continuous(breaks = c(1:n)/n*H+0.05, labels=ylabel) +
    scale_x_continuous(breaks = seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels=xlabel) +
    scale_colour_manual(values=cols)+
    xlab("")+
    theme_minimal()+ theme(axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.grid.major.x = element_blank(),
                           legend.position="top",
                           legend.margin=margin(0,0,0,0),
                           legend.box.margin=margin(5,5,-17,-17),
                           legend.title = element_blank())
  g

  ggsave("OriginalAndAdjustedBonferroni.pdf", width = 10, height = 6.5)

  # GGPLOT: Original and adjusted confidence sets, including bonferroni ------------------------------
  if (identificationapproach == 1) {
    Z2 <- X[ ,2]
    Z2_rescaled <- Z2

  # Uses Bonferroni

  W <- 7
  H <- 6.5
  ylabel <- studynames
  xlabel <- seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  segment_data_og <- data.frame(
    x = (Z1-1.96)/R*W,
    xend = (Z1+1.96)/R*W,
    y = c(1:n)/n*H+0.1,
    yend = c(1:n)/n*H+0.1
  )

  segment_data_adj <- data.frame(
    x = Z1_LB/R*W,
    xend = Z1_UB/R*W,
    y = c(1:n)/n*H,
    yend = c(1:n)/n*H
  )

  cols <- c("Adjusted Estimates"=color[5],"Original Estimates"=color[8], "Replication Estimates"=color[3])

  g <- ggplot()+
    geom_point(aes(x = Z1/R*W, y = c(1:n)/n*H+0.1, colour = "Original Estimates")) + # original
    geom_point(aes(x = Z1_M/R*W, y = c(1:n)/n*H, colour = "Adjusted Estimates")) + # adjusted
    geom_point(aes(x = Z2_rescaled/R*W, y = c(1:n)/n*H+0.2, colour = "Replication Estimates")) + # replication
    geom_segment(data = segment_data_og, aes(x = x, y = y, xend = xend, yend = yend), color = blue2) + # original
    geom_segment(data = segment_data_adj, aes(x = x, y = y, xend = xend, yend = yend), color = blue1) + # adjusted
    geom_line(aes(x = c(1.96/R*W,1.96/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    geom_line(aes(x = c(-1.96/R*W,-1.96/R*W),y = c(0,max(segment_data_og$yend)+0.1)), color = "grey") +
    scale_y_continuous(breaks = c(1:n)/n*H+0.1, labels=ylabel) +
    scale_x_continuous(breaks = seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels=xlabel) +
    scale_colour_manual(values=cols)+
    xlab("")+
    theme_minimal()+ theme(axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.grid.major.x = element_blank(),
                           legend.position="top",
                           legend.margin=margin(0,0,0,0),
                           legend.box.margin=margin(5,5,-17,-17),
                           legend.title = element_blank())
  g

  ggsave("OriginalReplicationAndAdjusted.pdf", width = 10, height = 6.5)

  }

  # GGPLOT: Correction Plot  ------------------------------

  xmin=min(xgrid)
  xmax=max(xgrid)

  ymin=min(min(Theta_LB_store,xmin-1.96))
  ymax=max(max(Theta_UB_store,xmax+1.96))

  cols <- c("95% confidence Bounds"=color[5],"Bonferroni Corrected 95% Bounds"=color[8], "Median Unbiased Estimator"=color[3])

  g <- ggplot()+
    geom_line(aes(x = xgrid, y = Theta_L_store,colour = "95% confidence Bounds"), lwd = 1) +
    geom_line(aes(x = xgrid, y = Theta_M_store,colour = "Median Unbiased Estimator"), lwd = 1) +
    geom_line(aes(x = xgrid, y = Theta_LB_store,colour = "Bonferroni Corrected 95% Bounds"), lwd = 1) +
    geom_line(aes(x = xgrid, y = Theta_UB_store,colour = "Bonferroni Corrected 95% Bounds"), lwd = 1) +
    geom_line(aes(x = xgrid, y = Theta_U_store,colour = "95% confidence Bounds"), lwd = 1) +
    geom_line(aes(x = xgrid, y = xgrid - 1.96), color = "grey") +
    geom_line(aes(x = xgrid, y = xgrid), color = "grey") +
    geom_line(aes(x = xgrid, y = xgrid + 1.96), color = "grey") +
    xlab(TeX('$X$')) +
    ylab(TeX('Estimation $\\theta$')) +
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

  ggsave("CorrectionPlot.pdf", width = 10, height = 6.5)

}
