

GA_Bfragilis_Scatter <- function(DfPlot, outfile = "Fig2_ScatterPlot.pdf", width_use = 1.75, height_use = 1.65, pseudo = -11.0, highlight=FALSE,plot_groups=FALSE) {
  library('ggplot2')
  library("extrafont")
  if (plot_groups) {
    plot <- ggplot(DfPlot, aes(x = Bfr, y = EI)) + geom_abline(intercept = 0) + geom_point(aes(fill = OrphanStatus), color = "black", pch = 21)
    
  } else {
    plot <- ggplot(DfPlot, aes(x = Bfr, y = EI)) + geom_abline(intercept = 0) + geom_point(data = subset(DfPlot, Bfr > pseudo), color = "black", pch = 21, fill = "#d3d3d3") 
    if (highlight) {
      #plot <- plot + geom_point(data = subset(DfPlot, (Bfr == pseudo) & (Highlight == "True")), color = "black", pch = 21, fill = "red")
      plot <- plot + geom_point(data = subset(DfPlot, (Bfr == pseudo) ), color = "black", pch = 21, fill = "red")
      print("Total number of points")
      print(dim(DfPlot))
      print("Highlight number of points")
      print(dim(subset(DfPlot, (Bfr == pseudo))))
    } else {
      plot <- plot + geom_point(data = subset(DfPlot, Bfr == pseudo), color = "black", pch = 21, fill = "#d3d3d3")
    } }
  plot <- plot + theme_classic() + scale_x_continuous(limits = c(pseudo,-7)) + scale_y_continuous(limits = c(pseudo,-7)) + geom_vline(xintercept=pseudo + 0.5, linetype = 2) + 
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      panel.grid.major.y = element_line(colour = "grey", size = 0.125),
      panel.grid.major.x = element_line(colour = "grey", size = 0.125),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25),
      panel.border = element_rect(color = "black", size = 0.5, fill = NA)
    ) + theme(legend.position = "none")
  
  if (!is.null(outfile)) {
    pdf(outfile, width = width_use, height = height_use,  useDingbats=FALSE)
    #pdf(outfile, width = 3, height = 3, family = "Arial")
    print(plot)
    dev.off()
  } else {
    print(plot)
  }
}

