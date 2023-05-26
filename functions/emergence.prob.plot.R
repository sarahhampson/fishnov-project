# Function to plot the results of GLM Intercepts and Confidence Intervals
# Requires a 4x4 dataframe with category, intercept, lower CI and upper CI

emergence.prob.plot <- function(prob.df, plotname, nov.type){
  
  # Verification variable
  #prob.df <- pilot.tax.GLM.df
  #plotname <- "trialtaxplot"
  #nov.type= "Tax"

  # Stop any plots
  graphics.off()
  
  # Change location of plot/name depending on novelty type
  if(nov.type == "Tax"){
    plot_filename <- paste0("./plots/Pilot Taxonomic Plots/", plotname, "_plot")
  }
  else if(nov.type == "Func"){
    plot_filename <- paste0("./plots/Pilot Functional Plots/", plotname, "_plot")
  }
  else{
    plot_filename <- paste0("./plots/", nov.type, "_", plotname, "_plot")
  }
  
  # Start writing PDF
  pdf(file=date.wrap(plot_filename, ".pdf"), width = 7, height = 6)
  
  # Cumulative colours
  cumul.col <- rgb(0.373,0.651,0.765)
  
  ## Plot Intercepts
  par(mar=c(4, 5, 0.1, 8), oma=c(1, 1, 3, 1), xpd=TRUE)
  plot(prob.df$Intercept[2:4] ~ rep(1:3, 1), 
       axes = F, 
       ylab = "Emergence Probability", 
       xlab = "Novelty Type", 
       xlim = c(0.5, 3.5), 
       ylim= c(0, max(prob.df$`CI Upper`[2:4]+0.01)), 
       type= "n"
  )
  
  # Draw lines for confidence intervals
  arrows(x0=1, y0=prob.df$`CI Lower`[2], 
         x1=1, y1=prob.df$`CI Upper`[2], 
         code=3, angle=90, length=0.05, col="black", lwd=1)
  arrows(x0=2, y0=prob.df$`CI Lower`[3], 
         x1=2, y1=prob.df$`CI Upper`[3], 
         code=3, angle=90, length=0.05, col="black", lwd=1)
  arrows(x0=3, y0=prob.df$`CI Lower`[4], 
         x1=3, y1=prob.df$`CI Upper`[4], 
         code=3, angle=90, length=0.05, col="black", lwd=1)
  
  # Include points with novelty colours
  points(prob.df$Intercept[2:4],
         col="black",
         pch=21,
         cex=1,
         bg= c(cumul.col, "red", "orange")
  )
  
  # Include axes
  axis(side = 1, # what axis are we plotting (same as mar, 1=bottom, 2=left, 3=top, 4=right)
       at = c(1:3), # where do we want to put tick marks?
       labels = c("C", "I", "N"), # what labels do we want to put beside those tick marks?
       tck = -0.02, #How big should the tick marks be?
       padj = -1.5, #How far away should the labels be from the tick marks
       las = 1, cex.axis = 0.8)
  axis(side = 2, 
       tck = -0.02, 
       hadj = 0.8, 
       las = 1, cex = 0.5,
       cex.axis=0.8)
  
  # Add legend
  legend("right", inset = c(-0.4, 0),
         legend=c("Cumulative novelty", "Instantaneous novelty", "Novel community"),
         col="black",
         pt.bg=c(cumul.col, "red", "orange"),
         pch=21,
         cex=0.75,
         bty="n"
  )
  
  # Box around plot
  box(which = "plot", lty = "solid")
  
  dev.off()

}

