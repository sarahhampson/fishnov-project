# A function to create venn diagrams of GLM intercepts and CI intervals
# Requires a 4x4 dataframe with category, intercept, lower CI and upper CI

emergence.prob.venn.diagram <- function(prob.df, vdname, nov.type){
  
  # Verification variable
  #prob.df <- pilot.tax.GLM.df
  #vdname <- "trialtaxvenndiagram"
  
  # Require libraries
  require(plotrix)
  require(raster)
  require(sp)
  require(rgeos)
  
  # Start pdf
  graphics.off()
  
  # Swap the two rows of instantaneous and cumulative
  prob.df <- prob.df[c(1,3,2,4),]
  
  # Set conditions for where plots are to go
  if(nov.type == "Tax"){
    vd_filename <- paste0("./plots/Pilot Taxonomic Plots/", vdname, "_venndiagram")
  }
  else if(nov.type == "Func"){
    vd_filename <- paste0("./plots/Pilot Functional Plots/", vdname, "_venndiagram")
  }
  else{
    plot_filename <- paste0("./plots/", nov.type, "_", plotname, "_plot")
  }
  
  # Start writing PDF
  pdf(file=date.wrap(vd_filename, ".pdf"),  height=2.5, width=2.5)
  
  # Cumulative colours
  cumul.col <- rgb(0.373,0.651,0.765)
  
  # Set area for plot
  par(mar=c(0,0,0,0), ps=8, tcl = -0.2) #mgp=c(3,0.5,0)) <- originally for space underneath
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # Define circle variables
  circle.cent<-c(0.35,0.65,0.35)
  circle.radius <- c(0.3, 0.3)
  circle.y <- 0.5
  
  # Draw circles relative to the circle variables
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius[1], nv=500, 
                 border="black", col="white", lwd=2)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius[2], nv=500, 
                 border="black", col="white", lwd=2)
  
  # Convert circles to spatial polygons to get overlap between
  I.sp <- Polygon(cbind(I$x, I$y))
  I.sp <- SpatialPolygons(list(Polygons(list(I.sp), ID = "a")))
  C.sp <- Polygon(cbind(C$x, C$y))
  C.sp <- SpatialPolygons(list(Polygons(list(C.sp), ID = "a")))
  
  # N.sp is the spatial overlap of the two polygons 
  N.sp <- raster::intersect(I.sp, C.sp) # Just give the overlap
  I.sub <- gDifference(I.sp, C.sp) # The leftover I not covered by C
  C.sub <- gDifference(C.sp, I.sp) # The leftover C not covered by I
  
  # Change colours of polygons
  plot(I.sub, col="red", add=TRUE, border="black", lwd=1)
  plot(C.sub, col=cumul.col, add=TRUE, border="black", lwd=1)
  plot(N.sp, col="orange", add=TRUE, border="black", lwd=1.5)
  
  # Plot the labels on the outside 
  text(x=0.3, y=circle.y+0.4,
       labels=c("Instantaneous\nnovelty (I)    "), 
       col="red", adj=1, cex=0.9)
  text(x=0.75, y=circle.y+0.4,
       labels=c("Cumulative\n   novelty (C)"), 
       col=cumul.col, adj=0, cex=0.9)
  text(x=0.65, y=circle.y-0.4,
       labels=c("True novelty (N)"), 
       col="orange", adj=1, cex=0.9)

  # Create df of the percentage variables to put inside circles
  overall.means <- data.frame(matrix(ncol=3, nrow=3))
  overall.means[,1] <- prob.df[2:4, 2]
  overall.means[,2] <- prob.df[2:4, 3]
  overall.means[,3] <- prob.df[2:4, 4]
  colnames(overall.means) <- c("fit", "lower", "upper")
  overall.means <- overall.means*100
  
  # Wherever the centre is, put the means (offset to the left)
  text(x=c(circle.cent + c(-0.14, 0.15, 0.16)),
       y= rep(circle.y, 3) + 0.025,
       # sprintf is ctype float coding (ensures only 1 decimal place)
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=1)
  
  # This puts the confidence intervals below the means
  text(x = c(circle.cent + c(-0.14, 0.15, 0.16)),
       y= rep(circle.y - 0.03, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), "-",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=0.9)
  
  # This portion of code include a figure label (i.e., 1A)
  # Note Relative axis point converts x y points into percentages 
  # and by changing lims, the code puts them in same place
  # text(x=relative.axis.point(0.1, "x"),
       #y=relative.axis.point(0.95, "y"),
       #labels="(A)", font=2, cex=0.9)
  
  # Close pdf plot environment
  dev.off()
  
}





