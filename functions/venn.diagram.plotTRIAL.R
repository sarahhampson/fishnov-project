# A function which takes the output of binomial regression model and plots the venn diagrams

venn.diagram.plotTRIAL <- function(GLM.prob.results, plot.name){
  
  # Verification variables
  plot.name <- "TRIAL PLOT"
  GLM.prob.results <- tax_prob_results
  
  # Load required packages
  require(plotrix)
  require(sp)
  require(raster)
  require(rgeos) 
  require(shape)
  
  # Retrieve probabilities of each category
  back.prob <- plogis(tax_prob_results[["back"]][["coefficients"]])
  instant.prob <- plogis(tax_prob_results[["instant"]][["coefficients"]])
  cumul.prob <- plogis(tax_prob_results[["cumul"]][["coefficients"]])
  novel.prob <- plogis(tax_prob_results[["novel"]][["coefficients"]])
  
  plot.new() # call an empty plot
  box("plot", col="red") # add boxes around the plot area
  box("figure", col="blue") # add box around figure area
  box("outer", col="green") # add box around outer area
}
  