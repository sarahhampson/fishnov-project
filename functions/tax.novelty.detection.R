# A function which will apply the novelty detection framework to a list of matrices
# This will be used for the taxonomic portion of the pilot analysis

# Define function
tax.novelty.detection <- function(matrix_list){
  
  # Required packages
  require(svMisc)
  
  # Verification variables
  #matrix_list <- taxonomic_matrix_list
  
  # Create empty results list
  tax_results <- list()
  
  # Loop through list to apply framework to each matrix
  for(i in 1:length(matrix_list)){
    
    # Identify matrix and TS ID to be used
    TS_mat <- matrix_list[[i]]
    TS_ID <- names(matrix_list[i])
  
    # Save plot to plots folder named according to TS ID
    plot_filename <- paste0("./plots/Taxonomic Plots/", TS_ID, "_plot.pdf")
    pdf(plot_filename, width = 10, height = 5)
    
    # Apply framework using identify.novelty.gam function
    novelty_output <- identify.novel.gam(TS_mat,
                                   #timeIds = rownames(TS_mat),
                                   alpha = 0.05,
                                   site = TS_ID,
                                   metric="bray", #method in new
                                   #dissWeight = NULL,
                                   #abundWeight = FALSE,
                                   #modelType = "beta", 
                                   plot = TRUE, 
                                   plot.data=FALSE,
                                   gam.max.k = -1) 
  
    # Close current plot device (i.e. the pdf being written)
    dev.off()
    # Remove first 5 time bin
    novelty_output <- novelty_output[-c(1:5),]
    # Make novelty output into list
    novelty_output <- list(novelty_output)
    # Name according to TS ID
    names(novelty_output) <- TS_ID
    # Append output to list
    tax_results <- append(tax_results, novelty_output)
  
    # Show progress
    progress(i, length(matrix_list))
    
  }
  return(tax_results)
}
