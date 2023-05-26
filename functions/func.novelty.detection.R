# A function which will apply the novelty detection framework to a list of matrices
# This will be used for the functional portion of the pilot analysis

# Define function
func.novelty.detection <- function(tax_A_matrix_list, diss_matrix_list){
  
  # Required packages
  require(svMisc)
  
  # Verification variables
  #tax_A_matrix_list <- functional_A_matrix_list
  #diss_matrix <- pilot_BL_dissimilarity.mat
  
  # Create empty results df
  func_results <- list()
  
  # Loop through list to apply framework to each matrix
  for(i in 1:length(tax_A_matrix_list)){
    
    # Identify matrix and TS ID to be used
    TS_ID <- names(tax_A_matrix_list[i])
    TS_mat <- tax_A_matrix_list[[i]]
    
    # Save plot to plots folder named according to TS ID
    plot_filename <- paste0("./plots/Functional Plots/", TS_ID, "_plot.pdf")
    pdf(plot_filename, width = 10, height = 5)
    
    # Apply framework using identify.novelty.gam function
    novelty_output <- novelTsComms(timeseries = TS_mat, 
                      timeIds = rownames(TS_mat),
                      dissWeight = diss_matrix_list,
                      alpha = 0.05, 
                      method = "weighted", # metric in old
                      abundWeight = TRUE,
                      modelType = "lognormal",
                      plot = TRUE, 
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
    func_results <- append(func_results, novelty_output)
      
    # Show progress
    progress(i, length(diss_matrix_list))
    }
    
  
  return(func_results)
}


