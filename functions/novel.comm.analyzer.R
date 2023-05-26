# # This functions creates summary information from the novelty results 
# and returns a data frame with community state information

novel.comm.summary <- function(results_list, ts_list){
  
  # Verification variables
  results_list <- pilot_tax_results
  ts_list <- pilot_ts_list
  
  # Loop through results list of lists and extract information
  for(i in 1:length(results_list)){
    # Extract "cat" vector from results list
    comm.states <- results_list[[i]][[1]][["cat"]]
    # Identify variables of interest and add to ts_list
    ts_list$novel_count[i] <- sum(comm.states == "novel")
    ts_list$back_count[i] <- sum(comm.states == "back")
    ts_list$instant_count[i] <- sum(comm.states == "instant")
    ts_list$cumul_count[i] <- sum(comm.states == "cumul")
    
    # Total number of states across all 
  }
  return(ts_list)
}
  
  
