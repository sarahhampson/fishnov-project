# Function which takes taxonomic matrix list and trait dataframes and 
# returns a dataframe of functional diversity indices and changes 
# in those indices between community time points.
# It also gives the community weighted mean trait value (CWM)

find.comm.FD.data <- function(comm, taxonomic_matrix_list, trait_data){
  
  # Require mFD package
  require("mFD")
  
  # Verification variables
  comm <- full_results_df %>% 
    dplyr::mutate(bins.before = as.numeric(bins)-as.numeric(bin.lag)) %>% 
    dplyr::select(sqID, bins, bins.before) %>% 
    filter(sqID=="G10003 Q4")
  taxonomic_matrix_list <- taxonomic_matrix_list
  trait_data <- trait_data
  
  # Define some subset variables
  ID <- comm$sqID[1]
  yearofcomm <- as.numeric(comm$bins[1])
  yearb4comm <- as.numeric(comm$bins.before[1])
  
  # Loop through matrix list
  for(i in 1:length(taxonomic_matrix_list)){
    
    mat_ID <- names(taxonomic_matrix_list[i])
    mat <- taxonomic_matrix_list[[i]]
    
    #mat_ID <- "G10003 Q4"
    #mat <- taxonomic_matrix_list[["G10003 Q4"]]
    
    # Check to see if matrix ID and comm ID match
    if(ID==mat_ID){
      
      # If they do, extract the year of community time point and year before data
      temp_df <- as.data.frame(mat) %>% mutate(Year=rownames(mat))
      year_of_comm <- subset(temp_df, Year==yearofcomm) %>% mutate(cat="yofComm")
      year_b4_comm <- subset(temp_df, Year==yearb4comm) %>% mutate(cat="yb4Comm")
      
      # Turn to long format
      year_of_comm_long <- pivot_longer(year_of_comm, 
                                        cols=-c("Year", "cat"),
                                        names_to="Species",
                                        values_to="Abundance")
      year_b4_comm_long <- pivot_longer(year_b4_comm, 
                                        cols=-c("Year", "cat"),
                                        names_to="Species",
                                        values_to="Abundance")
      
      # Merge these
      temp2_df <- rbind(year_of_comm_long, year_b4_comm_long)
      
      # Filter trait_data to only include species in the abundance matrix
      trait_data <- trait_data %>% dplyr::filter(Species %in% temp2_df$Species)
      
      # Use 
      
      # Summarise mean trait and variance for each year
      temp_df_sum <- temp2_df %>% group_by(Year) %>% 
        dplyr::summarise(meantrait = sum(tWAbund), sdtrait = sd(tWAbund))
      
      # Create new df with all this information and the variables to calculate
      output_df <- matrix(nrow=1, ncol=9)
      colnames(output_df) <- c("sqID", "YearofComm", "Yearb4Comm",
                               "meantrait", "sdtrait",
                               "traitchange", "sdchange",
                               "Abtraitchange", "Absdchange")
      
      # SqID
      output_df[1,1] <- ID
      # Year of comm and before comme
      output_df[1,2] <- as.numeric(yearofcomm)
      output_df[1,3] <- as.numeric(yearb4comm)
      # Mean trait and sd
      output_df[1,4] <- as.numeric(temp_df_sum$meantrait[2])
      output_df[1,5] <- as.numeric(temp_df_sum$sdtrait[2])
      # Change in trait and sd
      output_df[1,6] <- as.numeric(temp_df_sum$meantrait[2] - temp_df_sum$meantrait[1])
      output_df[1,7] <- as.numeric(temp_df_sum$sdtrait[2] - temp_df_sum$sdtrait[1])
      # Absolute differences in trait and sd
      output_df[1,8] <- as.numeric(abs(temp_df_sum$meantrait[2] - temp_df_sum$meantrait[1]))
      output_df[1,9] <- as.numeric(abs(temp_df_sum$sdtrait[2] - temp_df_sum$sdtrait[1]))
      
      
      # Make a list element
      output_list <- list(output_df)
      
      # Change name in list according to TS ID
      names(output_list) <- ID
    }}
  
  return(output_list)
}

# Make function to apply above function on multiple rows
find.commS.trait.data <- function(comms, taxonomic_matrix_list, trait_data, trait){
  
  # Verification variables
  
  # Create empty list for relative abundance community data
  comm_trait_data <- list()   
  
  # Loop through comm ids to get year before and year after abundance data
  for(i in 1:nrow(comms)){
    # Append size data list with output for one comm
    output <- find.comm.trait.data(comms[i,], 
                                   taxonomic_matrix_list,
                                   trait_data, trait)
    print(output)
    comm_trait_data <- append(comm_trait_data, output)
    
  }
  
  return(comm_trait_data)
}
