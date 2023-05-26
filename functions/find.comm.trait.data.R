# Functions for finding community trait data

# Go through abundance matrix list and see if the matrix list ID matches Comms ID
find.comm.trait.data <- function(comm, taxonomic_matrix_list, trait_data, trait){
  
  # Verification variables
  
  # Define some subset variables
  ID <- comm$sqID[1]
  yearofcomm <- as.numeric(comm$bins[1])
  yearb4comm <- as.numeric(comm$bins.before[1])
  
  # Loop through matrix list
  for(i in 1:length(taxonomic_matrix_list)){
    
    mat_ID <- names(taxonomic_matrix_list[i])
    mat <- taxonomic_matrix_list[[i]]
    
    # Testing for now
    #mat <- taxonomic_matrix_list[["G10021 Q4"]]
    #mat_ID = "G10021 Q4"
    
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
      
      # Filter trait_data to only include species and trait of interest
      trait_data <- trait_data %>% dplyr::select(Species, trait = grep(trait, names(trait_data)))
      
      # Merge with size data
      temp2_df <- merge(temp2_df, trait_data, by="Species")
      
      # Add trait weighted abundance column
      temp2_df <- temp2_df %>% 
        mutate(tWAbund = Abundance*(as.numeric(trait))) # actually community weighted means
      
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
