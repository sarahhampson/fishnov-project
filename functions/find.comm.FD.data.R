# Functions for finding community functional diversity data

# Go through abundance matrix list and see if the matrix list ID matches Comms ID
find.comm.FD.data <- function(matrix_list, trait_df){
  
  require(FD)
  require(tidyverse)
  
  # Create empty output list
  output_list <- list()
  
  # Loop through matrix list
  for (i in 1:length(matrix_list)){
  
    # Verification variables
    #matrix_list <- subset_matrix_list
    #trait_df <- impPCA_data
    
    # Select matrix
    matrix <- matrix_list[[i]]
    
    # Create ID name for this matrix
    ID <- names(matrix_list[i])
    print(ID)
    
    # Convert matrix to a dataframe
    relA_mat <- data.frame(matrix, check.names = FALSE)
    
    # Isolate species
    sp <- colnames(relA_mat)
    
    if (length(sp) > 2){
    
      # Drop species from the trait df that aren't in the matrix
      sub_trait_df <- trait_df %>% filter(Species %in% sp)
      rownames(sub_trait_df) <- sub_trait_df$Species
      sub_trait_df <- sub_trait_df %>% dplyr::select(-Species)
      
      # Use FD package to calculate FD indices
      FD_list <- dbFD(sub_trait_df,
                      relA_mat,
                      w.abun=T, # Weigh FD indices by abundance
                      stand.x=T, # Standardise variables
                      stand.FRic=T, # Standardises f richness, better for comparing FRic
                      calc.FRic=T,
                      m=5,
                      calc.CWM=T,
                      calc.FDiv = T)
                      #print.pco=T)
      
      # Create dataframe with raw outputs to start
      temp_df <- data.frame(cbind(FD_list$nbsp,
                                    FD_list$FRic, 
                                    FD_list$FEve, 
                                    FD_list$FDiv, 
                                    FD_list$FDis, 
                                    FD_list$RaoQ))
      colnames(temp_df) <- c("SpeciesRic", "FRic", "FEve", "FDiv", "FDis", "RaoQ")
      temp_df$Year <- rownames(temp_df)
      
      # Create empty variables for all FD changes
      temp_df$FRicChange <- "NA"
      temp_df$FEveChange <- "NA"
      temp_df$FDivChange <- "NA"
      temp_df$FDisChange <- "NA"
      temp_df$RaoQChange <- "NA"
      
      # Calculate the change in diversity measures excluding first five entries
      for (j in 6:nrow(temp_df)){
        temp_df$FRicChange[j] = temp_df$FRic[j] - temp_df$FRic[j-1]
        temp_df$FEveChange[j] = temp_df$FEve[j] - temp_df$FEve[j-1]
        temp_df$FDivChange[j] = temp_df$FDiv[j] - temp_df$FDiv[j-1]
        temp_df$FDisChange[j] = temp_df$FDis[j] - temp_df$FDis[j-1]
        temp_df$RaoQChange[j] = temp_df$RaoQ[j] - temp_df$RaoQ[j-1]
      }
      
      # Add sqID and qual of functional richness measures
      temp_df$sqID <- ID
      temp_df$qual.FRic <- FD_list$qual.FRic
    
    }else{
      temp_df <- ID
    }
    
    temp_df_list <- list(temp_df)
    names(temp_df_list) <- ID
    
    # Append this to the output list
    output_list <- append(output_list, temp_df_list)
  }
  return (output_list)
}

  