# A function which will transform a list of time series into list of matrices
# This will be used for the taxonomic portion of the pilot analysis

# Define function
taxonomic.matrix.maker <- function(surveydata, ts_data){
  
  # Verification variables
  # surveydata <- pilot_survey_data
  # ts_data <- pilot_ts_data
  
  # Select TS_ID as a list to iterate through
  TS_ID <- subset(ts_data, select = "TimeSeriesID")
    
  # Create empty list to store results
  matrix_list <- list()
  
  # Loop through survey data according to time series ID
  for(i in 1:length(TS_ID$TimeSeriesID)){
    
      # Select individual time series data
      TS_data <- surveydata %>%
        filter(surveydata$TimeSeriesID %in% TS_ID[i,])

      # Create species abundance matrix for timeseries ID
      TS_A_data.mat <- tapply(TS_data$Abundance, 
                            list(TS_data$Year, # This part will need to change 
                            TS_data$Species), 
                            mean, na.rm=TRUE)
      # Replace NA with 0
      TS_A_data.mat[is.na(TS_A_data.mat)] <- 0
      
      # Convert abundances into relative abundance # Note sqrt, accounts for rare taxa
      TS_data.mat <- make_relative(sqrt(TS_A_data.mat))
      
      # Convert matrix to list object
      TS_output <- list(TS_data.mat)
      
      # Add list object to full matrix list
      matrix_list <- append(matrix_list, TS_output)
      
      # Change name in list according to TS ID
      names(matrix_list)[i] <- TS_ID[i,]
  }
    
  # Return list of matrices
  return(matrix_list)
}
  
  
  