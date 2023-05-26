# Function to take a TS and convert into list of seasonal species abundance matrices
seasonal.matrix.maker <- function(surveydata, ts_data){
  
  # Required packages
  require(funrar)
  require(svMisc)
  
  # Verification variables
  #surveydata = pilot_survey_data
  #ts_data = pilot_ts_data
  
  # Create list for merge including number of quarters
  quarter_data <- surveydata %>% group_by(TimeSeriesID, Quarter) %>%
    dplyr::summarise(qcount = n_distinct(Year))
  
  # Merge with ts_data so each quarter is a unique TS entry
  ts_q_data <- merge(quarter_data, ts_data)
  
  # Create unique ID for each TS + quarter in both datasets
  ts_q_data$TS_qID <- str_c(ts_q_data$TimeSeriesID, ' Q', ts_q_data$Quarter)
  surveydata$TS_qID <- str_c(surveydata$TimeSeriesID, ' Q', surveydata$Quarter)
  
  
  # Create empty list to store matrices of all TS
  matrix_list <- list()
  
  # Loop through survey data according to time series ID
  for(i in 1:length(ts_q_data$TS_qID)){
    
    # Select individual time series data
    TS_data <- surveydata %>%
        dplyr::filter(surveydata$TS_qID %in% ts_q_data$TS_qID[i])
    
    # Create species abundance matrix for timeseries ID
    TS_A_data.mat <- tapply(TS_data$Abundance, 
                            list(TS_data$Year, 
                            TS_data$Species), 
                            mean, na.rm=TRUE) # Using mean to account for multiquarter surveys
    # Replace NA with 0
    TS_A_data.mat[is.na(TS_A_data.mat)] <- 0
    
    # Convert abundances into relative abundance # Note sqrt, accounts for rare taxa
    TS_data.mat <- make_relative(sqrt(TS_A_data.mat))
    
    # Convert matrix to list object
    TS_output <- list(TS_data.mat)
    
    # Add list object to full matrix list
    matrix_list <- append(matrix_list, TS_output)
    
    # Change name in list according to TS ID
    names(matrix_list)[i] <- ts_q_data$TS_qID[i]
    
    # Show progress
    progress(i, length(ts_q_data$TS_qID))
    #print(paste(round((i/length(ts_q_data$TS_qID)*100), 1), "% complete"))
    }

    # Return list of matrices
    return(matrix_list)
  }
