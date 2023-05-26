# Function to create a list of time series that qualify for novelty analysis
# Requirements are that there must be at least 10 unique years in the dataset, 
# and they must have greater than 2 species

# Function which gives filtered df for survey data and time series data
select.timeseries <- function(survey_data, time_series_data){
  
  # Remove extra spaces from the quarter column
  survey_data <- survey_data %>%
    mutate(Quarter = str_replace_all(Quarter, " ", ""))
  
  # If any 1/2 quarters appear replace the '/' with a '-'
  survey_data$Quarter <- gsub('/', '-', survey_data$Quarter)
  
  # Table of timeseries with min 10 consecutive yearly surveys and 2 species
  ts_table <- survey_data %>% group_by(TimeSeriesID) %>%
    dplyr::summarise(year_count = length(unique(Year)), first_year = min(Year), 
                     last_year = max(Year), first_to_last = last_year - first_year, 
                     full_data = (year_count == first_to_last+1), 
                     species_count = length(unique(Species)), 
                     quarter_count = length(unique(Quarter))) %>%
    filter(year_count >=10, species_count > 1) # removed filter full_data == TRUE
  
  # Modify survey data to only include timeseries from ts_table
  survey_data <- survey_data[survey_data$TimeSeriesID %in% 
                               ts_table$TimeSeriesID,]
  
  # Many TS are sampled in more than 1 quarter so I will filter out each quarter 
  # per time series (i.e. nested Quarter TS in one site)
  ts_q_table <- survey_data %>% group_by(TimeSeriesID, Quarter) %>%
    dplyr::summarise(year_count = length(unique(Year)),
                     species_count = length(unique(Species))) %>%
    filter(year_count >=10, species_count > 1)
  
  # Modify survey data to only include timeseries + quarters from ts_q_table
  survey_data <- merge(survey_data, ts_q_table, 
                       by.x = c("TimeSeriesID", "Quarter"))
  
  # Check quarter count frequency for the TS
  ts_data_qcount <- survey_data %>% group_by(TimeSeriesID) %>%
    dplyr::summarise(qcount = n_distinct(Quarter))
  
  # Check for frequency of TS with multiple surveys per quarter
  survey_quarter_table <- survey_data %>% group_by(TimeSeriesID, Year, Quarter) %>%
    dplyr::summarise(surveycount = length(unique(SurveyID)))
  # --> therefore I need to use the mean for abundance in quarterly time bins
  
  # Add site quarter index column
  survey_data <- survey_data %>% 
    mutate(sqID = paste0(TimeSeriesID, " Q", Quarter))
  
  return(survey_data)
}
