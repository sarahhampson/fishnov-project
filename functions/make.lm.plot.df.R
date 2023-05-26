# Function to create a plot df from a linear model
make.lm.plot.df <- function(model, log=NULL){
  
  # Verification vairables
  #model <- tax.lm.Absize 
  
  summary_df <- data.frame(
    Effect_Size = numeric(),
    Mean_Estimate = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    PI_Lower = numeric(),
    PI_Upper = numeric(),
    row.names = character()
  )
  
  # Set intercept
  intercept = coefficients(model)[1]
  
  # Loop through each fixed effect
  for (i in 1:length(coefficients(model))) {
    # Calculate the effect size
    effect_size <- abs(coefficients(model)[i])
    
    # Calculate the mean estimate
    mean_estimate <- coefficients(model)[i]
    
    # Calculate the confidence interval
    ci <- confint(model, level = 0.95)[i,]
    
    # Calculate the predictive interval
    pi <- predict(model, interval = "prediction")[,2:3][i,]
    
    # Add the results to the summary dataframe
    summary_df[i,] <- list(
      Effect_Size = effect_size,
      Mean_Estimate = mean_estimate,
      CI_Lower = ci[1],
      CI_Upper = ci[2],
      PI_Lower = pi[1],
      PI_Upper = pi[2]
    )
    
    # Set the row names to the fixed effect names
    row.names(summary_df)[i] <- names(coefficients(model))[i]
    
  }
  
  # Add the intercept to the mean estimates and CI for the binary fixed effect
  for(i in 2:4){
    summary_df[2,i] = summary_df[2,i]+intercept
  }
  
  # Print the summary dataframe
  return(summary_df)
}
