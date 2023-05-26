# Function to turn glm output into a plotable dataframe or just glm output

# Function to make dataframe of coefficients
make.glm.df <- function(glm, category, output){
  
  # Verification variable
  #glm <- tax.glm.tp
  #category <- "Tax"
  #output <- "plot"
  
  # Summarise glm and extract coefficients as dataframe
  model_df <- as.data.frame(coef(summary(glm)))
  
  # Add model name to df
  model_df$model_name <- deparse(substitute(glm))
  
  # Add effects row
  model_df$effect <- rownames(model_df)
  
  # Set intercept
  model_df$mean_est <- model_df$Estimate
  
  # If length of coeffs is greater than one, add them together to get means estimate
  if(nrow(model_df >1)){
    for(i in 2:nrow(model_df)){
      model_df$mean_est[i] = model_df$Estimate[i] + model_df$mean_est[i-1]
    }
  }
  
  # Create empty dataframe with names 
  plot.df <- data.frame(matrix(ncol=7, nrow=nrow(model_df)))
  colnames(plot.df) <- c("Model name", "Estimate",
                         "Effect", "Mean Estimate",
                         "CI Lower", "CI Upper", 
                         "Category")
  if(output=="plot"){
    
    for(i in 1:nrow(model_df)){
      
      # Input data from temp.df for cat and log transformed intercept
      plot.df$`Model name`[i] = model_df$model_name[i]
      plot.df$Category[i] = category
      plot.df$Effect[i] = model_df$effect[i]
      plot.df$Estimate[i] = model_df$Estimate[i]
      plot.df$`Mean Estimate`[i] = plogis(model_df$mean_est[i])
      plot.df$`CI Lower`[i] = plogis(model_df$mean_est[i]-model_df$`Std. Error`[i]*1.96)
      plot.df$`CI Upper`[i] = plogis(model_df$mean_est[i]+model_df$`Std. Error`[i]*1.96)
    }
    # Return the table dataframe
    return(plot.df)}
  
  if(output=="model"){
    # Change model coefs output to real probabilities
    return(model_df)
  }}
