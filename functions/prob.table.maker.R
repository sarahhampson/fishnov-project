# A function that takes probability results and converts into a table

prob.table.maker <- function(prob_results){
  
  # Verification variables
  # prob_results = func_prob_results
  
  # Extract coefficients data
  prob.coefs <- lapply(1:length(prob_results), function(x){
    summary(prob_results[[x]])$coefficients})
  names(prob.coefs) <- names(prob_results)
  
  # Convert into df
  temp.df <- data.frame(matrix(unlist(prob.coefs), nrow = length(prob.coefs), byrow = T))
  temp.df[, 5] <- names(prob_results)
  colnames(temp.df) <- c("Intercept", "SE", "Z score", "P value", "Category")
  
  # Create empty dataframe with names 
  plot.df <- data.frame(matrix(ncol=4, nrow=4))
  colnames(plot.df) <- c("Category", "Intercept",
                         "CI Lower", "CI Upper")
  
  # Input data from temp.df for cat and log transformed intercept
  plot.df$Category = temp.df$Category
  plot.df$Intercept = plogis(temp.df$Intercept)
  plot.df$`CI Lower` = plogis(temp.df$Intercept-temp.df$SE*1.96)
  plot.df$`CI Upper` = plogis(temp.df$Intercept+temp.df$SE*1.96)
  
  # Return the table dataframe
  return(plot.df)
  
}
