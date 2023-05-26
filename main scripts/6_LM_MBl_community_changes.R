# Title: LM models for MBl changes in community mean and sd
# Author: Sarah Hampson


# 0.  Load necessary dataframes -------------------------------------------

# Read results csvs
#tax_results_df <- read.csv("./detection results/tax_results_df.csv")
#func_results_df <- read.csv("./detection results/func_results_df.csv")
full_results_df <- read.csv("./detection results/full_results_df.csv")
trait_data <- read.csv("./inputs/data/trait_data.csv")

# Need taxonomic matrix list, if not already in environment run this
survey_data <- read.csv("./inputs/data/survey_data.csv")
ts_data <- read.csv("./inputs/data/ts_data.csv")
taxonomic_matrix_list <- seasonal.matrix.maker(survey_data, ts_data)

# 1.  Retrieve community size data  ---------------------------------------

# Includes:
# absolute change in community mean from one time point to the next
# raw change in community mean from one time point to the next
# absolute change in standard deviation from one time point to the next
# raw change in standard deviation from one time point to the next

# All of these variables are abundance weighted!
# Note that size (MBl) has been log transformed!

# Make list of communities and include year it emerges and the year of data beforehand
comms <- full_results_df %>% 
  mutate(bins.before = as.numeric(bins)-as.numeric(bin.lag)) %>% 
  dplyr::select(sqID, bins, bins.before)

# Use this function to extract the abundance weighted community mean trait
comms_size_data_list <- find.commS.trait.data(comms, taxonomic_matrix_list, trait_data, "MBl")

# 2.  Add and scale fixed effects -----------------------------------------

temp_df <- lapply(split(full_results_df, f=full_results_df$sqID),
                  function(x){
                    x$TS_Length <- nrow(x)
                    x$Yrs_passed <- x$bins - min(x$bins)
                    x$Yr_pos <- 1:nrow(x)
                    return(x)
                  })
full_results_df <- do.call("rbind", temp_df)

# Scale these new variables
full_results_df$bin.lagS <- scale(full_results_df$bin.lag)[,1]
full_results_df$Yrs_passedS <- scale(full_results_df$Yrs_passed)[,1]
full_results_df$logYr_posS <- scale(log(full_results_df$Yr_pos))[,1]

# 3.  Merge with detection results ----------------------------------------

# Row bind these into tax and func comm size data
comms_size_data <- do.call("rbind", comms_size_data_list)
comms_size_df <- data.frame(comms_size_data)

# Change character columns to numeric
comms_size_df[, 2:9] <- lapply(2:9, function(x){
  as.numeric(comms_size_df[[x]])
})

# Merge with full results df
full_results_df_size <- merge(size_results_df,comms_size_df, 
                              by.y=c("sqID", "YearofNov"),
                              by.x=c("sqID", "bins"))

# 4.  Save community size dataframe ---------------------------------------

# Save community size data into the model_data -> trait_models folder
write.csv(full_results_df_size, "./inputs/model_data/trait_models/MBl_df.csv")
comms_size_df <- read.csv("./inputs/model_data/trait_models/MBl_df.csv")

# 5.  Size changes ~ tax novelty models -----------------------------------

# Model trait change as a function of whether community is taxnovel or not
# Including same fixed/random effects as in other time point models
tax.lm.size <- lm(traitchange ~ taxnovel + bin.lagS,
                  data = comms_size_df)

# I have to remove data which have absolute zero change
comms_size_df_Absize <- comms_size_df %>% 
  dplyr::select(bin.lagS, taxnovel, funcnovel, Abtraitchange) %>% 
  filter(Abtraitchange>0)
# 0.0567 --> ~ 6% of datapoints had zero change in size
# Model absolute trait change as a function of whether community is
# tax novel or not
tax.lm.Absize <- lm(log(Abtraitchange) ~ taxnovel + bin.lagS,
                  data = comms_size_df_Absize)

# Model change in standard deviation of size as a function of whether
# community is tax novel or not
tax.lm.sdsize <- lm(sdchange ~ taxnovel + bin.lagS,
                    data = comms_size_df)
# See summary
summary(tax.lm.size)
summary(tax.lm.Absize)
summary(tax.lm.sdsize)

# 6.  Size changes ~ func novelty models ----------------------------------

# Model trait change as a function of whether community is funcnovel or not
size.lm.size <- lm(traitchange ~ funcnovel + bin.lagS,
                  data = comms_size_df)
size.lm.Absize <- lm(log(Abtraitchange) ~ funcnovel + bin.lagS,
                    data = comms_size_df_Absize)
size.lm.sdsize <- lm(sdchange ~ funcnovel + bin.lagS,
                    data = comms_size_df)
# See summary
summary(size.lm.size)
summary(size.lm.Absize)
summary(size.lm.sdsize)

# 7.  Create and save plot dataframes -------------------------------------

# Use function to create plot dataframes for each model
tax_size_plot_df <- make.lm.plot.df(tax.lm.size)
tax_Absize_plot_df <- make.lm.plot.df(tax.lm.Absize, log=T)
tax_sdsize_plot_df <- make.lm.plot.df(tax.lm.sdsize)
size_size_plot_df <- make.lm.plot.df(size.lm.size)
size_Absize_plot_df <- make.lm.plot.df(size.lm.Absize, log=T)
size_sdsize_plot_df <- make.lm.plot.df(size.lm.sdsize)

# Save the dataframes
write.csv(tax_size_plot_df, "./outputs/MBl/tax_lm_size.csv")
write.csv(tax_Absize_plot_df, "./outputs/MBl/tax_lm_Absize.csv")
write.csv(tax_sdsize_plot_df, "./outputs/MBl/tax_lm_sdsize.csv")
write.csv(size_size_plot_df, "./outputs/MBl/size_lm_size.csv")
write.csv(size_Absize_plot_df, "./outputs/MBl/size_lm_Absize.csv")
write.csv(size_sdsize_plot_df, "./outputs/MBl/size_lm_sdsize.csv")

# 8.  Cleanup ------------------------------------------------------------

rm(full_results_df, trait_data, taxonomic_matrix_list, survey_data, ts_data, 
   comms_size_data_list, comms_size_data, comms_size_df, tax.lm.size,
   tax.lm.Absize, tax.lm.sdsize, tax_size_plot_df, tax_Absize_plot_df,
   tax_sdsize_plot_df, func.lm.size, func.lm.Absize, func.lm.sdsize,
   func_size_plot_df, func_Absize_plot_df, func_sdsize_plot_df,
   full_results_df_size, full_results_df_Absize)






