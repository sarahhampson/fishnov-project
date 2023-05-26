# Title: Find community state changes in multi-trait space
# Author: Sarah Hampson

# 0.  Load necessary dataframes -------------------------------------------

# Load imputed PC trait data
PCtrait_data <- read.csv("./inputs/data/trait_and_PCA_data.csv")
PCtrait_data <- PCtrait_data %>% dplyr::select(-X)
#trait_data <- trait_data[order(trait_data$Species),]

# Subset taxonomic matrix list for multitrait only matrices
subset_matrix_list <- taxonomic_matrix_list[names(taxonomic_matrix_list) %in% multitrait_results_df$sqID]

# 1.  Find PCA axis trait changes ------------------------------------------

# Make list of communities and include year it emerges and the year of data beforehand
comms <- full_results_df %>% 
  mutate(bins.before = as.numeric(bins)-as.numeric(bin.lag)) %>% 
  dplyr::select(sqID, bins, bins.before)

# Abundance weighted PC1 and PC2 values for CWM and sd
comms_PC1_data_list <- find.commS.trait.data(comms, taxonomic_matrix_list, PCtrait_data, "Dim.1")
comms_PC2_data_list <- find.commS.trait.data(comms, taxonomic_matrix_list, PCtrait_data, "Dim.2")

# 2.  Merge with detection results ----------------------------------------

# Row bind these into tax and func comm size data
comms_PC1_data <- do.call("rbind", comms_PC1_data_list)
comms_PC1_df <- data.frame(comms_PC1_data)
comms_PC2_data <- do.call("rbind", comms_PC2_data_list)
comms_PC2_df <- data.frame(comms_PC2_data)

# Change character columns to numeric
comms_PC1_df[, 2:9] <- lapply(2:9, function(x){
  as.numeric(comms_PC1_df[[x]])
})
comms_PC2_df[, 2:9] <- lapply(2:9, function(x){
  as.numeric(comms_PC2_df[[x]])
})

# Merge with multitrait results df
tax_results_df <- full_results_df %>% 
  dplyr::select(Site, bins, bin.lag, sqID, Quarter, taxnovel)
func10_results_df <- func_results_df %>% 
  dplyr::select(Site, bins, bin.lag, sqID, Quarter, funcnovel, Country, HydroBasin)
multitrait_results_df <- merge(tax_results_df, func10_results_df, 
                               by=c("Site", "bins", "bin.lag", "sqID", "Quarter"))

multitrait_results_df_PC1 <- merge(multitrait_results_df, comms_PC1_df,
                               by.y=c("sqID", "YearofComm"),
                               by.x=c("sqID", "bins"))
multitrait_results_df_PC2 <- merge(full_results_df, comms_PC2_df,
                                   by.y=c("sqID", "YearofComm"),
                                   by.x=c("sqID", "bins"))

# 3.  Save community PC1 and PC2 dataframe ---------------------------------------

# Save community troph data into the model_data -> trait_models folder
write.csv(multitrait_results_df_PC1, "./inputs/model_data/trait_models/PC1_df.csv")
write.csv(multitrait_results_df_PC2, "./inputs/model_data/trait_models/PC2_df.csv")

# 4.  Add and scale fixed effects -----------------------------------------

temp_df <- lapply(split(multitrait_results_df_PC1, f=multitrait_results_df_PC1$sqID),
                  function(x){
                    x$TS_Length <- nrow(x)
                    x$Yrs_passed <- x$bins - min(x$bins)
                    x$Yr_pos <- 1:nrow(x)
                    return(x)
                  })
multitrait_results_df_PC1 <- do.call("rbind", temp_df)
temp_df <- lapply(split(multitrait_results_df_PC2, f=multitrait_results_df_PC2$sqID),
                  function(x){
                    x$TS_Length <- nrow(x)
                    x$Yrs_passed <- x$bins - min(x$bins)
                    x$Yr_pos <- 1:nrow(x)
                    return(x)
                  })
multitrait_results_df_PC2 <- do.call("rbind", temp_df)
rm(temp_df)

# Scale these new variables
multitrait_results_df_PC1$bin.lagS <- scale(multitrait_results_df_PC1$bin.lag)[,1]
multitrait_results_df_PC1$Yrs_passedS <- scale(multitrait_results_df_PC1$Yrs_passed)[,1]
multitrait_results_df_PC1$logYr_posS <- scale(log(multitrait_results_df_PC1$Yr_pos))[,1]
multitrait_results_df_PC2$bin.lagS <- scale(multitrait_results_df_PC2$bin.lag)[,1]
multitrait_results_df_PC2$Yrs_passedS <- scale(multitrait_results_df_PC2$Yrs_passed)[,1]
multitrait_results_df_PC2$logYr_posS <- scale(log(multitrait_results_df_PC2$Yr_pos))[,1]

# 5.  PC1 changes ~ tax novelty models  -----------------------------------

# Model trait change as a function of whether community is taxnovel or not
# Including same fixed/randomeffects as in other time point models
tax.lm.PC1 <- lm(traitchange ~ taxnovel + bin.lagS, 
                   data = multitrait_results_df_PC1)

# I have to remove data which have absolute zero change
multitrait_results_df_AbPC1 <- multitrait_results_df_PC1 %>% 
  filter(Abtraitchange>0)
# 0.0567 --> ~ 6% of datapoints had zero change in PC1 (i.e. community was exact same)
# Model absolute trait change as a function of whether community is
# tax novel or not
tax.lm.AbPC1 <- lm(log(Abtraitchange) ~ taxnovel + bin.lagS, 
                     data = multitrait_results_df_AbPC1)

# Model change in standard deviation of PC1 as a function of whether
# community is tax novel or not
tax.lm.sdPC1 <- lm(sdchange ~ taxnovel + bin.lagS, 
                     data = multitrait_results_df_PC1)
# See summary
summary(tax.lm.PC1)
summary(tax.lm.AbPC1)
summary(tax.lm.sdPC1)

# 6.  PC1 changes ~ multitrait novelty models ----------------------------------

# Model trait change as a function of whether community is funcnovel or not
multitrait.lm.PC1 <- lm(traitchange ~ funcnovel + bin.lagS, 
                     data = multitrait_results_df_PC1)
multitrait.lm.AbPC1 <- lm(log(Abtraitchange) ~ funcnovel + bin.lagS, 
                       data = multitrait_results_df_AbPC1)
multitrait.lm.sdPC1 <- lm(sdchange ~ funcnovel + bin.lagS, 
                       data = multitrait_results_df_PC1)
# See summary
summary(multitrait.lm.PC1)
summary(multitrait.lm.AbPC1)
summary(multitrait.lm.sdPC1)

# 7.  PC2 changes ~ tax novelty models  -----------------------------------

tax.lm.PC2 <- lm(traitchange ~ taxnovel + bin.lagS, 
                 data = multitrait_results_df_PC2)
multitrait_results_df_AbPC2 <- multitrait_results_df_PC2 %>% 
  filter(Abtraitchange>0)
tax.lm.AbPC2 <- lm(log(Abtraitchange) ~ taxnovel + bin.lagS, 
                   data = multitrait_results_df_AbPC2)
tax.lm.sdPC2 <- lm(sdchange ~ taxnovel + bin.lagS, 
                   data = multitrait_results_df_PC2)

summary(tax.lm.PC2)
summary(tax.lm.AbPC2)
summary(tax.lm.sdPC2)

# 8.  PC2 changes ~ multitrait novelty models ----------------------------------

# Model trait change as a function of whether community is funcnovel or not
multitrait.lm.PC2 <- lm(traitchange ~ funcnovel + bin.lagS, 
                        data = multitrait_results_df_PC2)
multitrait.lm.AbPC2 <- lm(log(Abtraitchange) ~ funcnovel + bin.lagS, 
                          data = multitrait_results_df_AbPC2)
multitrait.lm.sdPC2 <- lm(sdchange ~ funcnovel + bin.lagS, 
                          data = multitrait_results_df_PC2)
# See summary
summary(multitrait.lm.PC2)
summary(multitrait.lm.AbPC2)
summary(multitrait.lm.sdPC2)

# 9. Create and save plot dataframes -------------------------------------

# Use function to create plot dataframes for each model
#PC1
tax_PC1_plot_df <- make.lm.plot.df(tax.lm.PC1)
tax_AbPC1_plot_df <- make.lm.plot.df(tax.lm.AbPC1, log=T)
tax_sdPC1_plot_df <- make.lm.plot.df(tax.lm.sdPC1)
multitrait_PC1_plot_df <- make.lm.plot.df(multitrait.lm.PC1)
multitrait_AbPC1_plot_df <- make.lm.plot.df(multitrait.lm.AbPC1, log=T)
multitrait_sdPC1_plot_df <- make.lm.plot.df(multitrait.lm.sdPC1)
#PC2
tax_PC2_plot_df <- make.lm.plot.df(tax.lm.PC2)
tax_AbPC2_plot_df <- make.lm.plot.df(tax.lm.AbPC2, log=T)
tax_sdPC2_plot_df <- make.lm.plot.df(tax.lm.sdPC2)
multitrait_PC2_plot_df <- make.lm.plot.df(multitrait.lm.PC2)
multitrait_AbPC2_plot_df <- make.lm.plot.df(multitrait.lm.AbPC2, log=T)
multitrait_sdPC2_plot_df <- make.lm.plot.df(multitrait.lm.sdPC2)

# Save the dataframes
write.csv(tax_PC1_plot_df, "./outputs/multitrait/tax_lm_PC1.csv")
write.csv(tax_AbPC1_plot_df, "./outputs/multitrait/tax_lm_AbPC1.csv")
write.csv(tax_sdPC1_plot_df, "./outputs/multitrait/tax_lm_sdPC1.csv")
write.csv(multitrait_PC1_plot_df, "./outputs/multitrait/multitrait_lm_PC1.csv")
write.csv(multitrait_AbPC1_plot_df, "./outputs/multitrait/multitrait_lm_AbPC1.csv")
write.csv(multitrait_sdPC1_plot_df, "./outputs/multitrait/multitrait_lm_sdPC1.csv")
write.csv(tax_PC2_plot_df, "./outputs/multitrait/tax_lm_PC2.csv")
write.csv(tax_AbPC2_plot_df, "./outputs/multitrait/tax_lm_AbPC2.csv")
write.csv(tax_sdPC2_plot_df, "./outputs/multitrait/tax_lm_sdPC2.csv")
write.csv(multitrait_PC2_plot_df, "./outputs/multitrait/multitrait_lm_PC2.csv")
write.csv(multitrait_AbPC2_plot_df, "./outputs/multitrait/multitrait_lm_AbPC2.csv")
write.csv(multitrait_sdPC2_plot_df, "./outputs/multitrait/multitrait_lm_sdPC2.csv")

# 10.  Cleanup ------------------------------------------------------------

rm(comms, comms_PC1_data_list, comms_PC2_data_list, comms_PC1_data, comms_PC2_data,
   comms_PC1_df, comms_PC2_df, multitrait_results_df_PC1, multitrait_results_df_PC2,
   multitrait_results_df_AbPC1, multitrait_results_df_AbPC2, tax.lm.PC1, tax.lm.AbPC1,
   tax.lm.sdPC1, tax.lm.PC2, tax.lm.AbPC2, tax.lm.sdPC2, multitrait.lm.PC1, 
   mulitrait.lm.AbPC1, multitrait.lm.sdPC1, multitrait.lm.PC2, multitrait.lm.AbPC2,
   multitrait.lm.sdPC2, tax_PC1_plot_df, tax_AbPC2_plot_df, tax_sdPC1_plot_df,
   tax_PC2_plot_df, tax_AbPC2_plot_df, tax_sdPC2_plot_df, multitrait_PC1_plot_df, 
   multitrait_AbPC1_plot_df, multitrait_sdPC1_plot_df, multitrait_PC2_plot_df,
   multitrait_AbPC2_plot_df, multitrait_sdPC2_plot_df)

# OLD OPTION.  Find comm FD data  ---------------------------------------------------

# OLD WAY --> not good
FD_change_list <- find.comm.FD.data(subset_matrix_list, trait_data)

# Test below for one abundance matrix 
test_relA_df <- data.frame(subset_matrix_list[["G166 Q4"]], check.names = FALSE)
relA_species <- colnames(test_relA_df)
print(test_relA_df)

# Drop the species from the trait_data that aren't in the abundance data
test_trait_data <- trait_data %>% filter(Species %in% relA_species)
rownames(test_trait_data) <- test_trait_data$Species
test_trait_data <- test_trait_data %>% dplyr::select(-X, -Species)

test_FD <- dbFD(test_trait_data,
             test_relA_df,
             w.abun=T,
             stand.x=T, # Standardise variables
             stand.FRic=T, # Standardises f richness, better for comparing FRic
             calc.FRic=T,
             m=5,
             calc.CWM=T,
             calc.FDiv = T,
             print.pco=T
             )

test <- data.frame(cbind(test_FD$FRic, test_FD$FEve, test_FD$FDiv, test_FD$FDis, test_FD$RaoQ))
colnames(test) <- c("FRic", "FEve", "FDiv", "FDis", "RaoQ")
test$Year <- rownames(test)
test

rm(test_trait_data, test_relA, relA_species, test_FD)






