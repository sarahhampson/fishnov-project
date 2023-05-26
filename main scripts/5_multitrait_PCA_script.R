# Title: PCA for trait data
# Author: Sarah Hampson

# 0.  Load necessary dataframes and packages -------------------------------

survey_data <- read.csv("./inputs/data/survey_data.csv")
ts_data <- read.csv("./inputs/data/ts_data.csv")
trait_data <- read.csv("./inputs/data/trait_data.csv")
full_results_df <- read.csv("./detection results/full_results_df.csv")

# Create a list of relative abundance matrices using matrix maker function
taxonomic_matrix_list <- seasonal.matrix.maker(survey_data, ts_data)

# Subset continuous data and species data
PCA_data <- trait_data %>% 
  dplyr::select(-X, -Species, -troph_level)
sp_data <- trait_data %>% 
  dplyr::select(Species)

# Package for handling PCA with missing data
install.packages(c("missMDA", "FactoMineR"))
library(missMDA)
library(FactoMineR)

# 1.  OLD Normalise variable distributions ------------------------------------

# Because I am doing a non-parametric PCA I don't need to do transformations

# Data need to have a normal distribution for use in PCA
# Looking at all traits visually to inspect for normal distribution
histogram(log(trait_data$MBl), breaks=50) # --> log transform
histogram(log(trait_data$BEl), breaks=50) # --> log transform
histogram(trait_data$VEp, breaks=50) # --> normal
histogram(trait_data$REs, breaks=50) # --> normalish
histogram(trait_data$OGp, breaks=50) # --> normalish with some low values
histogram(log(trait_data$RMl), breaks=50) # --> log transform
histogram(trait_data$BLs, breaks=50) # --> normal
histogram(trait_data$PFv, breaks=50) # --> pretty normal
histogram(trait_data$PFs, breaks=50) # --> normal
histogram(trait_data$CPt, breaks=50) # --> normal with some low values
histogram(trait_data$troph_level, breaks=50) # --> normalish

#PCA_data <- PCA_data %>% 
 # mutate(logMBl=ifelse(is.na(MBl), NA, ifelse(MBl==0, 0, log(MBl))), 
     #    logBEl=ifelse(is.na(BEl), NA, ifelse(BEl==0, 0, log(BEl))), 
     #    logRMl=ifelse(is.na(RMl), NA, ifelse(MBl==0, 0, log(BEl)))) %>% 
 # dplyr::select(-MBl, -BEl, -RMl)

# 2.  Estimate number of principle components -----------------------------

# Use estim_ncPCA function and scale all variables
# Use gcv cross validation to impute NA
trait_PC <- estim_ncpPCA(PCA_data, ncp.min=5, ncp.max=7, scale=TRUE, 
                         method.cv="gcv", verbose=FALSE)

# Plot these according to mean squares estimate prediction
plot(trait_PC$criterion, xlab = "nb dim", ylab = "MSEP")

# Create imputed dataset using iterative PCA algorithm and cbind
impPCA_data_list <- imputePCA(PCA_data, ncp = trait_PC$ncp)
impPCA_data <- cbind.data.frame(sp_data, impPCA_data_list$completeObs)

# Check correlation between imputed variables
corr_mat <- as.matrix(cor(testimp_PCA[,2:11]))
corr_mat <- cor(matrix(rnorm(100), ncol = 10))
# Get all unique correlations
corr_vec <- as.vector(corr_mat) # convert correlation matrix to a vector
n_vars <- ncol(corr_mat)
corr_pairs <- combn(n_vars, 2) # get all pairs of indices
corr_df <- data.frame(Variable1 = rep(names(corr_mat), each = n_vars-1), 
                      Variable2 = names(corr_mat)[corr_pairs[2,]], 
                      Correlation = corr_vec[corr_pairs[1,] + (corr_pairs[2,]-1)*n_vars]) # combine pairs and correlations
corr_df <- corr_df[order(corr_df$Variable1, corr_df$Variable2),] # order by variable names
# order by variable names


# Print the resulting dataframe
corr_df


# 3.  Perform PCA ---------------------------------------------------------

# Log transform variables 
PCA_df <- impPCA_data %>% 
  mutate(logMBl=ifelse(is.na(MBl), NA, ifelse(MBl==0, 0, log(MBl))), 
         logBEl=ifelse(is.na(BEl), NA, ifelse(BEl==0, 0, log(BEl))), 
         logRMl=ifelse(is.na(RMl), NA, ifelse(MBl==0, 0, log(BEl)))) %>% 
  dplyr::select(-MBl, -BEl, -RMl)

# PCA using imputed data
PCA1 <- PCA(PCA_df, quali.sup=1, ncp=10, graph=T, axes=c(1,2))
summary(PCA1)

# See plots of variable loadings
plot(PCA1, choix="ind")
plot(PCA1, choix="var")

# Add PCA coordinates onto trait df
trait_and_PCA_data <- cbind(impPCA_data, PCA1$ind$coord[,1:2])

# 4.  Save data ---------------------------------------------

write.csv(impPCA_data, "./outputs/imputed_PCA_data.csv") # Note this is not scaled
write.csv(trait_and_PCA_data, "./inputs/data/trait_and_PCA_data.csv")

# 5.  Compute community weighted mean trait values ------------------------

comms <- full_results_df %>% 
  mutate(bins.before = as.numeric(bins)-as.numeric(bin.lag)) %>% 
  dplyr::select(sqID, bins, bins.before)

# Use this function to extract the abundance weighted community mean trait for all 11 traits
comms_logMBl_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "logMBl")
comms_logBEl_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "logBEl")
comms_logRMl_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "logRMl")
comms_VEp_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "VEp")
comms_REs_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "REs")
comms_OGp_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "OGp")
comms_BLs_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "BLs")
comms_PFv_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "PFv")
comms_PFs_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "PFs")
comms_CPt_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "CPt")
comms_troph_list <- find.commS.trait.data(comms, taxonomic_matrix_list, impPCA_data, "troph_level")

# 6.  Cleanup ------------------------------------------------------------

rm(full_results_df, taxonomic_matrix_list, comms, trait_data, survey_data,
   ts_data, PCA_data, sp_data, trait_PC, impPCA_data_list, impPCA_data,
   PCA1, PCA2, comms_logBEl_list, comms_logMBl_list, comms_logRMl_list,
   comms_VEp_list, comms_OGp_list, comms_REs_list, 






