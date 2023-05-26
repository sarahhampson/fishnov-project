# Title: Tax and Func Novelty Detection
# Author: Sarah Hampson


# 0.  Load in necessary dataframes ----------------------------------------

survey_data <- read.csv("./inputs/data/survey_data.csv")
ts_data <- read.csv("./inputs/data/ts_data.csv")
trait_data <- read.csv("./inputs/data/trait_data.csv")
#func_dissims_matrix <- as.matrix(read.csv("./inputs/data/func_dissims_matrix.csv", 
                                #row.names=1))
#colnames(func_dissims_matrix) <- rownames(func_dissims_matrix)

# 1.  Make relative abundance dataframes ----------------------------------

# Create a list of relative abundance matrices using matrix maker function
# Note these are split into quarter time series
taxonomic_matrix_list <- seasonal.matrix.maker(survey_data, ts_data)

# 2.  Taxonomic novelty detection -----------------------------------------

# Apply framework using taxonomiv novelty detection function 
# Note this removes the first 5 time bins, and also uses the old 
# identify.novel.gam function
tax_results_list <- tax.novelty.detection(taxonomic_matrix_list)

# Empty plots folder if wanted
unlink("./plots/Taxonomic Plots/*")

# 3.  Make func dissims matrix --------------------------------------------

# Need to make a multitrait dissimilarity matrix between all species using gowdis
# Only use morphological traits
func_dissims_matrix <- as.matrix(gowdis(trait_data[,3:12]), 
                                 dimnames = list(trait_data$Species))
dimnames(func_dissims_matrix) <- list(trait_data$Species, trait_data$Species)
attributes(func_dissims_matrix) # All traits should be numeric

# Change to numeric
class(func_dissims_matrix) <- "numeric"

# See first 5 rows and cols of diss matrix
func_dissims_matrix[1:5, 1:5]

# Save func dissims matrix as csv if needed later
filename_dissims <- "./inputs/data/func_dissims_matrix.csv"
write.csv(func_dissims_matrix, filename_dissims)

# 3.  Functional novelty detection ----------------------------------------

# Apply framework using functional novelty detection function
func_results_list <- func.novelty.detection(taxonomic_matrix_list, 
                                            func_dissims_matrix)
# Empty plots folder if wanted
unlink("./plots/Functional Plots/*")

# 4.  Clean up results ----------------------------------------------------

# Include cat.after column in tax and func results
tax_results_list <- lapply(tax_results_list[1:length(tax_results_list)], add.cat.after)
func_results_list <- lapply(func_results_list[1:length(func_results_list)], add.cat.after)

# Convert results lists into dataframes
tax_results_df <- do.call("rbind", tax_results_list)
func_results_df <- do.call("rbind", func_results_list)

# Change col names of func nov results columns to match tax nov results
func_results_df <- func_results_df %>% 
  rename(bins = timeIds, funcnovel = isNovelComm, bin.lag = timeLag)

# Change col names of tax nov results so novel is taxnovel
tax_results_df <- tax_results_df %>% 
  rename(taxnovel = novel)

# Add sqID index into each df
tax_results_df$sqID <- rep(names(tax_results_list), 
                           sapply(tax_results_list, nrow))
func_results_df$sqID <- rep(names(func_results_list), 
                            sapply(func_results_list, nrow))

# Add TS and quarter columns
tax_results_df[c('Site', 'Quarter')] <- 
  str_split_fixed(tax_results_df$site, " Q", 2)
func_results_df[c('Site', 'Quarter')] <- 
  str_split_fixed(func_results_df$sqID, " Q", 2)

# 5.  Create full nov detection dataframe ---------------------------------

# Create a full results df with other necessary TS information for analysis
full_results_df <- merge(tax_results_df, func_results_df,
                         by=c("bins", "bin.lag", "sqID", 
                              "Site", "Quarter"),
                         suffixes = c(".tax",".func"))
                         
# 6.  Merge with site data ------------------------------------------------

# Read in ts_data
ts_data <- read.csv("./inputs/data/ts_data.csv")

# Merge with ts data for each results df
tax_results_df <- merge(tax_results_df, ts_data, by.x="Site", by.y="TimeSeriesID")
func_results_df <- merge(func_results_df, ts_data, by.x="Site", by.y="TimeSeriesID")
full_results_df <- merge(full_results_df, ts_data, by.x="Site", by.y="TimeSeriesID")

# Conserve rownames
rownames(full_results_df) <- rownames(func_results_df)

# Filter tax results df to remove extra time series results
tax_results_df <- tax_results_df %>% 
  filter(tax_results_df$sqID %in% func_results_df$sqID)

# 7.  Save results dataframes ---------------------------------------------

write.csv(tax_results_df, "./detection results/tax_results_df.csv")
write.csv(func_results_df, "./detection results/func10_results_df.csv")
write.csv(full_results_df, "./detection results/full_results_df.csv")

# 8.  Cleanup -------------------------------------------------------------

# Note - do not remove tax_results_list or func_results_list from env

rm(taxonomic_matrix_list, func_dissims_matrix, trait_data, ts_data, 
   survey_data, tax_results_df, func_results_df, full_results_df)

