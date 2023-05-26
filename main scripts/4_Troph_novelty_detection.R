# Title: Trophic novelty detection


# 0.  Load necessary dataframes -------------------------------------------

# RivFishTIME data with standardised species
time_series_data <-read.csv("inputs/raw_data/1873_2_RivFishTIME_TimeseriesTable.csv")
raw_survey_data <- read.csv("inputs/raw_data/1873_2_RivFishTIME_SurveyTable.csv")

# 1.  Select time series --------------------------------------------------

# Select time series data using function
survey_data <- select.timeseries(survey_data, time_series_data)

# Standardise species names
survey_data <- standardise.species.names(survey_data)

# Filter out survey data again with less than 10 years of data 
# (needed again after standardising names)
survey_data <- survey_data %>% 
  group_by(sqID) %>% 
  filter(n_distinct(Year)>=10)

# Modify time series to only include time series from new ts_table
ts_data <- time_series_data[time_series_data$TimeSeriesID %in%
                              survey_data$TimeSeriesID,]

# Look at number of sites, time series, and species
n_distinct(survey_data$TimeSeriesID)
n_distinct(survey_data$sqID)
n_distinct(survey_data$Species)

# 2.  Make relative abundance matrices ------------------------------------

# Create a list of relative abundance matrices using matrix maker function
# Note these are split into quarter time series
taxonomic_matrix_list <- seasonal.matrix.maker(survey_data, ts_data)

# 3.  Find trophic data ---------------------------------------------------

species_list <- as.data.frame(unique(survey_data$Species))
colnames(species_list) <- "Species"

# Retrieve trophic level data from rfishbase
trophic_data_FB <- ecology(species_list$Species, 
                           fields = c("Species", "FoodTroph", "DietTroph"))

# Make list of species missing from troph data
missing_troph_sp <- species_list %>% 
  filter(!Species %in% trophic_data_FB$Species)

# Show more missing species from above troph list which have NAs and filter from df
missing_troph_sp2 <- trophic_data_FB %>% filter(is.na(FoodTroph), is.na(DietTroph)) %>% 
  dplyr::select(Species)
trophic_data_FB <- trophic_data_FB %>% filter(!Species %in% missing_troph_sp2$Species) 

# Make sure no double ups in troph data from FB by taking mean diet and food troph
trophic_data_FB <- trophic_data_FB %>% 
  group_by(Species) %>% 
  dplyr::summarise(FoodTroph = mean(FoodTroph, na.rm=T),
                   DietTroph = mean(DietTroph, na.rm=T))

# There are 227+48 so 275 species missing, which I found information on via online fishbase

missing_troph_sp <- rbind(missing_troph_sp, missing_troph_sp2)
write.csv(missing_troph_sp, "./outputs/missing_troph_species.csv")

# Read in csv of model based fishbase species without specific trophic information
troph_Est_data <- read_csv('./inputs/raw_data/trophic_level_nov_FBModels.csv', col_names = T)
troph_Est_data$Species <- gsub("_", " ", troph_Est_data$Species)
troph_Est_data <- troph_Est_data %>% 
  rename(troph_level = Trophic_level) %>% 
  mutate(troph_cat = "Est") %>% 
  dplyr::select(-Cat)

# Filter troph diet data first
troph_data_diet <- trophic_data_FB %>% 
  filter(!is.na(DietTroph)) %>% 
  rename(troph_level = DietTroph) %>% 
  dplyr::select(-FoodTroph) %>% 
  mutate(troph_cat = "Diet")

# Filter troph food data
troph_data_food <- trophic_data_FB %>% 
  filter(!Species %in% troph_data_diet$Species) %>% 
  rename(troph_level = FoodTroph) %>% 
  dplyr::select(-DietTroph) %>% 
  mutate(troph_cat = ifelse(troph_level>0, "Food", NA))

# Combine the data into one trophic df
troph_data <- rbind(troph_Est_data, troph_data_diet, troph_data_food)
write.csv(troph_data, "./inputs/data/trophi_data.csv")

rm(missing_troph_sp, missing_troph_sp2, troph_data_diet, troph_data_food,
   troph_data_FB, troph_Est_data, troph_data_FB)

# 4.  Subset time series and relative abundance matrices --------------------------------------------------

# Filter out time series and survey data which have missing troph species
missing_troph_sp <- troph_data %>% filter(is.na(troph_level)) %>% 
  dplyr::select(Species)
missing_troph_sqID <- survey_data %>% filter(Species %in% missing_troph_sp$Species) %>% 
  dplyr::summarise(sqID =unique(sqID))
troph_survey_data <- survey_data %>% filter(!sqID %in% missing_troph_sqID$sqID)

# Subset time series data
troph_ts_data <- ts_data %>% filter(TimeSeriesID %in% troph_survey_data$TimeSeriesID)

# Subset relative abundance matrices to exlude those with missing troph sp
subset_matrix_list <- taxonomic_matrix_list[!names(taxonomic_matrix_list) %in% missing_troph_sqID$sqID]

# 5.  Species trophic dissimilarity matrix --------------------------------

# Convert size dataframe to dissimilarity matrix
troph_diss_mat <- as.matrix(dist(troph_data$troph_level))
dimnames(troph_diss_mat) <- list(troph_data$Species, troph_data$Species)

# Change from character matrix to numeric matrix
class(troph_diss_mat) <- "numeric"

# Save troph dissims matrix as csv if needed later
filename_dissims <- "./inputs/data/troph_dissims_matrix.csv"
write.csv(troph_diss_mat, filename_dissims)








# 6.  Troph novelty detection ----------------------------------------

# Apply framework using functional novelty detection function
toph_results_list <- func.novelty.detection(subset_matrix_list, 
                                            troph_diss_mat)
# Empty plots folder if wanted
unlink("./plots/Functional Plots/*")

# 7.  Clean up results ----------------------------------------------------

# Convert results lists into dataframes
troph_results_df <- do.call("rbind", troph_results_list)

# Change col names of troph nov results columns to match tax nov results
troph_results_df <- troph_results_df %>% 
  rename(bins = timeIds, trophnovel = isNovelComm, bin.lag = timeLag)

# Add sqID index into each df
troph_results_df$sqID <- rep(names(troph_results_list), 
                            sapply(troph_results_list, nrow))

# Add TS and quarter columns
troph_results_df[c('Site', 'Quarter')] <- 
  str_split_fixed(troph_results_df$sqID, " Q", 2)

# Merge with ts data
troph_results_df <- merge(troph_results_df, ts_data, by.x="Site", by.y="TimeSeriesID")
rownames(troph_results_df) <- rownames(troph_results_df)

# 8.  Save results dataframe ---------------------------------------------

write.csv(troph_results_df, "./detection results/troph_results_df.csv")




# 9.  How many times did novelty occur?  ------------------------

# Count of times novelty occurred
sum(troph_results_df$trophnovel=="TRUE")

# Summary table 
troph_results_df %>%
  filter(trophnovel) %>%
  group_by(Site) %>%
  summarize(trophnovel_count = n()) %>%
  group_by(trophnovel_count) %>%
  summarize(site_count = n())

# 10. What country/biorealm did they occur in? ----------------------------

# Countries
troph_results_df %>%
  filter(trophnovel) %>%
  group_by(Country) %>%
  summarize(trophnovel_count = n())

# BioRealm
troph_results_df %>%
  filter(trophnovel) %>%
  group_by(BioRealm) %>%
  summarize(trophnovel_count = n())

# 11. Add and scale fixed effects -----------------------------------------

# Add time series length, years passed and year position 
troph_results_df$bins <- as.numeric(troph_results_df$bins)
troph_results_temp <- lapply(split(troph_results_df, f=troph_results_df$sqID),
                            function(x){
                              x$TS_Length <- nrow(x)
                              x$Yrs_passed <- x$bins - min(x$bins)
                              x$Yr_pos <- 1:nrow(x)
                              return(x)
                            })
troph_results_df <- do.call("rbind", troph_results_temp)

# Use scale function to be able to compare effect sizes
troph_results_df$bin.lagS <- scale(troph_results_df$bin.lag)[,1]
troph_results_df$Yrs_passedS <- scale(troph_results_df$Yrs_passed)[,1]
troph_results_df$logYr_posS <- scale(log(troph_results_df$Yr_pos))[,1]

# 12. Create dataframes for site level models -----------------------------

# Create new df for number of sites
# This includes adding mean time series fill and total bin count 
troph_results_df_site <- comms_troph_df %>% group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill = bincount/ts_length) %>% 
  ungroup() %>% group_by(Site, Country, HydroBasin) %>% 
  dplyr::summarise(total_bincount = sum(bincount), mean_ts_fill = mean(ts_fill),
                   novel.TF = ifelse(sum(trophnovel=="TRUE")>0,1,0))

# Site level df for co-occurrence models
tax_results_df <- read.csv("./archive/size_results_df.csv")

# Change classes of some variables
tax_results_df$Site <- tax_results_df$TimeSeriesID
troph_results_df$bins <- as.integer(troph_results_df$bins)
tax_results_df$HydroBasin <- as.character(tax_results_df$HydroBasin)
troph_results_df$HydroBasin <- as.character(troph_results_df$HydroBasin)
tax_results_df$bin.lag <- as.numeric(tax_results_df$bin.lag)

# Merge df
full_results_df <- tax_results_df %>% 
  dplyr::select(Site, bins, bin.lag, sqID, Quarter, HydroBasin, Country, taxnovel) %>% 
  filter(sqID %in% troph_results_df$sqID)
full_results_df <- merge(full_results_df, troph_results_df, by=c("Site", "Quarter", "sqID", 
                                                              "bins", "bin.lag", "Country",
                                                              "HydroBasin"))
# Make cooccurence df
cooccur_results_df_site <- comms_troph_df %>% 
  mutate(cooccur.TF = ifelse(taxnovel==TRUE & trophnovel==TRUE, 1, 0)) %>% 
  group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill=bincount/ts_length) %>% 
  ungroup() %>% group_by(Site, Country, HydroBasin) %>% 
  dplyr::summarise(taxnovel=ifelse(sum(taxnovel=="TRUE")>0, TRUE, FALSE), 
                   trophnovel=ifelse(sum(trophnovel=="TRUE")>0, TRUE, FALSE),
                   cooccur=sum(cooccur.TF),
                   total_bincount = sum(bincount), 
                   mean_ts_fill = mean(ts_fill))

# Scale effects
troph_results_df_site$mean_ts_fillS <- scale(as.numeric(troph_results_df_site$mean_ts_fill))[,1]
troph_results_df_site$bincountS <- scale(as.numeric(troph_results_df_site$total_bincount))[,1]
cooccur_results_df_site$mean_ts_fillS <- scale(as.numeric(cooccur_results_df_site$mean_ts_fill))[,1]
cooccur_results_df_site$bincountS <- scale(as.numeric(cooccur_results_df_site$total_bincount))[,1]

# 13. Save model dataframes -----------------------------------------------

# These will be needed potentially for making the same models later for
# Final figures

# Time point level df
write.csv(troph_results_df, "./inputs/model_data/troph_glm_tp_df.csv")
write.csv(full_results_df, "./inputs/model_data/trophcooccur_glm_tp_df.csv")

# Site level df
write.csv(troph_results_df_site, "./inputs/model_data/troph_glm_site_df.csv")
write.csv(cooccur_results_df_site, "./inputs/model_data/troph_cooccur_glm_site_df.csv")

# 14. Troph novelty models ----------------------------------

# Model of probability of tax novelty at least once for an individual site
# Fixed effects are total number of bins and mean time series fill (no. bins/ts length)
# Random effects are hydrobasin nested in country
troph.glm.site <- glmer(novel.TF ~ bincountS + mean_ts_fillS + (1|Country/HydroBasin), 
                      data = troph_results_df_site,
                      family = binomial)

# See summaries
summary(troph.glm.site)

# Create results df
make.glm.df(troph.glm.site, "Troph", "plot")

# Cooccurrence stuff
# At level of community states/time point
full_results_df %>% dplyr::summarise(sum(taxnovel==TRUE & trophnovel==TRUE)) #39

# At level of community as a whole
nov_comms <- full_results_df %>% 
  mutate(cooccur.TF = ifelse(taxnovel==TRUE & trophnovel==TRUE, 1, 0)) %>% 
  group_by(Site, Country, HydroBasin) %>% 
  dplyr::summarise(taxnovel=ifelse(sum(taxnovel=="TRUE")>0, TRUE, FALSE), 
                   trophnovel=ifelse(sum(trophnovel=="TRUE")>0, TRUE, FALSE),
                   cooccur=sum(cooccur.TF))
sum(nov_comms$taxnovel == TRUE & nov_comms$trophnovel == TRUE) 
rm(nov_comms)

# Site level model of novelty cooccurrence using full results df
cooccur_results_df_site$HydroBasin <- as.numeric(cooccur_results_df_site$HydroBasin)
cooccur.glm.site <- glmer(trophnovel ~ taxnovel + bincountS + mean_ts_fillS
                          + (1|Country/HydroBasin), 
                          data = cooccur_results_df_site,
                          family = binomial)

# See summaries
summary(cooccur.glm.site)

# Create results df
make.glm.df(cooccur.glm.site, "cooccur", "plot")

# 15. Test co-occurrence models -------------------------------------------

# Create subset site data for tax novelty == TRUE
sub_site_df <- subset(cooccur_results_df_site, cooccur_results_df_site$taxnovel==TRUE)

# Model
cooccur.model.site.test <- glmer(trophnovel ~ 1 + bincountS + mean_ts_fillS +
                                   (1|Country/HydroBasin),
                                 data = sub_site_df, 
                                 family = binomial)

summary(cooccur.model.site.test)

# Create results df
make.glm.df(cooccur.model.site.test, "cooccurtest", "plot")





