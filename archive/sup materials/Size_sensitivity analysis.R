# Sensitivity test removing time series with the 3 fish 
# that had standard length instead of MBl

# Species names
sp3 <- c("Hypophthalmus oremaculatus", "Percina nevisense", "Pterygoplichthys ambrosettii")

# Full survey data that I used
raw_survey_data <- read.csv("inputs/raw_data/1873_2_RivFishTIME_SurveyTable.csv")
time_series_data <-read.csv("inputs/raw_data/1873_2_RivFishTIME_TimeseriesTable.csv")

# Add sqID
taxsize_survey_data <- select.timeseries(raw_survey_data, time_series_data)

# Identify which sqID have the three species of SL
sp3_sqID <- taxsize_survey_data %>%
  filter(Species %in% sp3) %>% 
  summarise(sqID =unique(sqID))

# Remove these sqIDs from comms_size_df
comms_size_df_subsp3 <- comms_size_df %>% 
  filter(!sqID %in% sp3_sqID$sqID)

# Make site based size results
size_results_df_site <- comms_size_df_subsp3 %>% group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill=bincount/ts_length) %>% 
  ungroup() %>% group_by(Site=TimeSeriesID, Country, HydroBasin) %>% 
  dplyr::summarise(total_bincount = sum(bincount), mean_ts_fill = mean(ts_fill),
                   novel.TF = ifelse(sum(funcnovel=="TRUE")>0,1,0))

# Site level df for co-occurrence models
cooccur_results_df_site <- comms_size_df_subsp3 %>% 
  mutate(cooccur.TF = ifelse(taxnovel==TRUE & funcnovel==TRUE, 1, 0)) %>% 
  group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill=bincount/ts_length) %>% 
  ungroup() %>% group_by(Site=TimeSeriesID, Country, HydroBasin) %>% 
  dplyr::summarise(taxnovel=ifelse(sum(taxnovel=="TRUE")>0, TRUE, FALSE), 
                   funcnovel=ifelse(sum(funcnovel=="TRUE")>0, TRUE, FALSE),
                   cooccur=sum(cooccur.TF),
                   total_bincount = sum(bincount), 
                   mean_ts_fill = mean(ts_fill))

# Scale effects
size_results_df_site$mean_ts_fillS <- scale(as.numeric(size_results_df_site$mean_ts_fill))[,1]
size_results_df_site$bincountS <- scale(as.numeric(size_results_df_site$total_bincount))[,1]
cooccur_results_df_site$mean_ts_fillS <- scale(as.numeric(cooccur_results_df_site$mean_ts_fill))[,1]
cooccur_results_df_site$bincountS <- scale(as.numeric(cooccur_results_df_site$total_bincount))[,1]

# Make Models
size.glm.site.sens <- glmer(novel.TF ~ bincountS + mean_ts_fillS + (1|Country/HydroBasin), 
                       data = size_results_df_site,
                       family = binomial)
cooccur.glm.site.sens <- glmer(funcnovel ~ taxnovel + bincountS + mean_ts_fillS
                          + (1|Country/HydroBasin), 
                          data = cooccur_results_df_site,
                          family = binomial)

# See summaries
summary(size.glm.site.sens)
summary(cooccur.glm.site.sens)


# Create results df and compare to original
make.glm.df(size.glm.site.sens, "", "plot")
make.glm.df(size.glm.site, "", "plot")
make.glm.df(cooccur.glm.site.sens, "", "plot")
make.glm.df(cooccur.glm.site, "", "plot")


