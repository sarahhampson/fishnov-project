# Title: Tax and Size Analysis
# Author: Sarah Hampson


# 0.  Load necessary dataframes -------------------------------------------

# Load in the size+tax results and size data

size_results_df <- read.csv("~/Documents/UNI/2023/Honours/BIOL6502/Data/sarahfish/outputs/Full_resultsMARCH 2023-03-06.csv")
size_data <- read.csv("~/Documents/UNI/2023/Honours/BIOL6502/Data/sarahfish/outputs/minithesis.bodylengthdata.csv")

# Clean up results df
size_results_df <- size_results_df %>%
  rename(Site = TimeSeriesID) %>% 
  dplyr::select(Site, bins, bin.lag, sqID, Quarter,
                taxnovel, funcnovel, Country.x, 
                HydroBasin.x, BioRealm.x) %>% 
  rename(Country=Country.x,
         HydroBasin = HydroBasin.x,
         BioRealm = BioRealm.x)

# 1.  Add and scale fixed effects -----------------------------------------

# Add time series length, years passed and year position to results dfs
size_results_temp <- lapply(split(size_results_df, f=size_results_df$sqID),
                           function(x){
                             x$TS_Length <- nrow(x)
                             x$Yrs_passed <- x$bins - min(x$bins)
                             x$Yr_pos <- 1:nrow(x)
                             return(x)
                           })
size_results_df <- do.call("rbind", size_results_temp)

# Use scale function to be able to compare effect sizes
size_results_df$bin.lagS <- scale(size_results_df$bin.lag)[,1]
size_results_df$Yrs_passedS <- scale(size_results_df$Yrs_passed)[,1]
size_results_df$logYr_posS <- scale(log(size_results_df$Yr_pos))[,1]

# 2.  Create dataframes for site level models -----------------------------

# Create new df for number of sites
# This includes adding mean time series fill and total bin count 
tax_results_df_site <- size_results_df %>% group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill = bincount/ts_length) %>% 
  ungroup() %>% group_by(Site, Country, HydroBasin) %>% 
  dplyr::summarise(total_bincount = sum(bincount), mean_ts_fill = mean(ts_fill),
                   novel.TF = ifelse(sum(taxnovel=="TRUE")>0,1,0))
size_results_df_site <- size_results_df %>% group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill=bincount/ts_length) %>% 
  ungroup() %>% group_by(Site, Country, HydroBasin) %>% 
  dplyr::summarise(total_bincount = sum(bincount), mean_ts_fill = mean(ts_fill),
                   novel.TF = ifelse(sum(funcnovel=="TRUE")>0,1,0))

# Site level df for co-occurrence models
#cooccur_results_df_site <- full_results_df %>% 
#mutate(cooccur = ifelse(funcnovel=="TRUE" & taxnovel=="TRUE",1,0)) %>% 
#group_by(sqID) %>% 
#mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill=bincount/ts_length) %>% 
#ungroup() %>% group_by(Site, Country, HydroBasin) %>% 
#dplyr::summarise(total_bincount = sum(bincount), mean_ts_fill = mean(ts_fill),
# cooccur.TF =ifelse(sum(cooccur)>0, 1, 0), # There are no sites where co occur>1
#funcnov =ifelse(sum(funcnovel=="TRUE")>0, 1, 0))
cooccur_results_df_site <- size_results_df %>% 
  mutate(cooccur.TF = ifelse(taxnovel==TRUE & funcnovel==TRUE, 1, 0)) %>% 
  group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill=bincount/ts_length) %>% 
  ungroup() %>% group_by(Site, Country, HydroBasin) %>% 
  dplyr::summarise(taxnovel=ifelse(sum(taxnovel=="TRUE")>0, TRUE, FALSE), 
                   funcnovel=ifelse(sum(funcnovel=="TRUE")>0, TRUE, FALSE),
                   cooccur=sum(cooccur.TF),
                   total_bincount = sum(bincount), 
                   mean_ts_fill = mean(ts_fill))

# Scale effects
tax_results_df_site$mean_ts_fillS <- scale(as.numeric(tax_results_df_site$mean_ts_fill))[,1]
tax_results_df_site$bincountS <- scale(as.numeric(tax_results_df_site$total_bincount))[,1]
size_results_df_site$mean_ts_fillS <- scale(as.numeric(size_results_df_site$mean_ts_fill))[,1]
size_results_df_site$bincountS <- scale(as.numeric(size_results_df_site$total_bincount))[,1]
cooccur_results_df_site$mean_ts_fillS <- scale(as.numeric(cooccur_results_df_site$mean_ts_fill))[,1]
cooccur_results_df_site$bincountS <- scale(as.numeric(cooccur_results_df_site$total_bincount))[,1]

# 3.  Save model dataframes -----------------------------------------------

# These will be needed potentially for making the same models later for
# Final figures

# Time point level df
write.csv(size_results_df, "./archive/Size_Only_Stuff/size_glm_tp_df.csv")

# Site level df
write.csv(tax_results_df_site, "./archive/Size_Only_Stuff/tax_glm_site_df.csv")
write.csv(size_results_df_site, "./archive/Size_Only_Stuff/size_glm_site_df.csv")
write.csv(cooccur_results_df_site, "./archive/Size_Only_Stuff/cooccur_glm_site_df.csv")

# 4.  Taxonomic novelty occurrence models ----------------------------------

# Model of probability of tax novelty emergence across all time points (tp)
# Fixed effects are bin lag and ln bin position (scaled)
# Random effects are quarter time series nested in site and hydrobasin nested in country
tax.glm.tp <- glmer(taxnovel ~ 1 + bin.lagS + logYr_posS + 
                      (1|Site/Quarter) + (1|Country/HydroBasin), # Use this one
                    data = size_results_df,
                    family = binomial)

# Model of probability of tax novelty at least once for an individual site
# Fixed effects are total number of bins and mean time series fill (no. bins/ts length)
# Random effects are hydrobasin nested in country
tax.glm.site <- glmer(novel.TF ~ bincountS + mean_ts_fillS + (1|Country/HydroBasin), 
                      data = tax_results_df_site,
                      family = binomial)

# See summaries
summary(tax.glm.tp)
summary(tax.glm.site)

# Create results df
make.glm.df(tax.glm.tp, "Tax", "plot")
make.glm.df(tax.glm.site, "Tax", "plot")

# Ranef gives us variance predicted by random effect intercepts
# shows sites don't explain any more on their own than when compiling together
ranef(tax.glm.tp) 
ranef(tax.glm.site)

# Look at performance of the models
install.packages("performance")
library(performance)
performance(tax.glm.tp)
performance(tax.glm.site)

# 5.  Size novelty occurrence models ---------------------------------

# Model of probability of size novelty emergence across all time points (tp)
# Fixed effects are bin lag and ln bin position (scaled)
# Random effects are quarter time series nested in site and hydrobasin nested in country
size.glm.tp <- glmer(funcnovel ~ 1 + bin.lagS + logYr_posS + 
                       (1|Site/Quarter) + (1|Country/HydroBasin), # Use this one
                     data = size_results_df,
                     family = binomial)

# Model of probability of func novelty at least once for an individual site
# Fixed effects are total number of bins and mean time series fill (no. bins/ts length)
# Random effects are hydrobasin nested in country
size.glm.site <- glmer(novel.TF ~ bincountS + mean_ts_fillS + (1|Country/HydroBasin), 
                       data = size_results_df_site,
                       family = binomial)

# See summaries
summary(size.glm.tp)
summary(size.glm.site)

# Create results df
make.glm.df(size.glm.tp, "Size", "plot")
make.glm.df(size.glm.site, "Size", "plot")

# 6.  Novelty co-occurrence models ----------------------------------------

# Time point model of novelty cooccurrence using full results df
cooccur.glm.tp <- glmer(funcnovel ~ taxnovel + bin.lagS + logYr_posS
                        + (1|Site/Quarter) + (1|Country/HydroBasin),
                        data = size_results_df,
                        family=binomial)

# Site level model of novelty cooccurrence using full results df
cooccur.glm.site <- glmer(funcnovel ~ taxnovel + bincountS + mean_ts_fillS
                          + (1|Country/HydroBasin), 
                          data = cooccur_results_df_site,
                          family = binomial)

# See summaries
summary(cooccur.glm.tp)
summary(cooccur.glm.site)

# Create results df
make.glm.df(cooccur.glm.tp, "cooccur", "plot")
make.glm.df(cooccur.glm.site, "cooccur", "plot")

# 7.  Test co-occurrence models -------------------------------------------

# Create another model to test co-occurrence models
# Models will be functional novelty as response with covariates

# Create subset timepoint data for tax novelty == TRUE
sub_tp_df <- subset(size_results_df, size_results_df$taxnovel==TRUE)

# Create subset site data for tax novelty == TRUE
sub_site_df <- subset(cooccur_results_df_site, cooccur_results_df_site$taxnovel==TRUE)

# Create models for func novelty using subset tax novel data
cooccur.model.tp.test <- glmer(funcnovel ~ 1 + bin.lagS + logYr_posS +
                                 (1|Site/Quarter) + (1|Country/HydroBasin),
                               data = sub_tp_df, 
                               family = binomial)
cooccur.model.site.test <- glmer(funcnovel ~ 1 + bincountS + mean_ts_fillS +
                                   (1|Country/HydroBasin),
                                 data = sub_site_df, 
                                 family = binomial)

# See summaries
summary(cooccur.model.tp.test)
summary(cooccur.model.site.test)

# Create results df
make.glm.df(cooccur.model.tp.test, "cooccurtest", "plot")
make.glm.df(cooccur.model.site.test, "cooccurtest", "plot")



# 8.  Import time size change data ----------------------------------------

size_change_data <- read.csv("archive/Size_Only_Stuff/sizechange_data.csv")

# 9.  Plot size change data ----------------------------------------------------

#### Plot time series ####

ggplot(data=survey_data, aes(x=Year)) + 
  stat_bin(aes(y = ..density..), binwidth = 1, fill = "#50b6b9", alpha = 0.8) +
  geom_density(alpha = 0.1, color = "black") + theme_classic()

#### Plot size data

# Assuming your data is in a data.frame called `size_data` with columns `species` and `length`
ggplot(size_data, aes(x = Length)) +
  geom_histogram(binwidth = 10, color = "black", fill = "lightblue") +
  labs(x = "Maximum body length (cm)", y = "Number of Species", title = "Distribution of Fish Lengths") +
  theme_classic()

