# Title: GLM models for novelty occurrence
# Author: Sarah Hampson


# 0.  Load necessary dataframes -------------------------------------------

# Read results csvs
tax_results_df <- read.csv("./detection results/tax_results_df.csv")
func_results_df <- read.csv("./detection results/func10_results_df.csv")
full_results_df <- read.csv("./detection results/full_results_df.csv")

# 1.  Add and scale fixed effects -----------------------------------------

# Add time series length, years passed and year position to results dfs
tax_results_temp <- lapply(split(tax_results_df, f=tax_results_df$sqID),
                           function(x){
                             x$TS_Length <- nrow(x)
                             x$Yrs_passed <- x$bins - min(x$bins)
                             x$Yr_pos <- 1:nrow(x)
                             return(x)
                           })
tax_results_df <- do.call("rbind", tax_results_temp)
func_results_temp <- lapply(split(func_results_df, f=func_results_df$sqID),
                            function(x){
                              x$TS_Length <- nrow(x)
                              x$Yrs_passed <- x$bins - min(x$bins)
                              x$Yr_pos <- 1:nrow(x)
                              return(x)
                            })
func_results_df <- do.call("rbind", func_results_temp)
full_results_temp <- lapply(split(multitrait_results_df, f=multitrait_results_df$sqID),
                            function(x){
                              x$TS_Length <- nrow(x)
                              x$Yrs_passed <- x$bins - min(x$bins)
                              x$Yr_pos <- 1:nrow(x)
                              return(x)
                            })
full_results_df <- do.call("rbind", full_results_temp)

# Use scale function to be able to compare effect sizes
tax_results_df$bin.lagS <- scale(tax_results_df$bin.lag)[,1]
tax_results_df$Yrs_passedS <- scale(tax_results_df$Yrs_passed)[,1]
tax_results_df$logYr_posS <- scale(log(tax_results_df$Yr_pos))[,1] # Using log for bin position
func_results_df$bin.lagS <- scale(func_results_df$bin.lag)[,1]
func_results_df$Yrs_passedS <- scale(func_results_df$Yrs_passed)[,1]
func_results_df$logYr_posS <- scale(log(func_results_df$Yr_pos))[,1]
full_results_df$bin.lagS <- scale(full_results_df$bin.lag)[,1]
full_results_df$Yrs_passedS <- scale(full_results_df$Yrs_passed)[,1]
full_results_df$logYr_posS <- scale(log(full_results_df$Yr_pos))[,1]

# 2.  Create dataframes for site level models -----------------------------

# Create new df for number of sites
# This includes adding mean time series fill and total bin count 
tax_results_df_site <- comms_size_df %>% group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill = bincount/ts_length) %>% 
  ungroup() %>% group_by(Site=TimeSeriesID, Country, HydroBasin) %>% 
  dplyr::summarise(total_bincount = sum(bincount), mean_ts_fill = mean(ts_fill),
                   novel.TF = ifelse(sum(taxnovel=="TRUE")>0,1,0))
func_results_df_site <- full_results_df %>% group_by(sqID) %>% 
  mutate(bincount = n(), ts_length = max(bins) - min(bins)+1, ts_fill=bincount/ts_length) %>% 
  ungroup() %>% group_by(Site, Country, HydroBasin) %>% 
  dplyr::summarise(total_bincount = sum(bincount), mean_ts_fill = mean(ts_fill),
                   novel.TF = ifelse(sum(funcnovel=="TRUE")>0,1,0))

# Site level df for co-occurrence models
cooccur_results_df_site <- comms_size_df %>% 
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
tax_results_df_site$mean_ts_fillS <- scale(as.numeric(tax_results_df_site$mean_ts_fill))[,1]
tax_results_df_site$bincountS <- scale(as.numeric(tax_results_df_site$total_bincount))[,1]
func_results_df_site$mean_ts_fillS <- scale(as.numeric(func_results_df_site$mean_ts_fill))[,1]
func_results_df_site$bincountS <- scale(as.numeric(func_results_df_site$total_bincount))[,1]
cooccur_results_df_site$mean_ts_fillS <- scale(as.numeric(cooccur_results_df_site$mean_ts_fill))[,1]
cooccur_results_df_site$bincountS <- scale(as.numeric(cooccur_results_df_site$total_bincount))[,1]

# 3.  Save model dataframes -----------------------------------------------

# These will be needed potentially for making the same models later for
# Final figures

# Time point level df
write.csv(tax_results_df, "./inputs/model_data/tax_glm_tp_df.csv")
write.csv(func_results_df, "./inputs/model_data/func_glm_tp_df.csv")
write.csv(full_results_df, "./inputs/model_data/cooccur_glm_tp_df.csv")

# Site level df
write.csv(tax_results_df_site, "./inputs/model_data/tax_glm_site_df.csv")
write.csv(func_results_df_site, "./inputs/model_data/func_glm_site_df.csv")
write.csv(cooccur_results_df_site, "./inputs/model_data/cooccur_glm_site_df.csv")

# 4.  Taxonomic novelty occurrence models ----------------------------------

# Model of probability of tax novelty emergence across all time points (tp)
# Fixed effects are bin lag and ln bin position (scaled)
# Random effects are quarter time series nested in site and hydrobasin nested in country
tax.glm.tp <- glmer(taxnovel ~ 1 + bin.lagS + logYr_posS + 
                       (1|Site/Quarter) + (1|Country/HydroBasin), # Use this one
                     data = tax_results_df,
                     family = binomial)

# Model of probability of tax novelty at least once for an individual site
# Fixed effects are total number of bins and mean time series fill (no. bins/ts length)
# Random effects are hydrobasin nested in country
tax_results_df_site <- read.csv("./archive/Size_Only_Stuff/tax_glm_site_df.csv")
size_results_df_site <- read.csv("./archive/Size_Only_Stuff/size_glm_site_df.csv")
tax.glm.site <- glmer(novel.TF ~ bincountS + mean_ts_fillS + (1|Country/HydroBasin), 
                          data = tax_results_df_site,
                          family = binomial)
size.glm.site <- glmer(novel.TF ~ bincountS + mean_ts_fillS + (1|Country/HydroBasin), 
                      data = size_results_df_site,
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

# 5.  Multi-trait novelty occurrence models ---------------------------------

# Model of probability of func novelty emergence across all time points (tp)
# Fixed effects are bin lag and ln bin position (scaled)
# Random effects are quarter time series nested in site and hydrobasin nested in country
func.glm.tp <- glmer(funcnovel ~ 1 + bin.lagS + logYr_posS + 
                      (1|Site/Quarter) + (1|Country/HydroBasin), # Use this one
                    data = func_results_df,
                    family = binomial)

# Model of probability of func novelty at least once for an individual site
# Fixed effects are total number of bins and mean time series fill (no. bins/ts length)
# Random effects are hydrobasin nested in country
func.glm.site <- glmer(novel.TF ~ bincountS + mean_ts_fillS + (1|Country/HydroBasin), 
                      data = func_results_df_site,
                      family = binomial)

# See summaries
summary(func.glm.tp)
summary(func.glm.site)

# Create results df
make.glm.df(func.glm.tp, "Func", "plot")
make.glm.df(func.glm.site, "Func", "plot")

# Ranef gives us variance predicted by random effect intercepts
# shows sites don't explain any more on their own than when compiling together
ranef(func.glm.tp) 
ranef(func.glm.site)

# Look at performance of the models
performance(func.glm.tp)
performance(func.glm.site)

# 6.  Novelty co-occurrence models ----------------------------------------

# Time point model of novelty cooccurrence using full results df
cooccur.glm.tp <- glmer(funcnovel ~ taxnovel + bin.lagS + logYr_posS
                        + (1|Site/Quarter) + (1|Country/HydroBasin),
                        data = full_results_df,
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

# Ranef gives us variance predicted by random effect intercepts
# shows sites don't explain any more on their own than when compiling together
ranef(cooccur.glm.tp) 
ranef(cooccur.glm.site)

# Look at performance of the models
performance(cooccur.glm.tp)
performance(cooccur.glm.site)

# 7.  Test co-occurrence models -------------------------------------------

# Create another model to test co-occurrence models
# Models will be functional novelty as response with covariates

# Create subset timepoint data for tax novelty == TRUE
sub_tp_df <- subset(full_results_df, full_results_df$taxnovel==TRUE)

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

# 8.  Make/Save GLM tables -----------------------------------------------------

library(sjPlot)

tab_model


# Time point level models
tax_model_df_tp <- make.glm.df(tax.glm.tp, "Tax", "plot")
func_model_df_tp <- make.glm.df(func.glm.tp, "Func", "plot")
cooccur_model_df_tp <- make.glm.df(cooccur.glm.tp, "Co-occur", "plot")
write.csv(tax_model_df_tp, "./outputs/GLM/tax_glm_tp.csv")
write.csv(func_model_df_tp, "./outputs/GLM/func_glm_tp.csv")
write.csv(cooccur_model_df_tp, "./outputs/GLM/cooccur_glm_tp.csv")

# Site level models
tax_model_df_site <- make.glm.df(tax.glm.site, "Tax", "plot")
func_model_df_site <- make.glm.df(func.glm.site, "Func", "plot")
cooccur_model_df_site <- make.glm.df(cooccur.glm.site, "Co-occur", "plot")
write.csv(tax_model_df_site, "./outputs/GLM/tax_glm_site.csv")
write.csv(func_model_df_site, "./outputs/GLM/func_glm_site.csv")
write.csv(cooccur_model_df_site, "./outputs/GLM/cooccur_glm_site.csv")

# 9.  Cleanup ------------------------------------------------------------

rm(tax_results_df, func_results_df, full_results_df, tax_results_temp,
   func_results_temp, tax_results_df_site, func_results_df_site, tax.glm.tp,
   tax.glm.site, func.glm.tp, func.glm.site, cooccur_results_df_site,
   cooccur.glm.site, cooccur.glm.tp, full_results_temp, tax_model_df_tp,
   func_model_df_tp, cooccur_model_df_tp, tax_model_df_site, func_model_df_site,
   cooccur_model_df_site, cooccur.model.site.test, cooccur.model.tp.test,
   sub_site_df, sub_tp_df)
