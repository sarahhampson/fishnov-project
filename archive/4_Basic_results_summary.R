# Title: Basic results diagnostics/summaries
# Author: Sarah Hampson

# 0.  Load necessary dataframes -------------------------------------------

# Load full results csv
full_results_df <- read.csv("./detection results/full_results_df.csv")

# 1.  How many times did novelty occur? -----------------------------------

# Count of times novelty occurred
sum(full_results_df$taxnovel=="TRUE")
sum(full_results_df$funcnovel=="TRUE")

# Summary table for tax novelty
full_results_df %>%
  filter(taxnovel) %>%
  group_by(Site) %>%
  summarize(taxnovel_count = n()) %>%
  group_by(taxnovel_count) %>%
  summarize(site_count = n())
# Summary table for func novelty
full_results_df %>%
  filter(funcnovel) %>%
  group_by(Site) %>%
  summarize(funcnovel_count = n()) %>%
  group_by(funcnovel_count) %>%
  summarize(site_count = n())

# 2.  What country/biorealm did they occur in? ----------------------------

# Countries
# Summary table for tax novelty
full_results_df %>%
  filter(taxnovel) %>%
  group_by(Country) %>%
  summarize(taxnovel_count = n())
# Summary table for func novelty
full_results_df %>%
  filter(funcnovel) %>%
  group_by(Country) %>%
  summarize(funcnovel_count = n())

# BioRealm
# Summary table for tax novelty
full_results_df %>%
  filter(taxnovel) %>%
  group_by(BioRealm) %>%
  summarize(taxnovel_count = n())
# Summary table for func novelty
full_results_df %>%
  filter(funcnovel) %>%
  group_by(BioRealm) %>%
  summarize(funcnovel_count = n())

# 3. How often did novelty co-occur? --------------------------------------

# At level of community states/time point
full_results_df %>% dplyr::summarise(sum(taxnovel==TRUE & funcnovel==TRUE)) #40

# At level of community as a whole
nov_comms <- full_results_df %>% 
  mutate(cooccur.TF = ifelse(taxnovel==TRUE & funcnovel==TRUE, 1, 0)) %>% 
  group_by(Site, Country, HydroBasin) %>% 
  dplyr::summarise(taxnovel=ifelse(sum(taxnovel=="TRUE")>0, TRUE, FALSE), 
                   funcnovel=ifelse(sum(funcnovel=="TRUE")>0, TRUE, FALSE),
                   cooccur=sum(cooccur.TF))
sum(nov_comms$taxnovel == TRUE & nov_comms$funcnovel == TRUE) # 53


# 4. Cleanup --------------------------------------------------------------

rm(full_results_df, nov_comms)
